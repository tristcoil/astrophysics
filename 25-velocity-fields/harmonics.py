"""Non-circular motions: harmonic decomposition of the velocity field.

Pure circular rotation puts all the velocity signal around a ring into one
term, cos(theta). Anything else -- gas streaming along spiral arms, flowing
along a bar, a disc that is oval rather than round -- adds sines and higher
harmonics. Following Schoenmakers (1997) and Trachternach et al. (2008), the
velocities along each ring (geometry fixed from the tilted-ring fit, no
down-weighting of the minor axis) are expanded as

    v = c0 + sum_{m=1..3} [ c_m cos(m theta) + s_m sin(m theta) ],

and the non-circular amplitude is A_r = sqrt(s1^2 + c2^2 + s2^2 + c3^2 + s3^2),
divided by sin(i) to express it in the plane of the disc. Trachternach et al.
summarize a galaxy by the median of A_r over its rings, and over its inner
kiloparsec, where the cusp-core question is decided.

A global elongation of the potential leaves a characteristic signature,
eps_pot sin(2 phi_2) = (s3 - s1) (1 + 2q^2 + 5q^4) / (c1 (1 - q^4)), q = cos i
(Schoenmakers et al. 1997).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from config import HARMONIC_ORDER
from tiltedring import TiltedRingResult, _select, disc_coordinates, fit_ring


@dataclass
class Harmonics:
    radius: np.ndarray  # arcsec
    coeff: np.ndarray  # (n_rings, 2N+1): c0, c1, s1, c2, s2, ... line of sight, km/s
    err: np.ndarray
    sin_i: np.ndarray
    incl: np.ndarray
    npix: np.ndarray
    residual_median_abs: np.ndarray  # median |v - harmonic model| per ring

    def term(self, kind: str, m: int, deproject=True) -> np.ndarray:
        j = 2 * m - 1 if kind == "c" else 2 * m
        if m == 0:
            j = 0
        out = self.coeff[:, j]
        return out / self.sin_i if deproject and m > 0 else out

    @property
    def a_r(self) -> np.ndarray:
        """Quadratic sum of the non-circular terms, deprojected (km/s)."""
        terms = [self.term("s", 1)]
        for m in range(2, HARMONIC_ORDER + 1):
            terms += [self.term("c", m), self.term("s", m)]
        return np.sqrt(np.sum(np.square(terms), axis=0))

    def a_m(self, m: int) -> np.ndarray:
        if m == 1:
            return np.abs(self.term("s", 1))
        return np.hypot(self.term("c", m), self.term("s", m))

    @property
    def eps_pot(self) -> np.ndarray:
        q = np.cos(np.radians(self.incl))
        c1 = self.term("c", 1, deproject=False)
        s1, s3 = self.term("s", 1, deproject=False), self.term("s", 3, deproject=False)
        return (s3 - s1) * (1 + 2 * q ** 2 + 5 * q ** 4) / (c1 * (1 - q ** 4))


def decompose(x, y, v, tr: TiltedRingResult, order: int = HARMONIC_ORDER, min_pixels: int = 30,
              geometry: str = "adopted") -> Harmonics:
    """Harmonic expansion ring by ring.

    ``geometry='adopted'`` expands around the smoothed rings of the final
    rotation curve; ``'free'`` around each ring's own stage-2 position angle
    and inclination (fitted with |cos theta| weights); ``'reswri'`` refits
    each ring with no weighting first, as GIPSY's RESWRI does -- the setting
    of Trachternach et al. (2008). An inclination that does not suit all
    azimuths equally leaves a cos(3 theta) term (Schoenmakers 1999), which
    the unweighted refit removes.
    """
    n = len(tr.radius)
    pa_use, incl_use = tr.pa.copy(), tr.incl.copy()
    if geometry == "free":
        ok = np.isfinite(tr.free_pa) & np.isfinite(tr.free_incl)
        pa_use[ok], incl_use[ok] = tr.free_pa[ok], tr.free_incl[ok]
    elif geometry == "reswri":
        # refit each ring's PA, inclination and V_rot with no minor-axis
        # down-weighting (centre and V_sys fixed), as RESWRI does
        for k in range(n):
            ring = fit_ring(x, y, v, tr.radius[k], tr.width, tr.ring_params(k), ("pa", "incl", "vrot"),
                            min_pixels=min_pixels, weight_power=0)
            if ring.ok and 6 < ring.p["incl"] < 87:
                pa_use[k] = tr.pa[k] + (ring.p["pa"] - tr.pa[k] + 180) % 360 - 180
                incl_use[k] = ring.p["incl"]
    ncoef = 2 * order + 1
    coeff = np.full((n, ncoef), np.nan)
    err = np.full((n, ncoef), np.nan)
    npix = np.zeros(n, dtype=int)
    resid = np.full(n, np.nan)
    for k in range(n):
        geo = {**tr.ring_params(k), "pa": pa_use[k], "incl": incl_use[k]}
        sel = _select(x, y, v, geo, tr.radius[k], tr.width, "both")
        if sel.sum() < min_pixels:
            continue
        _, cos_t, sin_t = disc_coordinates(x[sel], y[sel], geo["x0"], geo["y0"], geo["pa"], geo["incl"])
        theta = np.arctan2(sin_t, cos_t)
        cols = [np.ones_like(theta)]
        for m in range(1, order + 1):
            cols += [np.cos(m * theta), np.sin(m * theta)]
        a = np.column_stack(cols)
        # the azimuthal coverage must be good enough to separate the terms
        if np.linalg.cond(a.T @ a) > 1e4:
            continue
        sol, *_ = np.linalg.lstsq(a, v[sel], rcond=None)
        r = v[sel] - a @ sol
        s2 = np.sum(r ** 2) / max(sel.sum() - ncoef, 1)
        cov = np.linalg.inv(a.T @ a) * s2
        coeff[k], err[k], npix[k] = sol, np.sqrt(np.diag(cov)), sel.sum()
        resid[k] = np.median(np.abs(r))
    return Harmonics(tr.radius, coeff, err, np.sin(np.radians(incl_use)), incl_use, npix, resid)


def summary(h: Harmonics, kpc_per_arcsec: float, r_max_arcsec: float | None = None) -> dict:
    """Trachternach et al.'s summary numbers for one galaxy."""
    ok = np.isfinite(h.a_r)
    if r_max_arcsec:
        ok &= h.radius <= r_max_arcsec
    inner = ok & (h.radius * kpc_per_arcsec <= 1.0)
    a = h.a_r[ok]
    eps = h.eps_pot[ok]
    eps = eps[np.isfinite(eps)]
    out = {
        "A_r_median": float(np.median(a)),
        "A_r_quartiles": [float(q) for q in np.percentile(a, [25, 75])],
        "A_r_1kpc_median": float(np.median(h.a_r[inner])) if inner.any() else None,
        "n_rings": int(ok.sum()),
        "n_rings_1kpc": int(inner.sum()),
        "A_m_median": {m: float(np.nanmedian(h.a_m(m)[ok])) for m in range(1, HARMONIC_ORDER + 1)},
        "eps_pot_mean": float(np.mean(eps)) if len(eps) else None,
        "eps_pot_sem": float(np.std(eps) / np.sqrt(len(eps))) if len(eps) > 1 else None,
        "residual_median_abs": float(np.nanmedian(h.residual_median_abs[ok])),
    }
    return out
