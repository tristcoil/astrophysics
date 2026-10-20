"""Eclipse timing of contact binaries: time systems, light-curve templates, O - C.

A contact binary varies all the time, so the best timing uses the whole light
curve, not just the bottom of each eclipse (Molnar et al. 2017; Kovacs et al.
2019). We fit a template light curve, shifted in phase, scaled and offset, to
any group of photometry (one night, one season, one TESS sector) and read off
the shift. Applied to every data set with the same template and the same
reference ephemeris, the shifts form the O - C diagram.

Times are converted to BJD_TDB first: survey archives use JD_UTC, HJD_UTC or
MJD in various conventions, and at the minute level (TDB - UTC is ~65 s) these
differences matter.
"""

from __future__ import annotations

import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time

import config as C


# ---------------------------------------------------------------------------
# Time systems
# ---------------------------------------------------------------------------


#: Light-travel corrections are computed for the geocentre; the observatory's
#: offset from it changes them by at most 21 ms.
GEOCENTRE = EarthLocation.from_geocentric(0.0, 0.0, 0.0, unit=u.m)


def _coord(ra, dec):
    return SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")


def jd_utc_to_bjd_tdb(jd_utc, ra, dec):
    """Geocentric JD(UTC) -> BJD(TDB). The observatory's offset from the
    geocentre (< 21 ms) is ignored."""
    t = Time(np.asarray(jd_utc, dtype=float), format="jd", scale="utc", location=GEOCENTRE)
    ltt = t.light_travel_time(_coord(ra, dec), kind="barycentric")
    return (t.tdb + ltt).jd


def hjd_utc_to_bjd_tdb(hjd_utc, ra, dec):
    """HJD(UTC) -> BJD(TDB): undo the heliocentric correction (evaluated at the
    HJD itself, an error far below a millisecond), then apply the barycentric one."""
    t = Time(np.asarray(hjd_utc, dtype=float), format="jd", scale="utc", location=GEOCENTRE)
    c = _coord(ra, dec)
    jd = t - t.light_travel_time(c, kind="heliocentric")
    return (jd.tdb + jd.light_travel_time(c, kind="barycentric")).jd


# ---------------------------------------------------------------------------
# Templates: Fourier series in orbital phase
# ---------------------------------------------------------------------------


class Template:
    """A periodic light curve f(phase) as a Fourier series, phase zero at the
    deeper minimum."""

    def __init__(self, coef: np.ndarray):
        self.coef = np.asarray(coef, dtype=float)  # [a0, a1, b1, a2, b2, ...]
        self.order = (self.coef.size - 1) // 2

    @staticmethod
    def design(phase, order):
        cols = [np.ones_like(phase)]
        for k in range(1, order + 1):
            cols += [np.cos(2 * np.pi * k * phase), np.sin(2 * np.pi * k * phase)]
        return np.column_stack(cols)

    @classmethod
    def fit(cls, phase, flux, order=12, weights=None):
        A = cls.design(phase, order)
        w = np.ones_like(flux) if weights is None else weights
        coef, *_ = np.linalg.lstsq(A * w[:, None], flux * w, rcond=None)
        tpl = cls(coef)
        # move phase zero to the deeper minimum, and normalise to mean 1
        ph = np.linspace(0, 1, 20001, endpoint=False)
        shift = ph[np.argmin(tpl(ph))]
        tpl = tpl.shifted(shift)
        tpl.coef = tpl.coef / tpl.coef[0]
        return tpl, shift

    def __call__(self, phase):
        return self.design(np.asarray(phase, dtype=float), self.order) @ self.coef

    def shifted(self, dphi):
        """Template g(phase) = f(phase + dphi)."""
        c = self.coef.copy()
        for k in range(1, self.order + 1):
            a, b = c[2 * k - 1], c[2 * k]
            ck, sk = np.cos(2 * np.pi * k * dphi), np.sin(2 * np.pi * k * dphi)
            # cos(k(x+d)) = cos kx cos kd - sin kx sin kd ; sin(k(x+d)) = sin kx cos kd + cos kx sin kd
            c[2 * k - 1] = a * ck + b * sk
            c[2 * k] = -a * sk + b * ck
        return Template(c)

    def depth_ratio(self):
        ph = np.linspace(0, 1, 4001, endpoint=False)
        f = self(ph)
        fmax = f.max()
        return (fmax - f[(ph < 0.25) | (ph > 0.75)].min()) / (fmax - f[(ph > 0.25) & (ph < 0.75)].min())


# ---------------------------------------------------------------------------
# Fitting a phase shift
# ---------------------------------------------------------------------------


def fit_shift(t, flux, err, tpl: Template, p_ref, t_ref, *, window=0.25, n_grid=501,
              trend=False, groups=None):
    """Least-squares phase shift of ``tpl`` for photometry (t, flux, err).

    Model: flux = offset_g + scale_g * tpl(phase - delta), with
    phase = (t - t_ref)/p_ref and a separate offset and scale for each label
    in ``groups`` (e.g. survey cameras, or AAVSO observer-band pairs, each its
    own photometric system). For a given delta the groups decouple, so each
    one's offset and scale come from its own 2x2 normal equations: the cost
    per trial shift is one pass over the data. ``delta`` is searched on a grid
    of +-``window`` cycles and refined with parabolas. Scales must be
    positive: a light curve close to a pure second harmonic (a contact binary
    at moderate inclination) looks the same flipped and shifted by a quarter
    cycle. ``trend`` adds a linear slope in time per group.
    Returns delta (cycles), its 1-sigma error (chi^2 rescaled to dof), the
    reduced chi^2, the number of points and the best-fit model.
    """
    t = np.asarray(t, float)
    y = np.asarray(flux, float)
    w2 = 1.0 / np.asarray(err, float) ** 2
    phase = (t - t_ref) / p_ref
    lab = np.zeros(t.size, int) if groups is None else np.unique(groups, return_inverse=True)[1]
    n_lab = int(lab.max()) + 1
    tm = t.mean()
    tt = (t - tm) / max(1.0, np.ptp(t))

    def sums(v):
        return np.bincount(lab, weights=v, minlength=n_lab)

    S = sums(w2)
    Sy = sums(w2 * y)
    Syy = sums(w2 * y * y)
    if trend:
        St, Stt, Sty = sums(w2 * tt), sums(w2 * tt * tt), sums(w2 * tt * y)

    basis = Template.design(np.mod(phase, 1.0), tpl.order)  # f(phase - delta) = basis @ rotated coefficients

    def solve(delta, want_model=False):
        f = basis @ tpl.shifted(-delta).coef
        Sf, Sff, Sfy = sums(w2 * f), sums(w2 * f * f), sums(w2 * f * y)
        if not trend:
            det = S * Sff - Sf * Sf
            b = (S * Sfy - Sf * Sy) / det
            a = (Sy - b * Sf) / S
            chi = float(np.sum(Syy - a * Sy - b * Sfy))
            if np.any(b <= 0):
                chi = np.inf
            if want_model:
                return chi, a[lab] + b[lab] * f
            return chi, None
        # with a slope: 3x3 normal equations per group
        Sft = sums(w2 * f * tt)
        M = np.empty((n_lab, 3, 3))
        M[:, 0, 0], M[:, 0, 1], M[:, 0, 2] = S, Sf, St
        M[:, 1, 0], M[:, 1, 1], M[:, 1, 2] = Sf, Sff, Sft
        M[:, 2, 0], M[:, 2, 1], M[:, 2, 2] = St, Sft, Stt
        rhs = np.stack([Sy, Sfy, Sty], axis=1)
        c = np.linalg.solve(M, rhs[..., None])[..., 0]
        chi = float(np.sum(Syy - np.sum(c * rhs, axis=1)))
        if np.any(c[:, 1] <= 0):
            chi = np.inf
        if want_model:
            return chi, c[lab, 0] + c[lab, 1] * f + c[lab, 2] * tt
        return chi, None

    grid = np.linspace(-window, window, n_grid)
    chi = np.array([solve(d)[0] for d in grid])
    j = int(np.argmin(chi))
    lo, hi = max(j - 3, 0), min(j + 4, grid.size)
    ok = np.isfinite(chi[lo:hi])
    pp = np.polyfit(grid[lo:hi][ok] - grid[j], chi[lo:hi][ok], 2) if ok.sum() >= 3 else np.zeros(3)
    d0 = grid[j] - pp[1] / (2 * pp[0]) if pp[0] > 0 else grid[j]
    # final evaluation and a local curvature on a fine grid
    fine = np.linspace(d0 - 3 * (grid[1] - grid[0]), d0 + 3 * (grid[1] - grid[0]), 41)
    cf = np.array([solve(d)[0] for d in fine])
    okf = np.isfinite(cf)
    pf = np.polyfit(fine[okf] - d0, cf[okf], 2) if okf.sum() >= 3 else np.zeros(3)
    d_best = d0 - pf[1] / (2 * pf[0]) if pf[0] > 0 else d0
    chi_best, model = solve(d_best, want_model=True)
    dof = max(1, t.size - (2 + int(trend)) * n_lab - 1)
    red = chi_best / dof
    sig = np.sqrt(red / pf[0]) if pf[0] > 0 else np.nan
    return {"delta": float(d_best), "delta_err": float(sig), "chi2_red": float(red), "n": int(t.size),
            "model": model, "t_mean": float(tm), "t_min": float(t.min()), "t_max": float(t.max())}


def sigma_clip_fit(t, flux, err, tpl, p_ref, t_ref, nsig=4.0, iters=4, **kw):
    """fit_shift with iterative rejection of outliers (|resid| > nsig * robust rms)."""
    keep = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    res = None
    for _ in range(iters):
        res = fit_shift(t[keep], flux[keep], err[keep], tpl, p_ref, t_ref,
                        groups=None if kw.get("groups") is None else kw["groups"][keep],
                        **{k: v for k, v in kw.items() if k != "groups"})
        r = flux[keep] - res["model"]
        mad = 1.4826 * np.median(np.abs(r - np.median(r)))
        bad = np.abs(r) > nsig * mad
        if not bad.any():
            break
        idx = np.nonzero(keep)[0]
        keep[idx[bad]] = False
    res["keep"] = keep
    return res


def minima_times(t, flux, err, tpl, p_ref, t_ref, half_width=0.12, min_points=8):
    """Individual primary (phase 0) and secondary (phase 0.5) minimum times.

    Each minimum is fitted with the template restricted to +-half_width in
    phase, with its own offset and scale. Returns arrays of (cycle number E
    with .5 for secondaries, time of minimum, error, reduced chi^2).
    """
    phase = (t - t_ref) / p_ref
    e_lo, e_hi = np.floor(phase.min()), np.ceil(phase.max())
    out = []
    for e in np.arange(e_lo, e_hi + 0.5, 0.5):
        s = np.abs(phase - e) < half_width
        if s.sum() < min_points:
            continue
        # require points on both sides of the minimum
        if (phase[s] < e).sum() < 3 or (phase[s] > e).sum() < 3:
            continue
        tp = tpl if e % 1 == 0 else tpl.shifted(0.5)
        r = fit_shift(t[s], flux[s], err[s], tp, p_ref, t_ref + e * p_ref, window=0.05, n_grid=101)
        if not np.isfinite(r["delta_err"]):
            continue
        out.append((e, t_ref + (e + r["delta"]) * p_ref, r["delta_err"] * p_ref, r["chi2_red"]))
    out = np.array(out)
    return {"E": out[:, 0], "t_min": out[:, 1], "err": out[:, 2], "chi2_red": out[:, 3]} if out.size else None
