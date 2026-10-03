"""TESS photometry of 55 Cnc e: the transit that settled an alias, and a radius.

A planet crossing its star blocks a fraction ~(Rp/R*)^2 of the light, but a
real star is brighter at the centre than at the limb ("limb darkening"), so
the dip is deeper than the area ratio when the planet crosses the middle.
``transit_flux`` integrates a quadratic limb-darkening law exactly (to
numerical precision) over the part of the disc the planet covers, ring by
ring, rather than using the uniform-disc approximation of project A1: for a
0.03%-deep transit that difference is a 10% error in radius.

Limb-darkening coefficients for the TESS band come from Claret (2017),
interpolated to the star's temperature, gravity and metallicity.
"""

from __future__ import annotations

import gzip
import warnings

import numpy as np
from scipy import optimize

from config import CLARET, TEFF_LOGG_FEH_55CNC
from fetch_data import tess_files

R_SUN_M = 6.957e8
R_EARTH_M = 6.3781e6


# ---------------------------------------------------------------------------
# Limb darkening and the transit model
# ---------------------------------------------------------------------------

def claret_quadratic(teff: float, logg: float, feh: float, xi: float = 2.0) -> tuple[float, float]:
    """Quadratic TESS limb-darkening (u1, u2) from Claret (2017) table 25
    (ATLAS models, least-squares method), trilinearly interpolated in
    (Teff, log g, [M/H]) at microturbulence ``xi`` (the only value tabulated
    at non-solar metallicity is 2 km/s)."""
    rows = []
    with gzip.open(CLARET.local("table25.dat.gz"), "rt") as fh:
        for line in fh:
            try:
                rows.append((float(line[0:5]), float(line[6:12]), float(line[13:17]),
                             float(line[19:22]), float(line[23:31]), float(line[33:41])))
            except ValueError:
                continue
    tab = np.array(rows)
    tab = tab[np.isclose(tab[:, 3], xi)]
    from scipy.interpolate import LinearNDInterpolator
    sel = (np.abs(tab[:, 1] - teff) <= 500) & (np.abs(tab[:, 0] - logg) <= 0.6) & (np.abs(tab[:, 2] - feh) <= 0.6)
    pts = tab[sel][:, [1, 0, 2]]
    u1 = LinearNDInterpolator(pts, tab[sel][:, 4])(teff, logg, feh)
    u2 = LinearNDInterpolator(pts, tab[sel][:, 5])(teff, logg, feh)
    return float(u1), float(u2)


def _intensity(r, u1, u2):
    mu = np.sqrt(np.clip(1.0 - r * r, 0.0, 1.0))
    return 1.0 - u1 * (1.0 - mu) - u2 * (1.0 - mu) ** 2


def transit_flux(z, p, u1, u2, *, n_rings: int = 400):
    """Relative flux of a quadratically limb-darkened star occulted by a dark
    disc of radius ratio ``p`` at projected separation ``z`` (units of R*).

    For each ring of the stellar disc at radius r, the planet covers an arc
    of angle 2 arccos((r^2 + z^2 - p^2) / (2 r z)); the blocked flux is the
    intensity-weighted integral of that arc over r, divided by the total,
    pi (1 - u1/3 - u2/6). Integrated with Gauss-Legendre nodes across the
    radii the planet actually spans, so a small planet is resolved finely.
    """
    z = np.atleast_1d(np.asarray(z, float))
    out = np.ones_like(z)
    hit = z < 1.0 + p
    if not np.any(hit):
        return out
    zh = z[hit]
    r_lo = np.clip(zh - p, 0.0, 1.0)
    r_hi = np.clip(zh + p, 0.0, 1.0)
    x, w = np.polynomial.legendre.leggauss(n_rings)
    r = 0.5 * (r_hi - r_lo)[:, None] * (x[None, :] + 1.0) + r_lo[:, None]
    wr = 0.5 * (r_hi - r_lo)[:, None] * w[None, :]
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (r * r + zh[:, None] ** 2 - p * p) / (2.0 * r * zh[:, None])
    angle = np.where(c <= -1.0, 2 * np.pi, np.where(c >= 1.0, 0.0, 2.0 * np.arccos(np.clip(c, -1, 1))))
    # a ring entirely inside the planet (possible only when z < p)
    angle = np.where(r < p - zh[:, None], 2 * np.pi, angle)
    blocked = np.sum(wr * r * angle * _intensity(r, u1, u2), axis=1)
    total = np.pi * (1.0 - u1 / 3.0 - u2 / 6.0)
    out[hit] = 1.0 - blocked / total
    return out


def projected_separation(t, tc, P, a_over_r, inc_rad):
    """Sky-projected star-planet distance (stellar radii), circular orbit.
    Only the half-orbit in front of the star can transit."""
    phase = 2 * np.pi * (t - tc) / P
    x = a_over_r * np.sin(phase)
    y = a_over_r * np.cos(phase) * np.cos(inc_rad)
    z = np.hypot(x, y)
    front = np.cos(phase) > 0
    return np.where(front, z, 10.0 + a_over_r)


def transit_model(t, tc, P, p, b, a_over_r, u1, u2):
    inc = np.arccos(np.clip(b / a_over_r, 0.0, 1.0))
    return transit_flux(projected_separation(t, tc, P, a_over_r, inc), p, u1, u2)


def a_over_r_from_density(P_days: float, rho_star_solar: float) -> float:
    """Kepler's third law in stellar units: (a/R*)^3 = G rho* P^2 / (3 pi)."""
    G = 6.674e-11
    rho = rho_star_solar * 1408.0  # mean solar density, kg/m^3
    P = P_days * 86400.0
    return float(np.cbrt(G * rho * P ** 2 / (3 * np.pi)))


# ---------------------------------------------------------------------------
# Light curves
# ---------------------------------------------------------------------------

def load_light_curves(quiet: bool = True):
    """Stitch the cached SPOC 2-min light curves: PDCSAP flux, good quality
    only, each sector normalised to its median. Returns (t_bjd, flux, err,
    sector)."""
    import lightkurve as lk

    ts, fs, es, ss = [], [], [], []
    for path in tess_files():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            lc = lk.read(str(path), flux_column="pdcsap_flux", quality_bitmask="default")
        lc = lc.remove_nans()
        f = np.asarray(lc.flux.value, float)
        e = np.asarray(lc.flux_err.value, float)
        med = np.nanmedian(f)
        ts.append(np.asarray(lc.time.value, float) + 2457000.0)  # BTJD -> BJD(TDB)
        fs.append(f / med)
        es.append(e / med)
        ss.append(np.full(f.size, int(lc.meta.get("SECTOR", -1))))
    t = np.concatenate(ts); f = np.concatenate(fs); e = np.concatenate(es); s = np.concatenate(ss)
    order = np.argsort(t)
    return t[order], f[order], e[order], s[order]


def flatten(t, f, *, window_d: float = 0.5, mask=None):
    """Divide out slow variability with a running median of width
    ``window_d`` (in-transit points masked so they do not drag it down).
    The transit lasts 1.6 h, so a 12-h window leaves it untouched."""
    trend = np.empty_like(f)
    use = np.ones(f.size, bool) if mask is None else ~mask
    tu, fu = t[use], f[use]
    half = window_d / 2
    lo = np.searchsorted(tu, t - half)
    hi = np.searchsorted(tu, t + half)
    for i in range(t.size):
        trend[i] = np.median(fu[lo[i]:hi[i]]) if hi[i] > lo[i] else np.nan
    return f / trend


def bls_search(t, f, period_min=0.5, period_max=5.0, durations=(0.04, 0.06, 0.08)):
    """Box least squares, as in project A1: the most transit-like period."""
    from astropy.timeseries import BoxLeastSquares

    bls = BoxLeastSquares(t, f)
    periods = np.exp(np.linspace(np.log(period_min), np.log(period_max), 60000))
    res = bls.power(periods, durations, objective="snr")
    i = int(np.argmax(res.power))
    return {"periods": periods, "power": np.asarray(res.power), "best_period": float(periods[i]),
            "t0": float(res.transit_time[i]), "depth": float(res.depth[i]),
            "duration": float(res.duration[i]), "snr": float(res.depth_snr[i])}


def fold_binned(t, f, P, tc, *, bin_minutes: float = 10.0, span=0.5):
    """Fold on (P, tc) and average in bins; returns phase (days from mid-
    transit), mean flux, standard error, count."""
    dt = (np.mod(t - tc + 0.5 * P, P) - 0.5 * P)
    sel = np.abs(dt) < span * P
    edges = np.arange(-span * P, span * P + 1e-9, bin_minutes / 1440.0)
    idx = np.digitize(dt[sel], edges) - 1
    n = np.bincount(idx, minlength=edges.size - 1)[:edges.size - 1]
    s = np.bincount(idx, weights=f[sel], minlength=edges.size - 1)[:edges.size - 1]
    s2 = np.bincount(idx, weights=f[sel] ** 2, minlength=edges.size - 1)[:edges.size - 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s / n
        sd = np.sqrt(np.maximum(s2 / n - mean ** 2, 0)) / np.sqrt(n)
    centres = 0.5 * (edges[1:] + edges[:-1])
    ok = n > 5
    return centres[ok], mean[ok], sd[ok], n[ok]


def fit_transit(t, f, err, P0, tc0, a_over_r0, u1, u2, *, window_d: float = 0.25):
    """Least-squares fit of (tc, P, Rp/R*, b, a/R*) to the flattened
    photometry within ``window_d`` of each transit. Errors from the scatter
    of the out-of-transit data (not the pipeline errors) and the Jacobian."""
    dt = np.mod(t - tc0 + 0.5 * P0, P0) - 0.5 * P0
    sel = np.abs(dt) < window_d
    ts, fs = t[sel], f[sel]
    oot = np.abs(dt[sel]) > 0.06
    sigma = 1.4826 * np.median(np.abs(fs[oot] - np.median(fs[oot])))

    def resid(x):
        tc, P, p, b, aR, f0 = x
        return (fs - f0 * transit_model(ts, tc, P, p, b, aR, u1, u2)) / sigma

    x0 = [tc0, P0, 0.0185, 0.4, a_over_r0, 1.0]
    lo = [tc0 - 0.02, P0 * (1 - 1e-4), 0.005, 0.0, 1.5, 0.99]
    hi = [tc0 + 0.02, P0 * (1 + 1e-4), 0.05, 1.0, 8.0, 1.01]
    sol = optimize.least_squares(resid, x0, bounds=(lo, hi), x_scale="jac")
    J = sol.jac
    cov = np.linalg.inv(J.T @ J) * max(1.0, np.sum(sol.fun ** 2) / (sol.fun.size - len(x0)))
    err_x = np.sqrt(np.diag(cov))
    names = ["tc", "P", "rp_rs", "b", "a_rs", "f0"]
    return {"params": dict(zip(names, map(float, sol.x))), "errors": dict(zip(names, map(float, err_x))),
            "cov": cov, "sigma_per_point": float(sigma), "n_points": int(sel.sum())}
