"""White et al. (2024)'s two time-dilation measurements, replicated.

The idea: Type Ia supernovae are similar enough to serve as standard clocks.
For a *target* supernova at redshift z observed in some band, collect the
light curves of all other supernovae seen at the same *rest-frame* wavelength
(a redder band for a farther supernova, a bluer one for a nearer). That pile
of points is the target's *reference light curve*.

Method 1 (``dispersion_curves``) stretches every reference point's time by
1/(1 + z_j)^b, using each point's own redshift z_j, and asks which b makes the
pile tightest. No light-curve model and no time-dilation assumption enter.

Method 2 (``fit_width``) first squeezes the reference by (1 + z_j), then finds
the width w by which the target's own light curve must be compressed to match
it. Expansion predicts w = 1 + z; tired light predicts w = 1. White et al.
fitted w by Nelder-Mead minimisation started at w = 1 + z, the expected
answer; ``fit_width`` also reports the global minimum on a grid and the
Nelder-Mead result started at w = 1, as a check that the start does not steer
the answer.

Everything here follows White et al.'s published code (Methods.py,
StretchMethod.py in github.com/ryanwhite1/DES-Time-Dilation), including its
choices: the target's chi-square is divided by (number of usable points - 1),
and a point is usable only when the median reference flux at its time is > 0.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import optimize

from config import (BANDS, DES_FILTERS, MIN_REFERENCE_POINTS, MIN_TARGET_POINTS, REF_DELTA,
                    REF_TAU_DAYS, WIDTH_MC_ITER)

#: Resolution (days) of the precomputed running-median profile.
PROFILE_STEP = 0.01
#: Width grid for the global-minimum check.
W_GRID = np.arange(0.40, 4.00, 0.002)


class MedianProfile:
    """Running median of a reference light curve on a fine time grid.

    ``value(x)`` is the median reference flux within +/- tau of x (0 where the
    window is empty), as in White et al.'s chi-square; precomputing it makes a
    chi-square evaluation cost a few array look-ups.
    """

    def __init__(self, t, flux, tau=REF_TAU_DAYS, step=PROFILE_STEP):
        order = np.argsort(t)
        self.t = np.asarray(t)[order]
        self.f = np.asarray(flux)[order]
        self.step = step
        self.start = self.t[0] - tau
        n = int(np.ceil((self.t[-1] + tau - self.start) / step)) + 1
        centres = self.start + step * np.arange(n)
        lo = np.searchsorted(self.t, centres - tau, side="left")
        hi = np.searchsorted(self.t, centres + tau, side="right")
        med = np.zeros(n)
        for j in np.nonzero(hi > lo)[0]:
            med[j] = np.median(self.f[lo[j]:hi[j]])
        self.med = med

    def value(self, x):
        j = np.rint((np.asarray(x) - self.start) / self.step).astype(np.int64)
        inside = (j >= 0) & (j < self.med.size)
        out = np.zeros(np.shape(x))
        out[inside] = self.med[j[inside]]
        return out


def reduced_chi2(profile: MedianProfile, t, flux, err, w):
    """White et al.'s objective for one trial width w (their stretch_chi_square)."""
    ref = profile.value(t / w)
    ok = ref > 0
    return np.sum(((flux[ok] - ref[ok]) / err[ok]) ** 2) / (ok.sum() - 1)


def chi2_on_grid(profile: MedianProfile, t, flux, err, w_grid, min_points=MIN_TARGET_POINTS):
    """The same objective for many widths at once; inf where too few points are usable.

    White et al.'s objective drops target points that fall where the reference
    has no data, so a width that pushes almost all points off the reference
    scores well by default. The grid search therefore requires at least
    ``min_points`` usable points before a width can count as the minimum.
    """
    ref = profile.value(t[None, :] / w_grid[:, None])
    ok = ref > 0
    n = ok.sum(axis=1)
    chi2 = np.where(ok, ((flux[None, :] - ref) / err[None, :]) ** 2, 0.0).sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = chi2 / (n - 1)
    return np.where(n >= min_points, out, np.inf)


def binned_dispersion(t, flux, n_bins):
    """Median of the flux standard deviations in ``n_bins`` equal-count time bins.

    White et al.'s ``binned_dispersions``: sort by time, cut into chunks of
    ceil(N / n_bins) points, take each chunk's standard deviation, return the
    median (chunks left empty at the end are ignored).
    """
    order = np.argsort(t)
    f = np.asarray(flux)[order]
    size = len(f) / n_bins
    size = int(size + 1) if round(size) != size else int(size)
    sig = [np.std(f[j * size:(j + 1) * size]) for j in range(n_bins) if j * size < len(f)]
    sig = np.array([s for s in sig if np.isfinite(s)])
    return float(np.median(sig))


class ReferenceSet:
    """All normalised light curves, grouped for fast reference-curve assembly."""

    def __init__(self, points: pd.DataFrame, sn: pd.DataFrame):
        ia = sn[sn.is_ia]
        self.cids = ia["cid"].to_numpy()
        self.zs = ia["z"].to_numpy()
        self.z_of = dict(zip(self.cids, self.zs))
        self.curves = {}
        for (cid, band), g in points.groupby(["cid", "band"], sort=False):
            self.curves[(cid, band)] = (g["t"].to_numpy(), g["flux"].to_numpy(),
                                        g["err"].to_numpy())

    def members(self, cid0, z0, fit_band, delta=REF_DELTA):
        """(cid, band, z) of every light curve that samples the target's rest wavelength.

        A light curve in band r at redshift z_r samples the same rest-frame
        wavelength as the target's band f at z when
        lambda_r (1 + z) / lambda_f = 1 + z_r, give or take delta times the
        target band's width (White et al.'s equation 3).
        """
        lam_f, fwhm_f = DES_FILTERS[fit_band]
        wiggle = delta * fwhm_f
        out = []
        for band in BANDS:
            lam_b = DES_FILTERS[band][0]
            upper = (lam_b * (1 + z0) + wiggle) / lam_f - 1
            lower = (lam_b * (1 + z0) - wiggle) / lam_f - 1
            sel = (self.zs >= lower) & (self.zs <= upper) & (self.cids != cid0)
            for cid, z in zip(self.cids[sel], self.zs[sel]):
                if (cid, band) in self.curves:
                    out.append((cid, band, z))
        return out

    def reference(self, cid0, z0, fit_band, *, power=1.0, delta=REF_DELTA):
        """Reference points: times divided by (1 + z_j)^power, fluxes, and z_j."""
        T, F, Z = [], [], []
        for cid, band, z in self.members(cid0, z0, fit_band, delta):
            t, f, _ = self.curves[(cid, band)]
            T.append(t / (1.0 + z) ** power)
            F.append(f)
            Z.append(np.full(t.size, z))
        if not T:
            return np.array([]), np.array([]), np.array([])
        return np.concatenate(T), np.concatenate(F), np.concatenate(Z)


def _nelder_mead(fun, x0):
    res = optimize.minimize(lambda x: fun(float(x[0])), [x0], method="Nelder-Mead")
    return float(res.x[0])


def fit_width(refset: ReferenceSet, cid, band, *, power=1.0, n_mc=WIDTH_MC_ITER, seed=0,
              require_prepeak=False):
    """One target light curve's width relative to its reference (White et al. method 2).

    Returns None when White et al.'s quality criteria fail, else a dict with
    the Nelder-Mead width started at 1 + z (theirs) and at 1, the global grid
    minimum, the Monte Carlo error, the error floor, and the adopted error.
    """
    if (cid, band) not in refset.curves:
        return None
    z = refset.z_of[cid]
    t, f, e = refset.curves[(cid, band)]
    if t.size < MIN_TARGET_POINTS:
        return None
    if require_prepeak and not np.any(t <= 0):
        return None
    rt, rf, _ = refset.reference(cid, z, band, power=power)
    if rt.size <= MIN_REFERENCE_POINTS:
        return None
    use = (t <= rt.max()) & (t >= rt.min())
    if use.sum() < MIN_TARGET_POINTS:
        return None
    t, f, e = t[use], f[use], e[use]

    profile = MedianProfile(rt, rf)
    with np.errstate(divide="ignore", invalid="ignore"):
        w_nm = _nelder_mead(lambda w: reduced_chi2(profile, t, f, e, w), 1.0 + z)
        w_nm1 = _nelder_mead(lambda w: reduced_chi2(profile, t, f, e, w), 1.0)
    grid = chi2_on_grid(profile, t, f, e, W_GRID)
    w_grid = float(W_GRID[np.argmin(grid)]) if np.isfinite(grid).any() else np.nan

    rng = np.random.default_rng(seed)
    draws = rng.normal(f, e, size=(n_mc, f.size))
    with np.errstate(divide="ignore", invalid="ignore"):
        mc = np.array([_nelder_mead(lambda w, d=d: reduced_chi2(profile, t, d, e, w), 1.0 + z)
                       for d in draws])
    err_mc = float(np.std(mc))
    err_floor = binned_dispersion(rt, rf, 30) * (1.0 + z)
    return {
        "cid": cid, "band": band, "z": z, "n_points": int(t.size), "n_ref": int(rt.size),
        "width": w_nm, "width_from_1": w_nm1, "width_grid": w_grid,
        "err": max(err_floor, err_mc), "err_floor": err_floor, "err_mc": err_mc,
    }


def dispersion_curve(refset: ReferenceSet, cid, band, b_values, n_bins=15):
    """Method 1 for one target: reference scatter as a function of b.

    The reference is assembled *without* any time scaling; each point's time is
    then divided by (1 + z_j)^b for every trial b. Returns None if the
    reference is too small (same criterion as the width fits).
    """
    if (cid, band) not in refset.curves:
        return None
    z = refset.z_of[cid]
    rt, rf, rz = refset.reference(cid, z, band, power=0.0)
    if rt.size <= MIN_REFERENCE_POINTS:
        return None
    return np.array([binned_dispersion(rt / (1.0 + rz) ** b, rf, n_bins) for b in b_values])


def average_widths(per_band: dict[str, pd.DataFrame], column="width") -> pd.DataFrame:
    """Per-supernova mean width over bands, as in White et al.'s final figure.

    Widths whose error exceeds the width itself are dropped (their cut). The
    error of the mean is sqrt(sum sigma^2) / n, i.e. the bands are treated as
    independent, as they did.
    """
    rows = []
    for band, df in per_band.items():
        good = df[df["err"] < df[column]]
        rows.append(good[["cid", "z", column, "err"]].assign(band=band))
    allw = pd.concat(rows)
    agg = allw.groupby("cid").agg(z=("z", "first"), width=(column, "mean"),
                                  n_bands=(column, "size"),
                                  err=("err", lambda s: np.sqrt(np.sum(s ** 2)) / len(s)))
    return agg.reset_index()


def fit_power(z, w, err):
    """Least-squares b in w = (1 + z)^b with absolute errors (as White et al.)."""
    popt, pcov = optimize.curve_fit(lambda x, b: x ** b, 1 + np.asarray(z), np.asarray(w),
                                    p0=(1.0,), sigma=np.asarray(err), absolute_sigma=True)
    b = float(popt[0])
    chi2 = float(np.sum((((1 + np.asarray(z)) ** b - w) / err) ** 2))
    return {"b": b, "err": float(np.sqrt(pcov[0, 0])), "chi2": chi2,
            "chi2_nu": chi2 / (len(w) - 1), "n": int(len(w))}


def fit_linear(z, w, err):
    """w = m (1 + z) + c, the linear form White et al. also quote."""
    popt, pcov = optimize.curve_fit(lambda x, m, c: m * x + c, 1 + np.asarray(z), np.asarray(w),
                                    p0=(1.0, 0.0), sigma=np.asarray(err), absolute_sigma=True)
    return {"slope": float(popt[0]), "slope_err": float(np.sqrt(pcov[0, 0])),
            "intercept": float(popt[1]), "intercept_err": float(np.sqrt(pcov[1, 1]))}
