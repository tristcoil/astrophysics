"""Straight-line fits of the Leavitt law, as OGLE does them.

    m = a + b (log P - x0)

fitted by ordinary least squares with iterative sigma clipping. Unweighted
fits are standard here: the photometric errors of OGLE's mean magnitudes
(~0.01 mag) are far smaller than the law's own scatter (0.08-0.2 mag), which
comes from the width of the instability strip, dust and depth. Clipping
removes blends and misidentified modes. Errors on a and b come from the scatter
of the kept points.
"""

from __future__ import annotations

import numpy as np
from scipy import stats

import config as C


def fit_line(x, y, *, slope=None, x0=C.LOGP_PIVOT, nsigma=C.CLIP_SIGMA, iters=50, keep=None):
    """Iteratively clipped least-squares line. Returns a dict with ``mask`` of kept points.

    With ``slope`` given, only the zero point is fitted (the standard way to
    compare galaxies: one slope, different distances).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if keep is not None:
        ok &= keep
    mask = ok.copy()
    for _ in range(iters):
        xx, yy = x[mask] - x0, y[mask]
        if slope is None:
            A = np.column_stack([np.ones_like(xx), xx])
            coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
            a, b = coef
        else:
            a, b = np.mean(yy - slope * xx), slope
        r = y - (a + b * (x - x0))
        s = np.std(r[mask], ddof=2 if slope is None else 1)
        new = ok & (np.abs(r) < nsigma * s)
        if np.array_equal(new, mask):
            break
        mask = new
    xx = x[mask] - x0
    n = int(mask.sum())
    if slope is None:
        A = np.column_stack([np.ones_like(xx), xx])
        cov = s ** 2 * np.linalg.inv(A.T @ A)
        ea, eb = np.sqrt(np.diag(cov))
    else:
        cov = np.array([[s ** 2 / n, 0.0], [0.0, 0.0]])
        ea, eb = np.sqrt(cov[0, 0]), 0.0
    return {"a": float(a), "b": float(b), "ea": float(ea), "eb": float(eb), "sigma": float(s), "n": n,
            "n_clipped": int(ok.sum() - n), "x0": x0, "cov": cov, "mask": mask,
            "resid": y - (a + b * (x - x0))}


def predict(fit, x):
    return fit["a"] + fit["b"] * (np.asarray(x) - fit["x0"])


def fit_broken(x, y, x_break=1.0, *, keep=None, nsigma=C.CLIP_SIGMA):
    """Two independent lines either side of ``x_break`` vs one line, on the same clipped stars.

    Returns both slopes and an F-test of whether the extra two parameters are
    warranted (Kanbur & Ngeow's test for the 10-day break).
    """
    one = fit_line(x, y, keep=keep, nsigma=nsigma)
    m = one["mask"]
    x, y = np.asarray(x, float), np.asarray(y, float)
    lo, hi = m & (x < x_break), m & (x >= x_break)
    f_lo = fit_line(x, y, keep=lo, nsigma=1e9)
    f_hi = fit_line(x, y, keep=hi, nsigma=1e9)
    rss1 = np.sum(one["resid"][m] ** 2)
    rss2 = np.sum(f_lo["resid"][lo] ** 2) + np.sum(f_hi["resid"][hi] ** 2)
    n = int(m.sum())
    F = ((rss1 - rss2) / 2) / (rss2 / (n - 4))
    p = float(stats.f.sf(F, 2, n - 4))
    return {"single": one, "short": f_lo, "long": f_hi, "F": float(F), "p": p, "n": n,
            "n_long": int(hi.sum()), "slope_diff": f_hi["b"] - f_lo["b"],
            "slope_diff_err": float(np.hypot(f_hi["eb"], f_lo["eb"]))}


def fit_plane_pl(logp, color, mag, *, keep=None, nsigma=C.CLIP_SIGMA, x0=C.LOGP_PIVOT):
    """Period-luminosity-colour relation m = a + b (log P - x0) + c * colour, clipped."""
    logp, color, mag = (np.asarray(v, float) for v in (logp, color, mag))
    ok = np.isfinite(logp) & np.isfinite(color) & np.isfinite(mag)
    if keep is not None:
        ok &= keep
    mask = ok.copy()
    for _ in range(50):
        A = np.column_stack([np.ones(mask.sum()), logp[mask] - x0, color[mask]])
        coef, *_ = np.linalg.lstsq(A, mag[mask], rcond=None)
        r = mag - (coef[0] + coef[1] * (logp - x0) + coef[2] * color)
        s = np.std(r[mask], ddof=3)
        new = ok & (np.abs(r) < nsigma * s)
        if np.array_equal(new, mask):
            break
        mask = new
    A = np.column_stack([np.ones(mask.sum()), logp[mask] - x0, color[mask]])
    cov = s ** 2 * np.linalg.inv(A.T @ A)
    return {"a": float(coef[0]), "b": float(coef[1]), "c": float(coef[2]),
            "err": np.sqrt(np.diag(cov)).tolist(), "sigma": float(s), "n": int(mask.sum()), "mask": mask}
