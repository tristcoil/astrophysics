"""Delays by cross-correlation, the reverberation-mapping way (project I2).

The interpolated cross-correlation function (ICCF; Gaskell & Sparke 1986;
White & Peterson 1994) and its flux-randomization / random-subset-selection
error analysis (Peterson et al. 1998, 2004) are copied from I2's
``33-reverberation-mapping/iccf.py``, where they reproduce PyCCF to 1e-16.

What changes for lensed quasars is microlensing. In reverberation mapping the
two light curves differ by a lag and a smooth response; here each image also
drifts by tenths of a magnitude over the years as stars in the lens galaxy
move. Left in, that drift dominates the correlation. So before correlating,
each image is *high-pass filtered*: a stiff spline (knots every ``knot_ml``
days) fitted to that image alone, plus a zero point per telescope, is
subtracted. What survives are the quasar's own fast wiggles (weeks to a few
months), which all images share up to the delay.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from lightcurves import LightCurves
import delays as D

THRESHOLD = 0.8  # centroid uses r > THRESHOLD * r_max
MIN_RMAX = 0.2


def _one_way(ta, ya, tb, yb, lags):
    """r(tau) from interpolating curve b onto the shifted times of a (I2's code)."""
    shifted = ta[None, :] + lags[:, None]
    inside = (shifted >= tb[0]) & (shifted <= tb[-1])
    w = inside.astype(float)
    yb_i = np.interp(shifted.ravel(), tb, yb).reshape(shifted.shape)
    ya2 = np.broadcast_to(ya, shifted.shape)
    n = w.sum(axis=1)
    sa = (w * ya2).sum(axis=1)
    sb = (w * yb_i).sum(axis=1)
    saa = (w * ya2 * ya2).sum(axis=1)
    sbb = (w * yb_i * yb_i).sum(axis=1)
    sab = (w * ya2 * yb_i).sum(axis=1)
    da = n * saa - sa * sa
    db = n * sbb - sb * sb
    den = np.sqrt(np.clip(da, 0.0, None) * np.clip(db, 0.0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.where(den > 0.0, (n * sab - sa * sb) / den, 0.0)
    return r, n


def iccf(t1, y1, t2, y2, lags):
    """ICCF of curve 2 against curve 1 (average of both interpolation directions).

    A positive lag means curve 2 follows curve 1.
    """
    t1, y1, t2, y2 = (np.asarray(a, float) for a in (t1, y1, t2, y2))
    lags = np.asarray(lags, float)
    r21, n21 = _one_way(t2, y2, t1, y1, -lags)
    r12, n12 = _one_way(t1, y1, t2, y2, lags)
    return 0.5 * (r12 + r21), 0.5 * (n12 + n21)


def peak_and_centroid(lags, r, threshold=THRESHOLD, min_rmax=MIN_RMAX):
    """(ok, r_max, peak lag, centroid lag) with PyCCF's rules (I2's code)."""
    i = int(np.argmax(r))
    r_max = float(r[i])
    if r_max < min_rmax or i == 0 or i == len(r) - 1:
        return False, r_max, np.nan, np.nan
    below = r < threshold * r_max
    left = np.nonzero(below[:i])[0]
    right = np.nonzero(below[i + 1:])[0]
    if left.size == 0 or right.size == 0:
        return False, r_max, float(lags[i]), np.nan
    sel = slice(left[-1] + 1, i + 1 + right[0])
    return True, r_max, float(lags[i]), float(np.sum(r[sel] * lags[sel]) / np.sum(r[sel]))


def _rss(rng, t, y, e):
    """Random subset selection: N draws with replacement, duplicates merged."""
    n = len(t)
    idx = rng.integers(0, n, n)
    uniq, counts = np.unique(idx, return_counts=True)
    return t[uniq], y[uniq], e[uniq] / np.sqrt(counts)


# ---------------------------------------------------------------------------
# Microlensing removal
# ---------------------------------------------------------------------------


def highpass(lc: LightCurves, im: str, knot_ml: float = 200.0, curv: float = 0.3) -> np.ndarray:
    """Image ``im`` minus its own stiff spline and telescope zero points."""
    t = lc.t
    kn = D.uniform_knots(t.min() - 1.0, t.max() + 1.0, knot_ml)
    B = D.basis(t, kn).toarray()
    tels = lc.telescopes
    O = np.zeros((lc.n, max(len(tels) - 1, 0)))
    for j, tel in enumerate(tels[1:]):
        O[lc.tel == tel, j] = 1.0
    A = np.hstack([B, O])
    w = 1.0 / lc.err[im] ** 2
    P = np.zeros((A.shape[1], A.shape[1]))
    d2 = D.second_difference(B.shape[1])
    P[: B.shape[1], : B.shape[1]] = d2.T @ d2 / curv**2
    P[B.shape[1]:, B.shape[1]:] += np.eye(O.shape[1])
    coef = np.linalg.solve(A.T @ (A * w[:, None]) + P, A.T @ (w * lc.mag[im]))
    return lc.mag[im] - A @ coef


@dataclass
class CCFDelay:
    pair: tuple  # (reference image, image)
    lags: np.ndarray
    r: np.ndarray
    ok: bool
    r_max: float
    tau_peak: float
    tau_cent: float
    cent: tuple  # (median, minus, plus) of the FR/RSS centroid distribution
    success: float
    overlap: np.ndarray = field(repr=False, default=None)  # fraction of points overlapping at each lag
    cccd: np.ndarray = field(repr=False, default=None)


def ccf_delay(lc: LightCurves, ref: str, im: str, lag_lo: float, lag_hi: float, step: float = 0.25,
              knot_ml: float = 200.0, n_mc: int = 500, seed: int = 1, min_overlap: float = 0.0) -> CCFDelay:
    """Delay of image ``im`` behind ``ref`` (t_im - t_ref) from the high-passed ICCF.

    The lag search is restricted to [lag_lo, lag_hi]; the FR/RSS centroid
    distribution gives the error. Lags at which less than ``min_overlap`` of
    the points overlap (shifted past the end of a single season) get r = 0.
    """
    lags = np.arange(lag_lo, lag_hi + step / 2, step)
    y1 = highpass(lc, ref, knot_ml)
    y2 = highpass(lc, im, knot_ml)
    e1, e2 = lc.err[ref], lc.err[im]
    t = lc.t
    r, npts = iccf(t, y1, t, y2, lags)
    overlap = npts / len(t)
    if min_overlap > 0:
        r = np.where(overlap >= min_overlap, r, 0.0)
    ok, rmax, tp, tc = peak_and_centroid(lags, r)
    rng = np.random.default_rng(seed)
    cents = []
    for _ in range(n_mc):
        a = _rss(rng, t, y1, e1)
        b = _rss(rng, t, y2, e2)
        rr, nn = iccf(a[0], rng.normal(a[1], a[2]), b[0], rng.normal(b[1], b[2]), lags)
        if min_overlap > 0:
            rr = np.where(nn / len(a[0]) >= min_overlap, rr, 0.0)
        good, _, _, c = peak_and_centroid(lags, rr)
        if good:
            cents.append(c)
    cents = np.array(cents)
    if cents.size:
        lo, med, hi = np.percentile(cents, [15.865525, 50.0, 84.134475])
        summ = (float(med), float(med - lo), float(hi - med))
    else:
        summ = (np.nan, np.nan, np.nan)
    return CCFDelay((ref, im), lags, r, ok, rmax, tp, tc, summ, cents.size / max(n_mc, 1), cents, overlap)
