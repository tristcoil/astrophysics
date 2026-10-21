"""The interpolated cross-correlation function and its error analysis.

This is our own implementation of the standard reverberation-mapping recipe:

* **ICCF** (Gaskell & Sparke 1986; Gaskell & Peterson 1987; White & Peterson
  1994). Two unevenly sampled light curves cannot be correlated point by point,
  so for each trial lag tau one curve is linearly interpolated onto the other's
  shifted times. Only the points that land inside the other curve's time span
  are used, and the Pearson coefficient is computed from *their* means and
  variances. This is done both ways (interpolating the line, then the
  continuum) and the two coefficients are averaged.
* **Peak and centroid.** tau_peak is where r is largest; tau_cent is the
  r-weighted mean lag over the contiguous region around the peak where
  r > 0.8 r_max. The centroid is the lag used for black hole masses
  (Koratkar & Gaskell 1991; Peterson et al. 2004).
* **FR/RSS errors** (Peterson et al. 1998; 2004). Each Monte Carlo realization
  draws N points with replacement from each light curve (random subset
  selection; a point drawn m times keeps one copy with its error divided by
  sqrt(m)) and adds Gaussian noise of each point's error (flux randomization).
  The lag is the median of the resulting centroid (or peak) distribution, the
  error its central 68%.
* **Alias weighting** for multi-year survey light curves (Grier et al. 2017,
  2019; Shen et al. 2024): lags at which most shifted line points fall into
  the continuum's seasonal gaps are down-weighted by P(tau) = [N(tau)/N(0)]^2
  convolved with the continuum's autocorrelation, and only the primary peak of
  the weighted lag distribution is kept.

Conventions match PyCCF (Sun, Grier & Peterson 2018), the code behind most
published lags since 2018: a positive lag means light curve 2 (the line)
follows light curve 1 (the continuum); a realization "fails" if its r_max is
below 0.2, if its peak sits on the edge of the lag range, or if r never drops
below 0.8 r_max on both sides of the peak.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

THRESHOLD = 0.8  # centroid uses r > THRESHOLD * r_max
MIN_RMAX = 0.2  # PyCCF's "sigmode": realizations with r_max below this fail


# ---------------------------------------------------------------------------
# The cross-correlation function
# ---------------------------------------------------------------------------


def _one_way(ta, ya, tb, yb, lags):
    """r(tau) from interpolating light curve b onto the shifted times of a.

    For each lag: times ta + tau inside [tb.min, tb.max] get yb interpolated
    there; r is the Pearson coefficient of those pairs only.
    Returns (r, n) arrays over ``lags``.
    """
    shifted = ta[None, :] + lags[:, None]  # (nlag, na)
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


def ccf(t1, y1, t2, y2, lags, mode: str = "both"):
    """Interpolated cross-correlation function of light curve 2 against 1.

    ``mode``: "both" (average of the two interpolation directions, the
    standard choice), "interp2" (only light curve 2 interpolated) or
    "interp1" (only light curve 1 interpolated). Times must be increasing.
    Returns (r, npts): r(tau) on ``lags`` and the mean number of overlapping
    points. Lags with no overlap at all get r = 0, as in PyCCF.
    """
    t1, y1, t2, y2 = (np.asarray(a, float) for a in (t1, y1, t2, y2))
    lags = np.asarray(lags, float)
    if mode == "interp2":
        return _one_way(t1, y1, t2, y2, lags)
    r21, n21 = _one_way(t2, y2, t1, y1, -lags)  # interpolate 1 at t2 - tau
    if mode == "interp1":
        return r21, n21
    r12, n12 = _one_way(t1, y1, t2, y2, lags)  # interpolate 2 at t1 + tau
    return 0.5 * (r12 + r21), 0.5 * (n12 + n21)


@dataclass
class PeakCentroid:
    ok: bool
    r_max: float = np.nan
    tau_peak: float = np.nan
    tau_cent: float = np.nan
    window: tuple = (np.nan, np.nan)  # lags bounding the r > 0.8 r_max region


def peak_and_centroid(lags, r, threshold: float = THRESHOLD, min_rmax: float = MIN_RMAX) -> PeakCentroid:
    """Peak lag and the centroid of the r > threshold * r_max region (PyCCF rules)."""
    i = int(np.argmax(r))
    r_max = float(r[i])
    if r_max < min_rmax or i == 0 or i == len(r) - 1:
        return PeakCentroid(False, r_max)
    below = r < threshold * r_max
    left = np.nonzero(below[:i])[0]
    right = np.nonzero(below[i + 1:])[0]
    if left.size == 0 or right.size == 0:
        return PeakCentroid(False, r_max, float(lags[i]))
    lo, hi = left[-1], i + 1 + right[0]  # first crossing on each side (excluded)
    sel = slice(lo + 1, hi)
    cent = float(np.sum(r[sel] * lags[sel]) / np.sum(r[sel]))
    return PeakCentroid(True, r_max, float(lags[i]), cent, (float(lags[lo]), float(lags[hi])))


# ---------------------------------------------------------------------------
# Monte Carlo errors (flux randomization + random subset selection)
# ---------------------------------------------------------------------------


def _rss(rng, t, y, e):
    """Random subset selection: N draws with replacement, duplicates merged."""
    n = len(t)
    while True:
        idx = rng.integers(0, n, n)
        uniq, counts = np.unique(idx, return_counts=True)
        if uniq.size > 1:
            return t[uniq], y[uniq], e[uniq] / np.sqrt(counts)


@dataclass
class LagResult:
    """Distributions and summaries from an FR/RSS run."""
    lags: np.ndarray
    r: np.ndarray  # ICCF of the observed light curves
    observed: PeakCentroid
    cccd: np.ndarray  # centroids of successful realizations
    ccpd: np.ndarray  # peaks of successful realizations
    rmax: np.ndarray
    n_mc: int
    extra: dict = field(default_factory=dict)

    @property
    def success(self) -> float:
        return self.cccd.size / self.n_mc

    @staticmethod
    def _summ(x):
        if x.size == 0:
            return (np.nan, np.nan, np.nan)
        lo, med, hi = np.percentile(x, [15.865525, 50.0, 84.134475])
        return float(med), float(med - lo), float(hi - med)

    @property
    def cent(self):
        """(median, minus error, plus error) of the centroid distribution."""
        return self._summ(self.cccd)

    @property
    def peak(self):
        return self._summ(self.ccpd)


def fr_rss(t1, y1, e1, t2, y2, e2, lags, n_mc: int = 1000, seed: int = 1,
           fr: bool = True, rss: bool = True, threshold: float = THRESHOLD,
           min_rmax: float = MIN_RMAX) -> LagResult:
    """ICCF of the data plus ``n_mc`` FR/RSS realizations (Peterson et al. 1998, 2004)."""
    t1, y1, e1, t2, y2, e2 = (np.asarray(a, float) for a in (t1, y1, e1, t2, y2, e2))
    lags = np.asarray(lags, float)
    rng = np.random.default_rng(seed)
    r_obs, _ = ccf(t1, y1, t2, y2, lags)
    obs = peak_and_centroid(lags, r_obs, threshold, min_rmax)
    cents, peaks, rmaxs = [], [], []
    for _ in range(n_mc):
        a = _rss(rng, t1, y1, e1) if rss else (t1, y1, e1)
        b = _rss(rng, t2, y2, e2) if rss else (t2, y2, e2)
        ya = rng.normal(a[1], a[2]) if fr else a[1]
        yb = rng.normal(b[1], b[2]) if fr else b[1]
        r, _ = ccf(a[0], ya, b[0], yb, lags)
        pc = peak_and_centroid(lags, r, threshold, min_rmax)
        if pc.ok:
            cents.append(pc.tau_cent)
            peaks.append(pc.tau_peak)
            rmaxs.append(pc.r_max)
    return LagResult(lags, r_obs, obs, np.array(cents), np.array(peaks), np.array(rmaxs), n_mc)


# ---------------------------------------------------------------------------
# Light-curve helpers
# ---------------------------------------------------------------------------


def detrend(t, y):
    """Subtract a straight-line least-squares fit (Welsh 1999), keeping the mean."""
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    p = np.polyfit(t - t.mean(), y, 1)
    return y - np.polyval(p, t - t.mean()) + y.mean()


def lag_grid(lo: float, hi: float, step: float) -> np.ndarray:
    n = int(round((hi - lo) / step))
    return lo + step * np.arange(n + 1)


# ---------------------------------------------------------------------------
# Alias mitigation for survey light curves (Grier et al. 2017, 2019)
# ---------------------------------------------------------------------------


def seasonal_gaps(t, gap: float):
    """(start, end) pairs of gaps longer than ``gap`` days in a sorted time array."""
    t = np.sort(np.asarray(t, float))
    i = np.nonzero(np.diff(t) > gap)[0]
    return np.column_stack([t[i], t[i + 1]])


def overlap_counts(t_cont, t_line, lags, gap: float):
    """N(tau): line epochs that, shifted back by tau, land inside the continuum
    span and outside its seasonal gaps."""
    t_cont = np.sort(np.asarray(t_cont, float))
    gaps = seasonal_gaps(t_cont, gap)
    shifted = np.asarray(t_line, float)[None, :] - np.asarray(lags, float)[:, None]
    lo = np.maximum(t_cont[0], shifted.min(axis=1))[:, None]
    hi = np.minimum(t_cont[-1], shifted.max(axis=1))[:, None]
    good = (shifted >= lo) & (shifted <= hi)
    for g0, g1 in gaps:
        good &= ~((shifted > g0) & (shifted < g1))
    return good.sum(axis=1)


def alias_weights(t_cont, y_cont, t_line, lags, k: int = 2, gap: float = 30.0):
    """Weight w(tau) = P(tau) convolved with the continuum ACF, normalized to 1.

    P(tau) = [N(tau)/N(0)]^k. The ACF is truncated at its first negative value
    on each side of zero lag, so only the central lobe smears P(tau).
    """
    lags = np.asarray(lags, float)
    n = overlap_counts(t_cont, t_line, lags, gap)
    n0 = overlap_counts(t_cont, t_line, np.array([0.0]), gap)[0]
    p = (n / max(n0, 1)) ** k
    step = lags[1] - lags[0]
    half = max(abs(lags[0]), abs(lags[-1]))
    acf_lags = lag_grid(-half, half, step)
    acf, _ = ccf(t_cont, y_cont, t_cont, y_cont, acf_lags)
    mid = len(acf_lags) // 2
    neg_r = np.nonzero(acf[mid:] < 0)[0]
    neg_l = np.nonzero(acf[:mid][::-1] < 0)[0]
    kern = acf.copy()
    if neg_r.size:
        kern[mid + neg_r[0]:] = 0.0
    if neg_l.size:
        kern[:mid - neg_l[0]] = 0.0
    kern = np.clip(kern, 0.0, None)
    w = np.convolve(p, kern, mode="same") if len(kern) <= len(p) else _conv_centered(p, kern)
    return w / w.max(), n, n0


def _conv_centered(p, kern):
    full = np.convolve(p, kern, mode="full")
    start = (len(kern) - 1) // 2
    return full[start:start + len(p)]


def primary_peak(samples, weights, lags, smooth: float = 15.0):
    """Bounds of the tallest peak of the weighted, smoothed lag distribution.

    The samples are binned on ``lags``, multiplied by ``weights``, smoothed with
    a Gaussian of ``smooth`` days (centred on zero lag), and the peak is bounded
    by the nearest local minima on either side.
    Returns (lo, peak, hi, smoothed weighted histogram).
    """
    lags = np.asarray(lags, float)
    step = lags[1] - lags[0]
    edges = np.concatenate([lags - step / 2, [lags[-1] + step / 2]])
    hist, _ = np.histogram(samples, bins=edges)
    x = np.arange(-int(4 * smooth / step) - 1, int(4 * smooth / step) + 2) * step
    g = np.exp(-0.5 * (x / smooth) ** 2)
    sm = np.convolve(hist * weights, g / g.sum(), mode="same")
    ip = int(np.argmax(sm))
    lo = ip
    while lo > 0 and sm[lo - 1] <= sm[lo]:
        lo -= 1
    hi = ip
    while hi < len(sm) - 1 and sm[hi + 1] <= sm[hi]:
        hi += 1
    return float(lags[lo]), float(lags[ip]), float(lags[hi]), sm


def weighted_lag(res: LagResult, t_cont, y_cont, t_line, k: int = 2, gap: float = 30.0,
                 smooth: float = 15.0):
    """Alias-mitigated centroid lag: median of the CCCD inside the primary peak."""
    w, n, n0 = alias_weights(t_cont, y_cont, t_line, res.lags, k, gap)
    if res.cccd.size == 0:
        return dict(lag=np.nan, err=(np.nan, np.nan), bounds=(np.nan, np.nan), frac_rej=1.0,
                    weights=w, ntau=n, n0=n0)
    lo, pk, hi, sm = primary_peak(res.cccd, w, res.lags, smooth)
    inside = res.cccd[(res.cccd >= lo) & (res.cccd <= hi)]
    if inside.size == 0:
        return dict(lag=np.nan, err=(np.nan, np.nan), bounds=(lo, hi), frac_rej=1.0,
                    weights=w, ntau=n, n0=n0)
    q16, q50, q84 = np.percentile(inside, [15.865525, 50.0, 84.134475])
    # r_max within +/-1 sigma of the lag, from the ICCF of the observed curves
    near = (res.lags >= q16) & (res.lags <= q84)
    rmax_near = float(res.r[near].max()) if near.any() else float(np.interp(q50, res.lags, res.r))
    return dict(lag=float(q50), err=(float(q50 - q16), float(q84 - q50)), bounds=(lo, hi),
                frac_rej=1.0 - inside.size / res.cccd.size, weights=w, ntau=n, n0=n0,
                rmax_near=rmax_near, smoothed=sm)
