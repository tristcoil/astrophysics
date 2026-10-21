"""Finding the length of the day from the dot's flicker.

Three tools: a Lomb-Scargle periodogram (unevenly sampled data), an
autocorrelation on an hourly grid that tolerates gaps (the method Jiang et al.
2018 used), and Jiang et al.'s experiment of thinning the data to ask how
rarely one could observe and still find the rotation.
"""

from __future__ import annotations

import numpy as np
from astropy.timeseries import LombScargle

SOLAR_DAY_H = 24.0
#: Sidereal day: 86164.0905 s.
SIDEREAL_DAY_H = 86164.0905 / 3600.0


def periodogram(t_d: np.ndarray, y: np.ndarray, f_min: float, f_max: float,
                oversample: float = 10.0) -> tuple[np.ndarray, np.ndarray]:
    """Lomb-Scargle power on a grid fine enough to resolve the peaks."""
    span = t_d.max() - t_d.min()
    n = int(np.ceil((f_max - f_min) * span * oversample)) + 1
    f = np.linspace(f_min, f_max, n)
    return f, LombScargle(t_d, y).power(f)


def peak_near(f: np.ndarray, p: np.ndarray, f0: float, half_width: float) -> tuple[float, float]:
    """Frequency and power of the highest peak within f0 +- half_width,
    refined by a parabola through the three highest grid points."""
    sel = np.nonzero((f > f0 - half_width) & (f < f0 + half_width))[0]
    i = sel[np.argmax(p[sel])]
    if 0 < i < len(f) - 1:
        y0, y1, y2 = p[i - 1], p[i], p[i + 1]
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
        return float(f[i] + shift * (f[1] - f[0])), float(y1)
    return float(f[i]), float(p[i])


def bootstrap_peak(t_d: np.ndarray, y: np.ndarray, f0: float, half_width: float,
                   rng: np.random.Generator, n_boot: int = 200, oversample: float = 20.0) -> dict:
    """Peak frequency with an uncertainty from resampling whole days.

    Weather makes neighbouring images within a day correlated, so days (not
    images) are the independent units."""
    f_grid = np.linspace(f0 - half_width, f0 + half_width,
                         int(2 * half_width * (t_d.max() - t_d.min()) * oversample) + 1)
    ls = LombScargle(t_d, y)
    f_best, _ = peak_near(f_grid, ls.power(f_grid), f0, half_width)
    day = np.floor(t_d).astype(int)
    days = np.unique(day)
    groups = {d: np.nonzero(day == d)[0] for d in days}
    peaks = []
    for _ in range(n_boot):
        pick = rng.choice(days, size=len(days), replace=True)
        idx = np.concatenate([groups[d] for d in pick])
        # Keep each resampled day at its own time: the phases stay coherent.
        fb, _ = peak_near(f_grid, LombScargle(t_d[idx], y[idx]).power(f_grid), f0, half_width)
        peaks.append(fb)
    peaks = np.array(peaks)
    return {"f": f_best, "f_err": float(np.std(peaks)), "period_h": 24.0 / f_best,
            "period_err_h": float(24.0 / f_best**2 * np.std(peaks))}


def hourly_grid(t_d: np.ndarray, y: np.ndarray, bin_h: float = 1.0):
    """Average onto a regular grid; empty bins are NaN."""
    edges = np.arange(t_d.min(), t_d.max() + bin_h / 24, bin_h / 24)
    idx = np.digitize(t_d, edges) - 1
    s = np.bincount(idx, weights=y, minlength=len(edges))
    n = np.bincount(idx, minlength=len(edges))
    with np.errstate(invalid="ignore"):
        g = s / n
    g[n == 0] = np.nan
    return edges[: len(g)], g


def acf(grid: np.ndarray, max_lag: int) -> np.ndarray:
    """Autocorrelation of a gappy series: each lag uses all pairs where both ends exist."""
    x = grid - np.nanmean(grid)
    out = np.full(max_lag + 1, np.nan)
    for k in range(max_lag + 1):
        a, b = x[: len(x) - k], x[k:]
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() > 10:
            out[k] = np.sum(a[m] * b[m]) / np.sqrt(np.sum(a[m] ** 2) * np.sum(b[m] ** 2))
    return out


def acf_peak(lags_h: np.ndarray, r: np.ndarray, lo: float = 16.0, hi: float = 32.0) -> tuple[float, float]:
    """First ACF maximum between lo and hi hours (parabolic refinement)."""
    sel = np.nonzero((lags_h >= lo) & (lags_h <= hi))[0]
    i = sel[np.nanargmax(r[sel])]
    y0, y1, y2 = r[i - 1], r[i], r[i + 1]
    denom = y0 - 2 * y1 + y2
    shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    return float(lags_h[i] + shift * (lags_h[1] - lags_h[0])), float(y1)


def jackknife_peak(t_d: np.ndarray, y: np.ndarray, groups: np.ndarray, f0: float,
                   half_width: float, oversample: float = 20.0) -> dict:
    """Peak frequency with a delete-one-group jackknife error (groups = months).

    Clouds and the seasonal shift of the land/ocean pattern correlate whole
    weeks, so a bootstrap over single days understates the error; leaving out
    one month at a time does not."""
    f_grid = np.linspace(f0 - half_width, f0 + half_width,
                         int(2 * half_width * (t_d.max() - t_d.min()) * oversample) + 1)
    f_best, _ = peak_near(f_grid, LombScargle(t_d, y).power(f_grid), f0, half_width)
    labels = np.unique(groups)
    est = np.array([peak_near(f_grid, LombScargle(t_d[groups != g], y[groups != g]).power(f_grid),
                              f0, half_width)[0] for g in labels])
    n = len(labels)
    f_err = float(np.sqrt((n - 1) / n * np.sum((est - est.mean()) ** 2)))
    return {"f": f_best, "f_err": f_err, "period_h": 24.0 / f_best,
            "period_err_h": 24.0 / f_best**2 * f_err, "n_groups": n}


def _thin_regular(t_d, start, seg_days, dt, rng):
    phase = rng.uniform(0, dt)
    targets = np.arange(start + phase, start + seg_days, dt)
    j = np.clip(np.searchsorted(t_d, targets), 1, len(t_d) - 1)
    nearest = np.where(np.abs(t_d[j - 1] - targets) < np.abs(t_d[j] - targets), j - 1, j)
    good = np.abs(t_d[nearest] - targets) < min(dt / 2, 1.0 / 24)
    return np.unique(nearest[good])


def _thin_random(t_d, start, seg_days, dt, rng):
    inside = np.nonzero((t_d >= start) & (t_d < start + seg_days))[0]
    n = int(round(seg_days / dt))
    if len(inside) < n:
        return inside
    return np.sort(rng.choice(inside, size=n, replace=False))


def sampling_test(t_d: np.ndarray, y: np.ndarray, base_interval_s: float, multiples,
                  rng: np.random.Generator, *, mode: str = "random", n_segments: int = 100,
                  seg_days: float = 7.0, band=(0.95, 1.05), f_max: float = 6.0) -> list[dict]:
    """Jiang et al. (2018) Table 2, redone: keep on average one image every
    k x base_interval in random week-long segments and count how often the
    strongest periodogram peak lies at the rotation frequency (0.95-1.05 /day).

    mode="random": the kept images are a random subset of that size (what
    Jiang et al. describe: "we randomly create a sampling"); irregular times
    do not alias. mode="regular": one image per interval, the nearest real one
    to an evenly spaced grid, as a telescope on a fixed schedule would get;
    below two samples per day the rotation then aliases to other frequencies.
    The search runs from one cycle per segment to ``f_max``."""
    thin = _thin_random if mode == "random" else _thin_regular
    out = []
    t0, t1 = t_d.min(), t_d.max() - seg_days
    f = np.linspace(1.0 / seg_days, f_max, 6000)
    for k in multiples:
        dt = k * base_interval_s / 86400.0
        hits, used = 0, 0
        for _ in range(n_segments * 3):  # extra tries for segments that land in gaps
            if used == n_segments:
                break
            idx = thin(t_d, rng.uniform(t0, t1), seg_days, dt, rng)
            if len(idx) < max(4, 0.6 * seg_days / dt):
                continue
            used += 1
            fbest = f[np.argmax(LombScargle(t_d[idx], y[idx]).power(f))]
            hits += band[0] <= fbest <= band[1]
        out.append({"interval_h": dt * 24, "per_day": 1 / dt, "n_segments": used,
                    "p_detect": hits / used if used else np.nan, "mode": mode})
    return out
