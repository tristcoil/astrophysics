"""The large frequency separation Dnu, the echelle diagram and the radial modes.

Modes of the same angular degree l and consecutive radial order n are spaced
by almost exactly Dnu, so the oscillation spectrum is a comb. Three steps:

1. **Autocorrelation.** Divide the spectrum by the fitted background (a
   signal-to-noise spectrum), take the part within a few envelope widths of
   numax, and autocorrelate it. The comb makes the autocorrelation peak at a
   lag of Dnu. Red giants also correlate at Dnu/2 (the l = 1 modes sit
   half-way between the l = 0 ones) and at Dnu +- dnu02, so the peak is looked
   for near the value numax predicts (Stello et al. 2009), with the Gaussian
   weighting of Yu et al. (2018). The peak's position is refined with a
   parabola. Its uncertainty comes from repeating the measurement on
   spectra perturbed with chi^2_2 noise (Huber et al. 2009).
2. **Echelle diagram.** Cut the spectrum into Dnu-wide strips and stack them.
   Modes of each degree line up in near-vertical ridges. Where the ridges sit
   (the phase epsilon of the l = 0 modes) is found by matching the stack to the
   red giant "universal pattern" (Mosser et al. 2011, A&A 525, L9): l = 2
   modes just below l = 0, l = 1 half an order away.
3. **Radial modes.** In each order the strongest peak near the l = 0 ridge is
   taken as the radial mode; a straight line (plus the small curvature of
   Mosser et al. 2013) through them gives Dnu from the radial modes alone.
"""

from __future__ import annotations

import numpy as np

import config as C
from spectrum import smooth_boxcar


def expected_dnu(numax):
    a, b = C.DNU_NUMAX_PRIOR
    return a * np.asarray(numax) ** b


def universal_pattern(dnu: float) -> dict:
    """Mosser et al. 2011 (A&A 525, L9, their Table 1), in units of Dnu.

    nu(n, l) = (n + l/2 + eps - d0l + alpha/2 (n - nmax)^2) Dnu
    with eps = 0.634 + 0.546 log10 Dnu, d01 = -0.056 - 0.002 log10 Dnu,
    d02 = 0.131 - 0.033 log10 Dnu, d03 = 0.280; curvature alpha =
    0.015 Dnu^-0.32 (Mosser et al. 2013, A&A 550, A126).
    """
    lg = np.log10(dnu)
    return {"eps": 0.634 + 0.546 * lg, "d01": -0.056 - 0.002 * lg, "d02": 0.131 - 0.033 * lg,
            "d03": 0.280, "alpha": 0.015 * dnu ** -0.32}


def _acf(x: np.ndarray) -> np.ndarray:
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfft(x, n=2 * n)
    a = np.fft.irfft(np.abs(f) ** 2)[:n]
    return a / a[0]


def _parabola_peak(lags, acf, i):
    if i <= 0 or i >= len(acf) - 1:
        return lags[i]
    y0, y1, y2 = acf[i - 1], acf[i], acf[i + 1]
    den = y0 - 2 * y1 + y2
    if den >= 0:
        return lags[i]
    return lags[i] + 0.5 * (y0 - y2) / den * (lags[1] - lags[0])


def acf_window(numax: float, sigma_env: float, nu_lo: float, nu_hi: float):
    dn = float(expected_dnu(numax))
    w = max(3.0 * sigma_env, 4.0 * dn)
    return max(nu_lo, numax - w), min(nu_hi, numax + w)


def _peak_in_acf(lags, acf, dnu_exp, frac, smooth_bins):
    a = smooth_boxcar(acf, smooth_bins)
    weight = np.exp(-0.5 * ((lags - dnu_exp) / (frac * dnu_exp)) ** 2)
    rng_ = (lags > 0.4 * dnu_exp) & (lags < 1.8 * dnu_exp)
    # local maxima within the search range
    i = np.arange(1, len(a) - 1)
    is_max = (a[i] > a[i - 1]) & (a[i] >= a[i + 1]) & rng_[i]
    cand = i[is_max]
    if len(cand) == 0:
        return np.nan, np.nan
    best = cand[np.argmax(a[cand] * weight[cand])]
    return _parabola_peak(lags, a, best), a[best]


def measure_dnu(nu, psd, bg, numax: float, sigma_env: float, *, seed: int = 0,
                n_perturb: int = C.DNU_PERTURBATIONS) -> dict:
    """Dnu from the autocorrelation of the background-divided spectrum."""
    lo, hi = acf_window(numax, sigma_env, nu[0], nu[-1])
    m = (nu >= lo) & (nu <= hi)
    snr = psd[m] / bg[m]
    df = nu[1] - nu[0]
    lags = np.arange(m.sum()) * df
    dnu_exp = float(expected_dnu(numax))
    frac = C.DNU_WINDOW_FRAC_CLUMP if C.DNU_WINDOW_CLUMP_RANGE[0] < numax < C.DNU_WINDOW_CLUMP_RANGE[1] \
        else C.DNU_WINDOW_FRAC
    # smooth the ACF over ~1% of Dnu so that mixed-mode fine structure does not split the peak
    smooth_bins = max(1, int(round(0.01 * dnu_exp / df)))
    acf = _acf(snr)
    dnu, height = _peak_in_acf(lags, acf, dnu_exp, frac, smooth_bins)
    # chi^2_2 perturbations of the observed spectrum (Huber et al. 2009)
    rng = np.random.default_rng(seed)
    pert = []
    near = (lags > 0.85 * dnu) & (lags < 1.15 * dnu)
    for _ in range(n_perturb):
        s = snr * rng.exponential(1.0, size=snr.size)
        a = smooth_boxcar(_acf(s), smooth_bins)
        j = np.flatnonzero(near)
        k = j[np.argmax(a[j])]
        pert.append(_parabola_peak(lags, a, k))
    pert = np.asarray(pert)
    return {"dnu": float(dnu), "dnu_err": float(np.std(pert)), "acf_height": float(height),
            "dnu_expected": dnu_exp, "window": [float(lo), float(hi)], "lags": lags,
            "acf": smooth_boxcar(acf, smooth_bins)}


def _folded(nu, s, dnu, numax, alpha, nbins):
    """Folded (collapsed) spectrum on Dnu, straightening the ridges' curvature."""
    n = nu / dnu
    x = (n - 0.5 * alpha * (n - numax / dnu) ** 2) % 1.0
    h, edges = np.histogram(x, bins=nbins, range=(0, 1), weights=s)
    c, _ = np.histogram(x, bins=nbins, range=(0, 1))
    m = h / np.maximum(c, 1)
    return 0.5 * (edges[1:] + edges[:-1]), m - m.mean()


def _pattern_score(nu, s, dnu, numax, eps_width=0.15, nbins=100):
    """Best match of the folded spectrum to the universal pattern at this Dnu."""
    up = universal_pattern(dnu)
    x, c = _folded(nu, s, dnu, numax, up["alpha"], nbins)
    eps_grid = (up["eps"] + np.arange(-eps_width, eps_width + 1e-9, 0.005)) % 1.0
    scores = [float(np.sum(c * _template(x, e, up))) for e in eps_grid]
    j = int(np.argmax(scores))
    return scores[j], float(eps_grid[j])


def refine_dnu(nu, psd, bg, numax: float, sigma_env: float, dnu0: float, *, seed: int = 0,
               n_perturb: int = 50, span: float = 0.05, step: float = 0.0005,
               eps_width: float = 0.15) -> dict:
    """Refine Dnu and epsilon together by fitting the universal pattern.

    The autocorrelation estimate is a good start but, in red-clump stars
    whose mixed modes are broad and whose spectra are less regular, it is
    easily pulled by 1-2% (Yu et al. 2018, Section 3.2). Here the spectrum
    within two envelope widths of numax is folded on trial values of Dnu
    (within +-5%) and matched to the red-giant universal pattern (Mosser et
    al. 2011): l = 0 and l = 2 ridges, a broad l = 1 ridge half an order away,
    a faint l = 3 ridge. The pattern's phase epsilon is free within +-0.15 of
    Mosser et al.'s epsilon(Dnu) relation, which breaks the degeneracy between
    a slightly wrong Dnu and a shifted epsilon. The uncertainty comes from
    chi^2_2-perturbed spectra, as for the autocorrelation.
    """
    lo, hi = numax - 2.0 * sigma_env, numax + 2.0 * sigma_env
    m = (nu >= lo) & (nu <= hi)
    nu_m = nu[m]
    snr = psd[m] / bg[m]
    df = nu[1] - nu[0]
    w = max(1, int(round(0.01 * dnu0 / df)))
    snr_s = smooth_boxcar(snr, w)
    trials = dnu0 * (1.0 + np.arange(-span, span + 1e-12, step))
    sc = np.array([_pattern_score(nu_m, snr_s, t, numax, eps_width)[0] for t in trials])
    k = int(np.argmax(sc))
    best = float(_parabola_peak(trials, sc, k))
    eps = _pattern_score(nu_m, snr_s, best, numax, eps_width)[1]
    rng = np.random.default_rng(seed + 1)
    pert = []
    near = np.abs(trials - best) < 0.015 * best
    for _ in range(n_perturb):
        s = smooth_boxcar(snr * rng.exponential(1.0, size=snr.size), w)
        ee = np.array([_pattern_score(nu_m, s, t, numax, eps_width)[0] for t in trials[near]])
        j = int(np.argmax(ee))
        pert.append(_parabola_peak(trials[near], ee, j))
    return {"dnu": best, "dnu_err": float(np.std(pert)), "epsilon": eps,
            "trials": trials, "score": sc}


# ---------------------------------------------------------------------------
# Echelle diagram and epsilon
# ---------------------------------------------------------------------------


def collapsed_echelle(nu, snr, dnu, lo, hi, nbins=200):
    """Sum of the signal-to-noise spectrum folded on Dnu (x = nu/Dnu mod 1)."""
    m = (nu >= lo) & (nu <= hi)
    x = (nu[m] / dnu) % 1.0
    h, edges = np.histogram(x, bins=nbins, range=(0, 1), weights=snr[m])
    c, _ = np.histogram(x, bins=nbins, range=(0, 1))
    return 0.5 * (edges[1:] + edges[:-1]), h / np.maximum(c, 1)


def _template(x, eps, up):
    """Universal-pattern template on [0, 1): l=0, l=2, l=1 (broad), l=3."""
    def g(x0, w):
        d = (x - x0 + 0.5) % 1.0 - 0.5
        return np.exp(-0.5 * (d / w) ** 2)
    return (1.0 * g(eps, 0.025) + 0.62 * g(eps - up["d02"], 0.025)
            + 1.0 * g(eps + 0.5 - up["d01"], 0.07) + 0.07 * g(eps + 1.5 - up["d03"], 0.03))


def find_epsilon(nu, snr, dnu, numax, sigma_env) -> dict:
    """Phase of the l = 0 ridge, by matching the folded spectrum to the pattern."""
    up = universal_pattern(dnu)
    lo, hi = numax - 2.0 * sigma_env, numax + 2.0 * sigma_env
    x, ce = collapsed_echelle(nu, snr, dnu, lo, hi)
    ce = ce - np.median(ce)
    grid = np.linspace(0, 1, 400, endpoint=False)
    score = np.array([np.sum(ce * _template(x, e, up)) for e in grid])
    eps = float(grid[np.argmax(score)])
    return {"epsilon": eps, "collapsed_x": x, "collapsed": ce, "score_grid": grid, "score": score,
            "pattern": up}


def radial_modes(nu, snr, dnu, eps, numax, sigma_env, *, min_snr: float = 8.0) -> dict:
    """Strongest peak near the l = 0 ridge in each order, and a fit through them."""
    up = universal_pattern(dnu)
    df = nu[1] - nu[0]
    sm = smooth_boxcar(snr, max(1, int(round(0.004 * dnu / df))))
    n_max = numax / dnu - eps
    orders = np.arange(int(np.floor((numax - 3.5 * sigma_env) / dnu - eps)),
                       int(np.ceil((numax + 3.5 * sigma_env) / dnu - eps)) + 1)
    rows = []
    for n in orders:
        centre = (n + eps + 0.5 * up["alpha"] * (n - n_max) ** 2) * dnu
        m = np.abs(nu - centre) < 0.06 * dnu
        if m.sum() < 3:
            continue
        j = np.flatnonzero(m)[np.argmax(sm[m])]
        if sm[j] >= min_snr:
            rows.append((int(n), float(nu[j]), float(sm[j])))
    out = {"modes": rows}
    if len(rows) >= 4:
        n_arr = np.array([r[0] for r in rows], float)
        f_arr = np.array([r[1] for r in rows])
        # nu = Dnu (n + eps) + (alpha Dnu / 2) (n - n_max)^2
        A = np.vstack([n_arr, np.ones_like(n_arr), (n_arr - n_max) ** 2]).T
        coef, res, *_ = np.linalg.lstsq(A, f_arr, rcond=None)
        resid = f_arr - A @ coef
        dof = max(1, len(f_arr) - 3)
        cov = np.linalg.inv(A.T @ A) * (resid @ resid) / dof
        out.update({"dnu_radial": float(coef[0]), "dnu_radial_err": float(np.sqrt(cov[0, 0])),
                    "eps_radial": float((coef[1] / coef[0]) % 1.0),
                    "alpha": float(2 * coef[2] / coef[0]), "rms_resid": float(np.std(resid))})
    return out


def echelle_image(nu, snr, dnu, lo, hi, *, x_offset: float = 0.0, nx: int = 240):
    """Echelle image: one row per order, ``nx`` columns across one Dnu.

    ``x_offset`` (in units of Dnu) shifts where each strip starts, so the
    ridges are not cut at the edge. Returns (x edges in uHz, row start
    frequencies, image).
    """
    start0 = (np.floor(lo / dnu - x_offset) + x_offset) * dnu
    n_rows = int(np.ceil((hi - start0) / dnu))
    img = np.zeros((n_rows, nx))
    cnt = np.zeros((n_rows, nx))
    m = (nu >= start0) & (nu < start0 + n_rows * dnu)
    rel = nu[m] - start0
    r = np.floor(rel / dnu).astype(int)
    c = np.floor((rel - r * dnu) / dnu * nx).astype(int)
    c = np.clip(c, 0, nx - 1)
    np.add.at(img, (r, c), snr[m])
    np.add.at(cnt, (r, c), 1)
    img = img / np.maximum(cnt, 1)
    rows = start0 + np.arange(n_rows) * dnu
    xedges = np.linspace(0, dnu, nx + 1) + x_offset * dnu
    return xedges, rows, img
