"""Multiple stellar populations: the "chromosome map" of each cluster's red giants.

A globular cluster is not one simple population. Most of its stars ("second
generation", 2G) carry the ashes of hot hydrogen burning: more nitrogen and
sodium, less carbon and oxygen than the "first generation" (1G), which looks
like normal halo stars of the same iron content. The ultraviolet filters were
chosen to see this. F275W contains OH molecular bands and F336W an NH band,
so O-poor/N-rich 2G giants are bluer in m_F275W - m_F814W and lower in the
pseudo-colour C = (m_F275W - m_F336W) - (m_F336W - m_F438W).

Recipe of Milone et al. (2017, Sect. 3), followed step by step:
1. red giants between 1.0 and 3.5 F814W magnitudes above the turnoff, well
   measured in all five bands, members, corrected for differential reddening;
2. in F814W bins, the 4th and 96th percentiles of each colour define a blue and
   a red envelope (boxcar-smoothed over three bins);
3. each star's position between the envelopes, scaled by the RGB width W 2 mag
   above the turnoff, gives Delta_F275W,F814W and Delta_C (their Eqs. 1-2);
4. rotating the map counter-clockwise by 18 degrees makes the 1G sequence
   horizontal (Delta_2 ~ 0); the weight of the lowest Gaussian in a mixture fit
   to the Delta_2 distribution is the 1G fraction (Milone et al. fit a Gaussian
   to the 1G peak alone), with a bootstrap error.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.stats import gaussian_kde

import catalog as K_
import literature as L
from isochrones import Extinction

THETA_DEG = 18.0
RGB_RANGE = (1.0, 3.5)  # F814W mag above the turnoff
WIDTH_LEVEL = 2.0  # W is measured this far above the turnoff
N_BOOT = 100


def _envelopes(mag: np.ndarray, col: np.ndarray, lo: float, hi: float, dm: float = 0.2):
    """4th/96th-percentile envelopes of col in magnitude bins (step dm/3), smoothed."""
    centres = np.arange(lo + dm / 2, hi - dm / 2 + 1e-9, dm / 3)
    blue, red, mc = [], [], []
    for cen in centres:
        s = abs(mag - cen) < dm / 2
        if s.sum() < 12:
            continue
        blue.append(np.percentile(col[s], 4))
        red.append(np.percentile(col[s], 96))
        mc.append(np.mean(mag[s]))
    blue = uniform_filter1d(np.array(blue), 3, mode="nearest")
    red = uniform_filter1d(np.array(red), 3, mode="nearest")
    return np.array(mc), blue, red


def chromosome_map(prep: dict, inputs: dict, rng: np.random.Generator) -> dict:
    cat = K_.read_hugs(prep["meta"]["hugs"])
    idx = prep["sel_index"]
    uv_ok = K_.quality_mask(cat, ("F275W", "F336W", "F438W"))[idx]

    # differential reddening in the UV bands, with ratios for a red giant
    ext_acs = Extinction(inputs["feh"], 0.4, "HST_ACS_WFC")
    ext_uv = Extinction(inputs["feh"], 0.4, "HST_WFC3")
    lt, lg = 3.70, 2.5
    a = {b: ext_uv.ratio(b, lt, lg) for b in ("F275W", "F336W", "F438W")}
    a["F814W"] = ext_acs.ratio("F814W", lt, lg)
    dE = prep["dE"]
    m275 = prep["m_F275W_raw"] - a["F275W"] * dE
    m336 = prep["m_F336W_raw"] - a["F336W"] * dE
    m438 = prep["m_F438W_raw"] - a["F438W"] * dE
    m814 = prep["m814"]  # already corrected

    ridge = prep["ridge"]
    m814_to = ridge.m_to - ridge.c_to
    lo, hi = m814_to - RGB_RANGE[1], m814_to - RGB_RANGE[0]
    x_col = m275 - m814
    c_idx = (m275 - m336) - (m336 - m438)
    # giants: near the F606W - F814W ridge on the post-turnoff side
    from fiducial import K, project

    post = np.arange(len(ridge.m)) > ridge.i_to
    s_r, t_r = project(K * ridge.c[post], ridge.m[post], K * prep["c"], prep["m606"])
    rgb = (uv_ok & np.isfinite(m275) & np.isfinite(m336) & np.isfinite(m438)
           & (m814 > lo) & (m814 < hi) & (abs(t_r) < 0.25))
    mag, X, Y = m814[rgb], x_col[rgb], c_idx[rgb]

    mx, xb, xr = _envelopes(mag, X, lo, hi)
    my, yb, yr = _envelopes(mag, Y, lo, hi)
    xb_s, xr_s = np.interp(mag, mx, xb), np.interp(mag, mx, xr)
    yb_s, yr_s = np.interp(mag, my, yb), np.interp(mag, my, yr)
    level = m814_to - WIDTH_LEVEL
    w_x = float(np.interp(level, mx, xr) - np.interp(level, mx, xb))
    w_c = float(np.interp(level, my, yr) - np.interp(level, my, yb))
    d_x = w_x * (X - xr_s) / (xr_s - xb_s)
    d_c = w_c * (yr_s - Y) / (yr_s - yb_s)

    th = np.radians(THETA_DEG)
    d1 = d_x * np.cos(th) - d_c * np.sin(th)
    d2 = d_x * np.sin(th) + d_c * np.cos(th)
    f1g, peak, sigma, is_1g = first_generation_fraction(d2, return_labels=True)
    boot = []
    for _ in range(N_BOOT):
        k = rng.integers(0, len(d2), len(d2))
        boot.append(first_generation_fraction(d2[k])[0])
    lit = L.milone2017().get(prep["meta"]["hugs"], {})
    return {
        "n_rgb": int(rgb.sum()), "W_275_814": w_x, "W_C": w_c,
        "f1G": f1g, "f1G_err": float(np.std(boot)), "d2_peak": peak, "d2_sigma": sigma,
        "extinction_ratios": a,
        "milone2017": {k: list(v) if isinstance(v, tuple) else v for k, v in lit.items()},
    }, {"d_x": d_x, "d_c": d_c, "d1": d1, "d2": d2, "mag": mag, "X": X, "Y": Y,
        "env_x": (mx, xb, xr), "env_y": (my, yb, yr), "level": level, "is_1g": is_1g}


def first_generation_fraction(d2: np.ndarray, max_components: int = 4, return_labels: bool = False):
    """Weight of the lowest-Delta_2 Gaussian in a 1D mixture fit (components chosen by BIC).

    Milone et al. fit a Gaussian to the 1G peak of the Delta_2 histogram and take
    its area over the total; a mixture fit does the same with the 2G stars
    modelled too, instead of assuming the 1G peak is uncontaminated. With
    ``return_labels`` each star is also assigned to its most probable component
    (True = first generation), consistent with the fraction.
    """
    from sklearn.mixture import GaussianMixture
    from threadpoolctl import threadpool_limits

    x = d2.reshape(-1, 1)
    best, best_bic = None, np.inf
    # One thread, set *after* scikit-learn has loaded its OpenMP runtime: the
    # k-means start of each fit otherwise spawns a thread per core, which on a
    # busy machine made the bootstrap ~30x slower.
    with threadpool_limits(limits=1):
        for k in range(2, max_components + 1):
            if len(d2) < 10 * k:
                break
            gm = GaussianMixture(k, n_init=3, random_state=0).fit(x)
            bic = gm.bic(x)
            if bic < best_bic:
                best, best_bic = gm, bic
    if best is None:
        out = (np.nan, np.nan, np.nan)
        return (*out, np.zeros(len(d2), bool)) if return_labels else out
    i = int(np.argmin(best.means_[:, 0]))
    out = (float(best.weights_[i]), float(best.means_[i, 0]), float(np.sqrt(best.covariances_[i].ravel()[0])))
    if return_labels:
        return (*out, best.predict(x) == i)
    return out


def run_all(prepared: dict, inputs: dict) -> tuple[dict, dict]:
    rng = np.random.default_rng(7)
    res, art = {}, {}
    for name, prep in prepared.items():
        res[name], art[name] = chromosome_map(prep, inputs[name], rng)
    return res, art
