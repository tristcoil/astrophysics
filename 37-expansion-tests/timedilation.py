"""Test 1: do distant supernovae unfold in slow motion?

Runs both measurements on the DES 5-year Type Ia supernovae and caches them:

* the replication of White et al. (2024): per-band light-curve widths against
  reference curves (``reference_method``), and the scatter-versus-b curves;
* our own global template fit with the dilation law free (``template_fit``):
  the best b with a bootstrap error, the chi-square of every b from 0 (tired
  light) to 1.6, one dilation factor per redshift bin, and the same fit on
  subsamples and with changed settings.

Heavy steps run in parallel worker processes, one BLAS thread each (the
machine is shared with other long jobs).
"""

from __future__ import annotations

import multiprocessing as mp
import os
import pickle
import time

import numpy as np
import pandas as pd

from config import BANDS, DERIVED_DIR, WHITE24
import fetch_data
import reference_method as rm
import supernovae
import template_fit as tf

N_WORKERS = int(os.environ.get("D7_WORKERS", "6"))
B_DISPERSION = np.linspace(-1.0, 2.0, 30)       # White et al.'s grid
B_SCAN = np.round(np.arange(0.0, 1.61, 0.1), 2)
N_BOOTSTRAP = int(os.environ.get("D7_BOOTSTRAP", "200"))
#: Redshift bins for the model-free dilation factors; ratios are taken
#: relative to the most populated bin.
Z_EDGES = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.15])
Z_REF_BIN = 4  # 0.45 <= z < 0.55

_G = {}  # data shared with forked workers


def _limit_threads():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def _pool(n):
    return mp.get_context("fork").Pool(n, initializer=_limit_threads)


# ---------------------------------------------------------------------------
# White et al. replication
# ---------------------------------------------------------------------------

def _width_task(args):
    cid, band, power, n_mc = args
    seed = int(cid) * 10 + BANDS.index(band) + (5 if power == 0 else 0)
    return rm.fit_width(_G["refset"], cid, band, power=power, n_mc=n_mc, seed=seed)


def _dispersion_task(args):
    cid, band = args
    return cid, band, rm.dispersion_curve(_G["refset"], cid, band, B_DISPERSION)


def white_replication(sn, *, refresh=False):
    """Widths in every band (de-redshifted reference and the null test) + dispersion curves."""
    cache = DERIVED_DIR / "white_replication.pkl"
    if cache.exists() and not refresh:
        return pickle.loads(cache.read_bytes())
    points, _ = supernovae.salt3_view(sn)
    _G["refset"] = rm.ReferenceSet(points, sn)
    ia = sn[sn.is_ia]
    out = {}
    with _pool(N_WORKERS) as pool:
        for label, power, n_mc in (("widths", 1.0, rm.WIDTH_MC_ITER), ("null_test", 0.0, 50)):
            t0 = time.time()
            tasks = [(cid, band, power, n_mc) for cid in ia["cid"] for band in BANDS]
            res = [r for r in pool.map(_width_task, tasks, chunksize=8) if r is not None]
            df = pd.DataFrame(res)
            out[label] = {b: df[df["band"] == b].reset_index(drop=True) for b in BANDS}
            print(f"  {label}: {len(df)} light curves in {time.time() - t0:.0f} s")
        t0 = time.time()
        tasks = [(cid, band) for cid in ia["cid"] for band in BANDS]
        res = pool.map(_dispersion_task, tasks, chunksize=16)
        disp = {b: np.array([d for c, bb, d in res if bb == b and d is not None
                             and np.all(np.isfinite(d))]) for b in BANDS}
        out["dispersion"] = disp
        print(f"  dispersion curves in {time.time() - t0:.0f} s")
    cache.write_bytes(pickle.dumps(out))
    return out


def summarise_white(rep) -> dict:
    """Fits to the replicated widths, compared with White et al.'s own widths."""
    res = {"per_band": {}, "average": {}, "comparison": {}}
    widths = rep["widths"]
    for b in BANDS:
        df = widths[b]
        good = df[df["err"] < df["width"]]
        res["per_band"][b] = {"n": int(len(good)),
                              **{k: v for k, v in rm.fit_power(good.z, good.width, good.err).items()
                                 if k in ("b", "err", "chi2_nu")}}
    # The three ways of finding each width: Nelder-Mead from 1+z (White et
    # al.), Nelder-Mead from 1, and the global minimum on a grid.
    for col in ("width", "width_from_1", "width_grid"):
        avg = rm.average_widths(widths, column=col)
        fit = rm.fit_power(avg.z, avg.width, avg.err)
        lin = rm.fit_linear(avg.z, avg.width, avg.err)
        boot = _bootstrap_power(avg)
        res["average"][col] = {**fit, "bootstrap_err": boot, "linear": lin,
                               "scatter": float(np.std(avg.width - (1 + avg.z) ** fit["b"]))}
    res["averaged_table"] = rm.average_widths(widths)
    # Null test: reference curves NOT divided by (1 + z).
    null = rm.average_widths(rep["null_test"])
    res["null_test"] = {**rm.fit_power(null.z, null.width, null.err),
                        "linear": rm.fit_linear(null.z, null.width, null.err)}
    res["null_table"] = null
    # Their published widths.
    white = {b: fetch_data.load_white_widths(b) for b in BANDS}
    for b in BANDS:
        ours = widths[b].set_index("cid")
        theirs = white[b].set_index("CID")
        common = ours.index.intersection(theirs.index)
        ratio = ours.loc[common, "width"] / theirs.loc[common, "Width"]
        res["comparison"][b] = {"n_ours": int(len(ours)), "n_white": int(len(theirs)),
                                "n_common": int(len(common)),
                                "median_ratio": float(np.median(ratio)),
                                "mad_ratio": float(1.4826 * np.median(np.abs(ratio - np.median(ratio)))),
                                "frac_within_2pct": float(np.mean(np.abs(ratio - 1) < 0.02))}
    wd = {b: white[b].rename(columns={"CID": "cid", "Width": "width", "Width_err": "err"})
          for b in BANDS}
    wavg = rm.average_widths(wd)
    res["white_refit"] = {**rm.fit_power(wavg.z, wavg.width, wavg.err),
                          "bootstrap_err": _bootstrap_power(wavg)}
    res["white_average_table"] = wavg
    # Method 1: dispersion against b.
    disp = rep["dispersion"]
    res["dispersion"] = {}
    for b in BANDS:
        if len(disp[b]) == 0:
            continue
        med = np.median(disp[b], axis=0)
        res["dispersion"][b] = {"n": int(len(disp[b])), "median_curve": med,
                                "b_min": _parabolic_min(B_DISPERSION, med)}
    riz = np.mean([res["dispersion"][b]["median_curve"] for b in "riz"], axis=0)
    res["dispersion"]["riz_mean_curve"] = riz
    res["dispersion"]["riz_b_min"] = _parabolic_min(B_DISPERSION, riz)
    wdisp = {b: fetch_data.load_white_pickle(f"power_dispersion_pickle_{b}") for b in BANDS}
    wriz = np.mean([np.median(wdisp[b], axis=0) for b in "riz"], axis=0)
    res["dispersion"]["white_riz_mean_curve"] = wriz
    res["dispersion"]["white_riz_b_min"] = _parabolic_min(B_DISPERSION, wriz)
    return res


def _parabolic_min(x, y):
    i = int(np.clip(np.argmin(y), 1, len(y) - 2))
    c = np.polyfit(x[i - 1:i + 2], y[i - 1:i + 2], 2)
    return float(-c[1] / (2 * c[0]))


def _bootstrap_power(avg, n=1000, seed=7):
    rng = np.random.default_rng(seed)
    z, w, e = avg.z.to_numpy(), avg.width.to_numpy(), avg.err.to_numpy()
    bs = []
    for _ in range(n):
        k = rng.integers(0, len(z), len(z))
        bs.append(rm.fit_power(z[k], w[k], e[k])["b"])
    return float(np.std(bs))


# ---------------------------------------------------------------------------
# Our global template fit
# ---------------------------------------------------------------------------

def _bootstrap_task(seed):
    F = _G["fitter"]
    rng = np.random.default_rng(seed)
    sn_ids = _G["sn_ids"]
    draw = rng.integers(0, len(sn_ids), len(sn_ids))
    mult = np.bincount(draw, minlength=len(sn_ids))
    weights = mult[_G["sn_index"]]
    F.u = _G["u_best"].copy()
    r = F.fit_power(_G["b_best"], weights=weights)
    return r["b"]


def _variant_task(i):
    label, kwargs, select = _G["variants"][i]
    pts = _G["points"]
    if select is not None:
        pts = pts[select(pts)]
    F = tf.TemplateFitter(pts, **kwargs)
    r = F.fit_power(1.0)
    return label, {"b": r["b"], "b_err_curvature": r["b_err_curvature"],
                   "n_sne": int(pts["cid"].nunique()), "n_points": r["n_points"],
                   "chi2": r["chi2"], "chi2_b0": F.chi2_power(0.0)}


def _variants(sn):
    """The same fit on subsamples and with changed settings."""
    fields = _G["points"].groupby("cid")["field"].agg(lambda s: s.mode().iat[0])
    deep = set(fields[fields.isin(supernovae.DEEP_FIELDS)].index)
    mass = sn.set_index("cid")["host_logmass"]
    prob = sn.set_index("cid")["prob_ia"]
    hi_mass = set(mass[mass >= 10].index)
    sure = set(prob[prob > 0.99].index)
    return [
        ("baseline", {}, None),
        ("error floor 2% of peak", {"eps": 0.02}, None),
        ("error floor 10% of peak", {"eps": 0.10}, None),
        ("knots every 2 days", {"tau_step": 2.0}, None),
        ("knots every 4 days", {"tau_step": 4.0}, None),
        ("rest wavelength > 350 nm", {"lambda_range": (350.0, 860.0)}, None),
        ("window -20 to +60 days", {"fit_window": (-20.0, 60.0)}, None),
        ("no g band", {}, lambda p: p["band"] != "g"),
        ("deep fields only", {}, lambda p: p["cid"].isin(deep)),
        ("shallow fields only", {}, lambda p: ~p["cid"].isin(deep)),
        ("host mass >= 10^10", {}, lambda p: p["cid"].isin(hi_mass)),
        ("host mass < 10^10", {}, lambda p: ~p["cid"].isin(hi_mass)),
        ("z < 0.6", {}, lambda p: p["z"] < 0.6),
        ("P(Ia) > 0.99", {}, lambda p: p["cid"].isin(sure)),
    ]


def dilation_bins(sn, *, refresh=False):
    """One free dilation factor per redshift bin (no power law), ratios to the z ~ 0.5 bin.

    The lowest bin is held fixed while fitting, so every other bin's factor is
    >= 1 and no point is pushed beyond the template's time range; ratios to the
    reference bin and their errors then follow from the fitted covariance.
    Only ratios are measurable: the template's own time axis absorbs an
    overall factor.
    """
    cache = DERIVED_DIR / "dilation_bins.pkl"
    if cache.exists() and not refresh:
        return pickle.loads(cache.read_bytes())
    F = tf.TemplateFitter(supernovae.photometry_view(sn))
    F.fit_power(1.0)
    k = np.digitize(F.z, Z_EDGES) - 1
    nb = len(Z_EDGES) - 1
    phi = np.stack([(k == j).astype(float) for j in range(1, nb)], axis=1)
    zmed = np.array([np.median(F.z[k == j]) for j in range(nb)])
    theta0 = np.log((1 + zmed[1:]) / (1 + zmed[0]))
    r = F.fit(phi=phi, theta0=theta0, free_theta=True)
    theta = np.concatenate([[0.0], r["theta"]])
    cov = np.zeros((nb, nb))
    cov[1:, 1:] = r["theta_cov"]
    ref = Z_REF_BIN
    ln_ratio = theta - theta[ref]
    var = np.diag(cov) + cov[ref, ref] - 2 * cov[:, ref]
    ratio = np.exp(ln_ratio)
    out = {"edges": Z_EDGES, "ref": ref, "D_ratio": ratio,
           "D_ratio_err": ratio * np.sqrt(np.clip(var, 0, None)),
           # Effective (1 + z) of each bin: the mean over its points, so that
           # both sides of the comparison average the same way.
           "one_plus_z": np.array([np.mean(1 + F.z[k == j]) for j in range(nb)]),
           "n_sne": np.array([np.unique(F.cids[k == j]).size for j in range(nb)]),
           "chi2": r["chi2"], "theta": theta, "cov": cov}
    cache.write_bytes(pickle.dumps(out))
    return out


def template_analysis(sn, *, refresh=False):
    cache = DERIVED_DIR / "template_analysis.pkl"
    if cache.exists() and not refresh:
        return pickle.loads(cache.read_bytes())
    t_start = time.time()
    points = supernovae.photometry_view(sn)
    F = tf.TemplateFitter(points)
    out = {"n_points": int(len(F.f)), "n_curves": int(F.n_curves),
           "n_sne": int(np.unique(F.cids).size)}

    best = F.fit_power(1.0)
    out["b"] = best["b"]
    out["b_err_curvature"] = best["b_err_curvature"]
    out["chi2_best"] = best["chi2"]
    u_best = F.u.copy()
    # Keep what the figures need from the best fit.
    out["fit_points"] = F.points.assign(model=best["model"], tau=best["tau"],
                                         amp=best["amp"][F.curve_codes])
    tau_grid = np.linspace(-20, 80, 201)
    lam_grid = np.geomspace(300, 860, 120)
    out["template"] = {"tau": tau_grid, "lam": lam_grid,
                       "S": F.template_grid(tau_grid, lam_grid)}
    print(f"  global fit: b = {out['b']:.4f} +/- {out['b_err_curvature']:.4f} (curvature)")

    out["scan"] = {"b": B_SCAN, "chi2": np.array([F.chi2_power(b) for b in B_SCAN])}
    F.u = u_best.copy()
    print(f"  chi2 scan done ({time.time() - t_start:.0f} s)")

    # Bootstrap over supernovae.
    sn_ids, sn_index = np.unique(F.cids, return_inverse=True)
    _G.update(fitter=F, sn_ids=sn_ids, sn_index=sn_index, u_best=u_best, b_best=out["b"],
              points=points)
    with _pool(N_WORKERS) as pool:
        boot = np.array(pool.map(_bootstrap_task, range(N_BOOTSTRAP)))
    out["bootstrap"] = boot
    out["b_err_bootstrap"] = float(np.std(boot))
    print(f"  bootstrap: sigma_b = {out['b_err_bootstrap']:.4f} from {len(boot)} resamples "
          f"({time.time() - t_start:.0f} s)")

    _G["variants"] = _variants(sn)
    with _pool(N_WORKERS) as pool:
        out["variants"] = dict(pool.map(_variant_task, range(len(_G["variants"]))))
    print(f"  {len(out['variants'])} variants done ({time.time() - t_start:.0f} s)")
    out["runtime_s"] = time.time() - t_start
    cache.write_bytes(pickle.dumps(out))
    return out


def run(*, refresh=False):
    sn = supernovae.supernova_table()
    print("Time dilation: White et al. replication")
    rep = white_replication(sn, refresh=refresh)
    white = summarise_white(rep)
    print("Time dilation: global template fit")
    tmpl = template_analysis(sn, refresh=refresh)
    print("Time dilation: one dilation factor per redshift bin")
    bins = dilation_bins(sn, refresh=refresh)
    # The bins' errors come from the chi-square curvature, which ignores the
    # supernova-to-supernova stretch scatter; scale them by how much the
    # bootstrap widened b's error over its curvature error.
    scale = max(1.0, tmpl["b_err_bootstrap"] / tmpl["b_err_curvature"])
    bins["D_ratio_err_formal"] = bins["D_ratio_err"]
    bins["D_ratio_err"] = bins["D_ratio_err"] * scale
    bins["error_scale"] = scale
    tmpl["bins"] = bins
    print("Time dilation: injection tests of the global fit")
    tmpl["injections"] = injection_tests(sn, refresh=refresh)
    sample = {"n_des": int(len(sn)), "n_ia": int(sn.is_ia.sum()),
              "z_median": float(sn.loc[sn.is_ia, "z"].median()),
              "x1_trend": _x1_trend(sn)}
    return {"sample": sample, "white": white, "template": tmpl, "published": WHITE24}


def _x1_trend(sn):
    """Slope of SALT3 stretch x1 with redshift (a selection-bias diagnostic)."""
    ia = sn[sn.is_ia & np.isfinite(sn.x1)]
    c, cov = np.polyfit(ia.z, ia.x1, 1, cov=True)
    return {"slope": float(c[0]), "slope_err": float(np.sqrt(cov[0, 0]))}


# ---------------------------------------------------------------------------
# Does the global fit recover a known answer? Injection tests.
# ---------------------------------------------------------------------------

INJECTIONS = [(b, seed) for b in (0.0, 0.5, 1.0) for seed in (11, 12, 13)]
STRETCH_SIGMA = 0.08


def _injection_task(args):
    """Simulate the DES data with dilation (1 + z)^b_true and refit it blind.

    The real epochs, bands, redshifts and flux errors are kept; the fluxes are
    replaced by the best-fit template (with the best-fit amplitudes), each
    supernova stretched by its own random factor (8% scatter, like real Type Ia
    diversity), plus Gaussian noise with the real errors.
    """
    b_true, seed = args
    F0, amp = _G["inj_fitter"], _G["inj_amp"]
    rng = np.random.default_rng(seed)
    sn_ids, idx = np.unique(F0.cids, return_inverse=True)
    stretch = 1.0 + STRETCH_SIGMA * rng.standard_normal(sn_ids.size)
    tau = F0.t / ((1.0 + F0.z) ** b_true * stretch[idx])
    flux = amp[F0.curve_codes] * (F0.design(tau) @ F0.coef) + rng.normal(0.0, F0.e)
    inside = (tau >= tf.TAU_RANGE[0]) & (tau <= tf.TAU_RANGE[1])   # no extrapolated template
    pts = F0.points[["cid", "band", "z", "t", "err"]].assign(flux=flux)[inside]
    F = tf.TemplateFitter(pts)
    r = F.fit_power(0.5)
    return {"b_true": b_true, "seed": seed, "b": r["b"], "b_err": r["b_err_curvature"],
            "chi2_b0": F.chi2_power(0.0), "chi2_b1": F.chi2_power(1.0)}


def injection_tests(sn, *, refresh=False):
    cache = DERIVED_DIR / "injection_tests.pkl"
    if cache.exists() and not refresh:
        return pickle.loads(cache.read_bytes())
    t0 = time.time()
    F = tf.TemplateFitter(supernovae.photometry_view(sn))
    best = F.fit_power(1.0)
    _G.update(inj_fitter=F, inj_amp=best["amp"])
    with _pool(N_WORKERS) as pool:
        out = pd.DataFrame(pool.map(_injection_task, INJECTIONS))
    print(f"  {len(out)} injection tests in {time.time() - t0:.0f} s")
    cache.write_bytes(pickle.dumps(out))
    return out
