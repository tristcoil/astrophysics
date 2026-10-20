"""The whole analysis, in order; returns everything the figures and README need."""

from __future__ import annotations

import time
import warnings

import numpy as np

warnings.filterwarnings("ignore", module="astropy")

import catalog as K  # noqa: E402
import cmb_link  # noqa: E402
import config as C  # noqa: E402
import cosmology as CO  # noqa: E402
import dark_energy as DE  # noqa: E402
import footprint as F  # noqa: E402
import inference as I  # noqa: E402
import literature as L  # noqa: E402
import masses as MS  # noqa: E402
import model as MD  # noqa: E402
import selection as S  # noqa: E402


def _log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def survey_section(fp) -> dict:
    clean, total = F.digitize_fig2b()
    out = {"area_primary": float(fp["area_primary"]), "area_cosmo": float(fp["area_cosmo"]),
           "area_base": float(fp["area_base"]), "radius_px": int(fp["radius_px"]),
           "trial_radii": fp["trial_radii"].tolist(), "trial_areas": fp["trial_areas"].tolist()}
    checks = []
    for name, mask, curve in (("clean", fp["mask_primary"], clean), ("total", fp["mask_de"], total)):
        e = np.sort(fp["exposure"][mask & np.isfinite(fp["exposure"])])[::-1].astype(float)
        cum = np.arange(1, e.size + 1) * F.PIX_AREA
        for tq in (80, 100, 150, 200, 300, 500, 1000):
            j = np.argmin(np.abs(np.log(curve[:, 0] / tq)))
            ours = float(cum[np.searchsorted(-e, -tq) - 1])
            checks.append({"curve": name, "T": tq, "ours": ours, "bulbul": float(curve[j, 1]),
                           "ratio": ours / float(curve[j, 1])})
        out[f"cum_{name}"] = (e[:: max(1, e.size // 3000)], cum[:: max(1, e.size // 3000)])
    out["fig2b_checks"] = checks
    out["fig2b_clean"], out["fig2b_total"] = clean, total
    m = fp["mask_cosmo"]
    for key in ("exposure", "nh", "zvlim"):
        v = np.asarray(fp[key], float)[m]
        out[f"{key}_pct"] = np.nanpercentile(v, [1, 5, 25, 50, 75, 95, 99]).tolist()
    return out


#: Global mass scale (M_true / M_assumed) for the sensitivity test: 1.0 is the
#: lensing calibration; 0.58-0.8 span what hydrostatic masses would imply.
MASS_SCALES = (0.58, 0.69, 0.76, 0.85, 1.0, 1.15)


def run(processes: int | None = None, mcmc_steps: int = C.MCMC_STEPS, skip_mcmc: bool = False) -> dict:
    processes = I.default_processes() if processes is None else processes
    res, summ = {}, {}

    _log("footprint")
    fp = F.build()
    res["survey"] = survey_section(fp)
    classes = F.survey_classes(fp)

    _log("catalogues")
    df = K.load_cosmology_sample(fp)
    primary = K.load_primary()
    psz2 = K.load_psz2()
    res["sample"] = {"n_published": df.attrs["n_published"], "n_used": len(df),
                     "n_zero_rate": df.attrs["n_zero_rate"], "n_missing_primary": df.attrs["n_missing_primary"],
                     "sum_weights": float(df.weight.sum()), "purity": float(df.weight.mean()),
                     "in_mask_frac": float(df.in_mask.mean()),
                     "z_median": float(df.z.median()), "cr_median": float(df.cr.median()),
                     "exp_median": float(df.exp.median())}
    res["meas_error"] = K.measurement_error_fit(df)
    res["absorption"] = K.absorption_fit(df)

    _log("selection function")
    sel = S.SelectionModel()
    res["selection"] = {"nodes": sel.nodes["params"].tolist(), "g": sel.nodes["g"], "rms": sel.nodes["rms"],
                        "validate": S.validate(sel),
                        "validate_no_clerc": S.validate(S.SelectionModel(S.fit_nodes(use_clerc=False))),
                        "cr50": {f"{z}": {f"{t}": float(sel.cr50(z, t)[0]) for t in (50, 150, 500, 2000)}
                                 for z in (0.1, 0.2, 0.3, 0.5, 0.65, 0.8)}}
    res["fig4_curves"] = S.extract_fig4_curves()
    res["clerc_points"] = S.digitize_clerc()

    _log("counts model at the published best fit")
    kernel = MD.build_kernel(sel, classes)
    cm = MD.CountsModel(kernel)
    data = MD.data_counts(df)
    data_unw = MD.data_counts(df, weighted=False)
    pub = CO.fiducial()
    N_pub = cm.expected(pub)
    res["counts"] = {"data": data, "data_unweighted": data_unw, "pred_published": N_pub,
                     "cr_edges": cm.cr_edges, "z_edges": cm.z_edges}
    summ["published_point"] = {"N_pred": float(N_pub.sum()), "N_obs_weighted": float(data.sum()),
                               "N_obs": int(len(df))}

    _log("consistency scan (published scaling relation)")
    grid = I.grid_scan(I.BASELINE, data, processes=processes)
    res["grid"] = grid
    best = CO.Cosmology(Om=grid["best"]["Om"], sigma8=grid["best"]["sigma8"])
    N_best = cm.expected(best)
    res["counts"]["pred_best"] = N_best
    res["counts"]["gof"] = goodness_of_fit(data, N_best)
    res["counts"]["gof_published"] = goodness_of_fit(data, N_pub)

    _log("how clumpy: predictions for three sigma_8 values")
    res["clumpy"] = {}
    for s8 in (0.76, 0.81, 0.86):
        om = 0.3
        res["clumpy"][s8] = cm.expected(CO.Cosmology(Om=om, sigma8=s8))

    if not skip_mcmc:
        _log(f"MCMC with the scaling relation free ({mcmc_steps} steps, {processes} processes)")
        res["mcmc"] = I.run_mcmc(I.BASELINE, data, nsteps=mcmc_steps,
                                 burn=min(C.MCMC_BURN, mcmc_steps // 3), processes=processes)

    if "mcmc" in res:
        med = res["mcmc"]["summary"]
        cos_m = CO.Cosmology(Om=med["Om"][0], sigma8=med["sigma8"][0])
        sr_m = MD.Scaling(*(med[k][0] for k in ("A_X", "B_X", "F_X", "G_X", "sigma_X")))
        res["counts"]["pred_mcmc"] = cm.expected(cos_m, sr_m)
        res["counts"]["gof_mcmc"] = goodness_of_fit(data, res["counts"]["pred_mcmc"])

    _log("robustness variants")
    res["variants"] = I.run_variants(df, processes=processes)

    _log("cluster masses and Planck's hydrostatic bias")
    lnm, P = MS.cluster_mass_posteriors(df, pub)
    med, lo, hi = MS.summarize_posteriors(lnm, P)
    res["masses"] = {"lnm": lnm, "sum_pdf": P.sum(0), "median": med, "lo": lo, "hi": hi,
                     "z": df.z.values, "vs_catalogue": MS.compare_with_catalogue(df, med),
                     "cat_m500": df.m500.values}
    hb = MS.hydrostatic_bias(primary, psz2, df=df, cosmo=pub, our_masses=med)
    res["hydrostatic"] = hb

    _log("mass-scale sensitivity: what if the masses were biased?")
    res["mass_scale"] = {}
    for mscale in MASS_SCALES:
        spec = dict(I.BASELINE, mass_scale=mscale, label=f"mass scale {mscale:.2f}")
        res["mass_scale"][mscale] = I.grid_scan(spec, data, np.linspace(0.12, 0.48, 37),
                                                np.linspace(0.52, 1.12, 61), processes)

    _log("D3's CMB, grown to today (CAMB)")
    s = cmb_link.sample(processes=processes)
    res["d3"] = {"samples": s, "summary": cmb_link.summary(s), "source": cmb_link.d3_fit()["source"],
                 "planck_check_sigma8": cmb_link.planck_check()}

    _log("dark energy's brake")
    res["brake"] = DE.brake()

    summ.update(build_summary(res))
    res["summary"] = summ
    return res


def goodness_of_fit(data, model) -> dict:
    """Pearson chi^2 of the marginal distributions and of the 2D bins with N > 5."""
    out = {}
    for name, ax in (("z", 0), ("cr", 1)):
        d, m = data.sum(ax), model.sum(ax)
        ok = m > 5
        out[f"chi2_{name}"] = float(np.sum((d[ok] - m[ok]) ** 2 / m[ok]))
        out[f"nbins_{name}"] = int(ok.sum())
    ok = model > 5
    out["chi2_2d"] = float(np.sum((data[ok] - model[ok]) ** 2 / model[ok]))
    out["nbins_2d"] = int(ok.sum())
    return out


def build_summary(res) -> dict:
    g = res["grid"]
    s = {"fixed_scaling": {k: g[k] for k in ("best", "Om", "sigma8", "S8")}}
    if "mcmc" in res:
        s["mcmc"] = res["mcmc"]["summary"]
        s["mcmc_meta"] = {k: res["mcmc"][k] for k in ("acceptance", "autocorr", "seconds")}
    s["variants"] = [{"label": v["label"], "Om": v["Om"], "sigma8": v["sigma8"], "S8": v["S8"]}
                     for v in res["variants"]]
    s["mass_scale"] = {f"{k:.2f}": {"S8": v["S8"], "sigma8": v["sigma8"], "Om": v["Om"],
                                    "touches_edge": v["touches_edge"]} for k, v in res["mass_scale"].items()}
    s["hydrostatic"] = {"table": res["hydrostatic"]["table"], "n_matched": res["hydrostatic"]["n_matched"],
                        "slope": res["hydrostatic"]["slope_dlog_ratio_dlogM"]}
    s["d3"] = res["d3"]["summary"]
    s["d3"]["planck_check_sigma8"] = res["d3"]["planck_check_sigma8"]
    b = res["brake"]
    s["brake"] = {k: b[k] for k in ("ratio_since_cmb", "sigma8_lcdm", "sigma8_noL_phys", "h_noL",
                                    "boost_gt_1e+14", "boost_gt_1e+15")}
    s["masses_vs_catalogue"] = res["masses"]["vs_catalogue"]
    s["gof"] = res["counts"]["gof"]
    s["gof_published"] = res["counts"]["gof_published"]
    if "gof_mcmc" in res["counts"]:
        s["gof_mcmc"] = res["counts"]["gof_mcmc"]
    s["survey"] = {k: res["survey"][k] for k in ("area_primary", "area_cosmo", "radius_px", "fig2b_checks",
                                                 "exposure_pct", "nh_pct", "zvlim_pct")}
    s["sample"] = res["sample"]
    s["meas_error"] = {k: res["meas_error"][k] for k in ("a", "poisson_ratio", "n_used")}
    s["absorption"] = res["absorption"]
    s["selection"] = {"g": res["selection"]["g"], "validate": res["selection"]["validate"],
                      "validate_no_clerc": res["selection"]["validate_no_clerc"], "cr50": res["selection"]["cr50"]}
    # tensions with Planck 2018 (Gaussian, on S8)
    p = L.PLANCK2018["S8"]
    for key, val in (("fixed_scaling", s["fixed_scaling"]["S8"]), ("mcmc", s.get("mcmc", {}).get("S8"))):
        if val is None:
            continue
        med, lo, hi = val
        sig = (hi - lo) / 2
        s[f"S8_vs_planck_{key}"] = float((med - p[0]) / np.hypot(sig, p[1]))
    return s
