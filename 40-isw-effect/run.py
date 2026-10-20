"""Run the whole D10 analysis from cached data and write data/results.json.

    python fetch_data.py   # once: downloads (~0.4 GB new; Planck maps linked from D3)
    python run.py          # analysis (first run ~1.5-2 h, mostly simulations; then minutes)
    python figures.py      # every figure from the cached results

Steps, in order:
 1. CMB maps (Planck SMICA, WMAP ILC) at the working resolutions.
 2. The filter bank and 1,000 simulated skies through it.
 3. BOSS: footprint, 3D grid, bias, Wiener filter, predicted ISW skies (+ variants).
 4. Unblinded stacks (all choices were frozen in data/preregistration.json first):
    Granett et al. replication, Nadathur & Hotchkiss, Mao et al.; matched-filter
    amplitudes; look-elsewhere test; random-position errors; stacked images.
 5. Cross-correlations with NVSS and BOSS; the BOSS template fit; rotated maps.
 6. Joint amplitude from all probes, with their covariance from shared simulations.
 7. The Cold Spot and the Eridanus supervoid.
"""

from __future__ import annotations

import json
import time
from datetime import datetime, timezone

import numpy as np
from scipy.stats import norm

import boss
import catalogs
import cmb
import coldspot as CS
import crosscorr as X
import simpass
import stacking as S
import superstructures as SS
from config import (LITERATURE, N_SIMS, PREREG_PATH, RESULTS_PATH, SCAN_ALPHA, SCAN_RADII_DEG,
                    SMHW_SCALE_DEG)

TAGS = X.TEMPLATE_TAGS


def _log(msg, t0):
    print(f"[{time.monotonic() - t0:7.0f} s] {msg}", flush=True)


def _clean(o):
    """Make numpy objects JSON-serialisable (arrays become lists)."""
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


def stacking_results(t0) -> tuple[dict, dict]:
    layout = SS.build_layout()
    sims = simpass.bank_sims()
    data = {w: SS.real_values(w, layout) for w in ("smica", "wmap")}
    preds = {t: SS.predicted_values(t, layout) for t in TAGS}
    res: dict = {"stacks": {}, "matched_filter": {}, "scans": {}}
    for key, setting in SS.PRIMARY.items():
        cat = layout.catalogs[key]
        rows = layout.rows(key, setting)
        entry = {"n": len(cat), "n_valid": int(layout.bank.valid[rows].sum()),
                 "setting": list(setting), "predicted_uK": {t: float(SS.stack(preds[t], layout, key,
                                                                               setting))
                                                           for t in TAGS}}
        for w, v in data.items():
            entry[w] = S.summarize(SS.stack(v, layout, key, setting), SS.stack(sims, layout, key,
                                                                                setting), cat.sign)
        res["stacks"][key] = entry
        # exploratory scans
        settings = SS.settings_for(key)
        res["scans"][key] = {
            "x": [s[1] for s in settings], "kind": settings[0][0],
            "data": [float(SS.stack(data["smica"], layout, key, s)) for s in settings],
            "sigma": [float(np.std(SS.stack(sims, layout, key, s), ddof=1)) for s in settings],
            "pred": [float(SS.stack(preds["fiducial"], layout, key, s)) for s in settings]}
    comb = {w: S.summarize(SS.granett_combined(v, layout, 4.0),
                           SS.granett_combined(sims, layout, 4.0), 1) for w, v in data.items()}
    comb["predicted_uK"] = {t: float(SS.granett_combined(preds[t], layout, 4.0)) for t in TAGS}
    res["stacks"]["granett_combined"] = comb
    res["scans"]["granett_combined"] = {
        "x": list(SCAN_RADII_DEG),
        "data": [float(SS.granett_combined(data["smica"], layout, r)) for r in SCAN_RADII_DEG],
        "data_wmap": [float(SS.granett_combined(data["wmap"], layout, r)) for r in SCAN_RADII_DEG],
        "sigma": [float(np.std(SS.granett_combined(sims, layout, r), ddof=1))
                  for r in SCAN_RADII_DEG],
        "pred": [float(SS.granett_combined(preds["fiducial"], layout, r)) for r in SCAN_RADII_DEG]}
    _log("stacks measured", t0)

    groups = {"granett": list(SS.GRANETT_KEYS), "nh14": list(SS.NH14_KEYS), "mao": ["mao_voids"]}
    a_sims = {}
    for g, keys in groups.items():
        rows = np.concatenate([layout.rows(k, SS.PRIMARY[k]) for k in keys])
        out = {}
        for t in TAGS:
            r = SS.matched_filter(data["smica"], preds[t], sims, rows)
            if t == "fiducial":
                a_sims[g] = r["A_sims"]
            out[t] = {k: v for k, v in r.items() if k != "A_sims"}
        out["wmap"] = {k: v for k, v in SS.matched_filter(data["wmap"], preds["fiducial"], sims,
                                                           rows).items() if k != "A_sims"}
        res["matched_filter"][g] = out
    _log("matched filters done", t0)

    le = SS.look_elsewhere(data["smica"], sims, layout)
    res["look_elsewhere"] = {k: v for k, v in le.items() if k != "max_snr_sims"}
    res["look_elsewhere"]["max_snr_sims_percentiles"] = {
        p: float(np.percentile(le["max_snr_sims"], p)) for p in (50, 84, 95, 99)}
    res["look_elsewhere"]["global_sigma"] = float(norm.isf(le["p_global"]))
    res["look_elsewhere"]["fixed_sigma_gaussian"] = float(le["fixed_snr"])
    np.save(SS.DERIVED_DIR / "lookelsewhere_max_snr_sims.npy", le["max_snr_sims"])

    res["random_positions"] = {}
    for key in SS.PRIMARY:
        r = SS.random_position_errors(key, SS.PRIMARY[key], layout)
        res["random_positions"][key] = {k: v for k, v in r.items() if k != "set_means"}
    _log("look-elsewhere and random positions done", t0)
    SS.stacked_images()
    return res, a_sims


def crosscorr_results(t0) -> tuple[dict, dict]:
    models = X.models()
    tr = X._prepared_tracers()
    tf = X.template_filters()
    d = cmb.prepare_maps()
    data = X.measure(d["smica64"], d["mask64"].astype(float), tr, tf)
    data_wmap = X.measure(d["wmap64"], d["mask64"].astype(float), tr, tf)
    sims = X.run_cc_sims(N_SIMS)
    rot = X.rotated_measurements()
    res: dict = {"tracers": {}, "template": {}}
    a_sims = {}
    for name in tr:
        m = models[name]
        fit = X.fit_amplitude(data[name], m["gT_model"], sims[name])
        fit_w = X.fit_amplitude(data_wmap[name], m["gT_model"], sims[name])
        # rotated-map errors: scatter of the same estimator on 35 rotations
        n, p = sims[name].shape
        icov = np.linalg.inv(np.cov(sims[name], rowvar=False)) * (n - p - 2) / (n - 1)
        w = icov @ m["gT_model"] / (m["gT_model"] @ icov @ m["gT_model"])
        a_rot = rot[name] @ w
        a_sims[name] = fit["A_sims"]
        res["tracers"][name] = {
            "bias": m["bias_label"], "bias_amp": m["bias_amp"],
            "ell": X.ELL_BIN, "cl_gT_data": data[name], "cl_gT_wmap": data_wmap[name],
            "cl_gT_model": m["gT_model"], "cl_gT_sigma": sims[name].std(axis=0, ddof=1),
            "cl_gT_null_mean": sims[name].mean(axis=0),
            # gg_model was computed with the fitted bias already in the window
            "cl_gg_data": m["gg_data"], "cl_gg_model": m["gg_model"],
            "cl_gg_err": m["gg_err"],
            "A": fit["A"], "sigma": fit["sigma"], "sigma_sims": fit["sigma_sims"], "snr": fit["snr"],
            "chi2": fit["chi2"], "dof": fit["dof"], "A_wmap": fit_w["A"],
            "sigma_rotations": float(np.std(a_rot, ddof=1)), "A_rotations_mean": float(a_rot.mean()),
            "n_rotations": int(len(a_rot)),
            "meta": {k: v for k, v in tracers_meta(tr[name]).items()}}
    for t in TAGS:
        k = f"template_{t}"
        s = sims[k]
        res["template"][t] = {"A": data[k], "sigma": float(np.std(s, ddof=1)),
                              "null_mean": float(np.mean(s)), "snr": data[k] / float(np.std(s, ddof=1)),
                              "A_wmap": data_wmap[k], "sigma_rotations": float(np.std(rot[k], ddof=1))}
    a_sims["template"] = sims["template_fiducial"] - sims["template_fiducial"].mean()
    _log("cross-correlations done", t0)
    return res, a_sims


def tracers_meta(t):
    return {k: t[k] for k in ("nbar_sr", "fsky", "n_sources", "n_gal") if k in t}


def combine(measurements: dict, a_sims: dict, keys) -> dict:
    """Generalised least-squares average of several A_ISW with their joint covariance."""
    keys = [k for k in keys if k in measurements]
    a = np.array([measurements[k] for k in keys])
    m = np.stack([a_sims[k] for k in keys], axis=1)
    cov = np.cov(m, rowvar=False)
    icov = np.linalg.inv(cov)
    one = np.ones(len(keys))
    sig2 = 1.0 / (one @ icov @ one)
    a_comb = float(sig2 * one @ icov @ a)
    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    chi2 = float((a - a_comb) @ icov @ (a - a_comb))
    return {"probes": keys, "A": a_comb, "sigma": float(np.sqrt(sig2)),
            "snr": a_comb / float(np.sqrt(sig2)), "p_zero_one_sided": float(norm.sf(a_comb / np.sqrt(sig2))),
            "chi2_consistency": chi2, "dof": len(keys) - 1, "correlation": corr}


def coldspot_results(t0) -> dict:
    an = CS.Analyzer()
    obs = an.score(an.map128)
    sims = CS.run_sims(N_SIMS)
    res = {"filters": {}}
    for k, v in obs.items():
        s = sims[k]["min_sigma"]
        res["filters"][f"{k[0]}_{k[1]:g}"] = {
            **v, "at_literature_position": CS.coldspot_value(an, an.map128, k),
            "p_min": float((1 + np.sum(s <= v["min_sigma"])) / (len(s) + 1)),
            "sims_min_percentiles": {p: float(np.percentile(s, p)) for p in (1, 5, 50)},
            "p_kurtosis": float((1 + np.sum(sims[k]["kurtosis"] >= v["kurtosis"])) / (len(s) + 1))}
    # the supervoid
    import healpy as hp

    key = ("smhw", SMHW_SCALE_DEG)
    centre = (obs[key]["l"], obs[key]["b"])
    d = cmb.prepare_maps()
    edges = np.arange(0, 26, 1.0)
    prof = CS.radial_profile(d["smica"], d["mask"], 512, centre, edges)
    fdat = an.filtered(an.map128)[key]
    pix = hp.ang2pix(CS.NSIDE_CS, *centre, lonlat=True)
    th = np.linspace(0, 60, 601)
    models = {}
    for name, (d0, r0, z, q) in CS.SUPERVOID_MODELS.items():
        p = CS.supervoid_profile(d0, r0, z, th, q)
        mm = CS.profile_to_map(th, p, centre)
        f = hp.alm2map(hp.almxfl(hp.map2alm(mm, lmax=CS.LMAX_CS), an.filters[key]), CS.NSIDE_CS,
                       lmax=CS.LMAX_CS)
        inner = float(np.mean(p[th < 5]))
        models[name] = {"delta0": d0, "r0": r0, "z": z, "quoted": q, "centre_uK": float(p[0]),
                        "mean_within_5deg_uK": inner,
                        "fraction_of_central_decrement": inner / float(np.nanmean(prof[:5])),
                        "fraction_of_smhw_coefficient": float(f[pix] / fdat[pix]),
                        "theta": th, "profile": p}
    # null scatter of a ring profile around a random point, from simulations
    prof_sims = []
    cl = cmb.theory_cl()
    rng = np.random.default_rng(3)
    for i in range(200):
        m = cmb.alm_to_map(cmb.sim_alm(i, cl), 512)
        lon, lat = rng.uniform(0, 360), np.degrees(np.arcsin(rng.uniform(-1, -0.4)))
        prof_sims.append(CS.radial_profile(m, d["mask"], 512, (lon, lat), edges))
    res["centre"] = centre
    res["profile"] = {"edges": edges, "data": prof, "null_sigma": np.nanstd(prof_sims, axis=0)}
    res["smhw_coefficient_uK"] = float(fdat[pix])
    res["supervoid_models"] = models
    _log("cold spot done", t0)
    return res


def main() -> int:
    t0 = time.monotonic()
    cmb.prepare_maps()
    _log("CMB maps ready", t0)
    SS.build_layout()
    simpass.bank_sims()
    _log("filter bank and simulations ready", t0)
    bias = boss.measure_bias()
    for t in TAGS:
        boss.predicted_isw(t)
    _log("BOSS reconstruction ready", t0)

    st, a_st = stacking_results(t0)
    cc, a_cc = crosscorr_results(t0)
    a_meas = {"granett": st["matched_filter"]["granett"]["fiducial"]["A"],
              "nh14": st["matched_filter"]["nh14"]["fiducial"]["A"],
              "mao": st["matched_filter"]["mao"]["fiducial"]["A"],
              "template": cc["template"]["fiducial"]["A"],
              **{k: v["A"] for k, v in cc["tracers"].items()}}
    a_sims = {**a_st, **a_cc}
    combos = {
        "nvss_template": combine(a_meas, a_sims, ["nvss", "template"]),
        "nvss_template_voids": combine(a_meas, a_sims, ["nvss", "template", "mao", "nh14"]),
        "all_crosscorr": combine(a_meas, a_sims, ["nvss", "boss_0.15_0.43", "boss_0.43_0.75"]),
        "voids_only": combine(a_meas, a_sims, ["mao", "nh14"]),
    }
    _log("amplitudes combined", t0)
    cs = coldspot_results(t0)

    gal = boss.galaxies()
    results = {
        "_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "preregistration": json.loads(PREREG_PATH.read_text()),
        "n_sims": N_SIMS,
        "boss": {"n_galaxies": int(len(gal["z"])), "bias": {k: {kk: v[kk] for kk in
                                                                ("zeff", "b", "b_err", "b_sigma8",
                                                                 "sigma8_0", "chi2", "nbins")}
                                                            for k, v in bias.items()}},
        "stacking": st, "crosscorr": cc, "amplitudes": a_meas, "combined": combos,
        "coldspot": cs, "literature": {k: {kk: (str(vv) if isinstance(vv, dict) else vv)
                                           for kk, vv in v.items()} for k, v in LITERATURE.items()},
    }
    RESULTS_PATH.write_text(json.dumps(_clean(results), indent=1))
    _log(f"wrote {RESULTS_PATH}", t0)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
