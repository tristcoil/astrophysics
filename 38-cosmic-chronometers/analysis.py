"""Every fit and derived number. ``run_all()`` returns (results, artefacts):
``results`` is JSON-serialisable (written to data/results.json); ``artefacts``
holds the arrays the figures need (pickled to data/analysis_cache.pkl).
"""

from __future__ import annotations

import time

import numpy as np
from scipy import linalg

import clock
import config as cfg
import cosmology as cosmo
import fetch_data
import gp
import likelihood as lk

H0_AXIS = np.linspace(30.0, 120.0, 9001)
ZS = np.linspace(0.0, cfg.Z_MAX_PLOT, 251)


def _s(x, w=None):
    return lk.summarize(x, w)


def _flat_lcdm_age(H0, Om):
    """Age in Gyr of a flat matter + Lambda universe (no radiation; analytic)."""
    OL = 1.0 - Om
    return (2.0 / 3.0) * cosmo.T_H / (H0 * np.sqrt(OL)) * np.arcsinh(np.sqrt(OL / Om))


def _flat_summary(fit: lk.GridFit, rng, n=200_000) -> dict:
    s = fit.sample(n, rng)
    zt = cosmo.z_transition_flat(s["Om"])
    return {
        "H0": lk.grid_summary(H0_AXIS, fit.h0_density(H0_AXIS)),
        "Om": lk.grid_summary(*fit.marginal("Om")),
        "z_transition": _s(zt),
        "frac_accelerating": float(np.mean(s["Om"] < 2.0 / 3.0)),
        "q0": _s(1.5 * s["Om"] - 1.0),
        "age_gyr": _s(_flat_lcdm_age(s["H0"], s["Om"])),
        "best": fit.best(),
    }


def _gp_run(z, H, cov, *, kind="se", anchor=None, rng, n_curves=cfg.GP_N_FUNCTION_SAMPLES,
            marginalise=False, seed=cfg.RNG_SEED, mean_sd=(0.0, 0.0)) -> dict:
    """GP reconstruction; optionally with an H0 anchor added at z = 0."""
    if anchor is not None:
        z = np.concatenate([[0.0], z])
        H = np.concatenate([[anchor[0]], H])
        cov = linalg.block_diag([[anchor[1] ** 2]], cov)
    g = gp.GaussianProcess(z, H, cov, kind=kind, mean_sd=mean_sd)
    theta = g.optimize()
    if marginalise:
        hyper = g.sample_hyperparameters(cfg.GP_N_HYPER_SAMPLES, seed)
        curves = np.vstack([g.sample(ZS, th, 1, rng) for th in hyper])
    else:
        hyper = None
        curves = g.sample(ZS, theta, n_curves, rng)
    d = gp.derived_from_curves(ZS, curves)
    out = {
        "theta": {"sigma_f": float(np.exp(theta[0])), "ell": float(np.exp(theta[1])),
                  "ln_evidence": g.log_marginal(theta)},
        "H0": _s(d["H0"]), "q0": _s(d["q0"]),
        "p_accelerating_today": float(np.mean(d["q0"] < 0)),
        "z_transition": _s(d["zt"]),
        "frac_transition_in_data": d["frac_transition"],
        "frac_decelerating_today": d["frac_decelerating_today"],
        "frac_no_turn_in_data": d["frac_no_turn_in_data"],
    }
    if hyper is not None:
        out["hyper"] = {"sigma_f": _s(np.exp(hyper[:, 0])), "ell": _s(np.exp(hyper[:, 1]))}
    bands = {}
    for key, arr in (("H", curves), ("adot", d["adot"]), ("q", d["q"])):
        bands[key] = np.nanpercentile(arr, [2.275, 15.865, 50, 84.135, 97.725], axis=0)
    art = {"bands": bands, "zt": d["zt"], "H0": d["H0"], "q0": d["q0"]}
    if anchor is None:
        mu, cv = g.predict(ZS, theta)
        art["predict"] = (mu, cv)
        art["gp"] = g
        art["theta"] = theta
    return out, art


def _desi_hubble(desi, rd, rd_err):
    """H(z) = c / (D_H/r_d * r_d) from DESI's D_H/r_d, with r_d from the CMB."""
    sel = desi["quantity"] == "DH_over_rs"
    z, dh = desi["z"][sel], desi["value"][sel]
    cov_dh = desi["cov"][np.ix_(sel, sel)]
    H = cfg.C_KM_S / (dh * rd)
    J = np.diag(-H / dh)
    cov = J @ cov_dh @ J.T + np.outer(H, H) * (rd_err / rd) ** 2
    order = np.argsort(z)
    return {"z": z[order], "H": H[order], "cov": cov[np.ix_(order, order)],
            "err": np.sqrt(np.diag(cov))[order]}


def _pull(value, err_lo, err_hi, model, model_err):
    """Signed distance from a prediction, using the error bar on its side."""
    err = err_lo if value > model else err_hi
    return (value - model) / np.hypot(err, model_err)


class _Clock:
    """Prints how long each stage took (the GP stages dominate)."""

    def __init__(self):
        self.t = time.time()

    def __call__(self, label):
        now = time.time()
        print(f"  {label:48s} {now - self.t:6.1f} s", flush=True)
        self.t = now


def run_all(seed: int = cfg.RNG_SEED):
    rng = np.random.default_rng(seed)
    results: dict = {}
    art: dict = {}
    tick = _Clock()

    d1, d3 = fetch_data.load_d1(), fetch_data.load_d3()
    results["inputs"] = {"D1": d1, "D3": d3}
    cc = lk.compile_chronometers()
    syst = fetch_data.load_systematics()
    z, H, err = cc["z"], cc["H"], cc["err"]
    # Baseline: Moresco et al.'s covariance with the shared systematic scaled by
    # a smooth reconstruction of H(z); "published recipe" scales it by the
    # measured values, which biases fits low (see smooth_scaled_covariance).
    lib_terms = cfg.SYSTEMATIC_TERMS + ("stlib",)
    cov_smooth, h_smooth = lk.smooth_scaled_covariance(z, H, err, syst=syst)
    cov_smooth_lib, _ = lk.smooth_scaled_covariance(z, H, err, syst=syst, terms=lib_terms)
    covs = {
        "full": cov_smooth,
        "full, published recipe": lk.chronometer_covariance(z, H, err, syst=syst),
        "full + stellar library": cov_smooth_lib,
        "full + stellar library, published recipe": lk.chronometer_covariance(
            z, H, err, syst=syst, terms=lib_terms),
        "systematics on Moresco's 15 only": lk.chronometer_covariance(
            z, h_smooth, err, syst=syst, apply_to=cc["is_moresco"]),
        "diagonal only": np.diag(err**2),
    }
    cov = covs["full"]
    art.update(cc=cc, syst=syst, covs=covs, h_smooth=h_smooth)
    tick("covariances")

    # -- the data set ------------------------------------------------------
    diag_sd = np.sqrt(np.diag(cov))
    results["data"] = {
        "n": int(z.size), "z_range": [float(z.min()), float(z.max())],
        "by_method": {m: int(np.sum(cc["method"] == m)) for m in cfg.METHODS},
        "n_moresco": int(cc["is_moresco"].sum()),
        "table_vs_gitlab": [
            {"z_table": float(a), "z": float(b), "H_table": float(c), "H": float(d),
             "err_table": float(e), "err": float(f)}
            for a, b, c, d, e, f in zip(cc["z_table"], z, cc["H_table"], H, cc["err_table"], err)
            if (abs(a - b) > 1e-9 or abs(c - d) > 1e-9 or abs(e - f) > 1e-9)],
        "median_frac_error_diag": float(np.median(err / H)),
        "median_frac_error_total": float(np.median(diag_sd / H)),
        "median_syst_frac": float(np.median(np.sqrt(np.diag(cov) - err**2) / H)),
        "max_smooth_vs_measured": float(np.max(np.abs(h_smooth / H - 1))),
    }
    corr = cov / np.outer(diag_sd, diag_sd)
    results["data"]["mean_offdiag_correlation"] = float(
        corr[~np.eye(z.size, dtype=bool)].mean())

    # -- flat LCDM: the covariance treatments ----------------------------------
    flat = {name: lk.fit_flat_lcdm(z, H, C) for name, C in covs.items()}
    results["flat_lcdm"] = {name: _flat_summary(f, rng) for name, f in flat.items()}
    art["flat_fits"] = flat
    om_axis = np.linspace(0.02, 0.98, 481)
    h_axis = np.linspace(45, 95, 501)
    art["planes"] = {name: (h_axis, om_axis, lk.flat_lcdm_plane(z, H, covs[name], h_axis, om_axis))
                     for name in ("full", "full, published recipe", "diagonal only")}
    art["h0_densities"] = {name: f.h0_density(H0_AXIS) for name, f in flat.items()}
    art["H0_AXIS"] = H0_AXIS

    tick("flat LCDM, six covariance treatments")

    # -- validation against published fits ------------------------------------
    mm = fetch_data.load_moresco_table("BC03")
    cov15 = lk.chronometer_covariance(mm["z"], mm["H"], mm["err"], syst=syst)
    val = {"primer_cov_determinant": float(np.linalg.det(cov15))}
    for label, C in (("full", cov15), ("diag", np.diag(mm["err"] ** 2))):
        f = lk.fit_flat_lcdm(mm["z"], mm["H"], C)
        val[f"primer_{label}"] = {"H0": lk.grid_summary(H0_AXIS, f.h0_density(H0_AXIS)),
                                  "Om": lk.grid_summary(*f.marginal("Om"))}
    rec = results["flat_lcdm"]["full, published recipe"]
    val["lrr_32"] = {"H0": rec["H0"], "Om": rec["Om"],
                     "chi2_red": rec["best"]["chi2"] / rec["best"]["dof"]}
    results["validation"] = val
    tick("validation fits")

    # -- how much does scaling by the measured values bias a fit? (mock data) --
    th = gp.GaussianProcess(z, H, cov).optimize()
    results["scaling_bias"], art["scaling_bias"] = scaling_bias_simulation(
        z, err, syst, (d3["H0"][0], d3["Om"][0]), th, rng)
    tick("scaling-bias simulation")

    # -- open LCDM (curvature free) -----------------------------------------
    olcdm = lk.fit_lcdm(z, H, cov)
    s = olcdm.sample(200_000, rng)
    results["open_lcdm"] = {
        "Om": lk.grid_summary(*olcdm.marginal("Om")),
        "OL": lk.grid_summary(*olcdm.marginal("OL")),
        "Ok": _s(1 - s["Om"] - s["OL"]),
        "H0": lk.grid_summary(H0_AXIS, olcdm.h0_density(H0_AXIS)),
        "q0": _s(0.5 * s["Om"] - s["OL"]),
        "p_accelerating_today": float(np.mean(0.5 * s["Om"] - s["OL"] < 0)),
        "z_transition": _s(cosmo.z_transition_lcdm(s["Om"], s["OL"])),
        "best": olcdm.best(),
    }
    art["olcdm"] = olcdm

    tick("open LCDM")

    # -- Moresco et al. 2016, replicated ----------------------------------------
    m16 = fetch_data.load_m16_table()
    cov16 = np.diag(m16["err"] ** 2)
    prior16 = cfg.PUBLISHED["m16_H0_prior"]
    f16 = lk.fit_lcdm_zt(m16["z"], m16["H"], cov16, h0_prior=prior16)
    zt_ax = f16.axes[1]
    prof = np.nanmin(f16.chi2_profile, axis=0)
    y16, cy16 = lk.adot_data(m16["z"], m16["H"], cov16)
    pw16 = lk.piecewise_transition(m16["z"], y16, cy16, anchor=prior16)

    def _profile_intervals(axis, chi2):
        d = chi2 - np.nanmin(chi2)
        out = {"best": float(axis[np.nanargmin(chi2)])}
        for lev, name in ((1.0, "68"), (4.0, "95"), (9.0, "99.7")):
            inside = axis[d < lev]
            out[name] = [float(inside.min()), float(inside.max())]
        return out

    results["moresco2016_replication"] = {
        "olcdm_profile": _profile_intervals(zt_ax, prof),
        "olcdm_marginal": lk.grid_summary(*f16.marginal("zt")),
        "piecewise_profile": _profile_intervals(pw16["zt_axis"], pw16["chi2"]),
        "piecewise_vs_line": {"piecewise": pw16["piecewise_ic"], "line": pw16["line_ic"]},
    }
    art["m16"] = {"olcdm_axis": zt_ax, "olcdm_chi2": prof, "piecewise": pw16}

    tick("Moresco 2016 replication")

    # -- the transition, model-independently: two straight lines ----------------
    ydot, cydot = lk.adot_data(z, H, cov)
    anchors = {
        "none": None,
        "D1 H0 (Cepheids + supernovae)": d1["H0"],
        "D3 H0 (Planck CMB in LCDM)": d3["H0"],
    }
    pw = {}
    for name, anc in anchors.items():
        r = lk.piecewise_transition(z, ydot, cydot, anchor=anc)
        pw[name] = r
        results.setdefault("piecewise", {})[name] = {
            "profile": _profile_intervals(r["zt_axis"], r["chi2"]),
            "best": r["best"], "line": r["line"],
            "delta_chi2_line_minus_piecewise": r["delta_chi2_line"],
            "AICc": {"piecewise": r["piecewise_ic"]["AICc"], "line": r["line_ic"]["AICc"]},
            "BIC": {"piecewise": r["piecewise_ic"]["BIC"], "line": r["line_ic"]["BIC"]},
        }
    art["piecewise"] = pw

    tick("two-line transition fits")

    # -- Gaussian process ----------------------------------------------------
    gp_res, gp_art = {}, {}
    runs = {
        "SE kernel, full covariance": dict(cov=cov, kind="se"),
        "SE kernel, published-recipe covariance": dict(cov=covs["full, published recipe"],
                                                       kind="se"),
        "SE kernel, hyperparameters marginalised": dict(cov=cov, kind="se", marginalise=True),
        "Matern-9/2 kernel, full covariance": dict(cov=cov, kind="m92"),
        "SE kernel, constant mean marginalised": dict(cov=cov, kind="se", mean_sd=(300.0, 0.0)),
        "SE kernel, straight-line mean marginalised": dict(cov=cov, kind="se",
                                                           mean_sd=(300.0, 300.0)),
        "SE kernel, diagonal errors only": dict(cov=covs["diagonal only"], kind="se"),
        "SE kernel, + D1 H0 anchor": dict(cov=cov, kind="se", anchor=d1["H0"]),
        "SE kernel, + D3 H0 anchor": dict(cov=cov, kind="se", anchor=d3["H0"]),
    }
    for name, kw in runs.items():
        C = kw.pop("cov")
        gp_res[name], gp_art[name] = _gp_run(z, H, C, rng=rng, seed=seed, **kw)
        tick(f"GP: {name}")
    results["gp"] = gp_res
    art["gp"] = gp_art
    art["ZS"] = ZS

    # -- method checks on an exact LCDM curve (no noise) --------------------
    checks = {}
    for label, (h0, om) in {"D3 LCDM": (d3["H0"][0], d3["Om"][0]),
                            "Planck 2018 LCDM": (cfg.PLANCK18["H0"][0], cfg.PLANCK18["Om"][0])}.items():
        h_true = cosmo.hubble(z, h0, om)
        entry = {"true_zt": float(cosmo.z_transition_flat(om)), "true_H0": float(h0),
                 "true_q0": float(1.5 * om - 1)}
        for mname, msd in (("zero mean", (0.0, 0.0)), ("constant mean marginalised", (300.0, 0.0)),
                           ("straight-line mean marginalised", (300.0, 300.0))):
            g = gp.GaussianProcess(z, h_true, cov, mean_sd=msd)
            th = g.optimize()
            mu, _ = g.predict(ZS, th)
            dd = gp.derived_from_curves(ZS, mu[None, :])
            entry[f"gp ({mname})"] = {"zt": float(dd["zt"][0]), "H0": float(mu[0]),
                                      "q0": float(dd["q0"][0])}
        yt, cyt = lk.adot_data(z, h_true, cov)
        entry["piecewise_zt"] = lk.piecewise_transition(z, yt, cyt)["best"]["zt"]
        entry["piecewise_zt_with_exact_H0"] = lk.piecewise_transition(
            z, yt, cyt, anchor=(h0, d1["H0"][1]))["best"]["zt"]
        flat_exact = lk.fit_flat_lcdm(z, h_true, cov)
        entry["flat_lcdm_fit"] = {"H0": float(flat_exact.best()["H0"]),
                                  "Om": float(flat_exact.best()["Om"])}
        checks[label] = entry
    results["method_checks"] = checks

    tick("method checks on exact LCDM")

    # -- DESI DR2: H(z) from the BAO scale, calibrated by the CMB ---------------
    desi = fetch_data.load_desi_dr2()
    desi_h = _desi_hubble(desi, *d3["rdrag"])
    desi_h_planck = _desi_hubble(desi, *cfg.PLANCK18["rdrag"])
    in_range = desi_h["z"] <= cfg.Z_MAX_DATA
    g_main = gp_art["SE kernel, full covariance"]["gp"]
    th_main = gp_art["SE kernel, full covariance"]["theta"]
    mu_d, cov_gp_d = g_main.predict(desi_h["z"][in_range], th_main)
    resid = desi_h["H"][in_range] - mu_d
    best_flat = results["flat_lcdm"]["full"]["best"]
    def _chi2_model(dz, dH, dcov, h0, om):
        return lk.chi2(dH - cosmo.hubble(dz, h0, om), dcov)
    results["desi"] = {
        "rd_used": d3["rdrag"],
        "points": [{"z": float(a), "H": float(b), "err": float(c), "adot": float(b / (1 + a))}
                   for a, b, c in zip(desi_h["z"], desi_h["H"], desi_h["err"])],
        "chi2_vs_cc_gp": lk.chi2(resid, desi_h["cov"][np.ix_(in_range, in_range)] + cov_gp_d),
        "n_vs_cc_gp": int(in_range.sum()),
        "chi2_vs_cc_flat_best": _chi2_model(desi_h["z"], desi_h["H"], desi_h["cov"],
                                            best_flat["H0"], best_flat["Om"]),
        "chi2_vs_D1": _chi2_model(desi_h["z"], desi_h["H"], desi_h["cov"],
                                  d1["H0"][0], d1["Om_calibrated"][0]),
        "chi2_vs_D3": _chi2_model(desi_h["z"], desi_h["H"], desi_h["cov"],
                                  d3["H0"][0], d3["Om"][0]),
        "n": int(desi_h["z"].size),
        "ratio_D1_over_DESI_at_z051": float(cosmo.hubble(0.51, d1["H0"][0], d1["Om_calibrated"][0])
                                            / desi_h["H"][0]),
        "rd_planck_shift_percent": float(100 * (desi_h_planck["H"][0] / desi_h["H"][0] - 1)),
    }
    art["desi"] = desi_h

    # -- measurements published since the compilation -----------------------
    mu_all, cov_all = g_main.predict(ZS, th_main)
    sd_all = np.sqrt(np.diag(cov_all))
    since = []
    for zz, hh, lo, hi, label, note in cfg.SINCE_2022:
        m = float(np.interp(zz, ZS, mu_all))
        sm = float(np.interp(zz, ZS, sd_all))
        m_flat = float(cosmo.hubble(zz, best_flat["H0"], best_flat["Om"]))
        since.append({"z": zz, "H": hh, "err_minus": lo, "err_plus": hi, "label": label,
                      "note": note, "gp_mean": m, "gp_sd": sm,
                      "pull_vs_gp": float(_pull(hh, lo, hi, m, sm)),
                      "pull_vs_flat_best": float(_pull(hh, lo, hi, m_flat, 0.0))})
    results["since_2022"] = since

    # -- inside a chronometer: the 4000 Å break ------------------------------
    m12 = clock.reread_moresco2012()
    simple = clock.sdss_simple_cut_medians()
    simple_h = {model: clock.h_from_rows(simple, rows) for model, rows in m12["A_rows"].items()} \
        if simple else {}
    results["clock"] = {
        "moresco2012_reread": m12["combined"],
        "moresco2012_published": {str(k): v for k, v in cfg.PUBLISHED["m12_sdss"].items()},
        "moresco2012_published_M11": {"0.1791": 80.52, "0.1993": 80.61},
        "simple_cuts_medians": simple,
        "simple_cuts_H": {model: {m: [{"z": p["z"], "H": p["H"], "stat": p["stat"]}
                                      for p in prs] for m, prs in per.items()}
                          for model, per in simple_h.items()},
    }
    art["clock"] = {"m12": m12, "simple": simple, "simple_h": simple_h}

    # -- the speedometer reading -------------------------------------------
    h0_3, om_3 = d3["H0"][0], d3["Om"][0]
    zt_3 = float(cosmo.z_transition_flat(om_3))
    adot_min = float(h0_3 * np.sqrt(cosmo.E2(zt_3, om_3)) / (1 + zt_3))
    results["speedometer"] = {
        "model": "D3 LCDM", "adot_today": h0_3, "z_transition": zt_3,
        "adot_at_transition": adot_min,
        "lookback_to_transition_gyr": float(cosmo.lookback_time(zt_3, h0_3, om_3)),
        "adot_at_z2": float(h0_3 * np.sqrt(cosmo.E2(2.0, om_3)) / 3.0),
        "d1_lookback_to_transition_gyr": float(cosmo.lookback_time(
            d1["z_transition"][0], d1["H0"][0], d1["Om_calibrated"][0])),
    }

    tick("DESI, new measurements, clock, speedometer")

    # -- the future -----------------------------------------------------------
    results["future"], art["future"] = future_scenarios(d1, d3)
    tick("future scenarios")
    return results, art


def scaling_bias_simulation(z, err, syst, truth, theta_gp, rng, n_mock=300):
    """Mock chronometer data drawn from a known flat LCDM with the full
    covariance, fitted three ways: the shared systematic scaled by each mock's
    own measured values (the published recipe), by the smooth GP mean (our
    baseline; hyperparameters fixed at the real data's optimum), and by the
    true curve (the unattainable ideal). Compares the H0 they return."""
    hm = cosmo.hubble(z, *truth)
    out, art = {}, {}
    om_ax = np.linspace(0.01, 0.99, 197)
    h_ax = np.linspace(30.0, 120.0, 1801)
    for label, terms in (("standard terms", cfg.SYSTEMATIC_TERMS),
                         ("with stellar library", cfg.SYSTEMATIC_TERMS + ("stlib",))):
        cov_true = lk.chronometer_covariance(z, hm, err, syst=syst, terms=terms)
        L = np.linalg.cholesky(cov_true)
        est = {"measured values (published recipe)": [], "smooth reconstruction (ours)": [],
               "true curve (ideal)": []}
        for _ in range(n_mock):
            y = hm + L @ rng.standard_normal(z.size)
            scale = y
            for _ in range(4):
                c = lk.chronometer_covariance(z, scale, err, syst=syst, terms=terms)
                scale, _ = gp.GaussianProcess(z, y, c).predict(z, theta_gp)
            for key, C in (("measured values (published recipe)",
                            lk.chronometer_covariance(z, y, err, syst=syst, terms=terms)),
                           ("smooth reconstruction (ours)",
                            lk.chronometer_covariance(z, scale, err, syst=syst, terms=terms)),
                           ("true curve (ideal)", cov_true)):
                f = lk.fit_flat_lcdm(z, y, C, om_axis=om_ax)
                est[key].append(lk.grid_summary(h_ax, f.h0_density(h_ax))["median"])
        est = {k: np.array(v) for k, v in est.items()}
        out[label] = {"truth_H0": truth[0], "n_mock": n_mock,
                      **{k: {"mean": float(v.mean()), "std": float(v.std()),
                             "bias": float(v.mean() - truth[0])} for k, v in est.items()},
                      "expected_factor_1_over_1_plus_nf2": float(1.0 / (1.0 + z.size * np.mean(
                          sum((np.interp(z, syst["z"], syst[t]) / 100) ** 2 for t in terms))))}
        art[label] = est
    return out, art


def future_scenarios(d1, d3):
    """a(t) from the Big Bang to +100 Gyr for several dark energies (projections)."""
    h0, om = d3["H0"][0], d3["Om"][0]
    desi = cfg.DESI_W0WA["DESI+CMB+Pantheon+"]
    scen = {
        "LCDM (D3's CMB fit)": dict(H0=h0, Om=om),
        "Evolving dark energy (DESI DR2+CMB+Pantheon+, extrapolated)":
            dict(H0=desi["H0"], Om=desi["Om"], w0=desi["w0"], wa=desi["wa"]),
        f"Phantom, w = {cfg.D1_W_95_LOW} (D1's 95% limit)": dict(H0=h0, Om=om, w0=cfg.D1_W_95_LOW),
        "Phantom, w = -1.5 (Caldwell et al. 2003)": dict(H0=h0, Om=om, w0=-1.5),
        "Matter only (no dark energy)": dict(H0=h0, Om=1.0),
    }
    res, art = {}, {}
    for name, p in scen.items():
        hist = cosmo.expansion_history(p["H0"], p["Om"], None, p.get("w0", -1.0),
                                       p.get("wa", 0.0), ln_a_max=14.0, n=120001)
        t_now = float(np.interp(0.0, hist["ln_a"], hist["t"]))
        entry = {"params": p, "age_gyr": t_now}
        if p.get("w0", -1.0) < -1 and p.get("wa", 0.0) == 0.0:
            entry["big_rip_from_now_gyr"] = cosmo.big_rip_time(p["H0"], p["Om"], p["w0"])
            entry["big_rip_approx_gyr"] = cosmo.big_rip_time_approx(p["H0"], p["Om"], p["w0"])
        if "wa" in p:
            entry["acceleration_end"] = cosmo.acceleration_end(p["H0"], p["Om"], p["w0"], p["wa"])
            entry["w_today"] = p["w0"]
            entry["w_at_a2"] = float(cosmo.w_of_a(2.0, p["w0"], p["wa"]))
        sel = hist["t"] - t_now <= 100.0
        entry["a_at_plus100gyr"] = float(np.interp(100.0, hist["t"] - t_now, hist["a"])) \
            if (hist["t"] - t_now).max() > 100 else None
        zt = cosmo.z_transition(p["Om"], None, p.get("w0", -1.0), p.get("wa", 0.0))
        entry["z_transition"] = None if not np.isfinite(zt) else float(zt)
        res[name] = entry
        art[name] = {"t": hist["t"][sel] - t_now, "a": hist["a"][sel]}
    # the other DESI supernova combinations: does the acceleration end, and when?
    res["DESI w0wa variants"] = {
        k: {"params": v, "acceleration_end": cosmo.acceleration_end(v["H0"], v["Om"], v["w0"], v["wa"])}
        for k, v in cfg.DESI_W0WA.items()}
    # the most phantom constant w allowed by DESI+CMB+Pantheon+ at 95%
    w_lo = cfg.DESI_WCDM["w"][0] - 2 * cfg.DESI_WCDM["w"][1]
    res["DESI wCDM 95% phantom limit"] = {
        "w": w_lo, "big_rip_from_now_gyr": cosmo.big_rip_time(cfg.DESI_WCDM["H0"],
                                                              cfg.DESI_WCDM["Om"], w_lo)}
    c03 = cfg.CALDWELL03
    res["Caldwell 2003 check"] = {
        "published_gyr": c03["t_rip_from_now_gyr"],
        "approx_formula_gyr": cosmo.big_rip_time_approx(c03["H0"], c03["Om"], c03["w"]),
        "exact_gyr": cosmo.big_rip_time(c03["H0"], c03["Om"], c03["w"]),
    }
    res["timeline"] = rip_timeline()
    return res, art


def rip_timeline() -> dict:
    """How long before the Big Rip each bound system comes apart."""
    G = 6.67430e-11
    msun, pc = 1.98847e30, 3.0856775814913673e16
    periods = {}
    for label, P in cfg.BOUND_SYSTEMS:
        if P is None and "cluster" in label:
            M, R = 1e15 * msun, 2e6 * pc
            P = 2 * np.pi * np.sqrt(R**3 / (G * M))
        elif P is None and "Milky Way" in label:
            P = 2 * np.pi * 8.2e3 * pc / 233e3
        elif P is None and "Earth itself" in label:
            P = 2 * np.pi * np.sqrt(6.371e6**3 / 3.986004e14)
        periods[label] = P
    bohr = 2 * np.pi * 5.29177e-11 / 2.18769e6  # the hydrogen electron's orbital period
    periods["A hydrogen atom (electron's orbital period)"] = bohr
    out = {}
    for w in (-1.5, cfg.D1_W_95_LOW):
        out[str(w)] = {label: float(cosmo.rip_lead_time(P, w)) for label, P in periods.items()}
    out["periods_s"] = {k: float(v) for k, v in periods.items()}
    return out
