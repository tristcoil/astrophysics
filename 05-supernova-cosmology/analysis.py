"""The analysis: every fit, every derived number, and the arrays the figures need.

Run via ``run.py``. Nothing here downloads anything; it reads the cache that
``fetch_data.py`` fills.

The fits, in the order a cosmologist would do them:

1. **Flat LambdaCDM from supernovae alone.** One shape parameter, Omega_m.
2. **The (Omega_m, Omega_Lambda) plane, curvature left free.** The 1998 figure.
   How confidently is Omega_Lambda > 0, and is the universe accelerating?
3. **Textbook universes compared.** Empty, Einstein-de Sitter, open matter-only,
   de Sitter, Planck's CMB fit, and our best fit, each allowed its own best zero
   point, scored by chi^2.
4. **Flat wCDM.** Is dark energy a cosmological constant (w = -1)?
5. **H0 from the distance ladder.** Adding the Cepheid-calibrated supernovae
   fixes the absolute magnitude M, and with it H0. Then the Hubble tension.
6. **Cross-checks.** The same H0 from only the 277 supernovae SH0ES used with
   a model-free expansion; the effect of statistical vs systematic covariance;
   and our posteriors against the Pantheon+ team's own chains.
"""

from __future__ import annotations

import json
import time

import numpy as np

import cosmology as cosmo
import likelihood as lk
from config import (
    CHAIN_FILES,
    CHAIN_PRIOR_OM,
    CHAIN_PRIOR_W,
    DISCOVERY_1998,
    GRID_FLAT_OM,
    GRID_OL,
    GRID_OL_COARSE,
    GRID_OM,
    GRID_OM_COARSE,
    GRID_W_OM,
    GRID_W_W,
    MANIFEST_PATH,
    PLANCK18_H0,
    PLANCK18_OM,
    PUBLISHED_PANTHEON,
    RESULTS_PATH,
    SH0ES_R22_H0,
    SURVEY_GROUPS,
    SURVEYS,
    Z_MIN_COSMOLOGY,
    Z_SH0ES_HUBBLE_FLOW,
)
from fetch_data import load_chain, load_covariance, load_table

#: Seed for the Monte Carlo propagation of uncertainties into derived numbers.
RNG_SEED = 20260929


def _grid(spec) -> np.ndarray:
    lo, hi, n = spec
    return np.linspace(lo, hi, int(n))


def _log(msg: str, t0: float) -> None:
    print(f"  [{time.time() - t0:5.1f}s] {msg}")


def survey_group(ids: np.ndarray) -> np.ndarray:
    """Index into ``SURVEY_GROUPS`` for each IDSURVEY code."""
    out = np.full(ids.shape, -1)
    for k, members in enumerate(SURVEY_GROUPS.values()):
        out[np.isin(ids, members)] = k
    if (out < 0).any():
        raise ValueError(f"unknown survey IDs: {sorted(set(ids[out < 0]))}")
    return out


def _bin_edges(z: np.ndarray, lo: float, hi: float, n: int, min_count: int) -> np.ndarray:
    """Log-spaced redshift bins, merged from the top down until each is populated."""
    edges = list(np.geomspace(lo, hi, n + 1))
    k = len(edges) - 2
    while k > 0:
        if np.sum((z >= edges[k]) & (z < edges[k + 1])) < min_count:
            del edges[k]  # merge bin k into the one below it
        k -= 1
    return np.array(edges)


def dataset_summary(table) -> dict:
    ids = table["IDSURVEY"].to_numpy()
    z = table["zHD"].to_numpy()
    cal = table["IS_CALIBRATOR"].to_numpy().astype(bool)
    groups = survey_group(ids)
    return {
        "n_light_curves": int(len(table)),
        "n_unique_names": int(table["CID"].nunique()),
        "n_cosmology": int(np.sum(z > Z_MIN_COSMOLOGY)),
        "n_calibrator_light_curves": int(cal.sum()),
        "n_calibrator_supernovae": int(table.loc[cal, "CID"].nunique()),
        "n_excluded_low_z": int(np.sum((z <= Z_MIN_COSMOLOGY) & ~cal)),
        "n_sh0es_hubble_flow": int(table["USED_IN_SH0ES_HF"].sum()),
        "z_min": float(z.min()),
        "z_max": float(z.max()),
        "n_above_z1": int(np.sum(z > 1)),
        "surveys": {SURVEYS[i]: int(np.sum(ids == i)) for i in sorted(SURVEYS)},
        "survey_groups": {label: int(np.sum(groups == k))
                          for k, label in enumerate(SURVEY_GROUPS)},
    }


def run(*, validate: bool = True) -> tuple[dict, dict]:
    """Every fit. Returns (results: JSON-safe numbers, art: arrays for figures)."""
    t0 = time.time()
    print("Loading data")
    table = load_table()
    cov = {"statsys": load_covariance("statsys"), "stat": load_covariance("stat")}
    results: dict = {"dataset": dataset_summary(table)}
    art: dict = {"table": table}
    rng = np.random.default_rng(RNG_SEED)

    sn = lk.build_sample(table, cov["statsys"], name="SN, stat+sys")
    sn_stat = lk.build_sample(table, cov["stat"], name="SN, stat only")
    sn_diag = lk.build_sample(table, cov["statsys"], name="SN, diagonal only",
                              diagonal_only=True)
    cal = lk.build_sample(table, cov["statsys"], name="SN + Cepheid calibrators",
                          calibrators=True)
    _log(f"samples built: {sn.n} Hubble-flow light curves, "
         f"{int(cal.is_calibrator.sum())} calibrator light curves", t0)

    # --- 1. Flat LambdaCDM, supernovae alone -----------------------------------
    om1 = _grid(GRID_FLAT_OM)
    flat = {}
    for label, sample in (("statsys", sn), ("stat", sn_stat), ("diag", sn_diag)):
        chi2, offset = lk.scan(sample, om1, 1.0 - om1)
        post = lk.grid_posterior(chi2)
        i = int(np.argmin(chi2))
        flat[label] = {"Om": lk.summarize_1d(om1, post), "best_Om": float(om1[i]),
                       "chi2_min": float(chi2[i]), "dof": int(sample.n - 2)}
        if label == "statsys":
            art["flat"] = {"Om": om1, "chi2": chi2, "post": post, "offset": offset}
    best_om = flat["statsys"]["best_Om"]
    om_draws = rng.choice(om1, size=200_000, p=art["flat"]["post"])
    flat["q0"] = lk.summarize_samples(cosmo.q0(om_draws, 1.0 - om_draws))
    flat["z_transition"] = lk.summarize_samples(cosmo.transition_redshift(om_draws, 1.0 - om_draws))
    # Goodness of fit with the table's own per-supernova error column instead of
    # the covariance -- to see which one describes the real scatter.
    col_err = table["m_b_corr_err_DIAG"].to_numpy()[sn.rows]
    shape_best = lk.shape_magnitudes(sn, best_om, 1.0 - best_om)[0]
    r = sn.m - shape_best - art["flat"]["offset"][np.argmin(art["flat"]["chi2"])]
    flat["chi2_with_table_error_column"] = float(np.sum((r / col_err) ** 2))
    results["flat_lcdm"] = flat
    _log(f"flat LCDM: Om = {flat['statsys']['Om']['mean']:.4f} "
         f"+/- {flat['statsys']['Om']['std']:.4f}", t0)

    # --- 2. The (Omega_m, Omega_Lambda) plane -----------------------------------
    om2, ol2 = _grid(GRID_OM), _grid(GRID_OL)
    OM, OL = np.meshgrid(om2, ol2, indexing="ij")
    chi2_plane, offset_plane = lk.scan(sn, OM, OL)
    # Same flat priors as the published chain: Om in [0.1, 0.9], Ok in [-1, 1].
    chain_box = ((OM >= CHAIN_PRIOR_OM[0]) & (OM <= CHAIN_PRIOR_OM[1])
                 & (OM + OL >= 0.0) & (OM + OL <= 2.0))
    post_plane = lk.grid_posterior(chi2_plane, chain_box)
    post_plane_physical = lk.grid_posterior(chi2_plane)  # only Om >= 0 (the grid edge)
    ib = np.unravel_index(np.argmin(chi2_plane), chi2_plane.shape)
    plane = {
        "best": {"Om": float(om2[ib[0]]), "OL": float(ol2[ib[1]]),
                 "chi2": float(chi2_plane[ib])},
        "Om": lk.summarize_1d(om2, post_plane.sum(1)),
        "OL": lk.summarize_1d(ol2, post_plane.sum(0)),
        "corr_Om_OL": float(np.sum(post_plane * (OM - np.sum(post_plane * OM))
                                   * (OL - np.sum(post_plane * OL)))
                            / np.sqrt(np.sum(post_plane * (OM - np.sum(post_plane * OM)) ** 2)
                                      * np.sum(post_plane * (OL - np.sum(post_plane * OL)) ** 2))),
    }
    # q0 = Om/2 - OL and Ok = 1 - Om - OL, marginalised over the plane.
    flat_idx = rng.choice(post_plane.size, size=400_000, p=post_plane.ravel())
    om_s, ol_s = OM.ravel()[flat_idx], OL.ravel()[flat_idx]
    om_s = om_s + rng.uniform(-0.5, 0.5, om_s.size) * (om2[1] - om2[0])
    ol_s = ol_s + rng.uniform(-0.5, 0.5, ol_s.size) * (ol2[1] - ol2[0])
    plane["q0"] = lk.summarize_samples(cosmo.q0(om_s, ol_s))
    plane["Ok"] = lk.summarize_samples(1.0 - om_s - ol_s)

    # How strongly is a universe without acceleration / without Lambda excluded?
    # Profile along the boundary lines; the only prior is Om >= 0.
    line = np.linspace(0.0, 1.2, 1201)
    chi2_no_L, _ = lk.scan(sn, line, np.zeros_like(line))
    chi2_q0, _ = lk.scan(sn, line, 0.5 * line)
    in_region_L = np.min(chi2_plane[OL <= 0.0])
    in_region_q = np.min(chi2_plane[OL <= 0.5 * OM])
    best_no_L = min(chi2_no_L.min(), in_region_L)
    best_q0 = min(chi2_q0.min(), in_region_q)
    d_L = best_no_L - plane["best"]["chi2"]
    d_q = best_q0 - plane["best"]["chi2"]
    plane["exclusion"] = {
        "no_lambda": {"delta_chi2": float(d_L), "sigma": lk.sigma_from_delta_chi2(d_L),
                      "best_Om_on_line": float(line[np.argmin(chi2_no_L)])},
        "no_acceleration": {"delta_chi2": float(d_q), "sigma": lk.sigma_from_delta_chi2(d_q),
                            "best_Om_on_line": float(line[np.argmin(chi2_q0)])},
        "posterior_mass_OL_le_0": float(post_plane_physical[OL <= 0].sum()),
        "posterior_mass_q0_ge_0": float(post_plane_physical[OL <= 0.5 * OM].sum()),
        "riess98_sigma": DISCOVERY_1998["riess98_sigma_OL_positive"],
    }
    results["lcdm_plane"] = plane
    art["plane"] = {"Om": om2, "OL": ol2, "chi2": chi2_plane, "post": post_plane,
                    "post_physical": post_plane_physical, "offset": offset_plane}
    _log(f"(Om, OL) plane: {OM.size:,} universes; Lambda > 0 at "
         f"{plane['exclusion']['no_lambda']['sigma']:.1f} sigma", t0)

    # Same plane with statistical-only and diagonal-only errors (coarser grid).
    omc, olc = _grid(GRID_OM_COARSE), _grid(GRID_OL_COARSE)
    OMc, OLc = np.meshgrid(omc, olc, indexing="ij")
    box_c = ((OMc >= CHAIN_PRIOR_OM[0]) & (OMc <= CHAIN_PRIOR_OM[1])
             & (OMc + OLc >= 0.0) & (OMc + OLc <= 2.0))
    variants = {}
    art["plane_variants"] = {"Om": omc, "OL": olc}
    for label, sample in (("statsys", sn), ("stat", sn_stat), ("diag", sn_diag)):
        c2, _ = lk.scan(sample, OMc, OLc)
        p = lk.grid_posterior(c2, box_c)
        art["plane_variants"][label] = p
        variants[label] = {"Om": lk.summarize_1d(omc, p.sum(1)),
                           "OL": lk.summarize_1d(olc, p.sum(0))}
    results["covariance_variants"] = {"flat_Om": {k: flat[k]["Om"] for k in ("statsys", "stat", "diag")},
                                      "plane": variants}
    _log("covariance variants done", t0)

    # --- 3. Textbook universes, each at its own best zero point ---------------
    models = {
        "empty": ("Empty, coasting (Omega_m = 0, Omega_L = 0)", 0.0, 0.0),
        "eds": ("Einstein-de Sitter (Omega_m = 1)", 1.0, 0.0),
        "open": ("Open, matter only (Omega_m = 0.3)", 0.3, 0.0),
        "de_sitter": ("Lambda only (Omega_L = 1)", 0.0, 1.0),
        "planck": ("Planck 2018 flat LambdaCDM", PLANCK18_OM[0], 1.0 - PLANCK18_OM[0]),
        "best_flat": ("Our best-fit flat LambdaCDM", best_om, 1.0 - best_om),
    }
    comparison = {}
    for key, (label, om, ol) in models.items():
        c2, off = lk.scan(sn, np.array([om]), np.array([ol]))
        comparison[key] = {"label": label, "Om": om, "OL": ol, "chi2": float(c2[0]),
                           "offset": float(off[0])}
    ref = comparison["best_flat"]["chi2"]
    for key in comparison:
        comparison[key]["delta_chi2"] = comparison[key]["chi2"] - ref
    results["model_comparison"] = comparison

    # --- 4. Flat wCDM ------------------------------------------------------------
    omw, ww = _grid(GRID_W_OM), _grid(GRID_W_W)
    WOM, WW = np.meshgrid(omw, ww, indexing="ij")
    chi2_w, _ = lk.scan(sn, WOM, 1.0 - WOM, WW)
    box_w = ((WOM >= CHAIN_PRIOR_OM[0]) & (WOM <= CHAIN_PRIOR_OM[1])
             & (WW >= CHAIN_PRIOR_W[0]) & (WW <= CHAIN_PRIOR_W[1]))
    post_w = lk.grid_posterior(chi2_w, box_w)
    jb = np.unravel_index(np.argmin(chi2_w), chi2_w.shape)
    wres = {
        "best": {"Om": float(omw[jb[0]]), "w": float(ww[jb[1]]), "chi2": float(chi2_w[jb])},
        "w": lk.summarize_1d(ww, post_w.sum(0)),
        "Om": lk.summarize_1d(omw, post_w.sum(1)),
    }
    # How far is w = -1 from the best fit, in the 2D sense?
    chi2_lambda = float(np.min(chi2_w[:, np.argmin(np.abs(ww + 1.0))]))
    wres["delta_chi2_w_minus_1"] = chi2_lambda - wres["best"]["chi2"]
    results["flat_wcdm"] = wres
    art["wcdm"] = {"Om": omw, "w": ww, "chi2": chi2_w, "post": post_w,
                   "post_unboxed": lk.grid_posterior(chi2_w)}
    _log(f"flat wCDM: w = {wres['w']['mean']:.3f} +/- {wres['w']['std']:.3f}", t0)

    # --- 5. H0 from the Cepheid-calibrated supernovae ---------------------------
    fit = lk.fit_calibrated(cal, om1)
    h0_grid = np.linspace(66.0, 81.0, 3001)
    h0_dens, om_wts = lk.h0_posterior(fit, h0_grid)
    ic = int(np.argmin(fit.chi2))
    # Monte Carlo: draw Om from its posterior, then (M, 5 log10 H0) from their
    # Gaussian at that Om, and push every draw through the derived quantities.
    n_draw = 200_000
    om_c = rng.choice(om1, size=n_draw, p=om_wts)
    idx = np.searchsorted(om1, om_c)
    beta = np.column_stack([fit.M[idx], fit.h5[idx]]) + rng.multivariate_normal(
        [0.0, 0.0], fit.cov_beta, size=n_draw)
    h0_draws = 10.0 ** (beta[:, 1] / 5.0)
    ol_c = 1.0 - om_c
    age_draws = (cosmo.hubble_time_gyr(h0_draws) * 2.0 / (3.0 * np.sqrt(ol_c))
                 * np.arcsinh(np.sqrt(ol_c / om_c)))
    h0_mean = float(np.trapezoid(h0_grid * h0_dens, h0_grid))
    h0_std = float(np.sqrt(np.trapezoid((h0_grid - h0_mean) ** 2 * h0_dens, h0_grid)))
    tension = (h0_mean - PLANCK18_H0[0]) / np.hypot(h0_std, PLANCK18_H0[1])
    calib = {
        "n_used": int(cal.n),
        "H0": lk.summarize_1d(h0_grid, h0_dens),
        "M": lk.summarize_samples(beta[:, 0]),
        "Om": lk.summarize_1d(om1, om_wts),
        "corr_M_H0": float(np.corrcoef(beta[:, 0], beta[:, 1])[0, 1]),
        "best": {"Om": float(om1[ic]), "M": float(fit.M[ic]),
                 "H0": float(10 ** (fit.h5[ic] / 5)), "chi2": float(fit.chi2[ic])},
        "age_gyr": lk.summarize_samples(age_draws),
        "tension_with_planck_sigma": float(tension),
        "planck_H0": PLANCK18_H0,
        "sh0es_r22_H0": SH0ES_R22_H0,
    }
    best_H0 = calib["best"]["H0"]
    best_M = calib["best"]["M"]
    zt_best = float(cosmo.transition_redshift(calib["best"]["Om"], 1 - calib["best"]["Om"]))
    calib["transition"] = {
        "z": zt_best,
        "lookback_gyr": cosmo.lookback_time_gyr(zt_best, calib["best"]["Om"],
                                                1 - calib["best"]["Om"], H0=best_H0),
    }
    calib["age_best_gyr"] = cosmo.age_gyr(calib["best"]["Om"], 1 - calib["best"]["Om"],
                                          H0=best_H0)
    calib["age_eds_same_H0_gyr"] = cosmo.age_gyr(1.0, 0.0, H0=best_H0)
    calib["age_empty_same_H0_gyr"] = cosmo.age_gyr(0.0, 0.0, H0=best_H0)
    calib["lookback_farthest_sn_gyr"] = cosmo.lookback_time_gyr(
        results["dataset"]["z_max"], calib["best"]["Om"], 1 - calib["best"]["Om"], H0=best_H0)

    # Per-calibrator absolute magnitudes, M_i = m_i - mu_Cepheid,i: the second
    # rung of the ladder made visible.
    cm = cal.is_calibrator
    cal_err = np.sqrt(np.diag(cal.cov))[cm]
    calib["calibrator_M_scatter"] = float(np.std(cal.m[cm] - cal.ceph_mu[cm]))
    art["calibrated"] = {"fit": fit, "h0_grid": h0_grid, "h0_dens": h0_dens,
                         "om_wts": om_wts, "h0_draws": h0_draws,
                         "cal_M": cal.m[cm] - cal.ceph_mu[cm], "cal_err": cal_err,
                         "cal_mu": cal.ceph_mu[cm], "cal_cid": table["CID"].to_numpy()[cal.rows][cm],
                         "cal_z": cal.z_hd[cm]}

    # Stat-only covariance for the same fit: how much of the H0 error is systematic?
    cal_stat = lk.build_sample(table, cov["stat"], name="calibrated, stat only",
                               calibrators=True)
    fit_stat = lk.fit_calibrated(cal_stat, om1)
    dens_stat, _ = lk.h0_posterior(fit_stat, h0_grid)
    calib["H0_stat_only"] = lk.summarize_1d(h0_grid, dens_stat)
    results["calibrated"] = calib
    _log(f"H0 = {calib['H0']['mean']:.2f} +/- {calib['H0']['std']:.2f} km/s/Mpc, "
         f"{tension:.1f} sigma from Planck", t0)

    # --- 6a. Cross-check: SH0ES's own Hubble-flow sample, model-free distances --
    hf_rows = np.flatnonzero(table["IS_CALIBRATOR"].to_numpy().astype(bool)
                             | (table["USED_IN_SH0ES_HF"].to_numpy() == 1))
    shoes = lk.build_sample(table, cov["statsys"], name="SH0ES Hubble flow",
                            calibrators=True, rows=hf_rows)
    kin = lk.fit_calibrated(shoes, np.array([np.nan]), dl_shape_fn=lk.kinematic_dl)
    kin_H0 = 10 ** (kin.h5[0] / 5)
    kin_err = kin_H0 * np.log(10) / 5 * np.sqrt(kin.cov_beta[1, 1])
    z_hf = shoes.z_hd[shoes.hubble_flow]
    results["cross_check_sh0es_hf"] = {
        "n_hubble_flow": int(shoes.hubble_flow.sum()),
        "z_range": [float(z_hf.min()), float(z_hf.max())],
        "window": list(Z_SH0ES_HUBBLE_FLOW),
        "H0": float(kin_H0), "H0_err": float(kin_err), "M": float(kin.M[0]),
    }
    _log(f"cross-check, SH0ES Hubble flow only: H0 = {kin_H0:.2f} +/- {kin_err:.2f}", t0)

    # --- Arrays for the Hubble diagram figures -----------------------------------
    t_all = table
    z_all = t_all["zHD"].to_numpy()
    zhel_all = t_all["zHEL"].to_numpy()
    mu_obs = t_all["m_b_corr"].to_numpy() - best_M
    err_cov = np.sqrt(np.diag(cov["statsys"]))
    ref_shape = cosmo.dl_milne(z_all) * (1 + zhel_all) / (1 + z_all)
    mu_milne = cosmo.distance_modulus(ref_shape * cosmo.hubble_distance_mpc(best_H0))
    art["hubble"] = {
        "z": z_all, "z_hel": zhel_all, "mu": mu_obs, "err": err_cov,
        "resid_milne": mu_obs - mu_milne,
        "group": survey_group(t_all["IDSURVEY"].to_numpy()),
        "is_cal": t_all["IS_CALIBRATOR"].to_numpy().astype(bool),
        "used": z_all > Z_MIN_COSMOLOGY,
        "H0": best_H0, "M": best_M, "Om": calib["best"]["Om"],
    }
    # Binned residuals relative to the empty universe, compressed with the full
    # covariance. Display only.
    sn_resid = art["hubble"]["resid_milne"][sn.rows]
    edges = _bin_edges(sn.z_hd, Z_MIN_COSMOLOGY, 2.3, 18, min_count=8)
    art["binned"] = lk.binned_residuals(sn, sn_resid, edges)
    art["binned"]["edges"] = edges
    # Zero points: A = M + 5 log10(c/H0) + 25 of the data's reference vs each model's own.
    a_ref = best_M + 5 * np.log10(cosmo.hubble_distance_mpc(best_H0)) + 25
    art["model_offsets"] = {k: v["offset"] - a_ref for k, v in comparison.items()}
    # Deceleration-to-acceleration: how much fainter than "no dark energy" at z = 0.5?
    zz = np.array([0.5])
    dl_best = cosmo.luminosity_distance_shape(zz, zz, calib["best"]["Om"], 1 - calib["best"]["Om"])[0, 0]
    dl_open = cosmo.luminosity_distance_shape(zz, zz, 0.3, 0.0)[0, 0]
    results["farther_than_open_at_z05_percent"] = float(100 * (dl_best / dl_open - 1))

    # --- 6b. Robustness: do the conclusions hinge on any one redshift range? -----
    results["robustness"] = robustness(table, cov["statsys"], om1)
    _log("robustness cuts done", t0)

    # --- 6c. Validation against the published chains -----------------------------
    if validate:
        results["validation"] = validate_against_chains(results, art)
        _log("validated against the Pantheon+ chains", t0)

    results["published"] = {k: list(v) for k, v in PUBLISHED_PANTHEON.items()}
    results["planck18"] = {"H0": list(PLANCK18_H0), "Om": list(PLANCK18_OM)}
    if MANIFEST_PATH.exists():
        results["provenance"] = json.loads(MANIFEST_PATH.read_text())
    results["runtime_s"] = round(time.time() - t0, 1)
    return results, art


def robustness(table, cov: np.ndarray, om_grid: np.ndarray) -> dict:
    """Refit flat LambdaCDM (and H0) with slices of the data removed.

    The two binned points furthest from the best fit in figure 02 are the nearest
    bin (0.010 < z < 0.0135) and the z ~ 1 bin (0.93 < z < 1.26), 2.7 and 2.9
    sigma high. If either were driving the answer, dropping it would move the
    result by a large fraction of its error bar.
    """
    z = table["zHD"].to_numpy()
    is_cal = table["IS_CALIBRATOR"].to_numpy().astype(bool)
    hubble_flow = z > Z_MIN_COSMOLOGY
    cuts = {
        "baseline": hubble_flow,
        "drop 0.010 < z < 0.0135 (nearest bin)": hubble_flow & (z >= 0.0135),
        "drop 0.93 < z < 1.26 (the z ~ 1 bin)": hubble_flow & ~((z > 0.93) & (z < 1.26)),
        "drop everything beyond z = 1": hubble_flow & (z <= 1.0),
        "only z > 0.023 (the SH0ES Hubble-flow cut)": z > 0.023,
    }
    out = {}
    for label, keep in cuts.items():
        # As in the published supernova-only likelihood, calibrators above the
        # cut stay in as ordinary Hubble-flow points.
        sn = lk.build_sample(table, cov, name=label, rows=np.flatnonzero(keep))
        chi2, _ = lk.scan(sn, om_grid, 1.0 - om_grid)
        with_cal = lk.build_sample(table, cov, name=label, calibrators=True,
                                   rows=np.flatnonzero(keep | is_cal))
        fit = lk.fit_calibrated(with_cal, om_grid)
        h0_grid = np.linspace(66.0, 81.0, 1501)
        dens, _ = lk.h0_posterior(fit, h0_grid)
        out[label] = {"n": int(sn.n), "Om": lk.summarize_1d(om_grid, lk.grid_posterior(chi2)),
                      "H0": lk.summarize_1d(h0_grid, dens)}
    return out


def validate_against_chains(results: dict, art: dict) -> dict:
    """Our grid posteriors vs the Pantheon+ team's MCMC chains, parameter by parameter."""
    out = {}
    pairs = [
        ("flat_lcdm", "Om", results["flat_lcdm"]["statsys"]["Om"]),
        ("lcdm", "Om", results["lcdm_plane"]["Om"]),
        ("lcdm", "OL", results["lcdm_plane"]["OL"]),
        ("flat_wcdm", "Om", results["flat_wcdm"]["Om"]),
        ("flat_wcdm", "w", results["flat_wcdm"]["w"]),
        ("shoes_flat_lcdm", "H0", results["calibrated"]["H0"]),
        ("shoes_flat_lcdm", "M", results["calibrated"]["M"]),
        ("shoes_flat_lcdm", "Om", results["calibrated"]["Om"]),
    ]
    chains = {name: load_chain(name) for name in CHAIN_FILES}
    for chain_name, par, ours in pairs:
        theirs = lk.summarize_samples(chains[chain_name][par])
        diff = ours["mean"] - theirs["mean"]
        out[f"{chain_name}:{par}"] = {
            "ours": ours, "chain": theirs,
            "diff_in_sigma": float(diff / theirs["std"]),
            "width_ratio": float(ours["std"] / theirs["std"]),
        }
    art["chains"] = chains
    return out


def save_results(results: dict) -> None:
    def _default(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))

    RESULTS_PATH.write_text(json.dumps(results, indent=2, default=_default))


if __name__ == "__main__":
    res, _ = run()
    save_results(res)
    print(json.dumps({k: res[k] for k in ("flat_lcdm", "calibrated")}, indent=1, default=float)[:3000])
