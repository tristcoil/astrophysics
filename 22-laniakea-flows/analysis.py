"""Every measurement of the project, in order, with results gathered for figures and JSON.

Expensive steps cache their products in ``data/derived/``: the Lambda-CDM
covariance of the bulk-flow fits (85 s per catalogue), the 2MRS reconstruction
(1 min), the beta* profile likelihood (5 min), and the basins of attraction.
Delete that folder to recompute everything from the raw catalogues.
"""

from __future__ import annotations

import json
import pickle

import numpy as np
import pandas as pd

import basins
import bulkflow
import catalogs
import cf4fix
import config
import cosmology as cg
import density
import flows
import gravity
import hubble


def _cached(name: str, fn, *, refresh: bool = False):
    path = config.DERIVED_DIR / f"{name}.pkl"
    if path.exists() and not refresh:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    obj = fn()
    with open(path, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return obj


def frames() -> dict:
    lg = catalogs.local_group_motion()
    return {"sun_cmb": vars(config.SUN_CMB), "sun_lg": vars(config.SUN_LG),
            "lg_cmb_ours": {"speed": lg["speed"], "l": lg["l"], "b": lg["b"]},
            "lg_cmb_planck_table": vars(config.LG_CMB),
            "sigma_v_1d_lcdm": cg.sigma_v_1d(), "growth_rate_f": cg.growth_rate(config.PLANCK18.Om0)}


def hubble_section() -> dict:
    return {v: hubble.hubble_constant(v) for v in ("published", "corrected")}


def estimator_section() -> dict:
    return flows.estimator_check()


def sdss_section() -> dict:
    return {"evidence": cf4fix.evidence(), "zero_point": cf4fix.zero_point(),
            "n_groups_affected": int(cf4fix.corrected_groups()["sdss_corrected"].sum())}


def bulk_section(verbose: bool = True) -> dict:
    def run():
        out = {}
        for version in ("corrected", "published"):
            if verbose:
                print(f"  bulk flows, {version} catalogue")
            out[version] = bulkflow.analyse("all", version=version, verbose=verbose)
        for which in ("TF", "SN", "FP north", "FP south"):
            if verbose:
                print(f"  bulk flows, {which} only")
            out[which] = bulkflow.analyse(which, version="corrected", verbose=verbose)
        out["ideal"] = bulkflow.ideal_sphere_expectation()
        return out

    return _cached("bulk_flows", run)


def comparison_with_nusser_tully(bulk: dict) -> dict:
    nt = config.PUBLISHED["nusser_tully2026"]
    out = {}
    for version in ("corrected", "published"):
        rows = []
        for k, ref in nt["per_bin_sigma"].items():
            d = bulk[version]["bins"].get(k)
            if d is None:
                continue
            rows.append({"bin": k, "ours": d["sigma_dev"], "nt": list(ref),
                         "diff": (np.array(d["sigma_dev"]) - np.array(ref)).tolist()})
        diffs = np.concatenate([r["diff"] for r in rows])
        out[version] = {"rows": rows, "rms_diff_sigma": float(np.sqrt(np.mean(diffs**2))), "n": int(len(diffs))}
    return out


def reconstruction(verbose: bool = True) -> dict:
    return density.reconstruct(verbose=verbose)


def dipole_section(rec: dict) -> dict:
    lf = rec["lf"]
    out = {"lg": gravity.flux_dipole(lf, frame="lg"), "cmb": gravity.flux_dipole(lf, frame="cmb")}
    for fr in ("lg", "cmb"):
        out[f"beta_{fr}_60"] = gravity.beta_from_dipole(out[fr], 60.0)
        out[f"beta_{fr}_130"] = gravity.beta_from_dipole(out[fr], 130.0)
    out["rho_L_isophotal"] = out["lg"]["rho_L"]
    # Isophotal K20 magnitudes miss ~0.2 mag of light (Kochanek et al. 2001).
    out["rho_L_total_equiv"] = out["lg"]["rho_L"] * 10 ** (0.4 * 0.20)
    out["published"] = config.PUBLISHED["erdogdu2006"]
    return out


def velocity_section(rec: dict, verbose: bool = True) -> dict:
    def run():
        res = gravity.velocity_comparison(rec, verbose=verbose)
        per = res.pop("_per_object")
        best = res.pop("_best")
        return {"res": res, "per_object": per, "vext": best["vext"], "beta_grid": best["beta"]}

    return _cached("velocity_comparison", run)


def lg_decomposition(rec: dict, vc: dict) -> list:
    return gravity.shell_decomposition(rec, vc["beta_grid"], np.asarray(vc["vext"]))


def basin_section(rec: dict, vc: dict) -> dict:
    def run():
        b = basins.all_basins(rec, vc["beta_grid"], np.asarray(vc["vext"]))
        return {k: {"report": v["report"], "labels": v["labels"].astype(np.int32),
                    "psi": v["psi"].astype(np.float32)} for k, v in b.items()}

    return _cached("basins", run)


def run_all(verbose: bool = True) -> tuple[dict, dict]:
    """Every result (JSON-able) and the arrays the figures need."""
    res, art = {}, {}
    res["frames"] = frames()
    if verbose:
        print("Hubble diagram")
    res["hubble"] = hubble_section()
    art["hubble_points"] = hubble.hubble_parameters(version="corrected")
    if verbose:
        print("Peculiar-velocity estimators and CF4's published columns")
    res["estimators"] = estimator_section()
    if verbose:
        print("The SDSS distances in CF4")
    res["sdss_sign"] = sdss_section()
    if verbose:
        print("Maps")
    g = flows.group_velocities("corrected")
    art["groups"] = g
    art["groups_published"] = flows.group_velocities("published")
    art["sky_maps"] = flows.sky_maps(g)
    res["sky_extremes"] = flows.sky_dipole_check(g)
    art["landmarks"] = flows.landmarks()
    if verbose:
        print("Bulk flows")
    bulk = bulk_section(verbose)
    res["bulk_flow"] = bulk
    res["nusser_tully_comparison"] = comparison_with_nusser_tully(bulk)
    if verbose:
        print("2MRS reconstruction")
    rec = reconstruction(verbose)
    art["rec"] = rec
    res["luminosity_function"] = rec["lf"]
    res["reconstruction"] = {"betas": rec["betas"].tolist(), "v_lg_internal_by_beta": np.linalg.norm(rec["v_lg"], axis=1).tolist(),
                             "n_galaxies": int(rec["r"].shape[1]), "n_clones": int(rec["n_clones"])}
    if verbose:
        print("Flux dipole")
    res["flux_dipole"] = dipole_section(rec)
    if verbose:
        print("Gravity against the measured flows")
    vc = velocity_section(rec, verbose)
    res["velocity_comparison"] = vc["res"]
    art["vc"] = vc
    res["lg_decomposition"] = lg_decomposition(rec, vc)
    res["sdss_sign"]["gravity_test"] = _cached("sdss_gravity_test",
                                               lambda: gravity.sdss_sign_test(rec, vc["beta_grid"], np.asarray(vc["vext"])))
    if verbose:
        print("Basins of attraction")
    bas = basin_section(rec, vc)
    res["basins"] = {k: v["report"] for k, v in bas.items()}
    art["basins"] = bas
    res["published"] = config.PUBLISHED
    return res, art


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, pd.DataFrame):
        return None
    return x


def save_results(res: dict) -> None:
    config.RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=1))
