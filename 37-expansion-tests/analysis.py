"""Run the three tests and collect every number the README quotes.

``run()`` returns (summary, artifacts): the summary is plain numbers (saved as
data/results.json), the artifacts are the arrays and tables the figures need.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

import tcmb
import timedilation
import tolman
from config import (BLONDIN08_B, GOLDHABER01_B, KLIMENKO20_T0, KOTANI26_BETA, NOTERDAEME11_BETA,
                    RESULTS_PATH, RIECHERS22_BETA, T0_FIRAS, WHITE24)


def _sigma_from(b, err, ref):
    return abs(b - ref) / err


def summarise_time_dilation(td):
    w, t = td["white"], td["template"]
    avg = w["average"]
    s = {
        "sample": td["sample"],
        "white_replication": {
            "b": avg["width"]["b"], "b_err_formal": avg["width"]["err"],
            "b_err_bootstrap": avg["width"]["bootstrap_err"], "chi2_nu": avg["width"]["chi2_nu"],
            "n_sne": avg["width"]["n"], "linear": avg["width"]["linear"],
            "scatter": avg["width"]["scatter"],
            "b_nelder_mead_from_1": avg["width_from_1"]["b"],
            "b_global_grid_minimum": avg["width_grid"]["b"],
            "per_band": w["per_band"], "comparison_with_white_widths": w["comparison"],
            "white_widths_refit": {k: w["white_refit"][k] for k in ("b", "err", "bootstrap_err",
                                                                    "chi2_nu", "n")},
            "null_test_b": w["null_test"]["b"], "null_test_b_err": w["null_test"]["err"],
            "dispersion_b_min_riz": w["dispersion"]["riz_b_min"],
            "dispersion_b_min_riz_white": w["dispersion"]["white_riz_b_min"],
            "dispersion_b_min_per_band": {b: v["b_min"] for b, v in w["dispersion"].items()
                                          if isinstance(v, dict)},
        },
        "global_fit": {
            "b": t["b"], "b_err_curvature": t["b_err_curvature"],
            "b_err_bootstrap": t["b_err_bootstrap"], "n_bootstrap": int(len(t["bootstrap"])),
            "n_sne": t["n_sne"], "n_points": t["n_points"], "n_curves": t["n_curves"],
            "delta_chi2_tired_light": float(t["scan"]["chi2"][0] - t["chi2_best"]),
            "delta_chi2_b1": float(np.interp(1.0, t["scan"]["b"], t["scan"]["chi2"]) - t["chi2_best"]),
            "sigma_from_b1": _sigma_from(t["b"], t["b_err_bootstrap"], 1.0),
            "sigma_from_b0": _sigma_from(t["b"], t["b_err_bootstrap"], 0.0),
            "bins": {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                     for k, v in t["bins"].items()},
            "variants": t["variants"],
        },
        "published": {"white24": WHITE24, "goldhaber01": GOLDHABER01_B, "blondin08": BLONDIN08_B},
    }
    inj = t.get("injections")
    if inj is not None:
        g = inj.groupby("b_true")["b"]
        s["global_fit"]["injections"] = {float(k): {"mean_b": float(v.mean()), "std_b": float(v.std()),
                                                    "n": int(v.size)} for k, v in g}
    return s


def summarise_tcmb(tc):
    keep = ("n_baseline", "n_upper_limits", "baseline", "baseline_models", "t0_highz", "by_method",
            "with_muller13", "noterdaeme11", "riechers22_all_rows", "riechers22_one_per_absorber")
    s = {k: tc[k] for k in keep}
    s["published"] = {"noterdaeme11_all": NOTERDAEME11_BETA,
                      "noterdaeme11_table": {"SZ": (0.040, 0.079), "SZ+C": (0.029, 0.053),
                                             "SZ+CO": (-0.012, 0.029), "all": (-0.007, 0.027)},
                      "riechers22": RIECHERS22_BETA, "kotani26": KOTANI26_BETA,
                      "klimenko20_T0": KLIMENKO20_T0, "firas_T0": T0_FIRAS}
    b = tc["baseline"]
    err = max(b["plus"], b["minus"])
    s["tired_light_sigma"] = (1.0 - b["beta"]) / err
    return s


def summarise_tolman(to):
    sig = to["signal"]
    return {"n_local": to["n_local"], "n_distant": to["n_distant"],
            "results": to["results"], "published": to["published"],
            "local_fit_check": {str(k): v for k, v in to["local_fit_check"].items()},
            "lookback_gyr": to["lookback_gyr"],
            "per_cluster_q0_half": sig[sig["model"] == "q0=0.5"].round(4).to_dict("records")}


def run(*, refresh=False):
    t0 = time.time()
    td = timedilation.run(refresh=refresh)
    print("CMB temperature")
    tc = tcmb.run()
    print("Tolman test")
    to = tolman.run()
    summary = {"time_dilation": summarise_time_dilation(td), "cmb_temperature": summarise_tcmb(tc),
               "tolman": summarise_tolman(to), "runtime_s": time.time() - t0}
    return summary, (td, tc, to)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        return None if not np.isfinite(x) else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, pd.DataFrame):
        return _jsonable(x.to_dict("records"))
    return x


def save_results(summary):
    RESULTS_PATH.write_text(json.dumps(_jsonable(summary), indent=2))
