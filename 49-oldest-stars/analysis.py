"""Every number in the README: ages, their error budgets, robustness tests, comparisons.

    import analysis; results, art = analysis.run_all()

``results`` is JSON-ready (written to data/results.json by run.py); ``art``
holds the arrays the figures need. The Monte Carlo fits run in parallel worker
processes (fork), each with one BLAS thread.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import time

import numpy as np

import agefit as A
import catalog as K_
import config as C
import cosmology as CO
import photometry as P
from isochrones import load_grid

#: Monte Carlo sizes. The baseline has every input varying at once; the error
#: budget varies one input at a time; variants change one analysis choice.
N_BASE = C.N_MC
N_BUDGET = 300
N_VARIANT = 300
N_WORKERS = int(os.environ.get("B6_WORKERS", "6"))

MODEL_SETS = list(C.MODEL_SETS)
BUDGET_TERMS = {
    "distance": ("dm",),
    "reddening": ("ebv",),
    "[Fe/H]": ("feh",),
    "[alpha/Fe]": ("alpha",),
    "ridge statistics": ("ridge",),
}

# Filled before the workers fork, so they share the grids and prepared data.
_GRIDS: dict = {}
_PREPARED: dict = {}
_INPUTS: dict = {}


def _fitter(cluster: str, model: str, **kw) -> A.ClusterFitter:
    return A.ClusterFitter(_GRIDS[model], _PREPARED[cluster]["ridge"], feh0=_INPUTS[cluster]["feh"], **kw)


def _job(task: dict) -> dict:
    """One Monte Carlo run (executed in a worker)."""
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass
    rng = np.random.default_rng(task["seed"])
    kw = task.get("fitter", {})
    f = _fitter(task["cluster"], task["model"], **kw)
    if task.get("method") == "sgb_point":
        draws = A.draw_inputs(_INPUTS[task["cluster"]], task["n"], rng, vary=task["vary"])
        ages = np.array([f.sgb_point_age(d) for d in draws])
        return {**task, "ages": ages, "summary": A.summarize(ages)}
    out = A.monte_carlo(f, _INPUTS[task["cluster"]], task["n"], rng, vary=task["vary"],
                        dm_key=task.get("dm_key", "dm"))
    return {**task, **out}


def _tasks() -> list[dict]:
    seeds = np.random.SeedSequence(C.RANDOM_SEED)
    tasks = []
    all_vary = ("dm", "ebv", "feh", "alpha", "ridge")
    for cluster in C.CORE_CLUSTERS:
        for model in MODEL_SETS:
            tasks.append({"kind": "baseline", "cluster": cluster, "model": model, "n": N_BASE,
                          "vary": all_vary})
            for term, vary in BUDGET_TERMS.items():
                tasks.append({"kind": "budget", "term": term, "cluster": cluster, "model": model,
                              "n": N_BUDGET, "vary": vary})
            if "dm_geo" in _INPUTS[cluster]:
                tasks.append({"kind": "geometric_distance", "cluster": cluster, "model": model,
                              "n": N_VARIANT, "vary": all_vary, "dm_key": "dm_geo"})
            tasks.append({"kind": "sgb_point", "cluster": cluster, "model": model, "n": N_VARIANT,
                          "vary": ("dm", "ebv", "feh", "alpha"), "method": "sgb_point"})
            tasks.append({"kind": "no_colour_shift", "cluster": cluster, "model": model,
                          "n": N_VARIANT, "vary": all_vary, "fitter": {"colour_shift": False}})
        # helium: BaSTI with Y = 0.275 instead of 0.247 (alpha-enhanced grid only)
        tasks.append({"kind": "helium", "cluster": cluster, "model": "BaSTI_Y275", "n": N_VARIANT,
                      "vary": ("dm", "ebv", "feh", "ridge")})
        tasks.append({"kind": "helium_ref", "cluster": cluster, "model": "BaSTI", "n": N_VARIANT,
                      "vary": ("dm", "ebv", "feh", "ridge")})
    children = seeds.spawn(len(tasks))
    for t, s in zip(tasks, children):
        t["seed"] = s
    return tasks


def _nominal(cluster: str, model: str, **kw) -> dict:
    inp = _INPUTS[cluster]
    d = A.Draw(inp["feh"], inp["alpha"], inp["dm"], inp["ebv"])
    f = _fitter(cluster, model, **kw)
    b = f.best(d)
    chi2, dc, r = f.chi2_curve(d, ages=f.ages)
    b["chi2_curve"] = chi2
    b["ages"] = f.ages
    b["residuals"] = r[int(np.argmin(chi2))]
    return b


def _region_tests(cluster: str, model: str) -> dict:
    out = {}
    for label, kw in {
        "baseline (MS 2.0, RGB 1.5 mag)": {},
        "main sequence to 1.0 mag below the turnoff": {"ms_depth": 1.0},
        "main sequence to 3.0 mag below the turnoff": {"ms_depth": 3.0},
        "giant branch to 0.5 mag above the subgiant point": {"rgb_height": 0.5},
        "giant branch to 2.5 mag above the subgiant point": {"rgb_height": 2.5},
        "model floor halved": {"floor_colour": 0.005, "floor_mag": 0.01},
        "model floor doubled": {"floor_colour": 0.02, "floor_mag": 0.04},
    }.items():
        out[label] = _nominal(cluster, model, **kw)["age"]
    return out


def run_ages() -> tuple[dict, dict]:
    t0 = time.time()
    for name in MODEL_SETS + ["BaSTI_Y275"]:
        _GRIDS[name] = load_grid(name)
    for cluster in C.CORE_CLUSTERS:
        _INPUTS[cluster] = K_.cluster_inputs(cluster)
        _PREPARED[cluster] = P.prepare(cluster)
    print(f"  grids and photometry ready ({time.time() - t0:.0f} s)")

    tasks = _tasks()
    print(f"  {len(tasks)} Monte Carlo runs on {N_WORKERS} workers")
    ctx = mp.get_context("fork")
    done = []
    with ctx.Pool(N_WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(_job, tasks)):
            done.append(res)
            if (i + 1) % 10 == 0 or i + 1 == len(tasks):
                print(f"    {i + 1}/{len(tasks)} runs ({time.time() - t0:.0f} s)", flush=True)

    results, art = {"clusters": {}}, {"mc": {}, "nominal": {}}
    for cluster in C.CORE_CLUSTERS:
        inp = _INPUTS[cluster]
        cres = {"inputs": {k: v for k, v in inp.items() if k not in ("bv21", "harris")},
                "photometry": _PREPARED[cluster]["meta"], "models": {}}
        for model in MODEL_SETS:
            mine = [r for r in done if r["cluster"] == cluster and r["model"] == model]
            base = next(r for r in mine if r["kind"] == "baseline")
            budget = {r["term"]: r["summary"] for r in mine if r["kind"] == "budget"}
            nominal = _nominal(cluster, model)
            entry = {
                "nominal_age": nominal["age"], "nominal_dc": nominal["dc"],
                "nominal_chi2": nominal["chi2"], "n_ridge_points": nominal["n"],
                "mc": base["summary"], "mc_edge_fraction": base["n_edge"] / base["n"],
                "mc_dc": A.summarize(base["dc"]),
                "budget": {k: 0.5 * (v["plus"] + v["minus"]) for k, v in budget.items()},
                "variants": {},
            }
            for r in mine:
                if r["kind"] in ("geometric_distance", "sgb_point", "no_colour_shift"):
                    entry["variants"][r["kind"]] = r["summary"]
            entry["region_tests"] = _region_tests(cluster, model)
            cres["models"][model] = entry
            art["mc"][(cluster, model)] = {"ages": base["ages"], "draws": base["draws"], "dc": base["dc"]}
            art["nominal"][(cluster, model)] = nominal
        he = next(r for r in done if r["cluster"] == cluster and r["kind"] == "helium")
        he_ref = next(r for r in done if r["cluster"] == cluster and r["kind"] == "helium_ref")
        cres["helium_test"] = {"Y0.247": he_ref["summary"], "Y0.275": he["summary"],
                               "shift": he["summary"]["median"] - he_ref["summary"]["median"]}
        cres["combined"] = combine(cres["models"])
        results["clusters"][cluster] = cres
    results["timing_s"] = time.time() - t0
    return results, art


def combine(models: dict) -> dict:
    """One age per cluster: the mean of the model sets, with their spread as a systematic."""
    med = np.array([m["mc"]["median"] for m in models.values()])
    stat = np.array([0.5 * (m["mc"]["plus"] + m["mc"]["minus"]) for m in models.values()])
    age = float(np.mean(med))
    spread = float(np.std(med, ddof=1)) if len(med) > 1 else 0.0
    inputs_err = float(np.sqrt(np.mean(stat**2)))
    return {"age": age, "inputs_err": inputs_err, "models_err": spread,
            "total_err": float(np.hypot(inputs_err, spread)),
            "range": [float(med.min()), float(med.max())]}


def run_all() -> tuple[dict, dict]:
    results, art = run_ages()
    print("  ages done; cosmology and multiple populations")
    import chromosome as CH

    results["chromosome_maps"], art["chromosome_maps"] = CH.run_all({c: _PREPARED[c] for c in C.CORE_CLUSTERS},
                                                                     _INPUTS)
    art["prepared"] = _PREPARED
    art["inputs"] = _INPUTS
    post_process(results, art)
    return results, art


def post_process(results: dict, art: dict) -> None:
    """Summaries computed from the stored fits in seconds (run.py --figures-only redoes them)."""
    results["literature"] = literature_comparison(results)
    results["cosmology"], art["cosmology"] = CO.age_crisis(results)
    results["ying_comparison"] = ying_comparison(art)
    results["relative_ages"] = relative_ages(art)


#: Distance moduli and reddenings for the like-for-like comparison with Ying et
#: al. (2023): their own distance modulus with their preferred reddening range.
YING_SETTINGS = [("B&V21 distance, Harris E(B-V)", None, None), ("Ying distance, E(B-V) = 0.02", 14.72, 0.02),
                 ("Ying distance, E(B-V) = 0.01", 14.72, 0.01), ("Ying distance, E(B-V) = 0.00", 14.72, 0.0)]


def ying_comparison(art: dict) -> dict:
    """M92 refitted at Ying et al.'s distance modulus (14.72) for each isochrone set (central inputs)."""
    inp, ridge = art["inputs"]["M92"], art["prepared"]["M92"]["ridge"]
    out = {}
    for model in MODEL_SETS:
        f = A.ClusterFitter(load_grid(model), ridge, feh0=inp["feh"])
        row = {}
        for label, dm, ebv in YING_SETTINGS:
            d = A.Draw(inp["feh"], inp["alpha"], inp["dm"] if dm is None else dm, inp["ebv"] if ebv is None else ebv)
            row[label] = f.best(d)["age"]
        out[model] = row
    return out


def relative_ages(art: dict, seed: int = 0) -> dict:
    """How much younger 47 Tuc is than the mean of the three metal-poor clusters, per isochrone set.

    The Monte Carlo draws of different clusters are independent, so randomly
    paired draws give the distribution of the difference.
    """
    rng = np.random.default_rng(seed)
    out = {}
    for model in MODEL_SETS:
        mp = np.mean([rng.permutation(art["mc"][(c, model)]["ages"]) for c in ("M92", "M15", "NGC 6397")], axis=0)
        tuc = rng.permutation(art["mc"][("47 Tuc", model)]["ages"])
        d = mp - tuc
        out[model] = {"median": float(np.median(d)),
                      "sigma": float(0.5 * (np.percentile(d, 84.135) - np.percentile(d, 15.865))),
                      "p_47tuc_younger": float(np.mean(d > 0))}
    meds = [v["median"] for v in out.values()]
    out["mean_over_sets"] = float(np.mean(meds))
    out["set_spread"] = float(np.std(meds, ddof=1))
    return out


def literature_comparison(results: dict) -> dict:
    import literature as L

    out = {}
    vdb = L.vandenberg2013()
    val = L.valcin2025()
    for cluster in C.CORE_CLUSTERS:
        hid = C.CLUSTERS[cluster]["hugs"]
        comb = results["clusters"][cluster]["combined"]
        rows = {"this work": (comb["age"], comb["total_err"])}
        if hid in vdb:
            rows["VandenBerg et al. 2013"] = vdb[hid]
        if hid in val:
            rows["Valcin et al. 2025"] = val[hid]
        if cluster == "M92":
            rows["Ying et al. 2023"] = (C.YING2023_M92["age"], C.YING2023_M92["err"])
        for clock in C.OTHER_CLOCKS.get(cluster, []):
            rows[clock["label"]] = (clock["age"], clock["err"])
        out[cluster] = {k: {"age": v[0], "err": v[1]} for k, v in rows.items()}
    return out
