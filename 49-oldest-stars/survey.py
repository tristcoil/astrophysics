"""Extension: the same age machinery on every globular cluster in HUGS.

    python survey.py            # needs: python fetch_data.py --all (+3 GB)

For each of the 56 globular clusters in the survey (the 57th target, NGC 6791,
is an old *open* cluster), with the same inputs as the core analysis (Harris
[Fe/H] and E(B-V), Baumgardt & Vasiliev distance + 0.025 mag floor, [alpha/Fe] =
+0.4 +- 0.1):
* the cleaned, dereddened CMD and its ridge line (photometry.py);
* a MIST age with a 200-draw Monte Carlo over the inputs;
* single fits at the central inputs with BaSTI and DSED, for the model spread.

Clusters are flagged when the result should not be trusted at face value:
* E(B-V) > 0.35: the reddening (and its patchiness) is large, and the
  extinction law's uncertainty starts to matter;
* "Type II" clusters of Milone et al. (2017), whose subgiant branch is split
  into populations with different iron or C+N+O (a single isochrone cannot
  describe them);
* fits that hit the edge of the age grid, need a colour shift above 0.06 mag,
  or fit badly (chi^2 per ridge point > 3).

Results go to data/survey.json; run.py --figures-only draws them if present.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
import zlib

import numpy as np

import agefit as A
import catalog as K_
import config as C
import literature as L
import photometry as P
from isochrones import load_grid

SURVEY_PATH = C.DATA_DIR / "survey.json"
N_MC = 200
N_WORKERS = int(os.environ.get("B6_WORKERS", "6"))
#: Above this reddening, the extinction law matters (M4, E(B-V) = 0.35, has R_V ~ 3.6-4
#: rather than 3.1; Hendricks et al. 2012) and differential reddening is large.
EBV_FLAG = 0.30
_GRIDS: dict = {}


def survey_ids() -> list[str]:
    ids = sorted(p.name.split("_")[4] for p in C.HUGS_DIR.glob("hlsp_hugs_*_catalog-meth1.txt"))
    return [i for i in ids if i != "ngc6791"]


def fit_one(hid: str) -> dict:
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass
    t0 = time.time()
    try:
        inputs = K_.cluster_inputs(hid)
        prep = P.prepare_hugs(hid, inputs)
        rng = np.random.default_rng([C.RANDOM_SEED, zlib.crc32(hid.encode())])  # stable across runs
        d0 = A.Draw(inputs["feh"], inputs["alpha"], inputs["dm"], inputs["ebv"])
        out = {"id": hid, "inputs": {k: v for k, v in inputs.items() if k not in ("bv21", "harris")},
               "photometry": prep["meta"], "nominal": {}}
        for model in ("MIST", "BaSTI", "DSED"):
            f = A.ClusterFitter(_GRIDS[model], prep["ridge"], feh0=inputs["feh"])
            out["nominal"][model] = f.best(d0)
        f = A.ClusterFitter(_GRIDS["MIST"], prep["ridge"], feh0=inputs["feh"])
        mc = A.monte_carlo(f, inputs, N_MC, rng)
        out["mc_mist"] = mc["summary"]
        out["mc_edge_fraction"] = mc["n_edge"] / N_MC
        out["seconds"] = time.time() - t0
        return out
    except Exception as exc:  # report and move on: one odd cluster must not stop the survey
        return {"id": hid, "error": f"{type(exc).__name__}: {exc}", "seconds": time.time() - t0}


def flags(entry: dict, milone: dict) -> list[str]:
    out = []
    if "error" in entry:
        return ["failed: " + entry["error"]]
    if entry["inputs"]["ebv"] >= EBV_FLAG:
        out.append(f"high reddening E(B-V) = {entry['inputs']['ebv']:.2f}")
    if milone.get(entry["id"], {}).get("type") == "II":
        out.append("Type II (split subgiant branch)")
    m = entry["nominal"]["MIST"]
    if m["edge"] or entry["mc_edge_fraction"] > 0.1:
        out.append("age at the edge of the grid")
    if abs(m["dc"]) > 0.06:
        out.append(f"colour shift {m['dc']:+.3f}")
    if m["chi2"] / max(m["n"], 1) > 3:
        out.append(f"poor fit (chi2/point {m['chi2'] / m['n']:.1f})")
    return out


def run() -> dict:
    for model in ("MIST", "BaSTI", "DSED"):
        _GRIDS[model] = load_grid(model)
    ids = survey_ids()
    print(f"survey: {len(ids)} clusters on {N_WORKERS} workers")
    t0 = time.time()
    entries = []
    with mp.get_context("fork").Pool(N_WORKERS) as pool:
        for k, e in enumerate(pool.imap_unordered(fit_one, ids)):
            entries.append(e)
            msg = e.get("error") or f"MIST {e['nominal']['MIST']['age']:.2f} Gyr"
            print(f"  {k + 1:2d}/{len(ids)} {e['id']}: {msg} ({e['seconds']:.0f} s)", flush=True)
    out = finish(entries)
    print(f"survey done in {time.time() - t0:.0f} s -> {SURVEY_PATH}")
    return out


def finish(entries: list) -> dict:
    """Flags, origins, combined ages and the summary; writes data/survey.json."""
    milone = L.milone2017()
    prog = L.massari2019()
    vdb, val = L.vandenberg2013(), L.valcin2025()
    out = {}
    for e in sorted(entries, key=lambda x: x["id"]):
        e["flags"] = flags(e, milone)
        e["progenitor"] = prog.get(e["id"])
        e["origin"] = L.origin(e["progenitor"])
        if "nominal" in e:
            ages = np.array([e["nominal"][m]["age"] for m in ("MIST", "BaSTI", "DSED")])
            e["model_spread"] = float(np.std(ages, ddof=1))
            e["age"] = float(np.mean(ages))
            stat = 0.5 * (e["mc_mist"]["plus"] + e["mc_mist"]["minus"])
            e["age_err"] = float(np.hypot(stat, e["model_spread"]))
        e["vandenberg2013"] = vdb.get(e["id"])
        e["valcin2025"] = val.get(e["id"])
        out[e["id"]] = e
    summary = summarize(out)
    SURVEY_PATH.write_text(json.dumps({"clusters": out, "summary": summary}, indent=1, default=_np))
    return {"clusters": out, "summary": summary}


def refresh() -> dict:
    """Recompute flags and summary from the stored fits (no refitting)."""
    stored = json.loads(SURVEY_PATH.read_text())["clusters"]
    return finish(list(stored.values()))


def _np(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


def summarize(out: dict) -> dict:
    good = [e for e in out.values() if "age" in e and not e["flags"]]
    res = {"n_total": len(out), "n_clean": len(good)}

    def branch(sel, lo, hi):
        pts = [(e["inputs"]["feh"], e["age"], e["age_err"]) for e in good if sel(e) and lo <= e["inputs"]["feh"] < hi]
        if len(pts) < 3:
            return None
        f, a, s = map(np.array, zip(*pts))
        w = 1 / s**2
        return {"n": len(pts), "mean_age": float(np.sum(w * a) / np.sum(w)), "err": float(np.sqrt(1 / np.sum(w)))}

    for org in ("in situ", "accreted", "unassociated"):
        for label, lo, hi in (("metal-poor", -3, -1.7), ("intermediate", -1.7, -1.0), ("metal-rich", -1.0, 0)):
            res[f"{org} {label}"] = branch(lambda e, o=org: e["origin"] == o, lo, hi)
    for ref in ("vandenberg2013", "valcin2025"):
        pairs = [(e["age"], e[ref][0]) for e in good if e.get(ref)]
        if pairs:
            a, b = map(np.array, zip(*pairs))
            res[f"vs_{ref}"] = {"n": len(pairs), "median_diff": float(np.median(a - b)),
                                "rms_diff": float(np.sqrt(np.mean((a - b - np.median(a - b)) ** 2))),
                                "pearson_r": float(np.corrcoef(a, b)[0, 1])}
    # In situ against accreted at intermediate metallicity, separately for each
    # isochrone set (central-input fits), so the gap is not an artefact of one set.
    res["in_situ_minus_accreted_intermediate"] = {}
    for model in ("MIST", "BaSTI", "DSED"):
        ages = {o: [e["nominal"][model]["age"] for e in good if e["origin"] == o and -1.7 <= e["inputs"]["feh"] < -1.0]
                for o in ("in situ", "accreted")}
        if all(len(v) >= 3 for v in ages.values()):
            res["in_situ_minus_accreted_intermediate"][model] = float(np.mean(ages["in situ"]) - np.mean(ages["accreted"]))
    # The fitted colour shifts against reddening: a slope means the reddened model
    # colours are off, i.e. the extinction law's E(F606W-F814W)/E(B-V) (MIST's is
    # 1.17 for a turnoff star) or the catalogue reddenings are biased. Baseline:
    # clean clusters (flagged ones include shifts clipped at the +-0.08 limit);
    # the other selections show how much the slope depends on that choice.
    def dc_fit(sel):
        pts = [(e["inputs"]["ebv"], e["nominal"]["MIST"]["dc"]) for e in out.values() if "nominal" in e and sel(e)]
        if len(pts) <= 5:
            return None
        x, y = map(np.array, zip(*pts))
        coef, cov = np.polyfit(x, y, 1, cov=True)
        return {"n": len(pts), "slope": float(coef[0]), "slope_err": float(np.sqrt(cov[0, 0])),
                "intercept": float(coef[1]), "intercept_err": float(np.sqrt(cov[1, 1]))}

    not_type2 = lambda e: "Type II (split subgiant branch)" not in e["flags"]  # noqa: E731
    unclipped = lambda e: abs(e["nominal"]["MIST"]["dc"]) < 0.079  # noqa: E731
    res["dc_vs_ebv"] = dc_fit(lambda e: not e["flags"])
    res["dc_vs_ebv_variants"] = {
        "all but Type II": dc_fit(not_type2),
        "all but Type II, unclipped shifts": dc_fit(lambda e: not_type2(e) and unclipped(e)),
        "same, E(B-V) < 0.30": dc_fit(lambda e: not_type2(e) and unclipped(e) and e["inputs"]["ebv"] < 0.30),
    }
    oldest = sorted(good, key=lambda e: -e["age"])[:5]
    res["oldest"] = [(e["id"], e["age"], e["age_err"]) for e in oldest]
    mp = [e for e in good if e["inputs"]["feh"] < -1.5]
    if mp:
        res["metal_poor_n"] = len(mp)
        res["metal_poor_mean"] = float(np.mean([e["age"] for e in mp]))
        res["metal_poor_std"] = float(np.std([e["age"] for e in mp], ddof=1))
        res["metal_poor_median_err"] = float(np.median([e["age_err"] for e in mp]))
    res["flag_counts"] = {}
    for e in out.values():
        for f in e["flags"]:
            key = f.split(" E(B-V)")[0].split(" (")[0].split(" -")[0].split(" +")[0]
            res["flag_counts"][key] = res["flag_counts"].get(key, 0) + 1
    res["n_flagged"] = sum(1 for e in out.values() if e["flags"])
    res["clean_with_dsed_at_15_gyr"] = [e["id"] for e in good if e["nominal"]["DSED"]["age"] > 14.9]
    return res


if __name__ == "__main__":
    import sys
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    if "--refresh" in sys.argv:
        refresh()
    else:
        run()
