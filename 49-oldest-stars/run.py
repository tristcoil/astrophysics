"""Run every fit, write data/results.json, and draw every figure.

    python run.py                 # analysis + figures (~35 minutes on 6 cores)
    python run.py --figures-only  # redraw from the cached analysis
    python run.py --no-figures    # numbers only

Needs the cache from ``python fetch_data.py``; after that it runs offline. The
first run also builds the isochrone grids (MIST streams its 1.1 GB tarball,
~2-3 minutes) and the cleaned photometry; both are cached in data/derived/.
Set B6_WORKERS to change the number of worker processes (default 6).
"""

from __future__ import annotations

import argparse
import json
import pickle
import time
import warnings

import numpy as np

import config as C

CACHE = C.DATA_DIR / "analysis_cache.pkl"


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (np.floating, float)):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def _pm(s: dict, nd: int = 2) -> str:
    return f"{s['median']:.{nd}f} +{s['plus']:.{nd}f} -{s['minus']:.{nd}f}"


def summary(r: dict) -> None:
    print("\n=== Ages of the oldest stars: HUGS turnoffs vs three isochrone sets ===")
    for name, c in r["clusters"].items():
        inp = c["inputs"]
        print(f"\n{name}: [Fe/H] = {inp['feh']:+.2f}, E(B-V) = {inp['ebv']:.2f}, "
              f"(m-M)0 = {inp['dm']:.3f} +- {inp['dm_err']:.3f} ({inp['dist_kpc']:.2f} kpc)")
        for model, m in c["models"].items():
            b = ", ".join(f"{k} {v:.2f}" for k, v in m["budget"].items())
            print(f"  {model:5s} age {_pm(m['mc'])} Gyr  (nominal {m['nominal_age']:.2f}, "
                  f"dc {m['nominal_dc']:+.3f}, chi2 {m['nominal_chi2']:.0f}/{m['n_ridge_points']})")
            print(f"        budget: {b}")
            for k, v in m["variants"].items():
                print(f"        {k}: {_pm(v)}")
        comb = c["combined"]
        print(f"  combined: {comb['age']:.2f} +- {comb['total_err']:.2f} Gyr "
              f"(inputs {comb['inputs_err']:.2f}, models {comb['models_err']:.2f})")
        he = c["helium_test"]
        print(f"  helium Y = 0.275 vs 0.247 (BaSTI): {he['shift']:+.2f} Gyr")
    if "cosmology" in r:
        co = r["cosmology"]
        print(f"\nmetal-poor clusters: t_GC = {co['t_gc']:.2f} +- {co['t_gc_err']:.2f} Gyr;"
              f" + formation delay {co['delay_gyr']['median']:.2f} Gyr -> universe at least "
              f"{co['t_universe_min']['median']:.2f} Gyr (5th percentile {co['t_universe_min']['p5']:.2f})")
        for k, t in co["tensions"].items():
            print(f"  {k}: universe {t['universe_age']:.2f} Gyr, older than the clusters by "
                  f"{t['difference']:+.2f} Gyr ({t['sigma']:+.1f} sigma)")
        for k, lim in co["h0_limits"].items():
            print(f"  H0 < {lim['H0_max_95']:.1f} km/s/Mpc (95%) for {k} = {lim['Om']:.3f}")
    if "ying_comparison" in r:
        print("M92 at Ying et al.'s distance modulus (14.72), central inputs:")
        for model, row in r["ying_comparison"].items():
            print(f"  {model:5s} " + "  ".join(f"{k}: {v:.2f}" for k, v in row.items()))
    if "relative_ages" in r:
        ra = r["relative_ages"]
        print("47 Tuc younger than the metal-poor mean by: " + ", ".join(
            f"{m} {v['median']:.2f} +- {v['sigma']:.2f} (P = {v['p_47tuc_younger']:.3f})"
            for m, v in ra.items() if isinstance(v, dict))
              + f"; mean {ra['mean_over_sets']:.2f}, spread {ra['set_spread']:.2f} Gyr")
    if "chromosome_maps" in r:
        for name, ch in r["chromosome_maps"].items():
            lit = ch["milone2017"].get("f1G", [np.nan, np.nan])
            print(f"chromosome map {name}: 1G fraction {ch['f1G']:.3f} +- {ch['f1G_err']:.3f}"
                  f" (Milone+17 {lit[0]:.3f} +- {lit[1]:.3f}), {ch['n_rgb']} giants")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figures-only", action="store_true")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    t0 = time.time()
    if args.figures_only and CACHE.exists():
        results, art = pickle.loads(CACHE.read_bytes())
        # cheap summaries are recomputed, so edits to the reference values in
        # config.py show up without rerunning the Monte Carlo
        import analysis

        analysis.post_process(results, art)
        C.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
        CACHE.write_bytes(pickle.dumps((results, art)))
    else:
        import analysis

        results, art = analysis.run_all()
        C.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
        CACHE.write_bytes(pickle.dumps((results, art)))
        print(f"analysis: {time.time() - t0:.0f} s")
    summary(results)
    if not args.no_figures:
        import figures

        written = figures.make_all(results, art)
        print(f"\n{len(written)} figure files written to figures/")
    print(f"total: {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
