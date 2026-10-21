"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run

The first run downloads ~3 MB (0.5 MB is kept). Once the cache is warm the
analysis runs offline in about six minutes, most of it Markov chains.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import config as C
import figures
import interactive
from fetch_data import fetch_all

CACHE = C.CACHE_DIR / "analysis.pkl"


def _fit(f: dict, nd=2) -> str:
    return (f"log M = {f['a']:.{nd}f} ± {f['ea']:.{nd}f} + ({f['b']:.{nd}f} ± {f['eb']:.{nd}f}) log(σ/200), "
            f"scatter {f['eps']:.{nd}f} dex, N = {f['n']}" + (f" ({f['n_upper']} limits)" if f.get("n_upper") else ""))


def print_summary(s: dict) -> None:
    f = s["fits"]
    print("\n" + "=" * 78)
    print("RESULTS: the M-σ relation")
    print("=" * 78)
    print("\nOur headline fit (van den Bosch's 230 galaxies, Bayesian, censored limits):")
    print("  " + _fit(f["vdb16"]))
    print("  detections only:      " + _fit(f["vdb16_detections"]))
    print("  limits made sharp:    " + _fit(f["vdb16_sharp_limits"]))
    print("  measured σ as given:  " + _fit(f["vdb16_conditional"]))
    print("\nEllipticals + classical bulges (Kormendy & Ho's 45):")
    print("  " + _fit(f["kh13"]))
    print("\nReproducing the published fits from each paper's own sample and method:")
    for r in s["reproduction"]:
        print(f"  {r['label']:<36} slope {r['b']:.3f} ± {r['eb']:.3f}  (published {r['pub_b']} ± {r['pub_eb']}, "
              f"{r['db_sigma']:+.2f}σ); zero point {r['a']:.3f} ({r['pub_a']:.3f}); scatter {r['eps']:.2f} ({r['pub_eps']})")
    print("\nWhy the slopes differ (McConnell & Ma -> Kormendy & Ho, cumulative):")
    for st in s["slope_budget"]["steps"]:
        print(f"  {st['label']:<46} N = {st['n']:<3} slope {st['b']:.2f} ± {st['eb']:.2f}")
    b = s["bulges"]
    print(f"\nPseudobulges sit {b['pseudo']['mean']:+.2f} ± {b['pseudo']['sem']:.2f} dex off the classical relation "
          f"({100 * b['pseudo']['frac_below']:.0f}% below; KS p = {b['pseudo']['ks_pvalue']:.1e}); "
          f"their own correlation r = {b['pseudo']['pearson_r']:.2f} (p = {b['pseudo']['pearson_p']:.2f})")
    a = s["anchors"]
    print(f"Milky Way (H1 mass {C.MILKY_WAY['M']:.3e}): {a['milky_way']['kh13']['res']:+.2f} dex from the classical "
          f"relation, {a['milky_way']['vdb16']['res']:+.2f} from the all-galaxy one")
    print(f"M87* (EHT 6.5e9): {a['m87']['kh13']['res']:+.2f} dex from the classical relation, "
          f"{a['m87']['vdb16']['res']:+.2f} from the all-galaxy one")
    bm = s["bulge_mass"]
    print(f"\nBlack hole / bulge mass at 1e11 Msun: {bm['ratio_at_1e11_percent_ours']:.2f}% "
          f"(Kormendy & Ho's method on their data: {bm['ratio_at_1e11_percent_eq10']:.2f}%; published 0.49%)")
    sc = s["scales"]
    print(f"Sphere of influence: median {sc['r_infl_pc_median']:.0f} pc = {100 * sc['r_infl_over_re_median']:.2f}% "
          f"of the half-light radius")
    hb = s["host_bias"]["hosts"]
    print(f"ATLAS3D: black hole hosts have σ {hb['mean_dlogsig']:+.3f} ± {hb['sem']:.3f} dex higher than other "
          f"galaxies of the same luminosity (permutation p = {hb['perm_p']:.4f})")
    fb = s["feedback"]
    print(f"King (2003) predicts {fb['king_M_at_200']:.2e} Msun at 200 km/s; ellipticals + classical bulges have "
          f"{fb['ratio_observed_to_king_kh13']:.2f}× that. Energy released / bulge binding energy: median "
          f"{fb['energy_ratio_median']:.0f}×")
    print("Virial factor <f> for reverberation masses, matched to:")
    for c in s["agn"]["calibrations"]:
        print(f"  {c['relation']:<66} f = {c['f']:.1f} ({c['f_lo']:.1f}-{c['f_hi']:.1f})")
    print()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        res = analysis.run()
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(res["summary"], indent=1, default=str))
    paths = figures.make_all(res)
    paths.append(interactive.write(res))
    print_summary(res["summary"])
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html")):
            print("  ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
