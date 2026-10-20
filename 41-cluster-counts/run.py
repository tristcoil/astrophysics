"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run
    python run.py --no-figures     # analysis only (results cached for --figures-only)

The first run downloads ~47 MB (77 MB kept after unpacking). Once the cache is
warm the analysis runs offline in about half an hour on six cores, most of it
the Markov chain; the robustness scans take a few minutes more.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import config as C  # noqa: E402

CACHE = C.CACHE_DIR / "analysis.pkl"


def _rng(v, nd=3):
    med, lo, hi = v
    return f"{med:.{nd}f} (-{med - lo:.{nd}f} / +{hi - med:.{nd}f})"


def print_summary(s: dict) -> None:
    print("\n" + "=" * 78)
    print("RESULTS: counting eRASS1's galaxy clusters")
    print("=" * 78)
    sp = s["published_point"]
    print(f"\nSample: {s['sample']['n_used']} clusters ({s['sample']['sum_weights']:.0f} after contamination weights)")
    print(f"At Ghirardini et al.'s best fit our model predicts {sp['N_pred']:.0f} clusters; "
          f"observed {sp['N_obs_weighted']:.0f}")
    f = s["fixed_scaling"]
    print("\nPublished count rate-mass relation held fixed (consistency scan):")
    print(f"  Omega_m = {_rng(f['Om'])}   sigma_8 = {_rng(f['sigma8'])}   S_8 = {_rng(f['S8'])}")
    if "mcmc" in s:
        m = s["mcmc"]
        print("Relation free within its published uncertainties (MCMC):")
        print(f"  Omega_m = {_rng(m['Om'])}   sigma_8 = {_rng(m['sigma8'])}   S_8 = {_rng(m['S8'])}")
        print(f"  (acceptance {s['mcmc_meta']['acceptance']:.2f}, {s['mcmc_meta']['seconds'] / 60:.0f} min)")
    print("Ghirardini et al. 2024: Omega_m = 0.29 (+0.01/-0.02), sigma_8 = 0.88 +- 0.02, S_8 = 0.86 +- 0.01")
    d3 = s["d3"]
    print(f"\nD3's CMB fit grown to today: sigma_8 = {d3['sigma8'][0]:.3f} +- {d3['sigma8'][1]:.3f}, "
          f"S_8 = {d3['S8'][0]:.3f} +- {d3['S8'][1]:.3f} (Planck best-fit check: {d3['planck_check_sigma8']:.4f} vs 0.8120)")
    for k in ("S8_vs_planck_fixed_scaling", "S8_vs_planck_mcmc"):
        if k in s:
            print(f"  S_8 difference from Planck 2018 ({k.split('_')[-1]}): {s[k]:+.1f} sigma")
    print("\nPlanck SZ masses / eRASS1 masses for the same clusters (S/N > 6), by mass estimator:")
    for r in s["hydrostatic"]["table"]:
        if r["snr_min"] == 6.0:
            print(f"  {r['estimator']:<48s} {r['one_minus_b']:.3f} +- {r['err']:.3f}  ({r['n']} clusters)")
    print("Counts re-read with a global mass bias (M_true / M_lensing):")
    for k, v in s["mass_scale"].items():
        print(f"  mass scale {k}: S_8 = {_rng(v['S8'])}  sigma_8 = {v['sigma8'][0]:.3f}")
    b = s["brake"]
    print(f"\nDark energy's brake: growth since z = 1100 is {b['ratio_since_cmb'][0]:.3f} of a universe without it;"
          f" sigma_8 would be {b['sigma8_noL_phys']:.2f} instead of {b['sigma8_lcdm']:.3f};"
          f" {b['boost_gt_1e+15'][0]:.1f}x as many clusters above 1e15 Msun today")
    print("\nRobustness (consistency scan, one change at a time):")
    for v in s["variants"]:
        print(f"  {v['label']:<45s} sigma_8 {v['sigma8'][0]:.3f}  S_8 {v['S8'][0]:.3f}")
    print()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    ap.add_argument("--no-figures", action="store_true", help="run the analysis, skip the figures")
    ap.add_argument("--processes", type=int, default=None, help="worker processes (default: up to 6)")
    ap.add_argument("--steps", type=int, default=C.MCMC_STEPS, help="MCMC steps per walker")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        from fetch_data import fetch_all

        import analysis

        fetch_all()
        res = analysis.run(processes=args.processes, mcmc_steps=args.steps)
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(res["summary"], indent=1, default=str))
    print_summary(res["summary"])
    if args.no_figures:
        return 0
    import figures
    import interactive

    paths = figures.make_all(res)
    paths.append(interactive.write(res))
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html")):
            print("  ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
