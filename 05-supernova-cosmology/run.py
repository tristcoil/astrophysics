"""Run the whole project end to end.

    python run.py

Downloads whatever is not already cached, runs every fit, writes the figures,
and prints a summary comparing our numbers to the published ones. Once the cache
is warm this runs offline, in about two minutes.

    python run.py --figures-only    # redraw figures from the last analysis run
"""

from __future__ import annotations

import argparse
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
from config import DATA_DIR, FIGURE_DIR
from fetch_data import fetch_all

CACHE = DATA_DIR / "analysis_cache.pkl"


def _pm(s: dict, spec: str = ".3f") -> str:
    return f"{s['mean']:{spec}} ± {s['std']:{spec}}"


def print_summary(r: dict) -> None:
    d = r["dataset"]
    print("\n" + "=" * 78)
    print("RESULTS — Type Ia supernovae and the expansion of the universe")
    print("=" * 78)
    print(f"\n{d['n_light_curves']} light curves, z = {d['z_min']:.4f} to {d['z_max']:.2f}; "
          f"{d['n_cosmology']} at z > 0.01 used for cosmology, "
          f"{d['n_calibrator_light_curves']} in Cepheid-calibrated hosts.")

    pub = r["published"]
    rows = [
        ("Flat ΛCDM: Ωm", r["flat_lcdm"]["statsys"]["Om"], pub["flat_lcdm_Om"]),
        ("Curvature free: Ωm", r["lcdm_plane"]["Om"], pub["lcdm_Om"]),
        ("Curvature free: ΩΛ", r["lcdm_plane"]["OL"], pub["lcdm_OL"]),
        ("Flat wCDM: w", r["flat_wcdm"]["w"], pub["flat_wcdm_w"]),
        ("With Cepheids: H0", r["calibrated"]["H0"], pub["shoes_flat_lcdm_H0"]),
    ]
    print(f"\n  {'quantity':<24}{'ours':>18}{'Brout+2022':>18}{'chain check':>16}")
    keys = ["flat_lcdm:Om", "lcdm:Om", "lcdm:OL", "flat_wcdm:w", "shoes_flat_lcdm:H0"]
    for (label, ours, published), key in zip(rows, keys):
        spec = ".2f" if "H0" in label else ".3f"
        v = r.get("validation", {}).get(key)
        check = f"{v['diff_in_sigma']:+.2f}σ" if v else "—"
        print(f"  {label:<24}{_pm(ours, spec):>18}"
              f"{published[0]:>11{spec}} ± {published[1]:<5{spec}}{check:>12}")

    ex = r["lcdm_plane"]["exclusion"]
    print("\n  The discovery:")
    print(f"    ΩΛ > 0 at {ex['no_lambda']['sigma']:.1f}σ (Δχ² = {ex['no_lambda']['delta_chi2']:.1f}); "
          f"Riess et al. 1998 had 3.0–4.0σ")
    print(f"    accelerating today (q0 < 0) at {ex['no_acceleration']['sigma']:.1f}σ; "
          f"q0 = {_pm(r['lcdm_plane']['q0'], '.2f')} (curvature free)")
    print("    model comparison (Δχ² vs our best flat ΛCDM, each at its own best zero point):")
    for k, v in r["model_comparison"].items():
        print(f"      {v['label']:<46} Δχ² = {v['delta_chi2']:+8.1f}")

    c = r["calibrated"]
    print("\n  The Hubble constant:")
    print(f"    H0 = {_pm(c['H0'], '.2f')} km/s/Mpc   (stat-only covariance: {_pm(c['H0_stat_only'], '.2f')})")
    k = r["cross_check_sh0es_hf"]
    print(f"    cross-check, {k['n_hubble_flow']} SH0ES Hubble-flow SNe, no cosmology: "
          f"{k['H0']:.2f} ± {k['H0_err']:.2f}")
    print(f"    supernova absolute magnitude M = {_pm(c['M'], '.3f')}")
    print(f"    Planck 2018 (CMB): {c['planck_H0'][0]} ± {c['planck_H0'][1]}  ->  "
          f"tension {c['tension_with_planck_sigma']:.1f}σ")

    print("\n  Robustness (flat ΛCDM Ωm and H0 with slices of the data removed):")
    for label, v in r.get("robustness", {}).items():
        print(f"    {label:<44} n = {v['n']:4d}   Ωm = {_pm(v['Om'])}   H0 = {_pm(v['H0'], '.2f')}")

    print("\n  What it implies (flat ΛCDM, our H0 and Ωm):")
    print(f"    age of the universe {_pm(c['age_gyr'], '.2f')} Gyr "
          f"(matter-only at the same H0: {c['age_eds_same_H0_gyr']:.1f} Gyr)")
    print(f"    acceleration began at z = {c['transition']['z']:.2f}, "
          f"{c['transition']['lookback_gyr']:.1f} Gyr ago")
    print(f"\n  Runtime {r['runtime_s']:.0f} s.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--figures-only", action="store_true",
                        help="redraw figures from the cached analysis, without refitting")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    if args.figures_only:
        if not CACHE.exists():
            print("no cached analysis yet; run `python run.py` first")
            return 1
        results, art = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        results, art = analysis.run()
        analysis.save_results(results)
        CACHE.write_bytes(pickle.dumps((results, art)))

    if not args.no_figures:
        print("Figures")
        written = figures.make_all(results, art)
        print(f"  wrote {len(written)} files to {FIGURE_DIR}")

    print_summary(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
