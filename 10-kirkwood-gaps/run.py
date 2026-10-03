"""Run the whole project end to end.

    python run.py

Parses the cached catalogues (fetching them first if needed), measures
everything, writes the figures, and prints a summary. Once the cache is warm
this runs offline. The N-body simulation is a separate, slower step
(``python simulation.py``, ~15 minutes on 6 cores); if its output exists, the
simulation figure is drawn from it.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import astdys
import figures
import mpcorb
import simulation
from config import FIGURE_DIR, MPCORB_PATH, RESULTS_PATH


def print_summary(res: dict) -> None:
    c = res["census"]
    print("\n" + "=" * 78)
    print("E1 -- KIRKWOOD GAPS: RESULTS")
    print("=" * 78)
    print(f"\n{c['total']:,} orbits ({c['numbered']:,} numbered); "
          f"{c['quality_sample']:,} with orbit quality U <= 2")
    for name, n in c["regions"].items():
        print(f"  {name:<30} {n:>9,}")

    g = res["gaps"]
    print("\nGaps found blind in the osculating semi-major axes, then identified:")
    print(f"  {'centre (AU)':>16} {'width':>7} {'depth':>6} {'sigma':>6}   "
          f"{'resonance':<10} {'predicted':>9}  {'offset':>12}")
    for x in g["detections"]:
        ident = x["label"] if x["credible"] else "unidentified"
        kind = "  (edge)" if x["kind"] == "edge" else ""
        print(f"  {x['center']:>8.4f} ±{x['center_err']:.4f} {x['fwhm']:>7.4f} "
              f"{x['depth']:>6.2f} {x['sigma']:>6.1f}   {ident:<10} "
              f"{x['a_resonance']:>9.4f}  {1000 * x['delta_a']:>+7.1f} mAU{kind}")
    ch = g["chance_alignment"]
    if ch["hits"] == 0:
        print(f"  chance alignment of the 4 strongest interior gaps "
              f"({', '.join(g['classic_four'])}): < {3 / ch['trials']:.1e} "
              f"(0 hits in {ch['trials']:,} trials)")
    else:
        print(f"  chance alignment of the 4 strongest interior gaps: {ch['probability']:.1e}")

    print("\nNamed resonances, osculating vs proper elements (centre - predicted, mAU):")
    for r in g["named"]:
        o, pr = r["osculating"], r.get("proper", {})
        edge = " edge" if o.get("walls", {}).get("is_edge") else ""
        print(f"  {r['label']:>5}  predicted {r['a_predicted']:.4f}   "
              f"osc {1000 * (o['center'] - r['a_predicted']):+6.1f} ± {1000 * o['center_err']:.1f}"
              f"   proper {1000 * (pr.get('center', float('nan')) - r['a_predicted']):+6.1f} ± "
              f"{1000 * pr.get('center_err', float('nan')):.1f}   depth {o['depth']:.2f}{edge}")

    re_ = res["reality"]
    print("\nGap depth by asteroid size (1 = empty):")
    for row in re_["by_H"]:
        lo, hi = row["H_range"]
        print(f"  H {lo:>5.1f}-{hi:<5.1f} n={row['n']:>8,}   " + "  ".join(
            f"{k} {row[k]:.2f}" for k in ("3:1", "5:2", "7:3", "9:4")))
    e0 = re_["by_era"][0]
    print(f"  Kirkwood's era (first 91): 3:1 holds {e0['3:1_inside']} "
          f"(expected {e0['3:1_expected']:.1f}), 5:2 holds {e0['5:2_inside']} "
          f"(expected {e0['5:2_expected']:.1f})")

    h = res["resonant"]["hildas"]
    print(f"\nHildas: {h['n']:,}; {h['frac_sigma_within_90']:.1%} have |sigma| < 90 deg "
          f"(control {h['control_frac_sigma_within_90']:.1%}); triangle corners at "
          + ", ".join(f"{x:+.0f}" for x in h["theta_peaks_deg"]) + " deg from Jupiter")
    t = res["resonant"]["trojans"]
    print(f"Trojans: {t['n']:,}. L4/L5 ratio by brightness limit "
          f"(literature {t['literature_ratio']} +- {t['literature_ratio_err']}):")
    for r in t["L4_L5_by_H"]:
        lim = "all" if r["H_max"] > 50 else f"H<{r['H_max']:g}"
        print(f"  {lim:>6}: {r['L4']:>5,} / {r['L5']:>5,} = {r['ratio']:.3f} +- {r['ratio_err']:.3f}")

    if "families" in res:
        f = res["families"]
        print(f"\nFamilies: {f['fraction_in_families']:.1%} of main-belt asteroids with "
              f"proper elements belong to one. Sharpening, osculating -> proper:")
        for row in f["top"]:
            print(f"  {row['name']:<12} {row['n_members']:>6,} members  "
                  f"sin i x{row['sini_sharpening']:.1f}  e x{row['e_sharpening']:.1f}")

    s = res.get("simulation")
    if s:
        print(f"\nSIMULATION ({s['n_particles']} particles, {s['years'] / 1e6:g} Myr, "
              f"energy drift {s['energy_drift_max']:.0e}):")
        for z in s["zones"]:
            med = f"{z['median_removal_kyr']:.0f} kyr" if z["median_removal_kyr"] else "-"
            print(f"  {z['label']:>4} zone: {z['removed_frac']:.0%} removed "
                  f"(median time {med}) of {z['n']}")
        print(f"  inner edge (a < 2.3 AU): {s['inner_edge']['removed_frac']:.0%} removed")
        print(f"  elsewhere (a > 2.3, away from resonances): "
              f"{s['away_from_resonances']['removed_frac']:.0%} removed")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    parser.add_argument("--reparse", action="store_true",
                        help="re-parse the raw catalogues instead of using the .npz caches")
    args = parser.parse_args(argv)

    if not MPCORB_PATH.exists():
        import fetch_data
        fetch_data.main([])

    print("loading catalogues")
    d = mpcorb.load_mpcorb(refresh=args.reparse)
    proper = astdys.load_proper(refresh=args.reparse)
    famtab = astdys.family_table()
    sim = simulation.load_simulation()
    if sim is None:
        print("  (no simulation output yet: run `python simulation.py` for figure 08)")

    print("analysing")
    res, arr = analysis.analyse(d, proper, famtab, sim)
    analysis.save_results(res)
    print(f"  results -> {RESULTS_PATH}")

    if not args.no_figures:
        print("drawing figures")
        written = figures.make_all(res, arr, d, proper, sim)
        print(f"  wrote {len(written)} files to {FIGURE_DIR}")

    print_summary(res)
    return 0


if __name__ == "__main__":
    sys.exit(main())
