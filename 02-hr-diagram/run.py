"""Run the whole project end to end.

    python run.py

Loads the cached Gaia data, applies the quality cuts, measures the features of
the diagram, writes every figure, and prints a summary. Run ``fetch_data.py``
once first; after that this runs entirely offline.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
from config import FIGURE_DIR
from fetch_data import (
    load_candidates,
    load_dwarf_sequence,
    load_named_stars,
    load_sky_counts,
)


def print_summary(res: dict) -> None:
    line = "=" * 78
    print(f"\n{line}\nRESULTS — the Gaia DR3 HR diagram of the solar neighbourhood\n{line}")

    print("\nSample (each star is attributed to the first cut it fails):")
    for step in res["cut_flow"]:
        removed = f"−{step['removed']:,}" if step["removed"] else ""
        print(f"  {step['cut']:<58}{removed:>10}{step['remaining']:>11,}")

    c = res["composition"]
    print(f"\nWhat is in it ({res['n_final']:,} stars within 100 pc):")
    for key, label in (("main_sequence", "main sequence"), ("white_dwarfs", "white dwarfs"),
                       ("giants", "giants and subgiants"), ("other", "between sequences")):
        n = c[key]
        print(f"  {label:<24}{n:>9,}   {n / res['n_final']:6.2%}")
    s = res["sun"]
    print(f"  the Sun is brighter than {s['fraction_fainter']:.1%} of them "
          f"(a lower limit: Gaia misses the faintest)")

    print("\nMeasured features vs published values:")
    for key in ("binary_sequence", "red_clump", "jao_gap", "white_dwarf_branches"):
        m = res["features"].get(key, {})
        if m.get("summary"):
            print(f"  {m['summary']}")

    k = res["kinematics"]
    print("\nKinematics (tangential velocity relative to the Sun):")
    for name, cls in k["classes"].items():
        print(f"  {name:<30}{cls['n']:>9,}   {cls['fraction']:6.2%}")
    p = k.get("parenago", {})
    if p.get("summary"):
        print(f"  {p['summary']}")

    sky = res.get("sky")
    if sky:
        print("\nAll-sky star counts:")
        print(f"  {sky['n_catalogue']:,} Gaia DR3 sources, mapped from a "
              f"{sky['sample_fraction']:.1%} random subsample")
        print(f"  half of them lie within {sky['half_of_stars_within_abs_b_deg']:.1f}° "
              f"of the Galactic plane, a band covering "
              f"{sky['sky_fraction_within_that_band']:.0%} of the sky")

    print("\nFamous stars:")
    for name, rec in res["named_stars"].items():
        where = rec["status"]
        if rec.get("kept"):
            where = (f"M_G = {rec['abs_g']:5.2f}, BP−RP = {rec['bp_rp']:5.2f}, "
                     f"{rec['distance_pc']:5.1f} pc, v_tan = {rec['v_tan_km_s']:4.0f} km/s")
        print(f"  {name:<20}{where}")
    absent = ", ".join(res["absent_bright_stars"])
    print(f"  not in Gaia DR3 at all (too bright): {absent}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    candidates = load_candidates()
    dwarf = load_dwarf_sequence()
    named = load_named_stars()
    try:
        sky = load_sky_counts()
    except FileNotFoundError:
        print("no cached sky counts; the sky map will be skipped "
              "(run `python fetch_data.py` to fetch them)")
        sky = None

    res, art = analysis.analyse(candidates, dwarf, named, sky)
    analysis.save_results(res)

    if not args.no_figures:
        written = figures.make_all(res, art)
        print(f"wrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
