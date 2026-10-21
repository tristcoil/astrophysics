"""Run the whole project end to end.

    python run.py

Downloads whatever is not already cached, builds the catalogue, makes every
measurement, writes the nine figures and the interactive 3D map, and prints a
summary. Once the cache is warm this runs offline.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
import map3d
from config import FIGURE_DIR, MAP_HTML_PATH, NAMED_RADIUS_PC
from fetch_data import fetch_all


def print_summary(res: dict) -> None:
    b, s = res["build"], res["sample"]
    print("\n" + "=" * 78)
    print("THE SOLAR NEIGHBOURHOOD: every star within 100 pc")
    print("=" * 78)

    g, h = b["gaia"], b["hipparcos"]
    print(f"\nCatalogue: {s['n_total']:,} stars = {s['n_gaia']:,} Gaia DR3 + {s['n_hipparcos']:,} Hipparcos")
    print(f"  Gaia candidates (parallax > 10 mas, 10% precision): {g['n_candidates']:,}")
    print(f"  rejected as spurious (fidelity < 0.5): {g['n_rejected_fidelity']:,}; "
          f"{100 * g['inner_galaxy_share_rejected']:.0f}% of them toward the crowded inner Galaxy, "
          f"a patch holding {100 * g['inner_galaxy_sky_share']:.1f}% of the sky")
    print(f"  kept sample in that patch: {100 * g['inner_galaxy_share_kept']:.1f}% (isotropic would be "
          f"{100 * g['inner_galaxy_sky_share']:.1f}%)")
    print(f"  a RUWE < 1.4 cut would also have dropped {100 * g['kept_fraction_ruwe_above_1p4']:.0f}% of the real stars")
    print(f"  our Galactic coordinates vs Gaia's own: max difference {g['lb_max_disagreement_arcsec']:.3f} arcsec")
    print(f"  colours contaminated by neighbours, reclassified by brightness: {g['n_colour_contaminated']:,}")
    print(f"  Hipparcos additions: {h['added_from_hipparcos']:,} ({h['added_brighter_than_v4']} brighter than V = 4; "
          f"{h['added_in_hipparcos_multiple_systems']} in multiple systems)")
    print(f"  Hipparcos-based G vs Gaia's measured G: bias {h['g_estimate_bias_mag']:+.3f}, "
          f"scatter {h['g_estimate_scatter_mag']:.3f} mag ({h['g_estimate_n_stars']:,} stars)")

    c = res["census"]
    r20 = c["by_radius"][str(int(NAMED_RADIUS_PC))]
    eye = c["naked_eye"][str(int(NAMED_RADIUS_PC))]
    print(f"\nCensus within {NAMED_RADIUS_PC:g} pc ({r20['total']:,} objects):")
    print(f"    {'class':<22}{'all':>8}{'naked eye':>12}")
    for cls in analysis.CLASS_ORDER:
        if r20["counts"][cls] or eye["counts"][cls]:
            print(f"    {analysis.CLASS_LABELS[cls]:<22}{r20['counts'][cls]:>8,}{eye['counts'][cls]:>12,}")
    print(f"  red dwarfs are {100 * c['m_dwarf_share_of_stars_20pc']:.1f}% of the stars; none is visible to the naked eye")
    sky = res["naked_eye_sky"]
    print(f"  of the {sky['n_visible']:,} stars visible to the naked eye anywhere in the sky, "
          f"{100 * sky['fraction_inside_map']:.0f}% lie within this map (median distance {sky['median_distance_pc']:.0f} pc)")

    st = res["structure"]
    print(f"\nStructure (neighbours within {st['search_radius_pc']:g} pc; best clump by chance: "
          f"{max(st['null_max_neighbours'])} neighbours)")
    for p in st["peaks"]:
        m = p.get("match") or {}
        verdict = "REAL" if p["beats_null"] else "not significant"
        print(f"  [{verdict}] {m.get('name', 'unidentified')}: {p['peak_neighbours']} neighbours; centre at "
              f"{p['distance_pc']:.1f} pc (l = {p['l_deg']:.1f}, b = {p['b_deg']:.1f})")
        if m:
            print(f"      literature {m['literature_distance_pc']:.1f} pc ({m['reference']}): "
                  f"{m['distance_difference_percent']:+.1f}%; ~{p['excess_within_mask_radius']:.0f} stars above the "
                  f"field within 10 pc")

    ns = res["north_south"]
    print(f"\nAbove or below the disc: {100 * (ns['pooled_ratio_10_90pc'] - 1):.1f} +- {100 * ns['pooled_ratio_err']:.1f}% "
          f"more stars below the Sun than above ({ns['pooled_significance_sigma']:.1f} sigma), |z| = 10-90 pc")

    v = res["validation"]
    gc = v["gcns"]
    print("\nValidation")
    print(f"  vs GCNS: {100 * gc['ours_in_gcns_fraction']:.2f}% of our Gaia stars are in it; we recover "
          f"{100 * gc['gcns_within_100pc_found_fraction']:.1f}% of its {gc['gcns_within_100pc']:,} objects within 100 pc")
    print(f"           white dwarfs: {gc['wd_both']:,} agree, {gc['wd_ours_only']} ours only, "
          f"{gc['wd_gcns_only']} GCNS only")
    st_t = v["simbad_types"]
    print(f"  vs SIMBAD spectral types (within {NAMED_RADIUS_PC:g} pc): dwarf letter agrees for "
          f"{100 * st_t['dwarf_letter_agreement']:.1f}% of {st_t['n_dwarfs_compared']:,} "
          f"({100 * st_t['dwarf_letter_within_one_class']:.1f}% within one class)")
    t10 = v["census_10pc"]["by_kind"]
    for kind in ("star", "white dwarf", "brown dwarf"):
        k = t10[kind]
        reasons = ", ".join(f"{n} {r}" for r, n in k["missed_reasons"].items() if n)
        print(f"  vs 10 pc census, {kind}s: {k['gaia'] + k['hipparcos']} of {k['total']} found"
              + (f"; missed: {reasons}" if reasons else ""))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download every dataset")
    parser.add_argument("--no-figures", action="store_true", help="skip the static figures")
    parser.add_argument("--no-map", action="store_true", help="skip the interactive 3D map")
    args = parser.parse_args(argv)

    fetch_all(refresh=args.refresh)
    print("Analysing")
    results, artefacts = analysis.analyse()
    analysis.save_results(results, artefacts["catalogue"])

    if not args.no_figures:
        written = figures.make_all(results, artefacts)
        print(f"  wrote {len(written)} figure files to {FIGURE_DIR}")
    if not args.no_map:
        path = map3d.make_map(artefacts["catalogue"], results)
        print(f"  wrote {path} ({MAP_HTML_PATH.stat().st_size / 1e6:.1f} MB)")

    print_summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
