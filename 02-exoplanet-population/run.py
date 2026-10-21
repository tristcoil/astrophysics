"""Run the whole project end to end.

    python run.py

Fetches whatever is not already cached, builds the samples, runs the four
analyses, writes every figure, and prints a summary. Once the cache is warm
this runs offline in a few seconds.
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

import analysis
import catalog
from config import FIGURE_DIR, RANDOM_SEED, RESULTS_PATH
from fetch_data import ensure_data, fetch_catalog, load_catalog_meta


def _rule(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


def print_summary(results: dict) -> None:
    meta = results["catalogue"]
    fill = results["fillers"]
    print("=" * 78)
    print(f"EXOPLANET POPULATION ATLAS — {meta['rows']:,} planets, "
          f"retrieved {meta['retrieved_utc'][:10]}")
    print("=" * 78)

    _rule("What the archive measured, and what it filled in")
    print(f"  masses from the empirical M-R relation : {fill['mass_from_relation']:>5,} "
          f"({fill['fraction_mass_from_relation']:.0%} of all planets) — excluded")
    print(f"  radii calculated from masses           : {fill['radius_calculated']:>5,} — excluded")
    print(f"  genuine true masses / minimum masses   : {fill['mass_true']:>5,} / {fill['mass_minimum_only']:,}")
    print(f"  planets with measured mass AND radius  : {fill['both_measured']:>5,}")
    print(f"  flagged controversial                  : {fill['controversial']:>5,} — census only")

    c = results["census"]
    _rule("Census")
    for method, n in c["by_method"].items():
        print(f"  {method:<16} {n:>5,}")
    print(f"  first planets: {c['earliest']['planet']} ({c['earliest']['year']}, "
          f"{c['earliest']['method'].lower()})")
    milestones = ", ".join(f"{int(k):,} by {v}" for k, v in c["cumulative_milestones"].items())
    print(f"  running total: {milestones}")
    for s in c["spikes"]:
        print(f"  {s['year']}: {s['planets']:,} planets, {s['reference_planets']:,} of them "
              f"from {s['reference']}")
    an = results["solar_analogs"]
    look = " · ".join(f"{k} {v['planets']}" for k, v in an["analogs"].items())
    print(f"  known planets within x{an['factor']:.0f} of a Solar System planet's mass and orbit:")
    print(f"    {look}")

    v = results["valley"]
    _rule(f"Radius valley ({v['sample_size']:,} planets, radii to < 10%)")
    lo, mid, hi = v["valley_rearth_16_50_84"]
    print(f"  valley          {v['valley_rearth']:.2f} R_e   (bootstrap 68%: {lo:.2f}–{hi:.2f})")
    print(f"  peaks           {v['peak_super_earth_rearth']:.2f} R_e (super-Earths), "
          f"{v['peak_sub_neptune_rearth']:.2f} R_e (sub-Neptunes)")
    print(f"  contrast        {v['contrast']:.2f}   (valley density / mean peak density)")
    sens = ", ".join(f"{float(bw):.3f}: {val:.2f}" for bw, val in v["bandwidth_sensitivity"].items())
    print(f"  vs kernel width {sens}")
    t = v["tilt"]
    ve = v["van_eylen_2018"]
    print(f"  tilt            slope {t['slope']:+.2f}  (68%: {t['slope_16_50_84'][0]:+.2f} to "
          f"{t['slope_16_50_84'][2]:+.2f}; 95%: {t['slope_2p5_97p5'][0]:+.2f} to "
          f"{t['slope_2p5_97p5'][1]:+.2f}; negative in {t['fraction_negative']:.1%})")
    print(f"                  Van Eylen et al. 2018: {ve['slope']:+.2f} (+0.02/-0.04)")
    print(f"  at P = 10 d     {t['valley_at_pivot_rearth']:.2f} R_e (ours) vs "
          f"{ve['valley_at_pivot_rearth']:.2f} R_e (Van Eylen)")
    print("  precision tiers (same cuts, split into thirds by radius precision):")
    for tier in v["precision_tiers"]:
        refs = ", ".join(f"{r} {f:.0%}" for r, f in tier["top_radius_references"].items())
        print(f"    {tier['label']:>8}: contrast {tier['contrast']:.2f}, "
              f"{'valley' if tier['valley_found'] else 'NO valley'}; "
              f"{tier['kepler_fraction']:.0%} Kepler; radii mostly from {refs}")
    b = v["blur_simulation"]
    print(f"  SIMULATION: best third blurred to worst third's quoted errors -> contrast "
          f"{b['contrast_mean']:.2f} ± {b['contrast_std']:.2f}")

    comp = results["composition"]
    _rule(f"Composition ({comp['small_planets']} planets < 4 R_e with mass and radius)")
    print(f"  below the valley ({comp['below_valley']}): {comp['below_consistent_with_rock']:.0%} "
          f"consistent with bare rock (within {comp['excess_threshold_sigma']:.0f} sigma)")
    print(f"  above the valley ({comp['above_valley']}): {comp['above_too_big_for_rock']:.0%} too big "
          f"for rock, {comp['above_too_big_for_water']:.0%} too big even for pure water")
    env = comp["envelope_illustration"]
    print(f"  a {env['core_mass_mearth']:.0f} M_e Earth-like core: {env['bare_rearth']:.2f} R_e bare, "
          f"{env['with_1pct_h2_rearth']:.2f} R_e with 1% H2 ({env['temperature_k']} K)")

    hz = results["habitable_zone"]
    _rule("Habitable zone (Kopparapu et al. 2013/2014)")
    print(f"  optimistic zone   : {hz['in_optimistic_hz']} planets "
          f"({hz['rocky_size_in_hz']} with radius <= {hz['rocky_size_threshold_rearth']} R_e)")
    print(f"  conservative zone : {hz['in_conservative_hz']} planets "
          f"({hz['rocky_size_in_conservative_hz']} rocky-size)")
    print(f"  red-dwarf hosts   : {hz['m_dwarf_fraction_rocky_hz']:.0%} of rocky-size HZ planets, "
          f"{hz['m_dwarf_fraction_all']:.0%} of all planets")
    names = ", ".join(p["planet"] for p in hz["rocky_size_planets"])
    print(f"  rocky-size, hottest first: {names}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download the archive table first (it changes weekly)")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    if args.refresh:
        fetch_catalog(refresh=True)
    ensure_data(quiet=True)
    meta = load_catalog_meta()

    df = catalog.load_catalog()
    samples = catalog.all_samples(df)
    rng = np.random.default_rng(RANDOM_SEED)

    census, census_art = analysis.analyse_census(samples["census"], meta["retrieved_utc"])
    analogs = analysis.analyse_solar_analogs(samples["mass_separation"])
    valley, valley_art = analysis.analyse_valley(samples["valley"], samples["valley_all"], rng)
    composition, comp_art = analysis.analyse_composition(samples["mass_radius"], valley)
    hz, hz_art = analysis.analyse_habitable_zone(samples["habitable_zone"])

    results = {
        "catalogue": meta,
        "fillers": catalog.summarise_fillers(df),
        "samples": {name: s.cut_table() for name, s in samples.items()},
        "census": census,
        "solar_analogs": analogs,
        "valley": valley,
        "composition": composition,
        "habitable_zone": hz,
    }
    analysis.save_results(results)

    if not args.no_figures:
        import figures

        artefacts = {"census": census_art, "valley": valley_art,
                     "composition": comp_art, "habitable_zone": hz_art}
        written = figures.make_all(samples, results, artefacts, meta)
        print(f"wrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    print(f"\nNumbers saved to {RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
