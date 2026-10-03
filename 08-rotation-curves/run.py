"""Run the whole project end to end.

    python run.py

Fetches whatever is not already cached, analyses the showcase galaxies and the
full SPARC sample, writes every figure, and prints a summary that grades our
numbers against the published ones. Once the cache is warm this runs offline.

    python run.py --galaxy NGC6503 --galaxy F568-3

additionally analyses (and draws) any other SPARC galaxies by name.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import sparc
from config import FIGURE_DIR, FLAGSHIP, SHOWCASE, UPSILON_BULGE, UPSILON_DISK
from fetch_data import fetch_all, missing_files


def _pct(x: float) -> str:
    return f"{100 * x:.0f}%"


def print_galaxies(results: list[analysis.GalaxyResult]) -> None:
    print(f"\nWhat the visible matter predicts, at the outermost measured point "
          f"(stellar M/L: disc {UPSILON_DISK}, bulge {UPSILON_BULGE})")
    print(f"  {'galaxy':<10}{'type':<5}{'R_last':>8}{'V_obs':>7}{'V_vis':>7}"
          f"{'pull/visible':>14}{'unseen':>8}{'mass needed/':>14}   dark matter dominates")
    print(f"  {'':<10}{'':<5}{'':>8}{'km/s':>7}{'km/s':>7}{'':>14}{'share':>8}"
          f"{'all visible':>14}")
    for r in results:
        where = ("at every radius" if r.dark_dominates_everywhere
                 else f"beyond {r.dark_dominates_from_kpc:.1f} kpc" if r.dark_dominates_from_kpc
                 else "nowhere measured")
        print(f"  {r.label:<10}{r.hubble_type:<5}{r.r_last_kpc:>5.1f} kpc{r.v_obs_last:>7.0f}"
              f"{r.v_bar_last:>7.0f}{r.mass_ratio_last:>13.1f}x{_pct(r.dark_fraction_last):>8}"
              f"{r.dynamical_over_all_visible_last:>13.1f}x   {where}")
    print("  (pull/visible: the observed V^2 over the visible matter's, at the outermost point;"
          "\n   mass needed/all visible: mass inside that radius over every star and gas cloud "
          "in the galaxy)")

    print("\nCould heavier stars do it? The maximum-disc test "
          f"({results[0].max_disk['tolerance_sigma']:.0f} sigma of slack at both ends)")
    print(f"  {'galaxy':<10}{'M/L allowed by centre':>23}{'M/L needed at edge':>20}"
          f"{'ratio':>8}   with the needed M/L, the centre would spin")
    for r in results:
        md = r.max_disk
        if not md:
            continue
        print(f"  {r.label:<10}{md['upsilon_disk_max']:>23.2f}{md['upsilon_disk_needed']:>20.2f}"
              f"{md['overshoot']:>8.1f}   {md['inner_excess']:.1f}x faster than observed")

    print("\nDark haloes that fill the gap (stars fixed at their population value)")
    print(f"  {'':<10}{'--- cored (ISO) ---':>21}{'--------- cuspy (NFW) ---------':>39}")
    print(f"  {'galaxy':<10}{'chi2/dof':>10}{'rms km/s':>11}{'chi2/dof':>10}{'rms km/s':>11}"
          f"{'c':>7}{'LCDM c':>8}{'offset':>8}")
    for r in results:
        iso, nfw = r.halo_fits["ISO"], r.halo_fits["NFW"]
        d = nfw["derived"]
        flag = "  (a parameter hit its bound)" if iso["at_bound"] or nfw["at_bound"] else ""
        print(f"  {r.label:<10}{iso['chi2_nu']:>10.2f}{iso['rms_residual_kms']:>11.1f}"
              f"{nfw['chi2_nu']:>10.2f}{nfw['rms_residual_kms']:>11.1f}"
              f"{nfw['params']['concentration']:>7.1f}{d['concentration_lcdm']:>8.1f}"
              f"{d['concentration_offset_sigma']:>+7.1f}s{flag}")
    print("  (offset: fitted NFW concentration vs the simulations' median for that halo mass, "
          "in units of their 0.11 dex scatter)")


def print_population(pop: dict, cmp: dict) -> None:
    s, rar, btfr, md = pop["sample"], pop["rar"], pop["btfr"], pop["max_disk"]
    c_rar, c_btfr, c_md = cmp["rar"], cmp["btfr"], cmp["max_disk"]
    print(f"\nThe whole sample: {s['n_galaxies']} of {s['n_catalogue']} galaxies pass the cuts "
          f"({s['n_rejected_quality']} too disturbed, {s['n_rejected_inclination']} too face-on)")

    print("\n  Radial acceleration relation")
    print(f"    {rar['n_points']:,} points measured to 10%, from {rar['n_galaxies']} of the "
          f"{s['n_galaxies']} galaxies (published: {c_rar['n_points_published']:,} points)")
    print(f"    reproduces the published points one-to-one: {c_rar['points_matched_within_rounding']:,} "
          f"within rounding, {c_rar['points_matched_within_0p015']:,} within 0.015 dex")
    g = rar["g_dagger"] / 1e-10
    print(f"    g_dagger = {g:.3f} x 1e-10 m/s^2   +/- {rar['g_dagger_err_formal'] / 1e-10:.3f} formal, "
          f"+/- {rar['g_dagger_err_bootstrap'] / 1e-10:.3f} resampling galaxies")
    print(f"               {rar['g_dagger_upsilon_plus20'] / 1e-10:.2f} to "
          f"{rar['g_dagger_upsilon_minus20'] / 1e-10:.2f} if the stars are 20% heavier/lighter")
    print(f"    published  1.20 +/- 0.02 (random) +/- 0.24 (systematic); the published points "
          f"refitted with our method give {c_rar['g_dagger_published_data_refit'] / 1e-10:.3f}")
    print(f"    scatter: rms {rar['residual_rms_dex']:.3f} dex, core width "
          f"{rar['residual_gaussian_width_dex']:.3f} dex (published 0.13 and 0.11)")
    print(f"    the observed pull exceeds the visible matter's by 10% below g_bar = "
          f"{10 ** rar['log_g_bar_where_discrepancy_10pct']:.1e} m/s^2, and is double it below "
          f"{10 ** rar['log_g_bar_where_discrepancy_2x']:.1e} m/s^2")

    print("\n  Maximum-disc test, every galaxy")
    print(f"    no single mass-to-light ratio works in {md['n_no_single_upsilon_works']} of "
          f"{md['n_galaxies']} galaxies ({md['tolerance_sigma']:.0f} sigma slack; "
          f"{pop['max_disk_1sigma']['n_no_single_upsilon_works']} with 1 sigma)")
    exc = md["exceptions"]
    plausible = md["n_explained_with_plausible_upsilon"]
    print(f"    of the {len(exc)} where one does, {len(exc) - plausible} need M/L above what stars can "
          f"have; {plausible} work with a plausible M/L -- all with short curves "
          f"(<= {max(e['n_points'] for e in exc if 0.2 <= e['upsilon_disk_needed'] <= 0.8) if plausible else 0} points)")
    print(f"    our maximum-disc M/L vs Starkman et al. (2018): Spearman {c_md['spearman']:.2f}, "
          f"median ratio {c_md['median_ratio_ours_over_published']:.2f}")

    print("\n  Baryonic Tully-Fisher relation")
    print(f"    slope {btfr['slope']:.2f} +/- {btfr['slope_err_bootstrap']:.2f} from "
          f"{btfr['n_galaxies']} galaxies (published {c_btfr['slope_published']} +/- "
          f"{c_btfr['slope_err_published']}; their table through our code: "
          f"{c_btfr['slope_published_data_refit']:.2f})")
    st = btfr["stellar_only"]
    print(f"    intrinsic scatter {btfr['sigma_int_dex']:.3f} dex perpendicular to the line "
          f"(published fit to their table: {c_btfr['sigma_int_published_data_refit']:.3f})")
    print(f"    scatter in mass about the line: {btfr['rms_vertical_dex']:.2f} dex for stars + gas, "
          f"{st['rms_vertical_dex']:.2f} dex for stars alone (slope {st['slope']:.2f})")
    print(f"    fitting with vertical instead of perpendicular scatter gives slope "
          f"{btfr['vertical_scatter_fit']['slope']:.2f}")
    print(f"    our baryonic masses vs Lelli et al. (2019): median difference "
          f"{c_btfr['log_mbar_median_diff_dex']:+.3f} dex, rms {c_btfr['log_mbar_rms_diff_dex']:.3f} dex")

    fc = pop["freeman_check"]
    print(f"\n  Units check: SPARC disc curves vs ideal exponential discs of the catalogued "
          f"luminosity: median ratio {fc['median_ratio']:.3f} over {fc['n_galaxies']} galaxies "
          f"(1.000 means the M/L = 1 normalisation is as documented)")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--galaxy", action="append", default=[],
                        help="also analyse this SPARC galaxy (repeatable), e.g. NGC6503")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    if missing_files():
        print("Some data are not cached yet; fetching")
        fetch_all()

    galaxies = sparc.load_sparc()
    unknown = [n for n in args.galaxy if n not in galaxies]
    if unknown:
        print(f"unknown galaxy name(s): {unknown}. Names are as in SPARC, e.g. NGC6503, "
              f"UGC00128, DDO154, F568-3.", file=sys.stderr)
        return 1

    print("Analysing the showcase galaxies")
    results, arts = [], {}
    for s in SHOWCASE:
        r, art = analysis.analyse_galaxy(galaxies[s.name], s)
        results.append(r)
        arts[s.name] = art

    extra, extra_arts = [], {}
    for name in args.galaxy:
        r, art = analysis.analyse_galaxy(galaxies[name])
        extra.append(r)
        extra_arts[name] = art

    print("Analysing the whole sample (includes 3,000 bootstrap refits; ~10 s)")
    population, pop_art = analysis.analyse_population(galaxies)
    comparison = analysis.compare_with_literature(population, pop_art, results)
    analysis.save_results(results + extra, population, comparison)

    if not args.no_figures:
        import figures
        import interactive

        written = figures.make_all(results, arts, population, pop_art, comparison, FLAGSHIP)
        for r in extra:
            written += figures.figure_single_galaxy(r, extra_arts[r.name])
        written.append(interactive.write_explorer(galaxies))
        print(f"Wrote {len(written)} figure files to {FIGURE_DIR}")

    print("\n" + "=" * 94)
    print("D4 -- Dark matter from galaxy rotation curves (SPARC)")
    print("=" * 94)
    print_galaxies(results + extra)
    print_population(population, comparison)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
