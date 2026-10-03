"""Run the whole project end to end.

    python run.py            # everything (a few minutes, most of it the Macquart fits)
    python run.py --quick    # baseline Macquart fit only, no systematics variants
    python run.py --figures-only   # redraw figures from the last run, seconds

Downloads whatever is not already cached, runs every analysis, writes
``data/results.json`` and the figures, and prints the headline numbers. Once the
cache is warm this runs offline.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
from config import FIGURE_DIR, RESULTS_PATH
from fetch_data import fetch_all


def _q(d: dict, spec=".3f") -> str:
    """'median (+up / -down)' from a quantile dict."""
    return (f"{format(d['q500'], spec)} (+{format(d['q840'] - d['q500'], spec)}"
            f" / -{format(d['q500'] - d['q160'], spec)})")


def print_summary(res: dict) -> None:
    line = "=" * 78
    print(f"\n{line}\nRESULTS\n{line}")

    print("\nCHIME/FRB Catalog 2 selection")
    for label, n in res["chime_cut_flow"]:
        print(f"  {label:<44}{n:>7,}")

    chk = res["ne2001_map_check"]
    print("\nNE2001 map frame check (CHIME's NE2001 / map value, median):")
    print(f"  Galactic   {chk['galactic']['median_ratio_catalogue_over_map']:.3f}   "
          f"equatorial {chk['equatorial']['median_ratio_catalogue_over_map']:.3f}   -> map read as Galactic")

    if "dispersion" in res:
        print("\n1. Dispersion law, measured on single bursts")
        for tns, d in res["dispersion"].items():
            print(f"  {tns}: index alpha = {d['alpha']:.5f} +/- {d['alpha_err']:.5f} (cold plasma: -2), "
                  f"{d['n_subbands']} sub-bands, sweep {d['sweep_s']:.2f} s")
            print(f"    DM (alpha fixed) {d['dm_alpha_fixed']:.3f} +/- {d['dm_alpha_fixed_err']:.3f} "
                  f"vs CHIME {d['dm_chime']:.3f}; residual rms {d['residual_rms_ms']:.3f} ms")

    s = res["sky"]
    print("\n2. Sky distribution (|b| < 10 deg: observed vs isotropic, same exposure)")
    for key, label in (("frb", "FRB sources"), ("pulsar_chime_sky", "pulsars (Dec > -11)")):
        v = s[key]
        print(f"  {label:<22} {v['frac_abs_b_lt_10_observed']:6.1%} vs {v['frac_abs_b_lt_10_isotropic']:6.1%}"
              f"   KS p = {v['ks_vs_isotropic_pvalue']:.1e}")
    e = s["exposure"]
    print(f"  exposure: median {e['median_hours_dec_above_70']:.0f} h above Dec 70 vs "
          f"{e['median_hours_dec_below_70']:.0f} h below; {e['fraction_sources_dec_above_70']:.0%} of sources "
          f"in {e['fraction_sky_dec_above_70']:.0%} of the visible sky")

    d = res["dm"]["frb"]
    print("\n3. Dispersion measure vs the Milky Way")
    print(f"  {d['n']:,} FRB sources: median DM {d['median_dm']:.0f}, median Milky Way disk (NE2001) "
          f"{d['median_dm_mw_ne2001']:.0f}, median ratio {d['median_ratio']:.1f}x")
    bl = d["below_mw_plus_max_halo"]
    print(f"  above the disk: {d['frac_ratio_gt_1']:.1%};  above disk + 80 halo: "
          f"{d['frac_above_mw_plus_max_halo']:.1%}")
    if bl["n"]:
        rng = bl["high_latitude_dm_range"]
        print(f"  the {bl['n']} below: {bl['n_within_10deg_of_plane']} within 10 deg of the plane; "
              f"{bl['n_high_latitude']} at high latitude with low DM"
              + (f" ({rng[0]:.0f}-{rng[1]:.0f})" if rng else ""))
    print(f"  median extragalactic DM {d['median_dm_eg']:.0f} -> z <~ {d['median_z_max']:.2f} "
          f"(90th percentile z <~ {d['p90_z_max']:.2f}); highest DM {d['max_dm']:.0f} ({d['max_dm_source']})")
    for grp, v in res["dm"]["pulsar"].items():
        print(f"  pulsars, {grp:<18} n={v['n']:<5} above NE2001: {v['frac_above_ne2001']:6.1%}   "
              f"above NE2001+80: {v['frac_above_ne2001_plus_halo']:6.1%}")

    r = res["repeaters"]
    print("\n4. Repeaters vs one-offs "
          f"({r['n_repeater_sources']} repeaters' first bursts vs {r['n_one_off']:,} one-offs)")
    for q in ("width_ms", "bandwidth_mhz", "peak_freq", "dm_exc_ne2001"):
        v = r["quantities"][q]
        t = v["test_one_off_vs_repeater_first"]
        print(f"  {v['label']:<40} {v['repeater_first']['median']:8.2f} vs {v['one_off']['median']:8.2f}"
              f"   KS p = {t['ks_pvalue']:.1e}")
    bf = r["band_filling_fraction"]
    mc = r["multi_component_fraction"]
    print(f"  fill the whole band: {bf['repeater_first']:.0%} vs {bf['one_off']:.0%};  "
          f"multi-component: {mc['repeater_first']:.0%} vs {mc['one_off']:.0%}")

    m = res["macquart"]
    print("\n5. Macquart relation")
    for label, n in m["sample_cut_flow"]:
        print(f"  {label:<48}{n:>5}")
    print(f"  {'variant':<52}{'f_d':>22}{'Omega_b h70 (f_d=0.844)':>28}")
    for k, v in m["fits"].items():
        print(f"  {v['label'][:50]:<52}{_q(v['f_d']):>22}{_q(v['omega_b_h70_if_fd_fiducial'], '.4f'):>28}")
        for w in v["prior_edge_warnings"]:
            print(f"      note: {w}")
    b = m["fits"]["baseline"]
    print(f"  baseline host DM median {b['host_median_dm']['q500']:.0f} pc cm^-3, sigma_host "
          f"{b['sigma_host']['q500']:.2f}, log10 F {b['log10_F']['q500']:+.2f}")
    print(f"  Planck 2018: Omega_b h70 = {m['omega_b_h70_planck']:.4f}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--quick", action="store_true", help="baseline Macquart fit only")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    parser.add_argument("--no-fetch", action="store_true", help="never touch the network")
    parser.add_argument("--figures-only", action="store_true",
                        help="redraw figures from the last run's cached results")
    args = parser.parse_args(argv)

    if args.figures_only:
        results, art = analysis.load_cached()
    else:
        if not args.no_fetch:
            fetch_all(quiet=True)
        print("Running analyses")
        results, art = analysis.run_all(quick=args.quick)
        analysis.save_results(results)
        analysis.save_artefacts(art)

    if not args.no_figures:
        print("Drawing figures")
        written = figures.make_all(results, art)
        print(f"  wrote {len(written)} files to {FIGURE_DIR}")

    print_summary(results)
    print(f"\nFull results in {RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
