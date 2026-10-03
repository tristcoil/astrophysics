"""Run the whole project end to end.

    python run.py                 # fetch if needed, analyse, draw, summarise
    python run.py --no-figures    # numbers only

Downloads the data if it is not already cached, finds every cycle, measures the
butterfly diagram and the Sun's rotation, writes the figures and the interactive
page, and prints every number the README quotes. Once the cache is warm this runs
offline.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
from config import FIGURE_DIR, LIT_PANEL_2019
from fetch_data import fetch_all


def _r(c: dict) -> str:
    return (f"r = {c['r']:+.2f} (95% {c['r_ci95'][0]:+.2f} to {c['r_ci95'][1]:+.2f}; "
            f"p = {c['p']:.1g}; n = {c['n']})")


def print_summary(res: dict, manifest: dict) -> None:
    rec = res["record"]
    print("\n" + "=" * 80)
    print(f"WDC-SILSO sunspot number v2.0, {rec['first_month']} to {rec['last_month']} "
          f"(retrieved {manifest['retrieved_utc'][:10]}); provisional: "
          f"{', '.join(rec['provisional_months'])}")
    print("=" * 80)
    sc = rec["smoothing_check"]
    print(f"  our 13-month smoothing vs SILSO's: max difference {sc['max_abs_diff']:.3f} over "
          f"{sc['n_months']:,} months (SILSO rounds to 0.1)")
    oc = rec["official_comparison"]
    print(f"  {rec['n_cycles']} cycles found; vs SILSO's official table: maxima {oc['maxima_exact']}/"
          f"{len(oc['rows'])} exact, minima {oc['minima_exact']}/{len(oc['rows'])} exact")
    for row in oc["rows"]:
        if row["min_offset_months"]:
            print(f"    Cycle {row['cycle']} minimum: ours {row['min_ours']} ({row['min_sn_ours']:.2f}), "
                  f"SILSO {row['min_official']} ({row['min_sn_official']})")
    L = rec["lengths"]
    a = rec["amplitudes"]
    print(f"  length {L['mean_yr']:.2f} ± {L['sd_yr']:.2f} yr ({L['n']} complete cycles); shortest "
          f"Cycle {L['shortest']['cycle']} ({L['shortest']['yr']:.2f} yr), longest Cycle "
          f"{L['longest']['cycle']} ({L['longest']['yr']:.2f} yr)")
    print(f"  mean rise {L['mean_rise_yr']:.1f} yr, mean decline {L['mean_decline_yr']:.1f} yr")
    print(f"  strongest: Cycle {a['strongest']['cycle']} ({a['strongest']['max_sn']:.1f}, "
          f"{a['strongest']['date']}); weakest: Cycle {a['weakest']['cycle']} ({a['weakest']['max_sn']:.1f}, "
          f"{a['weakest']['date']})")
    print(f"  Dalton cycles 5-7 average peak {a['dalton_cycles_mean_max']:.0f} vs "
          f"{a['other_cycles_mean_max']:.0f} for the rest")
    pg = rec["periodogram"]
    print(f"  power spectrum peak: {pg['peak_period_yr']:.2f} yr (half power "
          f"{pg['half_power_periods_yr'][0]:.2f}-{pg['half_power_periods_yr'][1]:.2f} yr)")
    sh = rec["shape_start_times_vs_hathaway"]
    print(f"  fitted cycle start times vs Hathaway (2011): {sh['mean_offset_months']:+.1f} ± "
          f"{sh['sd_offset_months']:.1f} months ({sh['n']} cycles)")
    sd = sorted(rec["spotless_days"], key=lambda s: -s["spotless_days"])
    print("  most spotless days around a minimum: " + ", ".join(
        f"{s['minimum'][:4]} ({s['spotless_days']})" for s in sd[:5]))

    ru = res["rules"]
    print("\nRULES OF THE CYCLE")
    print(f"  Waldmeier (rise time vs amplitude): {_r(ru['waldmeier'])}")
    print(f"     from Cycle 10 on: {_r(ru['waldmeier_from_cycle10'])}")
    print(f"     rise timed from 20% to 80% of the peak: {_r(ru['waldmeier_20_80'])}")
    print(f"  activity 2 yr after minimum vs amplitude: {_r(ru['early_rise'])}")
    print(f"  length of cycle n vs amplitude of n+1: {_r(ru['length_vs_next_amplitude'])}")
    print(f"  level at minimum vs amplitude: {_r(ru['minimum_vs_amplitude'])}")
    go = ru["gnevyshev_ohl"]
    last = go["pairs"][-1]
    print(f"  Gnevyshev-Ohl: holds for {go['n_hold']} of {go['n_complete']} complete pairs; broken by "
          f"{', '.join(go['violations'])}; Cycle 25 so far {last['odd_total']:.0f} vs Cycle 24's "
          f"{last['even_total']:.0f}")

    c = res["cycle25"]
    print("\nCYCLE 25")
    print(f"  minimum {c['observed_min_date']} ({c['observed_min']:.1f}); panel forecast "
          f"{LIT_PANEL_2019['min_date']} ± {LIT_PANEL_2019['min_date_pm_months']} months -> "
          f"{c['min_months_early']} months early")
    print(f"  maximum {c['observed_max']:.1f} in {c['observed_max_date']} (SILSO announced "
          f"{c['silso_announced']['value']} in {c['silso_announced']['date']}); panel forecast "
          f"{c['panel_max']:.0f} (105-125) in {c['panel_max_date']}")
    print(f"     -> {c['max_ratio_to_panel']:.2f}x the forecast, {c['max_months_early']} months early; "
          f"{c['vs_cycle24']:.2f}x Cycle 24, {c['vs_cycle23']:.2f}x Cycle 23")
    print(f"     ranks {c['rank']} of 25; mean peak of Cycles 10-24 {c['mean_max_calibration_cycles']:.0f} "
          f"(Cycle 25 is {c['vs_mean_max']:.2f}x that)")
    print("  what the early rise predicted (calibrated on Cycles 10-24):")
    for h in c["hindcast"]:
        print(f"    by {h['available']} (+{h['lead_months']:2d} mo): {h['prediction']:5.0f}  "
              f"68% {h['lo68']:.0f}-{h['hi68']:.0f}  95% {h['lo95']:.0f}-{h['hi95']:.0f}  (r = {h['r']:.2f})")
    fe = c["panel_first_excluded_68"]
    print(f"  the panel's 115 first fell outside the 68% range with data to {fe['available']}")

    hm = res["hemispheres"]
    print("\nTWO HEMISPHERES")
    for row in hm["silso"]:
        print(f"  Cycle {row['cycle']}: north peaked {row['north_peak']} ({row['north_sn']:.0f}), south "
              f"{row['south_peak']} ({row['south_sn']:.0f}): {row['south_minus_north_yr']:+.1f} yr")
    print(f"  from spot areas, Cycles 12-25: mean |lag| {hm['mean_abs_lag_yr']:.1f} yr (max "
          f"{hm['max_abs_lag_yr']:.1f}); south later in {hm['n_south_later']} of {hm['n_cycles']}")

    b = res["butterfly"]
    v, cen = b["validation"], b["census"]
    print("\nTHE BUTTERFLY DIAGRAM (RGO/USAF/NOAA sunspot groups)")
    print(f"  {cen['records']:,} daily records of {cen['groups']:,} groups, {cen['first_day']} to "
          f"{cen['last_day']}, {cen['days_observed']:,} observed days")
    print(f"  daily areas match Hathaway's derived file on {v['days_matching']:,} of {v['days_compared']:,} "
          f"days; our {v['missing_days_ours']} unobserved days = his list: {v['missing_days_agree']}")
    print(f"  spot area within 5 deg of the equator {cen['area_fraction_within_5deg']:.1%}, beyond 30 deg "
          f"{cen['area_fraction_beyond_30deg']:.1%}, beyond 40 deg {cen['area_fraction_beyond_40deg']:.2%}; "
          f"99% within {cen['lat_containing_99pct_area']:.0f} deg")
    s = b["seam"]
    print(f"  1977 seam: USAF/NOAA areas {s['factor']:.2f}x smaller (68% {s['factor_ci68'][0]:.2f}-"
          f"{s['factor_ci68'][1]:.2f}); 1977-81 {s['factor_usaf_1977_1981']:.2f}x, 1982 on "
          f"{s['factor_noaa_1982_on']:.2f}x (published 1.4, 1.43-1.49)")
    d = b["drift_law"]
    f = d["fit"]
    print(f"  drift law: {f['lambda0_deg']:.2f} deg x exp(-t / {f['tau_months']:.1f} months) "
          f"(68%: {f['lambda0_ci68'][0]:.1f}-{f['lambda0_ci68'][1]:.1f} deg, {f['tau_ci68'][0]:.0f}-"
          f"{f['tau_ci68'][1]:.0f} months); scatter {f['rms_deg']:.1f} deg; Hathaway 28 deg, 90 months")
    e = d["fit_first_108_months"]
    h = d["hathaway_t0_cycles_12_23"]
    m = d["fit_minimum_clock"]
    print(f"     first 9 years only: {e['lambda0_deg']:.1f} deg, {e['tau_months']:.1f} months; with Hathaway's own "
          f"start times: {h['lambda0_deg']:.1f} deg, {h['tau_months']:.1f} months; clock at minimum: "
          f"{m['lambda0_deg']:.1f} deg, {m['tau_months']:.1f} months")
    for k, lab in (("first_3_years", "mean latitude, first 3 yr"), ("t0", "latitude scale, t0 clock"),
                   ("minimum", "latitude scale, minimum clock"), ("hathaway_t0", "with Hathaway's t0")):
        print(f"  stronger cycles higher? {lab}: {_r(d['strength_dependence'][k]['correlation'])}")
    ov = b["overlap"]
    print(f"  near minimum: new cycle at {ov['median_lat_new_near_min']:.0f} deg, old at "
          f"{ov['median_lat_old_near_min']:.0f} deg; first new-cycle group {ov['median_first_new_months_before_min']:.0f}"
          f" months before minimum, last old one {ov['median_last_old_months_after_min']:.0f} after (medians)")

    ro = res["rotation"]
    fr = ro["fit"]
    print("\nDIFFERENTIAL ROTATION (RGO 1874-1976)")
    print(f"  {ro['n_pairs']:,} day-to-day moves of {ro['n_groups']:,} groups: Omega = {fr['A']:.3f} "
          f"± {fr['A_err_boot']:.3f} {fr['B']:+.3f} ± {fr['B_err_boot']:.3f} sin^2(lat) deg/day")
    print(f"     Balthasar et al. (1986): {ro['literature']['A']} {ro['literature']['B']:+.2f} -> differences "
          f"{ro['A_diff']:+.3f}, {ro['B_diff']:+.3f}")
    for k, p in ro["periods"].items():
        print(f"     {k:>5}: sidereal {p['sidereal_days']:.2f} d, synodic {p['synodic_days']:.2f} d")
    rb = ro["robustness"]
    print("     robustness: " + "; ".join(f"{k} A={v['A']:.3f} B={v['B']:.2f}" for k, v in rb.items()))

    lr = res["long_record"]
    print("\nBEFORE 1749")
    for k, v in lr["series"].items():
        extra = f"; 1749-1849 activity {v['early_ratio_1749_1849']:.2f}x the sunspot number's" \
            if "early_ratio_1749_1849" in v else ""
        print(f"  {v['label']}: x{v['scale_to_sn_v2']:.2f} to the sunspot number{extra}")
        if "maunder_mean_sn_equiv" in v:
            print(f"     Maunder minimum mean {v['maunder_mean_sn_equiv']:.1f} "
                  f"({v['maunder_fraction_of_modern_mean']:.1%} of the post-1749 mean)")
        if "maunder_years_with_no_spot_seen" in v:
            print(f"     no spot seen in {v['maunder_years_with_no_spot_seen']} of {v['maunder_years']} years; "
                  f"median {v['maunder_median_obs_in_blank_years']:.0f} observations in those years")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figures and the interactive page")
    args = parser.parse_args(argv)

    manifest = fetch_all(quiet=True)
    tables, results = analysis.analyse()
    analysis.save_results(results)

    if not args.no_figures:
        import figures
        import interactive

        written = figures.make_all(tables, results)
        written.append(interactive.write_page(tables, results))
        print(f"wrote {len(written)} files to {FIGURE_DIR}")

    print_summary(results, manifest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
