"""Run the whole project end to end.

    python fetch_data.py     # once: ~75 MB from JPL Horizons, about four minutes
    python run.py            # every measurement, eight figures and the viewer

Once the cache is warm this runs offline, in about half a minute. It prints
every headline measurement next to the published value it is graded against.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
import viewer
from config import FIGURE_DIR, SNAPSHOT_DATE


def _rule(title: str) -> None:
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def print_summary(res: dict) -> None:
    snap = res["snapshot"]["bodies"]
    _rule(f"WHERE EVERYTHING IS — {SNAPSHOT_DATE}, 00:00 TDB")
    for key, s in snap.items():
        print(f"  {s['name']:<15} {s['dist_au']:>9.3f} au from the Sun   {s['speed_kms']:>6.2f} km/s"
              f"   {'bound' if s['bound'] else 'UNBOUND'}")

    k = res["kepler3"]
    _rule("KEPLER'S THIRD LAW — P² ∝ a³, and what the Sun weighs")
    print(f"    {'body':<9}{'orbits':>7}{'P (days)':>13}{'a (au)':>12}"
          f"{'GM excess':>12}{'model':>9}{'meas-model':>12}   (parts per million)")
    for key, r in k["bodies"].items():
        err = f"±{r['excess_err_ppm']:.1f}" if r["excess_err_ppm"] is not None else ""
        print(f"    {r['name']:<9}{r['n_orbits']:>7}{r['period_days']:>13.4f}{r['a_au']:>12.7f}"
              f"{r['excess_ppm']:>9.1f}{err:>8}{r['predicted_ppm']:>8.1f}{r['residual_vs_model_ppm']:>10.1f}")
    print(f"  log P vs log a slope: {k['loglog_slope']:.5f} ± {k['loglog_slope_err']:.5f}  (Kepler: 1.5)")
    print(f"  inner planets agree on GM☉ to within {k['gm_sun_raw_inner_spread_ppm']:.1f} ppm "
          "before any correction")
    print(f"  GM☉ = {k['gm_sun_measured_si']:.8e} m³/s²  — JPL DE440 "
          f"{k['gm_sun_vs_de440_ppm']:+.2f} ± {k['gm_sun_measured_err_ppm']:.2f} ppm")
    print(f"  M☉ = GM☉ / G = {k['sun_mass_kg']:.5e} kg  (±{k['sun_mass_rel_err'] * 1e6:.0f} ppm, "
          "limited by G)")
    ratio = k["jupiter_sun_over_jupiter_measured"]
    err = k["jupiter_mass_ratio_err"] / k["jupiter_mass_ratio_measured"] * ratio
    print(f"  Jupiter's mass from its own orbit: Sun/Jupiter = {ratio:.1f} ± {err:.1f}  "
          f"(JPL: {k['jupiter_sun_over_jupiter_jpl']:.2f})")
    h = k["halley"]
    print(f"  Halley's comet, one orbit: GM excess {h['excess_ppm']:+.0f} ppm "
          "(planetary tugs change its period by years)")

    k2 = res["kepler2"]
    _rule("KEPLER'S SECOND LAW — Halley's comet")
    print(f"  perihelion {k2['perihelion_1_utc']} UTC at {k2['q1_au']:.4f} au; next "
          f"{k2['perihelion_2_utc']} UTC ({k2['perihelion_interval_years']:.2f} yr)")
    print(f"  aphelion {k2['aphelion_utc']} UTC at {k2['Q_au']:.3f} au")
    print(f"  speed {k2['speed_min_kms']:.3f} – {k2['speed_max_kms']:.2f} km/s ({k2['speed_ratio']:.0f}×); "
          f"area swept per day constant to ±{k2['areal_velocity_max_dev_percent']:.2f}%")
    print(f"  one year at perihelion sweeps {k2['year_at_perihelion']['angle_deg']:.1f}° and "
          f"{k2['year_at_perihelion']['area_au2']:.4f} au²; at aphelion "
          f"{k2['year_at_aphelion']['angle_deg']:.2f}° and {k2['year_at_aphelion']['area_au2']:.4f} au²")
    print(f"  ripple period {k2['ripple_period_years']:.2f} yr (Jupiter: 11.86); Sun's reflex "
          f"{k2['sun_reflex_speed_ms']:.1f} m/s predicts ±{k2['expected_ripple_percent']:.2f}%")

    iso = res["interstellar"]
    _rule("INTERSTELLAR VISITORS")
    for key, o in iso["objects"].items():
        c = o["comparison"]
        print(f"  {o['name']}  (orbit {o['orbit_solution']}, {o['n_obs']} observations)")
        print(f"    perihelion {o['perihelion_utc']} UTC, q = {o['q_au']:.4f} au, e = {o['e']:.4f}, "
              f"{o['speed_at_perihelion_kms']:.1f} km/s at perihelion")
        print(f"    v∞ in {o['v_inf_in_kms']:.4f} km/s (published {c['v_inf_kms']['published']}), "
              f"out {o['v_inf_out_kms']:.4f} (change {o['v_inf_change_ms']:+.1f} m/s)")
        print(f"    radiant RA {o['radiant_ra_deg']:.3f}° Dec {o['radiant_dec_deg']:+.3f}° in "
              f"{o['radiant_constellation']}: {c['radiant_separation_deg'] * 60:.1f}′ from the "
              "published radiant")
        uvw = ", ".join(f"{x:+.2f}" for x in o["uvw_kms"])
        extra = ""
        if "uvw_diff_kms" in c:
            extra = "  (minus published: " + ", ".join(f"{x:+.2f}" for x in c["uvw_diff_kms"]) + ")"
        print(f"    (U, V, W) = ({uvw}) km/s{extra}")
        print(f"    {o['speed_vs_lsr_kms']:.1f} km/s relative to the local stars; radiant "
              f"{o['radiant_from_apex_deg']:.1f}° from the solar apex")
    a = iso["solar_apex"]
    print(f"  solar apex: RA {a['ra_deg']:.2f}° Dec {a['dec_deg']:+.2f}° ({a['constellation']})")

    sc = res["spacecraft"]
    _rule("THE ESCAPING SPACECRAFT")
    for key, c in sc["craft"].items():
        how = "launched unbound" if c["launched_unbound"] else f"unbound after {c['escape_by'].capitalize()}"
        print(f"  {c['name']:<13} {c['today_r_au']:6.1f} au, {c['today_au_per_year']:.2f} au/yr, "
              f"{c['today_light_hours']:5.2f} light-h from Earth; heading for "
              f"{c['heading_constellation']} at v∞ {c['v_inf_kms']:.2f} km/s ({how})")
    print(f"  New Horizons' heading is {sc['new_horizons_vs_atlas_radiant_deg']:.2f}° from "
          "3I/ATLAS's radiant")
    for key, h in sc["heliopause"].items():
        print(f"  {key}: {h['ours_au']:.2f} au from the Sun on {h['date']}, when its instruments saw "
              f"it leave the heliosphere (published {h['published_au']} au)")
    print("  flybys (closest approach from the planet's centre; speed relative to the Sun, "
          "±40 days):")
    for f in sc["flybys"]:
        print(f"    {f['craft']:<13}{f['planet']:<8}{f['closest_approach_utc']}  "
              f"{f['closest_approach_km']:>10,.0f} km   v∞ {f['v_inf_in_kms']:.2f}→{f['v_inf_out_kms']:.2f}"
              f"   turn {f['turn_observed_deg']:6.2f}° (formula {f['turn_predicted_deg']:6.2f}°)"
              f"   {f['helio_speed_before_kms']:5.2f}→{f['helio_speed_after_kms']:5.2f} km/s")

    ld = res["light_day"]
    _rule("VOYAGER 1: ONE LIGHT-DAY FROM EARTH")
    print(f"  crossing {ld['crossing_utc']} UTC; NASA {ld['published_utc']} UTC; "
          f"difference {ld['difference_s']:+.3f} s")
    print(f"  never comes back: next minimum {ld['next_minimum_utc']} at {ld['next_minimum_au']:.3f} au "
          f"({ld['margin_light_minutes']:.2f} light-minutes beyond)")
    print(f"  24-hour command uplink from {ld['uplink_utc']} UTC; 24-hour downlink from "
          f"{ld['downlink_utc']} UTC; one light-day from the Sun {ld['sun_light_day_utc']} UTC")

    sky = res["sky"]
    _rule(f"THE SKY TONIGHT FROM {sky['site']['name'].upper()} ({sky['site']['tz_label']})")
    tw = sky["twilight_local"]
    print(f"  sunset {tw['sunset']['evening']}, dark {tw['astronomical']['evening']}–"
          f"{tw['astronomical']['morning']}, sunrise {tw['sunset']['morning']}")
    for key, b in sky["bodies"].items():
        if b["visible"]:
            print(f"  {b['name']:<8} mag {b['magnitude']:5.1f} in {b['constellation']:<10} highest in "
                  f"the dark {b['max_alt_dark_deg']:4.0f}° at {b['max_alt_time_local']}")
        else:
            print(f"  {b['name']:<8} below the horizon while it is dark")

    _rule("SELF-CHECKS ON REAL DATA")
    x = res["elements_crosscheck"]
    print(f"  our elements vs Horizons' own, 10 orbits: worst a {x['a_rel']:.1e} (relative), "
          f"e {x['e_abs']:.1e}, inclination {x['inc_deg']:.1e}°")
    chk = sky["astropy_check"]
    print(f"  Saturn's position vs astropy's independent ephemeris: altitude "
          f"{chk['alt_diff_arcsec']:+.1f}″, azimuth {chk['az_diff_arcsec']:+.1f}″")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip the static figures")
    parser.add_argument("--no-viewer", action="store_true", help="skip the interactive viewer")
    args = parser.parse_args(argv)

    results, art = analysis.run_all()
    analysis.save_results(results)
    if not args.no_figures:
        written = figures.make_all(results, art)
        print(f"wrote {len(written)} figure files")
    if not args.no_viewer:
        path, n = viewer.build_viewer(results)
        print(f"wrote {path} ({n / 1e6:.2f} MB of embedded data)")
    print_summary(results)
    print(f"\nFigures in {FIGURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
