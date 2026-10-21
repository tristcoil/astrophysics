"""Run the whole project end to end.

    python run.py

Fetches whatever is not already cached, runs the five analyses, writes every
figure, and prints a summary. Once the cache is warm this runs offline.
"""

from __future__ import annotations

import argparse
import warnings

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)
# Converting 2029 times to UTC: leap seconds after today are not yet known, so
# ERFA flags the year as "dubious". The effect is at most a second or two.
warnings.filterwarnings("ignore", message=".*dubious year.*")

import analysis
import catalog
from config import FIGURE_DIR, RANDOM_SEED, RESULTS_PATH
from fetch_data import ensure_data, load_manifest


def _rule(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


def _years(interval: float) -> str:
    return f"{interval:,.0f} yr" if interval < 1e5 else f"{interval / 1e6:.2f} Myr"


def print_summary(r: dict) -> None:
    s = r["snapshot"]
    print("=" * 78)
    print(f"NEAR-EARTH ASTEROID CLOSE APPROACHES — JPL data retrieved {s['now_utc'][:10]}")
    print("=" * 78)

    _rule(f"The snapshot: {s['window'][0]} to {s['window'][1]}")
    p, f = s["past"], s["future"]
    print(f"  last {s['days_before']} days : {p['count']:>4} approaches within 0.05 au "
          f"({p['per_day']:.1f}/day), {p['within_1ld']} inside the Moon's distance")
    fo = p["first_observed"]
    print(f"                   first observed before: {fo['before']}, same day: {fo['same_day']}, "
          f"only after: {fo['after']}; {p['first_observed_in_window']} first seen this month")
    print(f"  next {s['days_after']} days  : {f['count']:>4} predicted ({f['per_day']:.2f}/day), "
          f"{f['within_1ld']} inside the Moon's distance; {f['found_over_a_year_ahead']} found "
          f"over a year ahead")
    print(f"  past rate / future rate: {s['rate_ratio_past_to_future']:.0f}x")
    for key in ("closest", "largest"):
        o = s[key]
        print(f"  {key:<8}: {o['name']}, {o['time'][:10]}, {o['dist_ld']:.3f} LD "
              f"({o['dist_km']:,.0f} km), ~{o['diameter_m']:.0f} m")

    e = r["encounters"]
    y0, y1 = e["survey_years"]
    _rule(f"Encounters, {y0}–{y1}")
    print(f"  inside 1 LD: {e['inside_1ld_total']:,} passes ({e['inside_1ld_per_year']:.0f}/yr); "
          f"inside geostationary distance: {e['inside_geo_total']}")
    for label, c in e["per_class"].items():
        print(f"    {label:<14} {c['inside_1ld_per_year']:>6.1f}/yr inside 1 LD   "
              f"{c['within_005au_per_year']:>7.1f}/yr within 0.05 au")
    big = e["largest_inside_1ld"]
    print(f"  largest inside 1 LD: {big['name']} ({big['time'][:10]}, ~{big['diameter_m']:.0f} m, "
          f"{big['dist_ld']:.2f} LD)")
    c0 = e["closest_known_passes"][0]
    print(f"  closest known pass (no impact): {c0['name']}, {c0['time'][:10]}, "
          f"{c0['dist_km']:,.0f} km from Earth's centre = {c0['altitude_km']:,.0f} km up, ~{c0['diameter_m']:.1f} m")
    ry0, ry1 = e["recent_years"]
    print(f"  {ry0}–{ry1}: {e['recent_inside_1ld']} passes inside 1 LD, first observed only after "
          f"closest approach: {e['recent_fraction_after']:.0%}, on the day: "
          f"{e['recent_fraction_same_day']:.0%}; within ±3 days: {e['recent_within_3_days']:.0%}")

    ir = r["impact_rates"]
    _rule("From near misses to impact rates")
    print("  distance index k in N(<d) ∝ d^k, approaches within 0.05 au, 1900–2200 (d² law: k = 2)")
    for label, v in ir["distance_index"].items():
        lo, _, hi = v["k_16_50_84"]
        print(f"    {label:<14} k = {v['k']:.2f} ({lo:.2f}–{hi:.2f})  from {v['approaches']:,} "
              f"approaches by {v['asteroids']:,} asteroids")
    for name, v in ir["large"].items():
        lo, hi = v["interval_16_84"]
        print(f"  ≥ {name:<6}: one impact per {_years(v['interval_years'])} "
              f"(68%: {_years(lo)}–{_years(hi)}), rate {v['rate_per_year']:.2e}/yr "
              f"[{v['approaches']:,} approaches, {v['asteroids']:,} asteroids]")
    print(f"    published ≥ 1 km: one per {_years(ir['large']['1 km']['published']['interval_years'])} "
          f"(Stuart & Binzel 2004)")
    pub = ", ".join(f"{k} {v:.1e}" for k, v in ir["large"]["140 m"]["published"]["rates_per_year"].items())
    print(f"    published ≥ 140 m: {pub} per yr")
    a = ir["apophis_like_passes"]
    print(f"  ≥ 340 m passing within {a['distance_km']:,.0f} km: one per {_years(a['interval_years'])} "
          f"(68%: {_years(a['interval_16_84'][0])}–{_years(a['interval_16_84'][1])})")
    fb = ir["fireballs"]
    print(f"  fireballs {fb['years'][0]}–{fb['years'][1]} ({fb['events']} events) per year ≥ E kt vs "
          f"Brown et al. 2002:")
    for k_, v in fb["per_year_at_least_kt"].items():
        print(f"    ≥ {k_:>5} kt: {v:7.3f}   Brown: {fb['brown2002_per_year_at_least_kt'][k_]:7.3f}")
    for k_, v in ir["seen_vs_fireballs"].items():
        print(f"  {k_:>7}: known passes imply ≥ {v['close_approach_lower_limit']:.3f}/yr; fireballs "
              f"{v['fireballs']:.3f}/yr -> surveys see ~{v['fraction_seen']:.0%} of them")

    k = r["risk"]
    _rule(f"The risk list: {k['objects']:,} objects, {k['possible_impacts']:,} possible impacts")
    print(f"  Torino (max per object): {k['torino_max_counts']}")
    tc = k["torino_check"]
    print(f"  our Torino for {tc['rated_dates']:,} dates within 100 yr: all zero = {tc['ours_all_zero']}, "
          f"agree with Sentry = {tc['agree_with_sentry']}; beyond 100 yr would rate "
          f"{tc['what_if_within_100yr']}")
    pc = k["palermo_check"]
    print(f"  our Palermo within 0.02 of Sentry's: {pc['within_0p02']:.1%} "
          f"(impacts > 2 yr ahead: {pc['within_0p02_more_than_2yr_ahead']:.1%}, max diff "
          f"{pc['max_diff_more_than_2yr_ahead']:.3f})")
    for o in k["above_minus2"]:
        print(f"  Palermo > -2: {o['name']} {o['ps_cum']:+.2f}, p = {o['ip']:.1e}, "
              f"{o['diameter_m']:,.0f} m, {o['years']}")
    mp = k["most_probable"][0]
    print(f"  most probable: {mp['name']} p = {mp['ip']:.1%}, {mp['diameter_m']:.0f} m, "
          f"{mp['years']}, last seen {mp['last_observed']}")
    print(f"  median size {k['median_diameter_m']:.0f} m; ≥140 m: {k['at_least_140m']}, ≥1 km: "
          f"{k['at_least_1km']}; {k['below_1mt_fraction']:.0%} below 1 Mt; "
          f"{k['not_seen_for_5_years_fraction']:.0%} unobserved for 5+ years")
    print(f"  removed since {k['removed_first']}: {k['removed_total']:,} (ever listed: "
          f"{k['ever_listed']:,}); Apophis removed {k['famous_removed']['99942']}, "
          f"2024 YR4 {k['famous_removed']['2024 YR4']}")

    ap = r["apophis"]
    c = ap["closest_approach"]
    _rule("Apophis, 13 April 2029")
    print(f"  closest approach {c['time_utc']} UTC: {c['distance_km']:,.1f} km from Earth's centre "
          f"({c['altitude_km']:,.0f} km up), 3σ ±{c['horizons_3sigma_range_km']:.1f} km, "
          f"{c['speed_kms']:.3f} km/s, declination {c['declination_deg']:+.1f}°")
    tb = ap["two_body"]
    print(f"  path bent by {tb['measured_turn_deg_pm1day']:.1f}° (two-body formula: "
          f"{tb['turning_angle_deg']:.1f}°), e = {tb['eccentricity']:.2f}")
    g = ap["geostationary"]
    print(f"  closer than the geostationary distance: {g['inside_geo_distance']}; nearest point "
          f"of the ring {g['min_distance_to_ring_km']:,.0f} km away; crosses the equator plane at "
          f"{', '.join(f'{x:,.0f}' for x in g['equator_crossing_radius_km'])} km")
    b, a_ = ap["orbit_before"], ap["orbit_after"]
    print(f"  orbit: {b['class']} a = {b['a_au']:.3f} au, P = {b['period_days']:.1f} d  ->  "
          f"{a_['class']} a = {a_['a_au']:.3f} au, P = {a_['period_days']:.1f} d")
    rec = ap["record"]
    print(f"  orbit from {rec['optical_observations']:,} optical + {rec['radar_delay']} radar-range + "
          f"{rec['radar_doppler']} radar-Doppler observations, {rec['first_obs']} to {rec['last_obs']}")
    ca = ap["approaches"]
    print(f"  JPL lists {ca['listed']} approaches within 0.5 au, {ca['first']} to {ca['last']}; "
          f"3σ width before 2029 ≤ {ca['max_3sigma_width_km_before_2029']:,.0f} km, after up to "
          f"{ca['max_3sigma_width_km_after_2029']:,.0f} km")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download everything first (JPL's tables change daily)")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    ensure_data(refresh=args.refresh, quiet=True)
    manifest = load_manifest()
    rng = np.random.default_rng(RANDOM_SEED)

    ap = catalog.load_approaches()
    now = catalog.snapshot_time()
    apophis_cad = catalog.load_apophis_approaches()
    apophis_q_km = float(apophis_cad.loc[apophis_cad["dist_au"].idxmin(), "dist_km"])

    snapshot, snap_art = analysis.analyse_snapshot(ap, now)
    encounters, enc_art = analysis.analyse_encounters(ap, now)
    rates, rate_art = analysis.analyse_impact_rates(ap, catalog.load_fireballs(), rng, apophis_q_km)
    now_year = catalog.decimal_year(manifest["sentry_impacts"]["retrieved_utc"][:19])
    risk, risk_art = analysis.analyse_risk(catalog.load_sentry_objects(), catalog.load_sentry_impacts(),
                                           catalog.load_sentry_removed(), now_year)
    apophis, apo_art = analysis.analyse_apophis(
        apophis_cad, catalog.load_apophis_record(),
        catalog.load_horizons_vectors("apophis_vectors_0"),
        catalog.load_horizons_vectors("apophis_vectors_1"),
        catalog.load_apophis_range_sigma(), catalog.load_apophis_elements(),
        catalog.load_horizons_vectors("moon_vectors"))

    results = {
        "data": {k: {kk: v[kk] for kk in ("url", "retrieved_utc", "rows", "bytes")}
                 for k, v in manifest.items()},
        "snapshot": snapshot, "encounters": encounters, "impact_rates": rates,
        "risk": risk, "apophis": apophis,
    }
    analysis.save_results(results)

    if not args.no_figures:
        import figures

        artefacts = {"snapshot": snap_art, "encounters": enc_art, "impact_rates": rate_art,
                     "risk": risk_art, "apophis": apo_art}
        written = figures.make_all(results, artefacts)
        print(f"wrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    print(f"\nNumbers saved to {RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
