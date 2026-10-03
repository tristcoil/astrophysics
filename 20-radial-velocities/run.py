"""Run the whole project end to end.

    python run.py                 # everything
    python run.py --only 51peg    # one star (51peg, 55cnc, proxima, gj581)
    python run.py --no-figures

Fetches whatever is not cached, runs the analyses, writes every figure and
data/results.json, and prints a summary. MCMC chains are cached in
data/cache/, so only the first run is slow (about an hour, most of it the
33-parameter 55 Cnc model); later runs take a few minutes, all offline.
"""

from __future__ import annotations

import argparse
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

import analysis
from config import FIGURE_DIR, RANDOM_SEED, RESULTS_PATH
from fetch_data import ensure_data

STARS = [
    ("51peg", analysis.analyse_51peg),
    ("55cnc", analysis.analyse_55cnc),
    ("proxima", analysis.analyse_proxima),
    ("gj581", analysis.analyse_gj581),
]


def _rule(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


def _pm(s, fmt=".2f") -> str:
    return f"{s['median']:{fmt}} ± {0.5 * (s['plus'] + s['minus']):{fmt}}"


def print_summary(results: dict) -> None:
    print("=" * 78)
    print("RADIAL VELOCITIES — weighing planets by their stars' wobble")
    print("=" * 78)
    if "51peg" in results:
        r = results["51peg"]
        c = r["fits"]["circular"]["planet"]
        e = r["fits"]["eccentric"]["planet"]["e"]
        ds = r["discovery_season"]
        _rule(f"51 Pegasi b — {r['n']} velocities, 4 spectrographs, {r['baseline_days'] / 365.25:.0f} years")
        print(f"  discovery season ({ds['n']} ELODIE spectra, {ds['first']} to {ds['last']}): "
              f"P = {ds['P']:.4f} ± {ds['P_err']:.4f} d, K = {ds['K']:.0f} ± {ds['K_err']:.0f} m/s")
        print(f"  all data:  P = {c['P']['median']:.7f} ± {c['P']['plus']:.7f} d   "
              f"K = {_pm(c['K'])} m/s   e < {e['upper_95']:.3f} (95%)")
        print(f"  m sin i = {_pm(c['msini_mjup'], '.3f')} M_Jup = {_pm(c['msini_mearth'], '.0f')} M_E,  "
              f"a = {c['a_au']['median']:.4f} AU")
        for ref, tm in r["true_mass"].items():
            print(f"  true mass with K_p from {ref.split(' (')[0]}: {_pm(tm['mass_mjup'], '.3f')} M_Jup, "
                  f"i > {tm['incl_deg_95_lower']:.0f} deg (95%)")
        print(f"  alias of the 4.2308-d peak at {r['aliases_predicted']['sidereal_day'][0]:.4f} d "
              f"(1 sidereal day - f_b); nothing else significant (next peak "
              f"{r['next_signal']['period']:.2f} d, FAP {r['next_signal']['fap']:.0%})")
        b = r["elodie_berv_term"]
        print(f"  ELODIE archive systematic: {b['slope_ms_per_kms']:.2f} ± {b['slope_err']:.2f} m/s per km/s of "
              f"barycentric correction (jitter {b['elodie_jitter_without']:.1f} -> {b['elodie_jitter_with']:.1f} m/s)")
    if "55cnc" in results:
        r = results["55cnc"]
        _rule(f"55 Cancri — {r['n']} Lick + Keck velocities, {r['baseline_days'] / 365.25:.0f} years")
        order = " -> ".join(f"{s['identified_as'] or '?'} ({s['peak_period']:.4g} d)" for s in r["search"][:5])
        print(f"  found in order: {order}")
        es = r["e_step"]
        print(f"  planet e: {es['best_period_search_below_1d']:.4f} d; a search stopping at 1 d finds its alias "
              f"{es['best_period_if_search_stops_at_1d']:.3f} d (Butler et al. 2017: 2.817 d)")
        at = r["alias_test"]
        for k, v in at["hypotheses"].items():
            print(f"    if {k} d were the truth, a peak ratio this large occurs in {v['fraction_as_extreme']:.1%} "
                  f"of simulations (same ordering in {v['fraction_same_order']:.0%})")
        for p in r["fiducial"]["planets"][:5]:
            ecc = f"  e = {_pm(p['e'], '.3f')}" if "e" in p else ""
            print(f"  {p['label']:9s} P = {p['P']['median']:.5f} d  K = {_pm(p['K'])} m/s  "
                  f"m sin i = {_pm(p['msini_mearth'], '.1f')} M_E{ecc}")
        t = r["tess"]
        print(f"  TESS: Rp = {_pm(t['radius_earth'])} R_E, true mass {_pm(t['true_mass_earth'])} M_E, "
              f"density {_pm(t['density_gcc'], '.1f')} g/cm3; RV predicts the transit to "
              f"{t['rv_predicted_transit_offset_minutes']:+.0f} ± {t['rv_predicted_transit_sigma_minutes']:.0f} min")
        print(f"  HARPS 2012 (never fitted): scatter {r['harps_check']['rms_with_e']:.2f} m/s with e, "
              f"{r['harps_check']['rms_without_e']:.2f} without")
    if "proxima" in results:
        r = results["proxima"]
        p = r["fiducial"]["planet"]
        _rule(f"Proxima Centauri b — {r['n']} HARPS spectra")
        print(f"  strongest peak {r['search'][0]['peak_period']:.4f} d (FAP {r['search'][0]['fap']:.0e}); "
              f"K = {_pm(p['K'])} m/s, m sin i = {_pm(p['msini_mearth'])} M_E, a = {p['a_au']['median']:.4f} AU")
        ks = [v["K_b"] for v in r["activity_models"].values()]
        print(f"  K under five treatments of activity: {min(ks):.2f}-{max(ks):.2f} m/s")
        worst = max(v["power_at_planets"][next(iter(v["power_at_planets"]))] for v in r["indicators"].values())
        print(f"  activity indicators at 11.19 d: power <= {worst:.3f} in all {len(r['indicators'])}")
    if "gj581" in results:
        r = results["gj581"]
        _rule(f"Gliese 581 — {r['n']} HARPS spectra")
        for p in r["fiducial"]["planets"][:3]:
            print(f"  {p['label']:9s} P = {p['P']['median']:.4f} d  K = {_pm(p['K'])} m/s  "
                  f"m sin i = {_pm(p['msini_mearth'])} M_E")
        d = r["signal_d"]
        h = r["halpha_rotation"]
        print(f"  'planet d' found blind at {d['P']:.2f} d (K {d['K']:.2f} m/s); H-alpha rotation {h['period']:.0f} d; "
              f"RV-H-alpha r = {r['rv_halpha_correlation']['all']:+.2f}")
        m = r["models"]
        print(f"  power left at 66.8 d: {m['b, c, e']['power_near_66.8d']:.2f} -> "
              f"{m['b, c, e + H-alpha']['power_near_66.8d']:.2f} with an H-alpha term; "
              f"epochs agree p = {r['epochs']['p_consistent']:.2f}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", choices=[k for k, _ in STARS], help="run one star only")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    ensure_data(quiet=True)
    results = {"stellar_masses": analysis.check_stellar_masses()}
    artefacts = {}
    for i, (key, fn) in enumerate(STARS):
        if args.only and key != args.only:
            continue
        t0 = time.time()
        rng = np.random.default_rng([RANDOM_SEED, i])
        results[key], artefacts[key] = fn(rng)
        print(f"  ({time.time() - t0:.0f} s)")
    if args.only and RESULTS_PATH.exists():
        # keep the other stars' numbers from the last full run
        import json
        old = json.loads(RESULTS_PATH.read_text())
        results = {**old, **analysis._jsonable(results)}
    analysis.save_results(results)

    if not args.no_figures:
        import figures
        written = []
        makers = {"51peg": figures.make_51peg, "55cnc": figures.make_55cnc,
                  "proxima": figures.fig_proxima, "gj581": figures.fig_gj581}
        for key, art in artefacts.items():
            written += makers[key](results[key], art)
        if all(k in results for k, _ in STARS):
            written += figures.fig_mass_period(results)
            written += figures.fig_comparison(results)
        if all(k in artefacts for k, _ in STARS):
            import explorer
            written += explorer.build(results, artefacts)
        print(f"\nwrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    print(f"\nNumbers saved to {RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
