"""Run the whole project end to end.

    python run.py                # everything
    python run.py --no-figures

Fetches whatever is not cached (6 MB), runs every analysis, writes the
figures, data/results.json and data/blackcat_adopted.csv, and prints a
summary. The 15 population fits (nested sampling) take about half an hour
the first time and are cached in data/cache/; later runs take a minute or
two, all offline.
"""

from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)

import analysis
from config import DATA_DIR, FIGURE_DIR, RANDOM_SEED, RESULTS_PATH
from fetch_data import ensure_data


def _rule(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


def _pm(s, fmt=".2f") -> str:
    return f"{s['median']:{fmt}} +{s['plus']:{fmt}} -{s['minus']:{fmt}}"


def print_summary(r: dict) -> None:
    print("=" * 78)
    print("STELLAR-MASS BLACK HOLES IN X-RAY BINARIES — weighed by their companions")
    print("=" * 78)
    w = r["weighing"]
    _rule(f"BlackCAT online: {w['n_objects_online']} transients, {w['n_dynamical_online']} dynamically "
          f"confirmed (2016 table: {w['n_dynamical_2016']})")
    print(f"  new since 2016: {', '.join(w['new_since_2016'])}; no longer flagged as dynamical: "
          f"{', '.join(w['dropped_since_2016'])}")
    print(f"  mass function above 3.2 M_sun (no neutron star can be this heavy): "
          f"{len(w['fm_above_3_2'])} of {w['n_dynamical_online']}; above 2.9 M_sun: {len(w['fm_above_2_9'])}")
    v = w["v404_check"]
    print(f"  V404 Cyg f(M) = {_pm(v['fm_ours'], '.3f')} M_sun  (Casares & Charles 1994: "
          f"{v['published'][0]} +- {v['published'][1]})")
    c = w["cygx1_check"]
    print(f"  Cyg X-1 from its orbit: {_pm(c['m_bh'], '.1f')} M_sun  (Miller-Jones et al. 2021: "
          f"{c['published'][0]} +- {c['published'][1]})")
    for x in w["fm_inconsistent"]:
        print(f"  catalogue inconsistency: {x['name']} f(M) from its P, K2 = {x['ours']:.2f} vs listed "
              f"{x['published']} ({x['pull']:+.1f} sigma)")
    for s in w["systems"]:
        m = s.get("m1_ours") or s.get("m1_min_ours")
        tag = "" if s["in_population"] else "  (minimum)"
        print(f"    {s['name']:20s} f = {s['fm_ours']['median']:5.2f}   M = {_pm(m, '.1f'):18s}{tag:11s} "
              f"BlackCAT: {s['m1_published']}")
    g = r["mass_gap"]["fits"]
    _rule("Is there a mass gap? M_1% (lightest black hole), 10/50/90% quantiles [M_sun]")
    for label, fits in g.items():
        txt = "   ".join(f"{m[:5]} {f['m1pct_q10_q50_q90'][0]:.1f}/{f['m1pct_q10_q50_q90'][1]:.1f}/"
                         f"{f['m1pct_q10_q50_q90'][2]:.1f}" for m, f in fits.items())
        best = max(fits, key=lambda m: fits[m]["ln_evidence"])
        print(f"  {label:18s} n={fits[best]['n_systems']:2d}  {txt}   (preferred: {best})")
    print("  Farr et al. 2011 Table 8: power 4.3/6.1/6.6  expon 5.3/6.0/6.5  gauss 2.9/4.4/5.5")
    print("  Rocha et al. 2026: M_1% > 4.06 (combined, power law), > 2.56 (LMXBs, Gaussian) at 90%")
    inc = r["inclination"]
    _rule("Inclination, the dominant systematic")
    for who, h in inc["a0620"]["history"].items():
        print(f"  A0620-00 at i = {h['i']:.0f} deg ({who}): published {h['m_quoted']:g}, ours {h['m_ours']:.1f} M_sun")
    close = [d for d in inc["distance_to_gap"] if d["disc_fraction_needed"] is not None][:4]
    print("  closest to the gap: " + "; ".join(
        f"{d['common_name'].split(' (')[0] or d['name']} needs {100 * d['disc_fraction_needed']:.0f}% hidden disc light"
        for d in close))
    gv = r["graveyard"]
    _rule(f"Light vs gravitational waves: {gv['n_bbh']} binary black holes of {gv['n_gw_events']} GW events")
    x = gv["xrb_stats"]["quantiles"]
    o, wv, wc = gv["gw_stats"]["observed"]["all"], gv["gw_stats"]["w_vmax"]["all"], gv["gw_stats"]["w_chirp"]["all"]
    print(f"  X-ray binaries        10/50/90%: {x[0]:.1f} / {x[1]:.1f} / {x[2]:.1f} M_sun")
    print(f"  GW, as detected       10/50/90%: {o[0]:.1f} / {o[1]:.1f} / {o[2]:.1f} M_sun")
    print(f"  GW, 1/V_max weighted  10/50/90%: {wv[0]:.1f} / {wv[1]:.1f} / {wv[2]:.1f} M_sun")
    print(f"  GW, Mc^-5/2 weighted  10/50/90%: {wc[0]:.1f} / {wc[1]:.1f} / {wc[2]:.1f} M_sun")
    print(f"  neutron stars from pulsar timing: {gv['n_ns_pulsar_timing']} (heaviest {gv['ns_max']['name']} "
          f"{gv['ns_max']['mass']:.2f} M_sun)")
    _rule("Gaia's dormant black holes, from the DR3 orbits alone")
    for name, q in r["gaia"].items():
        chk = q["check"]
        print(f"  {name}: a0 = {_pm(q['a0_mas'])} mas, i = {_pm(q['incl_deg'], '.1f')} deg, "
              f"f_ast = {_pm(q['f_ast'])} M_sun, M2 = {_pm(q['m2'], '.1f')} M_sun")
        if "m2" in chk:
            print(f"     El-Badry et al. (astrometry only): a0 {chk['a0_mas'][0]} +- {chk['a0_mas'][1]}, "
                  f"i {chk['incl_deg'][0]} +- {chk['incl_deg'][1]}, f_ast {chk['f_ast'][0]} +- {chk['f_ast'][1]}, "
                  f"M2 {chk['m2'][0]} +- {chk['m2'][1]}; joint {chk['joint_m2'][0]} +- {chk['joint_m2'][1]}")
        else:
            m1 = q["m2_for_mstar_1"]
            print(f"     M2 for M* = 1: {_pm(m1, '.2f')} (El-Badry et al.: {chk['m2_for_mstar_1'][0]} +- "
                  f"{chk['m2_for_mstar_1'][1]}; f_ast {chk['f_ast'][0]} +- {chk['f_ast'][1]}); joint "
                  f"{chk['joint_m2'][0]} +- {chk['joint_m2'][1]}")
    ob = r["outbursts"]
    j = ob["j1820"]
    tr = j["transitions"]
    _rule("Outbursts: MAXI J1820+070 (2018) and GX 339-4")
    print(f"  hard -> soft: MJD {tr['hard_to_soft_mjd']:.1f} (window {tr['hard_to_soft_window'][0]:.1f}-"
          f"{tr['hard_to_soft_window'][1]:.1f}); Homan et al. 2020: MJD 58302-58308")
    print(f"  soft -> hard: MJD {tr['soft_to_hard_mjd']:.1f} (window {tr['soft_to_hard_window'][0]:.1f}-"
          f"{tr['soft_to_hard_window'][1]:.1f})")
    print(f"  peak: {j['peak_bat']['crab']:.1f} Crab (BAT) on MJD {j['peak_bat']['mjd']:.1f}; "
          f"L(2-50 keV) = {100 * j['peak_2_50_kev']['l_ledd']['median']:.0f}% L_Edd "
          f"(Atri et al. 2020, bolometric: 15 +- 3%)")
    print(f"  soft -> hard at L(2-50 keV) = {100 * j['soft_to_hard']['l_ledd']['median']:.1f}% L_Edd "
          f"(Atri et al. 2020, bolometric: 3 +- 1%; Maccarone 2003: 1-4%)")
    full = [o for o in ob["gx339"] if o["kind"] == "full"]
    print(f"  GX 339-4: {len(ob['gx339'])} outbursts since 2009, {len(full)} reached the soft state")
    for o in full:
        print(f"    MJD {o['start_mjd']:.0f}-{o['end_mjd']:.0f}: hard->soft at "
              f"{o.get('hard_to_soft_crab_2_20', np.nan):.3f} Crab, soft->hard at "
              f"{o.get('soft_to_hard_crab_2_20', np.nan):.3f} Crab")
    d = ob["gx339_distance_kpc"]
    print(f"  GX 339-4 distance if it switches back at MAXI J1820+070's luminosity: {_pm(d, '.1f')} kpc")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    ensure_data(quiet=True)
    rng = np.random.default_rng(RANDOM_SEED)
    results, art = {}, {}
    steps = [
        ("weighing", lambda: analysis.analyse_weighing(rng)),
        ("mass_gap", lambda: analysis.analyse_mass_gap(art["weighing"], rng)),
        ("inclination", lambda: (analysis.analyse_inclination(results["weighing"], art["weighing"], rng), None)),
        ("graveyard", lambda: analysis.analyse_graveyard(art["weighing"], rng)),
        ("gaia", lambda: analysis.analyse_gaia(rng)),
        ("outbursts", lambda: analysis.analyse_outbursts(art["weighing"], rng)),
    ]
    for key, fn in steps:
        t0 = time.time()
        results[key], art[key] = fn()
        print(f"  {key}: {time.time() - t0:.0f} s")
    analysis.save_results(results)
    adopted = pd.DataFrame([{k: v for k, v in s.items() if not isinstance(v, (dict, list))}
                            | {"fm": s["fm_ours"]["median"],
                               "m1": (s.get("m1_ours") or s.get("m1_min_ours"))["median"]}
                            for s in results["weighing"]["systems"]])
    adopted.to_csv(DATA_DIR / "blackcat_adopted.csv", index=False)

    if not args.no_figures:
        import figures
        written = []
        written += figures.fig_mass_functions(results, art)
        written += figures.fig_inclination(results, art, np.random.default_rng(RANDOM_SEED + 1))
        written += figures.fig_masses(results, art)
        written += figures.fig_population(results, art)
        written += figures.fig_inclination_bias(results, art)
        written += figures.fig_graveyard(results, art)
        written += figures.fig_light_vs_gravity(results, art)
        written += figures.fig_outburst(results, art)
        written += figures.fig_hid(results, art)
        written += figures.fig_gaia(results, art)
        import explorer
        written += explorer.build(results, art)
        print(f"\nwrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    print(f"\nNumbers saved to {RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
