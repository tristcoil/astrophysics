"""Run the whole project end to end.

    python run.py                # fetch what is missing, analyse, draw everything
    python run.py --no-figures

Fetches whatever is not cached (the two whole-catalogue Gaia scans take
~15-30 min each on the first run), runs every analysis, writes
data/results.json and the figures, and prints a summary. Later runs work
offline and take a few minutes.
"""

from __future__ import annotations

import argparse
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import analysis as A
from config import RANDOM_SEED
from fetch_data import ensure_data


def _rule(title: str) -> None:
    print(f"\n{title}\n" + "-" * len(title))


def _pm(s, fmt=".0f", unit="") -> str:
    return f"{s['median']:{fmt}} (+{s['plus']:{fmt}} −{s['minus']:{fmt}}){unit}"


def print_summary(r: dict) -> None:
    print("=" * 78)
    print("L5 — RUNAWAY SURVIVORS OF TYPE Ia SUPERNOVAE")
    print("=" * 78)
    s = r["search"]
    _rule("1. The search")
    print(f"  Gaia DR3 sources with face-value v_tan > 600 km/s: {s['n_positive']:,} with parallax > +3σ, "
          f"{s['n_negative']:,} with parallax < −3σ (pure errors)")
    for key in ("no cut", "all of these"):
        rows = {x["threshold"]: x for x in s["mirror"][key]}
        print(f"  [{key}] " + "  ".join(f">{k}σ: {rows[k]['n_pos']}/{rows[k]['n_neg']}"
                                        for k in (3, 5, 10, 20)) + "   (positive/negative)")
    st = s["shift_test"]
    print(f"  parallax shift needed to fall below 600 km/s: {st['frac_within_1sigma']:.0%} within 1σ (quoted errors)")
    for key, c in st["counts"].items():
        print(f"    [{key}] " + "  ".join(f">{k}σ: {n}" for k, n in c.items()))
    print(f"  robust (> 3σ, zero point + external errors, clean): {st['n_robust']}; unbound for any RV: "
          f"{st['n_unbound_any_rv']}")
    for u in st["unbound"]:
        print(f"    {u['source_id']} {u['known_as'] or '(new)'}: G {u['G']:.1f}, BP-RP {u['bp_rp']:.2f}, "
              f"M_G {u['M_G']:.1f}, {u['z']:.1f}σ, v_gc >= {u['v_gc_min']:.0f} (escape {u['v_esc']:.0f})")
    e = s["elbadry"]
    print(f"  El-Badry et al.'s query today: {e['n_returned']} sources; identical to their Table 1: "
          f"{e['identical_to_table1']} (missing {e['missing']}, extra {e['extra']})")
    print(f"  of {e['n_observed']} with spectra, {e['n_runaways']} are runaways; top 5 by v_tan,lower all runaways: "
          f"{e['top5_lower_all_runaways']}")
    print(f"  P(v_tan > 1000 | parallax): runaways {np.round(e['p_fast1000_runaways'], 2).tolist()}")
    print(f"                              impostors max {max(e['p_fast1000_false']):.2f}, median "
          f"{np.median(e['p_fast1000_false']):.2f}")
    print(f"  union selection (their cuts without plx/err < 5): {s['union']['n']} sources "
          f"({s['union']['n_precise']} with plx/err >= 5)")
    for x in s["shen"]:
        print(f"  {x['name']:5s} DR2 ϖ {x['plx_dr2']:.2f}±{x['err_dr2']:.2f} -> DR3 {x['plx_dr3']:.2f}±{x['err_dr3']:.2f} "
              f"({x['shift_sigma']:+.1f}σ), v_tan {x['vtan_dr2']:.0f} -> {x['vtan_dr3']:.0f} km/s, RUWE {x['ruwe_dr3']:.2f}")

    _rule("2. Kinematics (fiducial: zero-point-corrected parallaxes, McMillan17 rotation)")
    for name, k in r["kinematics"].items():
        eb = k.get("eb23") or {}
        ebv = f" | EB23 {eb['vej'][0]:.0f}" if "vej" in eb else ""
        print(f"  {name:11s} d {_pm(k['d'], '.2f')} kpc  v_tot {k['v_tot']['median']:5.0f}  "
              f"v_ej {_pm(k['v_ej'])}{ebv}  (EB conventions: {k['eb_conv_v_ej']['median']:.0f})  "
              f"t {k['t_disk']['median']:+.2f} Myr  P(unbound) {k['p_unbound']:.2f}")

    _rule("3. Orbits in McMillan (2017)'s potential")
    for name, o in r["orbits"].items():
        print(f"  {name:11s} v_ej with potential {o['v_ej_potential']['median']:5.0f} vs straight line "
              f"{o['v_ej_straight']['median']:5.0f} km/s")

    _rule("4. Remnants")
    rem = r["remnants"]
    b = rem["D6-2"]["best"]
    sc = rem["D6-2"]["shen_check"]
    print(f"  D6-2 -> {b['name']}: closest {b['sep_deg']:.2f}° from centre (radius {b['radius_deg']:.2f}°) "
          f"{b['t_yr'] / 1e3:.0f} kyr ago; inside in {100 * b['frac_inside']:.0f}% of draws; "
          f"at 90 kyr {sc['sep_at_90kyr_deg']:.2f}° (Shen et al.: {sc['published_deg']}°)")
    print(f"  chance: D6-2 alone {rem['D6-2']['chance']['p_chance']:.3f}; any D6 star "
          f"{rem['look_elsewhere']['p_any_d6']:.3f}; any of 17 {rem['look_elsewhere']['p_any_all']:.3f}")
    others = [n for n, e in rem.items() if isinstance(e, dict) and e.get("matches") and n != "D6-2"]
    print(f"  other stars with a remnant on their track: {others or 'none'}")

    _rule("5. Donor masses")
    c = r["binary"]["check"]
    print(f"  0.6 on 1.0 Msun: {c['v_06_on_10_bedard']:.0f} (Bedard) / {c['v_06_on_10_eggleton']:.0f} (cold) km/s; "
          f"1.0+1.0: {c['v_10_on_10_bedard']:.0f} / {c['v_10_on_10_eggleton']:.0f}; Shen+18: {c['shen18']}")
    print(f"  two 1.05: {c['vmax_two_105']:.0f} km/s (El-Badry+23: {c['eb23_vmax']:.0f}); two 1.20: {c['vmax_two_120']:.0f}")
    for name, x in r["binary"]["stars"].items():
        f = lambda v: "impossible" if v is None else f"> {v:.2f}"
        print(f"  {name:11s} v_ej > {x['v_ej_2sigma_lower']:.0f} (2σ): donor {f(x['mdon_min_105'])} Msun "
              f"(accretor ≤ 1.05), {f(x['mdon_min_120'])} (≤ 1.20)")

    _rule("6. Survivors vs the supernova rate (fractions are upper bounds at the oldest age seen)")
    for g, c in r["population"]["groups"].items():
        fo = c.get("fraction_at_oldest")
        if fo:
            lo, hi = c["fraction_at_oldest_range"]
            print(f"  {g}: {c['n_found']} found ({', '.join(c['names'])}); tau >= {c['oldest_age_myr']:.2f} Myr "
                  f"-> at most {100 * fo:.1f}% (68%: {100 * lo:.1f}-{100 * hi:.1f}%) of the parent rate")
    for n, p in r["population"]["per_star"].items():
        if p.get("fraction_of_parent_at_age"):
            print(f"    {n:11s} age {p['age_myr']:.2f} Myr: stars like it need <= "
                  f"{100 * p['fraction_of_parent_at_age']:.1f}% of the parent rate")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    ensure_data(quiet=True)
    rng = np.random.default_rng(RANDOM_SEED)
    results, t0 = {}, time.time()

    print("search ...", flush=True)
    results["search"], search_art = A.analyse_search()
    print("kinematics ...", flush=True)
    results["kinematics"], samples = A.analyse_kinematics(np.random.default_rng([RANDOM_SEED, 1]))
    print("orbits ...", flush=True)
    results["orbits"], orbit_art = A.analyse_orbits(samples, np.random.default_rng([RANDOM_SEED, 2]))
    print("remnants ...", flush=True)
    results["remnants"], rem_art = A.analyse_remnants(samples, np.random.default_rng([RANDOM_SEED, 3]))
    print("binaries ...", flush=True)
    results["binary"] = A.analyse_binary(results["kinematics"])
    print("population ...", flush=True)
    remnant_ages = {n: e["best"]["t_yr"] / 1e6 for n, e in results["remnants"].items()
                    if isinstance(e, dict) and "best" in e}
    results["population"], pop_art = A.analyse_population(
        results["kinematics"], results["search"]["union"]["ids"], np.random.default_rng([RANDOM_SEED, 4]),
        remnant_ages=remnant_ages)
    A.save_results(results)
    print(f"analysis done in {time.time() - t0:.0f} s")

    if not args.no_figures:
        import json
        import figures as F
        import interactive
        from config import RESULTS_PATH
        res = json.loads(RESULTS_PATH.read_text())
        written = []
        written += F.fig_search(res, search_art)
        written += F.fig_parallax_errors(res)
        written += F.fig_hr(res)
        written += F.fig_velocities(res)
        written += F.fig_trajectories(res, orbit_art)
        written += F.fig_remnant(res, rem_art)
        written += F.fig_donors(res)
        written += F.fig_population(res)
        written += F.fig_sky(res, orbit_art, rem_art)
        written += interactive.build(res, orbit_art, rem_art)
        print(f"{len(written)} figure files written")
    print_summary(A._jsonable(results))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
