"""Run the whole project end to end.

    python run.py                 # fetch (if needed), analysis, figures, explorer, summary
    python run.py --figures-only  # redraw figures and the explorer from data/results.json
    python run.py --no-fetch      # skip the download step (everything already cached)

The first run downloads ~3 MB. After that it runs offline in about a minute.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import config as C  # noqa: E402


def print_summary(r: dict) -> None:
    ck, cmp_, nu, cal = r["checks"], r["comparisons"], r["neutrinos"], r["calendar"]
    ev = {e["key"]: e for e in cal["events"]}
    print("\n" + "=" * 78)
    print("P2 — the thermal history of the universe, from quarks to atoms")
    print("=" * 78)
    print(f"g* at 10 TeV: ideal gases {ck['g_at_10TeV']['ideal']:.2f}, this work {ck['g_at_10TeV']['lattice']:.2f} "
          f"(Saikawa & Shirai {ck['g_at_10TeV']['SS2018']:.2f} ± {ck['g_at_10TeV']['SS2018_err']:.2f})")
    print(f"inside Saikawa & Shirai's error band at {100 * cmp_['fraction_of_grid_within_SS_error']:.0f}% of "
          f"temperatures; the ideal gas misses by up to {100 * (cmp_['ideal_max_ratio_to_SS']['ratio'] - 1):.0f}% "
          f"at {1000 * cmp_['ideal_max_ratio_to_SS']['T_GeV']:.0f} MeV (the QCD crossover)")
    print(f"today: g*rho = {ck['today_full']['g_rho']:.4f}, g*s = {ck['today_full']['g_s']:.4f} "
          f"(textbook {ck['today_instantaneous']['g_rho']:.4f}, {ck['today_instantaneous']['g_s']:.4f})")
    print(f"N_eff = {nu['N_eff_ours']:.4f} (Bennett et al. 2021: 3.0440 ± 0.0002); neutrinos stop exchanging "
          f"energy at {nu['decoupling_T_MeV']['nu_mu']:.2f} (nu_mu,tau) and {nu['decoupling_T_MeV']['nu_e']:.2f} MeV (nu_e)")
    print(f"relic neutrinos today: {nu['T_nu_today_K']['instantaneous']:.4f} K, "
          f"{nu['n_nu_today_cm3']['all_instantaneous']:.1f} per cm^3 (textbook)")
    for key in ("ew_crossover", "qcd_crossover", "nu_mu_decoupling", "np_freeze_out", "nuclei_form", "equality",
                "last_scattering", "today"):
        e = ev[key]
        print(f"  {e['label'][:62]:<62} t = {e['t_s']:.4g} s")
    q = r["eos_softening"]["qcd"]
    print(f"QCD softening: w = {q['w_min_ours']:.3f} at {1000 * q['T_w_min_ours_GeV']:.0f} MeV, "
          f"{q['M_hubble_at_w_min_Msun']:.1f} solar masses inside the horizon")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--no-fetch", action="store_true")
    args = ap.parse_args(argv)
    if args.figures_only:
        results = json.loads(C.RESULTS_PATH.read_text())
    else:
        if not args.no_fetch:
            from fetch_data import fetch_all

            fetch_all()
        import analysis

        results = analysis.run_core()
    print_summary(results)
    import explorer
    import figures

    for p in figures.make_all(results):
        print("  wrote", Path(p).relative_to(C.HERE))
    print("  wrote", explorer.build(results).relative_to(C.HERE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
