"""Build the whole project with one command.

    python run.py            # fetch what is missing, analysis (cached steps reused),
                             # every figure, the interactive page, and a summary
    python run.py --force    # recompute every analysis step (~10-25 min)
    python run.py --check    # also run the self-checks in test_physics.py

Once ``data/`` is filled this runs offline: nothing cached is fetched again.
"""

from __future__ import annotations

import argparse
import sys
import warnings


def summary(r: dict) -> None:
    b, bb, d = r["baryons"], r["bbn"], r["domains"]["modern"]
    cr, cl, sk = r["cosmic_rays"], r["clusters"], r["sakharov"]
    p = bb["rates"]["primat"]["pdg2026_12"]
    n = bb["rates"]["nacreii"]["pdg2026_12"]
    print("\n" + "=" * 78)
    print("P5 — Where did the antimatter go?")
    print("=" * 78)
    print(f"Today: {b['planck']['n_gamma_per_cm3']:.1f} photons per cm^3, {b['planck']['n_b_per_m3']:.3f} "
          f"baryons per m^3: one baryon per {b['planck']['photons_per_baryon'] / 1e9:.2f} billion photons")
    print(f"Above 1 GeV: one extra quark per {b['quark_excess']['pairs_per_extra_quark_1GeV'] / 1e6:.0f} "
          "million quark-antiquark pairs")
    c = b["catastrophe"]
    print(f"A symmetric universe: freeze-out at {c['T_freezeout_MeV']:.0f} MeV leaves eta = "
          f"{c['eta_symmetric']:.1e} ({c['observed_over_symmetric'] / 1e9:.1f} billion times too little)")
    print(f"eta from deuterium: {p['eta_median'] * 1e10:.3f} +- {p['eta_err'] * 1e10:.3f} (PRIMAT rates), "
          f"{n['eta_median'] * 1e10:.3f} +- {n['eta_err'] * 1e10:.3f} (NACRE-II) x 1e-10")
    print(f"eta from the CMB:   {b['planck']['eta'] * 1e10:.3f} (Planck 2018), "
          f"{b['d3']['eta'] * 1e10:.3f} +- {b['d3']['eta_err'] * 1e10:.3f} (D3)")
    print(f"Antimatter domains must exceed d0 = {d['headline_d0_min_gpc']:.0f} Gpc (COMPTEL, conservative heating); "
          f"distance to the CMB {d['comoving_distance_y1100_gpc']:.1f} Gpc")
    s = cr["pbar_over_p"]["2011-2015/II"]
    print(f"Antiprotons: secondary prediction / AMS-02 above 20 GV = {s['mean_model_over_data_above_20GV']:.2f}")
    print(f"Bullet Cluster antimatter fraction < {cl['f_limit_steigman_FX']:.1e} (Fermi-LAT; Steigman 2008: 3e-6)")
    jm = sk["jarlskog_measure"]["pdg2026"]
    print(f"Standard Model CP measure at 100 GeV: {jm['T100']:.1e} (needed ~1e-10); sphalerons freeze out at "
          f"{sk['T_star_ours']:.1f} GeV; v/T at m_H = 125 GeV: {sk['v_over_T_at_125']:.2f} (needed ~1)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="recompute every cached analysis step")
    ap.add_argument("--check", action="store_true", help="also run test_physics.py")
    args = ap.parse_args(argv)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    import analysis
    import fetch_data
    import figures

    fetch_data.fetch_all()
    results = analysis.run(force=args.force)
    print("Figures")
    figures.main([])
    try:
        import explorer

        print("Interactive page")
        explorer.main()
    except ImportError:
        pass
    summary(results)
    if args.check:
        import test_physics

        return test_physics.main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
