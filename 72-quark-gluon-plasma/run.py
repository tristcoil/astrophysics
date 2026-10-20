"""Build the whole project with one command.

    python run.py            # fetch what is missing, analysis (cached steps reused),
                             # every figure, the interactive page, and a summary
    python run.py --force    # recompute every analysis step (~4-5 min)
    python run.py --check    # also run the self-checks in test_physics.py (~15 s)

Once ``data/`` is filled this runs offline: nothing cached is fetched again.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def summary(r: dict) -> None:
    f = r["results"]["freeze_out"]
    a = r["results"]["andronic_2018"]
    t = r["results"]["fit_with_their_data"]
    nuc = r["results"]["nuclei_only"]
    u = r["universe"]
    ph = r["photons"]
    print("\n" + "=" * 78)
    print("P6 - Recreating the first microseconds: the quark-gluon plasma")
    print("=" * 78)
    print(f"ALICE Pb-Pb 2.76 TeV, 0-10%: T = {f['T_MeV']:.1f} +- {f['T_err']:.1f} (fit) "
          f"+- {f['T_model_sys']:.1f} (model) MeV, mu_B = {f['muB_MeV']:.1f} +- {f['muB_err']:.1f} MeV")
    print(f"   V = {f['V_fm3']:.0f} +- {f['V_err']:.0f} fm^3 (point-like); "
          f"{f['V_excluded_volume_fm3']:.0f} fm^3 with r0 = 0.3 fm; chi2 = {f['chi2']:.1f}/{f['ndf']}")
    print(f"   Andronic et al. 2018: T = {a['T']} +- {a['T_err']}, mu_B = {a['muB']} +- {a['muB_err']}, "
          f"V = {a['V']:.0f} +- {a['V_err']:.0f}")
    print(f"   our code on their data: T = {t['T']:.1f} +- {t['T_err']:.1f}, mu_B = {t['muB']:.1f}, "
          f"V = {t['V_EV']:.0f}")
    print(f"   nuclei alone: T = {nuc['T']:.1f} +- {nuc['T_err']:.1f} MeV (Andronic: 159 +- 5)")
    print(f"Lattice QCD crossover (HotQCD 2019): 156.5 +- 1.5 MeV")
    ag = u["ages"]["T_CF (our fit)"]
    age = ag["age_s_P2"] or ag["age_s_own"]
    print(f"The universe was at our T_CF {age * 1e6:.1f} +- {ag['age_err_s'] * 1e6:.1f} us after the "
          f"Big Bang, with mu_B/T = {u['mu_at_Tpc']['zero']['muB_over_T']:.1e}")
    print(f"   antibaryons annihilated (mu_B = T) at {u['T_antibaryons_gone_MeV']:.0f} MeV, "
          f"{u['age_antibaryons_gone_s'] * 1e3:.2f} ms")
    print(f"Direct photons: T_eff = {ph['T_eff_MeV']:.0f} +- {ph['stat_MeV']:.0f} (stat) +- "
          f"{ph['syst_MeV']:.0f} (syst) MeV (ALICE: 304 +- 11 +- 40 without pQCD subtraction)")
    for e, rec in r["star"]["ours"].items():
        q = rec["eq"]
        print(f"STAR {float(e):5.1f} GeV: T = {q['T']:.1f} +- {q['T_err']:.1f}, "
              f"mu_B = {q['muB']:.0f} +- {q['muB_err']:.0f} MeV")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="recompute every cached analysis step")
    ap.add_argument("--check", action="store_true", help="also run test_physics.py")
    args = ap.parse_args(argv)
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    import analysis
    import explorer
    import fetch_data
    import figures

    fetch_data.fetch_all()
    results = analysis.run(force=args.force)
    print("Figures")
    figures.main([])
    print("Interactive page")
    explorer.main()
    summary(results)
    if args.check:
        import test_physics

        return test_physics.main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
