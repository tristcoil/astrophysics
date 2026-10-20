"""Build the whole project with one command.

    python run.py            # fetch what is missing, every analysis step (cached steps
                             # reused), every figure, the interactive page, a summary
    python run.py --force    # recompute every cached analysis step (~2-3 h of CAMB fits)
    python run.py --check    # also run the self-checks in test_physics.py (under a minute)

Once ``data/`` is filled this runs offline: nothing cached is fetched again.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))


def summary(r: dict) -> None:
    import config as cfg

    print("\n" + "=" * 78)
    print("P3 — The cosmic neutrino background")
    print("=" * 78)
    rl = r["relic"]
    print(f"Relic neutrinos today: {rl['n_nu_total_cm3']:.1f} per cm³ "
          f"({rl['n_nu_total_cm3_dec']:.1f} with decoupling heat), T = {rl['T_nu_K']:.4f} K; "
          f"CMB photons {rl['n_gamma_cm3']:.1f} per cm³")
    d = r["decoupling"]["runs"][cfg.DECOUPLING_HEADLINE]
    print(f"Our N_eff from neutrino decoupling: {d['neff']:.4f} "
          f"(Bennett et al. 2021: 3.0440 ± 0.0002)")
    for o in ("NO", "IO"):
        m = r["masses"][o]
        print(f"{o}: minimum sum {m['min_sum_eV'] * 1e3:.1f} meV; KATRIN's m_β limit allows "
              f"Σm < {m['katrin_sum_limit_eV']:.2f} eV")
    fits = r.get("fits", {})
    for key, s in fits.get("neff", {}).items():
        iv = s["interval68"]
        cmp_ = s.get("compare") or {}
        ref = cmp_.get("value")
        ref_s = f"  vs {cmp_['label']}: {ref[0]:.2f} ± {ref[1]:.2f}" if ref else ""
        print(f"N_eff {s['label']:46s} {iv['best']:.2f} -{iv['minus']:.2f} +{iv['plus']:.2f}"
              + ref_s)
    for key, s in fits.get("mnu", {}).items():
        lim = s.get("limit95_flat_ge_0")
        cmp_ = s.get("compare") or {}
        ref_s = f"  vs {cmp_['label']}: < {cmp_['value']}" if cmp_.get("value") else ""
        print(f"Σm_ν  {s['label']:46s} < {lim:.3f} eV (95%)" + ref_s if lim else
              f"Σm_ν  {s['label']:46s} (limit not bounded by the grid)")


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
    explorer.build()
    summary(results)
    if args.check:
        import test_physics

        return test_physics.main([])
    return 0


if __name__ == "__main__":
    sys.exit(main())
