"""Build the whole project with one command.

    python run.py            # fetch what is missing, analysis (cached steps reused),
                             # every figure, the interactive page, and a summary
    python run.py --force    # recompute every analysis step (~10 min)
    python run.py --check    # also run the self-checks in test_physics.py (~5 s)

Once ``data/`` is filled this runs offline: nothing cached is fetched again.
"""

from __future__ import annotations

import argparse
import sys
import warnings


def summary(r: dict) -> None:
    b, c, m, h = r["bk18"], r["constraints"], r["models"]["models"], r["horizons"]
    o = c["ours"]
    print("\n" + "=" * 78)
    print("P1 — Inflation on trial")
    print("=" * 78)
    print(f"Our BK18 likelihood: r < {b['r95']:.4f} (95%; BK18 published < 0.036), "
          f"chi2 matches Cobaya's to {r['bk18_validation']['rms_diff']:.3f} rms")
    print(f"Planck 2018 + BK18 (ours): n_s = {o['ns']['mean']:.4f} ± {o['ns']['std']:.4f}, "
          f"r < {o['r95']:.4f}  (ACT's joint run: r < "
          f"{c['chains']['act_planck-bk_r_camb']['r95']:.4f})")
    for key, lab in (("act_p-actlite-bk-l-b_r_camb", "P-ACT-LB + BK18"),
                     ("cmbspa_desi", "CMB-SPA + DESI"), ("planck18_lcdm", "Planck 2018 LCDM")):
        s = c["chains"][key]
        line = f"{lab:>18}: n_s = {s['ns']['mean']:.4f} ± {s['ns']['std']:.4f}"
        if "r95" in s:
            line += f", r < {s['r95']:.4f}"
        print(line)
    print(f"D3's own n_s (TT only): {c['d3']['value']:.4f} ± {c['d3']['error']:.4f}")
    s = m["starobinsky"]
    print(f"Starobinsky: n_s = {s['table']['N55']['ns']:.4f}, r = {s['table']['N55']['r']:.4f} "
          f"at N* = 55; instant reheating N* = {s['reheating']['instant']['n_star']:.1f}")
    print("Verdicts (best case in N* = 50-60; sigma-equivalent of the 2D region):")
    for name, v in c["verdicts"].items():
        cells = []
        for key, lab in (("planck18_bk18_ours", "Planck+BK18"), ("pactlb_bk18", "P-ACT-LB+BK18"),
                         ("cmbspa_desi_x_bk18", "CMB-SPA+DESI")):
            x = v[key]["best_in_band"]
            cells.append(f"{lab} {'>' if x['beyond_resolution'] else ''}{x['sigma']:.1f}σ")
        print(f"  {name:12s} " + " | ".join(cells))
    hz = h["horizon"]
    print(f"Horizon problem: theta_H = {hz['theta_H_deg']:.2f} deg, "
          f"~{hz['patches_exchange']:,.0f} patches; inflation needs "
          f"{h['history']['n_min_horizon']:.1f} e-folds (Starobinsky, instant reheating)")
    print(f"Flatness: |Omega - 1| < {h['omega_dev_planck_temperature']:.0e} at the Planck "
          "temperature without inflation")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true", help="recompute every cached analysis step")
    ap.add_argument("--check", action="store_true", help="also run test_physics.py")
    args = ap.parse_args(argv)
    warnings.filterwarnings("ignore")

    import analysis
    import fetch_data
    import figures

    fetch_data.fetch_all()
    results = analysis.run(force=args.force)
    print("Figures")
    figures.main([])
    try:
        import explorer
    except ImportError:
        explorer = None
    if explorer is not None:
        print("Interactive page")
        explorer.main()
    summary(results)
    if args.check:
        import test_physics

        return test_physics.main()
    return 0


if __name__ == "__main__":
    sys.exit(main())
