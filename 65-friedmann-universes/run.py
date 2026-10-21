"""Build the whole project with one command.

    python run.py            # fetch what is missing, analysis (cached steps reused),
                             # every figure, the interactive page, and a summary
    python run.py --force    # recompute every analysis step (~8 min)
    python run.py --check    # also run the self-checks in test_physics.py (~1 min)

Once ``data/`` is filled this runs offline: nothing cached is fetched again.
"""

from __future__ import annotations

import argparse
import sys
import warnings


def summary(r: dict) -> None:
    c, f, o = r["combos"], r["flat"], r["ours"]
    ep = o["epochs"]
    print("\n" + "=" * 78)
    print("O1 — Build a universe: the Friedmann equations")
    print("=" * 78)
    print(f"Supernovae alone (D1's Pantheon+): Om = {r['sn']['Om']['mean']:.3f} ± "
          f"{r['sn']['Om']['std']:.3f}, OL = {r['sn']['OL']['mean']:.3f} ± {r['sn']['OL']['std']:.3f}")
    b = r["bao"]
    print(f"DESI BAO alone: Ok = {b['Ok']['mean']:+.4f} ± {b['Ok']['std']:.4f} "
          f"(DESI's chain {b['desi_chain']['Ok']['mean']:+.4f} ± {b['desi_chain']['Ok']['std']:.4f})")
    for key in ("bao+cmb", "sn+bao+cmb"):
        s = c[key]
        print(f"{s['label']:>16}: Ok = {s['Ok']['mean']:+.4f} ± {s['Ok']['std']:.4f}, "
              f"H0 = {s['H0']['mean']:.2f} ± {s['H0']['std']:.2f}, age = "
              f"{s['age_gyr']['mean']:.3f} ± {s['age_gyr']['std']:.3f} Gyr")
    s = f["sn+bao+cmb"]
    print(f"   flat, all three: age = {s['age_gyr']['mean']:.3f} ± {s['age_gyr']['std']:.3f} Gyr "
          "(Planck 2018 + BAO: 13.787 ± 0.020)")
    print(f"Epochs: matter-radiation equality z = {ep['z_eq']:.0f} ({ep['t_eq_yr'] / 1e3:.1f} kyr), "
          f"CMB z = {ep['z_star']:.1f} ({ep['t_star_yr'] / 1e3:.0f} kyr), acceleration from "
          f"z = {ep['z_acc']:.2f} ({ep['lookback_acc_gyr']:.2f} Gyr ago)")
    m = r["map"]
    print(f"No bouncing universe can show anything beyond z = {m['max_bounce_z_radiation']:.2f} "
          f"(radiation alone), {m['max_bounce_z_atoms']:.2f} (with atoms' Om >= 0.049)")
    a = r["ages"]
    print(f"Age crisis: Einstein-de Sitter with H0 = {a['cepheids']['H0']:.1f} is "
          f"{a['cepheids']['eds_age_gyr']:.2f} Gyr, {a['cepheids']['eds_vs_stars_sigma']:.1f} sigma "
          f"younger than the oldest stars ({a['oldest_stars']['t_min_gyr']:.1f} ± "
          f"{a['oldest_stars']['t_min_err']:.1f} Gyr, B6)")


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
