"""Run the whole project end to end.

    python run.py                 # fetch (if needed), analysis, figures, explorer, summary
    python run.py --figures-only  # redraw figures and the explorer from data/results.json
    python run.py --no-fetch      # skip the download step (everything already cached)
    python run.py --force         # recompute the cached pieces (dwarf profiles, Higgs portal)

The first run downloads ~26 MB (less when P2's and P5's caches exist) and takes ~35 min, most of
it computing per-dwarf likelihood profiles and the Higgs-portal line; after that ~5 min offline.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import config as C  # noqa: E402


def print_summary(r: dict) -> None:
    h = r["thermal"]["headline"]
    v = r["thermal"]["validation_ss2020"]["majorana_s"]
    print("\n" + "=" * 78)
    print("P4 - freeze-out: why there is dark matter (if it is a WIMP)")
    print("=" * 78)
    print(f"thermal <sigma v> for Omega_c h^2 = 0.120 (Majorana, s-wave): {h['sigv_100GeV']:.3e} cm^3/s at 100 GeV, "
          f"{h['range_10GeV_10TeV'][0]:.2e}-{h['range_10GeV_10TeV'][1]:.2e} for 10 GeV-10 TeV; "
          f"peak {h['low_mass_peak']['sigv']:.2e} at {h['low_mass_peak']['mass']:.2f} GeV")
    print(f"  vs Saikawa & Shirai 2020: median ratio {v['median_ratio']:.4f} ({v['min_ratio']:.4f}-{v['max_ratio']:.4f})"
          f"; freeze-out x_f = {h['x_f_100GeV']:.1f} at 100 GeV")
    d = r["thermal"]["decomposition_100GeV"]
    print(f"  SDB12 set-up {d[0]['sigv_100GeV']:.3e} -> today's {d[-1]['sigv_100GeV']:.3e} (Omega h^2 0.11 -> 0.120 "
          f"dominates)")
    a = r["analytic"]
    print(f"Kolb & Turner's estimate: thermal <sigma v> {a['thermal_ratio_kt'][0]:.3f}-{a['thermal_ratio_kt'][1]:.3f} "
          f"x numerical")
    n = r["neutrons"]["runs"]["escudero"]
    pr = r["neutrons"].get("prymordial", {})
    print(f"neutrons: freeze-out at T = {n['T_freeze_MeV']:.2f} MeV (t = {n['t_freeze_s']:.2f} s), n/p(1 s) = "
          f"{n['np_1s']:.3f}; helium forms at T = {n['T_nuc_keV']:.0f} keV (t = {n['t_nuc_s']:.0f} s)")
    line = f"  Y_p = {n['Y_p']:.4f}, D/H = {n['D_H']:.3e} (Born weak rates, 8 reactions)"
    if "born_small_instant_a" in pr:
        line += f"; PRyMordial same physics {pr['born_small_instant_a']['Y_p']:.4f}"
    if "full" in pr:
        line += f", full {pr['full']['Y_p']:.4f}"
    print(line)
    rc = r["recombination"]
    print(f"recombination: electrons freeze out at z = {rc['z_freeze']:.0f}; x_e(z=200) = {rc['x_e_200']:.3e} "
          f"(CAMB {rc['x_e_200_camb']:.3e}; Saha would give {rc['saha_at_200']:.1e})")
    L = r["dwarfs"]["limits"]

    def masses(key):
        return ", ".join(f"{m:.0f}" for m in L[key]["crossing_GeV"]) or "none"

    for ch in ("bb", "tau"):
        print(f"Fermi-LAT dwarfs ({ch}): thermal WIMPs excluded below {masses(f'{ch}_pppc_benchmark')} GeV "
              f"(42 dwarfs; their own profiles {masses(f'{ch}_mcdsum_benchmark')} GeV, "
              f"Ackermann+15 {masses(f'{ch}_ackermann2015')} GeV)")
    hp = r["direct"]["higgs_portal"]
    windows = ", ".join(f"{lo:g}-{hi:g}" for lo, hi in hp["allowed_by_lz2024"]) or "none"
    print(f"Higgs-portal singlet: allowed by LZ 2024 only for m = {windows} GeV; "
          f"lambda(1 TeV) = {hp['lambda_1TeV']:.3f}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--no-fetch", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args(argv)
    if args.figures_only:
        results = json.loads(C.RESULTS_PATH.read_text())
    else:
        if not args.no_fetch:
            from fetch_data import fetch_all

            fetch_all()
        import analysis

        results = analysis.run_core(force=args.force)
    print_summary(results)
    import explorer
    import figures

    for p in figures.make_all(results):
        print("  wrote", Path(p).relative_to(C.HERE))
    print("  wrote", explorer.build(results).relative_to(C.HERE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
