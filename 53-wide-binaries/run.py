"""Run the whole project end to end.

    python run.py                 # fetch what is missing, analyse, draw everything
    python run.py --no-figures
    python run.py --quick         # skip the forward-model variants (the slowest part)

First run: ~2.1 GB of file downloads (El-Badry et al.'s catalogues at ~6 MB/s
by parallel byte ranges, Hwang et al.'s eccentricities, Chae's samples, SFD
dust maps) plus 1.5 GB of Gaia DR3 query results for 8.9 M nearby stars in 63
sky patches, then roughly 8-10 hours of computing on a busy 16-core machine
(orbit libraries ~1 min each, each model library ~2-5 min, each forward fit +
MCMC ~15-25 min, each run of Chae's method on a mock ~4 min). Every step is
cached under data/derived/, so later runs take minutes.
"""

from __future__ import annotations

import argparse
import time

import numpy as np


def _rule(t):
    print(f"\n{t}\n" + "-" * len(t), flush=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--quick", action="store_true", help="skip the forward-model variants")
    args = ap.parse_args(argv)
    t0 = time.time()

    _rule("1. Data")
    import fetch_data
    fetch_data.ensure_data()
    import gaia_queries as GQ
    GQ.nearby_lean()
    import pairs as PR
    real, flow = PR.cached_pairs(False)
    sh, flow_s = PR.cached_pairs(True)
    ids = np.concatenate([real["source_id_a"][real["keep"]], real["source_id_b"][real["keep"]],
                          sh["source_id_a"][sh["keep"]], sh["source_id_b"][sh["keep"]]])
    GQ.full_rows("ps23_stars_full", ids)
    import elbadry
    elbadry.ensure_slim()
    print(f"  data ready ({time.time() - t0:.0f} s)")

    _rule("2. Physics: MOND force and orbit libraries")
    import qumond as Q
    import orbits as OB
    gf = Q.galactic_field()
    Q.cached_field(gf["gNe_a0"])
    for a in OB.ALPHA_GRID:
        OB.build_library(float(a))
    print(f"  g_e = {gf['ge_a0']:.3f} a0, g_Ne = {gf['gNe_a0']:.3f} a0, eta_inf = {gf['eta_inf']:.4f}")

    _rule("3. Analyses")
    import analysis as A
    r = A.forward_banik2024()
    a = r["posterior"]["alpha"]
    print(f"  forward-model test (Banik-style sample, {r['n_pairs_in_region']:,} pairs): alpha = {a['median']:.3f} "
          f"(+{a['hi'] - a['median']:.3f} -{a['median'] - a['lo']:.3f}); Newton over MOND {r['sigma_newton_over_mond']:.1f} sigma")
    inj = A.injection_chae()
    for k, v in inj["runs"].items():
        print(f"  Chae's method on {k:34s}: g_obs/g_pred = {np.round(v['gamma_g'], 3)}")
    # Chae's method: its choices on the real sky, control mocks, the mechanism, subsets
    A.chae_variants()
    A.injection_controls()
    A.injection_mechanism()
    A.injection_mechanism_opik()
    A.injection_subsets()
    A.noise_inflation()
    # the eccentricity law: fixed-gamma fits and the v-r angle check
    A.gamma_profile()
    A.vr_angle_test()
    A.vr_angle_noise()
    A.alpha_at_fixed_gamma()
    A.median_profiles()
    A.profile_tests()
    # 3D velocities
    A.three_d()
    A.three_d_robustness()
    A.three_d_saad_ting_like()
    # the forward model on Chae's sample, and mocks drawn from that fit
    A.forward_chae_all()
    A.alpha_at_fixed_gamma_chae()
    A.injection_chae_v2()
    if not args.quick:
        for name in A.VARIANTS:
            if name != "nominal":
                A.run_variant(name)
    A.write_summary()

    if not args.no_figures:
        _rule("4. Figures")
        import figures
        for name, fn in figures.FIGURES.items():
            out = fn()
            print(f"  {name}: {out[0] if out else 'skipped (inputs missing)'}")
        import interactive
        print(f"  explorer: {interactive.build()}")
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min. Results: data/derived/results/, summary data/results.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
