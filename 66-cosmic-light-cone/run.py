"""Run the whole project end to end.

    python run.py                 # fetch (if needed), analysis, figures, explorer, summary
    python run.py --figures-only  # redraw figures and the explorer from data/results.json
    python run.py --no-fetch      # skip the download step (everything already cached)

The first run downloads ~5 MB when the sibling projects (D1, D2, D5, I3) have
their caches, ~90 MB otherwise. After that it runs offline in about a minute.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import config as C  # noqa: E402


def print_summary(r: dict) -> None:
    s = r["headline"]["summary"]
    u = r["uncertainties"]["planck18"]
    print("\n" + "=" * 78)
    print("O2 — the cosmic light cone (Planck 2018)")
    print("=" * 78)
    print(f"age {s['age_gyr']:.3f} Gyr; Hubble distance c/H0 {s['hubble_distance_gly']:.2f} billion ly")
    print(f"particle horizon {s['particle_horizon_gly']:.2f} ± {u['particle_horizon_gly']['sigma']:.2f} "
          f"billion ly (spec 46.2); event horizon {s['event_horizon_gly']:.2f} ± "
          f"{u['event_horizon_gly']['sigma']:.2f} (spec 16.6)")
    print(f"recession = c today at z = {s['z_hubble_sphere']:.3f} (1.48); smallest angular size at "
          f"z = {s['z_min_angular_size']:.3f} (1.59); event horizon at z = {s['z_event_horizon']:.3f} (1.85)")
    print(f"reachable share of the observable universe {100 * s['reachable_fraction']:.2f}% (about 5%)")
    c = r["headline"]["cmb"]
    print(f"CMB: emitted {c['d_emit_gly'] * 1000:.0f} million ly away at {c['t_emit_gyr'] * 1e6:.0f} thousand "
          f"years, now {c['chi_gly']:.1f} billion ly; receding {c['v_then_c']:.1f}c then, {c['v_now_c']:.2f}c now")
    d = r["desi"]["rows"][-1]
    print(f"DESI Ly-alpha z = 2.33: receding {d['v_then_c']:.3f} ± {d['v_then_sigma']:.3f} c when its light "
          f"left ({d['sigma_from_c']:.1f} sigma above c), no ruler size or H0 needed")
    b = r["sn_velocities"]["bins"][-1]
    print(f"Pantheon+ beyond z = 1.48 ({b['n']} SNe): {b['v_now_c']:.3f} ± {b['v_now_sigma']:.3f} c today")
    rep = r["sr_test"]["dl2004_reproduction"]["rows"]
    pm = r["sr_test"]["pantheon"]["models"]
    print(f"Doppler reading of redshift: D&L's drawn curve {rep['SR_as_plotted']['sigma_fixedM_scaled_err']:.1f} "
          f"sigma (they quote 23); Pantheon+ {pm['SR']['sigma']:.1f} sigma; empty universe "
          f"{pm['milne']['sigma']:.1f} sigma")


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
