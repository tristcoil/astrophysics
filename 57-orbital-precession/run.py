"""Rebuild the whole project with one command.

    python run.py              # fetch (if needed), simulate (if needed), analyse, draw
    python run.py --no-movie   # skip the 3-minute animation
    python run.py --workers 6  # parallel N-body runs (default 4)

Stages, each cached so a re-run only redoes what is missing:

1. ``fetch_data``  -- DE440 (120 MB), Horizons, literature tables, asteroid orbits
2. ``nbody``       -- 32 REBOUND/REBOUNDx integrations (~25 min on 4 cores, mostly the
                      six-millennium full-physics run; ~5 min without it)
3. ``analysis``    -- every number in the README -> data/results.json (~30 s)
4. ``figures``     -- figures/*.png and *.svg (~30 s)
5. ``animate``     -- figures/mercury_rosette.{mp4,gif,png} (~3 min)
"""

from __future__ import annotations

import argparse
import json
import time

import config as C


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-movie", action="store_true")
    ap.add_argument("--no-fetch", action="store_true", help="assume data/raw is complete (offline)")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args(argv)
    t0 = time.time()

    if not args.no_fetch:
        print("1. inputs", flush=True)
        import fetch_data
        fetch_data.main([])

    print("2. N-body integrations (cached in data/derived/runs/)", flush=True)
    import nbody
    names = [r.name for r in C.RUNS + C.LOO_RUNS + C.WHFAST_RUNS + C.MASS_SCALE_RUNS
             + C.MILLENNIA_RUNS]
    nbody.run_all(names, workers=args.workers)

    print("3. analysis", flush=True)
    import analysis
    res = analysis.run()

    print("4. figures", flush=True)
    import figures
    figures.main()

    if not args.no_movie:
        print("5. animation", flush=True)
        import animate
        print("  ", animate.make())

    c = res["mercury"]["contributions"]
    print(f"\nMercury, in-plane perihelion rate at J2000 (arcsec/century):")
    print(f"  planets (Newton)      {c['all_newton']['rate']:10.4f}")
    print(f"  general relativity    {c['gr_in_system']['gr_full']['rate']:10.4f}")
    print(f"  Sun's J2 + spin       {c['j2_plus_lt_in_system']['rate']:10.4f}")
    print(f"  our total             {c['full_model_total']:10.4f}   (+ asteroid estimate: "
          f"{c['full_model_plus_asteroids']:.4f})")
    print(f"  DE440                 {c['de440_total']:10.4f}")
    print(f"  MESSENGER (Park+2017) {c['park2017_total']:10.4f}")
    print(f"\ndone in {(time.time() - t0) / 60:.1f} min; results in {C.RESULTS_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
