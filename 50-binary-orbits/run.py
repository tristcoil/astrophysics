"""Run the whole project end to end.

    python run.py                  # download what is missing, fit, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run

The first run downloads ~0.6 GB (570 MB of it Gaia's DR3 scanning law, used
only for the Sirius B parallax test) and takes ~20 minutes: the Sirius MCMC
(~10 min), the Gaia scan-law search (~15 s), 20,000 simulated Gaia solutions
(~1 min), the mass-luminosity MCMC (~2 min) and the animation (~3 min). Every
expensive step is cached in data/cache/, so later runs take a few minutes and
work offline.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import config as C
import figures
import report
from fetch_data import fetch_all

CACHE = C.CACHE_DIR / "analysis.pkl"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    ap.add_argument("--no-movie", action="store_true", help="skip the (slow) MP4 animation")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        res = analysis.run()
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(report.summarize(res), indent=1, default=float))
    paths = figures.make_all(res, movie=not args.no_movie)
    report.print_summary(res)
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html", ".mp4")):
            print("  ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
