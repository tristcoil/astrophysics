"""Run the whole project end to end.

    python run.py                  # download what is missing, fit light curves, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run

The first run downloads ~50 MB (210 MB once the light curves are unpacked) and
fits all 9,667 OGLE light curves (about 20-30 minutes on 6 cores; checkpointed,
so an interrupted run resumes). After that the analysis takes about two minutes
offline.
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
import lightcurves
from fetch_data import fetch_all

CACHE = C.CACHE_DIR / "analysis.pkl"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    ap.add_argument("--workers", type=int, default=6, help="processes for the light-curve fits")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        lightcurves.run_all(args.workers)
        print("Analysis ...", flush=True)
        res = analysis.run()
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(analysis.summarize(res), indent=1, default=float))
    paths = figures.make_all(res)
    analysis.print_summary(res)
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html")):
            print("  ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
