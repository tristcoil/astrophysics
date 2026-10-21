"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run

The first run downloads ~0.7 GB (Gaia DR3 tables 0.32 GB, gaiamock's scan files
0.29 GB zipped / 0.99 GB unpacked, VizieR tables 0.05 GB), builds GSL for
gaiamock's C kernel, and simulates ~50,000 mock binaries through Gaia's
detection pipelines (about an hour on 6 cores on a busy machine; cached in
data/cache/ so later runs take a few minutes, offline).
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
from fetch_data import fetch_all

CACHE = C.CACHE_DIR / "analysis.pkl"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        print("Analysis ...", flush=True)
        res = analysis.run()
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(analysis.summarize(res), indent=1, default=_jsonable))
    paths = figures.make_all(res)
    analysis.print_summary(res)
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html")):
            print("  ", p)
    return 0


def _jsonable(o):
    import numpy as np
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, tuple):
        return list(o)
    return str(o)


if __name__ == "__main__":
    sys.exit(main())
