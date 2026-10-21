"""Run the three analyses and save every number to data/results.json.

Each part returns (results, artefacts): the results are plain numbers that go
into results.json and the README; the artefacts are the arrays the figures
need. Nothing here downloads anything.
"""

from __future__ import annotations

import json
import time

import numpy as np

from config import GD1_G_MAX, RESULTS_PATH


def _clean(obj):
    """Make results JSON-serialisable (numpy scalars, tuples, private keys)."""
    if isinstance(obj, dict):
        return {str(k): _clean(v) for k, v in obj.items() if not str(k).startswith("_")}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _clean(obj.tolist())
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, str):
        return str(obj)
    return obj


def run(parts=("gd1", "gse", "spiral"), n_boot: int = 200, quiet: bool = False):
    """Run the requested parts; returns ({part: results}, {part: artefacts})."""
    results, artefacts = {}, {}
    for part in parts:
        t0 = time.time()
        if part == "gd1":
            import gd1

            results[part], artefacts[part] = gd1.analyse(GD1_G_MAX, n_boot=n_boot)
        elif part == "gse":
            import halo

            results[part], artefacts[part] = halo.analyse()
        elif part == "spiral":
            import spiral

            results[part], artefacts[part] = spiral.analyse(n_boot=min(n_boot, 100))
        else:
            raise ValueError(part)
        if not quiet:
            print(f"  {part}: analysed in {time.time() - t0:.0f} s", flush=True)
    return results, artefacts


def save_results(results: dict) -> None:
    """Merge into results.json, so running one part keeps the others."""
    old = json.loads(RESULTS_PATH.read_text()) if RESULTS_PATH.exists() else {}
    old.update(_clean(results))
    RESULTS_PATH.write_text(json.dumps(old, indent=1) + "\n")
