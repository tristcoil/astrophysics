"""Use project H1's own code (17-galactic-centre) for S2's data and orbit model.

H1 fitted S2's orbit; N6 tests the conserved quantities along it. Re-using
H1's modules guarantees the same data selection, sign conventions, reference
frames and orbit model (Kepler + light-travel delay + relativistic redshift +
Schwarzschild precession). H1's modules are imported under their own names
(``config``, ``orbit_physics``, ...), which would clash with ours, so they are
loaded with H1's folder first on the path and then moved out of
``sys.modules`` again; the module objects keep working afterwards.
"""

from __future__ import annotations

import importlib
import sys
from functools import lru_cache
from types import SimpleNamespace

import numpy as np

from config import H1_CHAIN, H1_DIR

_H1_NAMES = ("config", "orbit_physics", "s2data", "fitting", "analysis", "black_hole_case")


@lru_cache(maxsize=1)
def modules() -> SimpleNamespace:
    if not H1_DIR.exists():
        raise FileNotFoundError(f"project H1 not found at {H1_DIR}")
    saved = {m: sys.modules.pop(m) for m in _H1_NAMES if m in sys.modules}
    sys.path.insert(0, str(H1_DIR))
    try:
        ns = SimpleNamespace(op=importlib.import_module("orbit_physics"),
                             sd=importlib.import_module("s2data"),
                             ft=importlib.import_module("fitting"),
                             an=importlib.import_module("analysis"))
    finally:
        sys.path.remove(str(H1_DIR))
        for m in _H1_NAMES:
            sys.modules.pop(m, None)
        sys.modules.update(saved)
    return ns


@lru_cache(maxsize=1)
def chain() -> dict:
    """H1's joint Keck + VLT GR posterior: {name: samples}, plus the best fit."""
    z = np.load(H1_CHAIN)
    names = [str(n) for n in z["names"]]
    out = {n: z["chain"][:, k] for k, n in enumerate(names)}
    out["_best"] = dict(zip(names, z["best_theta"]))
    out["_names"] = names
    out["_n"] = z["chain"].shape[0]
    return out


@lru_cache(maxsize=1)
def model():
    """H1's joint GR model object (frames, data, likelihood)."""
    h = modules()
    return h.an.build_model("joint_gr")
