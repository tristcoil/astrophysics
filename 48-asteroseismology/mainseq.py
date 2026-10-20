"""16 Cyg A and B: solar-like oscillations of two Sun-like stars.

16 Cygni A and B are a wide pair of near-solar twins 21 pc away, bright
enough (V ~ 6) that Kepler watched them with custom pixel masks at 1-minute
cadence. Their oscillations are the Sun's five-minute oscillations shifted to
slightly lower frequencies (the stars are a little larger than the Sun), and
because these are pure pressure modes, with none of the red giants' mixed
modes, their echelle diagrams show four clean ridges, l = 0, 1, 2 and 3.

The same pipeline as for the red giants runs on one year (Quarters 7-10) of
short-cadence data. Checks: Lund et al. (2017)'s peak-bagged frequencies and
global parameters (the Kepler LEGACY sample), and interferometric radii from
the CHARA array (White et al. 2013), which do not use seismology at all.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import lightcurves as L
from stars import global_parameters

MS_RESULTS = C.CACHE_DIR / "sixteen_cyg.json"
MS_ARRAYS = C.CACHE_DIR / "sixteen_cyg.npz"
MS_SEARCH = (500.0, 8000.0, 24)
MS_FIT_MIN = 30.0
STARS = {"16 Cyg A": 12069424, "16 Cyg B": 12069449}


def analyse_16cyg(*, force: bool = False, quiet: bool = False) -> dict:
    if MS_RESULTS.exists() and not force:
        return json.loads(MS_RESULTS.read_text())
    res, arrays = {}, {}
    for name, kic in STARS.items():
        ts = L.kepler_short_cadence(name)
        out, arr = global_parameters(ts, seed=kic, search=MS_SEARCH, fit_nu_min=MS_FIT_MIN,
                                     period_spacing=False, eps_width=0.5, n_perturb=100,
                                     dnu_method="acf")
        out["kic"] = kic
        res[name] = out
        key = name.replace(" ", "_").lower()
        for k, v in arr.items():
            arrays[f"{key}_{k}"] = v
        if not quiet:
            print(f"    {name}: numax {out['numax']:.1f} +- {out['numax_err']:.1f}  "
                  f"Dnu {out['dnu']:.3f} +- {out['dnu_err']:.3f} uHz", flush=True)
    np.savez_compressed(MS_ARRAYS, **arrays)
    MS_RESULTS.write_text(json.dumps(res, indent=1, default=float))
    return res
