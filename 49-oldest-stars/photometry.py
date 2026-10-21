"""From a HUGS catalogue to a clean, dereddened CMD and its ridge line.

    prepare(name) -> dict with the selected stars, their differential-reddening
                     corrections, and the ridge line of F606W vs F606W - F814W.

The steps, in order:
1. well-measured stars in F606W and F814W (catalog.quality_mask) that are
   proper-motion members (P >= 90%);
2. a first ridge line, started from a MIST isochrone slid onto the data;
3. differential reddening from main-sequence reference stars (catalog.
   differential_reddening), using the extinction ratios of a turnoff-like star
   from the MIST tables;
4. the final ridge line of the corrected CMD.
The results are cached in data/derived/prepared_<hugs id>.npz.
"""

from __future__ import annotations

import json

import numpy as np

import catalog as K_
import config as C
import fiducial as F
from isochrones import Extinction, load_grid


def _rough_turnoff(m: np.ndarray, c: np.ndarray) -> float:
    """Magnitude of the bluest median colour along the upper main sequence."""
    lo, hi = np.percentile(m, [5, 95])
    edges = np.arange(lo, hi, 0.1)
    best, m_to = np.inf, np.nan
    for a, b in zip(edges[:-1], edges[1:]):
        s = (m >= a) & (m < b)
        if s.sum() < 50:
            continue
        med = np.median(c[s])
        if med < best:
            best, m_to = med, 0.5 * (a + b)
    return float(m_to)


def prepare(name: str, *, use_cache: bool = True) -> dict:
    inputs = K_.cluster_inputs(name) if name in C.CLUSTERS else None
    return prepare_hugs(inputs["hugs"], inputs, use_cache=use_cache)


def prepare_hugs(hugs_id: str, inputs: dict, *, use_cache: bool = True) -> dict:
    cache = C.DERIVED_DIR / f"prepared_{hugs_id}.npz"
    if use_cache and cache.exists():
        z = np.load(cache, allow_pickle=False)
        out = {k: z[k] for k in z.files if k != "meta"}
        out["meta"] = json.loads(str(z["meta"]))
        out["ridge"] = F.Ridge(*(out.pop(f"ridge_{f}") for f in ("s", "c", "m", "sig", "nx", "ny", "n")),
                               i_to=int(out["meta"]["ridge_i_to"]), width=out.pop("ridge_width"))
        return out

    cat = K_.read_hugs(hugs_id)
    q606 = K_.quality_mask(cat, ("F606W", "F814W"))
    mem = K_.members(cat)
    sel = q606 & mem
    m606, m814 = cat["m_F606W"][sel], cat["m_F814W"][sel]
    x, y = cat["x"][sel], cat["y"][sel]
    c = m606 - m814

    # extinction ratios of a turnoff-like star (log Teff 3.79, log g 4.2)
    ext = Extinction(inputs["feh"], 0.4)
    a606 = ext.ratio("F606W", 3.79, 4.2)
    a814 = ext.ratio("F814W", 3.79, 4.2)

    # first guide: MIST at 12.5 Gyr placed at the nominal distance and reddening
    grid = load_grid("MIST")
    iso = grid.interpolate(inputs["feh"], 0.4, [12.5])
    gm = iso["F606W"][0] + inputs["dm"] + a606 * inputs["ebv"]
    gc = iso["F606W"][0] - iso["F814W"][0] + (a606 - a814) * inputs["ebv"]
    ok = np.isfinite(gm) & np.isfinite(gc)
    m_to0 = _rough_turnoff(m606, c)
    win = ok & (gm < m_to0 + 3.6) & (gm > m_to0 - 4.0)
    gm, gc = F.slide_guide(m606, c, gm[win], gc[win])

    # differential reddening from main-sequence reference stars
    r0 = F.ridge_line(m606, c, gm, gc)
    s_star, t_star = F.project(F.K * r0.c, r0.m, F.K * c, m606)
    ms = (m606 > r0.m_to + 0.5) & (m606 < r0.m_to + 3.0) & (abs(t_star) < 0.10)
    dE = K_.differential_reddening(m606, c, x, y, ms, a606, a606 - a814)
    m606c, m814c = m606 - a606 * dE, m814 - a814 * dE
    cc = m606c - m814c

    ridge = F.ridge_line(m606c, cc, r0.m, r0.c)
    c_to, m_to, m_sgb = F.turnoff_and_sgb_point(ridge.c, ridge.m)
    meta = {
        "hugs": hugs_id, "n_catalogue": int(len(cat["x"])), "n_quality": int(q606.sum()),
        "n_members": int(mem.sum()), "n_selected": int(sel.sum()), "n_reference_dr": int(ms.sum()),
        "a606_per_ebv": a606, "a814_per_ebv": a814,
        "dE_rms": float(np.std(dE)), "dE_p95": float(np.percentile(abs(dE), 95)),
        "c_to": c_to, "m_to": m_to, "m_sgb": m_sgb, "ridge_i_to": int(ridge.i_to),
    }
    out = {
        "m606": m606c, "m814": m814c, "c": cc, "m606_raw": m606, "c_raw": c, "x": x, "y": y, "dE": dE,
        "ms_ref": ms, "ra": cat["ra"][sel], "dec": cat["dec"][sel], "sel_index": np.flatnonzero(sel),
    }
    for band in ("F275W", "F336W", "F438W"):
        out[f"m_{band}_raw"] = cat[f"m_{band}"][sel]
    np.savez_compressed(cache, **out, meta=json.dumps(meta),
                        **{f"ridge_{f}": getattr(ridge, f) for f in ("s", "c", "m", "sig", "nx", "ny", "n",
                                                                       "width")})
    out["meta"] = meta
    out["ridge"] = ridge
    return out
