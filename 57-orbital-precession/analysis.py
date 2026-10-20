"""Every number in the README, computed from the cached runs and data.

Writes ``data/results.json`` (all results, with their sources) and
``data/derived/series.npz`` (the time series the figures draw), so that
``figures.py`` never has to touch the raw runs.
"""

from __future__ import annotations

import json

import numpy as np

import binaries
import budget
import ephemeris as eph
import horizons_check
import ladder
import nbody
import orbit_physics as op
import pulsars
import s2
import config as C
from config import (ARCSEC_PER_RAD, CENTURY_DAYS, DERIVED_DIR, J2000_JD, PUBLISHED, RESULTS_PATH,
                    RUN, RUNS_DIR)

SERIES_PATH = DERIVED_DIR / "series.npz"


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


def perihelion_passages(t_days: np.ndarray, states: np.ndarray) -> dict:
    """Every perihelion passage in a 2-day-sampled series of heliocentric states.

    The minimum of r is located by a parabola through the three samples around
    each local minimum; the perihelion direction there is the osculating
    eccentricity vector, interpolated linearly.
    """
    r = np.linalg.norm(states[:, :3], axis=1)
    i = np.flatnonzero((r[1:-1] < r[:-2]) & (r[1:-1] <= r[2:])) + 1
    r0, r1, r2 = r[i - 1], r[i], r[i + 1]
    denom = r0 - 2 * r1 + r2
    frac = 0.5 * (r0 - r2) / denom  # offset of the minimum in samples
    dt = t_days[1] - t_days[0]
    t_peri = t_days[i] + frac * dt
    el = op.elements(states[:, :3].T, states[:, 3:].T, eph.gm("sun") + eph.gm("mercury"))
    eh = el["ehat"]
    j = np.where(frac >= 0, i, i - 1)
    w = np.where(frac >= 0, frac, frac + 1)
    e_at = eh[:, j] * (1 - w) + eh[:, j + 1] * w
    e_at /= np.linalg.norm(e_at, axis=0)
    q = (1 - el["e"][i]) * el["a"][i]
    return {"t_days": t_peri, "ehat": e_at.T, "q": q}


def mercury_series() -> dict:
    """Time series for the Mercury figures (decimated to 20 days)."""
    t = nbody.load("all_newton")["t"]
    de = budget.de440_series(t)
    out = {"t_cy": t[::10] / CENTURY_DAYS}
    for name in ("all_newton", "all_gr_full", "full_model"):
        s = budget.run_series(name)
        out[f"offset_{name}"] = budget.direction_offset(s, de)[::10] * ARCSEC_PER_RAD
        rec = nbody.load(name)
        r_de, _ = eph.heliocentric_state("mercury", J2000_JD + rec["t"][::10])
        out[f"km_{name}"] = np.linalg.norm(rec["mercury"][::10, :3].T - r_de, axis=0) * eph.au_km()
    trend = 575.0 / ARCSEC_PER_RAD * de["t_cy"]
    for k in ("psi", "varpi_ecl", "varpi_icrf"):
        y = de[k] - trend
        out[f"frame_{k}"] = (y - y[np.argmin(np.abs(t))])[::10] * ARCSEC_PER_RAD
    # Perihelion passages: DE440 and the Newton-only run.
    r_de, v_de = eph.heliocentric_state("mercury", J2000_JD + t)
    for label, states in (("de440", np.vstack([r_de, v_de]).T),
                          ("newton", nbody.load("all_newton")["mercury"]),
                          ("full", nbody.load("full_model")["mercury"])):
        pp = perihelion_passages(t, states)
        out[f"peri_t_{label}"] = pp["t_days"]
        out[f"peri_e_{label}"] = pp["ehat"]
        out[f"peri_q_{label}"] = pp["q"]
    # A reference orbit (DE440 at J2000) for drawing.
    out["orbit_j2000"] = np.vstack(eph.heliocentric_state("mercury", J2000_JD))
    return out


def millennia() -> dict | None:
    if not all((RUNS_DIR / f"{n}.npz").exists() for n in ("millennia_full", "millennia_newton")):
        return None
    res = {name: horizons_check.millennia_vs_de441(name)
           for name in ("millennia_full", "millennia_newton")}
    # How the rate itself changes: in-plane rate in 300-year windows.
    s = budget.run_series("millennia_full")
    t, psi = s["t_cy"], s["psi"]
    centres = np.arange(-28.5, 28.6, 1.0)
    rates = []
    for c in centres:
        m = np.abs(t - c) <= 1.5
        rates.append(op.fit_secular(t[m] - c, psi[m], deg=1)["rate"])
    res["rate_vs_epoch"] = {"year": (2000 + 100 * centres).tolist(), "rate": rates}
    res["full_span_fit"] = op.fit_secular(t, psi, deg=2)
    return res


def run() -> dict:
    print("  Mercury's budget (all runs) ...", flush=True)
    res = {"mercury": budget.analyse()}
    res["mercury"]["asteroids"] = budget.asteroid_belt_estimate()
    c = res["mercury"]["contributions"]
    c["full_model_plus_asteroids"] = c["full_model_total"] + res["mercury"]["asteroids"]["main_belt_rate"] \
        + res["mercury"]["asteroids"]["kuiper_belt_rate_at_44au"]
    c["residual_vs_de440_after_asteroids"] = c["full_model_plus_asteroids"] - c["de440_total"]
    pk = PUBLISHED["park2017"]
    c["park2017_total"] = pk["total_arcsec_cy"]
    c["park2017_minus_de440"] = pk["total_arcsec_cy"] - c["de440_total"]
    res["mercury"]["whfast"] = budget.whfast_test()
    ms = []
    for f in (1.0, 0.1, 0.01):
        name = "pair_jupiter" if f == 1.0 else f"pair_jupiter_x{f:g}"
        if (RUNS_DIR / f"{name}.npz").exists():
            rate = budget.difference_rate(budget.run_series(name),
                                          budget.run_series("two_body_newton"))["rate"]
            ms.append({"factor": f, "rate": rate, "rate_scaled_back": rate / f})
    res["mercury"]["mass_scaling"] = ms
    # Each planet's share defined as Park et al. (2017) do: leave it out of the full system.
    full = budget.run_series("all_newton")
    loo = {}
    for spec in C.LOO_RUNS:
        k = spec.name.replace("all_but_", "")
        if (RUNS_DIR / f"{spec.name}.npz").exists():
            loo[k] = budget.difference_rate(full, budget.run_series(spec.name))["rate"]
    res["mercury"]["leave_one_out"] = loo
    res["mercury"]["park_table3"] = C.PARK2017_TABLE3
    print("  Park et al.'s rate recipe on the six-millennium run ...", flush=True)
    if (RUNS_DIR / "millennia_full.npz").exists():
        import park_check
        res["mercury"]["park_recipe"] = park_check.analyse()
    print("  round trip 2000 -> 2650 -> 2000 ...", flush=True)
    res["mercury"]["reversibility"] = nbody.reversibility_test(649.5)
    print("  Horizons / DE441 checks ...", flush=True)
    res["horizons"] = {"elements": horizons_check.elements_check(),
                       "de440_vs_de441": horizons_check.de440_vs_de441(),
                       "millennia": millennia()}
    print("  binaries, pulsars, S2 ...", flush=True)
    res["binaries"] = binaries.analyse()
    res["di_her"] = binaries.di_her()
    res["pulsars"] = pulsars.analyse()
    sig = s2.precession_signal()
    res["s2"] = {k: v for k, v in sig.items() if not isinstance(v, np.ndarray)}
    res["s2"]["h1_public_data"] = s2.h1_results()
    res["s2"]["gravity"] = PUBLISHED["gravity_s2"]
    res["ladder"] = {"rungs": ladder.rungs(res["binaries"], res["pulsars"]),
                     "planets_predicted": ladder.predictions_only()}
    RESULTS_PATH.write_text(json.dumps(_jsonable(res), indent=1))

    print("  time series for the figures ...", flush=True)
    series = mercury_series()
    series.update({f"s2_{k}": v for k, v in sig.items() if isinstance(v, np.ndarray)})
    ros = s2.rosette()
    series["s2_rosette_sky"] = ros["sky"]
    series["s2_rosette_boost"] = np.array(ros["boost"])
    np.savez_compressed(SERIES_PATH, **series)
    return res


if __name__ == "__main__":
    run()
