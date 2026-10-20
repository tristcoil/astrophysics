"""Regular versus chaotic orbits: what Kerr's hidden symmetry buys.

Kerr and the Zipoy-Voorhees (ZV) spacetime are both exact vacuum solutions of
Einstein's equations, both static/stationary and axially symmetric, so both
have the same two Killing vectors and conserve E and L. What is left of a
geodesic is motion in the meridional plane (r and theta): two degrees of
freedom held by one constraint (the mass shell). A third conserved quantity
makes the motion *integrable*: every orbit then winds regularly around a
torus, and a Poincare section (the orbit's position and radial momentum each
time it crosses the equator upwards) shows smooth closed curves. Kerr has the
third quantity, Carter's constant from its Killing tensor. ZV does not (``jets.py``
finds only the four reducible Killing tensors), and orbits that come close to
the unstable circular orbit fill parts of the section chaotically
(Lukes-Gerakopoulos 2012, Phys. Rev. D 86, 044013).

For a fair comparison each spacetime gets an angular momentum chosen so that
the top of its centrifugal barrier sits at the same energy (0.995 c^2 per unit
mass), and orbits with an energy at the same fraction of the way from the
bottom of the potential well to the barrier top, close to the top, where they
graze the unstable circular orbit.

A second, quantitative test: two orbits started 1e-8 apart. In a regular
system they drift apart roughly linearly; in a chaotic one exponentially.
"""

from __future__ import annotations

import multiprocessing as mp
import pickle

import numpy as np
import sympy as sp
from scipy import optimize

import config as C
import geometry as geo
import orbits as ob

#: height of the centrifugal barrier, per unit rest mass, in both spacetimes
V_TOP = 0.995
#: energy position between well bottom (0) and barrier top (1)
ENERGY_FRACTION = 0.93
N_ORBITS = 12
CROSSINGS = 160

CASES = {
    "Kerr": {"label": "Kerr, a = 0.9 M", "mass": 1.0},
    "ZV": {"label": "Zipoy-Voorhees, delta = 2", "mass": 2.0},  # sigma = 1 -> M = delta sigma
}
_H = {}


def _metric(name):
    if name == "Kerr":
        return geo.kerr(), {geo.M: 1, geo.a: sp.Rational(9, 10)}
    return geo.zipoy_voorhees(2), {}


def hamiltonian(name):
    if name not in _H:
        met, sub = _metric(name)
        _H[name] = ob.AxisymmetricHamiltonian(met, sub)
    return _H[name]


def equatorial_potential(name, q1, L):
    """E for which a body at rest in the equatorial plane (p = 0) is on the mass shell."""
    H = hamiltonian(name)
    return np.array([H.energy_at_rest(q, 0.0, L) for q in np.atleast_1d(q1)])


def _barrier_and_well(name, L):
    mass = CASES[name]["mass"]
    q_in = 1.05 if name == "ZV" else 1.44   # just outside the horizon (r+ = 1.436 M for a = 0.9)
    q = np.geomspace(q_in, 200 * mass, 6000)
    V = equatorial_potential(name, q, L)
    ok = np.isfinite(V)
    q, V = q[ok], V[ok]
    d = np.diff(V)
    ext = np.where(np.sign(d[1:]) != np.sign(d[:-1]))[0] + 1
    maxima = [i for i in ext if d[i - 1] > 0]
    minima = [i for i in ext if d[i - 1] < 0]
    if not maxima or not minima:
        return None
    i_max = maxima[0]
    i_min = [i for i in minima if i > i_max][0]
    return q, V, i_max, i_min


def well(name) -> dict:
    """Pick L so the barrier top is at V_TOP, then E at ENERGY_FRACTION of the well depth."""
    mass = CASES[name]["mass"]

    def top(L):
        bw = _barrier_and_well(name, L)
        return (bw[1][bw[2]] - V_TOP) if bw else -1.0

    L = optimize.brentq(top, 3.0 * mass if name == "ZV" else 2.15, 6.0 * mass if name == "ZV" else 3.5, xtol=1e-10)
    q, V, i_max, i_min = _barrier_and_well(name, L)
    vmax, vmin = V[i_max], V[i_min]
    E = vmin + ENERGY_FRACTION * (vmax - vmin)
    f = lambda s_: equatorial_potential(name, s_, L)[0] - E
    inner = optimize.brentq(f, q[i_max], q[i_min], xtol=1e-12)
    outer = optimize.brentq(f, q[i_min], q[-1], xtol=1e-12)
    return {"L": float(L), "L_over_M": float(L / mass), "E": float(E), "q_barrier": float(q[i_max]),
            "V_barrier": float(vmax), "q_well": float(q[i_min]), "V_well": float(vmin),
            "inner": float(inner), "outer": float(outer), "q": q[::20].tolist(), "V": V[::20].tolist()}


def _section(args):
    name, q1, E, L, lam_max = args
    H = hamiltonian(name)
    p2 = H.p2_for(q1, 0.0, 0.0, E, L)
    if p2 is None:
        return name, q1, None
    sol = H.integrate(np.array([q1, 0.0, 0.0, p2]), E, L, lam_max)
    return name, q1, sol.y_events[0][:, [0, 2]]


def _divergence(args):
    name, q1, E, L, lam_max, n_out = args
    H = hamiltonian(name)
    out = []
    lam = np.linspace(0, lam_max, n_out)
    sols = []
    for dq in (0.0, 1e-8):
        p2 = H.p2_for(q1 + dq, 0.0, 0.0, E, L)
        sol = H.integrate(np.array([q1 + dq, 0.0, 0.0, p2]), E, L, lam_max, dense=True)
        sols.append(sol.sol(lam))
    d = np.sqrt(np.sum((sols[0] - sols[1]) ** 2, axis=0))
    out = {"lam": lam.tolist(), "dist": d.tolist()}
    return name, q1, out


def run(force: bool = False, workers: int = 4) -> dict:
    cache = C.CACHE_DIR / "chaos.pkl"
    if cache.exists() and not force:
        return pickle.loads(cache.read_bytes())
    wells = {n_: well(n_) for n_ in CASES}
    tasks, dtasks = [], []
    for n_, w in wells.items():
        # radial period estimate from a short test run, to size the integrations
        H = hamiltonian(n_)
        q0 = 0.5 * (w["inner"] + w["outer"])
        test = H.integrate(np.array([q0, 0.0, 0.0, H.p2_for(q0, 0.0, 0.0, w["E"], w["L"])]), w["E"], w["L"],
                           4000.0 * CASES[n_]["mass"])
        per = 4000.0 * CASES[n_]["mass"] / max(len(test.t_events[0]), 1)
        lam_max = per * CROSSINGS
        span = w["outer"] - w["inner"]
        for q1 in np.linspace(w["inner"] + 0.01 * span, w["outer"] - 0.01 * span, N_ORBITS):
            tasks.append((n_, float(q1), w["E"], w["L"], lam_max))
        for q1 in (w["inner"] + 0.05 * span, w["inner"] + 0.5 * span):
            dtasks.append((n_, float(q1), w["E"], w["L"], per * 60, 300))
        w["lam_per_crossing"] = per
    for n_ in CASES:   # build in the parent so forked workers inherit them
        hamiltonian(n_)
    with mp.get_context("fork").Pool(workers) as pool:
        secs = pool.map(_section, tasks)
        divs = pool.map(_divergence, dtasks)
    out = {"wells": wells, "sections": {n_: [] for n_ in CASES}, "divergence": {n_: [] for n_ in CASES},
           "settings": {"V_top": V_TOP, "energy_fraction": ENERGY_FRACTION, "n_orbits": N_ORBITS,
                        "crossings": CROSSINGS}}
    for n_, q1, pts in secs:
        if pts is not None:
            out["sections"][n_].append({"q1": q1, "points": pts.tolist()})
    for n_, q1, d in divs:
        dist = np.array(d["dist"])
        lam = np.array(d["lam"])
        good = dist > 0
        # growth: fit log(dist) vs lam over the second half (exponential) and compare with
        # log(dist) vs log(lam) (power law, as for regular orbits)
        half = good & (lam > lam.max() / 2)
        slope_exp = float(np.polyfit(lam[half], np.log(dist[half]), 1)[0]) if half.sum() > 5 else float("nan")
        d.update({"q1": q1, "final_over_initial": float(dist[-1] / 1e-8), "exp_rate_per_lam": slope_exp})
        out["divergence"][n_].append(d)
    cache.write_bytes(pickle.dumps(out))
    return out
