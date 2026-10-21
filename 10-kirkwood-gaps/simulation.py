"""EXTENSION -- A LABELLED SIMULATION: watch Jupiter carve the gaps.

The measured histogram shows *that* asteroids are missing at the resonances.
This N-body experiment shows *how*. It is a simulation, not data, and every
figure made from it says so.

Set-up
------
* **The planets are real.** The Sun, Jupiter and Saturn start from their JPL
  Horizons positions and velocities on the snapshot date. Saturn matters:
  its tug makes Jupiter's orbit precess, and that precession is what turns
  the resonances chaotic (the nu5 and nu6 secular resonances overlap them).
* **The test particles are real asteroids, moved.** Each one takes the
  eccentricity, inclination and orbital angles of a randomly chosen real
  main-belt asteroid, but its semi-major axis is drawn *uniformly* from
  2.0-3.4 AU. The starting belt therefore has no gaps at all. That uniform
  starting distribution is the one invented ingredient: it is the null
  hypothesis the experiment tests.
* Particles are massless, so they feel the planets but not each other.

Removal
-------
A particle is *removed* once its perihelion reaches inside Mars' aphelion
(1.666 AU) or its aphelion reaches Jupiter's perihelion (4.95 AU). Mars itself
is not in the simulation; from a Mars-crossing orbit, close encounters with
Mars and Earth (and ultimately collision with the Sun or a planet, or
ejection) clear an asteroid out on million-year timescales (Gladman et al.
1997). Counting the crossing as removal is the standard shortcut.

Integrator: WHFast (Wisdom-Holman symplectic, Rein & Tamayo 2015) in
democratic heliocentric coordinates, dt = 0.1 yr, via REBOUND.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import mpcorb
from config import (
    GAUSS_K,
    HORIZONS_PATH,
    JUPITER_PERIHELION_AU,
    MASS_RATIO,
    MAX_UNCERTAINTY_U,
    SIM_A_RANGE,
    SIM_DT_YEARS,
    SIM_N_PARTICLES,
    SIM_PATH,
    SIM_REMOVE_Q_AU,
    SIM_SEED,
    SIM_SNAPSHOTS,
    SIM_YEARS,
)

#: G*M_sun in AU^3/yr^2 with the Julian year: (k * 365.25)^2 = 39.4769. Not
#: 4 pi^2 = 39.4784, which would make the "year" 1 AU's own orbital period.
#: REBOUND's ('yr', 'AU', 'Msun') units use this value too (checked in _build).
GM_SUN = (GAUSS_K * 365.25) ** 2

CAUSE_NONE, CAUSE_MARS, CAUSE_JUPITER, CAUSE_UNBOUND = 0, 1, 2, 3


def xyz_to_aei(pos, vel, mu: float = GM_SUN):
    """Heliocentric osculating a, e, i (deg) from positions and velocities.

    Vectorised over rows. Vis-viva gives a; the eccentricity vector gives e;
    the angular momentum vector gives i. Unbound orbits get a < 0, e >= 1.
    """
    pos = np.atleast_2d(pos)
    vel = np.atleast_2d(vel)
    r = np.linalg.norm(pos, axis=1)
    v2 = np.einsum("ij,ij->i", vel, vel)
    a = 1.0 / (2.0 / r - v2 / mu)
    h = np.cross(pos, vel)
    hn = np.linalg.norm(h, axis=1)
    evec = np.cross(vel, h) / mu - pos / r[:, None]
    e = np.linalg.norm(evec, axis=1)
    inc = np.degrees(np.arccos(np.clip(h[:, 2] / hn, -1, 1)))
    return a, e, inc


def draw_initial_conditions(n: int, seed: int = SIM_SEED) -> dict:
    """Uniform semi-major axes; everything else borrowed from real asteroids.

    Donors are well-determined main-belt orbits (U <= 2, 2.1-3.3 AU). A donor
    is rejected if its eccentricity would make the particle planet-crossing
    at its new semi-major axis from the start.
    """
    d = mpcorb.load_mpcorb(quiet=True)
    donor_ok = ((d["U"] >= 0) & (d["U"] <= MAX_UNCERTAINTY_U)
                & (d["a"] > 2.1) & (d["a"] < 3.3))
    donors = np.flatnonzero(donor_ok)
    rng = np.random.default_rng(seed)
    a = np.sort(rng.uniform(*SIM_A_RANGE, n))
    picked = np.empty(n, dtype=np.int64)
    for k in range(n):
        while True:
            j = donors[rng.integers(len(donors))]
            e = d["e"][j]
            if a[k] * (1 - e) > SIM_REMOVE_Q_AU and a[k] * (1 + e) < JUPITER_PERIHELION_AU:
                picked[k] = j
                break
    return {
        "a": a,
        "e": d["e"][picked],
        "inc": d["incl"][picked],
        "node": d["node"][picked],
        "peri": d["peri"][picked],
        "M": d["M"][picked],
        "donor_number": d["number"][picked],
        "donor_readable": d["readable"][picked],
    }


def _build(ic: dict, idx: np.ndarray):
    import rebound

    planets = json.loads(HORIZONS_PATH.read_text())["bodies"]
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    assert abs(sim.G / GM_SUN - 1) < 1e-9, (sim.G, GM_SUN)
    sim.add(m=1.0, name="sun")
    for name in ("Jupiter", "Saturn"):
        b = planets[name]
        x, y, z = b["xyz_au"]
        vx, vy, vz = (np.array(b["vxyz_au_per_day"]) * 365.25).tolist()
        sim.add(m=MASS_RATIO[name], x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, name=name)
    sim.N_active = 3
    sim.testparticle_type = 0
    for k in idx:
        # Re-fetch the Sun each time: adding particles can reallocate the
        # particle array, which would leave a stored reference dangling.
        sim.add(primary=sim.particles[0], a=ic["a"][k], e=ic["e"][k], inc=np.radians(ic["inc"][k]),
                Omega=np.radians(ic["node"][k]), omega=np.radians(ic["peri"][k]),
                M=np.radians(ic["M"][k]), name=str(int(k)))
    sim.move_to_com()
    sim.integrator = "whfast"
    sim.integrator.coordinates = "democraticheliocentric"
    sim.integrator.safe_mode = 0
    sim.dt = SIM_DT_YEARS
    return sim


def _ids(sim) -> np.ndarray:
    """Test-particle ids (their names), in the simulation's current order."""
    return np.array([int(sim.particles[j].name) for j in range(3, sim.N)], dtype=np.int64)


def _state(sim):
    """Heliocentric positions and velocities of the test particles."""
    sim.synchronize()
    n = sim.N
    xyz = np.zeros((n, 3))
    vxyz = np.zeros((n, 3))
    sim.serialize_particle_data(xyz=xyz, vxvyvz=vxyz)
    return xyz[3:] - xyz[0], vxyz[3:] - vxyz[0]


def run_chunk(args) -> dict:
    """Integrate one chunk of test particles; returns sampled a, e, i."""
    ic, idx, times = args
    sim = _build(ic, idx)
    n = len(idx)
    pos_of = {int(k): j for j, k in enumerate(idx)}
    a_out = np.full((len(times), n), np.nan, dtype=np.float32)
    e_out = np.full((len(times), n), np.nan, dtype=np.float32)
    i_out = np.full((len(times), n), np.nan, dtype=np.float32)
    t_removed = np.full(n, np.nan)
    cause = np.zeros(n, dtype=np.int8)
    e0 = sim.energy()
    t_start = time.time()
    ids = _ids(sim)
    for it, t in enumerate(times):
        if len(ids) == 0:
            break
        if t > 0:
            sim.integrate(t, exact_finish_time=0)
        x, v = _state(sim)
        a, e, inc = xyz_to_aei(x, v)
        cols = np.array([pos_of[int(k)] for k in ids])
        a_out[it, cols] = a
        e_out[it, cols] = e
        i_out[it, cols] = inc
        q = a * (1 - e)
        Q = a * (1 + e)
        unbound = (a <= 0) | (e >= 1)
        gone = unbound | (q < SIM_REMOVE_Q_AU) | (Q > JUPITER_PERIHELION_AU)
        for j in np.flatnonzero(gone):
            c = cols[j]
            t_removed[c] = sim.t
            cause[c] = (CAUSE_UNBOUND if unbound[j] else
                        CAUSE_MARS if q[j] < SIM_REMOVE_Q_AU else CAUSE_JUPITER)
            sim.remove(str(int(ids[j])))
        if gone.any():
            ids = _ids(sim)
    energy_drift = abs((sim.energy() - e0) / e0)
    return {"idx": idx, "a": a_out, "e": e_out, "inc": i_out, "t_removed": t_removed,
            "cause": cause, "energy_drift": energy_drift,
            "wall_seconds": time.time() - t_start}


def run(n_particles: int = SIM_N_PARTICLES, years: float = SIM_YEARS,
        n_snapshots: int = SIM_SNAPSHOTS, workers: int = 6, quiet: bool = False) -> dict:
    ic = draw_initial_conditions(n_particles)
    times = np.linspace(0.0, years, n_snapshots)
    # Interleave particles across workers so each gets the full a range.
    chunks = [np.arange(w, n_particles, workers) for w in range(workers)]
    if not quiet:
        print(f"  integrating {n_particles} test particles + Jupiter + Saturn for "
              f"{years / 1e6:g} Myr on {workers} processes")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(run_chunk, [(ic, c, times) for c in chunks]))
    out = {
        "t_years": times,
        "a": np.full((len(times), n_particles), np.nan, dtype=np.float32),
        "e": np.full((len(times), n_particles), np.nan, dtype=np.float32),
        "inc": np.full((len(times), n_particles), np.nan, dtype=np.float32),
        "t_removed": np.full(n_particles, np.nan),
        "cause": np.zeros(n_particles, dtype=np.int8),
    }
    for r in results:
        for k in ("a", "e", "inc"):
            out[k][:, r["idx"]] = r[k]
        out["t_removed"][r["idx"]] = r["t_removed"]
        out["cause"][r["idx"]] = r["cause"]
    for k, v in ic.items():
        out["ic_" + k] = v
    out["energy_drift_max"] = np.array(max(r["energy_drift"] for r in results))
    out["wall_seconds"] = np.array(time.time() - t0)
    out["dt_years"] = np.array(SIM_DT_YEARS)
    out["seed"] = np.array(SIM_SEED)
    np.savez(SIM_PATH, **out)
    if not quiet:
        removed = np.isfinite(out["t_removed"]).sum()
        print(f"  done in {out['wall_seconds'] / 60:.1f} min: {removed} of {n_particles} "
              f"removed; planets' energy drift {float(out['energy_drift_max']):.1e}")
    return out


def load_simulation() -> dict | None:
    if not SIM_PATH.exists():
        return None
    with np.load(SIM_PATH) as z:
        return {k: z[k] for k in z.files}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Run the (labelled) N-body simulation.")
    parser.add_argument("--particles", type=int, default=SIM_N_PARTICLES)
    parser.add_argument("--years", type=float, default=SIM_YEARS)
    parser.add_argument("--snapshots", type=int, default=SIM_SNAPSHOTS)
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args(argv)
    run(args.particles, args.years, args.snapshots, args.workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
