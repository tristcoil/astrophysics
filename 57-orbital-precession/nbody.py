"""Our own integrations of the Solar System with REBOUND and REBOUNDx.

**These are simulations** (atlas ground rule 1): every one starts from DE440's
positions and velocities at J2000.0 and DE440's masses, and is then left to
evolve under the forces we choose. Comparing runs with and without a force is
how each force's share of Mercury's precession is measured.

Units: G = 1 with masses given as GM in au^3/day^2, so lengths are in au and
times in days, exactly as in DE440. The frame is DE440's (the ICRF, origin at
the barycentre). The integrator is IAS15, an adaptive 15th-order scheme whose
errors stay at the level of floating-point round-off (Rein & Spiegel 2015);
WHFast, a fixed-step symplectic scheme, is used only to show what a too-long
step does.

Each run integrates from J2000 back to 1550 and forward to 2650 and records
Mercury's heliocentric state every 2 days and every other planet's every 20
days. The records are cached in ``data/derived/runs/<name>.npz``.
"""

from __future__ import annotations

import json
import time

import numpy as np

import ephemeris as eph
from config import (J2000_JD, RUN, RUNS_DIR, SAMPLE_DAYS_OUTER, SUN_POLE_DEC_DEG,
                    SUN_POLE_RA_DEG, YEAR_DAYS, RunSpec)


def sun_pole() -> np.ndarray:
    return eph.unit_vector_radec(SUN_POLE_RA_DEG, SUN_POLE_DEC_DEG)


def build(spec: RunSpec, *, jd0: float = J2000_JD, gr_boost: float | None = None):
    """A REBOUND simulation for ``spec``, started from DE440 at ``jd0``.

    Returns (sim, rebx); keep a reference to ``rebx`` for as long as the
    simulation runs (REBOUNDx's forces live in it).
    """
    import rebound
    import reboundx

    sim = rebound.Simulation()
    sim.G = 1.0
    scale = dict(spec.mass_scale)
    for key in spec.bodies:
        r, v = eph.barycentric_state(key, jd0)
        sim.add(m=eph.gm(key) * scale.get(key, 1.0), x=r[0], y=r[1], z=r[2], vx=v[0], vy=v[1],
                vz=v[2], name=key)
    sim.move_to_com()  # DE440's barycentre includes asteroids we leave out
    if spec.integrator == "ias15":
        sim.integrator = "ias15"
    elif spec.integrator == "whfast":
        sim.integrator = "whfast"
        sim.dt = spec.dt_days
    else:
        raise ValueError(spec.integrator)

    rebx = None
    if spec.effects:
        rebx = reboundx.Extras(sim)
        boost = spec.gr_boost if gr_boost is None else gr_boost
        c = eph.c_au_per_day() / np.sqrt(boost)
        sun = sim.particles["sun"]
        phys = eph.sun_physical()
        spin = sun_pole() * phys["omega_rad_per_day"]
        for eff in spec.effects:
            if eff in ("gr", "gr_full", "gr_potential"):
                force = rebx.load_force(eff)
                rebx.add_force(force)
                force.params["c"] = c
            elif eff == "j2":
                force = rebx.load_force("gravitational_harmonics")
                rebx.add_force(force)
                sun.params["J2"] = phys["j2"]
                sun.params["R_eq"] = phys["radius_au"]
                sun.params["Omega"] = rebound.Vec3d(*spin)
            elif eff == "lt":
                force = rebx.load_force("lense_thirring")
                rebx.add_force(force)
                force.params["lt_c"] = c
                sun.params["I"] = phys["moment_of_inertia_factor"] * sun.m * phys["radius_au"] ** 2
                sun.params["Omega"] = rebound.Vec3d(*spin)
            else:
                raise ValueError(f"unknown effect {eff!r}")
    return sim, rebx


def _sample_times(spec: RunSpec, direction: int) -> tuple[np.ndarray, np.ndarray]:
    """Output times (days from J2000) for one direction, and which are 'outer' samples."""
    span = (spec.forward_yr if direction > 0 else spec.back_yr) * YEAR_DAYS
    step = spec.sample_days if spec.integrator == "ias15" else max(spec.sample_days, spec.dt_days)
    t = np.arange(int(np.floor(span / step + 1e-9)) + 1) * step
    if direction < 0:
        t = -t
    every = spec.outer_every if spec.integrator == "ias15" else max(
        1, int(round(SAMPLE_DAYS_OUTER / step)))
    outer = np.zeros(len(t), bool)
    outer[::max(every, 1)] = True
    return t, outer


def _integrate_direction(spec: RunSpec, direction: int) -> dict:
    sim, rebx = build(spec)
    keys = list(spec.bodies)
    i_sun, i_me = keys.index("sun"), keys.index("mercury")
    t_out, outer = _sample_times(spec, direction)
    n, n_outer = len(t_out), int(outer.sum())
    buf = np.zeros((len(keys), 6))
    t_actual = np.empty(n)
    merc = np.empty((n, 6))
    others = np.empty((n_outer, len(keys), 6))
    t_outer = np.empty(n_outer)
    exact = 1 if spec.integrator == "ias15" else 0
    j = 0
    e0 = sim.energy() if not spec.effects else None
    wall = time.time()
    for i, t in enumerate(t_out):
        if t != 0.0:
            sim.integrate(t, exact_finish_time=exact)
        sim.serialize_particle_data(xyzvxvyvz=buf)
        t_actual[i] = sim.t
        merc[i] = buf[i_me] - buf[i_sun]
        if outer[i]:
            others[j] = buf - buf[i_sun]
            t_outer[j] = sim.t
            j += 1
    wall = time.time() - wall
    info = {"steps": int(sim.steps_done), "wall_s": wall}
    if e0 is not None:
        info["rel_energy_error"] = float(abs((sim.energy() - e0) / e0))
    return {"t": t_actual, "mercury": merc, "t_outer": t_outer, "outer": others, "info": info,
            "keys": keys}


def run(spec: RunSpec, *, force: bool = False) -> dict:
    """Integrate ``spec`` both ways from J2000 (cached). Returns the merged record.

    ``t`` is in days from J2000 (TDB); ``mercury`` is (N, 6) heliocentric
    position+velocity; ``outer`` is (M, n_bodies, 6) heliocentric states of all
    bodies at the coarser cadence (the Sun's row is zero).
    """
    path = RUNS_DIR / f"{spec.name}.npz"
    if path.exists() and not force:
        return load(spec.name)
    parts = []
    for direction in (-1, +1):
        span = spec.back_yr if direction < 0 else spec.forward_yr
        if span > 0:
            parts.append((direction, _integrate_direction(spec, direction)))
    # Merge: backward part reversed (dropping its t=0 duplicate), then forward.
    t, merc, t_o, outer, info = [], [], [], [], {}
    for direction, p in parts:
        sl = slice(None, None, -1) if direction < 0 else slice(None)
        tt, mm, to, oo = p["t"][sl], p["mercury"][sl], p["t_outer"][sl], p["outer"][sl]
        if direction > 0 and t:
            tt, mm, to, oo = tt[1:], mm[1:], to[1:], oo[1:]
        t.append(tt); merc.append(mm); t_o.append(to); outer.append(oo)
        info["backward" if direction < 0 else "forward"] = p["info"]
    rec = {"t": np.concatenate(t), "mercury": np.concatenate(merc),
           "t_outer": np.concatenate(t_o), "outer": np.concatenate(outer)}
    meta = {"name": spec.name, "bodies": list(spec.bodies), "effects": list(spec.effects),
            "integrator": spec.integrator, "dt_days": spec.dt_days, "gr_boost": spec.gr_boost,
            "description": spec.description, "info": info, "mass_scale": list(spec.mass_scale),
            "c_au_per_day": eph.c_au_per_day(), "epoch_jd": J2000_JD}
    np.savez_compressed(path, meta=json.dumps(meta), **rec)
    return load(spec.name)


def load(name: str) -> dict:
    path = RUNS_DIR / f"{name}.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} not computed yet -- run `python run.py`")
    with np.load(path) as z:
        out = {k: z[k] for k in z.files if k != "meta"}
        out["meta"] = json.loads(str(z["meta"]))
    return out


def run_named(name: str, force: bool = False) -> str:
    """Pool-friendly wrapper: run one spec by name, return a one-line summary."""
    import threadpoolctl

    with threadpoolctl.threadpool_limits(1):
        t0 = time.time()
        rec = run(RUN[name], force=force)
    info = rec["meta"]["info"]
    steps = sum(v["steps"] for v in info.values())
    return f"{name}: {len(rec['t'])} samples, {steps} steps, {time.time() - t0:.0f} s"


# --------------------------------------------------------------------------
# Numerical checks
# --------------------------------------------------------------------------


def reversibility_test(years: float, spec_name: str = "all_newton") -> dict:
    """Integrate forward ``years``, then back to J2000; how far is Mercury off?

    IAS15 is not exactly time-reversible, so the residual measures the
    accumulated integration error of the round trip (twice the run length).
    """
    from orbit_physics import elements
    from config import ARCSEC_PER_RAD

    spec = RUN[spec_name]
    sim, rebx = build(spec)
    keys = list(spec.bodies)
    buf0 = np.zeros((len(keys), 6))
    sim.serialize_particle_data(xyzvxvyvz=buf0)
    t0 = time.time()
    sim.integrate(years * YEAR_DAYS)
    sim.integrate(0.0)
    wall = time.time() - t0
    buf1 = np.zeros_like(buf0)
    sim.serialize_particle_data(xyzvxvyvz=buf1)
    i_s, i_m = keys.index("sun"), keys.index("mercury")
    mu = eph.gm("sun") + eph.gm("mercury")
    s0, s1 = buf0[i_m] - buf0[i_s], buf1[i_m] - buf1[i_s]
    e0 = elements(s0[:3], s0[3:], mu)
    e1 = elements(s1[:3], s1[3:], mu)
    ang = np.arctan2(np.linalg.norm(np.cross(e0["ehat"][:, 0], e1["ehat"][:, 0])),
                     float(e0["ehat"][:, 0] @ e1["ehat"][:, 0]))
    return {"years_each_way": years, "position_error_km": float(np.linalg.norm(s1[:3] - s0[:3]) * eph.au_km()),
            "perihelion_direction_error_arcsec": float(ang * ARCSEC_PER_RAD),
            "time_offset_days": float(sim.t), "wall_s": wall}


def run_all(names: list[str], workers: int = 4, force: bool = False) -> None:
    """Run several specs in parallel worker processes (each single-threaded)."""
    import multiprocessing as mp

    todo = [n for n in names if force or not (RUNS_DIR / f"{n}.npz").exists()]
    if not todo:
        return
    # Longest first, so the pool finishes together.
    cost = {"gr_full": 10.0, "gr": 2.0, "gr_potential": 1.5, "j2": 1.5, "lt": 1.5}
    def weight(n):
        s = RUN[n]
        w = len(s.bodies) ** 2 * (s.back_yr + s.forward_yr)
        for e in s.effects:
            w *= cost.get(e, 1.0)
        return -w
    todo.sort(key=weight)
    ctx = mp.get_context("fork")
    with ctx.Pool(min(workers, len(todo)), maxtasksperchild=1) as pool:
        for line in pool.imap_unordered(_run_named_force if force else run_named, todo):
            print("  " + line, flush=True)


def _run_named_force(name: str) -> str:
    return run_named(name, force=True)
