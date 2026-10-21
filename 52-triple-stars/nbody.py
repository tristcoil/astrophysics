"""REBOUND helpers shared by every simulation in this project (all labelled as
simulations; their starting conditions come from published or fitted orbits).

Frames and conventions
----------------------
The sky frame has z toward the observer and x, y in the plane of the sky.
REBOUND then measures inclinations from the sky plane, and its ascending node
is where a body crosses the sky plane moving *toward* us. Astronomers' radial
velocities count positive *away* from us, so an argument of periastron from
spectroscopy relates to REBOUND's by omega_rebound = omega_RV - 180 deg for the
same body. For a single-lined outer orbit, omega_RV belongs to the inner
binary's barycentre; the third star's omega is 180 deg away, so the third
star's REBOUND omega equals the spectroscopic omega of the binary.

Units: au, days, solar masses (G = C.G_AU_DAY) unless stated.
"""

from __future__ import annotations

import numpy as np
import rebound

import config as C


def new_sim(units="day") -> rebound.Simulation:
    sim = rebound.Simulation()
    sim.G = C.G_AU_DAY if units == "day" else C.G_AU_YR
    return sim


def unit_normal(inc, Omega):
    """Orbit normal (unit angular momentum) for inclination and node in the
    sky frame (radians)."""
    return np.array([np.sin(inc) * np.sin(Omega), -np.sin(inc) * np.cos(Omega), np.cos(inc)])


def inc_node_from_normal(n):
    n = np.asarray(n, float) / np.linalg.norm(n)
    inc = np.arccos(np.clip(n[2], -1, 1))
    Omega = np.arctan2(n[0], -n[1])
    return inc, Omega


def add_inner_from_state(sim, m1, m2, r_rel, v_rel):
    """Primary at the origin, secondary at relative position/velocity (au, au/day)."""
    sim.add(m=m1)
    sim.add(m=m2, x=r_rel[0], y=r_rel[1], z=r_rel[2], vx=v_rel[0], vy=v_rel[1], vz=v_rel[2])


def add_third(sim, m3, P, e, inc, Omega, omega, M):
    """Third body on a Keplerian orbit about the inner pair's barycentre
    (Jacobi coordinates: REBOUND's default primary is the centre of mass of
    all particles added so far)."""
    sim.add(m=m3, P=P, e=e, inc=inc, Omega=Omega, omega=omega, M=M)


def relative_state(sim, i=1, j=0):
    p, q = sim.particles[i], sim.particles[j]
    return (np.array([p.x - q.x, p.y - q.y, p.z - q.z]),
            np.array([p.vx - q.vx, p.vy - q.vy, p.vz - q.vz]))


def inner_normal(sim):
    r, v = relative_state(sim, 1, 0)
    h = np.cross(r, v)
    return h / np.linalg.norm(h)


def outer_normal(sim):
    ps = sim.particles
    m12 = ps[0].m + ps[1].m
    rc = (ps[0].m * np.array([ps[0].x, ps[0].y, ps[0].z]) + ps[1].m * np.array([ps[1].x, ps[1].y, ps[1].z])) / m12
    vc = (ps[0].m * np.array([ps[0].vx, ps[0].vy, ps[0].vz]) + ps[1].m * np.array([ps[1].vx, ps[1].vy, ps[1].vz])) / m12
    r = np.array([ps[2].x, ps[2].y, ps[2].z]) - rc
    v = np.array([ps[2].vx, ps[2].vy, ps[2].vz]) - vc
    h = np.cross(r, v)
    return h / np.linalg.norm(h)


def inner_elements(sim):
    """Osculating inner-orbit elements (a, e, inc, Omega, omega) in the sky frame."""
    o = sim.particles[1].orbit(primary=sim.particles[0])
    return o.a, o.e, o.inc, o.Omega, o.omega


def setup_whfast(sim, dt):
    sim.integrator = "whfast"
    sim.dt = dt
    sim.integrator.safe_mode = 0
    sim.integrator.keep_unsynchronized = 1
    sim.move_to_com()


def record(sim, times, funcs: dict) -> dict:
    """Integrate through `times` (increasing), evaluating funcs[name](sim)."""
    out = {k: [] for k in funcs}
    for t in times:
        sim.integrate(t, exact_finish_time=0)
        for k, f in funcs.items():
            out[k].append(f(sim))
    return {k: np.array(v) for k, v in out.items()}
