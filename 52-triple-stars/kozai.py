"""Kozai-Lidov cycles and what stops them: Algol as the test case.

Algol's inner pair (A, a B8 dwarf, and B, a K subgiant filling its Roche
lobe; 2.87 d) and its third star C (680 d) orbit in planes 90.2 deg apart
(Baron et al. 2012). For point masses that is the worst possible geometry:
von Zeipel-Kozai-Lidov cycles would pump the inner eccentricity toward 1 in
a few centuries and the two stars would collide. Yet Algol has eclipsed
steadily since antiquity.

Labelled simulations, all started from Baron et al.'s orbits:
  1. Newtonian point masses;
  2. plus general relativity (REBOUNDx 'gr_full');
  3. plus tidal and rotational bulges of A and B with synchronised spins
     (REBOUNDx 'tides_spin', conservative: no tidal friction).
Each is followed until the stars touch (periastron < R_A + R_B) or for
3,000 years.

Physical inputs for the bulges (assumptions, varied in the README):
apsidal-motion constants k_2 of 0.005 (A, a radiative B-type star; Claret &
Gimenez tables give log k2 ~ -2.3) and 0.05 (B, a mostly convective,
Roche-lobe-filling subgiant), REBOUNDx's Love number k2 = 2 x k_2; squared
gyration radii 0.06 and 0.2.
"""

from __future__ import annotations

import json
import time

import numpy as np
import rebound
import reboundx

import config as C
import nbody as NB

SYSA = C.SYSTEMS["Algol"]
DEG = np.pi / 180
OUT = C.CACHE / "kozai"
OUT.mkdir(parents=True, exist_ok=True)

K2_APSIDAL = {"A": 0.005, "B": 0.05}
GYRATION2 = {"A": 0.06, "B": 0.2}


def vb_to_rebound(i_vb, Omega_vb, omega_vb=None):
    """Visual-binary angles (node = position angle where the companion
    recedes, measured north through east; i < 90 for counter-clockwise
    motion) -> REBOUND angles in the frame (x north, y east, z toward us)."""
    inc = i_vb * DEG
    Om = (Omega_vb + 180.0) * DEG
    om = None if omega_vb is None else (omega_vb + 180.0) * DEG
    return inc, Om, om


def build_algol(t0_jd=2455000.0):
    """Three-body Algol from Baron et al. (2012). Inner orbit circular; the
    phases are set from their epochs (primary minimum for the inner pair,
    periastron for the outer orbit)."""
    a = SYSA
    sim = NB.new_sim("day")
    sim.t = t0_jd
    sim.add(m=a["m_A"], r=a["R_A"] * C.R_SUN_AU)
    inc, Om, _ = vb_to_rebound(a["inner"]["i"], a["inner"]["Omega"])
    # primary eclipse = B in front of A: B on the near side, closest to the
    # line of sight -> true longitude 90 deg from the REBOUND node (z max)
    n_cyc = (t0_jd - a["inner"]["T0"]) / a["inner"]["P"]
    theta = 0.5 * np.pi + 2 * np.pi * (n_cyc % 1.0)
    sim.add(m=a["m_B"], r=a["R_B"] * C.R_SUN_AU, P=a["inner"]["P"], e=0.0, inc=inc, Omega=Om, theta=theta)
    inc, Om, om = vb_to_rebound(a["outer"]["i"], a["outer"]["Omega"], a["outer"]["omega"])
    M = 2 * np.pi * ((t0_jd - a["outer"]["T0"]) / a["outer"]["P"] % 1.0)
    sim.add(m=a["m_C"], r=a["R_C"] * C.R_SUN_AU, P=a["outer"]["P"], e=a["outer"]["e"], inc=inc,
            Omega=Om, omega=om, M=M)
    sim.move_to_com()
    return sim


def mutual_inclination(sim):
    return np.degrees(np.arccos(np.clip(np.dot(NB.inner_normal(sim), NB.outer_normal(sim)), -1, 1)))


#: variants of the bulge physics: which star's tidal bulge is kept, and how soft
VARIANTS = {
    "tides": {"A": K2_APSIDAL["A"], "B": K2_APSIDAL["B"]},
    "tides_Bsoft": {"A": K2_APSIDAL["A"], "B": K2_APSIDAL["B"] / 10},
    "tides_Aonly": {"A": K2_APSIDAL["A"], "B": 0.0},
}


def add_physics(sim, gr=False, tides=False, k2=None):
    rebx = reboundx.Extras(sim)
    if gr:
        f = rebx.load_force("gr_full")
        rebx.add_force(f)
        f.params["c"] = C.C_AU_DAY
    if tides:
        f = rebx.load_force("tides_spin")
        rebx.add_force(f)
        ps = sim.particles
        o = ps[1].orbit(primary=ps[0])
        n_orb = 2 * np.pi / o.P
        hn = NB.inner_normal(sim)
        for k, name in ((0, "A"), (1, "B")):
            p = ps[k]
            p.params["k2"] = 2 * (k2 or K2_APSIDAL)[name]
            p.params["tau"] = 0.0
            p.params["I"] = GYRATION2[name] * p.m * p.r**2
            p.params["Omega"] = rebound.Vec3d(*(n_orb * hn))
        rebx.initialize_spin_ode(f)
    return rebx


def run(label, gr=False, tides=False, years=3000.0, record_every_days=50.0, k2=None):
    path = OUT / f"algol_{label}.json"
    if path.exists():
        return json.loads(path.read_text())
    sim = build_algol()
    sim.integrator = "ias15"
    rebx = add_physics(sim, gr=gr, tides=tides, k2=k2)
    sim.collision = "none"
    a = SYSA
    r_touch = (a["R_A"] + a["R_B"]) * C.R_SUN_AU
    t_start = sim.t
    ts = t_start + np.arange(0, years * 365.25, record_every_days)
    out = dict(t_yr=[], e_in=[], a_in=[], J=[], i_sky=[])
    touched = None
    t0 = time.time()
    for t in ts:
        sim.integrate(t)
        o = sim.particles[1].orbit(primary=sim.particles[0])
        out["t_yr"].append((sim.t - t_start) / 365.25)
        out["e_in"].append(o.e)
        out["a_in"].append(o.a)
        out["J"].append(mutual_inclination(sim))
        out["i_sky"].append(np.degrees(o.inc))
        if o.a * (1 - o.e) < r_touch:
            touched = (sim.t - t_start) / 365.25
            break
    out["touched_yr"] = touched
    out["label"] = label
    out["wall_s"] = time.time() - t0
    out["e_max"] = float(np.max(out["e_in"]))
    path.write_text(json.dumps(out))
    return out


def timescales(k2=None):
    """The competing clocks (years), from the standard formulae; k2: apsidal
    constants of A and B (default K2_APSIDAL; see VARIANTS)."""
    a = SYSA
    M_in = a["m_A"] + a["m_B"]
    M = M_in + a["m_C"]
    P_in, P_out, e_out = a["inner"]["P"], a["outer"]["P"], a["outer"]["e"]
    a_in = C.semimajor_axis_au(M_in, P_in)
    # Kozai-Lidov (quadrupole) timescale, e.g. Holman et al. 1997 / Antognini 2015
    t_kl = (2 / (3 * np.pi)) * (M / a["m_C"]) * P_out**2 / P_in * (1 - e_out**2) ** 1.5 / 365.25
    # GR periastron advance of a circular inner orbit
    dw_gr = 6 * np.pi * C.GM_SUN * M_in / (C.C_LIGHT**2 * a_in * C.AU)      # rad per orbit
    t_gr = 2 * np.pi / dw_gr * P_in / 365.25
    # tidal + rotational bulges (Sterne 1939; synchronous rotation), e = 0
    RA, RB = a["R_A"] * C.R_SUN_AU / a_in, a["R_B"] * C.R_SUN_AU / a_in
    kA, kB = (k2 or K2_APSIDAL)["A"], (k2 or K2_APSIDAL)["B"]
    rate = (15 * (kA * a["m_B"] / a["m_A"] * RA**5 + kB * a["m_A"] / a["m_B"] * RB**5)
            + kA * RA**5 * (1 + a["m_B"] / a["m_A"]) + kB * RB**5 * (1 + a["m_A"] / a["m_B"]))
    t_tide = 1.0 / rate * P_in / 365.25            # rate is (domega/dt)/n
    return dict(t_KL=float(t_kl), t_GR=float(t_gr), t_tide=float(t_tide), a_in_au=float(a_in),
                rA=float(RA), rB=float(RB))


if __name__ == "__main__":
    import sys
    print(json.dumps(timescales(), indent=1))
    which = sys.argv[1:] or ["newton", "gr", "tides"]
    for lab in which:
        r = run(lab, gr=lab != "newton", tides=lab.startswith("tides"), k2=VARIANTS.get(lab))
        print(lab, "touched at", r["touched_yr"], "e_max", round(r["e_max"], 4), f"{r['wall_s']:.0f}s")
