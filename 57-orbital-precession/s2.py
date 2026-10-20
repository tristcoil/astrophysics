"""S2 around Sgr A*: the strongest-gravity rung of the precession ladder.

S2 swings to within 120 au of the Galactic Centre's black hole every 16 years
at 2.7% of the speed of light. General relativity turns its orbit by 12.1
arcminutes per revolution. GRAVITY detected that in 2020 (f_SP = 1.10 +- 0.19,
where 1 is GR and 0 is Newton).

GRAVITY's astrometry is not public, and H1 (17-galactic-centre) showed that
the public Keck and NTT/VLT positions cannot see the precession. So here we
compute what the signal *looks like*: starting from GRAVITY's 2020 orbit
(osculating at the 2010.35 apocentre), we integrate S2 once with Newton's
gravity and once with first-post-Newtonian gravity (REBOUNDx ``gr``) and
project both onto the sky. The difference is the precession signal, and it
can be set against each instrument's measurement errors. GRAVITY 2020 quotes
the result we should get: a 0.78 mas shift at apocentre, mostly in RA.

Sky convention (visual binaries, as GRAVITY and H1): X = offset towards East
(RA), Y = towards North (Dec); Omega is the position angle of the ascending
node, i > 90 deg means clockwise motion on the sky.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import ephemeris as eph
from config import PROJECT_DIR, PUBLISHED, YEAR_DAYS

H1_DIR = PROJECT_DIR.parent / "17-galactic-centre"


def _rotation(inc, Omega, omega):
    """Columns: unit vectors to periapsis (P) and 90 deg ahead (Q), in (X=E, Y=N, Z)."""
    cO, sO, ci, si, cw, sw = (np.cos(Omega), np.sin(Omega), np.cos(inc), np.sin(inc),
                              np.cos(omega), np.sin(omega))
    # Thiele-Innes: Dec (N) = A x + F y, RA (E) = B x + G y.
    A = cw * cO - sw * sO * ci
    B = cw * sO + sw * cO * ci
    F = -sw * cO - cw * sO * ci
    G = -sw * sO + cw * cO * ci
    C = sw * si
    H = cw * si
    P = np.array([B, A, C])  # (E, N, line of sight)
    Q = np.array([G, F, H])
    return P, Q


def initial_state(orbit: dict) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Position (au) and velocity (au/day) of S2 at the 2010.35 apocentre epoch."""
    gm = orbit["mass_1e6"] * 1e6 * eph.gm("sun")  # au^3/day^2
    a = orbit["a_mas"] / 1000.0 * orbit["r0_pc"]  # au (1 arcsec at 1 pc = 1 au)
    e = orbit["e"]
    t_ref = 2010.35
    period_days = 2 * np.pi * np.sqrt(a**3 / gm)
    M = 2 * np.pi * ((t_ref - orbit["t_peri_yr"]) * YEAR_DAYS / period_days)
    E = M
    for _ in range(60):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    n = 2 * np.pi / period_days
    x, y = a * (np.cos(E) - e), a * np.sqrt(1 - e**2) * np.sin(E)
    edot = n / (1 - e * np.cos(E))
    vx, vy = -a * np.sin(E) * edot, a * np.sqrt(1 - e**2) * np.cos(E) * edot
    P, Q = _rotation(np.radians(orbit["inc_deg"]), np.radians(orbit["Omega_deg"]),
                     np.radians(orbit["omega_deg"]))
    return P * x + Q * y, P * vx + Q * vy, gm, t_ref


def _simulation(r0, v0, gm, gr: bool, boost: float):
    import rebound
    import reboundx

    sim = rebound.Simulation()
    sim.G = 1.0
    sim.add(m=gm)
    sim.add(m=0.0, x=r0[0], y=r0[1], z=r0[2], vx=v0[0], vy=v0[1], vz=v0[2])
    sim.integrator = "ias15"
    rebx = None
    if gr:
        rebx = reboundx.Extras(sim)
        f = rebx.load_force("gr")
        rebx.add_force(f)
        f.params["c"] = eph.c_au_per_day() / np.sqrt(boost)
    return sim, rebx


def integrate(orbit: dict, years: np.ndarray, *, gr: bool, boost: float = 1.0) -> np.ndarray:
    """Sky positions (mas, X east / Y north) and line-of-sight z (au) at decimal years.

    Integrates outwards from the 2010.35 reference epoch in both directions, with
    first-post-Newtonian gravity (``gr``, optionally ``boost``-ed) or Newton's.
    """
    r0, v0, gm, t_ref = initial_state(orbit)
    t_days = (np.asarray(years) - t_ref) * YEAR_DAYS
    out = np.empty((len(t_days), 3))
    for backwards in (True, False):
        idx = np.flatnonzero((t_days < 0) == backwards)
        if not len(idx):
            continue
        sim, rebx = _simulation(r0, v0, gm, gr, boost)  # keep rebx alive while integrating
        for i in idx[np.argsort(np.abs(t_days[idx]))]:
            sim.integrate(t_days[i])
            p = sim.particles[1]
            out[i] = (p.x, p.y, p.z)
    au_to_mas = 1000.0 / orbit["r0_pc"]
    return np.column_stack([out[:, 0] * au_to_mas, out[:, 1] * au_to_mas, out[:, 2]])


def precession_signal(orbit: dict | None = None) -> dict:
    orbit = orbit or PUBLISHED["s2_orbit_gravity2020"]
    years = np.arange(1992.0, 2040.0, 2.0 / 365.25)
    newton = integrate(orbit, years, gr=False)
    einstein = integrate(orbit, years, gr=True)
    d = einstein[:, :2] - newton[:, :2]
    sep = np.hypot(d[:, 0], d[:, 1])
    # Apocentre after the 2018 pericentre: largest 3D distance from the black hole.
    after = (years > 2019) & (years < 2034)
    mas_to_au = orbit["r0_pc"] / 1000.0
    r3 = np.sqrt((einstein[:, 0] * mas_to_au) ** 2 + (einstein[:, 1] * mas_to_au) ** 2
                 + einstein[:, 2] ** 2)
    i_apo = np.flatnonzero(after)[np.argmax(r3[after])]
    gm = orbit["mass_1e6"] * 1e6 * eph.gm("sun")
    a = orbit["a_mas"] / 1000.0 * orbit["r0_pc"]
    x = gm / (eph.c_au_per_day() ** 2 * a * (1 - orbit["e"] ** 2))
    return {"years": years, "newton": newton, "einstein": einstein, "diff_mas": d,
            "apocentre_year": float(years[i_apo]), "shift_at_apocentre_mas": float(sep[i_apo]),
            "shift_at_apocentre_ra_dec_mas": [float(d[i_apo, 0]), float(d[i_apo, 1])],
            "precession_per_orbit_arcmin": float(np.degrees(6 * np.pi * x) * 60),
            "compactness": float(x), "a_au": float(a),
            "pericentre_au": float(a * (1 - orbit["e"])),
            "max_shift_1992_2040_mas": float(sep.max()),
            "max_shift_year": float(years[np.argmax(sep)])}


def rosette(orbit: dict | None = None, boost: float = 30.0, n_orbits: float = 12.0) -> dict:
    """S2's orbit with GR turned up ``boost`` times, for a picture of the rosette."""
    orbit = orbit or PUBLISHED["s2_orbit_gravity2020"]
    years = 2010.35 + np.linspace(0, n_orbits * orbit["period_yr"], 6000)
    return {"years": years, "sky": integrate(orbit, years, gr=True, boost=boost), "boost": boost}


# --------------------------------------------------------------------------
# What the public data could do (H1's analysis)
# --------------------------------------------------------------------------


def h1_results() -> dict | None:
    """H1's own precession fits to the public Keck + NTT/VLT positions, if built."""
    path = H1_DIR / "data" / "results.json"
    if not path.exists():
        return None
    r = json.loads(path.read_text())
    prof = r.get("precession_profile", {})
    rel = r.get("relativity", {}).get("joint_precession", {})
    return {"f_sp": rel.get("f_sp_median"), "f_sp_err": rel.get("f_sp_sd"),
            "profile_all": prof.get("all_positions"), "profile_ao_only": prof.get("ao_only"),
            "source": str(path.relative_to(PROJECT_DIR.parent))}


#: Typical 1-sigma astrometric errors per epoch (GRAVITY Collaboration 2020, Sect. 2).
INSTRUMENT_ERRORS_MAS = {
    "NTT speckle (SHARP), 1992-2002": (1992.0, 2002.0, 3.8),
    "VLT adaptive optics (NACO), 2002-2019": (2002.0, 2019.7, 0.40),
    "VLTI interferometry (GRAVITY), 2016.7-": (2016.7, 2040.0, 0.065),
}
