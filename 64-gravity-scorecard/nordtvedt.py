"""The Nordtvedt effect, simulated: what a violation of the strong equivalence principle does to the Moon.

If gravitational self-energy fell differently from other mass-energy, the Earth
(whose binding energy is 4.6e-10 of its mass) and the Moon (0.19e-10) would fall
toward the Sun at slightly different rates: (m_G/m_I) = 1 - eta * E_g/(m c^2)
(Will 2014 eq. 35, eta = 4 beta - gamma - 3 in PPN; GR: eta = 0). The Moon's
orbit would then be pushed toward or away from the Sun, a range signal at the
synodic month, delta_r ~ C * eta * cos D (D = mean elongation of the Moon).

This is a labelled simulation with real initial conditions: the Sun, the eight
planets' systems, Pluto, the Earth and the Moon from DE440, integrated with
REBOUND/IAS15 three times -- with eta = 0, +1 and -1 -- where the extra force
multiplies only the Sun's pull on the Earth and on the Moon. The violation is
switched on smoothly over two years (2004-2006) so it does not ring the orbit,
then followed through the APOLLO years. (eta = +-1 is huge, but the response is
linear: the +1 and -1 runs differ only by sign; the test checks that.)

Output: the Earth-Moon distance change per unit eta at any epoch (a template for
the lunar laser ranging fit) and its cos D amplitude C, which Nordtvedt (1995)
and Will (2014) give as 13.1 m.
"""

from __future__ import annotations

import json

import erfa
import numpy as np
import rebound

import config as C
from ephem import default as default_ephemeris

#: Gravitational self-energies |E_g|/(m c^2) (Will 2014, sec. 4.4: Earth 4.6e-10,
#: Moon 0.2e-10; Williams et al. 2009/2012 use 4.64e-10 and 1.90e-11).
SELF_ENERGY = {"earth": 4.64e-10, "moon": 1.90e-11}

START = (2004, 1, 1)          # integration start (TDB)
RAMP_YEARS = 2.0              # smooth switch-on
END = (2021, 3, 1)
STEP_DAYS = 0.25
CACHE = C.DERIVED_DIR / "nordtvedt_sim.npz"

# NAIF ids and DE440 GMs (km^3/s^2, de440_tech-comments / gm_de440.tpc)
BODIES = [("sun", 10, 132712440041.279419), ("mercury", 1, 22031.868551), ("venus", 2, 324858.592000),
          ("earth", 399, 398600.435507), ("moon", 301, 4902.800118), ("mars", 4, 42828.375816),
          ("jupiter", 5, 126712764.100000), ("saturn", 6, 37940584.841800), ("uranus", 7, 5794556.400000),
          ("neptune", 8, 6836527.100580), ("pluto", 9, 975.500000)]


def _day_sec(y, m, d):
    jd0, mjd = erfa.cal2jd(y, m, d)
    return int(mjd - 51545), 43200.0


def run(eta: float, *, eph=None):
    """Integrate with a given eta; returns times (days since START) and Earth-Moon vectors (km)."""
    eph = eph or default_ephemeris()
    day0, sec0 = _day_sec(*START)
    sim = rebound.Simulation()
    sim.G = 1.0                                    # masses are GM in km^3/s^2, time unit = 1 s ... rescaled below
    # use units of km and days: GM in km^3/day^2
    for name, naif, gm in BODIES:
        p, v = eph.posvel(naif, day0, sec0)
        sim.add(m=gm * C.DAY**2, x=p[0, 0], y=p[1, 0], z=p[2, 0],
                vx=v[0, 0] * C.DAY, vy=v[1, 0] * C.DAY, vz=v[2, 0] * C.DAY, name=name)
    sim.integrator = "ias15"
    sim.move_to_com()
    ps = sim.particles
    i_sun, i_e, i_m = 0, 3, 4
    d_e = -eta * SELF_ENERGY["earth"]              # (m_G/m_I) - 1
    d_m = -eta * SELF_ENERGY["moon"]
    gm_sun = BODIES[0][2] * C.DAY**2
    ramp = RAMP_YEARS * 365.25

    def extra(simp):
        t = simp.contents.t
        s = 1.0 if t >= ramp else 0.5 * (1 - np.cos(np.pi * t / ramp))
        sun = ps[i_sun]
        for idx, dlt in ((i_e, d_e), (i_m, d_m)):
            p = ps[idx]
            dx, dy, dz = p.x - sun.x, p.y - sun.y, p.z - sun.z
            r3 = (dx * dx + dy * dy + dz * dz) ** 1.5
            f = -s * dlt * gm_sun / r3
            p.ax += f * dx
            p.ay += f * dy
            p.az += f * dz

    if eta != 0.0:
        sim.additional_forces = extra
    _, mjd_end = erfa.cal2jd(*END)
    t_end = (mjd_end - 51545 - day0) + 0.0
    times = np.arange(0.0, t_end, STEP_DAYS)
    em = np.empty((len(times), 3))
    for k, t in enumerate(times):
        sim.integrate(t)
        em[k] = (ps[i_m].x - ps[i_e].x, ps[i_m].y - ps[i_e].y, ps[i_m].z - ps[i_e].z)
    return times, em


def simulate(force: bool = False) -> dict:
    if CACHE.exists() and not force:
        z = np.load(CACHE)
        return {k: z[k] for k in z.files}
    out = {}
    for label, eta in (("nominal", 0.0), ("plus", 1.0), ("minus", -1.0)):
        t, em = run(eta)
        out["t"] = t
        out[label] = em
        print(f"  Nordtvedt run eta={eta:+.0f} done")
    np.savez(CACHE, **out)
    return out


def template(sim: dict | None = None):
    """Earth-Moon distance change per unit eta (m) vs days since START, and the linearity check."""
    sim = sim or simulate()
    r0 = np.linalg.norm(sim["nominal"], axis=1)
    dp = np.linalg.norm(sim["plus"], axis=1) - r0
    dm = np.linalg.norm(sim["minus"], axis=1) - r0
    per_eta = 0.5 * (dp - dm) * 1e3            # metres per unit eta
    nonlin = 0.5 * (dp + dm) * 1e3             # even part: zero if the response is linear
    return sim["t"], per_eta, nonlin


def mean_elongation(t_days):
    """Delaunay D (rad) at TDB days since START."""
    day0, sec0 = _day_sec(*START)
    t_cent = (day0 + sec0 / C.DAY + t_days) / 36525.0
    return np.mod(erfa.fad03(t_cent), 2 * np.pi), np.mod(erfa.fal03(t_cent), 2 * np.pi), np.mod(erfa.faf03(t_cent), 2 * np.pi)


def amplitude(sim: dict | None = None) -> dict:
    """Fit the switched-on part with cos/sin D plus free-oscillation and drift terms."""
    t, y, nonlin = template(sim)
    D, l, F = mean_elongation(t)
    sel = t > RAMP_YEARS * 365.25 + 30
    tt = (t[sel] - t[sel].mean()) / 365.25
    cols = {"cosD": np.cos(D), "sinD": np.sin(D), "cos2D": np.cos(2 * D), "sin2D": np.sin(2 * D),
            "cosl": np.cos(l), "sinl": np.sin(l), "cosD-l": np.cos(D - l), "sinD-l": np.sin(D - l),
            "cosD+l": np.cos(D + l), "sinD+l": np.sin(D + l)}
    names = list(cols)
    A = [cols[k][sel] for k in names] + [np.ones(sel.sum()), tt, np.cos(l[sel]) * tt, np.sin(l[sel]) * tt]
    A = np.column_stack(A)
    x, *_ = np.linalg.lstsq(A, y[sel], rcond=None)
    resid = y[sel] - A @ x
    out = {k: float(v) for k, v in zip(names, x[: len(names)])}
    out.update(rms_m=float(np.sqrt(np.mean(resid**2))), nonlinearity_max_m=float(np.abs(nonlin).max()),
               peak_m=float(np.abs(y[sel]).max()))
    return out


if __name__ == "__main__":
    s = simulate()
    a = amplitude(s)
    print(json.dumps(a, indent=1))
