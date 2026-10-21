"""Orbit mathematics: Kepler's equation, sky positions, state vectors.

Conventions (the same as orbitize! and Wang et al. 2018, checked against
orbitize.kepler.calc_orbit in test_physics.py):

* Sky frame: x = Delta RA (east), y = Delta Dec (north), z = away from the
  observer. This frame is right-handed.
* Omega is the position angle (east of north) of the ascending node, the node
  where the planet moves away from us; omega is the planet's argument of
  periastron; i < 90 deg means the planet moves anticlockwise on the sky
  (position angle increasing), as all four HR 8799 planets do.
* tau is the epoch of periastron in units of the period after TAU_REF_MJD.

Multi-planet geometry: each planet's orbit is a Keplerian ellipse in Jacobi
coordinates (planet relative to the barycentre of the star and the planets
inside it, with mass M_k = M_star + sum of masses up to and including k).
The star itself is displaced from that barycentre, so the position of planet
k relative to the star is, exactly,

    X_k = r'_k + sum_{j inside k} (m_j / M_{k-1}) X_j ,

which adds a wobble of up to ~10 mas to the outer planets. With all planet
masses zero this reduces to independent Keplerian orbits (Wang et al.'s model).
"""

from __future__ import annotations

import numpy as np

import config as C

TWO_PI = 2.0 * np.pi
#: sqrt(G Msun) in au^1.5 / day, so that P[days] = 2 pi sqrt(a^3 / M) / K_GAUSS
K_GAUSS = np.sqrt(C.GM_SUN) * C.DAY / C.AU**1.5


def period_days(a_au, mass_msun):
    return TWO_PI * np.sqrt(np.asarray(a_au) ** 3 / np.asarray(mass_msun)) / K_GAUSS


def sma_from_period(p_days, mass_msun):
    return (np.asarray(p_days) * K_GAUSS / TWO_PI) ** (2.0 / 3.0) * np.asarray(mass_msun) ** (1.0 / 3.0)


def solve_kepler(M, e, tol=1e-11, max_iter=60):
    """Eccentric anomaly E for mean anomaly M (radians) and eccentricity e (arrays)."""
    M = np.mod(M, TWO_PI)
    e = np.broadcast_to(e, M.shape)
    E = M + 0.85 * e * np.where(np.sin(M) >= 0, 1.0, -1.0)      # Danby's starter
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        fpp = e * np.sin(E)
        # Halley's method: cubic convergence, robust up to e = 0.999
        dE = -f / (fp - 0.5 * f * fpp / fp)
        E = E + dE
        if np.max(np.abs(dE)) < tol:
            break
    return E


def mean_anomaly(t_mjd, period_d, tau, tau_ref=C.TAU_REF_MJD):
    """orbitize!'s definition: M = 2 pi [frac((t - t_ref)/P) - tau]."""
    return TWO_PI * (np.mod((t_mjd - tau_ref) / period_d, 1.0) - tau)


def orbit_xyz(t_mjd, a, e, inc, aop, pan, tau, mass, tau_ref=C.TAU_REF_MJD, velocity=False):
    """Position (au) in the sky frame, and optionally velocity (au/yr).

    All element arguments broadcast against ``t_mjd``. ``mass`` is the mass
    (Msun) setting the period (M_k for a Jacobi orbit).
    """
    P = period_days(a, mass)
    M = mean_anomaly(t_mjd, P, tau, tau_ref)
    E = solve_kepler(M, e)
    cosE, sinE = np.cos(E), np.sin(E)
    sq = np.sqrt(1.0 - e * e)
    # position in the orbital plane (periastron along +X)
    X = a * (cosE - e)
    Y = a * sq * sinE
    co, so = np.cos(aop), np.sin(aop)
    cO, sO = np.cos(pan), np.sin(pan)
    ci, si = np.cos(inc), np.sin(inc)
    # rotation: orbital plane -> sky (x east, y north, z away)
    # node direction (u = 0) is (sin O, cos O, 0); u = 90 deg is (cos i cos O, -cos i sin O, sin i)
    Ax, Ay, Az = co * sO + so * ci * cO, co * cO - so * ci * sO, so * si      # unit vector P (periastron)
    Bx, By, Bz = -so * sO + co * ci * cO, -so * cO - co * ci * sO, co * si   # unit vector Q
    x = X * Ax + Y * Bx
    y = X * Ay + Y * By
    z = X * Az + Y * Bz
    if not velocity:
        return x, y, z
    n = TWO_PI / (P / 365.25)                 # rad / yr
    Edot = n / (1.0 - e * cosE)
    Xd = -a * sinE * Edot
    Yd = a * sq * cosE * Edot
    return (x, y, z), (Xd * Ax + Yd * Bx, Xd * Ay + Yd * By, Xd * Az + Yd * Bz)


def sky_offsets(t_mjd, a, e, inc, aop, pan, tau, plx, mass, tau_ref=C.TAU_REF_MJD):
    """(Delta RA, Delta Dec) in mas of a single Keplerian orbit."""
    x, y, _ = orbit_xyz(t_mjd, a, e, inc, aop, pan, tau, mass, tau_ref)
    return x * plx, y * plx


def system_offsets(t_mjd, elements, plx, mstar, masses=None, tau_ref=C.TAU_REF_MJD):
    """Positions of all planets relative to the star at times ``t_mjd``.

    ``elements``: dict planet -> (a, e, inc, aop, pan, tau), each a scalar or
    an array broadcastable with t (Jacobi elements). ``masses``: dict planet ->
    Msun (None = all zero, i.e. independent Keplerian orbits). Returns dict
    planet -> (ra, dec) in mas.
    """
    masses = masses or {}
    out = {}
    m_inside = mstar                       # M_{k-1}: star + planets inside k
    for p in C.INNER_TO_OUTER:
        if p not in elements:
            continue
        m_p = masses.get(p, 0.0)
        a, e, inc, aop, pan, tau = elements[p]
        ra, dec = sky_offsets(t_mjd, a, e, inc, aop, pan, tau, plx, m_inside + m_p, tau_ref)
        # the star sits off the inner barycentre: + sum_j (m_j / M_{k-1}) X_j
        for q, (ra_q, dec_q) in out.items():
            m_q = masses.get(q, 0.0)
            ra = ra + (m_q / m_inside) * ra_q
            dec = dec + (m_q / m_inside) * dec_q
        out[p] = (ra, dec)
        m_inside = m_inside + m_p
    return out


def state_vectors(t_mjd, elements, mstar, masses, tau_ref=C.TAU_REF_MJD):
    """Barycentric positions (au) and velocities (au/yr) of star + planets.

    Jacobi Keplerian elements -> inertial Cartesian state at one epoch, in the
    sky frame (x east, y north, z away), ready for REBOUND. Returns a list of
    (m, x, y, z, vx, vy, vz), the star first, then planets inner to outer.
    """
    astro_r, astro_v = {}, {}
    m_inside = mstar
    for p in C.INNER_TO_OUTER:
        a, e, inc, aop, pan, tau = elements[p]
        m_p = masses[p]
        (x, y, z), (vx, vy, vz) = orbit_xyz(t_mjd, a, e, inc, aop, pan, tau, m_inside + m_p,
                                            tau_ref, velocity=True)
        r = np.array([x, y, z], float)
        v = np.array([vx, vy, vz], float)
        # shift from the inner barycentre to the star: + sum_j m_j/M_{k-1} X_j
        for q in astro_r:
            r = r + masses[q] / m_inside * astro_r[q]
            v = v + masses[q] / m_inside * astro_v[q]
        astro_r[p], astro_v[p] = r, v
        m_inside += m_p
    mtot = m_inside
    r_star = -sum(masses[p] * astro_r[p] for p in astro_r) / mtot
    v_star = -sum(masses[p] * astro_v[p] for p in astro_v) / mtot
    out = [(mstar, *r_star, *v_star)]
    for p in C.INNER_TO_OUTER:
        r = r_star + astro_r[p]
        v = v_star + astro_v[p]
        out.append((masses[p], *r, *v))
    return out


def wrap180(angle_deg):
    """Angles folded to [0, 180), the convention of Wang et al.'s tables."""
    return np.mod(angle_deg, 180.0)


def elements_from_state(r, v, mu):
    """Sky-frame state -> (a, e, inc, aop, pan, M) - the inverse of orbit_xyz.

    r (au), v (au/yr): arrays (..., 3) with x east, y north, z away;
    mu = G M in au^3/yr^2. Returns a <= 0 for unbound states.
    """
    r = np.asarray(r, float)
    v = np.asarray(v, float)
    rn = np.linalg.norm(r, axis=-1)
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=-1)
    inc = np.arccos(np.clip(-h[..., 2] / hn, -1, 1))          # h points at us for i < 90 deg
    pan = np.mod(np.arctan2(-h[..., 1], h[..., 0]), TWO_PI)    # node vector N = z x h
    N = np.stack([-h[..., 1], h[..., 0], np.zeros_like(hn)], axis=-1)
    Nn = np.linalg.norm(N, axis=-1)
    evec = np.cross(v, h) / np.asarray(mu)[..., None] - r / rn[..., None]
    e = np.linalg.norm(evec, axis=-1)
    a = 1.0 / (2.0 / rn - np.sum(v * v, axis=-1) / mu)
    cos_w = np.sum(N * evec, axis=-1) / (Nn * e)
    aop = np.arccos(np.clip(cos_w, -1, 1))
    aop = np.where(evec[..., 2] < 0, TWO_PI - aop, aop)
    cos_nu = np.sum(evec * r, axis=-1) / (e * rn)
    nu = np.arccos(np.clip(cos_nu, -1, 1))
    nu = np.where(np.sum(r * v, axis=-1) < 0, TWO_PI - nu, nu)
    E = 2.0 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(0.5 * nu))
    M = np.mod(E - e * np.sin(E), TWO_PI)
    return a, e, inc, aop, pan, M
