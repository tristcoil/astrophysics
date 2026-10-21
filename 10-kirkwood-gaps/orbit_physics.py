"""The formulas: Kepler's third law, resonance locations, and orbits -> positions.

Conventions: distances in AU, times in days unless a name says otherwise,
angles in degrees at the interfaces and radians inside. The asteroid is
treated as massless, so Kepler's third law reads n = k / a^(3/2) with k the
Gaussian gravitational constant. That is exactly the convention the MPC uses
for the mean daily motion it publishes alongside each orbit.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np

from config import GAUSS_K, JUPITER_PERIOD_DAYS, MASS_RATIO, SATURN_PERIOD_DAYS

DEG = np.pi / 180.0
DAYS_PER_YEAR = 365.25


# --------------------------------------------------------------------------
# Kepler's third law
# --------------------------------------------------------------------------

def mean_motion(a_au, mass_ratio: float = 0.0):
    """Mean motion in degrees/day for semi-major axis ``a_au``.

    ``mass_ratio`` is the orbiting body's mass over the Sun's; zero for an
    asteroid. (Including it matters for Jupiter: its own mass speeds its
    orbit by 0.05%, which moves every resonance by ~0.001 AU.)
    """
    a = np.asarray(a_au, dtype=float)
    return GAUSS_K * np.sqrt(1.0 + mass_ratio) / a ** 1.5 / DEG


def period_days(a_au, mass_ratio: float = 0.0):
    return 360.0 / mean_motion(a_au, mass_ratio)


def semi_major_axis(period_days_, mass_ratio: float = 0.0):
    """Invert Kepler's third law: orbital period (days) -> a (AU)."""
    n_rad = 2 * np.pi / np.asarray(period_days_, dtype=float)
    return (GAUSS_K ** 2 * (1.0 + mass_ratio) / n_rad ** 2) ** (1.0 / 3.0)


def jupiter_semi_major_axis() -> float:
    """Jupiter's mean semi-major axis implied by its mean period (~5.2026 AU)."""
    return float(semi_major_axis(JUPITER_PERIOD_DAYS, MASS_RATIO["Jupiter"]))


def period_ratio(a_au):
    """How many times an asteroid at ``a_au`` orbits per orbit of Jupiter,
    P_Jupiter / P_asteroid = n_asteroid / n_Jupiter."""
    return JUPITER_PERIOD_DAYS / period_days(a_au)


def a_from_period_ratio(ratio):
    """Semi-major axis at which an asteroid makes ``ratio`` orbits per
    Jupiter orbit -- the inverse of :func:`period_ratio`."""
    return semi_major_axis(JUPITER_PERIOD_DAYS / np.asarray(ratio, dtype=float))


# --------------------------------------------------------------------------
# Resonances
# --------------------------------------------------------------------------

def resonance_a(p: int, q: int) -> float:
    """Location of the p:q mean-motion resonance with Jupiter, in AU.

    Written the way the Kirkwood gaps are named: in the 3:1 resonance the
    asteroid completes 3 orbits while Jupiter completes 1. So n/n_J = p/q.
    """
    return float(a_from_period_ratio(Fraction(p, q)))


def resonance_order(p: int, q: int) -> int:
    """|p - q|. Resonance strength scales roughly as e^order, so first-order
    resonances (2:1, 3:2) are strongest and 7:3 (order 4) is weak."""
    return abs(p - q)


def three_body_resonance_a(m_j: int, m_s: int, m: int) -> float:
    """Location of the three-body resonance m_j*n_J + m_s*n_S + m*n = 0,
    involving Jupiter, Saturn and the asteroid (Nesvorny & Morbidelli 1998).
    E.g. (5, -2, -2): the "5J-2S-2" resonance near 3.17 AU."""
    n_j = 360.0 / JUPITER_PERIOD_DAYS
    n_s = 360.0 / SATURN_PERIOD_DAYS
    n = -(m_j * n_j + m_s * n_s) / m
    return float(semi_major_axis(360.0 / n))


def resonant_angle(p: int, q: int, lam, lam_j, varpi):
    """Critical angle of the p:q resonance, degrees in [-180, 180).

    sigma = p*lambda_J - q*lambda - (p - q)*varpi, where lambda is mean
    longitude and varpi longitude of perihelion. Inside the resonance sigma
    *librates* -- rocks back and forth about a fixed value -- instead of
    circulating through all angles. For the Hildas (3:2), sigma = 0 means
    every conjunction with Jupiter happens at the asteroid's perihelion,
    as far from Jupiter's orbit as it ever gets.
    """
    sigma = p * np.asarray(lam_j) - q * np.asarray(lam) - (p - q) * np.asarray(varpi)
    return wrap180(sigma)


def wrap180(angle_deg):
    return (np.asarray(angle_deg) + 180.0) % 360.0 - 180.0


def wrap360(angle_deg):
    return np.asarray(angle_deg) % 360.0


# --------------------------------------------------------------------------
# Orbit geometry
# --------------------------------------------------------------------------

def solve_kepler(mean_anomaly_rad, e, *, tol: float = 1e-12, max_iter: int = 50):
    """Eccentric anomaly E from M = E - e sin E, vectorised Newton iteration.

    Starting from E = M + e sin M (or pi for very eccentric orbits) Newton
    converges in a handful of steps for every e < 1.
    """
    M = np.mod(np.asarray(mean_anomaly_rad, dtype=float), 2 * np.pi)
    e = np.asarray(e, dtype=float)
    E = np.where(e < 0.8, M + e * np.sin(M), np.pi * np.ones_like(M))
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        dE = f / (1.0 - e * np.cos(E))
        E = E - dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def elements_to_xyz(a, e, incl_deg, node_deg, peri_deg, mean_anomaly_deg):
    """Heliocentric position (AU) in the frame the elements are referred to
    (for MPC and Horizons: the J2000 ecliptic). Returns x, y, z arrays."""
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)
    i = np.asarray(incl_deg) * DEG
    O = np.asarray(node_deg) * DEG
    w = np.asarray(peri_deg) * DEG
    E = solve_kepler(np.asarray(mean_anomaly_deg) * DEG, e)
    # Position in the orbital plane, x toward perihelion.
    xp = a * (np.cos(E) - e)
    yp = a * np.sqrt(1.0 - e * e) * np.sin(E)
    cO, sO, cw, sw, ci, si = np.cos(O), np.sin(O), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    x = (cO * cw - sO * sw * ci) * xp + (-cO * sw - sO * cw * ci) * yp
    y = (sO * cw + cO * sw * ci) * xp + (-sO * sw + cO * cw * ci) * yp
    z = (sw * si) * xp + (cw * si) * yp
    return x, y, z


def propagate_mean_anomaly(mean_anomaly_deg, n_deg_per_day, dt_days):
    """Two-body propagation: M advances linearly at the mean motion."""
    return wrap360(np.asarray(mean_anomaly_deg) + np.asarray(n_deg_per_day) * dt_days)


def ellipse_xy(a, e, peri_longitude_deg, n: int = 361):
    """Points around an orbit projected on the ecliptic (inclination ignored),
    for drawing planetary orbits."""
    E = np.linspace(0, 2 * np.pi, n)
    xp = a * (np.cos(E) - e)
    yp = a * np.sqrt(1 - e * e) * np.sin(E)
    w = peri_longitude_deg * DEG
    return xp * np.cos(w) - yp * np.sin(w), xp * np.sin(w) + yp * np.cos(w)


def e_mars_crossing(a_au, q_limit: float):
    """Eccentricity above which perihelion a(1-e) falls inside ``q_limit``."""
    return 1.0 - q_limit / np.asarray(a_au, dtype=float)


def e_aphelion_crossing(a_au, Q_limit: float):
    """Eccentricity above which aphelion a(1+e) reaches beyond ``Q_limit``."""
    return Q_limit / np.asarray(a_au, dtype=float) - 1.0


# --------------------------------------------------------------------------
# Sizes
# --------------------------------------------------------------------------

def diameter_km(H, albedo):
    """Asteroid diameter from absolute magnitude: D = 1329 km / sqrt(p) 10^(-H/5)
    (e.g. Pravec & Harris 2007). An albedo guess of 0.05 vs 0.25 changes D by
    a factor of 2.2, so these are scale indications, not measurements."""
    return 1329.0 / np.sqrt(albedo) * 10 ** (-np.asarray(H, dtype=float) / 5.0)


def longitude_from_jupiter(lam, lam_j):
    """Mean longitude measured from Jupiter's, degrees in [-180, 180).
    Trojans cluster at +60 (L4, leading Jupiter) and -60 (L5, trailing)."""
    return wrap180(np.asarray(lam) - np.asarray(lam_j))
