"""Orbital physics of binary stars: Kepler's law and what each method sees.

Units throughout: periods in days, masses in solar masses, separations in AU,
velocities in km/s, angles in radians unless a name says otherwise.

Kepler's third law in solar units, a^3 = M_tot P^2 (a in AU, P in years, M in
Msun), sets every scale:

* the stars' relative orbit, a;
* the primary's radial-velocity semi-amplitude K1 (spectroscopy sees this);
* the wobble of the light centre, a0 (astrometry sees this). A companion that
  is as bright as it is massive (a twin) leaves the light centre still: an
  astrometric blind spot.
"""

from __future__ import annotations

import numpy as np

YEAR_DAYS = 365.25
#: Earth's mean orbital speed, 2 pi AU / yr, in km/s: the velocity unit of
#: Kepler's law in solar units.
V_EARTH_KMS = 2 * np.pi * 149_597_870.7 / (YEAR_DAYS * 86_400.0)


def semi_major_axis_au(period_days, m_tot):
    """Relative orbit's semi-major axis from Kepler's third law."""
    return (m_tot * (np.asarray(period_days) / YEAR_DAYS) ** 2) ** (1 / 3)


def period_days_from_a(a_au, m_tot):
    """Inverse of :func:`semi_major_axis_au`."""
    return YEAR_DAYS * np.sqrt(np.asarray(a_au) ** 3 / m_tot)


def k1_kms(period_days, m1, m2, ecc=0.0, sin_i=1.0):
    """Radial-velocity semi-amplitude of star 1 (km/s)."""
    m_tot = m1 + m2
    p_yr = np.asarray(period_days) / YEAR_DAYS
    v_rel = V_EARTH_KMS * m_tot ** (1 / 3) * p_yr ** (-1 / 3)  # mean relative orbital speed
    return v_rel * (m2 / m_tot) * sin_i / np.sqrt(1 - np.asarray(ecc) ** 2)


def mass_function(period_days, k1, ecc=0.0):
    """Spectroscopic mass function (m2 sin i)^3 / (m1 + m2)^2 in Msun.

    f = P K^3 (1-e^2)^(3/2) / (2 pi G); with G Msun = 4 pi^2 AU^3/yr^2 and K in
    units of 2 pi AU/yr this is simply P[yr] (K / V_earth)^3 (1-e^2)^(3/2).
    """
    p_yr = np.asarray(period_days) / YEAR_DAYS
    return p_yr * (np.asarray(k1) / V_EARTH_KMS) ** 3 * (1 - np.asarray(ecc) ** 2) ** 1.5


def photocentre_factor(q, flux_ratio):
    """|B - beta|: the light centre's orbit as a fraction of the relative orbit,
    with B = m2/(m1+m2) the mass fraction and beta = F2/(F1+F2) the light
    fraction of the companion. Zero for a twin."""
    q = np.asarray(q, dtype=float)
    f = np.asarray(flux_ratio, dtype=float)
    return np.abs(q / (1 + q) - f / (1 + f))


def a0_mas(period_days, m1, m2, flux_ratio, parallax_mas):
    """Angular semi-major axis of the light centre's orbit (mas)."""
    a = semi_major_axis_au(period_days, m1 + m2)
    return a * parallax_mas * photocentre_factor(m2 / m1, flux_ratio)


def light_weighted_k(k1, q, flux_ratio):
    """Semi-amplitude of the light-weighted mean velocity of a blended pair,
    (v1 + beta' v2) / (1 + beta') with v2 = -v1/q and beta' = F2/F1: what a
    single-lined analysis measures when the two stars' lines overlap."""
    f = np.asarray(flux_ratio, dtype=float)
    return np.abs(np.asarray(k1) * (1 - f / np.asarray(q)) / (1 + f))


def e_max_envelope(period_days, p_min_days=2.0):
    """Upper envelope of eccentricity, e_max = 1 - (P / 2 d)^(-2/3) (Moe & Di
    Stefano 2017, Eq. 3): the periastron distance a(1-e) cannot shrink below
    the separation of a 2-day circular orbit without the stars touching or
    being circularized. (Mazeh 2008 and Halbwachs et al. 2005 use other forms
    of the same idea.)"""
    p = np.asarray(period_days, dtype=float)
    return np.clip(1 - (p / p_min_days) ** (-2 / 3), 0, 1)


# ---------------------------------------------------------------------------
# Orbits in 3D: random orientations and phases
# ---------------------------------------------------------------------------

def solve_kepler(mean_anomaly, ecc, iterations: int = 30):
    """Eccentric anomaly E from M = E - e sin E (Newton's method, vectorised)."""
    m = np.mod(mean_anomaly, 2 * np.pi)
    e = np.asarray(ecc, dtype=float)
    big_e = np.where(e < 0.8, m, np.pi * np.ones_like(m))
    for _ in range(iterations):
        f = big_e - e * np.sin(big_e) - m
        big_e = big_e - f / (1 - e * np.cos(big_e))
    return big_e


def random_orientation(rng, n):
    """Isotropic orbit orientations: cos i uniform, omega and Omega uniform."""
    cos_i = rng.uniform(-1, 1, n)
    return np.arccos(cos_i), rng.uniform(0, 2 * np.pi, n), rng.uniform(0, 2 * np.pi, n)


def projected_separation_factor(rng, ecc, n=None):
    """s / a for a random orientation and a random moment in the orbit:
    the projected separation in units of the semi-major axis. Averaged over
    everything, <s/a> ~ 0.9-1.1 depending on eccentricity (time-averaged r/a is
    1 + e^2/2), with a long tail to small values."""
    ecc = np.asarray(ecc, dtype=float)
    if n is None:
        n = ecc.size
    inc, omega, _ = random_orientation(rng, n)
    mean_anom = rng.uniform(0, 2 * np.pi, n)
    big_e = solve_kepler(mean_anom, ecc)
    r = 1 - ecc * np.cos(big_e)
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(big_e / 2), np.sqrt(1 - ecc) * np.cos(big_e / 2))
    # position in the orbital plane rotated by omega, then tilted by i
    u = omega + nu
    x = r * np.cos(u)
    y = r * np.sin(u) * np.cos(inc)
    return np.hypot(x, y)


def radial_velocity(t_days, period, t_peri, ecc, omega, k):
    """Keplerian radial velocity curve (km/s) of one star."""
    mean_anom = 2 * np.pi * (np.asarray(t_days) - t_peri) / period
    big_e = solve_kepler(mean_anom, ecc)
    nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(big_e / 2), np.sqrt(1 - ecc) * np.cos(big_e / 2))
    return k * (np.cos(nu + omega) + ecc * np.cos(omega))


# ---------------------------------------------------------------------------
# Tidal circularization
# ---------------------------------------------------------------------------

def meibom_mathieu_e(period_days, p_prime, alpha=0.35, beta=0.14):
    """Meibom & Mathieu (2005) circularization function: e = 0 below P', then
    e = alpha (1 - exp(beta (P' - P))) above it."""
    p = np.asarray(period_days, dtype=float)
    return np.where(p > p_prime, alpha * (1 - np.exp(beta * (p_prime - p))), 0.0)


def meibom_mathieu_pcirc(p_prime, alpha=0.35, beta=0.14, e_circ=0.01):
    """The circularization period: where the fitted function reaches e = 0.01."""
    return p_prime - np.log(1 - e_circ / alpha) / beta
