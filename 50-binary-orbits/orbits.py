"""Visual-binary orbit mathematics: Kepler's equation, Campbell and Thiele-Innes elements.

Conventions (those of ORB6, the WDS and Bond et al. 2017):

* The companion's position relative to the primary is given by its separation
  rho (arcsec) and position angle theta, measured from north through east.
  In Cartesian form x = rho cos(theta) points north, y = rho sin(theta) east.
* Campbell elements: period P (Julian years), time of periastron T (decimal
  year), eccentricity e, angular semi-major axis a (arcsec), inclination i,
  position angle of the node Omega, and the argument of periastron omega,
  measured from the node in the direction of motion. i < 90 deg means the
  companion moves anticlockwise on the sky (theta increasing); i > 90 deg
  means clockwise, as for Sirius (i = 136 deg).
* Without radial velocities the ascending and descending nodes cannot be told
  apart: (Omega, omega) and (Omega + 180, omega + 180) give the same sky orbit.
  We report Omega in [0, 180), as ORB6 does for orbits with unidentified nodes.
* Thiele-Innes constants (A, B, F, G) make the position linear in the
  geometric elements once (P, T, e) are fixed:
      x = A X + F Y,   y = B X + G Y,
  with X = cos E - e, Y = sqrt(1 - e^2) sin E, E the eccentric anomaly.

Times are decimal years throughout. The Sirius data are Besselian years; ORB6
periastron epochs are Besselian years or (M)JD. ``jd_to_year`` and
``besselian_to_jd`` convert between them; the differences matter at the level
of 0.001 yr, far below any orbit's timing precision here.
"""

from __future__ import annotations

import numpy as np

import config as C

TWO_PI = 2.0 * np.pi
DEG = np.pi / 180.0

# ---------------------------------------------------------------------------
# Time scales
# ---------------------------------------------------------------------------

#: Besselian epoch B -> JD (Lieske 1979): JD = 2415020.31352 + (B - 1900) * 365.242198781
B1900_JD = 2415020.31352
TROPICAL_YEAR = 365.242198781


def besselian_to_jd(b):
    return B1900_JD + (np.asarray(b, float) - 1900.0) * TROPICAL_YEAR


def jd_to_besselian(jd):
    return 1900.0 + (np.asarray(jd, float) - B1900_JD) / TROPICAL_YEAR


def jd_to_year(jd):
    """JD -> Julian epoch (decimal year of 365.25 days, J2000.0 = JD 2451545.0)."""
    return 2000.0 + (np.asarray(jd, float) - 2451545.0) / C.JULIAN_YEAR


def year_to_jd(y):
    return 2451545.0 + (np.asarray(y, float) - 2000.0) * C.JULIAN_YEAR


def besselian_to_year(b):
    """Besselian epoch -> Julian epoch (they differ by < 0.002 yr over 1860-2030)."""
    return jd_to_year(besselian_to_jd(b))


# ---------------------------------------------------------------------------
# Kepler's equation
# ---------------------------------------------------------------------------


def solve_kepler(M, e, tol=1e-12, max_iter=50):
    """Eccentric anomaly E (radians) for mean anomaly M and eccentricity e (broadcasting)."""
    M = np.mod(np.asarray(M, float), TWO_PI)
    e = np.asarray(e, float)
    M, e = np.broadcast_arrays(M, e)
    # Danby's starter, then Halley iterations: cubic convergence up to e ~ 0.99
    E = M + 0.85 * e * np.where(np.sin(M) >= 0.0, 1.0, -1.0)
    for _ in range(max_iter):
        sE, cE = np.sin(E), np.cos(E)
        f = E - e * sE - M
        fp = 1.0 - e * cE
        fpp = e * sE
        dE = -f / (fp - 0.5 * f * fpp / fp)
        E = E + dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def mean_anomaly(t, P, T):
    return TWO_PI * (np.asarray(t, float) - T) / P


def unit_orbit(t, P, T, e):
    """(X, Y) of the unit orbit: X = cos E - e, Y = sqrt(1 - e^2) sin E."""
    E = solve_kepler(mean_anomaly(t, P, T), e)
    return np.cos(E) - e, np.sqrt(1.0 - e * e) * np.sin(E)


def unit_orbit_rates(t, P, T, e):
    """d(X, Y)/dt in units of 1/yr."""
    E = solve_kepler(mean_anomaly(t, P, T), e)
    Edot = (TWO_PI / P) / (1.0 - e * np.cos(E))
    return -np.sin(E) * Edot, np.sqrt(1.0 - e * e) * np.cos(E) * Edot


# ---------------------------------------------------------------------------
# Element conversions
# ---------------------------------------------------------------------------


def thiele_innes(a, i_deg, Omega_deg, omega_deg):
    """Campbell geometric elements -> Thiele-Innes constants (A, B, F, G), arcsec."""
    i, O, w = (np.asarray(v, float) * DEG for v in (i_deg, Omega_deg, omega_deg))
    ci, cO, sO, cw, sw = np.cos(i), np.cos(O), np.sin(O), np.cos(w), np.sin(w)
    A = a * (cw * cO - sw * sO * ci)
    B = a * (cw * sO + sw * cO * ci)
    F = a * (-sw * cO - cw * sO * ci)
    G = a * (-sw * sO + cw * cO * ci)
    return A, B, F, G


def campbell_from_thiele_innes(A, B, F, G):
    """Thiele-Innes constants -> (a, i, Omega, omega) in arcsec and degrees.

    A + G = a cos(w + O)(1 + cos i), B - F = a sin(w + O)(1 + cos i),
    A - G = a cos(w - O)(1 - cos i), -B - F = a sin(w - O)(1 - cos i);
    A^2 + B^2 + F^2 + G^2 = a^2 (1 + cos^2 i), AG - BF = a^2 cos i.
    Omega is folded into [0, 180) together with omega (node ambiguity).
    """
    A, B, F, G = (np.asarray(v, float) for v in (A, B, F, G))
    u = np.arctan2(B - F, A + G)          # omega + Omega
    v = np.arctan2(-B - F, A - G)         # omega - Omega
    w = 0.5 * (u + v)
    O = 0.5 * (u - v)
    k = 0.5 * (A * A + B * B + F * F + G * G)
    j = A * G - B * F
    a2 = k + np.sqrt(np.maximum(k * k - j * j, 0.0))
    a = np.sqrt(a2)
    cosi = np.clip(j / a2, -1.0, 1.0)
    i = np.arccos(cosi)
    O_deg, w_deg = np.degrees(O), np.degrees(w)
    O_deg, w_deg = fold_node(O_deg, w_deg)
    return a, np.degrees(i), O_deg, w_deg


def fold_node(Omega_deg, omega_deg):
    """Put Omega in [0, 180) (shifting omega by the same 180 deg), omega in [0, 360)."""
    O = np.asarray(Omega_deg, float)
    w = np.asarray(omega_deg, float)
    flip = np.mod(O, 360.0) >= 180.0
    O = np.mod(O, 360.0) - 180.0 * flip
    w = np.mod(w + 180.0 * flip, 360.0)
    return O, w


# ---------------------------------------------------------------------------
# Ephemerides
# ---------------------------------------------------------------------------


def relative_xy(t, P, T, e, a, i_deg, Omega_deg, omega_deg):
    """Companion position (x north, y east) in arcsec at times t."""
    A, B, F, G = thiele_innes(a, i_deg, Omega_deg, omega_deg)
    X, Y = unit_orbit(t, P, T, e)
    return A * X + F * Y, B * X + G * Y


def rho_theta(t, P, T, e, a, i_deg, Omega_deg, omega_deg):
    """Separation (arcsec) and position angle (deg, [0, 360)) at times t."""
    x, y = relative_xy(t, P, T, e, a, i_deg, Omega_deg, omega_deg)
    return np.hypot(x, y), np.mod(np.degrees(np.arctan2(y, x)), 360.0)


def relative_velocity(t, P, T, e, a, i_deg, Omega_deg, omega_deg):
    """Companion's sky velocity relative to the primary (dx/dt north, dy/dt east), arcsec/yr."""
    A, B, F, G = thiele_innes(a, i_deg, Omega_deg, omega_deg)
    Xd, Yd = unit_orbit_rates(t, P, T, e)
    return A * Xd + F * Yd, B * Xd + G * Yd


def relative_acceleration(t, P, T, e, a, i_deg, Omega_deg, omega_deg):
    """Sky-plane acceleration (north, east) of the companion relative to the primary, arcsec/yr^2.

    Kepler motion: d2r/dt2 = -GM r / |r|^3 in 3D; the sky projection follows
    from the 3D position (the line-of-sight coordinate z enters |r|).
    """
    x, y = relative_xy(t, P, T, e, a, i_deg, Omega_deg, omega_deg)
    z = line_of_sight(t, P, T, e, a, i_deg, Omega_deg, omega_deg)
    r3 = (x * x + y * y + z * z) ** 1.5
    gm = TWO_PI**2 * a**3 / P**2            # arcsec^3/yr^2 (angular units)
    return -gm * x / r3, -gm * y / r3


def line_of_sight(t, P, T, e, a, i_deg, Omega_deg, omega_deg):
    """Companion's offset along the line of sight (arcsec-equivalent).

    Sign: positive = companion farther than the primary, for the ascending node
    at Omega (the node where the companion recedes). Only |z| matters for the
    projected acceleration; the sign needs radial velocities to pin down.
    """
    i, w = np.asarray(i_deg) * DEG, np.asarray(omega_deg) * DEG
    X, Y = unit_orbit(t, P, T, e)
    # orbital-plane coordinates (periastron along +X) rotated by omega, then tilted by i
    return a * np.sin(i) * (np.sin(w) * X + np.cos(w) * Y)


def true_orbit_xy(t, P, T, e, a, omega_deg=0.0):
    """Position in the orbital plane seen face-on (periastron along +x, rotated by omega)."""
    X, Y = unit_orbit(t, P, T, e)
    w = np.asarray(omega_deg) * DEG
    return a * (np.cos(w) * X - np.sin(w) * Y), a * (np.sin(w) * X + np.cos(w) * Y)


def ellipse_xy(P, T, e, a, i_deg, Omega_deg, omega_deg, n=721):
    """The apparent (projected) orbit as a closed curve, for plotting."""
    t = T + P * np.linspace(0.0, 1.0, n)
    return relative_xy(t, P, T, e, a, i_deg, Omega_deg, omega_deg)


# ---------------------------------------------------------------------------
# Masses
# ---------------------------------------------------------------------------


def total_mass(a_arcsec, plx_arcsec, P_yr):
    """Kepler's third law: M1 + M2 (Msun) = (a/parallax)^3 / P^2, P in Julian years."""
    return C.KEPLER * (np.asarray(a_arcsec) / np.asarray(plx_arcsec)) ** 3 / np.asarray(P_yr) ** 2


def total_mass_frac_err(a, a_err, plx, plx_err, P, P_err):
    """Fractional error of the total mass from uncorrelated errors in a, parallax and P."""
    return np.sqrt((3 * a_err / a) ** 2 + (3 * plx_err / plx) ** 2 + (2 * P_err / P) ** 2)


def wrap_angle_deg(d):
    """Angle difference folded to (-180, 180]."""
    return -np.mod(-np.asarray(d, float) + 180.0, 360.0) + 180.0
