"""The physics: every formula the project uses, in one place.

Nothing here touches data files, so each function can be checked against a
known answer in ``test_physics.py``.

Units: distances in km unless a name says otherwise (``_au``, ``_ld``, ``_m``),
speeds in km/s, energies in megatons (MT) or kilotons (kt) of TNT as named.
"""

from __future__ import annotations

import numpy as np

from config import (
    ALBEDO_TYPICAL, AU_KM, BROWN2002_A, BROWN2002_B, BROWN2002_RAD_A, BROWN2002_RAD_B,
    DENSITY_KG_M3, EARTH_RADIUS_KM, GM_EARTH_KM3_S2, H_TO_D_KM, KT_TNT_J, LD_KM, MT_TNT_J,
    PALERMO_FB_A, PALERMO_FB_B,
)


# ---------------------------------------------------------------------------
# Distances
# ---------------------------------------------------------------------------

def au_to_km(d_au):
    return np.asarray(d_au, float) * AU_KM


def au_to_ld(d_au):
    """Astronomical units to lunar distances (1 LD = 384,400 km)."""
    return np.asarray(d_au, float) * (AU_KM / LD_KM)


def ld_to_au(d_ld):
    return np.asarray(d_ld, float) * (LD_KM / AU_KM)


# ---------------------------------------------------------------------------
# Size from brightness
# ---------------------------------------------------------------------------

def diameter_km(h, albedo=ALBEDO_TYPICAL):
    """Diameter from absolute magnitude H and geometric albedo p.

    H is how bright the asteroid would look 1 au from both Sun and observer,
    fully lit. Brightness is area x reflectivity, so a size needs an assumed
    albedo: D = 1329 km / sqrt(p) * 10^(-H/5).
    """
    return H_TO_D_KM / np.sqrt(albedo) * 10.0 ** (-np.asarray(h, float) / 5.0)


def absolute_magnitude(d_km, albedo=ALBEDO_TYPICAL):
    """Inverse of :func:`diameter_km`."""
    return 5.0 * np.log10(H_TO_D_KM / (np.sqrt(albedo) * np.asarray(d_km, float)))


def mass_kg(d_km, density=DENSITY_KG_M3):
    """Mass of a uniform sphere."""
    r_m = np.asarray(d_km, float) * 500.0
    return density * 4.0 / 3.0 * np.pi * r_m**3


# ---------------------------------------------------------------------------
# Encounter dynamics
# ---------------------------------------------------------------------------

def escape_speed(r_km=EARTH_RADIUS_KM):
    """Escape speed from Earth at distance r from its centre (11.19 km/s at the surface)."""
    return np.sqrt(2.0 * GM_EARTH_KM3_S2 / np.asarray(r_km, float))


def impact_speed(v_inf):
    """Speed at Earth's surface of an asteroid arriving at v_inf.

    Energy conservation: v^2 = v_inf^2 + v_esc^2. Even the slowest possible
    arrival hits at 11.2 km/s; Earth's own gravity supplies that much.
    """
    return np.sqrt(np.asarray(v_inf, float) ** 2 + escape_speed() ** 2)


def speed_at_distance(q_km, v_inf):
    """Speed at closest approach q for an arrival speed v_inf (vis-viva)."""
    return np.sqrt(np.asarray(v_inf, float) ** 2 + 2.0 * GM_EARTH_KM3_S2 / np.asarray(q_km, float))


def impact_energy_mt(d_km, v_kms, density=DENSITY_KG_M3):
    """Kinetic energy (1/2 m v^2) in megatons of TNT."""
    v_ms = np.asarray(v_kms, float) * 1e3
    return 0.5 * mass_kg(d_km, density) * v_ms**2 / MT_TNT_J


def diameter_for_energy_km(energy_mt, v_kms, density=DENSITY_KG_M3):
    """Inverse of :func:`impact_energy_mt`: the diameter carrying ``energy_mt``."""
    m = 2.0 * np.asarray(energy_mt, float) * MT_TNT_J / (np.asarray(v_kms, float) * 1e3) ** 2
    return 2.0 * (3.0 * m / (4.0 * np.pi * density)) ** (1.0 / 3.0) / 1e3


def impact_parameter_km(q_km, v_inf):
    """Aim distance b that bends into a closest approach q.

    Far from Earth an asteroid travels in a straight line, aimed to miss Earth's
    centre by b. Gravity pulls it in, so it actually passes closer, at q.
    Angular momentum (b v_inf = q v_q) and energy give
    b^2 = q^2 (1 + 2GM / (q v_inf^2)).
    """
    q = np.asarray(q_km, float)
    return q * np.sqrt(1.0 + 2.0 * GM_EARTH_KM3_S2 / (q * np.asarray(v_inf, float) ** 2))


def capture_radius_km(v_inf, radius_km=EARTH_RADIUS_KM):
    """Aim distance within which an asteroid hits: Earth's radius, enlarged by
    gravitational focusing. Slow asteroids are pulled in from farther away."""
    return impact_parameter_km(radius_km, v_inf)


def impact_weight(v_inf, ref_dist_km):
    """Ratio of Earth's capture cross-section to that of a sphere of radius
    ``ref_dist_km``, both measured in the plane of aim (the "b-plane").

    Asteroids crossing Earth's neighbourhood are spread evenly over that plane,
    so for every approach within ``ref_dist_km`` this fraction would have hit.
    Summing it over all approaches within the reference distance estimates the
    impact rate (Öpik's reasoning).
    """
    return (capture_radius_km(v_inf) / impact_parameter_km(ref_dist_km, v_inf)) ** 2


def hyperbolic_eccentricity(q_km, v_inf):
    """Eccentricity of a flyby: e = 1 + q v_inf^2 / GM (e > 1: not captured)."""
    return 1.0 + np.asarray(q_km, float) * np.asarray(v_inf, float) ** 2 / GM_EARTH_KM3_S2


def turning_angle_deg(q_km, v_inf):
    """How far Earth's gravity bends the path: delta = 2 arcsin(1/e)."""
    return np.degrees(2.0 * np.arcsin(1.0 / hyperbolic_eccentricity(q_km, v_inf)))


# ---------------------------------------------------------------------------
# Risk scales
# ---------------------------------------------------------------------------

def background_impact_rate(energy_mt):
    """Annual frequency of natural impacts at least this energetic, as defined
    for the Palermo scale (Chesley et al. 2002): f_B = 0.03 E^-4/5 per year."""
    return PALERMO_FB_A * np.asarray(energy_mt, float) ** (-PALERMO_FB_B)


def palermo_scale(impact_probability, energy_mt, years_to_impact):
    """Palermo Technical Scale: log10 of the risk relative to the background.

    PS = log10( p / (f_B * dT) ) compares the chance of this impact with the
    chance that *some* impact at least as energetic happens anyway before the
    same date. PS = 0 means as likely as the background; PS = -2, 1% of it.
    """
    return np.log10(np.asarray(impact_probability, float)
                    / (background_impact_rate(energy_mt) * np.asarray(years_to_impact, float)))


def _torino_line(log_e, log_p, k):
    """The diagonal region boundaries of the Torino chart: zero on the line
    through (p = 1%, E = 10^k MT) along which log E falls by 1.5 for every
    decade of rising probability."""
    return (log_e - k) / 3.0 + (log_p + 2.0) / 2.0


def torino_scale(impact_probability, energy_mt):
    """Torino Impact Hazard Scale (Binzel 2000, wording revised 2005), 0-10.

    Evaluated from impact probability p and kinetic energy E (MT) using the
    chart's region boundaries: diagonal lines, plus the thresholds p = 1%,
    p = 99%, and E = 1, 100 and 100,000 MT. Anything below 1 MT is 0 however
    likely: such objects burst in the atmosphere. Defined only for impacts
    within the next 100 years (the caller's job to check).
    """
    p = np.atleast_1d(np.asarray(impact_probability, float))
    e = np.atleast_1d(np.asarray(energy_mt, float))
    p, e = np.broadcast_arrays(p, e)
    with np.errstate(divide="ignore", invalid="ignore"):
        x, y = np.log10(e), np.log10(p)
    out = np.zeros(p.shape, dtype=int)
    valid = np.isfinite(x) & np.isfinite(y) & (x >= 0)

    certain = valid & (p >= 0.99)
    out[certain & (x < 2)] = 8
    out[certain & (x >= 2) & (x < 5)] = 9
    out[certain & (x >= 5)] = 10

    likely = valid & ~certain & (y >= -2)  # 1% or more, not certain
    out[likely & (x < 2)] = 3
    out[likely & (x >= 2) & (x < 5) & (_torino_line(x, y, 5) < 0)] = 4
    out[likely & (x >= 2) & (x < 5) & (_torino_line(x, y, 5) >= 0)] = 5
    out[likely & (x >= 5)] = 7

    unlikely = valid & (y < -2)
    out[unlikely & (_torino_line(x, y, -1) >= 0) & (_torino_line(x, y, 2) < 0)] = 1
    out[unlikely & (_torino_line(x, y, 2) >= 0) & (_torino_line(x, y, 5) < 0)] = 2
    out[unlikely & (_torino_line(x, y, 5) >= 0)] = 6
    return out if np.ndim(impact_probability) or np.ndim(energy_mt) else int(out[0])


TORINO_COLOURS = {  # the scale's own colour code, by level
    0: "white", 1: "green", 2: "yellow", 3: "yellow", 4: "yellow",
    5: "orange", 6: "orange", 7: "orange", 8: "red", 9: "red", 10: "red",
}


# ---------------------------------------------------------------------------
# Fireballs
# ---------------------------------------------------------------------------

def fireball_impact_energy_kt(radiated_energy_j):
    """Total impact energy from the energy radiated as light (Brown et al.
    2002): E = 8.2508 E_rad^0.885, both in kilotons. Most of a small asteroid's
    energy goes into heat and blast; only a few percent comes out as light."""
    e_rad_kt = np.asarray(radiated_energy_j, float) / KT_TNT_J
    return BROWN2002_RAD_A * e_rad_kt**BROWN2002_RAD_B


def brown2002_rate(energy_kt):
    """Impacts per year at least as energetic as E (kt), N = 3.7 E^-0.90,
    fitted by Brown et al. (2002) to satellite bolide data over 0.1-10 kt."""
    return BROWN2002_A * np.asarray(energy_kt, float) ** (-BROWN2002_B)


# ---------------------------------------------------------------------------
# Orbits around the Sun
# ---------------------------------------------------------------------------

GM_SUN_AU3_D2 = 0.01720209895**2  # Gaussian gravitational constant squared


def yarkovsky_drift_m_per_year(a2_au_d2, a_au, e):
    """Orbit-averaged drift of the semi-major axis from a transverse
    acceleration A2 (1 au / r)^2 -- JPL's standard Yarkovsky model.

    Sunlight absorbed by a spinning asteroid is re-radiated as heat a little
    later, from the afternoon side, and the photons' recoil pushes it along
    (or against) its orbit. Gauss's equation da/dt = (2a^2/h)(p/r)T, averaged
    over one orbit with <r^-3> = a^-3 (1-e^2)^-3/2, gives
    da/dt = 2 A2 / (sqrt(GM a) (1 - e^2)).
    """
    rate_au_per_day = 2.0 * a2_au_d2 / (np.sqrt(GM_SUN_AU3_D2 * a_au) * (1.0 - e**2))
    return rate_au_per_day * 365.25 * AU_KM * 1e3


def orbit_xyz(a_au, e, i_deg, node_deg, peri_deg, true_anomaly_deg):
    """Heliocentric position (au) on a Keplerian ellipse, in the frame the
    elements refer to (here the ecliptic). Vectorised over true anomaly."""
    nu = np.radians(np.asarray(true_anomaly_deg, float))
    i, node, peri = np.radians(i_deg), np.radians(node_deg), np.radians(peri_deg)
    r = a_au * (1 - e**2) / (1 + e * np.cos(nu))
    u = peri + nu  # argument of latitude
    x = r * (np.cos(node) * np.cos(u) - np.sin(node) * np.sin(u) * np.cos(i))
    y = r * (np.sin(node) * np.cos(u) + np.cos(node) * np.sin(u) * np.cos(i))
    z = r * np.sin(u) * np.sin(i)
    return x, y, z


# ---------------------------------------------------------------------------
# Small numerical helpers used by the analysis
# ---------------------------------------------------------------------------

def loglog_slope(d, n):
    """Least-squares slope of log n against log d."""
    d, n = np.asarray(d, float), np.asarray(n, float)
    ok = (d > 0) & (n > 0)
    return float(np.polyfit(np.log10(d[ok]), np.log10(n[ok]), 1)[0])
