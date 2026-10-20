"""Where everything is: DSCOVR, the Sun, the Moon, and what EPIC sees.

All positions are geocentric J2000 vectors in km (the EPIC API's frame).
Natural-colour images are rotated so that Earth's north is up; seen from
outside the globe, east is then to the right.
"""

from __future__ import annotations

import numpy as np

from config import AU_KM, R_EARTH_EQ_KM, R_EARTH_MEAN_KM, R_MOON_KM

#: J2000 mean pole. The true pole of date differs by < 0.3 deg over 2015-2026
#: (precession), a rotation of < 4 pixels at the limb of a 2048-pixel image.
POLE = np.array([0.0, 0.0, 1.0])


def norm(v: np.ndarray) -> np.ndarray:
    return np.linalg.norm(v, axis=-1)


def unit(v: np.ndarray) -> np.ndarray:
    return v / norm(v)[..., None]


def angle_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    c = np.sum(unit(a) * unit(b), axis=-1)
    return np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))


def phase_angle_deg(dsc: np.ndarray, sun: np.ndarray) -> np.ndarray:
    """Sun-Earth-DSCOVR angle: Earth's phase angle as EPIC sees it."""
    return angle_deg(dsc, sun)


def lambert_phase(alpha_deg: np.ndarray) -> np.ndarray:
    """Disc-integrated brightness of a Lambert sphere relative to full phase,
    Phi(a) = [sin a + (pi - a) cos a] / pi."""
    a = np.radians(alpha_deg)
    return (np.sin(a) + (np.pi - a) * np.cos(a)) / np.pi


def sunlit_fraction_of_disc(alpha_deg: np.ndarray) -> np.ndarray:
    """Fraction of the projected disc that is sunlit, (1 + cos a) / 2."""
    return 0.5 * (1.0 + np.cos(np.radians(alpha_deg)))


def angular_radius_rad(dist_km: np.ndarray, radius_km: float = R_EARTH_EQ_KM) -> np.ndarray:
    return np.arcsin(radius_km / dist_km)


def image_axes(dsc: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Line of sight (DSCOVR -> Earth), image up (north) and image right (east)."""
    los = unit(-dsc)
    up = POLE - np.sum(POLE * los, axis=-1)[..., None] * los
    up = unit(up)
    right = np.cross(los, up)
    return los, up, right


def sky_offset_rad(dsc: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Angular offset (right, up) of a target from Earth's centre, seen from DSCOVR."""
    los, up, right = image_axes(dsc)
    w = unit(target - dsc)
    along = np.sum(w * los, axis=-1)
    return (np.arctan2(np.sum(w * right, axis=-1), along),
            np.arctan2(np.sum(w * up, axis=-1), along))


def moon_in_view(dsc: np.ndarray, moon: np.ndarray, sun: np.ndarray) -> dict:
    """The Moon as seen from DSCOVR: offset from Earth's centre, sizes, phase."""
    x, y = sky_offset_rad(dsc, moon)
    d_moon = norm(moon - dsc)
    d_earth = norm(dsc)
    return {
        "x_rad": x, "y_rad": y, "sep_rad": np.hypot(x, y),
        "r_moon_rad": angular_radius_rad(d_moon, R_MOON_KM),
        "r_earth_rad": angular_radius_rad(d_earth, R_EARTH_EQ_KM),
        "moon_dist_km": d_moon,
        "in_front": d_moon < d_earth,  # transit (in front) or occultation (behind)
        # Sun-Moon-DSCOVR angle: the Moon's phase angle as EPIC sees it.
        "moon_phase_deg": angle_deg(sun - moon, dsc - moon),
    }


def sub_point(v: np.ndarray, gmst_rad: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Geocentric latitude and longitude (deg) under a J2000 direction, given
    Greenwich mean sidereal time. Ignores precession/nutation (< 0.3 deg)."""
    u = unit(v)
    lat = np.degrees(np.arcsin(np.clip(u[..., 2], -1, 1)))
    ra = np.arctan2(u[..., 1], u[..., 0])
    lon = np.degrees(ra - gmst_rad)
    return lat, (lon + 180.0) % 360.0 - 180.0


def gmst_rad(t) -> np.ndarray:
    """Greenwich mean sidereal time (IAU 1982), from UTC timestamps. UT1-UTC
    (< 0.9 s) is ignored: 0.004 deg of rotation."""
    import pandas as pd
    jd = pd.DatetimeIndex(t).to_julian_date().to_numpy()
    d = jd - 2451545.0
    tc = d / 36525.0
    gmst_deg = 280.46061837 + 360.98564736629 * d + 0.000387933 * tc**2 - tc**3 / 38710000.0
    return np.radians(gmst_deg % 360.0)


def sun_earth_distance_au(sun: np.ndarray) -> np.ndarray:
    return norm(sun) / AU_KM


def mean_radius_km() -> float:
    return R_EARTH_MEAN_KM
