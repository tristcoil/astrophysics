"""JPL's DE440 planetary ephemeris: its constants and its trajectories.

DE440 comes as two files:

* ``header.440`` -- a text file of the 645 constants the integration used:
  the speed of light, the astronomical unit, the GM of the Sun and every planet
  system, the Sun's radius, oblateness (J2) and spin, and so on;
* ``de440.bsp`` -- a binary SPICE kernel holding Chebyshev polynomials for the
  position of each body relative to the solar-system barycentre, 1550-2650.

``jplephem`` evaluates the polynomials (positions in km; their derivative gives
velocities in km/day). Everything here returns au and au/day, in the ICRF
(equatorial) frame that DE440 is integrated in, with times as TDB Julian dates.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from config import BODY, DE440_BSP, DE440_HEADER, OBLIQUITY_J2000_ARCSEC


class MissingData(FileNotFoundError):
    pass


@lru_cache(maxsize=1)
def constants() -> dict[str, float]:
    """All 645 named constants of DE440, from its header file (groups 1040/1041)."""
    if not DE440_HEADER.exists():
        raise MissingData(f"{DE440_HEADER} missing -- run `python fetch_data.py` first")
    text = DE440_HEADER.read_text()
    names_block = text.split("GROUP   1040")[1].split("GROUP   1041")[0].split()
    values_block = text.split("GROUP   1041")[1].split("GROUP   1050")[0].split()
    n_names, n_values = int(names_block[0]), int(values_block[0])
    names = names_block[1:1 + n_names]
    values = [float(v.replace("D", "E")) for v in values_block[1:1 + n_values]]
    if n_names != n_values or len(names) != n_names or len(values) != n_values:
        raise ValueError(f"header.440: {n_names} names but {n_values} values")
    return dict(zip(names, values))


def gm(key: str) -> float:
    """GM of a body in au^3/day^2 (DE440's own units)."""
    return constants()[BODY[key].gm_key]


def au_km() -> float:
    return constants()["AU"]


def c_au_per_day() -> float:
    """The speed of light in au/day, from DE440's CLIGHT (km/s) and AU (km)."""
    c = constants()
    return c["CLIGHT"] * 86400.0 / c["AU"]


def sun_physical() -> dict:
    """The Sun's radius, oblateness and spin as DE440 adopts them.

    MOISUN is the moment of inertia in units of M R^2 and OMGSUN the rotation
    rate in rad/day; together they give the angular momentum that drives the
    Lense-Thirring effect (~1.9e41 kg m^2/s).
    """
    c = constants()
    return {"radius_au": c["ASUN"] / c["AU"], "radius_km": c["ASUN"], "j2": c["J2SUN"],
            "moment_of_inertia_factor": c["MOISUN"], "omega_rad_per_day": c["OMGSUN"]}


@lru_cache(maxsize=1)
def kernel():
    if not DE440_BSP.exists():
        raise MissingData(f"{DE440_BSP} missing -- run `python fetch_data.py` first")
    from jplephem.spk import SPK

    return SPK.open(str(DE440_BSP))


def coverage_jd() -> tuple[float, float]:
    """First and last Julian dates every segment of the kernel covers."""
    segs = kernel().segments
    return max(s.start_jd for s in segs), min(s.end_jd for s in segs)


def barycentric_state(key: str, jd) -> tuple[np.ndarray, np.ndarray]:
    """Position (au) and velocity (au/day) of a body relative to the barycentre.

    ``jd`` may be a scalar or an array; arrays give (3, N) outputs.
    """
    pos, vel = kernel()[0, BODY[key].naif].compute_and_differentiate(jd)
    a = au_km()
    return np.asarray(pos) / a, np.asarray(vel) / a


def heliocentric_state(key: str, jd) -> tuple[np.ndarray, np.ndarray]:
    """Position and velocity relative to the Sun's centre (au, au/day)."""
    r, v = barycentric_state(key, jd)
    rs, vs = barycentric_state("sun", jd)
    return r - rs, v - vs


# --------------------------------------------------------------------------
# Frames
# --------------------------------------------------------------------------

def rot_icrf_to_ecliptic() -> np.ndarray:
    """Rotation from the ICRF to the (fixed) ecliptic of J2000, as Horizons uses."""
    eps = np.deg2rad(OBLIQUITY_J2000_ARCSEC / 3600.0)
    ce, se = np.cos(eps), np.sin(eps)
    return np.array([[1.0, 0.0, 0.0], [0.0, ce, se], [0.0, -se, ce]])


def unit_vector_radec(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra, dec = np.deg2rad(ra_deg), np.deg2rad(dec_deg)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
