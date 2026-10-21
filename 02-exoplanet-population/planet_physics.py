"""The physics used to put the planet population in context.

Nothing in this module touches the catalogue. It holds the formulas that turn
catalogue columns into physical statements -- how much starlight a planet gets,
how hard its star wobbles, how likely its orbit is to transit -- plus the
published habitable-zone fits and a few Solar System reference values.

Constants come from ``astropy.constants`` (IAU 2015 nominal values) rather than
hand-copied numbers. Functions take plain floats or numpy arrays in the units
named in their arguments (``_au``, ``_msun``, ``_rearth``...), because they are
applied column-wise to a table of thousands of planets; ``test_physics.py``
checks each one against a case with a known answer.
"""

from __future__ import annotations

import numpy as np
from astropy import constants as const

# ---------------------------------------------------------------------------
# Constants (SI), from astropy's IAU 2015 nominal values
# ---------------------------------------------------------------------------

G = const.G.value
C = const.c.value
AU = const.au.value
PC = const.pc.value
M_SUN = const.M_sun.value
R_SUN = const.R_sun.value
M_EARTH = const.M_earth.value
R_EARTH = const.R_earth.value  # equatorial, 6378.1 km -- the exoplanet convention
M_JUP = const.M_jup.value
R_JUP = const.R_jup.value      # equatorial, 71492 km
DAY = 86400.0

#: Nominal solar effective temperature (IAU 2015 Resolution B3).
T_SUN = 5772.0

#: Jupiter's mass in Earth masses (~317.8).
MJUP_IN_MEARTH = M_JUP / M_EARTH
#: The deuterium-burning limit, ~13 Jupiter masses: objects heavier than this
#: briefly fuse deuterium. By convention it separates planets from brown dwarfs,
#: though nature draws no sharp line there. The archive admits objects to 30 M_J.
DEUTERIUM_LIMIT_MEARTH = 13.0 * MJUP_IN_MEARTH


# ---------------------------------------------------------------------------
# Orbits and starlight
# ---------------------------------------------------------------------------

def semi_major_axis_au(period_days, mstar_msun):
    """Kepler's third law: orbital distance from period and stellar mass.

    ``a**3 = G M P**2 / (4 pi**2)``, neglecting the planet's own mass, which is
    below 0.1% of the star's for everything in this catalogue that is a planet.
    """
    p = np.asarray(period_days, dtype=float) * DAY
    m = np.asarray(mstar_msun, dtype=float) * M_SUN
    return np.cbrt(G * m * p**2 / (4 * np.pi**2)) / AU


def period_days(a_au, mstar_msun):
    """Kepler's third law inverted: period from orbital distance."""
    a = np.asarray(a_au, dtype=float) * AU
    m = np.asarray(mstar_msun, dtype=float) * M_SUN
    return 2 * np.pi * np.sqrt(a**3 / (G * m)) / DAY


def luminosity_lsun(radius_rsun, teff_k):
    """Stefan-Boltzmann: ``L = 4 pi R**2 sigma T**4``, in solar units."""
    r = np.asarray(radius_rsun, dtype=float)
    t = np.asarray(teff_k, dtype=float)
    return r**2 * (t / T_SUN) ** 4


def insolation(luminosity_lsun_, a_au):
    """Starlight reaching a planet, relative to what Earth receives.

    ``S = L / a**2`` in solar/AU units, so Earth gets exactly 1 and a planet at
    0.1 AU from a Sun-like star gets 100. For an eccentric orbit this is the
    flux at the semi-major axis; the orbit-averaged flux is larger by
    ``1/sqrt(1 - e**2)`` -- a few percent for most planets.
    """
    return np.asarray(luminosity_lsun_, dtype=float) / np.asarray(a_au, dtype=float) ** 2


# ---------------------------------------------------------------------------
# The habitable zone: Kopparapu et al. (2013, ApJ 765, 131; 2014, ApJL 787, L29)
# ---------------------------------------------------------------------------

#: Coefficients ``(S_eff_sun, a, b, c, d)`` for the effective-flux limits, from
#: Kopparapu's own published calculator (hzcalc.py, 2014 update, for an
#: Earth-mass planet). The 2014 revision keeps the 2013 model and definitions
#: and updates the fits.
#:
#: ``S_eff = S_eff_sun + a T + b T**2 + c T**3 + d T**4`` with ``T = Teff - 5780 K``.
#:
#: - *runaway greenhouse* / *maximum greenhouse*: the conservative inner and
#:   outer edges from 1D climate models -- where oceans would boil away, and
#:   where even a thick CO2 atmosphere can no longer keep water liquid.
#: - *recent Venus* / *early Mars*: the optimistic edges, set empirically by the
#:   evidence that Venus may have had surface water until ~1 Gyr ago and early
#:   Mars ~3.8 Gyr ago, when each received this much (relative) sunlight.
KOPPARAPU_2014 = {
    "recent_venus": (1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15),
    "runaway_greenhouse": (1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15),
    "maximum_greenhouse": (0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16),
    "early_mars": (0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16),
}
#: The reference temperature in Kopparapu's polynomial (not the IAU 5772 K).
KOPPARAPU_TREF = 5780.0
KOPPARAPU_TEFF_RANGE = (2600.0, 7200.0)


def hz_flux_limit(limit: str, teff_k):
    """Effective stellar flux (Earth = 1) at one habitable-zone boundary.

    The limits depend on the star's temperature, not only on how much light
    arrives. Cooler stars emit redder light, which water and CO2 absorb more
    readily and which Rayleigh scattering reflects less, so the same total flux
    heats a planet more. The habitable zone around a red dwarf therefore sits
    at *lower* insolation than around the Sun: at 3000 K its inner edge is at
    17% less flux and its outer edge at 35% less.

    ``teff_k`` is clipped to the fits' validity range (2600-7200 K).
    """
    s0, a, b, c, d = KOPPARAPU_2014[limit]
    t = np.clip(np.asarray(teff_k, dtype=float), *KOPPARAPU_TEFF_RANGE) - KOPPARAPU_TREF
    return s0 + a * t + b * t**2 + c * t**3 + d * t**4


def hz_distance_au(limit: str, teff_k, luminosity_lsun_):
    """Orbital distance of a habitable-zone boundary: ``d = sqrt(L / S_eff)``."""
    return np.sqrt(np.asarray(luminosity_lsun_, dtype=float) / hz_flux_limit(limit, teff_k))


def hz_zone(insol, teff_k):
    """Classify planets as 'conservative', 'optimistic' or '' (outside).

    Conservative: between runaway and maximum greenhouse. Optimistic: between
    recent Venus and early Mars, excluding the conservative zone.
    """
    s = np.asarray(insol, dtype=float)
    inner_c = hz_flux_limit("runaway_greenhouse", teff_k)
    outer_c = hz_flux_limit("maximum_greenhouse", teff_k)
    inner_o = hz_flux_limit("recent_venus", teff_k)
    outer_o = hz_flux_limit("early_mars", teff_k)
    zone = np.full(s.shape, "", dtype=object)
    zone[(s <= inner_o) & (s >= outer_o)] = "optimistic"
    zone[(s <= inner_c) & (s >= outer_c)] = "conservative"
    return zone


# ---------------------------------------------------------------------------
# Why each detection method sees what it sees
# ---------------------------------------------------------------------------

def rv_semi_amplitude(mp_mearth, period_days_, mstar_msun, ecc=0.0, sin_i=1.0):
    """Radial-velocity semi-amplitude K, in m/s: how fast the star wobbles.

    ``K = (2 pi G / P)**(1/3) * Mp sin(i) / (M* + Mp)**(2/3) / sqrt(1 - e**2)``

    Jupiter moves the Sun at about 12.5 m/s -- a brisk walk. Earth moves it at
    9 cm/s, which is why an Earth twin around a Sun-like star has still not
    been found by this method.
    """
    mp = np.asarray(mp_mearth, dtype=float) * M_EARTH
    p = np.asarray(period_days_, dtype=float) * DAY
    ms = np.asarray(mstar_msun, dtype=float) * M_SUN
    e = np.asarray(ecc, dtype=float)
    return ((2 * np.pi * G / p) ** (1 / 3) * mp * sin_i
            / (ms + mp) ** (2 / 3) / np.sqrt(1 - e**2))


def rv_minimum_mass_mearth(k_ms, a_au, mstar_msun):
    """The planet mass that produces a wobble of ``k_ms`` at distance ``a``.

    For a circular orbit and a planet much lighter than its star, the star's
    speed is the planet's orbital speed scaled by the mass ratio, which gives
    ``Mp sin(i) = K * sqrt(a M* / G)``. At fixed K the detectable mass grows
    only as sqrt(a), so the signal itself fades slowly with distance; the
    harder limit on distant planets is the years of monitoring needed to
    cover a whole orbit (twelve years, for a Jupiter).
    """
    a = np.asarray(a_au, dtype=float) * AU
    ms = np.asarray(mstar_msun, dtype=float) * M_SUN
    return np.asarray(k_ms, dtype=float) * np.sqrt(a * ms / G) / M_EARTH


def transit_probability(rstar_rsun, a_au):
    """Chance that a randomly oriented circular orbit transits: ``R*/a``.

    Earth's is 0.47%: of every 215 Sun-like stars with an Earth twin, about one
    would show a transit. A hot Jupiter at 0.05 AU transits about 10% of the time.
    """
    return np.asarray(rstar_rsun, dtype=float) * R_SUN / (np.asarray(a_au, dtype=float) * AU)


def einstein_radius_au(lens_mass_msun, lens_distance_pc, source_distance_pc):
    """Physical Einstein radius of a gravitational lens, in AU.

    ``R_E = sqrt(4 G M / c**2 * D_L (D_S - D_L) / D_S)``

    A planet is most easily detected by microlensing when it lies near this
    ring, which for a typical lens -- a red dwarf halfway to the Galactic bulge
    -- is 2-3 AU. Microlensing needs no light from the planet or star, only
    mass, so it reaches cold, Earth-mass planets no other method can.
    """
    m = np.asarray(lens_mass_msun, dtype=float) * M_SUN
    dl = np.asarray(lens_distance_pc, dtype=float) * PC
    ds = np.asarray(source_distance_pc, dtype=float) * PC
    return np.sqrt(4 * G * m / C**2 * dl * (ds - dl) / ds) / AU


# ---------------------------------------------------------------------------
# Bulk properties
# ---------------------------------------------------------------------------

#: Mean density of a body with Earth's mass and nominal equatorial radius.
RHO_EARTH_UNITS_GCC = M_EARTH / (4 / 3 * np.pi * R_EARTH**3) / 1000.0


def bulk_density_gcc(mass_mearth, radius_rearth):
    """Mean density in g/cm^3 from mass and radius in Earth units."""
    return (RHO_EARTH_UNITS_GCC * np.asarray(mass_mearth, dtype=float)
            / np.asarray(radius_rearth, dtype=float) ** 3)


def valley_radius(period_days_, slope, intercept):
    """Radius-valley location for a line ``log10 R = slope * log10 P + intercept``."""
    return 10 ** (intercept + slope * np.log10(np.asarray(period_days_, dtype=float)))


# ---------------------------------------------------------------------------
# The Solar System, for scale
# ---------------------------------------------------------------------------

#: NASA Planetary Fact Sheet values. Masses in kg, equatorial radii in km (to
#: match the equatorial R_Earth and R_Jup that exoplanet radii are quoted in),
#: semi-major axes in AU, sidereal periods in days.
_SOLAR_SYSTEM_RAW = {
    "Mercury": (0.330e24, 2440.5, 0.387, 87.969),
    "Venus": (4.87e24, 6051.8, 0.723, 224.701),
    "Earth": (M_EARTH, R_EARTH / 1e3, 1.000, 365.256),
    "Mars": (0.642e24, 3396.2, 1.524, 686.980),
    "Jupiter": (M_JUP, R_JUP / 1e3, 5.203, 4332.59),
    "Saturn": (568e24, 60268.0, 9.537, 10759.22),
    "Uranus": (86.8e24, 25559.0, 19.19, 30685.4),
    "Neptune": (102e24, 24764.0, 30.07, 60189.0),
}

#: name -> dict(mass [M_e], radius [R_e], a [AU], period [d], insolation [S_e])
SOLAR_SYSTEM = {
    name: {
        "mass": m / M_EARTH,
        "radius": r_km * 1e3 / R_EARTH,
        "a": a,
        "period": p,
        "insolation": 1.0 / a**2,
    }
    for name, (m, r_km, a, p) in _SOLAR_SYSTEM_RAW.items()
}
