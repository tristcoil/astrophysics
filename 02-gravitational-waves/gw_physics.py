"""The physics: what a chirp's frequency sweep says about the objects making it.

Two orbiting masses radiate gravitational waves at twice their orbital
frequency. The waves carry away orbital energy, so the orbit shrinks, the
orbital frequency rises, and the waves rise with it -- the "chirp". To leading
(Newtonian) order the rate of that rise depends on the masses through a single
combination, the **chirp mass**::

    Mc = (m1 * m2) ** (3/5) / (m1 + m2) ** (1/5)

and obeys::

    df/dt = (96/5) * pi**(8/3) * (G*Mc/c**3)**(5/3) * f**(11/3)

That one equation is how the discovery paper went from a squiggle to "these are
black holes": measure f and df/dt, solve for Mc, bound the total mass, and ask
what could orbit that fast without touching.

Masses here are in solar masses and frequencies in hertz unless a function says
otherwise. Constants come from ``astropy.constants``. The convenient unit
``T_SUN = G*Msun/c**3`` is the Sun's mass expressed as a time, 4.93
microseconds, which is how mass enters every formula below.

The second half of this module is detector geometry: where the two LIGO sites
are, how their arms point, and what that implies for the arrival-time
difference and relative sign of a signal seen in both.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import constants as const
from scipy.optimize import brentq

G = const.G.value
C = const.c.value
M_SUN = const.M_sun.value
R_SUN = const.R_sun.value
#: The Sun's mass as a time, G*Msun/c^3 ~ 4.93 microseconds.
T_SUN = G * M_SUN / C**3
#: The Sun's mass as a length, G*Msun/c^2 ~ 1.48 km.
L_SUN = G * M_SUN / C**2


# --------------------------------------------------------------------------
# Masses
# --------------------------------------------------------------------------

def chirp_mass(m1, m2):
    """Chirp mass of a binary (same units as the inputs)."""
    m1, m2 = np.asarray(m1, float), np.asarray(m2, float)
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def symmetric_mass_ratio(m1, m2):
    """eta = m1*m2/(m1+m2)^2; 1/4 for equal masses, smaller otherwise."""
    m1, m2 = np.asarray(m1, float), np.asarray(m2, float)
    return m1 * m2 / (m1 + m2) ** 2


def min_total_mass(mc):
    """Smallest total mass consistent with a chirp mass.

    ``M = Mc * eta**(-3/5)`` and eta is at most 1/4 (equal masses), so the total
    mass can never be below ``2**(6/5) * Mc`` ~ 2.3 Mc. Unequal masses only make
    it larger. This lower bound is the whole of the discovery paper's "total
    mass >~ 70 Msun" -- it needs no knowledge of the individual masses.
    """
    return 2.0 ** 1.2 * np.asarray(mc, float)


def partner_mass(mc, m_known):
    """Mass of the other body, given a chirp mass and one component mass.

    Used for the discovery paper's black hole-neutron star argument: pair a
    neutron star with whatever companion reproduces the measured chirp mass, and
    see what total mass -- and so what merger frequency -- that implies.
    """
    if chirp_mass(m_known, 1e9) < mc:
        raise ValueError(f"no companion reaches Mc = {mc} with m = {m_known}")
    return brentq(lambda m: chirp_mass(m_known, m) - mc, 1e-6, 1e9, xtol=1e-10)


# --------------------------------------------------------------------------
# The Newtonian chirp
# --------------------------------------------------------------------------

def fdot_newtonian(f, mc):
    """Rate of change of GW frequency (Hz/s) at frequency ``f`` for chirp mass ``mc``."""
    tc = np.asarray(mc, float) * T_SUN
    return 96.0 / 5.0 * np.pi ** (8.0 / 3.0) * tc ** (5.0 / 3.0) * np.asarray(f, float) ** (11.0 / 3.0)


def chirp_mass_from_fdot(f, fdot):
    """Invert ``fdot_newtonian``: the discovery paper's formula, in solar masses.

    Mc = (c^3/G) * [ (5/96) * pi^(-8/3) * f^(-11/3) * fdot ]^(3/5)
    """
    x = 5.0 / 96.0 * np.pi ** (-8.0 / 3.0) * np.asarray(f, float) ** (-11.0 / 3.0) * np.asarray(fdot, float)
    return x ** 0.6 / T_SUN


def time_to_coalescence(f, mc):
    """Seconds left until merger, to Newtonian order, from GW frequency ``f``.

    Integrating ``fdot_newtonian`` gives ``tau = (5/256) (G Mc/c^3)^(-5/3) (pi f)^(-8/3)``.
    """
    tc = np.asarray(mc, float) * T_SUN
    return 5.0 / 256.0 * tc ** (-5.0 / 3.0) * (np.pi * np.asarray(f, float)) ** (-8.0 / 3.0)


def frequency_at(tau, mc):
    """GW frequency a time ``tau`` (s) before coalescence -- inverse of the above."""
    tc = np.asarray(mc, float) * T_SUN
    tau = np.asarray(tau, float)
    return (5.0 / 256.0 / tau) ** 0.375 * tc ** (-0.625) / np.pi


def chirp_slope(mc):
    """Slope of f^(-8/3) against time, in Hz^(-8/3) per second.

    Integrating the chirp law shows that ``f**(-8/3)`` falls *linearly* in time,
    reaching zero at coalescence::

        f^(-8/3) = (256/5) * pi^(8/3) * (G Mc/c^3)^(5/3) * (t_c - t)

    So plotting f^(-8/3) against t turns the chirp into a straight line whose
    slope is set by the chirp mass alone. That is the fit this project makes.
    """
    tc = np.asarray(mc, float) * T_SUN
    return -256.0 / 5.0 * np.pi ** (8.0 / 3.0) * tc ** (5.0 / 3.0)


def chirp_mass_from_slope(slope):
    """Invert ``chirp_slope``: chirp mass (Msun) from a fitted f^(-8/3)-vs-t slope."""
    slope = np.asarray(slope, float)
    if np.any(slope >= 0):
        raise ValueError("f^(-8/3) must fall with time; slope must be negative")
    return (-slope * 5.0 / 256.0 * np.pi ** (-8.0 / 3.0)) ** 0.6 / T_SUN


# --------------------------------------------------------------------------
# Orbits and compactness
# --------------------------------------------------------------------------

def keplerian_separation(f_gw, m_total):
    """Separation (m) of two point masses radiating at GW frequency ``f_gw``.

    The orbital frequency is half the GW frequency (the mass distribution looks
    the same after half an orbit), so the angular frequency is ``pi * f_gw``
    and Kepler's third law gives ``r^3 = G M / (pi f_gw)^2``.
    """
    gm = G * np.asarray(m_total, float) * M_SUN
    return (gm / (np.pi * np.asarray(f_gw, float)) ** 2) ** (1.0 / 3.0)


def orbital_speed(f_gw, m_total):
    """Relative orbital speed as a fraction of c: v/c = (pi G M f / c^3)^(1/3).

    This is the post-Newtonian expansion parameter. Once it approaches ~0.5 the
    leading-order chirp formula is no longer a good description.
    """
    return (np.pi * np.asarray(m_total, float) * T_SUN * np.asarray(f_gw, float)) ** (1.0 / 3.0)


def schwarzschild_radius(m):
    """Event-horizon radius 2GM/c^2 of a non-spinning black hole, in metres."""
    return 2.0 * L_SUN * np.asarray(m, float)


def isco_frequency(m_total):
    """GW frequency at the innermost stable circular orbit of a test mass.

    ``f = c^3 / (6^(3/2) pi G M)``, about 4400 Hz / (M/Msun). A rough marker of
    where inspiral gives way to plunge and merger.
    """
    return 1.0 / (6.0 ** 1.5 * np.pi * np.asarray(m_total, float) * T_SUN)


def contact_frequency(m_total, r_sum):
    """Highest GW frequency two bodies can reach before their surfaces touch.

    Kepler's law at separation ``r_sum`` (the sum of the two radii, metres),
    converted to GW frequency. Newtonian, so only an order-of-magnitude guide for
    black holes -- but orders of magnitude are all the compactness argument
    needs.
    """
    gm = G * np.asarray(m_total, float) * M_SUN
    return np.sqrt(gm / np.asarray(r_sum, float) ** 3) / np.pi


# Mass-radius relations for the compactness argument. Each gives the radius (m)
# of one body of mass m (Msun). They are deliberately simple: the argument turns
# on differences of many orders of magnitude, not on percent-level structure.

def radius_main_sequence(m):
    """Main-sequence star: R ~ M^0.8 below a solar mass, ~ M^0.57 above."""
    m = np.asarray(m, float)
    return R_SUN * np.where(m < 1.0, m ** 0.8, m ** 0.57)


#: Chandrasekhar mass for a carbon-oxygen white dwarf (mu_e = 2).
M_CHANDRASEKHAR = 1.44


def radius_white_dwarf(m):
    """Nauenberg (1972) relation: white dwarfs *shrink* as mass is added.

    R = 0.0112 Rsun * sqrt((M/Mch)^(-2/3) - (M/Mch)^(2/3)), going to zero at the
    Chandrasekhar mass, beyond which no white dwarf exists.
    """
    x = np.asarray(m, float) / M_CHANDRASEKHAR
    return 0.0112 * R_SUN * np.sqrt(np.clip(x ** (-2.0 / 3.0) - x ** (2.0 / 3.0), 0.0, None))


#: Neutron stars are ~12 km whatever their mass (NICER X-ray timing gives
#: 12-13 km); the heaviest known are ~2.1 Msun and theory caps them near 2.3.
R_NEUTRON_STAR = 12e3
M_NEUTRON_STAR_MAX = 2.3


# --------------------------------------------------------------------------
# Energy
# --------------------------------------------------------------------------

#: One solar mass of energy, E = M c^2, in joules (~1.8e47 J).
MSUN_C2_JOULES = M_SUN * C**2


def erg_per_s_to_msun_c2_per_s(luminosity_erg_s):
    """Convert a luminosity in erg/s to solar masses of energy per second."""
    return np.asarray(luminosity_erg_s, float) * 1e-7 / MSUN_C2_JOULES


# --------------------------------------------------------------------------
# Detector geometry
# --------------------------------------------------------------------------

#: WGS-84 reference ellipsoid, the Earth model LIGO's coordinates are quoted on.
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563


@dataclass(frozen=True)
class Detector:
    """An L-shaped interferometer: where its corner is and where its arms point.

    Angles follow LALSuite's convention: azimuth measured clockwise from local
    North, altitude measured up from the local horizontal (the arms are 4 km
    long and straight, while the Earth curves away beneath them, so each arm
    tilts very slightly relative to the ground at the corner).
    """

    name: str
    latitude: float    # geodetic, radians
    longitude: float   # radians, east positive
    elevation: float   # metres above the ellipsoid
    x_azimuth: float
    y_azimuth: float
    x_altitude: float
    y_altitude: float

    def _local_axes(self):
        """East, North and Up unit vectors at the vertex, in Earth-fixed axes."""
        lat, lon = self.latitude, self.longitude
        east = np.array([-np.sin(lon), np.cos(lon), 0.0])
        north = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
        up = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
        return east, north, up

    def arm(self, which: str) -> np.ndarray:
        """Unit vector along the x or y arm, in Earth-centred Earth-fixed axes."""
        az = self.x_azimuth if which == "x" else self.y_azimuth
        alt = self.x_altitude if which == "x" else self.y_altitude
        east, north, up = self._local_axes()
        return np.cos(alt) * (np.sin(az) * east + np.cos(az) * north) + np.sin(alt) * up

    @property
    def vertex(self) -> np.ndarray:
        """Earth-centred Earth-fixed position of the corner station, metres."""
        e2 = WGS84_F * (2.0 - WGS84_F)
        lat, lon, h = self.latitude, self.longitude, self.elevation
        n = WGS84_A / np.sqrt(1.0 - e2 * np.sin(lat) ** 2)
        return np.array([
            (n + h) * np.cos(lat) * np.cos(lon),
            (n + h) * np.cos(lat) * np.sin(lon),
            (n * (1.0 - e2) + h) * np.sin(lat),
        ])

    @property
    def tensor(self) -> np.ndarray:
        """Detector response tensor D = (x x - y y) / 2.

        The detector's output is the contraction of this tensor with the wave's
        strain tensor: it measures how much the wave stretches one arm relative
        to the other. Rotating an L-shape by 90 degrees swaps x and y, which
        flips D's sign -- so two detectors rotated 90 degrees apart see the same
        wave with opposite signs.
        """
        x, y = self.arm("x"), self.arm("y")
        return 0.5 * (np.outer(x, x) - np.outer(y, y))


#: LIGO Hanford and Livingston, constants from LALSuite's LALDetectors.h
#: (``LAL_LHO_4K_*`` and ``LAL_LLO_4K_*``), the reference every LIGO analysis uses.
DETECTORS = {
    "H1": Detector(
        name="LIGO Hanford", latitude=0.81079526383, longitude=-2.08405676917,
        elevation=142.554, x_azimuth=5.65487724844, y_azimuth=4.08408092164,
        x_altitude=-0.00061950000, y_altitude=0.00001250000,
    ),
    "L1": Detector(
        name="LIGO Livingston", latitude=0.53342313506, longitude=-1.58430937078,
        elevation=-6.574, x_azimuth=4.40317772346, y_azimuth=2.83238139666,
        x_altitude=-0.00031210000, y_altitude=-0.00061070000,
    ),
}


def light_travel_time(a: Detector, b: Detector) -> float:
    """Straight-line light travel time between two detector vertices, seconds.

    A gravitational wave travels at c, so no real signal can arrive at the two
    sites further apart in time than this. It is the first test any candidate
    must pass.
    """
    return float(np.linalg.norm(a.vertex - b.vertex) / C)


def tensor_overlap(a: Detector, b: Detector) -> float:
    """Normalised overlap of two detector tensors, between -1 and +1.

    +1 means identical orientation, -1 means rotated 90 degrees in the same
    plane: every wave then produces equal and opposite outputs. The two LIGO
    sites were deliberately built nearly anti-aligned (as far as the Earth's
    curvature between them allows) so that they respond to the same polarisation.
    """
    da, db = a.tensor, b.tensor
    return float(np.sum(da * db) / np.sqrt(np.sum(da * da) * np.sum(db * db)))


def source_direction(ra, dec, gmst):
    """Unit vector from Earth's centre toward a source, in Earth-fixed axes.

    ``gmst`` is Greenwich mean sidereal time in radians. The Earth-fixed
    longitude of the point beneath the source is ``ra - gmst``.
    """
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    lon = ra - gmst
    return np.stack([np.cos(dec) * np.cos(lon), np.cos(dec) * np.sin(lon), np.sin(dec)], axis=-1)


def arrival_delay(a: Detector, b: Detector, ra, dec, gmst):
    """Arrival time at detector ``a`` minus arrival time at ``b``, seconds.

    A plane wave from direction n reaches a detector at position r at
    ``t0 - (n . r)/c``: the detector nearer the source hears it first.
    """
    n = source_direction(ra, dec, gmst)
    return (n @ (b.vertex - a.vertex)) / C
