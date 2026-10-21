"""The formulas that turn a pulsar's spin into physics.

A pulsar is a rotating neutron star whose beamed radio emission sweeps past us
once per turn. Two numbers are measured by timing its pulses: the spin period P
(seconds) and its rate of change Pdot (seconds per second, dimensionless). Almost
everything on a P-Pdot diagram is derived from those two by assuming the star
is a spinning magnetic dipole that loses rotational energy to magnetic dipole
radiation. This module holds those derivations, one function each, with the
assumptions stated where they are made.

Units are cgs throughout -- gauss, erg, grams, centimetres -- because the pulsar
literature (and the ATNF catalogue) uses them, and mixing conventions is how
factors of 10^4 sneak into magnetic fields. Every function accepts numpy arrays.
Non-physical inputs (a zero or negative Pdot, which the catalogue does contain)
return NaN rather than a misleading number.

The fiducial neutron star is the conventional one: radius 10 km, moment of
inertia 10^45 g cm^2. Both are uncertain by tens of percent, and every derived
field and luminosity inherits that. They are fixed here because they are what
the ATNF catalogue uses, so our numbers can be checked against its own.
"""

from __future__ import annotations

import numpy as np
from astropy import constants as const
from astropy import units as u

C_CGS = const.c.cgs.value  # cm/s
G_CGS = const.G.cgs.value  # cm^3 g^-1 s^-2
SECONDS_PER_YEAR = 365.25 * 86400.0  # Julian year, as the catalogue uses
KPC_CM = u.kpc.to(u.cm)
MAS_PER_YR_IN_RAD_PER_S = (1 * u.mas / u.yr).to(u.rad / u.s).value
AGE_OF_UNIVERSE_YR = 13.8e9  # Planck 2018: 13.787 +/- 0.020 Gyr

#: Canonical neutron star. (2/5) M R^2 for 1.4 solar masses and 10 km is
#: 1.1e45 g cm^2; the round number is the convention.
MOMENT_OF_INERTIA = 1.0e45  # g cm^2
NS_RADIUS_CM = 1.0e6  # 10 km

#: Surface-field coefficient as the ATNF catalogue defines it, B = 3.2e19
#: sqrt(P Pdot) gauss. `dipole_field_coefficient()` derives it from first
#: principles; the test suite checks the two agree.
B_COEFFICIENT = 3.2e19


def _positive(x):
    """Replace zero and negative values with NaN, so logs and roots stay honest."""
    x = np.asarray(x, dtype=float)
    return np.where(x > 0, x, np.nan)


# --- Spin-down: age, field, power ---------------------------------------------


def characteristic_age_yr(period_s, pdot):
    """tau = P / (2 Pdot), in years.

    If a pulsar spins down as P Pdot = const (a constant magnetic dipole) and was
    born spinning much faster than it does now, this is its true age. Neither
    condition is guaranteed -- see `true_age_yr` for what the birth period does to
    it -- so astronomers call it the *characteristic* age.
    """
    return np.asarray(period_s, dtype=float) / (2.0 * _positive(pdot)) / SECONDS_PER_YEAR


def spin_down_age_yr(period_s, pdot, braking_index):
    """tau_n = P / ((n - 1) Pdot): the characteristic age for a general braking index."""
    return (np.asarray(period_s, dtype=float)
            / ((braking_index - 1.0) * _positive(pdot)) / SECONDS_PER_YEAR)


def true_age_yr(period_s, pdot, birth_period_s, braking_index=3.0):
    """Age of a pulsar born with ``birth_period_s``, spinning down with constant n.

    Integrating Pdot = K P^(2-n) from birth gives
        t = P / ((n-1) Pdot) * [1 - (P0 / P)^(n-1)],
    which falls to the characteristic age only when P0 << P.
    """
    ratio = np.asarray(birth_period_s, dtype=float) / np.asarray(period_s, dtype=float)
    n = braking_index
    return spin_down_age_yr(period_s, pdot, n) * (1.0 - ratio ** (n - 1.0))


def birth_period_s(period_s, pdot, age_yr, braking_index=3.0):
    """Invert `true_age_yr`: the birth period implied by a known true age.

    Returns NaN when the true age exceeds the spin-down age, which would require
    a birth period below zero -- i.e. the spin-down law assumed is inconsistent
    with the age.
    """
    n = braking_index
    frac = 1.0 - np.asarray(age_yr, dtype=float) / spin_down_age_yr(period_s, pdot, n)
    frac = np.where(frac > 0, frac, np.nan)
    return np.asarray(period_s, dtype=float) * frac ** (1.0 / (n - 1.0))


def dipole_field_coefficient(moment_of_inertia=MOMENT_OF_INERTIA, radius_cm=NS_RADIUS_CM):
    """The C in B = C sqrt(P Pdot), from equating spin-down to dipole radiation.

    A magnetic dipole m rotating at angular speed Omega in vacuum radiates
    L = 2 m^2 Omega^4 sin^2(alpha) / (3 c^3). Setting that equal to the loss of
    rotational energy, -I Omega dOmega/dt, and writing m = B R^3 for the field B at
    the magnetic equator gives
        B = sqrt(3 c^3 I P Pdot / (8 pi^2 R^6 sin^2 alpha)).
    With sin(alpha) = 1 (an orthogonal rotator) this is ~3.2e19 sqrt(P Pdot) G.
    The field at the magnetic *pole* is twice this; some papers quote that.
    """
    return np.sqrt(3.0 * C_CGS**3 * moment_of_inertia / (8.0 * np.pi**2 * radius_cm**6))


def surface_field_gauss(period_s, pdot):
    """Equatorial surface dipole field, B = 3.2e19 sqrt(P Pdot) gauss."""
    return B_COEFFICIENT * np.sqrt(np.asarray(period_s, dtype=float) * _positive(pdot))


def spin_down_luminosity(period_s, pdot, moment_of_inertia=MOMENT_OF_INERTIA):
    """Edot = 4 pi^2 I Pdot / P^3, in erg/s: the rate rotational energy is lost.

    Rotational energy is E = (1/2) I Omega^2 = 2 pi^2 I / P^2; differentiate.
    This is the power budget for everything the pulsar does -- radio pulses are a
    tiny fraction of it; most goes into a relativistic wind.
    """
    p = np.asarray(period_s, dtype=float)
    return 4.0 * np.pi**2 * moment_of_inertia * _positive(pdot) / p**3


def rotational_energy_erg(period_s, moment_of_inertia=MOMENT_OF_INERTIA):
    """E = 2 pi^2 I / P^2."""
    return 2.0 * np.pi**2 * moment_of_inertia / np.asarray(period_s, dtype=float) ** 2


def light_cylinder_radius_cm(period_s):
    """R_LC = c P / (2 pi): where co-rotating field lines would have to move at c."""
    return C_CGS * np.asarray(period_s, dtype=float) / (2.0 * np.pi)


def braking_index(f0, f1, f2):
    """n = f0 f2 / f1^2, from spin frequency and its first two derivatives.

    For spin-down torque proportional to Omega^n, a pure magnetic dipole gives n = 3.
    Measuring f2 is hard: for most pulsars it is swamped by timing noise and
    glitch recoveries, and only young, fast-slowing pulsars give a physical value.
    """
    f1 = np.asarray(f1, dtype=float)
    return np.asarray(f0, dtype=float) * np.asarray(f2, dtype=float) / f1**2


# --- Guide lines for the P-Pdot diagram -----------------------------------------


def pdot_for_age(period_s, age_yr):
    """Pdot along a line of constant characteristic age (slope +1 in log-log)."""
    return np.asarray(period_s, dtype=float) / (2.0 * age_yr * SECONDS_PER_YEAR)


def pdot_for_field(period_s, field_gauss):
    """Pdot along a line of constant surface field (slope -1 in log-log)."""
    return (field_gauss / B_COEFFICIENT) ** 2 / np.asarray(period_s, dtype=float)


def pdot_for_edot(period_s, edot, moment_of_inertia=MOMENT_OF_INERTIA):
    """Pdot along a line of constant spin-down luminosity (slope +3 in log-log)."""
    return edot * np.asarray(period_s, dtype=float) ** 3 / (4.0 * np.pi**2 * moment_of_inertia)


#: Death lines from Zhang, Harding & Muslimov (2000, ApJ 531, L135), their eqs.
#: (3) and (4): the vacuum-gap model of Ruderman & Sutherland (1975) with
#: curvature-radiation pair production and frame dragging included.
#:   'dipole'     log Pdot = (11/4) log P - 14.62   pure dipole surface field
#:   'multipole'  log Pdot = (9/4)  log P - 16.58   field lines curved on the
#:                                                  stellar-radius scale (rho_6 = 1)
#: Chen & Ruderman (1993) called the band between the two the "death valley".
DEATH_LINES = {
    "dipole": (11.0 / 4.0, -14.62),
    "multipole": (9.0 / 4.0, -16.58),
}


def death_line_pdot(period_s, model="dipole"):
    """Pdot below which the named model predicts radio emission switches off.

    Physically: a pulsar shines in radio only while the voltage its rotation
    induces across the polar cap can accelerate particles hard enough to make
    electron-positron pairs. Slower spin or a weaker field means less voltage; the
    line is where pair production stops. The exact position depends on the
    unknown field geometry near the surface, which is why there are two lines.
    """
    slope, intercept = DEATH_LINES[model]
    return 10.0 ** (slope * np.log10(np.asarray(period_s, dtype=float)) + intercept)


# --- Kinematic contamination of Pdot -------------------------------------------


def shklovskii_pdot(period_s, proper_motion_mas_yr, distance_kpc):
    """Apparent Pdot from transverse motion alone: P mu^2 d / c (Shklovskii 1970).

    A pulsar moving across the sky is getting farther away at an accelerating
    rate, so successive pulses have ever longer to travel. That looks exactly like
    spin-down. It is negligible for ordinary pulsars and can dominate for
    millisecond pulsars, whose intrinsic Pdot is ~10^5 times smaller.
    """
    mu = np.asarray(proper_motion_mas_yr, dtype=float) * MAS_PER_YR_IN_RAD_PER_S
    d = np.asarray(distance_kpc, dtype=float) * KPC_CM
    return np.asarray(period_s, dtype=float) * mu**2 * d / C_CGS


# --- What spin alone proves ------------------------------------------------------


def minimum_mean_density(period_s):
    """rho > 3 pi / (G P^2), in g/cm^3.

    A self-gravitating body rotating once every P holds together only if gravity
    at its equator exceeds the centripetal requirement. For a sphere of mean
    density rho that is exactly the condition that a surface orbit is not faster
    than the spin. This is the argument that ruled out white dwarfs in 1968: they
    cannot be dense enough to spin in milliseconds.
    """
    return 3.0 * np.pi / (G_CGS * np.asarray(period_s, dtype=float) ** 2)


def max_spin_frequency_hz(mass_msun, radius_km):
    """Mass-shedding limit, nu_max = 1045 (M/Msun)^(1/2) (10 km / R)^(3/2) Hz.

    Lattimer & Prakash (2004), fitted across many equations of state, where R is
    the radius of the non-rotating star of the same mass. It sits below the
    naive Kepler frequency because a fast-spinning star bulges at the equator.
    """
    return 1045.0 * np.sqrt(mass_msun) * (10.0 / np.asarray(radius_km, dtype=float)) ** 1.5


def max_radius_km(spin_hz, mass_msun):
    """Largest non-rotating radius that could survive spinning at ``spin_hz``."""
    return 10.0 * (1045.0 * np.sqrt(mass_msun) / np.asarray(spin_hz, dtype=float)) ** (2.0 / 3.0)


def equatorial_speed_km_s(period_s, radius_km=10.0):
    """Speed of the equator, 2 pi R / P."""
    return 2.0 * np.pi * radius_km / np.asarray(period_s, dtype=float)
