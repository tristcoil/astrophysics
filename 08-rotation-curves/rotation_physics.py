"""Rotation-curve physics: from orbital speeds and visible matter to missing mass.

Everything here works in the units galaxy dynamicists use -- speeds in km/s,
radii in kpc, masses in solar masses -- because those are the units SPARC
publishes, and because in them the numbers stay human-sized. The constants are
derived once from ``astropy.constants`` (CODATA 2018 / IAU 2015) rather than
hand-copied, and ``test_physics.py`` checks every conversion against a case with
a known answer.

The core idea fits in one line. Something on a circular orbit of radius ``r``
at speed ``V`` needs a centripetal acceleration ``V**2 / r``, and gravity has to
supply it. For a spherical distribution of mass that gives::

    V**2 / r = G * M(<r) / r**2      so      V = sqrt(G * M(<r) / r)

Once you are outside most of the mass, ``M(<r)`` stops growing and ``V`` falls
as ``1/sqrt(r)``: the Keplerian decline the planets follow. A *flat* rotation
curve, ``V = constant``, instead needs ``M(<r)`` proportional to ``r`` -- mass
that keeps accumulating, linearly, far beyond where the light runs out.

Galaxies are discs, not spheres, so the real calculation solves Poisson's
equation for the observed flattened distributions of stars and gas. SPARC has
already done that for every galaxy (see ``baryonic_velocity_squared``); the
spherical formulas here are used for the Solar System, for dark matter haloes
(which are round to a good approximation), and for enclosed-mass estimates that
are labelled as such wherever they appear.
"""

from __future__ import annotations

import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy import special

# --------------------------------------------------------------------------
# Constants, in galaxy units
# --------------------------------------------------------------------------

#: Newton's G in kpc (km/s)^2 / Msun. With it, V**2 = G*M/r needs no further
#: conversion when V is in km/s, r in kpc and M in solar masses.
G_KPC = float((const.G * u.M_sun / u.kpc).to(u.km**2 / u.s**2).value)

#: One kiloparsec in metres, for converting V**2/r into SI accelerations.
KPC_M = float(u.kpc.to(u.m))

#: One astronomical unit in kpc, for putting the Solar System on the same
#: footing as a galaxy.
AU_KPC = float(u.au.to(u.kpc))

#: GM of the Sun in km^3/s^2. The product GM is known far more precisely than
#: G or M separately, which is why planetary dynamics is done with it.
GM_SUN_KM3_S2 = float((const.GM_sun).to(u.km**3 / u.s**2).value)

#: Mass of gas per unit mass of atomic hydrogen. 21-cm radio observations
#: measure only neutral hydrogen, but primordial gas is ~25% helium by mass, so
#: the total is 1/0.75 = 1.33 times the hydrogen. SPARC's Vgas column already
#: includes this factor; we need it again only to turn catalogue HI masses into
#: total gas masses.
HELIUM_FACTOR = 1.33


# --------------------------------------------------------------------------
# Adding up the visible matter
# --------------------------------------------------------------------------

def signed_square(v):
    """``v * |v|``: a velocity squared that remembers its sign.

    Contributions to a rotation curve add in quadrature -- it is ``V**2`` that is
    proportional to the gravitational pull, not ``V``. SPARC reports each
    component as a velocity, but some are *negative*: a gas disc with a central
    hole, for instance, has most of its mass outside small radii, and in a disc
    (unlike a sphere) mass outside you pulls you *outward*. The signed square
    keeps that outward pull subtracting from the total instead of adding to it.
    """
    v = np.asarray(v, dtype=float)
    return v * np.abs(v)


def signed_sqrt(v2):
    """Inverse of ``signed_square``: the velocity a signed ``V**2`` corresponds to."""
    v2 = np.asarray(v2, dtype=float)
    return np.sign(v2) * np.sqrt(np.abs(v2))


def baryonic_velocity_squared(v_gas, v_disk, v_bulge, upsilon_disk, upsilon_bulge):
    """Signed ``V_bar**2`` (km^2/s^2): everything the visible matter can pull.

    ``V_bar**2 = Vgas|Vgas| + Y_disk * Vdisk|Vdisk| + Y_bulge * Vbul|Vbul|``

    SPARC computes each component by solving Poisson's equation for the observed
    surface density of that component -- 21-cm maps for the gas, 3.6-micron
    Spitzer images for the stars -- and publishes the stellar ones for a
    mass-to-light ratio of exactly 1 Msun/Lsun. Because ``V**2`` is proportional
    to mass, a different mass-to-light ratio ``Y`` simply multiplies the stellar
    ``V**2`` terms by ``Y``. That one number is the only freedom in the visible
    matter's prediction, which is what makes the test sharp.
    """
    return (
        signed_square(v_gas)
        + upsilon_disk * signed_square(v_disk)
        + upsilon_bulge * signed_square(v_bulge)
    )


def baryonic_mass(l_disk, l_bulge, m_hi, upsilon_disk, upsilon_bulge):
    """Total visible (baryonic) mass in solar masses.

    ``M_bar = Y_disk * L_disk + Y_bulge * L_bulge + 1.33 * M_HI``

    Luminosities in Lsun (3.6 micron), HI mass in Msun. Molecular gas is not
    included: it is not in SPARC, and in these galaxies it is usually a small
    fraction of the atomic gas -- one of the reasons the visible mass here is a
    slight underestimate, never an overestimate.
    """
    return upsilon_disk * l_disk + upsilon_bulge * l_bulge + HELIUM_FACTOR * m_hi


# --------------------------------------------------------------------------
# Speeds, accelerations and masses
# --------------------------------------------------------------------------

def acceleration(v_kms, r_kpc):
    """Centripetal acceleration ``V**2 / r`` in m/s^2.

    For galaxies this comes out around 1e-10 m/s^2 -- about a hundred billionth
    of the acceleration you feel standing on Earth. That tininess matters: it is
    precisely at accelerations below ~1e-10 m/s^2 that the visible matter stops
    accounting for the observed motions (see ``rar_mcgaugh``).
    """
    v = np.asarray(v_kms, dtype=float) * 1e3
    r = np.asarray(r_kpc, dtype=float) * KPC_M
    return v**2 / r


def enclosed_mass(v_kms, r_kpc):
    """Mass needed inside radius ``r`` to hold an orbit at speed ``V``, in Msun.

    ``M(<r) = V**2 * r / G`` -- exact for a spherical mass distribution. For a
    flattened disc the true enclosed mass is somewhat smaller (a disc pulls
    harder in its own plane than a sphere of the same mass), so treat this as a
    spherical-equivalent mass. Comparing two such numbers at the same radius, as
    we do for observed versus visible mass, is much less sensitive to the
    geometry than either number alone.
    """
    v = np.asarray(v_kms, dtype=float)
    return v**2 * np.asarray(r_kpc, dtype=float) / G_KPC


def keplerian_velocity(mass_msun, r_kpc):
    """Circular speed around a point mass (or outside a sphere): ``sqrt(G M / r)``."""
    return np.sqrt(G_KPC * np.asarray(mass_msun, dtype=float) / np.asarray(r_kpc, dtype=float))


def mean_orbital_speed(a_km, eccentricity, gm_km3_s2=GM_SUN_KM3_S2):
    """Time-averaged speed of a body on an elliptical orbit, in km/s.

    The average speed is the orbit's circumference divided by its period. An
    ellipse of semi-major axis ``a`` has circumference ``4 a E(e)``, with ``E``
    the complete elliptic integral of the second kind, and Kepler's third law
    gives the period, so::

        <v> = (2/pi) * E(e) * sqrt(GM / a)

    For a circle ``E(0) = pi/2`` and this is the familiar ``sqrt(GM/a)``. For
    Mercury (e = 0.21) the average is about 1% below the circular value -- which
    is exactly the size of the difference between the NASA fact-sheet speed and
    the naive formula, and why ``test_physics.py`` uses this version.
    """
    e = np.asarray(eccentricity, dtype=float)
    # scipy's ellipe takes the parameter m = k**2, and for an ellipse k = e.
    return (2.0 / np.pi) * special.ellipe(e**2) * np.sqrt(gm_km3_s2 / np.asarray(a_km, dtype=float))


# --------------------------------------------------------------------------
# An idealised disc, for checking that the published components mean what we
# think they mean
# --------------------------------------------------------------------------

def freeman_disk_velocity(r_kpc, mass_msun, r_disk_kpc):
    """Rotation curve of a razor-thin exponential disc (Freeman 1970), in km/s.

    Surface density ``Sigma(R) = Sigma0 * exp(-R / R_d)`` with total mass
    ``M = 2 pi Sigma0 R_d**2`` has::

        V**2 = (2 G M / R_d) * y**2 * [I0(y) K0(y) - I1(y) K1(y)],   y = R / (2 R_d)

    where ``I`` and ``K`` are modified Bessel functions. The curve rises, peaks
    near ``R = 2.15 R_d``, then falls -- approaching the Keplerian curve of the
    same total mass from slightly above, because a flattened mass pulls harder
    in its own plane than the same mass concentrated at a point. Real galaxy
    discs are only roughly exponential, so this is a cross-check on units and
    normalisation, not a substitute for SPARC's numerical mass models.
    """
    r = np.asarray(r_kpc, dtype=float)
    y = r / (2.0 * r_disk_kpc)
    # The exponentially scaled Bessel functions keep I*K finite at large y:
    # i0e(y) * k0e(y) = I0(y) e^-y * K0(y) e^y = I0(y) K0(y).
    bessel = special.i0e(y) * special.k0e(y) - special.i1e(y) * special.k1e(y)
    v2 = 2.0 * G_KPC * mass_msun / r_disk_kpc * y**2 * bessel
    return np.sqrt(np.clip(v2, 0.0, None))


# --------------------------------------------------------------------------
# Dark matter haloes
# --------------------------------------------------------------------------

def r200_kpc(v200_kms, h0_kms_mpc):
    """Radius at which a halo's mean density is 200 times the critical density.

    For such a halo ``M200 = 100 H0**2 R200**3 / G`` and ``V200**2 = G M200 / R200``,
    which combine to the neat ``V200 = 10 * H0 * R200``. With V200 in km/s and
    H0 in km/s/Mpc that gives R200 in Mpc; we return kpc.
    """
    return np.asarray(v200_kms, dtype=float) / (10.0 * h0_kms_mpc) * 1e3


def m200_msun(v200_kms, h0_kms_mpc):
    """Halo mass inside R200, in Msun: ``M200 = V200**2 * R200 / G``."""
    v200 = np.asarray(v200_kms, dtype=float)
    return v200**2 * r200_kpc(v200, h0_kms_mpc) / G_KPC


def v_nfw(r_kpc, v200_kms, concentration, h0_kms_mpc):
    """Circular speed of an NFW halo (Navarro, Frenk & White 1996, 1997), km/s.

    NFW is the density profile that cold-dark-matter simulations produce::

        rho(r) = rho_s / [(r/r_s) * (1 + r/r_s)**2]

    It is steep ("cuspy", rho ~ 1/r) at the centre and falls as 1/r**3 far out.
    Parameterised by ``V200`` and the concentration ``c = R200 / r_s``, the
    rotation curve is::

        V**2 = V200**2 * [ln(1 + c x) - c x / (1 + c x)] / (x * [ln(1 + c) - c / (1 + c)])

    with ``x = r / R200``. It peaks at ``r = 2.163 r_s``.
    """
    r = np.asarray(r_kpc, dtype=float)
    x = r / r200_kpc(v200_kms, h0_kms_mpc)
    cx = concentration * x
    mu = lambda t: np.log1p(t) - t / (1.0 + t)  # noqa: E731 -- the NFW mass function
    with np.errstate(divide="ignore", invalid="ignore"):
        v2 = v200_kms**2 * mu(cx) / (x * mu(concentration))
    return np.sqrt(np.where(r > 0, v2, 0.0))


#: Hubble parameter h of the Planck cosmology Dutton & Maccio (2014) simulated.
#: Their halo masses are quoted in units of 1/h solar masses.
DUTTON_MACCIO_H = 0.671


def nfw_concentration_lcdm(m200_msun):
    """Median NFW concentration that cold-dark-matter simulations predict, z = 0.

    Dutton & Maccio (2014, MNRAS 441, 3359), eq. 8, for relaxed haloes in a
    Planck cosmology::

        log10 c200 = 0.905 - 0.101 * log10(M200 / [1e12 Msun/h])

    with an intrinsic scatter of 0.11 dex about it. Lighter haloes form earlier,
    when the universe was denser, so they are *more* concentrated -- which is
    why dwarf galaxies are the sharpest test of the predicted central cusp.
    """
    m_h = np.asarray(m200_msun, dtype=float) * DUTTON_MACCIO_H
    return 10.0 ** (0.905 - 0.101 * np.log10(m_h / 1e12))


def v_iso(r_kpc, v_inf_kms, r_core_kpc):
    """Circular speed of a pseudo-isothermal halo, km/s.

    The classic empirical halo of the rotation-curve literature::

        rho(r) = rho_0 / (1 + (r / r_c)**2)

    has a constant-density *core* inside ``r_c`` and falls as 1/r**2 outside --
    the density law that by itself produces a flat rotation curve. Its curve
    rises and flattens to ``V_inf``::

        V**2 = V_inf**2 * [1 - (r_c / r) * arctan(r / r_c)],   V_inf**2 = 4 pi G rho_0 r_c**2
    """
    r = np.asarray(r_kpc, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        v2 = v_inf_kms**2 * (1.0 - (r_core_kpc / r) * np.arctan(r / r_core_kpc))
    return np.sqrt(np.where(r > 0, np.clip(v2, 0.0, None), 0.0))


def iso_central_density(v_inf_kms, r_core_kpc):
    """Central density of a pseudo-isothermal halo in Msun/pc^3."""
    rho_kpc3 = v_inf_kms**2 / (4.0 * np.pi * G_KPC * r_core_kpc**2)
    return rho_kpc3 / 1e9


# --------------------------------------------------------------------------
# The radial acceleration relation
# --------------------------------------------------------------------------

def rar_mcgaugh(g_bar, g_dagger):
    """Observed acceleration predicted from the baryonic one (McGaugh+2016, eq. 4).

    ``g_obs = g_bar / (1 - exp(-sqrt(g_bar / g_dagger)))``

    A one-parameter description of the radial acceleration relation. Where
    ``g_bar >> g_dagger`` it gives ``g_obs = g_bar``: visible matter explains
    everything and no dark matter is needed. Where ``g_bar << g_dagger`` it gives
    ``g_obs = sqrt(g_bar * g_dagger)``: the observed pull exceeds the visible
    matter's, by more the weaker the pull gets. ``g_dagger`` is the acceleration
    scale where the switch happens.

    The function is a *description* of the data, not a theory. It happens to
    have the form MOND predicts, which is why it sits at the centre of the
    dark-matter-versus-modified-gravity debate.
    """
    g_bar = np.asarray(g_bar, dtype=float)
    x = np.sqrt(g_bar / g_dagger)
    # -expm1(-x) is 1 - exp(-x) without the catastrophic cancellation at small x.
    return g_bar / -np.expm1(-x)


def rar_log_slope(g_bar, g_dagger):
    """``d log g_obs / d log g_bar`` along ``rar_mcgaugh``.

    With ``u = sqrt(g_bar / g_dagger)`` this is ``1 - (u/2) / (e**u - 1)``: 1 in
    the high-acceleration (Newtonian, no dark matter) regime and 1/2 in the
    low-acceleration one. A fit with errors on both axes needs it to turn an
    uncertainty in ``g_bar`` into an equivalent uncertainty in ``g_obs``.
    """
    u = np.sqrt(np.asarray(g_bar, dtype=float) / g_dagger)
    return 1.0 - 0.5 * u / np.expm1(u)


# --------------------------------------------------------------------------
# Error propagation
# --------------------------------------------------------------------------

def log_gobs_error(v_kms, e_v_kms, inclination_deg, e_inclination_deg,
                   distance_mpc, e_distance_mpc):
    """Uncertainty in log10(g_obs), from velocity, inclination and distance.

    ``g_obs = V**2 / R``, and three things move it:

    * the velocity measurement itself: ``d ln g = 2 dV / V``;
    * the disc inclination ``i``. We see only the line-of-sight part of the
      rotation, ``V sin i``, so the deprojected ``V`` scales as ``1/sin i`` and
      an inclination error gives ``d ln V = di / tan i`` -- small for edge-on
      discs, large for nearly face-on ones (hence the i > 30 degree cut);
    * the distance ``D``. Angles on the sky become radii through ``R = D *
      angle``, so ``d ln g = -dD / D``.

    The three are added in quadrature. Note that the distance error is really
    shared by every point in a galaxy, so treating it as independent per point
    -- which is what a point-by-point fit does -- is an approximation.
    """
    v = np.asarray(v_kms, dtype=float)
    inc = np.radians(np.asarray(inclination_deg, dtype=float))
    d_inc = np.radians(np.asarray(e_inclination_deg, dtype=float))
    frac = np.sqrt(
        (2.0 * np.asarray(e_v_kms, dtype=float) / v) ** 2
        + (2.0 * d_inc / np.tan(inc)) ** 2
        + (np.asarray(e_distance_mpc, dtype=float) / np.asarray(distance_mpc, dtype=float)) ** 2
    )
    return frac / np.log(10.0)


def log_gbar_error(stellar_fraction, upsilon_scatter_dex=0.11, gas_mass_error=0.10):
    """Uncertainty in log10(g_bar), from the stellar mass-to-light ratio and gas mass.

    ``g_bar`` depends on distance not at all: masses scale as ``D**2`` and radii
    as ``D``, so ``V_bar**2 ~ M/R ~ D`` and ``g_bar = V_bar**2 / R`` cancels.
    What does move it is the assumed mass-to-light ratio -- stellar populations
    vary from galaxy to galaxy by about 0.11 dex -- which matters in proportion
    to how much of ``g_bar`` comes from stars, plus a ~10% calibration
    uncertainty on the 21-cm flux that sets the gas mass.
    """
    f_star = np.clip(np.asarray(stellar_fraction, dtype=float), 0.0, 1.0)
    d_log_gas = np.log10(1.0 + gas_mass_error)
    return np.sqrt((f_star * upsilon_scatter_dex) ** 2 + ((1.0 - f_star) * d_log_gas) ** 2)
