"""Transit physics: turning a measured light-curve dip into planetary parameters.

Every constant comes from ``astropy.constants`` (IAU 2015 nominal values) rather
than hand-copied numbers, and every public function takes and returns explicit
``astropy.units`` quantities so a unit slip raises instead of silently producing
a wrong answer.

The core relation is simple. A planet of radius ``Rp`` crossing a star of radius
``Rstar`` blocks a fraction of the stellar disc equal to the ratio of their
areas::

    depth = (Rp / Rstar) ** 2

so a measured depth gives the *ratio* of radii directly, and an absolute planet
radius only once you know the star's radius from somewhere else. That dependence
is why the star's radius uncertainty propagates straight into ours -- see
``planet_radius`` for the propagation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import constants as const
from astropy import units as u

# Nominal solar/planetary radii (IAU 2015 Resolution B3). Using the "nominal"
# constants rather than measured ones is the convention in exoplanet work: they
# are defined exactly, so results are reproducible across papers.
R_SUN = const.R_sun
R_EARTH = const.R_earth
R_JUP = const.R_jup
M_SUN = const.M_sun

#: How many Earth radii fit across the Sun -- handy for sanity checks (~109).
RSUN_PER_REARTH = float((R_SUN / R_EARTH).decompose())
#: How many Jupiter radii fit across the Sun (~9.73).
RSUN_PER_RJUP = float((R_SUN / R_JUP).decompose())


# --------------------------------------------------------------------------
# Depth -> radius
# --------------------------------------------------------------------------

def radius_ratio(depth) -> float:
    """Ratio ``Rp/Rstar`` implied by a fractional transit depth.

    ``depth`` is dimensionless and fractional (1e-4 means a 0.01% dip, i.e.
    100 ppm). This inverts ``depth = (Rp/Rstar)**2``.

    This is the *geometric* inversion, which assumes a uniformly bright stellar
    disc. Real stars are limb darkened -- brighter at the centre than the edge --
    so a planet crossing the middle blocks disproportionately bright light and
    the observed dip is deeper than the area ratio alone. Taking sqrt(depth) of
    the *maximum* depth therefore OVERESTIMATES Rp/Rstar, typically by ~5-10% in
    the Kepler bandpass for a Sun-like star. See ``limb_darkening_note``.
    """
    depth = np.asarray(depth, dtype=float)
    if np.any(depth < 0):
        raise ValueError("transit depth cannot be negative")
    return np.sqrt(depth)


def planet_radius(depth, stellar_radius, *, depth_err=None, stellar_radius_err=None):
    """Planet radius from transit depth and stellar radius.

    Returns ``(radius, radius_err)`` as quantities in solar radii; ``radius_err``
    is ``None`` when no input uncertainties are supplied.

    Error propagation for ``Rp = Rstar * sqrt(depth)``. The partial derivatives
    give a fractional error of::

        (dRp/Rp)**2 = (dRstar/Rstar)**2 + (0.5 * ddepth/depth)**2

    Note the factor of one half on the depth term: because the depth enters under
    a square root, a 10% depth error is only a 5% radius error. In practice the
    stellar radius usually dominates the total uncertainty, which is worth
    stating plainly whenever a planet radius is quoted.
    """
    stellar_radius = u.Quantity(stellar_radius, u.R_sun)
    ratio = radius_ratio(depth)
    radius = ratio * stellar_radius

    if depth_err is None and stellar_radius_err is None:
        return radius, None

    frac_sq = 0.0
    if stellar_radius_err is not None:
        stellar_radius_err = u.Quantity(stellar_radius_err, u.R_sun)
        frac_sq += float((stellar_radius_err / stellar_radius).decompose()) ** 2
    if depth_err is not None:
        frac_sq += (0.5 * float(depth_err) / float(depth)) ** 2

    return radius, radius * np.sqrt(frac_sq)


def to_earth_radii(radius) -> float:
    """Express a radius quantity in Earth radii (dimensionless float)."""
    return float(u.Quantity(radius).to(u.R_earth).value)


def to_jupiter_radii(radius) -> float:
    """Express a radius quantity in Jupiter radii (dimensionless float)."""
    return float(u.Quantity(radius).to(u.R_jup).value)


def limb_darkening_note(depth_max, u1: float = 0.4, u2: float = 0.26) -> dict:
    """Estimate how much limb darkening inflates a naive depth-based radius.

    For a quadratic limb-darkening law the surface brightness at normalised
    radial coordinate ``mu = cos(theta)`` is::

        I(mu)/I(0) = 1 - u1*(1 - mu) - u2*(1 - mu)**2

    The disc-averaged intensity relative to the disc centre is::

        <I>/I(0) = 1 - u1/3 - u2/6

    A small planet transiting the *centre* of the disc blocks light of intensity
    ``I(0)``, but the depth we measure is relative to the total (disc-averaged)
    flux. So the observed central depth exceeds the true area ratio by the factor
    ``I(0)/<I>``, and the true area ratio is::

        (Rp/Rstar)**2 = depth_max * (1 - u1/3 - u2/6)

    The defaults ``u1=0.4, u2=0.26`` are typical for a Sun-like star (Teff~5800 K,
    log g~4.4) in the broad Kepler bandpass; they come from theoretical stellar
    atmosphere grids (e.g. Claret & Bloemen 2011), not from our data.

    This is a *first-order correction for a central transit of a small planet*,
    not a substitute for fitting a proper limb-darkened model (Mandel & Agol
    2002). It is included so the size of the systematic is visible rather than
    hidden. For a grazing transit (large impact parameter) the correction has the
    opposite sign, because the planet then covers only the dimmer stellar limb.

    Returns a dict with the correction factor and both radius ratios.
    """
    disc_average = 1.0 - u1 / 3.0 - u2 / 6.0
    corrected_depth = float(depth_max) * disc_average
    naive = float(np.sqrt(depth_max))
    corrected = float(np.sqrt(corrected_depth))
    return {
        "u1": u1,
        "u2": u2,
        "disc_average_intensity": disc_average,
        "corrected_depth": corrected_depth,
        "ratio_naive": naive,
        "ratio_corrected": corrected,
        "overestimate_percent": 100.0 * (naive / corrected - 1.0),
    }


# --------------------------------------------------------------------------
# Orbit
# --------------------------------------------------------------------------

def semi_major_axis(period, stellar_mass):
    """Orbital separation from Kepler's third law.

    ``a**3 = G * M * P**2 / (4 * pi**2)``

    This uses the star's mass alone, i.e. it neglects the planet's mass. For a
    hot Jupiter around a Sun-like star that approximation is good to ~0.03% in
    ``a``, far below the stellar-mass uncertainty, so it is not the limiting term.
    """
    period = u.Quantity(period, u.day)
    stellar_mass = u.Quantity(stellar_mass, u.M_sun)
    a_cubed = const.G * stellar_mass * period**2 / (4 * np.pi**2)
    return (a_cubed ** (1 / 3)).to(u.au)


def scaled_semi_major_axis(period, stellar_mass, stellar_radius) -> float:
    """``a/Rstar`` -- the orbit size in units of the star's radius.

    This dimensionless number controls transit geometry: it sets how long the
    transit lasts relative to the period, and whether a transit happens at all.
    """
    a = semi_major_axis(period, stellar_mass)
    stellar_radius = u.Quantity(stellar_radius, u.R_sun)
    return float((a / stellar_radius).decompose())


def transit_duration_central(period, stellar_mass, stellar_radius, planet_radius_q=None):
    """Expected duration of a *central* (impact parameter b=0) transit.

    Geometrically the planet crosses a chord of length ``2*(Rstar + Rp)`` while
    travelling around its orbit, giving::

        T14 = (P / pi) * arcsin((Rstar + Rp) / a)

    A real transit with impact parameter ``b`` is shorter by a factor
    ``sqrt(1 - b**2)``, so this returns the MAXIMUM possible duration for the
    system. Comparing it to the observed duration therefore constrains ``b``:
    a measured duration much shorter than this implies a grazing transit.
    """
    period = u.Quantity(period, u.day)
    a = semi_major_axis(period, stellar_mass)
    stellar_radius = u.Quantity(stellar_radius, u.R_sun)
    chord = stellar_radius
    if planet_radius_q is not None:
        chord = chord + u.Quantity(planet_radius_q, u.R_sun)
    ratio = float((chord / a).decompose())
    ratio = min(ratio, 1.0)  # guard against unphysical geometry
    return (period / np.pi * np.arcsin(ratio)).to(u.hour)


def impact_parameter_from_duration(observed_duration, max_duration) -> float | None:
    """Impact parameter ``b`` implied by a transit shorter than the central case.

    From ``T_obs / T_max = sqrt(1 - b**2)``. Returns ``None`` when the observed
    duration exceeds the central-transit maximum, which means the assumed stellar
    mass/radius are inconsistent with the data rather than that ``b`` is imaginary.
    """
    obs = u.Quantity(observed_duration).to(u.hour).value
    mx = u.Quantity(max_duration).to(u.hour).value
    frac = obs / mx
    if frac > 1.0:
        return None
    return float(np.sqrt(max(0.0, 1.0 - frac**2)))


def equilibrium_temperature(stellar_teff, stellar_radius, period, stellar_mass, albedo: float = 0.3):
    """Planet equilibrium temperature.

    ``T_eq = Teff * sqrt(Rstar / (2a)) * (1 - A)**0.25``

    This balances absorbed stellar flux against thermal re-radiation for a body
    that redistributes heat evenly over its whole surface. It is NOT a surface
    temperature: it ignores any greenhouse warming, internal heat, and tidal
    locking (a tidally locked planet has a hot day side well above ``T_eq`` and a
    cold night side below it). The Bond albedo ``A`` is assumed, not measured --
    0.3 is Earth-like, and hot Jupiters are typically darker (~0.1). Quoting
    ``T_eq`` without these caveats is one of the more common ways exoplanet
    numbers get oversold.
    """
    stellar_teff = u.Quantity(stellar_teff, u.K)
    a = semi_major_axis(period, stellar_mass)
    stellar_radius = u.Quantity(stellar_radius, u.R_sun)
    ratio = float((stellar_radius / (2 * a)).decompose())
    return (stellar_teff * np.sqrt(ratio) * (1.0 - albedo) ** 0.25).to(u.K)


# --------------------------------------------------------------------------
# Detection significance
# --------------------------------------------------------------------------

def transit_snr(depth: float, scatter: float, n_in_transit: int) -> float:
    """Signal-to-noise of a folded transit detection.

    ``SNR = (depth / sigma) * sqrt(N_in_transit)``

    where ``sigma`` is the point-to-point scatter of the out-of-transit flux and
    ``N_in_transit`` is the number of individual cadences that fall inside
    transit across ALL the transits being folded together. Folding many transits
    is what makes shallow planets detectable: the signal adds coherently while
    the noise grows only as sqrt(N).

    This assumes white (uncorrelated) noise. Real light curves have correlated
    "red" noise from stellar granulation and instrumental drift, so a real
    detection threshold sits above the white-noise SNR this returns -- treat it
    as an optimistic upper bound.
    """
    if scatter <= 0 or n_in_transit <= 0:
        return float("nan")
    return float(depth / scatter * np.sqrt(n_in_transit))


@dataclass(frozen=True)
class DerivedParameters:
    """Everything we can infer about a planet from its light curve plus the star."""

    radius_ratio: float
    planet_radius_rearth: float
    planet_radius_rjup: float
    planet_radius_err_rearth: float | None
    #: Radius after the first-order limb-darkening correction, which assumes a
    #: central transit. Because that assumption over-corrects for a planet with a
    #: non-zero impact parameter, the true radius lies between this value and the
    #: uncorrected one -- the pair brackets the answer.
    planet_radius_rearth_ld_corrected: float
    planet_radius_rjup_ld_corrected: float
    semi_major_axis_au: float
    scaled_semi_major_axis: float
    equilibrium_temperature_k: float
    max_duration_hours: float
    impact_parameter: float | None
    limb_darkening: dict


def derive_all(
    *,
    depth: float,
    period_days: float,
    stellar_radius_rsun: float,
    stellar_mass_msun: float,
    stellar_teff_k: float,
    observed_duration_hours: float | None = None,
    depth_err: float | None = None,
    stellar_radius_err_rsun: float | None = None,
    albedo: float = 0.3,
) -> DerivedParameters:
    """Run the full chain from a measured depth and period to planet parameters."""
    period = period_days * u.day
    rstar = stellar_radius_rsun * u.R_sun

    rp, rp_err = planet_radius(
        depth, rstar, depth_err=depth_err, stellar_radius_err=stellar_radius_err_rsun
    )
    ld = limb_darkening_note(depth)
    rp_ld = ld["ratio_corrected"] * rstar
    a = semi_major_axis(period, stellar_mass_msun)
    max_dur = transit_duration_central(period, stellar_mass_msun, rstar, rp)

    b = None
    if observed_duration_hours is not None:
        b = impact_parameter_from_duration(observed_duration_hours * u.hour, max_dur)

    return DerivedParameters(
        radius_ratio=float(radius_ratio(depth)),
        planet_radius_rearth=to_earth_radii(rp),
        planet_radius_rjup=to_jupiter_radii(rp),
        planet_radius_err_rearth=(to_earth_radii(rp_err) if rp_err is not None else None),
        planet_radius_rearth_ld_corrected=to_earth_radii(rp_ld),
        planet_radius_rjup_ld_corrected=to_jupiter_radii(rp_ld),
        semi_major_axis_au=float(a.to(u.au).value),
        scaled_semi_major_axis=scaled_semi_major_axis(period, stellar_mass_msun, rstar),
        equilibrium_temperature_k=float(
            equilibrium_temperature(
                stellar_teff_k, rstar, period, stellar_mass_msun, albedo=albedo
            ).value
        ),
        max_duration_hours=float(max_dur.to(u.hour).value),
        impact_parameter=b,
        limb_darkening=limb_darkening_note(depth),
    )
