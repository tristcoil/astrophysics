"""The formulas: everything that turns Gaia's measurements into physics.

Kept apart from the pipeline so that each one can be checked against cases with
known answers (``test_physics.py``). Four groups:

* **Distances and magnitudes** -- parallax to distance, apparent to absolute
  magnitude, proper motion to velocity.
* **Motion in three dimensions** -- how one shared space velocity looks from
  Earth at different places on the sky. Stars that move *together* do not all
  show the same proper motion, because we see the group in perspective. This is
  the same geometry that lets the Hyades' motions point at a "convergent point"
  on the sky, and it matters for any cluster that spans several degrees.
* **Photometric quality** -- the corrected BP/RP flux excess of Riello et al.
  (2021), which flags stars whose colours are contaminated by a neighbour or by
  nebulosity.
* **Parallax zero point** -- Gaia DR3 parallaxes are offset by a few hundredths
  of a milliarcsecond, depending on magnitude, colour and position (Lindegren et
  al. 2021). Negligible for a nearby star; a 4% distance error for M67.
"""

from __future__ import annotations

import warnings

import numpy as np

#: One astronomical unit per Julian year, in km/s (IAU 2012 au). A star with
#: proper motion mu (mas/yr) at parallax p (mas) moves across our line of sight
#: at K * mu / p km/s.
K_KMS = 149_597_870.7 / (365.25 * 86_400)

#: Extinction per unit reddening E(B-V) in Gaia's bands, averaged over typical
#: stars (Casagrande & VandenBerg 2018, MNRAS 479, L102, Table 2). With R_V = 3.1
#: these give A_G = 0.884 A_V and E(BP-RP) = 0.432 A_V.
R_G = 2.740
R_BP = 3.374
R_RP = 2.035
R_V = 3.1


# --------------------------------------------------------------------------
# Distances and magnitudes
# --------------------------------------------------------------------------

def distance_pc(parallax_mas):
    """Distance from parallax. Safe only where the parallax is precise."""
    return 1000.0 / np.asarray(parallax_mas, dtype=float)


def distance_modulus(distance_pc_):
    """m - M for a star at this distance: 0 at 10 pc, 5 at 100 pc."""
    return 5.0 * np.log10(np.asarray(distance_pc_, dtype=float)) - 5.0


def absolute_magnitude(apparent_mag, parallax_mas):
    """M = m + 5 log10(parallax / mas) - 10 -- the B1 formula."""
    return np.asarray(apparent_mag, dtype=float) + 5.0 * np.log10(
        np.asarray(parallax_mas, dtype=float)) - 10.0


def tangential_velocity(pm_masyr, parallax_mas):
    """Velocity across the line of sight, km/s."""
    return K_KMS * np.asarray(pm_masyr, dtype=float) / np.asarray(parallax_mas, dtype=float)


def extinction_from_av(av):
    """A_G and E(BP-RP) for a V-band extinction A_V, with mean coefficients.

    The true coefficients depend on each star's temperature by ~10%, which for
    the small extinctions of the clusters here (A_V ~ 0.1) is ~0.01 mag.
    """
    ebv = av / R_V
    return R_G * ebv, (R_BP - R_RP) * ebv


# --------------------------------------------------------------------------
# Sky geometry
# --------------------------------------------------------------------------

def sky_basis(ra_deg, dec_deg):
    """Unit vectors at a point on the sky, in ICRS Cartesian coordinates.

    ``n`` points from us to the star; ``e_ra`` and ``e_dec`` point east and
    north along the sky. A velocity V splits into proper motion along e_ra and
    e_dec, and radial velocity along n.
    """
    a = np.radians(np.asarray(ra_deg, dtype=float))
    d = np.radians(np.asarray(dec_deg, dtype=float))
    n = np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], axis=-1)
    e_ra = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    e_dec = np.stack([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)], axis=-1)
    return n, e_ra, e_dec


def space_velocity(ra_deg, dec_deg, pmra, pmdec, parallax_mas, rv_kms):
    """3D velocity (km/s, ICRS Cartesian) from the six observables.

    Proper motion becomes transverse velocity through the distance; the
    radial velocity supplies the third dimension. Together they give the
    velocity vector relative to the Sun.
    """
    n, e_ra, e_dec = sky_basis(ra_deg, dec_deg)
    scale, pmra, pmdec, rv = (
        np.asarray(x, dtype=float)[..., None]
        for x in (K_KMS / np.asarray(parallax_mas, dtype=float), pmra, pmdec, rv_kms)
    )
    return scale * (pmra * e_ra + pmdec * e_dec) + rv * n


def observed_motion(velocity_kms, ra_deg, dec_deg, parallax_mas):
    """What a star with this 3D velocity would show: (pmra, pmdec, rv).

    The inverse of :func:`space_velocity`, evaluated anywhere on the sky.
    """
    n, e_ra, e_dec = sky_basis(ra_deg, dec_deg)
    v = np.asarray(velocity_kms, dtype=float)
    to_pm = np.asarray(parallax_mas, dtype=float) / K_KMS
    return (e_ra @ v) * to_pm, (e_dec @ v) * to_pm, n @ v


def perspective_offset(ra_deg, dec_deg, ra0_deg, dec0_deg, velocity_kms, parallax0_mas):
    """How much a co-moving star's proper motion differs from the centre's.

    Every member of a cluster shares (nearly) one space velocity, but we see
    that velocity from a slightly different angle at each position. Two
    effects follow. The local east/north directions rotate as you move across
    the sky, which rotates the proper-motion vector; and part of the radial
    velocity leaks into proper motion away from the centre, making a receding
    cluster appear to contract and an approaching one to expand.

    Subtracting this offset puts every star into the cluster's own frame, where
    a genuinely co-moving member has the centre's proper motion wherever it
    sits. For the Pleiades the offset reaches ~1.5 mas/yr at the edge of the
    field -- comparable to the cluster's internal velocity spread.
    """
    pa, pd, _ = observed_motion(velocity_kms, ra_deg, dec_deg, parallax0_mas)
    pa0, pd0, _ = observed_motion(velocity_kms, ra0_deg, dec0_deg, parallax0_mas)
    return pa - pa0, pd - pd0


def vector_to_radec(v):
    """Direction of a 3-vector as (RA, Dec) in degrees."""
    v = np.asarray(v, dtype=float)
    r = np.linalg.norm(v)
    return float(np.degrees(np.arctan2(v[1], v[0])) % 360.0), float(np.degrees(np.arcsin(v[2] / r)))


def convergent_point(velocity_kms):
    """The point on the sky the cluster is heading towards.

    Parallel motions converge by perspective, like railway tracks, on the sky
    position that the velocity vector points at.
    """
    return vector_to_radec(velocity_kms)


def galactic_axes() -> np.ndarray:
    """The Galactic U, V, W directions as ICRS unit vectors (one per row).

    U points to the Galactic centre, V along the direction of Galactic
    rotation, W to the north Galactic pole. Built from astropy's own frame
    definitions rather than typed-in rotation constants.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    axes = SkyCoord(l=[0, 90, 0] * u.deg, b=[0, 0, 90] * u.deg, frame="galactic").icrs
    return np.asarray(axes.cartesian.xyz.value).T


def galactic_velocity(velocity_icrs_kms):
    """(U, V, W) in km/s, heliocentric -- i.e. relative to the Sun."""
    return galactic_axes() @ np.asarray(velocity_icrs_kms, dtype=float)


def angular_separation(ra_deg, dec_deg, ra0_deg, dec0_deg):
    """Great-circle distance in degrees (haversine: stable at small angles)."""
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    a0, d0 = np.radians(ra0_deg), np.radians(dec0_deg)
    h = np.sin((d - d0) / 2) ** 2 + np.cos(d) * np.cos(d0) * np.sin((a - a0) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1))))


def position_angle(ra0_deg, dec0_deg, ra_deg, dec_deg):
    """Direction from point 0 to a second point, degrees from north through east."""
    a0, d0 = np.radians(ra0_deg), np.radians(dec0_deg)
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    y = np.sin(a - a0) * np.cos(d)
    x = np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)
    return float(np.degrees(np.arctan2(y, x)) % 360.0)


def perspective_terms(dec0_deg, pmra, pmdec, parallax_mas, rv_kms):
    """Sizes of the two perspective effects, in mas/yr per radian across the sky.

    *Shrinking*: a receding cluster's members drift towards its centre at a rate
    set by the radial velocity, v_r / (K d). *Fanning*: parallel motions seen
    from slightly different angles aim at one vanishing point, an effect of
    size |mu| sin(dec) from the rotation of the east/north directions across
    the sky. Which one dominates decides what the pattern looks like.
    """
    shrink = abs(rv_kms) * parallax_mas / K_KMS
    fan = float(np.hypot(pmra, pmdec)) * abs(np.sin(np.radians(dec0_deg)))
    return shrink, fan


def gnomonic(ra_deg, dec_deg, ra0_deg, dec0_deg):
    """Tangent-plane coordinates (xi east, eta north) in degrees.

    A flat map of a few degrees of sky without the east-west stretching that
    plotting raw RA against Dec would introduce away from the equator.
    """
    a, d = np.radians(ra_deg), np.radians(dec_deg)
    a0, d0 = np.radians(ra0_deg), np.radians(dec0_deg)
    cosc = np.sin(d0) * np.sin(d) + np.cos(d0) * np.cos(d) * np.cos(a - a0)
    xi = np.cos(d) * np.sin(a - a0) / cosc
    eta = (np.cos(d0) * np.sin(d) - np.sin(d0) * np.cos(d) * np.cos(a - a0)) / cosc
    return np.degrees(xi), np.degrees(eta)


# --------------------------------------------------------------------------
# Measurement uncertainties
# --------------------------------------------------------------------------

def astrometric_covariance(df, *, with_parallax: bool, inflate_by_ruwe: bool = True):
    """Per-star covariance of (pmra, pmdec[, parallax]) from Gaia's errors.

    Gaia publishes each uncertainty and the correlation between every pair, so
    each star's measurement error is an ellipse (or ellipsoid), not a circle.

    Where RUWE > 1 the single-star model fits worse than the errors allow, so
    the errors are scaled up by RUWE -- the standard "Birge ratio" correction
    for an under-fitting model. This mostly affects unresolved binaries and the
    saturated brightest stars, and it lets them take part in the membership
    model with appropriately reduced weight instead of being thrown away.
    """
    sa = df["pmra_error"].to_numpy(float)
    sd = df["pmdec_error"].to_numpy(float)
    rad = np.nan_to_num(df["pmra_pmdec_corr"].to_numpy(float))
    if with_parallax:
        sp = df["parallax_error"].to_numpy(float)
        rpa = np.nan_to_num(df["parallax_pmra_corr"].to_numpy(float))
        rpd = np.nan_to_num(df["parallax_pmdec_corr"].to_numpy(float))
        cov = np.empty((len(df), 3, 3))
        cov[:, 0, 0] = sa ** 2
        cov[:, 1, 1] = sd ** 2
        cov[:, 2, 2] = sp ** 2
        cov[:, 0, 1] = cov[:, 1, 0] = rad * sa * sd
        cov[:, 0, 2] = cov[:, 2, 0] = rpa * sa * sp
        cov[:, 1, 2] = cov[:, 2, 1] = rpd * sd * sp
    else:
        cov = np.empty((len(df), 2, 2))
        cov[:, 0, 0] = sa ** 2
        cov[:, 1, 1] = sd ** 2
        cov[:, 0, 1] = cov[:, 1, 0] = rad * sa * sd
    if inflate_by_ruwe:
        scale = np.maximum(df["ruwe"].to_numpy(float), 1.0)
        cov *= (scale ** 2)[:, None, None]
    return cov


# --------------------------------------------------------------------------
# Photometric quality (Riello et al. 2021, A&A 649, A3)
# --------------------------------------------------------------------------

def expected_excess(bp_rp):
    """Flux excess C expected for clean photometry at this colour (Table 2).

    BP and RP together collect a little more light than G, by an amount that
    depends on colour. A star exceeding this is contaminated by something the
    wider BP/RP windows caught and the narrower G window did not.
    """
    x = np.asarray(bp_rp, dtype=float)
    blue = 1.154360 + 0.033772 * x + 0.032277 * x ** 2
    mid = 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3
    red = 1.057572 + 0.140537 * x
    return np.where(x < 0.5, blue, np.where(x < 4.0, mid, red))


def corrected_excess(excess_factor, bp_rp):
    """C* = C - f(BP-RP): close to zero for clean photometry (Eq. 6)."""
    return np.asarray(excess_factor, dtype=float) - expected_excess(bp_rp)


def sigma_corrected_excess(g):
    """1-sigma scatter of C* for well-behaved isolated stars (Eq. 18)."""
    return 0.0059898 + 8.817481e-12 * np.asarray(g, dtype=float) ** 7.618399


def clean_photometry(excess_factor, bp_rp, g, *, nsigma: float = 3.0, g_min: float = 4.0):
    """True where BP, RP and G agree to within ``nsigma`` of the expected scatter.

    Not applied brighter than ``g_min``: saturation dominates there, and
    Riello et al. advise against it.
    """
    g = np.asarray(g, dtype=float)
    cstar = corrected_excess(excess_factor, bp_rp)
    ok = np.abs(cstar) < nsigma * sigma_corrected_excess(g)
    return np.where(g <= g_min, np.isfinite(bp_rp), ok & np.isfinite(cstar))


# --------------------------------------------------------------------------
# Radial velocities (Katz et al. 2023, A&A 674, A5)
# --------------------------------------------------------------------------

def rv_magnitude_correction(grvs_mag):
    """Known bias of Gaia DR3 radial velocities for fainter stars, km/s (Eq. 5).

    Charge traps in the RVS detectors delay part of a faint star's signal as
    it crosses the CCD, nudging its spectrum redwards. The published DR3
    velocities are not corrected; Katz et al. recommend subtracting this,
    which is zero brighter than G_RVS = 11 and ~0.4 km/s at G_RVS = 14. Valid
    for stars with rv_template_teff < 8500 K.
    """
    g = np.asarray(grvs_mag, dtype=float)
    return np.where(g >= 11.0, 0.02755 * g ** 2 - 0.55863 * g + 2.81129, 0.0)


# --------------------------------------------------------------------------
# Parallax zero point (Lindegren et al. 2021, A&A 649, A4)
# --------------------------------------------------------------------------

_ZPT_LOADED = False


def parallax_zero_point(df):
    """Gaia DR3 parallax offset for each star, in mas, and where it is extrapolated.

    Corrected parallax = parallax - zero_point. The offset is typically -0.02
    to -0.04 mas, so correcting makes parallaxes slightly *larger* and
    distances slightly *smaller*. Uses the official implementation
    (``gaiadr3-zeropoint``). Outside its calibrated range -- G < 6 or > 21, or
    extreme colours -- it returns the value at the edge of the range, which we
    flag.
    """
    global _ZPT_LOADED
    from zero_point import zpt

    if not _ZPT_LOADED:
        zpt.load_tables()
        _ZPT_LOADED = True

    g = df["phot_g_mean_mag"].to_numpy(float)
    nu = df["nu_eff_used_in_astrometry"].to_numpy(float)
    pc = df["pseudocolour"].to_numpy(float)
    lat = df["ecl_lat"].to_numpy(float)
    solved = df["astrometric_params_solved"].to_numpy(int)

    five = solved == 31
    extrapolated = (g <= 6) | (g >= 21) \
        | (five & ((nu <= 1.1) | (nu >= 1.9))) \
        | (~five & ((pc <= 1.24) | (pc >= 1.72)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zp = np.asarray(zpt.get_zpt(g, nu, pc, lat, solved, _warnings=True), dtype=float)
    return zp, extrapolated
