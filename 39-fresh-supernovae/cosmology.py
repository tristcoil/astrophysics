"""Distances in a flat LambdaCDM universe, and the redshift frame change.

Only what this project needs: D1 fitted the cosmology, we evaluate it. The
distance modulus of a supernova at cosmological redshift z_cmb, seen from the
Sun (heliocentric redshift z_hel), is

    mu = 5 log10( d_L / 10 pc ),   d_L = (1 + z_hel) * (c / H0) * integral_0^z_cmb dz' / E(z')

with E(z) = sqrt(Om (1+z)^3 + 1 - Om). Same convention as Pantheon+ and D1:
the expansion between the supernova and us is measured in the CMB frame, while
the one factor of (1 + z) for the photons' energy uses what the Sun sees.
Radiation is neglected (D1 showed it changes distances by < 0.0004 mag here).
"""

from __future__ import annotations

import numpy as np

C_KM_S = 299_792.458

#: The Sun's velocity with respect to the cosmic microwave background
#: (Planck 2018 I, Table 3): 369.82 +- 0.11 km/s towards Galactic
#: (l, b) = (264.021, 48.253) deg.
CMB_DIPOLE_V = 369.82
CMB_DIPOLE_L = 264.021
CMB_DIPOLE_B = 48.253


def comoving_distance(z, Om: float, H0: float, n_grid: int = 8001) -> np.ndarray:
    """Line-of-sight comoving distance [Mpc] in flat LambdaCDM, vectorised.

    Integrates 1/E(z) with the trapezoid rule on a grid that is logarithmic
    in z (so relative accuracy is the same at z = 0.001 and z = 2) and
    interpolates; agrees with astropy to < 1e-6 mag up to z = 2.2 (test_physics.py).
    """
    z = np.asarray(z, dtype=float)
    zmax = max(float(np.nanmax(z)) if z.size else 0.0, 1e-6) * 1.0001
    grid = np.concatenate([[0.0], np.geomspace(1e-7, zmax, n_grid)])
    inv_e = 1.0 / np.sqrt(Om * (1 + grid) ** 3 + (1 - Om))
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (inv_e[1:] + inv_e[:-1]) * np.diff(grid))])
    return (C_KM_S / H0) * np.interp(z, grid, cum)


def distance_modulus(z_cmb, z_hel, Om: float, H0: float) -> np.ndarray:
    """mu = 5 log10(d_L / 10 pc) with d_L = (1 + z_hel) D_C(z_cmb)."""
    d_l = (1 + np.asarray(z_hel, float)) * comoving_distance(z_cmb, Om, H0)
    return 5 * np.log10(d_l) + 25


def z_hel_to_cmb(z_hel, ra_deg, dec_deg) -> np.ndarray:
    """Heliocentric -> CMB-frame redshift: remove the Sun's own motion.

    Redshifts compose multiplicatively: (1 + z_hel) = (1 + z_cmb)(1 + z_sun),
    where z_sun = -v cos(theta) / c is the Doppler shift from the Sun's motion
    (theta = angle between the object and the dipole apex). Objects towards
    the apex are blueshifted by our motion, so their z_cmb is larger than z_hel.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    obj = SkyCoord(ra=np.asarray(ra_deg) * u.deg, dec=np.asarray(dec_deg) * u.deg, frame="icrs")
    apex = SkyCoord(l=CMB_DIPOLE_L * u.deg, b=CMB_DIPOLE_B * u.deg, frame="galactic")
    cos_theta = np.cos(obj.separation(apex).rad)
    z_sun = -CMB_DIPOLE_V * cos_theta / C_KM_S
    return (1 + np.asarray(z_hel, float)) / (1 + z_sun) - 1


def sigma_mu_from_velocity(z, sigma_v: float) -> np.ndarray:
    """Magnitude scatter caused by a random velocity sigma_v [km/s] at redshift z.

    A velocity error dv shifts the inferred distance by dv / (c z) (low-z
    Hubble law), i.e. by 5/ln(10) * dv / (c z) magnitudes. The exact
    expression for distances at finite z replaces 1/z with
    (1 + z) / (z (1 + z/2)) (e.g. Davis et al. 2011); used here.
    """
    z = np.asarray(z, float)
    return (5 / np.log(10)) * (sigma_v / C_KM_S) * (1 + z) / (z * (1 + z / 2))
