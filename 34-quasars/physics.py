"""Cosmology, Eddington-limited growth and magnitudes: the small formulas used everywhere.

Two cosmologies appear. ``cosmo()`` is project D3's own fit to the Planck CMB
spectrum (H0 = 66.66, Omega_m = 0.325), used for every age and every new
measurement. ``ref_cosmo()`` is the (H0, Omega_m) = (70, 0.3) model that
Richards et al. 2006, Shen et al. 2011 and Shen & Kelly 2012 used, so that we
can reproduce their numbers like for like.
"""

from __future__ import annotations

from functools import lru_cache

import astropy.constants as const
import astropy.units as u
import numpy as np
from astropy.cosmology import FlatLambdaCDM

from config import D3_H0, D3_OM0, RADIATIVE_EFFICIENCY, REF_H0, REF_OM0

SEC_PER_MYR = 1e6 * 365.25 * 86400.0
DEG2_PER_SR = (180.0 / np.pi) ** 2


@lru_cache(maxsize=None)
def cosmo() -> FlatLambdaCDM:
    """D3's cosmology, with CMB radiation and one 0.06 eV neutrino as in Planck's base model."""
    h = D3_H0 / 100.0
    return FlatLambdaCDM(H0=D3_H0, Om0=D3_OM0, Tcmb0=2.7255, Neff=3.046,
                         m_nu=[0.0, 0.0, 0.06] * u.eV, Ob0=0.02221 / h**2, name="D3 (this atlas)")


@lru_cache(maxsize=None)
def ref_cosmo() -> FlatLambdaCDM:
    """The cosmology of the SDSS quasar papers we compare against (no radiation term)."""
    return FlatLambdaCDM(H0=REF_H0, Om0=REF_OM0, Tcmb0=0.0, name="h=0.7, Om=0.3")


# ---------------------------------------------------------------------------
# Ages and redshifts
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _age_table(which: str):
    c = cosmo() if which == "d3" else ref_cosmo()
    z = np.concatenate([np.linspace(0, 10, 2001), np.geomspace(10.0, 200.0, 800)[1:]])
    t = c.age(z).to_value(u.Myr)
    return z, t


def age_myr(z, which: str = "d3"):
    """Age of the universe at redshift z, in Myr (interpolated from a dense table)."""
    zt, tt = _age_table(which)
    return np.interp(np.asarray(z, float), zt, tt)


def z_at_age(t_myr, which: str = "d3"):
    """Inverse of age_myr."""
    zt, tt = _age_table(which)
    return np.interp(np.asarray(t_myr, float), tt[::-1], zt[::-1])


# ---------------------------------------------------------------------------
# Distances, volumes, magnitudes
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def _dist_table(which: str):
    c = cosmo() if which == "d3" else ref_cosmo()
    z = np.linspace(0.0, 8.0, 8001)
    z[0] = 1e-6
    dl = c.luminosity_distance(z).to_value(u.Mpc)
    dvdz = c.differential_comoving_volume(z).to_value(u.Mpc**3 / u.sr)
    return z, dl, dvdz


def distance_modulus(z, which: str = "d3"):
    zt, dl, _ = _dist_table(which)
    return 5.0 * np.log10(np.interp(np.asarray(z, float), zt, dl) * 1e6 / 10.0)


def dvdz_per_deg2(z, which: str = "d3"):
    """Comoving volume per unit redshift per square degree, Mpc^3."""
    zt, _, dv = _dist_table(which)
    return np.interp(np.asarray(z, float), zt, dv) / DEG2_PER_SR


def absolute_mag(m, z, kcorr, which: str = "d3"):
    """M = m - DM(z) - K(z), Oke & Sandage's sign convention (Richards et al. 2006)."""
    return np.asarray(m, float) - distance_modulus(z, which) - kcorr(z)


def luminosity_2500(mi2):
    """log10 of L_nu at 2500 A (erg/s/Hz) from M_i(z=2): Richards et al. 2006, eq. 4."""
    d10pc = (10.0 * u.pc).to_value(u.cm)
    return -0.4 * (np.asarray(mi2, float) + 48.60 + 2.5 * np.log10(3.0)) + np.log10(4 * np.pi * d10pc**2)


# ---------------------------------------------------------------------------
# Eddington-limited growth
# ---------------------------------------------------------------------------

def eddington_per_msun() -> float:
    """L_Edd per solar mass for ionised hydrogen, erg/s: 4 pi G m_p c / sigma_T x M_sun."""
    val = 4 * np.pi * const.G * const.m_p * const.c / const.sigma_T * const.M_sun
    return float(val.to_value(u.erg / u.s))


def eddington_time_myr() -> float:
    """sigma_T c / (4 pi G m_p): the 'Eddington time', ~450 Myr."""
    val = const.sigma_T * const.c / (4 * np.pi * const.G * const.m_p)
    return float(val.to_value(u.Myr))


def efold_time_myr(eps: float = RADIATIVE_EFFICIENCY, lam: float = 1.0) -> float:
    """Mass e-folding (Salpeter) time for radiative efficiency eps at Eddington ratio lam.

    A black hole radiating L = lam L_Edd = eps Mdot_in c^2 keeps (1 - eps) of the
    inflow, so dM/dt = lam (1 - eps) / eps x M / t_Edd.
    """
    return eps / (1.0 - eps) * eddington_time_myr() / lam


def seed_mass(m_final, dt_myr, eps: float = RADIATIVE_EFFICIENCY, lam=1.0):
    """The seed that grows to m_final in dt_myr of uninterrupted accretion."""
    lam = np.asarray(lam, float)
    return np.asarray(m_final, float) * np.exp(-np.asarray(dt_myr, float) * lam
                                               / efold_time_myr(eps, 1.0))


def time_to_grow_myr(m_seed, m_final, eps: float = RADIATIVE_EFFICIENCY, lam=1.0):
    return np.log(np.asarray(m_final, float) / np.asarray(m_seed, float)) * efold_time_myr(eps, 1.0) / lam
