"""The physics of the ISW effect: expansion, growth, potentials and CAMB spectra.

The key equation
----------------
A CMB photon crossing a gravitational potential Phi gains energy falling in and
loses it climbing out. If Phi does not change while the photon is inside, the
two cancel exactly. If Phi changes, the photon keeps the difference. Summed
along the whole path (Sachs & Wolfe 1967)::

    dT/T (n) = 2 * Integral dPhi/d(eta) d(eta)        (eta = conformal time, c = 1)

In linear theory a density pattern only grows in amplitude, delta(x, a) =
D(a) delta_0(x), and Poisson's equation in comoving coordinates,

    laplacian Phi = (3/2) Omega_m (H0/c)^2 delta / a,

says the potential is Phi(x, a) = g(a) Phi_0(x) with g = D/a. So everything
hinges on whether D grows like a:

* In a flat universe of matter only, D = a exactly. g is constant, potentials
  are frozen, and the ISW effect is **exactly zero**.
* With dark energy, growth falls behind expansion at late times (D grows
  slower than a), g decays, wells get shallower and hills flatter. A photon
  leaving a supercluster's well climbs out of a shallower well than it fell
  into: it arrives hotter. Through a supervoid's hill it arrives colder.

Writing dPhi/d(eta) = Phi_0 dg/d(eta) and dg/d(eta) = (aH/c) g (f - 1), with
f = d ln D / d ln a the growth rate, gives the line-of-sight integral used
throughout this project::

    dT/T (n) = Integral K(chi) Phi_0(chi n) d chi,    K = 2 (a H / c) g (f - 1).

K is negative when f < 1 (decaying potentials) and vanishes when f = 1. This is
the same as equation 3 of Kovacs et al. (2022), dT/T = -2 Int a (1-f) Phi dz.

Units: distances in h^-1 Mpc, H/c in h Mpc^-1, potentials dimensionless
(Phi / c^2), so dT/T is dimensionless.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy import integrate, interpolate

from config import COSMO, COSMO_OMEGA_M

C_KMS = 299792.458
#: Hubble distance c/H0 in h^-1 Mpc (H0 = 100 h km/s/Mpc).
DH = C_KMS / 100.0
#: Omega_nu h^2 per eV of neutrino mass (Planck convention, 93.14 eV).
NU_EV = 93.14


@dataclass(frozen=True)
class Cosmology:
    """A flat (or curved) Lambda + matter background with linear growth.

    Radiation is neglected: at z < 5 it changes distances by < 0.2% and the
    growth rate by less. The default is D3's measured cosmology.
    """

    omega_m: float = COSMO_OMEGA_M
    omega_l: float = 1.0 - COSMO_OMEGA_M
    h: float = COSMO["H0"] / 100.0

    @property
    def omega_k(self) -> float:
        return 1.0 - self.omega_m - self.omega_l

    def E(self, z):
        z = np.asarray(z, float)
        return np.sqrt(self.omega_m * (1 + z) ** 3 + self.omega_k * (1 + z) ** 2 + self.omega_l)

    # --- distances -----------------------------------------------------------

    @cached_property
    def _zgrid(self):
        return np.concatenate([np.linspace(0, 3, 3001), np.linspace(3.01, 1200, 4000)])

    @cached_property
    def _chi_table(self):
        z = self._zgrid
        chi = integrate.cumulative_trapezoid(1.0 / self.E(z), z, initial=0.0) * DH
        return z, chi

    def chi(self, z):
        """Comoving (radial) distance in h^-1 Mpc."""
        z0, c0 = self._chi_table
        return np.interp(z, z0, c0)

    def z_of_chi(self, chi):
        z0, c0 = self._chi_table
        return np.interp(chi, c0, z0)

    # --- growth --------------------------------------------------------------

    def _growth_integral(self, a):
        """I(a) = Int_0^a da' / (a' E(a'))^3, the core of the exact Lambda-CDM growth."""
        a = np.atleast_1d(np.asarray(a, float))
        out = np.empty_like(a)
        for i, ai in enumerate(a):
            out[i] = integrate.quad(lambda x: (x * self.E(1 / x - 1)) ** -3, 0.0, ai,
                                    epsabs=0, epsrel=1e-10, limit=200)[0]
        return out

    @cached_property
    def _growth_table(self):
        """D(a) and f(a) on a grid; D normalised to 1 today.

        For matter plus a cosmological constant (and curvature) the growing mode
        is exactly D(a) proportional to E(a) I(a) (Heath 1977). Its logarithmic
        derivative splits into two terms that are both analytic.
        """
        a = np.concatenate([np.geomspace(1e-3, 0.05, 60), np.linspace(0.0505, 1.0, 400)])
        z = 1 / a - 1
        I = self._growth_integral(a)
        E = self.E(z)
        D = E * I
        D /= D[-1]
        dlnE = (-3 * self.omega_m * a**-3 - 2 * self.omega_k * a**-2) / (2 * E**2)
        f = dlnE + 1.0 / (a**2 * E**3 * I)
        return a, D, f

    def D(self, z):
        a0, D0, _ = self._growth_table
        return np.interp(1 / (1 + np.asarray(z, float)), a0, D0)

    def f(self, z):
        a0, _, f0 = self._growth_table
        return np.interp(1 / (1 + np.asarray(z, float)), a0, f0)

    def g(self, z):
        """Potential growth g = D/a = D (1+z), equal to 1 today."""
        return self.D(z) * (1 + np.asarray(z, float))

    def isw_kernel(self, z):
        """K(z) = 2 (aH/c) g (f - 1) in h Mpc^-1: weight of Phi_0 per unit distance."""
        z = np.asarray(z, float)
        return 2.0 * self.E(z) / DH / (1 + z) * self.g(z) * (self.f(z) - 1.0)

    @property
    def poisson_factor(self) -> float:
        """(3/2) Omega_m (H0/c)^2 in h^2 Mpc^-2: laplacian Phi_0 = this * delta_0."""
        return 1.5 * self.omega_m / DH**2


LCDM = Cosmology()
#: A flat universe of matter only (Einstein-de Sitter): no ISW at all.
EDS = Cosmology(omega_m=1.0, omega_l=0.0)


# --------------------------------------------------------------------------
# Spherical structures
# --------------------------------------------------------------------------


def spherical_potential(r, delta0, cosmo: Cosmology = LCDM):
    """Phi_0(r) of a spherically symmetric linear density profile delta_0(r).

    Gauss's law for the enclosed mass plus the outer shells, with Phi -> 0 far
    away::

        Phi_0(r) = -P [ (1/r) Int_0^r delta r'^2 dr' + Int_r^inf delta r' dr' ],

    P = (3/2) Omega_m (H0/c)^2. A void (delta < 0) gives Phi_0 > 0: a hill.
    ``r`` must be a fine, increasing grid starting near 0 and reaching well past
    the structure.
    """
    r = np.asarray(r, float)
    d = np.asarray(delta0, float)
    inner = integrate.cumulative_trapezoid(d * r**2, r, initial=0.0)
    outer_c = integrate.cumulative_trapezoid(d * r, r, initial=0.0)
    outer = outer_c[-1] - outer_c
    with np.errstate(divide="ignore", invalid="ignore"):
        term = np.where(r > 0, inner / r, 0.0)
    return -cosmo.poisson_factor * (term + outer)


def spherical_isw_profile(theta_deg, z_c, r_grid, delta0_r, cosmo: Cosmology = LCDM,
                          t_cmb_uk: float = COSMO["TCMB"] * 1e6, n_los: int = 2001):
    """Temperature profile dT(theta) in microkelvin of one spherical structure.

    The structure sits at redshift ``z_c``; ``delta0_r`` is its linear density
    profile extrapolated to today (delta(z) = D(z) delta_0). Each line of sight is
    integrated exactly through the sphere -- the time dependence of K across the
    structure is kept, not frozen at its centre.
    """
    theta = np.radians(np.atleast_1d(theta_deg))
    phi_r = spherical_potential(r_grid, delta0_r, cosmo)
    rmax = r_grid[-1]
    chi_c = float(cosmo.chi(z_c))
    out = np.empty(len(theta))
    for i, th in enumerate(theta):
        # distances along the ray where the ray is within rmax of the centre
        b = chi_c * np.sin(th)
        if b >= rmax:
            out[i] = 0.0
            continue
        mid = chi_c * np.cos(th)
        half = np.sqrt(rmax**2 - b**2)
        chi = np.linspace(max(mid - half, 1e-3), mid + half, n_los)
        rr = np.sqrt(chi**2 + chi_c**2 - 2 * chi * chi_c * np.cos(th))
        phi = np.interp(rr, r_grid, phi_r, right=0.0)
        k = cosmo.isw_kernel(cosmo.z_of_chi(chi))
        out[i] = np.trapezoid(k * phi, chi) * t_cmb_uk
    return out


def finelli_profile(r, delta0: float, r0: float):
    """Compensated void profile of Finelli et al. (2016), eq. 6, with delta0 > 0 = depth.

    delta(r) = -delta0 (1 - 2 r^2 / (3 r0^2)) exp(-r^2 / r0^2). The underdense core
    is surrounded by an overdense shell; the total mass deficit integrates to zero.
    """
    x2 = (np.asarray(r, float) / r0) ** 2
    return -delta0 * (1 - 2 * x2 / 3) * np.exp(-x2)


# --------------------------------------------------------------------------
# CAMB
# --------------------------------------------------------------------------


def camb_params(*, lmax: int = 200, cosmo_overrides: dict | None = None):
    """CAMBparams for D3's cosmology (flat, one 0.06 eV neutrino, tau fixed)."""
    import camb

    c = dict(COSMO)
    if cosmo_overrides:
        c.update(cosmo_overrides)
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=c["H0"], ombh2=c["ombh2"], omch2=c["omch2"], mnu=c["mnu"],
                       tau=c["tau"], TCMB=c["TCMB"], omk=c.get("omk", 0.0))
    pars.InitPower.set_params(As=np.exp(c["logA"]) * 1e-10, ns=c["ns"])
    pars.set_for_lmax(lmax, lens_potential_accuracy=0)
    return pars


def linear_power(zs, *, kmax: float = 5.0):
    """Linear matter power P(k, z) in (h^-1 Mpc)^3 at the given redshifts, plus sigma8(z=0)."""
    import camb

    pars = camb_params()
    pars.set_matter_power(redshifts=sorted(set([0.0, *zs]), reverse=True), kmax=kmax)
    pars.NonLinear = camb.model.NonLinear_none
    res = camb.get_results(pars)
    k, z, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=kmax, npoints=600)
    return k, np.asarray(z), pk, float(res.get_sigma8_0())


def source_cls(windows, *, lmax: int = 300, counts_redshift: bool = True,
               counts_lensing: bool = False, cosmo_overrides: dict | None = None):
    """Angular spectra of CMB temperature and galaxy-count windows from CAMB.

    ``windows`` is a list of dicts with keys ``z``, ``dndz`` and ``bias``
    (array over z or a constant). Returns CAMB's raw C_l dictionary ('TxT',
    'TxW1', 'W1xW1', ...), l = 0..lmax, dimensionless (T in units of T_CMB).
    The temperature includes every source (Sachs-Wolfe, Doppler, ISW); the
    cross-correlation with low-redshift galaxies is ISW for l >~ 5.
    """
    import camb
    from camb.sources import SplinedSourceWindow

    pars = camb_params(lmax=lmax, cosmo_overrides=cosmo_overrides)
    pars.Want_CMB = True
    pars.SourceTerms.counts_redshift = counts_redshift
    pars.SourceTerms.counts_lensing = counts_lensing
    pars.SourceTerms.limber_windows = True
    ws = []
    for w in windows:
        z = np.asarray(w["z"], float)
        b = w["bias"]
        if np.ndim(b) == 0:
            ws.append(SplinedSourceWindow(z=z, W=np.asarray(w["dndz"], float),
                                          bias_z=np.full_like(z, float(b))))
        else:
            ws.append(SplinedSourceWindow(z=z, W=np.asarray(w["dndz"], float),
                                          bias_z=np.asarray(b, float)))
    pars.SourceWindows = ws
    res = camb.get_results(pars)
    # get_source_cls_dict leaves out the CMB; the unlensed scalar array holds
    # every auto and cross spectrum including TxW (dimensionless, raw C_l).
    return res.get_cmb_unlensed_scalar_array_dict(lmax=lmax, raw_cl=True, CMB_unit=None)
