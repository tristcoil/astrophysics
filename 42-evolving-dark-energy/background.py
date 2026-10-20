"""The expansion history of a universe, and the distances BAO and supernovae measure.

A (flat) universe here is described by five numbers:

* ``h``      -- today's expansion rate, H0 = 100 h km/s/Mpc
* ``ombh2``  -- physical density of ordinary matter (baryons), omega_b = Omega_b h^2
* ``omch2``  -- physical density of cold dark matter, omega_c
* ``w0, wa`` -- dark energy's equation of state, w(a) = w0 + wa (1 - a)
  (Chevallier & Polarski 2001; Linder 2003). A cosmological constant is
  w0 = -1, wa = 0.

Everything else is fixed by physics we already know: the CMB temperature
(2.7255 K) sets the photon density; the standard model sets the neutrino
background (N_eff = 3.044, of which one species carries a mass of 0.06 eV); and
flatness gives dark energy whatever is left of the critical density.

The expansion rate at redshift z, in physical-density units, is

    (H(z) / 100 km/s/Mpc)^2 = omega_cb (1+z)^3                 matter (b + c)
                            + omega_r (1+z)^4                  photons + massless nu
                            + omega_nu(z)                      the massive neutrino
                            + omega_de f_de(z)                 dark energy

with f_de(z) = (1+z)^(3(1+w0+wa)) exp(-3 wa z/(1+z)) for the w0-wa law. The
massive neutrino behaves like radiation early and like matter late; its density
is the exact Fermi-Dirac integral. Radiation matters here because the CMB
measures a distance out to z ~ 1090, where it is a third of the total.

Distances (all comoving, in Mpc):

    D_C(z) = int_0^z c dz' / H(z')        line-of-sight comoving distance
    D_M(z) = D_C(z)                        transverse (flat space)
    D_H(z) = c / H(z)                      "Hubble distance" at z
    D_V(z) = [z D_M^2 D_H]^(1/3)           the isotropic BAO average

and the sound horizon, the ruler BAO measure in:

    r_s(z) = int_z^inf c_s(z') dz' / H(z'),   c_s = c / sqrt(3 (1 + R)),
    R = 3 rho_b / (4 rho_gamma)  (baryons slow the sound down)

evaluated at the drag epoch z_d (r_d, for BAO) and at photon decoupling z*
(r_*, for the CMB's acoustic angle theta_* = r_* / D_M(z*)).

Every function here takes *batches* of universes (arrays of parameters) and
returns arrays of shape (n_models, n_redshifts), so an MCMC can evaluate all of
its walkers at once. Accuracy against CAMB is checked in ``test_physics.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from config import C_KM_S, M_NU_EV, N_EFF, NU_MASS_TO_OMEGA, T_CMB

# --- Physical constants (CODATA, identical to CAMB's) -------------------------
_G = 6.67430e-11
_SIGMA_SB = 5.670374419e-8
_C = 2.99792458e8
_MPC_M = 3.085677581e22
_K_B_EV = 8.617333262e-5  # Boltzmann constant in eV/K

#: Photon density omega_gamma = Omega_gamma h^2 for a 2.7255 K blackbody.
_H100_SI = 100e3 / _MPC_M
OMEGA_GAMMA = (4.0 * _SIGMA_SB * T_CMB**4 / _C**3) * 8.0 * np.pi * _G / (3.0 * _H100_SI**2)

#: One neutrino species, relativistic, relative to photons: (7/8) (4/11)^(4/3).
NU_PER_SPECIES = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)

#: CAMB's convention (``share_delta_neff``): the massive eigenstate carries
#: N_eff/3 of a species, the remaining 2 N_eff/3 are massless.
G_MASSIVE = N_EFF / 3.0
N_MASSLESS = N_EFF - G_MASSIVE
OMEGA_R = OMEGA_GAMMA * (1.0 + N_MASSLESS * NU_PER_SPECIES)  # photons + massless nu

#: Today's (non-relativistic) density of the 0.06 eV neutrino, as CAMB 2.0.4
#: computes it. The textbook m / 94.06 eV * (N_eff/3)^(3/4) gives 6.44867e-4,
#: 5e-5 lower (3e-8 in omega: immaterial); we take CAMB's value to compare exactly.
OMEGA_NU = 6.448993754203254e-4 if M_NU_EV == 0.06 else \
    M_NU_EV / 94.06410581217612 * G_MASSIVE**0.75

#: The mapping DESI/Cobaya use between their sampled Omega_m and omega_c:
#: omch2 = Omega_m h^2 - m_nu / 93.14 - ombh2.
OMEGA_NU_MAPPING = M_NU_EV / NU_MASS_TO_OMEGA

#: Speed of light over 100 km/s/Mpc, in Mpc: distances are (C_H100 / E_phys) integrals.
C_H100 = C_KM_S / 100.0

# --- The massive neutrino's density ------------------------------------------

_FD_F0 = 7.0 * np.pi**4 / 120.0          # int x^3 / (e^x + 1) dx
_FD_NR = 1.5 * 1.2020569031595942        # int x^2 / (e^x + 1) dx = (3/2) zeta(3)
#: Gauss-Laguerre nodes for int_0^inf e^-x g(x) dx; the Fermi-Dirac integrand,
#: once e^-x is factored out, is smooth, so 80 nodes give ~1e-12 accuracy.
_GL_X, _GL_W = np.polynomial.laguerre.laggauss(80)


def _fd_energy_ratio(y):
    """rho / rho_relativistic for a thermal fermion with m / kT = y:
    F(y)/F(0), F(y) = int x^2 sqrt(x^2 + y^2) / (e^x + 1) dx."""
    y = np.atleast_1d(np.asarray(y, float))
    x = _GL_X[None, :]
    g = x**2 * np.sqrt(x**2 + y[:, None] ** 2) / (1.0 + np.exp(-x))
    return (g @ _GL_W) / _FD_F0


#: m / kT for the massive neutrino today, chosen so that its non-relativistic
#: density equals OMEGA_NU: omega = G * nu_rel * omega_gamma * y * F_NR / F0.
Y_NU0 = OMEGA_NU / (G_MASSIVE * NU_PER_SPECIES * OMEGA_GAMMA * _FD_NR / _FD_F0)


def omega_nu_of_z(z) -> np.ndarray:
    """Physical density of the massive neutrino at redshift z (omega units, x (1+z)^0).

    Radiation-like (∝ (1+z)^4) while m << kT, matter-like (∝ (1+z)^3) after.
    For 0.06 eV the switch happens around z ~ 100.
    """
    zp1 = 1.0 + np.asarray(z, float)
    rel = G_MASSIVE * NU_PER_SPECIES * OMEGA_GAMMA * zp1**4
    return rel * _fd_energy_ratio(Y_NU0 / zp1).reshape(zp1.shape)


# --- A batch of universes ------------------------------------------------------


@dataclass
class Cosmo:
    """Parameters of n universes (1-D arrays of equal length; scalars broadcast)."""

    h: np.ndarray
    ombh2: np.ndarray
    omch2: np.ndarray
    w0: np.ndarray
    wa: np.ndarray

    def __post_init__(self):
        arrs = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float)) for v in
                                     (self.h, self.ombh2, self.omch2, self.w0, self.wa)))
        self.h, self.ombh2, self.omch2, self.w0, self.wa = (a.copy() for a in arrs)

    @property
    def n(self) -> int:
        return self.h.size

    @property
    def omega_cb(self) -> np.ndarray:
        return self.ombh2 + self.omch2

    @property
    def omega_de(self) -> np.ndarray:
        """Dark energy density today, by flatness (CAMB's closure: the massive
        neutrino counted at its non-relativistic density)."""
        return self.h**2 - self.omega_cb - OMEGA_NU - OMEGA_R

    @property
    def Om(self) -> np.ndarray:
        """Omega_m including the massive neutrino (DESI's and CAMB's definition)."""
        return (self.omega_cb + OMEGA_NU) / self.h**2

    @classmethod
    def from_Om(cls, Om, h, ombh2, w0=-1.0, wa=0.0) -> "Cosmo":
        """Build from (Omega_m, h, omega_b) exactly as DESI's runs map them."""
        Om, h, ombh2 = (np.atleast_1d(np.asarray(v, float)) for v in (Om, h, ombh2))
        return cls(h=h, ombh2=ombh2, omch2=Om * h**2 - OMEGA_NU_MAPPING - ombh2, w0=w0, wa=wa)

    def take(self, idx) -> "Cosmo":
        return Cosmo(self.h[idx], self.ombh2[idx], self.omch2[idx], self.w0[idx], self.wa[idx])


def f_de(z, w0, wa) -> np.ndarray:
    """rho_de(z) / rho_de(0) for w(a) = w0 + wa (1 - a). Shape (n_models, n_z)."""
    zp1 = 1.0 + np.asarray(z, float)[None, :]
    w0 = np.atleast_1d(w0)[:, None]
    wa = np.atleast_1d(wa)[:, None]
    return zp1 ** (3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (zp1 - 1.0) / zp1)


def w_of_z(z, w0, wa) -> np.ndarray:
    """The equation of state itself, w(z) = w0 + wa z/(1+z)."""
    z = np.asarray(z, float)
    return np.asarray(w0)[..., None] + np.asarray(wa)[..., None] * z / (1.0 + z)


def hubble_phys2(z, cosmo: Cosmo, omega_nu_z=None) -> np.ndarray:
    """(H(z) / 100 km/s/Mpc)^2 = sum of physical densities. Shape (n_models, n_z)."""
    z = np.asarray(z, float)
    zp1 = (1.0 + z)[None, :]
    if omega_nu_z is None:
        omega_nu_z = omega_nu_of_z(z)
    return (cosmo.omega_cb[:, None] * zp1**3 + OMEGA_R * zp1**4 + omega_nu_z[None, :]
            + cosmo.omega_de[:, None] * f_de(z, cosmo.w0, cosmo.wa))


def hubble(z, cosmo: Cosmo) -> np.ndarray:
    """H(z) in km/s/Mpc."""
    return 100.0 * np.sqrt(hubble_phys2(z, cosmo))


# --- Distances by integration on a shared grid ---------------------------------


class DistanceGrid:
    """Integration nodes shared by all the redshifts a likelihood needs.

    The comoving distance is integrated in x = ln(1+z), where the integrand
    (1+z)/E is smooth (constant for an empty universe, a gentle exponential for
    matter). The nodes are a uniform grid in x with spacing <= ``step``, plus
    every target redshift exactly, so no interpolation is needed for targets.
    Distances to other redshifts (e.g. z*, which varies from model to model)
    are interpolated between nodes (``comoving_at``).
    """

    def __init__(self, z_targets, z_max: float, step: float = 0.002):
        z_targets = np.asarray(z_targets, float)
        x_top = np.log1p(max(z_max, float(z_targets.max(initial=0.0))))
        fill = np.linspace(0.0, x_top, int(np.ceil(x_top / step)) + 1)
        x_nodes, inv = np.unique(np.concatenate([fill, np.log1p(z_targets)]),
                                 return_inverse=True)
        self.x = x_nodes
        self.z = np.expm1(x_nodes)
        self.index = inv[fill.size:]  # node index of each target
        self.z_targets = z_targets
        self.dx = np.diff(x_nodes)
        self.omega_nu = omega_nu_of_z(self.z)  # fixed physics: computed once

    def integrand(self, cosmo: Cosmo) -> np.ndarray:
        """(1+z) c / H(z) in Mpc on the nodes: dD_C/dx."""
        return C_H100 * (1.0 + self.z)[None, :] / np.sqrt(
            hubble_phys2(self.z, cosmo, self.omega_nu))

    def comoving_nodes(self, cosmo: Cosmo, integrand=None) -> np.ndarray:
        f = self.integrand(cosmo) if integrand is None else integrand
        steps = 0.5 * (f[:, 1:] + f[:, :-1]) * self.dx[None, :]
        return np.concatenate([np.zeros((f.shape[0], 1)), np.cumsum(steps, axis=1)], axis=1)

    def comoving(self, cosmo: Cosmo) -> np.ndarray:
        """D_C(z) at the target redshifts, Mpc; shape (n_models, n_targets)."""
        return self.comoving_nodes(cosmo)[:, self.index]

    def comoving_at(self, z_each, cosmo: Cosmo, dc_nodes=None, integrand=None) -> np.ndarray:
        """D_C to one (model-dependent) redshift per model, by integrating to the
        node below it and adding the remaining sliver with a linear integrand."""
        f = self.integrand(cosmo) if integrand is None else integrand
        dc = self.comoving_nodes(cosmo, f) if dc_nodes is None else dc_nodes
        x = np.log1p(np.asarray(z_each, float))
        k = np.clip(np.searchsorted(self.x, x) - 1, 0, self.x.size - 2)
        rows = np.arange(x.size)
        t = (x - self.x[k]) / self.dx[k]
        f_x = f[rows, k] * (1 - t) + f[rows, k + 1] * t
        return dc[rows, k] + 0.5 * (f[rows, k] + f_x) * (x - self.x[k])


# --- The sound horizon -------------------------------------------------------------


class SoundHorizonGrid:
    """r_s(z) = int_z^inf c_s / H dz for z in the recombination era.

    Integrated in x = ln(1+z) from the top down, on nodes from z = 600 to 1e8;
    above 1e8 the universe is pure radiation with c_s = c/sqrt(3), and the tail
    is added analytically.
    """

    def __init__(self, z_lo: float = 600.0, z_hi: float = 1e8, step: float = 0.001):
        self.x = np.linspace(np.log1p(z_lo), np.log1p(z_hi), int(np.ceil(
            (np.log1p(z_hi) - np.log1p(z_lo)) / step)) + 1)
        self.z = np.expm1(self.x)
        self.dx = np.diff(self.x)
        self.omega_nu = omega_nu_of_z(self.z)
        self.z_hi = z_hi

    def r_s(self, z_each, cosmo: Cosmo) -> np.ndarray:
        zp1 = (1.0 + self.z)[None, :]
        R = 0.75 * (cosmo.ombh2[:, None] / OMEGA_GAMMA) / zp1
        e2 = hubble_phys2(self.z, cosmo, self.omega_nu)
        f = C_H100 * zp1 / (np.sqrt(3.0 * (1.0 + R)) * np.sqrt(e2))
        # tail above z_hi: radiation era, H ∝ (1+z)^2, c_s -> c/sqrt(3)
        omega_rad_hi = e2[:, -1] / (1.0 + self.z_hi) ** 4
        tail = C_H100 / (np.sqrt(3.0 * omega_rad_hi) * (1.0 + self.z_hi))
        steps = 0.5 * (f[:, 1:] + f[:, :-1]) * self.dx[None, :]
        above = np.concatenate([np.cumsum(steps[:, ::-1], axis=1)[:, ::-1],
                                np.zeros((f.shape[0], 1))], axis=1) + tail[:, None]
        x = np.log1p(np.asarray(z_each, float))
        k = np.clip(np.searchsorted(self.x, x) - 1, 0, self.x.size - 2)
        rows = np.arange(x.size)
        t = (x - self.x[k]) / self.dx[k]
        f_x = f[rows, k] * (1 - t) + f[rows, k + 1] * t
        return above[rows, k + 1] + 0.5 * (f_x + f[rows, k + 1]) * (self.x[k + 1] - x)


# --- BAO observables -----------------------------------------------------------------


def bao_predictions(kinds, z, dm, hz, rd) -> np.ndarray:
    """D_M/r_d, D_H/r_d or D_V/r_d for each entry of a BAO data vector.

    ``dm`` and ``hz`` are (n_models, n_points) arrays of D_M (Mpc) and H (km/s/Mpc)
    at the points' redshifts; ``rd`` is (n_models,) in Mpc.
    """
    dh = C_KM_S / hz
    dv = np.cbrt(np.asarray(z)[None, :] * dm**2 * dh)
    out = np.empty_like(dm)
    for j, kind in enumerate(kinds):
        out[:, j] = {"DM_over_rs": dm[:, j], "DH_over_rs": dh[:, j],
                     "DV_over_rs": dv[:, j]}[kind]
    return out / np.asarray(rd)[:, None]


# --- Fitting formulae (for comparison only) -------------------------------------


def rd_desi_formula(ombh2, ombch2, neff=N_EFF):
    """DESI DR2 eq. 2 (after Brieden et al. 2023): r_d in Mpc."""
    return 147.05 * (ombh2 / 0.02236) ** -0.13 * (ombch2 / 0.1432) ** -0.23 * (neff / 3.04) ** -0.1


def zdrag_eisenstein_hu(ombh2, ommh2):
    """Drag-epoch redshift from Eisenstein & Hu (1998, eq. 4)."""
    b1 = 0.313 * ommh2**-0.419 * (1 + 0.607 * ommh2**0.674)
    b2 = 0.238 * ommh2**0.223
    return 1291 * ommh2**0.251 / (1 + 0.659 * ommh2**0.828) * (1 + b1 * ombh2**b2)


def zstar_hu_sugiyama(ombh2, ommh2):
    """Photon-decoupling redshift from Hu & Sugiyama (1996, eq. E1)."""
    g1 = 0.0783 * ombh2**-0.238 / (1 + 39.5 * ombh2**0.763)
    g2 = 0.560 / (1 + 21.1 * ombh2**1.81)
    return 1048 * (1 + 0.00124 * ombh2**-0.738) * (1 + g1 * ommh2**g2)
