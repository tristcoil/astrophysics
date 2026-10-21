"""Background cosmology, linear structure growth, and the halo mass function.

Everything a cluster count needs, written out in numpy so that one likelihood
evaluation takes milliseconds (test_physics.py checks each piece against
colossus and CAMB):

* E(z) = H(z)/H0, distances and the comoving volume of a flat universe with
  matter, radiation and a cosmological constant (or none, for the "no dark
  energy" comparison).
* The linear matter power spectrum: the Eisenstein & Hu (1998) transfer
  function with baryon wiggles, a primordial tilt n_s, normalised to sigma_8,
  the rms density contrast in spheres of 8 Mpc/h today.
* sigma(M, z): the rms fluctuation on the scale that encloses mass M,
  growing in time by the linear growth factor D(z).
* The Tinker et al. (2008) halo mass function, the abundance of collapsed
  halos per unit mass and comoving volume, as a function of sigma. Clusters
  are the rare peaks: halos with sigma(M) well below 1, where the abundance
  falls off as exp(-c/sigma^2). That exponential is why counting the most
  massive clusters measures sigma_8 so precisely.

Masses are M500c (Msun): the mass inside the radius where the mean density
is 500 times the critical density at the cluster's redshift. Tinker's fits
are tabulated against overdensity relative to the *mean matter* density,
Delta_m = 500 / Omega_m(z), and interpolated linearly in log10(Delta_m), as
in the CCL library that Ghirardini et al. (2024) used.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy import integrate

import config as C

C_KMS = 299792.458
RHO_CRIT_H2 = 2.77536627e11  # critical density today / h^2, Msun / Mpc^3
OMEGA_GAMMA_H2 = 2.4728e-5 * (C.T_CMB / 2.7255) ** 4  # photons
OMEGA_R_H2 = OMEGA_GAMMA_H2 * (1 + 0.2271 * C.N_EFF)  # + massless neutrinos

# Tinker et al. (2008) Table 2.
_T08_DELTA = np.log10([200, 300, 400, 600, 800, 1200, 1600, 2400, 3200])
_T08_A = np.array([0.186, 0.200, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260])
_T08_a = np.array([1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66])
_T08_b = np.array([2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41])
_T08_c = np.array([1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44])

_K = np.geomspace(1e-5, 2e2, 1500)  # Mpc^-1, integration grid for sigma(R)
_LNK = np.log(_K)
#: sigma(R) is computed on this fixed grid of radii (Mpc) and splined, so the
#: top-hat window matrix W^2(kR) is built once and each cosmology costs one
#: matrix-vector product.
_RGRID = np.geomspace(0.15, 80.0, 360)
_W2 = None


def _window_matrix():
    global _W2
    if _W2 is None:
        x = np.maximum(np.outer(_RGRID, _K), 1e-6)
        _W2 = (3 * (np.sin(x) - x * np.cos(x)) / x**3) ** 2 * (_K**3 * np.gradient(_LNK))[None, :]
    return _W2


@dataclass(frozen=True)
class Cosmology:
    Om: float = 0.3
    sigma8: float = 0.8
    h: float = C.H0_FID / 100
    Ob: float = C.OMEGA_B_FID
    ns: float = C.NS_FID
    Ode: float | None = None  # None = flat (Omega_Lambda = 1 - Om - Or)

    @property
    def Or(self):
        return OMEGA_R_H2 / self.h**2

    @property
    def OL(self):
        return 1.0 - self.Om - self.Or if self.Ode is None else self.Ode

    @property
    def Ok(self):
        return 1.0 - self.Om - self.Or - self.OL

    def with_(self, **kw):
        return replace(self, **kw)

    # -- background -------------------------------------------------------

    def E(self, z):
        zp = 1.0 + np.asarray(z, float)
        return np.sqrt(self.Om * zp**3 + self.Or * zp**4 + self.Ok * zp**2 + self.OL)

    def Om_z(self, z):
        zp = 1.0 + np.asarray(z, float)
        return self.Om * zp**3 / self.E(z) ** 2

    def comoving_distance(self, z):
        """Line-of-sight comoving distance (Mpc), on a fine grid + interpolation."""
        z = np.asarray(z, float)
        grid = np.linspace(0.0, max(float(np.max(z)), 1e-3) * 1.001, 2001)
        dc = integrate.cumulative_trapezoid(1.0 / self.E(grid), grid, initial=0.0)
        return C_KMS / (100 * self.h) * np.interp(z, grid, dc)

    def transverse_distance(self, z):
        dc = self.comoving_distance(z)
        if abs(self.Ok) < 1e-8:
            return dc
        dh = C_KMS / (100 * self.h)
        sk = np.sqrt(abs(self.Ok))
        return dh / sk * (np.sinh(sk * dc / dh) if self.Ok > 0 else np.sin(sk * dc / dh))

    def luminosity_distance(self, z):
        return (1 + np.asarray(z, float)) * self.transverse_distance(z)

    def dV_dz_dOmega(self, z):
        """Comoving volume per unit redshift per steradian (Mpc^3 / sr)."""
        dm = self.transverse_distance(z)
        return C_KMS / (100 * self.h) * dm**2 / self.E(z)

    @property
    def rho_m0(self):
        """Mean matter density today (comoving), Msun / Mpc^3."""
        return self.Om * RHO_CRIT_H2 * self.h**2

    # -- growth -------------------------------------------------------------

    def growth(self, z):
        """Linear growth factor D(z)/D(0) for z of a few or less.

        Heath (1977): D(a) proportional to E(a) * integral_0^a da' / (a' E(a'))^3,
        exact for matter + curvature + Lambda. Radiation is left out here; it
        changes D(z)/D(0) at z < 1 by less than 0.1%. (growth_ode keeps it.)
        """
        z = np.atleast_1d(np.asarray(z, float))
        a = np.geomspace(1e-5, 1.0, 4000)
        Ea = np.sqrt(self.Om * a**-3 + self.Ok * a**-2 + self.OL)
        I = integrate.cumulative_trapezoid(1.0 / (a * Ea) ** 3, a, initial=0.0)
        I += 0.4 * a[0] ** 2.5 / self.Om**1.5  # the matter-dominated part below a[0]
        D = Ea * I
        return np.interp(np.log(1 / (1 + z)), np.log(a), D / D[-1])

    def growth_ode(self, z):
        """D(z)/D(0) from the growth equation, radiation included (for z ~ 1000).

        Solved in ln a from a = 1e-6, deep in radiation domination, where the
        growing mode of the matter perturbation is ~ 1 + 3a/(2 a_eq)
        (Meszaros); the overall normalisation cancels in D(z)/D(0).
        """
        z = np.atleast_1d(np.asarray(z, float))
        a = np.geomspace(1e-6, 1.0, 4000)
        lna = np.log(a)
        E2 = self.Om * a**-3 + self.Or * a**-4 + self.Ok * a**-2 + self.OL
        dlnE2 = (-3 * self.Om * a**-3 - 4 * self.Or * a**-4 - 2 * self.Ok * a**-2) / E2
        om = self.Om * a**-3 / E2

        def rhs(x, y):
            i = np.interp(x, lna, np.arange(lna.size))
            j = int(min(i, lna.size - 2))
            w = i - j
            d2 = (1 - w) * dlnE2[j] + w * dlnE2[j + 1]
            o = (1 - w) * om[j] + w * om[j + 1]
            D, dD = y
            return [dD, -(2 + 0.5 * d2) * dD + 1.5 * o * D]

        aeq = self.Or / self.Om
        y0 = [1 + 1.5 * a[0] / aeq, 1.5 * a[0] / aeq]
        sol = integrate.solve_ivp(rhs, (lna[0], 0.0), y0, t_eval=lna, rtol=1e-9, atol=1e-12, method="DOP853")
        D = sol.y[0] / sol.y[0][-1]
        return np.interp(np.log(1 / (1 + z)), lna, D)

    # -- linear power spectrum ----------------------------------------------

    def transfer_eh98(self, k):
        """Eisenstein & Hu (1998) transfer function with baryon features; k in 1/Mpc."""
        h, Om, Ob = self.h, self.Om, self.Ob
        om, ob = Om * h**2, Ob * h**2
        fb, fc = Ob / Om, (Om - Ob) / Om
        th = C.T_CMB / 2.7
        zeq = 2.50e4 * om * th**-4
        keq = 7.46e-2 * om * th**-2
        b1 = 0.313 * om**-0.419 * (1 + 0.607 * om**0.674)
        b2 = 0.238 * om**0.223
        zd = 1291 * om**0.251 / (1 + 0.659 * om**0.828) * (1 + b1 * ob**b2)
        Rd = 31.5 * ob * th**-4 * (1000 / zd)
        Req = 31.5 * ob * th**-4 * (1000 / zeq)
        s = 2 / (3 * keq) * np.sqrt(6 / Req) * np.log((np.sqrt(1 + Rd) + np.sqrt(Rd + Req)) / (1 + np.sqrt(Req)))
        ksilk = 1.6 * ob**0.52 * om**0.73 * (1 + (10.4 * om) ** -0.95)
        a1 = (46.9 * om) ** 0.670 * (1 + (32.1 * om) ** -0.532)
        a2 = (12.0 * om) ** 0.424 * (1 + (45.0 * om) ** -0.582)
        alpha_c = a1**-fb * a2 ** (-fb**3)
        bb1 = 0.944 / (1 + (458 * om) ** -0.708)
        bb2 = (0.395 * om) ** -0.0266
        beta_c = 1 / (1 + bb1 * (fc**bb2 - 1))
        q = k / (13.41 * keq)

        def T0(alpha, beta):
            L = np.log(np.e + 1.8 * beta * q)
            Cq = 14.2 / alpha + 386 / (1 + 69.9 * q**1.08)
            return L / (L + Cq * q**2)

        f = 1 / (1 + (k * s / 5.4) ** 4)
        Tc = f * T0(1, beta_c) + (1 - f) * T0(alpha_c, beta_c)
        y = (1 + zeq) / (1 + zd)
        G = y * (-6 * np.sqrt(1 + y) + (2 + 3 * y) * np.log((np.sqrt(1 + y) + 1) / (np.sqrt(1 + y) - 1)))
        alpha_b = 2.07 * keq * s * (1 + Rd) ** -0.75 * G
        beta_node = 8.41 * om**0.435
        s_tilde = s / (1 + (beta_node / (k * s)) ** 3) ** (1 / 3)
        beta_b = 0.5 + fb + (3 - 2 * fb) * np.sqrt((17.2 * om) ** 2 + 1)
        Tb = (T0(1, 1) / (1 + (k * s / 5.2) ** 2)
              + alpha_b / (1 + (beta_b / (k * s)) ** 3) * np.exp(-((k / ksilk) ** 1.4))) * np.sinc(k * s_tilde / np.pi)
        return fb * Tb + fc * Tc

    def _shape(self, k):
        """Unnormalised P(k) shape, k^n_s T^2."""
        return k**self.ns * self.transfer_eh98(k) ** 2

    def _sigma2_spline(self):
        """Spline of ln sigma^2 (unnormalised, z = 0) in ln R on the fixed grid."""
        from scipy.interpolate import CubicSpline

        s2 = _window_matrix() @ self._shape(_K) / (2 * np.pi**2)
        return CubicSpline(np.log(_RGRID), np.log(s2))

    def _sigma2_unnorm(self, R):
        """Unnormalised sigma^2 at comoving radii R (Mpc), z = 0."""
        return np.exp(self._sigma2_spline()(np.log(np.atleast_1d(R))))

    def pk_norm(self):
        """Amplitude A such that P(k) = A k^ns T^2 has the requested sigma_8."""
        return self.sigma8**2 / self._sigma2_unnorm(np.array([8.0 / self.h]))[0]

    def sigma_M(self, M, z=0.0):
        """rms linear fluctuation for masses M (Msun) at redshift(s) z.

        Returns shape (len(z), len(M)); also d ln sigma / d ln M.
        """
        M = np.atleast_1d(np.asarray(M, float))
        R = (3 * M / (4 * np.pi * self.rho_m0)) ** (1 / 3)
        spl = self._sigma2_spline()
        norm = self.sigma8**2 / np.exp(spl(np.log(8.0 / self.h)))
        s0 = np.sqrt(norm * np.exp(spl(np.log(R))))
        dlns = 0.5 * spl(np.log(R), 1) / 3.0  # d ln sigma / d ln M = (1/6) d ln sigma^2 / d ln R
        D = self.growth(z)
        return D[:, None] * s0[None, :], dlns

    # -- halo mass function -------------------------------------------------

    @staticmethod
    def tinker08_params(delta_m, z):
        ld = np.log10(delta_m)
        A0, a0, b0, c0 = (np.interp(ld, _T08_DELTA, t) for t in (_T08_A, _T08_a, _T08_b, _T08_c))
        alpha = 10 ** (-((0.75 / np.log10(delta_m / 75.0)) ** 1.2))
        zp = 1.0 + np.asarray(z, float)
        return A0 * zp**-0.14, a0 * zp**-0.06, b0 * zp**-alpha, c0

    def dndlnM(self, M, z, delta_c=500.0):
        """Halo abundance dn/dlnM (comoving Mpc^-3) for M_{delta_c, crit}.

        Shape (len(z), len(M)).
        """
        z = np.atleast_1d(np.asarray(z, float))
        sig, dlns = self.sigma_M(M, z)
        delta_m = delta_c / self.Om_z(z)
        A, a, b, c = self.tinker08_params(delta_m, z)
        f = A[:, None] * ((sig / b[:, None]) ** (-a[:, None]) + 1) * np.exp(-np.atleast_1d(c)[:, None] / sig**2)
        return f * self.rho_m0 / np.atleast_1d(M)[None, :] * np.abs(dlns)[None, :]

    def dndlnM_mean(self, M, z, delta_m=200.0):
        """Same, for masses defined at fixed overdensity relative to the mean."""
        z = np.atleast_1d(np.asarray(z, float))
        sig, dlns = self.sigma_M(M, z)
        A, a, b, c = self.tinker08_params(np.full(z.size, delta_m), z)
        f = A[:, None] * ((sig / b[:, None]) ** (-a[:, None]) + 1) * np.exp(-np.atleast_1d(c)[:, None] / sig**2)
        return f * self.rho_m0 / np.atleast_1d(M)[None, :] * np.abs(dlns)[None, :]

    def S8(self):
        return self.sigma8 * np.sqrt(self.Om / 0.3)


def fiducial(**kw) -> Cosmology:
    """Ghirardini et al.'s eRASS1 best fit (Omega_m 0.29, sigma_8 0.88)."""
    return Cosmology(Om=0.29, sigma8=0.88).with_(**kw)
