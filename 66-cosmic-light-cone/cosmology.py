"""The expansion of a homogeneous universe, and the horizons and distances it defines.

Everything here follows from one function, the expansion rate H(a) at scale
factor a (a = 1 today, a = 1/(1+z) when light we now see with redshift z left
its source). The Friedmann equation gives

    H(a)^2 = H0^2 [ Om a^-3 + Or a^-4 + Onu(a) + Ok a^-2 + Ode f_de(a) ]

with matter (Om), radiation (photons and massless neutrinos, Or), the massive
neutrino (Onu: radiation-like early, matter-like late; its density is the exact
Fermi-Dirac integral), curvature (Ok) and dark energy (Ode, f_de(1) = 1).

Two clocks run through every result:

    t   = int da / (a H)       cosmic time, what a clock on a galaxy reads
    eta = int da / (a^2 H)     conformal time, the *comoving* distance a light
                               ray covers (c = 1 Gly/Gyr in our units)

so with eta0 = eta(today) and eta_inf = eta(t -> infinity):

    comoving distance to redshift z     chi(z)  = eta0 - eta(z)
    particle horizon at time t          chi_ph  = eta(t)          (how far light has got since the Big Bang)
    event horizon at time t             chi_eh  = eta_inf - eta(t) (how far light sent at t will ever get)
    Hubble sphere at time t             D_HS    = c / H(t)        (proper; recession speed = c there)
    our past light cone                 chi_lc  = eta0 - eta(t)

Proper distance at time t is a(t) times comoving distance, and every comoving
galaxy recedes at v = H D (Hubble's law, exact at all distances in an FRW
universe). Integrals are done once on a fine grid in ln a with 8-point
Gauss-Legendre per cell (accurate to ~1e-12), then evaluated anywhere through
cubic Hermite splines that use the exact derivatives.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

from config import KMSMPC_TO_PER_GYR, Params

# --- Physical constants (CODATA 2018) -------------------------------------------
_G = 6.67430e-11
_SIGMA_SB = 5.670374419e-8
_C = 2.99792458e8
_MPC_M = 3.0856775814913673e22
_K_B_EV = 8.617333262e-5          # eV/K

#: One neutrino species, relativistic, relative to photons: (7/8)(4/11)^(4/3).
NU_REL = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)

_FD_F0 = 7.0 * np.pi**4 / 120.0                    # int x^3/(e^x+1) dx
_GL_X, _GL_W = np.polynomial.laguerre.laggauss(80)  # nodes for int e^-x g(x) dx
_GL8_X, _GL8_W = np.polynomial.legendre.leggauss(8)


def omega_gamma_h2(T_cmb: float) -> float:
    """Photon density Omega_gamma h^2 for a blackbody at T_cmb (2.4728e-5 at 2.7255 K)."""
    h100 = 100e3 / _MPC_M
    rho_crit_100 = 3.0 * h100**2 / (8.0 * np.pi * _G)
    return 4.0 * _SIGMA_SB * T_cmb**4 / _C**3 / rho_crit_100


def fd_energy_ratio(y) -> np.ndarray:
    """rho / rho_relativistic for a thermal fermion with m / kT = y (1 at y = 0, ~0.317 y at y >> 1).

    F(y)/F(0) with F(y) = int x^2 sqrt(x^2 + y^2) / (e^x + 1) dx, by 80-node
    Gauss-Laguerre after factoring out e^-x (smooth remainder; ~1e-12 accurate).
    """
    y = np.atleast_1d(np.asarray(y, float))
    out = np.empty_like(y)
    x = _GL_X[None, :]
    for s in range(0, y.size, 20000):
        yy = y[s:s + 20000, None]
        g = x**2 * np.sqrt(x**2 + yy**2) / (1.0 + np.exp(-x))
        out[s:s + 20000] = (g @ _GL_W) / _FD_F0
    return out


class Universe:
    """One universe, integrated once; every method then answers in Gyr and Gly.

    Arrays in, arrays out. Redshift z, scale factor a, cosmic time t (Gyr since
    the Big Bang) and conformal time eta (Gly) are interchangeable through the
    ``*_of_*`` methods.
    """

    def __init__(self, p: Params, *, a_min: float = 1e-10, a_max: float = 1e4,
                 n_cells: int = 8000):
        self.p = p
        self.h = p.H0 / 100.0
        self.H0 = p.H0 * KMSMPC_TO_PER_GYR            # 1/Gyr
        self.hubble_distance = 1.0 / self.H0          # c/H0 in Gly
        # --- the energy budget today -------------------------------------------
        if p.T_cmb > 0:
            self.Og = omega_gamma_h2(p.T_cmb) / self.h**2
        else:
            self.Og = 0.0
        g = p.N_eff / 3.0 if p.N_eff > 0 else 0.0
        n_massive = len(p.m_nu)
        n_massless = max(p.N_eff - n_massive * g, 0.0)
        self.Or = self.Og * (1.0 + n_massless * NU_REL)  # photons + massless neutrinos
        t_nu0 = (4.0 / 11.0) ** (1.0 / 3.0) * p.T_cmb
        self._nu_y0 = np.array([m / (_K_B_EV * t_nu0) for m in p.m_nu]) if n_massive else np.zeros(0)
        self._nu_rel = g * NU_REL * self.Og              # each massive species, while relativistic
        self.Onu0 = float(sum(self._nu_rel * fd_energy_ratio(y)[0] for y in self._nu_y0))
        self.Om = p.Om
        if p.Ode is None:
            self.Ok = p.Ok
            self.Ode = 1.0 - p.Om - self.Or - self.Onu0 - p.Ok
        else:
            self.Ode = p.Ode
            self.Ok = 1.0 - p.Om - self.Or - self.Onu0 - p.Ode
        if abs(self.Ok) < 1e-15:
            self.Ok = 0.0
        if self.Om + self.Or + self.Onu0 == 0:
            # No Big Bang (empty or Lambda-only): conformal time diverges as a -> 0, so
            # start late enough for double precision; only differences are used then.
            a_min = max(a_min, 1e-4)
        self._integrate(a_min, a_max, n_cells)

    # ------------------------------------------------------------------ H(a)
    def _f_de(self, a):
        p = self.p
        a = np.asarray(a, float)
        cpl = a ** (-3.0 * (1.0 + p.w0 + p.wa)) * np.exp(-3.0 * p.wa * (1.0 - a))
        if p.w_future == "frozen":
            return np.where(a > 1.0, a ** (-3.0 * (1.0 + p.w0)), cpl)
        return cpl

    def _onu(self, a):
        a = np.asarray(a, float)
        out = np.zeros_like(a)
        for y0 in self._nu_y0:
            out = out + self._nu_rel * fd_energy_ratio(y0 * a.ravel()).reshape(a.shape) / a**4
        return out

    def E2_of_a(self, a) -> np.ndarray:
        """(H/H0)^2 at scale factor a."""
        a = np.asarray(a, float)
        e2 = self.Om / a**3 + self.Or / a**4 + self.Ok / a**2 + self.Ode * self._f_de(a)
        if self._nu_y0.size:
            e2 = e2 + self._onu(a)
        return e2

    def H_of_a(self, a) -> np.ndarray:
        """Expansion rate in 1/Gyr."""
        return self.H0 * np.sqrt(self.E2_of_a(a))

    def q_of_a(self, a, step: float = 1e-5) -> np.ndarray:
        """Deceleration parameter q = -1 - dlnH/dlna (q < 0: accelerating)."""
        a = np.asarray(a, float)
        up, dn = self.E2_of_a(a * np.exp(step)), self.E2_of_a(a * np.exp(-step))
        return -1.0 - 0.25 * (np.log(up) - np.log(dn)) / step

    # ------------------------------------------------------------ integration
    def _integrate(self, a_min, a_max, n_cells):
        # Today (x = ln a = 0) is a grid edge, so t0 and eta0 are exact sums.
        n_past = int(round(n_cells * -np.log(a_min) / (np.log(a_max) - np.log(a_min))))
        edges = np.concatenate([np.linspace(np.log(a_min), 0.0, n_past + 1),
                                np.linspace(0.0, np.log(a_max), n_cells - n_past + 1)[1:]])
        half = 0.5 * np.diff(edges)
        mid = 0.5 * (edges[1:] + edges[:-1])
        x = mid[:, None] + half[:, None] * _GL8_X[None, :]
        a = np.exp(x)
        H = self.H_of_a(a)
        dt = (half[:, None] * _GL8_W[None, :] / H).sum(axis=1)            # dt = dx / H
        deta = (half[:, None] * _GL8_W[None, :] / (a * H)).sum(axis=1)    # deta = dx / (a H)
        # Start values at a_min: radiation + matter era, exact closed forms
        # (both stay cancellation-free as either density goes to zero).
        a0 = a_min
        Or_early = self.Or + sum(self._nu_rel for _ in self._nu_y0)  # massive nu relativistic then
        S = np.sqrt(Or_early + self.Om * a0)
        if Or_early + self.Om > 0:
            t0 = (2.0 * a0**2 / 3.0) * (S + 2.0 * np.sqrt(Or_early)) / (S + np.sqrt(Or_early)) ** 2
            eta0 = 2.0 * a0 / (S + np.sqrt(Or_early))
            self.has_big_bang = True
        else:  # empty or Lambda-only universes: only differences are meaningful
            t0 = eta0 = 0.0
            self.has_big_bang = False
        self._x = edges
        self._t = np.concatenate([[t0 / self.H0], t0 / self.H0 + np.cumsum(dt)])
        self._eta = np.concatenate([[eta0 / self.H0], eta0 / self.H0 + np.cumsum(deta)])
        a_e = np.exp(edges)
        H_e = self.H_of_a(a_e)
        self._t_of_x = CubicHermiteSpline(edges, self._t, 1.0 / H_e)
        self._eta_of_x = CubicHermiteSpline(edges, self._eta, 1.0 / (a_e * H_e))
        self._x_of_t = CubicHermiteSpline(self._t, edges, H_e)
        self._x_of_eta = CubicHermiteSpline(self._eta, edges, a_e * H_e)
        self.a_min, self.a_max = a_min, a_max
        self.age = float(self._t_of_x(0.0))
        self.eta0 = float(self._eta_of_x(0.0))
        self.eta_inf = self._eta_infinity()

    def _eta_infinity(self) -> float:
        """Conformal time at t -> infinity: finite only if the expansion ends up accelerating."""
        a_end = self.a_max
        if self.q_of_a(a_end) > -0.05 or not np.isfinite(self.q_of_a(a_end)):
            return np.inf  # decelerating at late times: light eventually reaches everything
        # eta_inf - eta(a_end) = int_{a_end}^inf da/(a^2 H) = int_0^{1/a_end} ds / H(1/s)
        # (quad never evaluates the s = 0 endpoint itself).
        tail, _ = quad(lambda s: 1.0 / float(self.H_of_a(1.0 / s)), 0.0, 1.0 / a_end,
                       epsabs=0, epsrel=1e-12, limit=200)
        return float(self._eta[-1] + tail)

    # ------------------------------------------------------------ conversions
    @staticmethod
    def a_of_z(z):
        return 1.0 / (1.0 + np.asarray(z, float))

    def t_of_a(self, a):
        return self._t_of_x(np.log(np.asarray(a, float)))

    def eta_of_a(self, a):
        return self._eta_of_x(np.log(np.asarray(a, float)))

    def t_of_z(self, z):
        """Cosmic time when light now seen at redshift z left its source (Gyr)."""
        return self._t_of_x(-np.log1p(np.asarray(z, float)))

    def eta_of_z(self, z):
        return self._eta_of_x(-np.log1p(np.asarray(z, float)))

    def a_of_t(self, t):
        return np.exp(self._x_of_t(np.asarray(t, float)))

    def a_of_eta(self, eta):
        return np.exp(self._x_of_eta(np.asarray(eta, float)))

    def eta_of_t(self, t):
        return self._eta_of_x(self._x_of_t(np.asarray(t, float)))

    def t_of_eta(self, eta):
        return self._t_of_x(self._x_of_eta(np.asarray(eta, float)))

    def H_of_z(self, z):
        return self.H_of_a(self.a_of_z(z))

    def H_of_t(self, t):
        return self.H_of_a(self.a_of_t(t))

    # -------------------------------------------------------------- distances
    def chi_of_z(self, z):
        """Comoving distance = proper distance *today* to the source (Gly)."""
        return self.eta0 - self.eta_of_z(z)

    def z_of_chi(self, chi):
        """Redshift of a source now at comoving distance chi (inverse of chi_of_z)."""
        return 1.0 / self.a_of_eta(self.eta0 - np.asarray(chi, float)) - 1.0

    def _sk(self, chi):
        """Transverse comoving distance for a comoving radial distance chi (curvature)."""
        chi = np.asarray(chi, float)
        if self.Ok == 0:
            return chi
        k = np.sqrt(abs(self.Ok)) * self.H0
        return np.sinh(k * chi) / k if self.Ok > 0 else np.sin(k * chi) / k

    def dm_of_z(self, z):
        """Transverse comoving distance D_M (Gly)."""
        return self._sk(self.chi_of_z(z))

    def da_of_z(self, z):
        """Angular-diameter distance D_A = D_M / (1+z) (Gly). In flat space this is also
        the proper distance to the source *when the light left it*."""
        return self.dm_of_z(z) / (1.0 + np.asarray(z, float))

    def dl_of_z(self, z):
        """Luminosity distance D_L = (1+z) D_M (Gly)."""
        return self.dm_of_z(z) * (1.0 + np.asarray(z, float))

    def emission_distance(self, z):
        """Proper (radial) distance between us and the source when the light left: a * chi."""
        return self.chi_of_z(z) / (1.0 + np.asarray(z, float))

    def light_travel_distance(self, z):
        """c x (time the light was in flight): the 'lookback distance' (Gly)."""
        return self.age - self.t_of_z(z)

    # ------------------------------------------------------- recession speeds
    def v_now(self, z):
        """Today's recession speed of the source, in units of c: H0 * chi."""
        return self.H0 * self.chi_of_z(z)

    def v_then(self, z):
        """Recession speed when the light was emitted, in units of c: H(z) * a * chi."""
        return self.H_of_z(z) * self.emission_distance(z)

    # --------------------------------------------------------------- horizons
    def particle_horizon_comoving(self, t):
        """Comoving distance light has covered since the Big Bang, at time t (Gly)."""
        return self.eta_of_t(t)

    def event_horizon_comoving(self, t):
        """Comoving distance light sent at time t can ever cover (Gly); inf if none."""
        return self.eta_inf - self.eta_of_t(t)

    def hubble_sphere_proper(self, t):
        """Proper radius c/H(t) where recession speed equals c (Gly)."""
        return 1.0 / self.H_of_t(t)

    def past_light_cone_comoving(self, t, t_obs=None):
        """Comoving distance of light now arriving (at t_obs) that was at our cone at time t."""
        t_obs = self.age if t_obs is None else t_obs
        return float(self.eta_of_t(t_obs)) - self.eta_of_t(t)

    # ---------------------------------------------------------------- landmarks
    def z_where(self, func, lo: float = 0.01, hi: float = 50.0) -> float:
        """Root of func(z) = 0 (bracketed)."""
        return float(brentq(func, lo, hi, xtol=1e-12, rtol=1e-12, maxiter=400))

    def z_hubble_sphere_now(self) -> float:
        """Redshift of galaxies receding at exactly c today (chi = c/H0)."""
        return self.z_where(lambda z: self.v_now(z) - 1.0)

    def z_event_horizon_now(self) -> float:
        """Redshift of galaxies now crossing our event horizon (chi = eta_inf - eta0)."""
        if not np.isfinite(self.eta_inf):
            return np.inf
        return self.z_where(lambda z: self.chi_of_z(z) - (self.eta_inf - self.eta0))

    def z_max_angular_diameter(self) -> float:
        """Redshift where D_A peaks (smallest angle for a fixed-size object); inf if
        D_A grows forever (empty and Lambda-only universes)."""
        z = np.geomspace(0.05, 2000.0, 4000)
        da = self.da_of_z(z)
        i = int(np.argmax(da))
        if i in (0, z.size - 1):
            return np.inf
        # dD_A/dz = 0, refined by bisection on a finite-difference slope
        slope = lambda zz: float(self.da_of_z(zz * (1 + 1e-6)) - self.da_of_z(zz * (1 - 1e-6)))
        return self.z_where(slope, z[i - 1], z[i + 1])

    def z_light_cone_meets_hubble_sphere(self) -> float:
        """Redshift where our past light cone crosses the Hubble sphere: v_then = c.
        Light emitted beyond it first moved *away* from us. In flat space this is
        also where D_A peaks (proved in the README)."""
        return self.z_where(lambda z: self.v_then(z) - 1.0, 0.05, 50.0)

    def z_equality(self) -> float:
        """Matter-radiation equality (massive neutrinos counted by their actual density)."""
        def f(lz):
            a = 1.0 / (1.0 + np.exp(lz))
            rad = self.Or / a**4
            matter = self.Om / a**3
            if self._nu_y0.size:  # split the neutrino into its relativistic and matter parts
                nu = self._onu(np.array([a]))[0]
                nu_rel = sum(self._nu_rel for _ in self._nu_y0) / a**4
                rad, matter = rad + min(nu, nu_rel), matter + max(nu - nu_rel, 0.0)
            return rad - matter
        return float(np.exp(brentq(f, np.log(100.0), np.log(1e5))))

    def z_acceleration_start(self) -> float:
        """Redshift where q changes sign (deceleration -> acceleration)."""
        return self.z_where(lambda z: self.q_of_a(1.0 / (1.0 + z)), 0.0, 5.0)

    def t_last_seen(self, chi):
        """Cosmic time of the last event on a comoving galaxy at distance chi that light
        will *ever* bring to us (we see it frozen at that age as t -> infinity)."""
        chi = np.asarray(chi, float)
        if not np.isfinite(self.eta_inf):
            return np.full(chi.shape, np.inf)  # no event horizon: we will see everything
        target = self.eta_inf - chi
        out = np.full(chi.shape, np.nan)        # nan: beyond even the future visibility limit
        ok = target > self._eta[0]
        out[ok] = self.t_of_eta(target[ok])
        return out

    def summary(self) -> dict:
        """The headline numbers for this universe."""
        z_hs = self.z_hubble_sphere_now()
        z_eh = self.z_event_horizon_now()
        eh = self.eta_inf - self.eta0
        out = {
            "name": self.p.name, "label": self.p.label, "H0": self.p.H0,
            "Om_cb": self.Om, "Onu0": self.Onu0, "Om_incl_nu": self.Om + self.Onu0,
            "Or": self.Or, "Ode": self.Ode, "Ok": self.Ok, "w0": self.p.w0, "wa": self.p.wa,
            "w_future": self.p.w_future,
            "age_gyr": self.age,
            "hubble_distance_gly": self.hubble_distance,
            "particle_horizon_gly": self.eta0 if self.has_big_bang else np.inf,
            "event_horizon_gly": eh,
            "eta_inf_gly": self.eta_inf,
            "z_hubble_sphere": z_hs,
            "z_event_horizon": z_eh,
            "z_min_angular_size": self.z_max_angular_diameter(),
            "z_light_cone_meets_hubble_sphere": self.z_light_cone_meets_hubble_sphere(),
            "reachable_fraction": (eh / self.eta0) ** 3 if np.isfinite(eh) else 1.0,
            "inside_hubble_sphere_fraction": (self.hubble_distance / self.eta0) ** 3,
        }
        return {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in out.items()}
