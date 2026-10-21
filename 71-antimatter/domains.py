"""Could the universe be half antimatter? Gamma rays from domain boundaries.

Cohen, De Rujula & Glashow (1998, "CDG") asked what a universe made of
separate matter and antimatter regions would look like. After recombination
the regions must touch (the CMB is too smooth for empty gaps between them), so
gas flows into every boundary and annihilates. Each proton-antiproton
annihilation makes about five pions: the neutral ones decay into gamma rays
(a spectrum peaking at ~70 MeV), the charged ones into electrons and
positrons that heat the gas and speed up the flow. Those gamma rays, emitted
between recombination and the first galaxies and redshifted by a factor
20-1100, would fill today's sky at around 1 MeV.

This module repeats their calculation with our own code:

1. ``Universe``: expansion rate, gas and photon densities (their fiducial
   Einstein-de Sitter universe, or Planck 2018 with helium and radiation).
2. ``ElectronHeating``: how far annihilation electrons travel (they lose
   energy by Compton scattering CMB photons) and how much heat they leave in
   the gas along the way (Coulomb collisions).
3. ``solve_boundary``: a 1D finite-volume solution of CDG's fluid equations
   (their eqs. 1-3) next to a boundary, in comoving coordinates: Euler
   equations with an exact vacuum-sink boundary (perfect annihilation,
   sigma_ann -> infinity), Compton drag and Compton cooling on the CMB
   (stiff, integrated exactly over each step), and electron heating that
   feeds back on the annihilation rate. Output: J(y), annihilations per unit
   area and time.
4. ``AnnihilationPhotons``: the gamma-ray spectrum of one annihilation at rest
   (n-pion phase space with measured multiplicities, pi0 -> 2 gamma).
5. ``diffuse_spectrum``: redshifting and absorbing those photons (Klein-
   Nishina scattering and pair production) into today's spectrum, for domains
   of a given size d0 (1/d0 = boundary area per volume).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import integrate, interpolate, optimize

import config as cfg

C = cfg.C_LIGHT
EV = cfg.EV
ME_C2 = cfg.M_E_EV * EV                     # J
MP_C2 = cfg.M_P_EV * EV
LAMBDA_C = cfg.HBAR / (cfg.M_E_KG * C)      # reduced Compton wavelength (m)
A_RAD = 4 * 5.670374419e-8 / C              # J m^-3 K^-4
GAMMA_AD = 5.0 / 3.0


# ---------------------------------------------------------------------------
# 1. The background universe
# ---------------------------------------------------------------------------


@dataclass
class Universe:
    """Homogeneous background; y = 1 + z throughout."""

    name: str
    h: float
    omega_m: float
    omega_r: float
    omega_l: float
    eta: float
    yp: float                  # helium mass fraction (0 = hydrogen only, as CDG)
    t_cmb: float = cfg.T_CMB
    _t_of_lny: object = field(default=None, repr=False)
    _y_of_lnt: object = field(default=None, repr=False)

    @property
    def omega_k(self) -> float:
        return 1.0 - self.omega_m - self.omega_r - self.omega_l

    def hubble(self, y):
        """Expansion rate H(y) in 1/s."""
        y = np.asarray(y, float)
        e2 = (self.omega_r * y**4 + self.omega_m * y**3 + self.omega_k * y**2 + self.omega_l)
        return cfg.h0_si(100 * self.h) * np.sqrt(e2)

    def age(self, y):
        """Cosmic time t(y) in seconds (tabulated once)."""
        if self._t_of_lny is None:
            lny = np.linspace(np.log(1e5), 0.0, 4001)
            yy = np.exp(lny)
            # t(y) = int_y^inf dy'/(y' H); start from the radiation/matter-era
            # closed form at y = 1e5
            ymax = yy[0]
            if self.omega_r > 0:
                t0 = 1.0 / (2 * cfg.h0_si(100 * self.h) * math.sqrt(self.omega_r) * ymax**2)
            else:
                t0 = 2.0 / (3 * cfg.h0_si(100 * self.h) * math.sqrt(self.omega_m) * ymax**1.5)
            integrand = 1.0 / self.hubble(yy)              # dt = -dlny / H
            t = t0 + integrate.cumulative_trapezoid(integrand, -lny, initial=0.0)
            self._t_of_lny = interpolate.CubicSpline(lny[::-1], t[::-1])
        return self._t_of_lny(np.log(np.asarray(y, float)))

    def y_of_t(self, t):
        """Inverse of age(): y at cosmic time t (s)."""
        if self._y_of_lnt is None:
            lny = np.linspace(np.log(1e5), 0.0, 4001)          # increasing t
            lnt = np.log(self.age(np.exp(lny)))
            self._y_of_lnt = interpolate.CubicSpline(lnt, lny)
        return np.exp(self._y_of_lnt(np.log(t)))

    # particle content (fully ionised gas; per m^3)
    @property
    def n_gamma0(self) -> float:
        return 2 * cfg.ZETA3 / math.pi**2 * (cfg.K_B * self.t_cmb / (cfg.HBAR * C)) ** 3

    def n_b(self, y):
        return self.eta * self.n_gamma0 * np.asarray(y, float) ** 3

    @property
    def f_e(self) -> float:            # electrons per baryon
        return 1.0 - self.yp / 2.0

    @property
    def f_tot(self) -> float:          # free particles per baryon (H+, He++, e-)
        return 2.0 - 1.25 * self.yp

    @property
    def m_per_baryon(self) -> float:   # kg, including electrons
        return cfg.M_U_KG * (cfg.M_H1_U * (1 - self.yp) + cfg.M_HE4_U / 4 * self.yp)

    def t_gamma(self, y):
        """CMB temperature in kelvin."""
        return self.t_cmb * np.asarray(y, float)

    def t_gamma_ev(self, y):
        return cfg.K_B * self.t_gamma(y) / EV

    def u_gamma(self, y):
        """CMB energy density in J/m^3."""
        return A_RAD * self.t_gamma(y) ** 4

    def comoving_distance(self, y):
        """Comoving distance to redshift y-1, in m."""
        yy = np.linspace(1.0, float(y), 20001)
        return C * np.trapezoid(1.0 / self.hubble(yy), yy)


def cdg_fiducial() -> Universe:
    """CDG's universe: Omega = 1, Lambda = 0, h = 0.75, eta = 2e-10, hydrogen
    only. Their T(y) curve gives T_gamma = 0.2592 eV at y = 1101.6, i.e.
    T0 = 2.731 K; we use FIRAS 1996's 2.728 K."""
    return Universe("CDG 1998 fiducial", h=0.75, omega_m=1.0, omega_r=0.0, omega_l=0.0,
                    eta=2e-10, yp=0.0, t_cmb=2.728)


def planck2018(eta: float | None = None) -> Universe:
    """Planck 2018 LCDM with radiation (photons + 3.044 neutrino species)
    and helium; eta from Planck's Omega_b h^2 unless given."""
    import baryons

    h = cfg.H0 / 100
    om_r = cfg.omega_r_h2() / h**2
    if eta is None:
        eta = baryons.eta_from_ombh2(cfg.OMBH2_PLANCK)
    return Universe("Planck 2018", h=h, omega_m=cfg.OMEGA_M, omega_r=om_r,
                    omega_l=1 - cfg.OMEGA_M - om_r, eta=eta, yp=cfg.YP)


# ---------------------------------------------------------------------------
# 2. Annihilation electrons: range and heating (CDG Appendix A)
# ---------------------------------------------------------------------------


def _tan2_over_x_table():
    th = np.linspace(0.0, 1.5, 3001)
    f = np.where(th > 0, np.tan(th) ** 2 / np.where(th > 0, th, 1.0), 0.0)
    return th, integrate.cumulative_trapezoid(f, th, initial=0.0)


_TH, _HTAB = _tan2_over_x_table()


@dataclass
class ElectronHeating:
    """Energy loss of a relativistic annihilation electron at epoch y.

    * Compton scattering on the CMB: dp/dt = K_gamma beta gamma^2, with
      K_gamma = (4/3) sigma_T u_gamma (CDG eq. A2). Integrating gives the
      speed along the path, beta(x) = cos(K_gamma x / m_e c^2), for an
      initially ultra-relativistic electron, and its range
      L = (pi/2) m_e c^2 / K_gamma (eq. A4), independent of energy.
    * Coulomb collisions with the plasma's electrons: dp/dt = K_M / beta^2,
      K_M = 2 pi r_e^2 m_e c^2 n_e ln Lambda (eq. A3). This is the part that
      heats the gas; it takes over at beta_eq, where
      beta^3 gamma^2 = K_M / K_gamma.
    * Averaged over emission angles, the energy left per unit distance at a
      distance l from the boundary is K_M int_l^L' dx / (x beta(x)^2)
      (eq. A11), with L' the path length at beta_eq.
    """

    universe: Universe
    y: float

    def __post_init__(self):
        u = self.universe
        self.k_gamma = 4.0 / 3.0 * cfg.SIGMA_T * u.u_gamma(self.y)            # J/m
        self.n_e = u.f_e * u.n_b(self.y)
        coul_log = math.log((1 / LAMBDA_C) ** 3 / (16 * math.pi * self.n_e * cfg.ALPHA))
        self.coulomb_log = coul_log
        self.k_m = 2 * math.pi * cfg.R_E**2 * ME_C2 * self.n_e * coul_log      # J/m
        ratio = self.k_m / self.k_gamma
        self.beta_eq = optimize.brentq(lambda b: b**3 / (1 - b * b) - ratio, 1e-9, 1 - 1e-12)
        self.theta_eq = math.acos(self.beta_eq)
        self.range_m = 0.5 * math.pi * ME_C2 / self.k_gamma                      # L
        self.l_prime = self.theta_eq * ME_C2 / self.k_gamma                     # L'

    def path_factor(self, l):
        """G(l) = int_l^L' dx / (x beta^2): ln(L'/l) plus the slowing-down part."""
        l = np.asarray(l, float)
        th_l = np.clip(self.k_gamma * l / ME_C2, 1e-300, None)
        g = (np.log(self.theta_eq / th_l)
             + np.interp(self.theta_eq, _TH, _HTAB) - np.interp(th_l, _TH, _HTAB))
        return np.where(th_l < self.theta_eq, g, 0.0)

    def dedl(self, l, n_e_local=None):
        """Mean energy deposited in the gas per unit distance per electron (J/m)."""
        k_m = self.k_m if n_e_local is None else self.k_m * n_e_local / self.n_e
        return k_m * self.path_factor(l)


# ---------------------------------------------------------------------------
# 3. The flow into a boundary: CDG's fluid equations, solved numerically
# ---------------------------------------------------------------------------


def _hllc(rl, ul, pl, rr, ur, pr, g=GAMMA_AD):
    """HLLC fluxes for the 1D Euler equations (Toro 2009, ch. 10)."""
    cl, cr = np.sqrt(g * pl / rl), np.sqrt(g * pr / rr)
    el = pl / (g - 1) + 0.5 * rl * ul**2
    er = pr / (g - 1) + 0.5 * rr * ur**2
    sl = np.minimum(ul - cl, ur - cr)
    sr = np.maximum(ul + cl, ur + cr)
    sm = (pr - pl + rl * ul * (sl - ul) - rr * ur * (sr - ur)) / (rl * (sl - ul) - rr * (sr - ur))
    fl = np.array([rl * ul, rl * ul**2 + pl, (el + pl) * ul])
    fr = np.array([rr * ur, rr * ur**2 + pr, (er + pr) * ur])

    def star(r, u, p, e, s):
        f = r * (s - u) / (s - sm)
        return np.array([f, f * sm, f * (e / r + (sm - u) * (sm + p / (r * (s - u))))])

    ul_ = np.array([rl, rl * ul, el])
    ur_ = np.array([rr, rr * ur, er])
    fsl = fl + sl * (star(rl, ul, pl, el, sl) - ul_)
    fsr = fr + sr * (star(rr, ur, pr, er, sr) - ur_)
    return np.where(sl >= 0, fl, np.where(sm >= 0, fsl, np.where(sr > 0, fsr, fr)))


def _sink_flux(r, u, p, g=GAMMA_AD):
    """Flux into a perfectly absorbing boundary on the left (x = 0).

    Matter meeting antimatter vanishes, so each fluid sees vacuum: the exact
    solution is a rarefaction fan. Along it the Riemann invariant
    u - 2c/(g-1) is constant; at x/t = 0 the flow is sonic, u = -c, so
    c* = (2c - (g-1)u)/(g+1). For gas at rest this gives the mass flux
    (3/4)^4 rho c for g = 5/3 (CDG eq. 7). Returns (mass, momentum, energy)
    fluxes in the +x direction (negative = into the boundary).
    """
    c = math.sqrt(g * p / r)
    if u + c <= 0:                                  # supersonic into the boundary
        e = p / (g - 1) + 0.5 * r * u * u
        return r * u, r * u * u + p, (e + p) * u
    cs = (2 * c - (g - 1) * u) / (g + 1)
    if cs <= 0:                                     # gas moving away: vacuum
        return 0.0, 0.0, 0.0
    rs = r * (cs / c) ** (2 / (g - 1))
    ps = p * (rs / r) ** g
    us = -cs
    es = ps / (g - 1) + 0.5 * rs * us * us
    return rs * us, rs * us * us + ps, (es + ps) * us


#: CDG's conservative heating: <dE/dl> = 10 y^3 eV/Mpc against their Coulomb
#: coefficient K_M = 5.5 y^3 eV/Mpc (eqs. A3, A12) -- an effective path
#: factor of 1.8, where the full angle-averaged integral gives ~5-8.
CDG_LOG_FACTOR = 10.0 / 5.5


@dataclass
class BoundarySolution:
    y: np.ndarray                   # epochs (1 + z), decreasing
    J: np.ndarray                   # annihilations per m^2 per s
    T_far_ev: np.ndarray            # gas temperature well inside the depletion zone
    depletion_m: np.ndarray         # proper width where density < half of far value
    profiles: dict                  # snapshots: y -> (l_m, n/n_inf, T_eV, v)
    settings: dict

    @property
    def J_cgs(self):
        return self.J * 1e-4


def solve_boundary(universe: Universe, *, y_start: float = 1100.0, y_end: float = 20.0,
                   gbar_e: float = 3.05, heating: str | bool = "full", heating_scale: float = 1.0,
                   dedl_override_a: float | None = None, n_cells: int = 320,
                   dchi0_kpc: float = 0.004, chi_max_mpc: float = 2.0, cfl: float = 0.4,
                   snapshots=(1000, 500, 200, 100, 50, 20), verbose: bool = False
                   ) -> BoundarySolution:
    """Integrate CDG's fluid equations next to one matter-antimatter boundary.

    Comoving coordinate chi > 0 on the matter side; by symmetry the
    antimatter side is its mirror image, and the boundary at chi = 0 is a
    perfect sink. Operator splitting per time step:

    1. conservative Godunov update (HLLC fluxes, exact sink at chi = 0);
    2. Compton drag on the CMB, d(rho v)/dt = -beta_d rho v, integrated
       exactly *together with* the hydrodynamic force over the step, so the
       stiff diffusive regime (y > 400, where drag wins) relaxes to the right
       terminal velocity instead of being damped to zero;
    3. Compton cooling towards T_gamma plus electron heating, exact
       exponential relaxation;
    4. expansion: rho ~ R^-3, rho v ~ R^-4, energy ~ R^-5 (CDG eqs. 1-3).

    Heating by annihilation electrons (``heating``):

    * ``"full"``: the angle-averaged deposition K_M int_l^L' dx/(x beta^2)
      at each distance l from the boundary (CDG eq. A11);
    * ``"cdg"``: CDG's conservative shortcut, <dE/dl> = a y^3 eV/Mpc with
      a = 10 for their eta = 2e-10 hydrogen gas (eq. A12), i.e. a fixed
      10/5.5 times their Coulomb loss coefficient K_M; for other gas
      densities we keep that ratio, <dE/dl> = (10/5.5) K_M(n_e);
    * ``False``: no heating (the gas stays at the CMB temperature).

    ``dedl_override_a`` instead fixes <dE/dl> = a y^3 eV/Mpc exactly, as in
    CDG's own run (used to reproduce their figures).
    """
    u = universe
    g = GAMMA_AD
    mb = u.m_per_baryon
    # comoving grid, geometric from the boundary
    dchi0 = dchi0_kpc * 1e-3 * cfg.MPC
    chi_max = chi_max_mpc * cfg.MPC
    r = optimize.brentq(lambda q: dchi0 * (q**n_cells - 1) / (q - 1) - chi_max, 1.0001, 1.5)
    edges = np.concatenate([[0.0], np.cumsum(dchi0 * r ** np.arange(n_cells))])
    dchi = np.diff(edges)
    centres = 0.5 * (edges[1:] + edges[:-1])

    y = y_start
    t = float(u.age(y))
    t_end = float(u.age(y_end))
    nb = u.n_b(y)
    rho = np.full(n_cells, nb * mb)
    mom = np.zeros(n_cells)
    T0 = u.t_gamma_ev(y) * EV
    eth = np.full(n_cells, u.f_tot * nb * T0 / (g - 1))

    out_y, out_j, out_t, out_d = [], [], [], []
    profiles = {}
    snap = sorted(snapshots, reverse=True)
    j_now = 0.0
    step = 0
    next_log = y_start
    while t < t_end:
        R = 1.0 / y
        nb_inf = u.n_b(y)
        n_tot = u.f_tot * rho / mb
        p = (g - 1) * eth
        v = mom / rho
        cs = np.sqrt(g * p / rho)
        dt_cfl = cfl * np.min(R * dchi / (np.abs(v) + cs))
        dt = min(dt_cfl, 0.002 / float(u.hubble(y)), t_end - t)

        # 1. fluxes
        fl = _hllc(rho[:-1], v[:-1], p[:-1], rho[1:], v[1:], p[1:], g)
        f0 = np.array(_sink_flux(rho[0], v[0], p[0], g))
        fN = np.array([rho[-1] * v[-1], rho[-1] * v[-1] ** 2 + p[-1],
                       (eth[-1] + 0.5 * rho[-1] * v[-1] ** 2 + p[-1]) * v[-1]])
        flux = np.column_stack([f0, fl, fN])                       # (3, n+1)
        etot = eth + 0.5 * mom**2 / rho
        dU = -dt / (R * dchi) * (flux[:, 1:] - flux[:, :-1])
        rho_n = rho + dU[0]
        mom_h = mom + dU[1]
        etot_n = etot + dU[2]
        if np.any(rho_n <= 0):
            raise FloatingPointError(f"negative density at y = {y:.1f}")
        eth_n = etot_n - 0.5 * mom_h**2 / rho_n
        j_now = -f0[0] / mb                                        # annihilations / m^2 / s

        # 2. drag (exact with the hydrodynamic force held fixed over the step)
        ne_loc = u.f_e * rho_n / mb
        beta_d = 4.0 / 3.0 * cfg.SIGMA_T * u.u_gamma(y) * ne_loc / (rho_n * C)
        x = beta_d * dt
        decay = np.exp(-x)
        gain = np.where(x > 1e-8, -np.expm1(-x) / np.where(x > 0, x, 1.0), 1.0 - 0.5 * x)
        mom_n = mom * decay + (mom_h - mom) * gain

        # 3. Compton cooling + electron heating, per particle
        ntot_n = u.f_tot * rho_n / mb
        T = np.maximum((g - 1) * eth_n / ntot_n, 1e-6 * EV)
        kappa = 8.0 / 3.0 * cfg.SIGMA_T * u.u_gamma(y) / (cfg.M_E_KG * C) * u.f_e / u.f_tot
        if heating and j_now > 0:
            l = R * centres
            eh = ElectronHeating(u, y)
            if dedl_override_a is not None:
                dedl = dedl_override_a * y**3 * EV / cfg.MPC * ne_loc / eh.n_e
                dedl = np.where(l < eh.l_prime, dedl, 0.0)
            elif heating == "cdg":
                dedl = CDG_LOG_FACTOR * eh.k_m * ne_loc / eh.n_e
                dedl = np.where(l < eh.l_prime, dedl, 0.0)
            else:
                dedl = eh.dedl(l, ne_loc)
            q_vol = heating_scale * 0.5 * gbar_e * j_now * dedl          # W/m^3
            q_part = (g - 1) * q_vol / ntot_n                            # dT/dt from heating
        else:
            q_part = 0.0
        t_gam = u.t_gamma_ev(y) * EV
        t_eq = t_gam + q_part / kappa
        T = t_eq + (T - t_eq) * np.exp(-kappa * dt)
        eth_n = ntot_n * T / (g - 1)

        # 4. expansion over dt
        y_new = float(u.y_of_t(t + dt))
        s = y_new / y                                    # R_old / R_new
        rho = rho_n * s**3
        mom = mom_n * s**4
        eth = eth_n * s**5
        t += dt
        y = y_new
        step += 1

        if step % 20 == 0 or t >= t_end:
            prof_n = rho / mb / u.n_b(y)
            l = centres / y
            k = np.argmax(prof_n > 0.5)
            out_y.append(y)
            out_j.append(j_now)
            far = np.searchsorted(prof_n, 0.95)
            Tg = (g - 1) * eth / (u.f_tot * rho / mb) / EV
            out_t.append(float(Tg[min(far, n_cells - 1)]))
            out_d.append(float(l[k]))
            if snap and y <= snap[0]:
                profiles[snap.pop(0)] = (l.copy(), prof_n.copy(), Tg.copy(), (mom / rho).copy())
            if verbose and y < next_log:
                print(f"    y = {y:7.1f}  J = {j_now * 1e-4:9.3e} cm^-2 s^-1  "
                      f"T = {out_t[-1]:.3f} eV  depletion {out_d[-1] / cfg.MPC * 1e3:.3g} kpc")
                next_log = y / 1.25
    settings = {"universe": u.name, "y_start": y_start, "y_end": y_end, "gbar_e": gbar_e,
                "heating": heating, "heating_scale": heating_scale,
                "dedl_override_a": dedl_override_a, "n_cells": n_cells,
                "dchi0_kpc": dchi0_kpc, "chi_max_mpc": chi_max_mpc, "cfl": cfl, "steps": step}
    return BoundarySolution(np.array(out_y), np.array(out_j), np.array(out_t),
                            np.array(out_d), profiles, settings)


def j_hydrodynamic(universe: Universe, y, T_ev, gamma=GAMMA_AD):
    """CDG's estimate J = (3/4)^4 n c_s for gas at temperature T (their eq. 7),
    here with the sound speed of the whole plasma (electrons included)."""
    u = universe
    cs = np.sqrt(gamma * u.f_tot * np.asarray(T_ev) * EV / u.m_per_baryon)
    return (0.75**4) * u.n_b(y) * cs


def j_diffusive(universe: Universe, y, T_ev=None):
    """CDG's diffusive estimate J = n sqrt(5 D / (3 pi t)) (their eq. 6), with
    the drag-limited diffusion constant of the whole plasma,
    D = 3 c T n_tot / (4 sigma_T u_gamma n_e), at T = T_gamma unless given."""
    u = universe
    y = np.asarray(y, float)
    T = u.t_gamma_ev(y) if T_ev is None else np.asarray(T_ev)
    D = 3 * C * T * EV * u.f_tot / (4 * cfg.SIGMA_T * u.u_gamma(y) * u.f_e)
    return u.n_b(y) * np.sqrt(5 * D / (3 * math.pi * u.age(y)))


# ---------------------------------------------------------------------------
# 4. Photons from one annihilation at rest
# ---------------------------------------------------------------------------

M_PI0 = 134.9768      # MeV (PDG)
M_PIC = 139.57039


def _phase_space_energies(masses, w_total, n_events, rng):
    """Energies (MeV) of particles from isotropic n-body decay at rest
    (Raubold-Lynch / GENBOD with weights), w_total = total CM energy (MeV).
    Returns energies (n_events, n) and normalised weights."""
    masses = np.asarray(masses, float)
    n = len(masses)
    tecm = w_total - masses.sum()
    if n == 2:
        p = math.sqrt((w_total**2 - (masses[0] + masses[1]) ** 2)
                      * (w_total**2 - (masses[0] - masses[1]) ** 2)) / (2 * w_total)
        e = np.sqrt(p**2 + masses**2)
        return np.tile(e, (n_events, 1)), np.ones(n_events) / n_events
    rnd = np.sort(rng.random((n_events, n - 2)), axis=1)
    rnd = np.concatenate([np.zeros((n_events, 1)), rnd, np.ones((n_events, 1))], axis=1)
    cum_m = np.cumsum(masses)
    invm = rnd * tecm + cum_m[None, :]                          # invariant masses M_1..M_n
    def pdk(a, b, c):
        x = (a - b - c) * (a + b + c) * (a - b + c) * (a + b - c)
        return np.sqrt(np.clip(x, 0, None)) / (2 * a)
    pd = pdk(invm[:, 1:], invm[:, :-1], masses[None, 1:])        # (ev, n-1)
    weight = np.prod(pd, axis=1)
    # build momenta: start with particle 0 at rest frame of M_1, boost successively
    p4 = np.zeros((n_events, n, 4))
    p4[:, 0, :] = 0.0
    # first two-body: M_2 -> M_1 + m_1
    def rand_dir(k):
        cos = 2 * rng.random(k) - 1
        phi = 2 * math.pi * rng.random(k)
        sin = np.sqrt(1 - cos**2)
        return np.column_stack([sin * np.cos(phi), sin * np.sin(phi), cos])
    d = rand_dir(n_events)
    p = pd[:, 0][:, None] * d
    p4[:, 0, 1:] = p
    p4[:, 0, 0] = np.sqrt(masses[0] ** 2 + pd[:, 0] ** 2)
    p4[:, 1, 1:] = -p
    p4[:, 1, 0] = np.sqrt(masses[1] ** 2 + pd[:, 0] ** 2)
    for i in range(2, n):
        # particles 0..i-1 form invariant mass M_i (invm[:, i-1]), recoiling
        # against particle i in the rest frame of M_{i+1}
        d = rand_dir(n_events)
        q = pd[:, i - 1]
        e_sys = np.sqrt(invm[:, i - 1] ** 2 + q**2)
        bvec = -(q / e_sys)[:, None] * d                        # velocity of the system
        b2 = np.sum(bvec**2, axis=1)
        gam = 1 / np.sqrt(1 - b2)
        for j in range(i):
            bp = np.sum(bvec * p4[:, j, 1:], axis=1)
            coef = np.where(b2 > 0, (gam - 1) * bp / np.where(b2 > 0, b2, 1), 0.0)
            e_new = gam * (p4[:, j, 0] + bp)
            p4[:, j, 1:] = p4[:, j, 1:] + (coef + gam * p4[:, j, 0])[:, None] * bvec
            p4[:, j, 0] = e_new
        p4[:, i, 1:] = q[:, None] * d
        p4[:, i, 0] = np.sqrt(masses[i] ** 2 + q**2)
    return p4[:, :, 0], weight / weight.sum()


@dataclass
class AnnihilationPhotons:
    """Inclusive gamma-ray spectrum of p pbar annihilation at rest.

    The annihilation is modelled as n-pion phase space (the "statistical
    model"), with the channel probabilities given as {(n_charged, n_neutral):
    probability}. Neutral pions decay to two photons, flat in energy between
    E(1-beta)/2 and E(1+beta)/2. Returns dN/dE per annihilation on a log grid.
    """

    channels: dict
    w_total: float = 2 * 938.272
    n_events: int = 60000
    seed: int = 7

    def __post_init__(self):
        rng = np.random.default_rng(self.seed)
        self.e_grid = np.geomspace(1.0, 1000.0, 601)              # MeV
        dnde = np.zeros_like(self.e_grid)
        e_pi0, w_pi0 = [], []
        n_pi0_mean = n_pic_mean = 0.0
        e_charged = []
        prob_sum = sum(self.channels.values())
        for (nc, n0), prob in self.channels.items():
            prob /= prob_sum
            masses = [M_PIC] * nc + [M_PI0] * n0
            en, w = _phase_space_energies(masses, self.w_total, self.n_events, rng)
            n_pi0_mean += prob * n0
            n_pic_mean += prob * nc
            for j in range(nc):
                e_charged.append((en[:, j], w * prob))
            for j in range(nc, nc + n0):
                e_pi0.append(en[:, j])
                w_pi0.append(w * prob)
        e_pi0 = np.concatenate(e_pi0)
        w_pi0 = np.concatenate(w_pi0)
        beta = np.sqrt(1 - (M_PI0 / e_pi0) ** 2)
        lo, hi = e_pi0 * (1 - beta) / 2, e_pi0 * (1 + beta) / 2
        # each pi0 gives 2 photons uniformly in [lo, hi]: accumulate exactly per bin
        edges = np.sqrt(self.e_grid[1:] * self.e_grid[:-1])
        edges = np.concatenate([[self.e_grid[0] ** 2 / edges[0]], edges,
                                [self.e_grid[-1] ** 2 / edges[-1]]])
        width = hi - lo
        counts = np.zeros(len(edges) - 1)
        for k in range(0, len(lo), 4000):            # chunks keep memory at ~20 MB
            sl = slice(k, k + 4000)
            cdf = np.clip((edges[None, :] - lo[sl, None]) / width[sl, None], 0, 1)
            counts += 2 * (w_pi0[sl, None] * np.diff(cdf, axis=1)).sum(axis=0)
        dnde = counts / np.diff(edges)
        self.dnde = dnde                                          # photons / MeV / annihilation
        self.n_gamma = float(np.trapezoid(dnde, self.e_grid))
        self.mean_e_gamma = float(np.trapezoid(dnde * self.e_grid, self.e_grid) / self.n_gamma)
        self.n_pi0 = n_pi0_mean
        self.n_pic = n_pic_mean
        e_c = np.concatenate([e for e, _ in e_charged])
        w_c = np.concatenate([w for _, w in e_charged])
        self.mean_e_charged_pion = float(np.sum(e_c * w_c) / np.sum(w_c))
        self.energy_fraction_gamma = float(np.trapezoid(dnde * self.e_grid, self.e_grid)
                                           / self.w_total)

    def __call__(self, e_mev):
        e = np.asarray(e_mev, float)
        val = np.exp(np.interp(np.log(e), np.log(self.e_grid), np.log(np.maximum(self.dnde, 1e-300))))
        return np.where((e >= self.e_grid[0]) & (e <= self.e_grid[-1]), val, 0.0)


# ---------------------------------------------------------------------------
# 5. Today's diffuse spectrum
# ---------------------------------------------------------------------------


def klein_nishina(e_mev):
    """Total Klein-Nishina cross-section per electron (m^2)."""
    k = np.asarray(e_mev, float) / cfg.M_E_EV * 1e6
    k = np.maximum(k, 1e-6)
    t1 = (1 + k) / k**3 * (2 * k * (1 + k) / (1 + 2 * k) - np.log1p(2 * k))
    t2 = np.log1p(2 * k) / (2 * k) - (1 + 3 * k) / (1 + 2 * k) ** 2
    s = 0.75 * cfg.SIGMA_T * (t1 + t2)
    small = k < 1e-3
    return np.where(small, cfg.SIGMA_T * (1 - 2 * k + 5.2 * k * k), s)


def bethe_heitler(e_mev, z=1.0):
    """Pair production on a bare nucleus of charge z, no screening (ionised
    gas), Bethe-Heitler with the high-energy logarithm and the Maximon (1968)
    near-threshold expansion; m^2 per nucleus. Accurate to ~10%, plenty for
    an absorption correction that is itself small."""
    k = np.asarray(e_mev, float) / (cfg.M_E_EV * 1e-6)                    # in m_e c^2
    pref = cfg.ALPHA * cfg.R_E**2 * z * z
    out = np.zeros_like(k)
    hi = k > 6.0
    out[hi] = pref * (28.0 / 9.0 * np.log(2 * k[hi]) - 218.0 / 27.0)
    lo = (k > 2.0) & ~hi
    eps = (2 * k[lo] - 4) / (2 + k[lo] + 2 * np.sqrt(2 * k[lo]))
    out[lo] = pref * 2 * math.pi / 3 * ((k[lo] - 2) / k[lo]) ** 3 * (
        1 + eps / 2 + 23 / 40 * eps**2 + 11 / 60 * eps**3 + 29 / 960 * eps**4)
    # smooth junction at k = 6 is within ~5%
    return np.maximum(out, 0.0)


def triplet(e_mev):
    """Pair production on a free electron (triplet production), high-energy
    form alpha r_e^2 (28/9 ln 2k - 100/9) above its 4 m_e c^2 threshold."""
    k = np.asarray(e_mev, float) / (cfg.M_E_EV * 1e-6)
    val = cfg.ALPHA * cfg.R_E**2 * (28.0 / 9.0 * np.log(2 * np.maximum(k, 1.0)) - 100.0 / 9.0)
    return np.where(k > 4.0, np.maximum(val, 0.0), 0.0)


def gamma_opacity_per_baryon(universe: Universe, e_mev):
    """sigma_gamma per baryon (m^2) in fully ionised H + He gas: Compton
    scattering and triplet pair production on electrons, Bethe-Heitler pair
    production on protons and helium nuclei."""
    u = universe
    xh, yhe = 1 - u.yp, u.yp / 4
    return (u.f_e * (klein_nishina(e_mev) + triplet(e_mev)) + xh * bethe_heitler(e_mev, 1.0)
            + yhe * bethe_heitler(e_mev, 2.0))


def diffuse_spectrum(universe: Universe, sol_y, sol_J, photons, e_today_mev, *,
                     d0_mpc: float = 1000.0, y_min: float = 20.0, y_max: float = 1100.0,
                     absorb: bool = True, plane_distance_gpc: float | None = None):
    """Today's photon intensity (photons cm^-2 s^-1 sr^-1 MeV^-1), CDG eq. 16:

        F(E) = c / (4 pi d0) int dy Phi(E y) J(y) exp(-tau(E, y)) / (H(y) y^2)

    with d(y) = d0/y the proper surface-to-volume scale, Phi the photons per
    annihilation per MeV, and tau the optical depth from emission to today,
    counting every scattered photon as lost (CDG's conservative choice).

    ``plane_distance_gpc``: instead of many domains, ONE flat boundary at that
    comoving distance from us. Seen from here, a plane's area between
    comoving distances r and r + dr is 2 pi r dr (for r beyond the plane's
    closest point), against 4 pi r^2 dr / d0 for domains; so the sky-averaged
    intensity follows by replacing 1/d0 with 1/(2 r)."""
    u = universe
    E = np.atleast_1d(np.asarray(e_today_mev, float))
    order = np.argsort(sol_y)
    lny_tab = np.log(np.asarray(sol_y)[order])
    lnj_tab = np.log(np.maximum(np.asarray(sol_J)[order], 1e-300))
    yy = np.geomspace(y_min, y_max, 1200)
    jj = np.exp(np.interp(np.log(yy), lny_tab, lnj_tab))
    hh = u.hubble(yy)
    if plane_distance_gpc is None:
        inv_d = np.full_like(yy, 1.0 / (d0_mpc * cfg.MPC))
    else:
        r = C * integrate.cumulative_trapezoid(1.0 / u.hubble(np.geomspace(1, y_max, 6000)),
                                               np.geomspace(1, y_max, 6000), initial=0.0)
        rr = np.interp(yy, np.geomspace(1, y_max, 6000), r)
        inv_d = np.where(rr >= plane_distance_gpc * cfg.GPC, 1.0 / (2 * rr), 0.0)
    jj = jj * inv_d
    d0 = 1.0
    out = np.zeros_like(E)
    for i, e in enumerate(E):
        # Phi is per MeV of *emitted* energy; the y from dE_emit = y dE_today,
        # the y^-3 dilution and 1/d = y/d0 combine into CDG's 1/y^2
        src = photons(e * yy) * jj / (hh * yy**2)
        if absorb:
            yg = np.geomspace(1.0, y_max, 1500)
            dtau = (C * gamma_opacity_per_baryon(u, e * yg) * u.n_b(yg) / (u.hubble(yg) * yg))
            tau_g = integrate.cumulative_trapezoid(dtau, yg, initial=0.0)
            tau = np.interp(yy, yg, tau_g)
            src = src * np.exp(-tau)
        out[i] = C / (4 * math.pi * d0) * np.trapezoid(src, yy)
    # dimension: per m^2 s sr MeV -> per cm^2
    return out * 1e-4
