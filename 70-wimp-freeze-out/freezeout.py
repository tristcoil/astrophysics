"""Freeze-out: the Boltzmann equation for a relic's abundance, and its analytic approximations.

A particle chi of mass m and g internal states annihilates in pairs into lighter Standard Model
particles. Its number density n obeys (Lee & Weinberg 1977; Kolb & Turner 1990, ch. 5)

    dn/dt + 3 H n = -<sigma v> (n^2 - n_eq^2).

Per unit entropy, Y = n/s, and in terms of x = m/T this becomes (Gondolo & Gelmini 1991)

    dY/dx = -lambda(x) (Y^2 - Y_eq^2),   lambda(x) = sqrt(pi/45) M_Pl g*^(1/2)(T) m <sigma v>(x) / x^2,

with the Maxwell-Boltzmann equilibrium yield Y_eq = 45 g x^2 K_2(x) / (4 pi^4 g*s(T)). Early on
lambda Y_eq is enormous, so Y is pinned to Y_eq (the equation is "stiff"); once the annihilation
rate n <sigma v> drops below the expansion rate H, Y stops following Y_eq, which falls as e^-x,
and freezes. The relic density today is Omega h^2 = m s_0 Y_0 / (rho_crit / h^2).

We integrate in u = ln x for y = ln Y with SciPy's implicit Radau method (5th order, error
controlled), from x = 3 to T = 1 eV, as Saikawa & Shirai (2020) do.

<sigma v> = a + b <v_rel^2> with <v_rel^2> = 6/x: ``a`` is s-wave, ``b`` p-wave annihilation.
Conventions: n counts all dark-matter particles. For a self-conjugate particle (Majorana fermion,
real scalar) <sigma v> is the cross section for chi chi -> SM. For a Dirac fermion with n = n_chi
+ n_chibar, the equation above holds with <sigma v> = <sigma v>_(chi chibar) / 2, so the physical
cross section is twice the value solved for (Saikawa & Shirai's S_4 table uses the first form).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy import integrate, optimize, special

import config as C
import thermo

SQRT_PI_45 = math.sqrt(math.pi / 45.0)
#: Entropy-to-number-density bookkeeping for today's universe.
_CM3_PER_GEV3 = 1.0 / C.HBARC_GEV_CM ** 3   # n [cm^-3] = n [GeV^3] * this


def s0_cm3(g_s0: float, T0_K: float = C.T0_K) -> float:
    """Entropy density today in cm^-3 (2907 for g*s0 = 3.931 and 2.7255 K)."""
    T0 = T0_K * C.K_B_GEV_PER_K
    return 2 * math.pi ** 2 / 45 * g_s0 * T0 ** 3 * _CM3_PER_GEV3


def omega_h2_from_Y(m: float, Y0: float, g_s0: float, T0_K: float = C.T0_K) -> float:
    return m * Y0 * s0_cm3(g_s0, T0_K) / C.RHO_CRIT_H2_GEV_CM3


def ln_kve2(x):
    """ln(K_2(x) e^x). SciPy 1.18's kve(2, x) returns NaN above x ~ 1.1e9, so large arguments use
    the asymptotic series sqrt(pi/2x) (1 + 15/(8x) + 105/(128 x^2) - 315/(1024 x^3))."""
    x = np.asarray(x, dtype=float)
    big = x > 1.0e3
    xs = np.where(big, 1.0, x)
    xb = np.where(big, x, 1.0e3)
    asym = 0.5 * np.log(math.pi / (2 * xb)) + np.log1p(15 / (8 * xb) + 105 / (128 * xb ** 2) - 315 / (1024 * xb ** 3))
    return np.where(big, asym, np.log(special.kve(2, xs)))


def ln_yeq(x, g, g_s, nonrel: bool = False):
    """ln Y_eq. Exact Maxwell-Boltzmann (Bessel K_2) or the non-relativistic x^(3/2) e^-x form."""
    x = np.asarray(x, dtype=float)
    if nonrel:
        return (math.log(45.0 / (2 * math.pi ** 4) * math.sqrt(math.pi / 8)) + np.log(g / g_s)
                + 1.5 * np.log(x) - x)
    return math.log(45.0 / (4 * math.pi ** 4)) + np.log(g / g_s) + 2 * np.log(x) + ln_kve2(x) - x


@dataclass
class Relic:
    """A relic's specification. ``a`` and ``b`` are in cm^3/s (s-wave and p-wave)."""

    m: float                        # GeV
    a: float = 0.0                  # s-wave <sigma v>, cm^3/s
    b: float = 0.0                  # p-wave coefficient: <sigma v> = a + 6 b / x
    g: float = 2.0                  # internal degrees of freedom
    eos: str = "ss2020"
    #: optional exact thermal average: (ln x grid, ln <sigma v> in GeV^-2), held constant outside
    sigv_table: tuple | None = None

    def sigv_gev(self, x):
        if self.sigv_table is not None:
            lx, ls = self.sigv_table
            return np.exp(np.interp(np.log(np.asarray(x, dtype=float)), lx, ls))
        return (self.a + 6.0 * self.b / np.asarray(x)) / C.GEV2_TO_CM3S


@dataclass
class Solution:
    relic: Relic
    x: np.ndarray
    Y: np.ndarray
    Y_eq: np.ndarray
    omega_h2: float
    Y0: float
    x_f: float                      # where Y = 2 Y_eq (departure Delta = 1)
    x_gamma_eq_h: float             # where n_eq <sigma v> = H
    nfev: int
    extra: dict = field(default_factory=dict)


#: Past x = X_TAIL the equilibrium yield (~e^-150) is irrelevant; the rest is a quadrature.
X_TAIL = 150.0


def solve(relic: Relic, *, x_start: float = C.X_START, T_end: float = C.T_END,
          x_end: float | None = None, nonrel_yeq: bool = False, g_s0: float | None = None,
          rtol: float = C.RTOL, dense: bool = False, T0_K: float = C.T0_K) -> Solution:
    """Integrate the Boltzmann equation for one relic. Returns Y today and Omega h^2.

    The stiff solver runs until x = 150; after that Y_eq^2 / Y^2 < 1e-100 and the equation
    dY/dx = -lambda Y^2 integrates exactly: 1/Y_end = 1/Y(150) + integral lambda dx.
    """
    eos = thermo.load(relic.eos)
    m, g = relic.m, relic.g
    pref = SQRT_PI_45 * C.M_PL * m
    u0 = math.log(x_start)
    u_final = math.log(x_end if x_end is not None else m / T_end)
    u1 = min(u_final, math.log(X_TAIL))

    def lam_x(x):
        T = m / x
        return pref * eos.gstar_half(T) * relic.sigv_gev(x) / (x * x) * x   # x * lambda(x)

    def rhs(u, y):
        x = math.exp(u)
        T = m / x
        ly_eq = float(ln_yeq(x, g, eos.g_s(T), nonrel_yeq))
        xl = float(lam_x(x))
        return [-xl * (math.exp(y[0]) - math.exp(2 * ly_eq - y[0]))]

    def jac(u, y):
        x = math.exp(u)
        T = m / x
        ly_eq = float(ln_yeq(x, g, eos.g_s(T), nonrel_yeq))
        xl = float(lam_x(x))
        return [[-xl * (math.exp(y[0]) + math.exp(2 * ly_eq - y[0]))]]

    y0 = float(ln_yeq(x_start, g, eos.g_s(m / x_start), nonrel_yeq))
    sol = integrate.solve_ivp(rhs, (u0, u1), [y0], method="Radau", jac=jac, rtol=rtol,
                              atol=1e-12, dense_output=dense)
    if not sol.success:
        raise RuntimeError(f"Boltzmann solver failed for m = {m}: {sol.message}")
    x = np.exp(sol.t)
    Y = np.exp(sol.y[0])
    if u_final > u1:
        # exact tail: d(1/Y)/du = x lambda(x)
        ut = np.linspace(u1, u_final, int(max(400, 400 * (u_final - u1))) | 1)
        xt = np.exp(ut)
        xl_t = np.asarray(lam_x(xt), dtype=float)
        inv = 1.0 / Y[-1] + integrate.cumulative_simpson(xl_t, x=ut, initial=0.0)
        keep = slice(1, None, max(1, len(ut) // 200))
        x = np.concatenate([x, xt[keep], xt[-1:]])
        Y = np.concatenate([Y, 1.0 / inv[keep], 1.0 / inv[-1:]])
    T = m / x
    Y_eq = np.exp(ln_yeq(x, g, eos.g_s(T), nonrel_yeq))
    gs0 = eos.g_s0 if g_s0 is None else g_s0
    om = omega_h2_from_Y(m, Y[-1], gs0, T0_K)

    # Diagnostics on a fine grid from the dense solution (or the solver's own steps), x < 100.
    if dense:
        xf_grid = np.exp(np.linspace(u0, min(u1, math.log(100.0)), 4000))
        Yf = np.exp(sol.sol(np.log(xf_grid))[0])
    else:
        sel = x < 100.0
        xf_grid, Yf = x[sel], Y[sel]
    Yeqf = np.exp(ln_yeq(xf_grid, g, eos.g_s(m / xf_grid), nonrel_yeq))
    with np.errstate(divide="ignore", over="ignore"):
        ratio = Yf / Yeqf - 1.0
        x_f = float(np.interp(1.0, ratio, xf_grid)) if np.any(ratio > 1.0) else float("nan")
        # Gamma = n_eq <sigma v>; Gamma/H = s Y_eq <sigma v> / H
        Tg = m / xf_grid
        gam_h = (eos.entropy_density(Tg) * Yeqf * relic.sigv_gev(xf_grid) / eos.hubble(Tg))
        x_geh = float(np.interp(0.0, -np.log(gam_h), xf_grid)) if np.any(gam_h < 1.0) else float("nan")
    out = Solution(relic, x, Y, Y_eq, om, float(Y[-1]), x_f, x_geh, int(sol.nfev))
    if dense:
        out.extra["dense"] = sol.sol
    return out


def omega_h2(relic: Relic, **kw) -> float:
    return solve(relic, **kw).omega_h2


def thermal_cross_section(m: float, *, target: float = C.OMEGA_C_H2, g: float = 2.0, wave: str = "s",
                          eos: str = "ss2020", lo: float = 1e-28, hi: float = 1e-23, **kw) -> dict:
    """The <sigma v> (cm^3/s; p-wave: the coefficient b) that gives Omega h^2 = target."""
    def make(sv):
        return Relic(m, a=sv if wave == "s" else 0.0, b=sv if wave == "p" else 0.0, g=g, eos=eos)

    def f(ln_sv):
        return math.log(omega_h2(make(math.exp(ln_sv)), **kw) / target)

    # Omega ~ 1/<sigma v>: a secant-like bracket from a first guess, then Brent.
    guess = 2.2e-26 if wave == "s" else 1.6e-25
    f0 = f(math.log(guess))
    step = 0.3 if f0 > 0 else -0.3
    a, fa = math.log(guess), f0
    for _ in range(40):
        b = a + step
        fb = f(b)
        if fa * fb <= 0:
            break
        a, fa = b, fb
    else:
        raise RuntimeError(f"no bracket for m = {m}")
    lo_, hi_ = sorted((a, b))
    root = optimize.brentq(f, lo_, hi_, xtol=1e-10, rtol=1e-12)
    sv = math.exp(root)
    sol = solve(make(sv), **kw)
    return {"m": m, "sigv": sv, "x_f": sol.x_f, "x_gamma_eq_h": sol.x_gamma_eq_h, "Y0": sol.Y0,
            "omega_h2": sol.omega_h2}


# --- Analytic approximations ---------------------------------------------------------------


def kolb_turner(m: float, sigv_cm3s: float, g: float = 2.0, eos: str = "ss2020",
                g_s0: float | None = None) -> dict:
    """The textbook estimate (Kolb & Turner 1990, eqs. 5.43-5.45) for s-wave annihilation.

    x_f = ln[0.038 (g/sqrt(g*)) M_Pl m <sigma v>] - (1/2) ln x_f      (c(c+2) = 1)
    Y_inf = x_f / (sqrt(pi/45) M_Pl sqrt(g*) m <sigma v>)
    with g* = g*rho at freeze-out and no derivative term.
    """
    e = thermo.load(eos)
    sv = sigv_cm3s / C.GEV2_TO_CM3S
    x = 20.0
    for _ in range(50):
        gst = float(e.g_rho(m / x))
        arg = 45.0 / (4 * math.pi ** 4) * math.sqrt(math.pi / 2) * SQRT_PI_45 * g / math.sqrt(gst) * C.M_PL * m * sv
        x_new = math.log(arg) - 0.5 * math.log(x)
        if abs(x_new - x) < 1e-12:
            break
        x = x_new
    gst = float(e.g_rho(m / x))
    Y = x / (SQRT_PI_45 * C.M_PL * math.sqrt(gst) * m * sv)
    gs0 = e.g_s0 if g_s0 is None else g_s0
    return {"x_f": x, "g_star": gst, "Y0": Y, "omega_h2": omega_h2_from_Y(m, Y, gs0)}


def improved_analytic(m: float, sigv_cm3s: float, g: float = 2.0, eos: str = "ss2020",
                      b_cm3s: float = 0.0, c: float = math.sqrt(2.0) - 1.0) -> dict:
    """Freeze-out in two steps (Griest & Seckel 1991; Gondolo & Gelmini 1991).

    1. x_f where the departure from equilibrium Delta = Y/Y_eq - 1 reaches c, using
       dY_eq/dx ~ -lambda Delta (Delta + 2) Y_eq^2 with the exact Y_eq and g*^(1/2)(T):
       solve  -dln Y_eq/dx = lambda c (c + 2) Y_eq  (c(c+2) = 1 for c = sqrt 2 - 1).
    2. After that, ignore Y_eq:  1/Y_0 = 1/Y_f + integral_{x_f}^{x_0} lambda(x) dx, with the
       full temperature dependence of g*^(1/2) through the QCD era and e+e- annihilation.
    """
    e = thermo.load(eos)
    relic = Relic(m, a=sigv_cm3s, b=b_cm3s, g=g, eos=eos)

    def lam(x):
        return SQRT_PI_45 * C.M_PL * m * e.gstar_half(m / x) * relic.sigv_gev(x) / (x * x)

    def lyeq(x):
        return float(ln_yeq(x, g, e.g_s(m / x)))

    def cond(x):
        h = 1e-4 * x
        dly = (lyeq(x + h) - lyeq(x - h)) / (2 * h)
        if dly >= 0:      # at small x, Y_eq can still rise as g*s falls: deep in equilibrium
            return -1e3
        return math.log(-dly) - math.log(float(lam(x)) * c * (c + 2)) - lyeq(x)

    x_f = optimize.brentq(cond, 2.0, 300.0, xtol=1e-12)
    Y_f = (1.0 + c) * math.exp(lyeq(x_f))
    x_end = m / C.T_END
    # integrate lambda in ln x
    u = np.linspace(math.log(x_f), math.log(x_end), 6001)
    xs = np.exp(u)
    integral = integrate.simpson(np.asarray(lam(xs)) * xs, x=u)
    Y0 = 1.0 / (1.0 / Y_f + integral)
    return {"x_f": x_f, "Y_f": Y_f, "Y0": Y0, "omega_h2": omega_h2_from_Y(m, Y0, e.g_s0)}


# --- The same equation with a fast, vectorised implicit stepper (also used in the explorer) ----


def yield_grid(masses, a_cm3s, b_cm3s=0.0, g: float = 2.0, eos: str = "ss2020",
               steps_per_efold: int = 800, x_start: float = C.X_START, T_end: float = C.T_END,
               x_end=None, richardson: bool = True, nonrel_yeq: bool = False) -> np.ndarray:
    """Y today for many relics at once (arrays of equal shape), by backward Euler in u = ln x.

    Each step is a quadratic in Y_(n+1), solved exactly:
        Y_(n+1) = Y_n - A (Y_(n+1)^2 - Y_eq^2),   A = h x lambda(x_(n+1))
        ->  Y_(n+1) = 2 c / (1 + sqrt(1 + 4 A c)),  c = Y_n + A Y_eq^2.
    It is first order, L-stable (no oscillation however stiff), and two resolutions combined by
    Richardson extrapolation make it second order. Up to x = 150 all relics share the u grid;
    after that the tail is the exact quadrature used in ``solve``. The interactive explorer runs
    the same scheme in the browser; the tests compare it with the error-controlled Radau solver.
    """
    e = thermo.load(eos)
    m = np.atleast_1d(np.asarray(masses, dtype=float))
    a = np.broadcast_to(np.asarray(a_cm3s, dtype=float), m.shape) / C.GEV2_TO_CM3S
    b = np.broadcast_to(np.asarray(b_cm3s, dtype=float), m.shape) / C.GEV2_TO_CM3S
    pref = SQRT_PI_45 * C.M_PL * m
    u0 = math.log(x_start)
    u_fin = np.log(m / T_end) if x_end is None else np.full(m.shape, math.log(x_end))
    u_mid = min(math.log(X_TAIL), float(u_fin.min()))

    def xlam(x):
        T = m / x
        return pref * e.gstar_half(T) * (a + 6.0 * b / x) / x

    def run(n_per):
        n = int(math.ceil((u_mid - u0) * n_per))
        u = np.linspace(u0, u_mid, n + 1)
        h = u[1] - u[0]
        Y = np.exp(ln_yeq(math.exp(u0), g, e.g_s(m / math.exp(u0)), nonrel_yeq))
        for k in range(1, n + 1):
            x = math.exp(u[k])
            A = h * xlam(x)
            yeq2 = np.exp(2 * ln_yeq(x, g, e.g_s(m / x), nonrel_yeq))
            c = Y + A * yeq2
            Y = 2 * c / (1.0 + np.sqrt(1.0 + 4 * A * c))
        return Y

    Y = run(steps_per_efold)
    if richardson:
        Y = 2 * run(2 * steps_per_efold) - Y
    # exact tail from x_mid to each relic's end, ~100 points per e-fold (sharp g* features need it)
    N = int(max(801, 100 * float((u_fin - u_mid).max()))) | 1
    t = np.linspace(0.0, 1.0, N)
    ut = u_mid + (u_fin - u_mid)[:, None] * t[None, :]
    xt = np.exp(ut)
    T = m[:, None] / xt
    xl = pref[:, None] * e.gstar_half(T) * (a[:, None] + 6.0 * b[:, None] / xt) / xt
    integral = integrate.simpson(xl, x=ut, axis=1)
    return 1.0 / (1.0 / Y + integral)


def omega_grid(masses, a_cm3s, b_cm3s=0.0, g: float = 2.0, eos: str = "ss2020",
               g_s0: float | None = None, **kw) -> np.ndarray:
    m = np.atleast_1d(np.asarray(masses, dtype=float))
    Y = yield_grid(m, a_cm3s, b_cm3s, g=g, eos=eos, **kw)
    gs0 = thermo.load(eos).g_s0 if g_s0 is None else g_s0
    return omega_h2_from_Y(m, Y, gs0)


def thermal_cross_section_grid(masses, *, target: float = C.OMEGA_C_H2, g: float = 2.0,
                               wave: str = "s", eos: str = "ss2020", g_s0: float | None = None,
                               tol: float = 1e-9, **kw) -> np.ndarray:
    """<sigma v> (s-wave) or b (p-wave), in cm^3/s, giving Omega h^2 = target at every mass.

    Newton iteration in ln <sigma v> with a secant slope (Omega ~ x_f / <sigma v>, so the slope
    is close to -1 and a handful of iterations reach 1e-9).
    """
    m = np.atleast_1d(np.asarray(masses, dtype=float))
    target = np.broadcast_to(np.asarray(target, dtype=float), m.shape)

    def om(sv):
        if wave == "s":
            return omega_grid(m, sv, 0.0, g=g, eos=eos, g_s0=g_s0, **kw)
        return omega_grid(m, 0.0, sv, g=g, eos=eos, g_s0=g_s0, **kw)

    ln_sv = np.full(m.shape, math.log(2.2e-26 if wave == "s" else 1.6e-25))
    f = np.log(om(np.exp(ln_sv)) / target)
    slope = np.full(m.shape, -0.96)
    for _ in range(30):
        step = -f / slope
        ln_new = ln_sv + np.clip(step, -3.0, 3.0)
        f_new = np.log(om(np.exp(ln_new)) / target)
        with np.errstate(divide="ignore", invalid="ignore"):
            s_new = (f_new - f) / (ln_new - ln_sv)
        slope = np.where(np.isfinite(s_new) & (s_new < -0.2), s_new, slope)
        ln_sv, f = ln_new, f_new
        if np.max(np.abs(f)) < tol:
            break
    return np.exp(ln_sv)
