"""Inflation, solved exactly (a labelled calculation, not data).

For each potential we integrate the inflaton's equation of motion together
with the Friedmann equation, find where inflation ends, and then follow
individual quantum fluctuations from deep inside the horizon until they freeze
far outside it with the Mukhanov-Sasaki equation. Their frozen amplitudes are
the primordial spectra, P_R(k) and P_t(k), from which n_s, r and the running
follow -- no slow-roll approximation anywhere. The slow-roll formulas are
computed alongside, to see how good they are.

Conventions: reduced Planck mass M_Pl = 1; the potential is V = V0 * v(phi)
and is solved with V0 = 1 (the overall scale cancels in n_s and r and is set at
the end from the measured amplitude A_s). The time variable is the number of
e-folds N = ln a.

Background (exact):
    dphi/dN = pi,   dpi/dN = -(3 - eps) (pi + v'/v),   eps = pi^2/2 = -dlnH/dN,
    H^2 = v / (3 - eps);   inflation ends where eps = 1.

Perturbations of wavenumber k (comoving curvature R, tensor amplitude h):
    R'' + (3 - eps + eps'/eps) R' + (k/aH)^2 R = 0,
    h'' + (3 - eps) h' + (k/aH)^2 h = 0,
starting in the Bunch-Davies vacuum at k/aH = K_START and stopped at
k/aH = K_STOP, long after the mode froze. Then
    P_R = k^3 |R|^2 / (2 pi^2),   P_t = 2 k^3 |h|^2 / (2 pi^2)  (two polarizations).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from numba import njit
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

# Potential kinds (integers, for numba)
POWER, NATURAL, STAROBINSKY, HILLTOP, TMODEL, EXPONENTIAL = range(6)
KIND = {"power": POWER, "natural": NATURAL, "starobinsky": STAROBINSKY,
        "hilltop": HILLTOP, "tmodel": TMODEL, "exponential": EXPONENTIAL}

K_START = 300.0   # start each mode this far inside the horizon (k/aH)
K_STOP = 1e-4     # ... and stop it this far outside
B_STARO = math.sqrt(2.0 / 3.0)


# ---------------------------------------------------------------------------
# Potentials: v(phi), v'(phi), v''(phi) with V0 = 1
# ---------------------------------------------------------------------------


@njit(cache=True)
def potential(kind, p1, p2, phi):
    """Return (v, dv/dphi, d2v/dphi2).

    power        v = phi^p                            (p1 = p)
    natural      v = 1 + cos(phi/f)                   (p1 = f)
    starobinsky  v = (1 - exp(-sqrt(2/3) phi))^2
    hilltop      v = 1 - (phi/mu)^p                   (p1 = p, p2 = mu)
    tmodel       v = tanh^2(phi / sqrt(6 alpha))      (p1 = alpha)
    exponential  v = exp(-lam phi)                    (p1 = lam; power-law inflation, a test)
    """
    if kind == 0:
        p = p1
        return phi ** p, p * phi ** (p - 1.0), p * (p - 1.0) * phi ** (p - 2.0)
    elif kind == 1:
        f = p1
        return 1.0 + math.cos(phi / f), -math.sin(phi / f) / f, -math.cos(phi / f) / f ** 2
    elif kind == 2:
        b = 0.816496580927726  # sqrt(2/3)
        e = math.exp(-b * phi)
        return (1.0 - e) ** 2, 2.0 * b * e * (1.0 - e), 2.0 * b * b * e * (2.0 * e - 1.0)
    elif kind == 3:
        p, mu = p1, p2
        x = phi / mu
        return 1.0 - x ** p, -p * x ** (p - 1.0) / mu, -p * (p - 1.0) * x ** (p - 2.0) / mu ** 2
    elif kind == 4:
        s = math.sqrt(6.0 * p1)
        t = math.tanh(phi / s)
        sech2 = 1.0 - t * t
        return t * t, 2.0 * t * sech2 / s, 2.0 * sech2 * (sech2 - 2.0 * t * t) / s ** 2
    else:
        lam = p1
        e = math.exp(-lam * phi)
        return e, -lam * e, lam * lam * e


# ---------------------------------------------------------------------------
# The mode integrator (numba): background + one scalar and one tensor mode
# ---------------------------------------------------------------------------


@njit(cache=True)
def _rhs(kind, p1, p2, lnk, N, y, out):
    phi, pi = y[0], y[1]
    v, dv, _ = potential(kind, p1, p2, phi)
    eps = 0.5 * pi * pi
    dpi = -(3.0 - eps) * (pi + dv / v)
    lnH = 0.5 * (math.log(v) - math.log(3.0 - eps))
    q = math.exp(lnk - N - lnH)          # k / (aH)
    q2 = q * q
    friction_s = 3.0 - eps + 2.0 * dpi / pi   # 3 - eps + eps'/eps
    friction_t = 3.0 - eps
    out[0] = pi
    out[1] = dpi
    # scalar: R = y2 + i y3, R' = y4 + i y5
    out[2] = y[4]
    out[3] = y[5]
    out[4] = -friction_s * y[4] - q2 * y[2]
    out[5] = -friction_s * y[5] - q2 * y[3]
    # tensor: h = y6 + i y7, h' = y8 + i y9
    out[6] = y[8]
    out[7] = y[9]
    out[8] = -friction_t * y[8] - q2 * y[6]
    out[9] = -friction_t * y[9] - q2 * y[7]


@njit(cache=True)
def _dopri5(kind, p1, p2, lnk, N0, N1, y0, rtol, atol, h0):
    """Adaptive Dormand-Prince 5(4) from N0 to N1. Returns (y, n_steps)."""
    a21 = 1.0 / 5.0
    a31, a32 = 3.0 / 40.0, 9.0 / 40.0
    a41, a42, a43 = 44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0
    a51, a52, a53, a54 = 19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0
    a61, a62, a63, a64, a65 = (9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0,
                               -5103.0 / 18656.0)
    b1, b3, b4, b5, b6 = 35.0 / 384.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0
    e1, e3, e4, e5, e6, e7 = (71.0 / 57600.0, -71.0 / 16695.0, 71.0 / 1920.0,
                              -17253.0 / 339200.0, 22.0 / 525.0, -1.0 / 40.0)
    c2, c3, c4, c5 = 0.2, 0.3, 0.8, 8.0 / 9.0
    n = y0.size
    y = y0.copy()
    k1 = np.empty(n); k2 = np.empty(n); k3 = np.empty(n); k4 = np.empty(n)
    k5 = np.empty(n); k6 = np.empty(n); k7 = np.empty(n); yt = np.empty(n); yn = np.empty(n)
    N = N0
    h = h0
    steps = 0
    _rhs(kind, p1, p2, lnk, N, y, k1)
    while N < N1:
        if N + h > N1:
            h = N1 - N
        for i in range(n):
            yt[i] = y[i] + h * a21 * k1[i]
        _rhs(kind, p1, p2, lnk, N + c2 * h, yt, k2)
        for i in range(n):
            yt[i] = y[i] + h * (a31 * k1[i] + a32 * k2[i])
        _rhs(kind, p1, p2, lnk, N + c3 * h, yt, k3)
        for i in range(n):
            yt[i] = y[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i])
        _rhs(kind, p1, p2, lnk, N + c4 * h, yt, k4)
        for i in range(n):
            yt[i] = y[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i])
        _rhs(kind, p1, p2, lnk, N + c5 * h, yt, k5)
        for i in range(n):
            yt[i] = y[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i]
                                + a65 * k5[i])
        _rhs(kind, p1, p2, lnk, N + h, yt, k6)
        for i in range(n):
            yn[i] = y[i] + h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i])
        _rhs(kind, p1, p2, lnk, N + h, yn, k7)
        err = 0.0
        for i in range(n):
            sc = atol[i] + rtol * max(abs(y[i]), abs(yn[i]))
            e = h * (e1 * k1[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i]
                     + e7 * k7[i]) / sc
            err = max(err, abs(e))
        if err <= 1.0:
            N += h
            for i in range(n):
                y[i] = yn[i]
                k1[i] = k7[i]
            steps += 1
        fac = 0.9 * err ** -0.2 if err > 0 else 5.0
        h *= min(5.0, max(0.2, fac))
    return y, steps


# ---------------------------------------------------------------------------
# Models and their background trajectories
# ---------------------------------------------------------------------------


@dataclass
class Model:
    """One potential with fixed parameters."""
    kind: str
    p1: float = 0.0
    p2: float = 0.0
    name: str = ""

    @property
    def code(self) -> int:
        return KIND[self.kind]

    def v(self, phi):
        return potential(self.code, self.p1, self.p2, float(phi))

    def eps_v(self, phi) -> float:
        v, dv, _ = self.v(phi)
        return 0.5 * (dv / v) ** 2

    def eta_v(self, phi) -> float:
        v, _, d2v = self.v(phi)
        return d2v / v

    def rolls_up(self) -> bool:
        """True when the field rolls toward larger phi (hilltop-type potentials)."""
        return self.kind in ("natural", "hilltop")

    def field_range(self) -> tuple[float, float]:
        """Where the slow-roll end and the start of inflation are searched for.

        The bounds stay clear of points where the slope vanishes in floating point
        (1 + cos -> 0 at pi f, tanh -> 1) yet allow far more than 85 e-folds: near a
        hilltop N grows like mu^4/(8 phi^2) (quartic) or 2 f^2 ln(1/phi) (natural).
        """
        if self.kind == "power":
            return 1e-6, 200.0
        if self.kind == "natural":
            return 1e-6 * self.p1, math.pi * self.p1 * (1 - 1e-6)
        if self.kind == "starobinsky":
            return 1e-6, 40.0
        if self.kind == "hilltop":
            return 1e-3 * self.p2, self.p2 * (1 - 1e-9)
        if self.kind == "tmodel":
            return 1e-6 * math.sqrt(6 * self.p1), 15.0 * math.sqrt(6 * self.p1)
        raise ValueError(self.kind)

    def phi_end_sr(self) -> float:
        """Slow-roll end of inflation: eps_V = 1."""
        lo, hi = self.field_range()
        f = lambda x: self.eps_v(x) - 1.0
        if self.rolls_up():
            # eps_V grows from ~0 at the hilltop toward the minimum
            return brentq(f, lo * 10, hi, xtol=1e-14, rtol=1e-14, maxiter=500)
        return brentq(f, lo, hi, xtol=1e-14, rtol=1e-14, maxiter=500)

    def n_sr(self, phi) -> float:
        """Slow-roll e-folds from phi to the slow-roll end: |int v/v' dphi|."""
        pe = self.phi_end_sr()
        g = lambda x: self.v(x)[0] / self.v(x)[1]
        val, _ = quad(g, pe, phi, limit=400, epsabs=1e-11, epsrel=1e-11)
        return abs(val)

    def phi_at_n_sr(self, n: float) -> float:
        pe = self.phi_end_sr()
        lo, hi = self.field_range()
        f = lambda x: self.n_sr(x) - n
        if self.rolls_up():
            return brentq(f, lo, pe, xtol=1e-15, rtol=1e-13, maxiter=500)
        return brentq(f, pe, hi, xtol=1e-14, rtol=1e-13, maxiter=500)

    def slow_roll(self, n_star: float) -> dict:
        """First-order potential slow-roll predictions N_* e-folds before eps_V = 1."""
        phi = self.phi_at_n_sr(n_star)
        ev, et = self.eps_v(phi), self.eta_v(phi)
        return {"phi": phi, "eps_v": ev, "eta_v": et, "ns": 1 - 6 * ev + 2 * et, "r": 16 * ev}


@dataclass
class Background:
    model: Model
    n_end: float
    sol: object = field(repr=False)

    def state(self, N: float) -> tuple[float, float]:
        phi, pi = self.sol.sol(N)
        return float(phi), float(pi)

    def eps(self, N):
        return 0.5 * self.sol.sol(N)[1] ** 2

    def ln_h(self, N) -> float:
        phi, pi = self.state(N)
        v = self.model.v(phi)[0]
        return 0.5 * (math.log(v) - math.log(3.0 - 0.5 * pi * pi))

    def ln_ah(self, N) -> float:
        return N + self.ln_h(N)


def background(model: Model, n_target: float = 85.0) -> Background:
    """Exact background from ~n_target slow-roll e-folds before the end until eps = 1."""
    phi0 = model.phi_at_n_sr(n_target)
    v, dv, _ = model.v(phi0)
    pi0 = -dv / v  # the slow-roll attractor's velocity; transients die within a few e-folds
    code, p1, p2 = model.code, model.p1, model.p2

    def rhs(N, y):
        v, dv, _ = potential(code, p1, p2, y[0])
        eps = 0.5 * y[1] ** 2
        return [y[1], -(3.0 - eps) * (y[1] + dv / v)]

    def end(N, y):
        return 0.5 * y[1] ** 2 - 1.0
    end.terminal, end.direction = True, 1

    sol = solve_ivp(rhs, (0.0, 4 * n_target), [phi0, pi0], method="DOP853", rtol=1e-12,
                    atol=1e-14, dense_output=True, events=end)
    if sol.status != 1:
        raise RuntimeError(f"{model}: inflation did not end within {4 * n_target} e-folds")
    return Background(model, float(sol.t_events[0][0]), sol)


# ---------------------------------------------------------------------------
# Spectra
# ---------------------------------------------------------------------------


def mode(bg: Background, lnk: float, *, k_start: float = K_START, k_stop: float = K_STOP,
         rtol: float = 1e-10, return_steps: bool = False):
    """Integrate one scalar + one tensor mode; return (ln P_R, ln P_t)."""
    f_start = lambda N: lnk - bg.ln_ah(N) - math.log(k_start)
    f_stop = lambda N: lnk - bg.ln_ah(N) - math.log(k_stop)
    if f_start(0.0) < 0:
        raise ValueError("mode is not deep enough inside the horizon at the start")
    N0 = brentq(f_start, 0.0, bg.n_end, xtol=1e-13)
    N1 = bg.n_end if f_stop(bg.n_end) > 0 else brentq(f_stop, N0, bg.n_end, xtol=1e-12)
    phi0, pi0 = bg.state(N0)
    m = bg.model
    v, dv, _ = m.v(phi0)
    eps0 = 0.5 * pi0 * pi0
    dpi0 = -(3 - eps0) * (pi0 + dv / v)
    q0 = math.exp(lnk - bg.ln_ah(N0))
    # Bunch-Davies, rescaled by sqrt(2k) e^{N0}: R = 1/|pi|, h = 2 (phase set to zero)
    R0 = 1.0 / abs(pi0)
    dR = complex(-(1.0 + dpi0 / pi0), -q0) * R0
    h0 = 2.0
    dh = complex(-1.0, -q0) * h0
    y0 = np.array([phi0, pi0, R0, 0.0, dR.real, dR.imag, h0, 0.0, dh.real, dh.imag])
    atol = np.array([1e-14, 1e-14] + [1e-12 * R0] * 4 + [1e-12 * h0] * 4)
    y, steps = _dopri5(m.code, m.p1, m.p2, lnk, N0, N1, y0, rtol, atol, 1e-4)
    base = 2 * lnk - 2 * N0
    ln_pr = base - math.log(4 * math.pi ** 2) + math.log(y[2] ** 2 + y[3] ** 2)
    ln_pt = base - math.log(2 * math.pi ** 2) + math.log(y[6] ** 2 + y[7] ** 2)
    if return_steps:
        return ln_pr, ln_pt, steps
    return ln_pr, ln_pt


def spectra_at(bg: Background, n_star: float, *, offsets=(-0.5, -0.25, 0.0, 0.25, 0.5),
               extra_offsets: tuple = ()) -> dict:
    """n_s, r, running and amplitudes at the pivot that exits n_star e-folds before the end.

    The pivot wavenumber is k_* = a H at N_end - n_star (in the solver's units);
    neighbouring modes at ln(k/k_*) = offsets give the slopes by a quadratic fit.
    """
    n_exit = bg.n_end - n_star
    lnk_star = bg.ln_ah(n_exit)
    d = np.asarray(offsets, float)
    out = np.array([mode(bg, lnk_star + x) for x in d])
    lpr, lpt = out[:, 0], out[:, 1]
    cs = np.polyfit(d, lpr, 2)
    ct = np.polyfit(d, lpt, 2)
    i0 = int(np.argmin(np.abs(d)))
    res = {
        "n_star": n_star, "lnk_star": lnk_star,
        "ns": 1 + cs[1], "alpha_s": 2 * cs[0], "nt": ct[1],
        "r": math.exp(lpt[i0] - lpr[i0]), "ln_pr": lpr[i0], "ln_pt": lpt[i0],
    }
    for x in extra_offsets:
        a, b = mode(bg, lnk_star + x)
        res[f"r_at_{x:+.3f}"] = math.exp(b - a)
        res[f"ln_pr_at_{x:+.3f}"] = a
    phi, pi = bg.state(n_exit)
    res.update(phi_star=phi, pi_star=pi, phi_end=bg.state(bg.n_end)[0],
               ln_h_star=bg.ln_h(n_exit), ln_h_end=bg.ln_h(bg.n_end),
               v_star=bg.model.v(phi)[0], v_end=bg.model.v(bg.state(bg.n_end)[0])[0])
    res.update(hubble_flow(bg, n_exit))
    return res


def hubble_flow(bg: Background, N: float, h: float = 1e-3) -> dict:
    """Hubble-flow functions eps1 = -dlnH/dN, eps2 = dln eps1/dN, eps3 = dln eps2/dN."""
    def e1(n):
        return bg.eps(n)

    def e2(n):
        return (math.log(e1(n + h)) - math.log(e1(n - h))) / (2 * h)

    eps1, eps2 = e1(N), e2(N)
    eps3 = (math.log(abs(e2(N + h))) - math.log(abs(e2(N - h)))) / (2 * h)
    return {"eps1": eps1, "eps2": eps2, "eps3": eps3}


C_SR = 0.5772156649015329 + math.log(2.0) - 2.0  # Euler-Mascheroni + ln 2 - 2


def second_order(eps1: float, eps2: float, eps3: float) -> dict:
    """Second-order slow-roll spectra in Hubble-flow functions at horizon crossing
    (Schwarz, Terrero-Escalante & Garcia 2001; Martin, Ringeval & Vennin 2014)."""
    ns = 1 - 2 * eps1 - eps2 - 2 * eps1 ** 2 - (2 * C_SR + 3) * eps1 * eps2 - C_SR * eps2 * eps3
    r = 16 * eps1 * (1 + C_SR * eps2)
    return {"ns": ns, "r": r}


# ---------------------------------------------------------------------------
# Physical scales: energy of inflation and e-folds from reheating
# ---------------------------------------------------------------------------


def energy_scale(res: dict, ln10_10_as: float, m_pl_gev: float) -> dict:
    """Fix V0 by the measured amplitude A_s and return physical scales in GeV."""
    a_s = math.exp(ln10_10_as) * 1e-10
    v0 = a_s / math.exp(res["ln_pr"])            # V0 in units of M_Pl^4
    h_star = math.sqrt(v0) * math.exp(res["ln_h_star"]) * m_pl_gev
    h_end = math.sqrt(v0) * math.exp(res["ln_h_end"]) * m_pl_gev
    return {"V0": v0,
            "V_star_quarter_gev": (v0 * res["v_star"]) ** 0.25 * m_pl_gev,
            "V_end_quarter_gev": (v0 * res["v_end"]) ** 0.25 * m_pl_gev,
            "H_star_gev": h_star, "H_end_gev": h_end,
            "rho_end_gev4": 3.0 * (h_end / m_pl_gev) ** 2 * m_pl_gev ** 4,
            "delta_phi": abs(res["phi_star"] - res["phi_end"])}


def n_star_from_reheating(h_star_gev: float, rho_end_gev4: float, *, k_star_mpc: float,
                          h0_km_s_mpc: float, t0_k: float, g_s0: float, w: float = 0.0,
                          rho_reh_gev4: float | None = None, g_reh: float = 106.75,
                          mpc_inv_gev: float, kb_gev_per_k: float, c_km_s: float) -> float:
    """e-folds between the pivot's horizon exit and the end of inflation, for a given
    reheating history: equation of state w from the end of inflation until the
    universe thermalises at density rho_reh (instant reheating: rho_reh = rho_end).

        N_* = -ln(k_*/a0H0) + ln(H_*/H0) - ln(rho_end/rho_reh)/(3(1+w))
              + ln(T0/T_reh) + (1/3) ln(g_s0/g_reh)
    """
    if rho_reh_gev4 is None:
        rho_reh_gev4 = rho_end_gev4
    h0_gev = h0_km_s_mpc / c_km_s * mpc_inv_gev
    t_reh = (30.0 * rho_reh_gev4 / (math.pi ** 2 * g_reh)) ** 0.25
    t0 = t0_k * kb_gev_per_k
    k_gev = k_star_mpc * mpc_inv_gev
    return (-math.log(k_gev / h0_gev) + math.log(h_star_gev / h0_gev)
            - math.log(rho_end_gev4 / rho_reh_gev4) / (3.0 * (1.0 + w))
            + math.log(t0 / t_reh) + math.log(g_s0 / g_reh) / 3.0)


# ---------------------------------------------------------------------------
# Running a model: predictions over N_* and the reheating e-fold count
# ---------------------------------------------------------------------------


def model_from_spec(spec: dict, param: float | None = None) -> Model:
    kind = spec["kind"]
    if kind == "power":
        return Model("power", spec["p"])
    if kind == "natural":
        return Model("natural", param)
    if kind == "starobinsky":
        return Model("starobinsky")
    if kind == "hilltop":
        return Model("hilltop", spec["p"], param)
    if kind == "tmodel":
        return Model("tmodel", param)
    raise ValueError(kind)


def amplitude_sr(bg: Background, n_star: float) -> tuple[float, float]:
    """Slow-roll amplitude H^2/(8 pi^2 eps) at horizon exit (V0 = 1) and ln H there;
    only used to normalise V0 inside the reheating iteration (1% accuracy suffices:
    it moves N_* by ~0.003)."""
    n = bg.n_end - n_star
    eps = bg.eps(n)
    ln_h = bg.ln_h(n)
    return math.exp(2 * ln_h) / (8 * math.pi ** 2 * eps), ln_h


def n_star_reheating(bg: Background, *, ln10_10_as: float, m_pl_gev: float, w: float,
                     t_reh_gev: float | None, consts: dict, g_reh: float = 106.75) -> float:
    """Self-consistent N_* for a reheating history (instant if t_reh_gev is None)."""
    n = 55.0
    a_s = math.exp(ln10_10_as) * 1e-10
    v_end = bg.model.v(bg.state(bg.n_end)[0])[0]
    for _ in range(50):
        p_hat, ln_h = amplitude_sr(bg, n)
        v0 = a_s / p_hat
        h_star = math.sqrt(v0) * math.exp(ln_h) * m_pl_gev
        h_end2 = v0 * v_end / 2.0          # H^2 = V/(3 - eps) with eps = 1 at the end
        rho_end = 3.0 * h_end2 * m_pl_gev ** 4
        rho_reh = None if t_reh_gev is None else math.pi ** 2 / 30 * g_reh * t_reh_gev ** 4
        if rho_reh is not None and rho_reh > rho_end:
            return float("nan")  # cannot reheat hotter than the end of inflation
        new = n_star_from_reheating(h_star, rho_end, w=w, rho_reh_gev4=rho_reh,
                                    g_reh=g_reh, **consts)
        if abs(new - n) < 1e-6:
            return new
        n = new
    return n


def predictions(model: Model, n_stars, *, ln10_10_as: float, m_pl_gev: float,
                extra_offsets: tuple = ()) -> list[dict]:
    """Exact spectra and physical scales for each N_*, plus the slow-roll estimates.

    The background is solved from far enough back that every mode (including
    the offsets, e.g. r at 0.002 Mpc^-1, 3.2 e-folds earlier) starts 300 times
    inside the horizon after the initial transients have died away."""
    bg = background(model, n_target=max(85.0, float(np.max(n_stars)) + 15.0))
    out = []
    for n in n_stars:
        res = spectra_at(bg, float(n), extra_offsets=extra_offsets)
        res.update(energy_scale(res, ln10_10_as, m_pl_gev))
        sr = model.slow_roll(float(n))
        res["ns_sr1"], res["r_sr1"] = sr["ns"], sr["r"]
        so = second_order(res["eps1"], res["eps2"], res["eps3"])
        res["ns_sr2"], res["r_sr2"] = so["ns"], so["r"]
        out.append(res)
    return out


def mode_history(bg: Background, lnk: float, n_out: int = 600, *, k_start: float = K_START,
                 k_stop: float = 1e-3) -> dict:
    """One scalar and one tensor mode recorded along the way (for the figures):
    e-folds N, k/aH, |R| and |h| in units where the final values give P_R and P_t."""
    f_start = lambda N: lnk - bg.ln_ah(N) - math.log(k_start)
    f_stop = lambda N: lnk - bg.ln_ah(N) - math.log(k_stop)
    N0 = brentq(f_start, 0.0, bg.n_end, xtol=1e-13)
    N1 = bg.n_end if f_stop(bg.n_end) > 0 else brentq(f_stop, N0, bg.n_end, xtol=1e-12)
    phi0, pi0 = bg.state(N0)
    m = bg.model
    v, dv, _ = m.v(phi0)
    eps0 = 0.5 * pi0 * pi0
    dpi0 = -(3 - eps0) * (pi0 + dv / v)
    q0 = math.exp(lnk - bg.ln_ah(N0))
    R0 = 1.0 / abs(pi0)
    dR = complex(-(1.0 + dpi0 / pi0), -q0) * R0
    dh = complex(-1.0, -q0) * 2.0
    y = np.array([phi0, pi0, R0, 0.0, dR.real, dR.imag, 2.0, 0.0, dh.real, dh.imag])
    atol = np.array([1e-14, 1e-14] + [1e-12 * R0] * 4 + [2e-12] * 4)
    grid = np.linspace(N0, N1, n_out)
    out = np.empty((n_out, 5))
    # P = k^2 |X|^2 / (c pi^2 e^{2 N0}) with c = 4 (scalar) or 2 (tensor): record sqrt(P)
    pref = math.exp(lnk - N0)
    for i, Nn in enumerate(grid):
        if i > 0:
            y, _ = _dopri5(m.code, m.p1, m.p2, lnk, grid[i - 1], Nn, y, 1e-10, atol, 1e-4)
        out[i] = (Nn, math.exp(lnk - bg.ln_ah(Nn)),
                  pref * math.hypot(y[2], y[3]) / (2 * math.pi),
                  pref * math.hypot(y[6], y[7]) / (math.sqrt(2) * math.pi),
                  pref * y[2] / (2 * math.pi))
    return {"N": out[:, 0], "k_over_aH": out[:, 1], "sqrt_PR": out[:, 2], "sqrt_Pt": out[:, 3],
            "re_sqrt_PR": out[:, 4]}
