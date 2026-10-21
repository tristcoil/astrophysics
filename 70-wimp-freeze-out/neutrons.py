"""The same physics one second after the Big Bang: neutrons freeze out, decay, and end in helium.

Above ~1 MeV the weak reactions

    n + nu_e <-> p + e-,    n + e+ <-> p + anti-nu_e,    n <-> p + e- + anti-nu_e

keep the neutron-to-proton ratio at its equilibrium value n/p = exp(-Q/T), Q = m_n - m_p =
1.293 MeV. Their rates fall as T^5 while the expansion rate falls as T^2, so near T ~ 0.8 MeV
they lose the race and n/p freezes at ~1/6, exactly as a WIMP's abundance freezes when
n <sigma v> drops below H. Afterwards free neutrons decay (lifetime 878 s) until deuterium
survives the photons (the "deuterium bottleneck", T ~ 0.08 MeV, t ~ 3 min); then almost every
neutron ends up in helium-4, so Y_p ~ 2 (n/p) / (1 + n/p) ~ 0.25.

What is computed here:

* the six weak rates in the Born approximation (V-A theory, infinitely heavy nucleons, Fermi-Dirac
  electrons and neutrinos at their own temperatures), normalised so that the free-neutron decay
  rate equals 1/tau_n (Weinberg 1972; Bernstein, Brown & Feinberg 1989);
* the background: photons, e+- (exact Fermi-Dirac thermodynamics) and three neutrino species that
  decouple instantly, so e+- annihilation heats the photons by (11/4)^(1/3); or, alternatively,
  Escudero's (2020) non-instantaneous decoupling history (N_eff = 3.045);
* the neutron's fate: the n <-> p rate equation, plus the eight nuclear reactions that make
  helium-4 out of n, p, d, t and helium-3, with rates from the PRIMAT compilation (Pitrou et al.
  2018) and reverse rates from detailed balance with masses and spins.

This is a minimal network on purpose (no lithium, Born-level weak rates). K1 will run full BBN; here
the point is the freeze-out. PRyMordial (Burns, Tait & Valli 2024) run with the same simplifications
checks it, and run in full shows what the neglected corrections are worth.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from scipy import integrate, interpolate, optimize

import config as C
import thermo

ME = C.M_E
QN = C.Q_NP
QE = QN / ME                                  # q = Q / m_e = 2.531

# --- Born weak rates ------------------------------------------------------------------------


def _fd(e: float) -> float:
    """Fermi-Dirac occupation 1/(e^E + 1), E in units of the temperature: e^-E/(1 + e^-E) for E > 0.

    (0.5 (1 - tanh(E/2)) is the same function but cancels catastrophically for large E, which made
    the ~e^-60 proton-to-neutron rates at low temperature mostly rounding noise.)"""
    if e > 0:
        x = math.exp(-e)
        return x / (1.0 + x)
    return 1.0 / (1.0 + math.exp(e))


@lru_cache(maxsize=None)
def lambda0() -> float:
    """The free-decay phase-space integral  int_1^q eps sqrt(eps^2-1) (q-eps)^2 d eps  (~1.636)."""
    return integrate.quad(lambda e: e * math.sqrt(e * e - 1) * (QE - e) ** 2, 1.0, QE, epsabs=0,
                          epsrel=1e-13)[0]


def weak_rates(T_gamma: float, T_nu: float, tau_n: float = C.TAU_N) -> tuple[float, float]:
    """(lambda_n->p, lambda_p->n) in 1/s at photon temperature T_gamma and neutrino T_nu (GeV).

    Energies in units of m_e; eps is the electron or positron energy (>= 1), p = sqrt(eps^2 - 1).
    n -> p:  n nu -> p e   [nu energy eps - q]      n e+ -> p nubar [nubar energy eps + q]
             n -> p e nubar [nubar energy q - eps]
    p -> n:  p e -> n nu   [nu energy eps - q]      p nubar -> n e+ [nubar energy eps + q]
             p e nubar -> n [nubar energy q - eps]
    """
    z, zn = ME / T_gamma, ME / T_nu
    q = QE
    fe = (lambda e: float(_fd(e * z)))
    fn = (lambda e: float(_fd(e * zn)))
    scale = 1.0 / min(z, zn)              # the thermal energy in units of m_e
    top = q + 80.0 * scale + 60.0         # far into the Boltzmann tail

    def qd(g, a, b):
        """integral of eps sqrt(eps^2 - 1) g(eps) d eps from a to b. With eps = cosh(th) the
        threshold's square root disappears: eps p d eps = cosh(th) sinh(th)^2 d th. Relative
        tolerance only (at low T the p -> n integrals are ~e^-60 and an absolute one would return
        noise); the range is split at a few thermal energies."""
        edges = [a] + [a + k * scale for k in (2.0, 8.0, 25.0) if a + k * scale < b] + [b]
        th = [math.acosh(e) for e in edges]

        def f(t):
            e = math.cosh(t)
            return e * math.sinh(t) ** 2 * g(e)

        total = 0.0
        for lo, hi in zip(th[:-1], th[1:]):
            # far-tail pieces add ~1e-11 of the total: 1e-12 of what is already in is plenty
            total += integrate.quad(f, lo, hi, limit=200, epsabs=1e-12 * abs(total), epsrel=1e-10)[0]
        return total

    np_ = (qd(lambda e: (e - q) ** 2 * fn(e - q) * (1 - fe(e)), q, top)
           + qd(lambda e: (e + q) ** 2 * fe(e) * (1 - fn(e + q)), 1.0, top)
           + qd(lambda e: (q - e) ** 2 * (1 - fe(e)) * (1 - fn(q - e)), 1.0, q))
    pn = (qd(lambda e: (e - q) ** 2 * fe(e) * (1 - fn(e - q)), q, top)
          + qd(lambda e: (e + q) ** 2 * fn(e + q) * (1 - fe(e)), 1.0, top)
          + qd(lambda e: (q - e) ** 2 * fe(e) * fn(q - e), 1.0, q))
    norm = 1.0 / (tau_n * lambda0())
    return np_ * norm, pn * norm


# --- Backgrounds --------------------------------------------------------------------------------


@dataclass
class Background:
    """Photon temperature, neutrino temperature, expansion rate and time on a common grid.

    ``nb_factor`` = n_b / (eta_0 n_gamma(T_gamma)), with eta_0 today's baryon-to-photon ratio. While
    e+e- pairs annihilate, photons are still being added to every comoving volume, so the baryon
    density is higher than eta_0 n_gamma: by (2 + g_s,e)/2 = 2.7 at 1 MeV and still 3% at 70 keV,
    when deuterium burns. Defaults to that photon + e+- entropy factor.
    """

    name: str
    t: np.ndarray            # s, increasing
    T_gamma: np.ndarray      # GeV
    T_nu: np.ndarray         # GeV
    H: np.ndarray            # 1/s
    N_eff: float
    extra: dict = field(default_factory=dict)
    nb_factor: np.ndarray | None = None

    def __post_init__(self):
        if self.nb_factor is None:
            self.nb_factor = (2.0 + _electron_g(self.T_gamma)[1]) / 2.0
        lt = np.log(self.t)
        self._Tg = interpolate.CubicSpline(lt, np.log(self.T_gamma))
        self._Tn = interpolate.CubicSpline(lt, np.log(self.T_nu))
        self._H = interpolate.CubicSpline(lt, np.log(self.H))
        self._nbf = interpolate.PchipInterpolator(lt, np.log(self.nb_factor))
        lT = np.log(self.T_gamma[::-1])
        self._t_of_T = interpolate.CubicSpline(lT, lt[::-1])

    def baryon_factor(self, t):
        return np.exp(self._nbf(np.log(t)))

    def Tg(self, t):
        return np.exp(self._Tg(np.log(t)))

    def Tn(self, t):
        return np.exp(self._Tn(np.log(t)))

    def hubble(self, t):
        return np.exp(self._H(np.log(t)))

    def t_of_Tg(self, T):
        return np.exp(self._t_of_T(np.log(T)))


def _electron_g(T):
    """(g_rho, g_s) of e+e- (4 states) at temperature T, exact Fermi-Dirac."""
    return thermo._species_g(ME, 4, "F", np.asarray(T, dtype=float))


@lru_cache(maxsize=None)
def instant_background(N_eff: float = 3.0, T_start: float = 10e-3, T_end: float = 5e-6,
                       n: int = 3000) -> Background:
    """Neutrinos decouple instantly at T_start (e+- still relativistic), then T_nu ~ 1/a, while the
    photon + e+- plasma keeps its own entropy: (a T_gamma)^3 (2 + g_s,e(T_gamma)) = const.
    N_eff != 3 scales the neutrino energy density (the classic 'speed-up' test)."""
    gs_start = 2.0 + _electron_g(T_start)[1]
    Tn = np.geomspace(T_start, T_end * 0.25, n)
    # solve (T_gamma/T_nu)^3 (2 + g_s,e(T_gamma)) = gs_start for each T_nu
    zgrid = np.geomspace(1e-3, 1e4, 6000)                        # m_e / T_gamma
    gse = 2.0 + _electron_g(ME / zgrid)[1]
    # along a monotonic path: x = T_gamma/T_nu = (gs_start/gse)^(1/3), T_nu = T_gamma/x
    x_of_z = (gs_start / gse) ** (1.0 / 3.0)
    Tn_of_z = (ME / zgrid) / x_of_z
    order = np.argsort(Tn_of_z)
    Tg = np.exp(np.interp(np.log(Tn), np.log(Tn_of_z[order]), np.log(ME / zgrid[order])))
    grho_e = _electron_g(Tg)[0]
    rho = math.pi ** 2 / 30 * ((2.0 + grho_e) * Tg ** 4 + N_eff * 2 * 7 / 8 * Tn ** 4)
    H = np.sqrt(rho / 3.0) / C.M_PL_RED / C.HBAR_GEV_S          # 1/s
    # t: dt = -dln T_nu / H, starting from t = 1/(2H) (radiation era)
    lnTn = np.log(Tn)
    t = 1.0 / (2.0 * H[0]) + integrate.cumulative_simpson(1.0 / H, x=-lnTn, initial=0.0)
    keep = Tg >= T_end
    return Background(f"instant_Neff{N_eff:g}", t[keep], Tg[keep], Tn[keep], H[keep], N_eff,
                      {"Tg_over_Tnu_final": float(Tg[-1] / Tn[-1])})


@lru_cache(maxsize=None)
def escudero_background() -> Background:
    """Escudero (2020): non-instantaneous neutrino decoupling, 10 MeV -> 1 keV (N_eff = 3.045)."""
    d = np.loadtxt(C.RAW_DIR / "eos" / "escudero2020_sm_thermo.dat")
    t, Tg_mev, ratio, gstar, z_gam = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 6]
    Tg = Tg_mev * 1e-3
    Tn = Tg / ratio
    H = np.sqrt(math.pi ** 2 / 90 * gstar) * Tg ** 2 / C.M_PL_RED / C.HBAR_GEV_S
    # n_b ~ a^-3 and z_gamma = a T_gamma / m_e, so n_b / (eta_0 n_gamma) = (z_gamma,end / z_gamma)^3
    return Background("escudero2020", t, Tg, Tn, H, 3.045, nb_factor=(z_gam[-1] / z_gam) ** 3)


# --- Nuclear network ------------------------------------------------------------------------------

#: NUBASE2020 atomic mass excesses (keV) and spins, as in PRyMordial.
MASS_EXCESS_KEV = {"n": 8071.3171, "p": 7288.9706, "d": 13135.722, "t": 14949.81, "h": 14931.218,
                   "a": 2424.9156}
SPIN = {"n": 0.5, "p": 0.5, "d": 1.0, "t": 0.5, "h": 0.5, "a": 0.0}
A_NUM = {"n": 1, "p": 1, "d": 2, "t": 3, "h": 3, "a": 4}
SPECIES = ("n", "p", "d", "t", "h", "a")      # h = helium-3, a = helium-4
N_AVOGADRO = 6.02214076e23
T9_PER_GEV = 1.0 / (C.K_B_GEV_PER_K * 1e9)    # T9 = T[GeV] * this


Z_NUM = {"n": 0, "p": 1, "d": 1, "t": 1, "h": 2, "a": 2}


def mass_u(sp: str) -> float:
    """Nuclear mass in u: atomic mass (A + mass excess) minus Z electrons."""
    return A_NUM[sp] + MASS_EXCESS_KEV[sp] / (C.M_U * 1e6) - Z_NUM[sp] * C.M_E / C.M_U


def q_value_mev(ins, outs) -> float:
    return sum(MASS_EXCESS_KEV[s] for s in ins) / 1e3 - sum(MASS_EXCESS_KEV[s] for s in outs) / 1e3


#: (key, inputs, outputs). Inputs/outputs are tuples of species; 'g' = photon (not tracked).
REACTIONS = [
    ("npdg", ("n", "p"), ("d",)),
    ("dpHe3g", ("d", "p"), ("h",)),
    ("ddHe3n", ("d", "d"), ("h", "n")),
    ("ddtp", ("d", "d"), ("t", "p")),
    ("He3ntp", ("h", "n"), ("t", "p")),
    ("tdan", ("t", "d"), ("a", "n")),
    ("He3dap", ("h", "d"), ("a", "p")),
    ("tpag", ("t", "p"), ("a",)),
]


def reverse_coefficients(ins, outs) -> tuple[float, float, float]:
    """Detailed balance: <sigma v>_reverse / <sigma v>_forward = alpha T9^beta exp(gamma / T9).

    Saha: (prod n_in / prod n_out)_eq = (prod g_in m_in^(3/2) / prod g_out m_out^(3/2))
          * (m_u T / 2 pi)^(3/2 (N_in - N_out)) e^(-Q/T), and with identical particles in a channel
    the event rate carries 1/2 (n_d^2 <sigma v> / 2), which multiplies the reverse coefficient.
    Radiative captures (2 -> 1 + photon) give the photodisintegration rate in 1/s per nucleus.
    """
    g = (lambda s: 2 * SPIN[s] + 1)
    num = np.prod([g(s) * mass_u(s) ** 1.5 for s in ins])
    den = np.prod([g(s) * mass_u(s) ** 1.5 for s in outs])
    ratio = num / den
    sym_in = 0.5 if len(ins) == 2 and ins[0] == ins[1] else 1.0
    sym_out = 0.5 if len(outs) == 2 and outs[0] == outs[1] else 1.0
    ratio *= sym_in / sym_out
    dn = len(ins) - len(outs)
    beta = 1.5 * dn
    if dn:
        # (m_u k T / 2 pi hbar^2)^(3/2) at T9 = 1, in cm^-3, divided by N_A (rates are N_A<sigma v>)
        kT = C.K_B_GEV_PER_K * 1e9
        x = (C.M_U * kT / (2 * math.pi)) ** 1.5 / C.HBARC_GEV_CM ** 3
        ratio *= x ** dn / N_AVOGADRO ** dn
    gamma = -q_value_mev(ins, outs) * 1e-3 * T9_PER_GEV
    return float(ratio), float(beta), float(gamma)


@lru_cache(maxsize=None)
def primat_rate(key: str):
    """N_A <sigma v> (cm^3 mol^-1 s^-1) as a function of T9, log-log spline of PRIMAT's table."""
    d = np.loadtxt(C.RAW_DIR / "bbn" / f"primat_{key}.txt")
    return interpolate.CubicSpline(np.log(d[:, 0]), np.log(d[:, 1]))


def eta_from_omega_b(omega_b_h2: float = C.OMEGA_B_H2, Y_p: float = 0.247) -> float:
    """eta = n_b / n_gamma today, with the mean mass per baryon of H + He-4 (as PRyMordial)."""
    m_b = ((1 - Y_p) * 1.00782503223 + Y_p * 4.0026032541 / 4.0) * C.M_U      # GeV
    n_gamma = 2 * C.ZETA3 / math.pi ** 2 * C.T0_GEV ** 3 / C.HBARC_GEV_CM ** 3  # cm^-3
    return omega_b_h2 * C.RHO_CRIT_H2_GEV_CM3 / m_b / n_gamma


@dataclass
class BBNResult:
    t: np.ndarray
    T_gamma: np.ndarray
    Y: dict                     # species -> abundance per baryon
    X_n_free: np.ndarray        # free-neutron fraction n/(n+p) among free nucleons
    Y_p: float                  # helium-4 mass fraction (4 Y_a)
    D_H: float
    He3_H: float
    np_ratio_1s: float
    T_freeze: float             # GeV, where (lambda_np + lambda_pn) = H
    t_freeze: float
    t_nuc: float                # s, when half of the final helium has formed
    T_nuc: float
    X_n_at_nuc: float
    meta: dict = field(default_factory=dict)


def run(eta: float | None = None, tau_n: float = C.TAU_N, background: Background | None = None,
        nuclei: bool = True, n_rate_points: int = 220, T_stop: float = 8e-6,
        T_net: float = 1.0e-3) -> BBNResult:
    """Integrate free neutrons from 10 MeV and, if ``nuclei``, the 8-reaction network from T_net."""
    bg = background or instant_background()
    eta = eta_from_omega_b() if eta is None else eta
    # weak rates tabulated along the background
    Tg_tab = np.geomspace(bg.T_gamma[0], max(bg.T_gamma[-1], T_stop), n_rate_points)
    t_tab = bg.t_of_Tg(Tg_tab)
    Tn_tab = bg.Tn(t_tab)
    rates = np.array([weak_rates(a, b, tau_n) for a, b in zip(Tg_tab, Tn_tab)])
    lt_tab = np.log(t_tab)
    # Monotone interpolation in log space: lambda_p->n falls by hundreds of e-folds at low T and
    # a cubic spline would overshoot (which once created neutrons out of protons at 20 keV).
    s_np = interpolate.PchipInterpolator(lt_tab, np.log(rates[:, 0]))
    s_pn = interpolate.PchipInterpolator(lt_tab, np.log(np.maximum(rates[:, 1], 1e-300)))

    def lam(t):
        lt = math.log(t)
        return math.exp(s_np(lt)), math.exp(s_pn(lt))

    rev = {k: reverse_coefficients(i, o) for k, i, o in REACTIONS}
    spl = {k: primat_rate(k) for k, _, _ in REACTIONS}
    idx = {s: i for i, s in enumerate(SPECIES)}
    # n_b = nb_pref * T_gamma^3 * (baryon factor of the background)
    nb_pref = eta * 2 * C.ZETA3 / math.pi ** 2 / C.HBARC_GEV_CM ** 3
    net_def = [(k, [idx[s] for s in i], [idx[s] for s in o],
                0.5 if len(i) == 2 and i[0] == i[1] else 1.0) for k, i, o in REACTIONS]

    def _prod_except(Y, ids, j):
        p = 1.0
        for m_, s in enumerate(ids):
            if m_ != j:
                p *= Y[s]
        return p

    def terms(u, Y):
        """dY/du and its Jacobian (analytic: finite differences fail on 1e-30 abundances)."""
        t = math.exp(u)
        Tg = float(bg.Tg(t))
        l_np, l_pn = lam(t)
        d = np.zeros(6)
        J = np.zeros((6, 6))
        wn = l_pn * Y[1] - l_np * Y[0]
        d[0], d[1] = wn, -wn
        J[0, 0], J[0, 1], J[1, 0], J[1, 1] = -l_np, l_pn, l_np, -l_pn
        if nuclei:
            T9 = Tg * T9_PER_GEV
            lT9 = math.log(T9)
            nb = nb_pref * Tg ** 3 * float(bg.baryon_factor(t))
            for k, ins, outs, sym in net_def:
                fwd = math.exp(spl[k](lT9)) / N_AVOGADRO          # <sigma v>, cm^3/s
                a, b, c = rev[k]
                arg = c / T9
                rv = a * T9 ** b * math.exp(arg) * fwd if arg > -700 else 0.0
                Kf = sym * nb * fwd                                # events per baryon per second
                Kr = rv * N_AVOGADRO if len(outs) == 1 else nb * rv  # photodisintegration per nucleus
                net = Kf * _prod_except(Y, ins, -1) - Kr * _prod_except(Y, outs, -1)
                dnet = np.zeros(6)
                for j, s in enumerate(ins):
                    dnet[s] += Kf * _prod_except(Y, ins, j)
                for j, s in enumerate(outs):
                    dnet[s] -= Kr * _prod_except(Y, outs, j)
                for s in ins:
                    d[s] -= net
                    J[s] -= dnet
                for s in outs:
                    d[s] += net
                    J[s] += dnet
        return d * t, J * t

    def terms_np(u, Y):
        t = math.exp(u)
        l_np, l_pn = lam(t)
        wn = l_pn * Y[1] - l_np * Y[0]
        return np.array([wn, -wn]) * t, np.array([[-l_np, l_pn], [l_np, -l_pn]]) * t

    def nse(u, Yn, Yp):
        """Light nuclei in statistical equilibrium with n and p (the network's own fixed point:
        each abundance follows from the detailed-balance coefficients of the chain that makes it)."""
        Tg = float(bg.Tg(math.exp(u)))
        T9 = Tg * T9_PER_GEV
        nb = nb_pref * Tg ** 3 * float(bg.baryon_factor(math.exp(u)))

        def photo(k):                    # (photodisintegration rate) / (n_b <sigma v>_fwd)
            a, b, c = rev[k]
            return a * T9 ** b * math.exp(c / T9) * N_AVOGADRO / nb

        def swap(k):                     # 2 -> 2: <sigma v>_rev / <sigma v>_fwd
            a, b, c = rev[k]
            return a * T9 ** b * math.exp(c / T9)

        Yd = Yn * Yp / photo("npdg")
        Yh = Yd * Yp / photo("dpHe3g")
        Yt = Yh * Yn / (Yp * swap("He3ntp"))
        Ya = Yt * Yp / photo("tpag")
        return np.array([Yn, Yp, Yd, Yt, Yh, Ya])

    t0 = float(bg.t[0])
    l_np0, l_pn0 = lam(t0)
    Xn0 = l_pn0 / (l_np0 + l_pn0)
    t_end = float(bg.t_of_Tg(T_stop))
    u0, u_end = math.log(t0), math.log(t_end)
    # Stage 1: free nucleons only (light nuclei are ~1e-12 and irrelevant above ~1 MeV).
    u_net = min(math.log(float(bg.t_of_Tg(T_net))), u_end) if nuclei else u_end
    s1 = integrate.solve_ivp(lambda u, y: terms_np(u, y)[0], (u0, u_net), [Xn0, 1 - Xn0],
                             method="Radau", jac=lambda u, y: terms_np(u, y)[1], rtol=1e-10,
                             atol=1e-14, dense_output=True)
    if not s1.success:
        raise RuntimeError(s1.message)
    segments = [(s1.sol, u0, u_net, 2)]
    if nuclei:
        # Stage 2: the network from T_net, starting from statistical equilibrium (as PRyMordial).
        Y_start = nse(u_net, *s1.y[:, -1])
        s2 = integrate.solve_ivp(lambda u, y: terms(u, y)[0], (u_net, u_end), Y_start,
                                 method="BDF", jac=lambda u, y: terms(u, y)[1], rtol=1e-8,
                                 atol=1e-16, dense_output=True)
        if not s2.success:
            raise RuntimeError(s2.message)
        segments.append((s2.sol, u_net, u_end, 6))
    u = np.linspace(u0, u_end, 1500)
    Ys = np.zeros((6, u.size))
    for f, a_, b_, nvar in segments:
        sel = (u >= a_) & (u <= b_)
        Ys[:nvar, sel] = f(u[sel])
    t = np.exp(u)
    Tg = bg.Tg(t)
    Y = {s: Ys[i] for i, s in enumerate(SPECIES)}
    Xfree = Y["n"] / (Y["n"] + Y["p"])
    # freeze-out: the first time the total weak rate drops below H (late on, free decay at
    # 1/tau_n overtakes the falling H again, which is not a freeze-out)
    tot = np.array([sum(lam(tt)) for tt in t])
    Hh = bg.hubble(t)
    k = np.flatnonzero(tot < Hh)
    if len(k):
        j = k[0]
        r0, r1 = math.log(tot[j - 1] / Hh[j - 1]), math.log(tot[j] / Hh[j])
        t_fr = float(math.exp(u[j - 1] + (u[j] - u[j - 1]) * r0 / (r0 - r1)))
    else:
        t_fr = float("nan")
    out_Yp = 4 * float(Ys[idx["a"], -1])
    i1 = np.searchsorted(t, 1.0)
    res = BBNResult(t=t, T_gamma=Tg, Y=Y, X_n_free=Xfree, Y_p=out_Yp,
                    D_H=float(Ys[idx["d"], -1] / Ys[idx["p"], -1]),
                    # tritium (half-life 12 yr) ends up as helium-3; PRyMordial reports the sum too
                    He3_H=float((Ys[idx["h"], -1] + Ys[idx["t"], -1]) / Ys[idx["p"], -1]),
                    np_ratio_1s=float(Y["n"][i1] / Y["p"][i1]), T_freeze=float(bg.Tg(t_fr)),
                    t_freeze=t_fr, t_nuc=float("nan"), T_nuc=float("nan"), X_n_at_nuc=float("nan"),
                    meta={"eta": eta, "tau_n": tau_n, "background": bg.name, "N_eff": bg.N_eff,
                          "nuclei": nuclei})
    if nuclei:
        ya = Ys[idx["a"]]
        half = 0.5 * ya[-1]
        j = np.searchsorted(ya, half)
        res.t_nuc = float(np.exp(np.interp(half, ya[j - 1:j + 1], u[j - 1:j + 1])))
        res.T_nuc = float(bg.Tg(res.t_nuc))
        # all neutrons (free + bound) per baryon at that time, i.e. what helium can still take
        res.X_n_at_nuc = float(np.interp(res.t_nuc, t, Y["n"] + Y["d"] + 2 * Y["t"] + Y["h"] + 2 * Y["a"]))
    return res


def free_neutrons_only(tau_n: float = C.TAU_N, background: Background | None = None) -> BBNResult:
    """The n/p ratio with no nuclei at all: freeze-out, then pure decay."""
    return run(tau_n=tau_n, background=background, nuclei=False)


# --- PRyMordial as the independent check ------------------------------------------------------

PRYM_CACHE = C.DERIVED_DIR / "prymordial_runs.json"


def prymordial(**cfg) -> dict | None:
    """Run PRyMordial (our pinned copy in data/external) once; results cached by configuration.

    Returns None if the code is not available (fetch_data.py downloads it)."""
    import json
    import subprocess
    import sys

    from fetch_data import PRYM_DIR

    key = json.dumps(cfg, sort_keys=True)
    cache = json.loads(PRYM_CACHE.read_text()) if PRYM_CACHE.exists() else {}
    if key in cache:
        return cache[key]
    if not PRYM_DIR.exists():
        return None
    cmd = [sys.executable, "-I", str(C.HERE / "prym_runner.py"), str(PRYM_DIR), json.dumps(cfg)]
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    lines = [ln for ln in res.stdout.splitlines() if ln.startswith("{")]
    if res.returncode != 0 or not lines:
        raise RuntimeError(f"PRyMordial failed: {res.stderr[-1500:]}")
    out = json.loads(lines[-1])
    cache[key] = out
    PRYM_CACHE.parent.mkdir(parents=True, exist_ok=True)
    PRYM_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    return out


def prymordial_background(eta: float) -> Background | None:
    """PRyMordial's own t, T_gamma, T_nu table (written by a run with save_bckg)."""
    from fetch_data import PRYM_DIR

    f = PRYM_DIR / "PRyMrates" / "thermo" / "Tgamma_Tnu.txt"
    if not f.exists():
        if prymordial(eta=eta, born=True, small=True, aTid=False, save_bckg=True) is None:
            return None
    d = np.loadtxt(f)
    t, Tg, Tn = d[:, 0], d[:, 1] * 1e-3, d[:, 2] * 1e-3
    return Background("prymordial", t, Tg, Tn, -np.gradient(np.log(Tn), t), 3.044)
