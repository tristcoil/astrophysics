"""The hot plasma that sets the expansion rate: g*rho(T), g*s(T) and what follows from them.

In the radiation era the energy density and the entropy density are

    rho(T) = (pi^2/30) g*rho(T) T^4,        s(T) = (2 pi^2/45) g*s(T) T^3,

and the expansion rate is H = sqrt(rho/3) / M_P (M_P the reduced Planck mass). A particle species
counts fully in g* while it is relativistic (2 for the photon, 7/8 per fermion state) and drops
out when the temperature falls below its mass, handing its entropy to what is left.

Freeze-out needs one more combination (Gondolo & Gelmini 1991; Saikawa & Shirai 2020, eq. 3.6):

    g*^(1/2)(T) = (g*s / sqrt(g*rho)) (1 + (1/3) dln g*s / dln T),

because the Boltzmann equation is written per unit entropy and the entropy changes as species
annihilate. Everything else in this module is bookkeeping for the tables that give g*(T):

* ``ss2020``  Saikawa & Shirai (2020): the table behind their relic cross sections (default)
* ``ss2020_hi`` / ``ss2020_lo``: their upper and lower uncertainty curves
* ``ss2018``  Saikawa & Shirai (2018): the original 10,000-point table
* ``lms2015`` Laine & Meyer (2015, Dec 2019 update); its heat capacity gives the
              derivative analytically, as the authors recommend
* ``lms2015_gs_eq_grho``: the same with g*s set equal to g*rho, the approximation Steigman,
              Dasgupta & Beacom (2012) made
* ``ideal``  a textbook ideal gas of Standard Model particles: free quarks and gluons above the
              lattice crossover temperature, pions below, neutrinos decoupled at 1 MeV
* ``constant_<g>``: g* fixed at one value (the naive estimate)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy import integrate
from scipy.interpolate import CubicSpline, PchipInterpolator

import config as C

EOS_DIR = C.RAW_DIR / "eos"
#: P2's derived curves, used when that project has been run (otherwise those variants are skipped)
P2_CURVES = C.ROOT / "68-thermal-history" / "data" / "derived" / "curves.npz"


def available(name: str) -> bool:
    return P2_CURVES.exists() if name.startswith("p2") else True


@dataclass
class EoS:
    """g*rho, g*s and g*^(1/2) as smooth functions of temperature (GeV).

    Interpolation is a cubic spline in ln T of ln g; outside the table g is held at the edge
    values (below 10 keV nothing changes any more; above the top of the table the Standard
    Model has nothing left to switch on).
    """

    name: str
    label: str
    lnT: np.ndarray
    ln_grho: np.ndarray
    ln_gs: np.ndarray
    #: optional exact 1 + (1/3) dln gs/dln T on the same grid (LMS give it via the heat capacity)
    entropy_factor: np.ndarray | None = None
    #: monotone (PCHIP) interpolation for tables with steps, where a cubic spline would overshoot
    monotone: bool = False

    def __post_init__(self):
        spline = PchipInterpolator if self.monotone else CubicSpline
        self._grho = spline(self.lnT, self.ln_grho)
        self._gs = spline(self.lnT, self.ln_gs)
        self._dgs = self._gs.derivative()
        self._ef = CubicSpline(self.lnT, self.entropy_factor) if self.entropy_factor is not None else None
        self.T_min, self.T_max = float(np.exp(self.lnT[0])), float(np.exp(self.lnT[-1]))

    def _clip(self, T):
        return np.clip(np.log(np.asarray(T, dtype=float)), self.lnT[0], self.lnT[-1])

    def g_rho(self, T):
        return np.exp(self._grho(self._clip(T)))

    def g_s(self, T):
        return np.exp(self._gs(self._clip(T)))

    def dlngs_dlnT(self, T):
        lt = np.log(np.asarray(T, dtype=float))
        inside = (lt >= self.lnT[0]) & (lt <= self.lnT[-1])
        return np.where(inside, self._dgs(self._clip(T)), 0.0)

    def entropy_derivative_factor(self, T):
        """1 + (1/3) dln g*s / dln T."""
        if self._ef is not None:
            lt = np.log(np.asarray(T, dtype=float))
            inside = (lt >= self.lnT[0]) & (lt <= self.lnT[-1])
            return np.where(inside, self._ef(self._clip(T)), 1.0)
        return 1.0 + self.dlngs_dlnT(T) / 3.0

    def gstar_half(self, T):
        """g*^(1/2) = (g*s/sqrt(g*rho)) (1 + (1/3) dln g*s/dln T)."""
        return self.g_s(T) / np.sqrt(self.g_rho(T)) * self.entropy_derivative_factor(T)

    @property
    def g_s0(self) -> float:
        """g*s today (photons + the colder neutrinos), the table's low-temperature end."""
        return float(np.exp(self.ln_gs[0]))

    # --- derived quantities -----------------------------------------------------------

    def entropy_density(self, T):
        return 2 * math.pi ** 2 / 45 * self.g_s(T) * np.asarray(T) ** 3

    def hubble(self, T):
        """H(T) in GeV for a radiation-dominated universe."""
        return np.sqrt(math.pi ** 2 * self.g_rho(T) / 90.0) * np.asarray(T) ** 2 / C.M_PL_RED

    def time_of(self, T_grid):
        """Cosmic time t(T) in seconds by integrating dt = -dT / (T H) * dln s/dln T / 3.

        From s a^3 = const: dt = -(1/(3H)) dln s, and dln s = (3 + dln gs/dln T) dln T.
        Starts at t = 1/(2H) at the highest temperature of the grid (pure radiation there).
        """
        T_grid = np.sort(np.asarray(T_grid, dtype=float))[::-1]
        lnT = np.log(T_grid)
        integrand = (3.0 + self.dlngs_dlnT(T_grid)) / (3.0 * self.hubble(T_grid))
        t = np.concatenate([[0.0], integrate.cumulative_trapezoid(integrand, -lnT)])
        t += 1.0 / (2.0 * self.hubble(T_grid[0]))
        return T_grid, t * C.HBAR_GEV_S


# --- loaders ------------------------------------------------------------------------------


def _from_columns(name, label, T, grho, gs, ef=None, monotone=False) -> EoS:
    order = np.argsort(T)
    T, grho, gs = T[order], grho[order], gs[order]
    keep = T > 0
    if ef is not None:
        ef = ef[order][keep]
    return EoS(name, label, np.log(T[keep]), np.log(grho[keep]), np.log(gs[keep]), ef, monotone)


def _ss2020(path, name, label) -> EoS:
    d = np.loadtxt(path)                     # T [GeV], g_s, g_rho
    return _from_columns(name, label, d[:, 0], d[:, 2], d[:, 1])


def _lms2015(gs_equal_grho: bool = False) -> EoS:
    d = np.loadtxt(EOS_DIR / "lms2015_standardmodel.dat")
    T = d[:, 0] * 1e-3                       # MeV -> GeV
    s_T3, c_T3 = d[:, 3], d[:, 4]
    g_eff, h_eff = d[:, 7], d[:, 8]
    # dln s/dln T = c/s exactly (c = de/dT = T ds/dT), so 1 + (1/3) dln h/dln T = c/(3 s).
    ef = c_T3 / (3.0 * s_T3)
    if gs_equal_grho:
        # SDB12: g = g*rho = g*s. The derivative factor then comes from g*rho itself.
        e = _from_columns("lms2015_gs_eq_grho", "Laine et al. 2015, g*s = g*rho (as SDB12)", T, g_eff, g_eff)
        return e
    return _from_columns("lms2015", "Laine & Meyer 2015", T, g_eff, h_eff, ef)


@lru_cache(maxsize=None)
def load(name: str = "ss2020") -> EoS:
    if name == "ss2020":
        return _ss2020(EOS_DIR / "ss2020_eos2020.dat", name, "Saikawa & Shirai 2020")
    if name == "ss2020_hi":
        return _ss2020(EOS_DIR / "ss2020_eos2020_err1.dat", name, "Saikawa & Shirai 2020 (upper)")
    if name == "ss2020_lo":
        return _ss2020(EOS_DIR / "ss2020_eos2020_err2.dat", name, "Saikawa & Shirai 2020 (lower)")
    if name == "ss2018":
        d = np.loadtxt(EOS_DIR / "ss2018_standardmodel.dat")  # T, g_rho, err, g_s, err
        return _from_columns(name, "Saikawa & Shirai 2018", d[:, 0], d[:, 1], d[:, 3])
    if name == "lms2015":
        return _lms2015(False)
    if name == "lms2015_gs_eq_grho":
        return _lms2015(True)
    if name in ("p2", "p2_ideal"):
        # P2 (68-thermal-history): its own interacting build (hadron gas -> lattice QCD ->
        # perturbative QCD, with neutrino decoupling) and its ideal gas of PDG particles
        d = np.load(P2_CURVES)
        if name == "p2":
            return _from_columns(name, "P2's own g*(T) (68-thermal-history)", d["T"], d["g_rho"], d["g_s"],
                                 monotone=True)
        # The ideal gas switches abruptly at the crossover. Its entropy release is a near-delta spike in
        # g*^(1/2) that no quadrature on a log-T grid resolves, so the step is smoothed over 2% in T
        # (as in our own ideal gas; the integrated release is unchanged), with monotone interpolation.
        T = d["T"]
        out = []
        for col in ("g_rho_ideal", "g_s_ideal"):
            arr = np.array(d[col], dtype=float)
            near = np.flatnonzero(np.abs(np.log(T / T_QCD)) < 0.2)
            win = np.flatnonzero(np.abs(np.log(T / T_QCD)) < 0.4)
            k = np.exp(-0.5 * (np.log(T[near][:, None] / T[win][None, :]) / 0.02) ** 2)
            arr[near] = (k @ arr[win]) / k.sum(axis=1)
            out.append(arr)
        return _from_columns(name, "P2's ideal gas (68-thermal-history), QCD step smoothed over 2%", T, out[0], out[1],
                             monotone=True)
    if name == "ideal":
        return ideal_gas()
    if name.startswith("constant_"):
        g = float(name.split("_", 1)[1])
        T = np.array([1e-12, 1e12])
        return _from_columns(name, f"constant g* = {g:g}", T, np.full(2, g), np.full(2, g))
    raise KeyError(name)


# --- a textbook ideal gas of Standard Model particles ------------------------------------

#: (name, mass in GeV, internal states, 'B' or 'F', phase): 'all' = always present,
#: 'qgp' = only above the QCD crossover, 'hadron' = only below it. PDG 2024 masses
#: (MS-bar quark masses for the light quarks, pole masses for the top).
SM_SPECIES = [
    ("photon", 0.0, 2, "B", "all"),
    ("gluons", 0.0, 16, "B", "qgp"),
    ("W", 80.369, 6, "B", "all"),
    ("Z", 91.188, 3, "B", "all"),
    ("Higgs", 125.20, 1, "B", "all"),
    ("top", 172.57, 12, "F", "all"),
    ("bottom", 4.183, 12, "F", "all"),
    ("charm", 1.273, 12, "F", "all"),
    ("strange", 0.0935, 12, "F", "qgp"),
    ("down", 0.00470, 12, "F", "qgp"),
    ("up", 0.00216, 12, "F", "qgp"),
    ("tau", 1.77693, 4, "F", "all"),
    ("muon", 0.1056584, 4, "F", "all"),
    ("electron", C.M_E, 4, "F", "all"),
    ("pi+-", 0.13957, 2, "B", "hadron"),
    ("pi0", 0.134977, 1, "B", "hadron"),
]
#: Quarks and gluons switch to pions at HotQCD's crossover temperature (Bazavov et al. 2019).
T_QCD = 0.1565
#: Neutrinos decouple instantaneously at this temperature (then e+e- annihilation heats photons).
T_NU_DEC = 1.0e-3


@lru_cache(maxsize=None)
def _gas_tables():
    """I_rho(z) and I_s(z) = I_rho + I_p/3 for bosons and fermions on a grid of z = m/T."""
    z = np.concatenate([[0.0], np.logspace(-4, np.log10(60.0), 600)])
    out = {}
    for stat, sign in (("B", -1.0), ("F", 1.0)):
        irho, ip = [], []
        for zz in z:
            f = (lambda u, zz=zz: u * u * math.sqrt(max(u * u - zz * zz, 0.0)) / (math.exp(u) + sign))
            fp = (lambda u, zz=zz: max(u * u - zz * zz, 0.0) ** 1.5 / (math.exp(u) + sign))
            a = integrate.quad(f, zz, zz + 80.0, limit=200, epsabs=0, epsrel=1e-10)[0]
            b = integrate.quad(fp, zz, zz + 80.0, limit=200, epsabs=0, epsrel=1e-10)[0]
            irho.append(a)
            ip.append(b)
        out[stat] = (z, np.array(irho), np.array(ip))
    return out


def _species_g(m, g, stat, T):
    """Contributions (g_rho, g_s) of one species at temperatures T (vectorised)."""
    z_grid, irho, ip = _gas_tables()[stat]
    z = m / np.asarray(T)
    Ir = np.interp(z, z_grid, irho, right=0.0)
    Is = np.interp(z, z_grid, irho + ip / 3.0, right=0.0)
    return g * 15 / math.pi ** 4 * Ir, g * 45 / (4 * math.pi ** 4) * Is


def ideal_gas() -> EoS:
    T = np.logspace(-5.0, 5.0, 4001)
    grho = np.zeros_like(T)
    gs = np.zeros_like(T)
    for name, m, g, stat, phase in SM_SPECIES:
        r, s = _species_g(m, g, stat, T)
        if phase == "qgp":
            r, s = np.where(T >= T_QCD, r, 0.0), np.where(T >= T_QCD, s, 0.0)
        elif phase == "hadron":
            r, s = np.where(T < T_QCD, r, 0.0), np.where(T < T_QCD, s, 0.0)
        grho += r
        gs += s
    # Three neutrino species, 2 states each, at T_nu: equal to T before decoupling; after it,
    # the e+e-photon plasma keeps its own entropy, so (T_nu/T)^3 = g_s,e-gamma(T)/g_s,e-gamma(T_dec).
    r_e, s_e = _species_g(C.M_E, 4, "F", T)
    r_e_dec, s_e_dec = _species_g(C.M_E, 4, "F", np.array([T_NU_DEC]))
    ratio3 = np.where(T >= T_NU_DEC, 1.0, (2.0 + s_e) / (2.0 + s_e_dec[0]))
    tnu = ratio3 ** (1.0 / 3.0)
    grho += 6 * 7 / 8 * tnu ** 4
    gs += 6 * 7 / 8 * tnu ** 3
    # The step at T_QCD is smoothed over 2% in temperature so the spline stays well behaved.
    near = np.flatnonzero(np.abs(np.log(T / T_QCD)) < 0.2)
    win = np.flatnonzero(np.abs(np.log(T / T_QCD)) < 0.4)
    k = np.exp(-0.5 * (np.log(T[near][:, None] / T[win][None, :]) / 0.02) ** 2)
    for arr in (grho, gs):
        arr[near] = (k @ arr[win]) / k.sum(axis=1)
    return _from_columns("ideal", "textbook ideal gas (sharp QCD switch at 156.5 MeV)", T, grho, gs)


def describe(eos: EoS, temps=(1e3, 100.0, 10.0, 1.0, 0.2, 0.1, 1e-2, 1e-3, 1e-4)) -> list:
    return [{"T_GeV": T, "g_rho": float(eos.g_rho(T)), "g_s": float(eos.g_s(T)),
             "gstar_half": float(eos.gstar_half(T))} for T in temps]
