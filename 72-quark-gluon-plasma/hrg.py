"""The hadron resonance gas: every hadron as an ideal quantum gas.

For one species of mass m, degeneracy g and chemical potential mu at
temperature T, the grand-canonical density, energy density and pressure are
momentum integrals over the Bose-Einstein or Fermi-Dirac occupancy
f = 1 / (exp((E - mu)/T) -+ 1), E = sqrt(p^2 + m^2):

    n   = g / (2 pi^2) * int p^2 f dp
    eps = g / (2 pi^2) * int p^2 E f dp
    P   = g / (6 pi^2) * int p^4 / E f dp

(the Boltzmann limit would give the familiar n = g m^2 T K2(m/T) / (2 pi^2)).
The integrals are done with 64-point Gauss-Laguerre quadrature in p/T,
which agrees with the exact Bessel-function series to 1e-9 for every hadron
from the pion to 4He (test_physics.py checks this) at a fiftieth of the cost.

Each hadron's chemical potential is mu = B mu_B + Q mu_Q + S mu_S. For a given
(T, mu_B) the strangeness and charge potentials are solved from the initial
state's conservation laws: zero net strangeness, and Q/B of the colliding
nuclei. An optional strangeness under-saturation factor gamma_S multiplies
each hadron's fugacity by gamma_S^(number of strange + antistrange quarks).

Resonances with a width are smeared over a relativistic Breit-Wigner mass
distribution between max(threshold, m - 2 Gamma) and m + 2 Gamma (10-point
Gauss-Legendre, as Thermal-FIST's "BWTwoGamma"), or kept at their pole mass
("zero width").

Units: GeV for T, masses and chemical potentials; densities in fm^-3,
pressure and energy density in GeV/fm^3.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.polynomial.laguerre import laggauss

import config as cfg
import hadrons

HC3 = cfg.HBARC ** 3
TWO_PI2 = 2.0 * np.pi ** 2
N_QUAD = 64


@dataclass
class State:
    """Primordial (before decays) thermodynamics of every species."""
    T: float
    muB: float
    muQ: float
    muS: float
    gammaS: float
    n: np.ndarray        # fm^-3, per species
    p: np.ndarray        # GeV/fm^3
    e: np.ndarray        # GeV/fm^3
    dn_dmu: np.ndarray   # fm^-3 GeV^-1


class HRG:
    """A hadron resonance gas built from one of the hadron lists."""

    def __init__(self, list_name: str = cfg.DEFAULT_LIST, widths: str = cfg.DEFAULT_WIDTHS,
                 stats: str = "quantum", with_nuclei: bool = True, n_bw: int = 10,
                 exclude: tuple = ()):
        self.list_name, self.widths, self.stats = list_name, widths, stats
        self.hl = hl = hadrons.load(list_name, with_nuclei=with_nuclei)
        self.with_nuclei = with_nuclei
        n = len(hl)
        g = hl.g.astype(float).copy()
        for p in exclude:                 # drop states (systematic variants)
            if p in hl.index:
                g[hl.index[p]] = 0.0
        self.g = g
        # Mass points per species: the pole mass, or Breit-Wigner nodes
        xs, ws = np.polynomial.legendre.leggauss(n_bw)
        sp, ms, wt = [], [], []
        for i in range(n):
            m0, G, thr = hl.mass[i], hl.width[i], hl.threshold[i]
            if widths == "bw" and G > 1e-5 and not hl.is_nucleus[i]:
                a, b = max(thr, m0 - 2 * G), m0 + 2 * G
                if b > a:
                    m = 0.5 * (b - a) * xs + 0.5 * (b + a)
                    bw = m0 * G * m / ((m * m - m0 * m0) ** 2 + m0 * m0 * G * G)
                    w = 0.5 * (b - a) * ws * bw
                    sp.extend([i] * n_bw), ms.extend(m), wt.extend(w / w.sum())
                    continue
            sp.append(i), ms.append(m0), wt.append(1.0)
        self.pt_species = np.array(sp)
        self.pt_mass = np.array(ms)
        self.pt_weight = np.array(wt)
        stat = hl.stat[self.pt_species].astype(float)
        self.pt_stat = np.zeros_like(stat) if stats == "boltzmann" else stat
        self.pt_pref = g[self.pt_species] * self.pt_weight / TWO_PI2
        self.charges = np.vstack([hl.B, hl.Q, hl.S]).T.astype(float)   # (n, 3)
        x, w = laggauss(N_QUAD)
        self._x, self._w = x, w * np.exp(x)       # int_0^inf g(x) dx = sum w_k g(x_k)
        self._cacheT = None
        self.n_species = n

    # -- kernels -------------------------------------------------------------

    def _kernels(self, T: float) -> tuple:
        """Momentum-space weights at temperature T (cached: the conservation-law
        solver and the fitter call many chemical potentials at one T)."""
        if self._cacheT is not None and self._cacheT[0] == T:
            return self._cacheT[1]
        p = T * self._x[None, :]
        m = self.pt_mass[:, None]
        E = np.sqrt(p * p + m * m)
        w = T * self._w[None, :] * self.pt_pref[:, None]
        k = (E, w * p * p, w * p * p * E, w * p ** 4 / (3 * E))
        self._cacheT = (T, k)
        return k

    def _state_mu(self, T: float, mu_i: np.ndarray) -> State:
        mu = mu_i[self.pt_species]
        if np.any((mu >= self.pt_mass) & (self.pt_stat < 0)):
            raise ValueError("boson chemical potential above its mass (Bose condensation)")
        E, wn, we, wp = self._kernels(T)
        ea = np.exp(-(E - mu[:, None]) / T)
        st = self.pt_stat[:, None]
        f = ea / (1.0 + st * ea)
        N = self.n_species
        agg = lambda v: np.bincount(self.pt_species, v, minlength=N) / HC3  # noqa: E731
        return State(T, np.nan, np.nan, np.nan, 1.0,
                     agg(np.sum(wn * f, axis=1)), agg(np.sum(wp * f, axis=1)),
                     agg(np.sum(we * f, axis=1)),
                     agg(np.sum(wn * f * (1.0 - st * f), axis=1) / T))

    # -- one temperature and set of chemical potentials ---------------------

    def species_mu(self, muB: float, muQ: float, muS: float, T: float,
                   gammaS: float = 1.0) -> np.ndarray:
        hl = self.hl
        mu = hl.B * muB + hl.Q * muQ + hl.S * muS
        if gammaS != 1.0:
            mu = mu + T * hl.absS * np.log(gammaS)
        return mu

    def state(self, T: float, muB: float = 0.0, muQ: float = 0.0, muS: float = 0.0,
              gammaS: float = 1.0) -> State:
        st = self._state_mu(T, self.species_mu(muB, muQ, muS, T, gammaS))
        st.muB, st.muQ, st.muS, st.gammaS = muB, muQ, muS, gammaS
        return st

    # -- conservation laws --------------------------------------------------

    def solve_mu(self, T: float, muB: float, q_over_b: float | None = cfg.Q_OVER_B["Pb"],
                 gammaS: float = 1.0, muS0: float = 0.0, muQ0: float = 0.0,
                 tol: float = 1e-13) -> State:
        """Solve mu_S (net strangeness zero) and mu_Q (Q = q_over_b * B).

        With ``q_over_b=None`` mu_Q is held at zero (STAR's choice). Newton
        iterations on the exact Jacobian (the charge susceptibilities).
        """
        muS, muQ = muS0, (muQ0 if q_over_b is not None else 0.0)
        Bv, Qv, Sv = self.hl.B, self.hl.Q, self.hl.S
        for _ in range(60):
            st = self.state(T, muB, muQ, muS, gammaS)
            nB, nQ, nS = st.n @ Bv, st.n @ Qv, st.n @ Sv
            d = st.dn_dmu
            if q_over_b is None:
                step = nS / (d @ (Sv * Sv))
                muS -= step
                if abs(step) < tol:
                    break
                continue
            f = np.array([nS, nQ - q_over_b * nB])
            J = np.array([[d @ (Sv * Sv), d @ (Sv * Qv)],
                          [d @ (Qv * Sv) - q_over_b * d @ (Bv * Sv),
                           d @ (Qv * Qv) - q_over_b * d @ (Bv * Qv)]])
            step = np.linalg.solve(J, f)
            muS, muQ = muS - step[0], muQ - step[1]
            if np.max(np.abs(step)) < tol:
                break
        else:
            raise RuntimeError(f"mu_S/mu_Q did not converge at T={T}, muB={muB}")
        return self.state(T, muB, muQ, muS, gammaS)

    # -- bulk thermodynamics -----------------------------------------------

    def totals(self, st: State) -> dict:
        """Pressure, energy and entropy densities of the whole gas."""
        n_charge = st.n @ self.charges            # nB, nQ, nS
        P, E = st.p.sum(), st.e.sum()
        mu = np.array([st.muB, st.muQ, st.muS])
        s = (E + P - mu @ n_charge) / st.T        # fm^-3
        return {"p": P, "e": E, "s": s, "n": st.n.sum(), "nB": n_charge[0],
                "nQ": n_charge[1], "nS": n_charge[2]}

    def susceptibilities(self, T: float) -> dict:
        """chi_XY = d^2(p/T^4)/d(mu_X/T)d(mu_Y/T) at mu = 0."""
        st = self.state(T)
        d = st.dn_dmu * HC3 / T ** 2              # dimensionless per species
        hl = self.hl
        out = {}
        for a, va in (("B", hl.B), ("Q", hl.Q), ("S", hl.S)):
            for b, vb in (("B", hl.B), ("Q", hl.Q), ("S", hl.S)):
                out[a + b] = float(d @ (va * vb))
        return out

    # -- excluded volume ----------------------------------------------------

    def excluded_volume(self, st: State, r0_fm: float = cfg.EXCLUDED_VOLUME_R0_FM) -> dict:
        """Common eigenvolume b = 16 pi r0^3 / 3 for every hadron (the diagonal
        excluded-volume model of Rischke et al. 1991).

        Solves P = sum_i P_i^ideal(T, mu_i - b P) and returns each species'
        density suppression n_EV / n_ideal at the same (T, mu). With one
        eigenvolume for all, the factor is the same for every species (up to
        quantum statistics): it changes the fitted volume, not T or mu_B.
        """
        b = 16 * np.pi * r0_fm ** 3 / 3            # fm^3
        P = st.p.sum()
        for _ in range(500):
            P_new = self._shifted(st, b * P).p.sum()   # b P in GeV
            if abs(P_new - P) < 1e-15:
                P = P_new
                break
            P = P_new
        s2 = self._shifted(st, b * P)
        n_ev = s2.n / (1 + b * s2.n.sum())
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(st.n > 0, n_ev / st.n, 0.0)
        return {"p": P, "b_fm3": b, "suppression": ratio, "n_total": n_ev.sum()}

    def _shifted(self, st: State, shift: float) -> State:
        """State with every species' chemical potential lowered by ``shift``."""
        mu_i = self.species_mu(st.muB, st.muQ, st.muS, st.T, st.gammaS) - shift
        return self._state_mu(st.T, mu_i)


@lru_cache(maxsize=None)
def model(list_name: str = cfg.DEFAULT_LIST, widths: str = cfg.DEFAULT_WIDTHS,
          stats: str = "quantum", with_nuclei: bool = True, exclude: tuple = ()) -> HRG:
    return HRG(list_name, widths, stats, with_nuclei, exclude=exclude)
