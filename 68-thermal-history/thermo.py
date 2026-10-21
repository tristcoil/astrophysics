"""Ideal relativistic gases: the energy density, pressure, entropy and number density of one species.

For a species of mass m with g internal states, in equilibrium at temperature T with no chemical
potential (natural units),

    rho = g/(2 pi^2) Int_m^inf  E^2 sqrt(E^2 - m^2) / (exp(E/T) -+ 1) dE
    p   = g/(6 pi^2) Int_m^inf  (E^2 - m^2)^(3/2)  / (exp(E/T) -+ 1) dE
    n   = g/(2 pi^2) Int_m^inf  E   sqrt(E^2 - m^2) / (exp(E/T) -+ 1) dE
    s   = (rho + p) / T

with - for bosons and + for fermions. In x = m/T each is T^4 (T^3 for n and s) times a function
of x alone. We integrate those functions once, by adaptive quadrature to a relative accuracy of
~1e-12, on 2,400 values of x between 1e-5 and 200, and interpolate with splines. Below x = 1e-5
the massless values are used (relative error < 1e-10); above x = 200 (Boltzmann factor e^-200) a
species contributes nothing.

The g* bookkeeping: a species at temperature T_i adds

    g_rho = rho_i / (pi^2 T^4 / 30)      g_s = s_i / (2 pi^2 T^3 / 45)

to the totals, where T is the photon temperature; a massless boson adds g (1 for each state), a
massless fermion 7/8 g, each scaled by (T_i/T)^4 or (T_i/T)^3.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

PI2 = np.pi ** 2
ZETA3 = 1.2020569031595942
X_MIN, X_MAX, N_X = 1e-5, 200.0, 2400

#: Massless limits of rho/(g T^4), p/(g T^4), n/(g T^3).
MASSLESS = {
    +1: (PI2 / 30, PI2 / 90, ZETA3 / PI2),
    -1: (7 / 8 * PI2 / 30, 7 / 8 * PI2 / 90, 3 / 4 * ZETA3 / PI2),
}


def _integrals(x: float, stat: int) -> tuple[float, float, float]:
    """rho/(gT^4), p/(gT^4), n/(gT^3) times e^x, with eps = x + u (u from 0 to infinity)."""

    def occ(u):  # e^x / (exp(x+u) -+ 1) = e^-u / (1 -+ e^-(x+u))
        return np.exp(-u) / (1.0 - stat * np.exp(-(x + u)))

    def k(u):  # momentum p/T = sqrt(eps^2 - x^2)
        return np.sqrt(u * (2 * x + u))

    opts = dict(epsabs=0.0, epsrel=1e-12, limit=400)
    r = quad(lambda u: (x + u) ** 2 * k(u) * occ(u), 0, np.inf, **opts)[0] / (2 * PI2)
    p = quad(lambda u: k(u) ** 3 * occ(u), 0, np.inf, **opts)[0] / (6 * PI2)
    n = quad(lambda u: (x + u) * k(u) * occ(u), 0, np.inf, **opts)[0] / (2 * PI2)
    return r, p, n


@lru_cache(maxsize=None)
def _tables(stat: int):
    lx = np.linspace(np.log(X_MIN), np.log(X_MAX), N_X)
    x = np.exp(lx)
    vals = np.array([_integrals(xi, stat) for xi in x])  # scaled by e^x
    # Spline log(f e^x) in log x: smooth and positive everywhere.
    return {name: CubicSpline(lx, np.log(vals[:, i])) for i, name in enumerate(("rho", "p", "n"))}


def functions(x, stat: int) -> dict:
    """rho/(gT^4), p/(gT^4), n/(gT^3), s/(gT^3) and x-derivatives d/dx of rho and p, for x = m/T."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    tab = _tables(stat)
    r0, p0, n0 = MASSLESS[stat]
    out = {k: np.zeros_like(x) for k in ("rho", "p", "n", "drho", "dp")}
    small = x < X_MIN
    out["rho"][small], out["p"][small], out["n"][small] = r0, p0, n0
    mid = (x >= X_MIN) & (x <= X_MAX)
    if mid.any():
        lx = np.log(x[mid])
        for name in ("rho", "p", "n"):
            f = np.exp(tab[name](lx) - x[mid])
            out[name][mid] = f
            if name in ("rho", "p"):
                # d f / dx = f * (d log(f e^x)/d log x / x - 1)
                out["d" + name][mid] = f * (tab[name](lx, 1) / x[mid] - 1.0)
    out["s"] = out["rho"] + out["p"]
    return out


@dataclass(frozen=True)
class Species:
    """One particle species (all its spin, colour, charge and antiparticle states together)."""

    name: str
    label: str
    group: str
    mass: float          # GeV
    g: float             # internal states
    stat: int            # +1 boson, -1 fermion
    pdgid: int = 0

    @property
    def g_eff_massless(self) -> float:
        return self.g * (1.0 if self.stat > 0 else 7 / 8)


def gas(sp: Species, T, T_sp=None) -> dict:
    """Thermodynamics of one species in GeV units: rho, p, s, n and d(rho)/dT_sp.

    ``T_sp`` is the species' own temperature (default: T).
    """
    T = np.atleast_1d(np.asarray(T, dtype=float))
    Ts = T if T_sp is None else np.atleast_1d(np.asarray(T_sp, dtype=float))
    x = sp.mass / Ts
    f = functions(x, sp.stat)
    g = sp.g
    return {
        "rho": g * Ts ** 4 * f["rho"],
        "p": g * Ts ** 4 * f["p"],
        "s": g * Ts ** 3 * f["s"],
        "n": g * Ts ** 3 * f["n"],
        # d rho/dT = g T^3 (4 f - x f'(x))
        "drho_dT": g * Ts ** 3 * (4 * f["rho"] - x * f["drho"]),
    }


def gstar(sp: Species, T, T_sp=None) -> dict:
    """The species' contributions to g*rho, g*s and g*p (normalised with the photon temperature T)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    d = gas(sp, T, T_sp)
    return {
        "rho": d["rho"] / (PI2 / 30 * T ** 4),
        "s": d["s"] / (2 * PI2 / 45 * T ** 3),
        "p": d["p"] / (PI2 / 90 * T ** 4),
    }


def bessel_series(x, stat: int, n_terms: int = 400) -> tuple[np.ndarray, np.ndarray]:
    """Independent check: rho/(gT^4) and p/(gT^4) from the Bessel-function series.

    p/(gT^4) = 1/(2 pi^2) sum_k (+-1)^(k+1) x^2/k^2 K2(k x);
    rho/(gT^4) = 1/(2 pi^2) sum_k (+-1)^(k+1) [3 x^2/k^2 K2(k x) + x^3/k K1(k x)].
    """
    from scipy.special import kn

    x = np.atleast_1d(np.asarray(x, dtype=float))
    k = np.arange(1, n_terms + 1)[:, None]
    sign = 1.0 if stat > 0 else (-1.0) ** (k + 1)
    kx = k * x[None, :]
    p = (sign * x ** 2 / k ** 2 * kn(2, kx)).sum(0) / (2 * PI2)
    r = (sign * (3 * x ** 2 / k ** 2 * kn(2, kx) + x ** 3 / k * kn(1, kx))).sum(0) / (2 * PI2)
    return r, p
