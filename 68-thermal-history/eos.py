"""The equation of state of the hot universe: g*rho(T), g*s(T), w(T) and c_s^2(T).

Two builds, both from the same PDG particle list:

* ``"ideal"``: every particle a free (ideal) gas with its PDG mass. Quarks and gluons are free
  down to the lattice crossover temperature, 156.5 MeV, where they are swapped abruptly for a free
  gas of every PDG hadron. Neutrinos decouple before any e+- annihilation (the textbook
  (4/11)^(1/3)). This is the counting that gives the famous 106.75, 10.75 and 3.36.
* ``"lattice"``: the same, except that quarks, gluons and hadrons are replaced by the strong sector
  of ``qcd.strong_sector`` (hadron gas -> lattice QCD -> perturbative QCD, smoothly joined), the
  e+- gamma plasma gets its finite-temperature QED correction, and neutrinos follow our decoupling
  run (``neutrinos.solve``), which gives N_eff = 3.044.

T is always the photon temperature. Electroweak interactions (and the electroweak crossover) are
not modelled: W, Z, H and t keep their zero-temperature masses.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.interpolate import CubicSpline

import config as C
import neutrinos
import particles
import qcd
import thermo

PI2 = np.pi ** 2


def grid(t_min: float = C.T_MIN, t_max: float = C.T_MAX, per_decade: int = C.N_PER_DECADE) -> np.ndarray:
    n = int(round(np.log10(t_max / t_min) * per_decade)) + 1
    return np.logspace(np.log10(t_min), np.log10(t_max), n)


@lru_cache(maxsize=None)
def decoupling_run(rates: str = "fd_me", separate: bool = True, qed_order: int = 3) -> dict:
    return neutrinos.solve(rates=rates, separate=separate, qed_order=qed_order)


def nu_temperatures(T, model: str, run: dict | None = None) -> tuple[np.ndarray, np.ndarray]:
    """T_nu_e/T and T_nu_mu/T (= T_nu_tau/T) at photon temperature T."""
    T = np.asarray(T, dtype=float)
    if model.startswith("ideal"):
        r = neutrinos.instantaneous_ratio(T)
        return r, r
    run = decoupling_run() if run is None else run
    Tg = run["T_gamma"][::-1]
    out = []
    for key in ("T_nu_e", "T_nu_mu"):
        ratio = (run[key] / run["T_gamma"])[::-1]
        r = np.interp(np.log(T), np.log(Tg), ratio, left=ratio[0], right=1.0)
        out.append(np.where(T >= run["T_gamma"][0], 1.0, r))
    return out[0], out[1]


def _nonstrong_species(sm):
    return [s for s in sm if s.group not in ("quark", "gluon", "neutrino")]


def build(T=None, model: str = "lattice", quark_masses: str = "pdg", run: dict | None = None,
          qcd_kwargs: dict | None = None, hadron_kwargs: dict | None = None) -> dict:
    """g*rho, g*s, g*p, w, c_s^2 and a per-species breakdown of rho, on the temperature grid T.

    model: "lattice" (default), "ideal" (free quarks/gluons -> PDG hadron gas at 156.5 MeV) or
    "ideal_pions" (free quarks/gluons -> pions only).
    """
    T = grid() if T is None else np.asarray(T, dtype=float)
    sm = particles.standard_model(quark_masses)
    hs = particles.hadrons(**(hadron_kwargs or {}))
    norm_r, norm_s, norm_p = PI2 / 30 * T ** 4, 2 * PI2 / 45 * T ** 3, PI2 / 90 * T ** 4
    parts_rho, parts_s, parts_p = {}, {}, {}

    # photons, charged leptons, W, Z, Higgs: free gases at T
    for sp in _nonstrong_species(sm):
        d = thermo.gas(sp, T)
        parts_rho[sp.name], parts_s[sp.name], parts_p[sp.name] = d["rho"], d["s"], d["p"]

    # neutrinos at their own temperature
    re_, rm_ = nu_temperatures(T, model, run)
    for sp in (s for s in sm if s.group == "neutrino"):
        Tn = T * (re_ if sp.name == "nu_e" else rm_)
        d = thermo.gas(sp, T, Tn)
        parts_rho[sp.name] = d["rho"]
        parts_s[sp.name] = d["s"]  # entropy of a decoupled species: s = (rho + p)/T_nu
        parts_p[sp.name] = d["p"]

    if model.startswith("ideal"):
        above = T >= C.T_SWITCH_IDEAL
        for sp in (s for s in sm if s.group in ("quark", "gluon")):
            d = thermo.gas(sp, T)
            parts_rho[sp.name] = np.where(above, d["rho"], 0.0)
            parts_s[sp.name] = np.where(above, d["s"], 0.0)
            parts_p[sp.name] = np.where(above, d["p"], 0.0)
        had = hs if model == "ideal" else [h for h in hs if h.pdgid in (111, 211, -211)]
        hr, hsum, hp = np.zeros_like(T), np.zeros_like(T), np.zeros_like(T)
        for h in had:
            d = thermo.gas(h, T)
            hr += d["rho"]; hsum += d["s"]; hp += d["p"]
        parts_rho["hadrons"] = np.where(above, 0.0, hr)
        parts_s["hadrons"] = np.where(above, 0.0, hsum)
        parts_p["hadrons"] = np.where(above, 0.0, hp)
    elif model == "lattice":
        strong = qcd.strong_sector(T, hs=hs, **(qcd_kwargs or {}))
        parts_rho["strong"] = strong["rho"] * T ** 4
        parts_s["strong"] = strong["s"] * T ** 3
        parts_p["strong"] = strong["p"] * T ** 4
        q = neutrinos.qed(T)
        parts_rho["qed"], parts_s["qed"], parts_p["qed"] = q["rho"], q["s"], q["P"]
    else:
        raise ValueError(model)

    rho = sum(parts_rho.values())
    s = sum(parts_s.values())
    p = sum(parts_p.values())
    lnT = np.log(T)
    # c_s^2 = dp/drho along the evolution, from smooth splines in ln T
    if model == "lattice":
        cs2 = CubicSpline(lnT, p)(lnT, 1) / CubicSpline(lnT, rho)(lnT, 1)
    else:
        cs2 = np.full_like(T, np.nan)  # the abrupt switch makes it undefined there
    return {
        "T": T, "model": model,
        "g_rho": rho / norm_r, "g_s": s / norm_s, "g_p": p / norm_p,
        "w": p / rho, "cs2": cs2,
        "rho": rho, "s": s, "p": p,
        "T_nu_e_over_T": re_, "T_nu_mu_over_T": rm_,
        "parts_g_rho": {k: v / norm_r for k, v in parts_rho.items()},
        "parts_g_s": {k: v / norm_s for k, v in parts_s.items()},
    }


# --- Published curves, on our grid ----------------------------------------------------------


@lru_cache(maxsize=None)
def saikawa_shirai() -> dict:
    a = np.loadtxt(C.RAW_DIR / "ss2018_standardmodel.dat")
    w = np.loadtxt(C.RAW_DIR / "ss2018_eos_w.dat")
    return {"T": a[:, 0], "g_rho": a[:, 1], "g_rho_err": a[:, 2], "g_s": a[:, 3], "g_s_err": a[:, 4],
            "w_T": w[:, 0], "w": w[:, 1], "w_err": w[:, 2]}


@lru_cache(maxsize=None)
def laine_meyer_schroeder() -> dict:
    a = np.loadtxt(C.RAW_DIR / "lms2015_standardmodel.dat")
    return {"T": a[:, 0] / 1000.0, "p": a[:, 1], "rho": a[:, 2], "s": a[:, 3], "c": a[:, 4],
            "w": a[:, 5], "cs2": a[:, 6], "g_rho": a[:, 7], "g_s": a[:, 8], "i_eff": a[:, 9]}


@lru_cache(maxsize=None)
def escudero2020() -> dict:
    a = np.loadtxt(C.RAW_DIR / "escudero2020_sm_thermo.dat")
    return {"t": a[:, 0], "T": a[:, 1] / 1000.0, "Tg_over_Tnu": a[:, 2], "g_rho": a[:, 3],
            "g_s": a[:, 4], "a": a[:, 5], "z_gamma": a[:, 6], "z_nu": a[:, 7]}


def on_grid(src: dict, key: str, T) -> np.ndarray:
    """A published curve interpolated (linearly in log T) onto T; NaN outside its range."""
    T = np.asarray(T, dtype=float)
    x, y = np.log(src["T"]), src[key]
    o = np.argsort(x)
    x, y = x[o], y[o]
    out = np.interp(np.log(T), x, y)
    return np.where((T >= np.exp(x[0])) & (T <= np.exp(x[-1])), out, np.nan)
