"""The Higgs potential, the masses it gives, and when it switched on.

**The potential.** The Higgs field is a doublet H with potential
V = -mu^2 |H|^2 + lambda |H|^4 (the "Mexican hat"). Its minimum is not at
H = 0 but on a circle of radius |H| = v / sqrt(2), with v = mu / sqrt(lambda).
The Fermi constant of muon decay fixes v = (sqrt(2) G_F)^(-1/2) = 246.22 GeV;
the Higgs mass, the curvature of the potential at the minimum, fixes
lambda = M_h^2 / (2 v^2) and mu = M_h / sqrt(2) (at tree level).

**The masses.** Particles get mass by interacting with the field's constant
value: m_W = g v / 2, m_Z = sqrt(g^2 + g'^2) v / 2 for the weak bosons and
m_f = y_f v / sqrt(2) for each fermion, so a particle's coupling to the Higgs
boson is proportional to its mass. ATLAS's measured couplings (Nature 607,
52, 2022) test exactly this.

**The vacuum energy.** At the minimum the potential sits lambda v^4 / 4 below
its value at H = 0: an energy density ~10^8 GeV^4, against the dark energy
density of ~(2.2 meV)^4 that cosmology measures -- the cosmological-constant
problem in one number.

**The calendar.** In the hot early universe the minimum moved back to H = 0.
Lattice simulations (D'Onofrio & Rummukainen 2016) show the Higgs condensate
switching on smoothly -- a crossover, not a phase transition -- at
T_c = 159.6 +- 1.5 GeV. Project P2 (``68-thermal-history``) converts
temperature to time; its time-temperature table is read when present.
"""

from __future__ import annotations

import csv
import json

import numpy as np
from particle import Particle

import config as cfg

P2_DIR = cfg.ROOT_DIR / "68-thermal-history"
P2_CURVES = P2_DIR / "data" / "derived" / "curves.npz"
P2_RESULTS = P2_DIR / "data" / "results.json"

#: D'Onofrio & Rummukainen 2016 (PRD 93, 025003): crossover temperature from
#: the maximum of the condensate susceptibility; D'Onofrio, Rummukainen &
#: Tranberg 2014 (PRL 113, 141602): sphaleron freeze-out.
T_CROSSOVER = (159.6, 1.5)          # GeV (stat 0.1 + effective-theory 1.5)
T_SPHALERON = (131.7, 2.3)          # GeV

# Planck 2018 (TT,TE,EE+lowE+lensing) for the dark-energy density
H0_KMS_MPC = 67.36
OMEGA_L = 0.6847
MPC_KM = 3.0856775814913673e19


def higgs_potential() -> dict:
    gf = cfg.G_FERMI
    v = (np.sqrt(2) * gf) ** -0.5
    mh = cfg.M_HIGGS[0]
    lam = mh**2 / (2 * v**2)
    mu = mh / np.sqrt(2)
    depth = lam * v**4 / 4                       # GeV^4
    h0 = H0_KMS_MPC / MPC_KM * cfg.HBAR_GEV_S    # GeV
    rho_crit = 3 * h0**2 * cfg.REDUCED_PLANCK**2
    rho_l = OMEGA_L * rho_crit
    mw, mz = cfg.M_W[0], cfg.M_Z[0]
    return {"v": v, "lambda_tree": lam, "mu": mu, "depth_GeV4": depth,
            "depth_fourth_root_GeV": depth**0.25,
            "rho_lambda_GeV4": rho_l, "rho_lambda_fourth_root_meV": rho_l**0.25 * 1e12,
            "log10_depth_over_rho_lambda": float(np.log10(depth / rho_l)),
            "g": 2 * mw / v, "g_z": 2 * mz / v, "cos_theta_w_on_shell": mw / mz,
            "sin2_theta_w_on_shell": 1 - (mw / mz) ** 2}


def yukawas() -> list[dict]:
    """Tree-level couplings to the Higgs boson from PDG 2026 masses:
    y_f = sqrt(2) m_f / v for fermions, g_V = sqrt(2) m_V / v ("reduced"
    as in ATLAS's plot: kappa_V^(1/2) m_V / v)."""
    v = (np.sqrt(2) * cfg.G_FERMI) ** -0.5
    rows = []
    for name, pid, kind in [("t", 6, "f"), ("b", 5, "f"), ("c", 4, "f"), ("τ", 15, "f"),
                            ("μ", 13, "f"), ("s", 3, "f"), ("e", 11, "f"),
                            ("W", 24, "V"), ("Z", 23, "V")]:
        m = Particle.from_pdgid(pid).mass / 1000
        rows.append({"particle": name, "mass_GeV": m, "kind": kind,
                     "reduced_coupling": m / v, "yukawa": np.sqrt(2) * m / v if kind == "f"
                     else None})
    return rows


def atlas_couplings() -> list[dict]:
    """ATLAS 2022 (Nature 607, 52), HEPData Figure 5: kappa_F m_F / v and
    sqrt(kappa_V) m_V / v with their 68% intervals."""
    path = cfg.LITERATURE_DIR / "atlas2022_reduced_couplings.csv"
    if not path.exists():
        return []
    blocks, cur = [], []
    for line in path.read_text().splitlines():
        if line.startswith("#:"):
            continue
        if not line.strip():
            if cur:
                blocks.append(cur)
            cur = []
            continue
        cur.append(next(csv.reader([line])))
    if cur:
        blocks.append(cur)
    masses = {r[0]: float(r[1]) for r in blocks[0][1:]}
    red = {r[0]: (float(r[1]), float(r[2]), float(r[3])) for r in blocks[1][1:]}
    kap = {r[0]: (float(r[1]), float(r[2]), float(r[3])) for r in blocks[2][1:]}
    names = {"$\\kappa_{t}$": "t", "$\\kappa_{Z}$": "Z", "$\\kappa_{W}$": "W",
             "$\\kappa_{b}$": "b", "$\\kappa_{\\tau}$": "τ", "$\\kappa_{\\mu}$": "μ"}
    out = []
    for key, nm in names.items():
        out.append({"particle": nm, "mass_GeV": masses[key], "reduced": red[key][0],
                    "reduced_up": red[key][1], "reduced_down": red[key][2],
                    "kappa": kap[key][0], "kappa_up": kap[key][1], "kappa_down": kap[key][2]})
    return out


def lattice_condensate() -> dict:
    """<phi^dagger phi>/T against T across the crossover (the authors' table)."""
    path = cfg.LITERATURE_DIR / "donofrio_rummukainen2016_sm_data.txt"
    a = np.loadtxt(path, comments="#")
    return {"T": a[:, 1], "phi2_over_T": a[:, 3], "err": a[:, 4], "chi": a[:, 5],
            "chi_err": a[:, 6]}


def p2_time_of(T_gev) -> dict:
    """Time after the Big Bang at temperature T, from P2's table if present,
    else the radiation-era formula with the full Standard Model g* = 106.75."""
    T = np.atleast_1d(np.asarray(T_gev, dtype=float))
    if P2_CURVES.exists():
        z = np.load(P2_CURVES)
        cT, ct = z["cal_T"], z["cal_t"]
        o = np.argsort(cT)
        t = np.exp(np.interp(np.log(T), np.log(cT[o]), np.log(ct[o])))
        return {"t_s": t, "source": "P2 (68-thermal-history) time-temperature table"}
    gstar = 106.75
    hbar = cfg.HBAR_GEV_S
    t = 0.301 * gstar**-0.5 * cfg.M_PLANCK / T**2 * hbar
    return {"t_s": t, "source": "radiation era, g* = 106.75 (P2 not available)"}


def p2_half_temperatures() -> dict:
    if not P2_RESULTS.exists():
        return {}
    r = json.loads(P2_RESULTS.read_text())
    return {k: v["T_half_GeV"] for k, v in r.get("half_temperatures", {}).items()}


def calendar() -> dict:
    """The electroweak epoch: when the Higgs field switched on and what followed."""
    ev = {"crossover": T_CROSSOVER[0], "crossover_hi": T_CROSSOVER[0] + T_CROSSOVER[1],
          "crossover_lo": T_CROSSOVER[0] - T_CROSSOVER[1], "sphaleron_freeze_out": T_SPHALERON[0]}
    half = p2_half_temperatures()
    for k, v in half.items():
        ev[f"{k}_half"] = v
    times = p2_time_of(list(ev.values()))
    out = {k: {"T_GeV": T, "t_s": float(t)} for (k, T), t in zip(ev.items(), times["t_s"])}
    out["source_of_times"] = times["source"]
    return out
