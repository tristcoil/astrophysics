"""Lattice QCD: the equation of state and the crossover line, read from the
papers' own sources.

* HotQCD 2014 (Bazavov et al., PRD 90, 094503): continuum-extrapolated
  pressure, energy and entropy densities and trace anomaly, T = 130-400 MeV,
  Table I (``tab_fin.tex``, errors written ``0.439(-44)(+65)``), and the
  smooth ansatz for p/T^4 with its parameters (Eq. 16, Table II).
* Wuppertal-Budapest 2014 (Borsanyi et al., PLB 730, 99): the same
  quantities, T = 110-510 MeV, from the paper's ancillary file.
* The crossover line T_pc(mu_B) = T_pc(0) [1 - kappa2 (mu_B/T)^2 - kappa4 (mu_B/T)^4]
  of HotQCD 2019 and Borsanyi et al. 2020 (numbers in config.py, checked in
  the TeX of both papers).

The free quark-gluon gas (Stefan-Boltzmann limit) for three massless flavours:
p/T^4 = (pi^2/90)(16 + 21 N_f / 2) = 95 pi^2 / 180 = 5.21.
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np

import config as cfg
import fetch_data

SB_P = 95 * np.pi ** 2 / 180          # massless u, d, s quarks and gluons


def _val_err(tok: str) -> tuple:
    """'0.439(-44)(+65)' -> (0.439, 0.044, 0.065); '16.5(-0.9)(+3.0)' too."""
    m = re.match(r"\s*([0-9.]+)\(-([0-9.]+)\)\(\+([0-9.]+)\)", tok)
    v, lo, hi = m.group(1), m.group(2), m.group(3)
    val = float(v)
    if "." in lo:                      # errors written in absolute units
        return val, float(lo), float(hi)
    dec = len(v.split(".")[1]) if "." in v else 0
    scale = 10.0 ** (-dec)
    return val, float(lo) * scale, float(hi) * scale


@lru_cache(maxsize=None)
def hotqcd_2014() -> dict:
    tex = (fetch_data.arxiv_dir("1407.6387") / "tab_fin.tex").read_text()
    rows = []
    for line in tex.splitlines():
        if not re.match(r"^\s*\d{3}\s*&", line):
            continue
        cells = [c.strip().rstrip("\\").strip() for c in line.split("&")]
        T = float(cells[0])
        vals = [_val_err(c) for c in cells[1:7]]
        rows.append((T, *vals))
    rows.sort(key=lambda r: r[0])
    out = {"T": np.array([r[0] for r in rows])}
    for j, name in enumerate(("I", "p", "e", "s", "cv", "cs2")):
        out[name] = np.array([r[1 + j][0] for r in rows])
        out[name + "_err"] = np.array([0.5 * (r[1 + j][1] + r[1 + j][2]) for r in rows])
    return out


@lru_cache(maxsize=None)
def hotqcd_2014_ansatz_params() -> dict:
    tex = (fetch_data.arxiv_dir("1407.6387") / "sec4_thermo.tex").read_text()
    m = re.search(r"\$c_t\$&\$a_n\$&\$b_n\$&\$c_n\$&\$d_n\$ \\\\\s*\n([^\n]+)\\\\.*?"
                  r"\$t_0\$&\$a_d\$&\$b_d\$&\$c_d\$&\$d_d\$ \\\\\s*\n([^\n]+)\\\\", tex, re.S)
    a = [float(x) for x in m.group(1).split("&")]
    b = [float(x) for x in m.group(2).split("&")]
    return dict(zip(("ct", "an", "bn", "cn", "dn"), a)) | dict(zip(("t0", "ad", "bd", "cd", "dd"), b))


def hotqcd_pressure(T_mev: np.ndarray) -> np.ndarray:
    """HotQCD's smooth ansatz for p/T^4 (Tc = 154 MeV normalisation), 100-400 MeV."""
    q = hotqcd_2014_ansatz_params()
    t = np.asarray(T_mev, float) / 154.0
    num = SB_P + q["an"] / t + q["bn"] / t ** 2 + q["cn"] / t ** 3 + q["dn"] / t ** 4
    den = 1 + q["ad"] / t + q["bd"] / t ** 2 + q["cd"] / t ** 3 + q["dd"] / t ** 4
    return 0.5 * (1 + np.tanh(q["ct"] * (t - q["t0"]))) * num / den


def hotqcd_ansatz_thermo(T_mev: np.ndarray) -> dict:
    """p, e, s, I over T^4 (T^3 for s) from the ansatz by exact differentiation."""
    T = np.asarray(T_mev, float)
    h = 1e-3 * T
    p = hotqcd_pressure(T)
    dp = (hotqcd_pressure(T + h) - hotqcd_pressure(T - h)) / (2 * h)
    I = T * dp                      # (e - 3p)/T^4 = T d(p/T^4)/dT
    return {"p": p, "I": I, "e": I + 3 * p, "s": I + 4 * p}


@lru_cache(maxsize=None)
def wb_2014() -> dict:
    path = fetch_data.arxiv_dir("1309.5258") / "anc" / "WB-EoS.dat"
    d = np.loadtxt(path)
    names = ("T", "I", "I_err", "p", "p_err", "e", "e_err", "s", "s_err", "cs2", "cs2_err")
    return {n: d[:, i] for i, n in enumerate(names)}


@lru_cache(maxsize=None)
def hotqcd_2012_susceptibilities() -> dict:
    """Continuum second-order charge susceptibilities chi/T^2, 150-250 MeV
    (HotQCD 2012, Table V): BB, QQ, SS, BS, BQ, QS with symmetric errors."""
    tex = (fetch_data.arxiv_dir("1203.0784") / "table_chi2_extrap_rev2.text").read_text()
    rows = []
    for line in tex.splitlines():
        if not re.match(r"^\s*\d{3}\s*&", line):
            continue
        cells = [c.strip().rstrip("\\").strip() for c in line.split("&")]
        vals = []
        for c in cells[1:7]:
            m = re.match(r"(-?[0-9.]+)\(\s*([0-9]+)\)", c)
            v = m.group(1)
            dec = len(v.split(".")[1])
            vals.append((float(v), float(m.group(2)) * 10.0 ** (-dec)))
        rows.append((float(cells[0]), vals))
    out = {"T": np.array([r[0] for r in rows])}
    for j, name in enumerate(("BB", "QQ", "SS", "BS", "BQ", "QS")):
        out[name] = np.array([r[1][j][0] for r in rows])
        out[name + "_err"] = np.array([r[1][j][1] for r in rows])
    return out


# ---------------------------------------------------------------------------
# The crossover line
# ---------------------------------------------------------------------------


def tpc_line(mu_b: np.ndarray, source: str = "hotqcd", n_samples: int = 4000,
             seed: int = 7) -> dict:
    """T_pc(mu_B) with a 68% band from the published errors.

    HotQCD 2019 expands in mu_B / T_pc(0); Borsanyi et al. 2020 in
    mu_B / T_c(mu_B) (solved by iteration). Errors are sampled as independent
    Gaussians: Borsanyi et al. note that kappa2 and kappa4 are correlated, so
    their band is conservative.
    """
    p = cfg.HOTQCD_2019 if source == "hotqcd" else cfg.BORSANYI_2020
    T0 = p["Tpc"] if source == "hotqcd" else p["Tc"]
    T0e = p["Tpc_err"] if source == "hotqcd" else p["Tc_err"]
    rng = np.random.default_rng(seed)
    t0 = rng.normal(T0, T0e, n_samples)
    k2 = rng.normal(p["kappa2"], p["kappa2_err"], n_samples)
    k4 = rng.normal(p["kappa4"], p["kappa4_err"], n_samples)
    mu = np.asarray(mu_b, float)

    def curve(t0_, k2_, k4_):
        if source == "hotqcd":
            x = mu / t0_
            return t0_ * (1 - k2_ * x ** 2 - k4_ * x ** 4)
        T = np.full_like(mu, t0_)
        for _ in range(100):
            x = mu / T
            T = t0_ * (1 - k2_ * x ** 2 - k4_ * x ** 4)
        return T

    central = curve(T0, p["kappa2"], p["kappa4"])
    samples = np.array([curve(a, b, c) for a, b, c in zip(t0, k2, k4)])
    lo, hi = np.percentile(samples, [16, 84], axis=0)
    return {"mu_B": mu, "T": central, "lo": lo, "hi": hi}
