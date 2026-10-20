"""Strongly interacting matter: hadron gas, lattice QCD and perturbative QCD.

Everything here is about the quarks, gluons and hadrons together ("the strong sector"). Three
descriptions, each valid in its own temperature range:

* **Hadron-resonance gas (HRG)**, T below ~120 MeV: an ideal gas of every hadron in the PDG table.
  Interactions are mostly accounted for by the resonances themselves (Dashen, Ma & Bernstein 1969).
* **Lattice QCD**, ~100 MeV to 1 GeV: the continuum equation of state from computer simulations.
  We use three published results: Borsanyi et al. 2014 (2+1 flavours, table), HotQCD 2014
  (2+1 flavours, table) and Borsanyi et al. 2016 (2+1+1 flavours, i.e. with charm, a fit).
* **Perturbative QCD**, T above ~1 GeV: the pressure as a series in the strong coupling to order
  g^6 (Kajantie et al. 2003, as written by Saikawa & Shirai 2018, eq. 2.27), with the coupling
  run at four loops and the charm, bottom and top masses brought in with the tree-level factors of
  Borsanyi et al. 2016.

The central quantity is the trace anomaly Delta = (rho - 3p)/T^4 = T d(p/T^4)/dT. The pieces are
blended in Delta and integrated for p/T^4, so the result is thermodynamically consistent
(s = dp/dT and rho = Ts - p hold by construction).
"""

from __future__ import annotations

import re
from functools import lru_cache

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

import config as C
import particles
import thermo

PI2 = np.pi ** 2
ZETA3 = 1.2020569031595942


# --- Hadron-resonance gas ---------------------------------------------------------------


def hrg(T, hs=None) -> dict:
    """p/T^4, rho/T^4, s/T^3 and Delta of an ideal gas of hadrons (GeV units for T)."""
    T = np.atleast_1d(np.asarray(T, dtype=float))
    hs = particles.hadrons() if hs is None else hs
    p = np.zeros_like(T)
    r = np.zeros_like(T)
    for sp in hs:
        f = thermo.functions(sp.mass / T, sp.stat)
        p += sp.g * f["p"]
        r += sp.g * f["rho"]
    return {"p": p, "rho": r, "s": r + p, "delta": r - 3 * p}


# --- Lattice QCD: published continuum results ------------------------------------------


@lru_cache(maxsize=None)
def wb2014() -> dict:
    """Borsanyi et al. 2014 (PLB 730, 99), 2+1 flavours: T [GeV], I/T^4, p/T^4, e/T^4, s/T^3 with errors."""
    a = np.loadtxt(C.RAW_DIR / "wb2014_eos.dat")
    keys = ["T", "delta", "delta_err", "p", "p_err", "rho", "rho_err", "s", "s_err", "cs2", "cs2_err"]
    d = {k: a[:, i] for i, k in enumerate(keys)}
    d["T"] = d["T"] / 1000.0
    return d


def _parse_asym(tok: str) -> tuple[float, float, float]:
    """'0.439(-44)(+65)' -> (0.439, 0.044, 0.065); '16.5(-0.9)(+3.0)' -> (16.5, 0.9, 3.0)."""
    m = re.match(r"\s*([0-9.]+)\(([-+0-9.]+)\)\(([-+0-9.]+)\)", tok)
    if not m:
        raise ValueError(tok)
    val, lo, hi = m.groups()
    dec = len(val.split(".")[1]) if "." in val else 0
    def scale(e):
        e = e.lstrip("+-")
        return float(e) if "." in e else float(e) * 10 ** (-dec)
    return float(val), scale(lo), scale(hi)


@lru_cache(maxsize=None)
def hotqcd2014() -> dict:
    """HotQCD 2014 (Bazavov et al., PRD 90, 094503), Table I: 2+1 flavours, continuum, 130-400 MeV."""
    text = (C.LIT_DIR / "1407.6387" / "tab_fin.tex").read_text()
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not re.match(r"^\d{3}\s*&", line):
            continue
        cells = [c.strip() for c in line.rstrip("\\").split("&")]
        T = float(cells[0]) / 1000.0
        vals = [_parse_asym(c) for c in cells[1:7]]
        rows.append([T] + [v for trip in vals for v in trip])
    a = np.array(rows)
    names = ["delta", "p", "rho", "s", "cv", "cs2"]
    d = {"T": a[:, 0]}
    for i, n in enumerate(names):
        d[n], d[n + "_lo"], d[n + "_hi"] = a[:, 1 + 3 * i], a[:, 2 + 3 * i], a[:, 3 + 3 * i]
    return d


@lru_cache(maxsize=None)
def borsanyi2016_params() -> dict:
    """Coefficients of Borsanyi et al. 2016's trace-anomaly fit (their Table S3), read from the TeX."""
    text = (C.LIT_DIR / "1606.07494" / "som_eos.tex").read_text()
    out = {}
    for key, pat in (("2+1+1", r"2\+1\+1 flavors\}\s*&(.*?)\\\\"), ("2+1", r"2\+1 flavors\}\s*&(.*?)\\\\")):
        m = re.search(pat, text)
        vals = [float(v) for v in m.group(1).split("&")]
        out[key] = dict(zip(["h0", "h1", "h2", "f0", "f1", "f2", "g1", "g2"], vals))
    return out


@lru_cache(maxsize=None)
def borsanyi2016_spline_table() -> dict:
    """Borsanyi et al. 2016's g_rho(T) and g_rho/g_s spline knots for the whole Standard Model."""
    text = (C.LIT_DIR / "1606.07494" / "som_eos.tex").read_text()
    block = text.split(r"\label{tab:spline}")[0].split(r"$\log_{10} T/\mathrm{MeV}$")[-1]
    rows = re.findall(r"^\s*([0-9.]+)\s*&\s*([0-9.]+)\s*&\s*([0-9.]+)\s*\\\\", block, flags=re.M)
    a = np.array(rows, dtype=float)
    return {"T": 10 ** a[:, 0] / 1000.0, "g_rho": a[:, 1], "ratio": a[:, 2], "g_s": a[:, 1] / a[:, 2]}


def borsanyi2016_delta(T, flavours: str = "2+1+1") -> np.ndarray:
    """Trace anomaly I/T^4 from Borsanyi et al. 2016's fit (valid 100-1000 MeV); t = T/200 MeV."""
    c = borsanyi2016_params()[flavours]
    t = np.asarray(T, dtype=float) / 0.200
    return np.exp(-c["h1"] / t - c["h2"] / t ** 2) * (
        c["h0"] + c["f0"] * (np.tanh(c["f1"] * t + c["f2"]) + 1) / (1 + c["g1"] * t + c["g2"] * t ** 2))


# --- Perturbative QCD ----------------------------------------------------------------------


def _beta(nf: int) -> tuple[float, float, float, float]:
    """MS-bar beta-function coefficients for a = alpha_s/(4 pi): mu^2 da/dmu^2 = -sum b_i a^(i+2)
    (four loops: van Ritbergen, Vermaseren & Larin 1997)."""
    b0 = 11 - 2 / 3 * nf
    b1 = 102 - 38 / 3 * nf
    b2 = 2857 / 2 - 5033 / 18 * nf + 325 / 54 * nf ** 2
    b3 = (149753 / 6 + 3564 * ZETA3) - (1078361 / 162 + 6508 / 27 * ZETA3) * nf \
        + (50065 / 162 + 6472 / 81 * ZETA3) * nf ** 2 + 1093 / 729 * nf ** 3
    return b0, b1, b2, b3


def _run(a0: float, mu0: float, mu1: float, nf: int) -> float:
    """Run a = alpha_s/(4pi) from mu0 to mu1 at fixed nf (four loops)."""
    if mu0 == mu1:
        return a0
    b = _beta(nf)

    def rhs(lnmu2, y):
        a = y[0]
        return [-(b[0] * a ** 2 + b[1] * a ** 3 + b[2] * a ** 4 + b[3] * a ** 5)]

    sol = solve_ivp(rhs, (2 * np.log(mu0), 2 * np.log(mu1)), [a0], rtol=1e-11, atol=1e-14)
    return float(sol.y[0, -1])


#: MS-bar quark masses m_q(m_q) used as flavour thresholds (PDG: 1.27, 4.18, 162.5 GeV).
THRESHOLDS = {4: 1.27, 5: 4.18, 6: 162.5}


@lru_cache(maxsize=None)
def _alpha_table(nf: int):
    """alpha_s in the nf-flavour theory as a spline in log(mu), mu = 0.5 GeV ... 1e6 GeV.

    Starts from alpha_s(M_Z) in the 5-flavour theory and matches continuously at the thresholds.
    """
    a_mz = C.ALPHA_S_MZ / (4 * np.pi)
    # coupling of the nf-flavour theory at a reference scale
    if nf == 5:
        ref_mu, ref_a = C.M_Z, a_mz
    elif nf == 4:
        ref_mu, ref_a = THRESHOLDS[5], _run(a_mz, C.M_Z, THRESHOLDS[5], 5)
    elif nf == 3:
        a_b = _run(a_mz, C.M_Z, THRESHOLDS[5], 5)
        ref_mu, ref_a = THRESHOLDS[4], _run(a_b, THRESHOLDS[5], THRESHOLDS[4], 4)
    elif nf == 6:
        ref_mu, ref_a = THRESHOLDS[6], _run(a_mz, C.M_Z, THRESHOLDS[6], 5)
    else:
        raise ValueError(nf)
    lmu = np.linspace(np.log(0.8), np.log(1e7), 600)
    b = _beta(nf)

    def rhs(lnmu, y):
        a = y[0]
        return [-2 * (b[0] * a ** 2 + b[1] * a ** 3 + b[2] * a ** 4 + b[3] * a ** 5)]

    lr = np.log(ref_mu)
    up = solve_ivp(rhs, (lr, lmu[-1]), [ref_a], dense_output=True, rtol=1e-11, atol=1e-15)
    dn = solve_ivp(rhs, (lr, lmu[0]), [ref_a], dense_output=True, rtol=1e-11, atol=1e-15)
    a = np.where(lmu >= lr, up.sol(np.maximum(lmu, lr))[0], dn.sol(np.minimum(lmu, lr))[0])
    return CubicSpline(lmu, 4 * np.pi * a)


def alpha_s(mu, nf: int) -> np.ndarray:
    """Strong coupling alpha_s(mu) of the nf-flavour theory."""
    return _alpha_table(nf)(np.log(np.asarray(mu, dtype=float)))


def p_pqcd_massless(T, nf: int, mu_factor: float = 1.0, qc: float = C.QC_PQCD) -> np.ndarray:
    """p/T^4 of QCD with nf massless quarks to order g^6 (Kajantie et al. 2003; S&S eq. 2.27).

    The renormalisation scale is mu = mu_factor x 2 pi T.
    """
    T = np.asarray(T, dtype=float)
    mu = mu_factor * 2 * np.pi * T
    a = alpha_s(mu, nf) / np.pi
    L = np.log(mu_factor)  # ln(mu / 2 pi T)
    n = nf
    p0 = 1 + 21 / 32 * n
    p2 = -15 / 4 * (1 + 5 / 12 * n)
    p3 = 30 * (1 + n / 6) ** 1.5
    p4 = (237.2 + 15.96 * n - 0.4150 * n ** 2 + 135 / 2 * (1 + n / 6) * np.log(a * (1 + n / 6))
          - 165 / 8 * (1 + 5 / 12 * n) * (1 - 2 / 33 * n) * L)
    p5 = (1 + n / 6) ** 0.5 * (-799.1 - 21.96 * n - 1.926 * n ** 2
                              + 495 / 2 * (1 + n / 6) * (1 - 2 / 33 * n) * L)
    p6 = ((-659.2 - 65.89 * n - 7.653 * n ** 2 + 1485 / 2 * (1 + n / 6) * (1 - 2 / 33 * n) * L)
          * np.log(a * (1 + n / 6)) - 475.6 * np.log(a)
          - 1815 / 16 * (1 + 5 / 12 * n) * (1 - 2 / 33 * n) ** 2 * L ** 2
          + (2932.9 + 42.83 * n - 16.48 * n ** 2 + 0.2767 * n ** 3) * L + qc)
    return 8 * PI2 / 45 * (p0 + p2 * a + p3 * a ** 1.5 + p4 * a ** 2 + p5 * a ** 2.5 + p6 * a ** 3)


def p_sb(nf: int) -> float:
    """Stefan-Boltzmann p/T^4 of nf massless quark flavours plus gluons."""
    return PI2 / 90 * (16 + 10.5 * nf)


def _free_quark_p(T, m) -> np.ndarray:
    """p/T^4 of one free quark flavour (12 states)."""
    return 12 * thermo.functions(m / np.asarray(T, dtype=float), -1)["p"]


def p_pqcd(T, mu_factor: float = 1.0, qc4: float = C.QC_PQCD, qc6: float = C.QC_PQCD,
           masses: dict | None = None) -> np.ndarray:
    """p/T^4 of the strong sector at high T: u, d, s massless, with c, b and t mass thresholds.

    Following Borsanyi et al. 2016 and Saikawa & Shirai 2018 (eqs. 2.31-2.34):
    p(udsc)  = [p_SB(3) + p_c^0] / p_SB(4) x p_QCD(nf=4)
    p(udscb) = [p_SB(4) + p_b^0] / p_SB(4) x p(udsc)
    p(all)   = c_t p_QCD(nf=6) + (1 - c_t) p(udscb), c_t = p_t^0 / p_t^0(massless).
    """
    T = np.asarray(T, dtype=float)
    m = {4: particles.pdg_mass(4), 5: particles.pdg_mass(5), 6: particles.pdg_mass(6)}
    if masses:
        m.update(masses)
    p4 = p_pqcd_massless(T, 4, mu_factor, qc4)
    p_udsc = (p_sb(3) + _free_quark_p(T, m[4])) / p_sb(4) * p4
    p_udscb = (p_sb(4) + _free_quark_p(T, m[5])) / p_sb(4) * p_udsc
    c_t = _free_quark_p(T, m[6]) / (12 * 7 / 8 * PI2 / 90)
    return c_t * p_pqcd_massless(T, 6, mu_factor, qc6) + (1 - c_t) * p_udscb


# --- The strong sector, assembled -----------------------------------------------------------


def _smoothstep(lt, l0, l1):
    """0 below l0, 1 above l1, a smooth (C2) ramp in between (in log T)."""
    u = np.clip((lt - l0) / (l1 - l0), 0, 1)
    return u ** 3 * (10 - 15 * u + 6 * u ** 2)


def lattice_delta(T, lattice: str = "borsanyi2016") -> np.ndarray:
    """Trace anomaly of one lattice result inside its published range (0 outside)."""
    T = np.asarray(T, dtype=float)
    if lattice in ("borsanyi2016", "borsanyi2016_21"):
        ok = (T >= 0.08) & (T <= 1.5)  # the fits are published for 100-1000 MeV
        d = np.zeros_like(T)
        d[ok] = borsanyi2016_delta(T[ok], "2+1+1" if lattice == "borsanyi2016" else "2+1")
        return d
    tab = wb2014() if lattice == "wb2014" else hotqcd2014()
    sp = CubicSpline(np.log(tab["T"]), tab["delta"])
    ok = (T >= tab["T"][0]) & (T <= tab["T"][-1])
    return np.where(ok, sp(np.log(np.clip(T, tab["T"][0], tab["T"][-1]))), 0.0)


def strong_sector(T, lattice: str = "borsanyi2016", mu_factor: float = 1.0, qc: float = C.QC_PQCD,
                  t_hrg: float = C.T_HRG_TO_LATTICE, t_pqcd: float = C.T_LATTICE_TO_PQCD,
                  hs=None, hrg_dex: float = 0.04, pqcd_dex: float = 0.15, masses: dict | None = None) -> dict:
    """p/T^4, rho/T^4, s/T^3 and Delta of quarks, gluons and hadrons over the whole range.

    Low side: the trace anomaly of the hadron gas, handed over to lattice QCD at ``t_hrg`` (a smooth
    ramp over +- ``hrg_dex`` in log10 T), is integrated for p/T^4 upward from the bottom of the
    grid, where the hadron gas holds. High side: the perturbative pressure itself. The two
    pressures are blended around ``t_pqcd`` (+- ``pqcd_dex``), and Delta = d(p/T^4)/d ln T is
    taken from the result, so s = dp/dT and rho = Ts - p hold everywhere. ``qc=None`` fits the
    unknown g^6 constant so that the two pressures agree at ``t_pqcd``.
    """
    T = np.asarray(T, dtype=float)
    if T[0] > 0.01:
        raise ValueError("strong_sector needs a grid starting at or below 10 MeV")
    lt, lnT = np.log10(T), np.log(T)
    h = hrg(T, hs)
    d_lat = lattice_delta(T, lattice)
    w1 = _smoothstep(lt, np.log10(t_hrg) - hrg_dex, np.log10(t_hrg) + hrg_dex)
    d_low = (1 - w1) * h["delta"] + w1 * d_lat
    anti = CubicSpline(lnT, d_low).antiderivative()
    p_low = h["p"][0] + anti(lnT) - anti(lnT[0])
    if qc is None:
        # fix the unknown g^6 constant so that perturbation theory meets the lattice at t_pqcd
        from scipy.optimize import brentq

        i = int(np.argmin(np.abs(lnT - np.log(t_pqcd))))
        target = p_low[i]
        qc = brentq(lambda q: float(p_pqcd(T[i:i + 1], mu_factor, q, q, masses)[0]) - target, -1e5, 1e5)
    hi = T > 0.2
    p_high = np.full_like(T, np.nan)
    p_high[hi] = p_pqcd(T[hi], mu_factor, qc, qc, masses)
    w2 = _smoothstep(lt, np.log10(t_pqcd) - pqcd_dex, np.log10(t_pqcd) + pqcd_dex)
    p = np.where(w2 > 0, (1 - w2) * p_low + w2 * np.nan_to_num(p_high), p_low)
    delta = CubicSpline(lnT, p)(lnT, 1)
    rho = delta + 3 * p
    return {"p": p, "rho": rho, "s": rho + p, "delta": delta,
            "delta_hrg": h["delta"], "p_hrg": h["p"], "delta_lattice": d_lat,
            "p_low": p_low, "p_pqcd": p_high, "qc": qc}
