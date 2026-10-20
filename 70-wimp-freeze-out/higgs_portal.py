"""A thermal relic on the direct-detection plane: the Higgs-portal scalar singlet.

"The WIMP miracle" fixes an annihilation rate, not a scattering rate. To put a thermal relic on
the direct-detection plane one needs a model that ties the two together. The simplest is a real
scalar S, odd under a Z2 symmetry (so it is stable), whose only contact with ordinary matter is the
Higgs boson (Silveira & Zee 1985; McDonald 1994; Burgess, Pospelov & ter Veldhuis 2001):

    V  ⊃  (1/2) mu_S^2 S^2 + (1/2) lambda_hs S^2 |H|^2.

One coupling, lambda_hs, controls both:

* annihilation SS -> h* -> (everything a Higgs of mass sqrt(s) decays to), plus SS -> hh, which
  sets the relic density through freeze-out; and
* scattering on nuclei by Higgs exchange (spin-independent),
  sigma_SI = lambda_hs^2 f_N^2 mu^2 m_n^2 / (4 pi m_h^4 m_S^2).

So demanding Omega h^2 = 0.120 picks lambda_hs(m_S), and that predicts sigma_SI(m_S): a line on
the LZ plot. We follow Cline, Kainulainen, Scott & Weniger (2013, PRD 88, 055025; appendix of the
2015 version): for sqrt(s) < 300 GeV the width of a Standard-Model Higgs of mass sqrt(s) from the
LHC Higgs Working Group (Yellow Report 4 "BSM width" table, 20 GeV - 1 TeV); above, tree-level
widths into W, Z and top; the full SS -> hh cross section; and the exact thermal average
(Gondolo & Gelmini 1991), which matters near the 4-MeV-wide Higgs resonance at m_S ~ m_h/2.
The relic density then comes from our own Boltzmann solver with the Saikawa-Shirai g*(T).

Caveat (labelled in the README): near the resonance, dark matter can lose kinetic equilibrium
with the plasma before chemical freeze-out, which changes the relic density by up to an order of
magnitude for m_S just below m_h/2 (Binder, Bringmann, Gustafsson & Hryczuk 2017). We assume
kinetic equilibrium throughout, as Cline et al. did.
"""

from __future__ import annotations

import math
import re
import zipfile
from functools import lru_cache

import numpy as np
from scipy import integrate, interpolate, optimize, special

import config as C
import freezeout as fo

V0 = 246.22            # Higgs vacuum expectation value, GeV (from G_F)
M_H = 125.20           # PDG 2024
M_W = 80.369
M_Z = 91.188
M_T = 172.57
M_N_SI = 0.946         # nucleon mass used by Cline et al. in sigma_SI
F_N = 0.30             # Higgs-nucleon coupling factor, 0.30 +- 0.03 (Cline et al. 2013)
F_N_ERR = 0.03
ALPHA_S = 0.12         # as Cline et al. for the QCD correction to quark final states
#: invisible Higgs decays: ATLAS Run 1+2 combination, BR(h -> invisible) < 0.107 (95% CL)
BR_INV_MAX = 0.107
#: fermions for the tree-level width (mass in GeV, colour factor); MS-bar masses for light quarks
FERMIONS = {"b": (4.18, 3), "c": (1.27, 3), "tau": (1.77693, 1), "t": (M_T, 3), "mu": (0.1056584, 1)}


@lru_cache(maxsize=None)
def yr4_widths() -> dict:
    """Partial widths (GeV) of a SM-like Higgs vs its mass from the YR4 'BSM Width' sheet."""
    z = zipfile.ZipFile(C.RAW_DIR / "higgs" / "yr4_Higgs_XSBR_update.xlsx")
    wb = z.read("xl/workbook.xml").decode()
    rels = z.read("xl/_rels/workbook.xml.rels").decode()
    sheets = dict(re.findall(r'<sheet [^>]*name="([^"]*)"[^>]*r:id="([^"]*)"', wb))
    relmap = {}
    for m in re.finditer(r"<Relationship ([^>]*)/>", rels):
        a = m.group(1)
        relmap[re.search(r'Id="([^"]*)"', a).group(1)] = re.search(r'Target="([^"]*)"', a).group(1)
    x = z.read("xl/" + relmap[sheets["YR4 BSM Width"]].split("xl/")[-1]).decode()
    cols = {"A": "M", "B": "bb", "E": "tautau", "H": "mumu", "K": "cc", "N": "ss", "Q": "tt",
            "V": "gg", "Y": "gamgam", "AB": "Zgam", "AE": "WW", "AH": "ZZ"}
    rows = []
    for r in re.findall(r"<row [^>]*>(.*?)</row>", x, re.S):
        vals = {}
        for attrs, inner in re.findall(r"<c ([^>]*?)(?:/>|>(.*?)</c>)", r, re.S):
            v = re.search(r"<v>(.*?)</v>", inner or "")
            col = re.search(r'r="([A-Z]+)\d+"', attrs).group(1)
            if v and col in cols and 't="s"' not in attrs:
                vals[cols[col]] = float(v.group(1))
        if "M" in vals and "bb" in vals:
            rows.append(vals)
    out = {k: np.array([r.get(k, 0.0) for r in rows]) for k in cols.values()}
    out["total"] = sum(out[k] for k in cols.values() if k != "M")
    return out


@lru_cache(maxsize=None)
def _gamma_table_spline():
    t = yr4_widths()
    return interpolate.PchipInterpolator(np.log(t["M"]), np.log(t["total"]))


def gamma_sm_total(M: float) -> float:
    """Total width (GeV) of a SM Higgs of mass M, interpolated in the YR4 table (20 GeV - 1 TeV)."""
    return float(np.exp(_gamma_table_spline()(math.log(M))))


def gamma_h() -> float:
    return gamma_sm_total(M_H)


def _D2(s: float, gam_prop: float) -> float:
    return 1.0 / ((s - M_H ** 2) ** 2 + M_H ** 2 * gam_prop ** 2)


def _xq(mq: float, s: float) -> float:
    return 3.0 * (1.0 + (1.5 * math.log(mq * mq / s) + 2.25) * 4 * ALPHA_S / (3 * math.pi))


def sigv_sm(s: float, lam: float, gam_prop: float) -> float:
    """v_rel sigma(SS -> SM final states except hh), GeV^-2."""
    rs = math.sqrt(s)
    D2 = _D2(s, gam_prop)
    if 20.0 <= rs < 300.0:
        return 2 * lam ** 2 * V0 ** 2 * D2 * gamma_sm_total(rs) / rs
    tot = 0.0
    for MV, dv in ((M_W, 1.0), (M_Z, 0.5)):
        x = MV * MV / s
        if x < 0.25:
            tot += lam ** 2 * s / (8 * math.pi) * dv * math.sqrt(1 - 4 * x) * D2 * (1 - 4 * x + 12 * x * x)
    for name, (mf, nc) in FERMIONS.items():
        if 4 * mf * mf < s:
            X = _xq(mf, s) if nc == 3 else 1.0
            tot += lam ** 2 * mf ** 2 / (4 * math.pi) * X * (1 - 4 * mf * mf / s) ** 1.5 * D2
    return tot


def sigv_hh(s: float, m: float, lam: float, gam_prop: float) -> float:
    """v_rel sigma(SS -> hh), GeV^-2: contact + s-channel Higgs + t/u-channel S (Cline et al. A.4)."""
    if s <= 4 * M_H ** 2:
        return 0.0
    vS = math.sqrt(max(1 - 4 * m * m / s, 1e-300))
    vh = math.sqrt(1 - 4 * M_H ** 2 / s)
    D2 = _D2(s, gam_prop)
    aR = 1 + 3 * M_H ** 2 * (s - M_H ** 2) * D2
    aI = 3 * M_H ** 2 * math.sqrt(s) * gam_prop * D2
    tp = m * m + M_H ** 2 - 0.5 * s * (1 - vS * vh)
    tm = m * m + M_H ** 2 - 0.5 * s * (1 + vS * vh)
    term = ((aR ** 2 + aI ** 2) * s * vS * vh
            + 4 * lam * V0 ** 2 * (aR - lam * V0 ** 2 / (s - 2 * M_H ** 2)) * math.log(abs((m * m - tp) / (m * m - tm)))
            + 2 * lam ** 2 * V0 ** 4 * s * vS * vh / ((m * m - tm) * (m * m - tp)))
    return lam ** 2 / (16 * math.pi * s * s * vS) * term


def propagator_width(m: float, lam: float) -> float:
    """Higgs width in the propagator: SM plus h -> SS when open."""
    g = gamma_h()
    if 2 * m < M_H:
        g += lam ** 2 * V0 ** 2 / (32 * math.pi * M_H) * math.sqrt(1 - 4 * m * m / M_H ** 2)
    return g


def sigv_total(s: float, m: float, lam: float) -> float:
    gp = propagator_width(m, lam)
    return sigv_sm(s, lam, gp) + sigv_hh(s, m, lam, gp)


def thermal_average(m: float, lam: float, x: float) -> float:
    """<sigma v> (GeV^-2) at T = m/x, Gondolo & Gelmini (1991):
    int_{4m^2} s sqrt(s - 4 m^2) K1(sqrt s / T) sigma v_rel ds / (16 T m^4 K2(m/T)^2).
    Written in E = sqrt(s) with scaled Bessel functions; the Higgs pole is a breakpoint."""
    T = m / x
    norm = 16 * T * m ** 4 * special.kve(2, x) ** 2

    def f(E):
        s = E * E
        return (2 * E * s * math.sqrt(max(s - 4 * m * m, 0.0)) * special.kve(1, E / T)
                * math.exp(-(E - 2 * m) / T) * sigv_total(s, m, lam))

    lo, hi = 2 * m, 2 * m + 60 * T
    pts = [p for p in (M_H, 2 * M_H, 2 * M_W, 2 * M_Z, 2 * M_T, 300.0) if lo < p < hi]
    gp = propagator_width(m, lam)
    if lo < M_H < hi:
        pts += [M_H - 20 * gp, M_H + 20 * gp, M_H - 2000 * gp, M_H + 2000 * gp]
        pts = [p for p in pts if lo < p < hi]
    val = integrate.quad(f, lo, hi, points=sorted(set(pts)) or None, limit=500, epsabs=0, epsrel=1e-7)[0]
    return val / norm


#: x grid of the tabulated thermal average (constant beyond its ends)
LNX = np.log(np.geomspace(3.0, 3000.0, 70))


def sigv_table(m: float, lam: float) -> tuple:
    vals = np.array([thermal_average(m, lam, float(np.exp(lx))) for lx in LNX])
    return LNX, np.log(np.maximum(vals, 1e-300))


def omega_h2(m: float, lam: float) -> float:
    r = fo.Relic(m, g=1.0, sigv_table=sigv_table(m, lam))
    return fo.solve(r, rtol=1e-8).omega_h2


def relic_coupling(m: float, target: float = C.OMEGA_C_H2) -> dict:
    """lambda_hs giving Omega h^2 = target (Newton in ln lambda, Omega ~ lambda^-2)."""
    def f(lnl):
        return math.log(omega_h2(m, math.exp(lnl)) / target)

    # a first guess from the s-wave thermal value at threshold
    s0 = 4 * m * m * (1 + 1e-6)
    g1 = sigv_total(s0, m, 1.0)
    lam0 = math.sqrt(2.2e-26 / C.GEV2_TO_CM3S / max(g1, 1e-300)) if g1 > 0 else 0.1
    lam0 = min(max(lam0, 1e-5), 20.0)
    a = math.log(lam0)
    fa = f(a)
    b = a - fa / 2.0
    fb = f(b)
    for _ in range(30):
        if abs(fb) < 1e-7:
            break
        slope = (fb - fa) / (b - a) if b != a else -2.0
        if not -6.0 < slope < -0.3:
            slope = -2.0
        a, fa = b, fb
        b = b - fb / slope
        fb = f(b)
    lam = math.exp(b)
    return {"m": m, "lambda": lam, "omega_h2": target * math.exp(fb), "sigma_si_cm2": sigma_si(m, lam)}


def sigma_si(m: float, lam: float, f_n: float = F_N) -> float:
    """Spin-independent S-nucleon cross section (cm^2)."""
    mu = M_N_SI * m / (M_N_SI + m)
    return lam ** 2 * f_n ** 2 * mu ** 2 * M_N_SI ** 2 / (4 * math.pi * M_H ** 4 * m ** 2) * C.GEV2_TO_CM2


def br_invisible(m: float, lam: float) -> float:
    if 2 * m >= M_H:
        return 0.0
    gi = lam ** 2 * V0 ** 2 / (32 * math.pi * M_H) * math.sqrt(1 - 4 * m * m / M_H ** 2)
    return gi / (gi + gamma_h())


def lambda_invisible_max(m: float) -> float:
    """Largest lambda allowed by BR(h -> invisible) < BR_INV_MAX (only for m < m_h/2)."""
    if 2 * m >= M_H:
        return float("inf")
    gi = BR_INV_MAX / (1 - BR_INV_MAX) * gamma_h()
    return math.sqrt(gi * 32 * math.pi * M_H / (V0 ** 2 * math.sqrt(1 - 4 * m * m / M_H ** 2)))
