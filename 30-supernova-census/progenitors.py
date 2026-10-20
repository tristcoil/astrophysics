"""The stars that exploded, and the ones that are missing: the red supergiant problem.

For about 30 nearby core-collapse supernovae, archival images taken years
before the explosion show the star that blew up (or set a limit on how bright
it could have been). Type II-P supernovae come from red supergiants. Smartt
(2009, 2015) noticed that none of them was more luminous than log L/L_sun ~ 5.1
(an initial mass of ~18 M_sun), although red supergiants shine up to
log L ~ 5.5 (~25-30 M_sun). If the heaviest red supergiants collapse to black
holes without a bright explosion, this is what we would see.

This module re-runs that argument three ways, from the papers' own tables
(parsed from their arXiv LaTeX sources):

1. **Smartt (2015)'s counting argument.** 30 progenitors (detections and
   limits) lie below log L = 5.1. For a Salpeter IMF, 70% of 8-100 M_sun stars
   are below 18 M_sun, so ~13 more should lie above; at most one does.
2. **The maximum-likelihood mass range** (Smartt et al. 2009's method): fit a
   Salpeter mass function between M_lo and M_hi to the progenitor masses,
   counting the upper limits properly.
3. **Luminosity space** (Davies & Beasor 2018, 2020): with bolometric
   corrections measured on red supergiants near death, fit the progenitors'
   luminosity distribution directly, and ask how unlikely the observed
   brightest progenitor is if red supergiants exploded all the way to the
   Humphreys-Davidson limit (log L = 5.5). Then add the extra 0.2 dex
   scatter that Beasor, Smith & Jencson (2024) found single-band
   luminosities carry.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.special import erf

import config as C

# ---------------------------------------------------------------------------
# Parsing the published tables
# ---------------------------------------------------------------------------


def _table_rows(tex: str, label: str) -> list[list[str]]:
    """The data rows (split on '&') of the table carrying \\label{label}."""
    i = tex.index(r"\label{" + label + "}")
    start = tex.rfind(r"\begin{tabular}", 0, i)
    end = tex.find(r"\end{tabular}", start)
    body = tex[start:end].split("\n")
    rows = []
    for line in body:
        line = line.split("%")[0].strip()
        if re.match(r"SN\s?\d", line) and "&" in line:
            rows.append([c.strip() for c in line.rstrip("\\").split("&")])
    return rows


def _asym(s: str) -> tuple[float, float, float]:
    """'$12^{+3}_{-3}$' -> (12, 3, 3); '12' -> (12, nan, nan)."""
    m = re.search(r"([\d.]+)\$?\^\{\+([\d.]+)\}_\{-([\d.]+)\}", s)
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(3))
    m = re.search(r"([\d.]+)", s)
    return (float(m.group(1)), np.nan, np.nan) if m else (np.nan, np.nan, np.nan)


def load_smartt2015() -> pd.DataFrame:
    """Smartt (2015) Tables 1-3: type II detections, type II upper limits and
    type IIb detections. Masses M(S,G) are from the STARS and rotating Geneva
    models; M(K) from KEPLER."""
    tex = next((C.PAPERS_DIR / "smartt2015").glob("*.tex")).read_text(errors="replace")
    rows = []
    for r in _table_rows(tex, "tab:SNII"):
        m, mp, mm = _asym(r[7])
        rows.append({"sn": r[0], "type": r[1], "kind": "detection", "logT": float(r[2]),
                     "logL": float(r[4]), "logL_err": float(r[5]), "M_SG": m, "M_SG_up": mp,
                     "M_SG_lo": mm, "M_K": _asym(r[8])[0]})
    for r in _table_rows(tex, "tab:II-lims"):
        rows.append({"sn": r[0], "type": "II", "kind": "limit",
                     "logL": float(re.search(r"([\d.]+)", r[1]).group(1)),
                     "M_SG": float(re.search(r"([\d.]+)", r[3]).group(1)),
                     "M_K": float(re.search(r"([\d.]+)", r[4]).group(1))})
    for r in _table_rows(tex, "tab:IIb"):
        rows.append({"sn": r[0], "type": "IIb", "kind": "detection", "logT": float(r[1]),
                     "logL": float(r[3]), "logL_err": float(r[4])})
    return pd.DataFrame(rows)


def load_davies_beasor2018() -> pd.DataFrame:
    """Davies & Beasor (2018) Table 4: the same II-P/II-L progenitors with
    late-type bolometric corrections; log L and initial mass (STARS models)."""
    tex = (C.PAPERS_DIR / "davies_beasor2018" / "prog_MC_tab.tex").read_text()
    rows = []
    for line in tex.splitlines():
        if not line.strip().startswith("SN"):
            continue
        c = [x.strip() for x in line.rstrip("\\ ").split("&")]
        sn, band, Lcol, Mcol = c[0], c[3], c[6], c[7]
        limit = "<" in Lcol
        nums = [float(x) for x in re.findall(r"[\d.]+", Lcol)]
        if limit:
            logL, err = nums[0], nums[1]
            mass = float(re.findall(r"[\d.]+", Mcol)[0])
            rows.append({"sn": sn, "band": band, "kind": "limit", "logL": logL, "logL_err": err,
                         "M_init": mass})
        else:
            logL, err = nums[0], nums[1]
            m, mp, mm = _asym(Mcol)
            rows.append({"sn": sn, "band": band, "kind": "detection", "logL": logL,
                         "logL_err": err, "M_init": m, "M_up": mp, "M_lo": mm})
    df = pd.DataFrame(rows)
    # single-band = everything but multi-band SED fits (SN 2009md) and K-band
    df["single_band_optical"] = ~df["band"].isin(["multi", "K"])
    return df


# Luminosities of recent, well-studied progenitors from full optical-to-mid-IR
# SEDs, as summarised by Beasor, Smith & Jencson (2024, Sect. 3.2.1).
RECENT_SED_PROGENITORS = pd.DataFrame([
    {"sn": "SN2017eaw", "logL_lo": 4.9, "logL_hi": 5.1,
     "source": "Kilpatrick & Foley 2018 (4.9); Van Dyk et al. 2019 (5.1)"},
    {"sn": "SN2023ixf", "logL_lo": 4.74, "logL_hi": 5.1,
     "source": "Kilpatrick et al. 2023 (4.74) to Jencson et al. 2023, Niu et al. 2023 (~5.1)"},
    {"sn": "SN2024ggi", "logL_lo": 4.85, "logL_hi": 4.95, "source": "Xiang et al. 2024 (4.9)"},
])


# ---------------------------------------------------------------------------
# 1. Smartt's counting argument
# ---------------------------------------------------------------------------

def salpeter_number_fraction(m1, m2, lo=C.SMARTT_M_MIN, hi=C.SMARTT_M_MAX,
                             slope=C.SALPETER_SLOPE) -> float:
    """Fraction of stars born between lo and hi that have m1 < M < m2."""
    g = slope + 1.0
    return (m2 ** g - m1 ** g) / (hi ** g - lo ** g)


def counting_argument(n_below: int = 30, n_above_found: int = 1) -> dict:
    f_low = salpeter_number_fraction(C.SMARTT_M_MIN, C.SMARTT_M_SPLIT)
    n_total = n_below / f_low
    expected = n_total - n_below
    from scipy.stats import poisson
    return {"f_8_18_of_8_100": f_low, "n_below": n_below, "n_total": n_total,
            "n_expected_above": expected,
            "p_at_most_found": float(poisson.cdf(n_above_found, expected)),
            "p_none": float(poisson.pmf(0, expected))}


# ---------------------------------------------------------------------------
# 2. Maximum-likelihood mass range (Smartt et al. 2009, Appendix)
# ---------------------------------------------------------------------------

def _imf_pdf(m, lo, hi, slope=C.SALPETER_SLOPE):
    g = slope + 1.0
    norm = (hi ** g - lo ** g) / g
    return np.where((m >= lo) & (m <= hi), m ** slope / norm, 0.0)


def _imf_cdf(m, lo, hi, slope=C.SALPETER_SLOPE):
    g = slope + 1.0
    m = np.clip(m, lo, hi)
    return (m ** g - lo ** g) / (hi ** g - lo ** g)


def mass_range_likelihood(det_m, det_up, det_lo, lim_m, lo_grid, hi_grid) -> np.ndarray:
    """ln L on a (M_lo, M_hi) grid. A detection contributes the IMF convolved
    with its (split-normal) mass error; a limit the probability that a star
    drawn from the IMF lies below it."""
    mm = np.linspace(4.0, 60.0, 2241)
    lnL = np.zeros((len(lo_grid), len(hi_grid)))
    # each detection's error distribution on the mass grid
    kern = []
    for m, up, dn in zip(det_m, det_up, det_lo):
        s = np.where(mm >= m, up, dn)
        k = np.exp(-0.5 * ((mm - m) / s) ** 2)
        kern.append(k / np.trapezoid(k, mm))
    for i, lo in enumerate(lo_grid):
        for j, hi in enumerate(hi_grid):
            if hi <= lo + 0.5:
                lnL[i, j] = -np.inf
                continue
            pdf = _imf_pdf(mm, lo, hi)
            tot = 0.0
            for k in kern:
                tot += np.log(max(np.trapezoid(pdf * k, mm), 1e-300))
            for m in lim_m:
                tot += np.log(max(_imf_cdf(m, lo, hi), 1e-300))
            lnL[i, j] = tot
    return lnL


def profile_interval(grid, lnL_profile, delta=1.92) -> tuple[float, float, float]:
    """Best value and the interval where ln L is within delta of the peak
    (delta = 0.5 for 68%, 1.92 for 95%, for one parameter)."""
    k = int(np.argmax(lnL_profile))
    ok = lnL_profile >= lnL_profile[k] - delta
    return float(grid[k]), float(grid[ok].min()), float(grid[ok].max())


def smartt_likelihood(det_m, det_up, det_lo, lim_m, lo_grid, hi_grid,
                      slope=C.SALPETER_SLOPE) -> np.ndarray:
    """ln L on a (M_lo, M_hi) grid with Smartt et al. (2009)'s own likelihood
    (their Sect. 7, kept in Smartt 2015). A detection at m_j (range m_j - lo
    to m_j + up) contributes the IMF probability over its error range: below
    m_j a ramp from zero at the lower limit times the IMF value at m_j, above
    m_j the IMF itself, truncated at [M_lo, M_hi]. A non-detection with mass
    limit m_i contributes the IMF fraction below it; 1 if M_hi < m_i, and 0.16
    if M_lo > m_i (the limits are 84% confidence limits). The IMF is a
    normalised Salpeter power law (their printed normalisation lacks the
    factor -Gamma = 1.35, which would penalise every M_hi above a limit by
    that factor)."""
    G = slope + 1.0                      # Gamma = -1.35
    LO, HI = np.meshgrid(lo_grid, hi_grid, indexing="ij")
    valid = HI > LO + 0.5
    D = (LO ** G - HI ** G) / (-G)       # integral of m^(G-1) over [LO, HI]
    lnL = np.zeros_like(LO)
    for m, up, dn in zip(det_m, det_up, det_lo):
        a, b, c = m - dn, m, m + up
        x1, x2 = np.maximum(a, LO), np.minimum(b, HI)
        ramp = np.where(x2 > x1, b ** (G - 1) / (b - a) * ((x2 - a) ** 2 - (x1 - a) ** 2) / 2, 0.0)
        y1, y2 = np.maximum(b, LO), np.minimum(c, HI)
        upper = np.where(y2 > y1, (y2 ** G - y1 ** G) / G, 0.0)
        lnL += np.log(np.maximum((ramp + upper) / D, 1e-300))
    for m in lim_m:
        cdf = (np.minimum(m, HI) ** G - LO ** G) / G / D
        p = np.where(HI < m, 1.0, np.where(LO > m, 0.16, cdf))
        lnL += np.log(np.maximum(p, 1e-300))
    return np.where(valid, lnL, -np.inf)


def fit_mass_range(det: pd.DataFrame, lim: pd.DataFrame, mcol="M_SG", upcol="M_SG_up",
                   locol="M_SG_lo", limcol="M_SG") -> dict:
    """Smartt et al. (2009)'s two steps: M_lo (and a first M_hi) from the
    detections alone, then M_hi from detections and limits with M_lo fixed,
    because the limits say little about the low end. Intervals from the
    one-parameter profile (Delta chi2 = 1 for 68%, 3.84 for 95%)."""
    lo_grid = np.arange(5.0, 12.01, 0.05)
    hi_grid = np.arange(10.0, 40.01, 0.05)
    dm, du, dl = det[mcol].to_numpy(), det[upcol].to_numpy(), det[locol].to_numpy()
    L1 = smartt_likelihood(dm, du, dl, [], lo_grid, hi_grid)
    i, j = np.unravel_index(np.argmax(L1), L1.shape)
    m_lo = float(lo_grid[i])
    lo_prof = 2 * (L1.max() - L1.max(axis=1))
    L2 = smartt_likelihood(dm, du, dl, lim[limcol].to_numpy(), np.array([m_lo]), hi_grid)[0]
    d2 = 2 * (L2.max() - L2)

    def rng_(grid, d, level):
        ok = grid[d <= level]
        return float(ok.min()), float(ok.max())

    out = {"M_lo": m_lo, "M_lo_95": rng_(lo_grid, lo_prof, 3.84),
           "M_hi_detections_only": float(hi_grid[j]),
           "M_hi": float(hi_grid[np.argmax(L2)]), "M_hi_68": rng_(hi_grid, d2, 1.0),
           "M_hi_95": rng_(hi_grid, d2, 3.84), "hi_grid": hi_grid, "delta_chi2_hi": d2}
    for m in (18.0, 25.0, 30.0):
        out[f"delta_chi2_Mhi_{int(m)}"] = float(d2[int(np.argmin(np.abs(hi_grid - m)))])
    return out


# ---------------------------------------------------------------------------
# 3. Luminosity space
# ---------------------------------------------------------------------------

def _pl_pdf_logL(x, xlo, xhi, gamma):
    """pdf in x = log L for dN/dL ~ L^gamma between 10^xlo and 10^xhi."""
    a = gamma + 1.0
    L, Llo, Lhi = 10 ** x, 10 ** xlo, 10 ** xhi
    norm = (Lhi ** a - Llo ** a) / a if abs(a) > 1e-9 else np.log(Lhi / Llo)
    return np.where((x >= xlo) & (x <= xhi), L ** a * np.log(10) / norm, 0.0)


def luminosity_likelihood(det_x, det_err, lim_x, lim_err, xlo, xhi, gamma,
                          extra_err: np.ndarray | None = None) -> float:
    """ln L for progenitor luminosities: detections convolved with a Gaussian
    error in log L, limits as the probability of lying below the limit
    (itself uncertain by lim_err)."""
    xx = np.linspace(3.5, 6.5, 1501)
    pdf = _pl_pdf_logL(xx, xlo, xhi, gamma)
    tot = 0.0
    e_det = np.asarray(det_err, float)
    if extra_err is not None:
        e_det = np.hypot(e_det, extra_err[: len(e_det)])
    for x, e in zip(det_x, e_det):
        k = np.exp(-0.5 * ((xx - x) / e) ** 2) / (e * np.sqrt(2 * np.pi))
        tot += np.log(max(np.trapezoid(pdf * k, xx), 1e-300))
    for x, e in zip(lim_x, lim_err):
        # P(true luminosity below the limit), the limit itself uncertain by e:
        # integral of pdf(L) * Phi((limit - L) / e)
        w = 0.5 * (1 + erf((x - xx) / (np.sqrt(2) * max(e, 1e-3))))
        tot += np.log(max(np.trapezoid(pdf * w, xx), 1e-300))
    return tot


def fit_luminosity_range(db: pd.DataFrame, gamma: float, extra_single_band: float = 0.0) -> dict:
    det = db[db["kind"] == "detection"]
    lim = db[db["kind"] == "limit"]
    extra = np.where(det["single_band_optical"], extra_single_band, 0.0)
    lim_err = np.hypot(lim["logL_err"].to_numpy(),
                       np.where(lim["single_band_optical"], extra_single_band, 0.0))
    xlo_grid = np.arange(3.8, 4.71, 0.02)
    xhi_grid = np.arange(4.9, 6.21, 0.01)
    lnL = np.full((len(xlo_grid), len(xhi_grid)), -np.inf)
    for i, a in enumerate(xlo_grid):
        for j, b in enumerate(xhi_grid):
            lnL[i, j] = luminosity_likelihood(det["logL"], det["logL_err"], lim["logL"], lim_err,
                                              a, b, gamma, extra_err=extra)
    prof = lnL.max(axis=0)
    best = profile_interval(xhi_grid, prof, 0.5)
    k55 = int(np.argmin(np.abs(xhi_grid - C.LOGL_HD_LIMIT)))
    ib, jb = np.unravel_index(np.argmax(lnL), lnL.shape)
    xlo_best, xlo_at_55 = float(xlo_grid[ib]), float(xlo_grid[int(np.argmax(lnL[:, k55]))])
    d2 = float(2 * (prof.max() - prof[k55]))
    # one-sided significance of excluding L_hi = 5.5 from the likelihood ratio
    from scipy.stats import norm, chi2
    p = float(chi2.sf(d2, 1) / 2) if xhi_grid[np.argmax(prof)] < C.LOGL_HD_LIMIT else 0.5
    return {"gamma": gamma, "extra_err": extra_single_band, "logL_hi": best[0],
            "logL_lo_best": xlo_best, "logL_lo_at_5p5": xlo_at_55,
            "logL_hi_68": best[1:], "logL_hi_95": profile_interval(xhi_grid, prof, 1.92)[1:],
            "delta2lnL_at_5p5": d2, "p_one_sided": p,
            "sigma": float(norm.isf(p)) if p < 0.5 else 0.0,
            "xhi_grid": xhi_grid, "profile": prof}


def brightest_of_n(db: pd.DataFrame, gamma: float, xhi: float = C.LOGL_HD_LIMIT,
                   xlo: float = 4.3, extra_single_band: float = 0.0, n_mc: int = C.N_MC_PROGENITORS,
                   rng=None) -> dict:
    """Davies & Beasor (2020)'s order-statistics test: draw as many progenitors
    as were detected from L^gamma between 10^xlo and 10^xhi, blur each by the
    measurement error, and count how often the brightest is no brighter than
    the brightest observed."""
    rng = rng or np.random.default_rng(C.RNG_SEED)
    det = db[db["kind"] == "detection"]
    n = len(det)
    err = np.hypot(det["logL_err"].to_numpy(),
                   np.where(det["single_band_optical"], extra_single_band, 0.0))
    a = gamma + 1.0
    u = rng.random((n_mc, n))
    Llo, Lhi = 10 ** xlo, 10 ** xhi
    x = np.log10((Llo ** a + u * (Lhi ** a - Llo ** a)) ** (1 / a))
    x += rng.normal(size=x.shape) * err[None, :]
    observed_max = float(det["logL"].max())
    frac = float(np.mean(x.max(axis=1) <= observed_max))
    from scipy.stats import norm
    return {"n_detections": n, "observed_max": observed_max, "p": frac,
            "sigma": float(norm.isf(frac)), "extra_err": extra_single_band,
            "simulated_max_median": float(np.median(x.max(axis=1)))}


def nothing_brighter(db: pd.DataFrame, gamma: float, xhi: float = C.LOGL_HD_LIMIT,
                     xlo: float = 4.3, extra_single_band: float = 0.0,
                     n_mc: int = C.N_MC_PROGENITORS, rng=None) -> dict:
    """The red supergiant problem as one probability, using every supernova.

    Draw a true luminosity for each of the N supernovae in the sample (14
    detections, 10 limits) from L^gamma between 10^xlo and 10^xhi, blur it by
    that supernova's own measurement error, and ask whether the sample would
    have shown nothing brighter than we actually saw: every detection at most
    as bright as the brightest real detection (SN 2009hd, log L = 5.24), and
    every limit respected (two limits, SN 1999an and SN 2002hh, are looser)."""
    rng = rng or np.random.default_rng(C.RNG_SEED)
    n = len(db)
    err = np.hypot(db["logL_err"].to_numpy(),
                   np.where(db["single_band_optical"], extra_single_band, 0.0))
    top = float(db.loc[db["kind"] == "detection", "logL"].max())
    ceiling = np.where(db["kind"] == "limit", np.maximum(db["logL"], top), top)
    a = gamma + 1.0
    Llo, Lhi = 10 ** xlo, 10 ** xhi
    x = np.log10((Llo ** a + rng.random((n_mc, n)) * (Lhi ** a - Llo ** a)) ** (1 / a))
    x_obs = x + rng.normal(size=x.shape) * err[None, :]
    ok = np.all(x_obs <= ceiling[None, :], axis=1)
    from scipy.stats import norm
    p = float(ok.mean())
    return {"n": n, "brightest_seen": top, "p": p, "sigma": float(norm.isf(p)),
            "extra_err": extra_single_band, "simulated_max": x_obs.max(axis=1)}


def salpeter_gamma_L(alpha_mlr: float = 2.0, slope: float = C.SALPETER_SLOPE) -> float:
    """dN/dL ~ L^Gamma_L for dN/dM ~ M^slope and L ~ M^alpha (Davies & Beasor 2020)."""
    return (1 - alpha_mlr + slope) / alpha_mlr


def mass_from_logL(x) -> np.ndarray:
    """Davies & Beasor (2018)'s STARS mass-luminosity relation (neon burning)."""
    a, b = C.DB18_MLR
    return 10 ** (a + b * np.asarray(x, float))
