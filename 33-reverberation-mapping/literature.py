"""Published values we compare with, and parsers for the published tables we re-fit.

Nothing here feeds a measurement: these are the answers we grade ourselves
against, plus three tables read straight from the papers' arXiv TeX sources so
that their fits can be redone with our code:

* Bentz et al. 2013 (ApJ 767, 149), Table 14: 71 H-beta lags and host-corrected
  luminosities, the canonical R-L relation.
* Woo et al. 2015 (ApJ 801, 38), Table A1: virial products and stellar velocity
  dispersions of 29 reverberation-mapped AGN (the M-sigma calibration of f).
* Shen et al. 2024 (ApJS 272, 26), Table 3: virial products and BLR
  dynamical-model masses of 30 AGN (the dynamical calibration of f).
"""

from __future__ import annotations

import re
import tarfile
from functools import lru_cache

import numpy as np
import pandas as pd

import config as C

ARXIV = C.DATA_DIR / "arxiv"

# ---------------------------------------------------------------------------
# NGC 5548 during AGN STORM (2014)
# ---------------------------------------------------------------------------

#: Pei et al. 2017, Table 5: rest-frame lags (days) as (value, minus, plus).
PEI17_LAGS = {
    ("Hb", "F1158", "cent"): (6.23, 0.44, 0.39),
    ("Hb", "F1158", "peak"): (6.14, 0.98, 0.74),
    ("Hb", "F1158", "cent_T1"): (7.62, 0.49, 0.49),
    ("Hb", "F1158", "cent_T2"): (5.99, 0.75, 0.71),
    ("Hb", "F1367", "cent"): (5.89, 0.37, 0.37),
    ("Hb", "F5100", "cent"): (4.17, 0.36, 0.36),
    ("Hb", "F5100", "peak"): (4.42, 0.25, 0.98),
    ("Hb", "F5100", "cent_T1"): (4.99, 0.47, 0.40),
    ("Hb", "F5100", "cent_T2"): (3.10, 0.80, 0.77),
    ("Hb", "F5100", "jav"): (3.84, 0.59, 0.57),
    ("Hb", "V", "cent"): (3.79, 0.34, 0.37),
    ("HeII", "F1158", "cent"): (2.69, 0.25, 0.24),
    ("HeII", "F5100", "cent"): (0.79, 0.34, 0.35),
    ("F5100", "F1158", "cent"): (2.23, 0.26, 0.31),
}

#: Pei et al. 2017, Tables 6-7: H-beta rms line dispersion, virial product and
#: mass for the full campaign (f = 4.47, Woo et al. 2015).
PEI17_SIGMA_RMS = (4278.0, 671.0)  # km/s
PEI17_VP = (1.49e7, 0.49e7)  # Msun
PEI17_MBH = (6.66e7, 2.17e7)  # Msun
PEI17_F_AGN = (7.44, 0.50)  # 1e-15 erg/s/cm^2/A, host subtracted
PEI17_HB_ETA = (0.15, 0.01)  # responsivity before the anomaly
PEI17_HB_FLOST = 0.06  # mean flux deficit during the anomaly

#: Goad et al. 2016 (AGN STORM IV): UV lines vs the 1157 A continuum.
#: Table 1: CCF centroid lags over the first 75 days (days).
GOAD16_LAGS = {"Lya": (6.69, 0.41), "SiIV": (5.80, 0.58), "CIV": (4.97, 0.39), "HeII": (2.42, 0.44)}
#: Table 2: responsivity eta (log F_line vs log F_cont slope), blue + green points.
GOAD16_ETA = {"Lya": (0.30, 0.01), "SiIV": (0.45, 0.01), "CIV": (0.25, 0.01), "HeII": (0.58, 0.04)}
#: Section 3.4: time-averaged fraction of line flux lost during the anomaly.
GOAD16_FLOST = {"Lya": 0.09, "SiIV": 0.23, "CIV": 0.18, "HeII": 0.21}
#: Figure 1 caption: the colour-coded segments, HJD - 2,400,000, of the C IV
#: light curve shifted back by its lag: normal | onset | anomaly | recovery | normal.
GOAD16_SEGMENTS = {
    "blue": (56691.2813, 56765.1133), "cyan": (56766.1094, 56776.3984),
    "red": (56777.4336, 56813.9688), "magenta": (56814.7891, 56829.8281),
    "green": (56830.8242, 56856.8438),
}

#: Pancoast et al. 2014: BLR dynamical model of NGC 5548 (2008 LAMP data).
PANCOAST14_MBH = (3.24e7, 0.90e7, 2.26e7)  # value, minus, plus

# ---------------------------------------------------------------------------
# R-L relation and masses
# ---------------------------------------------------------------------------

#: Bentz et al. 2013 (LINMIX_ERR, all 71 measurements): log(R/lt-day) =
#: K + alpha log(L/1e44), with variance sigma^2 of the scatter.
BENTZ13_RL = {"K": (1.527, 0.031, 0.031), "alpha": (0.533, 0.033, 0.035), "var": (0.036, 0.008, 0.009),
              "scatter_dex": 0.19, "N": 71}

#: Shen et al. 2024, Table 4: SDSS-RM R-L regressions (Kelly 2007), L in units
#: of 1e44 (5100) or 1e45 (3000, 1350) erg/s: (a, b, sigma_int, N).
SHEN24_RL = {
    "Hb": dict(a=(1.458, 0.038), b=(0.41, 0.07), sig=(0.32, 0.03), N=80, Lnorm=44.0,
               a_fixed=(1.445, 0.038), b_fixed=0.533),
    "MgII": dict(a=(2.086, 0.030), b=(0.31, 0.06), sig=(0.32, 0.02), N=124, Lnorm=45.0,
                 a_fixed=(2.095, 0.030), b_fixed=0.39),
    "CIV": dict(a=(1.840, 0.074), b=(0.32, 0.11), sig=(0.51, 0.04), N=89, Lnorm=45.0,
                a_fixed=(1.783, 0.055), b_fixed=0.45),
}
#: Detections reported in Shen et al. 2024, Table 2 (PyROA, > 2 sigma).
SHEN24_NDET = {"HA": 23, "HB": 81, "MG2": 125, "C4": 110}
#: Shen et al. 2024, section 4.5: mean and intrinsic scatter of log f(sigma_rms).
SHEN24_LOGF = {"mean": (0.62, 0.07), "scatter": (0.31, 0.07)}

#: Woo et al. 2015, Appendix: joint M-sigma + virial-factor fit (sigma_line):
#: log M = alpha + beta log(sigma/200 km/s), log f, intrinsic scatter.
WOO15_FIT = {"alpha": (8.34, 0.05), "beta": (4.97, 0.28), "logf": (0.65, 0.12), "scatter": (0.43, 0.03)}

#: NGC 5548's own R-L slope (H-beta lag vs host-corrected L5100), Bentz et al.
#: 2007 (ApJ 662, 205, section 4.2): 0.73 +/- 0.14 from AGN Watch years 1-13
#: plus 2005; 0.66 +/- 0.13 without the poorly sampled Year 12 (their
#: preferred value). With F_opt ~ F_UV^(0.84 +/- 0.05) this becomes
#: r ~ L_UV^(0.55 +/- 0.14). Their virial slope, sigma_line vs tau: -0.5 +/- 0.1.
BENTZ07_NGC5548_SLOPE = (0.73, 0.14)
BENTZ07_NGC5548_SLOPE_NO_Y12 = (0.66, 0.13)
BENTZ07_OPT_UV = (0.84, 0.05)
BENTZ07_VIRIAL_SLOPE = (-0.5, 0.1)


# ---------------------------------------------------------------------------
# Table parsers (arXiv TeX)
# ---------------------------------------------------------------------------


def _tex(label: str, arxiv_id: str, member_suffix: str) -> str:
    path = ARXIV / f"{label}_{arxiv_id}.tar.gz"
    with tarfile.open(path) as tar:
        for m in tar.getmembers():
            if m.name.endswith(member_suffix):
                return tar.extractfile(m).read().decode("latin-1")
    raise FileNotFoundError(f"{member_suffix} not in {path}")


def _strip_comments(tex: str) -> list[str]:
    return [ln for ln in tex.splitlines() if not ln.lstrip().startswith("%")]


_PM = r"\$?\s*([-\d.]+)\s*\^\{\s*\+?\s*([-\d.]+)\s*\}_\{\s*-?\s*([-\d.]+)\s*\}\s*\$?"
_PMS = r"\$?\s*([-\d.]+)\s*\\pm\s*([-\d.]+)\s*\$?"


@lru_cache(maxsize=None)
def bentz13_table() -> pd.DataFrame:
    """Bentz et al. 2013 Table 14: object, lag (+, -), f_AGN, log lambda L_lambda (5100)."""
    tex = _tex("bentz2013", "1303.1742", "ms.tex")
    start = tex.index(r"\tablecaption{Rest-frame H$\beta$ Time Lags and 5100")
    block = tex[start:tex.index(r"\enddata", start)]
    rows, name = [], None
    for ln in _strip_comments(block):
        if "&" not in ln or r"\colhead" in ln:
            continue
        cells = [c.strip() for c in ln.replace(r"\\", "").split("&")]
        if len(cells) < 4:
            continue
        if cells[0]:
            name = cells[0].replace(r"\,", " ").strip()
        lag = re.search(_PM, cells[1])
        fa = re.search(_PMS, cells[2])
        ll = re.search(_PMS, cells[3])
        if not (lag and fa and ll):
            continue
        rows.append(dict(object=name, tau=float(lag.group(1)), tau_hi=float(lag.group(2)),
                         tau_lo=float(lag.group(3)), f_agn=float(fa.group(1)), f_agn_err=float(fa.group(2)),
                         logL=float(ll.group(1)), logL_err=float(ll.group(2))))
    return pd.DataFrame(rows)


@lru_cache(maxsize=None)
def woo15_table() -> pd.DataFrame:
    """Woo et al. 2015 Table A1, one row per AGN: mean virial product (sigma_line
    based, 1e6 Msun, +/-) and stellar velocity dispersion (km/s)."""
    tex = _tex("woo2015", "1412.7225", "msr_astroph.tex")
    start = tex.index(r"\tablecaption{Virial Products and \sigmastar\ of the Reverberation-mapped AGNs}")
    block = tex[start:tex.index(r"\enddata", start)]
    block = block[block.index(r"\startdata") + len(r"\startdata"):]
    out, name = [], None
    vp_pat = r"\$\s*([\d.]+)\s*\^\{\s*\+?\s*([\d.]+)\s*\}_\{\s*-\s*([\d.]+)\s*\}\s*\$"
    for ln in _strip_comments(block):
        if "&" not in ln:
            continue
        cells = [c.strip() for c in ln.replace(r"\\", "").split("&")]
        if cells[0]:
            name = cells[0]
        sig = re.search(r"(\d+)\s*\$\\pm\$\s*(\d+)", cells[9]) if len(cells) > 9 else None
        vp = re.search(vp_pat, cells[7]) if len(cells) > 7 else None
        if sig and vp:
            out.append(dict(object=name, vp=float(vp.group(1)) * 1e6, vp_hi=float(vp.group(2)) * 1e6,
                            vp_lo=float(vp.group(3)) * 1e6, sigma_star=float(sig.group(1)),
                            sigma_star_err=float(sig.group(2)), mean_of_several=(cells[1] == "mean")))
    return pd.DataFrame(out)


@lru_cache(maxsize=None)
def shen24_f_table() -> pd.DataFrame:
    """Shen et al. 2024 Table 3: log VP, log M_dyn and log f with errors (30 rows)."""
    tex = _tex("shen2024", "2305.01014", "final_results_rv1_astroph.tex")
    start = tex.index(r"\label{tab:f_sample}")
    block = tex[start:tex.index(r"\end{tabular}", start)]
    out = []
    for ln in _strip_comments(block):
        cells = [c.strip() for c in ln.replace(r"\\", "").replace(r"\hline", "").split("&")]
        if len(cells) != 9 or not re.match(r"[-\d.$]", cells[7]):
            continue
        try:
            logf = float(cells[7].replace("$", ""))
            slogf = float(cells[8])
        except ValueError:
            continue
        num = lambda s: float(s) if re.match(r"^-?[\d.]+$", s) else np.nan  # noqa: E731
        out.append(dict(object=re.sub(r"\$\^.\$", "", cells[0]).replace("--", "-"), logVP=num(cells[1]),
                        logVP_err=num(cells[2]), logMdyn=num(cells[4]), logMdyn_err=num(cells[5]),
                        logf=logf, logf_err=slogf, ref_dyn=cells[6]))
    return pd.DataFrame(out)
