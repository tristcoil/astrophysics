"""Published tests of gravity, read from the checked tables in ``literature/checks/``.

Five tables were transcribed for this project and every number in them was
checked against the primary paper's own text (abstract, table, TeX source or page
scan), recording the exact sentence, the confidence convention and corrections
to our first drafts: ``wep.csv`` (equivalence principle), ``gamma.csv`` (light
bending and time delay), ``beta_lunar.csv`` (perihelion, PPN limits, lunar
ranging), ``pulsars.csv`` and ``strongfield.csv``. Each has a ``*_notes.md``.

Published results come with different conventions: 1 sigma, 68%, 90%, 95%,
probable errors, separate statistical and systematic errors, one-sided limits.
``one_sigma`` turns each into a symmetric 1-sigma equivalent (the Watch-out of
the atlas spec): 90% / 1.645, 95% / 1.96, probable error x 1.4826, statistical
and systematic errors in quadrature, asymmetric errors averaged. The lunar-ranging
"realistic" errors (several times the formal error, as their authors recommend)
are used as 1-sigma.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

import config as C

CHECKS = C.LITERATURE_DIR / "checks"

#: divide a published error of this type by the factor to get 1 sigma
TO_ONE_SIGMA = {"1sigma": 1.0, "68%": 1.0, "standard_error": 1.0, "formal": 1.0, "realistic": 1.0,
                "not_stated": 1.0, "3sigma": 1.0,      # INPOP17a's "3-sigma" = 3 x formal, a realistic 1-sigma
                "2sigma": 2.0, "90%": 1.645, "95%": 1.96, "probable_error": 1 / 1.4826,
                "exact": 1.0}


@lru_cache(maxsize=None)
def table(name: str) -> pd.DataFrame:
    df = pd.read_csv(CHECKS / f"{name}.csv")
    if df["id"].duplicated().any():
        raise ValueError(f"duplicate ids in {name}.csv")
    return df.set_index("id")


def row(name: str, rid: str) -> pd.Series:
    return table(name).loc[rid]


def num(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def one_sigma(name: str, rid: str, *, treat_as: str | None = None) -> tuple[float, float]:
    """(value, symmetric 1-sigma error) of a checked row, conventions converted."""
    r = row(name, rid)
    v = num(r["value"])
    ep, em = num(r["err_plus"]), num(r["err_minus"])
    if np.isnan(em):
        em = ep
    err = np.nanmean([ep, em]) if not (np.isnan(ep) and np.isnan(em)) else np.nan
    et = treat_as or str(r["err_type"]).strip()
    if et.startswith("upper_limit") or et.startswith("range"):
        raise ValueError(f"{name}:{rid} is a limit ({et}); use limit()")
    err = err / TO_ONE_SIGMA[et]
    if "err_syst" in r.index and not np.isnan(num(r["err_syst"])):
        err = float(np.hypot(err, num(r["err_syst"]) / TO_ONE_SIGMA[et]))
    return v, float(err)


def limit(name: str, rid: str) -> tuple[float, float]:
    """(published bound, 1-sigma-equivalent width) of a one-sided limit."""
    r = row(name, rid)
    v = num(r["value"])
    et = str(r["err_type"]).strip()
    fac = {"upper_limit_95": 1.96, "upper_limit_90": 1.645, "upper_limit_68": 1.0}.get(et, 1.0)
    return v, v / fac


# ---------------------------------------------------------------------------
# The equivalence-principle timeline (wep.csv; plotting advice in wep_notes.md)
# ---------------------------------------------------------------------------

G_STD = 980.665e6        # microGal per g

#: (id, series, label, how the plotted 1-sigma precision is formed)
WEP_ENTRIES = [
    ("newton1687", "lab", "Newton (pendulums)", "limit"),
    ("bessel1832", "lab", "Bessel (pendulums)", "limit"),
    ("eotvos1890", "lab", "Eötvös 1890", "limit"),
    ("epf1922", "lab", "Eötvös, Pekár & Fekete", "limit"),
    ("renner1935", "lab_doubtful", "Renner (errors too small)", "limit"),
    ("rkd1964", "lab", "Roll, Krotkov & Dicke", "1sigma_override"),
    ("bp1972", "lab", "Braginsky & Panov", "sigma"),
    ("niebauer1987", "lab", "Niebauer et al. (free fall)", "microgal"),
    ("adelberger1990_CuBe", "lab", "Eöt-Wash 1990", "sigma"),
    ("su1994_earth_combined", "lab", "Eöt-Wash 1994", "sigma"),
    ("schlamminger2008", "lab", "Eöt-Wash 2008", "sigma"),
    ("wagner2012_BeAl_earth", "lab", "Eöt-Wash 2012", "sigma"),
    ("microscope2017", "space", "MICROSCOPE 2017", "sigma"),
    ("microscope2022", "space", "MICROSCOPE 2022", "sigma"),
    ("peters1999", "atom", "Peters, Chung & Chu", "limit"),
    ("fray2004", "atom", "Fray et al.", "sigma"),
    ("merlet2010", "atom", "Merlet et al.", "microgal"),
    ("schlippert2014", "atom", "Schlippert et al.", "sigma"),
    ("zhou2015", "atom", "Zhou et al. 2015", "sigma"),
    ("asenbaum2020", "atom", "Asenbaum et al.", "sigma"),
    ("llr_williams2004", "llr", "LLR, Williams et al. 2004", "sigma"),
    ("llr_williams2012", "llr", "LLR, Williams et al. 2012", "sigma"),
    ("llr_hofmann2018", "llr", "LLR, Hofmann & Müller", "sigma"),
    ("llr_viswanathan2018", "llr", "LLR, INPOP17a", "sigma"),
    ("llr_biskupek2021", "llr", "LLR, Biskupek et al.", "sigma"),
]


def wep_timeline() -> pd.DataFrame:
    rows = []
    for rid, series, label, how in WEP_ENTRIES:
        r = row("wep", rid)
        if how == "limit":
            val, prec = np.nan, num(r["value"])
        elif how == "microgal":       # differential acceleration in microGal -> eta
            v, e = one_sigma("wep", rid)
            val, prec = v / G_STD, e / G_STD
        elif how == "1sigma_override":   # RKD's "probable error" behaves as 1 sigma (wep_notes.md)
            val, prec = one_sigma("wep", rid, treat_as="1sigma")
        else:
            val, prec = one_sigma("wep", rid)
        rows.append({"id": rid, "year": int(r["year"]), "series": series, "label": label,
                     "value": val, "precision": prec, "kind": how})
    # LLR 1976 from eta_N = -0.001 +- 0.015 times the Earth-Moon self-energy difference 4.45e-10
    v, e = one_sigma("wep", "llr_shapiro1976")
    rows.append({"id": "llr_shapiro1976", "year": 1976, "series": "llr", "label": "LLR, Shapiro et al. 1976",
                 "value": v * 4.45e-10, "precision": e * 4.45e-10, "kind": "eta_N x 4.45e-10"})
    return pd.DataFrame(rows).sort_values("year").reset_index(drop=True)


def antimatter() -> dict:
    v, e = one_sigma("wep", "alphag2023")
    return {"a_over_g": v, "sigma": e, "ref": row("wep", "alphag2023")["reference"]}


# ---------------------------------------------------------------------------
# beta, gamma and the Nordtvedt parameter (beta_lunar.csv, gamma.csv)
# ---------------------------------------------------------------------------

def ppn_constraints() -> dict:
    """The constraints drawn in the gamma-beta plane (1 sigma)."""
    g, sg = one_sigma("beta_lunar", "PPN-B03")                 # Cassini: gamma = 1 + (2.1 +- 2.3)e-5
    b17, sb17 = one_sigma("beta_lunar", "BETA-P17")            # MESSENGER perihelion, gamma at Cassini value
    eta, seta = one_sigma("beta_lunar", "NORD-HM18b")          # LLR eta_N = 4 beta - gamma - 3 (WEP assumed)
    eta_c, seta_c = one_sigma("beta_lunar", "NORD-HOF17")      # same, WEP part subtracted ("clean" SEP)
    pp_b, pp_sb = one_sigma("beta_lunar", "BETA-PP13")         # EPM2011 joint fit
    pp_g, pp_sg = 4e-5, 6e-5                                   # its gamma - 1 (same paper, beta_lunar notes)
    mer72, smer72 = one_sigma("beta_lunar", "MER-SH72")        # (2 + 2 gamma - beta)/3 from radar (1972)
    return {
        "cassini_gamma_minus_1": (g, sg),
        # perihelion: the published beta - 1 assumes Cassini's gamma, so the constraint on
        # 2(gamma - 1) - (beta - 1) is 2 * g - b17 with b17's error (gamma's is negligible)
        "perihelion_2g_minus_b": (2 * g - b17, float(np.hypot(sb17, 2 * sg))),
        "nordtvedt_eta": (eta, seta),
        "nordtvedt_eta_clean": (eta_c, seta_c),
        "epm2011": ((pp_g, pp_sg), (pp_b, pp_sb)),
        "mercury_radar_1972": (mer72 - 1, smer72),
    }


@lru_cache(maxsize=1)
def gamma_table() -> pd.DataFrame | None:
    path = CHECKS / "gamma.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


GR_LIMB_ARCSEC = 1.7512          # gr_predictions.light_deflection_limb()

#: one point per independent experiment (gamma_notes.md, section 6); flag: "" | "rereduction" | "secondary" | "flagged"
GAMMA_ENTRIES = [
    ("E1919-sobral-4inch", ""), ("E1919-principe", ""), ("E1919-sobral-combined-harvey1979", "rereduction"),
    ("E1922-wallal-lick-15ft", ""), ("E1922-wallal-lick-5ft", ""), ("E1922-cordillo-downs", ""),
    ("E1929-takengon", ""), ("E1936-ussr-mikhailov-1959", ""), ("E1947-bocaiuva", ""), ("E1952-khartoum", ""),
    ("E1973-chinguetti", ""), ("E2017-wyoming-bruns", ""),
    ("R1969-owens-valley", ""), ("R1969-goldstone", ""), ("R1970-goldstone-haystack-vlbi", "secondary"),
    ("R1970-nrao-sramek", ""), ("R1970-cambridge-hill", ""), ("R1971-nrao-sramek", "secondary"),
    ("R1972-cambridge-riley", ""), ("R1972-westerbork", ""), ("R1972-vlbi-counselman", ""),
    ("R1973-westerbork-dual", ""), ("R1974-nrao-fomalont", ""), ("R1975-nrao-fomalont", ""),
    ("V1984-robertson-carter", ""), ("V1987-lebach", ""), ("V1991-robertson-carter-dillinger", ""),
    ("V2004-shapiro", ""), ("V2005-vlba-fomalont", ""), ("V2011-lambert", ""), ("V2018-titov-global", "flagged"),
    ("A1997-hipparcos", ""),
    ("T1967-shapiro-radar", "secondary"), ("T1970-shapiro-radar", ""), ("T1970-mariner6-7", ""),
    ("T1972-mariner9", "secondary"), ("T1976-viking", ""), ("T1985-psr1937", "secondary"), ("T1985-voyager2", ""),
    ("T2002-cassini", ""),
]


def gamma_timeline() -> pd.DataFrame:
    """(1+gamma)/2 with a 1-sigma error for each experiment (probable errors x 1.483)."""
    g = gamma_table().set_index("id")
    rows = []
    for rid, flag in GAMMA_ENTRIES:
        r = g.loc[rid]
        q = str(r["quantity_as_published"])
        v = num(r["value"])
        ep, em = num(r["err_plus"]), num(r["err_minus"])
        e = np.nanmean([ep, em])
        e = e / TO_ONE_SIGMA.get(str(r["err_type"]).strip(), 1.0)
        if q == "deflection_limb_arcsec":
            val, sig = v / GR_LIMB_ARCSEC, e / GR_LIMB_ARCSEC
        elif q in ("deflection_ratio_to_GR", "(1+gamma)/2"):
            val, sig = v, e
        elif q == "gamma":
            val, sig = (1 + v) / 2, e / 2
        elif q == "gamma-1":
            val, sig = 1 + v / 2, e / 2
        else:
            raise ValueError(q)
        method = str(r["method"])
        rows.append({"id": rid, "year": int(r["year"]), "method": method, "value": val, "sigma": sig, "flag": flag,
                     "err_type": str(r["err_type"]), "experiment": str(r["experiment"])})
    return pd.DataFrame(rows).sort_values("year").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Binary pulsars (pulsars.csv)
# ---------------------------------------------------------------------------

def pulsar_value(rid: str) -> tuple[float, float]:
    r = row("pulsars", rid)
    v = num(str(r["value"]).replace(" ", ""))
    e = num(str(r["err_plus"]).replace(" ", ""))
    return v, e
