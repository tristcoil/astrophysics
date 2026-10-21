"""Readers for the published catalogues we compare with.

Each reader returns a pandas DataFrame indexed by KIC number, with the
catalogues' sentinel values (-9999, -999.99, 0 for missing parallaxes) turned
into NaN and the columns renamed to plain identifiers.

* Yu et al. (2018, ApJS 236, 42): numax, Dnu, mode amplitude, granulation
  power (table 1) and masses/radii with and without a Dnu correction (table 2)
  for 16,094 Kepler red giants, measured with the SYD pipeline.
* APOKASC-3 (Pinsonneault et al. 2025, ApJS 276, 69): recommended (table 4)
  and alternative (table 5) properties of 15,808 evolved stars, and every
  pipeline's numax and Dnu (table 9).
* Vrard et al. (2016, A&A 588, A87): gravity-mode period spacings Delta Pi_1
  and evolutionary states of 6,111 red giants.
* Lund et al. (2017, ApJ 835, 172), the Kepler LEGACY sample: global and
  individual-mode parameters of 66 main-sequence stars incl. 16 Cyg A and B.
"""

from __future__ import annotations

import warnings
from functools import lru_cache

import numpy as np
import pandas as pd

import config as C


_ID_COLUMNS = {"KIC", "GaiaDR3"}


def _read_cds(folder: str, table: str) -> pd.DataFrame:
    from astropy.io import ascii

    d = C.CATALOG_DIR / folder
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tab = ascii.read(d / table, format="cds", readme=str(d / "ReadMe"))
    df = tab.to_pandas()
    for c in df.columns:
        if c in _ID_COLUMNS:
            # identifiers stay integers: a float64 cannot hold a 19-digit Gaia
            # source_id exactly (it has 15-16 significant digits)
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
        elif df[c].dtype.kind in "iuf" or str(df[c].dtype) in ("Int64", "Int32", "Float64"):
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("float64")
    return df


def _nan_sentinels(df: pd.DataFrame, cols, bad=(-9999.0, -999.99, -9999.99)) -> pd.DataFrame:
    for c in cols:
        if c in df and df[c].dtype.kind == "f":
            v = df[c].to_numpy()
            mask = np.zeros(len(v), bool)
            for b in bad:
                mask |= np.isclose(v, b)
            df.loc[mask, c] = np.nan
    return df


@lru_cache(maxsize=1)
def yu2018() -> pd.DataFrame:
    """Yu et al. 2018, tables 1 and 2 joined on KIC."""
    t1 = _read_cds("yu2018", "table1.dat")
    t2 = _read_cds("yu2018", "table2.dat")
    t1 = t1.rename(columns={"Kpmag": "kepmag", "QLen": "n_quarters", "DLen": "length_days",
                            "numax": "numax", "e_numax": "numax_err", "Delnu": "dnu",
                            "e_Delnu": "dnu_err", "A": "amp", "e_A": "amp_err",
                            "Width": "width", "e_Width": "width_err", "Gran": "gran",
                            "e_Gran": "gran_err"})
    t2 = t2.rename(columns={"Teff": "teff", "e_Teff": "teff_err", "log(g)": "logg",
                            "e_log(g)": "logg_err", "[Fe/H]": "feh", "e_[Fe/H]": "feh_err",
                            "NoCorM": "mass_nocorr", "e_NoCorM": "mass_nocorr_err",
                            "NoCorR": "radius_nocorr", "e_NoCorR": "radius_nocorr_err",
                            "RGBcorM": "mass_rgbcorr", "e_RGBcorM": "mass_rgbcorr_err",
                            "RGBcorR": "radius_rgbcorr", "e_RGBcorR": "radius_rgbcorr_err",
                            "ClcorM": "mass_clcorr", "e_ClcorM": "mass_clcorr_err",
                            "ClcorR": "radius_clcorr", "e_ClcorR": "radius_clcorr_err",
                            "Phase": "phase"})
    df = t1.merge(t2, on="KIC", how="outer").rename(columns={"KIC": "kic"})
    df["kic"] = df["kic"].astype(int)
    # Phase: 2 = core helium burning, 1 = red giant branch, 0 = unclassified.
    df["state"] = df["phase"].map({1.0: "RGB", 2.0: "RC"}).fillna("unclassified")
    return df.set_index("kic")


@lru_cache(maxsize=1)
def apokasc3() -> pd.DataFrame:
    """APOKASC-3 table 4 (recommended values) with selected table 5 columns."""
    t4 = _read_cds("apokasc3", "table4.dat")
    t4 = t4.rename(columns={
        "KIC": "kic", "EvolSt": "state", "EvolSrc": "state_source", "CatTab": "category",
        "SeisSrc": "seis_weighting", "SpecSrc": "spec_release", "NNumax": "n_numax",
        "NDNu": "n_dnu", "Quar": "n_quarters", "Numax": "numax", "e_Numax": "numax_err",
        "DNu": "dnu", "e_DNu": "dnu_err", "fDNu": "f_dnu", "e_fDNu": "f_dnu_err",
        "fNumax": "f_numax", "Mass": "mass", "e_Mass": "mass_err", "Radius": "radius",
        "e_Radius": "radius_err", "loggSeis": "logg_seis", "e_loggSeis": "logg_seis_err",
        "Teff": "teff", "e_Teff": "teff_err", "loggSpec": "logg_spec",
        "e_loggSpec": "logg_spec_err", "[Fe/H]": "feh", "e_[Fe/H]": "feh_err",
        "[a/Fe]": "alpha_fe", "e_[a/Fe]": "alpha_fe_err", "InvRGaia": "inv_r_gaia",
        "e_InvRGaia": "inv_r_gaia_err", "GaiaDR3": "gaia_dr3", "2MASS": "tmass"})
    keep = [c for c in t4.columns if not c.startswith("[") and not c.startswith("e_[")
            and not c.startswith("Age") and not c.startswith("E_Age") and not c.startswith("e_Age")]
    t4 = t4[keep]
    t5 = _read_cds("apokasc3", "table5.dat").rename(columns={
        "KIC": "kic", "plx": "plx_apokasc", "e_plx": "plx_err_apokasc", "Ksmag": "ksmag",
        "e_Ksmag": "ksmag_err", "fDNuSH": "f_dnu_sharma", "MSH": "mass_sharma",
        "RSH": "radius_sharma", "fDNuWh": "f_dnu_white", "MWh": "mass_white",
        "RWh": "radius_white"})
    t5 = t5[["kic", "plx_apokasc", "plx_err_apokasc", "ksmag", "ksmag_err", "f_dnu_sharma",
             "mass_sharma", "radius_sharma", "f_dnu_white", "mass_white", "radius_white"]]
    df = t4.merge(t5, on="kic", how="left")
    df = _nan_sentinels(df, df.columns)
    for c in ("plx_apokasc", "plx_err_apokasc", "ksmag", "ksmag_err"):
        df.loc[df[c] == 0, c] = np.nan
    df["kic"] = df["kic"].astype(int)
    for c in ("state", "category", "state_source", "tmass"):
        df[c] = df[c].astype(str).str.strip()
    df["gaia_dr3"] = df["gaia_dr3"].astype("Int64")
    df["radius_gaia"] = 1.0 / df["inv_r_gaia"]
    df["radius_gaia_err"] = df["inv_r_gaia_err"] / df["inv_r_gaia"] ** 2
    return df.set_index("kic")


#: Pipelines in APOKASC-3 table 9 (their Appendix A).
APOKASC_PIPELINES = ("COR", "ELS", "GAU", "SYD", "A2Z", "DIA", "CAN", "CV", "HON", "FLI")


@lru_cache(maxsize=1)
def apokasc3_pipelines() -> pd.DataFrame:
    """Every pipeline's zero-point-adjusted numax and Dnu (APOKASC-3 table 9)."""
    t9 = _read_cds("apokasc3", "table9.dat")
    t9 = _nan_sentinels(t9, t9.columns)
    t9 = t9.rename(columns={"KIC": "kic"})
    t9["kic"] = t9["kic"].astype(int)
    return t9.set_index("kic")


@lru_cache(maxsize=1)
def vrard2016() -> pd.DataFrame:
    t = _read_cds("vrard2016", "table2.dat").rename(columns={
        "KIC": "kic", "Dnu": "dnu", "DPi1": "dpi1", "e_DPi1": "dpi1_err", "q": "q",
        "M": "mass", "e_M": "mass_err", "Alias": "alias", "Measure": "n_codes",
        "Status": "status"})
    t["kic"] = t["kic"].astype(int)
    # 0 = RGB, 1 = red clump, 2 = secondary clump
    t["state"] = t["status"].map({0.0: "RGB", 1.0: "RC", 2.0: "2CL"}).fillna("unknown")
    return t.set_index("kic")


@lru_cache(maxsize=1)
def lund2017() -> tuple[pd.DataFrame, pd.DataFrame]:
    """LEGACY global parameters (table 1) and peak-bagged modes (table 6)."""
    t1 = _read_cds("lund2017", "table1.dat")
    t6 = _read_cds("lund2017", "table6.dat")
    t1 = t1.rename(columns={"KIC": "kic", "Name": "name", "numax": "numax", "Dnu": "dnu",
                            "Teff": "teff", "e_Teff": "teff_err", "[Fe/H]": "feh"})
    t6 = t6.rename(columns={"KIC": "kic", "Freq": "freq", "e_Freq": "freq_err_lo",
                            "E_Freq": "freq_err_hi"})
    t1["kic"] = t1["kic"].astype(str).str.strip()
    t6["kic"] = t6["kic"].astype(str).str.strip()
    return t1, t6
