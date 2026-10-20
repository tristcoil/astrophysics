"""Load the catalogues into tidy tables with one naming convention.

Each published table names its columns its own way; here they all become
``ra``, ``dec`` (host galaxy, degrees), ``z``, ``las_arcmin`` (angular
length), ``lp_mpc`` (projected proper length) and so on. Missing values
become NaN, never a placeholder zero.
"""

from __future__ import annotations

import gzip
from functools import lru_cache

import numpy as np
import pandas as pd
from astropy.io import ascii

import config
import cosmology
from config import CATALOG_DIR, DERIVED_DIR


def _read_cds(data_file: str, readme_file: str, table_name: str) -> pd.DataFrame:
    """Read a CDS fixed-width table using its ReadMe for the column layout.

    Our local file names differ from the archive's, so the reader is told
    which table of the ReadMe this is.
    """
    path = CATALOG_DIR / data_file
    if not path.exists():
        raise FileNotFoundError(f"{path} missing; run fetch_data.py --only catalogs")
    opener = gzip.open if data_file.endswith(".gz") else open
    with opener(path, "rt", encoding="latin-1") as fh:
        lines = fh.read().splitlines()
    reader = ascii.get_reader(ascii.Cds, readme=str(CATALOG_DIR / readme_file))
    reader.data.table_name = table_name
    table = reader.read(lines)
    df = table.to_pandas()
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip().replace({"": np.nan, "--": np.nan, "nan": np.nan})
    return df


@lru_cache(maxsize=None)
def oei2023() -> pd.DataFrame:
    """The 2,060 giants Oei et al. (2023) found by eye in LoTSS DR2."""
    raw = _read_cds("oei2023_tablef1.dat", "oei2023_ReadMe", "tablef1.dat")
    cols = {c.lower(): c for c in raw.columns}

    def col(*names):
        for n in names:
            if n.lower() in cols:
                return raw[cols[n.lower()]]
        raise KeyError(names)

    df = pd.DataFrame({
        "host": col("SDSS-DR12", "SDSS"),
        "ra": col("RAdeg", "RAJ2000").astype(float),
        "dec": col("DEdeg", "DEJ2000").astype(float),
        "z": col("z").astype(float),
        "z_err": col("e_z").astype(float),
        "z_spec": col("zsp?", "zsp").astype(float) == 1,
        "host_is_nearest_candidate": col("HNC?", "HNC").astype(float) == 1,
        "las_arcmin": col("LengthAng").astype(float),
        "lp_mpc": col("LengthMpc").astype(float),
        "lp_err_mpc": col("e_LengthMpc").astype(float),
        "mstar_1e11": col("M*", "Mstar").astype(float),
        "mbh_1e9": col("MSMBH").astype(float),
        "mbh_err_1e9": col("e_MSMBH").astype(float),
        "host_qso": col("HostQSO?", "HostQSO").astype(float),
    })
    # The archive stores "no estimate" as a blank; make sure no zero slipped in.
    for c in ("mstar_1e11", "mbh_1e9", "mbh_err_1e9"):
        df.loc[df[c] <= 0, c] = np.nan
    return df


@lru_cache(maxsize=None)
def mostert2024() -> pd.DataFrame:
    """All 11,585 known giants compiled by Mostert et al. (2024)."""
    raw = _read_cds("mostert2024_catalog.dat.gz", "mostert2024_ReadMe", "catalog.dat")
    df = pd.DataFrame({
        "name": raw["Name"],
        "ra": raw["RAdeg"].astype(float),
        "dec": raw["DEdeg"].astype(float),
        "lp_mpc": raw["LengthMpc"].astype(float),
        "lp_err_mpc": raw["e_LengthMpc"].astype(float),
        "las_arcmin": raw["LengthAng"].astype(float),
        "z": raw["z"].astype(float),
        "z_err": raw["e_z"].astype(float),
        "ref": raw["Ref"],
        "first_found_by": raw["FistDisc"] if "FistDisc" in raw else raw.iloc[:, -1],
    })
    method = {"RGZ (Hardcastle et al. 2023)": "LOFAR Galaxy Zoo (citizen science)",
              "ML (this work)": "machine learning (Mostert et al.)",
              "Oei et al. 2023": "by eye (Oei et al. 2023)"}
    df["found_by"] = df["first_found_by"].map(method).fillna("earlier literature")
    return df


def all_giants() -> pd.DataFrame:
    """Every known giant (Mostert et al.'s compilation includes Oei et al.'s)."""
    return mostert2024()


@lru_cache(maxsize=None)
def oei2024() -> dict[str, pd.DataFrame]:
    """Cosmic Web density and dynamical class from BORG SDSS (Oei et al. 2024)."""
    out = {}
    for key, fname in (("giants", "tab1grs.dat"), ("radiogalaxies", "tab1rg.dat")):
        raw = _read_cds(f"oei2024_{fname}", "oei2024_ReadMe", fname)
        df = pd.DataFrame({
            "ra": raw["RAdeg"].astype(float), "dec": raw["DEdeg"].astype(float),
            "z": raw["zspHost"].astype(float),
            "mu_fixed": raw["muFix"].astype(float), "sigma2_fixed": raw["sigma2Fix"].astype(float),
            "mu_flex": raw["muFlex"].astype(float), "sigma2_flex": raw["sigma2Flex"].astype(float),
            "p_cluster": raw["Pcluster"].astype(float), "p_filament": raw["Pfilament"].astype(float),
            "p_sheet": raw["Psheet"].astype(float), "p_void": raw["Pvoid"].astype(float),
            "cluster_mass_1e14": raw["MassHost"].astype(float),
            "is_bcg": raw["PBCG"].astype(float) == 1,
        })
        out[key] = df
    return out


@lru_cache(maxsize=None)
def andernach2025() -> pd.DataFrame:
    """The 142 giants larger than 3 Mpc, re-examined by Andernach et al. (2025)."""
    raw = _read_cds("andernach2025_tablea1.dat", "andernach2025_ReadMe", "tablea1.dat")
    return pd.DataFrame({
        "jname": raw["JName"], "ra": raw["RAdeg"].astype(float), "dec": raw["DEdeg"].astype(float),
        "z": raw["z"].astype(float), "z_kind": raw["ztype"], "host_class": raw["Class"],
        "host": raw["HostName"], "host_uncertain": raw["hflag"].astype(str).str.contains(r"\?"),
        "las_arcmin": raw["LAS"].astype(float), "lls_mpc": raw["LLS"].astype(float),
        "survey": raw["survFL"], "note": raw["Note"],
    })


# --------------------------------------------------------------------------
# The LoTSS DR2 value-added catalogue: ordinary radio galaxies
# --------------------------------------------------------------------------

VAC_CACHE = DERIVED_DIR / "vac_resolved.pkl"


@lru_cache(maxsize=None)
def vac() -> pd.DataFrame:
    """Resolved LoTSS DR2 sources with a redshift (Hardcastle et al. 2023)."""
    if VAC_CACHE.exists():
        return pd.read_pickle(VAC_CACHE)
    frames = []
    for lo, hi in config.VAC_DEC_STRIPS:
        path = CATALOG_DIR / f"vac_resolved_dec{lo:+03d}_{hi:+03d}.csv"
        if not path.exists():
            raise FileNotFoundError(f"{path} missing; run fetch_data.py --only catalogs")
        frames.append(pd.read_csv(path, low_memory=False))
    df = pd.concat(frames, ignore_index=True)
    for c in ("Source", "SCode", "r_LAS", "r_zbest"):
        df[c] = df[c].astype(str).str.strip()
    df = df.drop_duplicates("Source").reset_index(drop=True)
    df.to_pickle(VAC_CACHE)
    return df


def radio_galaxies(cosmo: config.Cosmology = config.COSMOLOGY) -> pd.DataFrame:
    """The comparison sample in our own units: projected size and 144-MHz
    luminosity recomputed with our cosmology from the catalogue's angular size,
    flux density and redshift."""
    v = vac()
    z = v["zbest"].to_numpy(float)
    las = v["LAS"].to_numpy(float) / 60.0  # arcsec -> arcmin
    host_ra = v["RAIdeg"].where(np.isfinite(v["RAIdeg"]), v["RAJ2000"]).to_numpy(float)
    host_dec = v["DEIdeg"].where(np.isfinite(v["DEIdeg"]), v["DEJ2000"]).to_numpy(float)
    df = pd.DataFrame({
        "source": v["Source"], "ra": v["RAJ2000"].to_numpy(float), "dec": v["DEJ2000"].to_numpy(float),
        "host_ra": host_ra, "host_dec": host_dec, "has_host": np.isfinite(v["RAIdeg"].to_numpy(float)),
        "z": z, "z_kind": v["r_zbest"], "las_arcmin": las, "las_from": v["r_LAS"],
        "pa_deg": v["CompPA"].to_numpy(float), "n_components": v["Assoc"].to_numpy(float),
        "flux_mjy": v["Ftot"].to_numpy(float),
        "logmass": v["logMass"].where(v["q_logMass"] == 1).to_numpy(float),
        "has_sdss_z": np.isfinite(v["zspecSDSS"].to_numpy(float)),
        "has_desi_z": np.isfinite(v["zDESI"].to_numpy(float)),
        "size_kpc_catalogue": v["Size"].to_numpy(float),
        "L144_catalogue": v["L144"].to_numpy(float),
    })
    df["lp_mpc"] = cosmology.projected_length_mpc(df["las_arcmin"].to_numpy(), z, cosmo)
    df["L144"] = cosmology.radio_luminosity(df["flux_mjy"].to_numpy() * 1e-3, z, cosmo, alpha=-0.7)
    df["is_giant_candidate"] = df["lp_mpc"] >= config.GIANT_THRESHOLD_MPC
    return df
