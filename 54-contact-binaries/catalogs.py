"""Readers for the cached catalogues (ZTF, Gaia, Chen et al. 2018, AllWISE).

Each reader returns a pandas DataFrame with lower-case, project-wide column
names, so the analysis modules never see archive-specific spellings.
"""

from __future__ import annotations

import gzip
import io

import numpy as np
import pandas as pd

import config as C

_ZTF_RENAME = {"ID": "ztf_id", "RAJ2000": "ra", "DEJ2000": "dec", "Per": "period", "R21": "r21",
               "phi21": "phi21", "T0": "t0_mjd", "gmag": "gmag", "rmag": "rmag", "Ng": "ng", "Nr": "nr",
               "gAmp": "amp_g", "rAmp": "amp_r", "log(FAP-g)": "logfap_g", "log(FAP-r)": "logfap_r",
               "R2-g": "r2_g", "R2-r": "r2_r", "Type": "type", "Dmin-g": "dmin_g", "Dmin-r": "dmin_r"}


def _read_vizier_tsv(path) -> pd.DataFrame:
    """VizieR ASU-TSV with '#' comments and three header lines (names, units, dashes)."""
    lines = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]
    names = lines[0].split("\t")
    body = "\n".join(lines[3:])
    df = pd.read_csv(io.StringIO(body), sep="\t", names=names, dtype=str, keep_default_na=False)
    for c in df.columns:
        s = df[c].str.strip()
        num = pd.to_numeric(s, errors="coerce")
        # A column is numeric if every non-blank entry parses.
        if num.notna().sum() == (s != "").sum():
            df[c] = num
        else:
            df[c] = s
    return df


def read_ztf(types=None) -> pd.DataFrame:
    """The ZTF periodic-variable catalogue (Chen et al. 2020), all 24 RA slices.

    ``zid`` is our own integer key (row order after sorting by ZTF ID), used to
    join the cross-matches back.
    """
    cache = C.DERIVED_DIR / "ztf_periodic.parquet"
    if cache.exists():
        df = pd.read_parquet(cache)
    else:
        parts = [_read_vizier_tsv(p) for p in sorted(C.ZTF_DIR.glob("ztf_periodic_ra*.tsv"))]
        if len(parts) != 24:
            raise FileNotFoundError(f"expected 24 ZTF slices, found {len(parts)}: run fetch_data.py ztf")
        df = pd.concat(parts, ignore_index=True).rename(columns=_ZTF_RENAME)
        df["type"] = df["type"].astype(str).str.strip()
        df = df.sort_values("ztf_id").reset_index(drop=True)
        if df["ztf_id"].duplicated().any():
            raise RuntimeError("duplicate ZTF IDs across RA slices")
        df["zid"] = np.arange(len(df), dtype=np.int64)
        df.to_parquet(cache)
    if types is not None:
        df = df[df["type"].isin(types)].reset_index(drop=True)
    return df


_XM_RENAME = {"angDist": "sep_arcsec", "Source": "source_id", "RAdeg": "ra_gaia", "DEdeg": "dec_gaia",
              "Plx": "parallax", "e_Plx": "parallax_error", "pmRA": "pmra", "pmDE": "pmdec", "RUWE": "ruwe",
              "Gmag": "g", "e_Gmag": "g_error", "BPmag": "bp", "RPmag": "rp", "BP-RP": "bp_rp",
              "E(BP/RP)": "bp_rp_excess", "nueff": "nu_eff", "pscol": "pseudocolour", "Solved": "solved",
              "VarFlag": "var_flag", "Teff": "teff_gspphot", "AG": "ag_gspphot",
              "Dist": "dist_gspphot", "NBPblend": "n_bp_blend", "NRPblend": "n_rp_blend",
              "IPDfmp": "ipd_frac_multi_peak", "Dup": "duplicated_source"}


def read_xmatch() -> pd.DataFrame:
    """All Gaia DR3 sources within the match radius of each ZTF binary (CDS xMatch)."""
    parts = []
    for p in sorted(C.GAIA_DIR.glob("xmatch_ztf_gaiadr3_*.csv.gz")):
        with gzip.open(p) as fh:
            parts.append(pd.read_csv(fh, low_memory=False))
    if not parts:
        raise FileNotFoundError("no xMatch chunks: run fetch_data.py xmatch")
    df = pd.concat(parts, ignore_index=True).rename(columns=_XM_RENAME)
    # Keep 19-digit Gaia IDs exact (a float would round them).
    df["source_id"] = df["source_id"].astype(np.int64)
    return df


def read_chen2018() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Chen et al. (2018): the 183 TGAS calibrators (table 1) and 55,603 candidates (table 4)."""
    t1 = _read_vizier_tsv(C.ZTF_DIR / "J_ApJ_859_140_table1.tsv")
    t4 = _read_vizier_tsv(C.ZTF_DIR / "J_ApJ_859_140_table4.tsv")
    return t1, t4


def read_allwise() -> pd.DataFrame:
    """AllWISE match (nearest within 2") for each ZTF binary with a good parallax."""
    parts = []
    for p in sorted(C.GAIA_DIR.glob("xmatch_allwise_*.csv.gz")):
        with gzip.open(p) as fh:
            parts.append(pd.read_csv(fh, dtype={"ccf": str, "qph": str}))
    df = pd.concat(parts, ignore_index=True).sort_values(["zid", "angDist"])
    n = df.groupby("zid")["zid"].transform("size")
    df = df.groupby("zid").head(1).copy()
    df["n_allwise"] = n.loc[df.index].to_numpy()
    df = df.rename(columns={"angDist": "allwise_sep", "W1mag": "w1", "e_W1mag": "w1_err", "W2mag": "w2",
                            "e_W2mag": "w2_err", "AllWISE": "allwise_id"})
    df["w1_clean"] = (df["ccf"].fillna("").str[:2] == "00") & (df["ex"] == 0) & (df["qph"].fillna("").str[0] == "A") \
        & (df["n_allwise"] == 1)
    return df.reset_index(drop=True)


def read_gaia_veb() -> pd.DataFrame:
    """Gaia DR3 eclipsing-binary candidates with P <= 2 d (VizieR I/358/veb)."""
    cache = C.DERIVED_DIR / "gaia_veb.parquet"
    if cache.exists():
        return pd.read_parquet(cache)
    parts = [_read_vizier_tsv(p) for p in sorted(C.GAIA_DIR.glob("veb_ra*.tsv"))]
    if len(parts) != 24:
        raise FileNotFoundError(f"expected 24 VizieR I/358/veb slices, found {len(parts)}")
    df = pd.concat(parts, ignore_index=True)
    df = df.rename(columns={"Source": "source_id", "Freq": "frequency", "ModelType": "model_type", "Rank": "rank",
                            "RA_ICRS": "ra", "DE_ICRS": "dec"})
    df["source_id"] = df["source_id"].astype(np.int64)
    df["period"] = 1.0 / df["frequency"]
    df.to_parquet(cache)
    return df
