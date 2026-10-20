"""The eRASS1 cluster samples, joined into one table per cluster.

* The **cosmology sample** (Bulbul et al. 2024): 5,259 clusters with extent
  likelihood > 6, photometric redshift 0.1-0.8 and richness > 3, in the
  12,791 deg^2 overlap of eROSITA-DE and the Legacy Surveys DR10-South. Its
  count rate CR500_NH0 is the observable of Ghirardini et al.'s analysis: the
  0.2-2.3 keV photons per second within R500, corrected for Galactic absorption.
* From the **primary catalogue** we add each cluster's exposure, its
  probability of being a contaminant (an AGN or a background fluctuation
  rather than a cluster; PCONT, from the paper's mixture model), its mass and
  its matches to other catalogues (including Planck's PSZ2).
* From the **optical catalogue** (Kluge et al. 2024) we add richness, the
  HI4PI column density and the optical depth limit ZVLIM. Its NHI column is
  labelled 1e21 cm^-2 but holds units of 1e22 cm^-2: so read, it matches our
  own HI4PI map at the same positions to 0.2% (median).

63 cosmology-sample clusters are missing from the primary catalogue (no
counterpart within 60 arcsec; several sit 2-4 arcmin from a primary cluster at
the same redshift, i.e. re-centred detections). Their exposure, N_H and ZVLIM
come from the maps and their contamination probability is the sample mean.
14 clusters have a count rate of exactly zero (failed fits); they are dropped.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import healpy as hp
from astropy.table import Table
from scipy.spatial import cKDTree

import config as C
import footprint as F

NHI_UNIT = 1e22  # see module docstring


def _strip(col):
    return np.char.strip(np.asarray(col).astype(str))


def load_primary() -> pd.DataFrame:
    """All 12,247 clusters of the primary catalogue (for maps and masses)."""
    P = Table.read(C.RAW_DIR / "erass1cl_primary_v3.2.fits")
    df = pd.DataFrame({
        "name": _strip(P["NAME"]), "detuid": _strip(P["DETUID"]),
        "ra": np.asarray(P["RA"], float), "dec": np.asarray(P["DEC"], float),
        "z": np.asarray(P["BEST_Z"], float), "z_type": _strip(P["BEST_Z_TYPE"]),
        "ext_like": np.asarray(P["EXT_LIKE"], float), "exp": np.asarray(P["EXP"], float),
        "pcont": np.asarray(P["PCONT"], float),
        "m500": np.asarray(P["M500"], float) * 1e13, "m500_lo": np.asarray(P["M500_L"], float) * 1e13,
        "m500_hi": np.asarray(P["M500_H"], float) * 1e13,
        "match": _strip(P["MATCH_NAME"]),
    })
    df.loc[df.m500 <= 0, ["m500", "m500_lo", "m500_hi"]] = np.nan
    return df


def primary_mass_pdfs():
    """The catalogue's M500 probability densities: (masses grid, pdf rows, detuid)."""
    P = Table.read(C.RAW_DIR / "erass1cl_primary_v3.2.fits")
    return np.asarray(P["M500_PDF_array"], float), np.asarray(P["M500_PDF"], float), _strip(P["DETUID"])


def _nearest_zvlim(ra, dec):
    K = Table.read(C.RAW_DIR / "eRASS1_clusters_optical.fits")
    run = _strip(K["RUN"])
    zv = np.asarray(K["ZVLIM_02"], float)
    ok = np.isfinite(zv) & np.char.startswith(run, "legacy_dr10")
    tree = cKDTree(np.array(hp.ang2vec(np.asarray(K["RA"], float)[ok], np.asarray(K["DEC"], float)[ok], lonlat=True)))
    _, j = tree.query(np.array(hp.ang2vec(ra, dec, lonlat=True)))
    return zv[ok][j]


def load_cosmology_sample(fp: dict | None = None) -> pd.DataFrame:
    """The 5,245 cosmology-sample clusters with a measured count rate, joined."""
    Cc = Table.read(C.RAW_DIR / "erass1cl_cosmology_v1.1.fits")
    df = pd.DataFrame({
        "name": _strip(Cc["NAME"]), "detuid": _strip(Cc["DETUID"]),
        "ra": np.asarray(Cc["RA"], float), "dec": np.asarray(Cc["DEC"], float),
        "z": np.asarray(Cc["Z_LAMBDA"], float), "z_err": np.asarray(Cc["Z_LAMBDA_E"], float),
        "cr": np.asarray(Cc["CR500_NH0"], float), "cr_lo": np.asarray(Cc["CR500_L_NH0"], float),
        "cr_hi": np.asarray(Cc["CR500_H_NH0"], float), "cr_abs": np.asarray(Cc["CR500"], float),
        "l500": np.asarray(Cc["L500"], float) * 1e42, "kt": np.asarray(Cc["KT"], float),
    })
    n_all = len(df)
    prim = load_primary().set_index("detuid")
    K = Table.read(C.RAW_DIR / "eRASS1_clusters_optical.fits")
    opt = pd.DataFrame({"detuid": _strip(K["DETUID"]), "lambda_norm": np.asarray(K["LAMBDA_NORM"], float),
                        "nh": np.asarray(K["NHI"], float) * NHI_UNIT, "zvlim": np.asarray(K["ZVLIM_02"], float)}
                       ).set_index("detuid")
    df = df.join(prim[["exp", "pcont", "ext_like", "m500", "m500_lo", "m500_hi", "match"]], on="detuid")
    df = df.join(opt, on="detuid")
    df["in_primary"] = df["exp"].notna()
    # Clusters absent from the primary catalogue: conditions from the maps.
    fp = F.build() if fp is None else fp
    miss = ~df.in_primary
    if miss.any():
        t, nh, _, _ = F.conditions_at(fp, df.ra[miss], df.dec[miss])
        df.loc[miss, "exp"] = t
        df.loc[miss, "nh"] = F.nh_map(df.ra[miss].values, df.dec[miss].values)
        df.loc[miss, "zvlim"] = _nearest_zvlim(df.ra[miss].values, df.dec[miss].values)
        df.loc[miss, "pcont"] = df.pcont[~miss].mean()
        df.loc[miss, "match"] = ""
    df["weight"] = 1.0 - df.pcont  # expected number of real clusters per entry
    df["counts"] = df.cr_abs * df.exp  # photons detected (absorbed rate x exposure)
    _, _, _, inside = F.conditions_at(fp, df.ra, df.dec)
    df["in_mask"] = inside
    df = df[df.cr > 0].reset_index(drop=True)
    df.attrs.update({"n_published": n_all, "n_zero_rate": int(n_all - len(df)), "n_missing_primary": int(miss.sum())})
    return df


def measurement_error_fit(df: pd.DataFrame) -> dict:
    """Count-rate error vs photons detected: sigma_ln = sqrt(a / N)."""
    ok = (df.cr_lo > 0) & (df.cr_hi > df.cr_lo) & (df.counts > 0)
    sig = 0.5 * np.log(df.cr_hi[ok] / df.cr_lo[ok])
    n = df.counts[ok]
    a = float(np.median(sig**2 * n))
    return {"a": a, "n_used": int(ok.sum()), "poisson_ratio": float(np.sqrt(a)),
            "binned": [(lo, hi, float(np.median(sig[(n >= lo) & (n < hi)])), float(np.sqrt(a / np.sqrt(lo * hi))))
                       for lo, hi in ((10, 20), (20, 40), (40, 80), (80, 160), (160, 400), (400, 2000))]}


def absorption_fit(df: pd.DataFrame) -> dict:
    """Fraction of the count rate surviving Galactic absorption vs N_H."""
    from scipy.optimize import curve_fit

    ok = (df.cr_abs > 0) & np.isfinite(df.nh) & df.in_primary
    a = (df.cr_abs / df.cr)[ok]
    nh = df.nh[ok]
    popt, _ = curve_fit(lambda x, k, p: -k * (x / 1e20) ** p, nh, np.log(a), p0=[0.05, 0.8])
    resid = np.log(a) + popt[0] * (nh / 1e20) ** popt[1]
    return {"k": float(popt[0]), "p": float(popt[1]), "rms_ln": float(resid.std()), "n": int(ok.sum())}


def load_psz2() -> pd.DataFrame:
    P = Table.read(C.RAW_DIR / "HFI_PCCS_SZ-union_R2.08.fits")
    return pd.DataFrame({
        "name": _strip(P["NAME"]), "ra": np.asarray(P["RA"], float), "dec": np.asarray(P["DEC"], float),
        "glon": np.asarray(P["GLON"], float), "glat": np.asarray(P["GLAT"], float),
        "snr": np.asarray(P["SNR"], float), "z": np.asarray(P["REDSHIFT"], float),
        "msz": np.asarray(P["MSZ"], float) * 1e14, "msz_up": np.asarray(P["MSZ_ERR_UP"], float) * 1e14,
        "msz_low": np.asarray(P["MSZ_ERR_LOW"], float) * 1e14, "cosmo": np.asarray(P["COSMO"], bool),
    })


def match_psz2(primary: pd.DataFrame, psz2: pd.DataFrame) -> pd.DataFrame:
    """eRASS1 clusters that Bulbul et al. matched to a PSZ2 detection."""
    rows = []
    by_name = {n: i for i, n in enumerate(psz2.name)}
    for i, m in enumerate(primary.match):
        for tok in m.split(","):
            if tok.startswith("PSZ2_"):
                j = by_name.get(tok[5:].strip())
                if j is not None:
                    rows.append((i, j))
    i, j = (np.array(v) for v in zip(*rows))
    out = primary.iloc[i].reset_index(drop=True).add_prefix("x_")
    return pd.concat([out, psz2.iloc[j].reset_index(drop=True).add_prefix("sz_")], axis=1)
