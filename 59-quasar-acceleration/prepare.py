"""Turn the raw downloads into one analysis table (data/cache/crf3.parquet).

* Gaia's null markers become NaN (FITS integer nulls arrive as -32768 etc.).
* One colour column for every source: nu_eff (the colour the astrometry
  assumed) for five-parameter solutions, the fitted pseudocolour for
  six-parameter ones.
* SDSS DR16Q spectroscopic redshifts are attached by position (within 1").
* Nothing is cut here: every selection happens in the analysis, where it is
  recorded.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from astropy.io import fits

from config import CACHE_DIR, DR16Q_MATCH_ARCSEC, RAW_DIR

CRF3_PARQUET = CACHE_DIR / "crf3.parquet"

#: FITS integer null values used by the archive for left-joined, missing rows.
_INT_NULLS = (-32768, -2147483648, -9223372036854775808)


def _column(data, name: str) -> np.ndarray:
    col = np.asarray(data[name])
    if col.dtype.kind in "iu":
        col = col.astype(np.int64)
        if np.isin(col, _INT_NULLS).any():
            out = col.astype(float)
            out[np.isin(col, _INT_NULLS)] = np.nan
            return out
        return col
    if col.dtype.kind == "f":
        return col.astype(np.float64 if name in ("ra", "dec", "l", "b", "pmra", "pmdec",
                                                 "parallax", "ecl_lat") else np.float32)
    if col.dtype.kind in "SU":
        return np.char.strip(col.astype(str))
    return col


def build(force: bool = False) -> pd.DataFrame:
    if CRF3_PARQUET.exists() and not force:
        return pd.read_parquet(CRF3_PARQUET)
    with fits.open(RAW_DIR / "gaia_crf3.fits", memmap=True) as hdul:
        data = hdul[1].data
        df = pd.DataFrame({name.lower(): _column(data, name) for name in data.columns.names})
    df["catalogue_name"] = df["catalogue_name"].astype("category")
    df["vari_best_class_name"] = df["vari_best_class_name"].replace("", np.nan).astype("category")
    df["five_param"] = df["astrometric_params_solved"] == 31
    df["colour"] = np.where(df["five_param"], df["nu_eff_used_in_astrometry"], df["pseudocolour"])

    z, sep = match_dr16q(df["ra"].to_numpy(), df["dec"].to_numpy())
    df["z_dr16q"] = z.astype(np.float32)
    df["dr16q_sep_arcsec"] = sep.astype(np.float32)
    df.to_parquet(CRF3_PARQUET, index=False)
    summary = {
        "rows": len(df), "five_parameter": int(df["five_param"].sum()),
        "with_dr16q_z": int(np.isfinite(df["z_dr16q"]).sum()),
        "with_qsoc_z": int(np.isfinite(df["redshift_qsoc"]).sum()),
    }
    (CACHE_DIR / "crf3_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return df


def match_dr16q(ra: np.ndarray, dec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Nearest DR16Q quasar within DR16Q_MATCH_ARCSEC of each Gaia source."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    q = pd.read_parquet(RAW_DIR / "dr16q_radecz.parquet")
    gaia = SkyCoord(ra * u.deg, dec * u.deg)
    sdss = SkyCoord(q["ra"].to_numpy() * u.deg, q["dec"].to_numpy() * u.deg)
    idx, d2d, _ = gaia.match_to_catalog_sky(sdss)
    sep = d2d.to_value(u.arcsec)
    ok = sep < DR16Q_MATCH_ARCSEC
    z = np.full(ra.size, np.nan)
    z[ok] = q["z"].to_numpy()[idx[ok]]
    # A DR16Q quasar claimed by two Gaia sources keeps only the closer one.
    claimed = pd.DataFrame({"i": idx[ok], "sep": sep[ok], "row": np.flatnonzero(ok)})
    dup = claimed.sort_values("sep").duplicated("i", keep="first")
    z[claimed["row"].to_numpy()[dup.to_numpy()]] = np.nan
    sep_out = np.where(np.isfinite(z), sep, np.nan)
    return z, sep_out


if __name__ == "__main__":
    d = build(force=True)
    print(d.shape)
    print(json.loads((CACHE_DIR / "crf3_summary.json").read_text()))
