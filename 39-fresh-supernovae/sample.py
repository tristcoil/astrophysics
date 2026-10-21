"""The ZTF DR2 sample: one table joining the release's metadata, host galaxy
properties, the ZTF team's own fits, and what we derive per supernova
(Milky Way dust, CMB-frame redshift, overlap with Pantheon+).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config
from cosmology import z_hel_to_cmb

#: Redshift sources (Rigault et al. 2025, Sect. 4.5). "galaxy": catalogue or
#: host emission lines in a non-SEDm spectrum, sigma_z <~ 1e-4. "sn": template
#: matching to the supernova's own spectral features (snid), sigma_z ~ 3e-3,
#: i.e. ~900 km/s. Host lines in low-resolution SEDm spectra are in between.
Z_SOURCE_GROUP = {
    "z_gal": "galaxy", "z_nonSEDm": "galaxy", "z_override": "galaxy",
    "z_SEDm": "SEDm host lines", "z_snid": "SN features", "z_unknown": "unknown",
}

#: DR2 sub-types that cosmology analyses drop (Rigault et al. 2025; Dimitriadis
#: et al. 2025): everything the release labels "snia-pec". 91T-like and
#: 99aa-like stay in, as in Pantheon+.
PECULIAR_TYPE = "snia-pec"


def mw_ebv(ra, dec) -> np.ndarray:
    """Milky Way E(B-V): the SFD map scaled by 0.86 (Schlafly & Finkbeiner 2011)."""
    from sfdmap2 import sfdmap

    m = sfdmap.SFDMap(str(config.DUST_DIR), scaling=config.SFD_SCALING)
    return np.asarray(m.ebv(np.asarray(ra, float), np.asarray(dec, float)), dtype=float)


def load_pantheon() -> pd.DataFrame:
    return pd.read_csv(config.PANTHEON_DIR / "Pantheon+SH0ES.dat", sep=r"\s+")


def load_dr2() -> pd.DataFrame:
    """All 3,628 DR2 supernovae with metadata, hosts and the team's fits."""
    t = config.ZTF_DIR / "tables"
    sn = pd.read_csv(t / "snia_data.csv", index_col=0).set_index("ztfname")
    s3 = pd.read_csv(t / "extra/snia_data_phase-10to40_gri_salt3.csv", index_col=0) \
        .set_index("ztfname")
    s2gr = pd.read_csv(t / "extra/snia_data_phase-10to40_gr_salt2-T21.csv", index_col=0) \
        .set_index("ztfname")
    ghost = pd.read_csv(t / "globalhost_data.csv", index_col=0)
    lhost = pd.read_csv(t / "localhost_data.csv", index_col=0)

    fit_cols = ["t0", "x0", "x1", "c", "t0_err", "x0_err", "x1_err", "c_err",
                "cov_x0_x1", "cov_x0_c", "cov_x1_c", "fitprob",
                "lccoverage_flag", "fitquality_flag"]
    meta = sn[["redshift", "redshift_err", "source", "ra", "dec", "sn_type", "sub_type",
               "iau_name", "mwebv"]].rename(columns={"mwebv": "mwebv_dr2"})
    df = meta.join(sn[fit_cols].add_prefix("t21_")) \
             .join(s3[fit_cols].add_prefix("s3zt_")) \
             .join(s2gr[["x1", "c", "fitquality_flag"]].add_prefix("t21gr_")) \
             .join(ghost[["mass", "mass_err", "restframe_gz", "d_dlr"]].add_prefix("host_")) \
             .join(lhost[["mass", "restframe_gz"]].add_prefix("local_"))
    for prefix in ("t21_", "s3zt_"):
        df[prefix + "mB"] = -2.5 * np.log10(df[prefix + "x0"]) + 10.635
    df["z_group"] = df["source"].map(Z_SOURCE_GROUP).fillna("unknown")
    df["peculiar"] = df["sn_type"].eq(PECULIAR_TYPE)
    df["mwebv"] = mw_ebv(df["ra"], df["dec"])
    df["z_cmb"] = z_hel_to_cmb(df["redshift"], df["ra"], df["dec"])
    df["in_pantheon"] = _pantheon_overlap(df)
    df["has_lightcurve"] = [(config.ZTF_LC_DIR / f"{n}_lc.csv.gz").exists() for n in df.index]
    return df


def _pantheon_overlap(df: pd.DataFrame, radius_arcsec: float = 5.0) -> pd.Series:
    """Is the supernova also in Pantheon+ (D1's sample)? Matched by position."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    p = load_pantheon()
    cz = SkyCoord(df["ra"].to_numpy() * u.deg, df["dec"].to_numpy() * u.deg)
    cp = SkyCoord(p["RA"].to_numpy() * u.deg, p["DEC"].to_numpy() * u.deg)
    _, sep, _ = cz.match_to_catalog_sky(cp)
    return pd.Series(sep.arcsec < radius_arcsec, index=df.index)


def fit_tasks(df: pd.DataFrame) -> list[dict]:
    """What the SALT fitter needs per supernova: heliocentric z and MW dust."""
    return [{"name": n, "z": float(r.redshift), "mwebv": float(r.mwebv)}
            for n, r in df.iterrows()]
