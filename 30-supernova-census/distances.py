"""Distances: cosmology, the CMB frame, and Cosmicflows-4 for the nearest hosts.

For most supernovae the redshift is the distance: at z = 0.05 a peculiar
velocity of 300 km/s changes it by 2% (0.04 mag). For the nearest ones it is
not. Their hosts move by several hundred km/s on top of a Hubble flow of only
1,000-3,000 km/s, so a redshift distance can be off by a factor of two. These
few events matter more than their number suggests: in a 1/V_max luminosity
function the faintest, nearest supernovae carry the largest weights.

So below cz = 4,000 km/s each host is looked up in Cosmicflows-4 (Tully et al.
2023; 55,877 galaxies with distances from Cepheids, the tip of the red-giant
branch, surface-brightness fluctuations, Tully-Fisher, supernovae...). The
galaxy's own modulus is used when it is the more precise, else its group's.

Two traps, documented in project J4 (22-laniakea-flows): CF4's Fundamental
Plane moduli for SDSS galaxies carry the wrong sign of the size residual, so a
galaxy (or a group) whose only distance is the Fundamental Plane gets no CF4
distance here; and CF4 distances sit on the Cepheid/TRGB scale (H0 = 74.6),
while the BTS uses H0 = 70, so CF4 moduli are shifted onto the H0 = 70 scale.
"""

from __future__ import annotations

import gzip
from functools import lru_cache

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

import config as C

NON_FP_METHODS = ("DMsnIa", "DMtf", "DMsbf", "DMsnII", "DMtrgb", "DMceph", "DMmas")


@lru_cache(maxsize=1)
def cosmology() -> FlatLambdaCDM:
    return FlatLambdaCDM(H0=C.H0, Om0=C.OMEGA_M)


@lru_cache(maxsize=1)
def _mu_table():
    z = np.concatenate([np.geomspace(1e-5, 0.01, 400), np.linspace(0.01, 3.0, 3000)[1:]])
    return z, cosmology().distmod(z).value


def distance_modulus(z) -> np.ndarray:
    """Luminosity-distance modulus (H0 = 70, Om = 0.3); NaN where z is missing."""
    z = np.asarray(z, float)
    zt, mt = _mu_table()
    out = np.full(z.shape, np.nan)
    ok = np.isfinite(z) & (z > 0)
    # interpolate in log z: mu is close to linear in log z at low redshift
    out[ok] = np.interp(np.log(z[ok]), np.log(zt), mt)
    return out


def redshift_at_modulus(mu) -> np.ndarray:
    """Inverse of distance_modulus."""
    zt, mt = _mu_table()
    return np.exp(np.interp(np.asarray(mu, float), mt, np.log(zt), left=np.nan, right=np.nan))


def _dipole_cos(ra_deg, dec_deg) -> np.ndarray:
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    apex = SkyCoord(l=C.CMB_DIPOLE_L_DEG * u.deg, b=C.CMB_DIPOLE_B_DEG * u.deg,
                    frame="galactic").icrs
    ra, dec, ra0, dec0 = map(np.radians, (ra_deg, dec_deg, apex.ra.deg, apex.dec.deg))
    return np.sin(dec) * np.sin(dec0) + np.cos(dec) * np.cos(dec0) * np.cos(ra - ra0)


def cmb_frame_redshift(z_helio, ra_deg, dec_deg) -> np.ndarray:
    """1 + z_cmb = (1 + z_helio) / (1 - v cos(theta) / c): we move toward the
    dipole apex at 370 km/s, which blueshifts galaxies in that direction."""
    beta = C.CMB_DIPOLE_KMS / C.C_KMS * _dipole_cos(ra_deg, dec_deg)
    return (1.0 + np.asarray(z_helio, float)) / (1.0 - beta) - 1.0


# ---------------------------------------------------------------------------
# Cosmicflows-4
# ---------------------------------------------------------------------------

def _read_cds(name: str) -> pd.DataFrame:
    from astropy.io.ascii import Cds
    reader = Cds(readme=str(C.CF4_DIR / "ReadMe"))
    # astropy matches the table by file name; ours is the gzipped copy
    reader.data.table_name = name
    lines = gzip.open(C.CF4_DIR / f"{name}.gz", "rt").read().splitlines()
    return reader.read(lines).to_pandas()


@lru_cache(maxsize=1)
def load_cf4() -> pd.DataFrame:
    """One row per CF4 galaxy with the distance modulus we would adopt for it."""
    gal = _read_cds("table2.dat")
    grp = _read_cds("table4.dat")[["1PGC", "DMzp", "e_DMzp"]]
    gal["has_non_fp"] = gal[list(NON_FP_METHODS)].notna().any(axis=1)
    grp_ok = gal.groupby("1PGC")["has_non_fp"].any().rename("group_has_non_fp")
    gal = gal.merge(grp, on="1PGC", how="left").merge(grp_ok, on="1PGC", how="left")

    own = gal["has_non_fp"] & (gal["e_DM"] <= gal["e_DMzp"].fillna(np.inf))
    use_grp = ~own & gal["group_has_non_fp"] & gal["DMzp"].notna()
    gal["dm_adopted"] = np.where(own, gal["DM"], np.where(use_grp, gal["DMzp"], np.nan))
    gal["e_dm_adopted"] = np.where(own, gal["e_DM"], np.where(use_grp, gal["e_DMzp"], np.nan))
    gal["dm_source"] = np.where(own, "galaxy", np.where(use_grp, "group", "none (FP only)"))
    return gal


def cf4_distances(df: pd.DataFrame) -> pd.DataFrame:
    """For every transient below cz = 4,000 km/s (CMB frame), the nearest CF4
    galaxy within 40 projected kpc (at least 2') whose velocity agrees to
    600 km/s. Returns distance moduli shifted to the H0 = 70 scale."""
    out = pd.DataFrame(index=df.index, columns=["pgc", "sep_kpc", "mu_h70", "e_mu", "source"],
                       dtype=object)
    for col in ("sep_kpc", "mu_h70", "e_mu"):
        out[col] = np.nan
    cz = df["z_cmb"].to_numpy() * C.C_KMS
    near = np.isfinite(cz) & (cz < C.CF4_MAX_CZ_KMS)
    if not near.any():
        return out
    cf4 = load_cf4()
    cf4 = cf4[cf4["dm_adopted"].notna()].reset_index(drop=True)
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    gal_c = SkyCoord(cf4["RAdeg"].values * u.deg, cf4["DEdeg"].values * u.deg)
    shift = 5 * np.log10(C.CF4_H0 / C.H0)
    for i in np.flatnonzero(near):
        d_mpc = max(cz[i], 300.0) / C.H0
        r_arcmin = max(np.degrees(C.CF4_MATCH_KPC / 1000.0 / d_mpc) * 60, C.CF4_MATCH_MIN_ARCMIN)
        sn_c = SkyCoord(df["ra"].iat[i] * u.deg, df["dec"].iat[i] * u.deg)
        sep = sn_c.separation(gal_c).arcmin
        ok = (sep < r_arcmin) & (np.abs(cf4["Vcmb"].to_numpy() - cz[i]) < C.CF4_MATCH_DV_KMS)
        if not ok.any():
            continue
        j = np.flatnonzero(ok)[np.argmin(sep[ok])]
        row = cf4.iloc[j]
        idx = df.index[i]
        out.at[idx, "pgc"] = int(row["PGC"])
        out.at[idx, "sep_kpc"] = np.radians(sep[j] / 60) * 10 ** ((row["dm_adopted"] - 25) / 5) * 1e3
        out.at[idx, "mu_h70"] = row["dm_adopted"] + shift
        out.at[idx, "e_mu"] = row["e_dm_adopted"]
        out.at[idx, "source"] = row["dm_source"]
    return out
