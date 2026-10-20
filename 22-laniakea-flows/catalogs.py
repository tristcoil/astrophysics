"""Read the catalogues, and the coordinate and reference-frame conversions they need.

Directions are handled as unit vectors. Galactic Cartesian axes point to the
Galactic centre (x), the direction of Galactic rotation (y, l = 90 deg) and the
north Galactic pole (z). Supergalactic axes follow de Vaucouleurs: SGZ toward
the supergalactic pole (l, b) = (47.37, 6.32) deg and SGX toward
(l, b) = (137.37, 0) deg, where the supergalactic plane crosses the Galactic
plane. The rotation between the two is computed once from astropy.

Redshifts come in three frames. *Heliocentric* is what a spectrograph
measures, corrected to the Solar System's barycentre. The *CMB frame* removes
the Sun's 369.82 km/s motion through the cosmic microwave background; it is
the frame of the expanding universe itself. The *Local Group frame* removes the
Sun's 299 km/s motion relative to our group of galaxies, so nearby galaxies,
which share most of our group's motion, look closest to a pure Hubble flow.
"""

from __future__ import annotations

import functools
import warnings

import numpy as np
import pandas as pd

import config
from config import C_KMS, Motion

warnings.filterwarnings("ignore", module="astropy")


# --------------------------------------------------------------------------
# Directions
# --------------------------------------------------------------------------

def unit_vector(lon_deg, lat_deg) -> np.ndarray:
    """Unit vectors (N, 3) for longitudes and latitudes in degrees."""
    lon, lat = np.radians(np.asarray(lon_deg, float)), np.radians(np.asarray(lat_deg, float))
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def lonlat(vec) -> tuple[np.ndarray, np.ndarray]:
    """Longitude in [0, 360) and latitude, in degrees, of vectors (..., 3)."""
    v = np.asarray(vec, float)
    lon = np.degrees(np.arctan2(v[..., 1], v[..., 0])) % 360.0
    lat = np.degrees(np.arctan2(v[..., 2], np.hypot(v[..., 0], v[..., 1])))
    return lon, lat


def angle_between(a, b) -> np.ndarray:
    """Angle in degrees between vectors (broadcasting)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    cross = np.linalg.norm(np.cross(a, b), axis=-1)
    return np.degrees(np.arctan2(cross, np.sum(a * b, axis=-1)))


@functools.lru_cache(maxsize=1)
def gal_to_sg_matrix() -> np.ndarray:
    """Rotation matrix R with sg_vector = R @ gal_vector."""
    from astropy.coordinates import Galactic, SkyCoord, Supergalactic
    import astropy.units as u

    basis = SkyCoord(l=[0.0, 90.0, 0.0] * u.deg, b=[0.0, 0.0, 90.0] * u.deg, frame=Galactic)
    sg = basis.transform_to(Supergalactic())
    # Column j of R is the supergalactic vector of Galactic basis vector j.
    return unit_vector(sg.sgl.deg, sg.sgb.deg).T


def gal_to_sg(vec) -> np.ndarray:
    return np.asarray(vec, float) @ gal_to_sg_matrix().T


def sg_to_gal(vec) -> np.ndarray:
    return np.asarray(vec, float) @ gal_to_sg_matrix()


def equatorial_to_galactic(ra_deg, dec_deg) -> tuple[np.ndarray, np.ndarray]:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c = SkyCoord(ra=np.atleast_1d(ra_deg) * u.deg, dec=np.atleast_1d(dec_deg) * u.deg, frame="icrs").galactic
    return c.l.deg, c.b.deg


# --------------------------------------------------------------------------
# Reference frames
# --------------------------------------------------------------------------

def motion_vector(m: Motion) -> np.ndarray:
    """A Motion as a Galactic Cartesian velocity vector, km/s."""
    return m.speed * unit_vector(m.l, m.b)


def boost_redshift(z, nhat_gal, v_observer: np.ndarray):
    """Redshift seen by an observer moving with ``v_observer`` (km/s, Galactic
    Cartesian), given the redshift ``z`` seen by one at rest, for sources in
    directions ``nhat_gal`` (N, 3). Exact special-relativistic Doppler factor:
    an observer moving toward a source sees it blueshifted by gamma (1 - beta cos theta).
    """
    beta = np.asarray(v_observer) / C_KMS
    gamma = 1.0 / np.sqrt(1.0 - beta @ beta)
    cos_term = np.asarray(nhat_gal) @ beta
    return (1.0 + np.asarray(z)) * gamma * (1.0 - cos_term) - 1.0


def helio_to_frame(z_helio, nhat_gal, sun_motion: Motion):
    """Heliocentric redshift -> redshift in a frame in which the Sun moves with
    ``sun_motion`` (e.g. the CMB frame for SUN_CMB, the Local Group frame for SUN_LG)."""
    v = motion_vector(sun_motion)
    beta = v / C_KMS
    gamma = 1.0 / np.sqrt(1.0 - beta @ beta)
    return (1.0 + np.asarray(z_helio)) / (gamma * (1.0 - np.asarray(nhat_gal) @ beta)) - 1.0


def frame_to_helio(z_frame, nhat_gal, sun_motion: Motion):
    return boost_redshift(z_frame, nhat_gal, motion_vector(sun_motion))


def local_group_motion() -> dict:
    """The Local Group's velocity through the CMB: the Sun's CMB velocity minus
    its velocity relative to the Local Group (non-relativistic, as Planck)."""
    v = motion_vector(config.SUN_CMB) - motion_vector(config.SUN_LG)
    l, b = lonlat(v)
    return {"vector_gal": v, "speed": float(np.linalg.norm(v)), "l": float(l), "b": float(b)}


# --------------------------------------------------------------------------
# Cosmicflows-4
# --------------------------------------------------------------------------

def _read_cds(folder: str, table: str) -> pd.DataFrame:
    from astropy.io import ascii

    path = config.DATA_DIR / folder
    tab = ascii.read(str(path / table), readme=str(path / "ReadMe"), format="cds")
    df = tab.to_pandas()
    return df


@functools.lru_cache(maxsize=None)
def cf4_galaxies() -> pd.DataFrame:
    """CF4 Table 2: 55,877 galaxies with distance moduli (combined and per method).

    Adds Galactic and supergalactic unit vectors and the CMB-frame redshift.
    """
    df = _read_cds("cf4", "table2.dat")
    df["z_cmb"] = df["Vcmb"] / C_KMS
    df[["gx", "gy", "gz"]] = unit_vector(df["GLON"], df["GLAT"])
    df[["sgx", "sgy", "sgz"]] = unit_vector(df["SGL"], df["SGB"])
    df["z_helio"] = frame_to_helio(df["z_cmb"].to_numpy(), df[["gx", "gy", "gz"]].to_numpy(), config.SUN_CMB)
    return df


@functools.lru_cache(maxsize=None)
def cf4_groups() -> pd.DataFrame:
    """CF4 Tables 3 and 4 joined: 38,053 groups (most are single galaxies).

    Table 4 carries the zero-pointed group distance modulus DMzp, its error,
    the luminosity distance, velocities in the heliocentric, Local Sheet and
    CMB frames, and CF4's own peculiar velocities. Table 3 adds the per-method
    group moduli and how many members each method measured.
    """
    t4 = _read_cds("cf4", "table4.dat")
    t3 = _read_cds("cf4", "table3.dat")
    keep3 = ["1PGC", "DMav", "e_DMav", "o_DMcal", "DMcal", "o_DMsnIa", "DMsnIa", "e_DMsnIa", "o_DMfp", "DMfp",
             "e_DMfp", "o_DMtf", "DMtf", "e_DMtf", "o_DMsbfo", "DMsbfo", "e_DMsbfo", "o_DMsbfi", "DMsbfi",
             "e_DMsbfi", "DMsnII", "e_DMsnII"]
    df = t4.merge(t3[keep3], on="1PGC", how="left", validate="one_to_one")
    df["z_cmb"] = df["V3k"] / C_KMS
    df[["gx", "gy", "gz"]] = unit_vector(df["GLON"], df["GLAT"])
    df[["sgx", "sgy", "sgz"]] = unit_vector(df["SGL"], df["SGB"])
    df["method"] = dominant_method(df)
    return df


METHOD_ORDER = ["calibrator", "SN Ia", "SBF", "SN II", "TF", "FP"]


def dominant_method(df: pd.DataFrame) -> pd.Series:
    """Label each group by its most precise distance method.

    Groups are often measured several ways; the label is the method with the
    smallest error among those present (Cepheid/TRGB/maser calibrators first).
    """
    errs = pd.DataFrame({
        "SN Ia": df["e_DMsnIa"], "TF": df["e_DMtf"], "FP": df["e_DMfp"],
        "SBF": df[["e_DMsbfo", "e_DMsbfi"]].min(axis=1), "SN II": df["e_DMsnII"],
    })
    has = errs.notna().any(axis=1)
    label = pd.Series("other", index=df.index)
    label[has] = errs[has].idxmin(axis=1, skipna=True)
    label = label.where(df["o_DMcal"].isna(), "calibrator")
    return label


# --------------------------------------------------------------------------
# 2MRS and its groups
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def twomrs() -> pd.DataFrame:
    """The 2MRS catalogue with redshifts in the three frames.

    Keeps every galaxy (with or without a redshift); ``has_cz`` marks those
    with one. Velocities: ``cz_helio`` as published, ``cz_cmb`` and ``cz_lg``
    after removing the Sun's motion through the CMB / relative to the Local Group.
    """
    df = _read_cds("2mrs", "table3.dat")
    df = df.rename(columns={"E(B-V)": "EBV"})
    n = unit_vector(df["GLON"], df["GLAT"])
    df[["gx", "gy", "gz"]] = n
    sg = gal_to_sg(n)
    df[["sgx", "sgy", "sgz"]] = sg
    df["SGL"], df["SGB"] = lonlat(sg)
    df["has_cz"] = df["cz"].notna()
    z = df["cz"].to_numpy(float) / C_KMS
    df["cz_helio"] = df["cz"].astype(float)
    df["cz_cmb"] = C_KMS * helio_to_frame(z, n, config.SUN_CMB)
    df["cz_lg"] = C_KMS * helio_to_frame(z, n, config.SUN_LG)
    df["in_zoa"] = in_zone_of_avoidance(df["GLON"].to_numpy(), df["GLAT"].to_numpy())
    return df


def in_zone_of_avoidance(l_deg, b_deg) -> np.ndarray:
    """The 2MRS mask: |b| < 5 deg, widened to |b| < 8 deg toward the bulge (|l| < 30 deg)."""
    l = (np.asarray(l_deg) + 180.0) % 360.0 - 180.0
    bcut = np.where(np.abs(l) < config.ZOA_L_BULGE, config.ZOA_B_BULGE, config.ZOA_B)
    return np.abs(np.asarray(b_deg)) < bcut


@functools.lru_cache(maxsize=None)
def tully_groups() -> pd.DataFrame:
    """Tully (2015) group memberships of the 2MRS galaxies."""
    df = _read_cds("tully2015", "table5.dat")
    return df


@functools.lru_cache(maxsize=None)
def twomrs_with_groups() -> pd.DataFrame:
    """2MRS galaxies matched to Tully's (2015) groups by position and velocity.

    The two catalogues share the same galaxies but no common identifier (2MASS
    ID vs PGC number), so each 2MRS galaxy is matched to the nearest Tully
    galaxy within 6 arcsec whose heliocentric velocity agrees to 150 km/s.
    Matched galaxies get their group ('nest') number and the group's mean
    velocity; unmatched ones are treated as groups of one.
    """
    from scipy.spatial import cKDTree

    g = twomrs().copy()
    t = tully_groups()
    tn = unit_vector(t["GLON"], t["GLAT"])
    tree = cKDTree(tn)
    chord = 2 * np.sin(np.radians(6.0 / 3600) / 2)
    dist, idx = tree.query(g[["gx", "gy", "gz"]].to_numpy(), k=4, distance_upper_bound=chord)
    nest = np.full(len(g), -1, dtype=np.int64)
    for k in range(idx.shape[1]):
        ok = (nest < 0) & np.isfinite(dist[:, k])
        j = np.where(ok, idx[:, k], 0)
        dv = np.abs(g["cz_helio"].to_numpy() - t["HV"].to_numpy()[j])
        good = ok & (dv < 150.0)
        nest[good] = t["Nest"].to_numpy()[j[good]]
    g["nest"] = nest
    return g
