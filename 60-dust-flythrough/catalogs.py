"""Loaders for the catalogues the project compares the dust map with.

All positions are converted to the heliocentric Galactic Cartesian frame of
the map (x to the Galactic Centre, y to l = 90 deg, z to the North Galactic
Pole, in pc).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from astropy.io import ascii, fits

import config as C
from dustmap import lbd_to_xyz

#: Hunt & Reffert (2023) "highly reliable" clusters: astrometric S/N (CST) above
#: 5 sigma and a median CMD class above 0.5 (their sect. 6).
HR23_CST_MIN = 5.0
HR23_CMD_MIN = 0.5


# ---------------------------------------------------------------------------
# Molecular-cloud distances (Zucker et al. 2020; Alves et al. 2020)
# ---------------------------------------------------------------------------

@lru_cache
def zucker2020() -> pd.DataFrame:
    """Zucker et al. (2020) Table A.1: 326 sightlines toward local molecular clouds.

    ``d16/d50/d84`` are the cloud distance's percentiles (pc); we recompute
    x, y, z from (l, b, d50) and keep the published ones as ``x_pub`` etc.
    """
    df = pd.read_csv(C.CATALOG_DIR / "Handbook_Distances_Zucker2020.dat", sep=r"\s+")
    df = df.rename(columns={"x": "x_pub", "y": "y_pub", "z": "z_pub", "New_Sightline?": "new_sightline"})
    df["x"], df["y"], df["z"] = lbd_to_xyz(df["l"], df["b"], df["d50"])
    df["complex"] = df["name"].map(complex_of)
    return df


@lru_cache
def bridging() -> pd.DataFrame:
    """Alves et al. (2020): 54 sightlines along the 'bridges' between Wave clouds."""
    df = pd.read_csv(C.CATALOG_DIR / "Bridging_Sightlines_Alves_2020.dat", sep=r"\s+")
    df = df.rename(columns={"x": "x_pub", "y": "y_pub", "z": "z_pub"})
    df["x"], df["y"], df["z"] = lbd_to_xyz(df["l"], df["b"], df["d50"])
    return df


#: Zucker et al.'s sightline names -> the complex we label on the maps.
_COMPLEX = {
    "Aquila_Rift": "Aquila Rift", "Aquila_S": "Aquila South", "California": "California",
    "CMa_OB1": "CMa OB1", "Cam": "Camelopardalis", "Cepheus": "Cepheus", "Cepheus_Far": "Cepheus",
    "Cepheus_L1251": "Cepheus", "Cepheus_OB3": "Cepheus", "Cepheus_South": "Cepheus",
    "Chamaeleon": "Chamaeleon", "Cham_I": "Chamaeleon", "Cham_II": "Chamaeleon",
    "Cham_III": "Chamaeleon", "Musca": "Musca", "CrA": "Corona Australis",
    "Cygnus_X": "Cygnus X", "Cygnus": "Cygnus X", "Lupus": "Lupus", "Ophiuchus": "Ophiuchus",
    "Orion_A": "Orion A", "Orion_B": "Orion B", "Orion": "Orion", "Lambda_Ori": "lambda Ori",
    "Mon_R2": "Mon R2", "Mon_OB1": "Mon OB1", "Rosette": "Rosette", "NAN": "North America",
    "Perseus": "Perseus", "Pipe": "Pipe", "Serpens": "Serpens", "Taurus": "Taurus",
    "Vela": "Vela", "Polaris": "Polaris Flare", "Ursa_Major": "Ursa Major",
}


def complex_of(name: str) -> str:
    if name in _COMPLEX:
        return _COMPLEX[name]
    for key, value in _COMPLEX.items():
        if name.startswith(key):
            return value
    return name.replace("_", " ")


def cloud_labels(min_sightlines: int = 3) -> pd.DataFrame:
    """Mean (x, y, z) of each named complex, for labelling the maps."""
    df = zucker2020()
    g = df.groupby("complex")
    out = g[["x", "y", "z", "d50"]].median()
    out["n"] = g.size()
    return out[out["n"] >= min_sightlines].reset_index()


# ---------------------------------------------------------------------------
# Young star clusters (Hunt & Reffert 2023)
# ---------------------------------------------------------------------------

@lru_cache
def hunt_reffert() -> pd.DataFrame:
    """Hunt & Reffert (2023) clusters with heliocentric x, y, z (from l, b, dist50)."""
    folder = C.CATALOG_DIR / "hunt_reffert_2023"
    table = ascii.read(folder / "clusters.dat", readme=folder / "ReadMe", format="cds")
    keep = ["Name", "Type", "CST", "N", "GLON", "GLAT", "r50pc", "pmRA", "e_pmRA", "pmDE", "e_pmDE",
            "Plx", "e_Plx", "dist16", "dist50", "dist84", "RV", "e_RV", "n_RV", "CMDCl50",
            "logAge16", "logAge50", "logAge84", "AV16", "AV50", "AV84", "diffAV50", "RAdeg", "DEdeg"]
    df = table[keep].to_pandas()
    for col in ("RV", "e_RV"):
        df[col] = df[col].astype(float)
    df["x"], df["y"], df["z"] = lbd_to_xyz(df["GLON"], df["GLAT"], df["dist50"])
    df["age_myr"] = 10 ** (df["logAge50"] - 6)
    df["reliable"] = (df["CST"] > HR23_CST_MIN) & (df["CMDCl50"] > HR23_CMD_MIN)
    return df


# ---------------------------------------------------------------------------
# The Local Bubble (O'Neill et al. 2024; Zucker et al. 2022)
# ---------------------------------------------------------------------------

@lru_cache
def oneill_shell() -> pd.DataFrame:
    """O'Neill et al. (2024) shell properties, one row per Nside-256 sightline."""
    path = C.CATALOG_DIR / "ONeill2024_LocalBubble_ShellProperties_A0.5.fits"
    with fits.open(path) as hdul:
        data = hdul[1].data
        df = pd.DataFrame({name: np.asarray(data[name]).astype(np.float64)
                           if data[name].dtype.kind == "f" else np.asarray(data[name])
                           for name in data.names})
    return df


def zucker2022_tracebacks() -> dict[str, pd.DataFrame]:
    """Zucker et al. (2022) cluster tracebacks: t (Myr, <= 0) and x, y, z (pc, LSR frame)."""
    out = {}
    for path in sorted((C.CATALOG_DIR / "zucker2022_tracebacks").glob("*_traceback.dat")):
        df = pd.read_csv(path, sep=r"\s+")
        df.columns = ["t", "x", "y", "z"]
        out[path.name.removesuffix("_traceback.dat")] = df
    return out


# ---------------------------------------------------------------------------
# The second 3D map (Vergely, Lallement & Cox 2022)
# ---------------------------------------------------------------------------

@lru_cache
def vergely() -> dict:
    """The 10 pc-resolution cube: density (A0 at 550 nm, mag/pc) as (nx, ny, nz), axes in pc.

    The header has no WCS, only STEP = 5 pc and the Sun at SUN_POS = (300.5,
    300.5, 80.5). We read that as the centre of pixel (300, 300, 80) counted
    from 0, i.e. pixel centres at (i - 300) x 5 pc: a box of +-1500 x +-1500 x
    +-400 pc centred on the Sun. (Reading SUN_POS as a 1-based pixel coordinate
    instead would move everything by 2.5 pc; ``compare.vergely_offset``
    cross-correlates with the Edenhofer map to check.)
    """
    with fits.open(C.VERGELY_DIR / C.VERGELY_FILE) as hdul:
        hdr = hdul[0].header
        cube = np.asarray(hdul[0].data, dtype=np.float32)       # FITS order (z, y, x)
    step = float(hdr["STEP"])
    sun = [hdr["SUN_POSX"], hdr["SUN_POSY"], hdr["SUN_POSZ"]]
    axes = [(np.arange(n) + 0.5 - s) * step for n, s in zip(cube.shape[::-1], sun)]
    return {"cube": np.transpose(cube), "x": axes[0], "y": axes[1], "z": axes[2], "header": hdr,
            "unit": hdr["UNIT"]}
