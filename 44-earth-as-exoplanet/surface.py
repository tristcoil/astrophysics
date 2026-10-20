"""What part of the dot's light comes from land, ocean and the polar caps.

For a Lambert surface, a patch of area dA contributes to the planet's total
reflected flux in proportion to mu0 * mu * dA (mu0: cosine of the Sun's zenith
angle at the patch; mu: cosine of the angle to the observer). Weighting a
land-sea mask this way gives, for each EPIC image, the share of the dot's
light that a uniformly reflecting Earth would get from land. Real land is
darker or brighter than ocean depending on wavelength, which is exactly what
the comparison with the measured colours tests.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import LANDMASK_PATH

GRID_DEG = 1.0  # the 0.25-degree mask averaged to 1 degree
POLAR_LAT = 60.0  # "polar" = poleward of 60 deg, as in Yang et al. (2018)


def land_grid(step_deg: float = GRID_DEG):
    """Land fraction on a regular lat/lon grid (cell centres), from the NOAA mask."""
    import xarray as xr

    with xr.open_dataset(LANDMASK_PATH) as ds:
        ocean = ds["lsmask"].isel(time=0).values.astype(float)  # 1 = ocean
        lat = ds["lat"].values
    land = 1.0 - ocean
    k = int(round(step_deg / 0.25))
    ny, nx = land.shape
    land = land.reshape(ny // k, k, nx // k, k).mean(axis=(1, 3))
    lat_c = lat.reshape(-1, k).mean(axis=1)
    lon_c = (np.arange(nx // k) + 0.5) * step_deg  # 0..360 east
    return lat_c, lon_c, land


def _unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def disc_fractions(obs_lat, obs_lon, sun_lat, sun_lon, chunk: int = 256) -> pd.DataFrame:
    """Per image: land share and polar share of a Lambert Earth's reflected light,
    plus the land share of the sunlit visible area (Yang et al.'s pixel-count
    definition)."""
    lat_c, lon_c, land = land_grid()
    LAT, LON = np.meshgrid(lat_c, lon_c, indexing="ij")
    normals = _unit(LAT, LON).reshape(-1, 3).astype(np.float32)
    area = np.cos(np.radians(LAT)).ravel().astype(np.float32)
    landv = land.ravel().astype(np.float32)
    polar = (np.abs(LAT) > POLAR_LAT).ravel().astype(np.float32)
    o = _unit(np.asarray(obs_lat), np.asarray(obs_lon)).astype(np.float32)
    s = _unit(np.asarray(sun_lat), np.asarray(sun_lon)).astype(np.float32)
    out = {"land_flux": [], "polar_flux": [], "land_area": []}
    for i in range(0, len(o), chunk):
        mu = np.clip(o[i:i + chunk] @ normals.T, 0, None)
        mu0 = np.clip(s[i:i + chunk] @ normals.T, 0, None)
        w = mu * mu0 * area
        wa = mu * (mu0 > 0) * area  # projected sunlit area: what the pixels count
        out["land_flux"].append((w @ landv) / w.sum(axis=1))
        out["polar_flux"].append((w @ polar) / w.sum(axis=1))
        out["land_area"].append((wa @ landv) / wa.sum(axis=1))
    return pd.DataFrame({k: np.concatenate(v) for k, v in out.items()})


def add_to(df: pd.DataFrame) -> pd.DataFrame:
    """Attach land/polar shares to a light-curve table (uses its geometry columns)."""
    fr = disc_fractions(df["subobs_lat_calc"].to_numpy(), df["subobs_lon_calc"].to_numpy(),
                        df["subsolar_lat"].to_numpy(), df["subsolar_lon"].to_numpy())
    fr.index = df.index
    return pd.concat([df, fr], axis=1)
