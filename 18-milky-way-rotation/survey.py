"""Reading the HI4PI Galactic-plane cubes as one 360-degree survey.

The 18 cubes share one pixel grid: 5' (1/12 deg) pixels at exact multiples of
1/12 deg in longitude and latitude, and identical velocity channels. Each cube
covers 22.2 deg of longitude around its centre, overlapping its neighbours by
about a degree, so a mosaic takes from each cube only the 20 degrees nearest
its centre; no resampling is ever needed.

Products that take a pass over all 2.8 GB (the longitude-velocity diagrams,
the column-density map) are computed once and cached in ``data/products/``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from astropy.io import fits

import config as cfg
from fetch_data import tile_path

PIX_PER_DEG = 12
N_LON = 360 * PIX_PER_DEG
TILE_CENTRES = np.arange(10, 360, 20)


@dataclass
class Tile:
    name: str
    path: Path
    l_centre: float
    data: np.ndarray          # (v, b, l) memmap, big-endian float32
    lon: np.ndarray           # deg, per column (may be negative or > 360)
    lat: np.ndarray           # deg, per row
    vel: np.ndarray           # km/s, per channel


def open_tile(name: str) -> Tile:
    path = tile_path(name)
    hdul = fits.open(path, memmap=True)
    h = hdul[0].header
    data = hdul[0].data
    nv, ny, nx = data.shape
    lon = h["CRVAL1"] + (np.arange(1, nx + 1) - h["CRPIX1"]) * h["CDELT1"]
    lat = h["CRVAL2"] + (np.arange(1, ny + 1) - h["CRPIX2"]) * h["CDELT2"]
    vel = (h["CRVAL3"] + (np.arange(1, nv + 1) - h["CRPIX3"]) * h["CDELT3"]) / 1000.0
    return Tile(name, path, float(h["CRVAL1"]), data, lon, lat, vel)


def tiles(names=None):
    for name in (cfg.HI4PI_ROW_E if names is None else names):
        yield open_tile(name)


#: Highest latitude a row-E cube reaches (b = +11.083 deg: 266 rows, b = 0 at
#: row 133). A row-F cube contributes only the rows above it.
ROW_E_TOP = 133 / PIX_PER_DEG


def own_rows(tile: Tile) -> np.ndarray:
    """Rows this cube contributes when rows E and F are combined."""
    if abs(tile.lat.mean()) < 1:  # row E, centred on the plane
        return np.arange(tile.lat.size)
    return np.flatnonzero(tile.lat > ROW_E_TOP + 0.5 / PIX_PER_DEG)


def own_columns(tile: Tile) -> tuple[np.ndarray, np.ndarray]:
    """Columns this cube contributes to the mosaic, and their global index.

    Global longitude index k means l = k / 12 deg, k = 0 .. 4319.
    """
    k = np.rint(tile.lon * PIX_PER_DEG).astype(int)
    lc = int(round(tile.l_centre * PIX_PER_DEG))
    half = 10 * PIX_PER_DEG
    cols = np.flatnonzero((k >= lc - half) & (k < lc + half))
    return cols, k[cols] % N_LON


def rows_for(tile: Tile, b_lo: float, b_hi: float) -> np.ndarray:
    eps = 1e-6
    return np.flatnonzero((tile.lat >= b_lo - eps) & (tile.lat <= b_hi + eps))


def velocity_axis() -> np.ndarray:
    return open_tile(cfg.HI4PI_ROW_E[0]).vel


def longitudes() -> np.ndarray:
    return np.arange(N_LON) / PIX_PER_DEG


# --------------------------------------------------------------------------
# Cached products
# --------------------------------------------------------------------------

def _cached(name: str, build, refresh: bool = False) -> dict:
    path = cfg.PRODUCTS_DIR / f"{name}.npz"
    if path.exists() and not refresh:
        with np.load(path) as z:
            return {k: z[k] for k in z.files}
    out = build()
    np.savez_compressed(path, **out)
    return out


def lv_diagram(b_lo: float, b_hi: float, *, stat: str = "mean",
               refresh: bool = False) -> dict:
    """Brightness temperature over b_lo <= b <= b_hi at every (l, v), combined
    across latitude by ``stat`` ("mean" or "median").

    Returns ``l`` (4320,), ``v`` (543,), ``T`` (4320, 543) in K.
    """
    tag = f"lv_b{b_lo:+.3f}_{b_hi:+.3f}_{stat}".replace(".", "p")
    combine = {"mean": np.nanmean, "median": np.nanmedian}[stat]

    def build():
        v = velocity_axis()
        T = np.full((N_LON, v.size), np.nan, dtype=np.float32)
        for t in tiles():
            if not np.allclose(t.vel, v):
                raise RuntimeError(f"{t.name}: velocity axis differs")
            cols, k = own_columns(t)
            rows = rows_for(t, b_lo, b_hi)
            block = np.asarray(t.data[:, rows[0]:rows[-1] + 1, :], dtype=np.float32)
            T[k] = combine(block[:, :, cols], axis=1).T
        return {"l": longitudes(), "v": v, "T": T, "b_range": np.array([b_lo, b_hi])}

    return _cached(tag, build, refresh)


def column_density_map(*, refresh: bool = False) -> dict:
    """N_HI(l, b) over the whole strip |b| <= 11 deg, integrating all kept
    channels (|v| <= 350 km/s), optically thin (a lower limit where the gas is
    opaque). Returns ``l``, ``b`` and ``NHI`` (b, l) in cm^-2."""

    def build():
        first = open_tile(cfg.HI4PI_ROW_E[0])
        dv = abs(first.vel[1] - first.vel[0])
        nhi = np.full((first.lat.size, N_LON), np.nan, dtype=np.float32)
        for t in tiles():
            cols, k = own_columns(t)
            acc = np.zeros((t.lat.size, cols.size), dtype=np.float64)
            for i0 in range(0, t.vel.size, 64):
                chunk = np.asarray(t.data[i0:i0 + 64][:, :, cols], dtype=np.float64)
                acc += np.nansum(chunk, axis=0)
            nhi[:, k] = (cfg.C_HI * dv * acc).astype(np.float32)
        return {"l": longitudes(), "b": first.lat, "NHI": nhi}

    return _cached("column_density", build, refresh)
