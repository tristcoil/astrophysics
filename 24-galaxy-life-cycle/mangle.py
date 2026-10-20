"""The survey's own footprint: SDSS tiling polygons in mangle format.

The legacy survey placed its spectroscopic plates ("tiles", circles 1.49 deg
in radius) over a region bounded by the imaging. The overlaps of tiles and
boundaries cut the sky into *sectors*; the file ``sdss_tiling_polygons.fits``
(SDSS SAS, tiling/final) describes each piece as a mangle polygon: the
intersection of a few spherical caps. A cap is an axis c and a number cm;
a point x (unit vector) is inside it when 1 - x.c < cm, or, for cm < 0,
outside the complementary cap: 1 - x.c > -cm (Hamilton & Tegmark 2004).

With these polygons we know exactly which sky the survey covered and how big
it is, rather than inferring it from where galaxies happen to be.
"""

from __future__ import annotations

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree

import cosmo_physics as cp

TILE_RADIUS_DEG = 1.49


def load_polygons(path) -> dict:
    d = fits.open(path)[1].data
    return {"xcaps": np.array(d["XCAPS"], float), "cmcaps": np.array(d["CMCAPS"], float),
            "ncaps": np.array(d["NCAPS"], int), "use_caps": np.array(d["USE_CAPS"], np.int64),
            "str": np.array(d["STR"], float), "sector": np.array(d["SECTOR"], int),
            "weight": np.array(d["WEIGHT"], float)}


def load_sectors(path) -> dict:
    d = fits.open(path)[1].data
    return {"sector": np.array(d["SECTOR"], int), "ntiles": np.array(d["NTILES"], int),
            "area": np.array(d["AREA"], float)}


def _bounding_cap(poly: dict, k: int):
    """The smallest positive cap of polygon k: the polygon lies inside it."""
    n = poly["ncaps"][k]
    use = poly["use_caps"][k]
    best = None
    for j in range(n):
        if not (use >> j) & 1:
            continue
        cm = poly["cmcaps"][k, j]
        if 0 < cm < 1 and (best is None or cm < best[1]):
            best = (poly["xcaps"][k, j], cm)
    return best


def _inside(poly: dict, k: int, x: np.ndarray) -> np.ndarray:
    ok = np.ones(len(x), bool)
    use = poly["use_caps"][k]
    for j in range(poly["ncaps"][k]):
        if not (use >> j) & 1:
            continue
        c, cm = poly["xcaps"][k, j], poly["cmcaps"][k, j]
        d = 1.0 - x @ c
        ok &= (d < cm) if cm >= 0 else (d > -cm)
    return ok


def which_polygon(poly: dict, ra, dec) -> np.ndarray:
    """Index of the polygon containing each point, or -1 if none.

    Each tiled polygon lies inside at least one tile circle (a positive cap of
    radius 1.49 deg); candidates are found with a KD-tree inside that circle
    and then tested against every cap. Polygons without a bounding cap
    (untiled sky, sector 0) are skipped: points there count as outside.
    """
    x = cp.unit_vectors(ra, dec)
    tree = cKDTree(x)
    out = np.full(len(x), -1, dtype=np.int64)
    for k in range(len(poly["str"])):
        bc = _bounding_cap(poly, k)
        if bc is None:
            continue
        axis, cm = bc
        chord = np.sqrt(2.0 * cm) * 1.0001  # |x - c|^2 = 2 (1 - x.c)
        cand = np.asarray(tree.query_ball_point(axis, chord), dtype=np.int64)
        if cand.size == 0:
            continue
        cand = cand[out[cand] < 0]
        hit = cand[_inside(poly, k, x[cand])]
        out[hit] = k
    return out


def tiled_mask(poly: dict, nside: int) -> np.ndarray:
    """HEALPix map: True where the pixel centre lies in a tiled polygon."""
    import healpy as hp

    ra, dec = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    # only test pixels near the survey (northern sky and the equatorial stripes)
    k = which_polygon(poly, ra, dec)
    return k >= 0
