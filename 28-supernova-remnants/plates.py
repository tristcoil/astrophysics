"""Register every image of a remnant onto one Gaia-defined grid, and cache it.

Each DSS scan is one image. A Pan-STARRS cutout can straddle up to four
"skycells" (the Crab sits on such a corner), each a separate stack with its
own input exposures and therefore its own date; we solve each piece
separately, then mosaic them, keeping a map of which date each pixel carries.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

import astrometry as A
import config

warnings.filterwarnings("ignore", category=FITSFixedWarning)

#: Minimum number of stars for a usable plate solution (10 coefficients per axis).
MIN_STARS = 40


@dataclass
class Registered:
    """One image of one remnant on the common grid."""

    plate: config.Plate
    image: np.ndarray  # NaN outside the image's footprint
    epoch_map: np.ndarray  # Julian year of each pixel (NaN outside)
    epoch: float  # median over the footprint
    fwhm_arcsec: float  # measured on the grid from faint Gaia stars
    solutions: list  # one dict per input file: stars used, rms, pixel size

    @property
    def epoch_span(self) -> tuple[float, float]:
        e = self.epoch_map[np.isfinite(self.epoch_map)]
        return float(e.min()), float(e.max())


def grid_for(target: config.Target, *, step: float | None = None, size_arcmin: float | None = None) -> A.Grid:
    step = step or (config.CRAB_GRID_ARCSEC if target.key == "crab" else 0.5)
    size = size_arcmin or (config.CRAB_GRID_SIZE_ARCMIN if target.key == "crab" else target.box_arcmin - 1.0)
    return A.Grid(target.ra, target.dec, step, int(round(size * 60 / step)))


def _input_files(plate: config.Plate) -> list:
    if plate.survey.startswith("ps1-"):
        return sorted(config.PLATE_DIR.glob(f"{plate.key}_*.fits"))
    return [config.PLATE_DIR / f"{plate.key}.fits"]


def _star_radius_px(pixel_arcsec: float) -> float:
    """Gaussian weighting radius for centroids: ~1.2" (seeing sigma on plates)."""
    return max(1.0, 1.2 / pixel_arcsec)


def _solve_file(path, plate, gaia, target):
    header = fits.getheader(path)
    img = fits.getdata(path).astype(float)
    if plate.survey.startswith("ps1-"):
        epoch, first, last = A.ps1_epoch(header)
    else:
        epoch = A.dss_epoch(header)
        first = last = epoch
    wcs = WCS(header)
    pix = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))) * 3600.0)
    sol = A.solve_plate(f"{plate.key}:{path.name}", img, wcs, epoch, gaia, target.ra, target.dec,
                        star_radius_px=_star_radius_px(pix))
    info = {"file": path.name, "epoch": epoch, "first_exposure": first, "last_exposure": last,
            "n_stars": sol.n_stars, "rms_arcsec": sol.rms_arcsec,
            "pixel_arcsec": sol.pixel_scale_arcsec()}
    if not plate.survey.startswith("ps1-"):
        info.update({k.lower(): str(header.get(k, "")).strip()
                     for k in ("SURVEY", "PLATEID", "EMULSION", "FILTER", "DATE-OBS", "EXPOSURE")})
    return img, sol, info


def _measure_fwhm(image: np.ndarray, grid: A.Grid, gaia: A.GaiaStars, epoch: float) -> float:
    ra, dec = gaia.at_epoch(epoch)
    sel = np.isfinite(ra) & (gaia.g > 15.0) & (gaia.g < 18.0)
    xi, eta = A.to_tangent(ra[sel], dec[sel], grid.ra0, grid.dec0)
    col, row = grid.to_index(xi, eta)
    half = int(np.ceil(5.0 / grid.step))  # 5" half-width stamps
    return A.stamp_fwhm(image, col, row, half=half) * grid.step


def register(plate: config.Plate, grid: A.Grid, gaia: A.GaiaStars, *, refresh: bool = False) -> Registered:
    target = config.TARGETS[plate.target]
    cache = config.DERIVED_DIR / f"reg_{plate.key}.npz"
    meta_path = cache.with_suffix(".json")
    if cache.exists() and meta_path.exists() and not refresh:
        z = np.load(cache)
        meta = json.loads(meta_path.read_text())
        if meta.get("grid") == [grid.ra0, grid.dec0, grid.step, grid.n]:
            return Registered(plate, z["image"], z["epoch_map"], meta["epoch"], meta["fwhm_arcsec"],
                              meta["solutions"])

    pieces = []
    for path in _input_files(plate):
        img, sol, info = _solve_file(path, plate, gaia, target)
        if sol.n_stars < MIN_STARS:
            continue
        on_grid = A.resample(img, sol, grid)
        pieces.append((np.isfinite(on_grid).sum(), on_grid, sol.epoch, info))
    if not pieces:
        raise RuntimeError(f"{plate.key}: no image piece had {MIN_STARS} usable stars")

    # Mosaic: start from the piece covering most of the grid; fill holes from
    # the others after scaling them to match in the overlap (stacks of
    # different skycells have slightly different backgrounds and depths).
    pieces.sort(key=lambda p: -p[0])
    image = pieces[0][1].copy()
    epoch_map = np.where(np.isfinite(image), pieces[0][2], np.nan)
    for _, extra, ep, _ in pieces[1:]:
        both = np.isfinite(image) & np.isfinite(extra)
        holes = ~np.isfinite(image) & np.isfinite(extra)
        if not holes.any():
            continue
        if both.sum() > 1000:
            a, b = np.polyfit(extra[both], image[both], 1)
            extra = a * extra + b
        image[holes] = extra[holes]
        epoch_map[holes] = ep

    epoch = float(np.nanmedian(epoch_map))
    fwhm = _measure_fwhm(image, grid, gaia, epoch)
    solutions = [p[3] for p in pieces]
    np.savez_compressed(cache, image=image.astype(np.float32), epoch_map=epoch_map.astype(np.float64))
    meta_path.write_text(json.dumps({"grid": [grid.ra0, grid.dec0, grid.step, grid.n], "epoch": epoch,
                                     "fwhm_arcsec": fwhm, "solutions": solutions}, indent=1))
    return Registered(plate, image, epoch_map, epoch, fwhm, solutions)


def register_all(target_key: str, *, refresh: bool = False, grid: A.Grid | None = None) -> dict:
    target = config.TARGETS[target_key]
    grid = grid or grid_for(target)
    gaia = A.read_gaia(target_key)
    return {p.key: register(p, grid, gaia, refresh=refresh) for p in config.PLATES if p.target == target_key}
