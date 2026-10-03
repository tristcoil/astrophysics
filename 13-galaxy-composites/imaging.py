"""Calibrated images on one grid, the catalogues that go with them, and measuring them.

The viewer's tiles each come on their own grid (a tangent-plane projection
centred on the tile). ``stitch`` puts them onto the galaxy's single target grid
with ``reproject``, which works out, for every output pixel, where it falls in
each input image and interpolates the surface brightness there. Because the
tiles overlap and come from the same underlying SDSS coadd, the result has no
seams; how faithful the stitching is gets checked against two galaxies that
were also fetched in one piece (``analysis.validate_stitching``).

Pixel values stay as the survey serves them -- nanomaggies per *native* SDSS
pixel (0.396") -- in ``Frame.native``; ``Frame.sb`` holds the same images as
surface brightness in nanomaggies per square arcsecond, the unit every other
module works in.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import pandas as pd
from astropy.coordinates import Angle
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from PIL import Image

from config import (
    BANDS,
    DATA_DIR,
    NATIVE_PIXSCALE,
    Galaxy,
    central_tile,
    target_wcs,
    tile_layout,
)

DERIVED_DIR = DATA_DIR / "derived"
DERIVED_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# Images
# --------------------------------------------------------------------------

def read_cutout(path) -> tuple[dict[str, np.ndarray], WCS]:
    """A viewer FITS cutout: (band -> image, celestial WCS)."""
    with fits.open(path) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float64)
        hdr = hdul[0].header
    bands = [str(hdr[f"BAND{k}"]).strip() for k in range(data.shape[0])]
    return {b: data[k] for k, b in enumerate(bands)}, WCS(hdr).celestial


def stitch(g: Galaxy) -> dict[str, np.ndarray]:
    """Reproject every tile onto its share of the target grid."""
    from reproject import reproject_interp

    tw = target_wcs(g)
    n = g.size_px
    out = {b: np.full((n, n), np.nan) for b in BANDS}
    for t in tile_layout(g):
        imgs, w = read_cutout(DATA_DIR / t.filename)
        sub = tw[t.y0:t.y1, t.x0:t.x1]
        shape = (t.y1 - t.y0, t.x1 - t.x0)
        for b in BANDS:
            arr, _ = reproject_interp((imgs[b], w), sub, shape_out=shape, order="bicubic")
            out[b][t.y0:t.y1, t.x0:t.x1] = arr
    bad = {b: int(np.isnan(v).sum()) for b, v in out.items()}
    if any(bad.values()):
        raise RuntimeError(f"{g.key}: stitched mosaic has unfilled pixels {bad}")
    return out


def _mosaic_path(g: Galaxy):
    return DERIVED_DIR / f"{g.key}_sdss_mosaic.fits"


def mosaic(g: Galaxy, *, rebuild: bool = False) -> dict[str, np.ndarray]:
    """The stitched g, r, i images (nanomaggies per native pixel), cached on disk."""
    path = _mosaic_path(g)
    if path.exists() and not rebuild:
        with fits.open(path) as hdul:
            return {b: np.asarray(hdul[b].data, dtype=np.float64) for b in BANDS}
    imgs = stitch(g)
    hdr = target_wcs(g).to_header()
    hdr["BUNIT"] = ("nanomaggy/native_pixel", "per 0.396 arcsec SDSS pixel")
    hdul = fits.HDUList([fits.PrimaryHDU()] +
                        [fits.ImageHDU(imgs[b].astype(np.float32), header=hdr, name=b) for b in BANDS])
    hdul.writeto(path, overwrite=True)
    return imgs


@dataclass
class Frame:
    """One galaxy's images on its target grid."""

    galaxy: Galaxy
    wcs: WCS
    #: nanomaggies per native (0.396") pixel, as served
    native: dict[str, np.ndarray]
    #: nanomaggies per square arcsecond
    sb: dict[str, np.ndarray]

    @property
    def pixscale(self) -> float:
        return self.galaxy.pixscale

    @property
    def shape(self) -> tuple[int, int]:
        return self.native[BANDS[0]].shape

    def world_to_pixel(self, ra, dec):
        return self.wcs.world_to_pixel_values(ra, dec)


def load_frame(g: Galaxy) -> Frame:
    native = mosaic(g)
    sb = {b: v / NATIVE_PIXSCALE**2 for b, v in native.items()}
    return Frame(g, target_wcs(g), native, sb)


def load_single_cutout(g: Galaxy) -> dict[str, np.ndarray]:
    """The whole frame fetched as one piece (only M81 and M82 have one)."""
    imgs, _ = read_cutout(DATA_DIR / f"cutouts/{g.key}_sdss_single.fits")
    return imgs


def load_tile(g: Galaxy):
    """The central tile, exactly as served: (images, WCS, viewer JPEG in FITS row order)."""
    t = central_tile(g)
    imgs, w = read_cutout(DATA_DIR / t.filename)
    jpg = np.asarray(Image.open(DATA_DIR / f"cutouts/{g.key}/tile_{t.iy}_{t.ix}_viewer.jpg").convert("RGB"))
    return imgs, w, jpg[::-1]


# --------------------------------------------------------------------------
# Catalogues
# --------------------------------------------------------------------------

def load_sdss_stars(g: Galaxy) -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / f"catalogs/{g.key}_sdss_stars.csv", comment="#")


def load_gaia(g: Galaxy) -> pd.DataFrame:
    return pd.read_csv(DATA_DIR / f"catalogs/{g.key}_gaia.csv")


@lru_cache(maxsize=None)
def load_rc3(key: str) -> pd.DataFrame:
    """RC3 entries in the frame, with sizes turned into arcminutes.

    D25 is the diameter of the 25 B-mag/arcsec^2 isophote (stored as log10 of
    tenths of an arcminute); R25 is log10 of the major/minor axis ratio; PA runs
    from north through east.
    """
    path = DATA_DIR / f"catalogs/{key}_rc3.tsv"
    lines = [ln for ln in path.read_text().splitlines() if ln and not ln.startswith("#")]
    header = lines[0].split("\t")
    rows = [ln.split("\t") for ln in lines[3:]]  # skip units and dashes
    df = pd.DataFrame(rows, columns=header)
    df = df.apply(lambda col: col.str.strip())
    out = pd.DataFrame({
        "name": df["name"].str.replace(r"\s+", " ", regex=True),
        "ra": [Angle(v, unit="hourangle").degree for v in df["RA2000"]],
        "dec": [Angle(v, unit="deg").degree for v in df["DE2000"]],
        "d25_arcmin": 0.1 * 10.0 ** pd.to_numeric(df["D25"], errors="coerce"),
        "axis_ratio": 10.0 ** -pd.to_numeric(df["R25"], errors="coerce"),
        "pa_deg": pd.to_numeric(df["PA"], errors="coerce"),
        "BT": pd.to_numeric(df["BT"], errors="coerce"),
        "B-VT": pd.to_numeric(df["B-VT"], errors="coerce"),
    })
    return out


def rc3_entry(g: Galaxy, name: str | None = None) -> pd.Series:
    df = load_rc3(g.key)
    target = (name or g.ngc).replace("NGC 0", "NGC ")
    hit = df[df["name"] == target]
    if hit.empty:
        # The Mice are listed as NGC 4676A and B; take the brighter half's row
        hit = df[df["name"].str.startswith(target)]
    if hit.empty:
        raise KeyError(f"{target} not in RC3 around {g.name}: {list(df['name'])}")
    return hit.iloc[0]


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def offsets_arcsec(frame: Frame, ra: float, dec: float) -> tuple[np.ndarray, np.ndarray]:
    """(east, north) offset of every pixel from a sky position, in arcsec.

    Uses the target grid's tangent plane directly: north is up and east is to
    the left, so over a few arcminutes the offsets are exact to far below a
    pixel.
    """
    x0, y0 = frame.world_to_pixel(ra, dec)
    ny, nx = frame.shape
    yy, xx = np.mgrid[:ny, :nx]
    return -(xx - x0) * frame.pixscale, (yy - y0) * frame.pixscale


def elliptical_radius(frame: Frame, ra, dec, axis_ratio: float, pa_deg: float) -> np.ndarray:
    """Semi-major-axis radius (arcsec) of the ellipse through each pixel."""
    east, north = offsets_arcsec(frame, ra, dec)
    pa = np.radians(pa_deg if np.isfinite(pa_deg) else 0.0)
    along = east * np.sin(pa) + north * np.cos(pa)
    across = east * np.cos(pa) - north * np.sin(pa)
    return np.hypot(along, across / max(axis_ratio, 0.05))


def ellipse_mask(frame: Frame, ra, dec, semi_major_arcsec, axis_ratio, pa_deg) -> np.ndarray:
    return elliptical_radius(frame, ra, dec, axis_ratio, pa_deg) <= semi_major_arcsec


def rectangle_mask(frame: Frame, ra, dec, length_arcsec, width_arcsec, pa_deg) -> np.ndarray:
    """A rectangle with its long side ``length`` along position angle ``pa``."""
    east, north = offsets_arcsec(frame, ra, dec)
    pa = np.radians(pa_deg)
    along = east * np.sin(pa) + north * np.cos(pa)
    across = east * np.cos(pa) - north * np.sin(pa)
    return (np.abs(along) <= length_arcsec / 2) & (np.abs(across) <= width_arcsec / 2)


def d25_mask(frame: Frame, entry: pd.Series, scale: float = 1.0) -> np.ndarray:
    a = 0.5 * entry["d25_arcmin"] * 60.0 * scale
    q = entry["axis_ratio"] if np.isfinite(entry["axis_ratio"]) else 1.0
    return ellipse_mask(frame, entry["ra"], entry["dec"], a, q, entry["pa_deg"])


# --------------------------------------------------------------------------
# Foreground stars
# --------------------------------------------------------------------------

def foreground_stars(gaia: pd.DataFrame) -> pd.DataFrame:
    """Gaia sources that are certainly Milky Way stars, not knots in the galaxy.

    A star in our own Galaxy shows a measurable parallax or proper motion; a
    star cluster or HII region in a galaxy millions of light-years away shows
    neither. Require one of them at more than 5 sigma.
    """
    g = gaia.copy()
    pm_err = np.hypot(g["pmra"] * g["pmra_error"], g["pmdec"] * g["pmdec_error"]) / g["pm"]
    pm_sig = g["pm"] / pm_err
    secure = (g["parallax_over_error"] > 5) | (pm_sig > 5)
    return g[secure.fillna(False)].reset_index(drop=True)


def star_mask_radius_arcsec(gmag) -> np.ndarray:
    """How far out a star's light matters: the Legacy Surveys' star-mask radius,
    1630" x 1.396^-G (``mask_radius_for_mag`` in legacypipe's reference.py),
    floored at 4" for faint stars and capped at 90"."""
    return np.clip(1630.0 * 1.396 ** (-np.asarray(gmag, float)), 4.0, 90.0)


def star_mask(frame: Frame, stars: pd.DataFrame, *, scale: float = 1.0) -> np.ndarray:
    ny, nx = frame.shape
    mask = np.zeros((ny, nx), bool)
    x, y = frame.world_to_pixel(stars["ra"].to_numpy(), stars["dec"].to_numpy())
    radii = star_mask_radius_arcsec(stars["phot_g_mean_mag"].to_numpy()) * scale / frame.pixscale
    for xi, yi, ri in zip(x, y, radii):
        if not (-ri <= xi < nx + ri and -ri <= yi < ny + ri):
            continue
        x0, x1 = int(max(0, np.floor(xi - ri))), int(min(nx, np.ceil(xi + ri) + 1))
        y0, y1 = int(max(0, np.floor(yi - ri))), int(min(ny, np.ceil(yi + ri) + 1))
        yy, xx = np.mgrid[y0:y1, x0:x1]
        mask[y0:y1, x0:x1] |= (xx - xi) ** 2 + (yy - yi) ** 2 <= ri**2
    return mask


def fill_star_masks(frame: Frame, stars: pd.DataFrame, mask: np.ndarray) -> dict[str, np.ndarray]:
    """Surface-brightness images with every foreground star replaced by the
    galaxy light around it.

    Each star's disc (the same radius rule as ``star_mask``) is filled with the
    median of a ring just outside it, per band, using only pixels outside every
    star mask. Simply dropping masked pixels would also drop the galaxy light
    underneath them -- several per cent of a galaxy's flux for the generous
    radii used here. Brightest stars are filled first, and a ring that is
    mostly masked is widened until it has enough clean pixels.
    """
    filled = {b: frame.sb[b].copy() for b in BANDS}
    ny, nx = frame.shape
    x, y = frame.world_to_pixel(stars["ra"].to_numpy(), stars["dec"].to_numpy())
    radii = star_mask_radius_arcsec(stars["phot_g_mean_mag"].to_numpy()) / frame.pixscale
    for k in np.argsort(stars["phot_g_mean_mag"].to_numpy()):
        xi, yi, ri = x[k], y[k], radii[k]
        if not (-ri <= xi < nx + ri and -ri <= yi < ny + ri):
            continue
        for width in (0.5, 1.0, 2.0):
            ro = ri * (1.0 + width) + 2
            x0, x1 = int(max(0, np.floor(xi - ro))), int(min(nx, np.ceil(xi + ro) + 1))
            y0, y1 = int(max(0, np.floor(yi - ro))), int(min(ny, np.ceil(yi + ro) + 1))
            yy, xx = np.mgrid[y0:y1, x0:x1]
            rr = np.hypot(xx - xi, yy - yi)
            disc = rr <= ri
            ring = (rr > ri) & (rr <= ro) & ~mask[y0:y1, x0:x1]
            if ring.sum() >= 30:
                break
        if ring.sum() < 10 or not disc.any():
            continue
        for b in BANDS:
            cut = filled[b][y0:y1, x0:x1]
            cut[disc] = np.median(cut[ring])
    return filled


# --------------------------------------------------------------------------
# Measuring
# --------------------------------------------------------------------------

def aperture_flux(sb: np.ndarray, mask: np.ndarray, pixscale: float) -> float:
    """Total flux (nanomaggies) of a surface-brightness image inside a mask."""
    return float(np.nansum(sb[mask]) * pixscale**2)


def sky_statistics(frame: Frame, exclude: np.ndarray) -> dict[str, dict[str, float]]:
    """Residual sky level and pixel noise per band, from pixels not excluded."""
    out = {}
    for b in BANDS:
        vals = frame.sb[b][~exclude]
        mean, median, std = sigma_clipped_stats(vals, sigma=3.0, maxiters=10)
        out[b] = {"median": float(median), "std": float(std), "n_pixels": int(vals.size)}
    return out


def star_photometry(frame: Frame, ra: np.ndarray, dec: np.ndarray, *, r_ap: float = 7.43,
                    r_in: float = 12.0, r_out: float = 20.0) -> pd.DataFrame:
    """Aperture magnitudes of point sources, sky from a surrounding annulus.

    7.43" is SDSS's own aperture 7, the radius its PSF magnitudes are
    normalised to. Also returns the local background surface brightness, so
    stars sitting on bright galaxy light (where the annulus cannot see the
    true sky) can be excluded.
    """
    ny, nx = frame.shape
    x, y = frame.world_to_pixel(np.asarray(ra), np.asarray(dec))
    rows = []
    rpx_out = r_out / frame.pixscale
    for xi, yi in zip(x, y):
        margin = rpx_out + 2
        rec = {"x": xi, "y": yi, "inside": bool(margin <= xi < nx - margin and margin <= yi < ny - margin)}
        if rec["inside"]:
            x0, x1 = int(np.floor(xi - rpx_out)) - 1, int(np.ceil(xi + rpx_out)) + 2
            y0, y1 = int(np.floor(yi - rpx_out)) - 1, int(np.ceil(yi + rpx_out)) + 2
            yy, xx = np.mgrid[y0:y1, x0:x1]
            r = np.hypot(xx - xi, yy - yi) * frame.pixscale
            ap, ann = r <= r_ap, (r >= r_in) & (r <= r_out)
            for b in BANDS:
                cut = frame.sb[b][y0:y1, x0:x1]
                bkg = float(np.median(cut[ann]))
                flux = float((cut[ap] - bkg).sum() * frame.pixscale**2)
                rec[f"flux_{b}"] = flux
                rec[f"bkg_{b}"] = bkg
                rec[f"mag_{b}"] = 22.5 - 2.5 * np.log10(flux) if flux > 0 else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)
