"""Radio images: noise, smoothing, where a giant ends, and its spectrum.

A LoTSS image is in Jy per beam: the flux density inside one resolution
element. Three operations recur:

* **Noise.** The rms of the empty sky, estimated with iterative 3-sigma
  clipping, so that emission and artefacts do not inflate it.
* **Smoothing.** Convolving to a broader beam trades resolution for
  sensitivity to faint, extended emission: a lobe spread over many 20" beams
  gains roughly a factor (60/20) in signal-to-noise at 60". A Gaussian beam of
  FWHM theta_1 convolved with a Gaussian of FWHM sqrt(theta_2^2 - theta_1^2)
  gives theta_2, and Jy/beam values scale by the ratio of beam areas.
* **Where a giant ends.** Pixels above 3 sigma form islands; the islands
  containing the hand-placed lobe seeds (plus any extra pieces) are the giant.
  Its largest angular size (LAS) is the greatest separation between any two of
  its pixels; we also measure the distance between the brightest points at
  each end (hotspot to hotspot). The 3-sigma edge of a bright hotspot lies
  beyond the hotspot itself by up to about one beam, so edge-to-edge lengths
  grow with the beam: we measure at the finest resolution that still shows
  the whole giant, and report the others.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from scipy import ndimage
from scipy.spatial import ConvexHull

import config
from config import IMAGE_DIR, Giant

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


@dataclass
class RadioImage:
    data: np.ndarray  # Jy/beam
    wcs: WCS
    beam_arcsec: float
    pixel_arcsec: float
    freq_mhz: float
    survey: str
    header: fits.Header = field(repr=False, default=None)

    @property
    def beam_area_pix(self) -> float:
        """Area of the Gaussian beam in pixels: pi theta^2 / (4 ln 2)."""
        return np.pi * self.beam_arcsec**2 / (4 * np.log(2)) / self.pixel_arcsec**2

    def world_to_pixel(self, ra, dec) -> tuple[np.ndarray, np.ndarray]:
        x, y = self.wcs.world_to_pixel_values(np.atleast_1d(ra), np.atleast_1d(dec))
        return np.asarray(x), np.asarray(y)

    def pixel_to_world(self, x, y) -> tuple[np.ndarray, np.ndarray]:
        ra, dec = self.wcs.pixel_to_world_values(np.atleast_1d(x), np.atleast_1d(y))
        return np.asarray(ra), np.asarray(dec)


def load_lotss(key: str, resolution: int, *, check: bool = False) -> RadioImage:
    """A giant's LoTSS image: from its own release (DR2, or DR3 for giants
    outside the DR2 footprint), stitched from tiles if it needed a mosaic.
    ``check=True`` loads the DR3 20" image of a DR2 giant, kept to compare the
    two releases."""
    g = config.giant_by_key(key)
    if check:
        path = IMAGE_DIR / f"{key}_dr3lotss20.fits"
    elif g.tiles > 1:
        if resolution != 20:
            raise FileNotFoundError(f"{key}: only a 20\" mosaic is fetched")
        path = _mosaic(g)
    else:
        path = IMAGE_DIR / f"{key}_lotss{resolution}.fits"
    hdu = fits.open(path)[0]
    data = np.squeeze(hdu.data).astype(np.float64)
    wcs = WCS(hdu.header).celestial
    pix = abs(hdu.header["CDELT2"]) * 3600
    beam = hdu.header["BMAJ"] * 3600
    release = "DR3" if (check or g.survey == "dr3") else "DR2"
    return RadioImage(data, wcs, beam, pix, config.LOTSS_FREQ_MHZ, f"LoTSS {release} {resolution}\"", hdu.header)


def standard_wcs(ra0: float, dec0: float, pixel_arcsec: float, npix: int) -> WCS:
    """A north-up, east-left TAN grid centred on (ra0, dec0)."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [(npix + 1) / 2, (npix + 1) / 2]
    w.wcs.cdelt = [-pixel_arcsec / 3600, pixel_arcsec / 3600]
    return w


def _mosaic(g: config.Giant):
    """Stitch a giant's 2 x 2 tiles onto one TAN grid centred on its host
    (mean where tiles overlap); cached in data/derived/."""
    from reproject import reproject_interp
    from reproject.mosaicking import reproject_and_coadd

    out = config.DERIVED_DIR / f"{g.key}_lotss20_mosaic.fits"
    tiles = sorted(IMAGE_DIR.glob(f"{g.key}_lotss20_tile*.fits"))
    if out.exists() and all(out.stat().st_mtime >= t.stat().st_mtime for t in tiles):
        return out
    if len(tiles) != g.tiles**2:
        raise FileNotFoundError(f"{g.key}: {len(tiles)} of {g.tiles ** 2} mosaic tiles cached")
    hdus = [fits.open(t)[0] for t in tiles]
    pix = abs(hdus[0].header["CDELT2"]) * 3600
    npix = int(np.ceil(g.fov_arcmin * 60 / pix))
    w = standard_wcs(g.ra, g.dec, pix, npix)
    inputs = [(np.squeeze(h.data).astype(np.float64), WCS(h.header).celestial) for h in hdus]
    data, footprint = reproject_and_coadd(inputs, w, shape_out=(npix, npix), reproject_function=reproject_interp,
                                          combine_function="mean")
    data[footprint == 0] = np.nan
    hdr = w.to_header()
    for k in ("BMAJ", "BMIN", "BPA", "BUNIT", "RESTFRQ", "TELESCOP"):
        if k in hdus[0].header:
            hdr[k] = hdus[0].header[k]
    hdr["HISTORY"] = f"Mosaic of {len(tiles)} LoTSS {g.survey.upper()} cutouts stitched by I18/radiomaps.py"
    fits.PrimaryHDU(data.astype(np.float32), hdr).writeto(out, overwrite=True)
    return out


def load_nvss(key: str) -> RadioImage:
    path = IMAGE_DIR / f"{key}_nvss.fits"
    hdu = fits.open(path)[0]
    data = np.squeeze(hdu.data).astype(np.float64)
    data[data == 0] = np.nan  # SkyView marks pixels outside the survey with 0
    wcs = WCS(hdu.header).celestial
    pix = abs(hdu.header["CDELT2"]) * 3600
    return RadioImage(data, wcs, config.NVSS_BEAM_ARCSEC, pix, config.NVSS_FREQ_MHZ, "NVSS", hdu.header)


def noise(img: RadioImage) -> float:
    """rms of empty sky: iterative 3-sigma clipping of all finite pixels."""
    d = img.data[np.isfinite(img.data)]
    _, _, sd = sigma_clipped_stats(d, sigma=3.0, maxiters=20)
    return float(sd)


def smooth_to(img: RadioImage, beam_arcsec: float) -> RadioImage:
    """Convolve to a broader circular Gaussian beam, keeping Jy/beam units.

    NaN pixels (outside the mosaic) are handled by normalised convolution so
    the edge does not bleed in as zeros, and stay NaN afterwards.
    """
    if beam_arcsec <= img.beam_arcsec + 1e-6:
        return img
    k_fwhm = np.sqrt(beam_arcsec**2 - img.beam_arcsec**2)
    k_sigma_pix = k_fwhm / FWHM_PER_SIGMA / img.pixel_arcsec
    valid = np.isfinite(img.data)
    num = ndimage.gaussian_filter(np.where(valid, img.data, 0.0), k_sigma_pix, mode="constant")
    den = ndimage.gaussian_filter(valid.astype(float), k_sigma_pix, mode="constant")
    with np.errstate(invalid="ignore", divide="ignore"):
        out = num / den * (beam_arcsec / img.beam_arcsec) ** 2
    out[~valid | (den < 0.5)] = np.nan
    return RadioImage(out, img.wcs, beam_arcsec, img.pixel_arcsec, img.freq_mhz,
                      f"{img.survey} smoothed to {beam_arcsec:.0f}\"", img.header)


# --------------------------------------------------------------------------
# Where a giant ends
# --------------------------------------------------------------------------

@dataclass
class Extent:
    beam_arcsec: float
    sigma_jy: float
    threshold_sigma: float
    mask: np.ndarray
    las_arcmin: float
    las_endpoints: tuple  # ((ra, dec), (ra, dec))
    peak_to_peak_arcmin: float | None
    peaks: tuple
    pa_deg: float  # position angle of the LAS axis, north through east
    n_islands: int
    flux_jy: float  # integrated flux density inside the mask
    missing_seeds: list


def _sky_offsets_arcsec(img: RadioImage, x, y, ra0, dec0) -> tuple[np.ndarray, np.ndarray]:
    """Tangent-plane offsets (east, north) in arcsec of pixels from (ra0, dec0)."""
    ra, dec = img.pixel_to_world(x, y)
    ra, dec, ra0r, dec0r = map(np.radians, (ra, dec, ra0, dec0))
    cosc = np.sin(dec0r) * np.sin(dec) + np.cos(dec0r) * np.cos(dec) * np.cos(ra - ra0r)
    xi = np.cos(dec) * np.sin(ra - ra0r) / cosc
    eta = (np.cos(dec0r) * np.sin(dec) - np.sin(dec0r) * np.cos(dec) * np.cos(ra - ra0r)) / cosc
    return np.degrees(xi) * 3600, np.degrees(eta) * 3600


def separation_arcmin(ra1, dec1, ra2, dec2) -> np.ndarray:
    """Great-circle separation (haversine)."""
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    h = np.sin((dec2 - dec1) / 2) ** 2 + np.cos(dec1) * np.cos(dec2) * np.sin((ra2 - ra1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(h))) * 60


def position_angle(ra1, dec1, ra2, dec2) -> float:
    """Position angle (deg, north through east) of point 2 seen from point 1."""
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    pa = np.arctan2(np.sin(ra2 - ra1), np.cos(dec1) * np.tan(dec2) - np.sin(dec1) * np.cos(ra2 - ra1))
    return float(np.degrees(pa) % 360)


def islands(img: RadioImage, sigma: float, threshold: float) -> tuple[np.ndarray, int]:
    above = np.isfinite(img.data) & (img.data >= threshold * sigma)
    lab, n = ndimage.label(above, structure=np.ones((3, 3)))
    if n:
        sizes = ndimage.sum_labels(np.ones_like(lab), lab, index=np.arange(1, n + 1))
        small = np.where(sizes < config.MIN_ISLAND_BEAMS * img.beam_area_pix)[0] + 1
        lab[np.isin(lab, small)] = 0
    return lab, n


def measure_extent(img: RadioImage, g: Giant, *, threshold: float = config.LAS_SIGMA,
                   sigma: float | None = None) -> Extent:
    """The giant's 3-sigma extent on one image, from its lobe seeds."""
    sigma = noise(img) if sigma is None else sigma
    lab, n = islands(img, sigma, threshold)
    seeds = list(g.seeds) + list(g.extra_seeds)
    member_labels, missing = set(), []
    for ra, dec in seeds:
        x, y = img.world_to_pixel(ra, dec)
        xi, yi = int(round(x[0])), int(round(y[0]))
        if not (0 <= yi < lab.shape[0] and 0 <= xi < lab.shape[1]):
            missing.append((ra, dec))
            continue
        # tolerate a seed placed a few pixels off its island
        r = max(1, int(round(0.5 * img.beam_arcsec / img.pixel_arcsec)))
        win = lab[max(0, yi - r):yi + r + 1, max(0, xi - r):xi + r + 1]
        labs = set(win[win > 0].ravel().tolist())
        if not labs:
            missing.append((ra, dec))
            continue
        member_labels |= labs
    for ra, dec in g.exclude:
        x, y = img.world_to_pixel(ra, dec)
        xi, yi = int(round(x[0])), int(round(y[0]))
        if 0 <= yi < lab.shape[0] and 0 <= xi < lab.shape[1] and lab[yi, xi] > 0:
            member_labels.discard(int(lab[yi, xi]))
    mask = np.isin(lab, list(member_labels)) if member_labels else np.zeros_like(lab, bool)
    yy, xx = np.nonzero(mask)
    las, ends, pa = np.nan, ((np.nan, np.nan), (np.nan, np.nan)), np.nan
    if xx.size >= 3:
        ex, ny = _sky_offsets_arcsec(img, xx, yy, g.ra, g.dec)
        pts = np.column_stack([ex, ny])
        try:
            hull = pts[ConvexHull(pts).vertices]
        except Exception:
            hull = pts
        d2 = ((hull[:, None, :] - hull[None, :, :]) ** 2).sum(-1)
        i, j = np.unravel_index(np.argmax(d2), d2.shape)
        # pixel centres underestimate the extent by up to half a pixel at each end
        las = (np.sqrt(d2[i, j]) + img.pixel_arcsec) / 60.0
        hi = np.where((ex == hull[i, 0]) & (ny == hull[i, 1]))[0][0]
        hj = np.where((ex == hull[j, 0]) & (ny == hull[j, 1]))[0][0]
        e1 = img.pixel_to_world(xx[hi], yy[hi])
        e2 = img.pixel_to_world(xx[hj], yy[hj])
        ends = ((float(e1[0][0]), float(e1[1][0])), (float(e2[0][0]), float(e2[1][0])))
        pa = position_angle(*ends[0], *ends[1]) % 180
    # hotspot to hotspot: the brightest pixel of the giant near each end seed
    # (within a fifth of the giant's length or three beams, whichever is
    # larger), so that a bright core in the same island is not mistaken for a
    # hotspot.
    peaks, p2p = [], None
    if mask.any() and np.isfinite(las):
        radius_arcmin = max(3 * img.beam_arcsec / 60, 0.2 * las)
        for ra, dec in g.seeds:
            ex, ny = _sky_offsets_arcsec(img, xx, yy, ra, dec)
            near = np.hypot(ex, ny) <= radius_arcmin * 60
            if near.any():
                k = np.argmax(np.where(near, img.data[yy, xx], -np.inf))
                pr, pd = img.pixel_to_world(xx[k], yy[k])
                peaks.append((float(pr[0]), float(pd[0])))
    if len(peaks) == 2:
        p2p = float(separation_arcmin(*peaks[0], *peaks[1]))
    flux = float(np.nansum(np.where(mask, img.data, 0.0)) / img.beam_area_pix)
    return Extent(img.beam_arcsec, sigma, threshold, mask, float(las), ends, p2p, tuple(peaks),
                  float(pa), len(member_labels), flux, missing)


# --------------------------------------------------------------------------
# Regridding
# --------------------------------------------------------------------------

def tangent_wcs(ra0: float, dec0: float, pixel_arcsec: float, npix: int, *, rotation_deg: float = 0.0) -> WCS:
    """A TAN grid centred on (ra0, dec0); ``rotation_deg`` turns the sky so that
    position angle ``rotation_deg`` points to +x (used to lay giants flat)."""
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [(npix + 1) / 2, (npix + 1) / 2]
    s = pixel_arcsec / 3600
    # rotate so that the given PA (north through east) lies along +x
    th = np.radians(270.0 - rotation_deg)
    c, sn = np.cos(th), np.sin(th)
    # +x points to position angle ``rotation_deg``; +y is 90 deg further (east of
    # it), the same handedness as an ordinary north-up, east-left image.
    w.wcs.cd = np.array([[-s * c, -s * sn], [-s * sn, s * c]])
    return w


def reproject_to(img: RadioImage, wcs_out: WCS, shape: tuple[int, int]) -> np.ndarray:
    from reproject import reproject_interp

    out, _ = reproject_interp((img.data, img.wcs), wcs_out, shape_out=shape, order="bilinear")
    return out
