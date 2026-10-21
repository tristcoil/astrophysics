"""Reading SDO images and the geometry of the solar disc.

Two instruments on NASA's Solar Dynamics Observatory feed this project:

* **AIA** (Atmospheric Imaging Assembly) photographs the whole Sun every 12 s in
  seven extreme-ultraviolet and three ultraviolet/visible channels. Each EUV channel
  is tuned to emission lines of one or two ions, so each shows gas in a narrow
  temperature range.
* **HMI** (Helioseismic and Magnetic Imager) images the Sun in visible light and
  measures the magnetic field along the line of sight, from the Zeeman splitting of
  the iron line at 6173 A.

All images here are "level 1.5": the Sun is centred, solar north is up, and the
pixel scale is uniform, so a pixel's position maps straight onto the Sun through
the World Coordinate System in each file's header. ``sunpy.map`` does that mapping.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy import ndimage

warnings.filterwarnings("ignore", module="sunpy")
warnings.filterwarnings("ignore", module="astropy")

import sunpy.map  # noqa: E402  (after the warning filters: sunpy is chatty about NRT headers)
from sunpy.coordinates import frames  # noqa: E402


def load_map(path: Path | str) -> "sunpy.map.GenericMap":
    """A map exactly as stored: AIA in DN per 0.6" pixel, HMI in gauss or DN/s."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sunpy.map.Map(str(path))


def load_aia(path: Path | str) -> "sunpy.map.GenericMap":
    """An AIA map in DN/s: the files hold counts per exposure, and AIA's automatic
    exposure control shortens exposures during flares, so only the rate compares."""
    m = load_map(path)
    exptime = float(m.meta["exptime"])
    data = m.data.astype(np.float64) / exptime
    meta = m.meta.copy()
    meta["bunit"] = "DN/s"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return sunpy.map.Map(data, meta)


def radius_map(m) -> np.ndarray:
    """Every pixel's distance from disc centre in units of the solar radius.

    Off the disc this is also the line of sight's closest approach to the Sun's
    centre, so ``(r - 1) * R_sun`` is the height above the limb it samples.
    """
    ny, nx = m.data.shape
    y, x = np.mgrid[0:ny, 0:nx]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = m.pixel_to_world(x * u.pix, y * u.pix)
    rho = np.hypot(c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec))
    return rho / m.rsun_obs.to_value(u.arcsec)


def mu_from_radius(r: np.ndarray) -> np.ndarray:
    """cos(theta), theta = angle between the local vertical and the line of sight.

    Exact for an observer at any distance when r is the ratio of sines of the angle
    from disc centre (sin(theta) = r, from the law of sines); for the Sun at 1 AU the
    angles are small enough that the ratio of angles used here is the same to 1e-5.
    """
    return np.sqrt(np.clip(1.0 - np.asarray(r, dtype=float) ** 2, 0.0, None))


def position_angle_map(m) -> np.ndarray:
    """Position angle of each pixel, degrees counter-clockwise from solar north."""
    ny, nx = m.data.shape
    y, x = np.mgrid[0:ny, 0:nx]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c = m.pixel_to_world(x * u.pix, y * u.pix)
    tx, ty = c.Tx.to_value(u.arcsec), c.Ty.to_value(u.arcsec)
    # north = +y, east = -x (solar east is on the left in these images)
    return np.degrees(np.arctan2(-tx, ty)) % 360.0


def arcsec_to_mm(m) -> float:
    """Megametres per arcsecond at the Sun's distance."""
    return float(m.dsun.to_value(u.Mm) * np.deg2rad(1 / 3600))


def locate_flare(pre_path: Path, peak_path: Path, *, smooth_px: float = 2.0,
                 max_radius: float = 1.15) -> dict:
    """Find a flare as the spot that brightened most between two images.

    Both images are converted to DN/s, the pre-flare image is subtracted, and the
    difference is smoothed a little so a single noisy or cosmic-ray pixel cannot
    win. Returns the position in helioprojective arcseconds (and heliographic
    latitude/longitude when it lies on the disc).
    """
    pre, peak = load_aia(pre_path), load_aia(peak_path)
    diff = ndimage.gaussian_filter(peak.data - pre.data, smooth_px)
    r = radius_map(peak)
    diff[r > max_radius] = -np.inf
    iy, ix = np.unravel_index(np.argmax(diff), diff.shape)
    # refine with the brightening-weighted centroid in a small box
    y0, y1, x0, x1 = iy - 4, iy + 5, ix - 4, ix + 5
    box = np.clip(diff[y0:y1, x0:x1], 0, None)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    cy, cx = (box * yy).sum() / box.sum(), (box * xx).sum() / box.sum()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hpc = peak.pixel_to_world(cx * u.pix, cy * u.pix)
    out = {
        "hpc_x_arcsec": float(hpc.Tx.to_value(u.arcsec)),
        "hpc_y_arcsec": float(hpc.Ty.to_value(u.arcsec)),
        "radius_rsun": float(np.hypot(hpc.Tx.to_value(u.arcsec), hpc.Ty.to_value(u.arcsec))
                             / peak.rsun_obs.to_value(u.arcsec)),
        "brightening_dn_s": float(diff[iy, ix]),
        "channel": int(peak.wavelength.to_value(u.AA)),
        "pre_time": pre.date.isot,
        "peak_time": peak.date.isot,
    }
    if out["radius_rsun"] < 1.0:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            hgs = hpc.transform_to(frames.HeliographicStonyhurst(obstime=peak.date))
        out["hgs_lon_deg"] = float(hgs.lon.to_value(u.deg))
        out["hgs_lat_deg"] = float(hgs.lat.to_value(u.deg))
    return out


def save_cutout(full_path: Path, dest: Path, x_arcsec: float, y_arcsec: float,
                half_width_arcsec: float) -> None:
    """Cut a square around (x, y) out of a full frame and save it, header and all.

    The data are stored exactly as in the original file (DN per 0.6" pixel, before
    exposure normalisation), losslessly compressed.
    """
    m = load_map(full_path)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bl = SkyCoord((x_arcsec - half_width_arcsec) * u.arcsec,
                      (y_arcsec - half_width_arcsec) * u.arcsec, frame=m.coordinate_frame)
        tr = SkyCoord((x_arcsec + half_width_arcsec) * u.arcsec,
                      (y_arcsec + half_width_arcsec) * u.arcsec, frame=m.coordinate_frame)
        sub = m.submap(bl, top_right=tr)
        header = sub.fits_header
    hdu = fits.CompImageHDU(data=np.asarray(sub.data, dtype=np.float32), header=header,
                            compression_type="GZIP_2")
    tmp = dest.with_name(dest.name + ".part")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(tmp, overwrite=True)
    tmp.replace(dest)
