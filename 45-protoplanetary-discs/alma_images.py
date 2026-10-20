"""Reading the ALMA images: pixels, sky offsets, beam, noise, flux and brightness temperature.

An ALMA image is not a photograph. The interferometer measures the Fourier transform of
the sky on the baselines between its antennas, and an image is reconstructed from those
"visibilities" (here with CASA's CLEAN). Its pixels are in Jy/beam: the flux density in
one resolution element, the "synthesized beam", an elliptical Gaussian of a few tens of
milliarcseconds. Three things follow that matter throughout this project:

* Neighbouring pixels are not independent; independent samples are a beam apart.
* Everything finer than the beam is smeared: narrow gaps look shallower and wider.
* Converting Jy/beam to a temperature needs the beam's solid angle.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy.io import fits

from config import ARCSEC, C_CGS, H_CGS, HLTAU_DIR, HLTAU_IMAGE_NAME, IMAGE_DIR, JY_CGS, K_B_CGS


@dataclass
class Image:
    key: str
    data: np.ndarray  # (ny, nx), Jy/beam
    x: np.ndarray  # (nx,) RA offset of each column from the phase centre, arcsec, positive = east
    y: np.ndarray  # (ny,) Dec offset of each row, arcsec, positive = north
    pix: float  # pixel size, arcsec
    bmaj: float  # beam FWHM, arcsec
    bmin: float
    bpa: float  # deg east of north
    freq_hz: float
    ra0_deg: float
    dec0_deg: float
    header: fits.Header

    @property
    def beam_sr(self) -> float:
        """Solid angle of the Gaussian beam: pi bmaj bmin / (4 ln 2)."""
        return np.pi * (self.bmaj * ARCSEC) * (self.bmin * ARCSEC) / (4.0 * np.log(2.0))

    @property
    def beam_pix(self) -> float:
        """Beam area in pixels, to turn a sum of Jy/beam pixels into Jy."""
        return np.pi * self.bmaj * self.bmin / (4.0 * np.log(2.0)) / self.pix**2

    @property
    def beam_mean(self) -> float:
        return float(np.sqrt(self.bmaj * self.bmin))

    def grids(self) -> tuple[np.ndarray, np.ndarray]:
        return np.meshgrid(self.x, self.y)

    def cutout(self, half_width: float, x0: float = 0.0, y0: float = 0.0) -> "Image":
        """A square sub-image of +/- half_width arcsec around (x0, y0)."""
        ix = np.where(np.abs(self.x - x0) <= half_width)[0]
        iy = np.where(np.abs(self.y - y0) <= half_width)[0]
        return Image(self.key, self.data[iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1].copy(), self.x[ix[0]:ix[-1] + 1].copy(),
                     self.y[iy[0]:iy[-1] + 1].copy(), self.pix, self.bmaj, self.bmin, self.bpa, self.freq_hz,
                     self.ra0_deg, self.dec0_deg, self.header)


def _axes(h: fits.Header, nx: int, ny: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Offsets from the reference pixel. Over a few arcsec the SIN projection is linear to
    better than 1e-9, so offset = (pixel - CRPIX) * CDELT."""
    if abs(abs(h["CDELT1"]) - abs(h["CDELT2"])) > 1e-12:
        raise ValueError("non-square pixels")
    x = (np.arange(nx) + 1 - h["CRPIX1"]) * h["CDELT1"] * 3600.0  # CDELT1 < 0: east is to the left
    y = (np.arange(ny) + 1 - h["CRPIX2"]) * h["CDELT2"] * 3600.0
    return x, y, abs(h["CDELT2"]) * 3600.0


def _freq(h: fits.Header) -> float:
    for i in range(1, h["NAXIS"] + 1):
        if h.get(f"CTYPE{i}", "").startswith("FREQ"):
            return float(h[f"CRVAL{i}"])
    raise KeyError("no frequency axis")


def read_image(path, key: str) -> Image:
    with fits.open(path, memmap=False) as hdul:
        h = hdul[0].header
        data = np.asarray(hdul[0].data, dtype=np.float64)
    data = np.squeeze(data)
    if data.ndim != 2:
        raise ValueError(f"{path}: expected a 2-D image, got shape {data.shape}")
    ny, nx = data.shape
    x, y, pix = _axes(h, nx, ny)
    return Image(key, data, x, y, pix, h["BMAJ"] * 3600.0, h["BMIN"] * 3600.0, float(h["BPA"]), _freq(h),
                 float(h["CRVAL1"]), float(h["CRVAL2"]), h)


def load_continuum(key: str) -> Image:
    if key == "HLTau":
        return read_image(HLTAU_DIR / HLTAU_IMAGE_NAME, key)
    return read_image(IMAGE_DIR / f"{key}_continuum.fits", key)


# ---------------------------------------------------------------------------
# Noise and flux
# ---------------------------------------------------------------------------
def robust_std(a: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    return float(1.4826 * np.median(np.abs(a - np.median(a))))


def measure_rms(img: Image, exclude_radius: float, x0: float = 0.0, y0: float = 0.0,
                outer_radius: float | None = None) -> float:
    """Noise from pixels outside ``exclude_radius`` arcsec of the disc (median absolute
    deviation, so stray emission does not inflate it)."""
    xx, yy = img.grids()
    rr = np.hypot(xx - x0, yy - y0)
    sel = rr > exclude_radius
    if outer_radius is not None:
        sel &= rr < outer_radius
    return robust_std(img.data[sel])


def flux_in_ellipse(img: Image, x0: float, y0: float, a: float, inc_deg: float, pa_deg: float) -> float:
    """Total flux density (Jy) inside a sky ellipse of semi-major axis a (arcsec)."""
    from geometry import disc_polar  # local import: geometry imports this module

    xx, yy = img.grids()
    r, _ = disc_polar(xx, yy, x0, y0, inc_deg, pa_deg)
    return float(img.data[r <= a].sum() / img.beam_pix)


# ---------------------------------------------------------------------------
# Brightness temperature
# ---------------------------------------------------------------------------
def jy_per_beam_to_cgs(i_jyb: np.ndarray | float, img: Image) -> np.ndarray:
    """Jy/beam -> specific intensity in erg s^-1 cm^-2 Hz^-1 sr^-1."""
    return np.asarray(i_jyb) * JY_CGS / img.beam_sr


def planck(nu: float, t: np.ndarray | float) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    with np.errstate(over="ignore", divide="ignore"):
        return 2.0 * H_CGS * nu**3 / C_CGS**2 / np.expm1(H_CGS * nu / (K_B_CGS * t))


def brightness_temperature(i_cgs: np.ndarray | float, nu: float, rayleigh_jeans: bool = False) -> np.ndarray:
    """Invert the Planck function (or its Rayleigh-Jeans limit); NaN where I <= 0."""
    i = np.asarray(i_cgs, dtype=float)
    if rayleigh_jeans:
        return i * C_CGS**2 / (2.0 * nu**2 * K_B_CGS)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = H_CGS * nu / K_B_CGS / np.log1p(2.0 * H_CGS * nu**3 / (C_CGS**2 * i))
    return np.where(i > 0, out, np.nan)


# ---------------------------------------------------------------------------
# Spectral cubes (CO)
# ---------------------------------------------------------------------------
@dataclass
class Cube:
    key: str
    data: np.ndarray  # (nchan, ny, nx) Jy/beam, cropped
    velocity: np.ndarray  # (nchan,) km/s, radio convention, LSRK
    x: np.ndarray
    y: np.ndarray
    pix: float
    bmaj: np.ndarray  # per channel, arcsec
    bmin: np.ndarray
    bpa: np.ndarray
    rest_hz: float

    @property
    def dv(self) -> float:
        return float(np.median(np.diff(self.velocity)))


def read_cube(path, key: str, half_width: float) -> Cube:
    """Read the central +/- half_width arcsec of a CASA cube, channel by channel from a
    memory map (the full 2000x2000x105 cube is 1.7 GB)."""
    with fits.open(path, memmap=True) as hdul:
        h = hdul[0].header
        raw = hdul[0].data  # (stokes, chan, y, x)
        nx, ny, nchan = h["NAXIS1"], h["NAXIS2"], h["NAXIS3"]
        x, y, pix = _axes(h, nx, ny)
        ix = np.where(np.abs(x) <= half_width)[0]
        iy = np.where(np.abs(y) <= half_width)[0]
        sub = np.empty((nchan, iy.size, ix.size), dtype=np.float32)
        for k in range(nchan):
            sub[k] = raw[0, k, iy[0]:iy[-1] + 1, ix[0]:ix[-1] + 1]
        freq = h["CRVAL3"] + (np.arange(nchan) + 1 - h["CRPIX3"]) * h["CDELT3"]
        rest = float(h["RESTFRQ"])
        if "BEAMS" in hdul:
            beams = hdul["BEAMS"].data
            unit = hdul["BEAMS"].header.get("TUNIT1", "arcsec").strip().lower()
            scale = 3600.0 if unit.startswith("deg") else 1.0
            bmaj = np.asarray(beams["BMAJ"], float) * scale
            bmin = np.asarray(beams["BMIN"], float) * scale
            bpa = np.asarray(beams["BPA"], float)
            order = np.asarray(beams["CHAN"], int)
            bmaj, bmin, bpa = bmaj[np.argsort(order)], bmin[np.argsort(order)], bpa[np.argsort(order)]
        else:
            bmaj = np.full(nchan, h["BMAJ"] * 3600.0)
            bmin = np.full(nchan, h["BMIN"] * 3600.0)
            bpa = np.full(nchan, h["BPA"])
    velocity = C_CGS / 1e5 * (1.0 - freq / rest)  # radio velocity, km/s
    order = np.argsort(velocity)
    return Cube(key, sub[order], velocity[order], x[ix[0]:ix[-1] + 1], y[iy[0]:iy[-1] + 1], pix,
                bmaj[order], bmin[order], bpa[order], rest)
