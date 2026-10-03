"""Turning three calibrated images into one colour picture.

A picture has to squeeze a galaxy's brightness range -- the core can be ten
thousand times brighter than the outer disc -- into the 256 levels of a screen,
and three filters into three display channels. This module holds the ways of
doing that which the project compares:

* ``linear_rgb`` and ``log_rgb``: each channel stretched on its own, the
  naive way.
* ``per_channel_asinh``: the arcsinh curve, still applied to each channel
  separately.
* ``lupton_rgb``: the algorithm of Lupton et al. (2004), written for SDSS's own
  colour images. The stretch is applied to the *total* intensity and every
  channel is scaled by the same factor, so a pixel's colour -- the ratio of its
  three fluxes -- survives the stretch. Bright cores stay coloured instead of
  burning out to white.
* ``viewer_sdss_rgb``: the Legacy Surveys viewer's JPEG recipe, re-implemented
  only to prove that we read the pixels exactly as the survey does.
* ``true_color_rgb``: our estimate of what the human eye would see, built from
  the colour matrix of ``color_physics`` rather than by assigning filters to
  channels.

All inputs are surface-brightness images in nanomaggies per square arcsecond.
Outputs are float arrays in [0, 1], shape (H, W, 3), row 0 at the *bottom*
(FITS order); ``to_uint8`` flips them for writing image files.
"""

from __future__ import annotations

import io

import numpy as np
from PIL import Image

import color_physics as cp
from config import BANDS, CHANNEL_OF_BAND


# --------------------------------------------------------------------------
# The Lupton et al. (2004) stretch
# --------------------------------------------------------------------------

def softening_for(Q: float, saturation: float) -> float:
    """The knee of the arcsinh curve, below which the stretch is linear.

    Lupton's curve is f(I) = asinh(I / s) / asinh(I_sat / s): linear for I << s,
    logarithmic for I >> s, reaching 1 at I_sat. Their parameter Q fixes the
    ratio I_sat / s = sinh(Q): Q = 8 compresses a factor sinh(8) ~ 1500 (about 8
    magnitudes) between "just above the noise" and "full brightness".
    """
    return saturation / np.sinh(Q)


def asinh_curve(intensity: np.ndarray, softening: float, saturation: float) -> np.ndarray:
    return np.arcsinh(intensity / softening) / np.arcsinh(saturation / softening)


def lupton_rgb(red: np.ndarray, green: np.ndarray, blue: np.ndarray, *,
               softening: float, saturation: float, minimum: float = 0.0) -> np.ndarray:
    """Lupton et al. (2004) colour-preserving arcsinh composite.

    1. I = mean of the three channels (after subtracting ``minimum``).
    2. Every channel is multiplied by the same factor f(I) / I.
    3. A pixel whose brightest channel would exceed 1 is scaled down as a whole,
       so it saturates to *its own colour*, not to white.
    """
    rgb = np.stack([red, green, blue], axis=-1).astype(float) - minimum
    intensity = rgb.mean(axis=-1)
    f = asinh_curve(intensity, softening, saturation)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(intensity > 0, f / intensity, 0.0)
    rgb = np.clip(rgb * ratio[..., None], 0.0, None)
    peak = rgb.max(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rgb = np.where((peak > 1.0)[..., None], rgb / peak[..., None], rgb)
    return rgb


def per_channel_asinh(red, green, blue, *, softening: float, saturation: float,
                      minimum: float = 0.0) -> np.ndarray:
    """The same arcsinh curve applied to each channel on its own."""
    chans = [asinh_curve(np.clip(c - minimum, 0.0, None), softening, saturation)
             for c in (red, green, blue)]
    return np.clip(np.stack(chans, axis=-1), 0.0, 1.0)


def linear_rgb(red, green, blue, *, saturation: float, minimum: float = 0.0) -> np.ndarray:
    chans = [(c - minimum) / (saturation - minimum) for c in (red, green, blue)]
    return np.clip(np.stack(chans, axis=-1), 0.0, 1.0)


def log_rgb(red, green, blue, *, floor: float, saturation: float) -> np.ndarray:
    """log10 stretch from ``floor`` (black) to ``saturation`` (white), per channel."""
    chans = []
    for c in (red, green, blue):
        with np.errstate(divide="ignore", invalid="ignore"):
            v = (np.log10(np.clip(c, floor, None)) - np.log10(floor)) / (np.log10(saturation) - np.log10(floor))
        chans.append(v)
    return np.clip(np.stack(chans, axis=-1), 0.0, 1.0)


# --------------------------------------------------------------------------
# Filters -> channels, and the white balance
# --------------------------------------------------------------------------

def channels(images: dict[str, np.ndarray], weights: dict[str, float]) -> tuple[np.ndarray, ...]:
    """(R, G, B) display channels: each filter's image times its white-balance
    weight, placed in the channel ``CHANNEL_OF_BAND`` assigns it."""
    out = {}
    for band in BANDS:
        out[CHANNEL_OF_BAND[band]] = images[band] * weights[band]
    return out["R"], out["G"], out["B"]


def composite(images: dict[str, np.ndarray], weights: dict[str, float], *,
              Q: float, saturation: float, minimum: float = 0.0) -> np.ndarray:
    """The project's standard picture: white-balanced channels through Lupton."""
    R, G, B = channels(images, weights)
    return lupton_rgb(R, G, B, softening=softening_for(Q, saturation),
                      saturation=saturation, minimum=minimum)


def weights_that_whiten(color_mags: dict[str, float]) -> dict[str, float]:
    """Weights that render a source with these AB magnitudes as neutral grey
    (normalised to r). ``color_mags`` needs only relative values."""
    return {b: 10.0 ** (0.4 * (color_mags[b] - color_mags["r"])) for b in BANDS}


def neutral_color_of(weights: dict[str, float]) -> dict[str, float]:
    """The colours (g-r, r-i) a set of channel weights renders as neutral grey."""
    # neutral when w_b f_b is the same for every band: f_b / f_r = w_r / w_b
    return {"g-r": 2.5 * np.log10(weights["g"] / weights["r"]),
            "r-i": 2.5 * np.log10(weights["r"] / weights["i"])}


def displayed_color(rgb: np.ndarray, weights: dict[str, float]) -> dict[str, np.ndarray]:
    """Invert a displayed pixel back to the (g-r, r-i) its channel ratios imply.

    For a colour-preserving stretch this returns the pixel's true colours; for a
    per-channel stretch it shows how far the picture has drifted from them.
    """
    R, G, B = rgb[..., 0], rgb[..., 1], rgb[..., 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        gr = -2.5 * np.log10((B / weights["g"]) / (G / weights["r"]))
        ri = -2.5 * np.log10((G / weights["r"]) / (R / weights["i"]))
    return {"g-r": gr, "r-i": ri}


# --------------------------------------------------------------------------
# The survey's own recipe (for validation)
# --------------------------------------------------------------------------

def viewer_sdss_rgb(images_native: dict[str, np.ndarray], *, scales: dict, m: float, Q: float) -> np.ndarray:
    """legacypipe's ``sdss_rgb``, as the Legacy Surveys viewer renders its SDSS layer.

    Inputs are the raw cutout values (nanomaggies per native pixel). Each band
    goes to a fixed plane with a fixed scale; I is the mean of
    max(0, scale x image + m); every plane is multiplied by asinh(Q I) /
    (sqrt(Q) I); then everything is clipped to [0, 1] -- which saturates bright
    cores to white, unlike Lupton's original.
    """
    bands = list(images_native)
    intensity = 0.0
    for b in bands:
        plane, scale = scales[b]
        intensity = intensity + np.maximum(0.0, images_native[b] * scale + m)
    intensity = intensity / len(bands)
    f = np.arcsinh(Q * intensity) / np.sqrt(Q)
    intensity = intensity + (intensity == 0.0) * 1e-6
    rgb = np.zeros(intensity.shape + (3,), np.float32)
    for b in bands:
        plane, scale = scales[b]
        rgb[..., plane] = (images_native[b] * scale + m) * f / intensity
    return np.clip(rgb, 0.0, 1.0)


# --------------------------------------------------------------------------
# What the eye would see
# --------------------------------------------------------------------------

def true_color_rgb(images: dict[str, np.ndarray], matrix: cp.ColorMatrix, *,
                   white_xyz: np.ndarray, reference_xyz_per_intensity: float,
                   Q: float, saturation: float) -> tuple[np.ndarray, dict]:
    """An estimate of the colours a human eye would see, if the galaxy were bright
    enough to switch on colour vision.

    Per pixel: the three fluxes -> CIE XYZ through the colour matrix -> adapted
    so that sunlight is white -> linear sRGB. Brightness is then mapped with the
    *same* arcsinh curve as the survey-style composite, applied to the eye's
    luminance Y and to encoded (on-screen) lightness, so a neutral pixel lands on
    exactly the grey level the standard composite gives it: the two pictures
    differ only in colour, never in brightness.

    ``reference_xyz_per_intensity`` is the Y of a source with the Sun's colours
    and one nanomaggy per arcsec^2 in r: it puts Y on the same scale as the
    composite's intensity I.

    Returns the image and a few diagnostics (fraction of out-of-gamut pixels).
    """
    fluxes = np.stack([images[b] for b in matrix.bands], axis=-1) * 1e-9  # maggies/arcsec^2
    xyz = matrix.xyz(fluxes)
    Y = xyz[..., 1] / reference_xyz_per_intensity
    rgb_lin = cp.xyz_to_display_rgb(xyz / reference_xyz_per_intensity, adapt_from=white_xyz)

    # out-of-gamut (a negative channel): desaturate towards the grey of equal Y
    lo = rgb_lin.min(axis=-1)
    oog = (lo < 0) & (Y > 0)
    t = np.where(oog, Y / np.where(oog, Y - lo, 1.0), 1.0)
    rgb_lin = Y[..., None] + t[..., None] * (rgb_lin - Y[..., None])

    target = cp.srgb_decode(asinh_curve(np.clip(Y, 0, None), softening_for(Q, saturation), saturation))
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(Y > 0, target / Y, 0.0)
    rgb_lin = np.clip(rgb_lin * scale[..., None], 0.0, None)
    peak = rgb_lin.max(axis=-1)
    with np.errstate(divide="ignore", invalid="ignore"):
        rgb_lin = np.where((peak > 1)[..., None], rgb_lin / peak[..., None], rgb_lin)
    # out-of-gamut pixels are almost all faint, noise-dominated ones; report the
    # share among pixels brighter than 1 nanomaggy/arcsec^2 (mu ~ 22.5), where
    # the noise is ~10% and a pixel's colour means something
    bright = Y > 1.0
    diag = {"out_of_gamut_fraction_bright": float(oog[bright].mean()) if bright.any() else 0.0,
            "out_of_gamut_fraction_all": float(oog[Y > 0].mean()) if (Y > 0).any() else 0.0}
    return cp.srgb_encode(rgb_lin), diag


# --------------------------------------------------------------------------
# Image files
# --------------------------------------------------------------------------

def to_uint8(rgb: np.ndarray, *, flip: bool = True, mode: str = "round") -> np.ndarray:
    """Float [0, 1] -> 8-bit, flipped so north is up in an image file."""
    v = np.clip(rgb, 0.0, 1.0) * 255.0
    v = np.floor(v) if mode == "floor" else np.round(v)
    v = v.astype(np.uint8)
    return v[::-1] if flip else v


def save_png(rgb: np.ndarray, path, *, flip: bool = True) -> str:
    Image.fromarray(to_uint8(rgb, flip=flip)).save(path, optimize=True)
    return str(path)


def save_jpeg(rgb: np.ndarray, path, *, flip: bool = True, quality: int = 90) -> str:
    """Photographs compress poorly without loss because the sky noise is real;
    quality 90 JPEG keeps them near 1 MB with no visible artefacts."""
    Image.fromarray(to_uint8(rgb, flip=flip)).save(path, format="JPEG", quality=quality, optimize=True)
    return str(path)


def jpeg_roundtrip(img_uint8: np.ndarray, quality: int) -> np.ndarray:
    """Compress and decompress an 8-bit image: what JPEG alone does to it."""
    buf = io.BytesIO()
    Image.fromarray(img_uint8).save(buf, format="JPEG", quality=quality)
    buf.seek(0)
    return np.asarray(Image.open(buf).convert("RGB"))
