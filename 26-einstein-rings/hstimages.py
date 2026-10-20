"""Hubble images for display: colour composites, and the lensed galaxy un-lensed.

Colour. Each filter's count rate is converted to flux per unit frequency
(f_nu, from the header's PHOTFLAM and PHOTPLAM), so that a source with a flat
f_nu spectrum would look white. Red is F814W (I band), blue is the bluest deep
filter available, green their mean. Old stars in the elliptical then come out
yellow-orange and the young stars of the lensed galaxy blue: the colour
difference is physical, not a choice of palette. The stretch is the arcsinh
mapping of Lupton et al. (2004), identical for every lens.

Un-lensing. Every pixel of the lens-subtracted image is traced back through our
best-fit SIE (beta = theta - alpha(theta)) and its flux dropped into a grid in
the source plane. Where several images of the same spot exist, they are
averaged. This shows the background galaxy as the lens model says it is,
pixel data and all, rather than only its fitted Sersic profile.
"""

from __future__ import annotations

import json

import numpy as np
from astropy.io import fits

import lensing as L
from config import HST_DIR, HST_PIXEL_ARCSEC, MODEL_DIR

BLUE_PREFERENCE = ["wfc3_uvis_f438w", "wfc3_uvis_f390w", "acs_wfc_f435w", "acs_wfc_f555w",
                   "wfc3_uvis_f336w", "acs_wfc_f606w"]
C_AA = 2.99792458e18  # speed of light in Angstrom/s


def _fnu(path):
    """Image in f_nu units (arbitrary but common scale): e/s * PHOTFLAM * PHOTPLAM^2 / c,
    with the pivot wavelength returned too."""
    with fits.open(path) as hdul:
        sci = np.nan_to_num(hdul["SCI"].data.astype(float))
        h = hdul["SCI"].header
        return sci * h["PHOTFLAM"] * h["PHOTPLAM"] ** 2 / C_AA * 1e29, float(h["PHOTPLAM"])


def _clean(img):
    """Subtract the sky (median of the outer frame) and replace cosmic-ray hits
    (pixels standing far above their 3x3 neighbourhood, much sharper than a
    star) with the local median. Display only; the lens fits do their own masking."""
    from scipy import ndimage

    edge = np.concatenate([img[:12].ravel(), img[-12:].ravel(), img[:, :12].ravel(), img[:, -12:].ravel()])
    sky = np.median(edge)
    noise = 1.4826 * np.median(np.abs(edge - sky))
    img = img - sky
    med3 = ndimage.median_filter(img, size=3)
    med7 = ndimage.median_filter(img, size=7)
    hit = ((img - med3) > 6 * noise) & ((img - med3) > 0.6 * np.maximum(img - med7, 1e-30))
    hit = ndimage.binary_dilation(hit, iterations=1)
    return np.where(hit, med3, img)


def _sun_fnu(lam_angstrom):
    """f_nu of a 5772-K blackbody, standing in for the Sun's colour."""
    h_over_k = 1.4387769e8  # hc/k in Angstrom K
    nu3 = (1.0 / lam_angstrom) ** 3
    return nu3 / np.expm1(h_over_k / (lam_angstrom * 5772.0))


def blue_band(name: str) -> str | None:
    """The bluest filter built from at least two exposures, so that cosmic rays
    were rejected. SLACS's F435W images are mostly single 420-s snapshots, full
    of them; a lens with nothing better is shown in grey instead."""
    fetched = json.loads((HST_DIR / name / "_fetched.json").read_text())["filters"]
    for f in BLUE_PREFERENCE:
        if f in fetched:
            ncomb = fits.getheader(HST_DIR / name / f"{f}.fits", "SCI").get("NCOMBINE", 1)
            if ncomb >= 2:
                return f
    return None


def composite(name: str, centre_pix: tuple[int, int], half: int, *, q: float = 6.0,
              stretch_frac: float = 0.05):
    """RGB array (half*2 x half*2 x 3) around the lens, or None without a blue band.

    White is the Sun's colour: each band is divided by a solar spectrum's f_nu
    there, so a star like the Sun would come out white, an old elliptical
    yellow-orange, and young stars blue."""
    from astropy.visualization import make_lupton_rgb

    blue = blue_band(name)
    if blue is None:
        return None, None
    r, lam_r = _fnu(HST_DIR / name / "acs_wfc_f814w.fits")
    b, lam_b = _fnu(HST_DIR / name / f"{blue}.fits")
    if r.shape != b.shape:
        return None, None
    cx, cy = centre_pix
    pad = 12
    sl = (slice(cy - half - pad, cy + half + pad), slice(cx - half - pad, cx + half + pad))
    r = _clean(r[sl])[pad:-pad, pad:-pad] / _sun_fnu(lam_r)
    b = _clean(b[sl])[pad:-pad, pad:-pad] / _sun_fnu(lam_b)
    g = 0.5 * (r + b)
    peak = np.percentile(r, 99.9)
    rgb = make_lupton_rgb(r, g, b, minimum=0, stretch=stretch_frac * peak, Q=q)
    return rgb, blue


def grey(name: str, centre_pix, half: int):
    r, _ = _fnu(HST_DIR / name / "acs_wfc_f814w.fits")
    cx, cy = centre_pix
    pad = 12
    img = r[cy - half - pad:cy + half + pad, cx - half - pad:cx + half + pad]
    return _clean(img)[pad:-pad, pad:-pad]


def load_model(name: str, variant: str = ""):
    from lensmodel import model_path
    path = model_path(name, variant, ".npz")
    if not path.exists():
        return None, None
    return np.load(path), json.loads(model_path(name, variant).read_text())


def model_geometry(npz):
    """Lens centre, Einstein radius, mass axis ratio and angle (image frame) of a fit."""
    from lensmodel import MASS_E_MAX, e_to_qphi

    p = npz["params"]
    x0, y0 = npz["light_shape"][:2]
    q, phi = e_to_qphi(p[1], p[2], MASS_E_MAX)
    return {"x0": float(x0), "y0": float(y0), "b": float(p[0]), "q": float(q), "phi": float(phi),
            "xs": float(p[3]), "ys": float(p[4]), "r_s": float(10 ** p[5]), "n_s": float(p[6])}


def coordinates(shape):
    n = shape[0]
    idx = (np.arange(n) - n // 2) * HST_PIXEL_ARCSEC
    return np.meshgrid(idx, idx)


def unlens(npz, half_width: float = 0.6, n: int = 70, min_snr: float = 1.0, sub: int = 3):
    """Back-trace the lens-subtracted image to the source plane.

    Each image pixel is split into sub x sub points that carry its surface
    brightness (lensing conserves surface brightness), traced back through the
    lens and averaged into a grid centred on the lens axis; this fills the grid
    evenly where the lens stretches the image. Only pixels inside the fit
    region and above min_snr contribute; empty cells are NaN."""
    geo = model_geometry(npz)
    data, light, var = npz["data"], npz["light"], npz["var"]
    resid = data - light
    x, y = coordinates(data.shape)
    sel = npz["fit_mask"] & (resid / np.sqrt(var) > min_snr)
    off = ((np.arange(sub) + 0.5) / sub - 0.5) * HST_PIXEL_ARCSEC
    ox, oy = np.meshgrid(off, off)
    px = (x[sel][:, None] + ox.ravel()[None, :]).ravel()
    py = (y[sel][:, None] + oy.ravel()[None, :]).ravel()
    val = np.repeat(resid[sel], sub * sub)
    ax, ay = L.sie_deflection_sky(px, py, geo["b"], geo["q"], geo["phi"], geo["x0"], geo["y0"])
    bx, by = px - ax - geo["x0"], py - ay - geo["y0"]
    edges = np.linspace(-half_width, half_width, n + 1)
    tot, _, _ = np.histogram2d(by, bx, bins=[edges, edges], weights=val)
    cnt, _, _ = np.histogram2d(by, bx, bins=[edges, edges])
    with np.errstate(invalid="ignore"):
        img = tot / cnt
    return img, [-half_width, half_width, -half_width, half_width]


def magnification(npz, oversample: int = 4) -> float:
    """Total magnification of the fitted source: lensed flux / intrinsic flux.

    Surface brightness is conserved by lensing, so the lensed flux is the
    source's brightness at each traced-back position summed over the image
    plane; the intrinsic flux is the Sersic profile's integral,
    2 pi n R_e^2 e^b_n Gamma(2n) / b_n^(2n) (intermediate-axis radius, so the
    axis ratio drops out)."""
    from scipy.special import gamma

    from lensmodel import e_to_qphi, sersic_source
    geo = model_geometry(npz)
    p = npz["params"]
    q_s, phi_s = e_to_qphi(p[7], p[8])
    n = npz["data"].shape[0]
    step = HST_PIXEL_ARCSEC / oversample
    half = n // 2 * HST_PIXEL_ARCSEC
    g = np.arange(-half, half, step) + step / 2
    x, y = np.meshgrid(g, g)
    ax, ay = L.sie_deflection_sky(x, y, geo["b"], geo["q"], geo["phi"], geo["x0"], geo["y0"])
    lensed = sersic_source(x - ax, y - ay, geo["xs"], geo["ys"], geo["r_s"], geo["n_s"], q_s, phi_s).sum() * step**2
    bn = L.sersic_b(geo["n_s"])
    ns = geo["n_s"]
    intrinsic = 2 * np.pi * ns * geo["r_s"] ** 2 * np.exp(bn) * gamma(2 * ns) / bn ** (2 * ns)
    return float(lensed / intrinsic)
