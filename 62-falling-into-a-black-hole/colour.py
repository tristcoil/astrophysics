"""Colours of starlight, and what a Doppler or gravitational shift does to them.

Every source on our sky (a star, or the summed light of the faint stars in one
patch of sky) is rendered as a blackbody of some temperature T, scaled to its
measured brightness. The reason is physics, not convenience: the specific
intensity I_nu / nu^3 is the same for every observer along a ray, so light
whose frequencies are all multiplied by g has I_nu(nu) = g^3 I_nu,0(nu / g).
For a blackbody that is exactly a blackbody at temperature g T. A shift turns
a 5,800 K Sun-like star into an 11,600 K blue-white one (g = 2) or a 2,900 K
red one (g = 1/2), and its total (bolometric) brightness by g^4.

* T comes from each star's colour through Mamajek's dwarf sequence (BP-RP for
  Gaia, B-V for Hipparcos). The observed colour includes interstellar
  reddening, which is what we want: the render shows the sky as it looks.
* Brightness: V magnitudes (Hipparcos directly; Gaia G converted with the DR3
  G - V(BP - RP) polynomial).
* Colour on screen: the blackbody spectrum through the CIE 1931 2-degree
  colour-matching functions, to linear sRGB (D65 white, so a 6,500 K body is
  white and the Sun a faint yellow-white).

Real stellar spectra are not blackbodies (absorption lines, the Balmer jump,
ultraviolet line blanketing), so strongly shifted colours are approximate: a
caveat recorded in the README.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

import config as C

H = 6.62607015e-34
KB = 1.380649e-23
CL = 299_792_458.0

#: XYZ -> linear sRGB (IEC 61966-2-1, D65).
XYZ_TO_RGB = np.array([[3.2406, -1.5372, -0.4986],
                       [-0.9689, 1.8758, 0.0415],
                       [0.0557, -0.2040, 1.0570]])

LOGT_MIN, LOGT_MAX, LOGT_N = 1.8, 9.5, 4000


@lru_cache(maxsize=1)
def cmf():
    data = np.loadtxt(C.REF_DIR / "ciexyz31_1.csv", delimiter=",")
    return data[:, 0] * 1e-9, data[:, 1:]          # wavelength (m), xbar ybar zbar


def _log_planck_lambda(lam, T):
    """ln B_lambda(T) up to a constant, robust for small T (Wien) and huge T (Rayleigh-Jeans)."""
    x = H * CL / (lam[None, :] * KB * T[:, None])
    # ln(1/(e^x - 1)) = -x - ln(1 - e^-x); for small x use -ln(expm1(x))
    lnfac = np.where(x > 1e-3, -x - np.log1p(-np.exp(-np.minimum(x, 700.0))), -np.log(np.expm1(np.maximum(x, 1e-300))))
    return -5.0 * np.log(lam)[None, :] + lnfac


@lru_cache(maxsize=1)
def tables():
    """log10 T grid, ln Y(T) (luminance) and chromaticity rgb(T) / Y(T) in linear sRGB."""
    lam, xyz = cmf()
    logT = np.linspace(LOGT_MIN, LOGT_MAX, LOGT_N)
    T = 10.0**logT
    lnB = _log_planck_lambda(lam, T)
    # integrate in a numerically safe way: factor out the max of ln B for each T
    m = lnB.max(axis=1, keepdims=True)
    w = np.exp(lnB - m)
    dl = np.gradient(lam)
    XYZ = (w * dl[None, :]) @ xyz                   # times exp(m)
    lnY = np.log(XYZ[:, 1]) + m[:, 0]
    rgb = XYZ @ XYZ_TO_RGB.T
    rgb = np.clip(rgb, 0.0, None)                    # very red blackbodies fall outside sRGB
    chroma = rgb / XYZ[:, 1:2]
    return logT, lnY, chroma


def uinterp(x, x0: float, dx: float, y):
    """Linear interpolation on a uniform grid x0 + i dx (clamped at the ends): np.interp
    without the binary search, several times faster for millions of points."""
    n = y.shape[0]
    t = (np.asarray(x, float) - x0) / dx
    t = np.clip(t, 0.0, n - 1.000001)
    i = t.astype(np.int64)
    f = t - i
    if y.ndim == 1:
        return y[i] * (1.0 - f) + y[i + 1] * f
    return y[i] * (1.0 - f)[..., None] + y[i + 1] * f[..., None]


def shifted_rgb(lum, logT, g):
    """Linear sRGB of sources with unshifted luminance ``lum`` and temperature 10**logT, shifted by g.

    lum * [Y(gT) / Y(T)] * chroma(gT). Works on arrays of any (broadcastable) shape.
    """
    grid, lnY, chroma = tables()
    x0, dx = grid[0], grid[1] - grid[0]
    lt = np.asarray(logT, float)
    lts = lt + np.log10(np.asarray(g, float))
    ratio = np.exp(uinterp(lts, x0, dx, lnY) - uinterp(lt, x0, dx, lnY))
    ratio = np.where(lts < grid[0], 0.0, ratio)      # redshifted out of sight
    out = uinterp(lts, x0, dx, chroma)
    return out * (np.asarray(lum, float) * ratio)[..., None]


def luminance_ratio(logT, g):
    """Y(gT)/Y(T): how much brighter a source looks to the eye after a shift by g."""
    grid, lnY, _ = tables()
    x0, dx = grid[0], grid[1] - grid[0]
    lt = np.asarray(logT, float)
    return np.exp(uinterp(lt + np.log10(g), x0, dx, lnY) - uinterp(lt, x0, dx, lnY))


# ---------------------------------------------------------------------------
# Temperatures from colours (Mamajek's dwarf sequence)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def mamajek():
    """(Teff, B-V, BP-RP) from Mamajek's 'A Modern Mean Dwarf Stellar Color and Effective
    Temperature Sequence' (the first table block of EEM_dwarf_UBVIJHK_colors_Teff.txt)."""
    rows = []
    started = False
    for line in (C.REF_DIR / "EEM_dwarf_UBVIJHK_colors_Teff.txt").read_text().splitlines():
        if line.startswith("#SpT"):
            if started:
                break
            started = True
            continue
        if not started or not line.strip() or line.startswith("#"):
            continue
        parts = line.split()

        def num(i):
            try:
                return float(parts[i])
            except (ValueError, IndexError):
                return math.nan
        rows.append((num(1), num(8), num(11)))
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def _bprp_extended():
    """Mamajek's BP-RP column, extended to O and early-B stars (where it is empty).

    Over B9V-A9V both colours exist and BP-RP is linear in B-V; that fit
    converts the hot stars' B-V into BP-RP. Hot stars brighter than V = 5 use
    Hipparcos B-V directly, so this only matters for fainter Gaia stars.
    """
    teff, bv, bprp = mamajek()
    both = np.isfinite(bv) & np.isfinite(bprp) & (teff >= 7400) & (teff <= 10700)
    slope, icpt = np.polyfit(bv[both], bprp[both], 1)
    fill = ~np.isfinite(bprp) & np.isfinite(bv) & (teff > 10700)
    out = bprp.copy()
    out[fill] = icpt + slope * bv[fill]
    return out


def _colour_to_logT(colour, which: int):
    teff, bv, _ = mamajek()
    col = (bv, _bprp_extended())[which]
    ok = np.isfinite(col) & np.isfinite(teff)
    x, y = col[ok], np.log10(teff[ok])
    order = np.argsort(x)
    x, y = x[order], y[order]
    # enforce monotonic (the sequence is monotonic in colour except for ties)
    keep = np.concatenate([[True], np.diff(x) > 0])
    return np.interp(np.asarray(colour, float), x[keep], y[keep])


def logT_from_bprp(bprp):
    return _colour_to_logT(bprp, 1)


def logT_from_bv(bv):
    return _colour_to_logT(bv, 0)


def v_from_gaia(g, bprp):
    """V magnitude from Gaia G and BP-RP (DR3 polynomial, colour clamped to its valid range)."""
    x = np.clip(np.asarray(bprp, float), -0.5, 2.75)
    c0, c1, c2, c3 = C.G_MINUS_V
    return np.asarray(g, float) - (c0 + c1 * x + c2 * x**2 + c3 * x**3)


def blackbody_srgb_swatch(T):
    """Display sRGB (0-1, gamma-encoded) of a blackbody at temperature T, normalised to max channel 1."""
    _, _, chroma = tables()
    grid = tables()[0]
    lt = np.log10(T)
    rgb = np.array([np.interp(lt, grid, chroma[:, k]) for k in range(3)])
    rgb = rgb / rgb.max()
    return np.where(rgb <= 0.0031308, 12.92 * rgb, 1.055 * rgb ** (1 / 2.4) - 0.055)
