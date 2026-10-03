"""The physics of colour: from a spectrum, or three filter fluxes, to what an eye sees.

Two different questions hide inside "what colour is this galaxy?":

* **What does the camera record?** A survey camera measures how much light
  gets through each of its filters. That is *photometry*: a flux per filter,
  computed here from any spectrum by integrating it against the filter's
  transmission curve (``band_flux``).

* **What would a person see?** The eye has its own three "filters", the cone
  cells, whose combined response is standardised by the CIE (International
  Commission on Illumination) as the 1931 colour-matching functions. Integrating
  a spectrum against them gives tristimulus values X, Y, Z (``spectrum_xyz``):
  Y is the brightness the eye perceives, and the ratios X:Y:Z are its colour.
  Any two spectra with the same XYZ look identical to a normal human observer.

The SDSS filters are not the eye's filters -- i is almost invisible to us -- so
survey colours are not eye colours. But the two are related: for the smooth
spectra of stars and galaxies, three broad filter fluxes pin down XYZ well. We
learn the linear map between them (the same "colour matrix" every digital
camera uses) from real galaxy spectra (``fit_color_matrix``), and check how
well it predicts spectra it was not trained on.

All wavelengths are in Angstrom and flux densities in erg s^-1 cm^-2 A^-1
unless a name says otherwise.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from astropy.io import fits

from config import (
    K_M_PHOTOPIC,
    K_M_SCOTOPIC,
    NANOMAGGY_JY,
    NATIVE_PIXSCALE,
    PAPER_REFLECTANCE,
    PICKLES_STARS,
    REFERENCE_DIR,
)

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

C_ANGSTROM_S = 2.99792458e18  # speed of light, Angstrom/s
AB_ZERO_JY = 3631.0  # the AB magnitude system's zero point
JY_CGS = 1e-23  # erg s^-1 cm^-2 Hz^-1
AU_CM = 1.495978707e13
PC_CM = 3.0856775814913673e18
ARCSEC2_SR = (np.pi / 180.0 / 3600.0) ** 2  # one square arcsecond in steradians
ERG_S_CM2_IN_W_M2 = 1e-3


# --------------------------------------------------------------------------
# Reference data
# --------------------------------------------------------------------------

@lru_cache(maxsize=None)
def load_filter(name: str) -> tuple[np.ndarray, np.ndarray]:
    """A filter's total response curve: (wavelength A, response)."""
    d = np.loadtxt(REFERENCE_DIR / f"filter_{name}.dat")
    order = np.argsort(d[:, 0])
    return d[order, 0], np.clip(d[order, 1], 0.0, None)


@lru_cache(maxsize=None)
def load_cmf() -> tuple[np.ndarray, np.ndarray]:
    """CIE 1931 2-degree colour-matching functions: (wavelength A, (N, 3) xbar ybar zbar)."""
    d = np.genfromtxt(REFERENCE_DIR / "CIE_xyz_1931_2deg.csv", delimiter=",")
    return d[:, 0] * 10.0, d[:, 1:4]


@lru_cache(maxsize=None)
def load_scotopic() -> tuple[np.ndarray, np.ndarray]:
    """CIE 1951 scotopic (rod) luminous efficiency V'(lambda)."""
    d = np.genfromtxt(REFERENCE_DIR / "CIE_sle_scotopic.csv", delimiter=",")
    return d[:, 0] * 10.0, np.nan_to_num(d[:, 1])


@lru_cache(maxsize=None)
def load_d65() -> tuple[np.ndarray, np.ndarray]:
    """CIE standard illuminant D65, relative spectral power."""
    d = np.genfromtxt(REFERENCE_DIR / "CIE_std_illum_D65.csv", delimiter=",")
    return d[:, 0] * 10.0, d[:, 1]


def _read_spectrum_fits(path) -> tuple[np.ndarray, np.ndarray]:
    t = fits.getdata(path, 1)
    names = [n.upper() for n in t.columns.names]
    lam = np.asarray(t.field(names.index("WAVELENGTH")), float)
    flux = np.asarray(t.field(names.index("FLUX")), float)
    order = np.argsort(lam)
    return lam[order], flux[order]


@lru_cache(maxsize=None)
def load_sun() -> tuple[np.ndarray, np.ndarray]:
    """The Sun's spectrum at 1 AU, above the atmosphere (CALSPEC)."""
    return _read_spectrum_fits(REFERENCE_DIR / "sun_reference_stis_002.fits")


@lru_cache(maxsize=None)
def load_vega() -> tuple[np.ndarray, np.ndarray]:
    return _read_spectrum_fits(REFERENCE_DIR / "alpha_lyr_stis_012.fits")


@lru_cache(maxsize=None)
def load_pickles(index: int) -> tuple[np.ndarray, np.ndarray]:
    """One Pickles (1998) library star (arbitrary flux scale)."""
    return _read_spectrum_fits(REFERENCE_DIR / f"pickles_uk_{index}.fits")


def pickles_stars() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {sp: load_pickles(i) for i, sp in PICKLES_STARS.items()}


@dataclass
class BrownAtlas:
    """The Brown et al. (2014) galaxy spectra and their published photometry."""

    #: name -> (observed wavelength A, f_lambda, source flag per point:
    #: 0 = model fill, 1 = optical spectrum, 2 = Akari, 3 = Spitzer IRS)
    spectra: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]
    #: name -> {band: AB mag} matched-aperture photometry (SDSS u g r i z)
    photometry: dict[str, dict[str, float]]
    #: name -> {band: foreground Galactic extinction, mag}
    extinction: dict[str, dict[str, float]]
    #: name -> morphology string
    morphology: dict[str, str]
    #: name -> (RA, Dec) in degrees, the centre of their apertures
    position: dict[str, tuple[float, float]]


@lru_cache(maxsize=None)
def load_brown() -> BrownAtlas:
    path = REFERENCE_DIR / "brown2014_galsedatlas.fits"
    spectra, phot, ext, morph, pos = {}, {}, {}, {}, {}
    with fits.open(path) as hdul:
        for hdu in hdul[1:]:
            name = hdu.name
            if name.upper().endswith("_SPEC"):
                d = hdu.data
                key = name[: -len("_spec")]
                lam = np.asarray(d["observed_wavelength"], float)
                order = np.argsort(lam)
                spectra[key] = (lam[order], np.asarray(d["flux_density"], float)[order],
                                np.asarray(d["source"], int)[order])
        for row in hdul["photometry"].data:
            key = row["Name"].strip().replace(" ", "_")
            phot[key] = {b: float(row[f"SDSS_{b}"]) for b in "ugriz"}
        for row in hdul["extinction"].data:
            key = row["Name"].strip().replace(" ", "_")
            ext[key] = {b: float(row[f"SDSS_{b}_a"]) for b in "ugriz"}
        for row in hdul["summary"].data:
            key = row["Name"].strip().replace(" ", "_")
            morph[key] = row["Morphology"].strip()
            pos[key] = (float(row["RA(J2000)"]), float(row["Decl(J2000)"]))
    return BrownAtlas(spectra, phot, ext, morph, pos)


# --------------------------------------------------------------------------
# Photometry: what the camera measures
# --------------------------------------------------------------------------

def band_flux(lam: np.ndarray, flam: np.ndarray, filter_name: str) -> float:
    """Mean flux density through a filter, in maggies (1 maggie = 3631 Jy).

    The SDSS camera counts photons, so each wavelength is weighted by
    lambda x response; this is the standard AB synthetic-photometry integral
        <f_nu> = int f_lambda R lambda dlambda / int R (c / lambda) dlambda.
    """
    fl, fr = load_filter(filter_name)
    resp = np.interp(lam, fl, fr, left=0.0, right=0.0)
    if lam.min() > fl[fr > 0.001 * fr.max()].min() or lam.max() < fl[fr > 0.001 * fr.max()].max():
        raise ValueError(f"spectrum does not cover filter {filter_name}")
    num = np.trapezoid(flam * resp * lam, lam)
    den = np.trapezoid(resp * C_ANGSTROM_S / lam, lam)
    return num / den / (AB_ZERO_JY * JY_CGS)


def ab_mag(lam, flam, filter_name: str) -> float:
    return -2.5 * np.log10(band_flux(lam, flam, filter_name))


def pivot_wavelength(filter_name: str) -> float:
    """The wavelength a filter "is at", for a photon-counting detector."""
    fl, fr = load_filter(filter_name)
    return float(np.sqrt(np.trapezoid(fr * fl, fl) / np.trapezoid(fr / fl, fl)))


def sun_absolute_magnitudes(bands=("u", "g", "r", "i", "z")) -> dict[str, float]:
    """The Sun's absolute AB magnitude in each SDSS band, from its measured spectrum."""
    lam, flam = load_sun()
    at_10pc = flam * (AU_CM / (10.0 * PC_CM)) ** 2
    return {b: ab_mag(lam, at_10pc, f"sdss_{b}") for b in bands}


def solar_white_weights(bands=("g", "r", "i")) -> dict[str, float]:
    """Per-band multipliers that make the Sun's colour exactly neutral grey.

    A pixel's display channel = weight x flux. Choosing weight_b proportional to
    1 / (the Sun's flux in b) makes the three channels equal for any source with
    solar colours. Normalised so the middle band has weight 1.
    """
    lam, flam = load_sun()
    f = {b: band_flux(lam, flam, f"sdss_{b}") for b in bands}
    mid = bands[len(bands) // 2]
    return {b: f[mid] / f[b] for b in bands}


# --------------------------------------------------------------------------
# Colorimetry: what the eye sees
# --------------------------------------------------------------------------

def spectrum_xyz(lam: np.ndarray, flam: np.ndarray) -> np.ndarray:
    """CIE 1931 tristimulus values X, Y, Z of a spectrum (units of flam x A).

    The spectrum is interpolated onto the CIE table's 1 nm grid, the sampling
    the CIE recommends for the integral.
    """
    cl, cmf = load_cmf()
    if lam.min() > cl.min() + 200 or lam.max() < 7800:
        raise ValueError("spectrum does not cover the visible range")
    f = np.interp(cl, lam, flam, left=0.0, right=0.0)
    return np.trapezoid(f[:, None] * cmf, cl, axis=0)


def spectrum_scotopic(lam: np.ndarray, flam: np.ndarray) -> float:
    """The rod-weighted integral of a spectrum (units of flam x A)."""
    sl, sv = load_scotopic()
    f = np.interp(sl, lam, flam, left=0.0, right=0.0)
    return float(np.trapezoid(f * sv, sl))


def xyz_to_xy(xyz: np.ndarray) -> np.ndarray:
    """Chromaticity (x, y): colour with brightness divided out."""
    xyz = np.asarray(xyz, float)
    s = xyz.sum(axis=-1, keepdims=True)
    return (xyz / s)[..., :2]


def xyz_to_upvp(xyz: np.ndarray) -> np.ndarray:
    """CIE 1976 u', v' chromaticity: a near-uniform colour space in which a
    distance of ~0.004 is roughly one just-noticeable difference in hue."""
    xyz = np.asarray(xyz, float)
    X, Y, Z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    d = X + 15.0 * Y + 3.0 * Z
    return np.stack([4.0 * X / d, 9.0 * Y / d], axis=-1)


def delta_upvp(xyz_a, xyz_b) -> np.ndarray:
    return np.hypot(*np.moveaxis(xyz_to_upvp(xyz_a) - xyz_to_upvp(xyz_b), -1, 0))


def xy_to_xyz(x: float, y: float, Y: float = 1.0) -> np.ndarray:
    return np.array([x * Y / y, Y, (1.0 - x - y) * Y / y])


# --- blackbodies and colour temperature --------------------------------------

H_CGS = 6.62607015e-27
K_CGS = 1.380649e-16


def blackbody_flam(lam: np.ndarray, temperature: float) -> np.ndarray:
    """Planck's law, B_lambda (arbitrary overall scale; shape is what matters)."""
    lam_cm = lam * 1e-8
    x = H_CGS * (C_ANGSTROM_S * 1e-8) / (lam_cm * K_CGS * temperature)
    return 1.0 / (lam_cm**5 * np.expm1(np.clip(x, 1e-8, 700)))


@lru_cache(maxsize=None)
def planckian_locus(t_min: float = 1000.0, t_max: float = 100_000.0, n: int = 3000):
    """(temperatures, CIE 1960 (u, v)) of blackbodies, computed through the CMFs."""
    cl, _ = load_cmf()
    temps = np.geomspace(t_min, t_max, n)
    uv = np.empty((n, 2))
    for k, t in enumerate(temps):
        uv[k] = xyz_to_uv1960(spectrum_xyz(cl, blackbody_flam(cl, t)))
    return temps, uv


def xyz_to_uv1960(xyz) -> np.ndarray:
    xyz = np.asarray(xyz, float)
    X, Y, Z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    d = X + 15.0 * Y + 3.0 * Z
    return np.stack([4.0 * X / d, 6.0 * Y / d], axis=-1)


def correlated_color_temperature(xyz) -> tuple[float, float]:
    """(CCT in K, Duv): the temperature of the blackbody whose colour is closest.

    Found by direct search along the Planckian locus in the CIE 1960 (u, v)
    diagram -- the definition of CCT, with no approximation formula. Duv is the
    distance from the locus: positive above it (greenish), negative below
    (pinkish). A light source counts as "white" roughly for |Duv| < 0.02.
    """
    temps, luv = planckian_locus()
    p = xyz_to_uv1960(xyz)
    d = np.hypot(luv[:, 0] - p[0], luv[:, 1] - p[1])
    k = int(np.argmin(d))
    # refine with a parabola through the three nearest points (in log T)
    if 0 < k < len(temps) - 1:
        lt = np.log(temps[k - 1:k + 2])
        a, b, _ = np.polyfit(lt, d[k - 1:k + 2] ** 2, 2)
        t_best = float(np.exp(-b / (2 * a))) if a > 0 else float(temps[k])
    else:
        t_best = float(temps[k])
    # Sign convention (ANSI C78.377): positive above the locus, i.e. towards
    # larger v. At the nearest locus point the offset is perpendicular to the
    # locus, which for T > ~1500 K never runs vertically, so the sign of the v
    # offset is the side.
    side = 1.0 if p[1] >= luv[k, 1] else -1.0
    return t_best, float(side * d[k])


# --- sRGB --------------------------------------------------------------------

#: sRGB (IEC 61966-2-1) primaries and white point, as chromaticities.
SRGB_PRIMARIES = {"R": (0.64, 0.33), "G": (0.30, 0.60), "B": (0.15, 0.06)}
D65_XY = (0.3127, 0.3290)


def _rgb_to_xyz_matrix(primaries: dict, white_xy: tuple[float, float]) -> np.ndarray:
    """The 3x3 matrix taking linear RGB to XYZ, built from the primaries' and
    white point's chromaticities (so RGB = (1,1,1) is the white with Y = 1)."""
    P = np.column_stack([xy_to_xyz(*primaries[c]) for c in "RGB"])
    S = np.linalg.solve(P, xy_to_xyz(*white_xy))
    return P * S


M_SRGB_TO_XYZ = _rgb_to_xyz_matrix(SRGB_PRIMARIES, D65_XY)
M_XYZ_TO_SRGB = np.linalg.inv(M_SRGB_TO_XYZ)

#: The Bradford cone-response matrix for chromatic adaptation (Lam 1985), the
#: transform used by ICC colour management.
M_BRADFORD = np.array([[0.8951, 0.2664, -0.1614],
                       [-0.7502, 1.7135, 0.0367],
                       [0.0389, -0.0685, 1.0296]])


def chromatic_adaptation(src_white_xyz, dst_white_xyz) -> np.ndarray:
    """Matrix mapping XYZ seen under one white to how it looks adapted to another.

    The eye rescales its three cone sensitivities so that the prevailing light
    looks white. Bradford does that rescaling in an approximate cone space.
    """
    s = M_BRADFORD @ np.asarray(src_white_xyz, float)
    d = M_BRADFORD @ np.asarray(dst_white_xyz, float)
    return np.linalg.inv(M_BRADFORD) @ np.diag(d / s) @ M_BRADFORD


def srgb_encode(linear: np.ndarray) -> np.ndarray:
    """The sRGB transfer curve (display "gamma"), for linear values in [0, 1]."""
    v = np.clip(linear, 0.0, 1.0)
    return np.where(v <= 0.0031308, 12.92 * v, 1.055 * np.power(v, 1 / 2.4) - 0.055)


def srgb_decode(encoded: np.ndarray) -> np.ndarray:
    v = np.asarray(encoded, float)
    return np.where(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055) ** 2.4)


def sun_white_xyz() -> np.ndarray:
    """The Sun's colour (XYZ normalised to Y = 1): what an eye adapted to
    daylight above the atmosphere calls white."""
    xyz = spectrum_xyz(*load_sun())
    return xyz / xyz[1]


def d65_white_xyz() -> np.ndarray:
    return xy_to_xyz(*D65_XY)


def xyz_to_display_rgb(xyz: np.ndarray, *, adapt_from: np.ndarray | None = None) -> np.ndarray:
    """XYZ -> linear sRGB, optionally chromatically adapting a white first.

    ``adapt_from`` is the XYZ of the light the viewer is adapted to; it is mapped
    to the display's white (D65). With ``adapt_from = sun_white_xyz()`` the Sun
    renders as pure white, which is what "white" means in this project.
    """
    xyz = np.asarray(xyz, float)
    M = M_XYZ_TO_SRGB
    if adapt_from is not None:
        M = M @ chromatic_adaptation(adapt_from, d65_white_xyz())
    return xyz @ M.T


def swatch(xyz, *, adapt_from=None) -> tuple[np.ndarray, str, bool]:
    """A displayable colour for a light, brightness normalised away.

    Returns (sRGB in [0, 1], hex string, in_gamut). The linear RGB is scaled so
    its largest channel is 1: "this colour, at a comfortable brightness". A
    colour outside what an sRGB screen can show (a negative channel -- pure
    spectral colours, for instance) gets just enough white added to bring it
    inside, the usual way of drawing a spectrum on a screen. Every galaxy and
    star in this project is inside the gamut, so this only matters for the
    rainbow strip in figure 03.
    """
    xyz = np.asarray(xyz, float)
    rgb = xyz_to_display_rgb(xyz / xyz[1], adapt_from=adapt_from)
    in_gamut = bool(np.all(rgb >= -1e-9))
    if not in_gamut:
        rgb = rgb - rgb.min()
    rgb = rgb / rgb.max()
    enc = srgb_encode(rgb)
    hexs = "#" + "".join(f"{int(round(255 * c)):02X}" for c in enc)
    return enc, hexs, in_gamut


# --------------------------------------------------------------------------
# From filter fluxes to what the eye sees: the colour matrix
# --------------------------------------------------------------------------

@dataclass
class ColorMatrix:
    """XYZ = fluxes @ matrix, for fluxes in maggies (per unit solid angle, or not)."""

    bands: tuple[str, ...]
    matrix: np.ndarray  # (n_bands, 3)
    #: Leave-one-out chromaticity errors on the training spectra, Delta u'v'.
    loo_duv: np.ndarray
    #: Leave-one-out fractional errors in the eye's brightness Y.
    loo_dY: np.ndarray
    #: The training set's names.
    names: list[str]

    def xyz(self, fluxes: np.ndarray) -> np.ndarray:
        return np.asarray(fluxes, float) @ self.matrix


def _normalised(fluxes: np.ndarray, xyz: np.ndarray):
    y = xyz[:, 1:2]
    return fluxes / y, xyz / y


def fit_color_matrix(fluxes: np.ndarray, xyz: np.ndarray, names: list[str],
                     bands=("g", "r", "i")) -> ColorMatrix:
    """Least-squares linear map from filter fluxes to XYZ ("camera characterisation").

    Every training spectrum is first scaled to Y = 1, so each counts equally
    whatever its brightness and the fit is to *colour*. Accuracy is judged by
    leave-one-out: refit without each spectrum in turn and predict it.
    """
    F, X = _normalised(np.asarray(fluxes, float), np.asarray(xyz, float))
    M, *_ = np.linalg.lstsq(F, X, rcond=None)
    loo = np.empty(len(F))
    loo_y = np.empty(len(F))
    for k in range(len(F)):
        keep = np.arange(len(F)) != k
        Mk, *_ = np.linalg.lstsq(F[keep], X[keep], rcond=None)
        pred = F[k] @ Mk
        loo[k] = float(delta_upvp(pred, X[k]))
        loo_y[k] = float(pred[1] / X[k][1] - 1.0)
    return ColorMatrix(tuple(bands), M, loo, loo_y, list(names))


# --------------------------------------------------------------------------
# Luminance: is it bright enough to see in colour?
# --------------------------------------------------------------------------

def surface_brightness_nmgy_arcsec2(pixel_values: np.ndarray) -> np.ndarray:
    """Cutout pixel values (nanomaggies per native SDSS pixel) -> nanomaggies/arcsec^2."""
    return np.asarray(pixel_values, float) / NATIVE_PIXSCALE**2


def mag_arcsec2(nmgy_arcsec2: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        return 22.5 - 2.5 * np.log10(nmgy_arcsec2)


def luminance_cd_m2(Y_per_arcsec2: np.ndarray) -> np.ndarray:
    """Photopic luminance from the CIE Y integral of a surface brightness.

    Y here is int f_lambda ybar dlambda per square arcsecond, in erg s^-1 cm^-2
    arcsec^-2 (ybar is the eye's photopic sensitivity V(lambda)). Convert to W
    m^-2 sr^-1 and multiply by 683 lm/W: candela per square metre, the unit of
    the brightness of a surface.
    """
    return K_M_PHOTOPIC * np.asarray(Y_per_arcsec2) * ERG_S_CM2_IN_W_M2 / ARCSEC2_SR


def scotopic_luminance_cd_m2(Vprime_per_arcsec2: np.ndarray) -> np.ndarray:
    return K_M_SCOTOPIC * np.asarray(Vprime_per_arcsec2) * ERG_S_CM2_IN_W_M2 / ARCSEC2_SR


def illuminance_lux(lam: np.ndarray, flam: np.ndarray) -> float:
    """Illuminance (lux) of light arriving with spectral flux density flam."""
    Y = spectrum_xyz(lam, flam)[1]
    return K_M_PHOTOPIC * Y * ERG_S_CM2_IN_W_M2


def paper_luminance(illuminance: float, reflectance: float = PAPER_REFLECTANCE) -> float:
    """Luminance of a matte white surface: reflectance x illuminance / pi."""
    return reflectance * illuminance / np.pi


def vision_regime(luminance: np.ndarray, lower: float, upper: float) -> np.ndarray:
    """0 = scotopic (rods only, no colour), 1 = mesopic, 2 = photopic."""
    L = np.asarray(luminance)
    return np.where(L < lower, 0, np.where(L < upper, 1, 2))


# --------------------------------------------------------------------------
# Magnitudes and fluxes
# --------------------------------------------------------------------------

def nmgy_to_mag(f):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 22.5 - 2.5 * np.log10(f)


def mag_to_nmgy(m):
    return 10.0 ** (-0.4 * (np.asarray(m, float) - 22.5))


def nanomaggies_to_maggies(f):
    return np.asarray(f, float) * 1e-9


def nanomaggy_in_cgs() -> float:
    """One nanomaggy as f_nu in erg s^-1 cm^-2 Hz^-1."""
    return NANOMAGGY_JY * JY_CGS
