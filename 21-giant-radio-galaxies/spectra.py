"""Spectral indices and synchrotron ages of the lobes.

Electrons in a radio lobe radiate their energy away: the more energetic ones
fastest. A population injected with a power-law spectrum (S ~ nu^alpha_inj)
therefore keeps its low-frequency spectrum but loses its high-frequency one
above a *break frequency* that falls with time. Comparing the 144-MHz LoTSS
image with the 1.4-GHz NVSS image measures how steep the spectrum has become
between the two, and through the break frequency, how long ago the electrons
were last accelerated:

    t = 1590 * B^0.5 / (B^2 + B_IC^2) / sqrt(nu_break (1 + z))   [Myr]

with fields in microgauss, nu_break in GHz (Miley 1980), and B_IC = 3.25
(1+z)^2 microgauss the field whose energy density equals the cosmic microwave
background's: in giant lobes, scattering of CMB photons, not synchrotron
radiation, drains most of the electrons' energy.

The spectral shape is the Jaffe & Perola (1973) model: one burst of injection,
pitch angles continuously re-isotropised. We compute it numerically from the
single-electron synchrotron kernel. The lobe field B is taken at minimum energy
(classical equipartition; Miley 1980, Eq. 2). Every one of these assumptions is
standard and uncertain: the ages are model ages (Harwood 2017), useful for
comparing regions and for an order of magnitude.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy import optimize

import config
import cosmology
import radiomaps as rm

NU_LO, NU_HI = config.LOTSS_FREQ_MHZ, config.NVSS_FREQ_MHZ
LOG_RATIO = np.log(NU_HI / NU_LO)


# --------------------------------------------------------------------------
# The Jaffe-Perola spectrum
# --------------------------------------------------------------------------

def synchrotron_kernel(x):
    """F(x) = x * integral_x^inf K_5/3(t) dt, the spectrum of one electron
    (x = nu / nu_critical). Approximation of Aharonian, Kelner & Prosekin
    (2010, PRD 82, 043002, Eq. D7), accurate to better than 0.2%."""
    x = np.asarray(x, float)
    x13 = np.cbrt(x)
    x23 = x13 * x13
    x43 = x23 * x23
    with np.errstate(over="ignore", under="ignore"):
        return (2.15 * x13 * (1 + 3.06 * x) ** (1 / 6) * (1 + 0.884 * x23 + 0.471 * x43)
                / (1 + 1.64 * x23 + 0.974 * x43) * np.exp(-x))


@lru_cache(maxsize=64)
def _jp_grid(alpha_inj: float):
    """JP spectral shape S(x) for x = nu / nu_break on a log grid.

    Electrons: N(y) ~ y^-p (1 - y)^(p - 2) for y = E / E_T < 1 (Jaffe & Perola
    1973), p = 1 - 2 alpha_inj. An electron of energy E at pitch angle theta
    radiates with critical frequency nu_break y^2 sin(theta); pitch angles are
    isotropic.
    """
    p = 1 - 2 * alpha_inj
    x = np.geomspace(1e-4, 1e3, 600)
    y = np.geomspace(1e-5, 1, 1600)[:, None]
    yy = y[:, 0]
    theta = np.linspace(1e-3, np.pi / 2, 120)
    s = np.empty(x.size)
    ne = yy ** (-p) * np.clip(1 - yy, 0, None) ** (p - 2)
    for i, xv in enumerate(x):
        # integrate over pitch angle (weight sin^2 theta: sin from the solid
        # angle, sin from the emitted power) and energy
        arg = xv / (yy[:, None] ** 2 * np.sin(theta)[None, :])
        k = synchrotron_kernel(arg) * np.sin(theta)[None, :] ** 2
        inner = np.trapezoid(k, theta, axis=1)
        s[i] = np.trapezoid(ne * inner, yy)
    return x, s


def jp_spectrum(nu, nu_break, alpha_inj: float = config.ALPHA_INJ) -> np.ndarray:
    """JP spectrum (arbitrary normalisation) at frequencies ``nu`` for a given
    break frequency (same units)."""
    x, s = _jp_grid(round(alpha_inj, 3))
    return np.exp(np.interp(np.log(np.asarray(nu, float) / nu_break), np.log(x), np.log(np.maximum(s, 1e-300))))


@lru_cache(maxsize=64)
def _ci_grid(alpha_inj: float):
    """Continuous-injection (CI) spectral shape: electrons injected at a steady
    rate for a time t, each batch then ageing as in the JP model. A batch of
    age t' has break nu_b (t/t')^2, so the integrated spectrum is
    S(nu) = integral_0^1 S_JP(nu; nu_b / u^2) du, u = t'/t. Below nu_b the index
    is alpha_inj; well above it, alpha_inj - 0.5 (Kardashev 1962)."""
    x, s_jp = _jp_grid(alpha_inj)
    u = np.geomspace(1e-5, 1.0, 600)
    lx, ls = np.log(x), np.log(np.maximum(s_jp, 1e-300))
    s = np.empty(x.size)
    for i, xv in enumerate(x):
        xu = xv * u**2  # the younger batch's spectrum is sampled further below its break
        vals = np.exp(np.interp(np.log(xu), lx, ls, right=-700.0))
        low = xu < x[0]  # below the table: the injection power law
        vals[low] = np.exp(ls[0]) * (xu[low] / x[0]) ** alpha_inj
        # Each batch holds the same electrons, so its emission far below its
        # break is the same whatever its age: weight by nu_b'^alpha ~ u^(-2 alpha).
        s[i] = np.trapezoid(vals * u ** (-2 * alpha_inj), u)
    return x, s


def ci_spectrum(nu, nu_break, alpha_inj: float = config.ALPHA_INJ) -> np.ndarray:
    x, s = _ci_grid(round(alpha_inj, 3))
    return np.exp(np.interp(np.log(np.asarray(nu, float) / nu_break), np.log(x), np.log(np.maximum(s, 1e-300))))


def jp_two_point_index(nu_break_mhz, alpha_inj: float = config.ALPHA_INJ, model: str = "JP") -> np.ndarray:
    """Spectral index between 144 MHz and 1.4 GHz for a JP (or CI) spectrum."""
    nb = np.atleast_1d(np.asarray(nu_break_mhz, float))
    spec = jp_spectrum if model == "JP" else ci_spectrum
    s1 = spec(NU_LO, nb, alpha_inj)
    s2 = spec(NU_HI, nb, alpha_inj)
    return np.log(s2 / s1) / LOG_RATIO


def break_from_index(alpha_obs: float, alpha_inj: float = config.ALPHA_INJ, model: str = "JP") -> float:
    """Observed-frame break frequency (MHz) giving the measured 144-1400 MHz
    index; NaN if the spectrum is not steeper than injection (no ageing seen)."""
    if not np.isfinite(alpha_obs) or alpha_obs >= alpha_inj - 0.01:
        return np.nan
    if model == "CI" and alpha_obs <= alpha_inj - 0.5:
        return np.nan  # a CI spectrum cannot be steeper than alpha_inj - 0.5
    f = lambda lnb: jp_two_point_index(np.exp(lnb), alpha_inj, model)[0] - alpha_obs
    lo, hi = np.log(1.0), np.log(1e7)
    if f(lo) > 0:  # steeper than even a 1-MHz break can make it: outside the model
        return np.nan
    return float(np.exp(optimize.brentq(f, lo, hi)))


def spectral_age_myr(nu_break_mhz, b_ug, z) -> np.ndarray:
    """Radiative age (Myr) from the break frequency (observed frame, MHz)."""
    b_ic = config.B_CMB0_UG * (1 + np.asarray(z, float)) ** 2
    nb_ghz_rest = np.asarray(nu_break_mhz, float) / 1e3 * (1 + np.asarray(z, float))
    return 1590.0 * np.sqrt(b_ug) / (b_ug**2 + b_ic**2) / np.sqrt(nb_ghz_rest)


def max_age_field(z) -> float:
    """The field (microgauss) that maximises the age for a given break: B_IC/sqrt(3)."""
    return config.B_CMB0_UG * (1 + z) ** 2 / np.sqrt(3)


# --------------------------------------------------------------------------
# Minimum-energy magnetic field
# --------------------------------------------------------------------------

def equipartition_field_ug(flux_jy: float, nu_mhz: float, alpha: float, z: float,
                           theta_x_arcsec: float, theta_y_arcsec: float, depth_kpc: float) -> float:
    """Classical minimum-energy field (Miley 1980, ARA&A 18, 165, Eq. 2), in
    microgauss, for a region of angular size theta_x by theta_y and line-of-
    sight depth ``depth_kpc``, flux density ``flux_jy`` at ``nu_mhz``, spectral
    index alpha (S ~ nu^alpha) between nu_min and nu_max (config), proton to
    electron energy ratio k, filling factor eta, field perpendicular to the
    line of sight:

      B = 5.69e-5 [ (1+k)/eta (1+z)^(3-alpha) / (th_x th_y s)
                    F0 / nu0^alpha (nu2^(alpha+1/2) - nu1^(alpha+1/2)) / (alpha+1/2) ]^(2/7) G
    (angles in arcsec, s in kpc, F0 in Jy, frequencies in GHz).
    """
    eq = config.EQUIPARTITION
    nu0 = nu_mhz / 1e3
    nu1, nu2 = eq["nu_min_mhz"] / 1e3, eq["nu_max_mhz"] / 1e3
    a = alpha
    if abs(a + 0.5) < 1e-6:
        integ = np.log(nu2 / nu1)
    else:
        integ = (nu2 ** (a + 0.5) - nu1 ** (a + 0.5)) / (a + 0.5)
    inner = ((1 + eq["k"]) / eq["filling"] * (1 + z) ** (3 - a)
             / (theta_x_arcsec * theta_y_arcsec * depth_kpc) * flux_jy / nu0**a * integ)
    return float(5.69e-5 * inner ** (2 / 7) * 1e6)


# --------------------------------------------------------------------------
# Maps and regions
# --------------------------------------------------------------------------

@dataclass
class SpectralMaps:
    key: str
    lotss: np.ndarray  # LoTSS smoothed to 45", on the NVSS grid (Jy/beam)
    nvss: np.ndarray
    sigma_lotss: float
    sigma_nvss: float
    alpha: np.ndarray  # where both > 3 sigma
    alpha_err: np.ndarray
    alpha_upper: np.ndarray  # upper limit where LoTSS > 5 sigma but NVSS < 3 sigma
    wcs: object
    pixel_arcsec: float
    beam_arcsec: float

    @property
    def beam_area_pix(self) -> float:
        return np.pi * self.beam_arcsec**2 / (4 * np.log(2)) / self.pixel_arcsec**2


def spectral_maps(g: config.Giant) -> SpectralMaps:
    nv = rm.load_nvss(g.key)
    lo = rm.smooth_to(rm.load_lotss(g.key, 20), config.NVSS_BEAM_ARCSEC)
    lo_on_nv = rm.reproject_to(lo, nv.wcs, nv.data.shape)
    sl = rm.noise(rm.RadioImage(lo_on_nv, nv.wcs, nv.beam_arcsec, nv.pixel_arcsec, NU_LO, "LoTSS 45\""))
    sn = rm.noise(nv)
    with np.errstate(invalid="ignore", divide="ignore"):
        both = (lo_on_nv > 3 * sl) & (nv.data > 3 * sn)
        alpha = np.where(both, np.log(nv.data / lo_on_nv) / LOG_RATIO, np.nan)
        rel_l = np.hypot(sl / lo_on_nv, config.LOTSS_FLUX_SCALE_ERR)
        rel_n = np.hypot(sn / nv.data, config.NVSS_FLUX_SCALE_ERR)
        alpha_err = np.where(both, np.hypot(rel_l, rel_n) / LOG_RATIO, np.nan)
        upper = (lo_on_nv > 5 * sl) & ~(nv.data > 3 * sn) & np.isfinite(nv.data)
        alpha_up = np.where(upper, np.log(3 * sn / lo_on_nv) / LOG_RATIO, np.nan)
    return SpectralMaps(g.key, lo_on_nv, nv.data, sl, sn, alpha, alpha_err, alpha_up, nv.wcs,
                        nv.pixel_arcsec, nv.beam_arcsec)


def axis_coordinates(sm: SpectralMaps, g: config.Giant, end) -> tuple[np.ndarray, np.ndarray]:
    """Distance of each pixel from the host along the giant's axis (arcmin,
    positive toward ``end``) and across it."""
    yy, xx = np.mgrid[0:sm.nvss.shape[0], 0:sm.nvss.shape[1]]
    ra, dec = sm.wcs.pixel_to_world_values(xx.ravel(), yy.ravel())
    pa = np.radians(rm.position_angle(g.ra, g.dec, *end))
    d = rm.separation_arcmin(g.ra, g.dec, ra, dec)
    ang = _pa_vec(g.ra, g.dec, ra, dec)
    return (d * np.cos(ang - pa)).reshape(xx.shape), (d * np.sin(ang - pa)).reshape(xx.shape)


def _pa_vec(ra1, dec1, ra2, dec2) -> np.ndarray:
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    return np.arctan2(np.sin(ra2 - ra1), np.cos(dec1) * np.tan(dec2) - np.sin(dec1) * np.cos(ra2 - ra1))


def region_flux(sm: SpectralMaps, mask: np.ndarray) -> dict:
    """Integrated flux densities (Jy) and the two-point index of a region."""
    n_beams = mask.sum() / sm.beam_area_pix
    fl = float(np.nansum(np.where(mask, sm.lotss, 0)) / sm.beam_area_pix)
    fn = float(np.nansum(np.where(mask, sm.nvss, 0)) / sm.beam_area_pix)
    el = np.hypot(sm.sigma_lotss * np.sqrt(max(n_beams, 1)), config.LOTSS_FLUX_SCALE_ERR * fl)
    en = np.hypot(sm.sigma_nvss * np.sqrt(max(n_beams, 1)), config.NVSS_FLUX_SCALE_ERR * fn)
    out = {"flux_144_jy": fl, "flux_144_err": float(el), "flux_1400_jy": fn, "flux_1400_err": float(en),
           "n_beams": float(n_beams)}
    noise_n = sm.sigma_nvss * np.sqrt(max(n_beams, 1))
    out["alpha"], out["alpha_err"], out["alpha_upper"] = np.nan, np.nan, np.nan
    if fl > 0 and fn > 3 * noise_n:
        out["alpha"] = float(np.log(fn / fl) / LOG_RATIO)
        out["alpha_err"] = float(np.hypot(el / fl, en / fn) / LOG_RATIO)
    elif fl > 0:
        out["alpha_upper"] = float((np.log(max(fn, 0.0) + 3 * noise_n) - np.log(fl)) / LOG_RATIO)
    return out
