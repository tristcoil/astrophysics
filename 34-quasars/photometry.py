"""Synthetic photometry: SDSS filters, colours of a redshifted spectrum, K-corrections.

A quasar's observed colours change with redshift because its emission lines
move through the five SDSS filters and, beyond z ~ 2.2, intergalactic hydrogen
eats the light blueward of Lyman alpha. Shifting our own composite spectrum
through the filters reproduces both effects, explains why the colour-selected
SDSS sample loses quasars near z = 2.7 and 3.5, and gives the K-correction
that turns an observed i magnitude into M_i(z=2).

Conventions follow Richards et al. (2006): AB magnitudes of a photon-counting
detector, and an i-band K-correction "pegged" at z = 2 (zero for a power law at
z = 2), split into a continuum part for an alpha_nu = -0.5 power law and an
emission-line part measured from the composite.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from config import BANDS, FILTER_DIR

C_ANGSTROM_PER_S = 2.99792458e18
LYA, LYB, LYG, LYD, LYMAN_LIMIT = 1215.67, 1025.72, 972.54, 949.74, 911.75


@lru_cache(maxsize=None)
def load_filter(band: str) -> tuple[np.ndarray, np.ndarray]:
    """SDSS response for point sources through 1.3 airmasses (column 2 of the DR7 files)."""
    data = np.loadtxt(FILTER_DIR / f"{band}.dat", comments="#")
    wave, resp = data[:, 0], data[:, 1]
    return wave, resp


def filters() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {b: load_filter(b) for b in BANDS}


def pivot_wavelength(band: str) -> float:
    w, r = load_filter(band)
    return float(np.sqrt(np.trapezoid(r * w, w) / np.trapezoid(r / w, w)))


def band_mean_fnu(wave, flam, band: str, shift: float = 1.0) -> float:
    """Photon-weighted mean f_nu of a spectrum seen through ``band``.

    ``shift`` > 1 blueshifts the filter onto a rest-frame spectrum: the filter's
    wavelengths are divided by ``shift`` = 1 + z. Returns NaN when the spectrum
    does not cover the whole (shifted) band.
    """
    fw, fr = load_filter(band)
    fw = fw / shift
    good = fr > 0
    lo, hi = fw[good].min(), fw[good].max()
    if lo < wave[0] or hi > wave[-1]:
        return np.nan
    fl = np.interp(fw, wave, flam)
    fnu = fl * fw**2 / C_ANGSTROM_PER_S
    return float(np.trapezoid(fnu * fr / fw, fw) / np.trapezoid(fr / fw, fw))


def ab_mag(wave, flam, band: str) -> float:
    """AB magnitude of an observed-frame f_lambda (erg/s/cm2/A) spectrum."""
    return -2.5 * np.log10(band_mean_fnu(wave, flam, band)) - 48.60


# ---------------------------------------------------------------------------
# Intergalactic absorption (Madau 1995)
# ---------------------------------------------------------------------------

_MADAU_SERIES = ((LYA, 3.6e-3), (LYB, 1.7e-3), (LYG, 1.2e-3), (LYD, 9.3e-4))


def madau95_transmission(lam_obs, z_source):
    """Mean IGM transmission exp(-tau_eff) at observed wavelengths for a source at z.

    Madau (1995, ApJ 441, 18): Lyman-series line blanketing (his eq. 15, first
    four lines) plus the approximate photoelectric optical depth of his eq. 16
    below the Lyman limit. A model, not data: it sets the forest's mean depth
    when the composite is shifted to redshifts its own spectra did not have.
    """
    lam = np.asarray(lam_obs, float)
    tau = np.zeros_like(lam)
    zs1 = 1.0 + z_source
    for lam_j, a_j in _MADAU_SERIES:
        m = lam < lam_j * zs1
        tau[m] += a_j * (lam[m] / lam_j) ** 3.46
    m = lam < LYMAN_LIMIT * zs1
    if np.any(m):
        xc = lam[m] / LYMAN_LIMIT
        xe = zs1
        tau[m] += (0.25 * xc**3 * (xe**0.46 - xc**0.46)
                   + 9.4 * xc**1.5 * (xe**0.18 - xc**0.18)
                   - 0.7 * xc**3 * (xc**-1.32 - xe**-1.32)
                   - 0.023 * (xe**1.68 - xc**1.68))
    return np.exp(-np.clip(tau, 0, 50))


def intrinsic_blue_extension(wave, flam, alpha_cont: float, norm_wave: float,
                             euv_alpha: float = -1.76):
    """Replace the forest-absorbed part of a composite (below Ly alpha) with an
    unabsorbed model: a power law of index ``euv_alpha`` (Telfer et al. 2002 find
    alpha_nu = -1.76 shortward of ~1200 A) joined to the near-UV power law at
    Ly alpha, plus the red wing of Ly alpha mirrored onto its blue side.

    Needed only to shift the composite to *other* redshifts, where the forest
    is then put back with the mean IGM transmission.
    """
    w = np.asarray(wave, float)
    f = np.array(flam, float)
    pl_norm = np.interp(norm_wave, w, f) / (norm_wave ** (-(alpha_cont + 2)))
    cont_lya = pl_norm * LYA ** (-(alpha_cont + 2))
    blue = w < LYA
    euv = cont_lya * (w[blue] / LYA) ** (-(euv_alpha + 2))
    # mirror the red wing of the line (above the continuum) onto the blue side
    red_w = 2 * LYA - w[blue]
    red_f = np.interp(red_w, w, f)
    red_cont = pl_norm * red_w ** (-(alpha_cont + 2))
    excess = np.clip(red_f - red_cont, 0, None)
    excess *= (euv / (pl_norm * w[blue] ** (-(alpha_cont + 2))))  # scale to the bluer continuum
    f[blue] = euv + excess
    return f


# ---------------------------------------------------------------------------
# Colours and K-corrections
# ---------------------------------------------------------------------------

def observed_spectrum(rest_wave, rest_flam, z: float, igm: bool = True):
    """A rest-frame spectrum moved to redshift z (shape only; normalisation arbitrary)."""
    lam_obs = rest_wave * (1 + z)
    flam = rest_flam / (1 + z)
    if igm:
        flam = flam * madau95_transmission(lam_obs, z)
    return lam_obs, flam


def colour_track(rest_wave, rest_flam, zgrid, igm: bool = True) -> dict[str, np.ndarray]:
    """u-g, g-r, r-i, i-z (and all five magnitudes) of the spectrum at each redshift."""
    mags = {b: np.full(len(zgrid), np.nan) for b in BANDS}
    for k, z in enumerate(zgrid):
        lo, fl = observed_spectrum(rest_wave, rest_flam, z, igm)
        for b in BANDS:
            v = band_mean_fnu(lo, fl, b)
            mags[b][k] = -2.5 * np.log10(v) - 48.60 if v > 0 else np.nan
    out = dict(mags)
    for a, b in zip(BANDS[:-1], BANDS[1:]):
        out[f"{a}-{b}"] = mags[a] - mags[b]
    return out


def kcorr_continuum(z, alpha_nu: float = -0.5, z0: float = 2.0):
    """K_cont for a power law f_nu ~ nu^alpha, pegged at z0: -2.5 (1+alpha) log10((1+z)/(1+z0))."""
    return -2.5 * (1.0 + alpha_nu) * np.log10((1.0 + np.asarray(z, float)) / (1.0 + z0))


def kcorr_emission(rest_wave, total_flam, continuum_flam, zgrid, band: str = "i"):
    """K_em(z) = -2.5 log10(<F_total> / <F_continuum>) through the band at each z."""
    out = np.full(len(zgrid), np.nan)
    for k, z in enumerate(zgrid):
        tot = band_mean_fnu(rest_wave, total_flam, band, shift=1 + z)
        con = band_mean_fnu(rest_wave, continuum_flam, band, shift=1 + z)
        if np.isfinite(tot) and np.isfinite(con) and tot > 0 and con > 0:
            out[k] = -2.5 * np.log10(tot / con)
    return out


def kcorr_total(rest_wave, flam, zgrid, band: str = "i", igm: bool = True):
    """Full K-correction of a spectrum from band X at redshift z to i at z = 2.

    X - DM(z) - K = the absolute magnitude through the i band blueshifted to
    z = 2. The observed f_nu through X is (1+z) <L_nu> over the blueshifted band
    (the factor is the bandwidth stretch), so K = -2.5 log10[(1+z) <L_nu>_X,z /
    (3 <L_nu>_i,2)]; for a power law it reduces to kcorr_continuum. Note that
    Richards et al. define M_i(z=2) without the emission lines, which is what
    kcorr_continuum + kcorr_emission give; this function keeps them.
    """
    out = np.full(len(zgrid), np.nan)
    ref = 3.0 * band_mean_fnu(rest_wave, flam, "i", shift=3.0)
    for k, z in enumerate(zgrid):
        lo, fl = observed_spectrum(rest_wave, flam, z, igm)
        v = band_mean_fnu(lo, fl, band)  # = (1+z) <L_nu>, with the IGM applied
        if np.isfinite(v) and v > 0:
            out[k] = -2.5 * np.log10(v / ref)
    return out
