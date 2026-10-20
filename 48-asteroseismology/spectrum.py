"""Power density spectra, calibrated so that Parseval's theorem holds.

The convention is the one used throughout Kepler asteroseismology (e.g.
Chaplin et al. 2011): a *single-sided power density* in ppm^2 per uHz whose
integral from zero to the Nyquist frequency equals the variance of the light
curve. A pure sine wave of amplitude A then shows a peak of total power A^2/2.

For unevenly sampled data (Kepler times are barycentric, and there are gaps)
the spectrum is a Lomb-Scargle periodogram (astropy, fast O(N log N) method),
which for evenly sampled data equals the classical Fourier periodogram.
astropy's ``psd`` normalisation gives P(nu) ~ |sum_k y_k exp(-2 pi i nu t_k)|^2 / N,
so the power density is PSD = 2 dt P (dt the cadence): white noise of variance
s^2 then has PSD = 2 dt s^2 everywhere, which integrates to s^2 over the
0..1/(2 dt) band (checked in test_physics.py).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Spectrum:
    nu: np.ndarray  # uHz
    psd: np.ndarray  # ppm^2 / uHz
    df: float  # grid spacing, uHz
    resolution: float  # 1 / (time span), uHz
    nyquist: float  # uHz
    oversample: int

    def window(self, lo: float, hi: float):
        m = (self.nu >= lo) & (self.nu <= hi)
        return self.nu[m], self.psd[m]


def power_density(time_days, flux_ppm, *, cadence_s: float, oversample: int = 1,
                  nu_max: float | None = None, nu_min: float | None = None) -> Spectrum:
    """Lomb-Scargle power density on a regular frequency grid.

    The grid step is the natural resolution 1/T divided by ``oversample``.
    With ``oversample=1`` neighbouring bins are (nearly) independent, as the
    chi-squared likelihood of the background fit assumes.
    """
    from astropy.timeseries import LombScargle

    t = (np.asarray(time_days) - np.min(time_days)) * 86400.0  # s
    y = np.asarray(flux_ppm, float)
    span = t.max() - t.min()
    resolution = 1e6 / span  # uHz
    nyquist = 1e6 / (2.0 * cadence_s)
    df = resolution / oversample
    lo = df if nu_min is None else max(df, nu_min)
    hi = nyquist if nu_max is None else min(nu_max, nyquist)
    nu = np.arange(lo, hi, df)
    ls = LombScargle(t, y, fit_mean=False, center_data=True, normalization="psd")
    p = ls.power(nu * 1e-6, method="fast", assume_regular_frequency=True)
    psd = 2.0 * cadence_s * p * 1e-6  # ppm^2/Hz * 1e-6 Hz/uHz -> ppm^2/uHz
    return Spectrum(nu=nu, psd=np.asarray(psd), df=df, resolution=resolution,
                    nyquist=nyquist, oversample=oversample)


def power_density_regular(flux_ppm_on_grid, valid, *, cadence_s: float,
                          oversample: int = 1) -> Spectrum:
    """Classical periodogram of data on a regular time grid (gaps = invalid).

    Gaps are filled with zeros, which for a regular grid is exactly the
    classical periodogram of the valid points. Used for the Sun, whose
    VIRGO samples are evenly spaced; an FFT of 13 million points is faster
    and exact.
    """
    x = np.where(valid, flux_ppm_on_grid - np.mean(flux_ppm_on_grid[valid]), 0.0)
    n_valid = int(valid.sum())
    n = len(x) * oversample
    X = np.fft.rfft(x, n=n)
    p = (np.abs(X) ** 2) / n_valid
    nu = np.fft.rfftfreq(n, d=cadence_s) * 1e6
    psd = 2.0 * cadence_s * p * 1e-6
    span = len(x) * cadence_s
    return Spectrum(nu=nu[1:], psd=psd[1:], df=nu[1] - nu[0], resolution=1e6 / span,
                    nyquist=1e6 / (2 * cadence_s), oversample=oversample)


def apodization(nu_uhz, nyquist_uhz: float) -> np.ndarray:
    """Power attenuation from integrating over each cadence: sinc^2.

    Each sample is an average over (nearly) the whole cadence, which damps a
    signal of frequency nu by sinc(pi nu / (2 nu_Nyq)) in amplitude (Chaplin
    et al. 2011; Kallinger et al. 2014). numpy's sinc(x) is sin(pi x)/(pi x).
    """
    return np.sinc(np.asarray(nu_uhz) / (2.0 * nyquist_uhz)) ** 2


def log_bin(nu, psd, *, nu_min: float, nu_max: float, dlnnu: float):
    """Average the spectrum into bins of constant width in ln(nu).

    Returns bin-centre frequency (mean of the frequencies averaged), mean
    power, and the number of raw bins in each average. Bins narrower than the
    grid hold one raw bin each.
    """
    m = (nu >= nu_min) & (nu <= nu_max)
    nu, psd = nu[m], psd[m]
    k = np.floor((np.log(nu) - np.log(nu_min)) / dlnnu).astype(np.int64)
    # merge raw bins that share a log bin
    uniq, start, counts = np.unique(k, return_index=True, return_counts=True)
    sums_p = np.add.reduceat(psd, start)
    sums_nu = np.add.reduceat(nu, start)
    return sums_nu / counts, sums_p / counts, counts.astype(float)


def smooth_boxcar(psd, width_bins: int) -> np.ndarray:
    """Running mean over ``width_bins`` frequency bins (edge-corrected)."""
    width_bins = max(1, int(width_bins))
    if width_bins == 1:
        return np.asarray(psd, float).copy()
    k = np.ones(width_bins)
    num = np.convolve(psd, k, mode="same")
    den = np.convolve(np.ones_like(psd), k, mode="same")
    return num / den


def smooth_gaussian(psd, sigma_bins: float) -> np.ndarray:
    from scipy.ndimage import gaussian_filter1d

    return gaussian_filter1d(np.asarray(psd, float), sigma_bins, mode="nearest")
