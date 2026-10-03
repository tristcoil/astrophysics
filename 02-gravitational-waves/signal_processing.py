"""Signal processing: getting a 1e-21 signal out of noise a hundred times larger.

LIGO's raw strain is dominated by noise that is anything but white. Below ~30 Hz
the ground shakes the mirrors (seismic noise, orders of magnitude above
everything else); narrow lines stick out of the spectrum (the 60 Hz mains and
its harmonics, calibration lines injected on purpose, the "violin" resonances of
the fibres the mirrors hang from); and above a few hundred Hz the quantum
graininess of the laser light (shot noise) takes over. GW150914 is invisible in
the raw data. Three tools, each here, make it visible:

- **Whitening** divides the data by its own noise amplitude at each frequency,
  so every frequency contributes noise of equal strength. Loud frequencies are
  turned down, quiet ones up, and a signal living where the detector is most
  sensitive stands out.
- **Band-passing and notching** keep only the band where the signal is and cut
  the instrumental lines. Used for the view in physical strain units.
- **Matched filtering** correlates the data against the waveform general
  relativity predicts. It is the optimal linear detector for a known signal in
  Gaussian noise, and it is how GW150914 was actually found.

Plus the Q-transform, a time-frequency map whose tiles narrow in time as
frequency rises, which is what makes a chirp look like a chirp.

Conventions. Fourier transforms follow numpy (``X_k = sum_n x_n e^{-2 pi i k n/N}``);
continuous-transform quantities are ``dt * X_k``. PSDs are one-sided, in
strain^2/Hz. The noise-weighted inner product is::

    <a, b> = 4 Re sum_{f>0} a(f) b*(f) / S(f) df
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal as sig
from scipy.ndimage import median_filter

# --------------------------------------------------------------------------
# Noise spectrum
# --------------------------------------------------------------------------


def estimate_psd(x: np.ndarray, fs: float, segment_s: float, *, average: str = "median"):
    """One-sided power spectral density by Welch's method.

    The data are cut into overlapping segments, each is tapered (Hann window) and
    Fourier transformed, and the squared amplitudes are averaged. Averaging is
    what turns a single, wildly fluctuating periodogram into a usable noise
    estimate. The median is used instead of the mean so that a short loud event
    -- a glitch, or the signal itself -- in a few segments cannot bias the
    estimate; scipy corrects the median's known low bias for this distribution.

    Returns ``(frequencies, psd)``. The amplitude spectral density, the quantity
    usually plotted, is ``sqrt(psd)`` in strain per root hertz.
    """
    nper = int(round(segment_s * fs))
    return sig.welch(x, fs=fs, window="hann", nperseg=nper, noverlap=nper // 2,
                     average=average, detrend=False)


def find_lines(freqs: np.ndarray, asd: np.ndarray, band, threshold: float,
               *, smooth_hz: float = 4.0, min_separation_hz: float = 0.3):
    """Narrow spectral lines: bins standing ``threshold`` times above the local median.

    The local median over ``smooth_hz`` traces the broadband noise floor and
    ignores narrow features, so the ratio of the spectrum to it picks out lines
    regardless of how steeply the floor itself slopes. Returns
    ``(line_frequencies, ratios)``.
    """
    df = freqs[1] - freqs[0]
    floor = np.exp(median_filter(np.log(asd), size=int(smooth_hz / df) | 1, mode="nearest"))
    ratio = asd / floor
    inside = (freqs >= band[0]) & (freqs <= band[1])
    peaks, _ = sig.find_peaks(np.where(inside, ratio, 0.0), height=threshold,
                              distance=max(1, int(min_separation_hz / df)))
    return freqs[peaks], ratio[peaks]


# --------------------------------------------------------------------------
# Filtering
# --------------------------------------------------------------------------


def tukey(n: int, alpha: float) -> np.ndarray:
    """A window flat in the middle with cosine tapers covering ``alpha`` of its length."""
    return sig.windows.tukey(n, alpha)


def bandpass(x: np.ndarray, fs: float, flo: float, fhi: float, *, order: int = 4) -> np.ndarray:
    """Zero-phase Butterworth band-pass.

    Run forwards and then backwards (``sosfiltfilt``) so that phase delays cancel:
    a single-pass filter would shift each frequency by a different time and
    distort the chirp's shape, as well as moving it relative to the other
    detector, which would corrupt the arrival-time measurement.
    """
    sos = sig.butter(order, [flo, fhi], btype="bandpass", fs=fs, output="sos")
    return sig.sosfiltfilt(sos, x)


def notch(x: np.ndarray, fs: float, line_freqs, *, width_hz: float = 1.0) -> np.ndarray:
    """Remove narrow spectral lines with zero-phase IIR notch filters."""
    y = x
    for f0 in line_freqs:
        b, a = sig.iirnotch(f0, Q=f0 / width_hz, fs=fs)
        y = sig.filtfilt(b, a, y)
    return y


def whiten(x: np.ndarray, fs: float, psd_freqs: np.ndarray, psd: np.ndarray,
           *, taper_alpha: float) -> np.ndarray:
    """Divide the data by its noise amplitude, frequency by frequency.

    The normalisation makes stationary Gaussian noise described by ``psd`` come
    out as white noise with unit variance per sample: for such noise
    ``E|X_k|^2 = N * S(f_k) * fs / 2``, so dividing each Fourier coefficient by
    ``sqrt(S(f_k) * fs / 2)`` leaves every frequency with the same expected
    power. The noise model is interpolated onto the data's (usually much finer)
    frequency grid.

    The taper is applied first so that the jump between the last and first
    sample -- the FFT treats the segment as periodic -- does not ring through
    the result. The tapered ends are unusable; callers trim them.
    """
    n = len(x)
    xf = np.fft.rfft(x * tukey(n, taper_alpha))
    f = np.fft.rfftfreq(n, 1.0 / fs)
    s = np.interp(f, psd_freqs, psd)
    xf /= np.sqrt(s * fs / 2.0)
    xf[0] = 0.0  # the mean carries no information and is poorly estimated
    return np.fft.irfft(xf, n=n)


def robust_sigma(x: np.ndarray) -> float:
    """Standard deviation estimated from the median absolute deviation.

    For Gaussian noise this equals the ordinary standard deviation; unlike it,
    a handful of glitches cannot drag it upward.
    """
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


# --------------------------------------------------------------------------
# Matched filtering
# --------------------------------------------------------------------------


def template_at_origin(h: np.ndarray, peak_index: int, n: int) -> np.ndarray:
    """Place a template in an ``n``-sample array with its peak at index 0.

    The inspiral before the peak wraps around to the end of the array. With this
    convention a matched-filter output peaking at sample ``k`` means the
    template's peak lines up with data sample ``k`` -- the arrival time can be
    read straight off the output. Requires ``len(h) <= n``.
    """
    out = np.zeros(n)
    out[: len(h)] = h
    return np.roll(out, -peak_index)


def to_frequency_domain(x: np.ndarray, fs: float, *, taper_alpha: float = 0.0) -> np.ndarray:
    """Continuous-transform approximation ``dt * rfft(x)``, optionally tapered first."""
    if taper_alpha:
        x = x * tukey(len(x), taper_alpha)
    return np.fft.rfft(x) / fs


@dataclass
class FilterSetup:
    """Everything about a matched filter that does not depend on the data.

    Transforming the template and interpolating the noise spectrum onto the
    frequency grid are done once here, then reused for every stretch of data and
    every chi-squared band -- on an hour of data each is a 16-million-point
    operation.
    """

    fs: float
    n: int
    freqs: np.ndarray
    inv_psd: np.ndarray      # 1/S(f) on the rfft grid, zero below flow
    template_fd: np.ndarray  # dt * rfft(template)

    @property
    def df(self) -> float:
        return self.fs / self.n

    def sigma(self, band=None) -> float:
        """Template norm sqrt(<h, h>), optionally restricted to a frequency band."""
        w = self._band_weights(band)
        return float(np.sqrt(4.0 * self.df * np.sum(np.abs(self.template_fd) ** 2 * w)))

    def _band_weights(self, band):
        if band is None:
            return self.inv_psd
        return np.where((self.freqs >= band[0]) & (self.freqs < band[1]), self.inv_psd, 0.0)


def prepare_filter(template: np.ndarray, fs: float, psd_freqs, psd, *, flow: float) -> FilterSetup:
    """Set up a matched filter for a template laid out by ``template_at_origin``."""
    n = len(template)
    f = np.fft.rfftfreq(n, 1.0 / fs)
    s = np.interp(f, psd_freqs, psd)
    inv = np.where(f >= flow, 1.0 / s, 0.0)
    return FilterSetup(fs=fs, n=n, freqs=f, inv_psd=inv,
                       template_fd=to_frequency_domain(template, fs))


@dataclass
class MatchedFilter:
    """Complex SNR time series of one template against one stretch of data.

    ``z`` is the complex filter output; ``snr = |z| / sigma``. The real part of
    ``z`` is the correlation with the template as given, the imaginary part with
    the template shifted by 90 degrees in phase, so ``|z|`` has already been
    maximised over the unknown phase of the signal. ``sigma`` is the template's
    own noise-weighted norm ``sqrt(<h, h>)``: the SNR a signal exactly equal to
    the template (at the template's amplitude) would produce.
    """

    z: np.ndarray
    sigma: float
    fs: float

    @property
    def snr(self) -> np.ndarray:
        return np.abs(self.z) / self.sigma


def matched_filter(setup: FilterSetup, data_fd: np.ndarray, band=None) -> MatchedFilter:
    """Correlate data with the template, weighted by the inverse noise spectrum.

    ``data_fd`` comes from ``to_frequency_domain`` (tapered). The output at every
    time shift t is::

        z(t) = 4 sum_{f>0} d(f) h*(f) / S(f) e^{2 pi i f t} df

    Dividing by S(f) is what makes this optimal: frequencies where the detector
    is noisy are trusted less. Summing only positive frequencies, rather than
    taking a real part, is what supplies the 90-degree-shifted copy of the
    template for free.
    """
    w = setup._band_weights(band)
    # Lay the one-sided product on the full complex FFT grid (negative
    # frequencies left at zero) and inverse-transform. numpy's ifft carries a
    # 1/N, so multiplying back by N * df gives the sum over frequencies above.
    q = np.zeros(setup.n, dtype=complex)
    q[: len(setup.freqs)] = data_fd * np.conj(setup.template_fd) * w
    z = np.fft.ifft(q)
    z *= 4.0 * setup.n * setup.df
    return MatchedFilter(z=z, sigma=setup.sigma(band), fs=setup.fs)


def chisq_band_edges(setup: FilterSetup, nbins: int) -> np.ndarray:
    """Frequency edges splitting the template into ``nbins`` bands of equal expected SNR.

    For a real signal matching the template each band then carries an equal
    share of the total SNR, while a glitch usually piles its power into one or
    two of them.
    """
    cum = np.cumsum(np.abs(setup.template_fd) ** 2 * setup.inv_psd)
    cum /= cum[-1]
    idx = np.searchsorted(cum, np.linspace(0.0, 1.0, nbins + 1)[1:-1])
    flow = setup.freqs[np.argmax(setup.inv_psd > 0)]
    return np.concatenate([[flow], setup.freqs[idx], [np.inf]])


def chisq_at(setup: FilterSetup, data_fd: np.ndarray, edges: np.ndarray,
             indices: np.ndarray) -> np.ndarray:
    """Allen's (2005) chi-squared signal-consistency test at chosen sample indices.

    The matched filter is re-run separately in each band. If the data are the
    templated signal plus Gaussian noise, each band's complex SNR carries 1/p of
    the total, and::

        chi^2 = p * sum_l |rho_l - rho/p|^2

    follows a chi-squared distribution with ``2p - 2`` degrees of freedom. A
    glitch that happens to correlate with the template does so in the wrong
    bands, and its chi^2 comes out large. Evaluated only at ``indices``: each
    band costs a full-length inverse FFT.
    """
    p = len(edges) - 1
    sigma = setup.sigma()
    parts = [matched_filter(setup, data_fd, band=(lo, hi)).z[indices] / sigma
             for lo, hi in zip(edges[:-1], edges[1:])]
    rho = np.sum(parts, axis=0)
    return p * np.sum([np.abs(r - rho / p) ** 2 for r in parts], axis=0)


def reweighted_snr(snr, chisq, dof):
    """The discovery search's ranking statistic, rho-hat.

    ``rho_hat = rho / {[1 + (chi2_r)^3] / 2}^(1/6)`` when the reduced chi-squared
    ``chi2_r = chi2/dof`` exceeds 1, else ``rho`` unchanged (PRL 116, 061102).
    A signal matching the template keeps its SNR; a glitch with a large chi^2
    is pushed down hard.
    """
    snr = np.asarray(snr, float)
    chi_r = np.asarray(chisq, float) / dof
    factor = np.where(chi_r > 1.0, ((1.0 + chi_r ** 3) / 2.0) ** (1.0 / 6.0), 1.0)
    return snr / factor


# --------------------------------------------------------------------------
# Time-frequency: the Q-transform
# --------------------------------------------------------------------------


@dataclass
class QPlane:
    """One Q-transform plane: normalised energy on a (frequency, time) grid."""

    q: float
    frequencies: np.ndarray   # row centre frequencies, Hz
    times: np.ndarray         # seconds, relative to the segment's reference time
    energy: np.ndarray        # shape (n_freq, n_time), normalised


def _q_values(q_range, mismatch):
    """Q planes spaced so neighbouring planes lose at most ``mismatch`` energy."""
    dm = 2.0 * np.sqrt(mismatch / 3.0)
    span = np.log(q_range[1] / q_range[0]) / np.sqrt(2.0)
    n = max(int(np.ceil(span / dm)), 1)
    return q_range[0] * np.exp(np.sqrt(2.0) * span / n * (np.arange(n) + 0.5))


def _q_frequencies(q, f_range, mismatch, duration):
    """Row frequencies for one Q, logarithmically spaced by the mismatch rule."""
    dm = 2.0 * np.sqrt(mismatch / 3.0)
    span = np.log(f_range[1] / f_range[0]) * np.sqrt(2.0 + q ** 2) / 2.0
    n = max(int(np.ceil(span / dm)), 1)
    step = span / n
    f = f_range[0] * np.exp(2.0 / np.sqrt(2.0 + q ** 2) * (np.arange(n) + 0.5) * step)
    return np.unique(np.round(f * duration) / duration)


def q_transform(white: np.ndarray, fs: float, *, q: float, f_range, mismatch: float,
                time_step: float = 1e-3, t_ref: float = 0.0,
                norm_mask: np.ndarray | None = None) -> QPlane:
    """Q-transform of whitened data at a single quality factor.

    Each row is a band-pass around frequency ``f0`` whose width grows with
    ``f0``: the filter is a bisquare window of half-width ``f0 * sqrt(11) / q``
    (Chatterji 2005, the construction LIGO's Omega/gwpy tools use). A fixed Q
    means a fixed number of cycles per tile, so low frequencies get long,
    frequency-sharp tiles and high frequencies short, time-sharp ones -- a match
    for a chirp, whose cycles shorten as it climbs.

    The squared magnitude of each row's complex output is the tile energy. It is
    normalised by the row's typical noise energy (median / ln 2, the mean of an
    exponential distribution estimated robustly) so that pure Gaussian noise
    gives values of about 1 and a tile of normalised energy Z has chance
    probability e^-Z. ``norm_mask`` selects which output samples define the
    noise level (the tapered edges should be excluded).
    """
    n = len(white)
    duration = n / fs
    xf = np.fft.fft(white)
    freqs = _q_frequencies(q, f_range, mismatch, duration)
    q_prime = q / np.sqrt(11.0)

    # Every row shares one output time grid, at least ``time_step`` fine and wide
    # enough to hold the broadest row's window without wrapping onto itself.
    widest = 2 * int(freqs[-1] / q_prime * duration) + 1
    n_out = int(2 ** np.ceil(np.log2(max(duration / time_step, widest))))
    energy = np.empty((len(freqs), n_out))
    for i, f0 in enumerate(freqs):
        half = int(f0 / q_prime * duration)
        j = np.arange(-half, half + 1)
        window = (1.0 - (j / (duration * f0 / q_prime)) ** 2) ** 2
        k0 = int(round(f0 * duration))
        row = np.zeros(n_out, dtype=complex)
        row[j % n_out] = xf[(k0 + j) % n] * window
        energy[i] = np.abs(np.fft.ifft(row)) ** 2
        ref = energy[i] if norm_mask is None else energy[i][_resample_mask(norm_mask, n_out)]
        energy[i] /= np.median(ref) / np.log(2.0)
    times = np.arange(n_out) * duration / n_out - t_ref
    return QPlane(q=float(q), frequencies=freqs, times=times, energy=energy)


def _resample_mask(mask: np.ndarray, n_out: int) -> np.ndarray:
    idx = (np.arange(n_out) * len(mask) / n_out).astype(int)
    return mask[idx]


def best_q_plane(white, fs, *, q_range, f_range, mismatch, search_window,
                 t_ref=0.0, norm_mask=None, time_step=1e-3) -> tuple[QPlane, list[dict]]:
    """Search Q planes and keep the one with the loudest tile inside ``search_window``.

    Different signals are best matched by different Q: a short burst wants low
    Q, a long quasi-sinusoid high Q. Scanning a range and keeping the plane with
    the highest peak energy is how the Omega pipeline chooses, and it is what
    ``gwpy``'s ``q_transform`` does by default.
    """
    best, scan = None, []
    for q in _q_values(q_range, mismatch):
        plane = q_transform(white, fs, q=q, f_range=f_range, mismatch=mismatch,
                            time_step=time_step, t_ref=t_ref, norm_mask=norm_mask)
        sel = (plane.times >= search_window[0]) & (plane.times <= search_window[1])
        peak = float(plane.energy[:, sel].max())
        scan.append({"q": float(q), "peak_energy": peak})
        if best is None or peak > best[1]:
            best = (plane, peak)
    return best[0], scan


# --------------------------------------------------------------------------
# Measuring the signal
# --------------------------------------------------------------------------


def cross_correlation_lag(a: np.ndarray, b: np.ndarray, fs: float, max_lag_s: float):
    """Time shift of ``a`` relative to ``b`` that best lines them up, and the sign.

    Computes the normalised correlation ``C(tau) = sum a(t+tau) b(t)`` over integer
    lags within +/- ``max_lag_s`` and returns the lag of the largest ``|C|``,
    refined to sub-sample precision by fitting a parabola through the peak and
    its neighbours. A positive lag means ``a`` happens later than ``b``; a
    negative correlation means ``a`` is inverted relative to ``b``.

    Returns ``(lag_s, correlation_at_peak, lags_s, correlation_curve)``.
    """
    k = int(np.ceil(max_lag_s * fs))
    n = len(b)
    lags = np.arange(-k, k + 1)
    core = slice(k, n - k)
    bb = b[core]
    c = np.array([np.dot(a[k + L: n - k + L], bb) for L in lags])
    c /= np.sqrt(np.dot(a[core], a[core]) * np.dot(bb, bb))
    i = int(np.argmax(np.abs(c)))
    shift = 0.0
    if 0 < i < len(c) - 1:
        y0, y1, y2 = np.abs(c[i - 1: i + 2])
        denom = y0 - 2 * y1 + y2
        shift = 0.5 * (y0 - y2) / denom if denom != 0 else 0.0
    return (lags[i] + shift) / fs, float(c[i]), lags / fs, c


def zero_crossings(t: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Times where ``x`` changes sign, linearly interpolated between samples."""
    s = np.signbit(x)
    i = np.where(s[:-1] != s[1:])[0]
    return t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i])


def half_cycle_frequencies(crossings: np.ndarray):
    """Instantaneous frequency from successive zero crossings.

    Consecutive crossings are half a cycle apart, so each gap gives a frequency
    ``1 / (2 * gap)`` belonging to the gap's midpoint. Returns ``(times, freqs)``.
    """
    gaps = np.diff(crossings)
    return 0.5 * (crossings[1:] + crossings[:-1]), 0.5 / gaps


def peak_time(t: np.ndarray, y: np.ndarray) -> float:
    """Time of the maximum of ``y``, refined by a parabola through the top three points."""
    i = int(np.argmax(y))
    if 0 < i < len(y) - 1:
        y0, y1, y2 = y[i - 1: i + 2]
        denom = y0 - 2 * y1 + y2
        if denom != 0:
            return float(t[i] + 0.5 * (y0 - y2) / denom * (t[1] - t[0]))
    return float(t[i])


# --------------------------------------------------------------------------
# Audio
# --------------------------------------------------------------------------


def frequency_shift(x: np.ndarray, fs: float, shift_hz: float) -> np.ndarray:
    """Move every frequency component up by ``shift_hz`` (single-sideband shift).

    The analytic signal (the data plus i times its Hilbert transform) has only
    positive frequencies, so multiplying it by ``exp(2 pi i shift t)`` slides the
    whole spectrum up without creating a mirror image below. A 35 Hz -> 250 Hz
    chirp becomes 435 Hz -> 650 Hz: the same sweep, now in a range laptop
    speakers can reproduce. Unlike speeding the recording up, this keeps the
    timing exactly as it happened.
    """
    t = np.arange(len(x)) / fs
    return np.real(sig.hilbert(x) * np.exp(2j * np.pi * shift_hz * t))
