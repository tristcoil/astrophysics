"""Read CHIME dynamic spectra ("waterfalls") and prepare them for display.

A dynamic spectrum is intensity as a function of frequency and time. CHIME
records 16384 channels, each 24.4 kHz wide, across 400-800 MHz, in ~0.98 ms
samples. The public files hold a short window around each burst that has
already been *dedispersed*: every channel shifted in time by the dispersion
delay for the burst's measured DM, so the burst lines up vertically.

Each file also stores, for every sample, the *true* time it was recorded. That
is what lets figure 02 put the data back on the time axis it arrived on and show
the sweep itself -- real samples, at their real times, with nothing simulated.
The public window is only ~0.1-0.3 s wide per channel, so outside the band that
follows the sweep there is simply no data in the release.
"""

from __future__ import annotations

from dataclasses import dataclass

import h5py
import numpy as np

import frb_physics as fp
from config import DISPLAY_FREQ_DOWNSAMPLE
from fetch_data import dynspec_path


@dataclass
class DynamicSpectrum:
    tns: str
    repeater_name: str
    #: DM (pc cm^-3) the file was dedispersed with -- CHIME's fitted value.
    dm: float
    #: Reference frequency of the dedispersion: its delay is defined as zero.
    ref_freq: float
    #: Dispersion constant CHIME used, s MHz^2 pc^-1 cm^3 (the 1/2.41e-4 convention).
    k_dm: float
    #: Channel centres after downsampling, MHz, ascending.
    freqs: np.ndarray
    #: Sample times in the dedispersed frame, s, zero at the burst's peak.
    t: np.ndarray
    #: Intensity in units of each channel's off-pulse noise; NaN where flagged.
    data: np.ndarray
    #: CHIME's best-fit burst model, in the same units.
    model: np.ndarray
    #: True recording time of every sample, s, zero when the burst's peak reached
    #: the top of the band.
    t_true: np.ndarray
    dt: float
    #: Channel width after downsampling, MHz.
    df: float

    @property
    def on_pulse(self) -> np.ndarray:
        """Samples where the frequency-averaged burst is above 3 sigma."""
        prof = self.profile()
        return prof > 3.0

    def profile(self) -> np.ndarray:
        """Frequency-averaged time series in units of its own off-pulse noise."""
        with np.errstate(invalid="ignore"):
            p = np.nanmean(self.data, axis=0)
        return _standardise_1d(p)

    def spectrum(self) -> np.ndarray:
        """On-pulse spectrum: each channel averaged over the burst's samples."""
        on = self.on_pulse
        if not on.any():
            on = np.abs(self.t) < 2 * self.dt
        with np.errstate(invalid="ignore"):
            return np.nanmean(self.data[:, on], axis=1)


def _robust_stats(x, axis=None):
    med = np.nanmedian(x, axis=axis, keepdims=True)
    mad = 1.4826 * np.nanmedian(np.abs(x - med), axis=axis, keepdims=True)
    return med, mad


def _standardise_1d(p: np.ndarray) -> np.ndarray:
    """Scale a time series to units of its off-pulse noise (median / MAD, with
    the brightest 20% of samples excluded so the burst does not set the scale)."""
    cut = np.nanpercentile(p, 80)
    off = p <= cut
    med, mad = _robust_stats(p[off])
    return (p - med.item()) / mad.item()


def load(tns: str, *, fdown: int = DISPLAY_FREQ_DOWNSAMPLE) -> DynamicSpectrum:
    """Load one burst, mask bad data, downsample in frequency, standardise.

    Downsampling averages ``fdown`` adjacent channels (default 64: 24.4 kHz ->
    1.56 MHz). Flagged samples and RFI-excised channels are excluded from the
    average; a downsampled channel with fewer than 30% valid inputs is blanked.
    Each channel is then scaled to its own off-pulse noise, so the colour scale
    reads directly as signal-to-noise per 1.56 MHz x 0.98 ms pixel.
    """
    with h5py.File(dynspec_path(tns), "r") as f:
        a = f.attrs
        data = f["data"][:].astype(np.float64)
        model = f["model"][:].astype(np.float64)
        valid = f["flag"][:] & f["good_freq"][:][:, None]
        freqs = f["index_map/freqs"][:]
        t_dd = f["index_map/times"][:]
        t_true = f["times"][:]
        dm = float(a["dm_incoherent"])
        ref_freq = float(a["ref_freq"])
        k_dm = float(a["dm_const"])
        dt = float(a["res_time"])
        rep = str(a.get("repeater_name", ""))

    data[~valid] = np.nan
    model[~valid] = np.nan
    t_true = np.where(valid, t_true, np.nan)

    nf, nt = data.shape
    nb = nf // fdown
    shape = (nb, fdown, nt)
    count = np.sum(valid.reshape(shape), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.nansum(data.reshape(shape), axis=1) / count
        m = np.nansum(model.reshape(shape), axis=1) / count
        tt = np.nansum(t_true.reshape(shape), axis=1) / count
    blank = count < 0.3 * fdown
    d[blank] = m[blank] = tt[blank] = np.nan
    fb = freqs.reshape(nb, fdown).mean(axis=1)

    order = np.argsort(fb)
    fb, d, m, tt = fb[order], d[order], m[order], tt[order]

    # Standardise every channel by its own off-pulse noise.
    prof = np.nanmean(d, axis=0)
    off = prof <= np.nanpercentile(prof, 80)
    med, mad = _robust_stats(d[:, off], axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = (d - med) / mad
        m = m / mad

    # Time zero at the burst's peak (dedispersed frame).
    peak = int(np.nanargmax(_standardise_1d(np.nanmean(d, axis=0))))
    t = t_dd - t_dd[peak]

    # True times: zero when the peak reached the top of the band. The peak's
    # true time in channel c is t_true[c, peak]; extrapolate to 800.2 MHz with
    # the same dispersion law the file documents.
    ok = np.isfinite(tt[:, peak])
    top = np.argmax(np.where(ok, fb, -np.inf))
    t0 = tt[top, peak] + fp.dispersion_delay(dm, 800.2, fb[top], k=k_dm)
    tt = tt - t0

    return DynamicSpectrum(
        tns=tns, repeater_name=rep, dm=dm, ref_freq=ref_freq, k_dm=k_dm,
        freqs=fb, t=t, data=d, model=m, t_true=tt, dt=dt, df=float(np.median(np.diff(fb))),
    )


def redisperse(ds: DynamicSpectrum, dm_new: float) -> np.ndarray:
    """The same data dedispersed with a different DM.

    A sample at dedispersed time tau in channel f was recorded at
    t_ref + delay(DM, f) + tau. Dedispersing with DM' instead places it at
    tau - [delay(DM', f) - delay(DM, f)]. We apply that shift channel by channel
    with linear interpolation; samples shifted in from outside the window are NaN.
    """
    shift = fp.dispersion_delay(dm_new - ds.dm, ds.freqs, ds.ref_freq, k=ds.k_dm)
    out = np.full_like(ds.data, np.nan)
    for i, s in enumerate(shift):
        row = ds.data[i]
        if not np.isfinite(row).any():
            continue
        good = np.isfinite(row)
        out[i] = np.interp(ds.t + s, ds.t[good], row[good], left=np.nan, right=np.nan)
    return out


def sweep_time(ds: DynamicSpectrum, freqs_mhz) -> np.ndarray:
    """When the burst's peak reached each frequency, s after it reached 800.2 MHz."""
    return fp.dispersion_delay(ds.dm, freqs_mhz, 800.2, k=ds.k_dm)


def subband_arrival_times(tns: str, *, fdown: int = 256, snr_min: float = 6.0,
                          max_err_ms: float = 0.25) -> dict:
    """Measure when the burst arrived in each sub-band, on the true time axis.

    The band is split into sub-bands of ``fdown`` channels (256 -> 6.25 MHz).
    In each, a Gaussian is fitted to the dedispersed time series around the
    burst; its centre, converted back to the true recording time with the
    file's own per-sample timestamps, is the arrival time in that sub-band.

    Two details matter at the sub-millisecond level:

    * Timestamps are averaged over the sub-band's channels. Because delay goes as
      f^-2 -- curved -- the average delay over a sub-band is the delay at
      f_eff = <f^-2>^(-1/2), not at the mean frequency. At 400 MHz the
      difference is ~0.7 ms, so each point is assigned f_eff.
    * Sub-bands where the burst is fainter than ``snr_min`` are skipped rather
      than fitted to noise, as are fits whose timing error exceeds
      ``max_err_ms`` (a quarter of a sample): those are sub-bands where
      interference or a flagged stretch leaves the pulse shape undefined.

    Nothing here assumes the f^-2 law: the arrival times are measurements, so
    fitting them with a free power law tests it.
    """
    from scipy.optimize import curve_fit

    with h5py.File(dynspec_path(tns), "r") as f:
        a = f.attrs
        data = f["data"][:].astype(np.float64)
        valid = f["flag"][:] & f["good_freq"][:][:, None]
        freqs = f["index_map/freqs"][:]
        t_true = f["times"][:]
        dt = float(a["res_time"])
        dm = float(a["dm_incoherent"])
        k_dm = float(a["dm_const"])
    data[~valid] = np.nan
    t_true = np.where(valid, t_true, np.nan)

    nf, nt = data.shape
    nb = nf // fdown
    d3 = data.reshape(nb, fdown, nt)
    t3 = t_true.reshape(nb, fdown, nt)
    f2 = freqs.reshape(nb, fdown)
    count = np.sum(np.isfinite(d3), axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        sub = np.nansum(d3, axis=1) / count
    peak = int(np.nanargmax(np.nanmean(sub, axis=0)))
    idx = np.arange(nt)

    def gauss(x, amp, x0, sig):
        return amp * np.exp(-0.5 * ((x - x0) / sig) ** 2)

    rows = []
    for i in range(nb):
        s = sub[i]
        if np.isfinite(s).sum() < 0.8 * nt:
            continue
        off = np.abs(idx - peak) > 10
        med = np.nanmedian(s[off])
        mad = 1.4826 * np.nanmedian(np.abs(s[off] - med))
        s = (s - med) / mad
        win = slice(max(peak - 8, 0), peak + 9)
        x, y = idx[win].astype(float), s[win]
        if not np.isfinite(y).all() or y.max() < snr_min:
            continue
        try:
            p, cov = curve_fit(gauss, x, y, p0=[y.max(), peak, 1.0], sigma=np.ones_like(y),
                               absolute_sigma=True, maxfev=5000)
        except RuntimeError:
            continue
        if not (0.2 < p[2] < 6) or not np.all(np.isfinite(cov)):
            continue
        if np.sqrt(cov[1, 1]) * dt * 1e3 > max_err_ms:
            continue
        j = int(np.floor(p[1]))
        frac = p[1] - j
        tj = np.nanmean(t3[i, :, j])
        tj1 = np.nanmean(t3[i, :, min(j + 1, nt - 1)])
        chans = f2[i][np.isfinite(d3[i]).any(axis=1)]
        rows.append((np.mean(chans ** -2.0) ** -0.5, tj + frac * (tj1 - tj),
                     np.sqrt(cov[1, 1]) * dt, p[0]))
    f_eff, t_arr, t_err, snr = np.array(rows).T
    return {"tns": tns, "dm_chime": dm, "k_dm": k_dm, "freq_mhz": f_eff,
            "t_arrival_s": t_arr - t_arr.min(), "t_err_s": t_err, "snr": snr}


def fit_dispersion_law(arr: dict) -> dict:
    """Fit arrival times with t = t0 + A (f/600 MHz)^alpha, alpha free, and with
    alpha fixed at -2 (which returns a DM to compare with CHIME's).

    Formal errors are scaled by sqrt(chi^2/dof) when the scatter exceeds them
    -- the standard, conservative correction when the model is not the whole
    story (here: real burst structure changing slightly with frequency).
    """
    from scipy.optimize import curve_fit

    f, t, e = arr["freq_mhz"], arr["t_arrival_s"], arr["t_err_s"]
    k = arr["k_dm"]

    def power(ff, t0, amp, alpha):
        return t0 + amp * (ff / 600.0) ** alpha

    p, cov = curve_fit(power, f, t, p0=[t.min(), k * arr["dm_chime"] / 600 ** 2, -2.0],
                       sigma=e, absolute_sigma=True)
    res = t - power(f, *p)
    chi2, dof = float(np.sum((res / e) ** 2)), len(f) - 3
    scale = max(1.0, np.sqrt(chi2 / dof))

    def cold_plasma(ff, t0, dm):
        return t0 + k * dm * ff ** -2.0

    p2, cov2 = curve_fit(cold_plasma, f, t, p0=[t.min(), arr["dm_chime"]], sigma=e, absolute_sigma=True)
    res2 = t - cold_plasma(f, *p2)
    chi2_2 = float(np.sum((res2 / e) ** 2))
    scale2 = max(1.0, np.sqrt(chi2_2 / (len(f) - 2)))
    return {
        "n_subbands": int(len(f)),
        "alpha": float(p[2]), "alpha_err": float(np.sqrt(cov[2, 2]) * scale),
        "chi2_dof_free_alpha": [chi2, dof],
        "dm_alpha_fixed": float(p2[1]), "dm_alpha_fixed_err": float(np.sqrt(cov2[1, 1]) * scale2),
        "dm_chime": arr["dm_chime"],
        "residual_rms_ms": float(np.std(res2) * 1e3),
        "median_timing_error_ms": float(np.median(e) * 1e3),
        "sweep_s": float(fp.dispersion_delay(p2[1], 400.2, 800.2, k=k)),
        "residuals_s": res2,
        "t0_alpha_fixed": float(p2[0]),
    }
