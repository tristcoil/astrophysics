"""The analysis pipeline: raw strain in; a chirp, a chirp mass and a significance out.

Each step is a function below, run in this order by ``analyse``:

1. **Noise model.** Estimate each detector's noise spectrum from the hour of data
   around the event, and find its narrow instrumental lines.
2. **Whitening.** Divide the data by that noise spectrum. The chirp appears.
3. **The discovery paper's view.** Band-pass and notch the raw strain exactly as
   the paper's Figure 1 did, in physical units, and check the result against the
   paper's own published curves.
4. **Matched filter.** Correlate an hour of data against a general-relativity
   waveform: the signal-to-noise ratio, the best-fitting waveform, and every
   noise trigger in the hour.
5. **Arrival times.** When did the wave reach each site, and with what sign?
6. **Time-frequency.** The Q-transform: the chirp as a rising track.
7. **Chirp mass.** Zero crossings -> frequency sweep -> the chirp mass, by the
   discovery paper's own back-of-envelope method.
8. **Compactness.** What that chirp mass and frequency say the objects must be.
9. **Significance.** How often does noise alone produce something this loud?
   Answered with time slides, as the discovery paper did.

Two kinds of numbers come out. **Measurements** are made from the strain alone.
**Published values** (from GWOSC's catalogue record and the GWTC-1 posterior
samples) are read only in ``compare_published``, after everything is measured.

Where a measurement's uncertainty cannot be written down from a formula, it is
found by **injection**: the signal as fitted is added to a stretch of real,
signal-free detector noise from elsewhere in the hour, the identical
measurement is repeated, and this is done hundreds of times. The spread of the
answers is the error bar, and their offset from the injected truth is the
method's bias. The discovery paper's companion analysis estimated its error bars
the same way.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np
from astropy.time import Time
from astropy.utils import iers
from scipy.signal import find_peaks, hilbert
from scipy.stats import norm

import gw_physics as gp
import signal_processing as sp
from config import (
    BANDPASS_ORDER,
    CHIRP_GATE_SIGMA,
    CHISQ_BINS,
    COINCIDENCE_WINDOW_S,
    DISPLAY_BAND_HZ,
    EVENT,
    EVENT_EXCLUSION_S,
    HOUR_EDGE_S,
    HOUR_TAPER_ALPHA,
    INJECTION_STEP_S,
    LINE_PSD_SEGMENT_S,
    LINE_SEARCH_BAND_HZ,
    LINE_THRESHOLD,
    MATCHED_FILTER_FLOW_HZ,
    MAX_LAG_S,
    NOTCH_WIDTH_HZ,
    PRL,
    PRL_EDGE_S,
    PRL_SPAN_S,
    PRL_TAPER_ALPHA,
    PSD_AVERAGE,
    PSD_SEGMENT_S,
    Q_FREQ_RANGE_HZ,
    Q_MISMATCH,
    Q_RANGE,
    Q_SEGMENT_S,
    RESULTS_PATH,
    SLIDE_STEP_S,
    TRIGGER_CLUSTER_S,
    TRIGGER_SNR_THRESHOLD,
)
from fetch_data import Strain, load_posterior, load_prl_figure1, load_strain, load_template

# Everything here is from 2015; astropy's bundled Earth-orientation tables cover
# it, so never reach for the network.
iers.conf.auto_download = False


def reference_gps() -> float:
    """GPS time of 2015-09-14 09:50:45 UTC, zero of the discovery paper's time axis."""
    return float(Time(EVENT.reference_utc, scale="utc").gps)


def _index(strain: Strain, gps: float) -> int:
    return int(round((gps - strain.t0) * strain.sample_rate))


# --------------------------------------------------------------------------
# 1. Noise model
# --------------------------------------------------------------------------

def noise_model(strain: Strain) -> dict:
    """Noise spectrum of an hour of data, its narrow lines, and its steadiness.

    The whitening and matched-filter noise model is a median-averaged Welch
    estimate over the full hour. A second, finer-resolution estimate is used
    only to pin down line frequencies. Steadiness is checked by comparing the
    spectra of the hour's two halves: a single noise model is only valid if
    the noise does not drift.
    """
    fs = strain.sample_rate
    x = strain.values
    f, psd = sp.estimate_psd(x, fs, PSD_SEGMENT_S, average=PSD_AVERAGE)
    fh, psdh = sp.estimate_psd(x, fs, LINE_PSD_SEGMENT_S, average=PSD_AVERAGE)
    lines, ratios = sp.find_lines(fh, np.sqrt(psdh), LINE_SEARCH_BAND_HZ, LINE_THRESHOLD)

    half = len(x) // 2
    _, pa = sp.estimate_psd(x[:half], fs, PSD_SEGMENT_S, average=PSD_AVERAGE)
    _, pb = sp.estimate_psd(x[half:], fs, PSD_SEGMENT_S, average=PSD_AVERAGE)
    band = (f > DISPLAY_BAND_HZ[0]) & (f < DISPLAY_BAND_HZ[1])
    drift = np.sqrt(pb[band] / pa[band])

    return {
        "freqs": f, "psd": psd, "freqs_hr": fh, "psd_hr": psdh,
        "lines": lines, "line_ratios": ratios,
        "halves_asd_ratio_median": float(np.median(drift)),
        "halves_asd_ratio_p5_p95": [float(v) for v in np.percentile(drift, [5, 95])],
    }


def quality_summary(strain: Strain) -> dict:
    """Which data-quality categories and hardware injections cover the hour."""
    return {
        "seconds": int(len(strain.dq_mask)),
        "dq_failures": {n: int((~strain.flag(n)).sum()) for n in strain.dq_names},
        # Injection bits are set when there is NO injection of that kind.
        "injections_active_seconds": {
            n.removeprefix("NO_"): int((~strain.flag(n, injections=True)).sum())
            for n in strain.inj_names
        },
    }


# --------------------------------------------------------------------------
# 2. Whitening
# --------------------------------------------------------------------------

def whiten_hour(strain: Strain, noise: dict, event_gps: float) -> dict:
    """Whiten the hour, band-pass it, and express it in units of its own noise.

    Returns the whitened series (for the Q-transform), the band-passed version
    divided by its noise standard deviation (so the y axis reads directly in
    sigma), and the mask of "quiet" samples -- away from the file edges and
    the event -- that the noise level was measured on.
    """
    fs = strain.sample_rate
    white = sp.whiten(strain.values, fs, noise["freqs"], noise["psd"], taper_alpha=HOUR_TAPER_ALPHA)
    band = sp.bandpass(white, fs, *DISPLAY_BAND_HZ, order=BANDPASS_ORDER)
    quiet = _quiet_mask(strain, event_gps)
    sigma = sp.robust_sigma(band[quiet])
    band /= sigma
    return {"white": white, "band": band, "sigma": sigma, "quiet": quiet}


def _quiet_mask(strain: Strain, event_gps: float) -> np.ndarray:
    t = strain.times
    edge = HOUR_EDGE_S
    return ((t > strain.t0 + edge) & (t < strain.t0 + len(t) * strain.dt - edge)
            & (np.abs(t - event_gps) > EVENT_EXCLUSION_S))


def loudness_census(band: np.ndarray, quiet: np.ndarray, fs: float, peak_index: int) -> dict:
    """How rare is the event's peak amplitude within its own hour?

    Counts separate moments (samples above threshold, merged when within 0.1 s)
    in the quiet part of the hour whose whitened amplitude beats the event's.
    """
    w = int(0.1 * fs)
    peak = float(np.max(np.abs(band[peak_index - w: peak_index + w])))
    over = np.where(quiet & (np.abs(band) > peak))[0]
    moments = int(np.sum(np.diff(over) > w) + 1) if len(over) else 0
    return {"event_peak_sigma": peak, "louder_moments_in_hour": moments,
            "quiet_seconds": float(quiet.sum() / fs)}


def whiten_like(series: np.ndarray, fs: float, noise: dict, sigma: float) -> np.ndarray:
    """Put any series (here, the fitted waveform) through the data's exact whitening.

    Same noise model, same band-pass, same sigma, so the result sits on the
    whitened data's axis in the same units and can be overlaid on it directly.
    """
    white = sp.whiten(series, fs, noise["freqs"], noise["psd"], taper_alpha=HOUR_TAPER_ALPHA)
    return sp.bandpass(white, fs, *DISPLAY_BAND_HZ, order=BANDPASS_ORDER) / sigma


# --------------------------------------------------------------------------
# 3. The discovery paper's view, in physical units
# --------------------------------------------------------------------------

def prl_filter(x: np.ndarray, fs: float, lines) -> np.ndarray:
    """35-350 Hz band-pass plus line notches: the recipe behind the paper's Figure 1."""
    x = x * sp.tukey(len(x), PRL_TAPER_ALPHA)
    x = sp.bandpass(x, fs, *DISPLAY_BAND_HZ, order=BANDPASS_ORDER)
    return sp.notch(x, fs, lines, width_hz=NOTCH_WIDTH_HZ)


def prl_view(strain: Strain, noise: dict, ref: float, series: np.ndarray | None = None):
    """Filter a stretch of +/- ``PRL_SPAN_S`` around the reference time.

    ``series`` defaults to the strain itself; the same filter is applied to the
    best-fit waveform so the two can be compared like for like.
    """
    fs = strain.sample_rate
    i0 = _index(strain, ref - PRL_SPAN_S)
    i1 = _index(strain, ref + PRL_SPAN_S)
    x = strain.values if series is None else series
    t = (strain.t0 + np.arange(i0, i1) * strain.dt) - ref
    return t, prl_filter(x[i0:i1], fs, noise["lines"])


def compare_to_paper(t: np.ndarray, ours: np.ndarray, paper: dict) -> dict:
    """Agreement between our filtered strain and the paper's published curve.

    The paper's files are sampled at 16384 Hz; ours at 4096 Hz is interpolated
    onto their time stamps. A constant time offset is allowed for (and
    reported) in case the two time axes differ by a sample.
    """
    tp, hp = paper["t"], paper["strain"]
    lags = np.arange(-8, 9) / 16384.0
    corr = [np.corrcoef(np.interp(tp + lag, t, ours), hp)[0, 1] for lag in lags]
    best = int(np.argmax(corr))
    y = np.interp(tp + lags[best], t, ours)
    return {
        "correlation": float(corr[best]),
        "correlation_no_shift": float(corr[len(lags) // 2]),
        "best_shift_ms": float(lags[best] * 1e3),
        "rms_ratio_ours_over_paper": float(np.sqrt(np.mean(y ** 2) / np.mean(hp ** 2))),
        "residual_rms_fraction": float(np.sqrt(np.mean((y - hp) ** 2) / np.mean(hp ** 2))),
        "samples_compared": int(len(tp)),
    }


# --------------------------------------------------------------------------
# 4. Matched filter
# --------------------------------------------------------------------------

@dataclass
class TemplateInfo:
    plus: np.ndarray
    peak_index: int
    m1: float
    m2: float
    meta: dict

    @property
    def chirp_mass(self) -> float:
        return float(gp.chirp_mass(self.m1, self.m2))


def load_template_info() -> TemplateInfo:
    t = load_template()
    peak = int(np.argmax(np.abs(t["plus"] + 1j * t["cross"])))
    return TemplateInfo(plus=t["plus"], peak_index=peak, m1=float(t["meta"]["m1"]),
                        m2=float(t["meta"]["m2"]), meta=t["meta"])


def matched_filter_hour(strain: Strain, noise: dict, template: TemplateInfo,
                        event_gps: float) -> dict:
    """Filter the whole hour against the template; characterise event and noise.

    Returns the event's SNR, phase, chi-squared and arrival time; the
    best-fitting waveform in strain units (the template scaled, phase-shifted and
    time-shifted to match the data -- three numbers fitted, nothing else); every
    noise trigger in the hour with its chi-squared; and the complex filter
    output itself, which the injection tests reuse.
    """
    fs = strain.sample_rate
    n = len(strain.values)
    setup = sp.prepare_filter(sp.template_at_origin(template.plus, template.peak_index, n),
                              fs, noise["freqs"], noise["psd"], flow=MATCHED_FILTER_FLOW_HZ)
    data_fd = sp.to_frequency_domain(strain.values, fs, taper_alpha=HOUR_TAPER_ALPHA)
    mf = sp.matched_filter(setup, data_fd)
    snr = mf.snr

    # The event: the loudest sample within 0.2 s of the catalogue time.
    iev = _index(strain, event_gps)
    w = int(0.2 * fs)
    k = iev - w + int(np.argmax(snr[iev - w: iev + w]))
    z_k = mf.z[k]
    edges = sp.chisq_band_edges(setup, CHISQ_BINS)
    dof = 2 * CHISQ_BINS - 2
    chisq_event = float(sp.chisq_at(setup, data_fd, edges, np.array([k]))[0])

    # Best fit: s(f) = (z_k / sigma^2) h(f) e^{-2 pi i f t_k}. See signal_processing.
    shift = np.exp(-2j * np.pi * setup.freqs * k / fs)
    fit = np.fft.irfft((z_k / mf.sigma ** 2) * setup.template_fd * shift, n=n) * fs
    z_signal = sp.matched_filter(setup, sp.to_frequency_domain(fit, fs)).z

    # Noise triggers across the hour, each tested with chi-squared.
    usable = snr.copy()
    edge = int(HOUR_EDGE_S * fs)
    usable[:edge] = 0.0
    usable[-edge:] = 0.0
    idx, _ = find_peaks(usable, height=TRIGGER_SNR_THRESHOLD, distance=int(TRIGGER_CLUSTER_S * fs))
    chisq = sp.chisq_at(setup, data_fd, edges, idx)
    rho_hat = sp.reweighted_snr(snr[idx], chisq, dof)

    t_rel = np.arange(n) / fs
    t_peak = strain.t0 + sp.peak_time(t_rel[k - 3: k + 4], snr[k - 3: k + 4])

    # A thinned copy of the hour-long SNR for plotting: the maximum in each
    # 0.25 s bin, so no spike can be lost to the thinning.
    b = int(0.25 * fs)
    nb = n // b
    snr_thin = snr[: nb * b].reshape(nb, b).max(axis=1)

    return {
        "sigma": mf.sigma,
        "snr": float(snr[k]),
        "phase_rad": float(np.angle(z_k)),
        "index": int(k),
        "gps_peak": float(t_peak),
        "chisq": chisq_event,
        "chisq_dof": dof,
        "rho_hat": float(sp.reweighted_snr(snr[k], chisq_event, dof)),
        "amplitude_scale": float(np.abs(z_k) / mf.sigma ** 2),
        "fit": fit,
        "z": mf.z.astype(np.complex64),
        "z_signal": z_signal.astype(np.complex64),
        "triggers": {"gps": strain.t0 + idx / fs, "snr": snr[idx], "chisq": chisq,
                     "rho_hat": rho_hat},
        "snr_thin": {"t": strain.t0 + (np.arange(nb) + 0.5) * b / fs, "snr": snr_thin},
        "edge_s": HOUR_EDGE_S,
    }


# --------------------------------------------------------------------------
# 5. Arrival times
# --------------------------------------------------------------------------

def _mf_peak_time(z: np.ndarray, sigma: float, fs: float) -> float:
    a = np.abs(z) / sigma
    return sp.peak_time(np.arange(len(a)) / fs, a)


def arrival_times(det: dict, fs: float) -> dict:
    """Arrival-time difference and relative sign between H1 and L1.

    Measured two ways. The matched-filter way compares the times at which each
    detector's SNR peaks. The model-independent way cross-correlates the two
    whitened data streams directly and asks which shift lines them up best --
    no waveform model at all.

    Both are then repeated on injections -- the fitted signals added to
    signal-free stretches of both detectors -- to measure their scatter.
    Filtering is linear, so an injection's filter output is just the noise's
    output at the injection site plus the signal's own output, and no data
    needs re-filtering.
    """
    h, l = det["H1"], det["L1"]
    t_h, t_l = h["mf"]["gps_peak"], l["mf"]["gps_peak"]
    lag_mf = t_h - t_l
    phase_diff = float(np.angle(np.exp(1j * (h["mf"]["phase_rad"] - l["mf"]["phase_rad"]))))

    # Model-independent: cross-correlate whitened, band-passed data around the
    # signal (both detectors' peaks lie within this window).
    centre = 0.5 * (t_h + t_l)
    w0, w1 = int((centre - 0.20 - h["strain"].t0) * fs), int((centre + 0.05 - h["strain"].t0) * fs)
    lag_xc, corr_xc, lags, curve = sp.cross_correlation_lag(
        h["white"]["band"][w0:w1], l["white"]["band"][w0:w1], fs, MAX_LAG_S)

    # --- injections ------------------------------------------------------
    k_h, k_l = h["mf"]["index"], l["mf"]["index"]
    half = int(0.05 * fs)
    win = np.arange(-half, half + 1)
    zs_h = h["mf"]["z_signal"][k_h + win]
    zs_l = l["mf"]["z_signal"][k_l + win]
    lag_true_mf = ((_mf_peak_time(zs_h, h["mf"]["sigma"], fs) + k_h / fs)
                   - (_mf_peak_time(zs_l, l["mf"]["sigma"], fs) + k_l / fs))

    sig_h, sig_l = h["white_fit"], l["white_fit"]
    lag_true_xc = sp.cross_correlation_lag(sig_h[w0:w1], sig_l[w0:w1], fs, MAX_LAG_S)[0]

    sites = _injection_sites(h["strain"], centre, span=None, margin=1.0)
    err_mf, err_xc = [], []
    for o in sites:
        zh = h["mf"]["z"][k_h + o + win] + zs_h
        zl = l["mf"]["z"][k_l + o + win] + zs_l
        lag = ((_mf_peak_time(zh, h["mf"]["sigma"], fs) + k_h / fs)
               - (_mf_peak_time(zl, l["mf"]["sigma"], fs) + k_l / fs))
        err_mf.append(lag - lag_true_mf)
        a = h["white"]["band"][w0 + o: w1 + o] + sig_h[w0:w1]
        b = l["white"]["band"][w0 + o: w1 + o] + sig_l[w0:w1]
        err_xc.append(sp.cross_correlation_lag(a, b, fs, MAX_LAG_S)[0] - lag_true_xc)
    err_mf, err_xc = np.array(err_mf), np.array(err_xc)

    return {
        "lag_mf_ms": lag_mf * 1e3,
        "lag_mf_err_ms": _half_width_68(err_mf) * 1e3,
        "phase_difference_deg": float(np.degrees(phase_diff)),
        "lag_xcorr_ms": lag_xc * 1e3,
        "lag_xcorr_err_ms": _half_width_68(err_xc) * 1e3,
        "xcorr_peak": corr_xc,
        "xcorr_curve": {"lags_ms": lags * 1e3, "corr": curve},
        "injections": len(sites),
        "injection_errors_ms": {"mf": err_mf * 1e3, "xcorr": err_xc * 1e3},
    }


def _half_width_68(x: np.ndarray) -> float:
    lo, hi = np.percentile(x, [16, 84])
    return float(0.5 * (hi - lo))


def _injection_sites(strain: Strain, centre_gps: float, *, span: float | None,
                     margin: float) -> list[int]:
    """Sample offsets, one per ``INJECTION_STEP_S``, to signal-free injection sites.

    Sites stay ``EVENT_EXCLUSION_S`` clear of the event and ``margin`` s clear of
    the usable data's edges (the whole hour, or +/- ``span`` around the event).
    """
    fs = strain.sample_rate
    t_start, t_end = strain.t0 + HOUR_EDGE_S, strain.t0 + len(strain.values) * strain.dt - HOUR_EDGE_S
    if span is not None:
        t_start, t_end = max(t_start, centre_gps - span), min(t_end, centre_gps + span)
    offsets = np.arange(t_start + margin - centre_gps, t_end - margin - centre_gps, INJECTION_STEP_S)
    offsets = offsets[np.abs(offsets) > EVENT_EXCLUSION_S]
    return [int(round(o * fs)) for o in offsets]


# --------------------------------------------------------------------------
# 6. Time-frequency
# --------------------------------------------------------------------------

def time_frequency(strain: Strain, white: np.ndarray, centre_gps: float, ref: float) -> dict:
    """The Q-transform of whitened data around the event, best Q chosen by search.

    Computed on 64 s of whitened data so the energy normalisation has plenty of
    noise to calibrate against, then cropped to the second around the event.
    """
    fs = strain.sample_rate
    i0 = _index(strain, centre_gps - Q_SEGMENT_S / 2)
    seg = white[i0: i0 + int(Q_SEGMENT_S * fs)]
    seg_t0 = strain.t0 + i0 * strain.dt
    edge = int(4 * fs)
    mask = np.zeros(len(seg), bool)
    mask[edge:-edge] = True
    c = centre_gps - ref
    plane, scan = sp.best_q_plane(seg, fs, q_range=Q_RANGE, f_range=Q_FREQ_RANGE_HZ,
                                  mismatch=Q_MISMATCH, search_window=(c - 0.5, c + 0.5),
                                  t_ref=ref - seg_t0, norm_mask=mask)
    keep = (plane.times > c - 1.0) & (plane.times < c + 0.5)
    energy = plane.energy[:, keep]
    i, j = np.unravel_index(np.argmax(energy), energy.shape)
    return {
        "q": plane.q, "scan": scan, "frequencies": plane.frequencies,
        "times": plane.times[keep], "energy": energy,
        "peak_energy": float(energy[i, j]), "peak_frequency": float(plane.frequencies[i]),
        "peak_time": float(plane.times[keep][j]),
    }


# --------------------------------------------------------------------------
# 7. Chirp mass from the frequency sweep
# --------------------------------------------------------------------------

def fractional_shift(y: np.ndarray, dt: float, fs: float) -> np.ndarray:
    """Return y(t + dt) for any dt, by a Fourier phase ramp (no resampling loss)."""
    f = np.fft.rfftfreq(len(y), 1.0 / fs)
    return np.fft.irfft(np.fft.rfft(y) * np.exp(2j * np.pi * f * dt), n=len(y))


def combine_detectors(l1: np.ndarray, h1: np.ndarray, lag: float, fs: float) -> np.ndarray:
    """L1 plus H1 advanced by the measured lag and inverted: one coherent stream.

    The signal adds in phase, the two detectors' independent noise does not,
    so the combination is cleaner than either -- the discovery paper's
    companion analysis fitted this same sum.
    """
    return l1 - fractional_shift(h1, lag, fs)


def sweep_chirp_mass(t: np.ndarray, y: np.ndarray, sigma: float, centre: float) -> dict:
    """The zero-crossing chirp-mass measurement on one stream.

    1. Find the amplitude peak (the Hilbert envelope's maximum) within 0.1 s of
       ``centre``.
    2. Find zero crossings; each gap between two is half a wave cycle, giving a
       frequency f = 1 / (2 * gap).
    3. Walk back from the half-cycle holding the peak, keeping half-cycles while
       their peak amplitude exceeds ``CHIRP_GATE_SIGMA`` times the noise.
    4. Fit a straight line to f^(-8/3) against time. The Newtonian chirp law
       says it *is* a straight line, with slope set by the chirp mass alone.

    Returns the fit, the points, and the frequencies at the amplitude peak.
    """
    near = np.abs(t - centre) < 0.5
    tt, yy = t[near], y[near]
    env = np.abs(hilbert(yy))
    s = np.abs(tt - centre) < 0.1
    t_peak = float(tt[np.where(s)[0][np.argmax(env[s])]])

    zc = sp.zero_crossings(tt, yy)
    tmid, freq = sp.half_cycle_frequencies(zc)
    starts = np.searchsorted(tt, zc)
    amp = np.array([np.abs(yy[a:b]).max() if b > a else 0.0
                    for a, b in zip(starts[:-1], starts[1:])])

    j_peak = min(int(np.searchsorted(zc, t_peak)) - 1, len(freq) - 1)
    j = j_peak
    kept = []
    while j >= 0 and amp[j] >= CHIRP_GATE_SIGMA * sigma:
        kept.append(j)
        j -= 1
    kept = np.array(sorted(kept), dtype=int)

    out = {"t_peak": t_peak, "n_half_cycles": int(len(kept)),
           "f_at_peak": float(freq[j_peak]) if j_peak >= 0 else float("nan"),
           "f_before_peak": float(freq[j_peak - 1]) if j_peak >= 1 else float("nan"),
           "crossings": zc, "t_mid": tmid, "freq": freq, "kept": kept, "amp": amp}
    if len(kept) < 4:
        out.update(chirp_mass=float("nan"), slope=float("nan"), intercept=float("nan"),
                   r_squared=float("nan"))
        return out
    x, v = tmid[kept], freq[kept] ** (-8.0 / 3.0)
    slope, intercept = np.polyfit(x, v, 1)
    resid = v - (slope * x + intercept)
    out.update(
        slope=float(slope), intercept=float(intercept),
        chirp_mass=float(gp.chirp_mass_from_slope(slope)) if slope < 0 else float("nan"),
        r_squared=float(1.0 - np.sum(resid ** 2) / np.sum((v - v.mean()) ** 2)),
        t_coalescence=float(-intercept / slope) if slope < 0 else float("nan"),
        f_min=float(freq[kept].min()),
    )
    return out


def chirp_mass_measurement(det: dict, lag_s: float, fs: float, template: TemplateInfo,
                           event_rel: float) -> dict:
    """Measure the chirp mass from the data, and calibrate the method by injection.

    The data are measured once. The same measurement on the noise-free fitted
    signal shows what the method returns with no noise at all. Then the fitted
    signal is added to hundreds of signal-free stretches and measured again: the
    spread is the statistical error, and the typical answer versus the injected
    waveform's true chirp mass is the method's systematic bias.

    That bias is not small. The Newtonian formula ignores the extra speed-up
    of the last orbits before merger, so it *over*-estimates the chirp mass --
    which is why the discovery team's own version of this estimate gave ~37
    Msun and they quoted 30 as a conservative figure.
    """
    t = det["H1"]["prl"]["t"]
    comb = combine_detectors(det["L1"]["prl"]["data"], det["H1"]["prl"]["data"], lag_s, fs)
    comb_fit = combine_detectors(det["L1"]["prl"]["fit"], det["H1"]["prl"]["fit"], lag_s, fs)
    quiet = (np.abs(t) < PRL_SPAN_S - PRL_EDGE_S) & (np.abs(t - event_rel) > EVENT_EXCLUSION_S)
    sigma = sp.robust_sigma(comb[quiet])

    data = sweep_chirp_mass(t, comb, sigma, event_rel)
    fit_centre = float(t[np.argmax(np.abs(hilbert(comb_fit)) * (np.abs(t - event_rel) < 0.2))])
    noise_free = sweep_chirp_mass(t, comb_fit, sigma, fit_centre)

    span = PRL_SPAN_S - PRL_EDGE_S - 1.0
    offsets = np.arange(-span, span, INJECTION_STEP_S)
    offsets = offsets[np.abs(offsets) > EVENT_EXCLUSION_S]
    recovered, f_peak_rec = [], []
    for o in offsets:
        y = np.roll(comb, -int(round(o * fs))) + comb_fit
        r = sweep_chirp_mass(t, y, sigma, fit_centre)
        recovered.append(r["chirp_mass"])
        f_peak_rec.append(r["f_before_peak"])
    recovered, f_peak_rec = np.array(recovered), np.array(f_peak_rec)
    ok = np.isfinite(recovered)
    q05, q16, q50, q84, q95 = np.percentile(recovered[ok], [5, 16, 50, 84, 95])

    # Calibrate: scale our raw value by (true / typical recovered) for a signal
    # like this one. The recovered spread, scaled the same way, is the interval.
    truth = template.chirp_mass
    raw = data["chirp_mass"]
    cal = {k: truth * raw / v for k, v in (("median", q50), ("lo68", q84), ("hi68", q16),
                                           ("lo90", q95), ("hi90", q05))}
    return {
        "time": t, "combined": comb, "combined_fit": comb_fit, "noise_sigma": sigma,
        "data": data, "noise_free": noise_free,
        "injections": {"offsets_s": offsets, "recovered": recovered, "f_peak": f_peak_rec,
                       "failed": int((~ok).sum()), "median": q50,
                       "p16_p84": [q16, q84], "p5_p95": [q05, q95]},
        "template_chirp_mass": truth,
        "raw": raw,
        "calibrated": cal["median"],
        "calibrated_68": [cal["lo68"], cal["hi68"]],
        "calibrated_90": [cal["lo90"], cal["hi90"]],
        "bias_factor": q50 / truth,
    }


# --------------------------------------------------------------------------
# 8. Compactness: what could these objects be?
# --------------------------------------------------------------------------

def compactness(mc_low: float, f_gw: float) -> dict:
    """The discovery paper's argument that only black holes fit, with our numbers.

    ``mc_low`` is a conservative (low) chirp mass and ``f_gw`` a conservative
    (low) frequency at peak amplitude. The total mass is at least 2^(6/5) Mc.
    At that mass, Kepler's law gives the separation at the observed orbital
    frequency; any larger mass only makes the objects *more* compact relative
    to their size, since separation grows as M^(1/3) but Schwarzschild radius
    as M. Then each candidate kind of object is tested.
    """
    m_min = float(gp.min_total_mass(mc_low))
    sep = float(gp.keplerian_separation(f_gw, m_min))
    rs = float(gp.schwarzschild_radius(m_min))

    # Neutron star + black hole: pair the heaviest neutron star with whatever
    # reproduces the chirp mass, and see where that binary would merge.
    m_bh = gp.partner_mass(mc_low, gp.M_NEUTRON_STAR_MAX)
    m_tot_bhns = m_bh + gp.M_NEUTRON_STAR_MAX

    families = {
        "binary black hole": {
            "max_total_mass": float("inf"),
            "contact_f_at_min_mass": float(gp.contact_frequency(m_min, rs)),
        },
        "binary neutron star": {
            "max_total_mass": 2 * gp.M_NEUTRON_STAR_MAX,
            "contact_f_at_max_mass": float(gp.contact_frequency(
                2 * gp.M_NEUTRON_STAR_MAX, 2 * gp.R_NEUTRON_STAR)),
        },
        "binary white dwarf": {
            "max_total_mass": 2 * gp.M_CHANDRASEKHAR,
        },
        "main-sequence stars": {
            "contact_f_at_min_mass": float(gp.contact_frequency(
                m_min, 2 * gp.radius_main_sequence(m_min / 2))),
        },
        "black hole + neutron star": {
            "black_hole_mass": float(m_bh),
            "total_mass": float(m_tot_bhns),
            "isco_frequency": float(gp.isco_frequency(m_tot_bhns)),
        },
    }
    return {
        "chirp_mass_used": mc_low,
        "f_gw_used": f_gw,
        "orbital_frequency": f_gw / 2.0,
        "min_total_mass": m_min,
        "separation_km": sep / 1e3,
        "schwarzschild_sum_km": rs / 1e3,
        "separation_over_schwarzschild": sep / rs,
        "orbital_speed_c": float(gp.orbital_speed(f_gw, m_min)),
        "families": families,
    }


# --------------------------------------------------------------------------
# 9. Significance by time slides
# --------------------------------------------------------------------------

def _coincidences(t_a, r_a, t_b, r_b, window):
    """Pairs of triggers within ``window`` seconds; returns combined rho-hats and times."""
    order = np.argsort(t_b)
    tb, rb = t_b[order], r_b[order]
    j = np.searchsorted(tb, t_a)
    out, times = [], []
    for i, ta in enumerate(t_a):
        best = None
        for jj in (j[i] - 1, j[i]):
            if 0 <= jj < len(tb) and abs(tb[jj] - ta) < window:
                cand = np.hypot(r_a[i], rb[jj])
                best = cand if best is None or cand > best else best
        if best is not None:
            out.append(best)
            times.append(ta)
    return np.array(out), np.array(times)


def significance(trig_h: dict, trig_l: dict, t_start: float, t_end: float,
                 event_gps: float) -> dict:
    """How often does noise alone produce a coincidence as loud as the event?

    Real signals appear in both detectors within 15 ms. Shift one detector's
    triggers by a second or more and every coincidence left is between noise
    in one detector and noise in the other -- an honest sample of what noise
    alone produces. Each shift of the hour's triggers is another hour of
    background, so ~4000 shifts give ~half a year. The discovery paper did
    exactly this, 10^7 times over 16 days of data.

    The background is built twice: with every trigger (the paper's
    conservative choice for GW150914's own significance, in which the event's
    loud H1 trigger pairs with L1 noise) and with the event's triggers removed.
    """
    span = t_end - t_start
    th, rh = trig_h["gps"], trig_h["rho_hat"]
    tl, rl = trig_l["gps"], trig_l["rho_hat"]
    zero, zero_t = _coincidences(th, rh, tl, rl, COINCIDENCE_WINDOW_S)

    n_slides = int(span / SLIDE_STEP_S) - 1
    away_h = np.abs(th - event_gps) > EVENT_EXCLUSION_S
    away_l = np.abs(tl - event_gps) > EVENT_EXCLUSION_S
    backgrounds = {}
    for label, (mh, ml) in {"all_triggers": (slice(None), slice(None)),
                            "event_removed": (away_h, away_l)}.items():
        stats = []
        for k in range(1, n_slides + 1):
            shifted = t_start + (tl[ml] - t_start + k * SLIDE_STEP_S) % span
            stats.append(_coincidences(th[mh], rh[mh], shifted, rl[ml], COINCIDENCE_WINDOW_S)[0])
        backgrounds[label] = np.sort(np.concatenate(stats))[::-1] if stats else np.array([])

    livetime = n_slides * span
    i_event = int(np.argmin(np.abs(zero_t - event_gps))) if len(zero_t) else None
    event_stat = float(zero[i_event]) if i_event is not None else float("nan")
    louder = int(np.sum(backgrounds["all_triggers"] >= event_stat))
    # With no louder background event only a bound exists: FAR < 1 / livetime.
    far = (louder if louder else 1) / livetime
    p_false = 1.0 - np.exp(-far * span)
    return {
        "event_statistic": event_stat,
        "zero_lag": {"statistic": zero, "gps": zero_t},
        "background": backgrounds,
        "slides": n_slides,
        "analysed_seconds": span,
        "background_livetime_s": livetime,
        "background_livetime_days": livetime / 86400.0,
        "louder_background": louder,
        "loudest_background": float(backgrounds["all_triggers"][0]) if len(backgrounds["all_triggers"]) else float("nan"),
        "loudest_background_event_removed": float(backgrounds["event_removed"][0]) if len(backgrounds["event_removed"]) else float("nan"),
        "far_per_year_bound": far * 365.25 * 86400,
        "false_alarm_probability_bound": float(p_false),
        "sigma_bound": float(norm.isf(p_false)),
        "is_upper_bound": louder == 0,
        "triggers": {"H1": int(len(th)), "L1": int(len(tl))},
    }


# --------------------------------------------------------------------------
# Published values -- read only here, after every measurement is made
# --------------------------------------------------------------------------

def compare_published(event_record: dict, event_gps: float) -> dict:
    """GWTC-1 catalogue values and what the posterior samples imply for us.

    The posterior samples give the *detector-frame* chirp mass, the quantity
    our frequency-sweep method actually measures (cosmological redshift
    stretches every period by 1 + z, so a sweep measured on Earth weighs
    Mc * (1 + z)), and sky positions, which predict the H1-L1 arrival lag.
    """
    post = load_posterior()
    mc_det = gp.chirp_mass(post["m1_detector_frame_Msun"], post["m2_detector_frame_Msun"])
    m_det = post["m1_detector_frame_Msun"] + post["m2_detector_frame_Msun"]
    gmst = Time(event_gps, format="gps").sidereal_time("mean", "greenwich").radian
    delays = gp.arrival_delay(gp.DETECTORS["H1"], gp.DETECTORS["L1"],
                              post["right_ascension"], post["declination"], gmst)

    def summary(x):
        q5, q50, q95 = np.percentile(x, [5, 50, 95])
        return {"median": float(q50), "p5": float(q5), "p95": float(q95)}

    pe = next(p for p in event_record["parameters"].values() if p.get("pipeline_type") == "pe")
    keys = ["mass_1_source", "mass_2_source", "chirp_mass_source", "final_mass_source",
            "luminosity_distance", "redshift", "chi_eff", "a_final", "E_rad", "L_peak", "sky_area"]
    catalogue = {}
    for k in keys:
        if k in pe:
            catalogue[k] = {"value": pe[k], "lower": pe.get(f"{k}_lower"),
                            "upper": pe.get(f"{k}_upper"), "unit": pe.get(f"{k}_unit")}
    return {
        "catalogue": catalogue,
        "network_snr_by_search": {k.removeprefix("gwtc1_").removesuffix("_GW150914"): v["network_matched_filter_snr"]
                                  for k, v in event_record["parameters"].items()
                                  if v.get("pipeline_type") == "search"},
        "posterior_samples": int(len(post)),
        "chirp_mass_detector_frame": summary(mc_det),
        "total_mass_detector_frame": summary(m_det),
        "predicted_lag_ms": summary(delays * 1e3),
        "predicted_lag_samples_ms": delays * 1e3,
        "gps": event_record["GPS"],
    }


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

def analyse(event_record: dict, *, quiet: bool = False) -> tuple[dict, dict]:
    """Run every step. Returns ``(results, artefacts)``.

    ``results`` is plain JSON-able numbers; ``artefacts`` holds the arrays the
    figures and audio are drawn from, so plotting never redoes analysis.
    """
    say = (lambda *a: None) if quiet else print
    event_gps = float(event_record["GPS"])
    ref = reference_gps()
    template = load_template_info()
    paper = load_prl_figure1()

    det: dict[str, dict] = {}
    for name in EVENT.detectors:
        say(f"  {name}: loading, noise model, whitening")
        strain = load_strain(name)
        fs = strain.sample_rate
        d = {"strain": strain, "noise": noise_model(strain), "quality": quality_summary(strain)}
        d["white"] = whiten_hour(strain, d["noise"], event_gps)
        say(f"  {name}: matched filter over {len(strain.values) / fs:.0f} s")
        d["mf"] = matched_filter_hour(strain, d["noise"], template, event_gps)
        d["white_fit"] = whiten_like(d["mf"]["fit"], fs, d["noise"], d["white"]["sigma"])
        d["census"] = loudness_census(d["white"]["band"], d["white"]["quiet"], fs, d["mf"]["index"])
        t, filt = prl_view(strain, d["noise"], ref)
        _, filt_fit = prl_view(strain, d["noise"], ref, series=d["mf"]["fit"])
        d["prl"] = {"t": t, "data": filt, "fit": filt_fit}
        letter = name[0]
        d["paper_match"] = {
            "observed": compare_to_paper(t, filt, paper[f"observed-{letter}"]),
            "waveform": compare_to_paper(t, filt_fit, paper[f"waveform-{letter}"]),
        }
        det[name] = d

    fs = det["H1"]["strain"].sample_rate
    say("  arrival times (with injections)")
    arrival = arrival_times(det, fs)
    say("  Q-transforms")
    centre = 0.5 * (det["H1"]["mf"]["gps_peak"] + det["L1"]["mf"]["gps_peak"])
    for name in EVENT.detectors:
        det[name]["qplane"] = time_frequency(det[name]["strain"], det[name]["white"]["white"],
                                             centre, ref)
    say("  chirp mass (with injections)")
    chirp = chirp_mass_measurement(det, arrival["lag_mf_ms"] / 1e3, fs, template,
                                   event_gps - ref)
    mc_low = chirp["calibrated_90"][0]
    f_low = chirp["data"]["f_before_peak"]
    compact = compactness(mc_low, f_low)
    say("  time slides")
    t_start = det["H1"]["strain"].t0 + HOUR_EDGE_S
    t_end = det["H1"]["strain"].t0 + len(det["H1"]["strain"].values) * det["H1"]["strain"].dt - HOUR_EDGE_S
    sig = significance(det["H1"]["mf"]["triggers"], det["L1"]["mf"]["triggers"],
                       t_start, t_end, event_gps)
    published = compare_published(event_record, event_gps)
    h1, l1 = gp.DETECTORS["H1"], gp.DETECTORS["L1"]

    results = {
        "event": EVENT.name,
        "gps_catalogue": event_gps,
        "reference_utc": EVENT.reference_utc,
        "reference_gps": ref,
        "detectors": {
            name: {
                "quality": d["quality"],
                "noise": {
                    "asd_at_100hz": float(np.sqrt(np.interp(100.0, d["noise"]["freqs"], d["noise"]["psd"]))),
                    **_sweet_spot(d["noise"]),
                    "halves_asd_ratio_median": d["noise"]["halves_asd_ratio_median"],
                    "halves_asd_ratio_p5_p95": d["noise"]["halves_asd_ratio_p5_p95"],
                    "lines_notched_hz": [round(float(v), 3) for v in d["noise"]["lines"]],
                },
                "whitened": d["census"],
                "paper_figure1": d["paper_match"],
                "matched_filter": {
                    k: d["mf"][k] for k in ("snr", "phase_rad", "gps_peak", "chisq", "chisq_dof",
                                            "rho_hat", "amplitude_scale", "sigma")
                },
                "peak_strain_filtered": float(np.max(np.abs(d["prl"]["data"][np.abs(d["prl"]["t"] - (event_gps - ref)) < 0.1]))),
                "triggers_above_threshold": int(len(d["mf"]["triggers"]["gps"])),
                "q_transform": {k: d["qplane"][k] for k in ("q", "peak_energy", "peak_frequency", "peak_time")},
            }
            for name, d in det.items()
        },
        "network_snr": float(np.hypot(det["H1"]["mf"]["snr"], det["L1"]["mf"]["snr"])),
        "arrival": {k: v for k, v in arrival.items() if k not in ("xcorr_curve", "injection_errors_ms")},
        "geometry": {
            "light_travel_ms": gp.light_travel_time(h1, l1) * 1e3,
            "separation_km": float(np.linalg.norm(h1.vertex - l1.vertex) / 1e3),
            "tensor_overlap": gp.tensor_overlap(h1, l1),
        },
        "chirp": {
            "raw": chirp["raw"],
            "r_squared": chirp["data"]["r_squared"],
            "half_cycles": chirp["data"]["n_half_cycles"],
            "f_min": chirp["data"].get("f_min"),
            "f_at_peak": chirp["data"]["f_at_peak"],
            "f_before_peak": chirp["data"]["f_before_peak"],
            "t_peak": chirp["data"]["t_peak"],
            "noise_free_template": chirp["noise_free"]["chirp_mass"],
            "template_chirp_mass": chirp["template_chirp_mass"],
            "injections": {k: v for k, v in chirp["injections"].items()
                           if k not in ("offsets_s", "recovered", "f_peak")},
            "injection_count": int(len(chirp["injections"]["recovered"])),
            "bias_factor": chirp["bias_factor"],
            "calibrated": chirp["calibrated"],
            "calibrated_68": chirp["calibrated_68"],
            "calibrated_90": chirp["calibrated_90"],
            "gate_sigma": CHIRP_GATE_SIGMA,
        },
        "compactness": compact,
        "significance": {k: v for k, v in sig.items() if k not in ("background", "zero_lag")},
        "published": {k: v for k, v in published.items() if k != "predicted_lag_samples_ms"},
        "paper_prl": PRL,
    }
    art = {"det": det, "arrival": arrival, "chirp": chirp, "significance": sig,
           "published": published, "template": template, "ref": ref, "event_gps": event_gps}
    return results, art


def compact_artefacts(art: dict) -> dict:
    """Cut the analysis artefacts down to what the figures and audio draw.

    The analysis holds several hour-long, 16-million-sample arrays per detector;
    the figures need a few seconds of each. Keeping only those windows makes the
    plotting stage light (and cacheable) and guarantees it cannot redo analysis.
    """
    ref, event_gps = art["ref"], art["event_gps"]
    out = {"ref": ref, "event_gps": event_gps, "det": {},
           "arrival": {k: art["arrival"][k] for k in
                       ("lag_mf_ms", "lag_mf_err_ms", "lag_xcorr_ms", "lag_xcorr_err_ms",
                        "xcorr_peak", "xcorr_curve", "injection_errors_ms", "phase_difference_deg")},
           "published": art["published"], "template_chirp_mass": art["template"].chirp_mass}
    for name, d in art["det"].items():
        st = d["strain"]
        fs = st.sample_rate

        def window(x, half):
            i0, i1 = _index(st, event_gps - half), _index(st, event_gps + half)
            return x[i0:i1]

        t16 = window(st.times, 16.0) - ref
        t3 = window(st.times, 3.0) - ref
        fit32 = window(d["mf"]["fit"], 32.0)
        fit_fd = np.abs(np.fft.rfft(fit32)) / fs
        prl_keep = (d["prl"]["t"] > -1.0) & (d["prl"]["t"] < 1.5)
        trig = d["mf"]["triggers"]
        out["det"][name] = {
            "noise": {k: d["noise"][k] for k in ("freqs", "psd", "lines", "line_ratios")},
            "signal_spectrum": {"f": np.fft.rfftfreq(len(fit32), 1.0 / fs), "abs_h": fit_fd},
            "raw": {"t": t16, "h": window(st.values, 16.0)},
            "white16": {"t": t16, "x": window(d["white"]["band"], 16.0)},
            "white3": {"t": t3, "x": window(d["white"]["band"], 3.0),
                       "fit": window(d["white_fit"], 3.0)},
            "prl": {"t": d["prl"]["t"][prl_keep], "data": d["prl"]["data"][prl_keep],
                    "fit": d["prl"]["fit"][prl_keep]},
            "snr_thin": {"t": d["mf"]["snr_thin"]["t"] - ref, "snr": d["mf"]["snr_thin"]["snr"]},
            "triggers": {"t": trig["gps"] - ref, "snr": trig["snr"], "chisq": trig["chisq"],
                         "rho_hat": trig["rho_hat"]},
            "mf": {k: d["mf"][k] for k in ("snr", "chisq", "chisq_dof", "rho_hat", "gps_peak", "edge_s")},
            "census": d["census"],
            "qplane": d["qplane"],
            "sample_rate": fs,
        }
    ch = art["chirp"]
    keep = (ch["time"] > -1.0) & (ch["time"] < 1.5)
    out["chirp"] = {k: v for k, v in ch.items() if k not in ("time", "combined", "combined_fit")}
    out["chirp"].update(time=ch["time"][keep], combined=ch["combined"][keep],
                        combined_fit=ch["combined_fit"][keep])
    sig = art["significance"]
    out["significance"] = sig
    return out


def _sweet_spot(noise: dict) -> dict:
    """Lowest noise in the science band (the anti-aliasing roll-off near Nyquist excluded)."""
    band = (noise["freqs"] >= 20.0) & (noise["freqs"] <= 1000.0)
    i = int(np.argmin(noise["psd"][band]))
    return {"asd_min": float(np.sqrt(noise["psd"][band][i])),
            "asd_min_frequency": float(noise["freqs"][band][i])}


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return v if np.isfinite(v) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def save_results(results: dict) -> None:
    RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=2))
