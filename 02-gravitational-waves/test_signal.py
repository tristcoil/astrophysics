"""Self-checks for the signal-processing module, on signals with known answers.

Run with ``python test_signal.py``. No network and no LIGO data: every input
here is manufactured to have a known answer, which is the only way to test a
filter. (The real analysis never touches manufactured data -- see README.)

What these guard against is the kind of error that still produces a plausible
plot: a whitening normalisation off by sqrt(2), an SNR scale off by a factor
of two, a time lag with the wrong sign. Each would quietly corrupt a headline
number.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy import signal as sig

import signal_processing as sp

FAILURES: list[str] = []
FS = 4096.0
RNG = np.random.default_rng(20150914)


def check(label: str, got: float, expected: float, tol: float, unit: str = "",
          absolute: bool = False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if absolute else (
        abs(got - expected) / abs(expected) if expected else abs(got))
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    how = f"+/-{tol:g} {unit}" if absolute else f"tol {tol * 100:g}%"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<52} {got:>11.5g} vs {expected:>10.5g} {unit:<5} ({how})")


def truth(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def coloured_noise(n: int) -> np.ndarray:
    """Gaussian noise with a steep low-frequency wall and a narrow line, like LIGO's.

    The wall is capped below 8 Hz. An uncapped f^-6 amplitude wall spans ~17
    decades of power, and a Hann-windowed Welch estimate then leaks that power
    into every higher bin -- a limit of the test signal, not of real LIGO data,
    whose whitened in-band variance the analysis checks directly.
    """
    white = RNG.standard_normal(n)
    f = np.fft.rfftfreq(n, 1 / FS)
    shape = 1.0 + (30.0 / np.maximum(f, 8.0)) ** 4 + 50.0 * np.exp(-0.5 * ((f - 60.0) / 0.2) ** 2)
    return np.fft.irfft(np.fft.rfft(white) * shape, n=n)


def chirp_template(n_tmpl: int, f0=30.0, f1=250.0) -> np.ndarray:
    """A chirp that rises to a peak and rings down smoothly -- a stand-in shape.

    Both ends taper to zero. A waveform cut off at its loudest point would make
    its own Hilbert transform (used below to phase-shift it) ring at the cut.
    """
    t = np.arange(n_tmpl) / FS
    t_peak = 0.8 * t[-1]
    x = sig.chirp(t, f0=f0, t1=t_peak, f1=f1, method="quadratic")
    env = np.where(t < t_peak, np.exp(-0.5 * ((t - t_peak) / 0.15) ** 2),
                   np.exp(-(t - t_peak) / 0.01))
    return x * env


def main() -> int:
    print("Noise spectrum")
    # White noise of variance s^2 at rate fs has one-sided PSD 2 s^2 / fs.
    x = 3.0 * RNG.standard_normal(int(64 * FS))
    f, p = sp.estimate_psd(x, FS, 4.0)
    # (Averaged over frequency bins with a mean: each bin is chi-squared
    # distributed, so a median across bins would itself sit ~1% low.)
    check("one-sided PSD of white noise", np.mean(p[(f > 10) & (f < 1900)]), 2 * 9.0 / FS, 0.01)
    f, p = sp.estimate_psd(x, FS, 4.0, average="mean")
    check("mean and median averaging agree (bias-corrected)",
          np.mean(p[(f > 10) & (f < 1900)]), 2 * 9.0 / FS, 0.01)

    print("\nWhitening")
    n = int(256 * FS)
    y = coloured_noise(n)
    f, p = sp.estimate_psd(y, FS, 4.0)
    w = sp.whiten(y, FS, f, p, taper_alpha=0.1)
    core = w[int(32 * FS): -int(32 * FS)]
    check("whitened coloured noise has unit variance", np.var(core), 1.0, 0.02)
    fw, pw = sp.estimate_psd(core, FS, 4.0)
    # Average the whitened spectrum in 12 broad bands (each bin alone scatters
    # by ~10% from chi-squared statistics) and require the bands to agree.
    bands = [pw[(fw >= a) & (fw < a + 150)].mean() for a in np.arange(20, 1800, 150)]
    check("whitened spectrum is flat (max/min band power)", max(bands) / min(bands), 1.0, 0.05)
    lines, ratios = sp.find_lines(f, np.sqrt(p), (30, 400), 5.0)
    truth(f"the injected 60 Hz line is found ({', '.join(f'{v:.2f}' for v in lines)} Hz)",
          any(abs(lines - 60.0) < 0.3))
    check("robust sigma of Gaussian noise", sp.robust_sigma(RNG.standard_normal(200000)), 1.0, 0.01)

    print("\nMatched filter: a known signal in known noise")
    n = int(64 * FS)
    noise = coloured_noise(n)
    f, p = sp.estimate_psd(noise, FS, 4.0)
    h = chirp_template(int(1.0 * FS))
    peak = int(np.argmax(np.abs(sig.hilbert(h))))
    tmpl = sp.template_at_origin(h, peak, n)
    setup = sp.prepare_filter(tmpl, FS, f, p, flow=20.0)
    # Inject the template phase-shifted by 60 degrees, peak at a known sample.
    # At SNR 200 the noise is negligible, so normalisation, timing and phase
    # convention are tested to a fraction of a percent; at SNR 12 -- roughly
    # GW150914 in Livingston -- the answers must be right to within noise.
    where, phase = int(30.3 * FS), np.deg2rad(60)
    shifted = np.real(sig.hilbert(h) * np.exp(-1j * phase))

    def inject(target_snr):
        out = np.zeros(n)
        out[where - peak: where - peak + len(h)] = target_snr / setup.sigma() * shifted
        return out

    loud = sp.matched_filter(setup, sp.to_frequency_domain(noise + inject(200.0), FS, taper_alpha=0.1))
    k = int(np.argmax(loud.snr[int(8 * FS): -int(8 * FS)])) + int(8 * FS)
    check("SNR-200 injection: recovered SNR", loud.snr[k], 200.0, 0.01)
    check("SNR-200 injection: arrival sample", k, where, 0, "samples", absolute=True)
    # A signal lagging the template in phase by phi comes out as z ~ e^{-i phi}.
    check("SNR-200 injection: phase", np.rad2deg(np.angle(loud.z[k])), -60.0, 0.5, "deg",
          absolute=True)

    injection = inject(12.0)
    mf = sp.matched_filter(setup, sp.to_frequency_domain(noise + injection, FS, taper_alpha=0.1))
    k = int(np.argmax(mf.snr[int(8 * FS): -int(8 * FS)])) + int(8 * FS)
    # Noise adds ~N(0, 1) to the SNR and ~1/SNR radians to the phase: allow 3 sigma.
    check("SNR-12 injection: SNR at the true time", mf.snr[where], 12.0, 3.0, absolute=True)
    check("SNR-12 injection: phase at the true time", np.rad2deg(np.angle(mf.z[where])), -60.0,
          3 * np.rad2deg(1 / 12.0), "deg", absolute=True)
    check("SNR-12 injection: loudest sample", k, where, 2, "samples", absolute=True)
    mf0 = sp.matched_filter(setup, sp.to_frequency_domain(noise, FS, taper_alpha=0.1))
    core = mf0.snr[int(8 * FS): -int(8 * FS)]
    # In Gaussian noise |z|^2/sigma^2 is chi-squared with 2 degrees of freedom.
    check("noise-only <SNR^2> = 2 (two phases)", np.mean(core ** 2), 2.0, 0.05)

    print("\nChi-squared test")
    edges = sp.chisq_band_edges(setup, 16)
    chi = sp.chisq_at(setup, sp.to_frequency_domain(noise + injection, FS, taper_alpha=0.1),
                      edges, np.array([k]))[0]
    check("chi^2 per dof for a matching signal", chi / 30, 1.0, 0.6, absolute=True)
    # A sine-Gaussian "glitch" loud enough to ring the template anyway.
    t = np.arange(n) / FS
    glitch = 4e3 * np.exp(-0.5 * ((t - 40.0) / 0.004) ** 2) * np.sin(2 * np.pi * 150 * (t - 40.0))
    gfd = sp.to_frequency_domain(noise + glitch, FS, taper_alpha=0.1)
    gsnr = sp.matched_filter(setup, gfd).snr
    kg = int(40.0 * FS) - 400 + int(np.argmax(gsnr[int(40.0 * FS) - 400: int(40.0 * FS) + 400]))
    chig = sp.chisq_at(setup, gfd, edges, np.array([kg]))[0]
    truth(f"a glitch with SNR {gsnr[kg]:.1f} has chi^2/dof = {chig / 30:.0f} >> 1", chig / 30 > 10)
    rh = sp.reweighted_snr([gsnr[kg], mf.snr[k]], [chig, chi], 30)
    truth(f"reweighting demotes the glitch ({gsnr[kg]:.1f} -> {rh[0]:.1f}) "
          f"but spares the signal ({mf.snr[k]:.1f} -> {rh[1]:.1f})",
          rh[0] < 0.5 * gsnr[kg] and rh[1] > 0.9 * mf.snr[k])

    print("\nQ-transform")
    n = int(16 * FS)
    white = RNG.standard_normal(n)
    t = np.arange(n) / FS
    burst = 6.0 * np.exp(-((t - 8.0) / 0.02) ** 2) * np.cos(2 * np.pi * 120.0 * (t - 8.0))
    plane = sp.q_transform(white + burst, FS, q=8.0, f_range=(20, 500), mismatch=0.2)
    edge = (plane.times > 1) & (plane.times < 15) & (np.abs(plane.times - 8.0) > 0.5)
    check("noise-only tiles have mean normalised energy 1", plane.energy[:, edge].mean(), 1.0, 0.1)
    i, j = np.unravel_index(np.argmax(plane.energy), plane.energy.shape)
    check("sine-Gaussian found at its frequency", plane.frequencies[i], 120.0, 0.08, "Hz")
    check("sine-Gaussian found at its time", plane.times[j], 8.0, 0.01, "s", absolute=True)

    print("\nArrival-time lag and sign")
    base = sp.bandpass(RNG.standard_normal(int(2 * FS)), FS, 35, 350)
    lag_true = 3.3e-3
    # a(t) = -b(t - lag): a is b delayed by lag and inverted.
    fb = np.fft.rfft(base)
    ff = np.fft.rfftfreq(len(base), 1 / FS)
    delayed = -np.fft.irfft(fb * np.exp(-2j * np.pi * ff * lag_true), n=len(base))
    lag, corr, _, _ = sp.cross_correlation_lag(delayed, base, FS, 0.012)
    check("lag of a delayed copy (ms)", lag * 1e3, lag_true * 1e3, 0.05, "ms", absolute=True)
    truth(f"inverted copy gives correlation {corr:+.3f} (negative)", corr < -0.95)

    print("\nFrequency from zero crossings")
    t = np.arange(int(0.5 * FS)) / FS
    tm, fz = sp.half_cycle_frequencies(sp.zero_crossings(t, np.sin(2 * np.pi * 47.3 * t + 0.4)))
    check("zero-crossing frequency of a 47.3 Hz tone", np.mean(fz), 47.3, 1e-3, "Hz")
    check("... and its scatter", np.std(fz), 0.0, 0.05, "Hz", absolute=True)

    print("\nFrequency shift for audio")
    t = np.arange(int(4 * FS)) / FS
    shifted = sp.frequency_shift(np.sin(2 * np.pi * 100.0 * t), FS, 400.0)
    spec = np.abs(np.fft.rfft(shifted * np.hanning(len(t))))
    check("100 Hz tone shifted by +400 Hz", np.fft.rfftfreq(len(t), 1 / FS)[np.argmax(spec)],
          500.0, 0.25, "Hz", absolute=True)
    below = spec[np.fft.rfftfreq(len(t), 1 / FS) < 300].max() / spec.max()
    truth(f"no mirror image left behind (residual {below:.1e})", below < 1e-3)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All signal-processing checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
