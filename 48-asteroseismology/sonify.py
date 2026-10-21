"""Hear the stars ring: the light curves themselves, sped up into the audible.

A red giant's oscillations have periods of hours to days, so they are far
below hearing. Playing the real Kepler light curve twenty million times faster
moves numax = 30 uHz (a red-clump star, periods of ~9 hours) to 600 Hz, close
to the D above middle C. Every giant uses the *same* speed-up, so their pitches
compare honestly: the bigger the star, the lower its note.

What you hear is the measured brightness, not a synthesis. It has been
band-passed around the oscillation hump (a Gaussian window two envelope widths
wide) to remove the low rumble of granulation and the hiss of photon noise;
gaps in the data are silent. Each clip is normalised to the same peak
loudness, so loudness does *not* compare between stars. Modes live for days
to weeks, a few tenths of a second after the speed-up, which is why the notes
waver: stochastically excited oscillations, like a bell in a sandstorm.

The Sun (SOHO/VIRGO) and 16 Cyg A (Kepler short cadence) oscillate a hundred
times faster (periods ~5 minutes) and need a smaller speed-up (100,000x) to
land in the same audible range.
"""

from __future__ import annotations

import base64
import io

import numpy as np
from scipy.io import wavfile
from scipy.signal import resample

import config as C

SPEEDUP_GIANTS = 2.0e7  # 1 uHz -> 20 Hz; 4 years of data -> 6.3 s
SPEEDUP_DWARFS = 1.0e5  # 1 uHz -> 0.1 Hz; the Sun's 3150 uHz -> 315 Hz
AUDIO_RATE = 22050
CLIP_SECONDS = 5.0
FADE_S = 0.05


def _uniform(ts):
    idx = np.round((ts.time - ts.time[0]) * 86400.0 / ts.cadence_s).astype(np.int64)
    x = np.zeros(idx.max() + 1)
    w = np.zeros(idx.max() + 1)
    x[idx] = ts.flux
    w[idx] = 1.0
    return x, w


def ringing_audio(ts, numax: float, sigma_env: float, *, speedup: float,
                  seconds: float = CLIP_SECONDS) -> tuple[np.ndarray, dict]:
    """Band-passed, sped-up light curve as float audio in [-1, 1]."""
    x, w = _uniform(ts)
    n_need = int(seconds * speedup / ts.cadence_s)
    if n_need >= len(x):
        start = 0
    else:  # the stretch with the most data in it
        cs = np.concatenate([[0], np.cumsum(w)])
        cover = cs[n_need:] - cs[:-n_need]
        start = int(np.argmax(cover))
    seg = x[start:start + n_need]
    f = np.fft.rfftfreq(len(seg), ts.cadence_s) * 1e6
    window = np.exp(-0.5 * ((f - numax) / (1.0 * sigma_env)) ** 2)
    y = np.fft.irfft(np.fft.rfft(seg) * window, n=len(seg))
    native_rate = speedup / ts.cadence_s
    y = resample(y, int(round(len(y) * AUDIO_RATE / native_rate)))
    y = y / (np.max(np.abs(y)) or 1.0) * 0.89
    nf = int(FADE_S * AUDIO_RATE)
    ramp = np.linspace(0.0, 1.0, nf)
    y[:nf] *= ramp
    y[-nf:] *= ramp[::-1]
    info = {"speedup": speedup, "pitch_hz": numax * 1e-6 * speedup,
            "data_days": len(seg) * ts.cadence_s / 86400.0,
            "start_day": float(ts.time[0] + start * ts.cadence_s / 86400.0),
            "duty": float(w[start:start + n_need].mean())}
    return y, info


def write_wav(path, y) -> None:
    wavfile.write(path, AUDIO_RATE, (np.clip(y, -1, 1) * 32767).astype(np.int16))


def wav_base64(y) -> str:
    buf = io.BytesIO()
    wavfile.write(buf, AUDIO_RATE, (np.clip(y, -1, 1) * 32767).astype(np.int16))
    return base64.b64encode(buf.getvalue()).decode("ascii")


def note_name(freq_hz: float) -> str:
    """Nearest equal-tempered note (A4 = 440 Hz), e.g. 'D5 +12 cents'."""
    names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    n = 12 * np.log2(freq_hz / 440.0) + 69
    k = int(round(n))
    cents = int(round(100 * (n - k)))
    return f"{names[k % 12]}{k // 12 - 1} {cents:+d} cents"


def make_audio(kics, *, quiet: bool = False) -> dict:
    """WAV files for the given giants, the Sun and 16 Cyg A; and a 'ladder'."""
    import json

    import lightcurves as L
    import stars
    import sun as SUN

    out = {}
    ladder = []
    gap = np.zeros(int(0.6 * AUDIO_RATE))
    for kic in kics:
        js, _ = stars.load_star(int(kic))
        ts = L.kepler_long_cadence(int(kic))
        y, info = ringing_audio(ts, js["numax"], js["envelope_fwhm"] / 2.3548, speedup=SPEEDUP_GIANTS)
        path = C.AUDIO_DIR / f"kic{int(kic)}.wav"
        write_wav(path, y)
        out[str(int(kic))] = {**info, "file": path.name, "note": note_name(info["pitch_hz"])}
        ladder += [y, gap]
        if not quiet:
            print(f"    KIC {kic}: numax {js['numax']:.1f} uHz -> {info['pitch_hz']:.0f} Hz "
                  f"({out[str(int(kic))]['note']})")
    write_wav(C.AUDIO_DIR / "red_giant_ladder.wav", np.concatenate(ladder))
    # the Sun: 2008, solar minimum
    sun = json.loads((C.CACHE_DIR / "sun.json").read_text())
    t, f, _ = SUN.virgo_series()
    ts = SUN.year_segment(t, f, 2008)
    s08 = sun["years"]["2008"]
    y, info = ringing_audio(ts, s08["numax"], s08["envelope_fwhm"] / 2.3548, speedup=SPEEDUP_DWARFS,
                            seconds=CLIP_SECONDS)
    write_wav(C.AUDIO_DIR / "sun.wav", y)
    out["sun"] = {**info, "file": "sun.wav", "note": note_name(info["pitch_hz"])}
    ms = json.loads((C.CACHE_DIR / "sixteen_cyg.json").read_text())["16 Cyg A"]
    ts = L.kepler_short_cadence("16 Cyg A")
    y, info = ringing_audio(ts, ms["numax"], ms["envelope_fwhm"] / 2.3548, speedup=SPEEDUP_DWARFS)
    write_wav(C.AUDIO_DIR / "16_cyg_a.wav", y)
    out["16 Cyg A"] = {**info, "file": "16_cyg_a.wav", "note": note_name(info["pitch_hz"])}
    (C.AUDIO_DIR / "audio.json").write_text(json.dumps(out, indent=1))
    return out
