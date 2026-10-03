"""Turn the chirp into sound.

Gravitational waves are not sound, but GW150914's frequencies happen to lie in
the range of human hearing: the signal sweeps from ~35 Hz to ~250 Hz, the
pitch range of a bass guitar, in 0.2 seconds. Played through a speaker, the
detector's whitened output *is* the chirp -- no conversion beyond making the
numbers into air pressure.

Four files, each three seconds long (two before the merger, one after):

- ``gw150914_H1.wav``, ``gw150914_L1.wav`` -- whitened, band-passed detector
  data in real time. Listen for a short upward "whoomp" near the two-second mark
  under the hiss of the noise.
- ``gw150914_H1_shifted.wav`` -- the same Hanford data with every frequency
  moved up by 400 Hz. Small speakers barely reproduce 35-250 Hz; at 435-650 Hz
  the chirp is unmistakable. The timing is untouched.
- ``gw150914_template_shifted.wav`` -- the best-fit general-relativity waveform
  alone (no noise), whitened like the data and shifted the same way: the chirp
  as the theory says it should sound. This one is a model, not data.
"""

from __future__ import annotations

import numpy as np
from scipy.io import wavfile
from scipy.signal import resample_poly

import signal_processing as sp
from config import AUDIO_DIR, AUDIO_RATE, AUDIO_SHIFT_HZ, AUDIO_WINDOW_S

FADE_S = 0.1


def _clip(art: dict, name: str, key: str):
    """The audio window of one whitened stream, relative to the event."""
    d = art["det"][name]
    t = d["white3"]["t"] - (art["event_gps"] - art["ref"])
    keep = (t >= AUDIO_WINDOW_S[0]) & (t < AUDIO_WINDOW_S[1])
    return d["white3"][key][keep], d["sample_rate"]


def write_wav(path, x: np.ndarray, fs: float) -> str:
    """Taper, resample to a standard audio rate, normalise, and save as 16-bit WAV.

    The taper (a short fade at each end) stops the speaker clicking when the
    clip starts and stops. The detector data are sampled at 4096 Hz, which few
    audio players accept, so the clip is resampled to 44.1 kHz; 44100/4096 is
    exactly 11025/1024, so a polyphase resampler does it without approximation.
    """
    fade = sp.tukey(len(x), 2 * FADE_S * fs / len(x))
    y = x * fade
    assert AUDIO_RATE * 1024 == int(fs) * 11025, "resampling ratio assumes 4096 Hz -> 44.1 kHz"
    y = resample_poly(y, 11025, 1024)
    y = 0.9 * y / np.max(np.abs(y))
    wavfile.write(path, AUDIO_RATE, (y * 32767).astype(np.int16))
    return str(path)


def make_all(res, art) -> list[str]:
    written = []
    for name in ("H1", "L1"):
        x, fs = _clip(art, name, "x")
        written.append(write_wav(AUDIO_DIR / f"gw150914_{name}.wav", x, fs))
    x, fs = _clip(art, "H1", "x")
    written.append(write_wav(AUDIO_DIR / "gw150914_H1_shifted.wav",
                             sp.frequency_shift(x, fs, AUDIO_SHIFT_HZ), fs))
    x, fs = _clip(art, "H1", "fit")
    written.append(write_wav(AUDIO_DIR / "gw150914_template_shifted.wav",
                             sp.frequency_shift(x, fs, AUDIO_SHIFT_HZ), fs))
    return written
