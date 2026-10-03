"""Project configuration: paths, the event, and every tunable analysis setting.

GW150914 is the first gravitational wave ever detected: two black holes
spiralling together and merging about 1.3 billion years ago (redshift 0.09).
The signal reached the two LIGO detectors on 2015-09-14 at 09:50:45 UTC. Everything in this project
is derived from the detectors' public strain data; the published values below are
used only to grade our own measurements at the end, never as inputs to them.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
STRAIN_DIR = DATA_DIR / "strain"
PRL_DIR = DATA_DIR / "prl_figure1"
FIGURE_DIR = PROJECT_DIR / "figures"
AUDIO_DIR = PROJECT_DIR / "audio"

EVENT_JSON_PATH = DATA_DIR / "event.json"
TEMPLATE_PATH = DATA_DIR / "GW150914_4_template.hdf5"
POSTERIOR_PATH = DATA_DIR / "GW150914_GWTC-1.hdf5"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, STRAIN_DIR, PRL_DIR, FIGURE_DIR, AUDIO_DIR):
    _d.mkdir(parents=True, exist_ok=True)


@dataclass(frozen=True)
class Event:
    """One gravitational-wave event and where its data lives."""

    name: str
    #: GWOSC catalogue and release version the strain files come from.
    catalog: str
    version: int
    #: Detectors whose strain we analyse, in the order they are plotted.
    detectors: tuple[str, ...]
    #: Sample rate of the files we download. The signal lives below ~500 Hz, so
    #: the 4096 Hz release (Nyquist 2048 Hz) loses nothing; the 16 kHz release
    #: would only quadruple the download.
    sample_rate: int
    #: Length of the files we download, in seconds. GWOSC offers 32 s and 4096 s.
    #: The long files cost ~130 MB each but buy two things the short ones cannot:
    #: a noise spectrum averaged over an hour instead of half a minute, and an
    #: hour of signal-free data to measure how often noise alone mimics a signal.
    duration: int
    #: UTC reference second the discovery paper measures its time axis from.
    reference_utc: str


EVENT = Event(
    name="GW150914",
    catalog="GWTC-1-confident",
    version=3,
    detectors=("H1", "L1"),
    sample_rate=4096,
    duration=4096,
    reference_utc="2015-09-14T09:50:45",
)

# --------------------------------------------------------------------------
# Signal processing
# --------------------------------------------------------------------------

#: Noise spectrum estimate: Welch's method on segments of this length. 4 s gives
#: 0.25 Hz resolution. The hour of data yields ~2000 overlapping segments.
PSD_SEGMENT_S = 4.0
#: Median rather than mean averaging, so a loud glitch -- or the signal itself --
#: in a handful of segments cannot inflate the noise estimate.
PSD_AVERAGE = "median"
#: Taper applied to each end of the hour before it is Fourier transformed, as a
#: fraction of its length (1/128 = 16 s at each end). Without it the implied
#: wrap-around discontinuity rings through every filter.
HOUR_TAPER_ALPHA = 1.0 / 128.0
#: Seconds discarded from each end of the hour: taper plus filter settling.
HOUR_EDGE_S = 40.0

#: The band the discovery paper plots: below ~35 Hz the detectors are dominated
#: by seismic and control noise, above ~350 Hz the signal has nothing left.
DISPLAY_BAND_HZ = (35.0, 350.0)
#: Butterworth order of that band-pass (applied forwards and backwards, so the
#: effective order is twice this). Order 4 reproduces the paper's Figure 1 best.
BANDPASS_ORDER = 4

#: Instrumental lines are notched only in the un-whitened, physical-units view
#: (whitening suppresses them automatically). A line is a peak standing this many
#: times above the local median spectrum, in a 16 s (0.0625 Hz) resolution
#: spectrum. The search band reaches past 350 Hz because the L1 mirror-suspension
#: "violin" resonances near 500 Hz stand up to ~1700x above the noise floor and
#: leak through the band-pass skirt unless notched.
LINE_PSD_SEGMENT_S = 16.0
LINE_THRESHOLD = 10.0
LINE_SEARCH_BAND_HZ = (30.0, 600.0)
NOTCH_WIDTH_HZ = 1.0
#: Half-width of the stretch filtered in physical units around the event: long
#: enough to hold hundreds of signal-free windows for the injection tests.
PRL_SPAN_S = 256.0
#: Taper for that stretch (5 s at each end), and the margin kept clear of it.
PRL_TAPER_ALPHA = 0.02
PRL_EDGE_S = 16.0

#: Seconds around the event excluded whenever "signal-free" data is needed:
#: noise levels, background triggers, injection sites.
EVENT_EXCLUSION_S = 8.0

#: Lowest frequency used in matched filtering. Below it the noise rises so
#: steeply that it contributes almost no SNR.
MATCHED_FILTER_FLOW_HZ = 20.0

#: Q-transform: range of quality factors searched, frequency range shown, and
#: the maximum fractional energy loss between adjacent tiles.
Q_RANGE = (4.0, 64.0)
Q_FREQ_RANGE_HZ = (20.0, 500.0)
Q_MISMATCH = 0.2
Q_SEGMENT_S = 64.0

# --------------------------------------------------------------------------
# Arrival times, chirp mass: both are checked with injections
# --------------------------------------------------------------------------

#: Largest lag searched when cross-correlating the detectors: a little beyond
#: the 10 ms light travel time, so a peak at the edge would be noticed.
MAX_LAG_S = 0.012

#: The chirp-mass estimate keeps a half-cycle (the stretch between two zero
#: crossings) only while it stands this many noise sigmas tall, walking back
#: from the amplitude peak. Where the signal is weaker than that, noise adds
#: spurious crossings and the frequency estimate falls apart. The value was
#: chosen by injecting the GR waveform into signal-free data and taking the
#: gate that recovered the chirp mass most precisely -- never by looking at the
#: event itself.
CHIRP_GATE_SIGMA = 2.5

#: Spacing of injection sites. Windows a second apart share no noise that
#: matters: the signal lasts ~0.2 s.
INJECTION_STEP_S = 1.0

# --------------------------------------------------------------------------
# Search / significance
# --------------------------------------------------------------------------

#: The discovery paper: "event pairs that occur within a 15-ms window ...
#: determined by the 10-ms intersite propagation time plus 5 ms for uncertainty
#: in arrival time of weak signals."
COINCIDENCE_WINDOW_S = 0.015
#: Single-detector triggers are local maxima of SNR above this threshold, kept
#: only if loudest within this window. The real search used 5.5 across 250,000
#: templates; with one template an hour of this clean data yields only ~5
#: triggers above 5.5, too few to trace the background's shape, so we go lower.
TRIGGER_SNR_THRESHOLD = 4.0
TRIGGER_CLUSTER_S = 1.0
#: Frequency bands used by the chi-squared signal-consistency test.
CHISQ_BINS = 16
#: Time-slide step. Larger than the coincidence window plus the trigger
#: clustering scale, so no genuine coincidence survives any non-zero slide.
SLIDE_STEP_S = 1.0

# --------------------------------------------------------------------------
# Audio
# --------------------------------------------------------------------------

#: Clip written to WAV, relative to the event time.
AUDIO_WINDOW_S = (-2.0, 1.0)
#: The raw chirp sweeps ~35-250 Hz: audible, but barely on laptop speakers.
#: Shifting every frequency up by this much keeps the sweep's shape (a constant
#: shift, not a pitch scaling) and lands it where small speakers work.
AUDIO_SHIFT_HZ = 400.0
AUDIO_RATE = 44100

# --------------------------------------------------------------------------
# Published values -- used only to grade our results
# --------------------------------------------------------------------------

#: Numbers quoted in the text of the discovery paper: Abbott et al. 2016,
#: "Observation of Gravitational Waves from a Binary Black Hole Merger",
#: Phys. Rev. Lett. 116, 061102 (arXiv:1602.03837). Each is a direct quote,
#: checked against the paper. The GWTC-1 catalogue values (masses, distance)
#: are fetched from GWOSC instead of being typed in here.
PRL = {
    "citation": "Abbott et al. 2016, PRL 116, 061102",
    # Fig. 1 caption: "GW150914 arrived first at L1 and 6.9 +0.5/-0.4 ms later
    # at H1 ... the H1 data are also shown, shifted in time by this amount and
    # inverted (to account for the detectors' relative orientations)."
    "arrival_lag_ms": 6.9,
    "arrival_lag_ms_upper": 0.5,
    "arrival_lag_ms_lower": -0.4,
    # "Over 0.2 s, the signal increases in frequency and amplitude in about 8
    # cycles from 35 to 150 Hz, where the amplitude reaches a maximum."
    "frequency_start_hz": 35.0,
    "frequency_peak_hz": 150.0,
    "cycles": 8,
    # Abstract: "... a peak gravitational-wave strain of 1.0 x 10^-21."
    "peak_strain": 1.0e-21,
    # "Estimating f and fdot from the data in Fig. 1, we obtain a chirp mass of
    # M ~ 30 Msun, implying that the total mass M = m1 + m2 is >~ 70 Msun in the
    # detector frame. This bounds the sum of the Schwarzschild radii ... to
    # 2GM/c^2 >~ 210 km ... equal Newtonian point masses orbiting at this
    # frequency [75 Hz] would be only ~350 km apart."
    "chirp_mass_estimate_msun": 30.0,
    "total_mass_bound_msun": 70.0,
    "schwarzschild_sum_km": 210.0,
    "separation_km": 350.0,
    # "The 15-ms window is determined by the 10-ms intersite propagation time
    # plus 5 ms for uncertainty in arrival time of weak signals."
    "light_travel_ms": 10.0,
    # "The GW150914 detection-statistic value of rho_c = 23.6 is larger than any
    # background event ... 1 in 203 000 years ... corresponding to 5.1 sigma."
    "detection_statistic": 23.6,
    "far_bound_years": 203000,
    "significance_sigma": 5.1,
    "templates_in_bank": 250000,
    "background_years": 608000,
}

GWOSC_ACKNOWLEDGEMENT = (
    "This research has made use of data or software obtained from the "
    "Gravitational Wave Open Science Center (gwosc.org), a service of the LIGO "
    "Scientific Collaboration, the Virgo Collaboration, and KAGRA."
)
