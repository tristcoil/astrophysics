"""Project configuration: paths, physical constants and every adopted input.

Everything we *assume* lives here, each value with its source: the Kepler
timing constants, the solar reference values of the scaling relations, the
sample definition and the analysis settings. Published *results* that we grade
ourselves against live separately in ``literature.py`` and never enter a fit.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
AUDIO_DIR = FIGURE_DIR / "audio"
CACHE_DIR = DATA_DIR / "cache"
CATALOG_DIR = DATA_DIR / "catalogs"
KEPLER_DIR = DATA_DIR / "kepler"
VIRGO_DIR = DATA_DIR / "virgo"
GAIA_DIR = DATA_DIR / "gaia"
STAR_CACHE_DIR = CACHE_DIR / "stars"
RESULTS_PATH = DATA_DIR / "results.json"
SAMPLE_PATH = DATA_DIR / "sample.csv"

for _d in (DATA_DIR, FIGURE_DIR, AUDIO_DIR, CACHE_DIR, CATALOG_DIR, KEPLER_DIR,
           VIRGO_DIR, GAIA_DIR, STAR_CACHE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: Plain tool user agent. CDS answers browser-like agents with an "are you a
#: bot" HTML page *with HTTP 200*, so downloads are also checked for HTML.
USER_AGENT = "astrophysics-atlas-B5/1.0 (python-urllib; asteroseismology project)"

# ---------------------------------------------------------------------------
# Kepler timing (Kepler Instrument Handbook, KSCI-19033; Gilliland et al. 2010)
# ---------------------------------------------------------------------------

#: One long-cadence sample is 270 co-added frames of 6.5387 s (6.0198 s
#: exposure + 0.5189 s readout): 1765.46 s = 29.4244 min.
KEPLER_LC_CADENCE_S = 1765.4615
#: One short-cadence sample is 9 frames: 58.85 s.
KEPLER_SC_CADENCE_S = 58.84876
NYQUIST_LC_UHZ = 1e6 / (2 * KEPLER_LC_CADENCE_S)  # 283.21 uHz
NYQUIST_SC_UHZ = 1e6 / (2 * KEPLER_SC_CADENCE_S)  # 8496.4 uHz

#: Kepler times in the light-curve files are BJD - 2454833.
KEPLER_BJD_OFFSET = 2454833.0

#: VIRGO/SPM samples the Sun once a minute (1-minute averages).
VIRGO_CADENCE_S = 60.0
NYQUIST_VIRGO_UHZ = 1e6 / (2 * VIRGO_CADENCE_S)  # 8333.3 uHz

# ---------------------------------------------------------------------------
# Solar reference values for the scaling relations
# ---------------------------------------------------------------------------
# The scaling relations are anchored on the Sun, and each pipeline measures the
# Sun's numax a little differently, so every catalogue comes with its own pair.
# We compare with each catalogue on *its* reference values and also measure the
# Sun with our own code (VIRGO photometry, analysis.solar_reference).

#: Huber et al. (2011, ApJ 743, 143), SYD pipeline on VIRGO photometry; adopted
#: by Yu et al. (2018) for their 16,000 red giants (their Section 3.3).
NUMAX_SUN_SYD = 3090.0
NUMAX_SUN_SYD_ERR = 30.0
DNU_SUN_SYD = 135.1
DNU_SUN_SYD_ERR = 0.1
TEFF_SUN_YU = 5777.0

#: APOKASC-2/3 (Pinsonneault et al. 2018; 2025, ApJS 276, 69, Section 3.1):
#: numax_sun = 3076 uHz, Dnu_sun = 135.1416 uHz. APOKASC-3 then multiplies
#: numax by a calibration f_numax (their Eq. 4-5) tied to Gaia radii.
NUMAX_SUN_APOKASC = 3076.0
DNU_SUN_APOKASC = 135.1416

#: Nominal solar effective temperature (IAU 2015 Resolution B3; Prsa et al.
#: 2016). APOKASC-3's tabulated masses are reproduced with this value
#: (checked in analysis.reproduce_apokasc_masses).
TEFF_SUN_NOMINAL = 5772.0

# ---------------------------------------------------------------------------
# Sample
# ---------------------------------------------------------------------------

#: Stars are drawn from APOKASC-3's Gold sample (numax and Dnu from 5+
#: pipelines, APOGEE spectra, an evolutionary state) that Yu et al. (2018) also
#: measured, with at least this many quarters of long-cadence data ...
MIN_QUARTERS = 15
#: ... and at least this bright in the Kepler band (Yu's Kp), so that the
#: high-numax stars, whose oscillations are smallest, are still well detected.
MAX_KEPMAG = 12.5

#: Bins in (Yu et al.'s) numax, log-spaced from 5 to 230 uHz: the range
#: Yu et al. cover (5-275 uHz) minus the top, where long cadence's Nyquist
#: frequency (283 uHz) starts to reflect the oscillation hump.
NUMAX_BIN_EDGES = (5.0, 7.5, 11.0, 16.0, 23.0, 33.0, 48.0, 70.0, 100.0, 150.0, 230.0)
#: Stars drawn at random from each bin. Where red-clump stars exist (bins
#: centred 16-100 uHz), half the draw is red clump and half red giant branch.
STARS_PER_BIN = 10
SAMPLE_SEED = 20261006

#: Stars added by hand because the literature uses them as showcases.
SHOWCASE_STARS = {
    6779699: "Bedding et al. 2011 (Nature 471, 608), Fig. 2: hydrogen-shell burning, dP_obs = 53 s",
    4902641: "Bedding et al. 2011, Fig. 2: core helium burning (red clump), dP_obs = 96 s",
    6928997: "Bedding et al. 2011, Fig. 1: g-mode period spacing 77.1 s in an echelle diagram",
    9145955: "Vrard et al. 2016 (A&A 588, A87), Fig. 3: the stretched-period method",
}

# ---------------------------------------------------------------------------
# Light-curve preparation
# ---------------------------------------------------------------------------

#: Width (Gaussian sigma, days) of the smooth trend removed from each Kepler
#: quarter before the power spectrum. Signals faster than ~1 uHz (periods
#: < 12 d) pass almost untouched (power transfer 0.95 at 1 uHz), which keeps
#: the granulation background of the most luminous giants (numax ~ 5 uHz).
TREND_SIGMA_DAYS_LC = 5.0
#: Same for short-cadence data and the Sun, whose signals of interest are at
#: hundreds to thousands of uHz.
TREND_SIGMA_DAYS_SC = 1.0
#: Points further than this many robust standard deviations from zero after
#: detrending are dropped (Yu et al. 2018 clip at 4 sigma).
CLIP_SIGMA = 4.5

# ---------------------------------------------------------------------------
# Global oscillation parameters
# ---------------------------------------------------------------------------

#: The background is fitted from this frequency up to the Nyquist frequency.
#: Below ~1 uHz our detrending filter removes power.
FIT_NU_MIN_UHZ = 1.0
#: Width of the logarithmic bins the spectrum is averaged into for the fit
#: (in ln nu). The likelihood accounts for the number of raw bins averaged.
FIT_LOG_BIN = 0.005
#: Trial starting values of numax for the blind search (uHz).
NUMAX_SEARCH = (3.0, 270.0, 24)
#: emcee settings for the background + envelope posterior.
MCMC_WALKERS = 32
MCMC_STEPS = 2000
MCMC_BURN = 800
#: Number of chi-squared-perturbed spectra used for the Dnu uncertainty
#: (Huber et al. 2009 / Yu et al. 2018 use 200).
DNU_PERTURBATIONS = 200

#: Expected Dnu from numax (Stello et al. 2009, MNRAS 400, L80):
#: Dnu = 0.263 numax^0.772. Used only to centre the window in which the
#: autocorrelation peak is searched for, never as a result.
DNU_NUMAX_PRIOR = (0.263, 0.772)
#: Half-width of that window, as a fraction of the expected Dnu (Yu et al.
#: 2018 use Gaussian weights of 30%, and 70% for 15 < numax < 60 uHz).
DNU_WINDOW_FRAC = 0.30
DNU_WINDOW_FRAC_CLUMP = 0.70
DNU_WINDOW_CLUMP_RANGE = (15.0, 60.0)
