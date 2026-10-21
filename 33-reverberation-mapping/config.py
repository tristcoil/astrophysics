"""Project configuration: paths, physical constants and every adopted input.

Everything we *assume* lives here, each value with its source: NGC 5548's
redshift, distance, host-galaxy starlight and stellar velocity dispersion, and
the analysis settings of the papers we reproduce. Published *results* we grade
ourselves against live separately in ``literature.py`` and never enter a fit.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
CACHE_DIR = DATA_DIR / "cache"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR, CACHE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants (cgs)
# ---------------------------------------------------------------------------

C_CGS = 2.99792458e10  # speed of light, cm/s
G_CGS = 6.67430e-8  # Newton's constant, cm^3 g^-1 s^-2 (CODATA 2018)
MSUN_G = 1.98841e33  # solar mass, g
DAY_S = 86400.0
LIGHT_DAY_CM = C_CGS * DAY_S  # 2.59e15 cm
PC_CM = 3.0856775814913673e18
MPC_CM = 1e6 * PC_CM
AU_CM = 1.495978707e13
L_EDD_PER_MSUN = 1.26e38  # Eddington luminosity per solar mass, erg/s (pure hydrogen)

# ---------------------------------------------------------------------------
# NGC 5548
# ---------------------------------------------------------------------------

#: Redshift of NGC 5548 (de Vaucouleurs et al. 1991, as used by Bentz et al.
#: 2013, Pei et al. 2017 and the AGN Black Hole Mass Database).
Z_NGC5548 = 0.01718

#: Luminosity distance adopted by Bentz et al. (2013) and Pei et al. (2017)
#: for NGC 5548's AGN luminosities: 75 Mpc. (The AGN mass database quotes
#: 72.5 Mpc for H0 = 72; the published log L values use 75 Mpc.)
DL_NGC5548_MPC = 75.0

#: Milky Way reddening toward NGC 5548 (Schlafly & Finkbeiner 2011), as used by
#: Pei et al. 2017 when turning fluxes into luminosities. A_5100 = 3.1 * 1.04 E(B-V)
#: is computed in ``physics.galactic_extinction_5100``.
EBV_NGC5548 = 0.017

#: Host-galaxy starlight at 5100 A through the AGN STORM ground-based aperture
#: (5" x 15", PA 0), measured on Bentz et al. 2013's AGN-free HST image (Pei et
#: al. 2017, section 3.2), in 1e-15 erg/s/cm^2/A.
HOST_FLUX_STORM = 4.52
HOST_FLUX_STORM_ERR = 0.45

#: Bulge stellar velocity dispersion of NGC 5548: 195 +/- 13 km/s (Woo et al.
#: 2010, ApJ 716, 269, from Keck spectra of the Ca triplet).
SIGMA_STAR_NGC5548 = 195.0
SIGMA_STAR_NGC5548_ERR = 13.0

# ---------------------------------------------------------------------------
# AGN STORM (2014) analysis settings, as in Pei et al. (2017), section 4.1
# ---------------------------------------------------------------------------

#: Lag grid of the ICCF, observed frame (days).
STORM_LAG_RANGE = (-20.0, 40.0)
STORM_LAG_STEP = 0.25
#: The last 10 MDM epochs have inconsistent flux calibration; Pei et al. cut the
#: light curves at THJD = 6828.75 (THJD = HJD - 2,450,000), leaving 133 epochs.
STORM_TRUNCATE_THJD = 6828.75
#: The campaign is split in two halves at THJD 6747 (Pei et al.), because the
#: lines stopped following the continuum in the second half ("BLR holiday").
STORM_SPLIT_THJD = 6747.0
#: Number of flux-randomization / random-subset-selection realizations.
N_MC = 1000

#: Goad et al. (2016) measured each UV line's lag and its normalization over
#: the first 75 days of the HST campaign (THJD 6690.6 to 6765.6), before the
#: anomaly began.
HST_START_THJD = 6690.6
HOLIDAY_REFERENCE_DAYS = 75.0

#: Narrow-line flux under each broad line, removed before computing
#: responsivities (Goad et al. 2016, Table 1, from low-state spectra of
#: Crenshaw et al. 1993), in 1e-13 erg/s/cm^2; for H-beta the narrow component
#: from Pei et al.'s fit to the mean MDM spectrum (48.4e-15).
NARROW_FLUX = {"Lya": 8.9, "SiIV": 1.2, "CIV": 7.0, "HeII": 1.2, "Hb": 0.484}

#: Segments of the H-beta light curve defined by Pei et al. (2017, section
#: 4.2): normal | decorrelating | suppressed | recovering | normal again.
HB_SEGMENTS_THJD = (6743.0, 6772.0, 6812.0, 6827.0)

# ---------------------------------------------------------------------------
# AGN Watch: NGC 5548 1988-2008 (Peterson et al. 2013 revised light curves)
# ---------------------------------------------------------------------------

#: Seasonal ICCF settings: lags searched 0-60 d, 0.5-d steps, as the published
#: yearly analyses (Peterson et al. 2002, 2004) searched well past the ~20-d
#: lags. Negative lags are searched too, so a failed season is visible as such.
WATCH_LAG_RANGE = (-30.0, 80.0)
WATCH_LAG_STEP = 0.5

# ---------------------------------------------------------------------------
# SDSS-RM (Shen et al. 2024)
# ---------------------------------------------------------------------------

#: ICCF settings for the SDSS-RM light curves: pyPETAL's defaults as used for
#: the final release (2-day grid, threshold 0.8, minimum r 0.2), over each
#: quasar's own lag search range (``LAG_SEARCH`` in their summary table).
SDSSRM_LAG_STEP = 2.0
SDSSRM_N_MC = 500
#: Alias weighting (Grier et al. 2017; 2019): P(tau) = [N(tau)/N(0)]^k with
#: k = 2, gaps longer than 30 days count as seasonal gaps, smoothed with a
#: 15-day Gaussian before picking the primary peak.
ALIAS_K = 2
ALIAS_GAP_DAYS = 30.0
ALIAS_SMOOTH_DAYS = 15.0
#: Doublet splitting of Mg II 2796/2803 in velocity, removed in quadrature
#: from the rms line dispersion (Shen et al. 2024).
MGII_SPLIT_KMS = 750.0

# ---------------------------------------------------------------------------
# Black hole masses
# ---------------------------------------------------------------------------

#: Virial factors for sigma_line(rms) virial products (log10 f):
#:   0.65 +/- 0.12: Woo et al. 2015 (M-sigma, quiescent galaxies of Kormendy & Ho 2013)
#:   0.62 +/- 0.07: Shen et al. 2024 (30 BLR dynamical models)
#: These are what the papers used; our own calibrations are in ``masses.py``.
LOG_F_WOO15 = 0.65
LOG_F_SHEN24 = 0.62
