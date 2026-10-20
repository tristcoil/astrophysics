"""Project configuration: paths, data sources, constants and published values.

Three tests of whether the universe really expands, each with its own data:

1. **Time dilation** -- light curves of the Type Ia supernovae of the Dark
   Energy Survey's 5-year sample (DES-SN5YR). In an expanding universe a
   supernova at redshift z is seen in slow motion, by a factor (1 + z).
2. **The CMB thermometer** -- published measurements of the temperature of the
   cosmic microwave background at high redshift (Sunyaev-Zel'dovich spectra of
   galaxy clusters, molecules and atoms in quasar absorption lines, water at
   z = 6.34). Expansion predicts T = T0 (1 + z).
3. **The Tolman test** -- surface brightnesses of early-type galaxies in three
   clusters at z = 0.76-0.92 against local galaxies (Lubin & Sandage 2001).
   Expansion predicts dimming by (1 + z)^4.

"Tired light" (Zwicky 1929) predicts no time dilation, an unchanging CMB
temperature, and surface-brightness dimming of only (1 + z)^1.

Published values below are used only to grade our results, never as fit inputs.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
DES_DIR = DATA_DIR / "des"
WHITE_DIR = DATA_DIR / "white2024"
PAPER_DIR = DATA_DIR / "papers"
DERIVED_DIR = DATA_DIR / "derived"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"
LITERATURE_DIR = PROJECT_DIR / "literature"

for _d in (DATA_DIR, DES_DIR, WHITE_DIR, PAPER_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: Speed of light, km/s (exact).
C_KM_S = 299_792.458

# ---------------------------------------------------------------------------
# 1. Supernova time dilation: DES-SN5YR
# ---------------------------------------------------------------------------

#: The DES-SN5YR data release on GitHub. Tag 1.3 is the release that Vincenzi
#: et al. (2024) and White et al. (2024) used; ``main`` has since moved on to a
#: 2026 recalibration ("Dovekie"). The photometry and SALT3 light-curve files
#: are byte-identical in both; the metadata are taken from 1.3 so that our
#: selection matches White et al.'s.
DES_TAG = "1.3"
DES_RAW = f"https://raw.githubusercontent.com/des-science/DES-SN5YR/{DES_TAG}/"

DES_FILES = {
    # One row per transient (19,706 of them): coordinates, redshifts, host
    # galaxy, pointers into the photometry file.
    "DES-SN5YR_DES_HEAD.FITS.gz": "0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz",
    # Every flux measurement of every transient (scene-modelling photometry).
    "DES-SN5YR_DES_PHOT.FITS.gz": "0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_PHOT.FITS.gz",
    "DES-SN5YR_DES.README": "0_DATA/DES-SN5YR_DES/DES-SN5YR_DES.README",
    # The SALT3 light-curve fits of the 1,635 cosmology supernovae, as model
    # light curves sampled every 2 days in each band, plus the data they fit.
    "DES5YR_SALT3_LCFIT.LCPLOT.gz": "0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz",
    # SALT3 fit results (peak date, stretch x1, colour c) and the probability
    # that each supernova is a Type Ia, for the Hubble-diagram sample.
    "DES-SN5YR_HD+MetaData.csv": "4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv",
    "README_distances.md": "4_DISTANCES_COVMAT/README.md",
}

#: White et al. (2024)'s analysis repository: their per-supernova light-curve
#: widths in each band are the answer key for our replication. They are Python
#: pickles of pandas DataFrames; ``fetch_data.load_white_widths`` opens them
#: with an unpickler that refuses anything except pandas/numpy reconstruction.
WHITE_RAW = "https://raw.githubusercontent.com/ryanwhite1/DES-Time-Dilation/main/"
WHITE_FILES = {f"{b}BAND_pickled_data": f"{b}BAND_pickled_data" for b in "griz"}
WHITE_FILES.update({f"{b}BAND_pickled_data_non_deredshifted":
                    f"{b}BAND_pickled_data_non_deredshifted" for b in "griz"})
WHITE_FILES.update({f"power_dispersion_pickle_{b}": f"power_dispersion_pickle_{b}"
                    for b in "griz"})

#: DECam filters as White et al. describe them: (central wavelength, width) in
#: nm. Only their ratios matter, for matching rest-frame wavelengths.
DES_FILTERS = {"g": (473.0, 150.0), "r": (642.0, 148.0), "i": (784.0, 147.0),
               "z": (926.0, 152.0)}
BANDS = ("g", "r", "i", "z")

# White et al. (2024) selection and method settings, reproduced as published.
#: Minimum probability of being a Type Ia (their ``PROBIa_BEAMS > 0.5``).
PROB_IA_MIN = 0.5
#: Photometric points with a flux error above this (FLUXCAL units, zero point
#: 27.5) are dropped.
FLUXERR_MAX = 20.0
#: Fraction of the target band's width allowed between rest-frame wavelengths
#: when building a reference light curve (their delta = 2^-4).
REF_DELTA = 1.0 / 16.0
#: Half-width (rest-frame days) of the reference window for each target point.
REF_TAU_DAYS = 2.0
#: Minimum number of target points and of reference points.
MIN_TARGET_POINTS = 5
MIN_REFERENCE_POINTS = 100
#: Monte Carlo resamplings for each width's error.
WIDTH_MC_ITER = 200

# Our own global fit: a smooth template in (rest-frame phase, rest-frame
# wavelength) fitted to all supernovae at once, with the time-dilation power b
# free.
#: Observer-frame window around peak used by the global fit (days). Wide enough
#: that no b in the scanned range pushes usable data out of it; the window
#: itself is fixed in the observer frame so it cannot prefer any b.
GLOBAL_WINDOW_OBS = (-40.0, 140.0)
#: Rest-frame phase range of the template (days); points that map outside it
#: for a given b are excluded from that b's chi-square, and the comparison
#: between b values is made on a fixed set of points (see timedilation.py).
TEMPLATE_PHASE = (-25.0, 85.0)
#: Rest-frame wavelength range of the template (nm). DES griz at z = 0.1-1.2.
TEMPLATE_WAVE = (200.0, 860.0)

# ---------------------------------------------------------------------------
# 2. CMB temperature at high redshift
# ---------------------------------------------------------------------------

#: Present-day CMB temperature, Fixsen (2009, ApJ 707, 916): COBE/FIRAS
#: recalibrated with WMAP. Project D3 (09-planck-cmb) confirms the blackbody.
T0_FIRAS = (2.72548, 0.00057)

#: The curated compilation (literature/tcmb_measurements.csv) is checked
#: against these sources by test_tcmb_table.py.
ARXIV_SOURCES = {
    # Lubin & Sandage Tolman test, Papers I, III and IV (TeX tables).
    "sandage_lubin_2001_paper1": "astro-ph/0102213",
    "lubin_sandage_2001_paper3": "astro-ph/0106563",
    "lubin_sandage_2001_paper4": "astro-ph/0106566",
    # T_CMB compilations and new measurements.
    "riechers_2022": "2202.00693",       # Nature 602, 58 (PDF only)
    "klimenko_2020": "2106.00119",       # Astron. Lett. 46, 715
    "noterdaeme_2011": "1012.3164",      # A&A 526, L7
    "avgoustidis_2016": "1511.04335",    # PRD 93, 043521
    "kotani_2025": "2509.20760",         # PKS 1830-211, z = 0.89
    "kotani_2026": "2602.22399",         # B0218+357, z = 0.68
    "muller_2013": "1212.5456",          # A&A 551, A109
    "luzzi_2009": "0909.2815",           # ApJ 705, 1122 (13 SZ clusters)
    # Time dilation.
    "white_2024": "2406.05050",          # MNRAS 533, 3365
}

#: Leiden Atomic and Molecular Database (Schoier et al. 2005): CO energy levels,
#: Einstein A coefficients and collision rates with H2 (Yang et al. 2010), for
#: the demonstration of why molecular thermometers read slightly warm.
LAMDA_CO_URL = "https://home.strw.leidenuniv.nl/~moldata/datafiles/co.dat"

# ---------------------------------------------------------------------------
# 3. The Tolman test
# ---------------------------------------------------------------------------

#: Hubble constant that Lubin & Sandage used for every linear radius. Any
#: value works as long as local and distant galaxies use the same one.
H0_TOLMAN = 50.0

#: The three distant clusters (Lubin & Sandage 2001, Paper III Table 4).
TOLMAN_CLUSTERS = {
    "Cl 1324+3011": {"z": 0.7565, "band": "I", "K": (0.71, 0.03)},
    "Cl 1604+4304": {"z": 0.8967, "band": "I", "K": (0.80, 0.07)},
    "Cl 1604+4321": {"z": 0.9243, "band": "R", "K": (1.89, 0.06)},
}
#: Petrosian eta radii they tabulate, and the two they trust (the HST
#: point-spread function biases eta < 1.7, Paper II).
ETAS = (1.0, 1.3, 1.5, 1.7, 2.0)
ETAS_FINAL = (1.7, 2.0)
#: Local early-type R - I colour that moves the R-band local relations to I.
R_MINUS_I_LOCAL = 0.62

# ---------------------------------------------------------------------------
# Published values (for comparison only)
# ---------------------------------------------------------------------------

#: White et al. 2024, MNRAS 533, 3365: time dilation from 1,504 DES SNe Ia.
WHITE24 = {
    "b": (1.003, 0.005), "b_sys": 0.010, "n_sne": 1504,
    "n_per_band": {"g": 361, "r": 1380, "i": 1465, "z": 1381},
    "linear": {"slope": (0.988, 0.016), "intercept": (0.020, 0.024)},
    "chi2_nu_average": 1.441,
}
#: Earlier supernova time-dilation measurements.
GOLDHABER01_B = (1.07, 0.06)      # 35 SNe, 0.3 < z < 0.7 (photometry)
BLONDIN08_B = (0.97, 0.10)        # spectral ageing of 13 SNe

#: T_CMB(z) = T0 (1+z)^(1-beta) fits.
NOTERDAEME11_BETA = (-0.007, 0.027)       # SZ + atomic carbon + CO
RIECHERS22_BETA = (0.0034, 0.0081, 0.0073)  # value, +err, -err (their ED Table 1)
KOTANI26_BETA = (0.0039, 0.0074, 0.0082)  # ~80 measurements + B0218+357
KLIMENKO20_BETA = (0.010, 0.013)          # Mol + SZ + CI + CO
KLIMENKO20_T0 = (2.719, 0.009)            # T0 from high-z data alone, beta = 0

#: Lubin & Sandage 2001 Paper IV: exponent n of the surface-brightness
#: dimming 2.5 log10 (1+z)^n, averaged over eta = 1.7 and 2.0, q0 = 1/2.
LUBIN_SANDAGE01 = {
    "n_R": (2.59, 0.17), "n_I": (3.37, 0.13),
    "n_R_q0_0": 2.28, "n_I_q0_0": 3.06, "n_R_q0_1": 2.81, "n_I_q0_1": 3.55,
    "n_R_tired": (1.61, 0.13), "n_I_tired": (2.27, 0.12),
}
