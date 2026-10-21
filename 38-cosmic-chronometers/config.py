"""Project configuration: paths, data sources, analysis choices, and reference values.

Everything tunable lives here, together with the published numbers we grade our
own results against. Published values are used only for comparison -- never as
inputs to a fit -- with three labelled exceptions: the H0 priors of the
"CC + local H0" variants, the DESI w0-wa parameters that drive one of the
*projected* futures, and the stellar-population calibration used to re-read
Moresco et al.'s 2012 measurement.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ATLAS_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
SOURCE_DIR = DATA_DIR / "sources"  # LaTeX tables pulled from arXiv sources
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, SOURCE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

C_KM_S = 299_792.458
#: 1 / (1 km/s/Mpc) in Gyr: one Mpc in km over one Julian Gyr in seconds.
HUBBLE_TIME_GYR = 3.0856775814913673e19 / 3.15576e16  # = 977.79...
#: Radiation today (photons + three massless neutrinos), Omega_r h^2, for
#: T_CMB = 2.7255 K and N_eff = 3.046. Matters only near the Big Bang.
OMEGA_R_H2 = 2.4728e-5 * (1.0 + 0.22711 * 3.046)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: Moresco's public cosmic-chronometer covariance primer (GitLab). The H(z)
#: table of his own 15 measurements, and the redshift-dependent systematic
#: errors (per cent) from Moresco et al. (2020): initial mass function,
#: stellar library, and stellar-population-synthesis model (with and without
#: the most discordant model). The fit notebook is kept as an answer key.
CC_REPO_RAW = "https://gitlab.com/mmoresco/CCcovariance/-/raw/master/"
CC_REPO_FILES = {
    "HzTable_MM_BC03.dat": "data/HzTable_MM_BC03.dat",
    "HzTable_MM_M11.dat": "data/HzTable_MM_M11.dat",
    "data_MM20.dat": "data/data_MM20.dat",
    "CC_fit.ipynb": "examples/CC_fit.ipynb",
}

#: arXiv LaTeX sources we read tables from: {local name: (arXiv id, file in the tarball)}.
ARXIV_SOURCES = {
    # Moresco et al. 2022, Living Reviews in Relativity 25, 6: Table 1 is the
    # 32-point cosmic-chronometer compilation; Table 2 holds its fits.
    "moresco2022_lrr.tex": ("2201.07241", "main.tex"),
    # Moresco et al. 2016, JCAP 05, 014: the BOSS measurements, the 30-point
    # table used for the first model-independent transition redshift.
    "moresco2016_boss.tex": ("1601.01701", "BOSS_sub_v2.tex"),
    # Moresco et al. 2012, JCAP 08, 006: D4000-redshift medians and the
    # D4000-age calibration behind the two SDSS points at z = 0.18 and 0.20.
    "moresco2012_d4000.tex": ("1201.3609", "Hz_v2_arxiv.tex"),
}

#: DESI DR2 BAO (DESI Collaboration 2025, arXiv:2503.14738), as packaged for
#: Cobaya: D_M/r_d, D_H/r_d and D_V/r_d with their full covariance.
DESI_RAW = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/"
DESI_FILES = {
    "desi_dr2_bao_mean.txt": "desi_gaussian_bao_ALL_GCcomb_mean.txt",
    "desi_dr2_bao_cov.txt": "desi_gaussian_bao_ALL_GCcomb_cov.txt",
}

#: Earlier atlas projects this one reads from (their own cached outputs).
D1_RESULTS = ATLAS_DIR / "05-supernova-cosmology" / "data" / "results.json"
D3_FIT = ATLAS_DIR / "09-planck-cmb" / "data" / "derived" / "tt_fit.json"
#: J3's SDSS Main Galaxy Sample with MPA-JHU spectral measurements (D4000_n).
J3_GALAXIES = ATLAS_DIR / "24-galaxy-life-cycle" / "data" / "mgs_mpajhu_dr17.npz"

# ---------------------------------------------------------------------------
# The cosmic-chronometer data set
# ---------------------------------------------------------------------------

#: Moresco's own measurements appear in the 2022 table rounded (z = 0.179 for
#: 0.1791, and 0.445 for 0.4497). Each is matched to the precise GitLab value by
#: (rounded z, H); the tolerance covers the rounding.
MATCH_DZ = 0.006
MATCH_DH = 0.6

#: Short labels for the methods in the 2022 table.
METHODS = {
    "D": "calibrated 4000 Å break",
    "F": "full-spectrum fitting",
    "L": "Lick indices",
}

#: Which systematic terms of Moresco et al. (2020) enter the covariance. The
#: primer notebook and the 2020 paper's recommendation for modern stellar
#: libraries: the IMF term plus the SPS-model term computed without the most
#: discordant model ("mod_ooo"). The stellar-library term is left out (the 2022
#: review: "overly conservative ... driven by a library that has now been
#: superseded"); it is kept as a variant to show what it would do.
SYSTEMATIC_TERMS = ("imf", "mod_ooo")

# ---------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------

#: Flat priors, matching the primer notebook's fit (50 < H0 < 100,
#: 0.01 < Omega_m < 0.99). The H0 grid is wider so that truncation is explicit.
H0_PRIOR = (50.0, 100.0)
OM_PRIOR = (0.01, 0.99)
OL_RANGE = (-0.5, 1.6)  # open LCDM: Omega_Lambda grid

#: Redshift grid for reconstructions and plots.
Z_MAX_PLOT = 2.5
Z_MAX_DATA = 1.965

#: Gaussian-process settings: hyperparameter prior box (natural log).
GP_LN_SIGMA_F = (1.0, 7.0)   # ~3 to ~1100 km/s/Mpc
GP_LN_LENGTH = (-1.0, 3.0)   # ~0.37 to ~20 in redshift
GP_N_HYPER_SAMPLES = 4000     # hyperparameter draws (one function each)
GP_N_FUNCTION_SAMPLES = 20000 # functions drawn at the optimum hyperparameters
RNG_SEED = 20261004

#: Piecewise-linear transition fit (Moresco et al. 2016): grid of z_t.
ZT_GRID = (0.0, 1.5, 1501)

# ---------------------------------------------------------------------------
# Reference values (comparison only unless stated)
# ---------------------------------------------------------------------------

#: D1 (05-supernova-cosmology) and D3 (09-planck-cmb), as published in their
#: READMEs. Used if their cached result files are missing; checked against those
#: files when present.
D1_FALLBACK = {
    "H0": (73.55, 1.02), "Om_flat": (0.332, 0.018), "Om_calibrated": (0.332, 0.018),
    "z_transition": (0.591, 0.044), "q0_flat": (-0.502, 0.027),
    "w_flat": (-0.914, 0.149), "age_gyr": (12.46, 0.24),
}
#: The most negative w D1's supernovae allow at 95% (2.5th percentile of the
#: marginal posterior on its (Omega_m, w) grid, read from D1's cache on
#: 2026-10-04: -1.241). Used for the "phantom energy within D1's limits" future.
D1_W_95_LOW = -1.24
D3_FALLBACK = {
    "H0": (66.66, 0.90), "Om": (0.3247, 0.0132), "rdrag": (146.89, 0.49),
    "age_gyr": (13.825, 0.033),
}

PLANCK18 = {"H0": (67.36, 0.54), "Om": (0.3153, 0.0073), "rdrag": (147.09, 0.26)}
SH0ES22 = {"H0": (73.04, 1.04)}

#: Answer keys for the fits.
PUBLISHED = {
    # Moresco's primer notebook (CC_fit.ipynb), 15 points, flat LCDM, emcee
    # medians and 16/84 percentiles. (value, +err, -err)
    "primer_full_H0": (65.995, 5.545, 5.591),
    "primer_full_Om": (0.355, 0.113, 0.088),
    "primer_diag_H0": (66.171, 3.770, 3.956),
    "primer_diag_Om": (0.354, 0.110, 0.088),
    # Moresco et al. 2022 (Living Reviews), Table 2: 32 points, flat LCDM.
    "lrr_H0": (66.5, 5.4, 5.4),
    "lrr_Om": (0.34, 0.08, 0.06),
    "lrr_chi2_red": 0.5,
    # Moresco et al. 2016: transition redshift from 30 CC points + an H0 prior.
    "m16_H0_prior": (73.0, 2.4),
    "m16_zt_piecewise": (0.4, 0.1, 0.1),
    "m16_zt_piecewise_95": (0.2, 0.8),
    "m16_zt_olcdm": (0.64, 0.11, 0.07),
    # Yu, Ratra & Wang 2018 (ApJ 856, 3): Gaussian process on H(z) data.
    "yu18_gp_H0": (67.0, 4.0),
    "yu18_zda_range": (0.33, 1.0),
    # Moresco et al. 2012 BC03 values at the two SDSS redshifts (GitLab table:
    # H, statistical part).
    "m12_sdss": {0.1791: (74.91, 3.8), 0.1993: (74.96, 4.9)},
}

#: H0 priors for the "chronometers + a local anchor" variants. Labelled inputs.
H0_ANCHORS = {
    "D1 (Cepheid-calibrated supernovae)": None,  # filled from D1 at run time
    "D3 (Planck CMB, assuming LCDM)": None,      # filled from D3 at run time
}

#: DESI DR2 (arXiv:2503.14738) Table V, w0waCDM posterior means. The
#: Pantheon+ combination drives the projected future (same supernovae as D1);
#: the other two show how much the projection depends on the supernova sample.
#: These are *inputs* to a labelled speculation, not data we fit.
DESI_W0WA = {
    "DESI+CMB+Pantheon+": {"Om": 0.3114, "H0": 67.51, "w0": -0.838, "wa": -0.62},
    "DESI+CMB+Union3": {"Om": 0.3275, "H0": 65.91, "w0": -0.667, "wa": -1.09},
    "DESI+CMB+DESY5": {"Om": 0.3191, "H0": 66.74, "w0": -0.752, "wa": -0.86},
}
DESI_W0WA_ERR = {"w0": (0.055, 0.055), "wa": (0.22, 0.19)}  # Pantheon+ combination
#: DESI DR2 + CMB + Pantheon+ with constant w: how phantom can w be?
DESI_WCDM = {"Om": 0.3047, "H0": 67.97, "w": (-0.995, 0.023)}

#: Caldwell, Kamionkowski & Weinberg 2003 (PRL 91, 071301): w = -3/2,
#: H0 = 70, Omega_m = 0.3 -> Big Rip in 22 Gyr; their timeline before the end.
CALDWELL03 = {
    "w": -1.5, "H0": 70.0, "Om": 0.3, "t_rip_from_now_gyr": 22.0,
    "before_rip": {"Galaxy clusters": "1 Gyr", "Milky Way": "60 Myr",
                   "Solar System": "3 months", "Earth": "30 minutes",
                   "Atoms": "1e-19 s"},
}

#: Bound systems for the Big Rip timeline: (label, orbital period in seconds).
#: Periods: Coma-like cluster (M = 1e15 Msun, R = 2 Mpc), the Sun's Galactic
#: orbit (R = 8.2 kpc, v = 233 km/s), Neptune, Earth, the Moon, and a grazing
#: orbit just above Earth's surface.
_YEAR = 3.15576e7
BOUND_SYSTEMS = [
    ("A rich galaxy cluster (Coma)", None),     # computed from M and R
    ("The Milky Way (the Sun's orbit)", None),  # computed from R and v
    ("Neptune's orbit", 164.79 * _YEAR),
    ("Earth's orbit", 1.0 * _YEAR),
    ("The Moon's orbit", 27.3217 * 86400.0),
    ("Earth itself (an orbit skimming its surface)", None),  # from M_E, R_E
]

#: Cosmic-chronometer measurements published since the 2022 compilation, from
#: the papers' abstracts (verified 2026-10-04). (z, H, err_minus, err_plus,
#: label, note). Errors combine statistical and systematic in quadrature.
#: Shown and compared, not fitted: several share galaxies with each other or
#: with the compilation, and their systematic covariance is not published in a
#: common form.
SINCE_2022 = [
    (0.80, 113.1, (15.1**2 + 11.3**2) ** 0.5, (15.1**2 + 29.1**2) ** 0.5,
     "Jiao+2023, LEGA-C", "same galaxies as the z = 0.75 Lick-index point"),
    (1.26, 135.0, 65.0, 65.0, "Tomasetti+2023, VANDELS", ""),
    (0.46, 88.48, (0.57**2 + 12.32**2) ** 0.5, (0.57**2 + 12.32**2) ** 0.5,
     "Loubser 2025, DESI DR1", "D4000; one of three"),
    (0.67, 119.45, (6.39**2 + 16.64**2) ** 0.5, (6.39**2 + 16.64**2) ** 0.5,
     "Loubser 2025, DESI DR1", ""),
    (0.83, 108.28, (10.07**2 + 15.08**2) ** 0.5, (10.07**2 + 15.08**2) ** 0.5,
     "Loubser 2025, DESI DR1", ""),
    (0.542, 66.0, (29.0**2 + 13.0**2) ** 0.5, (81.0**2 + 13.0**2) ** 0.5,
     "Tomasetti+2025, three clusters", "MUSE"),
    (0.12, 71.33, 4.20, 4.20, "Wang, Lei & Fan 2026, DESI DR1", ""),
    (0.65, 93.68, (28.27**2 + 10.67**2) ** 0.5, (28.27**2 + 10.67**2) ** 0.5,
     "Pradhan+2026, VIPERS photometry", ""),
    (0.61, 88.5, (12.6**2 + 8.1**2) ** 0.5, (6.7**2 + 8.1**2) ** 0.5,
     "Álvarez+2026, DESI DR1 LRGs", "reddest envelope"),
]

# ---------------------------------------------------------------------------
# SDSS clock (Moresco et al. 2012 re-read, and simple cuts on J3's galaxies)
# ---------------------------------------------------------------------------

#: Moresco et al. 2012's SDSS bins: two mass ranges, Delta z = 0.02 from 0.15.
SDSS_MASS_BINS = ((11.0, 11.25), (11.25, 11.5))
SDSS_Z_EDGES = (0.15, 0.17, 0.19, 0.21, 0.23)
#: Simple passive cuts we apply to J3's table (not Moresco's selection, which
#: used SED-template fits with 2MASS photometry): no H-alpha emission
#: (equivalent width above -3 Å; MPA-JHU's sign makes emission negative), red
#: u - r colour, and a low specific star-formation rate.
SDSS_SIMPLE_CUTS = {"ha_eqw_min": -3.0, "u_minus_r_min": 2.6, "log_ssfr_max": -11.5}
