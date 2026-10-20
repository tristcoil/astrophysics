"""Paths, sample cuts and the published numbers we check against (L4).

Every number in REFERENCE is quoted from the paper named beside it; the
analysis never uses them as inputs, only prints them next to our own results.
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGURES = HERE / "figures"
GALAH_DIR = DATA / "galah"
EU_DIR = GALAH_DIR / "europium"
KN_DIR = DATA / "kilonova"
LIT_DIR = DATA / "literature"
APOGEE_DIR = DATA / "apogee"
CACHE_DIR = DATA / "cache"
RESULTS_PATH = DATA / "results.json"

# ---------------------------------------------------------------------------
# Remote sources (fetch_data.py)
# ---------------------------------------------------------------------------
GALAH_BASE = "https://cloud.datacentral.org.au/teamdata/GALAH/public/GALAH_DR4/catalogs"
GALAH_FILES = {
    "allstar": "galah_dr4_allstar_240705.fits",          # 723 MB, 917,588 stars
    "dynamics": "galah_dr4_vac_dynamics_240705.fits",    # 155 MB, orbits (McMillan 2017)
    "zeropoints": "galah_dr4_zeropoints_240705.fits",    # adopted solar zero points
}
#: Kane et al. (2026, MNRAS; arXiv:2512.02125): europium re-derived for DR4.
EU_VAC_BASE = GALAH_BASE + "/galah_dr4_vac_europium/galah_dr4_vac_europium_v1"
EU_VAC_FILES = ["cnn_predictions.fits", "golden_sample.fits",
                "korg_secondary_catalogue.fits", "training_and_validation_dataset.fits"]

#: Villar et al. (2017, ApJL 851, L21): the unified AT2017gfo photometry is
#: Table 3, shipped in the arXiv source as phot_table.dat.
VILLAR_ARXIV = "1710.11576"
#: Asplund et al. (2009, ARA&A 47, 481): solar abundances, Table 1 in sun.tex.
ASPLUND_ARXIV = "0909.0948"
#: JINAbase (Abohalima & Frebel 2018), updated compilation (Mardini & Frebel).
JINABASE_URL = ("https://raw.githubusercontent.com/Mohammad-Mardini/JINAbase-updated/"
                "main/JINAbase_2021.csv")
#: Bisterzo et al. (2014, ApJ 787, 10): solar s-process fractions by element.
BISTERZO_VIZIER = "J/ApJ/787/10"
#: APOGEE DR17 allStar (Abdurro'uf et al. 2022) via VizieR III/286.
APOGEE_VIZIER = "III/286/catalog"
#: SVO Filter Profile Service (Rodrigo et al. 2012), as Villar et al. used.
SVO_URL = "http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={}"

# ---------------------------------------------------------------------------
# GALAH DR4 sample (Buder et al. 2025 best practices: S/N, flag_sp, flag_X_fe)
# ---------------------------------------------------------------------------
#: Field programmes only. Globular- and open-cluster pointings are dropped:
#: cluster stars carry Na-O and Mg-Al anomalies and would weigh the trends
#: toward a handful of clusters.
FIELD_SURVEYS = ("galah_main", "galah_bright", "galah_faint", "galah_phase2",
                 "k2_hermes", "tess_hermes")
SNR_MIN = 30.0                    # snr_px_ccd3 > 30 (DR4 recommendation)

#: Elements GALAH DR4 measures, in the order of the periodic table.
ELEMENTS = ["C", "N", "O", "Na", "Mg", "Al", "Si", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
            "Fe", "Co", "Ni", "Cu", "Zn", "Rb", "Sr", "Y", "Zr", "Mo", "Ru", "Ba", "La",
            "Ce", "Nd", "Sm", "Eu"]

#: "Compare like with like": two stellar classes analysed separately.
DWARFS = dict(logg=(3.5, 4.5), teff=(4200.0, 6700.0))   # Griffith et al. 2022's cuts
GIANTS = dict(logg=(1.0, 3.2), teff=(4000.0, 5500.0))

# ---------------------------------------------------------------------------
# Gaia-Sausage-Enceladus
# ---------------------------------------------------------------------------
#: Feuillet et al. (2021) selection as used on GALAH by Buder et al. (2022):
#: -500 < L_Z < 500 kpc km/s and 30 < sqrt(J_R / kpc km/s) < 55.
GSE_LZ = (-500.0, 500.0)
GSE_SQRT_JR = (30.0, 55.0)
#: Project H3, whose boxes (Helmi et al. 2018; Koppelman et al. 2019b) we
#: re-apply to H3's own Gaia sample and cross-match into GALAH.
H3_GLOB = "*-galactic-archaeology"

# ---------------------------------------------------------------------------
# AT2017gfo (Villar et al. 2017's choices)
# ---------------------------------------------------------------------------
KN_T_MERGER_MJD = 57982.528524     # GW170817 at 12:41:04.4 UTC
KN_DISTANCE_MPC = 39.5             # NGC 4993, as adopted by Villar et al.
KN_REDSHIFT = 0.009727             # NGC 4993 (OKC entry used by MOSFiT)
KN_EBV_MW = 0.105                  # Schlafly & Finkbeiner 2011
KN_RV = 3.1
KN_KAPPA = {"blue": 0.5, "purple": 3.0, "red": 10.0}     # cm^2/g, fixed

# ---------------------------------------------------------------------------
# Published values we compare with
# ---------------------------------------------------------------------------
REFERENCE = {
    "villar2017": {   # Table 2, spherical three-component model
        "blue": dict(m=0.020, m_err=(0.001, 0.001), v=0.266, v_err=(0.008, 0.008), T=674.0),
        "purple": dict(m=0.047, m_err=(0.002, 0.001), v=0.152, v_err=(0.005, 0.005), T=1308.0),
        "red": dict(m=0.011, m_err=(0.001, 0.002), v=0.137, v_err=(0.021, 0.025), T=3745.0),
        "sigma": 0.242,
        "n_points": 647,
        # Their Milky Way production-rate argument (Sect. 5), with the 2017 BNS
        # rate 1540 Gpc^-3 yr^-1 -> ~150 Myr^-1 per Milky Way:
        "mw_rate_per_myr": 150.0,
        "needed_light_msun": 5e-3,     # A <~ 140 r-process per merger
        "needed_heavy_msun": 7e-4,     # A >~ 140
    },
    "griffith2022_fcc": {   # GALAH DR3, Table 1: f_cc on the high-Ia sequence at [Mg/H] = 0
        "O": 1.13, "Si": 0.73, "Ca": 0.60, "Ti": 0.78, "Na": 0.53, "Al": 0.91, "K": 0.66,
        "Sc": 0.63, "Cr": 0.50, "Mn": 0.37, "Fe": 0.51, "Ni": 0.53, "Cu": 0.60, "Zn": 0.73,
        "Y": 0.29, "Ba": 0.31,
    },
    "griffith2022_plateau": 0.30,      # their assumed low-Ia [Mg/Fe] plateau
    "buder2022_gse_dynamical": {       # GALAH DR3, Table 4, dynamical selection (16/50/84)
        "fe_h": (-1.16, 0.32, 0.41), "mg_fe": (0.12, 0.12, 0.10), "eu_fe": (0.48, 0.17, 0.13),
        "ba_fe": (0.25, 0.36, 0.29), "na_fe": (-0.22, 0.17, 0.13),
    },
    "gse_knee": {
        "mackereth2019": -1.2,         # APOGEE, as quoted by Horta et al. (2023)
        "literature_range": (-2.5, -1.5),   # Matsuno+19, Monty+20, Aguado+21 (Horta+23)
    },
    "gwtc4_bns_rate": (7.6, 250.0),    # Gpc^-3 yr^-1, 90%, union of two models (GWTC-4.0 population, 2025)
    "gwtc4_bns_rate_simple": (13.0, 170.0),   # their "Simple Uniform BNS" model
    "bns_rate_2017": 1500.0,           # as used by Villar et al. 2017 (1500 +3200 -1220)
    "kane2026_golden": 54068,          # Eu VAC golden-sample size
}

#: Milky Way stellar mass (Bland-Hawthorn & Gerhard 2016): 5 +- 1 x 10^10 Msun.
MW_STELLAR_MASS = (5.0e10, 1.0e10)
#: Interstellar gas (HI + H2) within 30 kpc: 8.0 x 10^9 Msun (Nakanishi &
#: Sofue 2016, PASJ 68, 5).
MW_GAS_MASS = 8.0e9
#: Local stellar mass density measured in project J3 (log10 Msun/Mpc^3).
J3_LOG_RHO_STAR = 8.588
#: How long the Galaxy has been forming stars (Gyr).
MW_AGE_GYR = 12.0
