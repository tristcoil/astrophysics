"""Project configuration: paths, data sources, and every number we *assume*.

Each assumed value carries its source. Published *results* that we grade
ourselves against (Omega_m, sigma_8, S_8 from eROSITA, Planck, KiDS, DES...)
live separately in ``literature.py`` and are never used as fit inputs, with
one deliberate exception that the README states up front: the eRASS1 count
rate-mass relation (Ghirardini et al. 2024, Table 4), which this project
cannot measure because the weak-lensing shear profiles behind it are not public.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources (all public, no account needed)
# ---------------------------------------------------------------------------

#: eROSITA-DE Data Release 1 (eRASS1, western Galactic hemisphere). The plain
#: .fits links on this page answer 403; the .tgz bundles are the published route.
EROSITA_CATALOGUES = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/"
EROSITA_FILES = {
    "cosmology": ("BulbulE_DR1/erass1cl_cosmology_v1.1.fits.tgz", "erass1cl_cosmology_v1.1.fits"),
    "primary": ("BulbulE_DR1/erass1cl_primary_v3.2.fits.tgz", "erass1cl_primary_v3.2.fits"),
    "optical": ("KlugeM_DR1/eRASS1_clusters_optical.fits.tgz", "eRASS1_clusters_optical.fits"),
}

#: eRASS1 0.2-2.3 keV vignetted exposure map as a HiPS (Sanders, MPE). Only
#: Norder3 has an Allsky file: 768 tiles of 64x64 = a HEALPix map at order 9.
EXPOSURE_HIPS = "https://erosita.mpe.mpg.de/erodat/static/hips/eRASS1_024_ExposureMap_c010/"

#: DESI Legacy Surveys DR10 coverage, per band, as MOCs from the CDS MocServer.
LS_MOC = ("https://alasky.cds.unistra.fr/MocServer/query?"
          "ID=CDS/P/DESI-Legacy-Surveys/DR10/{band}&get=moc&order=8&fmt=fits")

#: HI4PI neutral-hydrogen column density, resampled by CDS hips2fits onto a
#: 0.25 deg equatorial grid (4 MB instead of the 50 MB HEALPix original).
HI4PI_NHI = ("https://alasky.cds.unistra.fr/hips-image-services/hips2fits?"
             "hips=CDS%2FP%2FHI4PI%2FNHI&width=1440&height=720&projection=CAR"
             "&coordsys=icrs&ra=180&dec=0&fov=360&format=fits")

#: Planck PSZ2 union catalogue (1,653 SZ detections), IRSA mirror of the PLA file.
PSZ2_URL = "https://irsa.ipac.caltech.edu/data/Planck/release_2/catalogs/HFI_PCCS_SZ-union_R2.08.fits.gz"

#: Published figures we read numbers from. Ghirardini et al. (2024) Fig. 4 is a
#: vector PDF inside the arXiv source bundle; the others are served one by one
#: by arXiv's HTML rendering of the papers.
ARXIV_EPRINT = "https://arxiv.org/e-print/{id}"
GHIRARDINI_ARXIV = "2402.08458"
GHIRARDINI_FIG4 = "figures_n/fig4.pdf"
BULBUL_FIG2A = "https://arxiv.org/html/2402.08452v1/figures/fig2a.png"
BULBUL_FIG2B = "https://arxiv.org/html/2402.08452v1/figures/fig2b.png"
CLERC_FIG = ("https://arxiv.org/html/2402.08457v1/images/"
             "GP_RateFLwoNH0223_nH_simbkg_Texp_z_Cosmo_cut_6-0_SEED000.png")

#: Project D3's own CMB fit, which this project turns into sigma_8 (optional:
#: the numbers are also copied into literature.py so this folder stands alone).
D3_FIT_PATH = PROJECT_DIR.parent / "09-planck-cmb" / "data" / "derived" / "tt_fit.json"

# ---------------------------------------------------------------------------
# The survey (Bulbul et al. 2024; Ghirardini et al. 2024)
# ---------------------------------------------------------------------------

AREA_COSMO_DEG2 = 12791.0  # cosmology footprint: eRASS1-DE x LS DR10-South
AREA_PRIMARY_DEG2 = 13116.0  # primary-catalogue footprint (adds LS DR9-North)
DEC_MAX_DR10_SOUTH = 32.375  # LS DR9/DR10 north-south split (Bulbul: "DEC <~ 32.5")
Z_MIN, Z_MAX = 0.1, 0.8  # photometric-redshift window of the cosmology sample
LAMBDA_MIN = 3.0  # optical richness cut of the cosmology sample

#: Reference sky conditions of the published selection-function curves
#: (Ghirardini et al. 2024, Fig. 4 caption).
NH_REF = 3.5e20  # cm^-2
BKG_REF = 5.2e-15  # erg/s/cm^2/arcmin^2 (stated, not modelled: see README)

# ---------------------------------------------------------------------------
# Cosmology held fixed in the fits (Ghirardini et al. 2024, Table 3 priors)
# ---------------------------------------------------------------------------
# Cluster counts barely constrain these, so the paper bounds them with priors:
# H0 ~ N(67.77, 0.6), Omega_b ~ U(0.046, 0.052), n_s ~ U(0.92, 1.0). We fix
# them at the prior centres and show in the robustness table that moving them
# across those ranges shifts sigma_8 by less than its error.
H0_FID = 67.77
OMEGA_B_FID = 0.049
NS_FID = 0.96
T_CMB = 2.7255  # K (Fixsen 2009; D3 measured 2.7250)
N_EFF = 3.046

# ---------------------------------------------------------------------------
# The count rate - mass relation (Ghirardini et al. 2024, Eqs. 22-24, Table 4)
# ---------------------------------------------------------------------------
# <ln(CR / CR_p) | M, z> = ln A_X + b_X(z) ln(M / M_p) + e_X(z)
#   b_X(z) = B_X + F_X ln((1+z)/(1+z_p))
#   e_X(z) = D_X ln(d_L(z)/d_L(z_p)) + E_X ln(E(z)/E(z_p)) + G_X ln((1+z)/(1+z_p))
# with log-normal intrinsic scatter sigma_X. "log" is natural: the richness
# normalisation ln A_lambda = 3.08 means ~22 galaxies at M_p, as expected.
CR_PIVOT = 0.1  # cts/s, unabsorbed 0.2-2.3 keV count rate within R500c
M_PIVOT = 2e14  # Msun, M500c
Z_PIVOT = 0.35
D_X = -2.0  # self-similar (fixed in the paper)
E_X = 2.0  # self-similar (fixed in the paper)

#: Posterior median and the larger of the two 68% half-widths, per parameter.
SCALING_PUBLISHED = {
    "A_X": (0.64, 0.06),
    "B_X": (1.38, 0.04),
    "F_X": (-0.33, 0.12),
    "G_X": (0.29, 0.13),
    "sigma_X": (0.98, 0.05),
}

#: The richness - mass relation (same table), used only for the lambda > 3 cut.
RICHNESS_PUBLISHED = {"lnA_l": 3.08, "B_l": 1.01, "C_l": -0.49, "D_l": 0.30, "sigma_l": 0.32}

#: Photometric-redshift error model (Table 3 prior centres): z_hat ~ N(b_z z, sigma_z (1+z)).
PHOTOZ_SIGMA = 0.005
PHOTOZ_BIAS = 1.005

# ---------------------------------------------------------------------------
# Analysis grids and bins
# ---------------------------------------------------------------------------

#: Observed-quantity bins for the Poisson likelihood.
CR_HAT_EDGES_LOG10 = (-2.4, 1.6, 20)  # 0.2 dex bins from 0.004 to 40 cts/s
Z_HAT_EDGES = (0.10, 0.80, 14)  # 0.05-wide photometric-redshift bins

#: Model integration ranges (Ghirardini et al. Sect. 4.1).
LNM_RANGE = (5e12, 5e15, 100)
Z_RANGE = (0.06, 0.86, 81)
LNCR_RANGE = (1e-3, 100.0, 151)

#: Parameter grid for the fixed-scaling-relation (Omega_m, sigma_8) scan.
OM_GRID = (0.18, 0.42, 61)
S8GRID_SIGMA8 = (0.72, 1.04, 65)

#: MCMC settings (emcee). The seed is passed to the sampler's own RNG: emcee
#: copies numpy's global state at construction, so seeding afterwards does nothing.
MCMC_WALKERS = 40
MCMC_STEPS = 3000
MCMC_BURN = 1000
MCMC_SEED = 20261004
