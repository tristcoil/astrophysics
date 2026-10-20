"""Project configuration: paths, data sources, physics conventions, the inflation
models on trial, and the published values our results are graded against.

Published numbers live here only for comparison at the end, never as inputs to
a fit -- with labelled exceptions where a published *measurement* is the input
by design: the amplitude of the primordial spectrum A_s (it sets the energy
scale of each model, not its n_s or r), and the BK18 likelihood's own data and
foreground priors (we re-implement that likelihood; its data are measurements).
"""

from __future__ import annotations

import os

# One OpenMP thread per process, set before CAMB (and its bundled libgomp) is
# loaded. CAMB run in the parent with a team of OpenMP threads leaves a
# fork-started worker pool deadlocked: the children wait forever on a thread
# team that was not copied into them. `run.py --force` hung this way inside
# step_bk18_validation, after bk18_mcmc had called CAMB in the parent.
# threadpool_limits() below cannot fix it: it only reaches libraries that are
# already loaded, and CAMB is imported later.
os.environ["OMP_NUM_THREADS"] = "1"

from pathlib import Path

from threadpoolctl import threadpool_limits

# One BLAS thread per process. The matrices here are small, and on a shared,
# busy machine OpenBLAS's 16 threads made one 594x594 inverse take 4.5 s
# instead of 0.05 s. Parallel work uses processes instead.
threadpool_limits(1)

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
DERIVED_DIR = DATA_DIR / "derived"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

BK18_DIR = DATA_DIR / "bk18"
ACT_DIR = DATA_DIR / "act_dr6"
SPT_DIR = DATA_DIR / "spt3g_d1"
PLANCK_DIR = DATA_DIR / "planck2018"
SKY_DIR = DATA_DIR / "planck_sky"
LITERATURE_DIR = DATA_DIR / "literature"

for _d in (DATA_DIR, FIGURE_DIR, DERIVED_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-atlas-P1/1.0 (educational; python-urllib)"

# ---------------------------------------------------------------------------
# Earlier atlas projects this one builds on
# ---------------------------------------------------------------------------
#
# Raw files a sibling project already downloaded are hard-linked from its cache
# after a checksum check (no second download, no extra disk). Results a sibling
# *computed* are read from its results file; if that file is missing, the value
# quoted in its README is used instead (``*_FALLBACK``).

D3_DIR = ROOT_DIR / "09-planck-cmb"            # D3: Planck CMB, CAMB fit (our own n_s)
O1_DIR = ROOT_DIR / "65-friedmann-universes"   # O1: curvature from BAO + CMB

D3_FIT = D3_DIR / "data" / "derived" / "tt_fit.json"
#: D3's TT-only CAMB fit (tau fixed at 0.0522), as quoted in its README.
D3_NS_FALLBACK = {"value": 0.9629, "error": 0.0051}

P2_DIR = ROOT_DIR / "68-thermal-history"       # P2: the Standard Model's g*(T), properly
P2_CURVES = P2_DIR / "data" / "derived" / "curves.npz"

O1_RESULTS = O1_DIR / "data" / "results.json"
#: O1's DESI BAO + CMB curvature, as quoted in its README.
O1_OK_FALLBACK = {"mean": 0.0021, "std": 0.0012}

#: Planck 2018 SMICA half-mission CMB maps and the common temperature mask
#: (D3 already has them; hard-linked, so no second download).
SMICA_URL = "https://irsa.ipac.caltech.edu/data/Planck/release_3/all-sky-maps/maps/component-maps/cmb/"
SMICA_FILES = ("COM_CMB_IQU-smica_2048_R3.00_hm1.fits", "COM_CMB_IQU-smica_2048_R3.00_hm2.fits")
MASK_URL = "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/masks/"
MASK_FILE = "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits"
SIBLING_COPIES = {SKY_DIR / name: D3_DIR / "data" / "planck" / name
                  for name in SMICA_FILES + (MASK_FILE,)}

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: BICEP/Keck BK18 (BK-XIII: Ade et al. 2021, PRL 127, 151301), LAMBDA's mirror
#: of the October 2021 release (bicepkeck.org itself is slow and often times out).
BK18_URL = "https://lambda.gsfc.nasa.gov/data/suborbital/BICEPK_2021/"
BK18_FILES = (
    "BK18_bandpowers_20210607.txt",    # all BB bandpowers in the likelihood (Fig. 15)
    "BK18_r_likelihood_20210607.txt",  # BK18's own r likelihood (Fig. 4) -- answer key
    "BK18_components_20210607.txt",    # CMB / dust / synchrotron per ell bin (Fig. 16)
    "BK18_Nl_fsky_20210607.txt",       # noise bandpowers and effective sky fraction
    "BK18_cosmomc.tgz",                # the likelihood's data: what we re-implement
    "rns_code.tgz",                    # BK18's own r-n_s figure code (cross-check)
)

#: ACT DR6.02 MCMC chains on LAMBDA (Louis et al. 2025; Calabrese et al. 2025).
#: Tensor-mode runs ("r"): every one includes BK18. Names: p = Planck, act =
#: ACT DR6 (lite = foreground-marginalised likelihood), bk = BK18, l = CMB
#: lensing, b = DESI BAO.
ACT_URL = "https://lambda.gsfc.nasa.gov/data/act/chains/"
ACT_R_CHAINS = {
    "planck-bk_r_camb": "Planck + BK18",
    "planck-bk-l-b_r_camb": "Planck + BK18 + lensing + BAO",
    "p-actlite-bk_r_camb": "P-ACT + BK18",
    "p-actlite-bk-l-b_r_camb": "P-ACT-LB + BK18",
}
ACT_LCDM_CHAINS = {
    "actlite_lcdm_camb": "ACT DR6",
    "p-actlite_lcdm_camb": "P-ACT",
    "p-actlite-l-b_lcdm_camb": "P-ACT-LB",
    "planck_lcdm_camb": "Planck (ACT's run)",
}

#: SPT-3G D1 chains on LAMBDA (Camphuis et al. 2025, arXiv:2506.20707): zip
#: archives read member by member with HTTP range requests. The CMB-SPA runs
#: (SPT + Planck + ACT) come from the February 2026 bug-fix release.
SPT_URL = "https://lambda.gsfc.nasa.gov/data/suborbital/SPT/spt_3g_d1/"
SPT_MAIN_ZIP = "chains_ttteee_winter1920.zip"
SPT_MAIN_INNER = "Compressed_Chains/LCDM.zip"   # stored (uncompressed) inside the main zip
SPT_FIXED_ZIP = "CMB_SPA_DESI_Fixed.zip"
SPT_CHAINS = {  # key -> (archive, folder, label)
    "spt3g": ("main", "LCDM/S1920lite_MPP_No_OLE", "SPT-3G D1"),
    "spt3g_desi": ("main", "LCDM/S1920lite_MPP_BAODR2_No_OLE", "SPT-3G D1 + DESI"),
    "cmbspa": ("fixed", "LCDM/S1920lite_MPP_PACT_PR4lens_actdr6lens_No_OLE", "CMB-SPA"),
    "cmbspa_desi": ("fixed", "LCDM/S1920lite_MPP_PACT_PR4lens_actdr6lens_BAODR2_No_OLE",
                    "CMB-SPA + DESI"),
}

#: Planck 2018 cosmological-parameter grid (Planck 2018 VI and X), read member
#: by member out of the 11.8 GB zip on IRSA's mirror of the Planck Legacy Archive.
PLANCK_GRID_URL = ("https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/"
                   "cosmoparams/COM_CosmoParams_fullGrid_R3.01.zip")
PLANCK_CHAINS = {  # key -> (grid folder, label)
    "planck18_lcdm": ("base/plikHM_TTTEEE_lowl_lowE_lensing", "Planck 2018"),
    "planck18_tt_lcdm": ("base/plikHM_TT_lowl_lowE", "Planck 2018 TT"),
    "planck18_r": ("base_r/plikHM_TTTEEE_lowl_lowE_lensing", "Planck 2018"),
    "planck18_r_bk15": ("base_r/plikHM_TTTEEE_lowl_lowE_BK15_lensing", "Planck 2018 + BK15"),
}

# ---------------------------------------------------------------------------
# Units and physical constants
# ---------------------------------------------------------------------------

#: Reduced Planck mass M_Pl = (8 pi G)^(-1/2) in GeV (CODATA 2018 G, hbar, c).
M_PL_GEV = 2.435_32e18
#: kpc and Mpc in metres; Mpc^-1 in GeV (hbar c = 1.973 269 804e-16 GeV m).
MPC_M = 3.085_677_581_491_367e22
HBARC_GEV_M = 1.973_269_804e-16
MPC_INV_GEV = HBARC_GEV_M / MPC_M  # 1 Mpc^-1 = 6.39e-39 GeV
#: Planck time and the age-of-universe conversion.
T_PLANCK_S = 5.391_247e-44
K_B_GEV_PER_K = 8.617_333_262e-14
C_KM_S = 299_792.458

#: Pivot scale for n_s and r: 0.05 Mpc^-1 (Planck, BK18, ACT and SPT all use it
#: for n_s; BK18 and ACT for r. Planck 2018 X quotes r at 0.002 Mpc^-1 too).
K_PIVOT = 0.05
K_PIVOT_PLANCK_R = 0.002

#: The primordial amplitude, ln(10^10 A_s) = 3.044 +- 0.014 (Planck 2018 VI,
#: TT,TE,EE+lowE+lensing). Used only to set each model's energy scale V^(1/4).
LN10_10_AS = 3.044
#: Today's photon temperature (Fixsen 2009) and the entropy degrees of freedom
#: today (photons + three neutrino species, 2 + 7/8 * 6 * 4/11 = 3.909).
T_CMB_K = 2.7255
G_S_TODAY = 3.909
#: Planck 2018 (TT,TE,EE+lowE+lensing) background used for the horizon and
#: flatness problems and for the reheating e-fold count.
H0_FID = 67.36
OMBH2_FID = 0.02237
OMCH2_FID = 0.1200
TAU_FID = 0.0544
NS_FID = 0.9649
MNU_FID = 0.06

# ---------------------------------------------------------------------------
# The models on trial (M_Pl = 1 throughout; potentials up to their overall scale)
# ---------------------------------------------------------------------------
#
# Every model is solved exactly (background + Mukhanov-Sasaki modes) for each
# parameter value and N_* = 50-60 e-folds before the end of inflation.

N_STAR_RANGE = (50.0, 60.0)
N_STAR_GRID = (50.0, 52.5, 55.0, 57.5, 60.0)

MODELS = {
    "phi23": {"label": "φ^(2/3) (axion monodromy)", "kind": "power", "p": 2.0 / 3.0},
    "phi1": {"label": "φ (linear monodromy)", "kind": "power", "p": 1.0},
    "phi2": {"label": "m²φ² (chaotic)", "kind": "power", "p": 2.0},
    "phi4": {"label": "λφ⁴ (chaotic)", "kind": "power", "p": 4.0},
    "natural": {"label": "natural inflation", "kind": "natural",
                "param": "f", "range": (3.0, 100.0)},
    "starobinsky": {"label": "Starobinsky R² / Higgs", "kind": "starobinsky"},
    "hilltop4": {"label": "hilltop (quartic)", "kind": "hilltop", "p": 4.0,
                 "param": "mu", "range": (5.0, 100.0)},
    "tmodel": {"label": "α-attractor (T-model)", "kind": "tmodel",
               "param": "alpha", "range": (0.01, 30.0)},
}

# ---------------------------------------------------------------------------
# Published values: the answer key (never inputs)
# ---------------------------------------------------------------------------

PUBLISHED = {
    # Planck 2018 VI (arXiv:1807.06209) and X (arXiv:1807.06211)
    "planck18_ns": {"value": 0.9649, "error": 0.0042,
                    "ref": "Planck 2018 VI Eq. 19 / X Eq. 13: TT,TE,EE+lowE+lensing, LCDM"},
    "planck18_tt_ns": {"value": 0.9626, "error": 0.0057, "ref": "Planck 2018 VI Table 2, TT+lowE"},
    "planck18_r002": {"value": 0.10, "ref": "Planck 2018 X Eq. 29 (TT+lowE+lensing; unchanged "
                                           "with high-l polarization), r_0.002 < 0.10 (95%)"},
    "planck18_bk15_r002": {"value": 0.056, "ref": "Planck 2018 X Eq. 31: TT,TE,EE+lowE+lensing+BK15"},
    "planck18_v_quarter_gev": {"value": 1.6e16, "ref": "Planck 2018 X Eq. 32: V* < (1.6e16 GeV)^4 "
                                                       "(95%, + BK15)"},
    # BK18 = BK-XIII (Ade et al. 2021, PRL 127, 151301; arXiv:2110.00483)
    "bk18_r95": {"value": 0.036, "ref": "BK-XIII abstract: r_0.05 < 0.036 (95%)"},
    "bk18_r_peak": {"value": 0.014, "plus": 0.010, "minus": 0.011,
                    "ref": "BK-XIII: r_0.05 = 0.014 +0.010 -0.011"},
    "bk18_L0": {"value": 0.46, "ref": "BK-XIII: zero-to-peak likelihood ratio"},
    "bk18_sigma_r": {"value": 0.009, "ref": "BK-XIII: sigma(r) from r = 0 simulations"},
    "bk18_Ad": {"value": 4.4, "plus": 0.8, "minus": 0.7,
                "ref": "BK-XIII Fig. 4 text: A_d at 353 GHz, l = 80 [uK^2]"},
    "bk18_betad": {"value": 1.49, "plus": 0.13, "minus": 0.12, "ref": "BK-XIII Fig. 4 text"},
    "bk18_Async95": {"value": 1.4, "ref": "BK-XIII: A_sync < 1.4 uK^2 (95%, 23 GHz, l = 80)"},
    "bk18_alphad_ml": {"value": -0.66, "ref": "BK-XIII: maximum-likelihood alpha_d"},
    "bk18_planck_bao_r95": {"value": 0.035, "ref": "BK-XIII Fig. 5: Planck 2018 + BK18 + BAO"},
    "tristram_r95": {"value": 0.032, "ref": "Tristram et al. 2022 (PRD 105, 083524): PR4 + BK18 + "
                                            "BAO + lensing"},
    # ACT DR6: Louis et al. 2025 (arXiv:2503.14452) Table 5; Calabrese et al. 2025 (2503.14454)
    "act_ns": {"value": 0.9666, "error": 0.0077, "ref": "Louis et al. 2025 Table 5, ACT"},
    "pact_ns": {"value": 0.9709, "error": 0.0038, "ref": "Louis et al. 2025 Table 5, P-ACT"},
    "pactlb_ns": {"value": 0.9743, "error": 0.0034,
                  "ref": "Louis et al. 2025 Table 5, P-ACT-LB (DESI DR1 BAO)"},
    "act_planck_ns": {"value": 0.9651, "error": 0.0044, "ref": "Louis et al. 2025 Table 5, Planck"},
    "pactlb_bk18_r95": {"value": 0.038, "ref": "Calabrese et al. 2025 Sect. IV.4: P-ACT-LB-BK18"},
    # SPT-3G D1: Camphuis et al. 2025/2026 (arXiv:2506.20707 v2 = PRD 113, 083504)
    "spt_ns": {"value": 0.951, "error": 0.011, "ref": "Camphuis et al. Table 1, SPT-3G"},
    "spt_desi_ns": {"value": 0.949, "error": 0.012, "ref": "Camphuis et al. Table 6, SPT-3G + DESI"},
    "spt_cmbspa_ns": {"value": 0.9679, "error": 0.0033,
                      "ref": "Camphuis et al. v2 Table 1, CMB-SPA (v1: 0.9684 +- 0.0030)"},
    "spt_cmbspa_desi_ns": {"value": 0.9726, "error": 0.0028,
                           "ref": "Camphuis et al. v2 Table 6, CMB-SPA + DESI DR2 "
                                  "(v1: 0.9728 +- 0.0027)"},
    #: Slow-roll predictions quoted in the atlas spec (N = 60).
    "phi2_n60": {"ns": 0.967, "r": 0.13},
    "phi4_n60": {"ns": 0.95, "r": 0.27},
    "starobinsky_n60": {"ns": 0.967, "r": 0.003},
    #: Future: Simons Observatory (2019 goal; 2026 forecast to 2035) and LiteBIRD.
    "so_sigma_r_2019": {"value": 0.003, "ref": "SO 2019 (arXiv:1808.07445)"},
    "so_sigma_r_2026": {"value": 0.0012, "ref": "SO 2026 (arXiv:2512.15833), survey to 2035"},
    "litebird_delta_r": {"value": 0.001, "ref": "LiteBIRD (arXiv:2202.02773) Sect. 3.1"},
}
