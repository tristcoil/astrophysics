"""Project configuration: paths, data sources, physics conventions, priors, and
the published values we grade our own results against.

Published numbers live here only for comparison at the end. None of them is an
input to a fit, with one deliberate exception: the compressed CMB prior *is*
a published summary of the Planck data (that is what it is for), and we also
rebuild it ourselves from the underlying chains and from project D3.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
DERIVED_DIR = DATA_DIR / "derived"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, FIGURE_DIR, DERIVED_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: DESI DR2 BAO (DESI Collaboration 2025, arXiv:2503.14738 and 2503.14739), as
#: distributed for the Cobaya sampler -- the same files DESI's own likelihood reads.
BAO_REPO_RAW = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/"
BAO_FILES = {
    "desi_dr2_mean.txt": "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    "desi_dr2_cov.txt": "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
    # DESI DR1 (2024), for the DR1 -> DR2 comparison and our own BAO measurement.
    "desi_dr1_mean.txt": "desi_2024_gaussian_bao_ALL_GCcomb_mean.txt",
    "desi_dr1_cov.txt": "desi_2024_gaussian_bao_ALL_GCcomb_cov.txt",
    # SDSS (BOSS DR12 + eBOSS DR16) consensus BAO, for the Hubble-diagram figure.
    "sdss_dr12_lrg_mean.dat": "sdss_DR12_LRG_BAO_DMDH.dat",
    "sdss_dr12_lrg_cov.txt": "sdss_DR12_LRG_BAO_DMDH_covtot.txt",
    "sdss_dr16_lrg_mean.dat": "sdss_DR16_LRG_BAO_DMDH.dat",
    "sdss_dr16_lrg_cov.txt": "sdss_DR16_LRG_BAO_DMDH_covtot.txt",
    "sdss_dr16_qso_mean.txt": "sdss_DR16_QSO_BAO_DMDH.txt",
    "sdss_dr16_qso_cov.txt": "sdss_DR16_QSO_BAO_DMDH_covtot.txt",
    "sdss_mgs_prob.txt": "sdss_MGS_prob.txt",
}

#: Type Ia supernova samples, as distributed for Cobaya (the versions DESI used,
#: plus the 2025 DES recalibration). Local folder -> {local name: repo path}.
SN_REPO_RAW = "https://raw.githubusercontent.com/CobayaSampler/sn_data/master/"
SN_FILES = {
    "PantheonPlus": {
        "Pantheon+SH0ES.dat": "PantheonPlus/Pantheon+SH0ES.dat",
        "Pantheon+SH0ES_STAT+SYS.cov": "PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov",
    },
    "Union3": {
        "lcparam_full.txt": "Union3/lcparam_full.txt",
        "mag_covmat.txt": "Union3/mag_covmat.txt",
    },
    "DESY5": {
        "DES-SN5YR_HD.csv": "DESY5/DES-SN5YR_HD.csv",
        "covsys_000.txt": "DESY5/covsys_000.txt",
    },
    "DES-Dovekie": {
        "DES-Dovekie_HD.csv": "DES-Dovekie/DES-Dovekie_HD.csv",
        "covtot_inv_000.npz": "DES-Dovekie/covtot_inv_000.npz",
    },
}

#: Lemos & Lewis (2023, arXiv:2302.12911): Planck PR4 chains of the *early*
#: universe only (lensing and late-time ISW marginalised). DESI's compressed CMB
#: prior is the Gaussian summary of these chains; we recompute it from them.
PLANCK_EARLY_RAW = "https://raw.githubusercontent.com/cmbant/PlanckEarlyLCDM/main/"
PLANCK_EARLY_ROOT = "spline_planck_PR4_TTTEEE_lowE_lensing_ISW"
PLANCK_EARLY_FILES = [f"{PLANCK_EARLY_ROOT}.{i}.txt" for i in (1, 2, 3, 4)] + [
    f"{PLANCK_EARLY_ROOT}.margestats"]

#: DESI's own published cosmology results (released 6 Oct 2025): MCMC chains and
#: best fits (iminuit). Used only as an answer key for our independent likelihood.
DESI_RESULTS_URL = "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/"
DESI_CHECKSUMS = "dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum"

_CMB_FULL = ("planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_"
             "planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing")
_CMBC = "CMB-compressed-theta-ombh2-ombch2"

#: DESI dataset folder names -> our short labels.
DESI_COMBOS = {
    "desi-bao-all": "DESI",
    "desi-bao-all_pantheonplus": "DESI+Pantheon+",
    "desi-bao-all_union3": "DESI+Union3",
    "desi-bao-all_desy5sn": "DESI+DESY5",
    f"desi-bao-all_{_CMBC}": "DESI+CMBc",
    f"desi-bao-all_pantheonplus_{_CMBC}": "DESI+CMBc+Pantheon+",
    f"desi-bao-all_union3_{_CMBC}": "DESI+CMBc+Union3",
    f"desi-bao-all_desy5sn_{_CMBC}": "DESI+CMBc+DESY5",
    f"desi-bao-all_{_CMB_FULL}": "DESI+CMB(full)",
    f"desi-bao-all_pantheonplus_{_CMB_FULL}": "DESI+CMB(full)+Pantheon+",
    f"desi-bao-all_union3_{_CMB_FULL}": "DESI+CMB(full)+Union3",
    f"desi-bao-all_desy5sn_{_CMB_FULL}": "DESI+CMB(full)+DESY5",
}

#: Best-fit (MAP) files we fetch: every minimiser run for these models/combos.
DESI_BESTFIT_MODELS = ("base", "base_w_wa")

#: Chains we fetch (about 410 MB in total), for posterior-level comparison.
DESI_CHAINS = [
    "base/desi-bao-all",
    "base_w_wa/desi-bao-all",
    "base_w_wa/desi-bao-all_pantheonplus",
    "base_w_wa/desi-bao-all_union3",
    "base_w_wa/desi-bao-all_desy5sn",
    f"base_w_wa/desi-bao-all_{_CMBC}",
    f"base_w_wa/desi-bao-all_{_CMBC}_desy5sn",
    # supernovae alone (about 0.1 GB more)
    "base/pantheonplus", "base/union3", "base/desy5sn",
    "base_w_wa/pantheonplus", "base_w_wa/union3", "base_w_wa/desy5sn",
]

#: Flat-LCDM Omega_m from each supernova sample alone, as its own team published
#: it: Brout et al. 2022 (Pantheon+), Rubin et al. 2023 (Union3/UNITY1.5),
#: DES Collaboration 2024 (DES-SN5YR), Popovic et al. 2025 (DES-Dovekie).
SN_OM_PUBLISHED = {"PantheonPlus": (0.334, 0.018, 0.018), "Union3": (0.356, 0.028, 0.026),
                   "DESY5": (0.352, 0.017, 0.017), "DES-Dovekie": (0.330, 0.015, 0.015)}

#: DESI DR1 large-scale-structure catalogues (DESI 2024 II; v1.5), for our own
#: BAO measurement: luminous red galaxies (LRGs) with weights, plus random
#: points that trace the survey's footprint and depth. Each random file holds
#: ~6.4 randoms per galaxy; two per sky cap give ~13x. About 3.2 GB in total.
LSS_URL = "https://data.desi.lbl.gov/public/dr1/survey/catalogs/dr1/LSS/iron/LSScats/v1.5/"
LSS_FILES = ["LRG_NGC_clustering.dat.fits", "LRG_SGC_clustering.dat.fits"] + [
    f"LRG_{cap}_{i}_clustering.ran.fits" for cap in ("NGC", "SGC") for i in (0, 1)]

#: Project D3's own fit of LCDM to the Planck TT spectrum (for our home-made
#: CMB prior). Optional: the analysis skips that variant if it is missing.
D3_FIT = ROOT_DIR / "09-planck-cmb" / "data" / "derived" / "tt_fit.json"
#: Project D1's flat-LCDM Omega_m from Pantheon+ (for the comparison in step 1).
D1_RESULTS = ROOT_DIR / "05-supernova-cosmology" / "data" / "results.json"

# ---------------------------------------------------------------------------
# Physics conventions -- identical to DESI's (CAMB) settings
# ---------------------------------------------------------------------------

#: Speed of light, km/s (exact).
C_KM_S = 299_792.458
#: CMB temperature today, K (Fixsen 2009; CAMB default).
T_CMB = 2.7255
#: Effective number of neutrino species (standard model value).
N_EFF = 3.044
#: Sum of neutrino masses, eV: one massive eigenstate of 0.06 eV (the minimum
#: allowed by oscillations), as in DESI's baseline.
M_NU_EV = 0.06
#: Cobaya/DESI convert neutrino mass to density with this factor when mapping
#: Omega_m -> omega_c (omch2 = omm h^2 - mnu/93.14 - ombh2). CAMB's own density
#: of the massive neutrino is computed exactly from the Fermi-Dirac integral.
NU_MASS_TO_OMEGA = 93.14
#: BBN table that sets the helium fraction from omega_b (as in DESI's runs).
BBN_PREDICTOR = "PArthENoPE_880.2_standard.dat"

#: For BAO without an absolute calibration DESI fixes H0 and omega_b at these
#: values. They only set the tiny radiation density (and r_d is then replaced by
#: the free parameter h*r_d), so the choice is immaterial; we copy it exactly.
UNCALIBRATED_H0 = 67.36
UNCALIBRATED_OMBH2 = 0.02237

# ---------------------------------------------------------------------------
# Priors (DESI DR2 paper, sec. VII; Table 2 of DESI 2024 VI)
# ---------------------------------------------------------------------------

PRIORS = {
    "H0": (20.0, 100.0),
    "ombh2": (0.005, 0.1),
    "Om": (0.01, 0.99),
    "w0": (-3.0, 1.0),
    "wa": (-3.0, 2.0),
    "hrd": (10.0, 1000.0),  # h * r_d in Mpc, for uncalibrated BAO
}
#: Dark energy must not dominate at early times: w0 + wa < 0.
REQUIRE_W0_PLUS_WA_NEGATIVE = True

# ---------------------------------------------------------------------------
# The compressed CMB prior on (theta_*, omega_b, omega_b + omega_c)
# ---------------------------------------------------------------------------

#: Exactly as in DESI's released likelihood configuration
#: (cobaya/base_w_wa/desi-bao-all_CMB-compressed-theta-ombh2-ombch2/chain.updated.yaml).
#: The paper prints the same numbers rounded (its eqs. A1-A2).
#: NOTE: the third mean, 0.14207901, is the chain mean of omega_m h^2, which
#: *includes* the massive neutrino (0.00064); the chain mean of omega_b + omega_c
#: is 0.14143. DESI's likelihood compares it with omega_b + omega_c. See README.
CMB_PRIOR_DESI = {
    "mean": (0.01041027, 0.02223208, 0.14207901),
    "cov": ((6.6209942e-12, 1.24442058e-10, -1.19287532e-09),
            (1.24442058e-10, 2.13441666e-08, -9.40008323e-08),
            (-1.19287532e-09, -9.40008323e-08, 1.48841714e-06)),
}
#: Burn-in fraction used when we recompute the prior from the chains
#: (the Lemos & Lewis README's own recipe).
PLANCK_EARLY_BURN = 0.3

# ---------------------------------------------------------------------------
# Published values (for comparison only)
# ---------------------------------------------------------------------------

#: DESI DR2 Results II (arXiv:2503.14738), Table V (Delta chi^2_MAP of w0waCDM
#: relative to LCDM) and the significance quoted for 2 extra parameters.
DESI_DCHI2 = {
    "DESI": (-4.7, 1.7),
    "DESI+CMBc": (-8.0, 2.4),
    "DESI+CMB(no lensing)": (-9.7, 2.7),
    "DESI+CMB(full)": (-12.5, 3.1),
    "DESI+Pantheon+": (-4.9, 1.7),
    "DESI+Union3": (-10.1, 2.7),
    "DESI+DESY5": (-13.6, 3.3),
    "DESI+CMB(full)+Pantheon+": (-10.7, 2.8),
    "DESI+CMB(full)+Union3": (-17.4, 3.8),
    "DESI+CMB(full)+DESY5": (-21.0, 4.2),
}

#: DESI DR2 Table IV (marginalised means and 68% intervals; asymmetric as
#: (value, +err, -err); one-sided limits are noted in comments).
DESI_PARAMS = {
    ("lcdm", "DESI"): {"Om": (0.2975, 0.0086, 0.0086), "hrd": (101.54, 0.73, 0.73)},
    ("w0wa", "DESI"): {"Om": (0.352, 0.041, 0.018), "w0": (-0.48, 0.35, 0.17)},  # wa < -1.34
    ("w0wa", "DESI+Pantheon+"): {"Om": (0.298, 0.025, 0.011), "w0": (-0.888, 0.055, 0.064),
                                  "wa": (-0.17, 0.46, 0.46)},
    ("w0wa", "DESI+Union3"): {"Om": (0.328, 0.019, 0.014), "w0": (-0.70, 0.11, 0.11),
                               "wa": (-0.99, 0.57, 0.57)},
    ("w0wa", "DESI+DESY5"): {"Om": (0.319, 0.017, 0.011), "w0": (-0.781, 0.067, 0.076),
                              "wa": (-0.72, 0.47, 0.47)},
    ("w0wa", "DESI+CMBc"): {"Om": (0.353, 0.022, 0.022), "H0": (63.7, 1.7, 2.2),
                             "w0": (-0.43, 0.22, 0.22), "wa": (-1.72, 0.64, 0.64)},
    ("w0wa", "DESI+CMB(full)"): {"Om": (0.353, 0.021, 0.021), "H0": (63.6, 1.6, 2.1),
                                  "w0": (-0.42, 0.21, 0.21), "wa": (-1.75, 0.58, 0.58)},
    ("w0wa", "DESI+CMB(full)+Pantheon+"): {"Om": (0.3114, 0.0057, 0.0057),
                                            "H0": (67.51, 0.59, 0.59),
                                            "w0": (-0.838, 0.055, 0.055),
                                            "wa": (-0.62, 0.22, 0.19)},
    ("w0wa", "DESI+CMB(full)+Union3"): {"Om": (0.3275, 0.0086, 0.0086),
                                         "H0": (65.91, 0.84, 0.84),
                                         "w0": (-0.667, 0.088, 0.088),
                                         "wa": (-1.09, 0.31, 0.27)},
    ("w0wa", "DESI+CMB(full)+DESY5"): {"Om": (0.3191, 0.0056, 0.0056),
                                        "H0": (66.74, 0.56, 0.56),
                                        "w0": (-0.752, 0.057, 0.057),
                                        "wa": (-0.86, 0.23, 0.20)},
}

#: DESI DR2: LCDM fit to BAO alone, chi^2 / dof = 10.2 / (13 - 2).
DESI_LCDM_CHI2 = 10.2

#: Popovic et al. 2025 (arXiv:2511.07517), the DES-Dovekie recalibration of the
#: DES 5-year supernovae: with DESI DR2 + CMB (Planck+ACT+SPT), w0 = -0.803 +-
#: 0.054, wa = -0.72 +- 0.21, a 3.2 sigma preference (4.2 sigma before).
DOVEKIE_PUBLISHED = {"w0": (-0.803, 0.054), "wa": (-0.72, 0.21), "sigma": 3.2,
                     "sigma_before": 4.2, "Om_flat_lcdm": (0.330, 0.015)}

#: Planck 2018 (TT,TE,EE+lowE+lensing) LCDM, the usual reference point.
PLANCK18 = {"H0": 67.36, "Om": 0.3153, "ombh2": 0.02237, "omch2": 0.1200}
