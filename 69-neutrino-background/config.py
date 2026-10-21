"""Project configuration: paths, data sources, physics conventions, fit
settings, and the published values we grade our results against.

Published numbers live at the bottom and are used only for comparison, with
two deliberate exceptions, both priors standing in for likelihoods we do not
run: Planck's low-ℓ polarisation measurement of the optical depth τ (a
Gaussian prior), and the calibration prior of Planck's high-ℓ likelihood.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
DERIVED_DIR = DATA_DIR / "derived"
THEORY_CACHE = DATA_DIR / "cache" / "theory"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, DERIVED_DIR, THEORY_CACHE, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

USER_AGENT = "astrophysics-atlas-P3/1.0 (educational; python-urllib)"

#: Worker processes for CAMB (one OpenMP thread each). CAMB takes ~10 s per
#: universe on one core at the accuracy we use.
N_WORKERS = 6

# ---------------------------------------------------------------------------
# Earlier atlas projects whose downloads we reuse (hard-linked, not refetched)
# ---------------------------------------------------------------------------

D3_DIR = ROOT_DIR / "09-planck-cmb"
D5_DIR = ROOT_DIR / "42-evolving-dark-energy"
#: D3's fit of LCDM to Planck's TT spectrum (our starting point and our check).
D3_FIT = D3_DIR / "data" / "derived" / "tt_fit.json"
#: Project P2 (the thermal history) solved neutrino decoupling with the same
#: method in independent code; its results are read for a cross-check if present.
P2_RESULTS = ROOT_DIR / "68-thermal-history" / "data" / "results.json"

# ---------------------------------------------------------------------------
# Physics conventions
# ---------------------------------------------------------------------------

#: CMB temperature today (Fixsen 2009, ApJ 707, 916: 2.72548 ± 0.00057 K).
T_CMB_K = 2.7255
#: The Standard Model's effective number of neutrino species (Bennett et al.
#: 2021). Used as CAMB's N_eff wherever N_eff is not varied, as DESI does.
NEFF_SM = 3.044
#: BBN table that sets the helium fraction from ω_b and N_eff (as DESI's runs).
BBN_PREDICTOR = "PArthENoPE_880.2_standard.dat"
#: Effective redshifts of the DESI DR2 BAO measurements (checked against the
#: data file when it is loaded).
DESI_Z_UNIQUE = (0.295, 0.510, 0.706, 0.934, 1.321, 1.484, 2.330)
SPEED_OF_LIGHT_KMS = 299_792.458

# ---------------------------------------------------------------------------
# Fit settings
# ---------------------------------------------------------------------------

#: The six LCDM parameters every fit varies.
LCDM = ["ombh2", "omch2", "theta100", "logA", "ns", "tau"]

#: Finite-difference half-steps for the Jacobians, about a quarter of each
#: parameter's Planck TT,TE,EE uncertainty: large enough to rise above CAMB's
#: numerical noise, small enough for the model to be linear across the step.
STEPS = {
    "ombh2": 0.00005, "omch2": 0.0005, "theta100": 0.0001, "logA": 0.004, "ns": 0.0015,
    "tau": 0.002, "A_planck": 0.001, "nnu": 0.05, "mnu": 0.02, "w0": 0.03, "wa": 0.1,
}
#: Hard bounds (the priors' edges; Σm_ν ≥ 0 is physical).
PARAM_BOUNDS = {
    "ombh2": (0.005, 0.1), "omch2": (0.001, 0.99), "theta100": (0.5, 10.0),
    "logA": (1.61, 3.91), "ns": (0.8, 1.2), "tau": (0.01, 0.8), "A_planck": (0.9, 1.1),
    "nnu": (0.05, 10.0), "mnu": (0.0, 5.0), "w0": (-3.0, 1.0), "wa": (-3.0, 2.0),
}

# ---------------------------------------------------------------------------
# Priors standing in for likelihoods we do not run
# ---------------------------------------------------------------------------

#: Optical depth from Planck's low-ℓ polarisation alone ("lowE"; Planck 2018 VI,
#: eq. 4, quoting Planck 2018 V): SimAll EE with low-ℓ TT, other parameters at a
#: fiducial ΛCDM. Our Gaussian stand-in for Planck's low-ℓ likelihoods.
TAU_PRIOR = (0.0506, 0.0086)
#: Planck's overall calibration y_P (A_planck), Gaussian prior of the 2018 high-ℓ
#: likelihoods (Planck 2018 V).
A_PLANCK_PRIOR = (1.0, 0.0025)

# ---------------------------------------------------------------------------
# Published cosmological constraints (comparison only; all checked against the
# papers and Planck's PLA tables on 2026-10-19)
# ---------------------------------------------------------------------------

#: N_eff, mean and 68% error. Planck 2018 VI eqs. 66-67 and the PLA baseline
#: tables v2 (base_nnu chains; "+lensing" by importance sampling).
NEFF_PUBLISHED = {
    "Planck TT+lowE": (3.00, 0.28),
    "Planck TT,TE,EE+lowE": (2.92, 0.19),
    "Planck TT,TE,EE+lowE+lensing": (2.89, 0.19),
    "Planck TT,TE,EE+lowE+lensing+BAO": (2.99, 0.17),
    # Calabrese et al. 2025 (ACT DR6 extended models), eq. 30. "LB" = ACT DR6 +
    # Planck PR4 lensing + DESI DR1 BAO.
    "ACT DR6 P-ACT": (2.73, 0.14),
    "ACT DR6 P-ACT-LB": (2.86, 0.13),
    # Camphuis et al. 2025 (SPT-3G D1), eq. 67: Planck + SPT + ACT (CMB-SPA).
    "SPT-3G CMB-SPA": (2.82, 0.12),
    # Elbers et al. 2025 (DESI DR2 BAO + CMB, CamSpec), Table II.
    "DESI DR2 + CMB": (3.23, 0.18),
}

#: Σm_ν 95% upper limits (eV), three degenerate states, flat prior Σm_ν ≥ 0.
MNU_PUBLISHED = {
    "Planck TT+lowE": 0.54,                       # Planck 2018 VI eq. 58a
    "Planck TT+lowE+BAO": 0.16,                   # eq. 62a
    "Planck TT,TE,EE+lowE": 0.26,                 # eq. 58b
    "Planck TT,TE,EE+lowE+lensing": 0.24,         # eq. 59b
    "Planck TT,TE,EE+lowE+lensing+BAO": 0.12,     # eq. 63b
    # DESI 2024 VI printed 0.072; DESI 2024 VII (eq. 4.3) corrected it to 0.082
    # after a bug fix in ACT's lensing likelihood (v1.1 -> v1.2).
    "DESI DR1 + CMB": 0.082,
    "DESI DR1 + CMB (as first printed)": 0.072,
    "ACT DR6 P-ACT-LB": 0.089,                    # Calabrese et al. 2025, eq. 33
    "DESI DR2 + CMB (CamSpec)": 0.0642,           # DESI DR2 Results II, eq. 29
    "DESI DR2 + CMB (plik)": 0.0691,              # eq. 30
    "DESI DR2 + CMB (LoLLiPoP+HiLLiPoP)": 0.0774,  # eq. 31
    "DESI DR2 + CMB, w0wa": 0.163,                # eq. 33
    "DESI DR2 + CMB + DESY5, w0wa": 0.129,        # eq. 34
    "DESI DR2 + CMB + Pantheon+, w0wa": 0.117,
    "DESI DR2 + CMB + Union3, w0wa": 0.139,
    "SPT-3G CMB-SPA + DESI DR2": 0.055,           # Camphuis et al. 2025, eq. 96
}
#: Elbers et al. 2025 (DESI DR2 + CMB): the "effective" mass allowed to go
#: negative (eq. 31; 68%), and their profile-likelihood parabola (Table III).
MNU_EFF_ELBERS = (-0.101, 0.047, 0.056)
MNU_PROFILE_ELBERS = {"mu": -0.036, "sigma": 0.043, "feldman_cousins_95": 0.053}

# ---------------------------------------------------------------------------
# Neutrino decoupling: our ladder of approximations (label, rates, QED order)
# ---------------------------------------------------------------------------

DECOUPLING_RUNS = [
    ("Instantaneous decoupling", "off", 0),
    ("+ QED plasma, O(e²)", "off", 2),
    ("+ QED plasma, O(e³)", "off", 3),
    ("Weak rates, Maxwell–Boltzmann", "MB", 0),
    ("Weak rates, Fermi–Dirac", "FD", 0),
    ("Fermi–Dirac + electron mass", "FD+me", 0),
    ("… + QED O(e²)", "FD+me", 2),
    ("… + QED O(e³): our N_eff", "FD+me", 3),
]
DECOUPLING_HEADLINE = DECOUPLING_RUNS[-1][0]
#: Bars of the "one ingredient at a time" panel: (label, run label).
DECOUPLING_WATERFALL = [
    ("instantaneous", DECOUPLING_RUNS[0][0]),
    ("+ QED, instantaneous", DECOUPLING_RUNS[2][0]),
    ("weak rates (MB)", DECOUPLING_RUNS[3][0]),
    ("Fermi–Dirac", DECOUPLING_RUNS[4][0]),
    ("+ electron mass", DECOUPLING_RUNS[5][0]),
    ("+ QED O(e²)", DECOUPLING_RUNS[6][0]),
    ("+ QED O(e³)", DECOUPLING_RUNS[7][0]),
]
#: Escudero (2020, JCAP 05, 048) Table 1, separate ν_e and ν_μ,τ temperatures,
#: the same method we follow, keyed by our run labels; and his eqs. A.12-A.13.
ESCUDERO_2020_TABLE1 = {
    DECOUPLING_RUNS[0][0]: 3.000, DECOUPLING_RUNS[1][0]: 3.011, DECOUPLING_RUNS[2][0]: 3.010,
    DECOUPLING_RUNS[3][0]: 3.042, DECOUPLING_RUNS[4][0]: 3.038, DECOUPLING_RUNS[5][0]: 3.035,
    DECOUPLING_RUNS[6][0]: 3.045, DECOUPLING_RUNS[7][0]: 3.044,
}
ESCUDERO_2020_A12 = {"z_gamma": 1.39791, "z_nue": 1.00237, "z_numu": 1.00098, "neff": 3.0443}
#: Full (momentum-resolved) calculations of the Standard Model value.
NEFF_SM_LITERATURE = {
    "Bennett et al. 2021": (3.0440, 0.0002),        # JCAP 04, 073, eq. 5.1
    "Froustey, Pitrou & Volpe 2020": (3.0440, 0.0001),
    "Escudero, Jackson, Laine & Sandner 2026": (3.0443, 0.0),
    "Cielo, Escudero, Mangano & Pisanti 2023": (3.0432, 0.0002),
    "Mangano et al. 2005": (3.046, 0.002),
}

# ---------------------------------------------------------------------------
# Neutrino oscillations: NuFIT 6.0 (Esteban et al. 2024, JHEP 12, 216), Table 1,
# "IC24 with SK-atm" -- best fit, minus and plus 1-sigma. dm in eV^2.
# ---------------------------------------------------------------------------

NUFIT = {
    "NO": {"s12sq": (0.308, 0.011, 0.012), "s13sq": (0.02215, 0.00058, 0.00056),
           "dm21": (7.49e-5, 0.19e-5, 0.19e-5), "dm3l": (2.513e-3, 0.019e-3, 0.021e-3)},
    "IO": {"s12sq": (0.308, 0.011, 0.012), "s13sq": (0.02231, 0.00056, 0.00056),
           "dm21": (7.49e-5, 0.19e-5, 0.19e-5), "dm3l": (-2.484e-3, 0.020e-3, 0.020e-3)},
}
#: Δχ² of the inverted relative to the normal ordering (NuFIT 6.0): 6.1 with
#: Super-Kamiokande's atmospheric data, −0.6 without (IO then fits slightly better).
NUFIT_DCHI2_IO = {"TByes": 6.1, "TBoff": -0.6}
#: Minimum Σm_ν (lightest mass zero) from the best fits above (NuFIT's own
#: eq. 5.8 quotes 0.058 and 0.098 eV at 95% CL).
MIN_SUM_NO = (7.49e-5) ** 0.5 + (2.513e-3) ** 0.5
MIN_SUM_IO = (2.484e-3) ** 0.5 + (2.484e-3 - 7.49e-5) ** 0.5

# ---------------------------------------------------------------------------
# Laboratory and direct-detection numbers (comparison only)
# ---------------------------------------------------------------------------

#: KATRIN (Aker et al. 2025, Science 388, 180): m_β < 0.45 eV at 90% CL
#: (Lokhov-Tkachov; equals the sensitivity because the best fit m² = -0.14 eV²
#: is negative); Feldman-Cousins would give 0.31 eV. Final aim: < 0.3 eV.
KATRIN_MBETA_LIMIT = 0.45
KATRIN_MBETA_FC = 0.31
KATRIN_FINAL_SENSITIVITY = 0.3
#: KATRIN's limit on the local relic overdensity (Aker et al. 2022, PRL 129,
#: 011806): η < 9.7e10 at 90% CL, relative to 56 per cm³ per state, assuming the
#: Majorana capture rate.
KATRIN_ETA_LIMIT = 9.7e10
#: Capture of non-relativistic ν_e on tritium: σv = 3.834e-45 cm² × c (Long,
#: Lunardini & Sabancilar 2014, eq. 2.29); 4.06 (Dirac) and 8.12 (Majorana)
#: captures per year for 100 g of tritium (their eq. 3.1).
TRITIUM_SIGMA_V_CM2 = 3.834e-45
PTOLEMY_PUBLISHED = {"dirac": 4.06, "majorana": 8.12}
#: Local overdensity n/n̄ of relic neutrinos from the Milky Way's gravity, for
#: one mass state: Mertsch et al. 2020 (abstract) and Zimmer, Correa & Ando
#: 2023 (median over simulated Milky Way analogues, read from their figure).
CLUSTERING = {"Mertsch et al. 2020": {0.01: 1.0053, 0.05: 1.12, 0.10: 1.50, 0.30: 6.0},
              "Zimmer et al. 2023": {0.01: 1.0011, 0.05: 1.034, 0.10: 1.14, 0.30: 1.95}}

#: Reference total matter density for scale estimates (Planck 2018 TT,TE,EE+
#: lowE+lensing+BAO, Table 2).
OMEGA_M_REF = 0.3111

#: N_eff from Big Bang nucleosynthesis alone (helium + deuterium, no CMB):
#: Yeh, Shelton, Olive & Fields 2022 (JCAP 10, 046), Table 2, N_ν = 2.889 ± 0.229,
#: converted with their eq. 3.2, N_eff = 1.0147 N_ν.
NEFF_BBN_YEH2022 = (2.889 * 1.0147, 0.229 * 1.0147)
#: The local expansion rate from Cepheid-calibrated supernovae (SH0ES; Riess et
#: al. 2022, ApJL 934, L7), for the N_eff-H0 panel.
H0_SH0ES = (73.04, 1.04)
#: Escape speed from the Milky Way at the Sun's position (Deason et al. 2019,
#: MNRAS 485, 3514: 528 +24/-25 km/s), for the speed panel of figure 01.
V_ESC_MW_KMS = 528.0
