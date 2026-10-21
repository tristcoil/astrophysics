"""Project configuration: paths, data sources, physics conventions, the named
universes of the gallery, and the published values our results are graded against.

Published numbers live here only for comparison at the end, never as inputs to
a fit -- with two deliberate exceptions, both labelled where they are used:
the compressed CMB prior (a published summary of the Planck data, which is what
it is for; we also rebuild it from the underlying chains and from project D3),
and the redshifts of the most distant objects known (a redshift *is* a
measurement).
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
ROOT_DIR = PROJECT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
DERIVED_DIR = DATA_DIR / "derived"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

PANTHEON_DIR = DATA_DIR / "pantheon"
BAO_DIR = DATA_DIR / "bao"
DESI_DIR = DATA_DIR / "desi_results"
PLANCK_EARLY_DIR = DATA_DIR / "planck_early_lcdm"
LITERATURE_DIR = DATA_DIR / "literature"

for _d in (DATA_DIR, FIGURE_DIR, DERIVED_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Earlier atlas projects this one builds on (their cached files)
# ---------------------------------------------------------------------------
#
# Raw files that a sibling project already downloaded are hard-linked from its
# cache after a checksum check (no second download, no extra disk). Results
# that a sibling *computed* are read from its results file; if that file is
# missing, the value quoted in its README is used instead (``*_FALLBACK``).

D1_DIR = ROOT_DIR / "05-supernova-cosmology"       # D1: Pantheon+ supernovae
D3_DIR = ROOT_DIR / "09-planck-cmb"                # D3: Planck CMB, CAMB fit
D5_DIR = ROOT_DIR / "42-evolving-dark-energy"      # D5: DESI BAO, CMB prior
D8_DIR = ROOT_DIR / "38-cosmic-chronometers"       # D8: H(z) from old galaxies
B6_DIR = ROOT_DIR / "49-oldest-stars"              # B6: globular-cluster ages
I3_DIR = ROOT_DIR / "34-quasars"                   # I3: quasar catalogues

#: Local file -> where a sibling project keeps an identical copy.
SIBLING_COPIES = {
    PANTHEON_DIR / "Pantheon+SH0ES.dat": D1_DIR / "data" / "Pantheon+SH0ES.dat",
    PANTHEON_DIR / "Pantheon+SH0ES_STAT+SYS.cov": D1_DIR / "data" / "Pantheon+SH0ES_STAT+SYS.cov",
    BAO_DIR / "desi_dr2_mean.txt": D5_DIR / "data" / "bao" / "desi_dr2_mean.txt",
    BAO_DIR / "desi_dr2_cov.txt": D5_DIR / "data" / "bao" / "desi_dr2_cov.txt",
    LITERATURE_DIR / "moresco2022_lrr.tex": D8_DIR / "data" / "sources" / "moresco2022_lrr.tex",
    LITERATURE_DIR / "J_ApJ_923_262_table2.tsv":
        I3_DIR / "data" / "literature" / "J_ApJ_923_262_table2.tsv",
}
for _i in (1, 2, 3, 4):
    _name = f"spline_planck_PR4_TTTEEE_lowE_lensing_ISW.{_i}.txt"
    SIBLING_COPIES[PLANCK_EARLY_DIR / _name] = D5_DIR / "data" / "planck_early_lcdm" / _name

D3_FIT = D3_DIR / "data" / "derived" / "tt_fit.json"
D1_RESULTS = D1_DIR / "data" / "results.json"
B6_RESULTS = B6_DIR / "data" / "results.json"
D5_RESULTS = D5_DIR / "data" / "derived" / "results.json"

# ---------------------------------------------------------------------------
# Data sources (used when no sibling copy exists)
# ---------------------------------------------------------------------------

#: Pantheon+ (Scolnic et al. 2022; Brout et al. 2022): 1,701 Type Ia supernova
#: light curves with their full statistical + systematic covariance.
PANTHEON_RAW = ("https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/"
                "Pantheon%2B_Data/4_DISTANCES_AND_COVAR/")
PANTHEON_FILES = {"Pantheon+SH0ES.dat": "Pantheon%2BSH0ES.dat",
                  "Pantheon+SH0ES_STAT+SYS.cov": "Pantheon%2BSH0ES_STAT%2BSYS.cov"}

#: DESI DR2 BAO (DESI Collaboration 2025, arXiv:2503.14738), as distributed for
#: the Cobaya sampler -- the same files DESI's own likelihood reads.
BAO_RAW = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/desi_bao_dr2/"
BAO_FILES = {"desi_dr2_mean.txt": "desi_gaussian_bao_ALL_GCcomb_mean.txt",
             "desi_dr2_cov.txt": "desi_gaussian_bao_ALL_GCcomb_cov.txt"}

#: Lemos & Lewis (2023, arXiv:2302.12911): Planck PR4 chains of the *early*
#: universe only (CMB lensing and the late-time ISW effect marginalised), from
#: which the compressed CMB prior on (theta_*, omega_b, omega_b + omega_c) is
#: built. Because late-time physics is marginalised, the prior can be applied to
#: universes with curvature, which only change the distance the light travels.
PLANCK_EARLY_RAW = "https://raw.githubusercontent.com/cmbant/PlanckEarlyLCDM/main/"
PLANCK_EARLY_ROOT = "spline_planck_PR4_TTTEEE_lowE_lensing_ISW"
PLANCK_EARLY_BURN = 0.3  # the README's own recipe (as in D5)

#: DESI's released cosmology results (6 Oct 2025): MCMC chains and best fits of
#: LCDM with free curvature ("base_omegak"). Answer key only, never an input.
DESI_RESULTS_URL = "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/"
DESI_CHECKSUMS = "dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum"
CMB_FULL = ("planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_"
            "planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing")
#: Chains: folder under cobaya/base_omegak/ -> our label.
DESI_OMK_CHAINS = {
    "desi-bao-all": "DESI BAO",
    CMB_FULL: "CMB (full Planck + ACT lensing)",
    f"desi-bao-all_{CMB_FULL}": "DESI BAO + CMB (full)",
}
#: Best fits (iminuit): folder under iminuit/base_omegak/ -> our label.
DESI_OMK_BESTFITS = {
    "desi-bao-all": "DESI BAO",
    "desi-bao-all_pantheonplus": "DESI BAO + Pantheon+",
    f"desi-bao-all_{CMB_FULL}": "DESI BAO + CMB (full)",
}
DESI_CHAIN_BURN = 0.3  # as D5 used for the same release

#: Moresco et al. (2022, Living Reviews in Relativity 25, 6): the 32 cosmic-
#: chronometer H(z) measurements, from the review's LaTeX (project D8's parse).
MORESCO_LRR_EPRINT = "https://arxiv.org/e-print/2201.07241"

#: Yang et al. (2021, ApJ 923, 262), table 2: 37 quasars at z > 6.3 including
#: J0313-1806 at z = 7.642, the most distant quasar known (Wang et al. 2021).
YANG21_VIZIER = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/923/262/table2"
                 "&-out.max=unlimited&-out.all")

USER_AGENT = "astrophysics-atlas-O1/1.0 (educational; python-urllib)"

# ---------------------------------------------------------------------------
# Physics conventions (identical to DESI's CAMB settings, as in D5)
# ---------------------------------------------------------------------------

C_KM_S = 299_792.458          # speed of light, km/s (exact)
T_CMB = 2.7255                # CMB temperature today, K (Fixsen 2009)
N_EFF = 3.044                 # effective number of neutrino species
M_NU_EV = 0.06                # one massive neutrino of 0.06 eV (oscillation minimum)
NU_MASS_TO_OMEGA = 93.14      # DESI/Cobaya's omega_nu = m_nu / 93.14 eV mapping
BBN_PREDICTOR = "PArthENoPE_880.2_standard.dat"

#: Without an absolute calibration (supernovae alone, BAO alone) the expansion
#: rate H0 and omega_b only set the tiny radiation density; DESI fixes them at
#: these values, and so do we.
UNCALIBRATED_H0 = 67.36
UNCALIBRATED_OMBH2 = 0.02237

#: Seconds per Julian gigayear and kilometres per megaparsec: 1/H0 in Gyr is
#: KM_PER_MPC / SEC_PER_GYR / H0 = 977.79 / H0[km/s/Mpc].
SEC_PER_GYR = 1e9 * 365.25 * 86400.0
KM_PER_MPC = 3.0856775814913673e19
HUBBLE_TIME_GYR_X_H0 = KM_PER_MPC / SEC_PER_GYR  # 977.79 Gyr km/s/Mpc

# ---------------------------------------------------------------------------
# Priors and grids
# ---------------------------------------------------------------------------

#: DESI's priors for LCDM + curvature (from the released chain configuration):
#: flat in Omega_m, Omega_k and h*r_d, plus H0 within [20, 100] km/s/Mpc wherever
#: H0 is a parameter. Used for every comparison with DESI's chains.
PRIOR_OMK = (-0.3, 0.3)
PRIOR_OM = (0.01, 0.99)
PRIOR_HRD = (10.0, 1000.0)
PRIOR_H0 = (20.0, 100.0)

#: The wide (Omega_m, Omega_Lambda) plane of the "map of universes".
MAP_OM = (0.0, 3.0, 601)
MAP_OL = (-1.0, 3.0, 801)

#: The plane on which the supernova contours are computed (D1's grid).
SN_OM = (0.0, 1.2, 241)
SN_OL = (-0.4, 1.6, 401)

#: Fine grid around where BAO, CMB and their combinations live.
ZOOM_OM = (0.22, 0.40, 361)
ZOOM_OL = (0.56, 0.84, 561)

#: The CMB band: a sheared grid that follows the band's ridge, (Omega_m, offset).
CMB_BAND_OM = (0.145, 1.20, 423)
CMB_BAND_HALFWIDTH = 0.08
CMB_BAND_N = 321

#: Radiation density used for the "toy" universes of the gallery and the map:
#: photons + three neutrino species, all relativistic, at H0 = 67.4. The data
#: fits use the full model (one massive neutrino) instead.
H0_GALLERY = 67.4

# ---------------------------------------------------------------------------
# The most distant objects seen (redshifts are measurements; see README)
# ---------------------------------------------------------------------------

#: Each one sets a lower limit on how much smaller the universe once was:
#: a universe that "bounced" at a_min shows nothing beyond z = 1/a_min - 1.
DISTANT_OBJECTS = {
    "quasar": {"label": "quasar J0313-1806", "z": 7.6423,
               "ref": "Wang et al. 2021; [C II] redshift from Yang et al. 2021 (I3's table)"},
    "galaxy": {"label": "galaxy MoM-z14", "z": 14.44,
               "ref": "Naidu et al. 2025, OJAp 9 (2026) (JWST NIRSpec)"},
    "cmb": {"label": "the CMB", "z": 1089.9, "ref": "Planck 2018 VI (z*)"},
}
#: The quasar record moved in July 2026 (not in I3's catalogues): a Euclid
#: quasar at z ~ 7.77 from rest-UV lines (D. Yang, Hennawi et al. 2026, A&A 711,
#: A104), and EUCL J1253+7054 at z = 7.6980 from [C II] (Belladitta et al. 2026).
QUASAR_RECORD_2026 = {"label": "EUCL J172902.75+641018.1", "z": 7.77,
                      "ref": "Yang, Hennawi et al. 2026, A&A 711, A104"}

# ---------------------------------------------------------------------------
# Published values (comparison only)
# ---------------------------------------------------------------------------

#: Planck 2018 VI (A&A 641, A6). Values checked against the paper (Tables 1-2,
#: Sect. 7.3). Without BAO = TT,TE,EE+lowE+lensing; "_bao" = with BAO added.
PLANCK18 = {
    "H0": (67.36, 0.54), "Om": (0.3153, 0.0073), "ombh2": 0.02237, "omch2": 0.1200,
    "age_gyr": (13.797, 0.023), "z_eq": (3402, 26), "z_star": (1089.92, 0.25),
    "H0_bao": (67.66, 0.42), "Om_bao": (0.3111, 0.0056), "age_gyr_bao": (13.787, 0.020),
    "z_eq_bao": (3387, 21),
    "Ok_TTTEEE_lowE": (-0.044, 0.018, 0.015),        # (value, +err, -err), eq. 46b
    "Ok_lensing": (-0.0106, 0.0065, 0.0065),          # eq. 47a
    "Ok_lensing_BAO": (0.0007, 0.0019, 0.0019),       # eq. 47b
}

#: DESI DR2 Results II (arXiv:2503.14738, PRD 112, 083515), Table V, LCDM + Omega_k.
#: Their "CMB" is the full Planck PR4 CamSpec + 2018 low-l + PR4/ACT DR6 lensing.
DESI_DR2_OMK = {"DESI BAO": (0.025, 0.041), "DESI BAO + CMB (full)": (0.0023, 0.0011),
                "CMB (full Planck + ACT lensing)": (-0.0107, 0.0064, 0.0053)}
DESI_DR2_CMB_LCDM_OMK = {"Om": (0.3034, 0.0037), "H0": (68.50, 0.33)}

#: The 1990s age crisis: Bolte & Hogan (1995, Nature 376, 399): globular
#: clusters 15.8 +- 2.1 Gyr; Chaboyer et al. (1996, Science 271, 957): median
#: 14.56 Gyr, 95% lower bound 12.07 Gyr; Freedman et al. (1994, Nature 371,
#: 757): H0 = 80 +- 17 km/s/Mpc from Cepheids in M100.
AGE_CRISIS_1990S = {"gc_age_bolte_hogan95": (15.8, 2.1), "gc_age_chaboyer96": 14.56,
                    "gc_age_chaboyer96_lower95": 12.07, "H0_freedman94": (80.0, 17.0)}

#: Petrosian, Salpeter & Szekeres (1967, ApJ 147, 1222): closed Lemaitre models
#: loitering near z ~ 2 to explain the pile-up of quasar redshifts near 1.95
#: (rho_0 = 1e-30 g/cm^3 at H0 = 75, i.e. Omega_m ~ 0.1).
PSS67 = {"z_loiter": 2.0, "Om": 0.095}

#: Project D1 (05-supernova-cosmology), as published in its README: Pantheon+
#: alone with curvature free.
D1_FALLBACK = {"Om": (0.295, 0.054), "OL": (0.609, 0.081), "H0": (73.55, 1.02),
               "H0_kinematic": (73.53, 1.05), "age_gyr": (12.46, 0.24)}
#: Project D3 (09-planck-cmb): CAMB fit to the Planck TT spectrum.
D3_FALLBACK = {"H0": (66.66, 0.90), "Om": (0.3247, 0.0132), "age_gyr": (13.825, 0.033)}
#: Project B6 (49-oldest-stars): the oldest globular clusters plus the time it
#: took them to form need a universe at least this old (median, std, 5th pct).
B6_FALLBACK = {"t_min_gyr": 13.605, "t_min_err": 0.835, "t_min_p5": 12.24}

#: Analytic ages in units of 1/H0, for the textbook universes.
ANALYTIC_AGES = {"Einstein-de Sitter": 2.0 / 3.0, "Milne (empty)": 1.0,
                 "radiation only": 0.5}
