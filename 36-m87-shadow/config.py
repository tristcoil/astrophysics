"""Paths, data sources, physical constants and the published numbers we grade against.

Every literature value used anywhere in the project lives here, with its source,
so a reader can check each comparison against the paper it came from.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths (location-independent, so the folder can be renumbered freely)
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"  # posterior samples, images, intermediate products
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

# H1 (17-galactic-centre) wrote its posterior here; used for the Sgr A* test if present.
H1_RESULTS_PATH = PROJECT_DIR.parent / "17-galactic-centre" / "data" / "results.json"


# ---------------------------------------------------------------------------
# Data sources: the EHT's own GitHub data releases, pinned to exact commits
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Source:
    key: str
    url: str
    filename: str
    description: str
    reference: str
    unpack_to: str | None = None  # sub-directory of RAW_DIR for tarballs


_GH = "https://raw.githubusercontent.com/eventhorizontelescope"

SOURCES = [
    Source(
        key="m87_2017",
        url=f"{_GH}/2019-D01-01/27ae76aacb9067a83388395c3826012b97aac371/"
        "EHTC_FirstM87Results_Apr2019_uvfits.tgz",
        filename="EHTC_FirstM87Results_Apr2019_uvfits.tgz",
        description="M87, April 5/6/10/11 2017, both bands, EHT-HOPS network-calibrated Stokes I, "
        "10-s averages (data product 2019-D01-01)",
        reference="EHT Collaboration 2019c, ApJL 875, L3 (M87 Paper III)",
        unpack_to="m87_2017",
    ),
    Source(
        key="m87_2018",
        url=f"{_GH}/2024-D01-01/fda41da417419fd9b54505ea3fa70a8cf6fb60ec/EHTC_M872018_uvfits.tgz",
        filename="EHTC_M872018_uvfits.tgz",
        description="M87, April 21/22/25 2018, four bands, HOPS and CASA pipelines, network-calibrated "
        "Stokes I, 10-s averages (data product 2024-D01-01)",
        reference="EHT Collaboration 2024, A&A 681, A79",
        unpack_to="m87_2018",
    ),
    Source(
        key="sgra_2017",
        url=f"{_GH}/2022-D02-01/52640fe116ac20fd67c7d1937c2e09a0ef49caa3/"
        "EHTC_FirstSgrAResults_May2022_uvfits.tgz",
        filename="EHTC_FirstSgrAResults_May2022_uvfits.tgz",
        description="Sgr A*, April 6/7 2017, both bands, HOPS and CASA pipelines, network- and "
        "LMT-calibrated, plus light-curve-normalised copies (data product 2022-D02-01)",
        reference="EHT Collaboration 2022b, ApJL 930, L13 (Sgr A* Paper II)",
        unpack_to="sgra_2017",
    ),
    Source(
        key="ehtim_pipeline",
        url=f"{_GH}/2019-D01-02/80d230e76c08edd31548e7cfe61dc7cf94300780/eht-imaging/eht-imaging_pipeline.py",
        filename="eht-imaging_pipeline.py",
        description="The EHT's own eht-imaging script that made the 2019 M87 image (data product 2019-D01-02)",
        reference="EHT Collaboration 2019d, ApJL 875, L4 (M87 Paper IV)",
    ),
]

# ---------------------------------------------------------------------------
# Observations
# ---------------------------------------------------------------------------
M87_2017_DAYS = {"095": "April 5", "096": "April 6", "100": "April 10", "101": "April 11"}
M87_2017_BANDS = ("lo", "hi")  # 227.1 and 229.1 GHz
M87_2018_DAYS = {"111": "April 21", "112": "April 22", "115": "April 25"}
M87_2018_BANDS = ("b1", "b2", "b3", "b4")
SGRA_DAYS = {"096": "April 6", "097": "April 7"}

STATION_NAMES = {
    # 2017 codes
    "AA": "ALMA", "AP": "APEX", "AZ": "SMT", "JC": "JCMT", "LM": "LMT", "PV": "IRAM 30 m",
    "SM": "SMA", "SP": "SPT",
    # 2018 codes (different convention in the 2018 release)
    "AX": "APEX", "GL": "GLT", "MG": "SMT", "MM": "JCMT", "SW": "SMA",
}
# Stations sharing a site: their baselines are a few hundred metres to ~2 km long and
# see the large-scale jet, not the ring. Closure triangles containing them are trivial.
INTRA_SITE_PAIRS = {frozenset(("AA", "AP")), frozenset(("JC", "SM")), frozenset(("AA", "AX")),
                    frozenset(("MM", "SW"))}
INTRA_SITE_UVDIST = 0.1e9  # lambda; same cut as the EHT pipelines (uv_zblcut)

# A-priori amplitude calibration uncertainty per station: the SEFD error budget of
# M87 Paper III, as used in the EHT's own eht-imaging script. A visibility amplitude
# scales as sqrt(SEFD_i SEFD_j), so a station's gain error is about half its SEFD error.
SEFD_ERROR_BUDGET = {"AA": 0.10, "AP": 0.11, "AZ": 0.07, "LM": 0.22, "PV": 0.10, "SM": 0.15, "JC": 0.14,
                     "SP": 0.07}
# The LMT's 2017 gains were far less stable than its budget (pointing and focus problems;
# Paper IV self-calibrates it separately), so it is given a loose prior instead.
LMT_GAIN_SIGMA = 0.5

# Data conditioning, following M87 Paper VI section 2.2
SYSTEMATIC_FRACTION = 0.01  # non-closing error added in quadrature to each amplitude
MIN_SNR = 1.0  # visibilities below this S/N are not used to form closure quantities;
# with this cut our closure counts land within a few of Paper VI's Table 1 (see tests)

# ---------------------------------------------------------------------------
# Physical constants (SI) and unit conversions
# ---------------------------------------------------------------------------
G = 6.67430e-11
C = 2.99792458e8
M_SUN = 1.98840987e30  # IAU 2015 nominal G*M_sun / G
PC = 3.0856775814913673e16
MPC = 1e6 * PC
RAD_PER_UAS = np.pi / 180.0 / 3600.0 / 1e6
UAS_PER_RAD = 1.0 / RAD_PER_UAS
SHADOW_FACTOR_SCHWARZSCHILD = 2.0 * np.sqrt(27.0)  # shadow diameter / (GM/c^2 D) = 10.392


# ---------------------------------------------------------------------------
# Published values we compare against
# ---------------------------------------------------------------------------
# M87 Paper VI (EHT Collaboration 2019f, ApJL 875, L6), Table 3: closure-only
# "xs-ring" crescent fits (dynesty), median and 68% interval, per day and band.
# d = crescent diameter (uas), fw = fractional width, pa = position angle (deg E of N).
PAPER6_XSRING = {
    ("095", "hi"): dict(d=(42.9, 0.59, 0.54), fw=0.39, pa=(160.5, 3.4, 3.1)),
    ("095", "lo"): dict(d=(43.5, 0.44, 0.41), fw=0.20, pa=(157.9, 1.7, 1.8)),
    ("096", "hi"): dict(d=(43.3, 0.44, 0.43), fw=0.34, pa=(149.1, 1.5, 1.5)),
    ("096", "lo"): dict(d=(43.4, 0.27, 0.26), fw=0.20, pa=(148.5, 1.4, 1.2)),
    ("100", "hi"): dict(d=(43.3, 0.79, 0.65), fw=0.50, pa=(194.3, 4.1, 4.6)),
    ("100", "lo"): dict(d=(44.4, 0.84, 0.87), fw=0.41, pa=(204.0, 3.9, 4.2)),
    ("101", "hi"): dict(d=(43.4, 0.74, 0.48), fw=0.44, pa=(180.1, 2.3, 1.8)),
    ("101", "lo"): dict(d=(41.6, 0.51, 0.46), fw=0.50, pa=(175.9, 2.1, 2.0)),
}
# The same table, THEMIS "xs-ringauss" fits to amplitudes + closure phases.
PAPER6_XSRINGAUSS = {
    ("095", "hi"): dict(d=(43.1, 0.35, 0.36), pa=160.6), ("095", "lo"): dict(d=(43.5, 0.27, 0.28), pa=160.9),
    ("096", "hi"): dict(d=(44.1, 0.23, 0.20), pa=146.4), ("096", "lo"): dict(d=(43.5, 0.14, 0.14), pa=153.0),
    ("100", "hi"): dict(d=(42.9, 1.09, 0.86), pa=199.8), ("100", "lo"): dict(d=(43.6, 1.50, 1.92), pa=204.2),
    ("101", "hi"): dict(d=(41.8, 0.46, 0.43), pa=207.4), ("101", "lo"): dict(d=(42.2, 0.43, 0.41), pa=201.1),
}
# Paper VI Table 1: number of closure phases and log closure amplitudes (minimal sets).
PAPER6_COUNTS = {
    ("095", "hi"): (81, 78), ("095", "lo"): (81, 78), ("096", "hi"): (141, 150), ("096", "lo"): (141, 149),
    ("100", "hi"): (53, 56), ("100", "lo"): (47, 48), ("101", "hi"): (110, 117), ("101", "lo"): (107, 113),
}

PUBLISHED = {
    # --- M87 2017 ---
    "m87_ring_diameter_uas": (42.0, 3.0),  # all methods (Paper I, VI)
    "m87_xsring_mean_uas": 43.2,  # weighted mean of the 8 closure-only fits, rms 0.69 (Paper VI)
    "m87_xsringauss_mean_uas": 43.4,  # rms 0.64 (Paper VI)
    "m87_xsring_obs_uncertainty_frac": (0.0175, 0.0169),  # Paper VI Table 3, "observational"
    # GRMHD-calibrated alpha = d / theta_g (Paper VI Table 4, MAD+SANE)
    "alpha_xsring": 11.56,
    "alpha_xsringauss": 11.35,
    "alpha_theory_frac": 0.115,  # ~ (+0.42, -0.43)/3.74 uas: the GRMHD-library spread
    "m87_theta_g_uas": (3.8, 0.4),
    "m87_mass_1e9": (6.5, 0.7),  # 6.5 +/- 0.2 (stat) +/- 0.7 (sys)
    "m87_distance_mpc": (16.8, 0.8),
    # --- dynamical masses, as angles theta_dyn = GM/(c^2 D) (Paper VI Table 9) and masses ---
    "theta_dyn_stars_uas": (3.62, 0.60, 0.34),  # Gebhardt et al. 2011
    "theta_dyn_gas_uas": (2.05, 0.48, 0.16),  # Walsh et al. 2013
    "mass_stars_gebhardt2011_1e9": (6.14, 1.07, 0.62),  # rescaled by Paper VI to D posterior
    "mass_gas_walsh2013_1e9": (3.45, 0.85, 0.26),
    "mass_stars_liepold2023_1e9": (5.37, 0.37, 0.25),  # +/-0.22 distance sys., D = 16.8 Mpc
    "mass_stars_simon2024_1e9": (8.7, 1.8, 1.8),  # 1.2 random (+) 1.3 systematic, D = 16.8 Mpc
    # --- M87 in later years ---
    "m87_2018_diameter_uas": (43.3, 1.5, 3.1),  # EHT 2024 method average
    "m87_2017_2021_mean_diameter_uas": (43.9, 0.6),  # EHT 2025 (arXiv:2509.24593)
    "m87_pa_2017_deg": (150.0, 200.0),  # range across methods (EHT 2025 Table)
    "m87_pa_2018_deg": (209.0, 222.0),
    # --- Sgr A* (EHT 2022 Paper IV, ApJL 930, L15) ---
    "sgra_ring_diameter_uas": (51.8, 2.3),
    "sgra_theta_g_uas": (4.8, 1.4, 0.7),
    "sgra_mass_1e6": (4.0, 1.1, 0.6),
    "sgra_shadow_uas": (48.7, 7.0),
    # variability noise budget measured in Paper IV (Table 3): sigma_var(|u|) = a4 (|u|/4 Glambda)^(-b/2)
    "sgra_var_a4": 0.019,
    "sgra_var_b": 2.4,
}

# H1's own measurement of Sgr A* (17-galactic-centre, joint Keck+VLT GR fit, 30 Sep 2026).
# Used if H1's results file is not available on this machine.
H1_FALLBACK = {
    "mass_1e6": (3.9565, 0.0423),
    "r0_kpc": (7.9311, 0.0438),
    "theta_g_uas": (4.9241, 0.0268),
    "shadow_diameter_uas": (51.1725, 0.2788),
    "mass_distance_correlation": 0.9757,
}

# Scattering towards the Galactic Centre (Johnson et al. 2018, ApJ 865, 104),
# as implemented in ehtim.scattering.ScatteringModel (its defaults are these values).
SGRA_SCATTERING = dict(theta_maj_mas_ref=1.380, theta_min_mas_ref=0.703, pos_ang_deg=81.9,
                       wavelength_reference_cm=1.0)


def ensure_dirs() -> None:
    for d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
        d.mkdir(parents=True, exist_ok=True)
