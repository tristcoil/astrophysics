"""Paths, data sources, targets, physical constants and the published numbers we grade against.

Every literature value used anywhere in the project lives here (or is parsed at run
time from the papers' own LaTeX sources by ``literature.py``), with its source, so a
reader can check each comparison against the paper it came from.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Paths (location-independent, so the folder can be renumbered freely)
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
DSHARP_DIR = RAW_DIR / "dsharp"
IMAGE_DIR = DSHARP_DIR / "images"
PROFILE_DIR = DSHARP_DIR / "profiles"
HLTAU_DIR = RAW_DIR / "hltau"
PAPER_DIR = RAW_DIR / "papers"
CACHE_DIR = DATA_DIR / "cache"  # derived products (geometry fits, profiles, velocity maps)
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"
GAIA_PATH = RAW_DIR / "gaia_dr3_hosts.csv"
EXOPLANET_PATH = RAW_DIR / "nasa_exoplanet_archive_pscomppars.csv"
CO_CUBE_PATH = IMAGE_DIR / "HD163296_CO.fits"

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------
DSHARP_URL = "https://almascience.eso.org/almadata/lp/DSHARP"
# The same release is mirrored at the other ALMA regional centres; fetch_data.py falls
# back to them if ESO is unreachable.
DSHARP_MIRRORS = (
    "https://almascience.eso.org/almadata/lp/DSHARP",
    "https://almascience.nrao.edu/almadata/lp/DSHARP",
    "https://almascience.nao.ac.jp/almadata/lp/DSHARP",
)
HLTAU_URL = "https://almascience.eso.org/almadata/sciver/HLTauBand6/HLTau_Band6+7_ReferenceImages.tgz"
GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_TAP_FALLBACK = "https://gaia.ari.uni-heidelberg.de/tap"  # full DR3 mirror (ARI Heidelberg)
EXOPLANET_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

# arXiv LaTeX sources whose tables are parsed at run time (literature.py), so that the
# published numbers enter the comparison without hand transcription.
ARXIV_SOURCES = {
    "dsharp1": ("1812.04040v1", "Andrews et al. 2018, ApJL 869, L41 (DSHARP I): sample, stellar properties, images"),
    "dsharp2": ("1812.04041v3", "Huang et al. 2018, ApJL 869, L42 (DSHARP II): annular substructures"),
    "dsharp7": ("1812.04045v1", "Zhang et al. 2018, ApJL 869, L47 (DSHARP VII): planet-disc interactions"),
}


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Target:
    key: str  # DSHARP file prefix, e.g. "AS209"
    name: str  # name as printed in the DSHARP papers
    tmass: str  # 2MASS designation (DSHARP I Table 1)
    gaia_dr3: int | None  # Gaia DR3 source_id (from the 2MASS best-neighbour table or a cone search)
    region: str
    multiple: bool = False  # AS 205 and HT Lup are triple systems (DSHARP IV); no ring analysis
    spiral: bool = False  # continuum spirals analysed in DSHARP III
    wedge_deg: float | None = None  # profile averaged only within +/- this of the major axis
    exclude_theta: tuple[tuple[float, float], ...] = ()  # azimuth ranges left out of profiles (deg)
    ridge_exclude: tuple[tuple[float, float], ...] = ()  # azimuths skipped when tracing rings (deg)
    spiral_range_au: tuple[float, float] | None = None  # radii spanned by continuum spiral arms (DSHARP III)
    image_mb: int = 35


# Order: decreasing stellar luminosity, as in DSHARP II's figures.
TARGETS = [
    Target("HD163296", "HD 163296", "17562128-2157218", 4070355800839765248, "isolated",
           exclude_theta=((50.0, 130.0),)),
    Target("HD142666", "HD 142666", "15564002-2201400", 6243489715051060480, "Upper Sco", wedge_deg=20.0),
    Target("Elias24", "Elias 24", "16262407-2416134", 6049168008547148416, "Oph L1688"),
    Target("HTLup", "HT Lup", "15451286-3417305", 6014696841553696768, "Lup I", multiple=True, image_mb=62),
    Target("HD143006", "HD 143006", "15583692-2257153", 6237190612936085120, "Upper Sco",
           exclude_theta=((80.0, 144.0),), ridge_exclude=((80.0, 144.0),)),  # the bright crescent (DSHARP II)
    Target("WaOph6", "WaOph 6", "16484562-1416359", 4326511906643771904, "Oph N 3a", spiral=True,
           spiral_range_au=(25.0, 75.0)),
    Target("IMLup", "IM Lup", "15560921-3756057", 6010135758090335232, "Lup II", spiral=True,
           spiral_range_au=(25.0, 110.0), image_mb=62),
    Target("Elias20", "Elias 20", "16261886-2428196", 6049137600178735744, "Oph L1688", wedge_deg=20.0),
    Target("AS205", "AS 205", "16113134-1838259", 6245891976154382720, "Upper Sco", multiple=True),
    Target("DoAr33", "DoAr 33", "16273901-2358187", 6050679321639383680, "Oph L1688"),
    Target("RULup", "RU Lup", "15564230-3749154", 6010114558131195392, "Lup II"),
    Target("AS209", "AS 209", "16491530-1422087", 4326521359869384576, "Oph N 3a"),
    Target("SR4", "SR 4", "16255615-2420481", 6049142410542091648, "Oph L1688"),
    Target("WSB52", "WSB 52", "16273942-2439155", 6049095406419962112, "Oph L1688", wedge_deg=20.0),
    Target("DoAr25", "DoAr 25", "16262367-2443138", 6049118388788178816, "Oph L1688", wedge_deg=20.0),
    Target("Elias27", "Elias 27", "16264502-2423077", 6049161819498750208, "Oph L1688", spiral=True,
           ridge_exclude=((-117.0, -76.0), (68.0, 123.0)), spiral_range_au=(50.0, 230.0),
           image_mb=62),  # ridge_exclude: where the spiral arms cross the gap D69
    Target("MYLup", "MY Lup", "16004452-4155310", 5995177933191206016, "Lup IV", wedge_deg=20.0),
    Target("Sz129", "Sz 129", "15591647-4157102", 5995168724780802944, "Lup IV"),
    Target("GWLup", "GW Lup", "15464473-3430354", 6014722194741392512, "Lup I"),
    Target("Sz114", "Sz 114", "16090185-3905124", 5997410491550194816, "Lup III", wedge_deg=20.0),
]
TARGET_BY_KEY = {t.key: t for t in TARGETS}
SINGLE_DISCS = [t for t in TARGETS if not t.multiple]  # the 18 discs of DSHARP II
SPIRAL_DISCS = [t for t in TARGETS if t.spiral]

# HL Tau, from the ALMA 2014 Long Baseline Campaign (science verification data). DSHARP II
# folds it in using the combined Band 6+7 (1.0 mm) image; so do we.
HLTAU = Target("HLTau", "HL Tau", "04313843+1813576", 3314312226646863488, "Taurus L1551", wedge_deg=20.0)
HLTAU_IMAGE_NAME = "HLTau_Band6+7_ReferenceImages/HLTau_B6B7cont_mscale_n2_ap.image.tt0.fits"  # combined 1.0 mm (tt0)
# Gaia DR3 has only a two-parameter (position-only) solution for HL Tau, which is deeply
# embedded (A_V ~ 25-33 mag); its distance comes from VLBA parallaxes of the L1551 cloud.
HLTAU_DISTANCE_PC = (147.3, 0.5)  # Galli et al. 2018, ApJ 859, 33 (adopted as 147 pc by DSHARP II)

# Distances. Gaia DR3 parallaxes corrected for the Lindegren et al. (2021) zero point, with the
# formal error inflated by the RUWE (and a 1.3 floor factor plus 10.3 uas in quadrature, after
# Maiz Apellaniz et al. 2021) to cover unresolved multiplicity and disc-contaminated photocentres.
# A parallax known to worse than 5% is replaced by the mean distance of the same cloud's
# DSHARP members (only Elias 27 at G = 20.3 falls under this rule).
DISTANCE_MAX_FRAC_ERROR = 0.05
PARALLAX_ERROR_FLOOR_FACTOR = 1.3
PARALLAX_SYSTEMATIC_MAS = 0.0103


# ---------------------------------------------------------------------------
# Physical constants (cgs, CODATA 2018 / IAU 2015 nominal)
# ---------------------------------------------------------------------------
G_CGS = 6.67430e-8
C_CGS = 2.99792458e10
H_CGS = 6.62607015e-27
K_B_CGS = 1.380649e-16
SIGMA_SB_CGS = 5.670374419e-5
M_H_CGS = 1.6735575e-24  # hydrogen atom
M_SUN_CGS = 1.98840987e33  # IAU 2015 GM_sun / G
L_SUN_CGS = 3.828e33  # IAU 2015 nominal
M_JUP_CGS = 1.89812460e30  # IAU 2015 GM_J / G
M_EARTH_CGS = 5.97216787e27  # IAU 2015 GM_E / G
AU_CGS = 1.495978707e13
PC_CGS = 3.0856775814913673e18
ARCSEC = np.pi / 180.0 / 3600.0  # radians
JY_CGS = 1e-23  # erg s^-1 cm^-2 Hz^-1

# Solar System planets (semi-major axes, au) drawn in the gallery
SOLAR_SYSTEM_AU = {"Jupiter": 5.2029, "Saturn": 9.537, "Uranus": 19.189, "Neptune": 30.070}
KUIPER_BELT_AU = (30.0, 50.0)

# ---------------------------------------------------------------------------
# Analysis settings
# ---------------------------------------------------------------------------
PROFILE_BIN_AU = 1.0  # DSHARP II bins its profiles in 1-au annuli (at the DSHARP I distance)
PROFILE_RMAX_AU = 320.0
WIDTH_INTERP_AU = 0.1  # DSHARP II interpolates profiles onto a 0.1-au grid for widths
GEOMETRY_PRIOR_WIDTH_DEG = 20.0  # half-range of uniform priors around our own starting geometry

# Disc temperature (DSHARP II eq. 1, DSHARP VII eq. 13): T = (phi L / (8 pi r^2 sigma))^(1/4)
FLARING_PHI = 0.02
MEAN_MOLECULAR_WEIGHT = 2.35
DUST_OPACITY_CM2_G = 0.43  # absorption opacity at 1.27 mm for the DSHARP grains (Birnstiel et al. 2018)

# ---------------------------------------------------------------------------
# Published values we compare against (transcribed; parsed tables are in literature.py)
# ---------------------------------------------------------------------------
# Zhang et al. 2018 (DSHARP VII) Table 1: gap width Delta = A K'^B with
# K' = q (h/r)^-0.18 alpha^-0.31, for each maximum Stokes number of the grains.
# Keyed by St_max; values (A, B, -err, +err) with err the 68% spread in log10 K'.
ZHANG_WIDTH_FITS = {
    "gas": (1.05, 0.26, 0.12, 0.03),
    # dust size distribution DSD1 (n(s) ~ s^-3.5, s_max = 0.1 mm); also used for the "1 mm" models
    1.57e-1: (1.09, 0.07, 1.12, 0.86),
    5.23e-2: (1.73, 0.24, 0.63, 0.53),
    1.57e-2: (2.00, 0.36, 0.50, 0.21),
    5.23e-3: (1.25, 0.27, 0.16, 0.22),
    1.57e-3: (1.18, 0.29, 0.17, 0.14),
    5.23e-4: (0.98, 0.25, 0.14, 0.16),
    1.57e-4: (1.11, 0.29, 0.16, 0.13),
}
ZHANG_WIDTH_FITS_DSD2 = {  # n(s) ~ s^-2.5, s_max = 1 cm
    1.57: None,  # no fit (Delta ~ 1: grains drift away)
    5.23e-1: (1.10, 0.05, 1.06, 0.80),
    1.57e-1: (1.13, 0.09, 0.77, 0.70),
    5.23e-2: (1.55, 0.23, 0.59, 0.65),
    1.57e-2: (2.00, 0.36, 0.29, 0.28),
}
ZHANG_ST_UNIT = 5.23e-4  # Table 3 lists St_max in units of this
ZHANG_MIN_WIDTH = 0.15  # narrower gaps are not converted with the fit (their Sec. 4)
ZHANG_ALPHAS = (1e-4, 1e-3, 1e-2)
ZHANG_PREFERRED = ("DSD1", 1e-3)  # "we prefer M_p,am3 with the DSD1 size distribution"

# Huang et al. 2018b (DSHARP III) Table 2: logarithmic spirals R = R0 exp(b theta).
# (R0 au, R0 err, b, b err, pitch deg, pitch err)
HUANG_SPIRALS_LOG = {
    ("Elias27", "S1"): (110.9, 0.6, -0.282, 0.004, 15.7, 0.2),
    ("Elias27", "S2"): (41.3, 0.7, -0.295, 0.004, 16.4, 0.2),
    ("IMLup", "S1a"): (94.0, 7.0, -0.40, 0.04, 22.0, 2.0),
    ("IMLup", "S1b"): (74.6, 1.6, -0.178, 0.018, 10.0, 1.0),
    ("IMLup", "S2a"): (30.0, 2.0, -0.34, 0.04, 19.0, 2.0),
    ("IMLup", "S2b"): (43.0, 3.0, -0.181, 0.018, 10.0, 1.0),
    ("WaOph6", "S1"): (45.9, 0.9, 0.238, 0.016, 13.4, 0.9),
    ("WaOph6", "S2"): (21.3, 1.35, 0.244, 0.016, 13.7, 0.9),
}
# Radial ranges (au) in which DSHARP III traced the arms, and its azimuthal sampling (deg)
HUANG_SPIRAL_RANGES = {
    "Elias27": ((50.0, 110.0, 8.0), (110.0, 230.0, 6.0)),
    "IMLup": ((25.0, 75.0, 15.0), (75.0, 110.0, 6.0)),
    "WaOph6": ((25.0, 50.0, 18.0), (50.0, 80.0, 9.0)),
}

# Teague, Bae & Bergin 2019 (Nature 574, 378), Methods: fit of the same DSHARP 12CO J=2-1
# line-centre map of HD 163296 with a Keplerian disc on a flared emission surface
# z(r) = z0 (r/1")^phi - z1 (r/1")^psi, at d = 101 pc.
TEAGUE2019_HD163296 = {
    "distance_pc": 101.0,
    "mstar": (2.022, 0.004),
    "inc_deg": (46.8, 0.1),
    "pa_deg": (312.83, 0.02),  # position angle of the red-shifted major axis
    "vlsr_kms": (5.763, 0.001),
    "x0_mas": (-5.0, 2.0),
    "y0_mas": (3.0, 2.0),
    "z0_mas": (265.0, 3.0),
    "phi": (1.29, 0.035),
    "z1_mas": (6.0, 2.5),
    "psi": (3.8, 0.3),
    "z0_mas_z1zero": (280.0, 2.0),
    "phi_z1zero": (0.994, 0.007),
    "beam_arcsec": (0.104, 0.095),
}
# The stellar mass DSHARP I Table 1 lists for HD 163296 (from evolutionary tracks): log M = 0.31 (+0.05/-0.03)

# Pinte et al. 2020 (ApJL 890, L9) Table 1: candidate velocity kinks in the DSHARP 12CO cubes.
# (separation arcsec, PA deg, S/N, gap radius au)
PINTE2020_KINKS = {
    "HD163296 #1": (2.20, 3.0, 27, 260.0),
    "HD163296 #2": (0.67, -93.0, 36, 86.0),
    "HD143006": (0.14, -107.0, 10, 22.0),
    "Elias27": (0.32, -6.0, 12, 69.0),
    "IMLup": (0.70, 69.0, 14, 117.0),
    "DoAr25": (0.44, 36.0, 7, 98.0),
    "GWLup": (0.42, 89.0, 12, 74.0),
    "Sz129": (0.29, 148.0, 11, 64.0),
    "WaOph6": (0.51, 37.0, 13, 79.0),
}
# Pinte et al. 2018 (ApJL 860, L13): the HD 163296 kink at ~260 au, a ~2 M_Jup planet.
PINTE2018_HD163296 = {"radius_au": 260.0, "mass_mjup": (2.0, 1.0)}

# ALMA Partnership 2015 (ApJL 808, L3) Table 3: HL Tau geometry from ellipse fits to 8 rings
# (weighted mean) and the ring semi-major axes at their assumed 140 pc.
ALMA2015_HLTAU_GEOMETRY = {"inc_deg": (46.72, 0.05), "pa_deg": (138.02, 0.07), "distance_pc": 140.0}

# Kurtovic et al. 2018 (DSHARP IV) Table 1: 2-D Gaussian fits to the primary discs of the two
# multiple systems (only used to draw projected orbits in the gallery).
KURTOVIC_GEOMETRY = {"AS205": (20.1, 114.0), "HTLup": (48.1, 166.1)}  # (inc, PA) deg

# DSHARP II Table 3: freeze-out temperature ranges (K) of the main volatiles for midplane
# gas densities 1e8-1e12 cm^-3.
SNOWLINES_K = {
    "CO$_2$": (57.0, 72.0),
    "CO (on water ice)": (26.0, 32.0),
    "CO (pure ice)": (19.0, 24.0),
    "N$_2$ (pure ice)": (17.0, 21.0),
}

PUBLISHED = {
    # DSHARP II: the deepest/widest gaps have Delta r / r > 0.25 (HD 163296 D48, Elias 24 D57,
    # HD 143006 D22, AS 209 D61); most widths/r between 0.05 and 0.2.
    "huang_wide_gaps": ("HD163296:D48", "Elias24:D57", "HD143006:D22", "AS209:D61"),
    # Pinte 2018: the kink is at 260 au; Teague 2018 (ApJL 860, L12): pressure bumps imply planets at
    # 83 and 137 au (~1 M_Jup each)
    "teague2018_planets_au": (83.0, 137.0),
}


def ensure_dirs() -> None:
    for d in (DATA_DIR, RAW_DIR, DSHARP_DIR, IMAGE_DIR, PROFILE_DIR, HLTAU_DIR, PAPER_DIR, CACHE_DIR, FIGURE_DIR):
        d.mkdir(parents=True, exist_ok=True)
