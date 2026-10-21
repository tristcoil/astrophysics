"""Project configuration: paths, data sources, physical constants and analysis choices.

Every number the project *assumes* lives here, each with its source, so that any
figure or printed result traces back to one line. Published measurements we
compare against live in ``literature/`` (transcribed and checked against the
primary papers); see ``literature.py``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
APOLLO_DIR = RAW_DIR / "apollo"
NAIF_DIR = RAW_DIR / "naif"
IERS_DIR = RAW_DIR / "iers"
LIT_DATA_DIR = DATA_DIR / "literature"      # downloaded papers and memos (not committed)
DERIVED_DIR = DATA_DIR / "derived"           # caches written by the analysis
LITERATURE_DIR = PROJECT_DIR / "literature"  # transcribed literature tables (committed)
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, APOLLO_DIR, NAIF_DIR, IERS_DIR, LIT_DATA_DIR, DERIVED_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

ATLAS_DIR = PROJECT_DIR.parent

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: APOLLO lunar laser ranging normal points (T. Murphy, UC San Diego). The page
#: lists every release; these are the current versions of groups a-f, 2006-04-07
#: to 2020-12-27 (groups e and f were re-reduced and replaced on 2022-01-05).
APOLLO_URL = "https://tmurphy.physics.ucsd.edu/apollo/{name}"
APOLLO_PAGE = "https://tmurphy.physics.ucsd.edu/apollo/norm_pts.html"
APOLLO_FILES = {
    "a": "121106a_apollo_np.txt",
    "b": "121106b_apollo_np.txt",
    "c": "131219c_apollo_np.txt",
    "d": "160725d_apollo_np.txt",
    "e": "160912e_apollo_np2.txt",
    "f": "190823f_apollo_np2.txt",
}
APOLLO_NOTES = ("121106_notes.txt", "131219_notes.txt", "160725_notes.txt",
                "160912_notes.txt", "190823_notes.txt")

#: JPL ephemeris DE440 (Park et al. 2021, AJ 161, 105) from NAIF: planetary and
#: lunar positions (TDB), the Moon's principal-axis orientation (libration Euler
#: angles integrated with the orbit) and the frame kernel that documents it.
NAIF_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/{path}"
NAIF_FILES = {
    "de440s.bsp": "spk/planets/de440s.bsp",
    "moon_pa_de440_200625.bpc": "pck/moon_pa_de440_200625.bpc",
    "moon_de440_250416.tf": "fk/satellites/moon_de440_250416.tf",
    "gm_de440.tpc": "pck/gm_de440.tpc",
}

#: IERS EOP 20 C04 series (daily polar motion, UT1-UTC and celestial pole offsets).
IERS_EOP_URL = "https://datacenter.iers.org/data/latestVersion/EOP_20_C04_one_file_1962-now.txt"
IERS_EOP_FILE = "EOP_20_C04_one_file_1962-now.txt"

#: Documents whose models and tables we implement or check against.
DOCUMENTS = {
    # Williams & Boggs (2020), JPL IOM 335N-20-01: the JPL lunar laser range model
    "jpl/LLR_Model_2020_DR.pdf": "https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/LLR_Model_2020_DR.pdf",
    # Park et al. (2021), AJ 161, 105: DE440 and DE441 (Table 1: reflector coordinates)
    "jpl/de440.aj.pdf": "https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/de440.aj.pdf",
    # IERS Conventions (2010), updated chapters 7 (tides), 8 (sub-daily EOP), 9 (troposphere)
    "iers/icc7.pdf": "https://iers-conventions.obspm.fr/content/chapter7/icc7.pdf",
    "iers/icc8.pdf": "https://iers-conventions.obspm.fr/content/chapter8/icc8.pdf",
    "iers/icc9.pdf": "https://iers-conventions.obspm.fr/content/chapter9/icc9.pdf",
    "iers/sw/ORTHO_EOP.F": "https://iers-conventions.obspm.fr/content/chapter8/software/ORTHO_EOP.F",
    "iers/sw/CNMTX.F": "https://iers-conventions.obspm.fr/content/chapter8/software/CNMTX.F",
    "iers/sw/FCUL_A.F": "https://iers-conventions.obspm.fr/content/chapter9/software/FCUL_A.F",
    "iers/sw/FCUL_ZD_HPA.F": "https://iers-conventions.obspm.fr/content/chapter9/software/FCUL_ZD_HPA.F",
    # Will (2014), Living Rev. Relativ. 17, 4 (arXiv:1403.7377): TeX source with its tables
    "will2014/src.tar.gz": "https://arxiv.org/src/1403.7377v1",
}

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

C_LIGHT = 299_792_458.0                     # m/s (exact)
C_KM = C_LIGHT / 1e3
AU = 1.495978707e11                         # m (IAU 2012, exact)
DAY = 86_400.0
JULIAN_YEAR = 365.25 * DAY

#: DE440 GM values (gm_de440.tpc, TDB-compatible), km^3/s^2.
GM_SUN = 132_712_440_041.279419
GM_EARTH = 398_600.435507
GM_MOON = 4_902.800118
GM_EMB = GM_EARTH + GM_MOON

#: Time-scale rate constants (IAU 2006 B3, IERS Conventions 2010 ch. 10).
L_B = 1.550519768e-8          # TDB vs TCB
L_G = 6.969290134e-10         # TT vs TCG
L_C = 1.48082686741e-8        # mean rate of TCB vs TCG (geocentre)
L_M = 1.4825e-8               # Williams & Boggs 2020, eq. (17): the same for the Moon

#: Earth figure (IERS 2010 / GRS80 ellipsoid for geodetic latitude).
R_EARTH_EQ = 6_378_136.6      # m, equatorial radius used in the tide formulae
WGS_A = 6_378_137.0           # m, GRS80 semi-major axis
WGS_F = 1 / 298.257222101     # GRS80 flattening
OMEGA_EARTH = 7.292115e-5     # rad/s, nominal rotation rate (station velocity only)

#: Moon radius used in the lunar tide formula (JPL lunar gravity reference radius).
R_MOON = 1_738_000.0          # m

# ---------------------------------------------------------------------------
# Lunar laser ranging model (Williams & Boggs 2020 = JPL IOM 335N-20-01)
# ---------------------------------------------------------------------------

#: APOLLO's telescope, Apache Point Observatory: approximate *geocentric* spherical
#: coordinates from the APOLLO normal-point page ("careful: not geodetic"). The fit
#: corrects them; plate motion is fitted as a velocity.
APOLLO_GEOCENTRIC = {"radius_km": 6374.69213, "lat_deg": 32.6054889, "lon_deg": 254.1795778}
APOLLO_STATION_ID = 70610

#: Retroreflector array coordinates in the DE440 principal-axis (PA) frame, metres
#: (Park et al. 2021, Table 1). APOLLO reflector codes: 0=A11 1=L1 2=A14 3=A15 4=L2.
REFLECTORS_PA = {
    "Apollo 11": (1591967.049, 690698.573, 21004.461),
    "Lunokhod 1": (1114291.452, -781299.273, 1076059.049),
    "Apollo 14": (1652689.369, -520998.431, -109729.869),
    "Apollo 15": (1554678.104, 98094.498, 765005.863),
    "Lunokhod 2": (1339363.598, 801870.995, 756359.260),
}
APOLLO_REFLECTOR_CODE = {0: "Apollo 11", 1: "Lunokhod 1", 2: "Apollo 14", 3: "Apollo 15", 4: "Lunokhod 2"}

#: Corner-cube delay n*h/c added to each leg (W&B 2020 sec. 3.18).
CC_INDEX = 1.4555
CC_DEPTH = {"Apollo 11": 0.0275, "Apollo 14": 0.0275, "Apollo 15": 0.0275,
            "Lunokhod 1": 0.0435, "Lunokhod 2": 0.0435}

#: Love numbers (W&B 2020 sec. 3.9 and 3.12); "real" values, constant part included.
EARTH_LOVE = {"h2": 0.6078, "l2": 0.0847, "h3": 0.292, "l3": 0.015}
MOON_LOVE = {"h2": 0.0423, "l2": 0.0107, "h3": 0.0234, "l3": 0.0030}

#: Step-2 correction for the K1 tide and its nodal sideband (IERS 2010 Table 7.3a),
#: radial in-phase/out-of-phase amplitudes in metres, applied as eq. (12a).
K1_STEP2 = {"K1": (12.00e-3, -0.78e-3), "K1_side": (1.73e-3, -0.12e-3)}

#: Solid-Earth pole tide (W&B 2020 sec. 3.6): mean pole and their coefficients.
POLE_TIDE_MEAN = {"x0": 0.056, "xdot": 0.0021, "y0": 0.346, "ydot": 0.0037}   # arcsec, arcsec/yr from 2000
POLE_TIDE_COEF = {"radial": 32.4e-3, "horizontal": 9.0e-3}                    # m per arcsec

#: Troposphere (IERS 2010 ch. 9): CO2 content for Mendes-Pavlis.
CO2_PPM = 375.0

#: PPN gamma in the Shapiro delay (GR = 1).
GAMMA_PPN = 1.0

#: APOLLO's recommended error scalings per group (release notes): factor or RSS term.
#: Groups e and f were released with no recommendation; we treat them like d.
APOLLO_ERROR_SCALE = {"a": (2.0, 17e-12), "b": (6.0, 60e-12), "c": (2.5, 20e-12),
                      "d": (1.0, 0.0), "e": (1.0, 0.0), "f": (1.0, 0.0)}

#: Analysis choices for the residual fit.
LLR_MODEL_FLOOR_M = 0.015       # one-way model-error floor added in quadrature (m); tuned so chi2/dof ~ 1
LLR_OUTLIER_SIGMA = 5.0         # iterative rejection threshold (robust sigma of residuals)
LLR_REFERENCE_EPOCH = 2013.0    # station position epoch (decimal year)
