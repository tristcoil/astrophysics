"""Project configuration: paths, data sources, physical constants and analysis choices.

Every number the project *assumes* lives here, each with its source, so that any
figure or printed result traces back to one line. Published *measurements* we
compare against live in ``literature.py``.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
SP3_DIR = RAW_DIR / "ilrs_orbits"      # weekly ILRS analysis-centre orbits (SP3), one folder per satellite
EOP_DIR = RAW_DIR / "iers"
GRAVITY_DIR = RAW_DIR / "gravity"
TIDE_DIR = RAW_DIR / "tides"
EPHEM_DIR = RAW_DIR / "ephemeris"
TLE_DIR = RAW_DIR / "tle"
LIT_DIR = DATA_DIR / "literature"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, SP3_DIR, EOP_DIR, GRAVITY_DIR, TIDE_DIR, EPHEM_DIR,
           TLE_DIR, LIT_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: The EUROLAS Data Center (EDC, DGFI-TUM, Munich) mirrors the ILRS archive with
#: anonymous FTP. Each ILRS analysis centre (AC) fits a 7-day orbit to the laser
#: ranges every week and publishes it as an SP3 file in the Earth-fixed frame
#: (SLRF2020), 2-minute positions and velocities. Folder = end date of the week.
EDC_ORBITS = "ftp://edc.dgfi.tum.de/pub/slr/products/orbits/TRF/{sat}/{year}/{week}/"

#: Satellites. LAGEOS 1 (1976) and LARES 2 (13 July 2022) fly at nearly the same
#: height in "supplementary" orbits (inclinations i and 180 deg - i), so the drift
#: of their orbital planes caused by the Earth's equatorial bulge is equal and
#: opposite while frame dragging pushes both the same way.
#: Masses/diameters: ILRS satellite pages; LARES 2 from Paolozzi et al. (2023).
SATELLITES = {
    "lageos1": {"name": "LAGEOS", "sp3_id": "L51", "cospar": "7603901", "norad": 8820,
                "mass_kg": 406.965, "diameter_m": 0.600, "launch": dt.date(1976, 5, 4)},
    "lares2": {"name": "LARES 2", "sp3_id": "L67", "cospar": "2208001", "norad": 53105,
               "mass_kg": 294.8, "diameter_m": 0.4244, "launch": dt.date(2022, 7, 13)},
    "lageos2": {"name": "LAGEOS 2", "sp3_id": "L52", "cospar": "9207002", "norad": 22195,
                "mass_kg": 405.38, "diameter_m": 0.600, "launch": dt.date(1992, 10, 22)},
}

#: IERS Earth orientation parameters, the 2020 C04 series (consistent with ITRF2020,
#: hence with the orbits' SLRF2020): polar motion, UT1-UTC, celestial pole offsets.
#: The Paris Observatory copy is current (about a month behind real time; the IERS
#: data centre's "latestVersion" file stopped in May 2025). The last weeks come from
#: the rapid service (finals2000A, IERS Bulletin A).
EOP_C04_URL = "https://hpiers.obspm.fr/iers/eop/eopc04/eopc04.1962-now"
EOP_FINALS_URL = "https://datacenter.iers.org/data/9/finals2000A.all"

#: IERS Conventions (2010), chapter 6: FES2004 ocean tides as geopotential
#: coefficients (main waves, normalized, unit 1e-11).
FES2004_URL = ("https://iers-conventions.obspm.fr/content/chapter6/additional_info/"
               "tidemodels/fes2004_Cnm-Snm.dat")
IERS_CONVENTIONS_URL = "https://iers-conventions.obspm.fr/content/chapter{ch}/icc{ch}.pdf"

#: Static gravity field models from ICGEM (GFZ Potsdam). GOCO06s is the main one;
#: GGM05S and EIGEN-6S4 (v2) are alternatives for the systematics budget.
ICGEM_URL = "https://icgem.gfz-potsdam.de/getmodel/gfc/{key}/{name}.gfc"
GRAVITY_MODELS = {
    "GOCO06s": "32ec2884630a02670476f752d2a2bf1c395d8c8d6d768090ed95b4f04b0d5863",
    "GGM05S": "06a6faa24892df587d29c8a345e09e7031428cf97d4fcc9435b31ae8e4ccc021",
    "EIGEN-6S4 (v2)": "4cc119d62d3f83adce857914bcabdfa7ca2f91677ed2905934cf52698584b4c9",
}

#: JPL planetary ephemeris DE440 (Sun and Moon for third-body forces and tides).
#: N1 (57-orbital-precession) already holds the full file; fetch_data hard-links it.
DE440_SHARED = PROJECT_DIR.parent / "57-orbital-precession" / "data" / "raw" / "de440" / "de440.bsp"
DE440S_URL = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp"

#: UT/CSR monthly C20 (the Earth's flattening) from laser ranging to ~5 satellites, RL06
#: (Cheng & Ries); the independent check of the J2 our two satellites measure.
CSR_C20_URL = "https://download.csr.utexas.edu/pub/slr/degree_2/C20_RL06.txt"

#: CelesTrak's current general-perturbations element sets (no account needed).
CELESTRAK_GP = "https://celestrak.org/NORAD/elements/gp.php?CATNR={norad}&FORMAT=TLE"

# ---------------------------------------------------------------------------
# Physical constants (IERS Conventions 2010, Table 1.1; IAU 2015 B3; DE440)
# ---------------------------------------------------------------------------

C_LIGHT = 299_792_458.0                 # m/s
G_NEWTON = 6.67430e-11                  # m^3 kg^-1 s^-2 (CODATA 2018)
GM_EARTH = 3.986004415e14               # m^3/s^2, TT-compatible (gravity models use it)
R_EARTH = 6_378_136.3                   # m, reference radius of the gravity models
OMEGA_EARTH = 7.292115e-5               # rad/s, nominal mean rotation rate (IERS 2010)
GM_SUN = 1.32712440041e20               # m^3/s^2 (TDB-compatible, DE440)
GM_MOON = 4.902800118e12                # m^3/s^2 (DE440)
AU = 1.495978707e11                     # m
R_SUN = 6.957e8                         # m, IAU 2015 nominal solar radius
SECONDS_PER_DAY = 86_400.0
JULIAN_YEAR_S = 365.25 * SECONDS_PER_DAY
MAS = 3.141592653589793 / 180.0 / 3600.0 / 1000.0   # radians per milliarcsecond

#: Earth's spin angular momentum per unit mass, J/M = C*omega/M (m^2/s).
#: C/(M a^2) = J2/H = 1.0826e-3 / 3.27379e-3 = 0.33069 (dynamical ellipticity H
#: from the IAU 2006 precession), with a = 6378137 m, so
#: J/M = 0.33069 * a^2 * 7.292115e-5 = 9.810e8 m^2/s (IERS 2010 rounds it to 9.8e8).
EARTH_SPIN_PER_MASS = 9.810e8

#: Solar radiation pressure at 1 au: total solar irradiance 1361 W/m^2 (Kopp & Lean 2011) / c.
SOLAR_PRESSURE_1AU = 1361.0 / C_LIGHT   # N/m^2

#: Leap seconds in force over the analysed span (TAI - UTC); none added since 2017.
TAI_MINUS_UTC = ((dt.date(2017, 1, 1), 37.0),)
TT_MINUS_TAI = 32.184

# ---------------------------------------------------------------------------
# Orbit analysis choices
# ---------------------------------------------------------------------------

#: Analysis centres whose weekly orbits we use, in order of preference.
ACS = ("asi", "dgfi", "jcet", "gfz", "esa", "bkg", "nsgf", "cnes")
PRIMARY_AC = "asi"            # ASI is the ILRS primary combination centre

#: Span analysed: from LARES 2's first full weeks of tracking to the last weekly orbit.
ANALYSIS_START = dt.date(2022, 7, 20)
ANALYSIS_END = dt.date(2026, 10, 10)

#: Cache tag of the adopted fits: J2 corrected from the pair ("_c20"), the radiation-pressure
#: scale held at one value per satellite ("_cr"). See campaign.py.
BASE_TAG = "_c20_cr"

GRAVITY_MODEL = "GOCO06s"
GRAVITY_DEGREE = 30           # LAGEOS-type orbits: degree n is damped by (R/a)^n ~ 0.52^n
OCEAN_TIDE_DEGREE = 20
INTEGRATION_STEP_S = 60.0     # fixed-step 8th-order Runge-Kutta
FIT_SAMPLE_S = 600.0          # use every 5th SP3 epoch (2-min sampling) in the arc fits
