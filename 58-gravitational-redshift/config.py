"""Project configuration: paths, data sources, physical constants and analysis choices.

Every number the project *assumes* lives here (published values we compare
against live in ``literature.py``), each with its source, so that any figure
or printed result traces back to one line.
"""

from __future__ import annotations

import datetime as dt
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
SP3_DIR = RAW_DIR / "esa_sp3"          # filtered daily SP3 files (E14, E18 + controls)
SP3_FULL_DIR = RAW_DIR / "esa_sp3_full"  # complete SP3 files, first day of each month
SLR_DIR = RAW_DIR / "esa_slr"          # ESA weekly SLR residual summaries
EOP_DIR = RAW_DIR / "iers"
SOLAR_DIR = RAW_DIR / "solar"
LIT_DIR = DATA_DIR / "literature"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, SP3_DIR, SP3_FULL_DIR, SLR_DIR, EOP_DIR, SOLAR_DIR,
           LIT_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: M1 (50-binary-orbits) measured Sirius B's mass and radius; we reuse its results.
M1_RESULTS = PROJECT_DIR.parent / "50-binary-orbits" / "data" / "results.json"

# ---------------------------------------------------------------------------
# GNSS data: ESA/ESOC Navigation Support Office final products
# ---------------------------------------------------------------------------

#: Daily multi-GNSS final products, open access, one folder per GPS week.
#: ESA0MGNFIN (a 2022/23 NAPEOS 4.8 reprocessing for 2014-2022, operational
#: afterwards) runs to June 2026; ESA0OPSFIN overlaps it from March 2026 on.
ESA_PRODUCTS_URL = "https://navigation-office.esa.int/products/gnss-products/{week}/"
ESA_SERIES = ("ESA0MGNFIN", "ESA0OPSFIN")   # order of preference

#: The two Galileo satellites stranded in eccentric orbits on 22 Aug 2014.
#: PRN -> (SVN, name). E18 = GSAT0201 "Doresa", E14 = GSAT0202 "Milena".
ECCENTRIC = {"E18": ("E201", "GSAT-0201"), "E14": ("E202", "GSAT-0202")}

#: Nominal (near-circular) Galileo satellites kept as controls. They have no
#: once-per-orbit relativistic swing to speak of (e ~ 0.0002-0.0006), so any
#: once-per-orbit signal in their clocks measures the orbit/clock systematics.
CONTROLS = ("E01", "E02", "E11", "E19", "E24", "E26", "E30")

KEEP_SATS = tuple(ECCENTRIC) + CONTROLS

#: ESA's products include E14 and E18 from August 2016 (earlier days lack them).
SP3_START = dt.date(2016, 8, 1)
SP3_END = dt.date(2026, 10, 6)       # last daily file on the server on 2026-10-10

#: GPS-time epoch and leap seconds (GPS - UTC) for the time span used.
GPS_EPOCH = dt.date(1980, 1, 6)
LEAP_SECONDS_GPS = ((dt.date(2015, 7, 1), 17), (dt.date(2017, 1, 1), 18))

# ---------------------------------------------------------------------------
# Physical constants (IERS Conventions 2010, Table 1.1; IAU 2015 Resolution B3)
# ---------------------------------------------------------------------------

C_LIGHT = 299_792_458.0                 # m/s
GM_EARTH = 3.986004418e14               # m^3/s^2 (TCG-compatible)
R_EARTH_EQ = 6_378_136.6                # m, equatorial radius for J2
J2_EARTH = 1.0826359e-3                 # zero-tide (IERS 2010)
W0_GEOID = 62_636_856.0                 # m^2/s^2, geoid potential (IERS 2010)
L_G = 6.969290134e-10                   # 1 - d(TT)/d(TCG) = W0/c^2
OMEGA_EARTH = 7.292115e-5               # rad/s, nominal (only for the budget plot)
GM_SUN = 1.32712440041e20               # m^3/s^2 (TDB-compatible, DE440)
GM_MOON = 4.902800118e12                # m^3/s^2 (DE440)
R_SUN = 6.957e8                         # m, IAU 2015 nominal solar radius
AU = 1.495978707e11                     # m
SECONDS_PER_DAY = 86_400.0

#: GPS satellites' 10.23 MHz clocks are set to 10.22999999543 MHz before launch,
#: a fractional offset of -4.465e-10 that cancels the relativistic gain in orbit
#: (IS-GPS-200, section 3.3.1.1).
GPS_FACTORY_OFFSET = -4.4647e-10

# ---------------------------------------------------------------------------
# Clock analysis choices
# ---------------------------------------------------------------------------

EPOCH_STEP_S = 300.0          # SP3 sampling (5 min); clocks equal the 30-s CLK files at these epochs
LAGRANGE_POINTS = 11          # interpolation order for velocities (10th-degree polynomial)
OUTLIER_SIGMA = 6.0           # robust (MAD) cut on daily-fit residuals
MIN_EPOCHS_PER_DAY = 200      # of 288; days with fewer clean epochs are dropped
JUMP_PS = 300.0               # a step this large between 5-min epochs splits the day into arcs
MC_REALIZATIONS = 400         # Monte Carlo noise realizations for the statistical error
BOOTSTRAP_BLOCK_DAYS = 30     # block length (days) for the block bootstrap
RANDOM_SEED = 20261010
