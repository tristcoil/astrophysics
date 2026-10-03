"""Project configuration: paths, API queries, physical constants, and every cut.

Every selection criterion, conversion constant and published comparison value
used anywhere in the project is defined here, once, so that each number quoted
in the README or printed by ``run.py`` traces back to a single line.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
#: One JSON record per cached download: URL, parameters, retrieval time, size.
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data services: JPL Solar System Dynamics (SSD) and CNEOS. All free, no key.
# ---------------------------------------------------------------------------

CAD_URL = "https://ssd-api.jpl.nasa.gov/cad.api"  # close approaches
SENTRY_URL = "https://ssd-api.jpl.nasa.gov/sentry.api"  # impact monitoring
FIREBALL_URL = "https://ssd-api.jpl.nasa.gov/fireball.api"  # bolides
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"  # one object in detail
SBDB_QUERY_URL = "https://ssd-api.jpl.nasa.gov/sbdb_query.api"  # many objects
HORIZONS_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"  # ephemerides

#: Close approaches of near-Earth *asteroids* to Earth within 0.05 au, the
#: distance used to define "potentially hazardous". Comets are left out: they
#: rarely come this close, and their brightness does not measure their size.
#: The table covers 1900-2200, JPL's full range. It is fetched in chunks so no
#: single request is large; each chunk is cached exactly as returned.
CAD_DIST_MAX_AU = 0.05
CAD_PARAMS = {
    "dist-max": str(CAD_DIST_MAX_AU),
    "nea": "true",
    "diameter": "true",  # include measured diameters where JPL has one
    "fullname": "true",
}
CAD_CHUNK_EDGES = [
    "1900-01-01", "1950-01-01", "2000-01-01", "2010-01-01", "2020-01-01",
    "2030-01-01", "2060-01-01", "2100-01-01", "2150-01-01", "2200-01-01",
]

#: The near-Earth-asteroid catalogue: one row per known NEA, for discovery
#: dates (first observation) and orbit quality. Joined to approaches on the
#: primary designation.
SBDB_QUERY_FIELDS = [
    "spkid", "pdes", "full_name", "first_obs", "last_obs", "H", "diameter",
    "albedo", "class", "pha", "moid", "condition_code", "n_obs_used", "data_arc",
]

#: Apophis, for the 2029 flyby. SPK-ID 2099942 = numbered asteroid 99942.
APOPHIS_DES = "99942"
APOPHIS_HORIZONS_ID = "99942;"  # trailing ';' = small body, not a major-body ID
#: Every Earth approach within this distance, for the "why we know" panel.
APOPHIS_CAD_DIST_MAX_AU = 0.5
#: Geocentric position/velocity: a coarse week and a fine window around the
#: closest approach (13 April 2029, 21:46 UTC).
APOPHIS_VECTOR_WINDOWS = [
    ("2029-04-10 00:00", "2029-04-17 00:00", "10 m"),
    ("2029-04-13 18:00", "2029-04-14 02:00", "1 m"),
]
#: Heliocentric osculating orbit, daily, across the flyby: the orbit Earth's
#: gravity leaves behind is not the one Apophis arrived on.
APOPHIS_ELEMENTS_WINDOW = ("2029-03-01", "2029-05-31", "1 d")
#: The Moon, geocentric, across one full orbit centred on the flyby.
MOON_VECTOR_WINDOW = ("2029-03-30 00:00", "2029-04-28 00:00", "1 h")


# ---------------------------------------------------------------------------
# Physical constants and conversions
# ---------------------------------------------------------------------------

AU_KM = 149_597_870.7  # IAU 2012, exact
#: One lunar distance: the Moon's average distance from Earth, rounded as NASA
#: rounds it. The Moon's actual distance varies from ~356,500 to ~406,700 km.
LD_KM = 384_400.0
EARTH_RADIUS_KM = 6_371.0  # mean radius
GM_EARTH_KM3_S2 = 398_600.4418
#: Geostationary satellites orbit 35,786 km above the equator.
GEO_RADIUS_KM = 42_164.0
#: GPS satellites: semi-major axis 26,560 km.
GPS_RADIUS_KM = 26_560.0
KT_TNT_J = 4.184e12
MT_TNT_J = 4.184e15

#: Diameter from absolute magnitude H and geometric albedo p:
#: D = 1329 km / sqrt(p) * 10^(-H/5)   (e.g. Pravec & Harris 2007).
H_TO_D_KM = 1329.0
#: A typical NEA albedo. With it, H = 22 is the familiar 140 m and H = 17.75
#: is 1 km. Real albedos run from ~0.05 (dark) to ~0.25+ (bright), which
#: scales a diameter by 1.67x to 0.75x.
ALBEDO_TYPICAL = 0.14
ALBEDO_RANGE = (0.05, 0.25)
#: Bulk density for masses and impact energies: JPL Sentry's convention.
DENSITY_KG_M3 = 2_600.0

#: Size classes (diameter in metres), largest first. Chelyabinsk (2013) was
#: ~18 m; Tunguska (1908) roughly 50 m; 140 m is the survey goal set by the
#: US Congress in 2005; 1 km is the threshold for global effects.
SIZE_CLASS_EDGES_M = [1000.0, 140.0, 50.0, 20.0, 10.0]
SIZE_CLASS_LABELS = ["≥ 1 km", "140 m – 1 km", "50 – 140 m", "20 – 50 m", "10 – 20 m", "< 10 m"]


# ---------------------------------------------------------------------------
# Time windows
# ---------------------------------------------------------------------------

#: "What's passing by": the month before the data were retrieved and the two
#: months after. The asymmetry is deliberate -- the past side shows what the
#: surveys found only as, or after, it went by.
SNAPSHOT_DAYS_BEFORE = 30
SNAPSHOT_DAYS_AFTER = 60

#: The modern survey era, for "how often" statistics of *observed* approaches:
#: complete calendar years only.
SURVEY_YEARS = (2016, 2025)  # inclusive: ten full years
#: Longer run of years for the discovery-history panel.
HISTORY_YEARS = (2000, 2025)
#: Recent five years for the discovery-timing statistics.
RECENT_YEARS = (2021, 2025)


# ---------------------------------------------------------------------------
# From close approaches to impact rates
# ---------------------------------------------------------------------------

#: Encounter geometry (Öpik theory): for a steady stream of asteroids crossing
#: Earth's neighbourhood, the number passing within distance d grows as d^2 --
#: the area of a target disc -- apart from gravitational focusing. The d^2 law
#: is tested on the data and, where it holds, extrapolated down to Earth's
#: surface. Two samples:
#:
#: 1. Large asteroids, whose catalogue is nearly complete: every computed
#:    approach 1900-2200 within 0.05 au. Their orbits are known well enough
#:    that JPL computes approaches for centuries either way.
#: 2. Small asteroids, found mostly as they pass: only approaches observed
#:    during the survey years, within 1 lunar distance. These give *lower
#:    limits*, because many small passers are never seen.
LARGE_WINDOW = ("1900-01-01", "2200-01-01")
LARGE_MIN_DIAMETER_M = 140.0
LARGE_REF_DIST_AU = 0.05
SMALL_REF_DIST_LD = 1.0
#: Distance range over which the slope of log N(<d) vs log d is fitted.
SLOPE_FIT_RANGE_AU = (0.005, 0.05)
N_BOOTSTRAP = 2000
RANDOM_SEED = 20290413  # Apophis's flyby date


# ---------------------------------------------------------------------------
# Published results used for comparison (never as inputs)
# ---------------------------------------------------------------------------

#: Stuart & Binzel (2004, Icarus 170, 295): asteroids >= 1 km strike Earth once
#: every 0.60 +/- 0.1 Myr on average.
STUART_BINZEL_1KM_INTERVAL_YR = (0.60e6, 0.10e6)
#: Annual impact probability for NEAs >= 140 m, as compiled by Nugent et al.
#: (2025, PSJ, arXiv:2508.02418): Granvik et al. 2018 ~1e-5; Harris & Chodas
#: 2021 and Nesvorny et al. 2024 4e-5; Nugent et al. 9.1e-5.
PUBLISHED_140M_RATES = {
    "Granvik et al. 2018": 1.0e-5,
    "Harris & Chodas 2021": 4.0e-5,
    "Nesvorný et al. 2024": 4.0e-5,
    "Nugent et al. 2025": 9.1e-5,
}
#: Brown et al. (2002, Nature 420, 294), satellite-observed bolides:
#: N(>=E) = 3.7 E^-0.90 impacts per year, E in kilotons, fitted 0.1-10 kt.
BROWN2002_A, BROWN2002_B = 3.7, 0.90
BROWN2002_RANGE_KT = (0.1, 10.0)
#: The same paper's conversion from optical (radiated) to total impact energy,
#: E_impact = 8.2508 E_rad^0.885, both in kt. JPL applies it to every fireball.
BROWN2002_RAD_A, BROWN2002_RAD_B = 8.2508, 0.885
#: Palermo scale background impact frequency, Chesley et al. (2002, Icarus 159,
#: 423): f_B = 0.03 E^-4/5 per year, E in megatons.
PALERMO_FB_A, PALERMO_FB_B = 0.03, 0.8

#: The fireball record: JPL's table starts in 1988, but coverage was sparse
#: early on. Events above 1 kt were logged at 1.2 a year in 1995-1999 against
#: 3.4-4.8 a year in every five-year block since (analysis.py reports the
#: blocks), so rates are measured from 2000.
FIREBALL_YEARS = (2000, 2025)


# ---------------------------------------------------------------------------
# Named objects
# ---------------------------------------------------------------------------

#: Removed from the risk list after further observations -- the usual fate of
#: a scary headline. Highlighted on the risk-list history.
FAMOUS_REMOVED = {"99942": "Apophis", "2024 YR4": "2024 YR4"}
