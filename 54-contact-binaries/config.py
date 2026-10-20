"""Project configuration: paths, physical constants and every adopted input.

Everything we *assume* lives here, each value with its source: target names and
coordinates, time-system offsets, survey constants and analysis settings.
Published *results* that we grade ourselves against live separately in
``literature.py`` and never enter a fit, except where a module says so
explicitly (e.g. when we reproduce a published fit on purpose).
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
RAW_DIR = DATA_DIR / "raw"
DERIVED_DIR = DATA_DIR / "derived"
OGLE_DIR = RAW_DIR / "ogle"
AAVSO_DIR = RAW_DIR / "aavso"
KEPLER_DIR = RAW_DIR / "kepler"
TESS_DIR = RAW_DIR / "tess"
ZTF_DIR = RAW_DIR / "ztf"
GAIA_DIR = RAW_DIR / "gaia"
LIT_DIR = RAW_DIR / "literature"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR, RAW_DIR, DERIVED_DIR, OGLE_DIR, AAVSO_DIR, KEPLER_DIR,
           TESS_DIR, ZTF_DIR, GAIA_DIR, LIT_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: Plain tool user agent. CDS answers browser-like agents with an "are you a
#: bot" HTML page *with HTTP 200*, so downloads are also checked for HTML.
USER_AGENT = "astrophysics-atlas-M5/1.0 (python-urllib; contact binaries project)"

# ---------------------------------------------------------------------------
# Physical constants (IAU 2015 nominal values; CODATA 2018 for G)
# ---------------------------------------------------------------------------
G_SI = 6.67430e-11
GM_SUN = 1.3271244e20          # m^3 s^-2 (IAU 2015 B3 nominal)
R_SUN_M = 6.957e8               # m
AU_M = 1.495978707e11
DAY_S = 86400.0
YEAR_D = 365.25                 # Julian year in days
JD_J2000 = 2451545.0

# ---------------------------------------------------------------------------
# Archives
# ---------------------------------------------------------------------------
VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
CDS_XMATCH = "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
VIZIER_GAIA_DR3 = "vizier:I/355/gaiadr3"
GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
AAVSO_VSX_API = "https://aavso.org/vsx/index.php"
OGLE_V1309_URL = "https://www.astrouw.edu.pl/ogle/ogle3/V1309_SCO/"
KEPLER_LC_URL = "https://archive.stsci.edu/missions/kepler/lightcurves/"
ARXIV_EPRINT = "https://arxiv.org/e-print/"

# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

#: KIC 9832227 (= 1SWASP J192916.01+463719.9 = TIC 63205800 = Gaia DR3
#: 2128318232616057856), the contact binary predicted to merge in 2022.
KIC9832227 = {
    "kic": 9832227,
    "tic": 63205800,
    "gaia_dr3": 2128318232616057856,
    "ra": 292.316397, "dec": 46.622157,   # Gaia DR3 (J2016.0)
}

#: V1309 Sco = OGLE-III BLG182.7 145645 (README of the OGLE data directory).
V1309SCO = {"ra_hms": "17:57:32.94", "dec_dms": "-30:43:09.9",
            "ogle_files": {"OIII_I": "V1309_SCO_OIII.dat", "OIII_V": "V1309_SCO_OIII_V.dat",
                           "OIV_I": "V1309_SCO.dat", "OIV_V": "V1309_SCO_V.dat"}}

#: T CrB in the AAVSO International Database (VSX OID 10602, AUID 000-BBW-825).
TCRB = {"name": "T CrB", "auid": "000-BBW-825", "ra": 239.87567, "dec": 25.92017}

# ---------------------------------------------------------------------------
# Time systems
# ---------------------------------------------------------------------------

#: Kepler light-curve TIME is BJD_TDB - 2454833 (Kepler Archive Manual).
KEPLER_BJD_OFFSET = 2454833.0
#: TESS light-curve TIME is BTJD = BJD_TDB - 2457000.
TESS_BJD_OFFSET = 2457000.0
#: OGLE times are HJD (UTC-based, full JD). For V1309 Sco the
#: heliocentric/barycentric difference (< 4 s) is irrelevant at a 1.4-day period.

#: One Kepler long-cadence sample: 270 frames of 6.5387 s = 1765.46 s.
KEPLER_LC_CADENCE_D = 1765.4615 / DAY_S

# ---------------------------------------------------------------------------
# AAVSO download plan for T CrB
# ---------------------------------------------------------------------------
#: Windows (JD start, JD end, step in days). The VSX API returns every
#: observation in a window; windows are split automatically if one returns
#: more than AAVSO_MAX_PER_REQUEST rows. Dense recent years get short steps.
AAVSO_WINDOWS = [
    (2400000.0, 2431000.0, 5000.0),     # 1858-1943: the 1866 eruption and quiescence
    (2431000.0, 2433000.0, 200.0),      # 1943-1949: the 1946 eruption
    (2433000.0, 2457000.0, 1000.0),     # 1949-2014
    (2457000.0, 2461400.0, 60.0),       # 2014-2026: super-active state and the dip
]
AAVSO_MAX_PER_REQUEST = 20000
AAVSO_PAUSE_S = 1.0

# ---------------------------------------------------------------------------
# Analysis settings
# ---------------------------------------------------------------------------

#: Cross-match radius between ZTF positions (epoch ~2018-2019) and Gaia DR3
#: (J2016.0). Contact binaries in ZTF are mostly > 200 pc away (proper motions
#: of tens of mas/yr at most), so 1.5" is generous.
XMATCH_RADIUS_ARCSEC = 1.5

#: ZTF variable types kept from Chen et al. (2020) for the binary analysis.
ZTF_BINARY_TYPES = ("EW", "EA")

#: Random seed used everywhere a random number is drawn.
SEED = 20261009
