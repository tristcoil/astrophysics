"""Project configuration: paths, data sources, physical constants and adopted inputs.

Every number the project *assumes* lives here (published results we grade
ourselves against live in ``transcribed/``), each with its source, so that any
figure or printed result traces back to one line.
"""

from __future__ import annotations

import math
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
TRANSCRIBED_DIR = PROJECT_DIR / "transcribed"   # hand-transcribed tables (committed)
DATA_DIR = PROJECT_DIR / "data"                 # downloads and caches (not committed)
RAW_DIR = DATA_DIR / "raw"
DYSON_DIR = RAW_DIR / "dyson1920"
PAGES_DIR = DYSON_DIR / "pages"
WALLAL_DIR = RAW_DIR / "wallal"
HARVEY_DIR = RAW_DIR / "harvey1979"
GAIA_DIR = RAW_DIR / "gaia"
HORIZONS_DIR = RAW_DIR / "horizons"
SKY_DIR = RAW_DIR / "sky"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, DYSON_DIR, PAGES_DIR, WALLAL_DIR, HARVEY_DIR, GAIA_DIR,
           HORIZONS_DIR, SKY_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: Dyson, Eddington & Davidson 1920, Phil. Trans. R. Soc. A 220, 291-333, as
#: scanned by the Royal Society for the Internet Archive. The item's PDF
#: endpoint returns HTTP 500 (October 2026); the page images work. Leaf n is
#: printed page 291 + n (45 leaves).
DYSON_IA_ITEM = "philtrans06337895"
DYSON_PAGE_URL = "https://archive.org/download/philtrans06337895/page/n{leaf}_w2560.jpg"
DYSON_N_LEAVES = 45

#: Scanned papers from the NASA ADS article service (image-only PDFs).
ADS_PDF_URL = "https://articles.adsabs.harvard.edu/pdf/{bibcode}"
ADS_SCANS = {
    "campbell_trumpler_1923": ("1923LicOB..11...41C", WALLAL_DIR),  # Wallal, preliminary
    "campbell_trumpler_1928": ("1928LicOB..13..130C", WALLAL_DIR),  # Wallal, five-foot cameras
    "harvey_1979": ("1979Obs....99..195H", HARVEY_DIR),             # 1919 Sobral plates remeasured
}

#: Gaia DR3 via VizieR (I/355/gaiadr3); the Hipparcos new reduction (I/311).
VIZIER_GAIA = "I/355/gaiadr3"
VIZIER_HIP2 = "I/311/hip2"

#: JPL Horizons API (vectors for the Sun, the Earth and Cassini).
HORIZONS_API = "https://ssd.jpl.nasa.gov/api/horizons.api"
CASSINI_ID = "-82"

#: CDS hips2fits, for a real image of the 1919 field.
HIPS2FITS_URL = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"

# ---------------------------------------------------------------------------
# Physical constants (IAU 2012/2015 nominal values, CODATA)
# ---------------------------------------------------------------------------

GM_SUN = 1.32712440041e20   # m^3 s^-2 (TDB-compatible; IAU 2015 nominal 1.3271244e20)
C_LIGHT = 299792458.0       # m/s
R_SUN = 6.957e8             # m (IAU 2015 nominal photospheric radius)
AU = 1.495978707e11         # m
ARCSEC = math.pi / 180.0 / 3600.0
DAY = 86400.0

#: GR's deflection of a ray grazing the Sun: 4 GM / (c^2 R) = 1.7512".
#: (1 + gamma)/2 times this in PPN; Soldner's Newtonian value is half.
DEFLECTION_LIMB_GR = 4.0 * GM_SUN / (C_LIGHT**2 * R_SUN) / ARCSEC
#: The round number used in 1919 (and by Dyson, Eddington & Davidson).
DEFLECTION_LIMB_1919 = 1.75

# ---------------------------------------------------------------------------
# The 1919 eclipse (Dyson, Eddington & Davidson 1920)
# ---------------------------------------------------------------------------

#: Micrometer screw values: Sobral 4-inch plates measured at Greenwich, 6.25"
#: per revolution (p.302); Sobral astrographic plates, 12.3" (p.310, Table IX).
REV_4INCH = 6.25
REV_ASTRO = 12.3
#: Principe measures (Cambridge machine) are in units of the fifth decimal of
#: a 5 mm (~5') reseau interval, = 0.003" (p.318); kappa = 184 units at 50'
#: is the full 1.75" deflection (p.322).
UNIT_PRINCIPE = 0.003
KAPPA_FULL = 184.0
#: Table I's coordinates are in units of 50'; the Sun's radius during the
#: eclipse was 15'.8 (p.306) or 15'.78 (p.322).
UNIT_TABLE1_ARCMIN = 50.0
SUN_RADIUS_1919_ARCMIN = 15.8
#: The Sun's true apparent radius that day (JPL Horizons, topocentric, Sobral and
#: Principe alike: R_sun = 695,700 km at 1.0137 au). Our own fits convert the
#: deflection at 50' to the limb with this; replications of the 1920 arithmetic
#: keep the paper's 15'.8.
SUN_RADIUS_TRUE_ARCMIN = 15.771

#: 1919 times are astronomical G.M.T., reckoned from noon: Sobral's second
#: contact "May 28, 23h 58m 18s G.M.T." (p.300) is civil 1919-05-29 11:58:18 UT.
SOBRAL_SECOND_CONTACT_UT = "1919-05-29T11:58:18"

#: Table I stars identified in Gaia DR3 by position (our match, rms 1.6"
#: against Table I's 3"-rounded coordinates; see stars.py). The paper's star
#: 11 "66 Tauri" sits where 56 Tau is (66 Tau is 12 degrees away), and star 9
#: "Piazzi IV 53" is the K giant HR 1370 = HD 27639 A (its 8.6-mag companion
#: BU 87 B is 1.9" away).
TABLE1_GAIA = {
    1: (145373377272257664, "BD+21 641", "HD 27808"),
    2: (145391484855481344, "HR 1403", "HD 28226"),
    3: (145484634106814592, "kappa2 Tau", "HD 27946"),
    4: (145674574739951872, "kappa1 Tau", "HD 27934"),
    5: (145259990135756928, "HR 1375", "HD 27742"),
    6: (145727729255780736, "upsilon Tau", "HD 28024"),
    7: (49251390712484864, "BD+20 741", "HD 27533"),
    8: (49231668222673920, "BD+20 740", "HD 27524"),
    9: (49177620354353152, "HR 1370", "HD 27639 A"),
    10: (146100115802055168, "72 Tau", "HD 28149"),
    11: (49520461823548288, "56 Tau", "HD 27309"),
    12: (49415527182759296, "53 Tau", "HD 27295"),
    13: (145769682495584896, "BD+22 688", "HD 27832"),
}

#: Epoch of the eclipse plates (mid-totality, Sobral) as a Julian year, and the
#: equinox of Table I's standard coordinates (fitted: rotation 0.01 deg).
EQUINOX_1919 = 1919.41
