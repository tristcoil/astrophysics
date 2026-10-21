"""Project configuration: paths, the archive query, and every sample cut.

Every selection criterion used anywhere in the project is defined here, once,
so that each number quoted in the README or printed by ``run.py`` traces back to
a single line. Change a cut here and every figure and summary follows it.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
MODEL_DIR = DATA_DIR / "models"
FIGURE_DIR = PROJECT_DIR / "figures"
CATALOG_PATH = DATA_DIR / "pscomppars.csv"
CATALOG_META_PATH = DATA_DIR / "pscomppars_meta.json"
MODEL_META_PATH = MODEL_DIR / "models_meta.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, MODEL_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# The archive query
# ---------------------------------------------------------------------------

#: NASA Exoplanet Archive Table Access Protocol endpoint for synchronous queries.
TAP_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

#: The Planetary Systems Composite Parameters table: one row per confirmed
#: planet, each parameter filled from whichever reference the archive judged
#: best. Complete, but not self-consistent -- see catalog.py for what that costs.
TABLE = "pscomppars"
TABLE_DOI = "10.26133/NEA13"

#: Columns fetched. The ``*_reflink`` columns matter more than they look: they
#: are the only way to tell a *measured* radius or mass from one the archive
#: filled in with an empirical formula.
COLUMNS = [
    # identity and flags
    "pl_name", "hostname", "sy_snum", "sy_pnum", "pl_controv_flag",
    "tran_flag", "rv_flag",
    # discovery
    "discoverymethod", "disc_year", "disc_facility", "disc_refname",
    # orbit
    "pl_orbper", "pl_orbpererr1", "pl_orbpererr2", "pl_orbperlim",
    "pl_orbsmax", "pl_orbsmaxlim", "pl_orbeccen",
    # radius (Earth radii)
    "pl_rade", "pl_radeerr1", "pl_radeerr2", "pl_radelim", "pl_rade_reflink",
    # mass or minimum mass (Earth masses); pl_bmassprov says which
    "pl_bmasse", "pl_bmasseerr1", "pl_bmasseerr2", "pl_bmasselim",
    "pl_bmasse_reflink", "pl_bmassprov",
    # irradiation (Earth = 1) and equilibrium temperature
    "pl_insol", "pl_insolerr1", "pl_insolerr2", "pl_insollim", "pl_insol_reflink",
    "pl_eqt",
    # host star
    "st_teff", "st_tefferr1", "st_tefferr2", "st_tefflim",
    "st_rad", "st_raderr1", "st_raderr2", "st_mass", "st_lum", "st_logg",
    "st_met", "st_spectype",
    # position, kept for cross-project use (e.g. a 3D map of planet hosts)
    "sy_dist", "ra", "dec",
]


# ---------------------------------------------------------------------------
# Published planet-interior models
# ---------------------------------------------------------------------------

#: Mass-radius curves of Zeng et al. (2019, PNAS 116, 9723), as tabulated on
#: Li Zeng's website. Model curves, not data: they say what radius a planet of a
#: given mass and composition would have.
ZENG_BASE_URL = "https://lweb.cfa.harvard.edu/~lzeng/tables/"

#: Temperature-dependent curves (water, hydrogen envelope) are taken at 700 K,
#: the tabulated temperature nearest the ~790 K median equilibrium temperature
#: of the small planets in the mass-radius sample. Cooler planets of the same
#: composition would be somewhat smaller.
ZENG_TEMPERATURE_K = 700

#: key -> (file name, legend label)
ZENG_MODELS: dict[str, tuple[str, str]] = {
    "iron": ("massradiusFe.txt", "Pure iron"),
    "earthlike": ("massradiusEarthlikeRocky.txt", "Earth-like (⅓ iron, ⅔ rock)"),
    "rock": ("massradiusmgsio3.txt", "Pure rock (no iron)"),
    "water50": ("massradius_50percentH2O_700K_1mbar.txt", "Half rock, half water"),
    "water100": ("massradius_100percentH2O_700K_1mbar.txt", "Pure water"),
    "h2_1pct": ("Earthlike1h700K1mbar.txt", "Earth-like + 1% hydrogen"),
    "hydrogen": ("massradiushydrogen.txt", "Pure hydrogen, cold"),
}

#: Pure-water curves at every tabulated temperature. The "too big even for
#: water" test compares each planet with the curve nearest its own equilibrium
#: temperature, because a hot steam envelope inflates a water world a lot.
ZENG_WATER_BY_TEMPERATURE: dict[int, str] = {
    300: "massradius_100percentH2O_300K_1mbar.txt",
    500: "massradius_100percentH2O_500K_1mbar.txt",
    700: "massradius_100percentH2O_700K_1mbar.txt",
    1000: "massradius_100percentH2O_1000K_1mbar.txt",
}


# ---------------------------------------------------------------------------
# Sample definitions
# ---------------------------------------------------------------------------

# --- Radius valley, in the spirit of Fulton et al. (2017, AJ 154, 109) ------
#: Close-in planets only: the valley is a feature of strongly irradiated orbits.
VALLEY_MAX_PERIOD_DAYS = 100.0
#: FGK hosts, the range Fulton et al. used. M dwarfs are left out because their
#: valley sits at smaller radii and would blur the combined histogram.
VALLEY_TEFF_RANGE_K = (4700.0, 6500.0)
#: Main-sequence hosts only. Evolved stars have larger, less certain radii,
#: and every planet radius inherits its star's.
VALLEY_MIN_LOGG = 3.9
#: The headline sample keeps radii known to better than 10%.
VALLEY_MAX_RADIUS_RELERR = 0.10

#: Where to look for the valley's minimum. Stated as an input, like the period
#: range of a transit search: generous margins around the published 1.5-2.0 R_e.
VALLEY_SEARCH_RANGE_REARTH = (1.4, 2.4)
#: Peaks are searched for between the valley and these limits.
VALLEY_PEAK_RANGE_REARTH = (1.0, 4.0)
#: Gaussian kernel width, in dex of radius (0.03 dex = 7%). The valley location
#: is stable to +-0.03 R_e for bandwidths 0.02-0.04; wider kernels smear it.
VALLEY_KDE_BANDWIDTH_DEX = 0.03
#: Periods are measured relative to this pivot when fitting the valley's tilt,
#: which makes the fitted radius the valley's location at a 10-day orbit.
VALLEY_PIVOT_PERIOD_DAYS = 10.0
#: Trial slopes d(log R)/d(log P) for the tilt fit.
VALLEY_SLOPE_GRID = (-0.30, 0.20, 0.01)

#: Van Eylen et al. (2018, MNRAS 479, 4786), support-vector-machine fit to 117
#: planets with asteroseismic host radii: log10 R = m log10 P + a.
VAN_EYLEN_SLOPE = -0.09
VAN_EYLEN_SLOPE_ERR = (0.02, 0.04)  # (+, -)
VAN_EYLEN_INTERCEPT = 0.37

# --- Mass-radius --------------------------------------------------------------
#: Masses to 25% and radii to 10%: loose enough to keep a few hundred small
#: planets, tight enough that a point's position means something.
MR_MAX_MASS_RELERR = 0.25
MR_MAX_RADIUS_RELERR = 0.10
#: A planet counts as "too big for rock" (or water) when it sits more than this
#: many standard errors above the model curve.
MR_EXCESS_SIGMA = 2.0
#: Core mass used to illustrate what a thin hydrogen envelope does to a radius.
MR_ILLUSTRATIVE_CORE_MASS = 5.0

# --- Habitable zone -------------------------------------------------------------
#: Kopparapu et al. (2013, 2014) fits are valid over this stellar Teff range.
#: Cooler hosts are evaluated at the lower bound (TRAPPIST-1, at 2566 K, is the
#: only famous case, and sits 34 K outside).
HZ_TEFF_VALID_K = (2600.0, 7200.0)
#: Radius below which a planet is plausibly rocky. Rogers (2015, ApJ 801, 41)
#: found most planets larger than ~1.6 R_e are not.
ROCKY_RADIUS_MAX = 1.6
#: Host stars cooler than this are counted as M dwarfs.
M_DWARF_TEFF_MAX = 3900.0

# --- Reproducibility --------------------------------------------------------------
N_BOOTSTRAP = 1000
RANDOM_SEED = 20170901


# ---------------------------------------------------------------------------
# Discovery methods
# ---------------------------------------------------------------------------

#: Archive method name -> display group. The archive distinguishes eleven
#: methods; the seven rare ones (80 planets between them) share "Other".
METHOD_GROUPS = {
    "Transit": "Transit",
    "Radial Velocity": "Radial velocity",
    "Microlensing": "Microlensing",
    "Imaging": "Direct imaging",
}
OTHER_METHOD = "Other"
#: Display order, used for stacking and for colour assignment.
METHOD_ORDER = ["Transit", "Radial velocity", "Microlensing", "Direct imaging", OTHER_METHOD]


# ---------------------------------------------------------------------------
# Named planets
# ---------------------------------------------------------------------------

#: Labelled on the habitable-zone figure. Kepler-186 f and Kepler-452 b, both
#: famous as "Earth cousins", are absent on purpose: the archive now flags them
#: as controversial, and controversial planets are excluded from every sample
#: except the discovery census.
HZ_HIGHLIGHT = [
    "TRAPPIST-1 e", "TRAPPIST-1 f", "TRAPPIST-1 g", "Proxima Cen b",
    "TOI-700 d", "LHS 1140 b", "K2-18 b", "Teegarden's Star b",
    "Kepler-442 b", "Kepler-62 f",
]
CONTROVERSIAL_FAMOUS = ["Kepler-186 f", "Kepler-452 b"]

#: Labelled on the mass-radius figure. Kepler-10 b is the rocky planet this
#: repository detected for itself in project A1.
MR_HIGHLIGHT = ["Kepler-10 b", "55 Cnc e", "K2-18 b", "GJ 1214 b", "TRAPPIST-1 e", "LHS 1140 b"]
