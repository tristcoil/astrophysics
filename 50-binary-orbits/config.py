"""Project configuration: paths, data sources, physical constants and the cuts we adopt.

Every number the project *assumes* lives here (or, for published values we
adopt or grade ourselves against, in ``literature.py``), each with its source,
so that any figure or printed result traces back to one line.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
ORB6_DIR = RAW_DIR / "orb6"
LIT_DIR = RAW_DIR / "literature"
GAIA_DIR = RAW_DIR / "gaia"
HIP_DIR = RAW_DIR / "hipparcos"
MONTREAL_DIR = RAW_DIR / "montreal"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, ORB6_DIR, LIT_DIR, GAIA_DIR, HIP_DIR, MONTREAL_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: Sixth Catalog of Orbits of Visual Binary Stars (Hartkopf, Mason & Worley;
#: now maintained at Georgia State University with the WDS). The https address
#: works again (October 2026); the atlas spec had found only http working.
ORB6_URL = "https://www.astro.gsu.edu/wds/orb6/orb6orbits.txt"
ORB6_FORMAT_URL = "https://www.astro.gsu.edu/wds/orb6/orb6format.txt"

CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/"

#: VizieR catalogues fetched whole from the CDS FTP mirror: key -> (catalogue, files).
CDS_CATALOGUES = {
    # Bond et al. 2017 (ApJ 840, 70): 2,367 measures of Sirius B relative to A, 1862-2016.
    "bond2017": ("J/ApJ/840/70", ("ReadMe", "sirius.dat", "notes.dat", "refs.dat")),
    # Eker et al. 2018 (MNRAS 479, 5491): 509 (+77 rejected) main-sequence stars in
    # detached eclipsing binaries, with masses, radii and temperatures.
    "eker2018": ("J/MNRAS/479/5491", ("ReadMe", "table1.dat")),
    # Hipparcos, the new reduction (van Leeuwen 2007).
    "hip2": ("I/311", ("ReadMe", "hip2.dat.gz")),
}

#: arXiv e-prints whose TeX tables we read (only the .tex files are kept).
ARXIV_SOURCES = {
    "bond2017_src": "1703.10625",       # Sirius: HST + historical astrometry, masses (Tables 2-7)
    "bond2015_procyon_src": "1510.00485",  # Procyon: HST astrometry, masses
    "bond2017_40eri_src": "1709.00478",    # 40 Eri B, and Table 2: WDs in visual binaries
    "joyce2018_src": "1809.01240",       # Sirius B gravitational redshift
    "parsons2017_src": "1706.05016",     # eclipsing white dwarf masses and radii
    "eker2018_src": "1807.02568",        # the eclipsing-binary mass-luminosity relation
}

#: Mamajek's "Modern Mean Dwarf Stellar Color and Effective Temperature Sequence"
#: (Pecaut & Mamajek 2013, ApJS 208, 9; version 2022.04.16).
MAMAJEK_URL = "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"

#: Montreal white dwarf models (Bedard et al. 2020, ApJ 901, 93): synthetic
#: photometry of pure-H atmospheres (Table_DA) and the evolutionary sequences
#: with thick hydrogen layers (q_H = 1e-4), C/O cores.
MONTREAL_BASE = "https://www.astro.umontreal.ca/~bergeron/CoolingModels/"
MONTREAL_TABLE_DA = "Tables/Table_DA"
MONTREAL_MASSES = (20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95,
                   100, 105, 110, 115, 120, 125, 130)  # in units of 0.01 Msun

#: Gaia DR3 TAP services: ESA first, then the ARI Heidelberg full mirror.
GAIA_TAP_ESA = "https://gea.esac.esa.int/tap-server/tap"
GAIA_TAP_ARI = "https://gaia.ari.uni-heidelberg.de/tap"
VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"

#: Sirius B in Gaia DR3 (G = 8.52, 10.7" from Sirius A in 2016; Sirius A itself
#: is too bright to be in Gaia). Found by a cone search around Sirius.
SIRIUS_B_GAIA_DR3 = 2947050466531873024
SIRIUS_HIP = 32349

# ---------------------------------------------------------------------------
# Physical constants (IAU 2012/2015 nominal values)
# ---------------------------------------------------------------------------

GM_SUN = 1.3271244e20       # m^3 s^-2 (IAU 2015 nominal)
G_NEWTON = 6.67430e-11      # m^3 kg^-1 s^-2 (CODATA 2018)
AU = 1.495978707e11         # m (IAU 2012, exact)
PC = AU * 648000.0 / 3.141592653589793
DAY = 86400.0
JULIAN_YEAR = 365.25        # days; every period in the project is in these years
C_LIGHT = 299792458.0       # m/s
R_SUN = 6.957e8             # m (IAU 2015 nominal)
L_SUN = 3.828e26            # W (IAU 2015 nominal)
M_SUN = GM_SUN / G_NEWTON   # kg
SIGMA_SB = 5.670374419e-8   # W m^-2 K^-4
TEFF_SUN = 5772.0           # K (IAU 2015 nominal)
MBOL_SUN = 4.74             # IAU 2015 B2 zero point
MV_SUN = 4.81               # Willmer 2018 (ApJS 236, 47), Vega system
AU_PER_YR_KMS = AU / (JULIAN_YEAR * DAY) / 1e3   # 4.7405 km/s

#: Kepler's third law in the units used by every visual-binary paper:
#: M[Msun] = KEPLER * a[au]^3 / P[yr]^2 with P in Julian years. KEPLER is
#: 1.00003 (it would be exactly 1 for sidereal years); kept for honesty.
KEPLER = 4.0 * 3.141592653589793**2 * AU**3 / (GM_SUN * (JULIAN_YEAR * DAY) ** 2)

# ---------------------------------------------------------------------------
# Sirius orbit fit
# ---------------------------------------------------------------------------

#: Bond et al. 2017 rejection codes we honour in our own fit: R2 (a false image
#: measured by Searle in 1866-67) and R3 (Schroeder et al. 2000 HST measures,
#: superseded by Bond's re-reduction of the same images). R1 rejections were
#: their 4-sigma clipping; our fit keeps those points and lets an outlier
#: mixture model weigh them instead.
SIRIUS_EXCLUDE_NOTES = ("R2", "R3")

#: Reference aperture for the micrometer error model, sigma ~ (D / D_REF)^-beta.
MICROMETER_D_REF = 0.5      # m

#: Broad component of the outlier mixture, in units of each group's own sigma.
OUTLIER_WIDTH = 8.0

#: emcee settings for the Sirius posterior (seeded; see emcee seeding trap).
MCMC_WALKERS = 64
MCMC_STEPS = 12000
MCMC_BURN = 4000
MCMC_THIN = 10
SEED = 20261009

# ---------------------------------------------------------------------------
# ORB6 mass-luminosity sample
# ---------------------------------------------------------------------------

#: Orbit grades used for dynamical masses (1 = definitive, 2 = good).
MLR_GRADES = (1, 2)
#: Gaia DR3 quality cuts for a "good" parallax.
GAIA_RUWE_MAX = 1.4
GAIA_PLX_SNR_MIN = 20.0     # 5% parallax -> 15% mass
#: Search radius around each WDS position (J2000) for Gaia DR3 sources (2016.0).
#: Proper motions of nearby binaries move them by up to ~1'/16 yr.
GAIA_SEARCH_RADIUS_ARCSEC = 60.0
#: After propagating Gaia positions back to J2000, a component must lie within
#: this distance of the WDS primary position (plus the pair's separation).
GAIA_MATCH_TOL_ARCSEC = 3.0
#: Distance limit for the mass-luminosity fit: inside ~100 pc interstellar
#: extinction is negligible (the Local Bubble), so V magnitudes need no correction.
MLR_MAX_DIST_PC = 100.0
#: Mass precision required for the mass-luminosity fit.
MLR_MAX_MASS_FRAC_ERR = 0.20
#: Systems this far (mag) above the main sequence in the colour-magnitude
#: diagram are treated as evolved (an equal-mass pair sits 0.75 mag above).
EVOLVED_DMAG = 1.5

# ---------------------------------------------------------------------------
# Famous-orbit gallery (WDS id, ORB6 discoverer designation, display name)
# ---------------------------------------------------------------------------

GALLERY = (
    ("06451-1643", "AGC   1AB", "Sirius"),
    ("07393+0514", "SHB   1AB", "Procyon"),
    ("14396-6050", "RHD   1AB", "α Centauri"),
    ("11182+3132", "STF1523AB", "ξ Ursae Majoris"),
    ("07346+3153", "STF1110AB", "Castor"),
    ("21069+3845", "STF2758AB", "61 Cygni"),
)

#: Gaia (E)DR3 commanded scanning law (570 MB) and DPAC's astrometric data gaps.
GAIA_SCANLAW_URL = ("https://cdn.gea.esac.esa.int/Gaia/gedr3/auxiliary/commanded_scan_law/"
                    "CommandedScanLaw_001.csv.gz")
GAIA_GAPS_URL = ("https://www.cosmos.esa.int/documents/29201/1658422/Astrometry_EDR3.csv/"
                 "7ec45023-6fae-8685-3d55-36e47ad68367?t=1632219024144&download=true")
#: CDS cross-match service (positions against any VizieR table).
CDS_XMATCH = "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
