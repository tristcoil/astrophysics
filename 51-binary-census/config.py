"""Project configuration: paths, data sources and the conventions we adopt.

Everything this project *assumes* lives here, each value with its source. The
published results we grade ourselves against (Raghavan et al. 2010's period
distribution, Moe & Di Stefano 2017's companion frequencies, circularization
periods...) are kept separately in ``literature.py`` and are never inputs to a
measurement, except where a step says explicitly that it adopts one.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
GAIA_DIR = RAW_DIR / "gaia"
CHUNK_DIR = GAIA_DIR / "chunks"
VIZIER_DIR = RAW_DIR / "vizier"
LIT_DIR = RAW_DIR / "literature"
GAIAMOCK_DIR = RAW_DIR / "gaiamock"
BUILD_DIR = DATA_DIR / "build"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, GAIA_DIR, CHUNK_DIR, VIZIER_DIR, LIT_DIR, GAIAMOCK_DIR,
           BUILD_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

#: Random seed for every Monte Carlo in the project.
SEED = 20261009

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: ESA Gaia archive, synchronous TAP. Big queries are split by source_id into
#: HEALPix pixels (source_id // 2**35 is the level-12 nested pixel), starting at
#: level 1 (48 pixels) and splitting a pixel in four if it times out.
GAIA_SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
HEALPIX_BASE_LEVEL = 1
HEALPIX_MAX_LEVEL = 4

#: VizieR's ASU interface (its TAP service returned HTTP 503 on 2026-10-09, as
#: it has on several days before; the ASU CGI kept working).
VIZIER_ASU = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
#: CDS xMatch: positional cross-matches against any VizieR table.
CDS_XMATCH = "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"

#: SB9, the 9th Catalogue of Spectroscopic Binary Orbits (Pourbaix et al. 2004,
#: A&A 424, 727; maintained at sb9.astro.ulb.ac.be, which timed out on
#: 2026-10-09). VizieR's copy: 4,079 systems, 5,099 orbits.
SB9_TABLES = {"main": "B/sb9/main", "orbits": "B/sb9/orbits", "alias": "B/sb9/alias"}

#: Raghavan et al. (2010, ApJS 190, 1): 454 Sun-like stars within 25 pc. Table
#: 13 lists every companion with its period or separation; table 18 the
#: companion masses.
RAGHAVAN_TABLES = {"stars": "J/ApJS/190/1/stars", "companions": "J/ApJS/190/1/table13",
                   "companion_masses": "J/ApJS/190/1/table18"}

#: El-Badry, Rix & Heintz (2021, MNRAS 506, 2269): 1.8 million candidate wide
#: pairs within 1 kpc from Gaia EDR3 (source_ids identical in DR3), with a
#: chance-alignment probability R for each.
ELBADRY21_TABLE = "J/MNRAS/506/2269/catalog"

#: Hunt & Reffert (2023, A&A 673, A114): open clusters with ages, and members
#: keyed by Gaia DR3 source_id.
HUNT23_CLUSTERS = "J/A+A/673/A114/clusters"
HUNT23_MEMBERS = "J/A+A/673/A114/members"

#: Eric Mamajek's mean dwarf sequence (colours, absolute magnitudes, masses by
#: spectral type; Pecaut & Mamajek 2013, ApJS 208, 9, updated online).
MAMAJEK_URL = "https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt"

#: DEBCat, John Southworth's catalogue of detached eclipsing binaries with
#: masses and radii measured to 2% or better (Southworth 2015, ASPC 496, 164).
DEBCAT_URL = "https://www.astro.keele.ac.uk/jkt/debcat/debs.dat"

#: gaiamock (El-Badry et al. 2024, OJAp 7, 100): simulates Gaia's epoch
#: astrometry with the real scanning law and runs DR3's astrometric
#: non-single-star cascade. Pinned to one commit; its scan-law files (one per
#: HEALPix level-6 pixel, 49,152 files, 987 MB unzipped) come from the author's
#: Box folder linked in the README.
GAIAMOCK_REPO = "https://raw.githubusercontent.com/kareemelbadry/gaiamock"
GAIAMOCK_COMMIT = "dd30fdbf787eb96878734605ac077ac69bf28c84"
GAIAMOCK_FILES = ("gaiamock.py", "kepler_solve_astrometry.c", "LICENSE", "README.md")
GAIAMOCK_SCANS_URL = ("https://caltech.app.box.com/index.php?rm=box_download_shared_file"
                      "&shared_name=4f7q6qdh0bku881bzvzxc4cm5u0902cf&file_id=f_1611600757604")
GAIAMOCK_SCANS_BYTES = 286_007_203
GAIAMOCK_SCANS_SHA256 = "86eff185a2e1e135157b69703d285216bacd8a50119e5a1f0e8ccf940ff1dc21"
#: gaiamock's C kernel needs the GNU Scientific Library; we build it locally
#: (static, position-independent) because this machine has no system GSL.
GSL_URL = "https://ftp.gnu.org/gnu/gsl/gsl-2.8.tar.gz"
GSL_SHA256 = "6a99eeed15632c6354895b1dd542ed5a855c0f15d9ad1326c6fe2b2c9e423190"

# ---------------------------------------------------------------------------
# Gaia DR3 facts used by the models
# ---------------------------------------------------------------------------

#: DR3 covers 2014-07-25 to 2017-05-28 (34 months; JD 2456863.9-2457901.9).
GAIA_DR3_START_JD = 2456863.9375
GAIA_DR3_END_JD = 2457901.8125
GAIA_DR3_BASELINE_DAYS = GAIA_DR3_END_JD - GAIA_DR3_START_JD

# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

#: The 100 pc sample: parallax > 10 mas.
PARALLAX_MIN_100PC = 10.0
#: A wider volume for wide companions of A and B stars, which are rare nearby.
PARALLAX_MIN_BRIGHT = 2.0
#: ...limited to stars brighter than this absolute G (late-A and hotter).
MG_MAX_BRIGHT = 2.0

#: Wide pairs: El-Badry et al. recommend R < 0.1 for a sample that is >90%
#: pure; R is the ratio of the local density of chance alignments (from their
#: shifted catalogue) to the density of all pairs at that separation.
R_CHANCE_MAX = 0.1
