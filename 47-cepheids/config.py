"""Project configuration: paths, data sources and the conventions we adopt.

Everything we *assume* lives here, each value with its source. The published
results we grade ourselves against (distances, slopes, geometry) are kept
separately in ``literature.py`` and are never used as inputs to a fit, except
where a fit explicitly says it adopts one (e.g. the eclipsing-binary distances
used to *measure* the metallicity term).
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
OGLE_DIR = RAW_DIR / "ogle"
GAIA_DIR = RAW_DIR / "gaia"
LIT_DIR = RAW_DIR / "literature"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, OGLE_DIR, GAIA_DIR, LIT_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: OGLE Collection of Variable Stars (OCVS), classical Cepheids of the
#: Magellanic Clouds: Soszynski et al. 2015 (AcA 65, 297), extended 2017-2023.
OGLE_BASE = "https://www.astrouw.edu.pl/ogle/ogle4/OCVS/"
OGLE_CEP_FILES = (
    "README", "ident.dat", "remarks.txt", "phot.tar.gz",
    "cepF.dat", "cep1O.dat", "cep2O.dat",
    "cepF1O.dat", "cep1O2O.dat", "cep1O3O.dat", "cep2O3O.dat",
    "cepF1O2O.dat", "cep1O2O3O.dat",
)
GALAXIES = ("lmc", "smc")

#: Pietrukowicz et al. (2021, AcA 71, 205) list of all known Galactic classical
#: Cepheids, compiled from many surveys and vetted by the OGLE team; updated
#: March 2026, with Gaia DR3 source_ids.
OGLE_GALACTIC_LIST = "allGalCep.listID"

#: Gaia archive (ESA) synchronous TAP endpoint, and CDS xMatch for positions.
GAIA_TAP_SYNC = "https://gea.esac.esa.int/tap-server/tap/sync"
CDS_XMATCH = "http://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync"
VIZIER_GAIA_DR3 = "vizier:I/355/gaiadr3"
#: Match radius for OGLE positions against Gaia DR3 (both are J2000/ICRS at
#: their survey epochs; Magellanic proper motions are ~2 mas/yr, negligible).
XMATCH_RADIUS_ARCSEC = 2.0

# ---------------------------------------------------------------------------
# Photometric conventions
# ---------------------------------------------------------------------------

#: OGLE's reddening-free Wesenheit index W_I = I - R (V - I) with
#: R = A_I / E(V-I) = 1.55 (Udalski et al. 1999; standard Galactic extinction
#: law). Any star moved by dust slides along a line of slope R in (V-I, I),
#: so W_I is unchanged by extinction if the law is right.
R_WI = 1.55

#: Gaia's Wesenheit W_G = G - 1.90 (BP - RP) (Ripepi et al. 2019, A&A 625, A14;
#: used for all Gaia DR2/DR3 Cepheid work).
R_WG = 1.90

#: Pivot period for the Leavitt law: M = a + b (log P - LOGP_PIVOT). Ten days
#: sits in the middle of the long-period Cepheids used for distances, which
#: decorrelates zero point and slope.
LOGP_PIVOT = 1.0

#: Iterative sigma clipping as in OGLE's own PL fits (outliers are mostly blends
#: and misidentified modes). Clip at this many standard deviations.
CLIP_SIGMA = 3.0

#: First-overtone Cepheids pulsate faster than a fundamental-mode star of the
#: same luminosity. The standard conversion to the "fundamentalised" period
#: (Feast & Catchpole 1997) is P_F = P_1O / (0.720 - 0.027 log P_1O). We use
#: it only for plots that mix the modes; fits treat the modes separately.
def fundamentalise(p1o):
    import numpy as np
    return p1o / (0.720 - 0.027 * np.log10(p1o))


# ---------------------------------------------------------------------------
# Gaia astrometry conventions
# ---------------------------------------------------------------------------

#: Quality cuts for the Milky Way calibrators. RUWE < 1.4 is the standard
#: single-star astrometry cut (Lindegren et al. 2021); G > 6 is the bright end
#: of the Lindegren et al. (2021) parallax zero-point recipe's calibration.
RUWE_MAX = 1.4
G_MIN_L21 = 6.0
G_MAX = 19.0

#: Random seed for every Monte Carlo and MCMC in the project.
SEED = 20261006
