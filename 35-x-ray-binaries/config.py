"""Project configuration: paths, data sources, physical constants, cuts.

Every selection criterion and every adopted literature number lives in this
file or in ``literature.py``, once, so that each value printed by ``run.py``
or quoted in the README traces back to a single line.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"  # expensive intermediate results (MCMC chains)
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"  # provenance of every download
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 20261004

# ---------------------------------------------------------------------------
# Physical constants (SI). GM_sun is the IAU 2015 nominal value, so a mass
# function computed from P and K never needs G and M_sun separately.
# ---------------------------------------------------------------------------

GM_SUN = 1.3271244e20  # m^3 s^-2
AU_M = 1.495978707e11
PC_M = 3.0856775814913673e16
DAY_S = 86400.0
YEAR_DAYS = 365.25
KPC_CM = 3.0856775814913673e21
# Eddington luminosity for ionised hydrogen, per solar mass
L_EDD_PER_MSUN = 1.26e38  # erg/s

# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/"

#: BlackCAT as published (Corral-Santana et al. 2016), from the CDS.
BLACKCAT_VIZIER = {
    "catalog": "J/A+A/587/A61",
    "files": ("ReadMe", "tablea1.dat", "tablea2.dat", "tablea4.dat", "refs.dat"),
    "citation": "Corral-Santana J.M. et al. 2016, A&A 587, A61 (BlackCAT)",
}
#: The living online BlackCAT: one page per black hole transient.
BLACKCAT_ONLINE = "https://www.astro.puc.cl/BlackCAT/"
BLACKCAT_POLITE_DELAY_S = 0.5

#: CatNoir (Rocha et al. 2026, arXiv:2609.14154): an independent living
#: catalogue of 42 Galactic black holes with their measurement distributions.
CATNOIR_URL = "https://catnoir.iag.usp.br/js/BH_Catalog.json"

#: Paulo Freire's page of precise neutron star masses from pulsar timing.
FREIRE_URL = "https://www3.mpifr-bonn.mpg.de/staff/pfreire/NS_masses.html"

#: GWOSC cumulative event list (every GWTC release, newest version of each event).
GWOSC_URL = "https://gwosc.org/eventapi/json/GWTC/"

#: MAXI/GSC daily light curves (version 7.7L processing, updated daily).
MAXI_URL = "http://maxi.riken.jp/star_data/{name}/{name}_g_lc_1day_all.dat"
#: Swift/BAT hard X-ray transient monitor, daily averages, 15-50 keV.
BAT_URL = "https://swift.gsfc.nasa.gov/results/transients/{path}.lc.txt"

#: X-ray sources: MAXI name and BAT monitor path for each.
XRAY_SOURCES = {
    "MAXI J1820+070": {"maxi": "J1820+071", "bat": "weak/MAXIJ1820p070"},
    "GX 339-4": {"maxi": "J1702-487", "bat": "weak/GX339-4"},
    "Crab": {"maxi": "J0534+220", "bat": "Crab"},
}

#: Gaia DR3 TAP endpoints, tried in order (ESA, then the ARI Heidelberg mirror).
GAIA_TAP_URLS = (
    "https://gea.esac.esa.int/tap-server/tap/sync",
    "https://gaia.ari.uni-heidelberg.de/tap/sync",
)
#: Gaia DR3 source_ids of the dormant black holes with DR3 orbits.
GAIA_BH_SOURCE_IDS = {
    "Gaia BH1": 4373465352415301632,
    "Gaia BH2": 5870569352746779008,
}

# ---------------------------------------------------------------------------
# Analysis settings
# ---------------------------------------------------------------------------

#: Monte Carlo draws per system when turning (P, K, q, i) into a mass.
N_MASS_SAMPLES = 20000
#: Draws per system for the hierarchical (population) fits, histogrammed on a
#: fine log-mass grid (masses.PopulationLikelihood).
N_POP_SAMPLES = 200_000

#: Rhoades & Ruffini (1974): no neutron star can exceed 3.2 M_sun if general
#: relativity holds and the sound speed stays below c. Kalogera & Baym (1996)
#: give 2.9 M_sun for an equation of state matched to nuclear data.
NS_MAX_CAUSAL = 3.2
NS_MAX_KALOGERA_BAYM = 2.9

#: Population model priors (Farr et al. 2011 / Kreidberg et al. 2012):
#: 0 <= M_min <= M_max <= 40 M_sun and -15 <= alpha <= 13 for the power law.
POP_MASS_MAX = 40.0
POP_ALPHA_RANGE = (-15.0, 13.0)
#: Live points for the nested sampling of each population model (dynesty).
POP_NLIVE = 1000

#: Gravitational-wave events: keep candidates with p_astro >= this.
GW_PASTRO_MIN = 0.5
#: Network SNR threshold for the 1/V_max selection weights.
GW_SNR_THRESHOLD = 8.0

#: Crab normalisation: running median over this window (days).
CRAB_WINDOW_D = 365.0
#: The Crab's spectrum used to turn Crab units into energy flux:
#: N(E) = 9.7 E^-2.1 photons cm^-2 s^-1 keV^-1 (Toor & Seward 1974).
CRAB_NORM = 9.7
CRAB_INDEX = 2.1

#: MAXI J1820+070 for Eddington ratios: VLBI parallax distance (Atri et al.
#: 2020) and the mass for the radio-jet inclination (Torres et al. 2020).
J1820_DISTANCE_KPC = (2.96, 0.33)
J1820_MASS = (8.48, 0.79, 0.72)  # median, +, -
#: 2018 outburst window (MJD) and the literature dates of its states
#: (Homan et al. 2020, ApJL 891, L29: hard state MJD 58189-58290 and
#: 58400-58443; the rapid hard-to-soft transition 2018 July 3-9).
J1820_OUTBURST = (58185.0, 58450.0)
J1820_HOMAN_HARD1 = (58189.0, 58290.0)
J1820_HOMAN_HARD2 = (58400.0, 58443.0)
J1820_HOMAN_TRANSITION = (58302.0, 58308.0)

#: GX 339-4 outburst detection: MAXI 2-20 keV flux (Crab units), smoothed.
GX339_OUTBURST_THRESHOLD_CRAB = 0.02
GX339_MIN_OUTBURST_DAYS = 40
