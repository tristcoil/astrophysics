"""Paths, data sources and physical constants shared by every script.

Every number that steers the analysis is defined here once, with where it
comes from, so a reader can see (and change) the choices in one place.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"            # downloads exactly as served
GAIA_DIR = DATA_DIR / "gaia"          # Gaia DR3 query results (+ JSON provenance)
LIT_DIR = DATA_DIR / "literature"     # tables copied from papers
DERIVED_DIR = DATA_DIR / "derived"    # caches of slow intermediate results
FIG_DIR = PROJECT_DIR / "figures"
PROVENANCE = DATA_DIR / "provenance.json"

for _d in (RAW_DIR, GAIA_DIR, LIT_DIR, DERIVED_DIR, FIG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 20261009

# --------------------------------------------------------------------------
# El-Badry, Rix & Heintz (2021, MNRAS 506, 2269): a million binaries from
# Gaia EDR3. Zenodo record 4435257 (files and MD5s as listed by the Zenodo
# API on 9 Oct 2026).
# --------------------------------------------------------------------------
ZENODO_RECORD = 4435257
ZENODO_FILES = {
    "all_columns_catalog.fits.gz": dict(size=1_419_826_533, md5="9644d53f056758fb2ab1d3d04a45e25d"),
    "all_columns_catalog_shift.fits.gz": dict(size=354_073_286, md5="dd34343c5a3947f5db0fbb79b85f4622"),
    "find_binaries_edr3.py": dict(size=14_017, md5="d22d1e72713af4b143dedb161a65844a"),
    "num_neighbors_edr3.py": dict(size=3_616, md5="34d44081ac70e8abfb8f8581408a1ae0"),
}

# --------------------------------------------------------------------------
# Physical constants (SI) and units
# --------------------------------------------------------------------------
G_SI = 6.67430e-11            # CODATA 2018
MSUN_KG = 1.988409870698051e30  # IAU 2015 nominal GM_sun / G
GM_SUN = 1.3271244e20         # m^3 s^-2 (IAU 2015 nominal)
AU_M = 1.495978707e11
PC_M = 3.0856775814913673e16
KMS = 1.0e3
YR_S = 365.25 * 86400.0
#: 1 mas/yr at 1 kpc in km/s (= AU/yr in km/s)
K_PM = 4.740470463533348

# --------------------------------------------------------------------------
# MOND
# --------------------------------------------------------------------------
#: Milgrom's acceleration scale. 1.2e-10 m/s^2 is the conventional value
#: (Begeman, Broeils & Sanders 1991; McGaugh et al. 2016 fit g† = 1.20e-10)
#: and the one Banik et al. (2024) and Chae (2023) adopt. Project D4 measured
#: g† = 1.11e-10 from SPARC with the RAR's own interpolating function.
A0_SI = 1.2e-10

# --------------------------------------------------------------------------
# Gaia DR3 TAP services, tried in order (name, URL, gaia_source table name).
# On 9 Oct 2026 both ESA and ARI answered a level-2 HEALPix count in 2-4 s.
# --------------------------------------------------------------------------
GAIA_SERVICES = [
    ("ESA Gaia archive", "https://gea.esac.esa.int/tap-server/tap", "gaiadr3.gaia_source"),
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap", "gaiadr3.gaia_source"),
]
