"""Project configuration: paths, the catalogue version, and every tunable choice.

The whole project runs on one file: the ATNF Pulsar Catalogue, the standard
compilation of every published pulsar. The version is pinned so that the numbers
quoted in the README are reproducible; bump ``CATALOGUE_VERSION`` and re-run
``fetch_data.py`` to move to a newer release.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"

for _d in (DATA_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --- The catalogue --------------------------------------------------------------

#: ATNF Pulsar Catalogue release. v2.8.1 was released 2026-06-18 with 4,393
#: pulsars; the release history is at .../psrcat/catalogueHistory.html.
CATALOGUE_VERSION = "2.8.1"
CATALOGUE_URL = (
    "https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/"
    f"psrcat_pkg.v{CATALOGUE_VERSION}.tar.gz"
)
TARBALL_PATH = DATA_DIR / f"psrcat_pkg.v{CATALOGUE_VERSION}.tar.gz"
#: The plain-text database inside the tarball, extracted so it can be read directly.
CATALOGUE_DB = DATA_DIR / f"psrcat_v{CATALOGUE_VERSION}.db"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

# --- Sorting pulsars into families ------------------------------------------------
#
# The catalogue's TYPE field flags the non-standard objects (magnetars, X-ray
# isolated neutron stars, radio-quiet pulsars, rotating radio transients). Anything
# not flagged is an ordinary radio pulsar -- or a millisecond pulsar, which the
# catalogue does not flag, so we define it here.

#: A pulsar is "millisecond" (recycled) if it spins faster than this...
MSP_MAX_PERIOD_S = 0.030
#: ...and, where Pdot is measured, slows down more gently than this. The second
#: condition keeps out the two young pulsars that also spin faster than 30 ms
#: (J0537-6910 and J0058-7218, both in the Magellanic Clouds) but have the huge
#: Pdot of youth. It is not applied inside globular clusters, where Pdot measures
#: the cluster's gravity rather than the pulsar's youth.
MSP_MAX_PDOT = 1e-17

#: DM-based distances at or above this value mean the pulsar's dispersion measure
#: exceeds everything the electron-density model puts along that sightline, so
#: the distance is a placeholder, not a measurement. Excluded from heights.
DISTANCE_CAP_KPC = 25.0

#: Binary fraction is measured in these spin-period bins (seconds).
#: (Above 3 s the bins are merged: only three such pulsars have a timing solution
#: past 10 s, too few for a fraction to mean anything.)
BINARY_FRACTION_BIN_EDGES_S = (0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 100.0)

#: Characteristic-age bins (years) for the "born in the plane" height analysis.
AGE_BIN_EDGES_YR = (1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e10)

# --- The historical check on characteristic ages ----------------------------------

#: The Crab pulsar was born in the supernova Chinese astronomers recorded on
#: 1054 July 4. That date, not any model, fixes its true age.
CRAB_JNAME = "J0534+2200"
CRAB_BIRTH_YEAR = 1054.51

# --- Pulsars worth labelling on the diagram ---------------------------------------
#
# Each earns its place by illustrating something on the plot. Offsets are in
# points, chosen by hand to keep labels off the data; ``ha`` is the side of the
# label the leader line attaches to.

FAMOUS = {
    "J0534+2200": dict(label="Crab", note="born in SN 1054", offset=(-50, 16), ha="right"),
    "J0835-4510": dict(label="Vela", note="", offset=(-22, -16), ha="right"),
    "J1915+1606": dict(label="Hulse–Taylor binary", note="Nobel Prize 1993",
                       offset=(-10, 40), ha="right"),
    "J0737-3039A": dict(label="Double Pulsar A", note="", offset=(-10, 26), ha="right"),
    "J0737-3039B": dict(label="Double Pulsar B", note="same binary, not recycled",
                        offset=(26, 20), ha="left"),
    "J1808-2024": dict(label="SGR 1806−20", note="strongest field", offset=(12, 6),
                       ha="left"),
    "J1935+2154": dict(label="SGR 1935+2154", note="source of a fast radio burst",
                       offset=(-34, 40), ha="right"),
    "J0901-4046": dict(label="J0901−4046", note="76 s: the slowest", offset=(-4, -30),
                       ha="center"),
    "J2144-3933": dict(label="J2144−3933", note="8.5 s, beyond the valley",
                       offset=(14, -26), ha="left"),
}
