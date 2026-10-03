"""Project configuration: paths, data sources, physical constants and sample cuts.

Everything tunable lives here, so a change of cut or constant propagates to
every measurement and figure. The physical constants are the handful of numbers
the whole project rests on; each carries its source.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Data sources
# --------------------------------------------------------------------------

#: The Minor Planet Center's orbit file: osculating elements for every minor
#: planet with an orbit, numbered or not. Regenerated daily.
MPCORB_URL = "https://minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz"
MPCORB_PATH = DATA_DIR / "MPCORB.DAT.gz"
#: Parsed columns, cached so the 300 MB text file is only parsed once.
MPCORB_PARSED = DATA_DIR / "mpcorb_parsed.npz"

#: AstDyS (University of Pisa / SpaceDyS): synthetic proper elements and the
#: family classification built on them.
ASTDYS_BASE = "https://newton.spacedys.com/~astdys2/propsynth/"
ASTDYS_FILES = {
    # Proper a, e, sin i, plus the proper frequencies g and s, for the main
    # belt and Hungarias (numbered and multi-opposition asteroids).
    "all.syn": "proper elements, main belt + Hungarias",
    # The ~26,000 asteroids locked in secular resonances are kept out of
    # all.syn, because their proper eccentricity is not defined; this file
    # has them, with e replaced by the libration amplitude. Their a and sin i
    # are ordinary proper elements, so they are merged back for those two.
    "secres.syn": "proper elements, secular-resonant asteroids",
    # Family membership: only asteroids that belong to some family are listed.
    "all_tro.members": "family membership",
    # One row per family: size, velocity cutoff, extent in proper elements.
    "all_tro.famtab": "family summary table",
    # Family names, as the HTML option list the AstDyS plot tool uses.
    "family_list.txt": "family names",
}
ASTDYS_PARSED = DATA_DIR / "astdys_parsed.npz"

#: JPL Horizons, for where the planets are at the snapshot instant.
HORIZONS_PATH = DATA_DIR / "horizons_planets.json"
#: Horizons ids: system barycentres, so each "planet" carries its moons' mass.
HORIZONS_BODIES = {
    "Mercury": "1",
    "Venus": "2",
    "Earth": "3",
    "Mars": "4",
    "Jupiter": "5",
    "Saturn": "6",
}

PROVENANCE_PATH = DATA_DIR / "provenance.json"

# --------------------------------------------------------------------------
# Physical constants
# --------------------------------------------------------------------------

#: Gaussian gravitational constant, AU^(3/2) day^-1 Msun^(-1/2). With it,
#: Kepler's third law for a massless body reads n = k / a^(3/2) in rad/day.
GAUSS_K = 0.01720209895

#: Jupiter's mean sidereal orbital period, days (NASA Jupiter fact sheet; the
#: long-term mean, not the osculating value, which wobbles by ~0.1% as Saturn
#: tugs on Jupiter). Every resonance location follows from this one number.
JUPITER_PERIOD_DAYS = 4332.589
#: Saturn's mean sidereal orbital period, days (NASA Saturn fact sheet). Needed
#: only to place the three-body (Jupiter-Saturn-asteroid) resonances.
SATURN_PERIOD_DAYS = 10759.22

#: Planet/Sun mass ratios for the system barycentres, from the DE440 GM values
#: (Park et al. 2021): GM_sun = 132712440041.279 km^3 s^-2.
MASS_RATIO = {
    "Jupiter": 126712764.1 / 132712440041.279,  # = 1/1047.349
    "Saturn": 37940584.8 / 132712440041.279,  # = 1/3497.90
}

#: Planet-crossing limits (NASA fact sheets). An asteroid whose perihelion
#: falls inside Mars' aphelion can meet Mars; one whose aphelion reaches past
#: Jupiter's perihelion can meet Jupiter.
MARS_APHELION_AU = 1.6662
EARTH_APHELION_AU = 1.0167
JUPITER_PERIHELION_AU = 4.9506

#: Saturn's secular frequencies g6 and s6, arcsec/yr (Laskar et al. 2004,
#: La2004 solution). An asteroid whose perihelion precesses at g6 is in the
#: nu6 secular resonance -- the belt's inner boundary.
G6_ARCSEC_PER_YR = 28.2450
S6_ARCSEC_PER_YR = -26.347855

# --------------------------------------------------------------------------
# Samples and cuts
# --------------------------------------------------------------------------

#: The MPC's orbit-quality "uncertainty parameter" U runs 0 (excellent) to 9
#: (barely constrained). For measuring narrow features we keep U <= this.
MAX_UNCERTAINTY_U = 2

#: Semi-major axis range of the headline histogram, AU.
A_RANGE_FULL = (1.8, 5.5)
#: The main belt proper, AU, for the gap measurement.
A_RANGE_BELT = (2.0, 3.5)

#: Histogram bin width for gap detection, AU. Narrow enough to resolve the
#: weak gaps (~0.01 AU wide), wide enough for ~1,000+ asteroids per bin.
GAP_BIN_AU = 0.002
#: Width of the running-median "background" the gaps are measured against, AU.
#: Must be several times wider than the widest gap (~0.05 AU at 3:1).
GAP_BACKGROUND_AU = 0.12
#: A candidate gap must be this many standard deviations below background.
GAP_MIN_SIGMA = 5.0
#: ...and at least this deep, as a fractional deficit, to be reported.
GAP_MIN_DEPTH = 0.15

#: Population boxes in osculating elements (the MPC's own classification
#: uses similar cuts). Hildas: 3:2 resonance; Trojans: 1:1 resonance.
HILDA_A = (3.70, 4.20)
HILDA_E_MAX = 0.30
HILDA_I_MAX_DEG = 20.0
TROJAN_A = (4.95, 5.45)

#: Absolute magnitude below which the main-belt catalogue is essentially
#: complete (every asteroid that bright has been found), used to show the
#: gaps are not a discovery artefact. Surveys have been complete to H ~ 15
#: across the whole main belt for years; we take a conservative value.
H_COMPLETE = 14.0

#: Geometric albedo assumed when converting H to a diameter for scale bars
#: and captions only. Real main-belt albedos range 0.03-0.5.
ASSUMED_ALBEDO = 0.14

# --------------------------------------------------------------------------
# N-body simulation (extension, a labelled simulation)
# --------------------------------------------------------------------------

SIM_PATH = DATA_DIR / "simulation.npz"
#: Test particles, spread uniformly in semi-major axis over this range, AU.
SIM_A_RANGE = (2.0, 3.4)
SIM_N_PARTICLES = 3000
#: Integration length and number of saved snapshots (every 1,000 years).
SIM_YEARS = 1.0e6
SIM_SNAPSHOTS = 1001
#: WHFast timestep, years: ~1/30 of the innermost test-particle orbit.
SIM_DT_YEARS = 0.1
#: A particle whose perihelion reaches inside Mars' aphelion is counted as
#: removed: from there, close encounters with Mars and Earth (not modelled)
#: clear it out on ~Myr timescales.
SIM_REMOVE_Q_AU = MARS_APHELION_AU
SIM_SEED = 1866  # the year Kirkwood announced the gaps
