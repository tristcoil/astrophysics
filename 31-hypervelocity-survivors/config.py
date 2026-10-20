"""Project configuration: paths, the Galaxy model, every selection number.

L5 asks whether some Type Ia supernovae leave a survivor behind. If one white
dwarf in a very tight pair explodes, its companion is released at the orbital
speed it had at that instant, well over 1,000 km/s, and the fastest stars in
the Galaxy should be such survivors. The project:

1. repeats the Gaia search for them (search.py) and asks how many candidates
   are real and how many are parallax errors;
2. turns the known survivors' astrometry and published radial velocities into
   3D velocities, flight times and launch speeds (kinematics.py);
3. traces them back through a Milky Way potential and tests whether any of
   them passes through a known supernova remnant, with a chance-alignment
   test (orbits.py, remnants.py);
4. inverts the launch speeds into the masses of the binaries that exploded
   (binary.py);
5. compares how many survivors Gaia sees with how many Type Ia supernovae the
   Galaxy makes (population.py).

The literature numbers we compare against live in literature.py; nothing in
this file is an answer the analysis is meant to find.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
GAIA_DIR = DATA_DIR / "gaia"
AUX_DIR = DATA_DIR / "aux"
CACHE_DIR = DATA_DIR / "cache"
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, GAIA_DIR, AUX_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 20261004

# --------------------------------------------------------------------------
# Gaia DR3 access
# --------------------------------------------------------------------------
# ESA's own archive was overloaded on 4 Oct 2026 (a one-row lookup took 24 s,
# a 2-million-row count timed out), and the CDS VizieR TAP returned 503, as it
# had on 30 Sep. Full copies of gaiadr3.gaia_source are also served by the
# Astronomisches Rechen-Institut in Heidelberg and by NOIRLab's Astro Data Lab;
# both ran whole-catalogue scans. Each query records which service answered.

#: (name, TAP URL, table name of gaia_source on that service), tried in order.
GAIA_SERVICES = [
    ("ARI Heidelberg", "https://gaia.ari.uni-heidelberg.de/tap", "gaiadr3.gaia_source"),
    ("NOIRLab Astro Data Lab", "https://datalab.noirlab.edu/tap", "gaia_dr3.gaia_source"),
    ("ESA Gaia archive", "https://gea.esac.esa.int/tap-server/tap", "gaiadr3.gaia_source"),
]

#: 1 mas/yr at 1 kpc, in km/s (1 AU per year).
K_KMS = 4.740470463533348

# --------------------------------------------------------------------------
# The Sun's place in the Galaxy (astropy's v4.0 Galactocentric defaults)
# --------------------------------------------------------------------------

#: Sun - Galactic Centre distance (GRAVITY Collaboration 2018).
R0_KPC = 8.122
#: Sun's height above the midplane (Bennett & Bovy 2019).
Z_SUN_KPC = 0.0208
#: Sun's Galactocentric velocity (U toward the GC, V with rotation, W to the
#: North Galactic Pole), km/s (Drimmel & Poggio 2018).
V_SUN_KMS = (12.9, 245.6, 7.78)

# --------------------------------------------------------------------------
# The search (search.py)
# --------------------------------------------------------------------------

#: El-Badry et al. (2023), section 2.1.2, copied verbatim (only the table name
#: is substituted per service). Their paper lists the 25 sources it returns.
ELBADRY_QUERY = """select {cols} from {table}
where phot_g_mean_mag < 20
and ruwe < 1.4
and pm > 50
and (4.74*pm/(parallax + parallax_error) > 400 or (parallax + parallax_error) < 0)
and (4.74*pm/parallax > 600 or parallax < 0)
and bp_rp < 0.5
and parallax_over_error < 5"""

#: Our broad search: every source whose tangential velocity, taken at face
#: value, exceeds V_TAN_MIN, with a parallax at least PLX_SNR_MIN sigma, at any
#: colour and without astrometric quality cuts (they are applied, and
#: studied, afterwards). Its mirror image -- parallaxes that are significantly
#: *negative*, which can only be measurement errors -- estimates how many of
#: the positive ones are errors too.
V_TAN_MIN = 600.0
PLX_SNR_MIN = 3.0
PM_MIN = 20.0
G_MAX = 20.5

#: Cuts studied in search.py (each one alone, then together).
RUWE_MAX = 1.4
FIDELITY_MIN = 0.5          # Rybizki et al. 2022 classifier, where available
IPD_MULTI_PEAK_MAX = 5      # per cent of transits with a double image
IPD_HARMONIC_MAX = 0.1      # image-shape asymmetry across scan angles
VIS_PERIODS_MIN = 10

#: Background colour-magnitude diagram: stars within 100 pc with good
#: parallaxes (the same sample project B1 draws its HR diagram from).
BACKGROUND_PLX_MIN = 10.0
BACKGROUND_PLX_SNR = 10.0

#: Random all-sky subsample used to measure Gaia's parallax precision as a
#: function of magnitude (for the population model's selection function).
#: random_index runs from 0 to ~1.8e9 over the whole catalogue.
PRECISION_RANDOM_INDEX_MAX = 400_000

# --------------------------------------------------------------------------
# Kinematics (kinematics.py)
# --------------------------------------------------------------------------

#: Posterior samples per star.
N_SAMPLES = 4000
#: Ejection speeds above this are treated as impossible (El-Badry et al.
#: 2023, eq. 4: roughly two 1.3 Msun white dwarfs at contact).
V_EJECTION_MAX = 3000.0
#: Distance grid for the posteriors, kpc.
D_GRID_MIN, D_GRID_MAX, D_GRID_N = 0.02, 40.0, 20000
#: A_G / E(B-V), as used by El-Badry et al. (2023) with the Schlafly &
#: Finkbeiner (2011) recalibration of the SFD dust map.
A_G_PER_EBV = 2.67
#: Flat rotation speed used by El-Badry et al. to remove the Galactic
#: rotation at the launch point. Our own runs use the potential's circular
#: speed instead; this one only for the like-for-like comparison.
V_CIRC_FLAT = 240.0

# --------------------------------------------------------------------------
# Orbits and remnants (orbits.py, remnants.py)
# --------------------------------------------------------------------------

#: How far back to trace the D6 stars when looking for remnants: radiative
#: remnants fade into the interstellar medium within a few 1e5 yr.
REMNANT_LOOKBACK_YR = 3.0e5
#: Trajectory samples per star used for the remnant search.
N_TRACK_SAMPLES = 600
#: Random re-placements of the remnant catalogue for the chance-alignment
#: test (shifted in Galactic longitude, latitudes kept).
N_SHUFFLES = 20000

# --------------------------------------------------------------------------
# Population model (population.py)
# --------------------------------------------------------------------------

#: Milky Way Type Ia rate (Li et al. 2011, LOSS: 0.54 +- 0.12 per century).
SN_IA_RATE_PER_YR = 0.0054
SN_IA_RATE_ERR = 0.0012
#: Stellar disc used to place the explosions (Bland-Hawthorn & Gerhard 2016),
#: as in Shen et al. (2018) and El-Badry et al. (2023).
DISC_SCALE_LENGTH_KPC = 2.6
DISC_SCALE_HEIGHT_KPC = 0.3
#: Dust layer for extinction in the model: exponential, A_V = DUST_AV_PER_KPC
#: per kpc in the midplane at the solar circle.
DUST_SCALE_HEIGHT_KPC = 0.12
DUST_AV_PER_KPC = 1.0
A_G_PER_A_V = 0.85
N_POP_STARS = 400_000
