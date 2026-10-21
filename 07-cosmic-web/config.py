"""Project configuration: paths, the SkyServer query, and every sample cut.

The whole project is one catalogue -- the SDSS Main Galaxy Sample -- cut in
different ways for different questions. Every cut lives here, next to the
reason for it, so any figure can be traced back to exactly which galaxies it
shows and why those and not others.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
#: One CSV per right-ascension chunk, exactly as SkyServer served it.
RAW_DIR = DATA_DIR / "skyserver"
#: The chunks merged into one typed catalogue.
CATALOG_PATH = DATA_DIR / "mgs_dr17.npz"
#: Abell rich-cluster catalogue, used only to put names on clusters we find.
ABELL_PATH = DATA_DIR / "abell_clusters.csv"
#: Pair counts are the expensive step; they are cached here keyed by settings.
PAIRCOUNT_DIR = DATA_DIR / "paircounts"
RESULTS_PATH = DATA_DIR / "results.json"
FIGURE_DIR = PROJECT_DIR / "figures"

for _d in (DATA_DIR, RAW_DIR, PAIRCOUNT_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Where the data comes from
# ---------------------------------------------------------------------------

#: SkyServer's public SQL endpoint for Data Release 17, the final release of
#: SDSS-IV and the last to reprocess the original (2000-2008) galaxy survey.
SKYSERVER_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
DATA_RELEASE = "DR17"

#: The survey is fetched in right-ascension strips so that no single query
#: returns more than ~50,000 rows -- well inside SkyServer's per-query limits,
#: and resumable if the connection drops halfway.
RA_CHUNK_DEG = 15.0

#: Main Galaxy Sample target bits in ``legacy_target1``: GALAXY (0x40),
#: GALAXY_BIG (0x80) and GALAXY_BRIGHT_CORE (0x100). Checked against
#: SkyServer's own ``dbo.fPrimTarget()`` function. Selecting on the *target*
#: bits, rather than on "anything classified as a galaxy", matters: SDSS also
#: observed luminous red galaxies, quasar candidates and special programmes, each
#: with a different selection. Mixing them would make the density of the map
#: depend on which programme happened to observe where.
MGS_TARGET_BITS = 0x40 | 0x80 | 0x100

#: The Main Galaxy Sample's flux limit: Petrosian r magnitude, corrected for
#: Milky Way dust, brighter than 17.77 (Strauss et al. 2002, AJ 124, 1810).
MGS_R_LIMIT = 17.77

#: Redshift range fetched. Below z = 0.002 (600 km/s) a galaxy's own motion
#: swamps the cosmic expansion, so redshift says little about distance. Above
#: z = 0.3 the Main Galaxy Sample has essentially no galaxies left.
FETCH_Z_MIN = 0.002
FETCH_Z_MAX = 0.30


# ---------------------------------------------------------------------------
# The slice: a wedge through the universe
# ---------------------------------------------------------------------------

#: Half-width of the equatorial slice in declination. A 4-degree slice matches
#: the one in which Gott et al. (2005) measured the Sloan Great Wall. Thinner
#: and the structures break up into scattered dots; thicker and walls seen at a
#: slant smear across each other in projection.
SLICE_DEC_HALF_WIDTH = 2.0

#: The slice's northern (spring-sky) arc, in right ascension, degrees: the
#: survey covers the equator from RA 130 to 250 here, with a small margin
#: either side. Beyond it the equator runs toward the Milky Way's disc, whose
#: dust hides the universe behind it and which SDSS did not survey.
SLICE_RA_RANGE = (128.0, 252.0)

#: Outer redshift of the wedge plot.
WEDGE_Z_MAX = 0.20


# ---------------------------------------------------------------------------
# The volume-limited sample: a fair census for measuring
# ---------------------------------------------------------------------------
#
# A flux-limited survey sees faint galaxies only nearby, so its density falls
# with distance for reasons that have nothing to do with the universe. For any
# measurement we switch to a *volume-limited* sample: every galaxy brighter than
# a luminosity threshold, out to the distance at which the faintest of them
# would still be bright enough to make the flux limit. Inside that volume the
# census is complete, so density means density.
#
# The sample below mirrors Zehavi et al. (2011, ApJ 736, 59) threshold sample
# "Mr < -20" so our clustering can be graded against theirs.

#: Faint flux limit for the volume-limited sample. Stricter than 17.77 because
#: early SDSS regions were targeted with an r < 17.6 limit; 17.6 is complete
#: everywhere (Zehavi et al. 2011, section 2.1).
VL_R_LIMIT = 17.6

#: Luminosity threshold, in the literature's convention: r-band absolute
#: magnitude K-corrected to the z = 0.1 bandpass, evolved to z = 0.1, and
#: quoted for h = 1 (written ``M_r - 5 log h``).
VL_MAG_LIMIT_H = -20.0

#: Redshift range of Zehavi et al.'s Mr < -20 sample: cz from 6,000 to
#: 31,900 km/s. We adopt it, then verify with our own K-corrections that the
#: sample really is complete to z_max.
VL_Z_MIN = 0.02
VL_Z_MAX = 31_900 / 299_792.458

#: Luminosity evolution per unit redshift in the r band, Q = 1.62: galaxies
#: were brighter in the past. Blanton et al. (2003, ApJ 592, 819); used by
#: Zehavi et al. to evolve every galaxy to its z = 0.1 equivalent.
EVOLUTION_Q = 1.62


# ---------------------------------------------------------------------------
# Clustering measurement
# ---------------------------------------------------------------------------

#: HEALPix resolution of the angular mask: 0.84 square degrees per pixel,
#: ~65 main-sample galaxies each. Coarse enough that no pixel inside the survey
#: is empty by chance, fine enough to trace the survey's ragged edge.
MASK_NSIDE = 64

#: A pixel counts as surveyed if it holds at least this fraction of the median
#: pixel's galaxy count; the clustering sample then keeps only pixels whose
#: eight neighbours are surveyed too, which trims the partially covered edge.
MASK_MIN_FRACTION = 0.3

#: Random points per real galaxy in the clustering measurement.
RANDOM_FACTOR = 3

#: Jackknife regions for error bars: the footprint is cut into this many
#: compact, equal-area patches and the measurement is repeated leaving each out.
N_JACKKNIFE = 50

#: Projected-separation bins: 0.2 dex wide from 0.126 to 31.6 h^-1 Mpc, the
#: binning of Zehavi et al. (2011) so the two measurements line up point for
#: point.
RP_LOG_EDGES = (-0.9, 1.5, 0.2)

#: Line-of-sight integration limit for w_p, h^-1 Mpc, as used by Zehavi et al.
#: for this sample.
PI_MAX = 60.0

#: Power-law fit range for w_p, h^-1 Mpc. Zehavi et al. fit r_p < 20. The lower
#: end is raised to 0.5 because we do not correct for fibre collisions, which
#: remove close pairs below ~0.2 h^-1 Mpc (see README).
FIT_RP_RANGE = (0.5, 20.0)

#: Linear bin width for the 2D xi(r_p, pi) map, h^-1 Mpc, and its extent.
XI2D_BIN = 1.0
XI2D_MAX = 30.0


# ---------------------------------------------------------------------------
# Galaxy environment and groups
# ---------------------------------------------------------------------------

#: Local density: projected distance to the 5th-nearest neighbour among
#: galaxies within +/-1000 km/s along the line of sight, the standard
#: Sigma_5 estimator (e.g. Baldry et al. 2006, MNRAS 373, 469). The velocity
#: window is wide enough to keep cluster members whose own motions smear them
#: out along the line of sight -- the "fingers of God".
DENSITY_NTH = 5
DENSITY_DV_KMS = 1000.0

#: Friends-of-friends group finder linking lengths, as fractions of the mean
#: inter-galaxy separation: perpendicular and parallel to the line of sight.
#: The parallel length is ~5x longer because cluster members are smeared out
#: along the line of sight by their orbital velocities. Values from Berlind
#: et al. (2006, ApJS 167, 1).
FOF_B_PERP = 0.14
FOF_B_PAR = 0.75

#: Looking up an Abell cluster name for a group we found: the group's centre
#: must lie within one Abell radius, 1.5 h^-1 Mpc projected -- the aperture
#: Abell counted galaxies in -- and, where the Abell cluster has a measured
#: redshift, agree with it to within 0.01. A fixed angle would not do: nearby
#: Coma spans over a degree, a cluster at z = 0.1 a quarter of that.
ABELL_MATCH_RADIUS_H = 1.5
ABELL_MATCH_DZ = 0.01

#: Seed for every random number drawn in the project, so reruns reproduce.
SEED = 20260929
