"""Project configuration: where the data comes from, and every sample cut.

The project joins four public catalogues -- Galaxy Zoo 2's volunteer votes,
SDSS photometry and redshifts, MPA-JHU stellar masses, and the expert
classifications of Nair & Abraham (2010) -- and then cuts them in different
ways for different questions. Every cut lives here, next to the reason for
it, so each figure can be traced back to exactly which galaxies it shows.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
#: Galaxy Zoo tables, exactly as served.
GZ_DIR = DATA_DIR / "galaxyzoo"
#: One CSV per right-ascension strip: GZ2 galaxies joined with SDSS DR17.
GZ2_SDSS_DIR = DATA_DIR / "skyserver_gz2"
#: One CSV per strip: the Main Galaxy Sample, the density tracers of project D2.
MGS_DIR = DATA_DIR / "skyserver_mgs"
VIZIER_DIR = DATA_DIR / "vizier"
#: JPEG cutouts, one sub-folder per use (poster, explorer, cnn).
CUTOUT_DIR = DATA_DIR / "cutouts"
#: The rows of the Galaxy Zoo DESI catalogue that match our galaxies.
GZDESI_DIR = DATA_DIR / "gzdesi"
#: The merged, typed catalogues built from the raw downloads.
CATALOG_PATH = DATA_DIR / "gz2_catalog.npz"
MGS_PATH = DATA_DIR / "mgs_tracers.npz"
RESULTS_PATH = DATA_DIR / "results.json"
FIGURE_DIR = PROJECT_DIR / "figures"

for _d in (DATA_DIR, GZ_DIR, GZ2_SDSS_DIR, MGS_DIR, VIZIER_DIR, CUTOUT_DIR, GZDESI_DIR,
           FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Where the data comes from
# ---------------------------------------------------------------------------

#: Galaxy Zoo 2, re-debiased by Hart et al. (2016, MNRAS 461, 3663): 239,695
#: galaxies with spectroscopic redshifts and r <= 17.0. The Galaxy Zoo team
#: "strongly advise" this table over the original Willett et al. (2013) one,
#: because it corrects every question of the decision tree consistently.
HART16_URL = "https://gz2hart.s3.amazonaws.com/gz2_hart16.csv.gz"

#: SkyServer's SQL endpoint for SDSS Data Release 17. It holds the original
#: GZ2 table (Willett et al. 2013, "Table 5", as ``zoo2MainSpecz``) keyed by
#: the DR7, DR8 and spectrum IDs, so photometry, redshifts and stellar masses
#: can be joined to every GZ2 galaxy on the server.
SKYSERVER_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
DATA_RELEASE = "DR17"

#: SkyServer's image cutout service: the same service that drew the 424 x 424
#: pixel gri colour images the GZ2 volunteers classified (Willett et al. 2013,
#: section 2.2), now serving the DR17 reprocessing of the same exposures.
IMGCUTOUT_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/ImgCutout/getjpeg"

#: Expert classifications of 14,034 SDSS galaxies by Nair & Abraham (2010,
#: ApJS 186, 427), on VizieR: Hubble T-types and bar types, used to check the
#: volunteers against professionals.
NA10_CATALOG = "J/ApJS/186/427/table2"

#: The K-correction authors' own code, fetched only so the tests can check
#: our transcription of their coefficients against it.
KCOR_CODE_URL = "http://kcor.sai.msu.ru/calc_kcor.py"

#: Galaxy Zoo DESI (Walmsley et al. 2023, MNRAS 526, 4768): the Zoobot deep
#: learning model's predicted vote fractions for 8.7 million galaxies in the
#: deeper DESI Legacy Surveys images. Only a handful of columns are read, by
#: HTTP byte range, never the whole 659 MB file.
GZDESI_URL = ("https://zenodo.org/api/records/8360385/files/"
              "gz_desi_deep_learning_catalog_friendly.parquet/content")

#: GZ2 sample limit: Petrosian r magnitude, corrected for Milky Way dust,
#: measured in DR7 (Willett et al. 2013, section 2.1).
GZ2_R_LIMIT = 17.0

#: GZ2 image recipe: 424 x 424 pixels at 0.02 x petroR90_r arcsec per pixel,
#: so every galaxy fills the frame the same way whatever its size.
GZ2_IMAGE_PIXELS = 424
GZ2_IMAGE_SCALE_PER_R90 = 0.02

#: SkyServer is queried in right-ascension strips so no single query returns
#: more than ~25,000 rows, and an interrupted download resumes where it stopped.
GZ2_RA_CHUNK_DEG = 10.0
MGS_RA_CHUNK_DEG = 15.0


# ---------------------------------------------------------------------------
# Reading the votes
# ---------------------------------------------------------------------------
#
# Galaxy Zoo asks a tree of questions: "smooth, or features/disk?" first; only
# volunteers who answer "features" are asked whether it is edge-on; only those
# who say "not edge-on" are asked about bars and spiral arms, and so on. A
# vote fraction deeper in the tree is therefore *conditional* on the answers
# above it, and the probability that a galaxy is, say, a face-on spiral is the
# product p_features x p_not-edge-on x p_spiral.

#: The "clean" threshold of Willett et al. (2013): a debiased vote fraction
#: above 0.8 for the answer in question. Conservative -- many real spirals and
#: ellipticals fall below it -- but the samples it picks are pure.
CLEAN_P = 0.8

#: Hart et al. (2016) define spirals by p_features x p_not-edge-on x p_spiral
#: > 0.5, and require at least 5 votes that name a definite number of arms
#: ("can't tell" does not count) before assigning an arm number.
HART_SPIRAL_P = 0.5
HART_MIN_ARM_VOTES = 5

#: Whether a galaxy was shown the bar question often enough to answer it:
#: Willett et al. (2013) Table 3, for at least 20 votes on the bar question.
#: p_features > 0.430 and p_not-edge-on > 0.715.
BAR_MIN_P_FEATURES = 0.430
BAR_MIN_P_NOTEDGEON = 0.715
BAR_MIN_VOTES = 20

#: A galaxy counts as barred when at least half of the volunteers who were
#: asked said so -- the majority rule of Masters et al. (2011), which picks
#: out strong bars (Willett et al. 2013, section 5.2.3).
BAR_P = 0.5

#: For our own bar census a disk needs at least 5 votes on the bar question,
#: so that its bar fraction is not a guess from one or two volunteers.
BAR_MIN_VOTES_DEBIASED = 5

#: Bars cannot be seen in steeply inclined disks: keep disks whose exponential
#: fit has axis ratio b/a > 0.5 (inclination below ~60 degrees; Masters et al.
#: 2011 used log(a/b) < 0.3, the same cut).
BAR_MIN_AXIS_RATIO = 0.5


# ---------------------------------------------------------------------------
# A fair census: volume-limited samples
# ---------------------------------------------------------------------------
#
# GZ2 is flux-limited: far away, only luminous galaxies are bright enough to
# be in it, so any fraction measured over the whole sample mixes up "how
# common" with "how far". Every fraction quoted is therefore measured in a
# sample that is complete for the galaxies it includes.

#: The completeness limit we test against, in DR17 magnitudes. GZ2 was cut at
#: r = 17.0 in DR7's photometry; DR8 onwards re-measured the same images with
#: a different sky subtraction, which moves individual magnitudes by a few
#: hundredths. A 0.1 mag margin keeps the samples complete despite that.
COMPLETENESS_R = 16.9

#: The redshift range of the main volume-limited sample: the Galaxy Zoo
#: standard (Bamford et al. 2009; Masters et al. 2010; Hart et al. 2016). Below
#: z = 0.03 a galaxy's own motion is a sizeable part of its redshift and the
#: biggest galaxies overflow the SDSS image processing; above 0.085 the GZ2
#: luminosity limit leaves only the most luminous galaxies.
VL_Z_RANGE = (0.03, 0.085)

#: Stellar mass bins (log solar masses) for every "at fixed mass" result.
#: Each bin uses only the redshifts out to which it is complete (see
#: ``analysis.mass_completeness``).
MASS_BINS = (10.0, 10.3, 10.6, 10.9, 11.2, 11.5)

#: Fraction of galaxies a mass bin may lose to the flux limit at its far edge:
#: the bin's redshift limit is where the reddest 2% of its galaxies would
#: drop out. The lost 2% are not random -- they are the reddest -- so this
#: limit is tested in the README by halving it.
MASS_COMPLETENESS_LOSS = 0.02


# ---------------------------------------------------------------------------
# Environment: project D2's density measure
# ---------------------------------------------------------------------------
#
# Local density is measured exactly as in project D2 (07-cosmic-web): the
# surface density out to the 5th-nearest neighbour within +/-1000 km/s,
# Sigma_5 = 5 / (pi d_5^2), counting only *tracer* galaxies from a sample that
# is complete throughout the volume -- so the measure means the same thing at
# every redshift.

#: Tracers: the Main Galaxy Sample (Strauss et al. 2002) as D2 fetched it:
#: target bits GALAXY / GALAXY_BIG / GALAXY_BRIGHT_CORE, zWarning = 0,
#: r < 17.77 ...
MGS_TARGET_BITS = 0x40 | 0x80 | 0x100
MGS_R_LIMIT = 17.77
MGS_FETCH_Z = (0.002, 0.30)
#: ... cut to D2's volume-limited sample: M_r - 5 log h < -20 in Zehavi et al.
#: (2011)'s convention, r < 17.6, 0.02 < z < 0.1064.
TRACER_R_LIMIT = 17.6
TRACER_Z_RANGE = (0.02, 31_900 / 299_792.458)
TRACER_M_LIT = -20.0
EVOLUTION_Q = 1.62
DENSITY_NTH = 5
DENSITY_DV_KMS = 1000.0

#: Number of equal-count density bins for environment trends.
N_DENSITY_BINS = 5


# ---------------------------------------------------------------------------
# Pictures
# ---------------------------------------------------------------------------

#: Tuning-fork poster: galaxies are chosen by rule, never by eye. Within each
#: Hubble class (defined in ``poster.py`` from vote fractions), candidates must
#: be big enough on the sky to show their structure ...
POSTER_MIN_R90_ARCSEC = 20.0
#: ... nearby enough that the volunteers saw them well ...
POSTER_Z_RANGE = (0.005, 0.06)
#: ... and the best-voted few are shown.
POSTER_PER_CLASS = 3

#: Seed for every random number drawn in the project, so reruns reproduce.
SEED = 20261003
