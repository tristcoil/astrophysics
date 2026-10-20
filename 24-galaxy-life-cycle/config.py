"""Project configuration: paths, data sources, cosmology and every sample cut.

The project has two halves. The first is a census of ~650,000 nearby galaxies
from the Sloan Digital Sky Survey (SDSS): their stellar masses and
star-formation rates, which galaxies are still forming stars, which have
stopped, and what decides it. The second is the history of the whole
universe's star formation, from the published compilation of Madau &
Dickinson (2014), anchored at the near end by our own SDSS measurement.

Every cut lives here, next to the reason for it, so any number in the README
can be traced back to exactly which galaxies it came from.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
#: SDSS Main Galaxy Sample galaxies, one CSV per RA strip as served.
GALAXY_RAW_DIR = DATA_DIR / "skyserver_galaxies"
#: The four MPA-JHU tables, one CSV per range of plates, as served.
GALSPEC_RAW_DIR = DATA_DIR / "skyserver_galspec"
#: Every Main Galaxy Sample *target* (observed or not), one CSV per RA strip.
TARGET_RAW_DIR = DATA_DIR / "skyserver_targets"
#: What the spectra of observed targets that are not clean galaxies turned out to be.
TARGET_SPEC_RAW_DIR = DATA_DIR / "skyserver_target_spectra"
#: The strips merged into typed catalogues.
GALAXY_CATALOG = DATA_DIR / "mgs_mpajhu_dr17.npz"
TARGET_CATALOG = DATA_DIR / "mgs_targets_dr17.npz"
#: SDSS legacy tiling geometry: mangle polygons of the spectroscopic footprint.
TILING_DIR = DATA_DIR / "tiling"
#: GALEX-SDSS-WISE Legacy Catalog, version 2 (Salim et al. 2016, 2018).
GSWLC_DIR = DATA_DIR / "gswlc"
#: LaTeX source of Madau & Dickinson (2014), whose Tables 1 and 2 we parse.
MD14_DIR = DATA_DIR / "madau_dickinson_2014"
#: SDSS colour cutouts for the gallery figure.
CUTOUT_DIR = DATA_DIR / "cutouts"
RESULTS_PATH = DATA_DIR / "results.json"
FIGURE_DIR = PROJECT_DIR / "figures"

for _d in (DATA_DIR, GALAXY_RAW_DIR, GALSPEC_RAW_DIR, TARGET_RAW_DIR, TARGET_SPEC_RAW_DIR,
           TILING_DIR, GSWLC_DIR, MD14_DIR, CUTOUT_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Where the data comes from
# ---------------------------------------------------------------------------

#: SkyServer's public SQL endpoint for SDSS Data Release 17. DR17 still carries
#: the MPA-JHU value-added tables (``galSpec*``), computed for the DR8
#: reprocessing of the original 2000-2008 galaxy survey.
SKYSERVER_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
SKYSERVER_CUTOUT_URL = "https://skyserver.sdss.org/dr17/SkyServerWS/ImgCutout/getjpeg"
DATA_RELEASE = "DR17"

#: RA strip width for the downloads: ~45,000 galaxies (or ~55,000 targets) per
#: query, inside SkyServer's limits and resumable if the connection drops.
RA_CHUNK_DEG = 15.0

#: The legacy survey's final tiling products on the SDSS Science Archive
#: Server: the footprint as mangle polygons (each tagged with its sector) and
#: the sector list (area, number of overlapping plates).
TILING_BASE_URL = "https://data.sdss.org/sas/dr17/sdss/tiling/final/"
TILING_FILES = ("sdss_tiling_polygons.fits", "sdss_sectors.fits")

#: GSWLC-2. "X" takes, for each galaxy, the deepest of the three GALEX depths;
#: "M" is the medium-deep catalogue the authors recommend for galaxies off the
#: main sequence (their passive-galaxy SFRs need deeper UV). Both list SED-fit
#: masses and SFRs from UV + optical + mid-IR photometry, Chabrier IMF.
GSWLC_BASE_URL = "https://salims.pages.iu.edu/gswlc/"
GSWLC_FILES = ("GSWLC-X2.dat.gz", "GSWLC-M2.dat.gz")

#: arXiv source of Madau & Dickinson (2014, ARA&A 52, 415): arXiv:1403.0007.
MD14_ARXIV_SRC = "https://arxiv.org/src/1403.0007"

#: Main Galaxy Sample target bits in ``legacy_target1`` / ``primTarget``:
#: GALAXY (0x40), GALAXY_BIG (0x80), GALAXY_BRIGHT_CORE (0x100), as in D2.
MGS_TARGET_BITS = 0x40 | 0x80 | 0x100

#: Flux limit of the Main Galaxy Sample: dust-corrected Petrosian r < 17.77
#: (Strauss et al. 2002).
MGS_R_LIMIT = 17.77

#: Redshift range fetched (as in D2): below z = 0.002 peculiar motions swamp
#: the Hubble flow; above z = 0.3 the sample has no galaxies left.
FETCH_Z_MIN = 0.002
FETCH_Z_MAX = 0.30


# ---------------------------------------------------------------------------
# Cosmology and the IMF: the conventions everything is quoted in
# ---------------------------------------------------------------------------

#: MPA-JHU, GSWLC and Madau & Dickinson all assume the same "concordance"
#: cosmology: H0 = 70 km/s/Mpc, Omega_m = 0.3, Omega_Lambda = 0.7. We adopt it
#: for every SDSS quantity so that masses, SFRs and volumes are mutually
#: consistent, and quote Planck 2018 (project D3) where it changes a number.
H0 = 70.0
OMEGA_M = 0.3

#: Stellar masses and SFRs depend on the assumed mass spectrum of newborn stars
#: (the initial mass function, IMF). MPA-JHU uses Kroupa (2001), GSWLC
#: Chabrier (2003), Madau & Dickinson Salpeter (1955). Their conversion
#: factors to Salpeter (Madau & Dickinson 2014, sections 3.1.1 and 3.1.2):
#: divide SFRs by 0.67 (Kroupa) or 0.63 (Chabrier), masses by 0.66 or 0.61.
IMF_TO_SALPETER_SFR = {"kroupa": 1 / 0.67, "chabrier": 1 / 0.63, "salpeter": 1.0}
IMF_TO_SALPETER_MASS = {"kroupa": 1 / 0.66, "chabrier": 1 / 0.61, "salpeter": 1.0}

#: Kroupa -> Chabrier for MPA-JHU quantities when comparing with GSWLC:
#: -0.025 dex (Salim et al. 2007, as applied in Salim et al. 2016).
KROUPA_TO_CHABRIER_DEX = -0.025


# ---------------------------------------------------------------------------
# The galaxy census
# ---------------------------------------------------------------------------

#: The SFR-mass plane and the quenching analysis use 0.02 < z < 0.085, the
#: range of Peng et al. (2010) and Renzini & Peng (2015). Below 0.02 galaxies
#: are so large on the sky that the 3" fibre samples only their nucleus and
#: peculiar velocities distort distances; beyond 0.085 the sample loses
#: low-mass galaxies quickly.
LOCAL_Z_RANGE = (0.02, 0.085)

#: For 1/V_max statistics we use a flux limit complete everywhere: early SDSS
#: regions were targeted to r < 17.6 rather than 17.77 (Zehavi et al. 2011;
#: D2), and a galaxy counted near 17.77 there would be given a volume it was
#: never searched in. The bright limit 14.5 avoids the incompleteness of very
#: bright, large galaxies, which the photometric pipeline breaks into pieces
#: (Blanton et al. 2003 use the same 14.5 < r < 17.77 window).
VMAX_R_FAINT = 17.6
VMAX_R_BRIGHT = 14.5

#: Missing-value sentinel used by MPA-JHU and SkyServer.
BAD_VALUE = -999.0

#: MPA-JHU's emission-line flux errors are underestimated; the team's DR7
#: documentation gives these scale factors from duplicate observations
#: (wwwmpa.mpa-garching.mpg.de/SDSS/DR7/raw_data.html).
LINE_ERROR_SCALE = {"h_alpha": 2.473, "h_beta": 1.882, "oiii_5007": 1.566,
                    "nii_6584": 2.039}
#: S/N required in all four BPT lines to classify a galaxy (Brinchmann et al.
#: 2004 use 3 after the scaling above).
BPT_MIN_SN = 3.0

#: Main sequence: ridge line fitted over this mass range (log M_sun), in bins of
#: 0.2 dex in mass and 0.2 dex in SFR, the binning of Renzini & Peng (2015).
MS_FIT_LOGM = (9.0, 11.0)
MS_BIN_DEX = 0.2


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------

#: Local density: projected 5th-nearest-neighbour density within +/-1000 km/s
#: (as in D2 and Peng et al. 2010), traced by a volume-limited sample complete
#: to the far edge of LOCAL_Z_RANGE. The tracer luminosity threshold is set in
#: analysis from the data so the tracers are complete at z = 0.085.
DENSITY_NTH = 5
DENSITY_DV_KMS = 1000.0
#: Galaxies whose 5th-neighbour circle reaches beyond the survey edge would get
#: a falsely low density; Peng et al. keep galaxies with >= 90% of the circle
#: inside the survey. We use the same rule on a HEALPix map of the footprint.
EDGE_MIN_FRACTION = 0.9
MASK_NSIDE = 512

#: Fibre collisions: two fibres cannot sit closer than 55" on one plate, so
#: ~7% of targets in crowded places were never observed. Each galaxy is
#: weighted by 1 / (fraction of MGS targets within 55" that got a spectrum),
#: the "target sampling rate" weighting of Peng et al. (2010).
FIBRE_COLLISION_ARCSEC = 55.0

#: Bins for the quenched fraction: 0.2 dex in mass, and in density (1 + delta).
QUENCH_LOGM_EDGES = (9.0, 11.6, 0.2)
QUENCH_DENSITY_QUANTILES = 6


# ---------------------------------------------------------------------------
# Cosmic star-formation history
# ---------------------------------------------------------------------------

#: Return fraction: the share of the mass formed into stars that goes back to
#: the gas as winds and supernova ejecta. R = 0.27 for a Salpeter IMF
#: (Madau & Dickinson 2014, section 2.3).
RETURN_FRACTION = 0.27

#: Reionization: ended by z ~ 6 (Madau & Dickinson's "epoch of reionization").
Z_REIONIZATION_END = 6.0

#: Seed for every random number drawn in the project, so reruns reproduce.
SEED = 20261003
