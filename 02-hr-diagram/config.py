"""Project configuration: paths, the Gaia queries, and every quality threshold.

All the numbers that decide which stars end up on the diagram live here, each
with the reason it has the value it does. Changing a cut means changing it in
exactly one place, and the figures and README numbers follow.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"

#: Every Gaia DR3 source with parallax > 10 mas, *before* any quality cut. We
#: download the loose sample and cut locally so the effect of each cut can be
#: shown rather than silently applied on the archive's side.
CANDIDATES_PATH = DATA_DIR / "gaia_dr3_parallax_gt_10mas.fits"
#: Star counts per sky pixel, from a random subsample of the whole catalogue.
SKY_COUNTS_PATH = DATA_DIR / "gaia_dr3_sky_counts.fits"
#: Eric Mamajek's mean dwarf sequence: spectral type, temperature, Gaia colour.
DWARF_SEQUENCE_PATH = DATA_DIR / "mamajek_dwarf_sequence.txt"
#: Gaia DR3 source_ids for a few famous stars, resolved once through SIMBAD.
NAMED_STARS_PATH = DATA_DIR / "named_stars.json"
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------
# The 100 pc sample
# --------------------------------------------------------------------------

#: parallax > 10 mas means a naive distance 1/parallax under 100 parsecs.
PARALLAX_MIN_MAS = 10.0

#: Columns fetched for every candidate. Kept to what the analysis uses, plus the
#: quality indicators needed to justify each cut.
CANDIDATE_COLUMNS = [
    "source_id",
    "ra", "dec", "l", "b",
    "parallax", "parallax_error", "parallax_over_error",
    "pmra", "pmdec",
    "ruwe",
    "astrometric_params_solved",   # 31 = 5-parameter, 95 = 6-parameter solution
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "phot_g_mean_flux_over_error",
    "phot_bp_mean_flux_over_error",
    "phot_rp_mean_flux_over_error",
    "phot_bp_rp_excess_factor",
    "radial_velocity",
]

CANDIDATE_QUERY = (
    f"SELECT {', '.join(CANDIDATE_COLUMNS)}\n"
    f"FROM gaiadr3.gaia_source\n"
    f"WHERE parallax > {PARALLAX_MIN_MAS}"
)

# Quality cuts, applied locally and in this order by ``analysis.apply_cuts``.

#: Parallax at least 10x its own uncertainty. Below this, 1/parallax is a
#: biased distance estimate and absolute magnitudes scatter by > 0.2 mag.
MIN_PARALLAX_OVER_ERROR = 10.0

#: Renormalised unit weight error: how well a single-star model fits Gaia's
#: measurements of the source's position over time. 1.0 is a perfect fit; above
#: ~1.4 the source is usually an unresolved binary or otherwise problematic
#: (Lindegren et al. 2018, Gaia technical note GAIA-C3-TN-LU-LL-124).
MAX_RUWE = 1.4

#: Consistency of the G, BP and RP photometry, via the corrected flux excess C*
#: of Riello et al. (2021, A&A 649, A3, Sect. 9.4). Sources whose BP+RP flux
#: disagrees with G by more than N standard deviations are contaminated by a
#: neighbour's light. Riello et al. offer 3 and 5 sigma; on this sample 5 sigma
#: clears 80% of the junk between the white dwarfs and the main sequence (3
#: sigma: 83%) while rejecting 6 instead of 39 of the 170 red dwarfs brighter
#: than G = 10 -- stars far too bright to be contaminated, whose C* is offset
#: by a calibration effect rather than a neighbour.
EXCESS_NSIGMA = 5.0

#: Riello et al. (2021) note the C* filter is unreliable brighter than G = 4
#: because of saturation, so it is not applied to those few stars.
EXCESS_G_MIN = 4.0

#: Faint-end limit on the BP magnitude. Riello et al. (2021, Sect. 9.2) show BP
#: fluxes become systematically too bright for G_BP between about 20.3 and 20.9,
#: which drags the colours of the faintest red dwarfs artificially blue.
MAX_BP_MAG = 20.5


# --------------------------------------------------------------------------
# The all-sky star-count map
# --------------------------------------------------------------------------

#: HEALPix level (order) for the sky map. Level 7 has 196,608 equal-area
#: pixels, each 0.21 square degrees.
SKY_HEALPIX_LEVEL = 7

#: Counting all 1.8 billion sources server-side would take hours of archive
#: time. Gaia's ``random_index`` is a random permutation of the catalogue, so
#: ``random_index < N`` is an honest random subsample of N sources, and the
#: counts scale up by (catalogue size / N).
SKY_SAMPLE_SIZE = 36_000_000

#: Gaia DR3 source_id carries the level-12 HEALPix index in its upper bits, so
#: integer division by 2**(59 - 2*level) gives the pixel at any lower level.
SKY_QUERY = (
    f"SELECT source_id / {2 ** (59 - 2 * SKY_HEALPIX_LEVEL)} AS hpx, COUNT(*) AS n\n"
    f"FROM gaiadr3.gaia_source_lite\n"
    f"WHERE random_index < {SKY_SAMPLE_SIZE}\n"
    f"GROUP BY hpx"
)


# --------------------------------------------------------------------------
# Reference values
# --------------------------------------------------------------------------

#: The Sun in Gaia's bands (Casagrande & VandenBerg 2018, MNRAS 479, L102).
#: Gaia cannot observe the Sun, so these come from synthetic photometry.
SUN_ABS_G = 4.67
SUN_BP_RP = 0.82

#: Kinematic populations by heliocentric tangential velocity, following Gaia
#: Collaboration, Babusiaux et al. (2018, A&A 616, A10, Sect. 4.3). Two velocity
#: components out of three are known, so the cuts are deliberately loose, and
#: the gaps between them keep the classes from bleeding into each other.
VTAN_THIN_DISK_MAX = 40.0      # km/s
VTAN_THICK_DISK = (60.0, 150.0)
VTAN_HALO_MIN = 200.0

#: Famous stars to find on the diagram: display name -> SIMBAD identifier. Their
#: Gaia DR3 source_ids are resolved once through SIMBAD and cached, and a star
#: is only labelled if it survives the same quality cuts as everything else.
NAMED_STARS = {
    # main sequence
    "Proxima Centauri": "Proxima Centauri",   # the nearest star to the Sun
    "Barnard's Star": "Barnard's star",
    "Wolf 359": "Wolf 359",
    "TRAPPIST-1": "TRAPPIST-1",               # seven Earth-sized planets
    "ε Eridani": "eps Eri",
    "τ Ceti": "tau Cet",
    "61 Cygni A": "61 Cyg A",                 # first stellar parallax, Bessel 1838
    "51 Pegasi": "51 Peg",                    # first planet around a Sun-like star
    "β Pictoris": "bet Pic",
    "HR 8799": "HR 8799",
    # slightly evolved and giant stars
    "δ Pavonis": "del Pav",
    "ε Tauri": "eps Tau",                     # Hyades red-clump giant
    "γ Tauri": "gam Tau",                     # Hyades red-clump giant
    # halo stars passing through the solar neighbourhood
    "Kapteyn's Star": "Kapteyn's star",
    "Groombridge 1830": "Groombridge 1830",
    # white dwarfs
    "40 Eridani B": "40 Eri B",
    "Sirius B": "Sirius B",
    "van Maanen's Star": "vMa 2",
}

#: Bright stars within 100 pc that Gaia DR3 does not contain at all: its
#: detectors saturate on them. ``fetch_data`` verifies the absence with a cone
#: search around each rather than taking it on trust.
BRIGHT_STARS_ABSENT = ["Sirius", "Vega", "Arcturus", "Procyon", "Altair",
                       "Aldebaran", "Pollux"]
