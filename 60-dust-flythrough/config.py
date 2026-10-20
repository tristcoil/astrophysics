"""Project configuration: paths, data sources, unit conversions and analysis choices.

Every number the project *assumes* lives here, each with its source, so that
any figure or printed result traces back to one line. Published values we
compare against live in ``literature.py``.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
EDENHOFER_DIR = RAW_DIR / "edenhofer2024"    # the 3D dust map (Zenodo 10658339)
LIT_DIR = RAW_DIR / "papers"                 # arXiv sources of the papers we check against
CATALOG_DIR = RAW_DIR / "catalogs"           # cloud distances, clusters, bubble models
VERGELY_DIR = RAW_DIR / "vergely2022"        # the second 3D map (VizieR J/A+A/664/A174)
DUSTMAPS_DIR = RAW_DIR / "dustmaps"          # dustmaps' data_dir (SFD lives here)
CACHE_DIR = DATA_DIR / "cache"               # derived arrays (grids, sightlines, fits)
FIGURE_DIR = PROJECT_DIR / "figures"
RESULTS_PATH = DATA_DIR / "results.json"
MANIFEST_PATH = DATA_DIR / "manifest.json"

for _d in (DATA_DIR, RAW_DIR, EDENHOFER_DIR, LIT_DIR, CATALOG_DIR, VERGELY_DIR,
           DUSTMAPS_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# The Edenhofer et al. (2024, A&A 685, A82) 3D dust map
# ---------------------------------------------------------------------------

#: Zenodo record v1.0.2 (concept DOI 10.5281/zenodo.8187942).
ZENODO_RECORD = "10658339"
ZENODO_FILE_URL = "https://zenodo.org/api/records/{record}/files/{key}/content"

#: Posterior mean and standard deviation on 516 HEALPix spheres (Nside 256,
#: 69-1250 pc). The 2 kpc file is the paper's appendix-D "validation"
#: reconstruction (fewer, higher-quality stars, out to 2 kpc); the authors
#: advise using the main map wherever both exist. Sizes and MD5s are Zenodo's.
EDENHOFER_FILES = {
    "main": ("mean_and_std_healpix.fits", 3_252_715_200, "10c823a5fcf81b47b6e15530bcdf54dc"),
    "2kpc": ("validation_with_less_data_but_2kpc_mean_and_std_healpix.fits",
             4_120_937_280, "0826b2f4cdfccad69bdc5e2d85bb4c54"),
}

#: The 12 posterior samples (19.5 GB, one float32 array of shape 12 x 516 x
#: 786432 starting at byte 5760). We fetch only a few adjacent shell pairs by
#: byte range (3.1 MB per sample and shell), enough to measure how the
#: log-space interpolation of the posterior *mean* differs from the mean of the
#: interpolated *samples* (see validation.py). Shell indices: pairs near 150,
#: 400, 800 and 1150 pc.
SAMPLES_FILE = "samples_healpix.fits"
SAMPLES_DATA_OFFSET = 5760
SAMPLE_SHELL_PAIRS = (138, 312, 435, 503)

# ---------------------------------------------------------------------------
# Units. The map is in the unitless "ZGR23" extinction of Zhang, Green & Rix
# (2023) per parsec. Conversions:
# ---------------------------------------------------------------------------

#: A_V at 540 nm per unit ZGR23 extinction (Edenhofer et al. 2024, sect. 5;
#: from ZGR23's published extinction curve).
ZGR23_TO_AV = 2.8
#: Gaia G band (673 nm) per unit ZGR23 (O'Neill et al. 2024, eq. 3).
ZGR23_TO_AG = 2.0407
#: Hydrogen nuclei per cm^3 for 1 ZGR23 mag/pc, assuming A_G/N_H = 4.0e-22
#: cm^2 mag (Draine 2003, R_V = 3.1), as in O'Neill et al. (2024, eq. 4):
#: n = 2.0407 / (4.0e-22 * 3.0857e18) = 1653 cm^-3.
AG_PER_NH = 4.0e-22          # mag cm^2
PC_CM = 3.0857e18            # cm per parsec
ZGR23_TO_NH = ZGR23_TO_AG / (AG_PER_NH * PC_CM)
#: Mass per hydrogen nucleus including helium and metals (1.37 m_p, as O'Neill).
MU_H = 1.37
M_PROTON_G = 1.67262192e-24
M_SUN_G = 1.98892e33

# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

#: Heliocentric Galactic Cartesian frame: x to the Galactic Centre (l = 0),
#: y to l = 90 deg (the direction of rotation), z to the North Galactic Pole.
#: The Sun sits at the origin; the map is heliocentric, so no solar offset.
#: Vergely, Lallement & Cox (2022): extinction density at 550 nm on a
#: 601 x 601 x 161 grid (5 pc sampling, 10 pc correlation length) covering
#: 3 x 3 x 0.8 kpc around the Sun.
VERGELY_FILE = "explore_cube_density_values_010pc_v2.fits"
#: ZGR23's curve at 552 nm (the table entry nearest 550 nm): A_550 = 2.658 E.
ZGR23_TO_A550 = 2.658

R_INNER = 69.0               # pc, the map's first sphere
R_OUTER_MAIN = 1250.0        # pc, the main map's outer edge
R_OUTER_2KPC = 2000.0        # pc, the validation map's outer edge
