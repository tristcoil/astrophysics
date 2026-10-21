"""Project configuration: paths, data sources, the EPIC bands, physical constants.

Every number used anywhere in the project is defined here, once, so each value
quoted in the README or printed by ``run.py`` traces back to a single line.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
CACHE_DIR = DATA_DIR / "cache"  # derived tables: metadata table, photometry
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"  # provenance of every download
RESULTS_PATH = DATA_DIR / "results.json"

for _d in (DATA_DIR, RAW_DIR, CACHE_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data sources
# ---------------------------------------------------------------------------

#: The EPIC team's public site at NASA Goddard: a JSON API with the metadata of
#: every natural-colour image (time, Earth centroid, and the J2000 positions of
#: DSCOVR, the Moon and the Sun), and the images themselves. No key needed.
EPIC_SITE = "https://epic.gsfc.nasa.gov"
EPIC_API = f"{EPIC_SITE}/api"
EPIC_ARCHIVE = f"{EPIC_SITE}/archive"
EPIC_CREDIT = "NASA EPIC Team (DSCOVR/EPIC, epic.gsfc.nasa.gov)"

API_DIR = RAW_DIR / "api"
#: Per-day metadata from ``/api/natural/date/YYYY-MM-DD``.
API_NATURAL_DIR = API_DIR / "natural"
#: Requests in flight at once against the EPIC site, and the pause between
#: requests of one worker. The site is a small public service.
EPIC_WORKERS = 4
EPIC_POLITE_DELAY_S = 0.1

#: The calibrated ten-band images (level 1B, HDF5, ~385 MB each) live only at
#: NASA's Atmospheric Science Data Center, behind a free Earthdata Login. The
#: CMR catalogue that lists them is public.
CMR_GRANULES = "https://cmr.earthdata.nasa.gov/search/granules.json"
L1B_COLLECTION_V3 = "C3860110022-LARC_CLOUD"  # DSCOVR_EPIC_L1B version 3
L1B_COLLECTION_V4 = "C4025357058-LARC_CLOUD"  # DSCOVR_EPIC_L1B version 4
ASDC_DOCS = "https://asdc.larc.nasa.gov/documents/dscovr/"
DOCS = {
    "calibration_factors_v03": "DSCOVR_EPIC_Calibration_Factors_V03.pdf",
    "dfcb_v03": "DSCOVR_EPIC_Data_Format_Control_Book_V03.pdf",
    "dfcb_2025": "EPIC_Data_Format_Control_Book_2025-12-19.pdf",
}
DOCS_DIR = RAW_DIR / "docs"

#: Natural-colour images, three sizes of the same picture: a 2048-pixel PNG
#: (~3 MB), a 1080-pixel JPEG (~0.2 MB) and a 120-pixel JPEG thumbnail
#: (~6 kB). Linearised disc-averaged brightness agrees between the PNG and the
#: 1080-pixel JPEG to 0.2%, but the thumbnails average the non-linear display
#: values over ~17 x 17 pixels and come out 2-4% low once linearised (tested on
#: 2024-03-20), so thumbnails are used only where display units suffice.
NATURAL_DIR = RAW_DIR / "natural"
IMAGE_KINDS = {"png": "png", "jpg": "jpg", "thumbs": "jpg"}

#: Years downloaded at 1080 pixels. The first is the year Yang et al. (2018)
#: analysed with calibrated L1B data, so their published numbers can be
#: compared with ours; the second is the best-covered recent year (29-31 days
#: in every month), used for the movie and as an independent repeat.
JPG_WINDOWS = {
    "yang2018_year": ("2015-08-06", "2016-08-03"),
    "year2024": ("2024-01-01", "2024-12-31"),
}

#: The EPIC team's galleries of special events, with full-size PNGs.
GALLERY_DIR = RAW_DIR / "gallery"
GALLERIES = {
    "lunar_transit_2015": "galleries/2015/lunar_transit",
    "lunar_transit_2016": "galleries/2016/lunar_transit",
    "lunar_transit_2021": "galleries/2021/lunar_transit",
    "lunar_occultation_2020": "galleries/2020/lunar_occultation",
}

#: NOAA's 0.25-degree land-sea mask for the OISST v2 sea-surface temperature
#: analysis (1 = ocean, 0 = land, incl. Antarctica), from NOAA PSL.
LANDMASK_URL = "https://downloads.psl.noaa.gov/Datasets/noaa.oisst.v2.highres/lsmask.oisst.nc"
LANDMASK_PATH = RAW_DIR / "lsmask.oisst.nc"

for _d in (API_DIR, API_NATURAL_DIR, DOCS_DIR, NATURAL_DIR, GALLERY_DIR):
    _d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# The ten EPIC bands
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Band:
    name: str  # as in the L1B group names, e.g. "Band443nm"
    wavelength_nm: float  # centre wavelength (EPIC instrument page, Table 1)
    fwhm_nm: float
    cal_factor: float  # reflectance per (count/s), DSCOVR_EPIC_Calibration_Factors_V03
    role: str


#: Centre wavelengths and widths from epic.gsfc.nasa.gov/about/epic (Table 1);
#: calibration factors from the ASDC sheet (accuracy "about 1-3%").
BANDS: tuple[Band, ...] = (
    Band("Band317nm", 317.5, 1.0, 1.216e-4, "ozone (absorbed)"),
    Band("Band325nm", 325.0, 2.0, 1.111e-4, "ozone (partly absorbed)"),
    Band("Band340nm", 340.0, 3.0, 1.975e-5, "UV continuum: aerosols"),
    Band("Band388nm", 388.0, 3.0, 2.685e-5, "UV continuum: aerosols, clouds"),
    Band("Band443nm", 443.0, 3.0, 8.340e-6, "blue"),
    Band("Band551nm", 551.0, 3.0, 6.660e-6, "green"),
    Band("Band680nm", 680.0, 3.0, 9.300e-6, "red"),
    Band("Band688nm", 687.75, 0.8, 2.020e-5, "oxygen B band (absorbed)"),
    Band("Band764nm", 764.0, 1.0, 2.360e-5, "oxygen A band (absorbed)"),
    Band("Band780nm", 779.5, 2.0, 1.435e-5, "near infrared: land, vegetation"),
)
BAND_BY_NM = {int(round(b.wavelength_nm)): b for b in BANDS}


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

AU_KM = 149_597_870.7
#: Earth's equatorial and polar radii (WGS-84) and the mean radius.
R_EARTH_EQ_KM = 6378.137
R_EARTH_POL_KM = 6356.752
R_EARTH_MEAN_KM = 6371.0
R_MOON_KM = 1737.4
#: EPIC's field of view (instrument page): 0.62 degrees across 2048 pixels.
EPIC_FOV_DEG = 0.62
EPIC_NPIX = 2048
#: Plate scale of the natural-colour images, pixels per radian at 1080 pixels
#: (scale by size/1080 for the PNG and thumbnails). Measured from the disc's
#: east-west extent at phase angles below 6 deg (99,750) and above (99,580):
#: the unlit crescent trims the limb by up to ~0.5% at 11 deg. The images are
#: never rescaled: the disc grows and shrinks with DSCOVR's distance.
PLATE_SCALE_1080 = 99_700.0


# ---------------------------------------------------------------------------
# Selections
# ---------------------------------------------------------------------------

#: Days whose images are not ordinary views of Earth: the Moon crosses the
#: disc (2021-02-11) or its shadow does (the 2024-04-08 total eclipse).
#: Each is analysed on its own and kept out of the light curves.
SPECIAL_DAYS = {
    "lunar_transit_2021": "2021-02-11",
    "high_cadence_2022": "2022-03-21",
    "eclipse_2024": "2024-04-08",
}
EXCLUDE_FROM_LIGHTCURVES = ("2021-02-11", "2024-04-08")
#: An image is used only if its disc is where and as large as geometry says:
#: centre within this many pixels (at 1080) of the frame centre and lit area
#: within this fraction of the predicted disc area.
MAX_CENTRE_OFFSET_PX = 6.0
DISC_AREA_TOLERANCE = 0.04
#: Frame-to-frame tone adjustments: the natural-colour pipeline brightens or
#: darkens whole frames by ~1.5% rms (display units) relative to their
#: neighbours, with tails to +-10%. Frames whose grey display brightness departs
#: from the mean of the frames before and after (each within 3 h) by more than
#: this are left out of the light curves (~4 sigma).
MAX_DISPLAY_JUMP = 0.06
JUMP_NEIGHBOUR_H = 3.0

RANDOM_SEED = 20150716  # the first EPIC lunar transit, 16 July 2015
