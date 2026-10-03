"""Project configuration: paths, data sources, and every tunable choice.

Everything this project measures comes from three public archives:

* **WDC-SILSO** (Royal Observatory of Belgium): the international sunspot number,
  version 2.0 -- daily since 1818, monthly since 1749, yearly since 1700 -- plus
  its hemispheric split since 1992 and four published reconstructions of the
  sunspot *group* number that reach back to 1610.
* **The RGO/USAF/NOAA sunspot-group database** as compiled and maintained by David
  Hathaway (solarcyclescience.com): the position and area of every sunspot group
  on every day since May 1874. This is what the butterfly diagram is drawn from.
* **NOAA SWPC**: the official Solar Cycle 25 forecast issued by the NOAA/NASA
  prediction panel in 2019, and SWPC's current forecast for the decline.

The literature values each result is checked against are collected at the bottom,
with their sources, so a comparison in the README can always be traced.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"
AR_DIR = DATA_DIR / "sunspot_groups"

for _d in (DATA_DIR, FIGURE_DIR, AR_DIR):
    _d.mkdir(parents=True, exist_ok=True)

MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

# --- WDC-SILSO -----------------------------------------------------------------------

SILSO_DATA_URL = "https://www.sidc.be/SILSO/DATA/"
#: Cite as "Source: WDC-SILSO, Royal Observatory of Belgium, Brussels" + this DOI.
SILSO_DOI = "https://doi.org/10.24414/qnza-ac80"

#: The sunspot-number files (semicolon-separated CSV). SILSO publishes SHA-256
#: hashes for exactly these files, and the download is checked against them.
SILSO_SERIES = {
    "daily": "SN_d_tot_V2.0.csv",
    "monthly": "SN_m_tot_V2.0.csv",
    "smoothed": "SN_ms_tot_V2.0.csv",
    "yearly": "SN_y_tot_V2.0.csv",
    "hemispheric_monthly": "SN_m_hem_V2.0.csv",
    "hemispheric_smoothed": "SN_ms_hem_V2.0.csv",
}
SILSO_HASH_MANIFEST = "manifest-paired.sha256"

#: SILSO's official table of cycle minima and maxima -- the answer key our own
#: cycle finder is tested against. (The CSV lags the web page by one entry: at the
#: time of writing it still lacks Cycle 25's maximum, which SILSO announced as
#: 160.9 in October 2024; see ``LIT_SILSO_CYCLE25_MAX``.)
SILSO_CYCLE_TABLE = "Cycles/TableCyclesMiMa.csv"

#: Reconstructions of the sunspot *group* number, the only way to see before 1700.
#: They disagree about how strong the 18th-century cycles were, and SILSO's own
#: page groups them into two "families" -- that disagreement is part of the story.
GROUP_NUMBER_SERIES = {
    "HS98": dict(
        file="GroupNumber/GNhs_y.txt",
        label="Hoyt & Schatten 1998",
        citation="Hoyt, D. V. & Schatten, K. H. 1998, Sol. Phys. 179, 189",
        family="rising trend",
    ),
    "SS16": dict(
        file="GroupNumber/GNbb2_y.txt",
        label="Svalgaard & Schatten 2016 (backbone)",
        citation="Svalgaard, L. & Schatten, K. H. 2016, Sol. Phys. 291, 2653",
        family="no trend",
    ),
    "U16": dict(
        file="GroupNumber/GNiu_y2.txt",
        label="Usoskin et al. 2016 (active days)",
        citation="Usoskin, I. G. et al. 2016, Sol. Phys. 291, 2685",
        family="rising trend",
        # Their method only covers 1749-1899; SILSO appends Hoyt & Schatten after.
        own_until=1899.5,
    ),
    "CL16": dict(
        file="GroupNumber/GNec_y.txt",
        label="Cliver & Ling 2016",
        citation="Cliver, E. W. & Ling, A. G. 2016, Sol. Phys. 291, 2763",
        family="no trend",
    ),
}

# --- The RGO / USAF / NOAA sunspot-group database -------------------------------------

AR_DATABASE_URL = "https://solarcyclescience.com/AR_Database/"
#: One fixed-width text file per year, g1874.txt onwards. The yearly files are the
#: maintained copies: the site's bundled zip predates later corrections (its g2024
#: still has a group at latitude -76 that the yearly file fixes to -6).
AR_FIRST_YEAR = 1874
AR_FILE_PATTERN = "g{year}.txt"
#: Documentation and Hathaway's own derived products. We do not use the derived
#: files as inputs -- they are the answer key our own daily areas are checked against.
AR_EXTRA_FILES = ("format.txt", "MissingDays.txt", "daily_area.txt")

#: The Royal Greenwich Observatory measured the photographic plates until the end
#: of 1976; from 1977 the data come from the US Air Force's SOON network, via NOAA.
RGO_LAST_YEAR = 1976

# --- NOAA SWPC forecasts ----------------------------------------------------------------

SWPC_URL = "https://services.swpc.noaa.gov/json/solar-cycle/"
SWPC_FILES = {
    # The NOAA/NASA Solar Cycle 25 Prediction Panel's curve (December 2019).
    "panel_2019": "solar-cycle-25-predicted.json",
    # The same curve with the panel's stated uncertainty: peak +/-10%, timing +/-6 months.
    "panel_2019_range": "solar-cycle-25-ssn-predicted-high-low.json",
    # SWPC's current forecast for the rest of the cycle.
    "current_forecast": "predicted-solar-cycle.json",
}

# --- The solar cycle, as measured from the sunspot number -------------------------------

#: Cycles are found as minima of the smoothed series at least this far apart. The
#: shortest cycle on record lasted 9.0 years; 7 leaves room without admitting the
#: brief dips inside a double-peaked maximum.
MIN_CYCLE_SEPARATION_YR = 7.0

#: Wolf's numbering: Cycle 1 is the one that began at the minimum of 1755.
FIRST_NUMBERED_MINIMUM_YEAR = 1755

#: Named episodes, for annotation. The Maunder minimum is Eddy's (1976) definition;
#: the Dalton minimum spans the weak cycles 5-7.
MAUNDER_MINIMUM = (1645, 1715)
DALTON_MINIMUM = (1790, 1830)

#: Cycle 25 "hindcast": for each lead time after minimum, predict the cycle maximum
#: from the smoothed sunspot number reached by then, calibrated on earlier cycles.
PREDICTION_LEAD_MONTHS = (12, 18, 24, 30, 36, 42)
#: Calibrate the rise-to-peak relation on these cycles. Cycles 1-9 rest on sparse
#: 18th- and early-19th-century observations; from Cycle 10 (1855) the daily record
#: is essentially complete.
PREDICTION_CALIBRATION_CYCLES = tuple(range(10, 25))

#: The old (version 1) sunspot number was scaled by 0.6 relative to version 2. Hathaway's
#: cycle-shape width relation (``solar_physics.shape_width_months``) and Hoyt &
#: Schatten's group-number scale were both calibrated on version 1.
SN_V1_OVER_V2 = 0.6

# --- The butterfly diagram and the drift of the sunspot zones ---------------------------

#: Mean synodic Carrington rotation (days): the natural time bin for sunspot maps.
CARRINGTON_SYNODIC_DAYS = 27.2753
#: JD at which Carrington rotation 1 began (1853 Nov 9), as used by sunpy/Meeus.
CARRINGTON_EPOCH_JD = 2398167.4

#: Latitude binning for the butterfly image (degrees).
BUTTERFLY_LAT_BIN_DEG = 1.0
BUTTERFLY_LAT_MAX_DEG = 50.0

#: Near each minimum the dying cycle's spots (near the equator) and the new cycle's
#: (at high latitude) coexist. Each group goes to the cycle under which it is most
#: probable: that cycle's activity at the time (its fitted shape) times a Gaussian in
#: latitude about its drift-law latitude. A new cycle's shape is zero before its
#: fitted start, so it gets this small floor (in sunspot-number units) for the
#: preceding months, to let its first high-latitude spots be claimed.
CYCLE_PRIOR_FLOOR = 1.0
CYCLE_EARLY_MONTHS = 36.0
#: Rotation bins whose total spot area (in millionths of a hemisphere, summed over
#: the days of the rotation) falls below this are too sparse for a centroid.
CENTROID_MIN_AREA = 100.0

#: After 1976 the USAF/NOAA areas are systematically smaller than RGO's. Hathaway
#: scales them up by 1.4 in his derived products; we measure the factor ourselves
#: (``analysis.measure_area_seam``) and apply this literature value only where we
#: compare with his files.
USAF_AREA_FACTOR_HATHAWAY = 1.4
#: Years either side of the seam used to measure the area-to-sunspot-number ratio.
AREA_SEAM_WINDOW_YR = 12

# --- Differential rotation from sunspot tracking ----------------------------------------

#: Only RGO data (1874-1976) are used: its positions are given to 0.1 deg and its
#: observation times to 0.001 day. USAF/NOAA positions are rounded to whole degrees.
ROTATION_YEARS = (1874, 1976)
#: Keep positions within this central-meridian distance: near the limb, positions
#: are foreshortened and less accurate.
ROTATION_MAX_CMD_DEG = 60.0
#: Consecutive sightings of a group must be 0.5-1.5 days apart (daily pairs).
ROTATION_PAIR_DT_DAYS = (0.5, 1.5)
#: Daily longitude shifts beyond this are misidentifications, not rotation (the
#: real spread is ~14 deg/day +/- 3).
ROTATION_MAX_DRIFT_DEG_PER_DAY = 5.0
ROTATION_LAT_BIN_DEG = 2.5

# --- Literature values we compare against ------------------------------------------------

#: Hathaway (2011, Sol. Phys. 273, 221): the sunspot-zone centroid in each hemisphere
#: follows lambda = 28 deg * exp(-(t - t0) / 90 months) for every cycle, with t0 the
#: start time from fitting the cycle shape. Table 1 gives t0 per cycle (version-1
#: sunspot numbers), as year-month.
LIT_HATHAWAY_DRIFT = dict(lambda0_deg=28.0, tau_months=90.0)
LIT_HATHAWAY_T0 = {
    12: "1878-05", 13: "1889-05", 14: "1901-08", 15: "1913-03", 16: "1923-02",
    17: "1933-11", 18: "1944-03", 19: "1954-04", 20: "1964-11", 21: "1976-06",
    22: "1986-03", 23: "1996-08",
}

#: Balthasar, Vazquez & Woehl (1986, A&A 155, 87): sidereal rotation of sunspot
#: groups from the complete Greenwich record 1874-1976,
#: Omega = 14.551 - 2.87 sin^2(latitude) deg/day.
LIT_BALTHASAR_1986 = dict(A=14.551, B=-2.87)

#: The NOAA/NASA Solar Cycle 25 Prediction Panel (December 2019): minimum in April
#: 2020 (+/-6 months), maximum of 115 (range 105-125) in July 2025 (+/-8 months).
LIT_PANEL_2019 = dict(max_ssn=115.0, max_range=(105.0, 125.0), max_date="2025-07",
                      max_date_pm_months=8, min_date="2020-04", min_date_pm_months=6)

#: SILSO's announced Cycle 25 maximum: smoothed sunspot number 160.9, October 2024.
LIT_SILSO_CYCLE25_MAX = dict(value=160.9, date="2024-10")

#: The USAF/NOAA-to-RGO area scale factor: Hathaway uses 1.4; Balmaceda et al.
#: (2009, JGR 114, A07104) found 1.43 for projected and 1.49 for corrected areas.
LIT_AREA_FACTOR = dict(hathaway=1.4, balmaceda_projected=1.43, balmaceda_corrected=1.49)

#: Group number -> sunspot number: Hoyt & Schatten scaled G by 12.08 to match the
#: version-1 sunspot number; on the version-2 scale that is 12.08 / 0.6 = 20.1.
LIT_HS98_SCALE_V2 = 12.08 / SN_V1_OVER_V2
