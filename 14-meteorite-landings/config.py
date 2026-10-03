"""Project configuration: paths, pinned data sources, and every tunable choice.

The project runs on three downloads, all small and all free:

* NASA's copy of the Meteoritical Society's catalogue of meteorites (the main
  data set -- a frozen 2013 snapshot of the Meteoritical Bulletin Database);
* Natural Earth country outlines (public domain), for the base map and for
  working out which country each meteorite landed in;
* World Bank population figures, for the test of whether witnessed falls
  follow people.

Every literature number the analysis compares against also lives here, with
its source, so the comparison is auditable in one place.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
FIGURE_DIR = PROJECT_DIR / "figures"

for _d in (DATA_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"

# --- The meteorite catalogue --------------------------------------------------------

#: NASA's open-data portal was rebuilt in 2025; the old Socrata links
#: (data.nasa.gov/resource/gh4g-9sfh) are dead and the file now lives here.
METEORITE_URL = (
    "https://data.nasa.gov/docs/legacy/meteorite_landings/Meteorite_Landings.csv"
)
METEORITE_LANDING_PAGE = "https://data.nasa.gov/dataset/meteorite-landings"
METEORITE_CSV = DATA_DIR / "Meteorite_Landings.csv"
#: The file is a frozen snapshot, so its checksum should never change. If it
#: does, fetch_data.py says so loudly: every number in the README came from this one.
METEORITE_SHA256 = "5b4c6a3596f177b3fe8b722e8a7e5f98f05628fe2167ec3bd007acdeb77036b5"
METEORITE_ROWS = 45_716
#: The snapshot was taken in 2013 (it holds the February 2013 Chelyabinsk fall;
#: NASA's notes date the update to 14 May 2013). Later years are errors.
SNAPSHOT_YEAR = 2013

# --- Natural Earth country outlines -------------------------------------------------

#: Natural Earth is maintained on GitHub by its editors; we pin a release tag.
NATURAL_EARTH_TAG = "v5.1.2"
_NE = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
       f"{NATURAL_EARTH_TAG}/geojson/ne_{{res}}_admin_0_countries.geojson")
#: 1:50 million -- for locating meteorites and measuring areas (coastlines good
#: to a few km), and for the printed maps.
COUNTRIES_50M_URL = _NE.format(res="50m")
COUNTRIES_50M = DATA_DIR / "ne_50m_admin_0_countries.geojson"
#: 1:110 million -- a lighter outline for the interactive map.
COUNTRIES_110M_URL = _NE.format(res="110m")
COUNTRIES_110M = DATA_DIR / "ne_110m_admin_0_countries.geojson"

# --- World Bank population -----------------------------------------------------------

#: The witness test compares falls in POPULATION_YEARS with the population of the
#: same years. World Bank series start in 1960.
POPULATION_YEARS = (1960, 2009)
POPULATION_URL = (
    "https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL"
    f"?date={POPULATION_YEARS[0]}:{POPULATION_YEARS[1]}&format=json&per_page=20000"
)
POPULATION_JSON = DATA_DIR / "worldbank_population.json"

# --- The Earth ------------------------------------------------------------------------

#: Radius of the sphere with the same surface area as the WGS84 ellipsoid (the
#: "authalic" radius). With it, sphere areas are true areas: 4 pi R^2 = 510.07
#: million km^2, the Earth's surface.
EARTH_RADIUS_KM = 6371.0072

#: South of this latitude is Antarctica (the Antarctic Treaty boundary; there is
#: no other land between 60 S and the Antarctic coast).
ANTARCTIC_LATITUDE = -60.0

# --- Cleaning the catalogue ------------------------------------------------------------

#: Every year before the Gregorian calendar reform is one too early in NASA's
#: file (Ensisheim, which fell in November 1492, is listed as 1491). The pattern
#: is what happens when a 1 January timestamp is converted between the Gregorian
#: and Julian calendars: before October 1582 the two differ by 4-10 days, which
#: pushes 1 January back into the previous year. We add the year back.
CALENDAR_REFORM_YEAR = 1583

#: Found on Mars by the Opportunity rover in 2005 ("Heat Shield Rock"). Its
#: coordinates are Martian (longitude 354.5 E), so it is left off every Earth map.
FOUND_ON_MARS = ("Meridiani Planum",)

#: Names used for meteorites bought through the trade, with the find site known
#: only to a region. All four regions are hot deserts.
DESERT_TRADE_NAMES = ("Northwest Africa", "Northeast Africa", "Sahara", "Arabian Peninsula")

#: "Nova NNN" is the Meteoritical Society's name for a meteorite whose origin is
#: unknown -- no environment can be assigned.
UNKNOWN_ORIGIN_NAMES = ("Nova",)

#: A few Antarctic search areas have no coordinates at all in this file (Larkman
#: Nunatak, Roberts Massif, ...). Antarctic meteorites are numbered season year +
#: serial -- Allan Hills 84001 was collected in 1984 -- while desert areas
#: count 001, 002, ... So a 5- or 6-digit number that starts with the find year
#: (give or take one: field seasons straddle New Year) marks an Antarctic find.
#: The rule is only used where coordinates are missing, and analysis.py measures
#: how often it is right on the meteorites whose coordinates we do know.
ANTARCTIC_NUMBER_DIGITS = (5, 6)

# --- The witness test -------------------------------------------------------------

#: Falls whose recovered mass reaches this are compared with the published flux of
#: falls that drop at least this much.
WITNESS_MIN_MASS_G = 100.0

#: Published rates of meteorite falls reaching the ground, per million km^2 per
#: year, with at least ``mass_g`` of meteorites. Camera networks count the
#: fireballs directly; desert estimates date the stones that accumulated.
FLUX_REFERENCES = (
    dict(mass_g=10, rate=83, method="camera network",
         source="Halliday, Blackwell & Griffin 1989, Meteoritics 24, 173 "
                "(as quoted by Drouard et al. 2019)"),
    dict(mass_g=100, rate=58, method="camera network",
         source="Halliday, Blackwell & Griffin 1989 (as quoted by Borovicka et al. 2022)"),
    dict(mass_g=1000, rate=9, method="camera network",
         source="Halliday, Blackwell & Griffin 1989 (as quoted by Borovicka et al. 2022)"),
    dict(mass_g=100, rate=48, method="camera network",
         source="Borovicka et al. 2022, A&A 667, A158 (European Fireball Network)"),
    dict(mass_g=10, rate=222, method="desert accumulation",
         source="Drouard et al. 2019, Geology 47, 673 (Atacama, 2-million-year average)"),
)

#: The reference flux for the witness test: Halliday et al.'s falls of >= 100 g.
WITNESS_REFERENCE_RATE = 58.0

# --- A published table to reproduce --------------------------------------------------

#: Evatt et al. 2016, Nature Communications 7, 10679, Table 1: iron and stony-iron
#: meteorites as a share of each collection, from the Meteoritical Bulletin
#: Database as of 18 December 2015. (k = iron + stony-iron, n = all.)
EVATT_2016 = {
    "LaPaz Icefield": dict(k=5, n=1670),
    "Frontier Mountain": dict(k=0, n=798),
    "All Antarctic finds": dict(k=239, n=34927),
    "Falls (rest of world)": dict(k=60, n=1097),
    "Finds (rest of world)": dict(k=1145, n=16650),
}

# --- Presentation -----------------------------------------------------------------------

#: Collection areas worth naming on the maps, and milestone years for the timeline.
MILESTONES = (
    (1803, "L'Aigle: Biot shows\nstones fall from the sky"),
    (1969, "Yamato: meteorites\nfound concentrated on ice"),
    (1976, "US Antarctic searches\n(ANSMET) begin"),
)
