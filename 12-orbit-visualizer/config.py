"""Project configuration: paths, the bodies we follow, the Horizons queries, and
the published numbers we grade ourselves against.

Everything comes from one source, JPL Horizons -- the ephemeris service JPL's
own mission navigators use. The bodies are chosen to span every kind of orbit
the Sun has:

* **bound and nearly circular** -- the planets, Pluto and Ceres;
* **bound and extremely elongated** -- Halley's comet and Hale-Bopp;
* **unbound, arriving from outside** -- the three interstellar objects known;
* **unbound, leaving** -- the five spacecraft on escape trajectories.

"Today" is a fixed date (``SNAPSHOT_DATE``) rather than the clock, so the
project reproduces exactly. To move it, edit the date and re-run
``fetch_data.py``: only the queries that depend on it are re-fetched.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
DATA_DIR = PROJECT_DIR / "data"
HORIZONS_DIR = DATA_DIR / "horizons"
FIGURE_DIR = PROJECT_DIR / "figures"
MANIFEST_PATH = DATA_DIR / "manifest.json"
RESULTS_PATH = DATA_DIR / "results.json"
VIEWER_PATH = FIGURE_DIR / "orbit_viewer.html"

for _d in (DATA_DIR, HORIZONS_DIR, FIGURE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# "Now", and where we stand to look at the sky
# --------------------------------------------------------------------------

#: The date of the snapshot, 00:00 TDB. TDB (Barycentric Dynamical Time) is the
#: time scale Horizons integrates in; it runs 69.2 s ahead of UTC, which matters
#: for nothing here except the one prediction quoted to the second (Voyager 1's
#: light-day crossing), where the conversion is done explicitly.
SNAPSHOT_DATE = "2026-09-29"

#: Observing site for the "sky tonight" figure: the Airy transit circle at the
#: Royal Observatory Greenwich, the prime meridian. It is a default, not a
#: requirement: any latitude/longitude works (east positive). Elevation in km.
OBSERVER = {"name": "Greenwich", "lat": 51.4778, "lon": -0.0015, "elevation_km": 0.046,
            "utc_offset_hours": 1, "tz_label": "BST"}

# --------------------------------------------------------------------------
# The bodies
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Body:
    """One object Horizons can give us a trajectory for."""

    #: Short slug, used in file names and as the dictionary key everywhere.
    key: str
    #: Display name.
    name: str
    #: Horizons COMMAND identifier.
    horizons_id: str
    #: planet | dwarf | comet | interstellar | spacecraft
    group: str
    #: astroquery ``id_type``; None lets Horizons treat integers as major bodies.
    id_type: str | None = None
    #: For comets: pick the orbit solution for the apparition nearest the epoch.
    closest_apparition: bool = False
    #: Approximate orbital period in days, used only to choose a sampling step.
    period_days: float | None = None
    #: First and last dates Horizons has a trajectory for (spacecraft only).
    start: str | None = None
    stop: str | None = None
    note: str = ""


#: Planets are fetched as *system barycentres* (Horizons IDs 1-9), not planet
#: centres. For Earth this matters: Earth's centre loops around the Earth-Moon
#: barycentre every month by 4,670 km, a wobble that has nothing to do with the
#: orbit around the Sun. For Mercury and Venus, which have no moons, the
#: barycentre is the planet.
PLANETS: list[Body] = [
    Body("mercury", "Mercury", "1", "planet", period_days=87.969),
    Body("venus", "Venus", "2", "planet", period_days=224.70),
    Body("earth", "Earth", "3", "planet", period_days=365.256,
         note="Earth-Moon barycentre"),
    Body("mars", "Mars", "4", "planet", period_days=686.98),
    Body("jupiter", "Jupiter", "5", "planet", period_days=4332.6),
    Body("saturn", "Saturn", "6", "planet", period_days=10759.2),
    Body("uranus", "Uranus", "7", "planet", period_days=30685.0),
    Body("neptune", "Neptune", "8", "planet", period_days=60190.0),
]

DWARFS: list[Body] = [
    Body("pluto", "Pluto", "9", "dwarf", period_days=90560.0, note="Pluto-Charon barycentre"),
    Body("ceres", "Ceres", "1", "dwarf", id_type="smallbody", period_days=1681.6),
]

COMETS: list[Body] = [
    Body("halley", "1P/Halley", "1P", "comet", id_type="designation",
         closest_apparition=True, period_days=27_600.0),
    Body("hale_bopp", "Hale-Bopp", "C/1995 O1", "comet", id_type="designation",
         closest_apparition=True),
]

#: The only three objects ever seen passing through from interstellar space.
#: ``perihelion`` is approximate -- it only centres the high-resolution window;
#: the analysis finds the real moment of closest approach from the data.
INTERSTELLAR: list[Body] = [
    Body("oumuamua", "1I/ʻOumuamua", "1I", "interstellar", id_type="designation",
         note="perihelion 2017-09-09"),
    Body("borisov", "2I/Borisov", "2I", "interstellar", id_type="designation",
         note="perihelion 2019-12-08"),
    Body("atlas", "3I/ATLAS", "3I", "interstellar", id_type="designation",
         note="perihelion 2025-10-29"),
]
INTERSTELLAR_PERIHELION = {"oumuamua": "2017-09-09", "borisov": "2019-12-08",
                           "atlas": "2025-10-29"}

#: The five spacecraft on escape trajectories. Start dates are the day after
#: the first epoch Horizons covers; stop dates are the end of the published
#: predictions (2100 for the Voyagers, 2050 for the others).
SPACECRAFT: list[Body] = [
    Body("voyager1", "Voyager 1", "-31", "spacecraft", start="1977-09-06", stop="2100-01-01"),
    Body("voyager2", "Voyager 2", "-32", "spacecraft", start="1977-08-21", stop="2100-01-01"),
    Body("pioneer10", "Pioneer 10", "-23", "spacecraft", start="1972-03-04", stop="2050-01-01"),
    Body("pioneer11", "Pioneer 11", "-24", "spacecraft", start="1973-04-07", stop="2050-01-01"),
    Body("new_horizons", "New Horizons", "-98", "spacecraft", start="2006-01-20",
         stop="2050-01-01"),
]

ALL_BODIES: list[Body] = PLANETS + DWARFS + COMETS + INTERSTELLAR + SPACECRAFT
BODIES = {b.key: b for b in ALL_BODIES}

#: Planet centres (not barycentres) for flyby geometry, and the Earth's centre
#: for anything "as seen from Earth". Horizons CENTER codes: 500@<id> means the
#: body's centre of mass.
CENTER_SUN = "500@10"
CENTER_SSB = "500@0"  # solar-system barycentre
CENTER_GEO = "500@399"
PLANET_CENTER = {"jupiter": "500@599", "saturn": "500@699", "uranus": "500@799",
                 "neptune": "500@899"}

#: Gravity-assist flybys, with the approximate UTC date of closest approach.
#: The window is ``FLYBY_HALF_WINDOW_DAYS`` either side, sampled every
#: ``FLYBY_STEP``; the analysis measures the actual closest approach.
FLYBYS: list[tuple[str, str, str]] = [
    ("pioneer10", "jupiter", "1973-12-04"),
    ("pioneer11", "jupiter", "1974-12-03"),
    ("pioneer11", "saturn", "1979-09-01"),
    ("voyager1", "jupiter", "1979-03-05"),
    ("voyager1", "saturn", "1980-11-12"),
    ("voyager2", "jupiter", "1979-07-09"),
    ("voyager2", "saturn", "1981-08-26"),
    ("voyager2", "uranus", "1986-01-24"),
    ("voyager2", "neptune", "1989-08-25"),
    ("new_horizons", "jupiter", "2007-02-28"),
]
FLYBY_HALF_WINDOW_DAYS = 6
FLYBY_STEP = "10m"

# --------------------------------------------------------------------------
# Time spans and sampling
# --------------------------------------------------------------------------

#: Kepler's-third-law spans. The inner planets complete tens to hundreds of
#: orbits in a century; the outer ones need several centuries for even one or
#: two. ``SAMPLES_PER_ORBIT`` sets the step: with velocities in hand, cubic
#: Hermite interpolation between samples is accurate to far better than any
#: effect we measure.
ORBIT_SPAN_INNER = ("1950-01-01", "2050-01-01")
ORBIT_SPAN_OUTER = ("1750-01-01", "2250-01-01")
ORBIT_SPAN_CERES = ("1900-01-01", "2100-01-01")
SAMPLES_PER_ORBIT = 120

#: Halley: one full orbit, perihelion 1986 to perihelion 2061, daily -- starting
#: early enough to cover the interactive viewer's clock, which opens in 1972.
HALLEY_SPAN = ("1970-01-01", "2066-01-01")
HALE_BOPP_SPAN = ("1970-01-01", "2050-01-01")

#: Interstellar objects: a long, coarse heliocentric track; a fine track through
#: perihelion; and a barycentric track from which the arrival speed is measured
#: far outside the planetary system.
ISO_LONG_SPAN = ("1900-01-01", "2100-01-01")
ISO_LONG_STEP = "5d"
ISO_PASSAGE_HALF_WINDOW_DAYS = 150
ISO_PASSAGE_STEP = "2h"
ISO_BARY_STEP = "30d"

SPACECRAFT_STEP = "5d"

#: Voyager 1's distance from Earth, for the light-day crossing.
LIGHT_DAY_SPAN = ("2026-01-01", "2028-12-31")
LIGHT_DAY_STEP = "6h"

#: The window the interactive viewer's clock runs over. It opens a day after
#: Pioneer 10's launch, the first of the five escaping spacecraft.
VIEWER_SPAN = ("1972-03-04", "2050-01-01")

# --------------------------------------------------------------------------
# Constants that are conventions rather than measurements
# --------------------------------------------------------------------------

#: The Sun's velocity relative to the Local Standard of Rest -- the average
#: motion of stars in our neighbourhood of the Galaxy. (U, V, W) in km/s:
#: towards the Galactic centre, along Galactic rotation, towards the north
#: Galactic pole. Schoenrich, Binney & Dehnen (2010), MNRAS 403, 1829; the same
#: value astropy's LSR frame uses.
SUN_PECULIAR_UVW_KMS = (11.1, 12.24, 7.25)

#: G, CODATA 2018. Only needed to turn the Sun's measured GM into a mass in kg,
#: and it is the weakest link in that conversion: GM is known to 11 digits, G
#: to 5.
G_SI = 6.67430e-11
G_SI_REL_ERR = 2.2e-5

# --------------------------------------------------------------------------
# Published values we compare against, consulted only at the end
# --------------------------------------------------------------------------

PUBLISHED = {
    "gm_sun": {
        "value_km3_s2": 132_712_440_041.279419,
        "source": "JPL planetary ephemeris DE440 (Park et al. 2021, AJ 161, 105)",
    },
    #: JPL's mean semi-major axes, fitted to its ephemeris over 1800-2050.
    "planet_a_au": {
        "values": {"mercury": 0.38709927, "venus": 0.72333566, "earth": 1.00000261,
                   "mars": 1.52371034, "jupiter": 5.20288700, "saturn": 9.53667594,
                   "uranus": 19.18916464, "neptune": 30.06992276},
        "source": ("JPL Solar System Dynamics, 'Approximate Positions of the Planets', "
                   "Table 1 (Standish), valid 1800-2050"),
    },
    "oumuamua": {
        "v_inf_kms": 26.3209, "v_inf_err_kms": 0.0099,
        "radiant_ra_deg": 279.804, "radiant_dec_deg": 33.997,
        "uvw_kms": (-11.457, -22.395, -7.746),
        "e": 1.19936, "q_au": 0.255287,
        "source": "Mamajek (2017), RNAAS 1, 21 -- from JPL orbit solution 13",
        "note": "solution 13 has no non-gravitational acceleration; the current one does",
    },
    "borisov": {
        "v_inf_kms": 32.286894,
        "radiant_ra_deg": 32.79720, "radiant_dec_deg": 59.44014,
        "source": "Bailer-Jones et al. (2020), A&A 634, A14 -- JPL orbit solution 47",
    },
    "atlas": {
        "v_inf_kms": 58.00,
        "radiant_ra_deg": 295.050, "radiant_dec_deg": -19.0693,
        "uvw_kms": (-51.25, -19.466, 18.94),
        "e": 6.141351, "q_au": 1.3564811,
        "source": ("radiant and (U,V,W): Hopkins et al. (2025), ApJL, from an early orbit; "
                   "v_inf is |(U,V,W)|; e and q: JPL orbit solution 54"),
    },
    #: Where the Voyagers left the heliosphere, detected by their own plasma
    #: and cosmic-ray instruments. Used only to mark the spot on the maps.
    "heliopause": {
        "voyager1": {"date": "2012-08-25", "au": 121.6,
                     "source": "Stone et al. (2013), Science 341, 150"},
        "voyager2": {"date": "2018-11-05", "au": 119.0,
                     "source": "Stone et al. (2019), Nature Astronomy 3, 1013"},
    },
    "voyager1_light_day": {
        "utc": "2026-11-18T10:16:07",
        "distance_km": 25_902_068_356,
        "source": ("NASA Science, 'Where are Voyager 1 and Voyager 2 now?': "
                   "'On Wednesday, Nov. 18, 2026, at 2:16:07 a.m. PST'"),
        "note": ("some press coverage gives '12:16:07 a.m. CST', which is 4 hours off; "
                 "2:16:07 a.m. PST is 4:16:07 a.m. CST"),
    },
}


def body(key: str) -> Body:
    try:
        return BODIES[key]
    except KeyError:
        raise KeyError(f"unknown body {key!r}; known: {sorted(BODIES)}") from None
