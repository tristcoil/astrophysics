"""Load NASA's meteorite table, repair what can be repaired, flag what cannot, and
label every meteorite with its family and the kind of place it came from.

The raw file has one row per meteorite: name, id, nametype, recclass, mass (g),
fall, year, reclat, reclong. Four things in it need care before any map or
statistic means anything:

* **Missing locations wear disguises.** 7,315 rows have no coordinates, and
  another 6,214 sit at exactly (0, 0) -- a placeholder, not a meteorite in the
  Gulf of Guinea. A further 223 Yamato meteorites (Antarctica, 71.5 S) have had
  their latitude zeroed, which would put an Antarctic collection area in Kenya.
* **One meteorite was found on Mars** (Meridiani Planum, by the Opportunity rover).
* **Early years are one too small**: every year before the 1582 calendar reform
  is shifted (see ``config.CALENDAR_REFORM_YEAR``), and one year (2101) is a typo.
* **A name encodes how a meteorite was found.** In a *dense collection area* --
  an Antarctic ice field or a desert plain that yields many meteorites -- each
  one is named after the area plus a number (Allan Hills 84001, Dhofar 019).
  A meteorite found on its own gets a place name of its own (Willamette).
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

import geo
from config import (
    ANTARCTIC_LATITUDE,
    ANTARCTIC_NUMBER_DIGITS,
    CALENDAR_REFORM_YEAR,
    COUNTRIES_50M,
    DESERT_TRADE_NAMES,
    FOUND_ON_MARS,
    METEORITE_CSV,
    SNAPSHOT_YEAR,
    UNKNOWN_ORIGIN_NAMES,
)

# --- Families -------------------------------------------------------------------------

#: Broad families in drawing order. Chondrites are unmelted rock from asteroids
#: that never melted; achondrites, stony-irons and irons come from bodies that
#: melted and separated into crust, mantle and metal core.
FAMILIES = ("ordinary", "carbonaceous", "other_chondrite", "achondrite", "stony_iron", "iron")

FAMILY_NAMES = {
    "ordinary": "Ordinary chondrites",
    "carbonaceous": "Carbonaceous chondrites",
    "other_chondrite": "Enstatite & other chondrites",
    "achondrite": "Achondrites",
    "stony_iron": "Stony-irons",
    "iron": "Irons",
    "unclassified": "Unclassified",
}

_ACHONDRITE_GROUPS = (
    # (prefix, group) -- tested in order
    ("Martian", "Martian"),
    ("Lunar", "Lunar"),
    ("Eucrite", "HED"), ("Howardite", "HED"), ("Diogenite", "HED"),
    ("Ureilite", "Ureilite"),
    ("Aubrite", "Aubrite"), ("Enst achon", "Aubrite"),
    ("Angrite", "Angrite"),
    ("Acapulcoite", "Primitive achondrite"), ("Lodranite", "Primitive achondrite"),
    ("Winonaite", "Primitive achondrite"), ("Brachinite", "Primitive achondrite"),
    ("Achondrite-prim", "Primitive achondrite"),
    ("Achondrite-ung", "Ungrouped achondrite"),
)

#: Chondrite group letters must be followed by a petrologic type, a qualifier or
#: nothing -- never a lower-case letter, which would make it a word such as
#: "Howardite" or "Lunar".
_OC = re.compile(r"^(OC|LL|L|H)(?![a-z])")
_CC = re.compile(r"^C(?=[IMOVKRHB0-9/-]|$)")
_EC = re.compile(r"^E(H|L)?(?![a-z])")
_RK = re.compile(r"^(R|K)(?![a-z])")
_UNCLASSIFIED = ("Stone-uncl", "Stone-ung", "Unknown", "Fusion crust",
                 "Chondrite-fusion crust", "Impact melt breccia")


def classify_recclass(recclass: str) -> tuple[str, str]:
    """Map a Meteoritical Bulletin class string to (family, group).

    >>> classify_recclass("L6")
    ('ordinary', 'L')
    >>> classify_recclass("Iron, IIIAB")
    ('iron', 'IIIAB')
    """
    s = str(recclass).strip()
    if s.startswith("Relict"):
        return "unclassified", "Relict"
    if s in _UNCLASSIFIED:
        return "unclassified", s
    if s.startswith("Iron"):
        grp = s.split(",", 1)[1].strip() if "," in s else "Iron (ungrouped or unstudied)"
        return "iron", grp.rstrip("?")
    if s.startswith("Pallasite"):
        return "stony_iron", "Pallasite"
    if s.startswith("Mesosiderite"):
        return "stony_iron", "Mesosiderite"
    for prefix, grp in _ACHONDRITE_GROUPS:
        if s.startswith(prefix):
            return "achondrite", grp
    m = _OC.match(s)
    if m:
        # Mixed H/L and L/LL (and "L(LL)") classes keep their first letter's group
        # only when unambiguous; otherwise they are "H/L" or "L/LL".
        if s.startswith(("H/L", "H(L)", "L(H)")):
            return "ordinary", "H/L"
        if s.startswith(("L/LL", "L(LL)", "LL(L)")):
            return "ordinary", "L/LL"
        return "ordinary", m.group(1)
    if _CC.match(s):
        g = re.match(r"^C([IMOVKRHB])", s)
        return "carbonaceous", ("C" + g.group(1)) if g else "C (ungrouped)"
    m = _EC.match(s)
    if m:
        return "other_chondrite", "E" + (m.group(1) or "")
    m = _RK.match(s)
    if m:
        return "other_chondrite", m.group(1)
    if s.startswith("Chondrite-ung"):
        return "other_chondrite", "Chondrite (ungrouped)"
    raise ValueError(f"unrecognised meteorite class {recclass!r}")


# --- Names ---------------------------------------------------------------------------------

#: "<Area> <number>", where the number may carry a letter prefix (Allan Hills
#: A77005) and rarely a suffix. Everything else is an individual name.
_NUMBERED = re.compile(r"^(?P<area>.+?)\s+(?P<prefix>[A-Z]{0,4})\s?(?P<digits>\d{2,6})"
                       r"[a-zA-Z]?(?:[,/-]\d+)*$")


def parse_name(name: str) -> tuple[str | None, str | None]:
    """(collection area, digits) for a numbered name; (None, None) otherwise."""
    m = _NUMBERED.match(str(name).strip())
    if not m:
        return None, None
    return m.group("area"), m.group("digits")


def antarctic_style_number(digits: str | None, year: float) -> bool:
    """Does a number follow the Antarctic scheme: season year + serial?

    Antarctic meteorites are numbered YY + serial (ANSMET: 5 digits; the Japanese
    and Chinese programmes: 6), YY being the field season -- which straddles New
    Year, so the find year may be YY or the year either side.
    """
    if digits is None or not np.isfinite(year) or len(digits) not in ANTARCTIC_NUMBER_DIGITS:
        return False
    yy = int(digits[:2])
    return any(yy == (int(year) + d) % 100 for d in (-1, 0, 1))


# --- Loading --------------------------------------------------------------------------------

def load_raw(path=METEORITE_CSV) -> pd.DataFrame:
    df = pd.read_csv(path)
    return df.rename(columns={
        "mass (g)": "mass_g", "reclat": "lat", "reclong": "lon", "recclass": "recclass",
    }).drop(columns=["GeoLocation"])


def clean(raw: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Return the cleaned catalogue and a ledger of every repair and flag."""
    df = raw.copy()
    ledger: dict = {"rows": len(df)}

    df["valid"] = df["nametype"].eq("Valid")
    df["fell"] = df["fall"].eq("Fell")
    ledger["relict"] = int((~df["valid"]).sum())

    # Years.
    df["year_raw"] = df["year"]
    early = df["year"] < CALENDAR_REFORM_YEAR
    df.loc[early, "year"] = df.loc[early, "year"] + 1
    df["year_calendar_fixed"] = early
    future = df["year"] > SNAPSHOT_YEAR
    ledger["year_calendar_fixed"] = int(early.sum())
    ledger["year_after_snapshot"] = df.loc[future, ["name", "year_raw"]].to_dict("records")
    df.loc[future, "year"] = np.nan
    ledger["year_missing"] = int(df["year"].isna().sum())

    # Masses: a recorded mass of zero means "unknown".
    ledger["mass_zero"] = int((df["mass_g"] == 0).sum())
    ledger["mass_missing"] = int(df["mass_g"].isna().sum())
    df.loc[df["mass_g"] <= 0, "mass_g"] = np.nan

    # Families.
    fam = df["recclass"].map(classify_recclass)
    df["family"] = fam.str[0]
    df["group"] = fam.str[1]

    # Names -> collection areas.
    parsed = df["name"].map(parse_name)
    df["area"] = parsed.str[0]
    df["area_digits"] = parsed.str[1]

    # Coordinates, first pass.
    status = pd.Series("ok", index=df.index)
    status[df["lat"].isna() | df["lon"].isna()] = "missing"
    status[(df["lat"] == 0) & (df["lon"] == 0)] = "placeholder_0_0"
    status[df["name"].isin(FOUND_ON_MARS)] = "mars"
    df["coord_status"] = status

    # Where is each collection area? The median of its located members.
    ok = df["coord_status"].eq("ok")
    area_lat = df[ok & df["area"].notna()].groupby("area")["lat"].median()
    df["area_lat"] = df["area"].map(area_lat)

    # Environment: how each meteorite came to be in the catalogue.
    env = pd.Series("chance", index=df.index)
    numbered = df["area"].notna()
    env[numbered & df["area_lat"].ge(ANTARCTIC_LATITUDE)] = "desert"
    env[numbered & df["area_lat"].lt(ANTARCTIC_LATITUDE)] = "antarctic"
    # Areas with no located member at all: judge by their numbering, area by area
    # (an area is Antarctic if nearly all its numbers follow the season-year scheme).
    df["antarctic_number"] = [antarctic_style_number(d, y)
                              for d, y in zip(df["area_digits"], df["year"])]
    no_located = numbered & df["area_lat"].isna()
    share = df[no_located].groupby("area")["antarctic_number"].mean()
    yy = df["area"].map(share).ge(0.9).to_numpy()
    env[no_located & yy] = "antarctic"
    env[no_located & ~yy] = "unknown"
    env[df["area"].isin(DESERT_TRADE_NAMES)] = "desert"
    env[df["area"].isin(UNKNOWN_ORIGIN_NAMES)] = "unknown"
    env[~numbered & ok & df["lat"].lt(ANTARCTIC_LATITUDE)] = "antarctic"
    env[df["fell"]] = "fall"
    env[~df["valid"]] = "relict"
    env[df["coord_status"].eq("mars")] = "mars"
    df["environment"] = env
    ledger["environment_by_number_rule"] = int((no_located & yy & env.eq("antarctic")).sum())

    # Coordinates, second pass: members of Antarctic areas placed outside Antarctica.
    stray = df["environment"].eq("antarctic") & ok & df["lat"].ge(ANTARCTIC_LATITUDE)
    ledger["antarctic_area_outside_antarctica"] = (
        df[stray].groupby(["area", "lat", "lon"]).size().rename("n").reset_index()
        .to_dict("records"))
    df.loc[stray, "coord_status"] = "zeroed_latitude"
    df["located"] = df["coord_status"].eq("ok")

    for s in ("missing", "placeholder_0_0", "zeroed_latitude", "mars"):
        ledger[f"coords_{s}"] = int(df["coord_status"].eq(s).sum())
    ledger["located"] = int(df["located"].sum())
    return df, ledger


def attach_countries(df: pd.DataFrame, countries: list) -> pd.DataFrame:
    """Add the country each located meteorite lies in (-1 / NaN offshore)."""
    df = df.copy()
    df["country_idx"] = -1
    df["snap_km"] = np.nan
    loc = df["located"].to_numpy()
    idx, snap = geo.assign_countries(df.loc[loc, "lon"].to_numpy(), df.loc[loc, "lat"].to_numpy(),
                                     countries)
    df.loc[loc, "country_idx"] = idx
    df.loc[loc, "snap_km"] = snap
    names = np.array([c.name for c in countries] + [None], dtype=object)
    isos = np.array([c.iso3 for c in countries] + [None], dtype=object)
    df["country"] = names[df["country_idx"].to_numpy()]
    df["iso3"] = isos[df["country_idx"].to_numpy()]
    return df


def load_catalog(*, with_countries: bool = True):
    """The cleaned catalogue, the cleaning ledger, and the country list."""
    df, ledger = clean(load_raw())
    countries = geo.load_countries(COUNTRIES_50M)
    if with_countries:
        df = attach_countries(df, countries)
        loc = df["located"]
        ledger["offshore"] = int((loc & df["country_idx"].lt(0)).sum())
        ledger["snapped_to_coast"] = int((loc & df["snap_km"].gt(0)).sum())
    return df, ledger, countries
