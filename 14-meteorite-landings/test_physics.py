"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework needed, and no network:
the geometry, parsing and statistics are checked against analytic results and
textbook cases; the last section also runs over the cached catalogue, if present.

Every map, every "finds per million km^2" and every recovery rate passes through
these functions. An equal-area projection that is 3% off would not look wrong
on a world map -- which is exactly why it is checked here.
"""

from __future__ import annotations

import sys

import numpy as np

import analysis
import catalog
import geo
from config import COUNTRIES_50M, EARTH_RADIUS_KM, METEORITE_CSV

FAILURES: list[str] = []
N_CHECKS = 0


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    global N_CHECKS
    N_CHECKS += 1
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    status = "ok  " if rel <= tol else "FAIL"
    if rel > tol:
        FAILURES.append(label)
    print(f"  [{status}] {label:<52} {got:>12.6g} vs {expected:>12.6g} {unit:<7} "
          f"({rel * 100:.2g}% off, tol {tol * 100:.2g}%)")


def check_true(label: str, condition: bool) -> None:
    global N_CHECKS
    N_CHECKS += 1
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def jacobian_det(proj, lon, lat, h=1e-5):
    """|d(x,y)/d(lambda,phi)| by central differences, in radians."""
    x1, y1 = proj(lon + h, lat)
    x0, y0 = proj(lon - h, lat)
    dxl, dyl = (x1 - x0) / np.radians(2 * h), (y1 - y0) / np.radians(2 * h)
    x1, y1 = proj(lon, lat + h)
    x0, y0 = proj(lon, lat - h)
    dxp, dyp = (x1 - x0) / np.radians(2 * h), (y1 - y0) / np.radians(2 * h)
    return np.abs(dxl * dyp - dxp * dyl)


# --- Calendars, for the year repair ---------------------------------------------------

def jdn_gregorian(y: int, m: int, d: int) -> int:
    """Julian Day Number of a proleptic-Gregorian date (Fliegel & Van Flandern 1968)."""
    a = (14 - m) // 12
    yy, mm = y + 4800 - a, m + 12 * a - 3
    return d + (153 * mm + 2) // 5 + 365 * yy + yy // 4 - yy // 100 + yy // 400 - 32045


def julian_calendar_date(jdn: int) -> tuple[int, int, int]:
    """The Julian-calendar date of a Julian Day Number."""
    c = jdn + 32082
    d = (4 * c + 3) // 1461
    e = c - 1461 * d // 4
    m = (5 * e + 2) // 153
    return d - 4800 + m // 10, m + 3 - 12 * (m // 10), e - (153 * m + 2) // 5 + 1


def hybrid_year_of_jan1(year: int) -> int:
    """The year you read back if 1 January of ``year`` (proleptic Gregorian) is displayed
    in a hybrid calendar -- Julian before 15 October 1582, Gregorian after -- as Java's
    GregorianCalendar and many databases do by default."""
    jdn = jdn_gregorian(year, 1, 1)
    return julian_calendar_date(jdn)[0] if jdn < jdn_gregorian(1582, 10, 15) else year


def main() -> int:
    print("Equal Earth projection (Savric, Patterson & Jenny 2018)")
    a1 = 1.340264
    check("half-width x(180 E, 0)", float(geo.equal_earth(180, 0)[0]),
          2 * np.sqrt(3) * np.pi / (3 * a1), 1e-12)
    ymax = float(geo.equal_earth(0, 90)[1])
    th = np.pi / 3  # the parametric latitude of the pole: asin(sqrt(3)/2)
    check("half-height y(0, 90 N)", ymax,
          th * (a1 - 0.081106 * th**2 + th**6 * (0.000893 + 0.003796 * th**2)), 1e-12)
    check("width : height (paper: 2.05 : 1)", float(geo.equal_earth(180, 0)[0]) / ymax, 2.05,
          0.005)
    rng = np.random.default_rng(1)
    lon, lat = rng.uniform(-179, 179, 400), rng.uniform(-89, 89, 400)
    ratio = jacobian_det(geo.equal_earth, lon, lat) / np.cos(np.radians(lat))
    check("equal-area: |Jacobian| / cos(lat), worst of 400", float(np.max(np.abs(ratio - 1))) + 1,
          1.0, 1e-6)
    ox, oy = geo.equal_earth(*geo.world_outline(n=20001))
    check("area inside the map outline = 4 pi", abs(geo.shoelace(ox, oy)), 4 * np.pi, 1e-6)
    blon, blat = geo.equal_earth_inverse(*geo.equal_earth(lon, lat))
    check_true("inverse projection round trip < 1e-9 degree",
               bool(np.max(np.abs(blon - lon)) < 1e-9 and np.max(np.abs(blat - lat)) < 1e-9))

    print("\nSouth polar Lambert azimuthal equal-area")
    lat_s = rng.uniform(-89, -50, 300)
    ratio = jacobian_det(geo.south_polar, lon[:300], lat_s) / np.cos(np.radians(lat_s))
    check("equal-area: |Jacobian| / cos(lat), worst of 300", float(np.max(np.abs(ratio - 1))) + 1,
          1.0, 1e-6)
    ring = np.linspace(-180, 180, 20001)
    cap = abs(geo.shoelace(*geo.south_polar(ring, np.full_like(ring, -60.0))))
    check("cap south of 60 S = 2 pi (1 - sin 60)", cap, 2 * np.pi * (1 - np.sin(np.radians(60))),
          1e-6)
    check_true("Greenwich meridian points up, 90 E right",
               bool(geo.south_polar(0, -70)[1] > 0 and geo.south_polar(90, -70)[0] > 0))

    print("\nAreas on the sphere")
    check("Earth's surface, 4 pi R^2 (million km^2)", geo.band_area_km2(-90, 90) / 1e6, 510.066,
          1e-4, "M km2")
    cell = np.array([[10, 0], [11, 0], [11, 1], [10, 1], [10, 0]], dtype=float)
    check("1 x 1 degree cell at the equator", geo.ring_area_km2(*cell.T),
          geo.band_area_km2(0, 1, 1), 1e-4, "km2")
    cell = np.array([[10, 60], [11, 60], [11, 61], [10, 61], [10, 60]], dtype=float)
    check("1 x 1 degree cell at 60 N", geo.ring_area_km2(*cell.T),
          geo.band_area_km2(60, 61, 1), 1e-4, "km2")
    lo, la = geo.destination(0.0, 0.0, 0.0, np.radians(1) * EARTH_RADIUS_KM)
    check("1 degree of arc north from (0, 0)", float(la), 1.0, 1e-9, "deg")

    print("\nCountries (Natural Earth 1:50m, if cached)")
    if COUNTRIES_50M.exists():
        countries = geo.load_countries(COUNTRIES_50M)
        by = {c.name: c for c in countries}
        # Areas against published land areas: close, not exact -- Natural Earth
        # leaves out large lakes and ice shelves.
        check("Australia (published 7.69 M km^2)", by["Australia"].area_km2 / 1e6, 7.692, 0.02)
        check("Oman (published 309,500 km^2)", by["Oman"].area_km2, 309_500, 0.02)
        check("Germany (published 357,600 km^2)", by["Germany"].area_km2, 357_600, 0.02)
        places = {"Paris": (2.35, 48.86, "France"), "Muscat": (58.4, 23.6, "Oman"),
                  "Canberra": (149.13, -35.28, "Australia"),
                  "Honolulu": (-157.86, 21.31, "United States of America"),
                  "McMurdo Station": (166.67, -77.85, "Antarctica"),
                  "mid-Pacific": (-140.0, 0.0, None), "(0, 0) placeholder": (0.0, 0.0, None)}
        idx, _ = geo.assign_countries([p[0] for p in places.values()],
                                      [p[1] for p in places.values()], countries)
        for (label, (_, _, want)), i in zip(places.items(), idx):
            got = countries[i].name if i >= 0 else None
            check_true(f"{label} -> {want or 'the sea'}", got == want)
    else:
        print("  (skipped: run fetch_data.py first)")

    print("\nMeteorite classes -> families")
    cases = {
        "L6": ("ordinary", "L"), "H5": ("ordinary", "H"), "LL3.8": ("ordinary", "LL"),
        "H/L3.6": ("ordinary", "H/L"), "L(LL)3": ("ordinary", "L/LL"), "OC": ("ordinary", "OC"),
        "H-metal": ("ordinary", "H"), "Howardite": ("achondrite", "HED"),
        "Lodranite": ("achondrite", "Primitive achondrite"), "Lunar (anorth)": ("achondrite", "Lunar"),
        "Martian (shergottite)": ("achondrite", "Martian"), "Eucrite-pmict": ("achondrite", "HED"),
        "Enst achon": ("achondrite", "Aubrite"), "CM2": ("carbonaceous", "CM"),
        "CV3": ("carbonaceous", "CV"), "C2-ung": ("carbonaceous", "C (ungrouped)"),
        "C": ("carbonaceous", "C (ungrouped)"), "CH/CBb": ("carbonaceous", "CH"),
        "EH4": ("other_chondrite", "EH"), "E3": ("other_chondrite", "E"),
        "EL-melt rock": ("other_chondrite", "EL"), "R3.8": ("other_chondrite", "R"),
        "K3": ("other_chondrite", "K"), "Chondrite-ung": ("other_chondrite", "Chondrite (ungrouped)"),
        "Iron, IIIAB": ("iron", "IIIAB"), "Iron?": ("iron", "Iron (ungrouped or unstudied)"),
        "Pallasite, PMG": ("stony_iron", "Pallasite"), "Mesosiderite-A1": ("stony_iron", "Mesosiderite"),
        "Stone-uncl": ("unclassified", "Stone-uncl"), "Relict OC": ("unclassified", "Relict"),
    }
    bad = [k for k, v in cases.items() if catalog.classify_recclass(k) != v]
    check_true(f"{len(cases)} class strings classified as expected" + (f" (wrong: {bad})" if bad else ""),
               not bad)

    print("\nNames -> collection areas")
    names = {
        "Allan Hills 84001": ("Allan Hills", "84001"), "Allan Hills A77005": ("Allan Hills", "77005"),
        "Northwest Africa 7034": ("Northwest Africa", "7034"), "Dhofar 019": ("Dhofar", "019"),
        "Yamato 791197": ("Yamato", "791197"), "Queen Alexandra Range 93148": (
            "Queen Alexandra Range", "93148"),
        "Willamette": (None, None), "Campo del Cielo": (None, None), "Sutter's Mill": (None, None),
        "Niger (L6)": (None, None),
    }
    bad = [k for k, v in names.items() if catalog.parse_name(k) != v]
    check_true(f"{len(names)} names parsed as expected" + (f" (wrong: {bad})" if bad else ""),
               not bad)
    rules = {("84001", 1984): True, ("04315", 2004): True, ("04075", 2003): True,
             ("791197", 1979): True, ("7034", 2011): False, ("019", 1999): False,
             ("400", 1998): False, ("84001", 1990): False}
    bad = [k for k, v in rules.items() if catalog.antarctic_style_number(*k) != v]
    check_true("season-year numbering rule on 8 cases" + (f" (wrong: {bad})" if bad else ""),
               not bad)

    print("\nThe year repair: a Gregorian 1 January read back in the Julian calendar")
    # Documented dates vs NASA's file: Nogata fell 19 May 861; the Elbogen iron
    # fell about 1400; Ensisheim fell 7 November 1492; Campo del Cielo was found
    # by a Spanish expedition in 1576. NASA lists 860, 1399, 1491 and 1575.
    for name, true_year, file_year in [("Nogata", 861, 860), ("Elbogen", 1400, 1399),
                                       ("Ensisheim", 1492, 1491), ("Campo del Cielo", 1576, 1575)]:
        check_true(f"{name}: {true_year} reads back as {file_year}",
                   hybrid_year_of_jan1(true_year) == file_year)
    for name, year in [("Castrovillari", 1583), ("L'Aigle", 1803)]:
        check_true(f"{name}: {year}, after the reform, is unchanged",
                   hybrid_year_of_jan1(year) == year)
    check_true("shift ends exactly at the reform (1582 -> 1581, 1583 -> 1583)",
               hybrid_year_of_jan1(1582) == 1581 and hybrid_year_of_jan1(1583) == 1583)

    print("\nStatistics")
    p = analysis.proportion(0, 798)
    check("Wilson upper limit, 0 of 798", p["ci95"][1], 1.959964**2 / (798 + 1.959964**2), 1e-4)
    check("Poisson (Garwood) upper limit, k = 0", analysis.poisson_ci(0)[1], 3.6889, 1e-4)
    lo, hi = analysis.poisson_ci(10)
    check("Poisson (Garwood) 95% interval, k = 10: lower", lo, 4.7954, 1e-4)
    check("Poisson (Garwood) 95% interval, k = 10: upper", hi, 18.3904, 1e-4)
    expected = np.array([50.0, 200.0, 1000.0, 5000.0, 20000.0])
    density = np.array([2.0, 10.0, 60.0, 250.0, 900.0])
    n = expected * 0.002 * (density / 100) ** 0.7
    fit = analysis.fit_efficiency(n, expected, density)
    check("Poisson fit recovers gamma from exact data", fit["gamma"], 0.7, 1e-5)
    check("Poisson fit recovers eps0 from exact data", fit["eps0"], 0.002, 1e-5)

    print("\nPublished flux, cross-checked")
    # Halliday et al.: ~9 falls per 10^6 km^2 per year drop >= 1 kg. Over the whole
    # Earth that is the "about 4500 falls a year" quoted by Borovicka et al. (2015).
    check("9 per 10^6 km^2/yr x Earth's surface", 9 * geo.EARTH_AREA_KM2 / 1e6, 4500, 0.03, "/yr")

    if METEORITE_CSV.exists():
        print("\nThe cached catalogue")
        raw = catalog.load_raw()
        classified = raw["recclass"].map(catalog.classify_recclass)
        check_true(f"all {raw['recclass'].nunique()} class strings in the file are recognised",
                   classified.notna().all())
        df, ledger = catalog.clean(raw)
        early = df[df["year_calendar_fixed"]]
        check_true(f"every repaired year ({len(early)}) matches the calendar explanation",
                   all(hybrid_year_of_jan1(int(y)) == int(r)
                       for y, r in zip(early["year"], early["year_raw"])))
        check_true("no Antarctic-area meteorite left outside Antarctica",
                   not (df["located"] & df["environment"].eq("antarctic")
                        & df["lat"].ge(-60)).any())

    print()
    if FAILURES:
        print(f"{len(FAILURES)} of {N_CHECKS} CHECKS FAILED: {', '.join(FAILURES)}")
        return 1
    print(f"All {N_CHECKS} checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
