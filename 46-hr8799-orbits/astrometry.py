"""Relative astrometry of HR 8799 b, c, d, e, parsed from the papers' LaTeX.

Every measurement becomes one row of a tidy table (``data/derived/astrometry.csv``)
with the planet, epoch, the instrument/paper it comes from ("source"), and
either an offset from the star (Delta RA = east, Delta Dec = north, mas) or a
separation and position angle (mas, degrees east of north), exactly as the
paper gave it. Nothing is rounded or converted on the way in; the likelihood
uses each measurement in its own coordinates.

Conventions checked here once (and tested in test_physics.py):

* Konopacky et al. (2016) tabulate detector offsets Delta x, Delta y in arcsec,
  with Delta x positive towards the west: Delta RA = -Delta x. (Zurlo et al.
  2022 list the same points with the sign flipped, which confirms it.)
* Balmer et al. (2025) give arcsec; everyone else mas.
* The GRAVITY point of planet e has correlated errors (covariance -0.00668 mas^2).
"""

from __future__ import annotations

import re
from datetime import datetime

import numpy as np
import pandas as pd
from astropy.time import Time

import config as C

COLUMNS = ["source", "instrument", "reference", "planet", "date", "mjd", "kind",
           "ra", "ra_err", "dec", "dec_err", "corr", "sep", "sep_err", "pa", "pa_err", "note"]

_PM = re.compile(r"(-?\d+(?:\.\d+)?)\s*(?:\$?\s*\\pm\s*\$?)\s*(\d+(?:\.\d+)?)")


def _tex(key: str, filename: str | None = None) -> str:
    src = C.ARXIV[key]
    return (src.folder / (filename or src.files[0])).read_text(errors="replace")


def _between(text: str, start: str, end: str, after: str | None = None) -> str:
    i0 = text.index(after) if after else 0
    i = text.index(start, i0) + len(start)
    j = text.index(end, i)
    return text[i:j]


def _pm(s: str) -> tuple[float, float]:
    s = re.sub(r"\\phantom\{[^}]*\}", "", s)       # alignment padding in Nguyen+2026
    s = s.replace("$", "").replace("{", "").replace("}", "")
    m = _PM.search(s)
    if not m:
        raise ValueError(f"no value +- error in {s!r}")
    return float(m.group(1)), float(m.group(2))


def _mjd_from_date(d: str) -> float:
    for fmt in ("%Y %b %d", "%Y-%m-%d"):
        try:
            return Time(datetime.strptime(d.strip(), fmt)).mjd
        except ValueError:
            pass
    raise ValueError(d)


def _iso(d: str) -> str:
    for fmt in ("%Y %b %d", "%Y-%m-%d"):
        try:
            return datetime.strptime(d.strip(), fmt).strftime("%Y-%m-%d")
        except ValueError:
            pass
    raise ValueError(d)


def _row(source, instrument, reference, planet, date, mjd, kind, **kw) -> dict:
    row = dict.fromkeys(COLUMNS, np.nan)
    row.update(source=source, instrument=instrument, reference=reference, planet=planet,
               date=date, mjd=float(mjd), kind=kind, note="")
    if kind == "radec":
        row["corr"] = 0.0
    row.update(kw)
    return row


# ---------------------------------------------------------------------------
# One parser per paper
# ---------------------------------------------------------------------------

def konopacky2016() -> list[dict]:
    """Keck/NIRC2, Konopacky et al. 2016 Table 2 (arcsec, Delta x = -Delta RA)."""
    body = _between(_tex("Konopacky+2016", "tab2.tex"), r"\startdata", r"\enddata")
    rows = []
    for line in body.split(r"\\"):
        line = re.sub(r"\\tablenotemark\{\w\}", "", line).strip()
        if not line:
            continue
        cells = [c.strip() for c in line.split("&")]
        date = _iso(cells[0])
        # cells: date, b_x, b_y, '', c_x, c_y, '', d_x, d_y, '', e_x, e_y
        for planet, (ix, iy) in zip(C.PLANETS, [(1, 2), (4, 5), (7, 8), (10, 11)]):
            if "n/a" in cells[ix]:
                continue
            x, xe = _pm(cells[ix])
            y, ye = _pm(cells[iy])
            r = _row("NIRC2-K16", "Keck/NIRC2", "Konopacky+2016", planet, date,
                     _mjd_from_date(date), "radec",
                     ra=-1000 * x, ra_err=1000 * xe, dec=1000 * y, dec_err=1000 * ye)
            if (date, planet) in C.KONOPACKY_EXCLUDE:
                r["note"] = "excluded by Konopacky+2016 (table notes)"
            rows.append(r)
    return rows


def wang2018_gpi() -> list[dict]:
    """GPI separations and position angles, Wang et al. 2018 Table 2."""
    tex = _tex("Wang+2018")
    body = _between(tex, r"\startdata", r"\enddata", after=r"\label{table:astrom}")
    rows, date = [], None
    for line in body.split(r"\\"):
        cells = [c.strip() for c in line.split("&")]
        if len(cells) < 6:
            continue
        if cells[0]:
            date = _iso(cells[0])
        sep, sep_e = _pm(cells[4])
        pa, pa_e = _pm(cells[5])
        rows.append(_row("GPI-W18", "Gemini/GPI", "Wang+2018", cells[1], date,
                         _mjd_from_date(date), "seppa",
                         sep=sep, sep_err=sep_e, pa=pa, pa_err=pa_e,
                         note="original GPI calibration (north offset -0.10 deg)"))
    return rows


def derosa2020_gpi() -> list[dict]:
    """The same GPI points with the revised calibration, De Rosa et al. 2020 Table 7."""
    tex = _tex("DeRosa+2020")
    body = _between(tex, r"\hline", r"\end{tabular}", after=r"\label{tbl:revised_astrometry}")
    rows = []
    for line in body.split(r"\\"):
        cells = [c.strip() for c in line.split("&")]
        if len(cells) < 7 or "HR 8799" not in cells[0]:
            continue
        planet = cells[0].split()[-1]
        sep, sep_e = _pm(cells[5])
        pa, pa_e = _pm(cells[6])
        rows.append(_row("GPI-D20", "Gemini/GPI", "DeRosa+2020", planet, cells[1],
                         _mjd_from_date(cells[1]), "seppa",
                         sep=sep, sep_err=sep_e, pa=pa, pa_err=pa_e,
                         note="revised GPI calibration (north offset +0.17 to +0.45 deg)"))
    return rows


def _zurlo_table(label: str, planets: tuple[str, ...], source: str, instrument: str) -> list[dict]:
    tex = _tex("Zurlo+2022")
    body = _between(tex, r"\label{" + label + "}", r"\end{tabular}")
    rows = []
    for line in body.split(r"\\"):
        cells = [c.strip() for c in line.replace(r"\hline", "").split("&")]
        if len(cells) != len(planets) + 1 or not re.match(r"\d{4}-\d\d-\d\d", cells[0]):
            continue
        date = cells[0]
        for planet, cell in zip(planets, cells[1:]):
            ra_s, dec_s = cell.split(";")
            ra, ra_e = _pm(ra_s)
            dec, dec_e = _pm(dec_s)
            rows.append(_row(source, instrument, "Zurlo+2022", planet, date,
                             _mjd_from_date(date), "radec",
                             ra=ra, ra_err=ra_e, dec=dec, dec_err=dec_e,
                             note="values and errors rounded to 1 mas in the paper"))
    return rows


def zurlo2022() -> list[dict]:
    return (_zurlo_table("t:astroirdis", C.PLANETS, "IRDIS-Z22", "VLT/SPHERE-IRDIS")
            + _zurlo_table("t:astroifs", ("d", "e"), "IFS-Z22", "VLT/SPHERE-IFS")
            + _zurlo_table("t:luci", C.PLANETS, "LUCI-Z22", "LBT/LUCI"))


def soummer2011() -> list[dict]:
    """HST/NICMOS 1998; the 'sigma' column (incl. 0.5 deg north) is the error used."""
    tex = _tex("Soummer+2011")
    body = _between(tex, r"\tableline", r"\end{tabular}", after=r"\label{tableResults}")
    epoch = lafreniere2009_epoch()
    rows = []
    for line in body.split(r"\\"):
        cells = [c.strip() for c in line.split("&")]
        if len(cells) < 5 or not cells[0].startswith("This work"):
            continue
        planet = cells[1][-1]
        sep, _ = _pm(cells[2])
        pa, _ = _pm(cells[3])
        sigma = float(cells[4])
        rows.append(_row("NICMOS-S11", "HST/NICMOS", "Soummer+2011", planet,
                         f"{epoch:.4f}", Time(epoch, format="decimalyear").mjd, "seppa",
                         sep=sep, sep_err=sigma, pa=pa, pa_err=np.degrees(sigma / sep),
                         note="isotropic error incl. 0.5 deg north systematics"))
    return rows


def lafreniere2009_epoch() -> float:
    tex = _tex("Lafreniere+2009")
    m = re.search(r"\n(1998\.\d+)\s*&", tex)
    return float(m.group(1))


def gravity2019() -> list[dict]:
    tex = _tex("GRAVITY+2019")
    body = _between(tex, r"\hline", r"\end{tabular}", after=r"\label{table:2}")
    mjds = [float(m) for m in re.findall(r"\n\s*(58358\.\d+)\s*&", body)]
    g = re.search(r"Global\s*&\s*\$(-?[\d.]+)\\pm([\d.]+)\$\s*&\s*\$(-?[\d.]+)\\pm([\d.]+)\$\s*&\s*(-?[\d.]+)",
                  body)
    ra, ra_e, dec, dec_e, cov = map(float, g.groups())
    mjd = float(np.mean(mjds))
    return [_row("GRAVITY-L19", "VLTI/GRAVITY", "GRAVITY+2019", "e",
                 Time(mjd, format="mjd").iso[:10], mjd, "radec",
                 ra=ra, ra_err=ra_e, dec=dec, dec_err=dec_e, corr=cov / (ra_e * dec_e),
                 note=f"mean of {len(mjds)} exposures")]


def balmer2025() -> list[dict]:
    """JWST/NIRCam GTO 1194 (2023 Nov 5), the 'all filters' rows (arcsec)."""
    tex = _tex("Balmer+2025")
    body = _between(tex, r"\startdata", r"\enddata", after=r"\label{tab:bka}")
    rows = []
    for line in body.split(r"\\"):
        if r"\cdots & \cdots &" not in line or " all " not in line:
            continue
        m = re.match(r"\s*\\hr(\w)\s*&\s*all", line)
        if not m:
            continue
        cells = [c.strip() for c in line.split("&")]
        ra, ra_e = _pm(cells[4])
        dec, dec_e = _pm(cells[5])
        rows.append(_row("NIRCam-B25", "JWST/NIRCam", "Balmer+2025", m.group(1), "2023-11-05",
                         _mjd_from_date("2023-11-05"), "radec",
                         ra=1000 * ra, ra_err=1000 * ra_e, dec=1000 * dec, dec_err=1000 * dec_e,
                         note="median and scatter over filters"))
    return rows


def nguyen2026() -> list[dict]:
    tex = _tex("Nguyen+2026")
    rows = []
    for label, source, instrument, date in [
        ("tab:astrophot", "AMI-N26", "JWST/NIRISS-AMI", "2023-08-03"),
        ("tab:kecknirc2", "NIRC2-N26", "Keck/NIRC2", "2025-10-09"),
    ]:
        body = _between(tex, r"\startdata", r"\enddata", after=r"\label{" + label + "}")
        for line in body.split(r"\\"):
            cells = [c.strip() for c in line.replace(r"\midrule", "").split("&")]
            if len(cells) < 3 or cells[0] not in C.PLANETS + ("f",):
                continue
            ra, ra_e = _pm(cells[1])
            dec, dec_e = _pm(cells[2])
            rows.append(_row(source, instrument, "Nguyen+2026", cells[0], date,
                             _mjd_from_date(date), "radec",
                             ra=ra, ra_err=ra_e, dec=dec, dec_err=dec_e,
                             note="candidate fifth planet" if cells[0] == "f" else ""))
    return rows


# ---------------------------------------------------------------------------
# The assembled table
# ---------------------------------------------------------------------------

PARSERS = [konopacky2016, wang2018_gpi, derosa2020_gpi, zurlo2022, soummer2011,
           gravity2019, balmer2025, nguyen2026]


def build() -> pd.DataFrame:
    rows = [r for parse in PARSERS for r in parse()]
    df = pd.DataFrame(rows, columns=COLUMNS)
    # fill the other coordinate pair for plotting (not used by the likelihood)
    sp = df.kind == "seppa"
    pa = np.radians(df.loc[sp, "pa"])
    df.loc[sp, "ra"] = df.loc[sp, "sep"] * np.sin(pa)
    df.loc[sp, "dec"] = df.loc[sp, "sep"] * np.cos(pa)
    df.loc[sp, "ra_err"] = np.hypot(df.loc[sp, "sep_err"] * np.sin(pa),
                                    df.loc[sp, "sep"] * np.radians(df.loc[sp, "pa_err"]) * np.cos(pa))
    df.loc[sp, "dec_err"] = np.hypot(df.loc[sp, "sep_err"] * np.cos(pa),
                                     df.loc[sp, "sep"] * np.radians(df.loc[sp, "pa_err"]) * np.sin(pa))
    rd = df.kind == "radec"
    df.loc[rd, "sep"] = np.hypot(df.loc[rd, "ra"], df.loc[rd, "dec"])
    df.loc[rd, "pa"] = np.degrees(np.arctan2(df.loc[rd, "ra"], df.loc[rd, "dec"])) % 360
    df["year"] = Time(df.mjd.to_numpy(), format="mjd").decimalyear
    df = df.sort_values(["planet", "mjd", "source"]).reset_index(drop=True)
    _check(df)
    return df


def _check(df: pd.DataFrame) -> None:
    """Spot checks against the papers and against Zurlo et al.'s compilation."""
    n = df.groupby("source").size().to_dict()
    expected = {"NIRC2-K16": 46 + 4 + 1, "GPI-W18": 9, "GPI-D20": 9, "IRDIS-Z22": 17 * 4,
                "IFS-Z22": 14 * 2, "LUCI-Z22": 4, "NICMOS-S11": 3, "GRAVITY-L19": 1,
                "NIRCam-B25": 4, "AMI-N26": 5, "NIRC2-N26": 4}
    if n != expected:
        raise AssertionError(f"row counts {n} != {expected}")
    k = df[(df.source == "NIRC2-K16") & (df.date == "2004-07-14") & (df.planet == "b")].iloc[0]
    assert (k.ra, k.dec) == (1471.0, 884.0), "Konopacky sign convention"   # Zurlo: 1471, 884
    g = df[(df.source == "GPI-D20") & (df.planet == "c") & (df.date == "2013-11-17")].iloc[0]
    assert abs(g.ra - (-537)) < 1 and abs(g.dec - 782) < 1, "De Rosa revision vs Zurlo"
    h = df[df.source == "NICMOS-S11"].set_index("planet")
    assert abs(h.loc["b", "ra"] - 1418) < 1 and abs(h.loc["d", "dec"] + 533) < 1


def load(sources=None, *, include_excluded: bool = False, before_mjd: float | None = None,
         planets=C.PLANETS) -> pd.DataFrame:
    """The astrometry table, optionally restricted to some sources / epochs."""
    path = C.DERIVED / "astrometry.csv"
    if not path.exists():
        build().to_csv(path, index=False)
    df = pd.read_csv(path, keep_default_na=True)
    df["note"] = df["note"].fillna("")
    if sources is not None:
        df = df[df.source.isin(sources)]
    if not include_excluded:
        df = df[~df.note.str.startswith(("excluded", "discarded"))]
    if before_mjd is not None:
        df = df[df.mjd < before_mjd]
    df = df[df.planet.isin(planets)]
    return df.reset_index(drop=True)


def dataset(name: str) -> pd.DataFrame:
    if name == "W18":
        return load(C.W18_SOURCES)
    if name == "ALL":
        return load(C.ALL_SOURCES)
    raise KeyError(name)


if __name__ == "__main__":
    df = build()
    df.to_csv(C.DERIVED / "astrometry.csv", index=False)
    print(df.groupby(["source", "planet"]).size().unstack(fill_value=0))
    for name in ("W18", "ALL"):
        d = dataset(name)
        print(f"{name}: {len(d)} measurements ({2 * len(d)} numbers), "
              f"{d.year.min():.2f}-{d.year.max():.2f}")
