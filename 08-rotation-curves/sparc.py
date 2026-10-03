"""Readers for the SPARC tables and the NASA Planetary Fact Sheet.

The SPARC files use the AAS "machine-readable table" (MRT) format: a header
that describes every column by byte position, then fixed-width rows. Astropy
can read MRT files -- but not these. In Table 1 the data rows sit one byte to
the right of the positions the header declares, so a reader that trusts the
byte positions cuts every number in the wrong place (astropy fails on the
distance column). Every field in these tables is a single token with no
embedded spaces, so we split on whitespace instead, and check the number of
fields in every row: a malformed row raises rather than being silently misread.

Units are converted once, here, on the way in: SPARC gives luminosities and
masses in units of 10^9 solar, and everything downstream uses plain solar units.
"""

from __future__ import annotations

import html
import re
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

import numpy as np

from config import DATA_DIR

HUBBLE_TYPES = {
    0: "S0", 1: "Sa", 2: "Sab", 3: "Sb", 4: "Sbc", 5: "Sc", 6: "Scd",
    7: "Sd", 8: "Sdm", 9: "Sm", 10: "Im", 11: "BCD",
}

DISTANCE_METHODS = {
    1: "Hubble flow",
    2: "tip of the red giant branch",
    3: "Cepheid variables",
    4: "Ursa Major cluster membership",
    5: "supernova light curve",
}


def _mrt_rows(path: Path, n_fields: int) -> list[list[str]]:
    """Whitespace-split data rows of an MRT file, checking every row's width.

    The data start after the last line of dashes, which closes the header.
    """
    lines = path.read_text().splitlines()
    rules = [i for i, line in enumerate(lines) if line.startswith("-----")]
    if not rules:
        raise ValueError(f"{path.name}: no header rule found; is this an MRT file?")
    rows = [line.split() for line in lines[rules[-1] + 1:] if line.strip()]
    bad = [r for r in rows if len(r) != n_fields]
    if bad:
        raise ValueError(f"{path.name}: {len(bad)} rows lack {n_fields} fields, e.g. {bad[0]}")
    return rows


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    return path


def display_name(name: str) -> str:
    """'NGC3198' -> 'NGC 3198', 'UGC00128' -> 'UGC 128'; others unchanged."""
    m = re.fullmatch(r"(NGC|UGCA|UGC|DDO|IC|PGC)0*(\d.*)", name)
    return f"{m.group(1)} {m.group(2)}" if m else name


# --------------------------------------------------------------------------
# The galaxies
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class RotationCurve:
    """One galaxy's rotation curve and the visible-matter components, by radius.

    The stellar components are for a mass-to-light ratio of 1 Msun/Lsun at 3.6
    microns; the gas already includes the 1.33 helium correction.
    """

    r: np.ndarray        #: galactocentric radius, kpc
    v_obs: np.ndarray    #: observed rotation speed, km/s
    e_v_obs: np.ndarray  #: its random uncertainty, km/s (excludes inclination)
    v_gas: np.ndarray    #: gas contribution, km/s (can be negative)
    v_disk: np.ndarray   #: stellar-disc contribution at Y = 1, km/s
    v_bulge: np.ndarray  #: bulge contribution at Y = 1, km/s
    sb_disk: np.ndarray  #: disc surface brightness, Lsun/pc^2
    sb_bulge: np.ndarray  #: bulge surface brightness, Lsun/pc^2

    def __len__(self) -> int:
        return len(self.r)


@dataclass(frozen=True)
class Galaxy:
    """One SPARC galaxy: the Table 1 properties plus its rotation curve."""

    name: str
    hubble_type: int
    distance_mpc: float
    distance_err_mpc: float
    distance_method: int
    inclination_deg: float
    inclination_err_deg: float
    #: Total 3.6-micron luminosity, Lsun.
    luminosity: float
    luminosity_err: float
    #: Half-light radius, kpc, and the mean surface brightness inside it, Lsun/pc^2.
    r_eff_kpc: float
    sb_eff: float
    #: Exponential disc scale length, kpc, and central surface brightness, Lsun/pc^2.
    r_disk_kpc: float
    sb_disk0: float
    #: Atomic hydrogen mass, Msun (hydrogen only -- no helium correction).
    m_hi: float
    r_hi_kpc: float
    #: Speed along the flat outer part of the rotation curve, km/s. ``None``
    #: where SPARC could not identify a flat part (the curve is still rising).
    v_flat: float | None
    v_flat_err: float | None
    #: 1 = high, 2 = medium, 3 = low quality rotation curve.
    quality: int
    references: str
    #: Bulge luminosity at 3.6 microns, Lsun (zero for bulgeless galaxies).
    luminosity_bulge: float
    curve: RotationCurve

    @property
    def label(self) -> str:
        return display_name(self.name)

    @property
    def type_name(self) -> str:
        return HUBBLE_TYPES.get(self.hubble_type, "?")

    @property
    def luminosity_disk(self) -> float:
        return max(self.luminosity - self.luminosity_bulge, 0.0)

    @property
    def has_bulge(self) -> bool:
        return bool(np.any(self.curve.v_bulge > 0))

    @cached_property
    def distance_method_name(self) -> str:
        return DISTANCE_METHODS.get(self.distance_method, "?")


def _load_bulges() -> dict[str, float]:
    path = _require(DATA_DIR / "Bulges.mrt")
    out = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        name, value = line.split()
        out[name] = float(value) * 1e9
    return out


def load_sparc() -> dict[str, Galaxy]:
    """Every SPARC galaxy, keyed by name, in catalogue order.

    Cross-checks the two tables against each other on the way in: every galaxy
    in Table 1 must have a rotation curve in Table 2, at the same distance.
    """
    t1 = _mrt_rows(_require(DATA_DIR / "SPARC_Lelli2016c.mrt"), 19)
    t2 = _mrt_rows(_require(DATA_DIR / "MassModels_Lelli2016c.mrt"), 10)
    bulges = _load_bulges()

    curves: dict[str, list[list[float]]] = {}
    distances: dict[str, float] = {}
    for row in t2:
        curves.setdefault(row[0], []).append([float(x) for x in row[2:]])
        distances.setdefault(row[0], float(row[1]))

    galaxies: dict[str, Galaxy] = {}
    for row in t1:
        name = row[0]
        if name not in curves:
            raise ValueError(f"{name} is in Table 1 but has no rotation curve in Table 2")
        d = float(row[2])
        if abs(distances[name] - d) > 0.01 * d:
            raise ValueError(f"{name}: Table 1 distance {d} Mpc disagrees with Table 2 "
                             f"{distances[name]} Mpc")
        a = np.array(curves[name])
        order = np.argsort(a[:, 0], kind="stable")
        a = a[order]
        v_flat = float(row[15])
        galaxies[name] = Galaxy(
            name=name,
            hubble_type=int(row[1]),
            distance_mpc=d,
            distance_err_mpc=float(row[3]),
            distance_method=int(row[4]),
            inclination_deg=float(row[5]),
            inclination_err_deg=float(row[6]),
            luminosity=float(row[7]) * 1e9,
            luminosity_err=float(row[8]) * 1e9,
            r_eff_kpc=float(row[9]),
            sb_eff=float(row[10]),
            r_disk_kpc=float(row[11]),
            sb_disk0=float(row[12]),
            m_hi=float(row[13]) * 1e9,
            r_hi_kpc=float(row[14]),
            v_flat=v_flat if v_flat > 0 else None,
            v_flat_err=float(row[16]) if v_flat > 0 else None,
            quality=int(row[17]),
            references=row[18],
            luminosity_bulge=bulges.get(name, 0.0),
            curve=RotationCurve(*(a[:, k] for k in range(8))),
        )

    orphans = set(curves) - set(galaxies)
    if orphans:
        raise ValueError(f"rotation curves without a Table 1 entry: {sorted(orphans)}")
    return galaxies


# --------------------------------------------------------------------------
# Published results, for comparison only
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class PublishedRAR:
    """Lelli et al. (2017), data behind their Figure 2: log10 accelerations, m/s^2."""

    log_gbar: np.ndarray
    e_log_gbar: np.ndarray
    log_gobs: np.ndarray
    e_log_gobs: np.ndarray


@dataclass(frozen=True)
class PublishedRARBins:
    log_gbar: np.ndarray
    log_gobs: np.ndarray
    sd: np.ndarray
    n: np.ndarray


def load_published_rar() -> tuple[PublishedRAR, PublishedRARBins]:
    pts = np.array(_mrt_rows(_require(DATA_DIR / "RAR.mrt"), 4), dtype=float)
    bins = np.array(_mrt_rows(_require(DATA_DIR / "RARbins.mrt"), 4), dtype=float)
    return (
        PublishedRAR(*(pts[:, k] for k in range(4))),
        PublishedRARBins(bins[:, 0], bins[:, 1], bins[:, 2], bins[:, 3].astype(int)),
    )


def load_published_btfr() -> dict[str, dict]:
    """Lelli et al. (2019): baryonic mass and flat velocity per galaxy.

    Their masses assume Y_disk = 0.5 at 3.6 microns, as ours do.
    """
    out = {}
    for row in _mrt_rows(_require(DATA_DIR / "BTFR_Lelli2019.mrt"), 19):
        v_f = float(row[5])
        out[row[0]] = {
            "log_mbar": float(row[1]),
            "e_log_mbar": float(row[2]),
            "v_flat": v_f if v_f > 0 else None,
            "e_v_flat": float(row[6]) if v_f > 0 else None,
        }
    return out


def load_published_max_disk() -> dict[str, dict]:
    """Starkman et al. (2018): maximum-disc stellar mass-to-light ratios.

    This file's byte-position header is malformed, but it carries a line of
    column names, so we read by name from that line.
    """
    path = _require(DATA_DIR / "MaximumDiskData.mrt")
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("ID "))
    columns = lines[start].split()
    out = {}
    for line in lines[start + 1:]:
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != len(columns):
            raise ValueError(f"{path.name}: row has {len(fields)} fields, expected {len(columns)}")
        rec = dict(zip(columns, fields))
        out[rec["ID"]] = {
            "upsilon_disk": float(rec["M2L"]),
            "e_upsilon_disk": float(rec["e_M2L"]),
            "upsilon_bulge": float(rec["M2Lb"]),
        }
    return out


# --------------------------------------------------------------------------
# The Solar System
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Planet:
    name: str
    #: Semi-major axis ("distance from Sun" on the fact sheet), km.
    a_km: float
    eccentricity: float
    #: Mean orbital speed, km/s, as published (rounded to 0.1 km/s).
    v_mean_kms: float
    period_days: float
    #: Pluto is on the fact sheet but is a dwarf planet.
    dwarf: bool = False


def load_planets() -> list[Planet]:
    """Planets from the NASA Planetary Fact Sheet (the Moon is skipped).

    Parses the page's single HTML table by row label. The Moon's column holds
    Earth-relative values (flagged with an asterisk), so it is dropped.
    """
    text = _require(DATA_DIR / "planetary_factsheet.html").read_text(encoding="latin-1")
    table: dict[str, list[str]] = {}
    header: list[str] | None = None
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", text, re.S | re.I):
        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row, re.S | re.I)
        cells = [" ".join(html.unescape(re.sub(r"<[^>]+>", "", c)).split()) for c in cells]
        if not cells:
            continue
        if header is None and "MERCURY" in cells:
            header = cells[1:]
        elif cells[0]:
            table[cells[0]] = cells[1:]
    if header is None:
        raise ValueError("fact sheet: planet header row not found")

    def row(prefix: str) -> list[str]:
        key = next((k for k in table if k.startswith(prefix)), None)
        if key is None:
            raise ValueError(f"fact sheet: no row starting {prefix!r}")
        return table[key]

    dist, vel, ecc, per = (row(p) for p in (
        "Distance from Sun", "Orbital Velocity", "Orbital Eccentricity", "Orbital Period"))
    planets = []
    for i, name in enumerate(header):
        if name == "MOON":
            continue
        num = lambda s: float(s.replace(",", "").rstrip("*"))  # noqa: E731
        planets.append(Planet(
            name=name.title(),
            a_km=num(dist[i]) * 1e6,
            eccentricity=num(ecc[i]),
            v_mean_kms=num(vel[i]),
            period_days=num(per[i]),
            dwarf=(name == "PLUTO"),
        ))
    return planets
