"""Stellar properties along the main sequence: mass <-> magnitude <-> colour.

Everything here comes from Eric Mamajek's "Modern Mean Dwarf Stellar Color and
Effective Temperature Sequence" (Pecaut & Mamajek 2013, updated online; the
copy we use is version 2022.04.16). It gives, spectral type by spectral type,
the Gaia absolute magnitude M_G, the colours BP-RP and G-RP, and a typical mass.

We need it for three things:

* a mass for each star in the 100 pc sample (from its colour, which an
  unresolved companion changes far less than it changes the magnitude);
* the brightness of a companion of a given mass, which decides whether Gaia can
  see the binary (a faint companion barely moves the light centre; a twin
  moves it not at all);
* masses of both stars in a resolved wide pair (from each star's own M_G).

The table is a mean sequence for solar-metallicity dwarfs: individual stars
scatter about it by ~0.1-0.3 mag in M_G (age, metallicity), i.e. ~3-8% in mass.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

import config as C

MAMAJEK_PATH = C.RAW_DIR / "mamajek_dwarf_sequence.txt"


def _value(token: str) -> float:
    token = token.rstrip(":")  # a trailing colon marks an uncertain value
    if token in ("...", "....", ".....", "") or set(token) <= {"."}:
        return np.nan
    try:
        return float(token)
    except ValueError:
        return np.nan


@lru_cache(maxsize=1)
def mamajek_table() -> dict[str, np.ndarray]:
    """The dwarf sequence as a dict of arrays keyed by Mamajek's column names."""
    lines = MAMAJEK_PATH.read_text().splitlines()
    header_idx = next(i for i, l in enumerate(lines) if l.startswith("#SpT") and "Msun" in l)
    names = lines[header_idx].lstrip("#").split()
    rows = []
    for line in lines[header_idx + 1:]:
        if line.startswith("#"):
            break  # the table ends where its header is repeated
        parts = line.split()
        if len(parts) < len(names) - 1:
            continue
        rows.append(parts[: len(names)])
    out: dict[str, np.ndarray] = {"SpT": np.array([r[0] for r in rows])}
    for j, name in enumerate(names[1:], start=1):
        if name == "SpT":
            continue
        out[name] = np.array([_value(r[j]) if j < len(r) else np.nan for r in rows])
    return out


def _monotone(x: np.ndarray, y: np.ndarray):
    """Keep finite pairs, sorted by x, dropping repeated x values."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    order = np.argsort(x)
    x, y = x[order], y[order]
    keep = np.concatenate([[True], np.diff(x) > 0])
    return x[keep], y[keep]


@lru_cache(maxsize=1)
def _relations():
    t = mamajek_table()
    mass, mg, bprp, grp = t["Msun"], t["M_G"], t["Bp-Rp"], t["G-Rp"]
    # Brown dwarfs and the hottest stars lack Gaia colours; the relations below
    # cover B9 to M9.5, i.e. masses ~0.08 to ~2.8 Msun.
    rel = {
        "mass_to_mg": _monotone(mass, mg),
        "mass_to_bprp": _monotone(mass, bprp),
        "mass_to_grp": _monotone(mass, grp),
        "bprp_to_mg": _monotone(bprp, mg),
        "bprp_to_mass": _monotone(bprp, mass),
        "mg_to_mass": _monotone(mg, mass),
        "teff_to_bprp": _monotone(t["Teff"], bprp),
        "bprp_to_teff": _monotone(bprp, t["Teff"]),
    }
    return rel


def _interp(key: str, x, *, extrapolate: bool = False):
    xs, ys = _relations()[key]
    x = np.asarray(x, dtype=float)
    y = np.interp(x, xs, ys)
    if not extrapolate:
        y = np.where((x < xs[0]) | (x > xs[-1]), np.nan, y)
    return y


def mg_from_mass(mass):
    """Absolute G magnitude of a main-sequence star of this mass (Msun)."""
    return _interp("mass_to_mg", mass)


def bprp_from_mass(mass):
    return _interp("mass_to_bprp", mass)


def grp_from_mass(mass):
    """G - G_RP colour of a dwarf of this mass."""
    return _interp("mass_to_grp", mass)


def mass_from_bprp(bprp):
    """Mass of a dwarf with this (dereddened) BP-RP colour."""
    return _interp("bprp_to_mass", bprp)


def mass_from_mg(mg):
    """Mass of a dwarf with this absolute G magnitude."""
    return _interp("mg_to_mass", mg)


def ms_ridge_mg(bprp):
    """Absolute G of Mamajek's mean dwarf sequence at this colour."""
    return _interp("bprp_to_mg", bprp)


def teff_from_bprp(bprp):
    return _interp("bprp_to_teff", bprp)


def bprp_from_teff(teff):
    return _interp("teff_to_bprp", teff)


def spectral_type_edges() -> dict[str, float]:
    """BP-RP colour of each tabulated dwarf spectral type (for axis labels)."""
    t = mamajek_table()
    return {s: c for s, c in zip(t["SpT"], t["Bp-Rp"]) if np.isfinite(c)}


def flux_ratio_g(m1, m2):
    """G-band flux ratio F2/F1 of two dwarfs; companions below the table's
    mass range (brown dwarfs) are treated as dark (flux ratio 0)."""
    dm = mg_from_mass(m2) - mg_from_mass(m1)
    f = 10 ** (-0.4 * dm)
    return np.where(np.isfinite(f), f, 0.0)


def grvs_minus_grp(bprp, coeffs) -> np.ndarray:
    """G_RVS - G_RP as a polynomial in BP-RP (coefficients fitted to the 100 pc
    sample's own Gaia DR3 grvs_mag in ``samples.fit_grvs_colour``)."""
    return np.polyval(coeffs, np.clip(bprp, -0.2, 4.5))


def flux_ratio_rvs(m1, m2, coeffs):
    """RVS-band (845-872 nm) flux ratio of two dwarfs: what decides whether the
    companion's lines are visible and how much they dilute the primary's."""
    def m_rvs(m):
        mg = mg_from_mass(m)
        grp = grp_from_mass(m)
        bprp = bprp_from_mass(m)
        return mg - grp + grvs_minus_grp(bprp, coeffs)

    dm = m_rvs(m2) - m_rvs(m1)
    f = 10 ** (-0.4 * dm)
    return np.where(np.isfinite(f), f, 0.0)
