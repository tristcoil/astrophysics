"""Turn the archive table into analysis-ready samples.

The composite table is built for completeness: every confirmed planet gets a
row, and where the literature has no value the archive often supplies one.
Two of those fillers matter enormously for this project:

* **Radii calculated from masses.** For a planet with a measured mass but no
  measured radius -- almost every radial-velocity planet -- the archive fills
  in a radius from the Chen & Kipping (2017) empirical mass-radius relation.
  About 1,500 radii were made this way as of September 2026.
* **Masses calculated from radii.** The same relation, inverted, supplies the
  mass of nearly every transiting planet that has never been weighed: about
  3,000 planets, nearly half the table.

On a mass-radius diagram those rows would not be data. They would be the Chen &
Kipping curve drawn thousands of times over, and they would make the diagram
look far tidier than nature made it. They are recognisable by a reference
column that reads "Calculated Value", and for masses also by the provenance
column (``pl_bmassprov == "M-R relationship"``); both are excluded from every
measurement here. The archive documents what it calculates at
https://exoplanetarchive.ipac.caltech.edu/docs/pscp_calc.html.

Calculated *insolations* are kept. The archive derives those as ``L*/a**2``
from measured quantities -- a physical calculation, not an empirical fill.

Every sample records the rows remaining after each cut, so the README can say
exactly how it was built.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from config import (
    M_DWARF_TEFF_MAX, METHOD_GROUPS, MR_MAX_MASS_RELERR, MR_MAX_RADIUS_RELERR,
    OTHER_METHOD, VALLEY_MAX_PERIOD_DAYS, VALLEY_MAX_RADIUS_RELERR,
    VALLEY_MIN_LOGG, VALLEY_TEFF_RANGE_K,
)
from fetch_data import load_catalog_raw

#: How the archive marks a value it computed rather than took from a paper.
CALCULATED = "CALCULATED_VALUE"
#: Mass provenance for masses filled from the mass-radius relation.
MASS_FROM_RELATION = "M-R relationship"
#: Provenances that are genuine true-mass measurements. "Msin(i)/sin(i)" is a
#: minimum mass divided by an independently measured inclination.
TRUE_MASS = ("Mass", "Msin(i)/sin(i)")


def _plain(value) -> str:
    """Strip the HTML anchor the archive wraps around reference names."""
    if not isinstance(value, str):
        return ""
    return re.sub(r"<[^>]+>", "", value).strip()


def _is_calculated(reflink: pd.Series) -> pd.Series:
    return reflink.fillna("").str.contains(CALCULATED, regex=False)


def _no_limit(flag: pd.Series) -> pd.Series:
    """True where a value is a measurement rather than an upper/lower limit.

    Limit flags are 0 for a measurement, 1 for an upper limit and -1 for a
    lower limit. A mass *upper limit* plotted as a point is a classic way to
    misread a mass-radius diagram.
    """
    return flag.fillna(0).eq(0)


def _relative_error(value: pd.Series, err_up: pd.Series, err_down: pd.Series) -> pd.Series:
    """Symmetrised fractional uncertainty, ``(|err+| + |err-|) / 2 / value``.

    NaN when the value is missing or non-positive or either error is missing;
    such rows cannot pass a precision cut.
    """
    rel = (err_up.abs() + err_down.abs()) / 2.0 / value
    return rel.where(value.gt(0))


def load_catalog() -> pd.DataFrame:
    """The cached archive table plus the flags every sample is built from."""
    df = load_catalog_raw()

    df["method"] = df["discoverymethod"].map(METHOD_GROUPS).fillna(OTHER_METHOD)
    df["controversial"] = df["pl_controv_flag"].eq(1)
    df["transiting"] = df["tran_flag"].eq(1)

    df["radius_calculated"] = _is_calculated(df["pl_rade_reflink"])
    df["radius_measured"] = (
        df["pl_rade"].gt(0) & ~df["radius_calculated"] & _no_limit(df["pl_radelim"])
    )
    df["radius_relerr"] = _relative_error(df["pl_rade"], df["pl_radeerr1"], df["pl_radeerr2"])

    df["mass_from_relation"] = (
        df["pl_bmassprov"].eq(MASS_FROM_RELATION) | _is_calculated(df["pl_bmasse_reflink"])
    )
    # A true mass or a minimum mass (M sin i), but never one from the relation.
    df["mass_measured"] = (
        df["pl_bmasse"].gt(0) & ~df["mass_from_relation"] & _no_limit(df["pl_bmasselim"])
    )
    df["mass_is_true"] = df["mass_measured"] & df["pl_bmassprov"].isin(TRUE_MASS)
    df["mass_relerr"] = _relative_error(df["pl_bmasse"], df["pl_bmasseerr1"], df["pl_bmasseerr2"])

    df["period_ok"] = df["pl_orbper"].gt(0) & _no_limit(df["pl_orbperlim"])
    df["sma_ok"] = df["pl_orbsmax"].gt(0) & _no_limit(df["pl_orbsmaxlim"])
    df["insol_ok"] = df["pl_insol"].gt(0) & _no_limit(df["pl_insollim"])
    df["insol_calculated"] = _is_calculated(df["pl_insol_reflink"])
    df["teff_ok"] = df["st_teff"].gt(0) & _no_limit(df["st_tefflim"])
    df["m_dwarf_host"] = df["st_teff"].lt(M_DWARF_TEFF_MAX)

    df["disc_reference"] = df["disc_refname"].map(_plain)
    df["radius_reference"] = df["pl_rade_reflink"].map(_plain)
    return df


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

@dataclass
class Sample:
    """A named subset of the catalogue, with the cuts that produced it."""

    name: str
    description: str
    frame: pd.DataFrame
    #: (cut, planets remaining after it), in the order applied
    cuts: list[tuple[str, int]] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.frame)

    def cut_table(self) -> list[dict]:
        return [{"cut": c, "remaining": n} for c, n in self.cuts]


def _select(df: pd.DataFrame, name: str, description: str,
            steps: list[tuple[str, pd.Series]]) -> Sample:
    mask = pd.Series(True, index=df.index)
    cuts = [("all confirmed planets", len(df))]
    for label, condition in steps:
        mask &= condition.fillna(False).astype(bool)
        cuts.append((label, int(mask.sum())))
    return Sample(name, description, df[mask].copy(), cuts)


def census(df: pd.DataFrame) -> Sample:
    """Every confirmed planet. Controversial ones included: this is the census
    as the archive counts it."""
    return Sample("census", "every planet in the archive", df.copy(),
                  [("all confirmed planets", len(df))])


def radius_period(df: pd.DataFrame) -> Sample:
    return _select(df, "radius_period", "planets with a measured radius and period", [
        ("not flagged controversial", ~df["controversial"]),
        ("radius measured (not calculated, not a limit)", df["radius_measured"]),
        ("orbital period measured", df["period_ok"]),
    ])


def mass_separation(df: pd.DataFrame) -> Sample:
    return _select(df, "mass_separation", "planets with a measured mass or minimum mass", [
        ("not flagged controversial", ~df["controversial"]),
        ("mass or minimum mass measured (not from the M-R relation, not a limit)",
         df["mass_measured"]),
        ("semi-major axis available", df["sma_ok"]),
    ])


def _valley_steps(df: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    lo, hi = VALLEY_TEFF_RANGE_K
    return [
        ("not flagged controversial", ~df["controversial"]),
        ("transiting (radius from the transit depth)", df["transiting"]),
        ("radius measured (not calculated, not a limit)", df["radius_measured"]),
        ("radius uncertainty quoted", df["radius_relerr"].notna()),
        (f"period < {VALLEY_MAX_PERIOD_DAYS:.0f} d", df["period_ok"] & df["pl_orbper"].lt(VALLEY_MAX_PERIOD_DAYS)),
        (f"host Teff {lo:.0f}-{hi:.0f} K (FGK)", df["st_teff"].between(lo, hi)),
        (f"host log g > {VALLEY_MIN_LOGG} (main sequence)", df["st_logg"].gt(VALLEY_MIN_LOGG)),
    ]


def valley_all_precisions(df: pd.DataFrame) -> Sample:
    """The radius-valley sample before any precision cut."""
    return _select(df, "valley_all", "close-in planets of FGK dwarfs, any radius precision",
                   _valley_steps(df))


def valley(df: pd.DataFrame) -> Sample:
    """The headline radius-valley sample: radii known to better than 10%."""
    steps = _valley_steps(df) + [
        (f"radius uncertainty < {VALLEY_MAX_RADIUS_RELERR:.0%}",
         df["radius_relerr"].lt(VALLEY_MAX_RADIUS_RELERR)),
    ]
    return _select(df, "valley", "close-in planets of FGK dwarfs with precise radii", steps)


def mass_radius(df: pd.DataFrame) -> Sample:
    return _select(df, "mass_radius", "planets with a measured mass and radius", [
        ("not flagged controversial", ~df["controversial"]),
        ("transiting (radius from the transit depth)", df["transiting"]),
        ("radius measured (not calculated, not a limit)", df["radius_measured"]),
        ("true mass measured (not M sin i, not from the M-R relation, not a limit)",
         df["mass_is_true"]),
        (f"mass uncertainty < {MR_MAX_MASS_RELERR:.0%}", df["mass_relerr"].lt(MR_MAX_MASS_RELERR)),
        (f"radius uncertainty < {MR_MAX_RADIUS_RELERR:.0%}",
         df["radius_relerr"].lt(MR_MAX_RADIUS_RELERR)),
    ])


def habitable_zone(df: pd.DataFrame) -> Sample:
    return _select(df, "habitable_zone", "planets with a known insolation and host temperature", [
        ("not flagged controversial", ~df["controversial"]),
        ("insolation available (measured or L/a^2)", df["insol_ok"]),
        ("host effective temperature available", df["teff_ok"]),
    ])


def all_samples(df: pd.DataFrame) -> dict[str, Sample]:
    builders = (census, radius_period, mass_separation, valley_all_precisions,
                valley, mass_radius, habitable_zone)
    return {s.name: s for s in (build(df) for build in builders)}


def summarise_fillers(df: pd.DataFrame) -> dict:
    """How much of the table is filled in rather than measured. For the README."""
    return {
        "planets": int(len(df)),
        "radius_present": int(df["pl_rade"].notna().sum()),
        "radius_calculated": int(df["radius_calculated"].sum()),
        "radius_measured": int(df["radius_measured"].sum()),
        "mass_present": int(df["pl_bmasse"].notna().sum()),
        "mass_from_relation": int(df["mass_from_relation"].sum()),
        "mass_measured": int(df["mass_measured"].sum()),
        "mass_true": int(df["mass_is_true"].sum()),
        "mass_minimum_only": int((df["mass_measured"] & ~df["mass_is_true"]).sum()),
        "mass_limits": int(df["pl_bmasselim"].fillna(0).ne(0).sum()),
        "both_measured": int((df["radius_measured"] & df["mass_is_true"]).sum()),
        "insolation_calculated": int(df["insol_calculated"].sum()),
        "controversial": int(df["controversial"].sum()),
        "fraction_mass_from_relation": float(df["mass_from_relation"].mean()),
    }


def symmetric_error(value: pd.Series, relerr: pd.Series) -> np.ndarray:
    """Absolute 1-sigma error from a value and its symmetrised relative error."""
    return np.asarray(value, dtype=float) * np.asarray(relerr, dtype=float)
