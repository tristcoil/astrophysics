"""A century of gamma: every independent light-deflection and time-delay result.

The PPN parameter gamma measures how much space curvature a unit of mass
produces; general relativity has gamma = 1. Light passing the Sun is bent by
(1 + gamma)/2 x 1.7512" at the limb and delayed by (1 + gamma) GM/c^3 x
ln(...), so every such experiment measures (1 + gamma)/2.

``transcribed/gamma_history.csv`` is an unchanged copy of the table that
project N7 (64-gravity-scorecard) checked against the primary papers in
October 2026: 58 results, each exactly as published (a deflection in
arcseconds, a ratio to Einstein's value, gamma, or gamma - 1) with its error
type (probable errors before about 1950; one-sigma, "realistic" or formal
errors later). This module converts them to one sigma and to gamma, picks one
point per independent experiment (N7's recommended set), and tracks the best
precision reached by each year.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

#: Divide a published error by this to get one sigma (N7's convention: a
#: probable error is 0.6745 sigma; INPOP-style "3 sigma" are 3 x formal, which
#: their authors treat as a realistic 1 sigma).
TO_ONE_SIGMA = {"1sigma": 1.0, "68%": 1.0, "standard_error": 1.0, "formal": 1.0, "realistic": 1.0,
                "not_stated": 1.0, "3sigma": 1.0, "2sigma": 2.0, "90%": 1.645, "95%": 1.96,
                "probable_error": 0.6745}

#: One point per independent experiment, as recommended in N7's notes (section 6):
#: (row id, flag). Flags: "rereduction" (a later reduction of earlier plates),
#: "secondary" (number read only in a secondary source), "flagged" (the authors
#: suspect systematics). Planetary-ephemeris and galaxy-lensing rows are left out:
#: their "intervals" are not standard errors, and lensing tests gamma on other scales.
ENTRIES = [
    ("E1919-sobral-4inch", ""), ("E1919-principe", ""), ("E1919-sobral-combined-harvey1979", "rereduction"),
    ("E1922-wallal-lick-15ft", ""), ("E1922-wallal-lick-5ft", ""), ("E1922-cordillo-downs", ""),
    ("E1929-takengon", ""), ("E1936-ussr-mikhailov-1959", ""), ("E1947-bocaiuva", ""), ("E1952-khartoum", ""),
    ("E1973-chinguetti", ""), ("E2017-wyoming-bruns", ""),
    ("R1969-owens-valley", ""), ("R1969-goldstone", ""), ("R1970-goldstone-haystack-vlbi", "secondary"),
    ("R1970-nrao-sramek", ""), ("R1970-cambridge-hill", ""), ("R1971-nrao-sramek", "secondary"),
    ("R1972-cambridge-riley", ""), ("R1972-westerbork", ""), ("R1972-vlbi-counselman", ""),
    ("R1973-westerbork-dual", ""), ("R1974-nrao-fomalont", ""), ("R1975-nrao-fomalont", ""),
    ("V1984-robertson-carter", ""), ("V1987-lebach", ""), ("V1991-robertson-carter-dillinger", ""),
    ("V2004-shapiro", ""), ("V2005-vlba-fomalont", ""), ("V2011-lambert", ""), ("V2018-titov-global", "flagged"),
    ("A1997-hipparcos", ""),
    ("T1967-shapiro-radar", "secondary"), ("T1970-shapiro-radar", ""), ("T1970-mariner6-7", ""),
    ("T1972-mariner9", "secondary"), ("T1976-viking", ""), ("T1985-psr1937", "secondary"), ("T1985-voyager2", ""),
    ("T2002-cassini", ""),
]

#: Broad families of method, for the figure: starlight, radio waves, timing.
FAMILY = {"optical_eclipse": "starlight", "space_astrometry": "starlight",
          "radio_interferometry": "radio", "vlbi": "radio",
          "time_delay_radar": "delay", "time_delay_spacecraft": "delay", "time_delay_pulsar": "delay"}


def load() -> pd.DataFrame:
    return pd.read_csv(C.TRANSCRIBED_DIR / "gamma_history.csv")


def to_half_one_plus_gamma(quantity: str, value: float, err: float) -> tuple[float, float]:
    """(1 + gamma)/2 and its error from a published value of any kind."""
    if quantity == "deflection_limb_arcsec":
        return value / C.DEFLECTION_LIMB_GR, err / C.DEFLECTION_LIMB_GR
    if quantity in ("deflection_ratio_to_GR", "(1+gamma)/2"):
        return value, err
    if quantity == "gamma":
        return (1 + value) / 2, err / 2
    if quantity == "gamma-1":
        return 1 + value / 2, err / 2
    raise ValueError(f"unknown quantity {quantity!r}")


def timeline() -> pd.DataFrame:
    """One row per independent experiment: (1 + gamma)/2 and gamma with one-sigma errors."""
    g = load().set_index("id")
    rows = []
    for rid, flag in ENTRIES:
        r = g.loc[rid]
        errs = pd.to_numeric(pd.Series([r["err_plus"], r["err_minus"]]), errors="coerce")
        err = float(np.nanmean(errs)) / TO_ONE_SIGMA[str(r["err_type"]).strip()]
        half, sig = to_half_one_plus_gamma(str(r["quantity_as_published"]), float(r["value"]), err)
        rows.append({"id": rid, "year": int(r["year"]), "method": r["method"], "family": FAMILY[r["method"]],
                     "experiment": r["experiment"], "flag": flag, "err_type": r["err_type"],
                     "half_1pg": half, "sigma_half": sig, "gamma": 2 * half - 1, "sigma_gamma": 2 * sig,
                     "converted_from_pe": r["err_type"] == "probable_error"})
    return pd.DataFrame(rows).sort_values(["year", "sigma_gamma"]).reset_index(drop=True)


def frontier(t: pd.DataFrame | None = None) -> pd.DataFrame:
    """The best (smallest) one-sigma error on gamma reached by each year, counting
    only first-hand, unflagged results: the 'precision ladder'."""
    t = timeline() if t is None else t
    ok = t[t["flag"] == ""].sort_values("year")
    best, rows = np.inf, []
    for year, grp in ok.groupby("year"):
        i = grp["sigma_gamma"].idxmin()
        if grp.loc[i, "sigma_gamma"] < best:
            best = float(grp.loc[i, "sigma_gamma"])
            rows.append({"year": int(year), "sigma_gamma": best, "id": grp.loc[i, "id"],
                         "family": grp.loc[i, "family"]})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    t = timeline()
    cols = ["year", "id", "family", "half_1pg", "sigma_half", "sigma_gamma", "flag"]
    with pd.option_context("display.width", 160, "display.max_rows", 200):
        print(t[cols].to_string(index=False, float_format=lambda v: f"{v:.6g}"))
        print(frontier(t).to_string(index=False))
