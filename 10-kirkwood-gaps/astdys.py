"""Read the AstDyS proper-element and family files.

**Why proper elements.** The orbital elements the MPC publishes are
*osculating*: the ellipse an asteroid would follow if the planets vanished at
the epoch. Jupiter's pull makes them oscillate -- eccentricity and inclination
swing by several hundredths over thousands of years -- so two fragments of the
same collision, launched together, quickly look unrelated. *Proper* elements
average those oscillations away by integrating each orbit numerically for
millions of years (Knezevic & Milani 2000). What remains is nearly constant
over hundreds of millions of years, and in it the families snap into focus.

Files (all from https://newton.spacedys.com/~astdys2/propsynth/):

* ``all.syn`` -- proper a, e, sin i for 1.15 million main-belt and Hungaria
  asteroids, plus two proper *frequencies*: g, the rate at which the
  perihelion turns, and s, the rate at which the orbital plane turns.
* ``secres.syn`` -- 26,000 more asteroids locked in secular resonances, whose
  proper e is undefined (the column holds the libration amplitude instead).
* ``all_tro.members`` -- family membership for the 265,000 asteroids that
  belong to one (Milani et al. 2014 classification, as updated by AstDyS).
* ``all_tro.famtab`` -- one row per family.
* ``family_list.txt`` -- family names, as the AstDyS plot tool's HTML list.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from config import ASTDYS_PARSED, DATA_DIR

SYN_COLUMNS = ["name", "H", "a", "e", "sini", "n_deg_yr", "g", "s", "lce", "my"]

#: Membership status codes in ``all_tro.members``. The meaning follows from
#: the family table: summing its per-family counts over all families
#: reproduces the status counts exactly (core 26,200; small 16,596; added
#: 221,692 + 71 attributed to two families).
STATUS = {3: "core", 1: "small family", 2: "added", 4: "double"}


def _read_syn(path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", comment="%", header=None, names=SYN_COLUMNS,
                     dtype={"name": str})
    return df


def _number(names: pd.Series) -> np.ndarray:
    """Asteroid number, or 0 for a provisional designation like '1979MW5'."""
    digits = names.str.fullmatch(r"\d+")
    out = np.zeros(len(names), dtype=np.int64)
    out[digits.to_numpy()] = names[digits].astype(np.int64).to_numpy()
    return out


def family_names() -> dict[int, str]:
    """Family id (the number of its namesake asteroid) -> name."""
    text = (DATA_DIR / "family_list.txt").read_text()
    return {int(k): v.strip() for k, v in
            re.findall(r'<option value="(\d+)">\s*\(\d+\)\s*([^<]+?)\s*</option>', text)}


def family_table() -> pd.DataFrame:
    cols = ["family", "zone", "d_max", "n_core", "n_small", "n_added", "n_total",
            "a_min", "a_max", "e_min", "e_max", "sini_min", "sini_max"]
    df = pd.read_csv(DATA_DIR / "all_tro.famtab", sep=r"\s+", comment="%",
                     header=None, names=cols, usecols=range(len(cols)))
    names = family_names()
    df["name"] = df["family"].map(names)
    return df.sort_values("n_total", ascending=False).reset_index(drop=True)


def _members() -> pd.DataFrame:
    cols = ["name", "H", "status", "family", "dv", "near", "family2", "dv2", "near2",
            "rescod"]
    return pd.read_csv(DATA_DIR / "all_tro.members", sep=r"\s+", comment="%",
                       header=None, names=cols, dtype={"name": str, "rescod": str})


def parse_proper() -> dict[str, np.ndarray]:
    main = _read_syn(DATA_DIR / "all.syn")
    main["secres"] = False
    sec = _read_syn(DATA_DIR / "secres.syn")
    sec["secres"] = True
    # In secres.syn the "e" column is a libration amplitude, not a proper e.
    sec["e"] = np.nan
    df = pd.concat([main, sec], ignore_index=True)
    dup = df["name"].duplicated(keep="first")
    df = df[~dup].reset_index(drop=True)

    mem = _members()[["name", "family", "status"]]
    df = df.merge(mem, on="name", how="left")
    df["family"] = df["family"].fillna(0).astype(np.int64)
    df["status"] = df["status"].fillna(0).astype(np.int8)

    out = {c: df[c].to_numpy() for c in
           ("H", "a", "e", "sini", "n_deg_yr", "g", "s", "lce", "my", "family",
            "status", "secres")}
    out["number"] = _number(df["name"])
    out["name"] = df["name"].to_numpy().astype("S12")
    out["n_duplicates_dropped"] = np.array(int(dup.sum()))
    return out


def load_proper(*, refresh: bool = False, quiet: bool = False) -> dict[str, np.ndarray]:
    """Proper elements (+ family id and status), cached as ``.npz``."""
    src = DATA_DIR / "all.syn"
    if (not refresh and ASTDYS_PARSED.exists()
            and ASTDYS_PARSED.stat().st_mtime > src.stat().st_mtime):
        with np.load(ASTDYS_PARSED) as z:
            return {k: z[k] for k in z.files}
    if not quiet:
        print("  parsing AstDyS proper elements (once; cached afterwards)")
    data = parse_proper()
    np.savez(ASTDYS_PARSED, **data)
    return data
