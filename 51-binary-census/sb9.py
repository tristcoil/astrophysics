"""SB9, the 9th Catalogue of Spectroscopic Binary Orbits, made analysis-ready.

SB9 (Pourbaix et al. 2004) collects published orbits of spectroscopic
binaries: a century of radial-velocity work by many observers, so it is
heterogeneous and favours bright, short-period, large-amplitude systems. Some
systems have several published orbits; we keep one per system, preferring the
best grade, then the most recent paper.

For each system we also find its Gaia DR3 counterpart (CDS xMatch within 10
arcsec, then the candidate whose G magnitude best matches SB9's), which gives a
distance and an absolute magnitude, so dwarfs can be told from giants without
trusting the spectral type alone.
"""

from __future__ import annotations

import re

import numpy as np
from astropy.table import Table

import config as C
import stellar

LUMINOSITY_ORDER = ("Ia", "Iab", "Ib", "II", "III", "IV", "V")


def _s(x) -> str:
    return x.decode() if isinstance(x, bytes) else str(x)


def _f(col) -> np.ndarray:
    """A float array with masked values as NaN."""
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def parse_spectral_type(sp: str) -> tuple[str, float, str]:
    """('G', 2.0, 'V') from 'G2V'; luminosity class '' when absent.

    Composite or range types ('F5-G3Ib', 'G8III+A') are reduced to their first
    component; peculiarity suffixes (m, p, e, n, s...) are ignored.
    """
    sp = sp.strip()
    m = re.match(r"^(sd|d|g)?([OBAFGKM])\s*(\d+(?:\.\d+)?)?", sp)
    if not m:
        return "", np.nan, ""
    prefix, letter = m.group(1), m.group(2)
    sub = float(m.group(3)) if m.group(3) else np.nan
    rest = sp[m.end():].split("+")[0]
    lum = ""
    # the first Roman-numeral class after the type, e.g. 'III', 'IV-V', 'Ib-II'
    numeral = r"(?:VI|IV|V|III|II|I)(?:ab|a|b)?"
    r = re.search(rf"{numeral}(?:[-/]{numeral})?", rest)
    if r:
        lum = r.group(0)
    if prefix == "sd":
        lum = "VI"
    elif prefix == "d" and not lum:
        lum = "V"
    elif prefix == "g" and not lum:
        lum = "III"
    return letter, sub, lum


def luminosity_group(lum: str) -> str:
    """'dwarf' (V, IV-V), 'subgiant' (IV, IV/III mix), 'giant' (III and brighter) or ''."""
    if not lum:
        return ""
    first = re.sub(r"(ab|a|b)$", "", re.split(r"[-/]", lum)[0])
    if first in ("V", "VI") or lum in ("IV-V", "IV/V"):
        return "dwarf"
    if first == "IV":
        return "subgiant"
    return "giant"


def _ref_year(ref: str) -> int:
    m = re.match(r"(\d{4})", ref.strip())
    return int(m.group(1)) if m else 0


def load_sb9() -> Table:
    """One orbit per SB9 system, with parsed spectral types and Gaia data."""
    orbits = Table.read(C.VIZIER_DIR / "sb9_orbits.fits")
    main = Table.read(C.VIZIER_DIR / "sb9_main.fits")

    grade = _f(orbits["Grade"])
    grade = np.where(np.isfinite(grade), grade, -1)
    year = np.array([_ref_year(_s(r)) for r in orbits["Ref"]])
    n_rv = _f(orbits["o_K1"])
    n_rv = np.where(np.isfinite(n_rv), n_rv, 0)
    seq = np.asarray(orbits["Seq"])
    # best grade first, then most recent paper, then most velocities
    order = np.lexsort((-n_rv, -year, -grade, seq))
    first = np.concatenate([[True], np.diff(seq[order]) != 0])
    best = orbits[order[first]]

    t = Table()
    t["seq"] = np.asarray(best["Seq"])
    t["period"] = _f(best["Per"])
    t["ecc"] = _f(best["e"])
    t["ecc_err"] = _f(best["e_e"])
    t["ecc_fixed"] = np.array([_s(v).strip() == "a" for v in best["f_e"]])
    t["k1"] = _f(best["K1"])
    t["k2"] = _f(best["K2"])
    t["grade"] = np.where(np.isfinite(_f(best["Grade"])), _f(best["Grade"]), -1)
    t["ref"] = np.array([_s(r).strip() for r in best["Ref"]])
    t["n_orbits"] = np.bincount(seq)[t["seq"]]

    idx = {s: i for i, s in enumerate(np.asarray(main["Seq"]))}
    rows = np.array([idx[s] for s in t["seq"]])
    t["name"] = np.array([_s(n).strip() for n in main["Name"][rows]])
    t["sp1"] = np.array([_s(s).strip() for s in main["Sp1"][rows]])
    t["sp2"] = np.array([_s(s).strip() for s in main["Sp2"][rows]])
    t["mag1"] = _f(main["mag1"][rows])
    t["mag1_band"] = np.array([_s(b).strip() for b in main["n_mag1"][rows]])
    parsed = [parse_spectral_type(s) for s in t["sp1"]]
    t["sp_class"] = np.array([p[0] for p in parsed])
    t["sp_sub"] = np.array([p[1] for p in parsed])
    t["lum_class"] = np.array([p[2] for p in parsed])
    t["lum_group"] = np.array([luminosity_group(p[2]) for p in parsed])
    t["double_lined"] = np.isfinite(t["k2"])
    _attach_gaia(t)
    return t


def _attach_gaia(t: Table) -> None:
    """Best Gaia DR3 counterpart per SB9 system.

    Among xMatch candidates within 10", keep those within 1.5 mag of SB9's
    magnitude (V for most systems; G - V is -0.1 to -0.9 for FGKM dwarfs), and
    take the brightest -- a binary's companion can be a separate, fainter
    Gaia source nearby.
    """
    x = Table.read(C.VIZIER_DIR / "sb9_gaia_xmatch.fits")
    seq = np.asarray(x["seq"])
    gmag = _f(x["Gmag"])
    order = np.argsort(gmag)
    best: dict[int, int] = {}
    mag_by_seq = dict(zip(t["seq"], t["mag1"]))
    for j in order:
        s = int(seq[j])
        if s in best:
            continue
        m = mag_by_seq.get(s, np.nan)
        if np.isfinite(m) and np.isfinite(gmag[j]) and abs(gmag[j] - m) > 1.5:
            continue
        best[s] = j
    cols = {"gaia_source_id": "Source", "gaia_g": "Gmag", "gaia_bprp": "BP-RP", "gaia_plx": "Plx",
            "gaia_plx_err": "e_Plx", "gaia_ruwe": "RUWE", "gaia_sep_arcsec": "angDist"}
    for new, old in cols.items():
        vals = _f(x[old]) if old != "Source" else np.asarray(x[old]).astype(np.int64)
        fill = -1 if old == "Source" else np.nan
        t[new] = np.array([vals[best[s]] if s in best else fill for s in t["seq"]])
    plx = t["gaia_plx"]
    good = (plx > 0) & (plx / t["gaia_plx_err"] > 5)
    t["abs_g"] = np.where(good, t["gaia_g"] + 5 * np.log10(np.where(good, plx, 1) / 100.0), np.nan)
    # height above Mamajek's dwarf sequence at the same colour (negative = brighter)
    t["dg_ms"] = t["abs_g"] - stellar.ms_ridge_mg(t["gaia_bprp"])


def main_sequence_mask(t: Table, *, max_above: float = 1.2, max_below: float = 0.6) -> np.ndarray:
    """SB9 systems on the main sequence by Gaia photometry *and* spectral type.

    A binary can sit up to 0.75 mag above the single-star sequence (an equal
    pair), more for triples; anything brighter than ``max_above`` mag above the
    sequence is evolved. Spectral luminosity class, where given, must not be
    giant.
    """
    phot = (t["dg_ms"] > -max_above) & (t["dg_ms"] < max_below)
    spec_ok = t["lum_group"] != "giant"
    return phot & spec_ok
