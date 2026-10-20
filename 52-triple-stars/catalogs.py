"""Readers for the catalogues: Tokovinin's Multiple Star Catalog (MSC, 2026
release) and the two VizieR tables used for the statistics.

The MSC describes each multiple system as a tree of "elementary binaries"
(subsystems). Every subsystem row names its primary and secondary component
and its *parent*: the component label, in a higher-level subsystem, that this
pair makes up. For TIC 290061484, for example:

    AB , C  , parent *     (outermost: the triple AB and the fourth star C)
    A  , B  , parent AB    (the 24.5-d orbit: pair A and star B)
    Aa , Ab , parent A     (the 1.8-d eclipsing pair)

``inner_outer_pairs`` walks these trees and returns every (inner, outer)
couple: the building block of a hierarchical triple.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

MSC_DIR = C.RAW / "msc"

#: period units of the MSC: days, years, kiloyears, megayears
_P_UNIT_DAYS = {"d": 1.0, "y": 365.25, "k": 365.25e3, "M": 365.25e6}
#: separation units: arcsec, mas, arcmin
_SEP_UNIT_ARCSEC = {'"': 1.0, "m": 1e-3, "'": 60.0}


def _read_pipe(path, names):
    df = pd.read_csv(path, sep="|", header=None, names=names, dtype=str,
                     keep_default_na=False, quoting=3, engine="python")
    for c in df.columns:
        df[c] = df[c].str.strip()
    return df


def _num(s):
    return pd.to_numeric(s.str.replace("*", "", regex=False), errors="coerce")


def read_msc():
    """The four MSC tables as DataFrames, numbers converted, zeros = unknown."""
    comp = _read_pipe(MSC_DIR / "comp.tsv",
                      ["wds", "ra", "dec", "plx", "plx_ref", "pmra", "pmdec", "rv", "comp", "sep",
                       "sp", "hip", "hd", "B", "V", "G", "J", "H", "K", "ncomp", "grade", "ids"])
    for c in ("ra", "dec", "plx", "pmra", "pmdec", "rv", "sep", "B", "V", "G", "J", "H", "K"):
        comp[c] = _num(comp[c])
    for c in ("hip", "hd", "ncomp", "grade"):
        comp[c] = _num(comp[c]).fillna(0).astype(int)

    sys_ = _read_pipe(MSC_DIR / "sys.tsv",
                      ["wds", "prim", "sec", "parent", "type", "P", "P_unit", "sep", "sep_unit", "pa",
                       "V1", "sp1", "V2", "sp2", "m1", "m1_code", "m2", "m2_code", "comment"])
    for c in ("P", "sep", "pa", "V1", "V2", "m1", "m2"):
        sys_[c] = _num(sys_[c])
    sys_["P_days"] = sys_["P"] * sys_["P_unit"].map(_P_UNIT_DAYS)
    sys_["sep_arcsec"] = sys_["sep"] * sys_["sep_unit"].map(_SEP_UNIT_ARCSEC)
    sys_.loc[sys_["P_days"] <= 0, "P_days"] = np.nan
    for c in ("m1", "m2"):
        sys_.loc[sys_[c] <= 0, c] = np.nan
    sys_["system"] = sys_["prim"] + "," + sys_["sec"]

    orb = _read_pipe(MSC_DIR / "orb.tsv",
                     ["wds", "system", "P", "T", "e", "a", "Omega", "omega", "i", "K1", "K2", "V0",
                      "node_flag", "P_unit", "note"])
    for c in ("P", "T", "e", "a", "Omega", "omega", "i", "K1", "K2", "V0"):
        orb[c] = _num(orb[c])
    orb["P_days"] = orb["P"] * orb["P_unit"].map({"d": 1.0, "y": 365.25}).fillna(np.nan)
    # a few systems are written "Aa.Ab" instead of "Aa,Ab"
    orb["system"] = orb["system"].str.replace(".", ",", regex=False)
    # visual orbit: semimajor axis and inclination both given
    orb["visual"] = (orb["a"] > 0) & (orb["i"] > 0)
    orb["spectroscopic"] = (orb["K1"] > 0) | (orb["K2"] > 0)

    notes = _read_pipe(MSC_DIR / "notes.tsv", ["wds", "remark", "bib"])
    return comp, sys_, orb, notes


def primary_components(comp):
    """One row per system: the main component (zero separation) with its grade,
    parallax and identifiers."""
    main = comp[comp["sep"].fillna(0) == 0].drop_duplicates("wds")
    return main.set_index("wds")


def inner_outer_pairs(sys_, orb, comp):
    """Every (inner subsystem, outer subsystem) couple of the MSC.

    Columns: wds, inner/outer system labels and types, periods (days), the
    inner masses m1, m2 and the tertiary mass m3 (the *other* component of
    the outer pair, itself possibly a binary), outer eccentricity when an
    orbit is known, and whether the periods come from orbits or from
    separations ("estimated").
    """
    main = primary_components(comp)
    rows = []
    by_wds = {w: g for w, g in sys_.groupby("wds", sort=False)}
    orb_e = (orb[orb["e"].notna()].sort_values("visual", ascending=False)
             .drop_duplicates(["wds", "system"]).set_index(["wds", "system"]))
    for wds, g in by_wds.items():
        g = g.reset_index(drop=True)
        for _, inn in g.iterrows():
            par = inn["parent"]
            if par in ("*", "t", ""):
                continue
            outs = g[(g["prim"] == par) | (g["sec"] == par)]
            if len(outs) == 0:
                continue
            out = outs.iloc[0]
            inner_is_primary = out["prim"] == par
            m3 = out["m2"] if inner_is_primary else out["m1"]
            m_in = (inn["m1"] if np.isfinite(inn["m1"]) else np.nan) + \
                   (inn["m2"] if np.isfinite(inn["m2"]) else np.nan)
            e_out = np.nan
            key = (wds, out["system"])
            if key in orb_e.index:
                e_out = float(orb_e.loc[key, "e"])
            e_in = np.nan
            key_in = (wds, inn["system"])
            if key_in in orb_e.index:
                e_in = float(orb_e.loc[key_in, "e"])
            has_orbit = lambda t: any(x in t for x in ("V", "S1", "S2", "A", "E"))
            rows.append(dict(
                wds=wds, inner=inn["system"], outer=out["system"], inner_type=inn["type"],
                outer_type=out["type"], P_in=inn["P_days"], P_out=out["P_days"],
                m1=inn["m1"], m2=inn["m2"], m_in=m_in, m3=m3, e_in=e_in, e_out=e_out,
                inner_from_orbit=has_orbit(inn["type"]), outer_from_orbit=has_orbit(out["type"]),
                inner_eclipsing="E" in inn["type"], outer_type_wide=("C" in out["type"]),
                grade=int(main.loc[wds, "grade"]) if wds in main.index else 0,
                plx=float(main.loc[wds, "plx"]) if wds in main.index else np.nan,
                ids=main.loc[wds, "ids"] if wds in main.index else "",
            ))
    df = pd.DataFrame(rows)
    df["ratio"] = df["P_out"] / df["P_in"]
    return df


def visual_orbit_pairs(pairs, orb):
    """Inner/outer couples where both subsystems have a visual orbit with a
    node and an inclination: the two candidate mutual inclinations Phi1, Phi2
    (the 180-degree node ambiguity of visual orbits, Tokovinin 2017 eq. 2)."""
    vis = orb[orb["visual"] & (orb["Omega"].notna())].drop_duplicates(["wds", "system"])
    vis = vis.set_index(["wds", "system"])
    out = []
    for _, r in pairs.iterrows():
        ki, ko = (r["wds"], r["inner"]), (r["wds"], r["outer"])
        if ki in vis.index and ko in vis.index:
            a, b = vis.loc[ki], vis.loc[ko]
            phi1, phi2 = mutual_inclinations(a["i"], a["Omega"], b["i"], b["Omega"])
            out.append(dict(wds=r["wds"], inner=r["inner"], outer=r["outer"], P_in_yr=a["P_days"] / 365.25,
                            P_out_yr=b["P_days"] / 365.25, a_in=a["a"], a_out=b["a"], e_in=a["e"],
                            e_out=b["e"], i_in=a["i"], Omega_in=a["Omega"], i_out=b["i"],
                            Omega_out=b["Omega"], node_in=a["node_flag"], node_out=b["node_flag"],
                            phi1=phi1, phi2=phi2, plx=r["plx"], m_in=r["m_in"], m3=r["m3"]))
    return pd.DataFrame(out)


def mutual_inclinations(i1, O1, i2, O2):
    """Both solutions of cos Phi = cos i1 cos i2 +- sin i1 sin i2 cos(O1 - O2)."""
    i1, O1, i2, O2 = np.radians([i1, O1, i2, O2])
    a = np.cos(i1) * np.cos(i2)
    b = np.sin(i1) * np.sin(i2) * np.cos(O1 - O2)
    return (np.degrees(np.arccos(np.clip(a + b, -1, 1))),
            np.degrees(np.arccos(np.clip(a - b, -1, 1))))


# ---------------------------------------------------------------------------
# VizieR tables
# ---------------------------------------------------------------------------

def _fixed(path, colspecs, names):
    return pd.read_fwf(path, colspecs=colspecs, names=names, header=None, dtype=str,
                       keep_default_na=False)


def read_tokovinin2017():
    """Tokovinin (2017, ApJ 844, 103): Table 1 (443 triples with a resolved
    inner orbit, outer revolution sense) and Table 3 (two visual orbits)."""
    folder = C.RAW / "vizier" / "J_ApJ_844_103"
    t1 = _fixed(folder / "table1.dat",
                [(0, 10), (11, 16), (17, 22), (23, 29), (30, 35), (36, 42), (43, 48), (49, 54),
                 (55, 62), (63, 68), (69, 71), (71, 77)],
                ["wds", "plx", "M1", "P1", "e1", "a1", "Omega1", "i1", "rho", "theta", "revdir", "sign"])
    for c in t1.columns[1:]:
        t1[c] = pd.to_numeric(t1[c], errors="coerce")
    t3 = _fixed(folder / "table3.dat",
                [(0, 10), (11, 12), (13, 20), (21, 26), (27, 32), (33, 38), (39, 44), (45, 46),
                 (47, 55), (56, 61)],
                ["wds", "io", "P", "e", "a", "Omega", "i", "grade", "ref", "phi"])
    for c in ("P", "e", "a", "Omega", "i", "grade", "phi"):
        t3[c] = pd.to_numeric(t3[c], errors="coerce")
    return t1, t3


def read_tokovinin2006():
    """Tokovinin, Thomas, Sterzik & Udry (2006, A&A 450, 681): Table 5, the
    subsystems of the 161 multiple systems containing their 165 SBs."""
    folder = C.RAW / "vizier" / "J_AA_450_681"
    t5 = _fixed(folder / "table5.dat",
                [(0, 6), (8, 15), (16, 22), (23, 30), (31, 32), (33, 41), (42, 43), (44, 48), (49, 50),
                 (51, 55), (56, 57), (58, 100)],
                ["hip", "comp", "type", "sep", "sep_unit", "P", "P_unit", "m1", "m1_code", "m2",
                 "m2_code", "remark"])
    for c in ("hip", "sep", "P", "m1", "m2"):
        t5[c] = pd.to_numeric(t5[c], errors="coerce")
    t5["P_days"] = t5["P"] * t5["P_unit"].map({"d": 1.0, "y": 365.25, "k": 365.25e3})
    parts = t5["comp"].str.split(",", expand=True)
    t5["prim"], t5["sec"], t5["parent"] = parts[0], parts[1], parts[2]
    t2 = _fixed(folder / "table2.dat",
                [(0, 6), (7, 14), (16, 18), (19, 21), (22, 24), (25, 26), (26, 28), (29, 31), (32, 34),
                 (35, 40), (41, 46), (47, 52), (53, 58), (59, 64), (65, 71), (73, 74), (75, 76), (77, 78),
                 (79, 80)],
                ["hip", "alt", "rah", "ram", "ras", "de_sign", "ded", "dem", "des", "pmra", "pmde", "plx",
                 "V", "BV", "sp", "ncomp", "ao", "twomass", "speckle"])
    for c in ("hip", "plx", "V", "BV", "ncomp", "ao", "twomass", "speckle"):
        t2[c] = pd.to_numeric(t2[c], errors="coerce")
    return t2, t5
