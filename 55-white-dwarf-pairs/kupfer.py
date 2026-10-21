"""The known LISA binaries: Tables 1 and 2 of Kupfer et al. (2024, ApJ 963, 100).

Parsed from the LaTeX source of arXiv:2302.12719 (the tables are not on VizieR, and
the LISA Consortium GitLab mirror cited in the paper did not answer in October 2026).
Every value keeps its kind: measured (with errors), assumed ("[...]" in the paper:
set from evolutionary arguments), a range ("0.49-0.80"), a lower limit (">0.17") or
approximate ("~38").

Two features of the source tables:

* TIC 378898110 appears in both the "verification" and the "detectable" block of
  each table. The paper's own counts (18 verification + 22 detectable = 40) only
  work if it is counted once, as detectable, which is what we do.
* SDSS J0651's masses are listed in the opposite order to the other rows
  (m1 = 0.247 < m2 = 0.49); chirp masses do not care, the fates figure sorts them.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

import config as C

TEX = C.LIT_DIR / "arxiv_2302.12719" / "sample631.tex"
OUT = C.DERIVED_DIR / "kupfer2024_binaries.csv"


def _clean(s: str) -> str:
    s = s.replace("\\,", "").replace("$", "").replace("{", "").replace("}", "")
    s = s.replace("--", "-").replace("\\approx", "~").replace("\\Msun", "").strip()
    return s


def parse_number_with_paren(s: str):
    """'321.529129(10)' -> (321.529129, 1.0e-5); '939.0(7.2)' -> (939.0, 7.2)."""
    s = _clean(s)
    m = re.match(r"^([0-9.]+)\(([0-9.]+)\)$", s)
    if m:
        val, err = m.group(1), m.group(2)
        if "." in err:
            return float(val), float(err)
        ndec = len(val.split(".")[1]) if "." in val else 0
        return float(val), float(err) * 10.0**-ndec
    return float(s), np.nan


def parse_value(raw: str):
    """Return (value, err_plus, err_minus, kind) for a mass/inclination/distance cell."""
    s = raw.strip()
    kind = "measured"
    if s.startswith("[") and s.endswith("]"):
        kind = "assumed"
        s = s[1:-1]
    s = _clean(s)
    s = re.sub(r"\^[a-z]$", "", s)  # footnote marks like ^c
    if s in ("", "-"):
        return np.nan, np.nan, np.nan, "missing"
    if s.startswith(">"):
        s = s[1:]
        kind = "lower_limit" if kind == "measured" else kind
        m = re.match(r"^([0-9.]+)(?:\\pm([0-9.]+))?$", s)
        val = float(m.group(1))
        err = float(m.group(2)) if m.group(2) else np.nan
        return val, err, err, kind
    if s.startswith("~"):
        return float(s[1:]), np.nan, np.nan, "approx"
    m = re.match(r"^([0-9.]+)\^\+([0-9.]+)_-([0-9.]+)$", s)
    if m:
        return float(m.group(1)), float(m.group(2)), float(m.group(3)), kind
    m = re.match(r"^([0-9.,]+)\\pm([0-9.]+)$", s)
    if m:
        v = float(m.group(1).replace(",", ""))
        e = float(m.group(2))
        return v, e, e, kind
    m = re.match(r"^([0-9.,]+)\s*-\s*([0-9.,]+)$", s)
    if m:
        lo, hi = float(m.group(1).replace(",", "")), float(m.group(2).replace(",", ""))
        return 0.5 * (lo + hi), 0.5 * (hi - lo), 0.5 * (hi - lo), "range" if kind == "measured" else kind
    return float(s.replace(",", "")), np.nan, np.nan, kind


def _rows(block: str):
    out = []
    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("%") or "multicolumn" in line or line.startswith("\\"):
            continue
        for part in line.split("\\\\"):
            part = part.strip()
            if part.count("&") >= 6:
                out.append([c.strip() for c in part.split("&")])
    return out


#: The two tables spell a few names differently.
ALIASES = {"SDSS0634": "SDSSJ0634", "4U 1820-30": "4U1820-30"}


def _name(cell: str) -> str:
    s = re.sub(r"\$\^\{?[0-9,*]+\}?\$", "", cell)  # reference superscripts
    s = s.replace("$^{*7,8}$", "").replace("$-$", "-").replace("$^\\circ$", "d").replace("--", "-")
    s = s.replace("\\,", " ").replace("$", "").replace("~", " ")
    s = re.sub(r"\^\{?[0-9,*]+\}?", "", s)
    s = " ".join(s.split())
    return ALIASES.get(s, s)


def parse(tex_path=TEX) -> pd.DataFrame:
    text = tex_path.read_text()
    t1 = text[text.index("\\label{sec:sample}"):text.index("\\label{tab:system}")]
    t1 = t1[t1.index("Verification binaries"):]
    t2 = text[text.index("\\caption{Measured EM properties (parallax, distance) and derived GW parameters of the verification"):]
    t2 = t2[t2.index("Verification binaries"):t2.index("\\label{tab:system1}")]

    def split_blocks(t):
        i = t.index("Detectable binaries")
        return t[:i], t[i:]

    rec1 = []
    for status, block in zip(("verification", "detectable"), split_blocks(t1)):
        for c in _rows(block):
            rec1.append({"name": _name(c[0]), "status": status, "type_raw": _clean(c[1]),
                         "P_raw": c[2], "l": c[3], "b": c[4], "bp_rp": c[5], "MG": c[6],
                         "m1_raw": c[7], "m2_raw": c[8], "incl_raw": c[9].rstrip("\\").strip()})
    rec2 = []
    for status, block in zip(("verification", "detectable"), split_blocks(t2)):
        for c in _rows(block):
            rec2.append({"name": _name(c[0]), "status2": status, "f_mHz_table": c[1], "plx_dr3": c[2],
                         "L_pc": c[4], "d_raw": c[5], "sigA_over_A_pct": c[6],
                         "dincl_deg": c[7].rstrip("\\").strip()})
    d1 = pd.DataFrame(rec1)
    d2 = pd.DataFrame(rec2)
    # TIC 378898110 is listed in both blocks of both tables: keep it once, as detectable
    d1 = d1[~((d1["name"] == "TIC378898110") & (d1["status"] == "verification"))]
    d2 = d2[~((d2["name"] == "TIC378898110") & (d2["status2"] == "verification"))]
    df = d1.merge(d2, on="name", how="left", validate="one_to_one")
    assert (df["status"] == df["status2"]).all(), "table blocks disagree"
    rows = []
    for r in df.to_dict("records"):
        P, Perr = parse_number_with_paren(r["P_raw"])
        out = {"name": r["name"], "status": r["status"], "type": re.sub(r"\^?\*", "", r["type_raw"]),
               "star_flag": "*" in r["type_raw"], "P_s": P, "P_err_s": Perr,
               "l_deg": float(_clean(r["l"])), "b_deg": float(_clean(r["b"]).replace("-", "-")),
               "f_mHz_table": float(_clean(r["f_mHz_table"]))}
        for key, raw in (("m1", r["m1_raw"]), ("m2", r["m2_raw"]), ("incl", r["incl_raw"]),
                         ("d", r["d_raw"])):
            v, ep, em, kind = parse_value(raw)
            out[key], out[f"{key}_errp"], out[f"{key}_errm"], out[f"{key}_kind"] = v, ep, em, kind
        out["d_pc"] = out.pop("d")
        out["d_pc_errp"] = out.pop("d_errp")
        out["d_pc_errm"] = out.pop("d_errm")
        out["d_pc_kind"] = out.pop("d_kind")
        out["L_pc"] = _clean(r["L_pc"]) or np.nan
        out["sigA_over_A_pct"] = pd.to_numeric(_clean(r["sigA_over_A_pct"]), errors="coerce")
        out["dincl_deg"] = pd.to_numeric(_clean(r["dincl_deg"]), errors="coerce")
        rows.append(out)
    res = pd.DataFrame(rows)
    res["f_mHz"] = 2.0 / res["P_s"] * 1e3
    return res


def load(force: bool = False) -> pd.DataFrame:
    if OUT.exists() and not force:
        return pd.read_csv(OUT)
    df = parse()
    df.to_csv(OUT, index=False)
    return df


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    df = load(force=True)
    print(df[["name", "status", "type", "P_s", "f_mHz", "f_mHz_table", "m1", "m1_kind", "m2", "m2_kind",
              "incl", "incl_kind", "d_pc", "d_pc_errp", "d_pc_kind"]].to_string())
    print(df.groupby(["status", "type"]).size())
