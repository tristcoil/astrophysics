"""BlackCAT, the catalogue of black hole X-ray transients: as published in
2016 (VizieR J/A+A/587/A61) and as maintained online today.

The online pages hold each dynamical parameter as text ("124±4", "0.25 – 0.31",
"≤79", "9⁺⁰·²₋₀.₆", two values from two papers on separate lines, ...) with
links to the papers on ADS. This module

1. parses every page into (value, references) entries without losing any;
2. turns the dynamically confirmed black holes into measurement
   distributions for the period P, velocity semi-amplitude K2, mass ratio q
   and inclination i, by fixed rules (``adopted_inputs``).

The rules, applied to every system alike:

* Where several values are listed, K2, q and i take the most recent paper;
  P takes the value from the same paper as the adopted K2 (so that the mass
  function uses one consistent orbit), else the most precise.
* ``a±b`` is a Gaussian, ``a⁺ᵇ₋ᶜ`` a split Gaussian. A range ``a – b`` is
  uniform for q and isotropic (cos i uniform) for i. An upper limit on q is
  U(0, limit). ``≈x`` (no error) is U(0.5x, 1.5x).
* An inclination given only as an upper limit gets the lower bound Kreidberg
  et al. (2012) derived from a donor of normal mass for its spectral type; a
  lower limit gets the no-eclipse upper bound for each drawn q.
* A missing q is U(0, 0.5), Farr et al.'s (2011) convention. A system whose
  inclination has no lower bound (none given, or an upper limit that
  Kreidberg et al. did not bound from below) has no upper bound on its mass:
  it gets only a minimum mass and is left out of the population fits.
* Two overrides, where BlackCAT's mass ratio disagrees with the paper it
  cites or with the measurement behind it, documented in ``OVERRIDES``.
"""

from __future__ import annotations

import html
import re
from pathlib import Path

import numpy as np
import pandas as pd

from fetch_data import BLACKCAT_DIR, CDS_DIR, blackcat_object_ids, blackcat_page_path
from literature import COMMON_NAME, KREIDBERG2012, KREIDBERG_NAME

FIELDS = {
    "Orbital period": "porb_h",
    "Star's spectral type": "sptype",
    "Star's radial velocity amplitude": "k2",
    "Mass function": "fm",
    "BH mass": "m1",
    "Mass ratio": "q",
    "Inclination angle": "incl",
    "Rot. broadening": "vsini",
}

#: Explicit departures from the rules, each with its reason.
OVERRIDES = {
    "1H 1659-487": {
        "q": (("N", 0.18, 0.05),
              "Heida et al. 2017 (ApJ 846, 132) measure q = 0.18 +- 0.05 from v sin i; "
              "BlackCAT lists it as an upper limit because the K2 behind it is a lower limit"),
    },
    "SAX J1819.3-2525": {
        "q": (("N", 0.455, 0.041),
              "MacDonald et al. 2014 (ApJ 784, 2) find M_BH/M_2 = 2.2 +- 0.2, i.e. "
              "q = 0.45 +- 0.04; BlackCAT cites that paper for q = 0.63-0.70, which is the "
              "older solution (q = 0.67 +- 0.04, Orosz et al. 2001, as listed by Kreidberg et al. 2012)"),
    },
}

_NUM = r"\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"
_FLAG = r"(?P<flag>[<>≤≥~∼≈]?)\s*"
_PATTERNS = [
    ("range", re.compile(rf"^{_FLAG}(?P<a>{_NUM})\s*(?:[–—]|\s-\s)\s*(?P<b>{_NUM})(?P<note>.*)$")),
    ("asym", re.compile(rf"^{_FLAG}(?P<v>{_NUM})\s*\+\s*(?P<p>{_NUM})\s*-\s*(?P<m>{_NUM})(?P<note>.*)$")),
    ("pm", re.compile(rf"^{_FLAG}(?P<v>{_NUM})\s*±\s*(?P<e>{_NUM})(?P<note>.*)$")),
    ("plain", re.compile(rf"^{_FLAG}(?P<v>{_NUM})(?P<note>.*)$")),
]
_TAG = re.compile(r"</?(?:sup|sub|a|span|strong|b|i|em|font|p|div|center|h\d)\b[^>]*>", re.I)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def clean_cell(cell_html: str) -> list[str]:
    """Cell HTML -> text lines. Asymmetric errors (<sup>+a</sup><sub>-b</sub>)
    become ' +a -b'. Only known tags are stripped: BlackCAT writes upper
    limits as a literal '<70', which a generic tag stripper would delete."""
    s = cell_html.replace("&plusmn;<sup>", "<sup>")
    s = re.sub(r"<sup>\s*\+?\s*([^<]+?)\s*</sup>\s*<sub>\s*[-−]?\s*([^<]+?)\s*</sub>",
               r" +\1 -\2", s)
    lines = re.split(r"<br\s*/?>", s)
    out = []
    for line in lines:
        line = html.unescape(_TAG.sub("", line))
        line = line.replace("−", "-").replace("\xa0", " ")
        line = re.sub(r"\s+", " ", line).strip()
        if line:
            out.append(line)
    return out


def parse_value(text: str) -> dict | None:
    """One BlackCAT value string -> a dict (None if it is not numeric)."""
    for kind, pat in _PATTERNS:
        m = pat.match(text)
        if not m:
            continue
        d = m.groupdict()
        out = {"kind": kind, "flag": d["flag"].replace("∼", "~"), "note": d["note"].strip()}
        if kind == "range":
            out.update(low=float(d["a"]), high=float(d["b"]))
        elif kind == "asym":
            out.update(value=float(d["v"]), plus=float(d["p"]), minus=float(d["m"]))
        elif kind == "pm":
            out.update(value=float(d["v"]), plus=float(d["e"]), minus=float(d["e"]))
        else:
            out.update(value=float(d["v"]))
        return out
    return None


def _year(bibcode: str) -> int:
    m = re.match(r"(\d{4})", bibcode)
    return int(m.group(1)) if m else 0


def parse_page(path: Path) -> dict:
    s = path.read_text(encoding="utf-8", errors="replace")
    obj_id = int(re.search(r"ID#\s*(\d+)", s).group(1))
    name = html.unescape(re.search(r"<h2><strong>(.*?)</h2>", s, re.S).group(1)).strip()
    cp = re.search(r"</h2></strong>\s*<h4>(.*?)</h4>", s, re.S)
    counterpart = html.unescape(_TAG.sub("", cp.group(1))).strip() if cp else ""
    out = {"id": obj_id, "name": name, "counterpart": counterpart,
           "dynamical": "Dynamically confirmed BH" in s, "fields": {}}
    if "Dynamical parameters" not in s:
        return out
    block = s[s.find("Dynamical parameters"):s.find("Optical/NIR magnitudes")]
    for m in re.finditer(r"<th[^>]*>(.*?)</th>\s*<td>(.*?)</td>\s*<td>(.*?)</td>", block, re.S):
        label = html.unescape(_TAG.sub("", m.group(1)))
        key = next((v for k, v in FIELDS.items() if label.strip().startswith(k)), None)
        if key is None:
            continue
        lines = clean_cell(m.group(2))
        refs = [(urlbib.replace("%26", "&"), txt.strip())
                for urlbib, txt in re.findall(r"abs/([^'\"]+)['\"][^>]*>([^<]+)<", m.group(3))]
        entries = []
        for k, line in enumerate(lines):
            if len(lines) == len(refs):
                line_refs = [refs[k]]
            elif len(lines) == 1:
                line_refs = refs
            else:  # more refs than lines: first line <- first ref, rest <- the others
                line_refs = refs[:1] if k == 0 else refs[1:]
            entries.append({"raw": line, "parsed": parse_value(line) if key != "sptype" else None,
                            "refs": line_refs,
                            "year": max((_year(b) for b, _ in line_refs), default=0)})
        out["fields"][key] = entries
    return out


def load_online() -> list[dict]:
    """Every object of the online catalogue (73 as of October 2026)."""
    return [parse_page(blackcat_page_path(i)) for i in blackcat_object_ids()]


def load_vizier_a4() -> pd.DataFrame:
    """Table A4 of Corral-Santana et al. (2016): the 18 dynamical systems of 2016
    (byte positions from the CDS ReadMe)."""
    cols = [("Name", 0, 18), ("n_Name", 19, 20), ("SpType", 32, 39), ("Porb_h", 40, 51),
            ("e_Porb", 52, 61), ("l_K2", 62, 63), ("K2", 63, 68), ("e_K2", 69, 73),
            ("fM", 75, 80), ("e_fM", 81, 85), ("l_M1", 86, 88), ("M1", 88, 93),
            ("E_M1", 94, 97), ("e_M1", 98, 101), ("M1u", 102, 106)]
    rows = []
    for line in (CDS_DIR / "tablea4.dat").read_text().splitlines():
        if not line.strip():
            continue
        line = line.ljust(252)
        rec = {}
        for name, a, b in cols:
            v = line[a:b].strip()
            if name in ("Name", "n_Name", "SpType", "l_K2", "l_M1"):
                rec[name] = v
            else:
                try:
                    rec[name] = float(v)
                except ValueError:
                    rec[name] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# From catalogue entries to measurement distributions
# ---------------------------------------------------------------------------

def _latest(entries: list[dict]) -> dict | None:
    numeric = [e for e in entries if e["parsed"] is not None]
    if not numeric:
        return None
    best = max(e["year"] for e in numeric)
    return [e for e in numeric if e["year"] == best][-1]


def _gauss(p: dict):
    if "plus" not in p:
        return None
    if abs(p["plus"] - p["minus"]) < 1e-12:
        return ("N", p["value"], p["plus"])
    return ("SN", p["value"], p["plus"], p["minus"])


def _q_spec(p: dict | None):
    if p is None:
        return ("U", 0.0, 0.5), "missing: U(0, 0.5), the convention of Farr et al. 2011"
    if p["kind"] == "range":
        return ("U", p["low"], p["high"]), "range: uniform"
    if p["flag"] in ("≤", "<") and p["kind"] in ("plain", "pm"):
        return ("U", 0.0, p["value"]), "upper limit: U(0, limit)"
    if p["flag"] in ("≈", "~") and "plus" not in p:
        return ("U", 0.5 * p["value"], 1.5 * p["value"]), "approximate: U(0.5x, 1.5x)"
    g = _gauss(p)
    if g is not None:
        return g, "measured"
    raise ValueError(f"cannot interpret q {p}")


def _i_spec(p: dict | None, kreidberg_name: str | None):
    """Returns (spec, note). spec may be ('I_ecl', i_lo) meaning isotropic from
    i_lo up to the no-eclipse limit of each drawn q, or None (unconstrained)."""
    if p is None:
        return None, "no constraint: minimum mass only"
    if p["kind"] == "range":
        return ("I", p["low"], p["high"]), "range: isotropic"
    if p["flag"] in ("≤", "<"):
        if kreidberg_name in KREIDBERG2012 and KREIDBERG2012[kreidberg_name][3] > 0:
            i_min = KREIDBERG2012[kreidberg_name][3]
            return ("I", i_min, p["value"]), (f"upper limit; lower bound {i_min:g} deg "
                                              "(donor of normal mass, Kreidberg et al. 2012)")
        return ("max", p["value"]), "upper limit only: minimum mass only"
    if p["flag"] in ("≥", ">") and "plus" not in p:
        return ("I_ecl", p["value"]), "lower limit; upper bound = no-eclipse limit"
    g = _gauss(p)
    if g is not None:
        return g, "measured"
    raise ValueError(f"cannot interpret i {p}")


def adopted_inputs(systems: list[dict] | None = None) -> list[dict]:
    """Measurement distributions for every dynamically confirmed black hole."""
    systems = load_online() if systems is None else systems
    out = []
    for s in systems:
        if not s["dynamical"]:
            continue
        f = s["fields"]
        k2 = _latest(f.get("k2", []))
        k2_refs = {b for b, _ in k2["refs"]}
        p_entries = [e for e in f.get("porb_h", []) if e["parsed"] is not None]
        same = [e for e in p_entries if {b for b, _ in e["refs"]} & k2_refs]
        porb = same[0] if same else min(p_entries, key=lambda e: e["parsed"].get("plus", np.inf))
        q_entry, i_entry = _latest(f.get("q", [])), _latest(f.get("incl", []))
        notes = []
        if k2["parsed"]["flag"] in (">", "≥"):
            notes.append("K2 quoted as a lower limit; used as measured")
        q_spec, q_note = _q_spec(q_entry["parsed"] if q_entry else None)
        kname = KREIDBERG_NAME.get(s["name"])
        i_spec, i_note = _i_spec(i_entry["parsed"] if i_entry else None, kname)
        for field, (spec, reason) in OVERRIDES.get(s["name"], {}).items():
            if field == "q":
                q_spec, q_note = spec, f"override: {reason}"
        fm_pub, m1_pub = _latest(f.get("fm", [])), _latest(f.get("m1", []))
        pp, kp = porb["parsed"], k2["parsed"]
        out.append({
            "id": s["id"], "name": s["name"], "common_name": COMMON_NAME.get(s["name"], ""),
            "counterpart": s["counterpart"],
            "porb_h": (pp["value"], pp.get("plus", 0.0)), "porb_ref": porb["refs"],
            "k2": (kp["value"], kp.get("plus", 0.0)), "k2_ref": k2["refs"],
            "q_spec": q_spec, "q_note": q_note, "q_ref": q_entry["refs"] if q_entry else [],
            "i_spec": i_spec, "i_note": i_note, "i_ref": i_entry["refs"] if i_entry else [],
            "fm_published": fm_pub["parsed"] if fm_pub else None,
            "fm_published_raw": fm_pub["raw"] if fm_pub else "",
            "fm_ref": fm_pub["refs"] if fm_pub else [],
            "m1_published_raw": m1_pub["raw"] if m1_pub else "",
            "m1_published": m1_pub["parsed"] if m1_pub else None,
            "m1_ref": m1_pub["refs"] if m1_pub else [],
            "sptype": " / ".join(e["raw"] for e in f.get("sptype", [])),
            "notes": notes,
        })
    return out


def raw_table(systems: list[dict]) -> pd.DataFrame:
    """Every dynamical entry as text, with its references, for data/ output."""
    rows = []
    for s in systems:
        if not s["dynamical"]:
            continue
        for key, entries in s["fields"].items():
            for e in entries:
                rows.append({"id": s["id"], "name": s["name"], "field": key, "value": e["raw"],
                             "bibcodes": "; ".join(b for b, _ in e["refs"]),
                             "ref_keys": "; ".join(k for _, k in e["refs"])})
    return pd.DataFrame(rows)
