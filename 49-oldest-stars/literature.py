"""Published tables read from the arXiv LaTeX sources (comparison only, never fit inputs).

All return dicts keyed by the HUGS-style identifier ("ngc0104", "ngc6341", ...).
"""

from __future__ import annotations

import re

import config as C


def _ngc_id(num: str) -> str:
    return f"ngc{int(num):04d}"


def vandenberg2013() -> dict:
    """VandenBerg et al. (2013) Table 2: {id: (age, internal error)}.

    The quoted errors are *relative* (internal) uncertainties, typically
    +-0.25-0.5 Gyr; the paper puts the absolute uncertainty at ~+-1.5-2 Gyr.
    """
    text = (C.LIT_DIR / "vandenberg2013.tex").read_text()
    start = text.index("Ages and Other Properties of the Globular Cluster Sample")
    body = text[start:text.index("\\enddata", start)]
    out = {}
    for m in re.finditer(r"^\s*(?:\\phantom\{0\})?\s*(\d{3,4})\s*&[^&]*&\s*\$([-+\d.]+)\$\s*&\s*"
                         r"(?:\\phantom\{1\})?\$?([\d.]+)\\pm\s*([\d.]+)\$", body, flags=re.M):
        out[_ngc_id(m.group(1))] = (float(m.group(3)), float(m.group(4)))
    return out


def vandenberg2013_feh() -> dict:
    text = (C.LIT_DIR / "vandenberg2013.tex").read_text()
    start = text.index("Ages and Other Properties of the Globular Cluster Sample")
    body = text[start:text.index("\\enddata", start)]
    return {_ngc_id(m.group(1)): float(m.group(2)) for m in
            re.finditer(r"^\s*(?:\\phantom\{0\})?\s*(\d{3,4})\s*&[^&]*&\s*\$([-+\d.]+)\$", body, flags=re.M)}


def valcin2025() -> dict:
    """Valcin et al. (2025) Table 5: {id: (age, statistical error)} plus extras via valcin2025_full."""
    return {k: (v["age"], v["age_err"]) for k, v in valcin2025_full().items()}


def valcin2025_full() -> dict:
    text = (C.LIT_DIR / "valcin2025.tex").read_text()
    start = text.index("Cluster\\_name & Age (Gyr)")
    out = {}
    for ln in text[start:].splitlines()[1:]:
        if "\\end{tabular}" in ln:
            break
        cells = [c.strip() for c in ln.replace("\\hline", "").split("\\\\")[0].split("&")]
        if len(cells) < 17:
            continue
        try:
            vals = [float(c) for c in cells[1:17]]
        except ValueError:
            continue
        name = cells[0]
        out[name] = {"age": vals[0], "age_err": vals[1], "feh": vals[2], "dist_kpc": vals[4],
                     "ebv": vals[6], "alpha": vals[8], "Y": vals[12]}
    return out


def milone2017() -> dict:
    """Milone et al. (2017) Table 2: RGB widths and the fraction of first-generation stars."""
    text = (C.LIT_DIR / "milone2017.tex").read_text()
    out = {}
    pm = r"([\d.]+)\$\\pm\$([\d.]+)"
    for m in re.finditer(r"^NGC\\,(\d{4})\s*&\s*" + pm + r"\s*&\s*" + pm + r"\s*&\s*(?:" + pm + r")?\s*&\s*(?:"
                         + pm + r")?\s*&\s*" + pm + r"\s*&\s*(I{1,2})", text, flags=re.M):
        g = m.groups()
        out[_ngc_id(g[0])] = {"W_C": (float(g[1]), float(g[2])), "W_275_814": (float(g[3]), float(g[4])),
                              "f1G": (float(g[9]), float(g[10])), "type": g[11]}
    return out


def massari2019() -> dict:
    """Massari, Koppelman & Helmi (2019) Table 1: {id: progenitor code}."""
    text = (C.LIT_DIR / "massari2019.tex").read_text()
    start = text.index("List of Galactic GCs and associated progenitors")
    body = text[start:text.index("\\end{tabular}", start)]
    out = {}
    for m in re.finditer(r"([A-Za-z][\w~ ]*?\d+[A-Za-z]?)\s*&\s*([A-Za-z][A-Za-z0-9-]*)", body):
        name = m.group(1).replace("~", " ").strip()
        num = re.fullmatch(r"NGC\s*(\d+)", name)
        key = _ngc_id(num.group(1)) if num else name.lower().replace(" ", "")
        out[key] = m.group(2)
    return out


#: Massari et al. progenitor codes, grouped as they do: formed in the Milky Way's
#: main disc or bulge; brought in by a known dwarf galaxy; or on low/high-energy
#: orbits not tied to a known progenitor ("unassociated"; the low-energy group
#: has since been linked to an early accretion, Kraken/Heracles).
IN_SITU = {"M-D", "M-B"}
ACCRETED = {"G-E", "Sag", "H99", "Seq"}
UNASSOCIATED = {"L-E", "H-E"}


def origin(code: str | None) -> str:
    if code in IN_SITU:
        return "in situ"
    if code in ACCRETED:
        return "accreted"
    if code in UNASSOCIATED:
        return "unassociated"
    return "unknown"


PROGENITOR_NAMES = {
    "M-D": "main disc", "M-B": "main bulge", "G-E": "Gaia-Enceladus", "Sag": "Sagittarius",
    "H99": "Helmi streams", "Seq": "Sequoia", "L-E": "low-energy group", "H-E": "high-energy group",
    "XXX": "no kinematics",
}
