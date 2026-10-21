"""Published DSHARP tables, parsed from the papers' own LaTeX sources (arXiv).

Parsing instead of transcribing keeps typing errors out of every comparison. Four tables:

* DSHARP I (Andrews et al. 2018) Table 1, host-star properties: distance (Gaia DR2),
  log L*, log M* (from evolutionary tracks), age, accretion rate.
* DSHARP I Table 4, fiducial continuum images: beam, rms noise, peak, total flux.
* DSHARP II (Huang et al. 2018) Table 1, the annular substructures: every gap ("D") and
  ring ("B") with its radius (mas and au), the ellipse-fit geometry where measured, the
  method (E = ellipse fit, R = radial-profile extremum, V = visual), width and depth.
* DSHARP II Table 2, disc geometries (centre offset, inclination, PA, R_dust).
* DSHARP VII (Zhang et al. 2018) Table 3, planet masses inferred from 12 gap widths
  (plus 7 shallow gaps matched by eye), for three alphas and three grain-size models.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from config import PAPER_DIR

NAME_TO_KEY = {
    "AS 209": "AS209", "DoAr 25": "DoAr25", "DoAr 33": "DoAr33", "Elias 20": "Elias20", "Elias 24": "Elias24",
    "Elias 27": "Elias27", "GW Lup": "GWLup", "HD 142666": "HD142666", "HD 143006": "HD143006",
    "HD 163296": "HD163296", "IM Lup": "IMLup", "MY Lup": "MYLup", "RU Lup": "RULup", "SR 4": "SR4",
    "Sz 114": "Sz114", "Sz 129": "Sz129", "WaOph 6": "WaOph6", "WSB 52": "WSB52", "AS 205": "AS205",
    "HT Lup": "HTLup", "TW Hya": "TWHya", "HL Tau": "HLTau",
}


# ---------------------------------------------------------------------------
# LaTeX cell cleaning
# ---------------------------------------------------------------------------
@dataclass
class Value:
    """A number from a table cell: value, lower/upper 1-sigma errors and flags."""

    value: float | None
    err_lo: float | None = None
    err_hi: float | None = None
    approx: bool = False  # "~"
    upper_limit: bool = False  # "<"
    lower_limit: bool = False  # ">"
    dagger: bool = False  # DSHARP II: width narrower than the beam's minor axis
    raw: str = ""

    @property
    def err(self) -> float | None:
        if self.err_lo is None or self.err_hi is None:
            return None
        return 0.5 * (self.err_lo + self.err_hi)


def _strip(cell: str) -> str:
    s = cell
    s = re.sub(r"\\tablenotemark\{[^}]*\}", "", s)
    s = re.sub(r"\\(phn|phd)", "", s)
    s = s.replace("\\text{", "{").replace("\\,", " ")
    s = s.replace("{", "").replace("}", "") if "substack" not in s and "^" not in s else s
    return s.strip()


_NUM = r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"


def parse_value(cell: str) -> Value:
    raw = cell.strip()
    s = re.sub(r"\\tablenotemark\{[^}]*\}", "", raw)
    s = re.sub(r"\\(phn|phd)", "", s).replace("\\text{", "").replace("\\,", "")
    dagger = "dagger" in s
    s = s.replace("^\\dagger", "").replace("^{\\dagger}", "").replace("\\dagger", "")
    s = s.replace("$", "").strip()
    if s in ("", "-", "--", "\\nodata", "\\ldots"):
        return Value(None, raw=raw)
    approx = "\\sim" in s
    upper = s.startswith("<") or s.startswith("\\lt")
    lower = s.startswith(">")
    s = s.replace("\\sim", "").replace("<", "").replace(">", "").strip()
    # asymmetric errors: 35.6\substack{+0.7 \\ -0.8} or 116 ^{+19}_{-10}
    m = re.match(rf"^({_NUM})\s*\\substack\{{\s*\+?({_NUM})\s*\\\\\s*-?({_NUM})\s*\}}", s)
    if m:
        return Value(float(m.group(1)), abs(float(m.group(3))), abs(float(m.group(2))), approx, upper, lower,
                     dagger, raw)
    m = re.match(rf"^({_NUM})\s*\^\{{\s*\+?({_NUM})\s*\}}\s*_\{{\s*-?({_NUM})\s*\}}", s)
    if m:
        return Value(float(m.group(1)), abs(float(m.group(3))), abs(float(m.group(2))), approx, upper, lower,
                     dagger, raw)
    m = re.match(rf"^({_NUM})\s*\\pm\s*({_NUM})", s)
    if m:
        e = abs(float(m.group(2)))
        return Value(float(m.group(1)), e, e, approx, upper, lower, dagger, raw)
    m = re.match(rf"^({_NUM})", s.replace(" ", ""))
    if m:
        return Value(float(m.group(1)), None, None, approx, upper, lower, dagger, raw)
    return Value(None, raw=raw)


def _table_body(tex: str, label: str) -> str:
    """The rows between \\startdata and \\enddata of the deluxetable carrying ``label``."""
    i = tex.index(label)
    start = tex.rfind("\\begin{deluxetable", 0, i)
    body = tex[start:tex.index("\\end{deluxetable", i)]
    return body[body.index("\\startdata") + len("\\startdata"):body.index("\\enddata")]


def _rows(body: str) -> list[list[str]]:
    lines = [ln for ln in body.splitlines() if not ln.strip().startswith("%")]
    text = " ".join(lines).replace("\\hline", "")
    rows = []
    for row in re.split(r"\\\\(?![^{]*\})", text):  # split on row ends, not on \\ inside \substack{...}
        if row.strip():
            rows.append([c.strip() for c in row.split("&")])
    return rows


@lru_cache(maxsize=None)
def _tex(arxiv_base: str) -> str:
    folder = PAPER_DIR / arxiv_base
    files = sorted(folder.glob("*.tex"), key=lambda p: p.stat().st_size, reverse=True)
    if not files:
        raise FileNotFoundError(f"no LaTeX source in {folder}; run fetch_data.py")
    return files[0].read_text(errors="replace")


# ---------------------------------------------------------------------------
# DSHARP I
# ---------------------------------------------------------------------------
@dataclass
class Star:
    name: str
    key: str
    region: str
    tmass: str
    distance: Value
    spt: str
    log_teff: Value
    log_lstar: Value
    log_mstar: Value
    log_age: Value
    log_mdot: Value

    @property
    def lstar(self) -> float:
        return 10.0 ** self.log_lstar.value

    @property
    def mstar(self) -> float:
        return 10.0 ** self.log_mstar.value


def _signed(cell: str) -> Value:
    """DSHARP I writes negative logs as '-$0.48\\pm0.20$': move the sign inside."""
    c = re.sub(r"\\(phn|phd)", "", cell).strip()
    neg = c.startswith("-$")
    if neg:
        c = "$-" + c[2:]
    c = c.replace("$<$", "<")
    return parse_value(c)


@lru_cache(maxsize=None)
def dsharp1_stars() -> dict[str, Star]:
    tex = _tex("1812.04040")
    out = {}
    for c in _rows(_table_body(tex, "\\label{table:stars}")):
        if len(c) < 10:
            continue
        name = re.sub(r"\\tablenotemark\{[^}]*\}", "", c[0]).strip()
        out[NAME_TO_KEY[name]] = Star(
            name=name, key=NAME_TO_KEY[name], region=c[1].strip(), tmass=c[2].strip(),
            distance=parse_value(c[3]), spt=c[4].strip(), log_teff=parse_value(c[5]), log_lstar=_signed(c[6]),
            log_mstar=_signed(c[7]), log_age=_signed(c[8]), log_mdot=_signed(c[9]),
        )
    return out


@dataclass
class ImageProps:
    name: str
    key: str
    freq_ghz: float
    beam_mas: tuple[float, float]
    beam_pa: float
    rms_ujy: float
    rms_k: float
    peak_mjy: float
    peak_k: float
    flux_mjy: float
    robust: float
    taper: str


@lru_cache(maxsize=None)
def dsharp1_images() -> dict[str, ImageProps]:
    tex = _tex("1812.04040")
    out = {}
    for c in _rows(_table_body(tex, "\\label{table:images}")):
        if len(c) < 8:
            continue
        cl = [re.sub(r"\\(phn|phd)", "", x).replace("$", "").strip() for x in c]
        name = cl[0]
        bm = re.match(rf"({_NUM})\s*\\times\s*({_NUM}),\s*({_NUM})", cl[2])
        rms = [float(x) for x in re.findall(_NUM, cl[3])]
        peak = [float(x) for x in re.findall(_NUM, cl[4])]
        out[NAME_TO_KEY[name]] = ImageProps(
            name=name, key=NAME_TO_KEY[name], freq_ghz=float(cl[1]),
            beam_mas=(float(bm.group(1)), float(bm.group(2))), beam_pa=float(bm.group(3)),
            rms_ujy=rms[0], rms_k=rms[1], peak_mjy=peak[0], peak_k=peak[1], flux_mjy=float(cl[5]),
            robust=float(cl[6].replace("\\phn", "")), taper=cl[7],
        )
    return out


# ---------------------------------------------------------------------------
# DSHARP II
# ---------------------------------------------------------------------------
@dataclass
class Feature:
    key: str  # disc
    label: str  # e.g. "D61"
    kind: str  # "D" (gap) or "B" (ring)
    dx_mas: Value
    dy_mas: Value
    r_mas: Value
    r_au: Value
    inc: Value
    pa: Value
    method: str  # E, R or V
    width_au: Value
    depth: Value
    note: str = ""

    @property
    def name(self) -> str:
        return f"{self.key}:{self.label}"


@lru_cache(maxsize=None)
def dsharp2_features() -> list[Feature]:
    tex = _tex("1812.04041")
    feats = []
    current = None
    for c in _rows(_table_body(tex, "\\label{tab:ringpositions}")):
        if len(c) < 9:
            continue
        name = re.sub(r"\\tablenotemark\{[^}]*\}", "", c[0]).strip()
        if name:
            current = NAME_TO_KEY[name]
        label_raw = c[1].strip()
        note = "".join(re.findall(r"\\tablenotemark\{([^}]*)\}", label_raw))
        label = re.sub(r"\\tablenotemark\{[^}]*\}", "", label_raw).strip()
        method = c[8].replace("\\text{", "").replace("}", "").replace("$", "").strip()
        feats.append(Feature(
            key=current, label=label, kind=label[0], dx_mas=parse_value(c[2]), dy_mas=parse_value(c[3]),
            r_mas=parse_value(c[4]), r_au=parse_value(c[5]), inc=parse_value(c[6]), pa=parse_value(c[7]),
            method=method or "-", width_au=parse_value(c[9]) if len(c) > 9 else Value(None),
            depth=parse_value(c[10]) if len(c) > 10 else Value(None), note=note,
        ))
    return feats


@dataclass
class Geometry:
    key: str
    dx_mas: Value
    dy_mas: Value
    inc: Value
    pa: Value
    method: str  # E (ellipse fits to rings) or G (2-D Gaussian)
    r_dust_au: Value
    beam_min_au: float


@lru_cache(maxsize=None)
def dsharp2_geometry() -> dict[str, Geometry]:
    tex = _tex("1812.04041")
    out = {}
    for c in _rows(_table_body(tex, "\\label{tab:orientation}")):
        if len(c) < 8:
            continue
        key = NAME_TO_KEY[c[0].strip()]
        out[key] = Geometry(key, parse_value(c[1]), parse_value(c[2]), parse_value(c[3]), parse_value(c[4]),
                            c[5].strip(), parse_value(c[6]), float(c[7]))
    return out


# ---------------------------------------------------------------------------
# DSHARP VII
# ---------------------------------------------------------------------------
@dataclass
class ZhangPlanet:
    key: str
    mstar: float
    r_gap_au: float
    width: float  # Delta = (r_out - r_in)/r_out
    sigma_dust: float  # g/cm^2, averaged over 1.1-2 r_gap
    h_over_r: float
    st_used: tuple[float | None, float | None, float | None]  # St_max in units of 5.23e-4: DSD1, "1 mm", DSD2
    masses: dict = field(default_factory=dict)  # (model, alpha) -> M_p in M_Jup
    log_unc: dict = field(default_factory=dict)  # model -> (-err, +err) in log10 M_p
    from_fit: bool = True  # upper part of the table (fitting formula) vs. matched by eye
    note: str = ""


def _triple(cell: str) -> list[float | None]:
    s = cell.replace("$", "").replace("(", "").replace(")", "")
    s = re.sub(r"\\tablenotemark\{[^}]*\}", "", s)
    toks = re.findall(r">?\s*\d+(?:\.\d+)?|--", s)
    vals = []
    for t in toks:
        vals.append(None if t == "--" else float(t.replace(">", "").strip()))
    while len(vals) < 3:
        vals.append(None)
    return vals[:3]


@lru_cache(maxsize=None)
def dsharp7_planets() -> list[ZhangPlanet]:
    tex = _tex("1812.04045")
    body = _table_body(tex, "\\label{table:Mp}")
    upper, _, lower = body.partition("\\hline")
    out = []
    for part, from_fit in ((upper, True), (lower, False)):
        for c in _rows(part):
            if len(c) < 14:
                continue
            name_raw = c[0]
            note = "".join(re.findall(r"\\tablenotemark\{([^}]*)\}", name_raw))
            key = NAME_TO_KEY[re.sub(r"\\tablenotemark\{[^}]*\}", "", name_raw).strip()]
            p = ZhangPlanet(key=key, mstar=float(c[1]), r_gap_au=float(c[2]), width=float(c[3]),
                            sigma_dust=float(c[4]), h_over_r=float(c[5]), st_used=tuple(_triple(c[9])),
                            from_fit=from_fit, note=note)
            for alpha, col in ((1e-4, 10), (1e-3, 11), (1e-2, 12)):
                for model, m in zip(("DSD1", "1mm", "DSD2"), _triple(c[col])):
                    if m is not None:
                        p.masses[(model, alpha)] = m
            # one entry per model, "--" where there is none; some cells read "$^{+ 0.16}_{-0.14}$"
            for model, entry in zip(("DSD1", "1mm", "DSD2"), c[13].split(",")):
                m = re.search(rf"\^\{{\s*\+\s*({_NUM})\s*\}}\s*_\{{\s*-\s*({_NUM})\s*\}}", entry)
                if m:
                    p.log_unc[model] = (float(m.group(2)), float(m.group(1)))
            out.append(p)
    return out


def summary() -> str:
    lines = []
    stars = dsharp1_stars()
    lines.append(f"DSHARP I Table 1: {len(stars)} stars")
    imgs = dsharp1_images()
    lines.append(f"DSHARP I Table 4: {len(imgs)} images")
    feats = dsharp2_features()
    lines.append(f"DSHARP II Table 1: {len(feats)} features in {len({f.key for f in feats})} discs")
    geo = dsharp2_geometry()
    lines.append(f"DSHARP II Table 2: {len(geo)} geometries")
    pl = dsharp7_planets()
    lines.append(f"DSHARP VII Table 3: {sum(p.from_fit for p in pl)} fitted gaps + {sum(not p.from_fit for p in pl)} by eye")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
