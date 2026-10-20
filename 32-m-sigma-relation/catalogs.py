"""Read the black hole compilations and comparison samples into tidy tables.

Every loader returns a pandas DataFrame with the same core columns, so the
analysis never needs to know which paper a number came from:

    name        canonical galaxy name ("NGC 4486", "Milky Way", "Mrk 335")
    logM        log10 black hole mass (Msun)
    elogM       symmetric 1-sigma error on logM (dex), the mean of the upper
                and lower errors in log space
    upper       True if logM is an upper limit (then logM is the limit)
    logsig      log10 velocity dispersion (km/s)
    elogsig     its 1-sigma error (dex)
    method      "stars", "gas", "maser", "reverb" (how the mass was measured)
    dist        distance (Mpc)

plus source-specific columns. Parsed tables are cached as CSV in data/cache/.

Kormendy & Ho's Table 3 needs care: it exists only as a dvips PostScript file
in their arXiv source, and the authors encode which galaxies they fitted by
colour (red = classical bulge, blue = pseudobulge, cyan = mass not used in
fits). We render it with Ghostscript, read the text with poppler's pdftotext
and read each word's colour off the rendered pixels.
"""

from __future__ import annotations

import io
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import config as C

# ---------------------------------------------------------------------------
# Galaxy names
# ---------------------------------------------------------------------------

#: Messier numbers of the galaxies that appear in black hole compilations.
MESSIER_NGC = {
    31: 224, 32: 221, 33: 598, 49: 4472, 51: 5194, 59: 4621, 60: 4649, 64: 4826, 77: 1068,
    81: 3031, 82: 3034, 84: 4374, 85: 4382, 86: 4406, 87: 4486, 89: 4552, 94: 4736,
    96: 3368, 104: 4594, 105: 3379, 106: 4258, 110: 205,
}
_PREFIX = {"N": "NGC", "NGC": "NGC", "IC": "IC", "U": "UGC", "UGC": "UGC", "MRK": "Mrk",
           "ARK": "Ark", "FAIRALL": "Fairall", "PGC": "PGC", "ESO": "ESO"}
_ALIASES = {"GALAXY": "Milky Way", "MW": "Milky Way", "MILKYWAY": "Milky Way",
            "CYGNUSA": "Cygnus A", "CYGA": "Cygnus A", "CIRCINUS": "Circinus",
            "CENA": "NGC 5128", "3C120": "3C 120", "3C390.3": "3C 390.3", "ARP151": "Mrk 40"}


def canonical(name: str) -> str:
    """One spelling per galaxy: 'N4486 (M87)', 'M{\\ts}87', 'NGC4486' -> 'NGC 4486'."""
    s = re.sub(r"\$\^\{?[^$]*\}?\$", "", str(name))  # TeX footnote marks like $^{b,x}$
    s = s.replace("{\\ts}", " ").replace("\\ts", " ").replace("~", " ")
    s = re.sub(r"\(.*?\)", "", s).strip()  # 'N4486 (M87)' -> 'N4486'
    key = re.sub(r"[\s_\-]+", "", s).upper()
    if key in _ALIASES:
        return _ALIASES[key]
    m = re.fullmatch(r"(A\d{4})-?BCG", key)
    if m:
        return f"{m.group(1)} BCG"
    m = re.fullmatch(r"(NGC|IC|UGC|MRK|ARK|FAIRALL|PGC|N|U|M)0*(\d+)([A-Z]?)", key)
    if m:
        prefix, num, suffix = m.group(1), int(m.group(2)), m.group(3)
        if prefix == "M":
            if suffix == "" and num in MESSIER_NGC:
                return f"NGC {MESSIER_NGC[num]}"
            return f"M{num}{suffix}"
        return f"{_PREFIX[prefix]} {num}{suffix}"
    return s.strip()


def _log_errors(value, lo, hi):
    """Upper and lower errors in log10 for a value with 1-sigma bounds (lo, hi)."""
    up = np.log10(hi / value)
    with np.errstate(divide="ignore"):
        down = np.where(lo > 0, np.log10(value / np.where(lo > 0, lo, 1.0)), np.inf)
    return up, down


# ---------------------------------------------------------------------------
# van den Bosch (2016)
# ---------------------------------------------------------------------------

_VDB_METHOD = {"star": "stars", "stars": "stars", "gas": "gas", "CO": "gas", "maser": "maser",
               "reverb": "reverb", "GC": "cluster"}


def load_vdb16() -> pd.DataFrame:
    """The 294-object compilation as the author distributes it (FITS)."""
    from astropy.table import Table

    t = Table.read(C.RAW_DIR / "vdb16" / "BHcompilation.fits")
    df = t.to_pandas()
    for col in ("NAME", "HETNAME", "TYPE", "REF"):
        df[col] = df[col].str.decode("ascii") if df[col].dtype == object and isinstance(df[col].iloc[0], bytes) \
            else df[col].astype(str)
        df[col] = df[col].str.strip()
    missing = lambda v: np.where(np.asarray(v, float) <= -1000, np.nan, np.asarray(v, float))  # noqa: E731
    out = pd.DataFrame({
        "name": [canonical(n) for n in df["NAME"]],
        "name_orig": df["NAME"],
        "selected": df["SELECTED"].astype(int),
        "dist": missing(df["DDIST"]),
        "edist": missing(df["DISTERR"]),
        "logM": missing(df["MBH"]),
        "elogM": missing(df["DMBH"]),
        "upper": df["UPPERLIMIT"].astype(float) > 0.5,
        "logsig": missing(df["SIG"]),
        "elogsig": missing(df["DSIG"]),
        "logLk": missing(df["LK"]),
        "elogLk": missing(df["DLK"]),
        "logRe": missing(df["R50"]),  # half-light radius, log kpc
        "elogRe": missing(df["DR50"]),
        "method": [_VDB_METHOD.get(m, m) for m in df["TYPE"]],
        "type_orig": df["TYPE"],
        "ref": df["REF"],
        "het_sig": missing(df["HETMGSSIG"]),
    })
    # Globular clusters (selected = -4) are not galaxies; masses of 0 there
    # mean "no constraint". Keep them out of everything.
    out = out[out["selected"] != -4].reset_index(drop=True)
    # For upper limits the compilation stores the middle of the allowed range
    # and its half-width (VizieR lists the range as logM -logM / +dlogM), so the
    # limit itself is logM + dlogM: M33 (NGC 598) < 10^3.15 Msun, as published
    # (Gebhardt et al. 2001: < 1500 Msun); NGC 4382 < 10^8.35, as in Kormendy &
    # Ho. The fits use the compilation's own convention (centre +/- half-width,
    # censored), as van den Bosch did; ``logM_limit`` is for plots and checks.
    out["logM_limit"] = np.where(out["upper"], out["logM"] + out["elogM"], np.nan)
    out["source"] = "vdB16"
    return out


def load_vdb16_vizier(table: str = "table2") -> pd.DataFrame:
    """The journal's version of the compilation (VizieR), for a consistency check."""
    lines = [ln for ln in (C.RAW_DIR / "vdb16" / f"vizier_{table}.tsv").read_text().splitlines()
             if ln and not ln.startswith("#")]
    head = lines[0].split("\t")
    body = [ln.split("\t") for ln in lines[3:]]
    df = pd.DataFrame(body, columns=head)
    for c in df.columns:
        df[c] = df[c].str.strip()
        conv = pd.to_numeric(df[c], errors="coerce")
        if conv.notna().sum() >= 0.5 * df[c].ne("").sum() and c not in ("Name", "Meth", "Refs", "SimbadName"):
            df[c] = conv
    df["name"] = [canonical(n) for n in df["Name"]]
    return df


# ---------------------------------------------------------------------------
# Kormendy & Ho (2013), Table 2: ellipticals (TeX)
# ---------------------------------------------------------------------------

def _tex_clean(s: str) -> str:
    s = s.replace("\\0", "").replace("\\1", "").replace("\\x", "").replace("\\y", "")
    s = s.replace("{\\ts}", " ").replace("\\ts", " ").replace("~", " ").replace("$", "")
    s = re.sub(r"\\textBlack|\\textRed|\\ARDiscard", "", s)
    return re.sub(r"\s+", " ", s).strip()


def _num(s: str) -> float:
    s = s.strip()
    return np.nan if s in ("", "\\dots", "...", "…") else float(s)


def parse_kh13_table2(tex: str) -> pd.DataFrame:
    rows = []
    for line in tex.splitlines():
        s = line.strip()
        if not s.startswith("&") or "\\end" not in s or "\\times 10" not in s:
            continue
        cells = s.split("\\end")[0].split("&")[1:]
        if len(cells) != 13:
            raise ValueError(f"KH13 Table 2: expected 13 cells, got {len(cells)}: {s[:80]}")
        name_cell, typ, dist, ks, mkt, mvt, vk, bv, mbul, mass, sig, flags, src = cells
        row_cyan = name_cell.strip().startswith("\\y")
        merger = name_cell.strip().startswith("\\g")
        mass_cyan = "\\y" in mass
        name = _tex_clean(name_cell.replace("\\g", ""))
        d_tokens = _tex_clean(dist).split()
        m = re.search(r"([\d.]+)\s*\(\s*([\d.]+)\s*-\s*([\d.]+)\s*\)\s*\\times\s*10\^\{?(\d+)\}?",
                      _tex_clean(mass))
        if not m:
            raise ValueError(f"KH13 Table 2: cannot read mass {mass!r}")
        exp = 10.0 ** int(m.group(4))
        mb = re.match(r"([\d.]+)\s*\\pm\s*([\d.]+)", _tex_clean(mbul))
        sc = _tex_clean(sig)
        sm = re.match(r"([\d.]+)\s*\\pm\s*([\d.]+)", sc)
        sa = re.match(r"([\d.]+)\^\{\+([\d.]+)\}_\{-([\d.]+)\}", sc)
        if sm:
            sig_v, sig_hi, sig_lo = float(sm.group(1)), float(sm.group(2)), float(sm.group(2))
        elif sa:
            sig_v, sig_hi, sig_lo = float(sa.group(1)), float(sa.group(2)), float(sa.group(3))
        else:
            sig_v = sig_hi = sig_lo = np.nan
        f = flags.split()
        rows.append({
            "name": canonical(name), "name_orig": name, "table": 2, "hubble_type": _tex_clean(typ),
            "bulge": "elliptical", "dist": float(d_tokens[0].rstrip(".")),
            "Ks": _num(_tex_clean(ks)), "MKs_bulge": _num(_tex_clean(mkt)),
            "logMbulge": float(mb.group(1)) if mb else np.nan,
            "elogMbulge": float(mb.group(2)) if mb else np.nan,
            "M": float(m.group(1)) * exp, "M_lo": float(m.group(2)) * exp, "M_hi": float(m.group(3)) * exp,
            "sigma": sig_v, "sigma_hi": sig_hi, "sigma_lo": sig_lo,
            "method_flag": int(f[0]), "core": int(f[1]), "corrected": int(f[2]),
            "source_ref": _tex_clean(src),
            "row_color": "cyan" if row_cyan else ("green" if merger else "black"),
            "mass_color": "cyan" if (mass_cyan or row_cyan) else ("green" if merger else "black"),
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Kormendy & Ho (2013), Table 3: spirals and S0s (PostScript, colour-coded)
# ---------------------------------------------------------------------------

#: Ghostscript page set-up that keeps the whole landscape table on the page:
#: dvips placed part of it at negative coordinates, beyond the letter page.
_GS_MEDIA = ["-dDEVICEWIDTHPOINTS=900", "-dDEVICEHEIGHTPOINTS=1100", "-dFIXEDMEDIA"]
_GS_SHIFT = "<</Install {150 250 translate}>> setpagedevice"
_RENDER_DPI = 200

#: Reference colours of the rendered glyphs (sRGB), measured on the render.
_COLOURS = {
    "red": (238, 38, 46),  # classical bulges
    "blue": (44, 76, 145),  # pseudobulges
    "cyan": (14, 149, 198),  # "not included in the correlation fits"
    "black": (40, 38, 39),
    "green": (40, 150, 40),
}


def _require(tool: str) -> str:
    path = shutil.which(tool)
    if path is None:
        raise RuntimeError(f"`{tool}` is needed to read Kormendy & Ho's Table 3 "
                           "(install ghostscript and poppler-utils)")
    return path


def render_kh13_table3(ps_path: Path, workdir: Path) -> tuple[Path, Path, Path]:
    """PostScript -> rotated PDF -> (layout text, word boxes, PNG)."""
    workdir.mkdir(parents=True, exist_ok=True)
    pdf = workdir / "kh13_table3.pdf"
    subprocess.run([_require("gs"), "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=pdfwrite", *_GS_MEDIA,
                    "-dAutoRotatePages=/PageByPage", f"-sOutputFile={pdf}", "-c", _GS_SHIFT,
                    "-f", str(ps_path)], check=True)
    txt = workdir / "kh13_table3.txt"
    subprocess.run([_require("pdftotext"), "-layout", str(pdf), str(txt)], check=True)
    bbox = workdir / "kh13_table3_bbox.html"
    subprocess.run([_require("pdftotext"), "-bbox", str(pdf), str(bbox)], check=True)
    stem = workdir / "kh13_table3"
    subprocess.run([_require("pdftoppm"), "-r", str(_RENDER_DPI), "-png", "-singlefile", str(pdf),
                    str(stem)], check=True)
    return txt, bbox, stem.with_suffix(".png")


def _word_colour(img: np.ndarray, box) -> str:
    """Name of the reference colour nearest to the darkest third of a word's pixels."""
    s = _RENDER_DPI / 72.0
    x0, y0, x1, y1 = (int(round(v * s)) for v in box)
    px = img[y0:y1 + 1, x0:x1 + 1].reshape(-1, 3).astype(float)
    ink = px[px.min(axis=1) < 200]
    if len(ink) == 0:
        return "none"
    ink = ink[np.argsort(ink.sum(axis=1))[: max(1, len(ink) // 3)]].mean(axis=0)
    return min(_COLOURS, key=lambda k: np.sum((ink - np.array(_COLOURS[k])) ** 2))


_T3_ROW = re.compile(
    r"^\s*(?P<name>(?:M|NGC|IC|UGC)\s+\d+[A-Z]?|Galaxy|Circinus)\s+(?P<type>\S+)\s+"
    r"(?P<dist>[\d.]+)(?:\s+(?P<dsrc>\d+))?\s+(?P<Ks>[\d.]+|\.\.\.)\s+"
    r"(?P<MKsT>-[\d.]+)\s+(?P<MKsb>-[\d.]+)\s+(?P<MKsd>-[\d.]+)\s+(?P<BT>[\d.]+)\s+(?P<PBT>[\d.]+)\s+"
    r"(?P<MVT>-[\d.]+|\.\.\.)\s+(?P<MVb>-[\d.]+|\.\.\.)\s+(?P<MVd>-[\d.]+|\.\.\.)\s+"
    r"(?P<VK>[\d.]+|\.\.\.)\s+(?P<BV>[\d.]+|\.\.\.)\s+(?P<logMb>[\d.]+)\s*±\s*(?P<elogMb>[\d.]+)\s+"
    r"(?P<M>[\d.]+)\s*\(\s*(?P<lo>[\d.]+)\s*-\s*(?P<hi>[\d.]+)\s*\)\s*×\s*10(?P<exp>\d+)\s+"
    r"(?P<sig>\d+)\s*±\s*(?P<esig>\d+)\s+(?:(?P<vc>\d+)\s*±\s*(?P<evc>\d+)|\.\.\.)\s+"
    r"(?P<flags>[\d,]+(?:\s+\d)?)\s+(?P<src>.+?)\s*$")


def parse_kh13_table3(txt_path: Path, bbox_path: Path, png_path: Path) -> pd.DataFrame:
    from PIL import Image

    img = np.asarray(Image.open(png_path).convert("RGB"))
    words = [(float(a), float(b), float(c), float(d), w) for a, b, c, d, w in re.findall(
        r'xMin="([\d.-]+)" yMin="([\d.-]+)" xMax="([\d.-]+)" yMax="([\d.-]+)">([^<]*)<',
        bbox_path.read_text())]
    # Column anchors from the "(1)" ... "(20)" header row.
    head = [w for w in words if re.fullmatch(r"\(\d+\)", w[4]) and int(w[4][1:-1]) <= 20 and w[1] < 300]
    anchors = {int(w[4][1:-1]): 0.5 * (w[0] + w[2]) for w in head}
    x_left = next(w[0] for w in head if w[4] == "(1)")  # names are left-aligned with "(1)"
    y_head = next(w[1] for w in head if w[4] == "(1)")
    # The mass cell starts just right of the log M_bulge error and its first
    # word (the mantissa, e.g. "1.43(1.12" or "1.7") begins left of anchor 16.
    x_mass = (anchors[15] + 15.0, anchors[16])
    name_words = sorted([w for w in words if abs(w[0] - x_left) < 3 and w[1] > y_head + 5
                         and re.fullmatch(r"M|NGC|IC|UGC|Galaxy|Circinus", w[4])], key=lambda w: w[1])
    coloured = []
    for nw in name_words:
        row = sorted([w for w in words if abs(w[1] - nw[1]) < 2.5], key=lambda w: w[0])
        label = nw[4] if nw[4] in ("Galaxy", "Circinus") else f"{row[0][4]} {row[1][4]}"
        mass_words = [w for w in row if x_mass[0] < w[0] < x_mass[1]]
        coloured.append((label, _word_colour(img, nw[:4]),
                         _word_colour(img, mass_words[0][:4]) if mass_words else "none"))

    lines = [ln.replace("−", "-") for ln in txt_path.read_text().splitlines()]
    data_lines = [ln for ln in lines if re.match(r"^\s*(?:M|NGC|IC|UGC)\s+\d+|^\s*(Galaxy|Circinus)\s", ln)
                  and "±" in ln]
    if len(data_lines) != len(coloured):
        raise ValueError(f"KH13 Table 3: {len(data_lines)} text rows but {len(coloured)} coloured names")
    rows = []
    block = "classical"  # the table lists classical bulges first, then pseudobulges
    for ln, (label, c_name, c_mass) in zip(data_lines, coloured):
        m = _T3_ROW.match(ln)
        if not m:
            raise ValueError(f"KH13 Table 3: cannot parse row {ln!r}")
        g = m.groupdict()
        if re.sub(r"\s+", " ", g["name"]) != label:
            raise ValueError(f"KH13 Table 3: text row {g['name']!r} does not line up with {label!r}")
        # Red rows are classical bulges and blue rows pseudobulges; a cyan row
        # (a galaxy whose mass is not fitted) belongs to the block it sits in.
        if c_name in ("red", "blue"):
            block = {"red": "classical", "blue": "pseudo"}[c_name]
        exp = 10.0 ** int(g["exp"])
        flags = g["flags"].replace(" ", "")
        rows.append({
            "name": canonical(g["name"]), "name_orig": g["name"], "table": 3, "hubble_type": g["type"],
            "bulge": block,
            "dist": float(g["dist"].rstrip(".")),
            "Ks": _num(g["Ks"]), "MKs_bulge": float(g["MKsb"]),
            "BT": float(g["BT"]), "PBT": float(g["PBT"]),
            "logMbulge": float(g["logMb"]), "elogMbulge": float(g["elogMb"]),
            "M": float(g["M"]) * exp, "M_lo": float(g["lo"]) * exp, "M_hi": float(g["hi"]) * exp,
            "sigma": float(g["sig"]), "sigma_hi": float(g["esig"]), "sigma_lo": float(g["esig"]),
            "vcirc": float(g["vc"]) if g["vc"] else np.nan,
            "evcirc": float(g["evc"]) if g["evc"] else np.nan,
            "method_flag": int(flags[:-1].rstrip(",").split(",")[0]),
            "method_flags_all": flags[:-1].rstrip(","),
            "corrected": int(flags[-1]), "core": 0,
            "source_ref": g["src"],
            "row_color": c_name, "mass_color": c_mass,
        })
    return pd.DataFrame(rows)


_KH_METHOD = {1: "stars", 2: "gas", 3: "maser", 4: "gas"}


def load_kh13() -> pd.DataFrame:
    """Kormendy & Ho (2013) Tables 2 and 3 in one table, with their fit flags."""
    cache = C.CACHE_DIR / "kh13.csv"
    if cache.exists():
        return pd.read_csv(cache, keep_default_na=True)
    t2 = parse_kh13_table2((C.RAW_DIR / "kh13" / "annrev4K-E.tex").read_text())
    txt, bbox, png = render_kh13_table3(C.RAW_DIR / "kh13" / "annrev4K-S.ps", C.CACHE_DIR / "kh13_render")
    t3 = parse_kh13_table3(txt, bbox, png)
    df = pd.concat([t2, t3], ignore_index=True)
    df["method"] = df["method_flag"].map(_KH_METHOD)
    df["upper"] = df["M_lo"] <= 0  # KH13 print upper limits as a range from zero
    df["logM"] = np.where(df["upper"], np.log10(df["M_hi"]), np.log10(df["M"]))
    up, down = _log_errors(df["M"].to_numpy(), df["M_lo"].to_numpy(), df["M_hi"].to_numpy())
    df["elogM_up"], df["elogM_down"] = up, down
    df["elogM"] = np.where(df["upper"], np.nan, 0.5 * (up + down))
    df["logsig"] = np.log10(df["sigma"])
    df["elogsig"] = 0.5 * (np.log10(df["sigma"] + df["sigma_hi"]) - np.log10(df["sigma"] - df["sigma_lo"]))
    df["merger"] = df["row_color"] == "green"
    df["source"] = "KH13"
    # KH13's fit sample (their eqs. 2-3, 10): ellipticals and classical bulges
    # whose masses are printed in black or red, i.e. not cyan (ionised gas
    # without line widths, NGC 2778, NGC 3607, the monsters NGC 1277 and
    # NGC 4486B, NGC 4342, and NGC 3842 / NGC 4889) and not green (mergers).
    df["kh13_fit"] = (df["bulge"].isin(["elliptical", "classical"])
                      & df["mass_color"].isin(["black", "red"]) & ~df["merger"] & ~df["upper"])
    df.to_csv(cache, index=False)
    return df


# ---------------------------------------------------------------------------
# McConnell & Ma (2013), Table A1 (TeX)
# ---------------------------------------------------------------------------

def parse_mm13(tex: str) -> pd.DataFrame:
    start = tex.index("\\label{tab:sample}")
    body = tex[start: tex.index("\\end{document}")]
    rows = []
    for line in body.splitlines():
        s = line.strip()
        if "&" not in s or not s.endswith("\\\\") or "\\bf" in s:
            continue
        cells = [c.strip() for c in s[:-2].split("&")]
        if len(cells) != 11:
            continue
        name, mass, ref, sig, lv, mbul, mref, rinf, morph, dist, method = cells
        mm = re.match(r"([\d.]+)\s*\(\s*([\d.]+)\s*,\s*([\d.]+)\s*\)\s*e(\d+)", mass)
        if not mm:  # the units line under each header
            continue
        notes = re.findall(r"\$\^\{?([a-z,]+)\}?\$", name)
        m_val, m_up, m_dn, exp = float(mm.group(1)), float(mm.group(2)), float(mm.group(3)), 10.0 ** int(mm.group(4))
        sp = re.match(r"\$\s*([\d.]+)\s*\\pm\s*([\d.]+)\s*\$", sig)
        sa = re.match(r"\$\s*([\d.]+)\^\{\+([\d.]+)\}_\{-([\d.]+)\}\s*\$", sig)
        if sp:
            s_v, s_hi, s_lo = float(sp.group(1)), float(sp.group(2)), float(sp.group(2))
        else:
            s_v, s_hi, s_lo = float(sa.group(1)), float(sa.group(2)), float(sa.group(3))
        lvm = re.match(r"\$\s*([\d.]+)\s*\\pm\s*([\d.]+)\s*\$", lv)
        rows.append({
            "name": canonical(name), "name_orig": re.sub(r"\$.*?\$", "", name).strip(),
            "notes": ",".join(notes),
            "M": m_val * exp, "M_up": m_up * exp, "M_dn": m_dn * exp,
            "mass_ref": ref, "sigma": s_v, "sigma_hi": s_hi, "sigma_lo": s_lo,
            "logLV": float(lvm.group(1)) if lvm else np.nan,
            "elogLV": float(lvm.group(2)) if lvm else np.nan,
            "Mbulge": float(mbul) if mbul else np.nan,
            "rinf_arcsec": float(rinf), "morph": morph, "dist": float(dist),
            "method": {"stars": "stars", "gas": "gas", "masers": "maser"}[method],
        })
    df = pd.DataFrame(rows)
    df["logM"] = np.log10(df["M"])
    df["elogM_up"] = np.log10((df["M"] + df["M_up"]) / df["M"])
    df["elogM_down"] = np.log10(df["M"] / (df["M"] - df["M_dn"]))
    df["elogM"] = 0.5 * (df["elogM_up"] + df["elogM_down"])
    df["logsig"] = np.log10(df["sigma"])
    df["elogsig"] = 0.5 * (np.log10(df["sigma"] + df["sigma_hi"]) - np.log10(df["sigma"] - df["sigma_lo"]))
    df["early_type"] = df["morph"].str.match(r"^(E|S0)")
    df["upper"] = False
    # NGC 1399 has two independent measurements; MM13 weight each by 1/2.
    df["weight"] = 1.0 / df.groupby("name")["name"].transform("size")
    df["source"] = "MM13"
    return df


def load_mm13() -> pd.DataFrame:
    return parse_mm13((C.RAW_DIR / "mm13" / "arxiv_v3.tex").read_text())


# ---------------------------------------------------------------------------
# AGN Black Hole Mass Database (Bentz & Katz 2015)
# ---------------------------------------------------------------------------

def load_agn(label: str = "virial_products") -> pd.DataFrame:
    """Reverberation masses: log (c tau dV^2 / G) for f = 1, or log M for another <f>."""
    text = (C.RAW_DIR / "agn" / f"agnmass_{label}.csv").read_text()
    f_applied = float(re.search(r"<f> value applied = ([\d.]+)", text).group(1))
    lines = text.splitlines()
    df = pd.read_csv(io.StringIO("\n".join(lines[2:])), na_values=["NAN"])
    df.columns = ["object", "ra", "dec", "z", "logM", "elogM_up", "elogM_down"]
    df["name"] = [canonical(n) for n in df["object"]]
    df["elogM"] = 0.5 * (df["elogM_up"] + df["elogM_down"])
    df["f_applied"] = f_applied
    df["method"] = "reverb"
    return df


# ---------------------------------------------------------------------------
# ATLAS3D (Cappellari et al. 2011a, 2013a)
# ---------------------------------------------------------------------------

def _read_hash_table(path: Path, names: list[str]) -> pd.DataFrame:
    rows = [ln.split() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    return pd.DataFrame(rows, columns=names)


def load_atlas3d() -> pd.DataFrame:
    """The 260 ATLAS3D early-type galaxies: distance, M_K, sigma_e, R_e."""
    samp = _read_hash_table(C.RAW_DIR / "atlas3d" / C.ATLAS3D_FILES["sample"],
                            ["galaxy", "ra", "dec", "sbf", "nedd", "virgo", "vhel", "dist", "MK",
                             "AB", "ttype", "logRe_arcsec_p1"])
    scal = _read_hash_table(C.RAW_DIR / "atlas3d" / C.ATLAS3D_FILES["scaling"],
                            ["galaxy", "logsig", "rsig_re", "logsig_kpc", "inc", "logML_jam",
                             "logVcmax", "quality", "logRe_maj", "logRe_arcsec", "logR12", "logRg",
                             "conc", "eps_e", "logLum_r"])
    df = samp.merge(scal, on="galaxy", how="inner")
    for c in df.columns:
        if c != "galaxy":
            df[c] = pd.to_numeric(df[c].replace("--", np.nan))  # "--": not measured
    df["name"] = [canonical(n) for n in df["galaxy"]]
    df["logLk"] = (C.MK_SUN - df["MK"]) / 2.5
    # R_e in kpc from the angular size and distance (1 arcsec = D[Mpc] * 4.848 pc).
    df["logRe"] = df["logRe_arcsec"] + np.log10(df["dist"] * 1e3 / C.ARCSEC_PER_RAD)
    return df


# ---------------------------------------------------------------------------
# Convenience
# ---------------------------------------------------------------------------

def load_all() -> dict[str, pd.DataFrame]:
    return {
        "vdb16": load_vdb16(),
        "kh13": load_kh13(),
        "mm13": load_mm13(),
        "agn_vp": load_agn("virial_products"),
        "agn_grier": load_agn("grier13"),
        "atlas3d": load_atlas3d(),
    }
