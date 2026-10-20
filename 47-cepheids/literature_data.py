"""Published tables we need as *inputs* (not answer keys): Cepheid metallicities.

* C-MetaLL IX (Ripepi et al. 2026, A&A 708, A216), VizieR J/A+A/708/A216:
  290 Milky Way classical Cepheids with homogeneous high-resolution [Fe/H]
  from one team's spectra and one analysis, keyed by Gaia DR3 source_id.
* Trentin et al. 2024 (C-MetaLL IV, A&A 681, A65), VizieR J/A+A/681/A65:
  910 Cepheids compiled from the literature (635 with high-resolution
  spectroscopy, flag f_[Fe/H] = 0; 275 from Gaia RVS spectra), keyed by
  Gaia DR3 source_id. Heterogeneous, so used as a cross-check.

Note (verified against the VizieR ReadMe): in both tables the parallax
"correction" column holds the Lindegren et al. zero point itself, not a
corrected parallax. We do not use those columns; our zero points are our own.
"""

from __future__ import annotations

import io

import numpy as np

import config as C

VIZIER_TSV = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source={src}&-out.all&-out.max=unlimited"

TABLES = {
    "cmetall9_tablea1.tsv": ("J/A+A/708/A216/tablea1",
                             "C-MetaLL IX (Ripepi et al. 2026) Table A.1: homogeneous [Fe/H] of 290 MW Cepheids"),
    "trentin2024_table1.tsv": ("J/A+A/681/A65/table1",
                               "Trentin et al. 2024 (C-MetaLL IV) Table 1: [Fe/H] compilation of 910 MW Cepheids"),
}


def fetch(manifest: dict, download) -> None:
    print("Literature tables (VizieR) ...")
    for name, (src, note) in TABLES.items():
        download(VIZIER_TSV.format(src=src), C.LIT_DIR / name, manifest, note)


def read_vizier_tsv(name: str) -> dict:
    """Parse a VizieR ASU-TSV file: '#' comments, header, units line, dashes line, then rows."""
    text = (C.LIT_DIR / name).read_text(errors="replace").splitlines()
    rows = [ln for ln in text if ln and not ln.startswith("#")]
    header = rows[0].split("\t")
    body = rows[3:]  # skip units and dashes
    cols = {h: [] for h in header}
    for ln in body:
        parts = ln.split("\t")
        if len(parts) != len(header):
            continue
        for h, v in zip(header, parts):
            cols[h].append(v.strip())
    return cols


def _floats(vals):
    out = []
    for v in vals:
        try:
            out.append(float(v))
        except ValueError:
            out.append(np.nan)
    return np.array(out)


def metallicities() -> dict:
    """{'cmetall9': {source_id: ([Fe/H], err)}, 'trentin24': {source_id: ([Fe/H], err, flag)}}."""
    out = {}
    t = read_vizier_tsv("cmetall9_tablea1.tsv")
    ids = [int(s) if s.isdigit() else 0 for s in t["GaiaDR3"]]
    feh, err = _floats(t["[Fe/H]"]), _floats(t["e_[Fe/H]"])
    out["cmetall9"] = {i: (f, e) for i, f, e in zip(ids, feh, err) if i > 0 and np.isfinite(f)}
    t = read_vizier_tsv("trentin2024_table1.tsv")
    ids = [int(s) if s.isdigit() else 0 for s in t["GaiaDR3"]]
    feh = _floats(t["[Fe/H]"])
    err = _floats(t["e_[Fe/H]"]) if "e_[Fe/H]" in t else np.full(len(ids), np.nan)
    flag = [s for s in t.get("f_[Fe/H]", [""] * len(ids))]
    out["trentin24"] = {i: (f, e, g) for i, f, e, g in zip(ids, feh, err, flag) if i > 0 and np.isfinite(f)}
    return out
