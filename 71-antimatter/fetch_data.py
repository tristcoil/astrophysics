"""Download and cache everything the analysis needs, and load it back.

What this project reads, all public:

1. **Cosmic-ray measurements** from the Cosmic-Ray DataBase (CRDB, LPSC
   Grenoble): AMS-02's antiproton flux and antiproton-to-proton ratio
   (2016, 2021 and 2025 releases), the proton, helium, carbon, nitrogen and
   oxygen fluxes, the boron-to-carbon ratio, and every published limit on
   antihelium (BESS, PAMELA, AMS-01), with the other experiments' points.
2. **Saikawa & Shirai's (2018)** Standard-Model g*(T) table (hard-linked from
   P2's cache when that copy matches our checksum).
3. **arXiv sources** of papers whose figures or tables we read numerically
   (only the files we use are kept).
4. Small literature values are typed into ``config.py`` with their sources.

Every file's URL, size, SHA-256 and retrieval time goes into
``data/manifest.json``; nothing cached is fetched twice, so after one run the
project works offline.

    python fetch_data.py
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import os
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import config as cfg


# ---------------------------------------------------------------------------
# Manifest and downloads
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _load_manifest() -> dict:
    return json.loads(cfg.MANIFEST_PATH.read_text()) if cfg.MANIFEST_PATH.exists() else {}


def _save_manifest(manifest: dict) -> None:
    tmp = cfg.MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    tmp.replace(cfg.MANIFEST_PATH)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _get(url: str, retries: int = 4, timeout: int = 180) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last = err
            time.sleep(2 + 4 * attempt)
    raise RuntimeError(f"download failed: {url}: {last}")


def _record(manifest: dict, path: Path, url: str, what: str, **extra) -> None:
    key = str(path.relative_to(cfg.DATA_DIR))
    manifest[key] = {"url": url, "what": what, "bytes": path.stat().st_size,
                     "sha256": _sha256(path), "fetched": _now(), **extra}


def _write_bytes(path: Path, blob: bytes) -> None:
    tmp = path.with_suffix(path.suffix + ".part")
    tmp.write_bytes(blob)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# 1. CRDB
# ---------------------------------------------------------------------------


def crdb_url(num: str, den: str, combo_level: int, etype: str = "R") -> str:
    q = {"num": num, "energy_type": etype, "combo_level": combo_level,
         "energy_convert_level": 0, "format": "csv"}
    if den:
        q["den"] = den
    return cfg.CRDB_URL + "?" + urllib.parse.urlencode(q)


def fetch_crdb(manifest: dict) -> None:
    for stem, (num, den, combo, etype) in cfg.CRDB_QUERIES.items():
        path = cfg.CRDB_DIR / f"{stem}.csv"
        if path.exists() and str(path.relative_to(cfg.DATA_DIR)) in manifest:
            continue
        url = crdb_url(num, den, combo, etype)
        blob = _get(url)
        text = blob.decode("utf-8", "replace")
        if '"' not in text:                    # CRDB answers errors as bare text
            raise RuntimeError(f"CRDB query {stem} returned no data: {text[:200]}")
        _write_bytes(path, blob)
        _record(manifest, path, url, f"CRDB export: {num}{'/' + den if den else ''} vs {etype}")
        print(f"  CRDB {stem}: {path.stat().st_size / 1e3:.0f} kB")
        _save_manifest(manifest)
        time.sleep(1.0)


CRDB_COLUMNS = ["quantity", "subexp", "etype", "e_mean", "e_low", "e_high", "value",
                "stat_lo", "stat_hi", "syst_lo", "syst_hi", "ads", "phi_mv", "dist_au",
                "dates", "is_ul"]


def load_crdb(stem: str):
    """One CRDB export as a pandas DataFrame (numbers parsed, one row per bin)."""
    import csv

    import pandas as pd

    rows = []
    with open(cfg.CRDB_DIR / f"{stem}.csv", newline="") as f:
        for rec in csv.reader(line for line in f if not line.startswith("#")):
            if len(rec) == len(CRDB_COLUMNS):
                rows.append(rec)
    df = pd.DataFrame(rows, columns=CRDB_COLUMNS)
    for col in ("e_mean", "e_low", "e_high", "value", "stat_lo", "stat_hi", "syst_lo",
                "syst_hi", "phi_mv", "dist_au"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["is_ul"] = df["is_ul"].astype(int).astype(bool)
    df["err"] = np.hypot(0.5 * (df.stat_lo + df.stat_hi), 0.5 * (df.syst_lo + df.syst_hi))
    return df


# ---------------------------------------------------------------------------
# 2. Saikawa & Shirai g*(T) (shared with P2)
# ---------------------------------------------------------------------------


def _link_or_fetch(manifest: dict, path: Path, url: str, sha: str, what: str) -> None:
    key = str(path.relative_to(cfg.DATA_DIR))
    if path.exists() and key in manifest:
        return
    sibling = cfg.SIBLING_COPIES.get(path)
    if sibling is not None and sibling.exists() and _sha256(sibling) == sha:
        if path.exists():
            path.unlink()
        os.link(sibling, path)
        _record(manifest, path, url, what, hard_linked_from=str(sibling.relative_to(cfg.ROOT_DIR)))
        print(f"  {path.name}: hard-linked from {sibling.relative_to(cfg.ROOT_DIR)}")
    else:
        _write_bytes(path, _get(url))
        if _sha256(path) != sha:
            print(f"  warning: {path.name} checksum differs from P2's copy")
        _record(manifest, path, url, what)
        print(f"  {path.name}: downloaded")
    _save_manifest(manifest)


def load_gstar():
    """Saikawa & Shirai 2018: T [GeV], g*rho, g*s (arrays sorted by T)."""
    d = np.loadtxt(cfg.SS2018_FILE)
    return d[:, 0], d[:, 1], d[:, 3]


# ---------------------------------------------------------------------------
# 3. VizieR tables
# ---------------------------------------------------------------------------


def fetch_vizier(manifest: dict) -> None:
    for name, url in list(cfg.VIZIER_FILES.items()) + list(cfg.PDG_FILES.items()):
        path = (cfg.VIZIER_DIR if name in cfg.VIZIER_FILES else cfg.PDG_DIR) / name
        if path.exists() and str(path.relative_to(cfg.DATA_DIR)) in manifest:
            continue
        _write_bytes(path, _get(url))
        what = ("VizieR " + url.split("/ftp/")[1]) if "/ftp/" in url else "PDG " + name
        _record(manifest, path, url, what)
        print(f"  {what}")
        _save_manifest(manifest)


def load_fermi_igrb(model: str = "A") -> dict:
    """Fermi-LAT IGRB and total EGB (Ackermann et al. 2015, Table 3), one
    foreground model, as bin-averaged intensities per MeV."""
    cols = ["model", "emin", "emax", "igrb", "igrb_up", "igrb_lo", "igrb_fg_up", "igrb_fg_lo",
            "egb", "egb_up", "egb_lo", "egb_fg_up", "egb_fg_lo", "src", "src_up", "src_lo"]
    rows = [line.split() for line in (cfg.VIZIER_DIR / "J_ApJ_799_86_table3.dat").read_text()
            .splitlines() if line.strip()]
    rows = [r for r in rows if r[0] == model]
    a = np.array([[float(v) for v in r[1:]] for r in rows])
    d = {c: a[:, i] for i, c in enumerate(cols[1:])}
    width = d["emax"] - d["emin"]
    out = {"E_lo": d["emin"], "E_hi": d["emax"], "E_mid": np.sqrt(d["emin"] * d["emax"])}
    for key in ("igrb", "egb"):
        out[key] = d[key] / width
        out[key + "_err_up"] = np.hypot(d[key + "_up"], d[key + "_fg_up"]) / width
        out[key + "_err_lo"] = np.hypot(d[key + "_lo"], d[key + "_fg_lo"]) / width
        out[key + "_stat_up"] = d[key + "_up"] / width
    return out


# ---------------------------------------------------------------------------
# 4. arXiv sources
# ---------------------------------------------------------------------------


def fetch_arxiv(manifest: dict) -> None:
    for arxiv_id, files in cfg.ARXIV_SOURCES.items():
        outdir = cfg.ARXIV_DIR / arxiv_id.replace("/", "_")
        wanted = {outdir / local: inner for local, inner in files.items()}
        if all(p.exists() and str(p.relative_to(cfg.DATA_DIR)) in manifest for p in wanted):
            continue
        outdir.mkdir(parents=True, exist_ok=True)
        url = f"https://arxiv.org/e-print/{arxiv_id}"
        blob = _get(url)
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tar:
            members = {m.name: m for m in tar.getmembers() if m.isfile()}
            for path, inner in wanted.items():
                if inner not in members:
                    raise RuntimeError(f"{arxiv_id}: {inner} not in source tarball")
                _write_bytes(path, tar.extractfile(members[inner]).read())
                _record(manifest, path, url, f"arXiv:{arxiv_id} source file {inner}",
                        source_tarball_bytes=len(blob))
        print(f"  arXiv:{arxiv_id}: {len(wanted)} files")
        _save_manifest(manifest)
        time.sleep(3.0)                          # arXiv asks for gentle crawling


# ---------------------------------------------------------------------------
# 5. PRyMordial (code, pinned commit)
# ---------------------------------------------------------------------------


def fetch_prymordial(manifest: dict) -> None:
    import zipfile

    zpath = cfg.EXTERNAL_DIR / f"PRyMordial-{cfg.PRYMORDIAL_COMMIT[:12]}.zip"
    key = str(zpath.relative_to(cfg.DATA_DIR))
    if not (zpath.exists() and key in manifest):
        _write_bytes(zpath, _get(cfg.PRYMORDIAL_URL))
        _record(manifest, zpath, cfg.PRYMORDIAL_URL,
                "PRyMordial source at a pinned commit (github.com/vallima/PRyMordial)")
        _save_manifest(manifest)
        print(f"  PRyMordial {cfg.PRYMORDIAL_COMMIT[:12]}: {zpath.stat().st_size / 1e6:.1f} MB")
    if not cfg.PRYMORDIAL_DIR.exists():
        with zipfile.ZipFile(zpath) as z:
            for m in z.namelist():                         # refuse path traversal
                if m.startswith("/") or ".." in Path(m).parts:
                    raise RuntimeError(f"unsafe path in PRyMordial zip: {m}")
            z.extractall(cfg.EXTERNAL_DIR)


# ---------------------------------------------------------------------------


def fetch_all() -> None:
    manifest = _load_manifest()
    print("Fetching data (cached files are skipped)")
    fetch_crdb(manifest)
    _link_or_fetch(manifest, cfg.SS2018_FILE, cfg.SS2018_URL, cfg.SS2018_SHA256,
                   "Saikawa & Shirai 2018: T [GeV], g*rho, err, g*s, err")
    fetch_vizier(manifest)
    fetch_arxiv(manifest)
    fetch_prymordial(manifest)
    _save_manifest(manifest)


if __name__ == "__main__":
    fetch_all()


def load_pdg_xs(kind: str = "total"):
    """PDG p-bar p cross-section compilation: p_lab (GeV/c), sigma (mb),
    combined error (mb). Fixed-width rows after the header block."""
    rows = []
    for line in (cfg.PDG_DIR / f"rpp2022-pbarp_{kind}.dat").read_text().splitlines():
        parts = line.split()
        if len(parts) < 10 or not parts[0].isdigit():
            continue
        try:
            plab, sig = float(parts[1]), float(parts[4])
            st = 0.5 * (float(parts[5]) + float(parts[6]))
            sy = 0.5 * (float(parts[7]) + float(parts[8])) / 100 * sig
        except ValueError:
            continue
        rows.append((plab, sig, math.hypot(st, sy)))
    a = np.array(rows)
    return a[np.argsort(a[:, 0])]
