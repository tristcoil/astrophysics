"""Download and cache every real object and measurement the light cone is drawn with.

Nine small groups of files, all public:

1. **Davis & Lineweaver (2004)**, arXiv source: the LaTeX and the EPS figures
   whose curves we read back as the check (``dl2004.py``).
2. **DESI DR2 BAO** (13 numbers + covariance): distances at seven redshifts.
3. **Pantheon+** (1,701 supernova light curves + covariance): distances, and the
   test of the special-relativistic reading of redshift.
4. **GWOSC** cumulative event list (GWTC-1 to GWTC-5): gravitational-wave mergers.
5. **JADES DR4** NIRSpec redshifts (Curtis-Lake et al. 2025): only the first FITS
   extension (4.4 MB of the 90 MB file), fetched by byte range.
6. **SDSS DR16Q** quasar redshifts (Lyke et al. 2020), 750,414 quasars.
7. **Banados et al. (2016)** quasars at z > 5.6, from VizieR.
8. **SDSS Main Galaxy Sample** redshifts (D2's selection), 656,795 galaxies.

Sibling projects of this atlas already hold several of these (D5's BAO table,
D1's Pantheon+, I3's DR16Q columns, D2's galaxies); when they are present the
files are hard-linked or read in place instead of downloaded again, and the
manifest says so. Every file lands in ``data/`` with its URL, size, SHA-256 and
retrieval time; nothing cached is fetched twice, so after one run everything
works offline.

    python fetch_data.py   # ~5.5 MB when the sibling caches exist, ~85 MB when nothing is cached
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402

import config as C  # noqa: E402

# --- Plumbing ------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def load_manifest() -> dict:
    return json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {}


def save_manifest(manifest: dict) -> None:
    merged = load_manifest()
    merged.update(manifest)
    C.MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    tmp = C.MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(merged, indent=2, sort_keys=True))
    tmp.replace(C.MANIFEST_PATH)


def _record(manifest: dict, dest: Path, url: str, **extra) -> None:
    key = str(dest.relative_to(C.DATA_DIR))
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                     "retrieved_utc": extra.pop("retrieved_utc", _now()), **extra}


def http_get(url: str, *, headers: dict | None = None, retries: int = 5,
             timeout: float = 300.0) -> bytes:
    """GET with retries on the transient failures this machine sees (DNS, resets)."""
    hdr = {"User-Agent": C.USER_AGENT, **(headers or {})}
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=hdr),
                                        timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if isinstance(err, urllib.error.HTTPError) and err.code in (400, 403, 404):
                raise
            if attempt == retries - 1:
                raise
            print(f"    retrying after {err!r}", flush=True)
            time.sleep(2 + 4 * attempt)
    raise RuntimeError("unreachable")


def download(url: str, dest: Path, manifest: dict, *, headers: dict | None = None,
             refresh: bool = False) -> Path:
    """Fetch one file unless cached; atomic write; provenance in the manifest."""
    if dest.exists() and not refresh:
        key = str(dest.relative_to(C.DATA_DIR))
        if key not in manifest:
            mtime = datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc)
            _record(manifest, dest, url, retrieved_utc=mtime.isoformat(timespec="seconds"))
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    data = http_get(url, headers=headers)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)
    _record(manifest, dest, url)
    print(f"  fetched {dest.relative_to(C.DATA_DIR)}  ({len(data) / 1e6:.2f} MB)", flush=True)
    return dest


def link_or_download(rel: str, url: str, manifest: dict) -> Path:
    """Hard-link a sibling project's identical cached file if it exists, else download."""
    dest = C.DATA_DIR / rel
    sibling = C.ROOT / C.SIBLING_FILES[rel] if rel in C.SIBLING_FILES else None
    if not dest.exists() and sibling is not None and sibling.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        try:
            os.link(sibling, dest)
        except OSError:
            shutil.copy2(sibling, dest)
        _record(manifest, dest, url, linked_from=str(sibling.relative_to(C.ROOT)))
        print(f"  linked {rel} from {sibling.relative_to(C.ROOT)}", flush=True)
        return dest
    return download(url, dest, manifest)


# --- 1. Davis & Lineweaver 2004 ---------------------------------------------------------

DL2004_MEMBERS = ("Davis_Lineweaver03v2.tex", "fig1a.eps", "fig1b.eps", "fig1c.eps",
                  "fig2.eps", "fig5.eps")


def fetch_dl2004(manifest: dict) -> Path:
    tgz = download(C.DL2004_ARXIV_SRC, C.LIT_DIR / "dl2004_src.tar.gz", manifest)
    out = C.LIT_DIR / "dl2004_src"
    if all((out / m).exists() for m in DL2004_MEMBERS):
        return out
    out.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tgz, "r:gz") as tar:
        for m in tar.getmembers():  # extract only the plain files we need, by bare name
            if m.isfile() and m.name in DL2004_MEMBERS:
                with tar.extractfile(m) as src, open(out / m.name, "wb") as dst:
                    shutil.copyfileobj(src, dst)
    missing = [m for m in DL2004_MEMBERS if not (out / m).exists()]
    if missing:
        raise RuntimeError(f"arXiv source lacks {missing}")
    return out


# --- 5. JADES DR4: one FITS extension by byte range --------------------------------------


def _fits_header_end(buf: bytes, start: int) -> tuple[int, dict]:
    """End offset of the FITS header that starts at ``start``, and its key cards."""
    cards = {}
    pos = start
    while pos + 2880 <= len(buf):
        block = buf[pos:pos + 2880]
        pos += 2880
        for i in range(0, 2880, 80):
            card = block[i:i + 80].decode("ascii", "replace")
            key = card[:8].strip()
            if key == "END":
                return pos, cards
            if "=" in card[:10]:
                cards[key] = card[10:].split("/")[0].strip().strip("'").strip()
    raise ValueError("header continues beyond the bytes fetched")


def fetch_jades(manifest: dict) -> Path:
    dest = C.RAW_DIR / "jades" / "jades_dr4_obs_info.fits"
    if dest.exists():
        return download(C.JADES_DR4_URL, dest, manifest)  # records provenance if needed
    head = http_get(C.JADES_DR4_URL, headers={"Range": "bytes=0-57599"})
    end0, primary = _fits_header_end(head, 0)
    if int(primary.get("NAXIS", "0")) != 0:
        raise RuntimeError("JADES primary HDU unexpectedly has data")
    end1, ext = _fits_header_end(head, end0)
    if ext.get("EXTNAME") != "Obs_info":
        raise RuntimeError(f"first extension is {ext.get('EXTNAME')!r}, expected Obs_info")
    nbytes = int(ext["NAXIS1"]) * int(ext["NAXIS2"]) + int(ext.get("PCOUNT", "0"))
    stop = end1 + -(-nbytes // 2880) * 2880  # data padded to whole 2880-byte blocks
    body = http_get(C.JADES_DR4_URL, headers={"Range": f"bytes=0-{stop - 1}"})
    if len(body) != stop:
        raise RuntimeError(f"JADES: got {len(body)} bytes, expected {stop}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(".part")
    tmp.write_bytes(body)
    tmp.replace(dest)
    _record(manifest, dest, C.JADES_DR4_URL, byte_range=f"0-{stop - 1}",
            note="primary HDU + Obs_info extension only (5190 rows)")
    print(f"  fetched JADES DR4 Obs_info ({stop / 1e6:.2f} MB of the 94 MB file)", flush=True)
    download(C.JADES_README_URL, C.RAW_DIR / "jades" / "Readme_DR4_catalogues.md", manifest)
    return dest


# --- VizieR -----------------------------------------------------------------------------


def parse_vizier_tsv(text: str) -> list[dict]:
    """Rows of a VizieR ASU TSV answer. Anchors on the dashes row (the units row may be
    blank, and skipping a fixed number of lines then drops the first data row)."""
    if "QUERY_STATUS=OVERFLOW" in text or "#Truncated" in text:
        raise RuntimeError("VizieR truncated the answer")
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    header = [h.strip() for h in lines[i].split("\t")]
    j = i + 1
    while j < len(lines) and not re.fullmatch(r"[\t ]*-[-\t ]*", lines[j]):
        j += 1
    rows = []
    for ln in lines[j + 1:]:
        if ln.strip():
            vals = ln.split("\t")
            rows.append({h: (vals[k].strip() if k < len(vals) else "") for k, h in enumerate(header)})
    return rows


def vizier_tsv(source: str, columns: list[str], dest: Path, manifest: dict,
               constraints: dict | None = None) -> Path:
    params = {"-source": source, "-out.max": "unlimited", "-out": ",".join(columns)}
    params.update(constraints or {})
    url = C.VIZIER_ASU + "?" + urllib.parse.urlencode(params)
    path = download(url, dest, manifest)
    if not path.read_text(encoding="utf-8", errors="replace").lstrip().startswith("#"):
        path.unlink()
        raise RuntimeError(f"unexpected VizieR answer for {source}")
    return path


# --- 6. DR16Q --------------------------------------------------------------------------

DR16Q_N = 750_414


def fetch_dr16q(manifest: dict, strips: int = 12) -> Path:
    dest = C.DERIVED_DIR / "dr16q_z.npz"
    if dest.exists():
        return dest
    sibling = C.ROOT / C.SIBLING_DR16Q
    if sibling.exists():
        import pyarrow.parquet as pq

        t = pq.read_table(sibling, columns=["RAJ2000", "DEJ2000", "z"])
        ra, dec, z = (t.column(c).to_numpy() for c in ("RAJ2000", "DEJ2000", "z"))
        origin = f"read from {C.SIBLING_DR16Q} (I3's VizieR VII/289 download)"
    else:
        ra, dec, z = [], [], []
        edges = np.linspace(0.0, 360.0, strips + 1)
        for k in range(strips):
            p = vizier_tsv(C.DR16Q_VIZIER, ["RAJ2000", "DEJ2000", "z"],
                           C.RAW_DIR / "dr16q" / f"dr16q_ra{k:02d}.tsv", manifest,
                           {"RAJ2000": f"{edges[k]:.6f}..{edges[k + 1]:.6f}"})
            for r in parse_vizier_tsv(p.read_text(encoding="utf-8")):
                ra.append(float(r["RAJ2000"])), dec.append(float(r["DEJ2000"])), z.append(float(r["z"]))
        ra, dec, z = (np.array(v) for v in (ra, dec, z))
        key = np.round(ra, 6) * 1e3 + np.round(dec, 6)   # strips share boundary values
        _, first = np.unique(key, return_index=True)
        ra, dec, z = ra[np.sort(first)], dec[np.sort(first)], z[np.sort(first)]
        origin = "VizieR VII/289/dr16q in RA strips"
    if z.size != DR16Q_N:
        raise RuntimeError(f"DR16Q: {z.size} quasars, expected {DR16Q_N}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest, z=z.astype(np.float64), ra=ra.astype(np.float32),
                        dec=dec.astype(np.float32))
    manifest[str(dest.relative_to(C.DATA_DIR))] = {
        "url": C.VIZIER_ASU + "?-source=" + C.DR16Q_VIZIER, "derived": origin,
        "bytes": dest.stat().st_size, "sha256": _sha256(dest), "retrieved_utc": _now()}
    print(f"  DR16Q: {z.size:,} quasar redshifts ({origin})", flush=True)
    return dest


# --- 8. SDSS Main Galaxy Sample ------------------------------------------------------------

#: D2's selection (07-cosmic-web/config.py): MGS target bits, r < 17.77 after dust,
#: 0.002 < z < 0.30, zWarning = 0, best spectrum and primary photometry.
MGS_SQL = """
SELECT z, ra, dec FROM SpecPhoto
WHERE (legacy_target1 & 448) != 0 AND class = 'GALAXY' AND sciencePrimary = 1 AND mode = 1
  AND zWarning = 0 AND z BETWEEN 0.002 AND 0.30 AND petroMag_r - extinction_r < 17.77
  AND ra >= {lo} AND ra < {hi}
"""
MGS_N = 656_795


def fetch_mgs(manifest: dict, strips: int = 12) -> Path:
    dest = C.DERIVED_DIR / "sdss_mgs_z.npz"
    if dest.exists():
        return dest
    sibling = C.ROOT / C.SIBLING_MGS
    if sibling.exists():
        d = np.load(sibling)
        z, ra, dec = d["z"], d["ra"], d["dec"]
        origin = f"read from {C.SIBLING_MGS} (D2's SkyServer download)"
    else:
        zs, ras, decs = [], [], []
        edges = np.linspace(0.0, 360.0, strips + 1)
        for k in range(strips):
            sql = " ".join(MGS_SQL.format(lo=edges[k], hi=edges[k + 1]).split())
            url = C.SKYSERVER_SQL + "?" + urllib.parse.urlencode({"cmd": sql, "format": "csv"})
            p = download(url, C.RAW_DIR / "sdss_mgs" / f"mgs_ra{k:02d}.csv", manifest)
            text = p.read_text()
            if not text.startswith("#Table1"):
                p.unlink()
                raise RuntimeError(f"SkyServer error: {text[:200]!r}")
            arr = np.loadtxt(io.StringIO(text), delimiter=",", skiprows=2, ndmin=2)
            zs.append(arr[:, 0]), ras.append(arr[:, 1]), decs.append(arr[:, 2])
        z, ra, dec = (np.concatenate(v) for v in (zs, ras, decs))
        origin = "SkyServer DR17 SpecPhoto, D2's selection"
    if z.size != MGS_N:
        raise RuntimeError(f"MGS: {z.size} galaxies, expected {MGS_N}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dest, z=z.astype(np.float64), ra=ra.astype(np.float32),
                        dec=dec.astype(np.float32))
    manifest[str(dest.relative_to(C.DATA_DIR))] = {
        "url": C.SKYSERVER_SQL, "derived": origin, "bytes": dest.stat().st_size,
        "sha256": _sha256(dest), "retrieved_utc": _now()}
    print(f"  SDSS MGS: {z.size:,} galaxy redshifts ({origin})", flush=True)
    return dest


# --- Everything ---------------------------------------------------------------------------


def fetch_all() -> dict:
    manifest = load_manifest()
    try:
        print("Davis & Lineweaver 2004 (arXiv source)")
        fetch_dl2004(manifest)
        print("DESI DR2 BAO")
        for name, path in C.BAO_FILES.items():
            link_or_download(f"raw/bao/{name}", C.BAO_REPO_RAW + path, manifest)
        print("Pantheon+")
        for name, path in C.SN_FILES.items():
            link_or_download(f"raw/sn/{name}", C.SN_REPO_RAW + path, manifest)
        print("GWOSC event list")
        download(C.GWOSC_URL, C.RAW_DIR / "gwosc" / "gwosc_GWTC.json", manifest)
        print("JADES DR4 redshifts")
        fetch_jades(manifest)
        print("SDSS DR16Q quasars")
        fetch_dr16q(manifest)
        print("Banados et al. 2016 quasars at z > 5.6")
        vizier_tsv("J/ApJS/227/11/highzqso", ["Name", "RAJ2000", "DEJ2000", "z", "n_z", "r_z"],
                   C.RAW_DIR / "highz_quasars" / "banados2016_highzqso.tsv", manifest)
        print("SDSS Main Galaxy Sample")
        fetch_mgs(manifest)
    finally:
        save_manifest(manifest)
    return manifest


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.split("\n")[0]).parse_args(argv)
    fetch_all()
    total = sum(v.get("bytes", 0) for v in load_manifest().values())
    print(f"done: {total / 1e6:.1f} MB cached in {C.DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
