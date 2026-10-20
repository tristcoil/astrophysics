"""Download and cache everything the analysis needs (about 9 MB in total).

What this project reads, all public:

1. **HEPData** (hepdata.net): the measured yields from ALICE (Pb-Pb at
   2.76 TeV: pi, K, p, K0S, Lambda, Xi, Omega, K*, phi, nuclei, hypertriton,
   direct photons) and STAR (Au+Au at 7.7-39 GeV), as the full YAML
   submission of each paper.
2. **arXiv e-print sources** of eleven papers: tables that are not on HEPData
   (deuteron, 3He and 4He yields), lattice-QCD tables and numbers, and the
   published thermal fit of Andronic et al. 2018, whose vector figures we read
   back point by point.
3. **Thermal-FIST**'s hadron lists and decay tables (PDG 2025 and PDG 2014),
   pinned to one commit, plus two of its transcriptions used only as
   cross-checks.

Each file's URL, size, SHA-256 and retrieval time goes into
``data/manifest.json``; nothing cached is fetched twice, so after one run the
whole project works offline.

    python fetch_data.py
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402


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
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(cfg.MANIFEST_PATH)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _get(url: str, retries: int = 5, timeout: int = 120) -> bytes:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as err:
            # 429: HEPData allows 60 requests a minute
            if err.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                raise
            print(f"    HTTP {err.code}, retrying")
            time.sleep(5 + 10 * attempt)
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if attempt == retries - 1:
                raise
            print(f"    retrying after {err!r}")
            time.sleep(2 + 3 * attempt)
    raise RuntimeError("unreachable")


def download(url: str, dest: Path, manifest: dict) -> Path:
    """Fetch one file unless it is cached; record it in the manifest."""
    key = str(dest.relative_to(cfg.DATA_DIR))
    if dest.exists():
        if key not in manifest:
            manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                             "retrieved_utc": datetime.fromtimestamp(
                                 dest.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")}
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    blob = _get(url)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(blob)
    tmp.replace(dest)
    manifest[key] = {"url": url, "bytes": len(blob), "sha256": _sha256(dest),
                     "retrieved_utc": _now()}
    print(f"  fetched  {key}  ({len(blob) / 1e3:.0f} kB)")
    return dest


def _extract(archive: Path, into: Path) -> None:
    """Unpack a downloaded archive (tar.gz, or a single gzipped file) safely."""
    marker = into / ".extracted"
    if marker.exists() and marker.read_text() == _sha256(archive):
        return
    into.mkdir(parents=True, exist_ok=True)
    raw = archive.read_bytes()
    try:
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:*") as tar:
            tar.extractall(into, filter="data")
    except tarfile.ReadError:
        # arXiv serves single-file submissions as a gzipped TeX file
        (into / "main.tex").write_bytes(gzip.decompress(raw))
    marker.write_text(_sha256(archive))


def hepdata_dir(key: str) -> Path:
    """Folder holding the unpacked YAML tables of one HEPData record."""
    ins, version, _ = cfg.HEPDATA[key]
    base = cfg.HEPDATA_DIR / key
    sub = base / f"HEPData-ins{ins}-v{version}-yaml"
    return sub if sub.exists() else base


def arxiv_dir(arxiv_id: str) -> Path:
    return cfg.ARXIV_DIR / arxiv_id / "source"


def fetch_all(verbose: bool = True) -> dict:
    manifest = _load_manifest()
    if verbose:
        print("Fetching data (cached files are reused)")
    for key, (ins, version, _desc) in cfg.HEPDATA.items():
        url = f"https://www.hepdata.net/download/submission/ins{ins}/{version}/yaml"
        arch = download(url, cfg.HEPDATA_DIR / key / "submission.tar.gz", manifest)
        _extract(arch, cfg.HEPDATA_DIR / key)
    for arxiv_id in cfg.ARXIV:
        arch = download(f"https://arxiv.org/e-print/{arxiv_id}",
                        cfg.ARXIV_DIR / arxiv_id / "source.tar.gz", manifest)
        _extract(arch, arxiv_dir(arxiv_id))
    for rel in cfg.THERMAL_FIST_FILES:
        download(cfg.THERMAL_FIST_RAW + rel, cfg.LISTS_DIR / rel, manifest)
    _save_manifest(manifest)
    total = sum(v["bytes"] for v in manifest.values())
    if verbose:
        print(f"  {len(manifest)} files, {total / 1e6:.1f} MB in data/ (manifest.json)")
    return manifest


if __name__ == "__main__":
    fetch_all()
