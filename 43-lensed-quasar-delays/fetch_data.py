"""Download and cache every input this project uses.

    python fetch_data.py            # everything (~190 MB, a few minutes)
    python fetch_data.py --no-big   # skip the 84 MB PG 1115 distance chain

Files land in ``data/`` and are never re-downloaded once present, so after one
run the whole analysis works offline. Sources (see the README for citations):

* VizieR/CDS: COSMOGRAIL light curves of the lensed quasars.
* GitHub shsuyu/H0LiCOW-public: H0LiCOW's time-delay distance posteriors.
* GitHub TDCOSMO/TDCOSMO2025_public: TDCOSMO 2025's likelihood inputs and chains.
"""

from __future__ import annotations

import argparse
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

#: CDS answers browser-like User-Agents with an anti-bot page (HTTP 200), so
#: every request uses a plain tool User-Agent.
UA = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"


def download(url: str, dest: Path, *, retries: int = 4, min_bytes: int = 1) -> bool:
    """Fetch ``url`` to ``dest`` unless it is already cached. Returns True if fetched."""
    if dest.exists() and dest.stat().st_size >= min_bytes:
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r, open(tmp, "wb") as f:
                while chunk := r.read(1 << 20):
                    f.write(chunk)
            head = tmp.read_bytes()[:200].lower()
            if b"<!doctype html" in head or b"<html" in head:
                raise OSError(f"got an HTML page instead of data from {url}")
            tmp.replace(dest)
            print(f"  fetched {dest.relative_to(C.HERE)} ({dest.stat().st_size / 1e6:.2f} MB)")
            return True
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            print(f"  attempt {attempt}/{retries} failed for {url}: {exc}", file=sys.stderr)
            time.sleep(3 * attempt)
    raise SystemExit(f"could not download {url}")


def fetch_cds() -> None:
    print("COSMOGRAIL light curves (VizieR/CDS)")
    for cat, (folder, files) in C.CDS_CATALOGS.items():
        base = C.CDS_FTP + cat + "/"
        for name in files:
            download(base + name, C.DATA / folder / name)
    # COSMOGRAIL XIX: the 23 light-curve files are listed in stars.dat.
    folder = C.DATA / C.CDS_CATALOGS["J/A+A/640/A105"][0]
    for line in (folder / "stars.dat").read_text().splitlines():
        name = line[43:67].strip()  # bytes 44-67: FileName
        if name:
            download(C.CDS_FTP + "J/A+A/640/A105/" + name, folder / name)


def fetch_h0licow(skip_big: bool) -> None:
    print("H0LiCOW distance posteriors (GitHub shsuyu/H0LiCOW-public)")
    for rel in C.H0LICOW_FILES:
        if skip_big and "PG1115" in rel:
            print("  skipping", rel)
            continue
        download(C.H0LICOW_RAW + rel, C.DATA / "h0licow" / Path(rel).name)


def fetch_tdcosmo() -> None:
    print("TDCOSMO 2025 release (GitHub TDCOSMO/TDCOSMO2025_public)")
    for rel in C.TDCOSMO_FILES:
        download(C.TDCOSMO_RAW + rel, C.DATA / "tdcosmo2025" / rel)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-big", action="store_true", help="skip the 84 MB PG 1115 chain")
    args = ap.parse_args()
    C.DATA.mkdir(exist_ok=True)
    fetch_cds()
    fetch_h0licow(args.no_big)
    fetch_tdcosmo()
    total = sum(p.stat().st_size for p in C.DATA.rglob("*") if p.is_file())
    print(f"done: data/ holds {total / 1e6:.0f} MB")


if __name__ == "__main__":
    main()
