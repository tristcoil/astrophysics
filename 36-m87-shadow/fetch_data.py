"""Download and cache the Event Horizon Telescope data releases this project uses.

About 11 MB in total, all from the EHT Collaboration's own GitHub data releases,
pinned to exact commits so the bytes can never silently change:

1. **M87, 2017** (data product 2019-D01-01, 1.6 MB): the calibrated visibilities
   behind the first black hole image. Four nights (April 5, 6, 10, 11) in two
   2-GHz bands at 227 and 229 GHz, from the EHT-HOPS pipeline, averaged to 10 s.
2. **M87, 2018** (2024-D01-01, 2.7 MB): the same source a year later, four bands,
   two independent calibration pipelines.
3. **Sgr A\\*, 2017** (2022-D02-01, 7.1 MB): the Milky Way's black hole on April 6
   and 7, including copies normalised by the light curve (Sgr A* flickers within
   minutes; the normalised data are what the EHT itself modelled).
4. **The EHT's imaging script** (2019-D01-02): the eht-imaging pipeline that made
   the famous 2019 picture, re-run here unchanged except for one line.

Everything lands in ``data/raw/``; tarballs are unpacked next to themselves.
``data/manifest.json`` records each file's URL, retrieval time and SHA-256.
Files already present are never re-downloaded (``--refresh`` forces it).

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

from config import MANIFEST_PATH, RAW_DIR, SOURCES, Source, ensure_dirs

USER_AGENT = "astrophysics-project-atlas/I5 (hobby research; python-urllib)"


def _download(source: Source, *, attempts: int = 4) -> bytes:
    """Fetch one URL, retrying with backoff; raise if every attempt fails."""
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last_error = err
            if attempt < attempts:
                time.sleep(3.0 * attempt)
    raise RuntimeError(f"could not download {source.url}: {last_error}")


def _unpack(path, target_dir) -> list[str]:
    """Extract a tarball safely (no absolute paths, no '..') and list what came out."""
    target_dir.mkdir(parents=True, exist_ok=True)
    names = []
    with tarfile.open(path, "r:gz") as tar:
        members = [m for m in tar.getmembers() if m.isfile()]
        for m in members:
            if m.name.startswith("/") or ".." in m.name.split("/"):
                raise ValueError(f"refusing to extract suspicious path {m.name}")
        tar.extractall(target_dir, members=members, filter="data")
        names = [m.name for m in members]
    return names


def fetch(refresh: bool = False) -> dict:
    ensure_dirs()
    manifest = json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {"files": {}}
    for source in SOURCES:
        path = RAW_DIR / source.filename
        if path.exists() and not refresh and source.filename in manifest["files"]:
            print(f"  cached   {source.filename}")
        else:
            print(f"  fetching {source.url}")
            payload = _download(source)
            path.write_bytes(payload)
            manifest["files"][source.filename] = {
                "key": source.key,
                "url": source.url,
                "description": source.description,
                "reference": source.reference,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
            print(f"           {len(payload) / 1e6:.2f} MB")
        if source.unpack_to:
            target = RAW_DIR / source.unpack_to
            if refresh or not target.exists():
                names = _unpack(path, target)
                manifest["files"][source.filename]["unpacked"] = sorted(names)
                print(f"           unpacked {len(names)} files into data/raw/{source.unpack_to}/")
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=1))
    return manifest


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--refresh", action="store_true", help="download again even if cached")
    args = parser.parse_args(argv)
    fetch(refresh=args.refresh)
    print(f"manifest: {MANIFEST_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
