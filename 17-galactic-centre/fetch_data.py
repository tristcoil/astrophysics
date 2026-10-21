"""Download and cache everything the analysis needs.

About 18 MB in total, nearly all of it the UCLA group's posterior chains
(the S2 measurements themselves are ~140 kB):

1. **Do et al. (2019) data release** -- the Keck positions and radial
   velocities of S0-2 (= S2) from 1995 to 2019, and the authors' own MCMC
   chains. The group's page (astro.ucla.edu/~ghezgroup/gc/gr/) went offline
   around 2023, so we fetch the Internet Archive's copy of the original
   ``data.zip``, captured on 1 September 2019.

2. **Gillessen et al. (2017) tables from CDS VizieR** -- the VLT team's S2
   positions (1992-2016) and radial velocities, plus the orbital elements of
   40 S-stars.

Everything lands in ``data/raw/`` alongside ``data/manifest.json``, which
records each file's source URL, retrieval time and SHA-256 hash. Files already
present are never re-downloaded; ``--refresh`` forces it. The zip is unpacked
into ``data/raw/do2019/``.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone

from config import MANIFEST_PATH, RAW_DIR, SOURCES, Source

USER_AGENT = "astrophysics-project-atlas/H1 (hobby research; python-urllib)"
DO2019_DIR = RAW_DIR / "do2019"


def _download(source: Source, *, attempts: int = 4) -> bytes:
    """Fetch one URL, retrying with backoff; raise if every attempt fails."""
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                payload = response.read()
            if source.marker and source.marker not in payload[:200_000]:
                raise ValueError(f"{source.filename}: content does not look like the expected file")
            return payload
        except (urllib.error.URLError, TimeoutError, ValueError, ConnectionError) as err:
            last_error = err
            if attempt < attempts:
                time.sleep(3.0 * attempt)
    raise RuntimeError(f"could not download {source.url}: {last_error}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _unpack_do2019(zip_path) -> list[str]:
    """Unpack the UCLA zip, and the chains zip nested inside it."""
    DO2019_DIR.mkdir(parents=True, exist_ok=True)
    names = []
    with zipfile.ZipFile(zip_path) as outer:
        for member in outer.namelist():
            if member.startswith("__MACOSX") or member.endswith("/"):
                continue
            payload = outer.read(member)
            if member.endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(payload)) as inner:
                    for sub in inner.namelist():
                        if sub.startswith("__MACOSX") or sub.endswith("/"):
                            continue
                        (DO2019_DIR / sub).write_bytes(inner.read(sub))
                        names.append(sub)
            else:
                (DO2019_DIR / member).write_bytes(payload)
                names.append(member)
    return names


def fetch_all(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Make sure every source file is cached (and unpacked); return the manifest.

    Writes are atomic (download to a temporary name, then rename), so an
    interrupted download never leaves a truncated file that later looks cached.
    """
    manifest = _load_manifest()
    for source in SOURCES:
        path = RAW_DIR / source.filename
        if path.exists() and not refresh:
            if not quiet:
                print(f"  cached     {source.filename}")
        else:
            if not quiet:
                print(f"  fetching   {source.filename} ...", end="", flush=True)
            payload = _download(source)
            tmp = path.with_suffix(path.suffix + ".part")
            tmp.write_bytes(payload)
            tmp.replace(path)
            manifest["files"][source.filename] = {
                "url": source.url,
                "original_url": source.original_url or source.url,
                "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
                "what": source.what,
                "reference": source.reference,
            }
            MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))
            if not quiet:
                print(f" {len(payload) / 1024:.0f} kB")

        if source.filename == "do2019_data.zip" and (refresh or not (DO2019_DIR / "S0-2_rv.csv").exists()):
            names = _unpack_do2019(path)
            if not quiet:
                print(f"  unpacked   {', '.join(names)}")
    return manifest


def missing_files() -> list[str]:
    need = [RAW_DIR / s.filename for s in SOURCES]
    need += [DO2019_DIR / n for n in ("S0-2_astro.csv", "S0-2_rv.csv", "chains_Newton.txt",
                                      "chains_GR.txt", "chains_redshift.txt")]
    return [str(p.relative_to(RAW_DIR)) for p in need if not p.exists()]


def verify_cache() -> list[str]:
    """Names of cached files whose bytes no longer match the recorded hash."""
    manifest = _load_manifest()
    changed = []
    for name, record in manifest.get("files", {}).items():
        path = RAW_DIR / name
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            changed.append(name)
    return changed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download every file even if it is cached")
    args = parser.parse_args(argv)

    print("Fetching the S2 data of Do et al. 2019 (Keck) and Gillessen et al. 2017 (VLT)")
    try:
        fetch_all(refresh=args.refresh)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1

    changed = verify_cache()
    if changed:
        print(f"\nWARNING: cached files differ from their recorded hashes: {changed}")
    missing = missing_files()
    if missing:
        print(f"\nWARNING: still missing: {missing}")
        return 1
    print(f"\nDone. Everything is in {RAW_DIR.relative_to(RAW_DIR.parents[1])}/; "
          "the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
