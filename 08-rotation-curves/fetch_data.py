"""Download and cache everything the analysis needs.

About 450 kB in total, all of it public:

1. **The SPARC database** (Lelli, McGaugh & Schombert 2016), hosted at Case
   Western Reserve University: 175 disc galaxies, each with a rotation curve
   measured from the Doppler shift of hydrogen gas, and the rotation the
   *visible* matter alone would produce, computed from Spitzer 3.6-micron
   images of the stars and 21-cm radio maps of the gas.

2. **Published results built on SPARC** -- the radial acceleration relation,
   the baryonic Tully-Fisher relation, and maximum-disc mass-to-light ratios.
   The analysis never reads these; they are used only at the end, to grade our
   own numbers against the literature.

3. **NASA's Planetary Fact Sheet**, for the Solar System's orbital speeds: the
   textbook example of what a Keplerian rotation curve looks like.

Everything lands in ``data/`` alongside ``manifest.json``, which records each
file's source URL, retrieval time and SHA-256 hash, so a result can always be
traced to the exact bytes it came from. Files already present are never
re-downloaded; ``--refresh`` forces it.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

from config import DATA_DIR, MANIFEST_PATH, SOURCES, Source

USER_AGENT = "astrophysics-project-atlas/D4 (hobby research; python-urllib)"

#: A cheap check that what arrived is the file we asked for and not, say, an
#: HTML error page served with status 200. Each file must contain its marker.
CONTENT_MARKERS = {
    ".mrt": ("Title:", "# CAUTION"),
    ".html": ("Orbital Velocity",),
}


def _looks_right(source: Source, payload: bytes) -> bool:
    suffix = source.filename[source.filename.rfind("."):]
    text = payload[:20000].decode("latin-1")
    return any(marker in text for marker in CONTENT_MARKERS.get(suffix, ("",)))


def _download(source: Source, *, attempts: int = 3) -> bytes:
    """Fetch one URL, retrying with backoff; raise if every attempt fails."""
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read()
            if not _looks_right(source, payload):
                raise ValueError(f"{source.filename}: content does not look like the expected file")
            return payload
        except (urllib.error.URLError, TimeoutError, ValueError) as err:
            last_error = err
            if attempt < attempts:
                time.sleep(2.0 * attempt)
    raise RuntimeError(f"could not download {source.url}: {last_error}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def fetch_all(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Make sure every source file is cached; return the manifest.

    Writes are atomic (download to a temporary name, then rename), so an
    interrupted download never leaves a truncated file that later looks cached.
    """
    manifest = _load_manifest()
    for source in SOURCES:
        path = DATA_DIR / source.filename
        if path.exists() and not refresh:
            if not quiet:
                print(f"  cached     {source.filename}")
            continue

        if not quiet:
            print(f"  fetching   {source.filename} ...", end="", flush=True)
        payload = _download(source)
        tmp = path.with_suffix(path.suffix + ".part")
        tmp.write_bytes(payload)
        tmp.replace(path)

        manifest["files"][source.filename] = {
            "url": source.url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "what": source.what,
            "reference": source.reference,
            "used_as": "analysis input" if source.input else "published comparison only",
        }
        MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))
        if not quiet:
            print(f" {len(payload) / 1024:.0f} kB")
    return manifest


def missing_files() -> list[str]:
    return [s.filename for s in SOURCES if not (DATA_DIR / s.filename).exists()]


def verify_cache() -> list[str]:
    """Names of cached files whose bytes no longer match the recorded hash."""
    manifest = _load_manifest()
    changed = []
    for name, record in manifest.get("files", {}).items():
        path = DATA_DIR / name
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            changed.append(name)
    return changed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download every file even if it is cached")
    args = parser.parse_args(argv)

    print("Fetching SPARC tables and the NASA planetary fact sheet")
    try:
        fetch_all(refresh=args.refresh)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1

    changed = verify_cache()
    if changed:
        print(f"\nWARNING: cached files differ from their recorded hashes: {changed}")
    print(f"\nDone. {len(SOURCES)} files in {DATA_DIR.name}/; analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
