"""Download and cache everything the project needs.

About 400 MB in total, all public:

1. **The gallery**: calibrated SDSS g, r, i images of twelve galaxies, cut out
   on a grid we choose by the Legacy Surveys viewer, plus the viewer's own
   colour JPEG of exactly the same pixels (to check we read the data the way
   the survey does).
2. **Catalogues of the stars in each field**: SDSS PSF photometry and Gaia DR3
   (with synthetic SDSS magnitudes from Gaia's spectra), to check the image
   calibration and to mask foreground stars; and RC3 galaxy sizes.
3. **Colour-science reference data**: the CIE colour-matching functions (how
   the human eye responds to each wavelength), filter transmission curves, the
   measured spectra of the Sun, Vega and a library of stars, and the Brown et
   al. (2014) atlas of integrated galaxy spectra.
4. **M51 across the spectrum**: ultraviolet (GALEX), H-alpha, infrared
   (Spitzer), radio (VLA 21 cm), and Hubble's mosaic.

Everything lands in ``data/`` with ``manifest.json`` recording each file's
source URL, retrieval time and SHA-256 hash. Files already present are never
re-downloaded; ``--refresh`` forces it. Requests go one at a time, with pauses,
and back off when a server says it is busy.

    python fetch_data.py                 # everything
    python fetch_data.py --only gallery  # just one group
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import os
import http.client
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

from config import DATA_DIR, MANIFEST_PATH, SOURCES, Source

USER_AGENT = "astrophysics-project-atlas/G1 (hobby research; python-urllib)"

GROUPS = {
    "gallery": ("cutouts/",),
    "catalogs": ("catalogs/",),
    "reference": ("reference/",),
    "multiwavelength": ("multiwavelength/",),
}


class BadPayload(ValueError):
    """The server answered, but not with the file we asked for."""


def _check_payload(source: Source, head: bytes) -> None:
    """Cheap check that the bytes are the expected kind of file -- not, say, an
    HTML error page or a TAP error document served with status 200."""
    text = head[:4000].decode("latin-1", errors="replace")
    if source.kind == "fits" and not head.startswith(b"SIMPLE"):
        raise BadPayload(f"{source.filename}: not a FITS file ({text[:80]!r})")
    if source.kind == "jpeg" and not head.startswith(b"\xff\xd8"):
        raise BadPayload(f"{source.filename}: not a JPEG ({text[:80]!r})")
    if source.kind == "xml" and 'status="ok"' not in text:
        raise BadPayload(f"{source.filename}: service did not report success ({text[:200]!r})")
    if source.kind in ("csv", "tsv", "text"):
        stripped = text.lstrip()
        if stripped.startswith("<") or "QUERY_STATUS" in text or "ERROR" in text[:200]:
            raise BadPayload(f"{source.filename}: service returned an error ({stripped[:200]!r})")
        if not stripped:
            raise BadPayload(f"{source.filename}: empty response")


def _download(source: Source, dest, *, attempts: int = 8) -> tuple[int, str]:
    """Stream one URL to ``dest``; return (bytes, sha256). Retries with backoff.

    HTTP 429 (too many requests) and 5xx answers are retried after a pause that
    grows each time; other 4xx answers mean the request itself is wrong and
    fail at once. Dropped connections, failed DNS lookups and transfers cut off
    half-way are retried quickly.
    """
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            sha = hashlib.sha256()
            n = 0
            with urllib.request.urlopen(request, timeout=300) as response, open(dest, "wb") as out:
                first = True
                while chunk := response.read(1 << 20):
                    if first:
                        _check_payload(source, chunk)
                        first = False
                    out.write(chunk)
                    sha.update(chunk)
                    n += len(chunk)
            if n == 0:
                raise BadPayload(f"{source.filename}: empty response")
            return n, sha.hexdigest()
        except urllib.error.HTTPError as err:
            last_error = err
            if err.code != 429 and err.code < 500:
                break
            wait = float(err.headers.get("Retry-After") or 15 * attempt)
        except (urllib.error.URLError, http.client.HTTPException, OSError, BadPayload) as err:
            last_error = err
            wait = 3.0 * attempt
        if attempt < attempts:
            time.sleep(wait)
    raise RuntimeError(f"could not download {source.filename} from {source.url[:120]}...: {last_error}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(filename: str, entry: dict) -> None:
    """Add one entry to the manifest: lock, re-read, update, write.

    Re-reading under a lock means two fetch processes running at once (say,
    one per group of files) never overwrite each other's entries.
    """
    lock = MANIFEST_PATH.with_suffix(".lock")
    with open(lock, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        manifest = _load_manifest()
        manifest["files"][filename] = entry
        tmp = MANIFEST_PATH.with_suffix(f".tmp{os.getpid()}")
        tmp.write_text(json.dumps(manifest, indent=2))
        tmp.replace(MANIFEST_PATH)
        fcntl.flock(fh, fcntl.LOCK_UN)


def _sha256(path) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            sha.update(chunk)
    return sha.hexdigest()


def repair_manifest() -> list[str]:
    """Record any cached file that is missing from the manifest (its hash now,
    its file time as the retrieval time). Returns the names added."""
    manifest = _load_manifest()
    added = []
    for source in SOURCES:
        path = DATA_DIR / source.filename
        if path.exists() and source.filename not in manifest["files"]:
            _record(source.filename, {
                "url": source.url,
                "retrieved_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
                "what": source.what,
                "reference": source.reference,
                "used_as": "analysis input" if source.input else "comparison / validation only",
                "note": "recorded after download (retrieval time = file time)",
            })
            added.append(source.filename)
    return added


def selected(only: str | None = None) -> list[Source]:
    if not only:
        return list(SOURCES)
    prefixes = GROUPS[only]
    return [s for s in SOURCES if s.filename.startswith(prefixes)]


def fetch_all(*, refresh: bool = False, quiet: bool = False, only: str | None = None) -> dict:
    """Make sure every source file is cached; return the manifest.

    Writes are atomic (download to a temporary name unique to this process,
    then rename), so an interrupted download never leaves a truncated file that
    later looks cached, and two fetch processes can safely run side by side.
    """
    failures = []
    for source in selected(only):
        path = DATA_DIR / source.filename
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not refresh:
            continue
        if source.pause:
            time.sleep(source.pause)
        if not quiet:
            print(f"  fetching   {source.filename} ...", end="", flush=True)
        tmp = path.with_name(f"{path.name}.part{os.getpid()}")
        try:
            nbytes, sha = _download(source, tmp)
        except RuntimeError as err:
            tmp.unlink(missing_ok=True)
            failures.append(str(err))
            if not quiet:
                print(" FAILED")
            continue
        tmp.replace(path)
        _record(source.filename, {
            "url": source.url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes": nbytes,
            "sha256": sha,
            "what": source.what,
            "reference": source.reference,
            "used_as": "analysis input" if source.input else "comparison / validation only",
        })
        if not quiet:
            print(f" {nbytes / 1e6:.2f} MB" if nbytes > 1e6 else f" {nbytes / 1e3:.0f} kB")
    if failures:
        raise RuntimeError("some downloads failed:\n  " + "\n  ".join(failures))
    return _load_manifest()


def missing_files(only: str | None = None) -> list[str]:
    return [s.filename for s in selected(only) if not (DATA_DIR / s.filename).exists()]


def verify_cache() -> list[str]:
    """Names of cached files whose bytes no longer match the recorded hash."""
    manifest = _load_manifest()
    changed = []
    for name, record in manifest.get("files", {}).items():
        path = DATA_DIR / name
        if path.exists() and _sha256(path) != record["sha256"]:
            changed.append(name)
    return changed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download every file even if it is cached")
    parser.add_argument("--only", choices=sorted(GROUPS), help="fetch one group of files")
    parser.add_argument("--repair-manifest", action="store_true",
                        help="record cached files missing from the manifest, then exit")
    args = parser.parse_args(argv)
    if args.repair_manifest:
        added = repair_manifest()
        print(f"recorded {len(added)} file(s) missing from the manifest")
        return 0

    todo = missing_files(args.only) if not args.refresh else [s.filename for s in selected(args.only)]
    print(f"{len(selected(args.only))} files in scope, {len(todo)} to fetch")
    try:
        fetch_all(refresh=args.refresh, only=args.only)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1

    changed = verify_cache()
    if changed:
        print(f"\nWARNING: cached files differ from their recorded hashes: {changed}")
    total = sum((DATA_DIR / s.filename).stat().st_size for s in SOURCES
                if (DATA_DIR / s.filename).exists())
    print(f"\nDone. {total / 1e6:.0f} MB cached in {DATA_DIR.name}/; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
