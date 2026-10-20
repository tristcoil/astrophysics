"""Download and cache everything the project needs (about 13 MB compressed).

Three public catalogues from the CDS archive in Strasbourg:

1. **Cosmicflows-4** (Tully et al. 2023): distances to 55,877 galaxies in
   38,053 groups, measured without using their redshifts.
2. **The 2MASS Redshift Survey** (Huchra et al. 2012): 44,599 galaxies brighter
   than Ks = 11.75 over 91% of the sky -- the map of where the mass is.
3. **Tully's (2015) 2MASS groups**: which 2MRS galaxies belong together, so that
   the random motions inside clusters can be averaged away.

Each file lands in ``data/`` (gzipped tables are unpacked) and
``data/manifest.json`` records its source URL, retrieval time and the SHA-256
of the bytes as served. Files already present are never re-downloaded;
``--refresh`` forces it.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import http.client
import json
import os
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

from config import DATA_DIR, MANIFEST_PATH, SOURCES, USER_AGENT, Source


class BadPayload(ValueError):
    """The server answered, but not with the file we asked for."""


def _check_payload(source: Source, head: bytes) -> None:
    if source.gzipped:
        if not head.startswith(b"\x1f\x8b"):
            raise BadPayload(f"{source.url}: not gzip ({head[:120]!r})")
        return
    if source.local.endswith(".pdf"):
        if not head.startswith(b"%PDF"):
            raise BadPayload(f"{source.url}: not a PDF")
        return
    text = head[:2000].decode("latin-1", errors="replace").lstrip()
    if not text or text.startswith("<"):
        raise BadPayload(f"{source.url}: empty or HTML response ({text[:200]!r})")


def _download(source: Source, *, attempts: int = 6) -> bytes:
    """Fetch one URL into memory, retrying transient failures with backoff.

    HTTP 429 and 5xx answers and dropped connections are retried; other 4xx
    answers mean the request itself is wrong and fail at once.
    """
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                payload = response.read()
            if not payload:
                raise BadPayload(f"{source.url}: empty response")
            _check_payload(source, payload)
            return payload
        except urllib.error.HTTPError as err:
            last = err
            if err.code != 429 and err.code < 500:
                break
            wait = float(err.headers.get("Retry-After") or 10 * attempt)
        except (urllib.error.URLError, http.client.HTTPException, OSError, BadPayload) as err:
            last = err
            if isinstance(err, BadPayload) and attempt >= 2:
                break
            wait = 3.0 * attempt
        if attempt < attempts:
            time.sleep(wait)
    raise RuntimeError(f"could not download {source.url}: {last}")


def _record(entry_name: str, entry: dict) -> None:
    manifest = json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {"files": {}}
    manifest["files"][entry_name] = entry
    tmp = MANIFEST_PATH.with_suffix(f".tmp{os.getpid()}")
    tmp.write_text(json.dumps(manifest, indent=1))
    tmp.replace(MANIFEST_PATH)


def missing() -> list[Source]:
    return [s for s in SOURCES if not (DATA_DIR / s.local).exists()]


def fetch_all(*, refresh: bool = False, quiet: bool = False) -> list[str]:
    """Make sure every source is cached; return the list of failures.

    Writes are atomic (temporary file, then rename), so an interrupted run never
    leaves a truncated file that later looks cached.
    """
    failures = []
    for source in SOURCES:
        path = DATA_DIR / source.local
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not refresh:
            continue
        if not quiet:
            print(f"  fetching {source.local:<28}", end="", flush=True)
        try:
            raw = _download(source)
        except RuntimeError as err:
            failures.append(str(err))
            if not quiet:
                print(" FAILED")
            continue
        data = gzip.decompress(raw) if source.gzipped else raw
        tmp = path.with_name(f"{path.name}.part{os.getpid()}")
        tmp.write_bytes(data)
        tmp.replace(path)
        _record(source.local, {
            "url": source.url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes_served": len(raw),
            "sha256_served": hashlib.sha256(raw).hexdigest(),
            "bytes_unpacked": len(data),
            "what": source.what,
            "reference": source.reference,
        })
        if not quiet:
            print(f" {len(raw) / 1e6:6.2f} MB")
        time.sleep(0.5)
    return failures


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download even if cached")
    args = parser.parse_args(argv)
    failures = fetch_all(refresh=args.refresh)
    for f in failures:
        print("ERROR:", f, file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
