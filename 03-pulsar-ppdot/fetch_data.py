"""Download and cache the ATNF Pulsar Catalogue.

The ATNF Pulsar Catalogue (Manchester, Hobbs, Teoh & Hobbs 2005, AJ 129, 1993) is
maintained by CSIRO's Australia Telescope National Facility. It compiles every
published pulsar -- radio pulsars, but also X-ray and gamma-ray pulsars and
magnetars -- with its spin period, period derivative, sky position, dispersion
measure, binary orbit and literature references, one plain-text record per pulsar.

The whole catalogue is a single ~1.4 MB tarball. We download one pinned release,
extract its database file ``psrcat.db``, and record where and when it came from in
``data/manifest.json``. Nothing is re-downloaded once it is on disk, and the
analysis reads the extracted file directly, so it runs offline from then on.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

import requests

from config import (
    CATALOGUE_DB,
    CATALOGUE_URL,
    CATALOGUE_VERSION,
    MANIFEST_PATH,
    TARBALL_PATH,
)


def _download(url: str, dest: Path) -> None:
    """Stream ``url`` to ``dest`` via a temporary file, so an interrupted download
    never leaves a truncated file that later runs would mistake for the cache."""
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
    tmp.replace(dest)


def _extract_db(tarball: Path, dest: Path) -> None:
    """Pull ``psrcat.db`` out of the release tarball (it sits in a subdirectory)."""
    with tarfile.open(tarball, mode="r:gz") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith("/psrcat.db")]
        if len(members) != 1:
            raise RuntimeError(f"expected one psrcat.db in {tarball.name}, found {len(members)}")
        src = tar.extractfile(members[0])
        dest.write_bytes(src.read())


def catalogue_version(db_path: Path) -> str | None:
    """The release number, as stamped in the database's own header comment."""
    with open(db_path, encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#") and "CATALOGUE" in line:
                return line.split()[-1].strip()
            if not line.startswith("#"):
                break
    return None


def count_records(db_path: Path) -> int:
    """Number of pulsars: every record carries exactly one ``PSRJ`` line."""
    with open(db_path, encoding="latin-1") as fh:
        return sum(1 for line in fh if line.startswith("PSRJ "))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def fetch_catalogue(*, refresh: bool = False, quiet: bool = False) -> Path:
    """Make sure the pinned catalogue release is on disk; return the database path."""
    if CATALOGUE_DB.exists() and MANIFEST_PATH.exists() and not refresh:
        if not quiet:
            print(f"catalogue v{CATALOGUE_VERSION} already cached: {CATALOGUE_DB.name}")
        return CATALOGUE_DB

    if refresh or not TARBALL_PATH.exists():
        if not quiet:
            print(f"downloading {CATALOGUE_URL}")
        _download(CATALOGUE_URL, TARBALL_PATH)

    _extract_db(TARBALL_PATH, CATALOGUE_DB)

    stamped = catalogue_version(CATALOGUE_DB)
    if stamped != CATALOGUE_VERSION:
        raise RuntimeError(
            f"asked for catalogue v{CATALOGUE_VERSION} but the file says v{stamped}"
        )

    manifest = {
        "source": "ATNF Pulsar Catalogue (CSIRO)",
        "citation": "Manchester, R. N., Hobbs, G. B., Teoh, A. & Hobbs, M. 2005, AJ, 129, 1993",
        "url": CATALOGUE_URL,
        "version": stamped,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tarball_bytes": TARBALL_PATH.stat().st_size,
        "tarball_sha256": _sha256(TARBALL_PATH),
        "n_pulsars": count_records(CATALOGUE_DB),
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    if not quiet:
        print(f"cached catalogue v{stamped}: {manifest['n_pulsars']:,} pulsars "
              f"({manifest['tarball_bytes'] / 1e6:.2f} MB tarball)")
    return CATALOGUE_DB


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="download again even if the release is already cached")
    args = parser.parse_args(argv)
    fetch_catalogue(refresh=args.refresh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
