"""Download and cache everything the analysis needs.

Four sources, all free and public:

1. **The Second CHIME/FRB Catalog** -- one CSV row per burst (sub-burst) with
   position, dispersion measure, width, scattering, spectrum, fluence and
   exposure. From CANFAR, DOI 10.11570/25.0066.
2. **CHIME dynamic spectra** for a handful of example bursts: the calibrated
   intensity in 16384 frequency channels x ~1 ms samples around each burst.
   About 10-20 MB each; skip with ``--no-dynamic-spectra``.
3. **Localized FRBs and their host galaxies** from the FRB community's open
   repository (github.com/FRBs/FRB), pinned to one commit. It also supplies an
   all-sky map of the NE2001 Galactic electron model.
4. **The ATNF Pulsar Catalogue**, as a control sample of radio sources that are
   known to lie inside the Milky Way.

Every file is recorded in ``data/manifest.json`` with its source URL, retrieval
time, size and SHA-256, so a result can always be traced to the exact bytes it
came from. Nothing already on disk is fetched again unless ``--refresh`` is
given, and once the cache is warm the analysis runs offline.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from datetime import datetime, timezone
from html import unescape
from pathlib import Path

import requests

from config import (
    ATNF_PARAMS,
    ATNF_PATH,
    ATNF_URL,
    CHIME_CATALOG_PATH,
    CHIME_CATALOG_URL,
    CHIME_DOI,
    CHIME_DYNSPEC_URL,
    DYNSPEC_DIR,
    EXAMPLE_BURSTS,
    FRB_REPO_COMMIT,
    FRB_REPO_DIR,
    FRB_REPO_FILES,
    MANIFEST_PATH,
)

USER_AGENT = "astrophysics-atlas-frb-project/1.0 (hobby research; python-requests)"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, dest: Path, *, params: dict | None = None, attempts: int = 3,
              timeout: float = 120.0) -> None:
    """Stream ``url`` to ``dest`` through a temporary file.

    Writing to ``.part`` and renaming at the end means an interrupted download
    never leaves a truncated file that a later run would mistake for the cache.
    CANFAR answers file requests with two redirects (to a storage node, then to a
    signed URL); ``requests`` follows them.
    """
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with requests.get(url, params=params, stream=True, timeout=timeout,
                              headers={"User-Agent": USER_AGENT}) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=1 << 16):
                        fh.write(chunk)
            tmp.replace(dest)
            return
        except requests.RequestException as exc:  # network hiccup: back off, retry
            last_error = exc
            if attempt < attempts:
                time.sleep(5 * attempt)
    tmp.unlink(missing_ok=True)
    raise RuntimeError(f"download failed after {attempts} attempts: {url}") from last_error


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(manifest: dict, path: Path, url: str, **extra) -> None:
    manifest["files"][str(path.relative_to(MANIFEST_PATH.parent))] = {
        "source": url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
        **extra,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))


def _fetch(manifest: dict, url: str, dest: Path, *, refresh: bool, quiet: bool, **extra) -> Path:
    if dest.exists() and not refresh:
        if not quiet:
            print(f"  cached   {dest.relative_to(MANIFEST_PATH.parent)}")
        return dest
    if not quiet:
        print(f"  fetching {dest.name} ...", end="", flush=True)
    _download(url, dest)
    _record(manifest, dest, url, **extra)
    if not quiet:
        print(f" {dest.stat().st_size / 1e6:.1f} MB")
    return dest


# ---------------------------------------------------------------------------
# The four sources
# ---------------------------------------------------------------------------


def fetch_chime_catalog(manifest: dict, *, refresh: bool = False, quiet: bool = False) -> Path:
    return _fetch(manifest, CHIME_CATALOG_URL, CHIME_CATALOG_PATH, refresh=refresh, quiet=quiet,
                  doi=CHIME_DOI, citation="CHIME/FRB Collaboration 2026, arXiv:2601.09399")


def dynspec_path(tns: str) -> Path:
    return DYNSPEC_DIR / f"{tns}_stokesi_dynamic_spectrum.h5"


def fetch_dynamic_spectra(manifest: dict, *, refresh: bool = False, quiet: bool = False) -> list[Path]:
    out = []
    for burst in EXAMPLE_BURSTS:
        out.append(_fetch(manifest, CHIME_DYNSPEC_URL.format(tns=burst.tns), dynspec_path(burst.tns),
                          refresh=refresh, quiet=quiet, doi=CHIME_DOI, role=burst.role))
    return out


def fetch_frb_repo(manifest: dict, *, refresh: bool = False, quiet: bool = False) -> list[Path]:
    out = []
    for name, url in FRB_REPO_FILES.items():
        out.append(_fetch(manifest, url, FRB_REPO_DIR / name, refresh=refresh, quiet=quiet,
                          commit=FRB_REPO_COMMIT, repository="https://github.com/FRBs/FRB"))
    return out


#: Columns of the ATNF table in the order our query requests them.
ATNF_COLUMNS = ["jname", "ra_deg", "dec_deg", "gl", "gb", "dm", "assoc", "type"]


def parse_atnf_html(html: str) -> tuple[str, list[list[str]]]:
    """Pull the catalogue version and the data rows out of the ATNF web page.

    The page returns a fixed-format table inside ``<pre>``. Rows start with a
    running index; missing values are ``*``. No field contains spaces (the
    association and type columns use commas and brackets), so splitting on
    whitespace is safe -- and checked: every row must have the expected width.
    """
    version_match = re.search(r"Catalogue Version:\s*([0-9.]+)", html)
    if not version_match:
        raise RuntimeError("ATNF response has no catalogue version -- query format changed?")
    pre = re.search(r"<pre>(.*?)</pre>", html, flags=re.S)
    if not pre:
        raise RuntimeError("ATNF response has no <pre> table")
    rows = []
    for line in unescape(pre.group(1)).splitlines():
        parts = line.split()
        if not parts or not parts[0].isdigit():
            continue  # rule lines, header lines, blank separators
        fields = parts[1:]
        if len(fields) != len(ATNF_COLUMNS):
            raise RuntimeError(f"unexpected ATNF row width {len(fields)}: {line!r}")
        rows.append(fields)
    return version_match.group(1), rows


def fetch_atnf(manifest: dict, *, refresh: bool = False, quiet: bool = False) -> Path:
    if ATNF_PATH.exists() and not refresh:
        if not quiet:
            print(f"  cached   {ATNF_PATH.relative_to(MANIFEST_PATH.parent)}")
        return ATNF_PATH
    if not quiet:
        print("  querying the ATNF Pulsar Catalogue ...", end="", flush=True)
    r = requests.get(ATNF_URL, params=ATNF_PARAMS, timeout=120, headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    version, rows = parse_atnf_html(r.text)
    tmp = ATNF_PATH.with_suffix(".part")
    with open(tmp, "w") as fh:
        fh.write(f"# ATNF Pulsar Catalogue version {version}; retrieved "
                 f"{datetime.now(timezone.utc).date()}; missing values are '*'\n")
        fh.write(",".join(ATNF_COLUMNS) + "\n")
        for row in rows:
            fh.write(",".join(f'"{v}"' if "," in v else v for v in row) + "\n")
    tmp.replace(ATNF_PATH)
    _record(manifest, ATNF_PATH, r.url, catalogue_version=version, n_pulsars=len(rows),
            citation="Manchester, Hobbs, Teoh & Hobbs 2005, AJ 129, 1993")
    if not quiet:
        print(f" version {version}, {len(rows)} pulsars")
    return ATNF_PATH


def fetch_all(*, dynamic_spectra: bool = True, refresh: bool = False, quiet: bool = False) -> dict:
    manifest = _load_manifest()
    if not quiet:
        print("CHIME/FRB Catalog 2 (CANFAR)")
    fetch_chime_catalog(manifest, refresh=refresh, quiet=quiet)
    if dynamic_spectra:
        if not quiet:
            print("CHIME dynamic spectra for the example bursts")
        fetch_dynamic_spectra(manifest, refresh=refresh, quiet=quiet)
    if not quiet:
        print(f"Localized FRBs + NE2001 map (FRBs/FRB @ {FRB_REPO_COMMIT[:10]})")
    fetch_frb_repo(manifest, refresh=refresh, quiet=quiet)
    if not quiet:
        print("ATNF Pulsar Catalogue")
    fetch_atnf(manifest, refresh=refresh, quiet=quiet)
    return _load_manifest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-dynamic-spectra", action="store_true",
                        help="skip the ~70 MB of example-burst dynamic spectra")
    parser.add_argument("--refresh", action="store_true", help="re-download even if cached")
    args = parser.parse_args(argv)
    try:
        fetch_all(dynamic_spectra=not args.no_dynamic_spectra, refresh=args.refresh)
    except RuntimeError as exc:
        print(f"\nERROR: {exc}", file=sys.stderr)
        return 1
    print("\nDone. Cache is warm; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
