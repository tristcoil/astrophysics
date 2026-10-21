"""Download and cache the three data sets the project uses.

1. **Meteorite Landings** -- NASA's copy of the Meteoritical Society's catalogue:
   45,716 meteorites, each with its name, class, recovered mass, whether it was
   seen to fall or found later, a year and (usually) coordinates. It is a snapshot
   of the Meteoritical Bulletin Database taken in 2013. The live database
   (lpi.usra.edu/meteor) is larger today but sits behind a browser check that
   scripts cannot pass, so this snapshot is the freely downloadable version.
2. **Natural Earth admin-0 countries** at 1:50m and 1:110m (public domain) --
   the base map, and the polygons that say which country a meteorite landed in.
3. **World Bank population** (series SP.POP.TOTL, CC BY 4.0) for 1960-2009, to
   test whether witnessed falls follow people.

Everything is cached in ``data/`` and never fetched twice; ``data/manifest.json``
records where each file came from, when, and its SHA-256. Once the cache is warm
the whole project runs offline.

    python fetch_data.py            # populate the cache
    python fetch_data.py --refresh  # download everything again
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

from config import (
    COUNTRIES_50M,
    COUNTRIES_50M_URL,
    COUNTRIES_110M,
    COUNTRIES_110M_URL,
    MANIFEST_PATH,
    METEORITE_CSV,
    METEORITE_LANDING_PAGE,
    METEORITE_ROWS,
    METEORITE_SHA256,
    METEORITE_URL,
    NATURAL_EARTH_TAG,
    POPULATION_JSON,
    POPULATION_URL,
)

USER_AGENT = "astrophysics-atlas-G2/1.0 (hobby project; caches every download)"

SOURCES = {
    "meteorites": dict(
        url=METEORITE_URL, path=METEORITE_CSV,
        description="Meteorite Landings (Meteoritical Society catalogue, 2013 snapshot)",
        publisher=f"NASA Open Data Portal, {METEORITE_LANDING_PAGE}",
        citation="The Meteoritical Society, Meteoritical Bulletin Database "
                 "(https://www.lpi.usra.edu/meteor/)",
    ),
    "countries_50m": dict(
        url=COUNTRIES_50M_URL, path=COUNTRIES_50M,
        description=f"Natural Earth {NATURAL_EARTH_TAG} admin-0 countries, 1:50m",
        publisher="Natural Earth (public domain), naturalearthdata.com",
        citation="Made with Natural Earth. Free vector and raster map data @ naturalearthdata.com",
    ),
    "countries_110m": dict(
        url=COUNTRIES_110M_URL, path=COUNTRIES_110M,
        description=f"Natural Earth {NATURAL_EARTH_TAG} admin-0 countries, 1:110m",
        publisher="Natural Earth (public domain), naturalearthdata.com",
        citation="Made with Natural Earth. Free vector and raster map data @ naturalearthdata.com",
    ),
    "population": dict(
        url=POPULATION_URL, path=POPULATION_JSON,
        description="World Bank total population (SP.POP.TOTL), 1960-2009",
        publisher="World Bank Open Data (CC BY 4.0), data.worldbank.org",
        citation="World Bank, World Development Indicators: Population, total (SP.POP.TOTL)",
    ),
}


def _download(url: str, dest: Path) -> None:
    """Stream ``url`` to ``dest`` via a temporary file, so an interrupted download
    never leaves a truncated file that later runs would mistake for the cache."""
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120, headers={"User-Agent": USER_AGENT}) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
    tmp.replace(dest)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def _check(key: str, path: Path) -> dict:
    """Sanity-check a downloaded file; return facts worth recording in the manifest."""
    if key == "meteorites":
        with open(path, encoding="utf-8") as fh:
            header = fh.readline().strip()
            rows = sum(1 for _ in fh)
        if not header.startswith("name,id,nametype,recclass,mass (g),fall,year"):
            raise RuntimeError(f"unexpected meteorite CSV header: {header!r}")
        digest = _sha256(path)
        if digest != METEORITE_SHA256 or rows != METEORITE_ROWS:
            print(f"WARNING: the meteorite file differs from the pinned snapshot "
                  f"({rows:,} rows, sha256 {digest[:12]}...). NASA has changed it; "
                  "the numbers in the README may no longer match.", file=sys.stderr)
        return {"rows": rows}
    if key.startswith("countries"):
        features = json.loads(path.read_text())["features"]
        return {"features": len(features)}
    if key == "population":
        meta, rows = json.loads(path.read_text())
        if meta["pages"] != 1:
            raise RuntimeError("World Bank response was paginated; raise per_page")
        return {"records": len(rows), "worldbank_last_updated": meta.get("lastupdated")}
    return {}


def fetch_all(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Make sure every source is on disk; return the manifest."""
    manifest = load_manifest() if MANIFEST_PATH.exists() else {}
    changed = False
    for key, src in SOURCES.items():
        path: Path = src["path"]
        if path.exists() and key in manifest and not refresh:
            continue
        if not quiet:
            print(f"downloading {src['url']}")
        _download(src["url"], path)
        manifest[key] = {
            "description": src["description"],
            "publisher": src["publisher"],
            "citation": src["citation"],
            "url": src["url"],
            "file": path.name,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes": path.stat().st_size,
            "sha256": _sha256(path),
            **_check(key, path),
        }
        changed = True
    if changed:
        MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    if not quiet:
        for key, m in manifest.items():
            print(f"  {key:<15} {m['bytes'] / 1e6:6.2f} MB  retrieved {m['retrieved_utc'][:10]}"
                  f"  {m['file']}")
    return manifest


def load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="download again even if a file is already cached")
    args = parser.parse_args(argv)
    fetch_all(refresh=args.refresh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
