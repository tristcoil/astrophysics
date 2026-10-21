"""Download and cache everything the analysis needs.

Two things are fetched, both free and public:

1. **The NASA Exoplanet Archive's Planetary Systems Composite Parameters table**
   (``pscomppars``, Caltech/IPAC), queried over the archive's Table Access
   Protocol (TAP) service. This is the list of every confirmed exoplanet, one
   row per planet. It is cached exactly as the archive returns it, as CSV, next
   to a small JSON record of the query and the time it ran -- the archive
   changes weekly, so the retrieval date is part of the result.

2. **Planet-interior model curves** from Zeng et al. (2019, PNAS 116, 9723):
   the radius a planet of a given mass would have if it were pure iron, Earth-
   like rock, water, and so on. These are published models, used only as
   reference lines on the mass-radius diagram.

Everything is cached under ``data/`` and nothing already present is fetched
again, so once this has run the whole project works offline.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from config import (
    CATALOG_META_PATH, CATALOG_PATH, COLUMNS, MODEL_DIR, MODEL_META_PATH, TABLE,
    TABLE_DOI, TAP_URL, ZENG_BASE_URL, ZENG_MODELS, ZENG_WATER_BY_TEMPERATURE,
)

#: The archive's query service can take a minute for the full table.
TIMEOUT_S = 300


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def build_query(columns: list[str] = COLUMNS) -> str:
    """The ADQL query sent to the archive. Recorded verbatim in the cache."""
    return f"select {','.join(columns)} from {TABLE}"


def fetch_catalog(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Download the composite-parameters table, unless it is already cached.

    Returns the provenance record (query, retrieval time, row count).

    The response is checked before it is written: a TAP service reports errors
    as a normal HTTP 200 with an error document in the body, so a status check
    alone would happily cache an error message as if it were the catalogue.
    """
    if CATALOG_PATH.exists() and CATALOG_META_PATH.exists() and not refresh:
        meta = json.loads(CATALOG_META_PATH.read_text())
        if not quiet:
            print(f"  catalogue cached ({meta['rows']:,} planets, "
                  f"retrieved {meta['retrieved_utc']})")
        return meta

    query = build_query()
    if not quiet:
        print(f"  querying {TAP_URL}\n  table {TABLE}, {len(COLUMNS)} columns...")
    response = requests.get(
        TAP_URL, params={"query": query, "format": "csv"}, timeout=TIMEOUT_S,
    )
    response.raise_for_status()
    text = response.text

    header = next(csv.reader(io.StringIO(text)))
    if header != COLUMNS:
        raise RuntimeError(
            "archive response does not look like the requested table; first "
            f"200 characters:\n{text[:200]}"
        )
    rows = pd.read_csv(io.StringIO(text), low_memory=False)
    if len(rows) < 1000:
        raise RuntimeError(f"only {len(rows)} rows returned -- refusing to cache")

    # Write to a temporary name first so an interrupted download can never
    # leave a truncated file where the analysis expects a complete one.
    tmp = CATALOG_PATH.with_name(CATALOG_PATH.name + ".part")
    tmp.write_text(text)
    tmp.replace(CATALOG_PATH)

    meta = {
        "retrieved_utc": _now_utc(),
        "source": "NASA Exoplanet Archive, Planetary Systems Composite Parameters",
        "table": TABLE,
        "doi": TABLE_DOI,
        "service": TAP_URL,
        "query": query,
        "rows": int(len(rows)),
        "bytes": len(text.encode()),
    }
    CATALOG_META_PATH.write_text(json.dumps(meta, indent=2))
    if not quiet:
        print(f"  wrote {CATALOG_PATH.name}: {meta['rows']:,} planets, "
              f"{meta['bytes'] / 1e6:.1f} MB")
    return meta


def fetch_models(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Download the Zeng et al. (2019) mass-radius tables that are not cached.

    Each file is two whitespace-separated columns: mass and radius, both in
    Earth units. The shape is checked before anything is written.
    """
    meta = json.loads(MODEL_META_PATH.read_text()) if MODEL_META_PATH.exists() else {}
    meta.setdefault("source", "Zeng et al. 2019, PNAS 116, 9723 (tables from Li Zeng's website)")
    meta.setdefault("files", {})

    filenames = sorted({f for f, _ in ZENG_MODELS.values()} | set(ZENG_WATER_BY_TEMPERATURE.values()))
    for filename in filenames:
        path = MODEL_DIR / filename
        if path.exists() and not refresh:
            continue
        url = ZENG_BASE_URL + filename
        if not quiet:
            print(f"  downloading {filename}")
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        table = np.loadtxt(io.StringIO(response.text))
        if table.ndim != 2 or table.shape[1] != 2 or len(table) < 5:
            raise RuntimeError(f"{url} is not a two-column mass-radius table")
        path.write_text(response.text)
        meta["files"][filename] = {"url": url, "retrieved_utc": _now_utc()}

    MODEL_META_PATH.write_text(json.dumps(meta, indent=2))
    if not quiet:
        print(f"  {len(filenames)} model curves cached in {MODEL_DIR.name}/")
    return meta


# ---------------------------------------------------------------------------
# Loading from the cache
# ---------------------------------------------------------------------------

def load_catalog_raw() -> pd.DataFrame:
    """The cached archive table, exactly as downloaded."""
    if not CATALOG_PATH.exists():
        raise FileNotFoundError(f"{CATALOG_PATH} not found -- run `python fetch_data.py` first")
    return pd.read_csv(CATALOG_PATH, low_memory=False)


def load_catalog_meta() -> dict:
    return json.loads(CATALOG_META_PATH.read_text())


def load_model_file(filename: str) -> tuple[np.ndarray, np.ndarray]:
    """A Zeng et al. model table as ``(mass, radius)`` arrays, sorted by mass."""
    path = MODEL_DIR / filename
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    mass, radius = np.loadtxt(path, unpack=True)
    order = np.argsort(mass)
    return mass[order], radius[order]


def load_model(key: str) -> tuple[np.ndarray, np.ndarray]:
    """One of the named curves in ``config.ZENG_MODELS``."""
    return load_model_file(ZENG_MODELS[key][0])


def ensure_data(*, quiet: bool = False) -> None:
    """Fetch anything missing from the cache; a no-op when it is complete."""
    fetch_catalog(quiet=quiet)
    fetch_models(quiet=quiet)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download even if cached (the archive updates weekly)")
    args = parser.parse_args(argv)

    print("NASA Exoplanet Archive")
    fetch_catalog(refresh=args.refresh)
    print("\nPlanet-interior models (Zeng et al. 2019)")
    fetch_models(refresh=args.refresh)
    print("\nDone. Cache is warm; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
