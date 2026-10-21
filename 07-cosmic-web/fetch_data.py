"""Download and cache everything the analysis needs.

Two things are fetched, both free and public:

1. **The SDSS Main Galaxy Sample** from SkyServer (Data Release 17): position,
   redshift, brightness and colour for every galaxy that the Sloan Digital Sky
   Survey targeted for its flux-limited galaxy survey and got a reliable
   redshift for -- about 657,000 of them. This is the data the whole project
   is built on.

2. **The Abell catalogue of rich clusters** (Abell, Corwin & Olowin 1989), from
   VizieR. It is used for one thing only: once our own group finder has picked
   out the richest clusters, we look up what astronomers already call them. It
   plays no part in finding them.

Everything is cached under ``data/``. The survey is fetched in right-ascension
strips, one query each, and every strip is saved as the CSV SkyServer served,
so an interrupted download resumes where it stopped and nothing already on disk
is fetched twice. Once the cache is warm the analysis runs offline.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from config import (
    ABELL_PATH,
    CATALOG_PATH,
    DATA_RELEASE,
    FETCH_Z_MAX,
    FETCH_Z_MIN,
    MGS_R_LIMIT,
    MGS_TARGET_BITS,
    RA_CHUNK_DEG,
    RAW_DIR,
    SKYSERVER_URL,
)

#: The query, one right-ascension strip at a time.
#:
#: ``SpecPhoto`` is SkyServer's pre-joined view of spectra and the imaging
#: they were targeted from. ``sciencePrimary = 1`` keeps the single best
#: spectrum of each object (some galaxies were observed more than once), and
#: ``mode = 1`` its primary photometric detection.
#:
#: ``zWarning = 0`` keeps only redshifts the SDSS pipeline has no doubts about.
#:
#: Magnitudes: ``petroMag_r`` (Petrosian) is the flux the survey was selected on,
#: so the flux limit is applied to it. ``dered_*`` are model magnitudes with
#: Milky Way dust removed; model magnitudes are the SDSS recommendation for
#: colours, because they measure every band through the same aperture.
QUERY_TEMPLATE = """
SELECT specObjID, ra, dec, z, zErr,
       petroMag_r, extinction_r,
       dered_u, dered_g, dered_r,
       modelMagErr_u, modelMagErr_g, modelMagErr_r
FROM SpecPhoto
WHERE (legacy_target1 & {bits}) != 0
  AND class = 'GALAXY'
  AND sciencePrimary = 1
  AND mode = 1
  AND zWarning = 0
  AND z BETWEEN {zmin} AND {zmax}
  AND petroMag_r - extinction_r < {rlim}
  AND ra >= {ra_lo} AND ra < {ra_hi}
"""

#: Columns and the dtype each is stored as. float32 would do for magnitudes,
#: but not for coordinates: at float32 precision a right ascension near 360
#: degrees is only good to ~0.1 arcsec. Keep float64 where it matters.
COLUMNS = {
    "specObjID": np.int64,
    "ra": np.float64,
    "dec": np.float64,
    "z": np.float64,
    "zErr": np.float32,
    "petroMag_r": np.float32,
    "extinction_r": np.float32,
    "dered_u": np.float32,
    "dered_g": np.float32,
    "dered_r": np.float32,
    "modelMagErr_u": np.float32,
    "modelMagErr_g": np.float32,
    "modelMagErr_r": np.float32,
}


def skyserver_sql(sql: str, *, retries: int = 4, timeout: float = 600.0) -> str:
    """Run one SQL query on SkyServer and return the CSV text.

    SkyServer reports SQL errors inside a 200 response, so the body is checked,
    not just the status code. Transient failures are retried with backoff.
    """
    last_error = None
    for attempt in range(retries):
        try:
            r = requests.get(
                SKYSERVER_URL,
                params={"cmd": " ".join(sql.split()), "format": "csv"},
                timeout=timeout,
            )
            r.raise_for_status()
            text = r.text
            if not text.startswith("#Table1"):
                raise RuntimeError(f"SkyServer returned an error: {text[:300]!r}")
            return text
        except (requests.RequestException, RuntimeError) as exc:
            last_error = exc
            wait = 5 * 2**attempt
            print(f"    attempt {attempt + 1} failed ({exc}); retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"SkyServer query failed after {retries} attempts: {last_error}")


def _parse(text: str) -> pd.DataFrame:
    """SkyServer CSV: a ``#Table1`` line, then a normal header and rows."""
    return pd.read_csv(io.StringIO(text), skiprows=1)


def _chunks() -> list[tuple[float, float]]:
    edges = np.arange(0.0, 360.0 + RA_CHUNK_DEG / 2, RA_CHUNK_DEG)
    return [(float(lo), float(hi)) for lo, hi in zip(edges[:-1], edges[1:])]


def _chunk_path(lo: float, hi: float):
    return RAW_DIR / f"mgs_ra{lo:05.1f}-{hi:05.1f}.csv"


def fetch_survey(*, refresh: bool = False, quiet: bool = False) -> None:
    """Download every right-ascension strip not already on disk."""
    chunks = _chunks()
    if not quiet:
        print(f"SDSS {DATA_RELEASE} Main Galaxy Sample, {len(chunks)} strips of "
              f"{RA_CHUNK_DEG:g} deg in right ascension")
    for i, (lo, hi) in enumerate(chunks, 1):
        path = _chunk_path(lo, hi)
        if path.exists() and not refresh:
            if not quiet:
                print(f"  [{i:2d}/{len(chunks)}] RA {lo:5.1f}-{hi:5.1f}: cached")
            continue
        sql = QUERY_TEMPLATE.format(
            bits=MGS_TARGET_BITS, zmin=FETCH_Z_MIN, zmax=FETCH_Z_MAX,
            rlim=MGS_R_LIMIT, ra_lo=lo, ra_hi=hi,
        )
        t0 = time.time()
        text = skyserver_sql(sql)
        n = text.count("\n") - 2
        # Write to a temporary name first so a crash never leaves a truncated
        # file that a later run would mistake for a complete strip.
        tmp = path.with_suffix(".part")
        tmp.write_text(text)
        tmp.replace(path)
        if not quiet:
            print(f"  [{i:2d}/{len(chunks)}] RA {lo:5.1f}-{hi:5.1f}: "
                  f"{n:6,d} galaxies in {time.time() - t0:4.1f}s")
        time.sleep(1.0)  # be polite to a free public service


def build_catalog(*, quiet: bool = False) -> dict:
    """Merge the cached strips into one typed catalogue, with provenance."""
    frames = []
    for lo, hi in _chunks():
        path = _chunk_path(lo, hi)
        if not path.exists():
            raise FileNotFoundError(f"{path.name} missing -- run `python fetch_data.py`")
        frame = _parse(path.read_text())
        if len(frame):  # strips crossing the Milky Way's disc hold no galaxies
            frames.append(frame)
    df = pd.concat(frames, ignore_index=True)

    missing = set(COLUMNS) - set(df.columns)
    if missing:
        raise RuntimeError(f"columns missing from SkyServer output: {sorted(missing)}")

    n_raw = len(df)
    df = df.drop_duplicates("specObjID")
    n_dupes = n_raw - len(df)

    # SkyServer marks a missing measurement with -9999. A galaxy lacking any
    # magnitude we need is dropped, and counted, rather than silently kept.
    mags = ["petroMag_r", "dered_u", "dered_g", "dered_r"]
    bad = (df[mags] < -1000).any(axis=1) | df[mags].isna().any(axis=1)
    n_bad = int(bad.sum())
    df = df[~bad]

    df = df.sort_values(["ra", "dec"]).reset_index(drop=True)
    arrays = {col: df[col].to_numpy(dtype=dtype) for col, dtype in COLUMNS.items()}

    meta = {
        "source": f"SDSS {DATA_RELEASE} SkyServer, view SpecPhoto",
        "url": SKYSERVER_URL,
        "query": " ".join(QUERY_TEMPLATE.split()),
        "retrieved_utc": _retrieval_date(),
        "n_rows_served": n_raw,
        "n_duplicates_dropped": n_dupes,
        "n_missing_photometry_dropped": n_bad,
        "n_galaxies": len(df),
    }
    np.savez_compressed(CATALOG_PATH, meta=json.dumps(meta), **arrays)
    if not quiet:
        print(f"  merged {n_raw:,} rows -> {len(df):,} galaxies "
              f"({n_dupes} duplicates, {n_bad} with missing photometry dropped)")
        print(f"  wrote {CATALOG_PATH.relative_to(CATALOG_PATH.parents[1])}")
    return meta


def _retrieval_date() -> str:
    """When the strips were downloaded: the newest file's modification time."""
    newest = max(p.stat().st_mtime for p in RAW_DIR.glob("mgs_ra*.csv"))
    return datetime.fromtimestamp(newest, timezone.utc).isoformat(timespec="seconds")


def load_catalog() -> tuple[dict, dict]:
    """Read the merged catalogue: (column arrays, provenance metadata)."""
    if not CATALOG_PATH.exists():
        raise FileNotFoundError(f"{CATALOG_PATH} not found -- run `python fetch_data.py` first")
    with np.load(CATALOG_PATH) as f:
        meta = json.loads(str(f["meta"]))
        arrays = {k: f[k] for k in f.files if k != "meta"}
    return arrays, meta


def fetch_abell(*, refresh: bool = False, quiet: bool = False) -> pd.DataFrame:
    """The northern Abell catalogue (VizieR VII/110A, table 3): 2,712 clusters.

    Covers everything north of declination -27, i.e. all of the SDSS footprint.
    Positions are VizieR's J2000 conversions of the original B1950 ones.
    Redshifts exist only for clusters someone has measured; the rest are blank.
    """
    if ABELL_PATH.exists() and not refresh:
        return pd.read_csv(ABELL_PATH, comment="#")

    from astroquery.vizier import Vizier

    if not quiet:
        print("Abell, Corwin & Olowin (1989) cluster catalogue from VizieR")
    v = Vizier(columns=["ACO", "_RA.icrs", "_DE.icrs", "z", "Rich", "Count", "BMtype"],
               row_limit=-1)
    table = v.get_catalogs("VII/110A/table3")[0]

    from astropy.coordinates import SkyCoord
    import astropy.units as u

    coords = SkyCoord(table["_RA.icrs"], table["_DE.icrs"], unit=(u.hourangle, u.deg))
    z = np.ma.filled(np.ma.asarray(table["z"], dtype=float), np.nan)
    df = pd.DataFrame({
        "aco": np.asarray(table["ACO"], dtype=int),
        "ra": coords.ra.deg,
        "dec": coords.dec.deg,
        "z": z,
        "richness_class": np.asarray(table["Rich"], dtype=int),
        "count": np.asarray(table["Count"], dtype=int),
        "bm_type": [str(s).strip() for s in table["BMtype"]],
    })
    header = (f"# Abell, Corwin & Olowin 1989, ApJS 70, 1; VizieR VII/110A/table3; "
              f"retrieved {datetime.now(timezone.utc).date().isoformat()}\n")
    ABELL_PATH.write_text(header + df.to_csv(index=False))
    if not quiet:
        print(f"  {len(df):,} clusters, {int(np.isfinite(z).sum()):,} with a redshift; "
              f"wrote {ABELL_PATH.name}")
    return df


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download even what is already cached")
    parser.add_argument("--skip-abell", action="store_true",
                        help="skip the Abell cluster catalogue")
    args = parser.parse_args(argv)

    fetch_survey(refresh=args.refresh)
    print("\nMerging strips")
    meta = build_catalog()

    if not args.skip_abell:
        print()
        fetch_abell(refresh=args.refresh)

    print(f"\nDone: {meta['n_galaxies']:,} galaxies cached. Analysis can now run offline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
