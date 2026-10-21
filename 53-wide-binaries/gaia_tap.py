"""Run ADQL queries against Gaia DR3, trying several archives, and cache them.

Copied from project L5 (31-hypervelocity-survivors/gaia_tap.py), where it
fetched whole-sky Gaia DR3 cuts patch by patch on 4 Oct 2026.

Every query is a template with a ``{table}`` placeholder for gaia_source's
name on each service (``gaiadr3.gaia_source`` at ESA and ARI,
``gaia_dr3.gaia_source`` at NOIRLab) and an optional ``{cols}`` placeholder.
Results are written to data/gaia/<name>.fits with a JSON sidecar holding the
exact query, the service that answered, the row count and the date, and are
never fetched twice.

Two silent failure modes are guarded against (both seen in earlier projects,
see 19-galactic-archaeology/gaia_archive.py): a result truncated at the
service's row limit, and column names whose letter case varies between
services. All column names are lower-cased.
"""

from __future__ import annotations

import json
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyvo
from astropy.table import Table

from config import GAIA_DIR, GAIA_SERVICES

MAXREC = 2_000_000


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _normalise(table: Table) -> Table:
    for name in table.colnames:
        if name != name.lower():
            table.rename_column(name, name.lower())
    # Object-typed (variable-length string) columns cannot go to FITS; masked
    # float columns are filled with nan so FITS round-trips cleanly.
    for name in table.colnames:
        col = table[name]
        if col.dtype.kind == "O":
            table[name] = np.array(["" if v is None or v is np.ma.masked else str(v)
                                    for v in np.asarray(col, dtype=object)])
            continue
        if getattr(col, "mask", None) is not None and np.any(col.mask):
            if col.dtype.kind == "f":
                table[name] = col.filled(np.nan)
    return table


def _run_async(svc, adql: str, timeout: float):
    """Submit a UWS job with the service's longest allowed run time. A
    whole-catalogue scan takes 15-30 min, beyond some sync limits (ARI: 30 min
    sync, 4 h async; NOIRLab: a 10-min gateway timeout on sync)."""
    job = svc.submit_job(adql, maxrec=MAXREC)
    try:
        try:
            job.execution_duration = int(timeout)
        except Exception:  # noqa: BLE001 - not every service lets us change it
            pass
        job.run()
        job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=timeout)
        job.raise_if_error()
        return job.fetch_result()
    finally:
        try:
            job.delete()
        except Exception:  # noqa: BLE001
            pass


def run_query(adql_template: str, *, cols: str = "*", timeout: float = 3600.0,
              services=None, quiet: bool = False, scan: bool = False) -> tuple[Table, dict]:
    """Run the query on the first service that answers; return the table and
    a provenance record. ``scan=True`` (a query that must read the whole
    catalogue) goes through the asynchronous job queue."""
    errors = []
    for name, url, table_name in (services or GAIA_SERVICES):
        adql = adql_template.format(table=table_name, cols=cols)
        t0 = time.time()
        try:
            svc = pyvo.dal.TAPService(url)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if scan:
                    result = _run_async(svc, adql, max(timeout, 7200.0))
                else:
                    result = svc.run_sync(adql, maxrec=MAXREC, timeout=timeout)
            table = _normalise(result.to_table())
        except Exception as exc:  # noqa: BLE001 - any failure: try the next service
            msg = f"{name}: {type(exc).__name__}: {str(exc)[:200]}"
            errors.append(msg)
            if not quiet:
                print(f"    {msg} -> trying the next service", flush=True)
            continue
        if len(table) >= MAXREC:
            errors.append(f"{name}: hit the {MAXREC:,}-row limit")
            continue
        if getattr(result, "query_status", "OK") == "OVERFLOW":
            errors.append(f"{name}: OVERFLOW")
            continue
        prov = dict(service=name, url=url, adql=adql, rows=len(table),
                    seconds=round(time.time() - t0, 1), retrieved=now_utc())
        return table, prov
    raise RuntimeError("every Gaia service failed:\n  " + "\n  ".join(errors))


def cached_query(name: str, adql_template: str, *, cols: str = "*", quiet: bool = False,
                 **kw) -> Table:
    """run_query, cached as data/gaia/<name>.fits (+ .json provenance)."""
    path = GAIA_DIR / f"{name}.fits"
    meta = path.with_suffix(".json")
    if path.exists() and meta.exists():
        return Table.read(path)
    if not quiet:
        print(f"  Gaia query '{name}' ...", flush=True)
    table, prov = run_query(adql_template, cols=cols, quiet=quiet, **kw)
    tmp = path.with_suffix(".part.fits")
    table.write(tmp, format="fits", overwrite=True)
    tmp.rename(path)
    meta.write_text(json.dumps(prov, indent=1))
    if not quiet:
        print(f"    {prov['rows']:,} rows from {prov['service']} in {prov['seconds']:.0f} s", flush=True)
    return table


def provenance(name: str) -> dict:
    meta = (GAIA_DIR / f"{name}.json")
    return json.loads(meta.read_text()) if meta.exists() else {}


def id_list(ids) -> str:
    return ", ".join(str(int(i)) for i in ids)


# --------------------------------------------------------------------------
# Whole-catalogue scans in sky patches
# --------------------------------------------------------------------------
# Every Gaia source_id encodes its position: source_id // 2**35 is its
# HEALPix level-12 pixel (nested), so a source_id range is a patch of sky.
# Level 1 cuts the sky into 48 patches; a patch that fails is split in four.

def _patch_range(level: int, pixel: int) -> tuple[int, int]:
    width = 2**35 * 4 ** (12 - level)
    return pixel * width, (pixel + 1) * width


def chunked_scan(name: str, adql_template: str, *, cols: str, service, base_level: int = 1,
                 max_level: int = 4, workers: int = 3, quiet: bool = False) -> Table:
    """Run ``adql_template`` (which must leave ``{lo}``/``{hi}`` placeholders
    for a source_id range, plus ``{table}``/``{cols}``) patch by patch on one
    service, caching each patch, and save the stacked result like
    cached_query does."""
    from concurrent.futures import ThreadPoolExecutor
    from astropy.table import vstack
    sname, url, table_name = service
    path = GAIA_DIR / f"{name}.fits"
    meta = path.with_suffix(".json")
    if path.exists() and meta.exists():
        return Table.read(path)
    chunk_dir = GAIA_DIR / "chunks" / name
    chunk_dir.mkdir(parents=True, exist_ok=True)
    svc = pyvo.dal.TAPService(url)
    t_start = time.time()

    def one(level, pixel):
        cpath = chunk_dir / f"L{level}_{pixel:06d}.fits"
        if cpath.exists():
            return [cpath]
        lo, hi = _patch_range(level, pixel)
        adql = (adql_template.replace("{lo}", str(lo)).replace("{hi}", str(hi))
                .format(table=table_name, cols=cols))
        t0 = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                tab = _normalise(svc.run_sync(adql, maxrec=MAXREC, timeout=700).to_table())
        except Exception as exc:  # noqa: BLE001
            if level >= max_level:
                raise
            if not quiet:
                print(f"    {name} L{level}/{pixel}: {type(exc).__name__} -> splitting", flush=True)
            out = []
            for child in range(4 * pixel, 4 * pixel + 4):
                out += one(level + 1, child)
            return out
        tmp = cpath.with_suffix(".part.fits")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tab.write(tmp, format="fits", overwrite=True)
        tmp.rename(cpath)
        if not quiet:
            print(f"    {name} L{level}/{pixel:>3}: {len(tab):>7,} rows in {time.time() - t0:5.0f} s", flush=True)
        return [cpath]

    with ThreadPoolExecutor(workers) as pool:
        parts = list(pool.map(lambda p: one(base_level, p), range(12 * 4 ** base_level)))
    paths = [p for ps in parts for p in ps]
    tables = [Table.read(p) for p in paths]
    tables = [t for t in tables if len(t)]
    stacked = vstack(tables, metadata_conflicts="silent")
    if len(np.unique(stacked["source_id"])) != len(stacked):
        raise RuntimeError(f"{name}: duplicate source_ids across patches")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        stacked.write(path, format="fits", overwrite=True)
    prov = dict(service=sname, url=url, adql=adql_template.format(table=table_name, cols=cols,
                                                                   lo="{lo}", hi="{hi}"),
                rows=len(stacked), patches=len(paths), seconds=round(time.time() - t_start, 1),
                retrieved=now_utc(),
                note="run patch by patch: source_id >= lo and source_id < hi (HEALPix level-1 "
                     "patches, split further on failure)")
    meta.write_text(json.dumps(prov, indent=1))
    if not quiet:
        print(f"    {name}: {len(stacked):,} rows from {len(paths)} patches on {sname}", flush=True)
    return stacked


# --------------------------------------------------------------------------
# Two-step scans: find source_ids in a narrow table, then fetch full rows
# --------------------------------------------------------------------------
# A whole-catalogue scan of gaia_source (~150 columns per row) takes 30-60
# minutes on any mirror. ARI also serves gaiadr3.gaia_source_lite (51
# columns), so the cut runs there and only the matching source_ids come back;
# their full records are then fetched by source_id, which is indexed.

def lite_then_full(name: str, lite_where: str, *, cols: str, quiet: bool = False,
                   batch: int = 1500) -> Table:
    from astropy.table import vstack
    path = GAIA_DIR / f"{name}.fits"
    meta = path.with_suffix(".json")
    if path.exists() and meta.exists():
        return Table.read(path)
    ari = next(s for s in GAIA_SERVICES if s[0].startswith("ARI"))
    lite = (ari[0], ari[1], "gaiadr3.gaia_source_lite")
    lite_adql = f"select source_id from gaiadr3.gaia_source_lite where {lite_where}"
    # One whole-table query ran past ARI's 30-min sync limit on 4 Oct 2026,
    # while a 1/48 sky patch (a source_id range) took ~45 s: so patch by patch.
    ids = chunked_scan(f"{name}_ids", "select {cols} from {table} where source_id >= {lo} and "
                       "source_id < {hi} and " + lite_where, cols="source_id", service=lite,
                       quiet=quiet)
    lite_prov = json.loads((GAIA_DIR / f"{name}_ids.json").read_text())
    sids = np.asarray(ids["source_id"], dtype=np.int64)
    parts = []
    for k in range(0, len(sids), batch):
        q = f"select {cols} from {{table}} where source_id in ({id_list(sids[k:k + batch])})"
        t, prov = run_query(q, quiet=quiet)
        parts.append(t)
    full = vstack(parts, metadata_conflicts="silent") if parts else Table()
    if len(full) != len(sids):
        raise RuntimeError(f"{name}: {len(sids)} ids but {len(full)} full rows")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full.write(path, format="fits", overwrite=True)
    meta.write_text(json.dumps(dict(
        service=lite_prov["service"], url=lite_prov["url"], adql=lite_adql, rows=len(full),
        seconds=lite_prov["seconds"], retrieved=now_utc(),
        note=f"step 1 selects source_ids from gaiadr3.gaia_source_lite; step 2 fetches these "
             f"columns from gaia_source by source_id: {cols}"), indent=1))
    if not quiet:
        print(f"    {name}: {len(full):,} full rows", flush=True)
    return full
