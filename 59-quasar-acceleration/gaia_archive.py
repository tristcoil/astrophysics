"""A small, careful client for the ESA Gaia archive's synchronous TAP endpoint.

Adapted from project H3 (19-galactic-archaeology/gaia_archive.py), whose
lessons still hold in October 2026:

* Synchronous queries start at once but give up after a minute or less (HTTP
  408); even a COUNT(*) over the whole quasar table timed out on 19 Oct 2026.
  So the sky is fetched one HEALPix pixel at a time (a source_id range), and a
  pixel that times out is split into its four children and retried.
* Results come back as FITS (smaller than VOTable) and are cached chunk by
  chunk, so an interrupted download resumes where it stopped.
* Two silent failure modes are checked: the archive's 3,000,000-row output cap
  (FITS output truncates without warning) and column names whose letter case
  varies between responses (stacking such tables silently splits a column).
"""

from __future__ import annotations

import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from config import CHUNK_DIR, GAIA_SYNC_URL, HEALPIX_BASE_LEVEL, HEALPIX_MAX_LEVEL

#: ESA's hard output limit (from the service's /capabilities document).
MAXREC = 3_000_000
#: source_id // 2**35 is the HEALPix level-12 nested pixel of the source.
_LEVEL12_SHIFT = 2**35


class QueryError(RuntimeError):
    """The archive rejected the query itself (bad ADQL): retrying cannot help."""


class TransientError(RuntimeError):
    """Timeout, overload or network trouble: worth retrying, or splitting."""


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _error_text(body: bytes) -> str:
    text = body.decode("utf-8", "replace")
    start = text.find('value="ERROR">')
    if start >= 0:
        text = text[start + 14:]
        text = text[: text.find("</INFO>")]
    return " ".join(text.split())[:600]


def sync_query(adql: str, *, timeout: float = 180.0, attempts: int = 3) -> Table:
    """Run one synchronous ADQL query and return an astropy Table.

    Column names are lower-cased. A result with MAXREC rows is refused,
    because it has almost certainly been truncated.
    """
    data = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "fits",
        "MAXREC": str(MAXREC), "QUERY": adql,
    }).encode()
    last = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(GAIA_SYNC_URL, data=data)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read()
            break
        except urllib.error.HTTPError as exc:
            body = exc.read()
            if exc.code == 400:
                raise QueryError(_error_text(body)) from None
            last = f"HTTP {exc.code}"
            # 408 = the query ran past the endpoint's time limit: retrying the
            # same query is pointless, the caller should split it instead.
            if exc.code == 408:
                raise TransientError(last) from None
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = f"{type(exc).__name__}: {exc}"
        if attempt < attempts:
            time.sleep(5 * attempt)
    else:
        raise TransientError(last)

    if not body.startswith(b"SIMPLE"):
        # The archive reports some failures as a VOTable with status 200.
        raise TransientError(f"not a FITS response: {_error_text(body)}")
    with fits.open(io.BytesIO(body)) as hdul:
        table = Table(hdul[1].data)
    for name in table.colnames:
        if name != name.lower():
            table.rename_column(name, name.lower())
    if len(table) >= MAXREC:
        raise TransientError(f"result hit the {MAXREC:,}-row cap; it is truncated")
    return table


def healpix_source_id_range(level: int, pixel: int) -> tuple[int, int]:
    """[lo, hi) source_id bounds of one nested HEALPix pixel at ``level``."""
    width = _LEVEL12_SHIFT * 4 ** (12 - level)
    return pixel * width, (pixel + 1) * width


def _chunk_path(name: str, level: int, pixel: int) -> Path:
    return CHUNK_DIR / name / f"L{level}_{pixel:06d}.fits"


def _fetch_pixel(name: str, template: str, level: int, pixel: int, log: list,
                 quiet: bool) -> list[Path]:
    """Fetch one sky pixel, splitting it into its four children on timeout."""
    path = _chunk_path(name, level, pixel)
    if path.exists():
        log.append({"level": level, "pixel": pixel,
                    "rows": int(fits.getheader(path, 1)["NAXIS2"]), "seconds": 0.0})
        return [path]
    lo, hi = healpix_source_id_range(level, pixel)
    t0 = time.time()
    try:
        table = sync_query(template.format(lo=lo, hi=hi))
    except TransientError as exc:
        if level >= HEALPIX_MAX_LEVEL:
            raise
        if not quiet:
            print(f"    L{level} pixel {pixel}: {exc} -> splitting in four", flush=True)
        paths = []
        for child in range(4 * pixel, 4 * pixel + 4):
            paths += _fetch_pixel(name, template, level + 1, child, log, quiet)
        return paths
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    table.write(tmp, format="fits", overwrite=True)
    tmp.rename(path)
    dt = time.time() - t0
    log.append({"level": level, "pixel": pixel, "rows": len(table), "seconds": round(dt, 1)})
    if not quiet:
        print(f"    L{level} pixel {pixel:>4}: {len(table):>7,} rows in {dt:5.1f} s", flush=True)
    return [path]


def _sky_tasks(name: str, start_level: int) -> list[tuple[int, int]]:
    """Pixels to fetch: every pixel at ``start_level``, except where a chunk
    cached at a coarser level already covers it."""
    tasks = []
    for pixel in range(12 * 4**HEALPIX_BASE_LEVEL):
        if start_level == HEALPIX_BASE_LEVEL or _chunk_path(name, HEALPIX_BASE_LEVEL, pixel).exists():
            tasks.append((HEALPIX_BASE_LEVEL, pixel))
            continue
        k = 4 ** (start_level - HEALPIX_BASE_LEVEL)
        tasks += [(start_level, child) for child in range(k * pixel, k * (pixel + 1))]
    return tasks


def fetch_all_sky(name: str, template: str, *, start_level: int = HEALPIX_BASE_LEVEL,
                  workers: int = 1, quiet: bool = False) -> tuple[Table, list]:
    """Run ``template`` (with ``{lo}``/``{hi}`` source_id placeholders) over the
    whole sky, one HEALPix pixel at a time, and stack the results.

    ``start_level`` skips coarse pixels known to time out; ``workers`` runs that
    many queries at once (the archive's per-row lookup, ~1.7 ms on 19 Oct 2026,
    not the network, sets the pace)."""
    log: list = []
    tasks = _sky_tasks(name, start_level)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(lambda t: _fetch_pixel(name, template, t[0], t[1], log, quiet),
                                tasks))
    paths = [p for group in results for p in group]
    tables = [Table.read(p, format="fits") for p in paths]
    names = [set(t.colnames) for t in tables]
    if any(n != names[0] for n in names):
        raise RuntimeError(f"{name}: chunks disagree on their columns")
    return vstack(tables, metadata_conflicts="silent"), log


def _count_pixel(template: str, level: int, pixel: int) -> dict:
    lo, hi = healpix_source_id_range(level, pixel)
    try:
        table = sync_query(template.format(lo=lo, hi=hi))
    except TransientError:
        if level >= HEALPIX_MAX_LEVEL:
            raise
        total: dict = {}
        for c in range(4 * pixel, 4 * pixel + 4):
            for key, n in _count_pixel(template, level + 1, c).items():
                total[key] = total.get(key, 0) + n
        return total
    cols = table.colnames
    if len(cols) == 1:  # plain COUNT(*)
        return {"all": int(np.asarray(table[cols[0]])[0])}
    return {int(k): int(n) for k, n in zip(table[cols[0]], table[cols[1]])}


def count_all_sky(template: str, *, start_level: int = HEALPIX_BASE_LEVEL,
                  workers: int = 1) -> dict:
    """Server-side COUNT(*) of a chunked query, splitting sky pixels that time
    out exactly as the download does. The template either SELECTs COUNT(*)
    alone (returns {"all": n}) or a key column and COUNT(*) with GROUP BY
    (returns {key: n})."""
    k = 4 ** (start_level - HEALPIX_BASE_LEVEL)
    pixels = range(12 * 4**HEALPIX_BASE_LEVEL * k)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(lambda p: _count_pixel(template, start_level, p), pixels))
    total: dict = {}
    for part in parts:
        for key, n in part.items():
            total[key] = total.get(key, 0) + n
    return total


def write_with_provenance(table: Table, path: Path, provenance: dict) -> None:
    """Save a table as FITS plus a JSON sidecar recording how it was made."""
    tmp = path.with_suffix(".part")
    table.write(tmp, format="fits", overwrite=True)
    tmp.rename(path)
    provenance = {"rows": len(table), "written_utc": now_utc(), **provenance}
    path.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")


def async_query(adql: str, *, poll_s: float = 15.0, timeout_s: float = 7200.0,
                quiet: bool = False) -> Table:
    """Run one query through ESA's asynchronous TAP queue (jobs can wait a few
    minutes before starting, but are not cut off after ~15-60 s like sync)."""
    base = GAIA_SYNC_URL.rsplit("/", 1)[0] + "/async"
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "fits",
                                   "PHASE": "RUN", "MAXREC": str(MAXREC), "QUERY": adql}).encode()

    class _NoRedirect(urllib.request.HTTPRedirectHandler):
        def redirect_request(self, *args, **kwargs):
            return None

    opener = urllib.request.build_opener(_NoRedirect)
    try:
        resp = opener.open(urllib.request.Request(base, data=data), timeout=120)
        job_url = resp.headers.get("Location") or resp.geturl()
    except urllib.error.HTTPError as exc:
        if exc.code not in (301, 302, 303):
            raise TransientError(f"async submit: HTTP {exc.code} {_error_text(exc.read())}") from None
        job_url = exc.headers["Location"]
    if not quiet:
        print(f"    async job {job_url.rsplit('/', 1)[-1]} submitted", flush=True)
    t0 = time.time()
    phase = ""
    while time.time() - t0 < timeout_s:
        with urllib.request.urlopen(job_url + "/phase", timeout=60) as r:
            phase = r.read().decode().strip()
        if phase in ("COMPLETED", "ERROR", "ABORTED"):
            break
        time.sleep(poll_s)
    if phase != "COMPLETED":
        detail = ""
        try:
            with urllib.request.urlopen(job_url + "/error", timeout=60) as r:
                detail = _error_text(r.read())
        except Exception:  # noqa: BLE001
            pass
        raise TransientError(f"async job ended in phase {phase!r} {detail}")
    with urllib.request.urlopen(job_url + "/results/result", timeout=600) as r:
        body = r.read()
    with fits.open(io.BytesIO(body)) as hdul:
        table = Table(hdul[1].data)
    for name in table.colnames:
        if name != name.lower():
            table.rename_column(name, name.lower())
    if len(table) in (2_000, 50_000, 100_000, MAXREC):
        # an answer of exactly a default page size has been truncated
        raise TransientError(f"async result has exactly {len(table):,} rows: truncated")
    if not quiet:
        print(f"    async job done in {time.time() - t0:.0f} s: {len(table):,} rows", flush=True)
    return table
