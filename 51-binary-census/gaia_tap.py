"""A small, careful client for the ESA Gaia archive's synchronous TAP endpoint.

Adapted unchanged in substance from H3 (``19-galactic-archaeology/gaia_archive.py``).
Lessons it encodes, from earlier projects in this repository:

* ESA's *asynchronous* job queue works but can sit for minutes before a job
  starts, while *synchronous* queries run at once -- but give up after about a
  minute. So big samples are fetched as many small synchronous queries, one
  patch of sky at a time (a range of source_id, which encodes the HEALPix
  pixel), and a patch that times out is split into four and retried.
* Results come back as FITS (28% smaller than the default VOTable), are cached
  chunk by chunk (so an interrupted download resumes where it stopped), and are
  checked for two silent failure modes: the archive's 3,000,000-row output cap,
  which truncates without warning in FITS output, and column names whose letter
  case varies between responses (which makes stacking tables silently split one
  column into two).
"""

from __future__ import annotations

import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from config import CHUNK_DIR, GAIA_SYNC_URL, HEALPIX_BASE_LEVEL, HEALPIX_MAX_LEVEL

#: ESA's hard output limit (from the service's /capabilities document).
MAXREC = 3_000_000
#: source_id // 2**35 is the HEALPix level-12 nested pixel.
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


def sync_query(adql: str, *, timeout: float = 150.0, attempts: int = 3,
               fmt: str = "fits") -> Table:
    """Run one synchronous ADQL query and return an astropy Table.

    Column names are lower-cased. A result with exactly MAXREC rows is
    refused, because it has almost certainly been truncated.
    """
    data = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt,
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

    if fmt != "fits":
        raise ValueError("only FITS output is parsed here")
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


@dataclass
class ChunkLog:
    level: int
    pixel: int
    rows: int
    seconds: float


def _chunk_path(name: str, level: int, pixel: int) -> Path:
    return CHUNK_DIR / name / f"L{level}_{pixel:06d}.fits"


def _fetch_pixel(name: str, template: str, level: int, pixel: int, log: list[ChunkLog],
                 quiet: bool) -> list[Path]:
    """Fetch one sky pixel, splitting it into its four children on timeout."""
    path = _chunk_path(name, level, pixel)
    if path.exists():
        n = fits.getheader(path, 1)["NAXIS2"]
        log.append(ChunkLog(level, pixel, n, 0.0))
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
    log.append(ChunkLog(level, pixel, len(table), dt))
    if not quiet:
        print(f"    L{level} pixel {pixel:>4}: {len(table):>8,} rows in {dt:5.1f} s", flush=True)
    return [path]


def fetch_all_sky(name: str, template: str, *, quiet: bool = False) -> tuple[Table, list[ChunkLog]]:
    """Run ``template`` (with ``{lo}``/``{hi}`` source_id placeholders) over the
    whole sky, one HEALPix pixel at a time, and stack the results."""
    log: list[ChunkLog] = []
    paths: list[Path] = []
    npix = 12 * 4**HEALPIX_BASE_LEVEL
    for pixel in range(npix):
        paths += _fetch_pixel(name, template, HEALPIX_BASE_LEVEL, pixel, log, quiet)
    tables = [Table.read(p, format="fits") for p in paths]
    names = [set(t.colnames) for t in tables]
    if any(n != names[0] for n in names):
        raise RuntimeError(f"{name}: chunks disagree on their columns")
    return vstack(tables, metadata_conflicts="silent"), log


def _count_pixel(template: str, level: int, pixel: int) -> int:
    lo, hi = healpix_source_id_range(level, pixel)
    try:
        return int(np.asarray(sync_query(template.format(lo=lo, hi=hi)).columns[0])[0])
    except TransientError:
        if level >= HEALPIX_MAX_LEVEL:
            raise
        return sum(_count_pixel(template, level + 1, c) for c in range(4 * pixel, 4 * pixel + 4))


def count_all_sky(template: str) -> int:
    """Server-side COUNT(*) of a chunked query (template must SELECT COUNT(*)),
    splitting sky pixels that time out exactly as the download does."""
    return sum(_count_pixel(template, HEALPIX_BASE_LEVEL, p)
               for p in range(12 * 4**HEALPIX_BASE_LEVEL))


def write_with_provenance(table: Table, path: Path, provenance: dict) -> None:
    """Save a table as FITS plus a JSON sidecar recording how it was made."""
    tmp = path.with_suffix(".part")
    table.write(tmp, format="fits", overwrite=True)
    tmp.rename(path)
    provenance = {"rows": len(table), "written_utc": now_utc(), **provenance}
    path.with_suffix(".json").write_text(json.dumps(provenance, indent=2) + "\n")
