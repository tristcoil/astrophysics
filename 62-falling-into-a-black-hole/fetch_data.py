"""Download and cache every input this project uses (about 0.3 GB, ~10-20 minutes).

    python fetch_data.py               # everything
    python fetch_data.py stars flux    # only some sources

The black holes' masses come from earlier projects' results files (H1, I4, I5),
so the only downloads are the sky and a few reference tables:

1. **Gaia DR3 point stars** (``stars``): every source with G < 11 (1,247,240
   stars): position, G, BP and RP magnitudes. Each one is drawn as a point
   and lensed individually.
2. **Gaia DR3 summed-flux map** (``flux``): for every HEALPix level-10 pixel
   (3.4 arcmin, 12.6 million pixels) the number of fainter (G >= 11) sources
   and the sums of their G, BP and RP fluxes, aggregated on the archive's
   server. This is the glow of the Milky Way, made from 1.8 billion stars.
3. **Hipparcos** (``hipparcos``): the main catalogue (V, B-V) and Gaia DR3's
   official Hipparcos-2 cross-match, to replace Gaia's saturated photometry of
   the brightest naked-eye stars.
4. **Reference tables** (``reference``): the CIE 1931 colour-matching functions
   (CVRL, UCL) and Mamajek's dwarf colour-temperature sequence.
5. **Literature** (``literature``): arXiv PDFs of the papers we check against.

Queries go to ESA's Gaia archive as many small synchronous ADQL queries, one
HEALPix patch of sky at a time (source_id ranges); a patch that times out is
split into its four children (the recipe of H3, 19-galactic-archaeology).
Nothing is re-downloaded once present. ``data/manifest.json`` records each
file's source, size, SHA-256 and download date.
"""

from __future__ import annotations

import concurrent.futures as cf
import hashlib
import io
import json
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

import config as C

UA = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
MAXREC = 3_000_000                 # ESA's hard output cap (from /capabilities)
LEVEL12_WIDTH = 2**35              # source_id // 2**35 = HEALPix level-12 nested index
N_THREADS = 4                      # concurrent archive queries (be polite)

_manifest_lock = threading.Lock()
_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(msg, flush=True)


def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def record(path: Path, source: str, **extra) -> None:
    """Add a file to data/manifest.json."""
    with _manifest_lock:
        manifest = json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {}
        manifest[str(path.relative_to(C.DATA_DIR))] = {
            "source": source, "bytes": path.stat().st_size, "sha256": sha256(path),
            "retrieved": now_utc(), **extra}
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))


# ---------------------------------------------------------------------------
# Plain HTTP downloads
# ---------------------------------------------------------------------------

def fetch_file(url: str, dest: Path, *, timeout: float = 120.0, retries: int = 4) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = resp.read()
            dest.parent.mkdir(parents=True, exist_ok=True)
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(payload)
            tmp.replace(dest)
            record(dest, url)
            log(f"  {dest.name}: {len(payload) / 1e6:.2f} MB from {url}")
            return dest
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404, 410):
                raise RuntimeError(f"HTTP {exc.code} for {url}") from exc
            last = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
        time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


# ---------------------------------------------------------------------------
# Gaia archive: chunked synchronous TAP queries
# ---------------------------------------------------------------------------

class QueryError(RuntimeError):
    """The archive rejected the ADQL itself: retrying cannot help."""


class TransientError(RuntimeError):
    """Timeout, overload or network trouble: retry, or split the sky patch."""


def _error_text(body: bytes) -> str:
    text = body.decode("utf-8", "replace")
    start = text.find('value="ERROR">')
    if start >= 0:
        text = text[start + 14:]
        text = text[: text.find("</INFO>")]
    return " ".join(text.split())[:500]


def sync_query(adql: str, *, timeout: float = 180.0, attempts: int = 3) -> Table:
    """One synchronous ADQL query -> astropy Table (FITS transport, lower-case columns)."""
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "fits",
                                   "MAXREC": str(MAXREC), "QUERY": adql}).encode()
    last = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(C.GAIA_SYNC_URL, data=data, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                body = resp.read()
            break
        except urllib.error.HTTPError as exc:
            body = exc.read()
            if exc.code == 400:
                raise QueryError(_error_text(body)) from None
            if exc.code == 408:          # ran past the endpoint's ~60 s limit: split instead
                raise TransientError("HTTP 408 (query time limit)") from None
            last = f"HTTP {exc.code}: {_error_text(body)}"
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = f"{type(exc).__name__}: {exc}"
        if attempt < attempts:
            time.sleep(5 * attempt)
    else:
        raise TransientError(str(last))
    if not body.startswith(b"SIMPLE"):
        raise TransientError(f"not a FITS response: {_error_text(body)}")
    with fits.open(io.BytesIO(body)) as hdul:
        table = Table(hdul[1].data)
    for name in table.colnames:
        if name != name.lower():
            table.rename_column(name, name.lower())
    if len(table) >= MAXREC:
        raise TransientError("result hit the 3,000,000-row cap (truncated)")
    return table


def pixel_range(level: int, pixel: int) -> tuple[int, int]:
    width = LEVEL12_WIDTH * 4 ** (12 - level)
    return pixel * width, (pixel + 1) * width


def _chunk(name: str, template: str, level: int, pixel: int) -> list[Path]:
    """Fetch one sky patch (cached), splitting it in four on a timeout."""
    path = C.CHUNK_DIR / name / f"L{level}_{pixel:06d}.fits"
    if path.exists():
        return [path]
    lo, hi = pixel_range(level, pixel)
    t0 = time.time()
    try:
        table = sync_query(template.format(lo=lo, hi=hi))
    except TransientError as exc:
        if level >= C.QUERY_MAX_LEVEL:
            raise
        log(f"    {name} L{level} {pixel}: {exc} -> four smaller patches")
        out: list[Path] = []
        for child in range(4 * pixel, 4 * pixel + 4):
            out += _chunk(name, template, level + 1, child)
        return out
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    table.write(tmp, format="fits", overwrite=True)
    tmp.replace(path)
    log(f"    {name} L{level} {pixel:>5}: {len(table):>8,} rows in {time.time() - t0:5.1f} s")
    return [path]


def _cleanup(name: str) -> None:
    """Delete the per-patch chunk files once the assembled, verified file is written."""
    d = C.CHUNK_DIR / name
    if d.exists():
        for f in d.glob("*.fits"):
            f.unlink()
        d.rmdir()


def all_sky(name: str, template: str, base_level: int) -> list[Path]:
    pixels = range(12 * 4**base_level)
    paths: list[Path] = []
    with cf.ThreadPoolExecutor(N_THREADS) as pool:
        for result in pool.map(lambda p: _chunk(name, template, base_level, p), pixels):
            paths += result
    return paths


# ---------------------------------------------------------------------------
# 1. Point stars
# ---------------------------------------------------------------------------

STARS_PATH = C.GAIA_DIR / "gaia_dr3_g_lt_11.fits"
STARS_ADQL = (
    "SELECT source_id, ra, dec, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag "
    "FROM gaiadr3.gaia_source WHERE source_id >= {lo} AND source_id < {hi} "
    f"AND phot_g_mean_mag < {C.G_POINT_LIMIT!r}")


def fetch_stars() -> None:
    if STARS_PATH.exists():
        log(f"  {STARS_PATH.name} present")
        return
    log(f"== Gaia DR3 stars with G < {C.G_POINT_LIMIT}")
    paths = all_sky("stars", STARS_ADQL, C.QUERY_BASE_LEVEL)
    table = vstack([Table.read(p) for p in paths], metadata_conflicts="silent")
    expected = int(np.asarray(sync_query(
        f"SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE phot_g_mean_mag < {C.G_POINT_LIMIT!r}")["n"])[0])
    if len(table) != expected or len(np.unique(table["source_id"])) != expected:
        raise RuntimeError(f"star list has {len(table)} rows, the archive counts {expected}")
    table.meta["QUERY"] = STARS_ADQL[:60]
    table.write(STARS_PATH, overwrite=True)
    _cleanup("stars")
    record(STARS_PATH, "ESA Gaia archive, gaiadr3.gaia_source", query=STARS_ADQL, rows=len(table),
           server_count=expected)
    log(f"  {len(table):,} stars (archive count {expected:,})")


# ---------------------------------------------------------------------------
# 2. Summed-flux map of all fainter stars
# ---------------------------------------------------------------------------

FLUX_PATH = C.GAIA_DIR / f"gaia_dr3_flux_map_L{C.FLUX_MAP_LEVEL}.npz"
_DIV = LEVEL12_WIDTH * 4 ** (12 - C.FLUX_MAP_LEVEL)
FLUX_ADQL = (
    f"SELECT CAST(source_id/{_DIV} AS INTEGER) AS hpx, CAST(COUNT(*) AS INTEGER) AS n, "
    "CAST(COUNT(phot_bp_mean_flux) AS INTEGER) AS n_bp, "
    "CAST(SUM(phot_g_mean_flux) AS REAL) AS sum_g, "
    "CAST(SUM(phot_bp_mean_flux) AS REAL) AS sum_bp, "
    "CAST(SUM(phot_rp_mean_flux) AS REAL) AS sum_rp "
    "FROM gaiadr3.gaia_source WHERE source_id >= {lo} AND source_id < {hi} "
    f"AND phot_g_mean_mag >= {C.G_POINT_LIMIT!r} GROUP BY hpx")


def fetch_flux_map() -> None:
    if FLUX_PATH.exists():
        log(f"  {FLUX_PATH.name} present")
        return
    log(f"== Gaia DR3 summed fluxes of G >= {C.G_POINT_LIMIT} stars per level-{C.FLUX_MAP_LEVEL} pixel")
    paths = all_sky("flux", FLUX_ADQL, 2)
    npix = 12 * 4**C.FLUX_MAP_LEVEL
    n = np.zeros(npix, np.int64)
    n_bp = np.zeros(npix, np.int64)
    sums = np.zeros((3, npix), np.float64)
    for p in paths:
        t = Table.read(p)
        hpx = np.asarray(t["hpx"], np.int64)
        if np.any(n[hpx] > 0):
            raise RuntimeError(f"{p.name}: pixel seen twice")
        n[hpx] = t["n"]
        n_bp[hpx] = t["n_bp"]
        for i, col in enumerate(("sum_g", "sum_bp", "sum_rp")):
            sums[i, hpx] = np.ma.filled(t[col], 0)
    # Independent check: server-side counts for a few patches, including the densest.
    checks = {}
    for level, pixel in ((2, 112), (2, 0), (3, 450), (4, 1801)):
        lo, hi = pixel_range(level, pixel)
        try:
            cnt = int(np.asarray(sync_query(
                f"SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE source_id >= {lo} AND source_id < {hi} "
                f"AND phot_g_mean_mag >= {C.G_POINT_LIMIT!r}")["n"])[0])
        except TransientError as exc:
            log(f"  count check L{level} {pixel} skipped: {exc}")
            continue
        k = 4 ** (C.FLUX_MAP_LEVEL - level)
        ours = int(n[pixel * k:(pixel + 1) * k].sum())
        checks[f"L{level}_{pixel}"] = (ours, cnt)
        if ours != cnt:
            raise RuntimeError(f"count mismatch in L{level} pixel {pixel}: {ours} vs {cnt}")
    np.savez_compressed(FLUX_PATH, n=n.astype(np.int32), n_bp=n_bp.astype(np.int32),
                        sum_g=sums[0].astype(np.float32), sum_bp=sums[1].astype(np.float32),
                        sum_rp=sums[2].astype(np.float32), level=C.FLUX_MAP_LEVEL)
    record(FLUX_PATH, "ESA Gaia archive, gaiadr3.gaia_source (server-side aggregation)", query=FLUX_ADQL,
           sources=int(n.sum()), empty_pixels=int((n == 0).sum()), count_checks=checks)
    log(f"  {n.sum():,} sources in {npix:,} pixels ({(n == 0).sum():,} empty); checks {checks}")
    _cleanup("flux")


# ---------------------------------------------------------------------------
# 3. Hipparcos + Gaia's cross-match
# ---------------------------------------------------------------------------

HIP_PATH = C.GAIA_DIR / "hipparcos_main.fits"
HIPXM_PATH = C.GAIA_DIR / "gaiadr3_hipparcos2_best_neighbour.fits"


def fetch_hipparcos() -> None:
    if not HIP_PATH.exists():
        log("== Hipparcos main catalogue (ESA 1997, via the Gaia archive)")
        parts = [sync_query("SELECT hip, ra, de, vmag, b_v, hpmag, plx, sptype FROM public.hipparcos "
                            f"WHERE hip >= {lo} AND hip < {lo + 30000}") for lo in range(0, 120001, 30000)]
        t = vstack(parts)
        if len(t) != 118218:
            raise RuntimeError(f"Hipparcos: {len(t)} rows, expected 118,218")
        t["sptype"] = [str(s) for s in np.ma.filled(t["sptype"], "")]
        t.write(HIP_PATH, overwrite=True)
        record(HIP_PATH, "ESA Gaia archive, public.hipparcos (ESA 1997, I/239)", rows=len(t))
        log(f"  {len(t):,} Hipparcos stars")
    if not HIPXM_PATH.exists():
        t = sync_query("SELECT source_id, original_ext_source_id, angular_distance, number_of_neighbours, "
                       "xm_flag FROM gaiadr3.hipparcos2_best_neighbour")
        t.write(HIPXM_PATH, overwrite=True)
        record(HIPXM_PATH, "ESA Gaia archive, gaiadr3.hipparcos2_best_neighbour", rows=len(t))
        log(f"  {len(t):,} Hipparcos-2 / Gaia DR3 best neighbours")


# ---------------------------------------------------------------------------
# 4. Reference tables, 5. literature
# ---------------------------------------------------------------------------

CMF_PATH = C.REF_DIR / "ciexyz31_1.csv"
MAMAJEK_PATH = C.REF_DIR / "EEM_dwarf_UBVIJHK_colors_Teff.txt"


def fetch_reference() -> None:
    log("== reference tables")
    fetch_file("http://cvrl.ucl.ac.uk/database/data/cmfs/ciexyz31_1.csv", CMF_PATH)
    fetch_file("https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt", MAMAJEK_PATH)


LITERATURE = {
    "hamilton_polhemus_2010.pdf": "1012.4043",     # seeing deep inside a black hole
    "hamilton_lisle_2008.pdf": "gr-qc/0411060",    # the river model of black holes
    "lvk_gw250114_2025.pdf": "2509.08054",         # GW250114 area law
    "zhao_2021_cygx1_spin.pdf": "2102.09093",      # Cyg X-1 spin
}


def fetch_literature() -> None:
    log("== literature (for checking numbers; not used by the analysis)")
    for name, arxiv_id in LITERATURE.items():
        try:
            fetch_file(f"https://arxiv.org/pdf/{arxiv_id}", C.RAW_DIR / "literature" / name, timeout=90)
        except Exception as exc:     # reference only: never fatal
            log(f"  could not fetch {arxiv_id}: {exc}")


SOURCES = {"reference": fetch_reference, "hipparcos": fetch_hipparcos, "stars": fetch_stars,
           "flux": fetch_flux_map, "literature": fetch_literature}


def main(which: list[str]) -> None:
    t0 = time.time()
    for key in which or list(SOURCES):
        SOURCES[key]()
    log(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main(sys.argv[1:])
