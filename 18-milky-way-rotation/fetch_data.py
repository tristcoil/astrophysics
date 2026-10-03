"""Download and cache everything the analysis needs.

1. **The HI4PI 21-cm survey** (HI4PI Collaboration 2016), from the CDS archive
   in Strasbourg: the brightness of hydrogen's 21-cm line at every sky position
   and velocity. We take the 18 survey cubes that cover the Galactic plane
   (|b| <= 11 deg, all longitudes) plus 7 just north of it (b = +9 to +31 deg,
   l = 20 to 160, where the outer disc warps upward), and from each only the
   velocity channels with |v| <= 350 km/s. A cube stores one full sky image per velocity channel,
   one after another, so a velocity range is one contiguous stretch of bytes
   and an HTTP range request fetches exactly it: 3.8 GB instead of 6.6 GB.
   Each subset is saved as a standard FITS file whose header is the original
   with the velocity axis re-indexed, plus HISTORY lines saying which bytes of
   which file it holds.

2. **Small reference tables**: published terminal velocities
   (McClure-Griffiths & Dickey 2016), maser parallaxes (Reid et al. 2019), and
   the SPARC rotation curves of 175 other galaxies (Lelli et al. 2016). The
   analysis measures everything from HI4PI; these are only for comparison.

Everything lands in ``data/`` with ``manifest.json``, which records each file's
source URL, byte range, retrieval time and SHA-256 hash. Interrupted downloads
resume where they stopped; finished files are never fetched again.

    python fetch_data.py            # everything
    python fetch_data.py --tiles 3  # only the first 3 survey cubes (for testing)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits

from config import (DATA_DIR, HI4PI_DIR, HI4PI_PAPER, HI4PI_ROW_E, HI4PI_ROW_F_WARP,
                    HI4PI_URL,
                    MANIFEST_PATH, SOURCES, V_KEEP_KMS, Source)

USER_AGENT = "astrophysics-project-atlas/H2 (hobby research; python-urllib)"
CHUNK_BYTES = 8 * 1024 * 1024
PARALLEL = 3
FITS_BLOCK = 2880

_manifest_lock = threading.Lock()
_print_lock = threading.Lock()


def _say(*args, **kwargs) -> None:
    with _print_lock:
        print(*args, **kwargs, flush=True)


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------

def _get(url: str, byte_range: tuple[int, int] | None = None, *,
         attempts: int = 6, timeout: float = 120.0) -> bytes:
    """GET a URL (optionally a byte range, inclusive), retrying with backoff."""
    headers = {"User-Agent": USER_AGENT}
    if byte_range is not None:
        headers["Range"] = f"bytes={byte_range[0]}-{byte_range[1]}"
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                        timeout=timeout) as resp:
                payload = resp.read()
                if byte_range is not None:
                    want = byte_range[1] - byte_range[0] + 1
                    if resp.status != 206 or len(payload) != want:
                        raise ValueError(f"asked for {want} bytes, got status {resp.status}, "
                                         f"{len(payload)} bytes")
            return payload
        except (urllib.error.URLError, TimeoutError, ConnectionError, ValueError) as err:
            last = err
            if attempt < attempts:
                time.sleep(min(60.0, 3.0 * attempt ** 2))
    raise RuntimeError(f"could not download {url} {byte_range or ''}: {last}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(name: str, entry: dict) -> None:
    with _manifest_lock:
        manifest = _load_manifest()
        manifest["files"][name] = entry
        MANIFEST_PATH.write_text(json.dumps(manifest, indent=2))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


# --------------------------------------------------------------------------
# Small tables
# --------------------------------------------------------------------------

def fetch_tables(*, quiet: bool = False) -> None:
    for source in SOURCES:
        path = DATA_DIR / source.filename
        if path.exists():
            if not quiet:
                _say(f"  cached     {source.filename}")
            continue
        payload = _get(source.url)
        if payload[:200].lstrip().lower().startswith((b"<!doctype", b"<html")):
            raise RuntimeError(f"{source.url} returned an HTML page, not the table")
        tmp = path.with_suffix(path.suffix + ".part")
        tmp.write_bytes(payload)
        tmp.replace(path)
        _record(source.filename, {
            "url": source.url, "retrieved_utc": _now(), "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(), "what": source.what,
            "reference": source.reference,
            "used_as": "analysis input" if source.input else "published comparison only",
        })
        if not quiet:
            _say(f"  fetched    {source.filename} ({len(payload) / 1024:.0f} kB)")


EILERS19_ARXIV = "https://arxiv.org/e-print/1810.09466"
EILERS19_FILE = "eilers19_vc.txt"


def fetch_eilers19(*, quiet: bool = False) -> None:
    """Eilers et al. (2019) Table 1, taken from the paper's arXiv source.

    The circular-velocity table is not on VizieR; the authors' LaTeX source
    holds it as a deluxetable (its first nine rows, inside 5 kpc, are
    commented out in the source and were not published; we skip them too).
    """
    import io
    import re
    import tarfile

    path = DATA_DIR / EILERS19_FILE
    if path.exists():
        if not quiet:
            _say(f"  cached     {EILERS19_FILE}")
        return
    payload = _get(EILERS19_ARXIV)
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:*") as tar:
        tex = tar.extractfile("main.tex").read().decode("utf-8")
    block = tex[tex.index("label{tab:measurements}"):]
    block = block[block.index(r"\startdata") + len(r"\startdata"):block.index(r"\enddata")]
    rows = []
    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("%"):
            continue
        nums = [float(x) for x in re.findall(r"[-+]?\d+\.?\d*", line)]
        if len(nums) != 4:
            raise RuntimeError(f"unexpected Eilers et al. table row: {line!r}")
        rows.append(nums)
    header = ("# Eilers, Hogg, Rix & Ness 2019, ApJ 871, 120, Table 1 (arXiv:1810.09466)\n"
              "# circular velocity of the Milky Way from red giants; R0 = 8.122 kpc\n"
              "# R_kpc  vc_kms  err_minus_kms  err_plus_kms\n")
    path.write_text(header + "".join(f"{r[0]:.2f} {r[1]:.2f} {r[2]:.2f} {r[3]:.2f}\n"
                                     for r in rows))
    _record(EILERS19_FILE, {
        "url": EILERS19_ARXIV, "retrieved_utc": _now(),
        "source_sha256": hashlib.sha256(payload).hexdigest(), "rows": len(rows),
        "what": "Milky Way circular velocity 5-25 kpc from 23,000 red giants "
                "(Gaia DR2 + APOGEE + WISE), extracted from the LaTeX table",
        "reference": "Eilers et al. 2019, ApJ 871, 120", "used_as": "published comparison only",
    })
    if not quiet:
        _say(f"  fetched    {EILERS19_FILE} ({len(rows)} rows from the arXiv source)")


# --------------------------------------------------------------------------
# HI4PI cubes, velocity subset by byte range
# --------------------------------------------------------------------------

def tile_path(name: str) -> Path:
    return HI4PI_DIR / name.replace(".fits", f"_v{int(V_KEEP_KMS)}.fits")


def _plan_tile(name: str) -> dict:
    """Read a cube's header from the archive and work out which bytes to fetch."""
    url = HI4PI_URL + "CUBES/GAL/CAR/" + name
    head = _get(url, (0, 4 * FITS_BLOCK - 1))
    cards = [head[i:i + 80] for i in range(0, len(head), 80)]
    n_end = next(i for i, card in enumerate(cards) if card.startswith(b"END "))
    header = fits.Header.fromstring(b"".join(cards[:n_end + 1]).decode("ascii"))
    header_bytes = -(-(n_end + 1) * 80 // FITS_BLOCK) * FITS_BLOCK
    nx, ny, nv = header["NAXIS1"], header["NAXIS2"], header["NAXIS3"]
    if header["BITPIX"] != -32 or header["CTYPE3"].strip() != "VRAD":
        raise RuntimeError(f"{name}: unexpected cube layout")
    k = np.arange(1, nv + 1)
    v_kms = (header["CRVAL3"] + (k - header["CRPIX3"]) * header["CDELT3"]) / 1000.0
    keep = np.flatnonzero(np.abs(v_kms) <= V_KEEP_KMS) + 1  # 1-based channels
    k_lo, k_hi = int(keep[0]), int(keep[-1])
    plane = nx * ny * 4
    start = header_bytes + (k_lo - 1) * plane
    stop = header_bytes + k_hi * plane - 1

    sub = header.copy()
    for key in ("CHECKSUM", "DATASUM"):
        sub.remove(key, ignore_missing=True)
    sub["NAXIS3"] = k_hi - k_lo + 1
    sub["CRPIX3"] = header["CRPIX3"] - (k_lo - 1)
    sub.add_history(f"Velocity subset |v| <= {V_KEEP_KMS:g} km/s of HI4PI {name}: "
                    f"channels {k_lo}-{k_hi} of {nv}")
    sub.add_history(f"bytes {start}-{stop} of {url}")
    sub.add_history(f"original DATASUM {header.get('DATASUM', '?')} (whole cube)")
    return {"name": name, "url": url, "header": sub, "start": start, "stop": stop,
            "channels": (k_lo, k_hi), "total_channels": nv}


def fetch_tile(name: str) -> Path:
    """Fetch one cube's velocity subset, resuming a partial download."""
    final = tile_path(name)
    if final.exists():
        return final
    plan = _plan_tile(name)
    header_block = plan["header"].tostring().encode("ascii")
    assert len(header_block) % FITS_BLOCK == 0
    nbytes = plan["stop"] - plan["start"] + 1

    part = final.with_suffix(".fits.part")
    have = 0
    if part.exists():
        have = max(0, part.stat().st_size - len(header_block))
        have -= have % CHUNK_BYTES  # restart the last, possibly short, chunk
    mode = "r+b" if part.exists() else "wb"
    t0 = time.time()
    with part.open(mode) as fh:
        fh.seek(0)
        fh.write(header_block)
        fh.seek(len(header_block) + have)
        fh.truncate()
        offset = have
        while offset < nbytes:
            n = min(CHUNK_BYTES, nbytes - offset)
            a = plan["start"] + offset
            fh.write(_get(plan["url"], (a, a + n - 1)))
            offset += n
        pad = (-(len(header_block) + nbytes)) % FITS_BLOCK
        fh.write(b"\0" * pad)

    # Sanity check: the file opens as FITS with the promised shape, and holds
    # plausible brightness temperatures (K), not garbage.
    with fits.open(part, memmap=True) as hdul:
        data = hdul[0].data
        if data.shape != (plan["header"]["NAXIS3"], plan["header"]["NAXIS2"],
                          plan["header"]["NAXIS1"]):
            raise RuntimeError(f"{name}: wrong shape {data.shape}")
        mid = np.asarray(data[data.shape[0] // 2])
        finite = mid[np.isfinite(mid)]
        if finite.size < 0.5 * mid.size or not (-5 < np.median(finite) < 300):
            raise RuntimeError(f"{name}: implausible values in the middle channel")
    part.replace(final)

    _record(f"hi4pi/{final.name}", {
        "url": plan["url"], "byte_range": [plan["start"], plan["stop"]],
        "channels": list(plan["channels"]), "of_channels": plan["total_channels"],
        "retrieved_utc": _now(), "bytes": final.stat().st_size, "sha256": _sha256(final),
        "what": f"HI4PI 21-cm cube {name}, velocity subset |v_LSR| <= {V_KEEP_KMS:g} km/s",
        "reference": HI4PI_PAPER, "used_as": "analysis input",
    })
    rate = nbytes / max(time.time() - t0, 1e-3) / 1e6
    _say(f"  fetched    {final.name} ({nbytes / 1e6:.0f} MB at {rate:.1f} MB/s)")
    return final


def fetch_tiles(names: list[str], *, parallel: int = PARALLEL, quiet: bool = False) -> list[Path]:
    todo = [n for n in names if not tile_path(n).exists()]
    for n in names:
        if n not in todo and not quiet:
            _say(f"  cached     {tile_path(n).name}")
    if todo:
        _say(f"  {len(todo)} survey cubes to fetch, {parallel} at a time "
             f"(~154 MB each)")
    with ThreadPoolExecutor(max_workers=parallel) as pool:
        return list(pool.map(fetch_tile, names))


ALL_TILES = HI4PI_ROW_E + HI4PI_ROW_F_WARP


def missing_tiles(names: list[str] = ALL_TILES) -> list[str]:
    return [n for n in names if not tile_path(n).exists()]


def fetch_all(*, tiles: int | None = None, quiet: bool = False) -> None:
    fetch_tables(quiet=quiet)
    fetch_eilers19(quiet=quiet)
    fetch_tiles(ALL_TILES if tiles is None else ALL_TILES[:tiles], quiet=quiet)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--tiles", type=int, default=None,
                        help="fetch only the first N survey cubes")
    args = parser.parse_args(argv)
    print("Fetching reference tables and the HI4PI Galactic-plane cubes")
    try:
        fetch_all(tiles=args.tiles)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1
    print("\nDone; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
