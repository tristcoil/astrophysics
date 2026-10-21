"""Download and cache every input of the analysis.

Run directly to populate the cache::

    python fetch_data.py              # everything (~0.4 GB new, plus 1.6 GB of Planck maps)
    python fetch_data.py --only nvss  # one product

What comes from where:

* **Planck PR3 SMICA maps, masks and best-fit spectra** -- the ESA Planck Legacy
  Archive (IRSA mirror as fallback). Project D3 already downloaded the same
  files; if its cache exists they are hard-linked, costing no extra disk.
* **WMAP 9-year ILC map** -- NASA's LAMBDA archive.
* **Granett et al. (2008) supervoids and superclusters** -- their web page is
  gone, but the Internet Archive kept both text files.
* **Mao et al. (2017) BOSS voids; Nadathur & Hotchkiss (2014) DR7 voids and
  superclusters; the NVSS radio catalogue** -- the CDS archive in Strasbourg.
* **BOSS DR12 galaxies and survey footprints** -- the SDSS Science Archive Server.

Every file is written under a temporary name and renamed only when complete
(and, where the size is known, only when the size matches), so an interrupted
run never leaves a truncated file that looks good. Nothing cached is fetched
again; once the cache is warm the whole analysis runs offline.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import requests

from config import D3_PLANCK_DIR, PRODUCTS, Product, product

CHUNK = 1 << 20
TIMEOUT = (20, 120)


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.0f} B" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GB"


def _stream(url: str, dest: Path, expected: int | None, agent: str, quiet: bool) -> None:
    """Stream ``url`` to ``dest``, resuming a partial download when the server allows."""
    part = dest.with_name(dest.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    if expected is not None and have > expected:
        part.unlink()
        have = 0
    headers = {"User-Agent": agent}
    if have:
        headers["Range"] = f"bytes={have}-"
    with requests.get(url, headers=headers, stream=True, timeout=TIMEOUT,
                      allow_redirects=True) as r:
        if r.status_code == 206 and have:
            mode = "ab"
        elif r.status_code == 200:
            mode, have = "wb", 0
        else:
            raise RuntimeError(f"HTTP {r.status_code}")
        done, t0, next_report = have, time.monotonic(), 0.0
        with open(part, mode) as fh:
            for chunk in r.iter_content(CHUNK):
                fh.write(chunk)
                done += len(chunk)
                if (not quiet and expected and expected > 20 * CHUNK
                        and done / expected >= next_report):
                    rate = (done - have) / max(time.monotonic() - t0, 1e-6)
                    print(f"    {done / expected:6.1%}  {_human(done):>9} / {_human(expected)}"
                          f"  ({_human(rate)}/s)", flush=True)
                    next_report += 0.1
    size = part.stat().st_size
    if size == 0 or (expected is not None and size != expected):
        if expected is not None and size > expected:
            part.unlink()
        raise RuntimeError(f"got {size:,} bytes, expected {expected if expected else 'more'}")
    part.rename(dest)


def _link_from_d3(p: Product, quiet: bool) -> bool:
    """Hard-link one of D3's identical Planck files, if present and complete."""
    src = D3_PLANCK_DIR / p.filename
    if p.subdir != "planck" or not src.exists():
        return False
    if p.size_bytes is not None and src.stat().st_size != p.size_bytes:
        return False
    try:
        os.link(src, p.path)
    except OSError:
        return False
    if not quiet:
        print(f"  linked   {p.filename}  (from D3's cache, no download)")
    return True


def fetch_product(p: Product, *, refresh: bool = False, quiet: bool = False) -> Path:
    dest = p.path
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not refresh:
        if p.size_bytes is None or dest.stat().st_size == p.size_bytes:
            if not quiet:
                print(f"  cached   {p.filename}")
            return dest
        print(f"  {p.filename} has the wrong size -- fetching again")
    if dest.exists():
        dest.unlink()
    if _link_from_d3(p, quiet):
        return dest
    if not quiet:
        size = _human(p.size_bytes) if p.size_bytes else "size unknown"
        print(f"  fetching {p.filename}  ({size})")
    errors = []
    for url in p.urls:
        host = url.split("/")[2]
        for attempt in range(3):
            try:
                _stream(url, dest, p.size_bytes, p.agent, quiet)
                if not quiet:
                    print(f"    ok, {_human(dest.stat().st_size)} from {host}")
                return dest
            except Exception as exc:  # network errors, wrong sizes, HTTP errors
                errors.append(f"{host} (attempt {attempt + 1}): {exc}")
                if not quiet:
                    print(f"    {host} failed ({exc})")
                time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"could not download {p.filename}:\n  " + "\n  ".join(errors))


def fetch_all(*, refresh: bool = False, quiet: bool = False) -> None:
    for p in PRODUCTS:
        fetch_product(p, refresh=refresh, quiet=quiet)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", action="append", metavar="KEY",
                        help=f"fetch only this product (repeatable): "
                             f"{', '.join(p.key for p in PRODUCTS)}")
    parser.add_argument("--refresh", action="store_true", help="download again even if cached")
    args = parser.parse_args(argv)
    if args.only:
        for key in args.only:
            fetch_product(product(key), refresh=args.refresh)
    else:
        fetch_all(refresh=args.refresh)
        print("\nDone. The cache is warm; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
