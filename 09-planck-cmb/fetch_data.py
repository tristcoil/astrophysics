"""Download and cache everything the analysis needs.

Run directly to populate the cache::

    python fetch_data.py                 # everything, ~1.6 GB
    python fetch_data.py --spectra-only  # skip the sky maps, ~0.4 MB

What comes from where:

* **COBE/FIRAS** spectrum -- NASA's LAMBDA archive.
* **Planck 2018 spectra, maps and masks** -- the ESA Planck Legacy Archive, with
  NASA/IPAC's IRSA mirror as the fallback.
* **The SMICA beam** -- the last few kilobytes of a 2 GB file, fetched alone
  with an HTTP range request (see ``fetch_beam``).

Every file is checked against its known size before it is accepted, and written
under a temporary name until it is complete, so an interrupted run can never
leave a truncated file masquerading as a good one. Nothing already cached is
downloaded again; once the cache is warm the analysis runs offline.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import requests

from config import (BEAM_PATH, BEAM_SOURCE_SIZE, BEAM_SOURCE_URL, PRODUCTS, Product,
                    product)

CHUNK = 1 << 20  # 1 MiB
TIMEOUT = (20, 120)  # connect, read -- seconds
USER_AGENT = "astrophysics-atlas-D3/1.0 (educational; python-requests)"


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024 or unit == "GB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n:.0f} B"
        n /= 1024
    return f"{n:.1f} GB"


def _stream(url: str, dest: Path, expected: int, *, quiet: bool) -> None:
    """Stream ``url`` into ``dest``, resuming a partial file where the server allows.

    The partial file lives at ``dest.part`` and is renamed only once its size is
    exactly ``expected``. A server that ignores the Range header (the PLA does)
    answers 200 with the whole file, in which case we start again from zero.
    """
    part = dest.with_name(dest.name + ".part")
    have = part.stat().st_size if part.exists() else 0
    if have > expected:
        part.unlink()
        have = 0

    headers = {"User-Agent": USER_AGENT}
    if have:
        headers["Range"] = f"bytes={have}-"

    with requests.get(url, headers=headers, stream=True, timeout=TIMEOUT) as r:
        if r.status_code == 206 and have:
            mode = "ab"
            if not quiet:
                print(f"    resuming at {_human(have)}")
        elif r.status_code == 200:
            mode, have = "wb", 0
        else:
            raise RuntimeError(f"HTTP {r.status_code}")

        total, done = expected, have
        next_report, t0 = 0.0, time.monotonic()
        with open(part, mode) as fh:
            for chunk in r.iter_content(CHUNK):
                fh.write(chunk)
                done += len(chunk)
                frac = done / total
                if not quiet and total > 20 * CHUNK and frac >= next_report:
                    rate = (done - have) / max(time.monotonic() - t0, 1e-6)
                    print(f"    {frac:6.1%}  {_human(done):>9} / {_human(total)}"
                          f"  ({_human(rate)}/s)", flush=True)
                    next_report += 0.1

    size = part.stat().st_size
    if size != expected:
        # Keep a short file for resuming; a long one is garbage.
        if size > expected:
            part.unlink()
        raise RuntimeError(f"got {size:,} bytes, expected {expected:,}")
    part.rename(dest)


def fetch_product(p: Product, *, refresh: bool = False, quiet: bool = False) -> Path:
    """Make sure one product is on disk, downloading it if necessary."""
    dest = p.path
    if dest.exists() and not refresh:
        if dest.stat().st_size == p.size_bytes:
            if not quiet:
                print(f"  cached   {p.filename}")
            return dest
        print(f"  {p.filename} has the wrong size -- fetching again")
    if refresh and dest.exists():
        dest.unlink()

    if not quiet:
        print(f"  fetching {p.filename}  ({_human(p.size_bytes)})")
    errors = []
    for url in p.urls:
        host = url.split("/")[2]
        try:
            _stream(url, dest, p.size_bytes, quiet=quiet)
            if not quiet:
                print(f"    ok, from {host}")
            return dest
        except Exception as exc:  # network errors, bad sizes, HTTP errors alike
            errors.append(f"{host}: {exc}")
            if not quiet:
                print(f"    {host} failed ({exc}); trying the next mirror")
    raise RuntimeError(f"could not download {p.filename}:\n  " + "\n  ".join(errors))


# --------------------------------------------------------------------------
# The beam
# --------------------------------------------------------------------------

_BLOCK = 2880  # FITS files are made of 2880-byte blocks


def _minimal_primary_header() -> bytes:
    """A valid empty primary HDU, so a bare extension becomes a readable FITS file."""
    cards = ["SIMPLE  =                    T / conforms to FITS standard",
             "BITPIX  =                    8 / array data type",
             "NAXIS   =                    0 / number of array dimensions",
             "EXTEND  =                    T",
             "COMMENT   Primary HDU added by fetch_data.py; the extension that follows",
             "COMMENT   is byte-for-byte the last HDU of COM_CMB_IQU-smica_2048_R3.00_full",
             "END"]
    raw = "".join(c.ljust(80) for c in cards).encode("ascii")
    return raw + b" " * (-len(raw) % _BLOCK)


def fetch_beam(*, refresh: bool = False, quiet: bool = False) -> Path:
    """Fetch the SMICA effective beam transfer function B_l.

    The beam is not published as a file of its own: it is stored as the final
    HDU of the 2 GB full-mission SMICA map, which we otherwise do not need (we
    use the two half-mission maps instead). So we ask the server for the last
    ~58 KB only, find the start of that final HDU -- a FITS header always begins
    on a 2880-byte boundary with the keyword ``XTENSION`` -- and write it out
    behind a minimal primary header.

    Why bother, when the beam is "a 5 arcmin Gaussian"? Because the details
    matter at the percent level: Planck's B_l departs from a pure Gaussian by up
    to ~4% near l = 3000. (Whether the HEALPix pixel window must be applied on
    top of it is its own story -- see ``sky.read_beam``.)
    """
    from astropy.io import fits

    if BEAM_PATH.exists() and not refresh:
        if not quiet:
            print(f"  cached   {BEAM_PATH.name}")
        return BEAM_PATH

    tail_len = 20 * _BLOCK
    start = BEAM_SOURCE_SIZE - tail_len
    if not quiet:
        print(f"  fetching {BEAM_PATH.name}  (the last {_human(tail_len)} of a "
              f"{_human(BEAM_SOURCE_SIZE)} file, by range request)")
    r = requests.get(BEAM_SOURCE_URL, timeout=TIMEOUT,
                     headers={"Range": f"bytes={start}-{BEAM_SOURCE_SIZE - 1}",
                              "User-Agent": USER_AGENT})
    if r.status_code != 206 or len(r.content) != tail_len:
        raise RuntimeError(f"range request failed: HTTP {r.status_code}, "
                           f"{len(r.content)} bytes")
    tail = r.content

    # The last XTENSION card that sits on a block boundary starts the final HDU.
    offset = None
    for k in range(len(tail) // _BLOCK - 1, -1, -1):
        if (start + k * _BLOCK) % _BLOCK == 0 and tail[k * _BLOCK:].startswith(b"XTENSION"):
            offset = k * _BLOCK
            break
    if offset is None:
        raise RuntimeError("no extension header found in the tail of the file")

    BEAM_PATH.write_bytes(_minimal_primary_header() + tail[offset:])
    with fits.open(BEAM_PATH) as hdul:
        cols = hdul[1].columns.names
        if "INT_BEAM" not in cols:
            BEAM_PATH.unlink()
            raise RuntimeError(f"unexpected beam columns {cols}")
        if not quiet:
            print(f"    ok, {len(hdul[1].data)} multipoles, columns {cols}")
    return BEAM_PATH


# --------------------------------------------------------------------------

def have_maps() -> bool:
    """True when every product needed for the map analysis is cached."""
    maps = [p for p in PRODUCTS if p.group == "maps"]
    return BEAM_PATH.exists() and all(
        p.path.exists() and p.path.stat().st_size == p.size_bytes for p in maps)


def fetch_all(*, include_maps: bool = True, refresh: bool = False, quiet: bool = False):
    groups = ["firas", "spectra"] + (["maps"] if include_maps else [])
    for group in groups:
        if not quiet:
            print(f"\n[{group}]")
        for p in PRODUCTS:
            if p.group == group:
                fetch_product(p, refresh=refresh, quiet=quiet)
        if group == "maps":
            fetch_beam(refresh=refresh, quiet=quiet)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--spectra-only", action="store_true",
                        help="skip the ~1.6 GB of sky maps and masks")
    parser.add_argument("--only", action="append", metavar="KEY",
                        help=f"fetch just this product (repeatable): "
                             f"{', '.join(p.key for p in PRODUCTS)}, beam")
    parser.add_argument("--refresh", action="store_true",
                        help="download again even if cached")
    args = parser.parse_args(argv)

    if args.only:
        for key in args.only:
            if key == "beam":
                fetch_beam(refresh=args.refresh)
            else:
                fetch_product(product(key), refresh=args.refresh)
        return 0

    fetch_all(include_maps=not args.spectra_only, refresh=args.refresh)
    print("\nDone. Cache is warm; the analysis can now run offline.")
    if args.spectra_only:
        print("(Sky maps skipped: the map figures and the map power spectrum will be "
              "left out.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
