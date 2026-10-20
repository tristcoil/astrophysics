"""Read chosen columns of FITS binary tables from a gzipped file in one streaming pass.

astropy cannot memory-map a gzipped FITS file, so opening the 1.4 GB Wu & Shen
catalogue (2.6 GB of table once inflated) would pull the whole thing into RAM.
This reader inflates the file once, block by block, parses each header with
astropy, turns the TFORMn keywords into a numpy record type, and keeps only the
requested columns. Peak memory is one chunk of rows plus the output columns.

Supported column types are the fixed-width ones (L, B, I, J, K, E, D, A, with
repeat counts and TDIM shapes, TZERO offsets for unsigned integers). Variable-
length (P/Q) columns are rejected: the catalogue has none.
"""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
from astropy.io import fits

BLOCK = 2880
_CARD_END = b"END" + b" " * 77

_TFORM_CODES = {
    "L": "S1",   # logical 'T'/'F'
    "B": "u1",
    "I": ">i2",
    "J": ">i4",
    "K": ">i8",
    "E": ">f4",
    "D": ">f8",
    "C": ">c8",
    "M": ">c16",
}


def _read_header(stream) -> fits.Header | None:
    """Read 2880-byte blocks up to and including the one holding END."""
    chunks = []
    while True:
        block = stream.read(BLOCK)
        if not block:
            return None if not chunks else _fail("file ends inside a header")
        if len(block) != BLOCK:
            _fail("truncated header block")
        chunks.append(block)
        cards = [block[i:i + 80] for i in range(0, BLOCK, 80)]
        if _CARD_END in cards:
            break
    return fits.Header.fromstring(b"".join(chunks).decode("ascii"))


def _fail(message: str):
    raise ValueError(f"FITS stream: {message}")


def row_dtype(header: fits.Header) -> np.dtype:
    """numpy record type of one table row, from the TFORMn/TDIMn keywords."""
    names, formats, shapes = [], [], []
    for i in range(1, header["TFIELDS"] + 1):
        name = header[f"TTYPE{i}"].strip()
        tform = header[f"TFORM{i}"].strip()
        m = re.fullmatch(r"(\d*)([A-Z])(.*)", tform)
        if not m:
            _fail(f"cannot parse TFORM{i} = {tform!r}")
        repeat = int(m.group(1)) if m.group(1) else 1
        code = m.group(2)
        if code in "PQ":
            _fail(f"variable-length column {name} is not supported")
        if code == "A":
            fmt, shape = f"S{repeat}", ()
        elif code == "X":
            fmt, shape = "u1", ((repeat + 7) // 8,)
        else:
            fmt = _TFORM_CODES[code]
            tdim = header.get(f"TDIM{i}")
            if tdim:
                dims = tuple(int(v) for v in tdim.strip("() ").split(","))[::-1]
                shape = dims
            else:
                shape = (repeat,) if repeat != 1 else ()
        names.append(name)
        formats.append(fmt)
        shapes.append(shape)
    dtype = np.dtype([(n, f, s) if s else (n, f) for n, f, s in zip(names, formats, shapes)])
    if dtype.itemsize != header["NAXIS1"]:
        _fail(f"row size {dtype.itemsize} != NAXIS1 {header['NAXIS1']}")
    return dtype


def _finish(header: fits.Header, name: str, data: np.ndarray) -> np.ndarray:
    """Apply TZERO/TSCAL, make native-endian, decode strings."""
    idx = [header[f"TTYPE{i}"].strip() for i in range(1, header["TFIELDS"] + 1)].index(name) + 1
    zero = header.get(f"TZERO{idx}", 0)
    scale = header.get(f"TSCAL{idx}", 1)
    if data.dtype.kind == "S":
        return np.char.strip(np.char.decode(data, "ascii"))
    out = data.astype(data.dtype.newbyteorder("="), copy=False)
    if zero or scale != 1:
        bits = 8 * out.dtype.itemsize
        if out.dtype.kind == "i" and scale == 1 and zero == 2 ** (bits - 1):
            # FITS stores unsigned integers as signed ones offset by TZERO = 2^(n-1);
            # adding that offset is the same as flipping the sign bit.
            unsigned = np.dtype(f"u{out.dtype.itemsize}")
            out = out.view(unsigned) ^ np.array(1 << (bits - 1), dtype=unsigned)
        else:
            out = out * scale + zero
    return out


def read_columns(path: str | Path, wanted: dict[int, list[str]], *, chunk_rows: int = 20000,
                 progress=None) -> dict[int, dict[str, np.ndarray]]:
    """Read ``wanted[hdu_index]`` columns from each binary-table HDU of a gzipped FITS file.

    HDU indices count from 0 (the primary HDU). Returns {hdu: {column: array}}.
    ``progress`` is called with (hdu, rows_done, rows_total) after each chunk.
    """
    out: dict[int, dict[str, np.ndarray]] = {}
    with gzip.open(path, "rb") as stream:
        hdu = 0
        while True:
            header = _read_header(stream)
            if header is None:
                break
            naxis = header.get("NAXIS", 0)
            if hdu == 0 or header.get("XTENSION", "").strip() != "BINTABLE":
                size = 0
                if naxis:
                    size = abs(header["BITPIX"]) // 8
                    for k in range(1, naxis + 1):
                        size *= header[f"NAXIS{k}"]
                    size += header.get("PCOUNT", 0)
                _skip(stream, size)
                _skip_padding(stream, size)
                hdu += 1
                continue
            nrows, rowlen = header["NAXIS2"], header["NAXIS1"]
            heap = header.get("PCOUNT", 0)
            dtype = row_dtype(header)
            cols = wanted.get(hdu, [])
            missing = [c for c in cols if c not in dtype.names]
            if missing:
                _fail(f"HDU {hdu} lacks columns {missing}")
            if cols:
                parts = {c: [] for c in cols}
                done = 0
                while done < nrows:
                    n = min(chunk_rows, nrows - done)
                    buf = _read_exact(stream, n * rowlen)
                    rows = np.frombuffer(buf, dtype=dtype, count=n)
                    for c in cols:
                        parts[c].append(rows[c].copy())
                    done += n
                    if progress:
                        progress(hdu, done, nrows)
                out[hdu] = {c: _finish(header, c, np.concatenate(parts[c])) for c in cols}
                _skip(stream, heap)
                _skip_padding(stream, nrows * rowlen + heap)
            else:
                _skip(stream, nrows * rowlen + heap)
                _skip_padding(stream, nrows * rowlen + heap)
            hdu += 1
            if all(h in out for h in wanted):
                break
    return out


def _read_exact(stream, n: int) -> bytes:
    buf = stream.read(n)
    if len(buf) != n:
        _fail(f"expected {n} bytes, got {len(buf)}")
    return buf


def _skip(stream, n: int) -> None:
    while n > 0:
        step = min(n, 1 << 24)
        got = stream.read(step)
        if not got:
            _fail("file ends inside a data unit")
        n -= len(got)


def _skip_padding(stream, nbytes: int) -> None:
    pad = (-nbytes) % BLOCK
    if pad:
        _skip(stream, pad)
