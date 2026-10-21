"""Read the Minor Planet Center's MPCORB.DAT into numpy arrays.

MPCORB is a fixed-width text file, one 202-character line per object, format
documented at https://minorplanetcenter.net/iau/info/MPOrbitFormat.html. Two
fields use the MPC's *packed* encodings, which squeeze dates and large numbers
into five characters:

* **Packed dates** (the epoch): ``K2669`` is 2026 June 9. The first character
  is the century (I=18, J=19, K=20), then two year digits, then month and day
  written in base 32 -- 1-9 as digits, then A=10, B=11, ... V=31.

* **Packed numbers**: asteroid numbers above 99,999 do not fit in five
  digits, so the first digit becomes a letter: A0000 = 100,000, Z9999 =
  359,999, a0000 = 360,000, z9999 = 619,999. Beyond that the MPC switched to
  a tilde plus four base-62 digits: ~0000 = 620,000.

The file has three sections separated by blank lines: numbered objects,
unnumbered multi-opposition orbits, and recent one-opposition orbits. The
section matters, because one-opposition orbits are far less certain.

Parsing 1.5 million lines takes ~20 s, so the result is cached as ``.npz``.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np

from config import MPCORB_PATH, MPCORB_PARSED

LINE_WIDTH = 202

#: (name, first column, last column) in the MPC's 1-based inclusive
#: convention, copied from the format page so the two can be compared by eye.
FIELDS = [
    ("packed", 1, 7),
    ("H", 9, 13),
    ("G", 15, 19),
    ("epoch", 21, 25),
    ("M", 27, 35),  # mean anomaly at epoch, deg
    ("peri", 38, 46),  # argument of perihelion, J2000 ecliptic, deg
    ("node", 49, 57),  # longitude of ascending node, J2000 ecliptic, deg
    ("incl", 60, 68),  # inclination to the J2000 ecliptic, deg
    ("e", 71, 79),
    ("n", 81, 91),  # mean daily motion, deg/day
    ("a", 93, 103),  # semi-major axis, AU
    ("U", 106, 106),  # uncertainty parameter 0-9 (or E/D/F flags)
    ("n_obs", 118, 122),
    ("n_opp", 124, 126),
    ("arc", 128, 136),  # "1801-2026" or "  34 days"
    ("rms", 138, 141),
    ("flags", 162, 165),  # 4 hex digits: orbit type + NEO/PHA bits
    ("readable", 167, 194),
    ("last_obs", 195, 202),
]

FLOAT_FIELDS = ("H", "G", "M", "peri", "node", "incl", "e", "n", "a", "rms")

#: Orbit-type codes in the low 6 bits of the hex flags (MPC format page).
ORBIT_TYPES = {
    1: "Atira", 2: "Aten", 3: "Apollo", 4: "Amor", 5: "q < 1.665 AU",
    6: "Hungaria", 8: "Hilda", 9: "Jupiter Trojan", 10: "Distant object",
}
FLAG_NEO = 2048
FLAG_PHA = 32768

_BASE62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


# --------------------------------------------------------------------------
# Packed encodings
# --------------------------------------------------------------------------

def _base32_digit(c: str) -> int:
    """MPC month/day digit: '1'-'9' are 1-9, 'A'-'V' are 10-31."""
    if c.isdigit():
        return int(c)
    if "A" <= c <= "V":
        return ord(c) - ord("A") + 10
    raise ValueError(f"bad packed digit {c!r}")


def unpack_date(packed: str) -> tuple[int, int, float]:
    """``'K2669'`` -> ``(2026, 6, 9.0)``.

    Any characters after the fifth are a decimal fraction of the day, which
    the MPC simply appends: ``'K01AM138303'`` -> 2001 Oct 22.138303.
    """
    packed = packed.strip()
    century = {"I": 18, "J": 19, "K": 20}[packed[0]]
    year = century * 100 + int(packed[1:3])
    month = _base32_digit(packed[3])
    day = float(_base32_digit(packed[4]))
    if len(packed) > 5:
        day += float("0." + packed[5:])
    return year, month, day


def calendar_to_jd(year: int, month: int, day: float) -> float:
    """Julian Date of a Gregorian calendar date (Meeus, *Astronomical
    Algorithms*, ch. 7). Day may be fractional; day 9.0 is 0h."""
    if month <= 2:
        year -= 1
        month += 12
    a = year // 100
    b = 2 - a + a // 4
    return (int(365.25 * (year + 4716)) + int(30.6001 * (month + 1))
            + day + b - 1524.5)


def packed_date_to_jd(packed: str) -> float:
    return calendar_to_jd(*unpack_date(packed))


def unpack_number(packed: str) -> int:
    """Packed permanent number -> integer, e.g. ``'A0345'`` -> 100345.

    Returns 0 for a provisional designation (an unnumbered object).
    """
    s = packed.strip()
    if len(s) != 5:
        return 0  # 7-character provisional designation
    c = s[0]
    if c == "~":
        v = 0
        for ch in s[1:]:
            v = v * 62 + _BASE62.index(ch)
        return 620000 + v
    if c.isdigit():
        return int(s)
    if "A" <= c <= "Z":
        return (ord(c) - ord("A") + 10) * 10000 + int(s[1:])
    if "a" <= c <= "z":
        return (ord(c) - ord("a") + 36) * 10000 + int(s[1:])
    raise ValueError(f"bad packed number {packed!r}")


# --------------------------------------------------------------------------
# Fixed-width parsing
# --------------------------------------------------------------------------

def _column(block: np.ndarray, first: int, last: int) -> np.ndarray:
    """Slice columns first..last (1-based, inclusive) out of an (N, 202)
    byte matrix as an array of fixed-width byte strings."""
    sub = np.ascontiguousarray(block[:, first - 1:last])
    return sub.view(f"S{last - first + 1}").ravel()


def _to_float(raw: np.ndarray) -> np.ndarray:
    """Byte strings -> float64, with blank fields as NaN."""
    stripped = np.char.strip(raw)
    out = np.full(len(raw), np.nan)
    ok = stripped != b""
    out[ok] = stripped[ok].astype(np.float64)
    return out


def _parse_block(lines: list[bytes]) -> dict[str, np.ndarray]:
    buf = b"".join(line.ljust(LINE_WIDTH)[:LINE_WIDTH] for line in lines)
    block = np.frombuffer(buf, dtype="S1").reshape(-1, LINE_WIDTH)
    raw = {name: _column(block, a, b) for name, a, b in FIELDS}
    out: dict[str, np.ndarray] = {}
    for name in FLOAT_FIELDS:
        out[name] = _to_float(raw[name])
    for name in ("n_obs", "n_opp"):
        out[name] = np.nan_to_num(_to_float(raw[name]), nan=0).astype(np.int32)

    # U is a digit for most orbits; one-opposition orbits may carry E (the
    # eccentricity was assumed, not fitted), D or F instead.
    u = raw["U"]
    out["U"] = np.where(np.char.isdigit(u), np.char.strip(u), b"-1").astype(np.int8)
    out["e_assumed"] = (u == b"E") | (u == b"F")

    # Arc: "YYYY-YYYY" for multi-opposition orbits, "NNNN days" otherwise.
    arc = raw["arc"]
    multi = _column(block, 132, 132) == b"-"
    arc_days = np.zeros(len(arc))
    first_year = np.zeros(len(arc), dtype=np.int16)
    if multi.any():
        y1 = _column(block[multi], 128, 131).astype(np.int16)
        y2 = _column(block[multi], 133, 136).astype(np.int16)
        first_year[multi] = y1
        arc_days[multi] = (y2 - y1) * 365.25
    single = ~multi
    if single.any():
        arc_days[single] = _to_float(_column(block[single], 128, 131))
    out["arc_days"] = arc_days
    out["first_obs_year"] = first_year

    flags = np.char.strip(raw["flags"])
    out["flags"] = np.array([int(f, 16) if f else 0 for f in flags.tolist()],
                            dtype=np.int32)
    out["packed"] = raw["packed"]
    out["epoch_packed"] = raw["epoch"]
    out["readable"] = np.char.strip(raw["readable"])
    return out


def parse_mpcorb(path: Path = MPCORB_PATH, *, chunk: int = 200_000,
                 quiet: bool = False) -> dict[str, np.ndarray]:
    """Parse the whole file. Streams in chunks to keep memory modest."""
    parts: list[dict] = []
    sections: list[np.ndarray] = []
    section = 0
    in_section = False
    pending: list[bytes] = []

    def flush():
        if pending:
            parts.append(_parse_block(pending))
            sections.append(np.full(len(pending), section, dtype=np.int8))
            pending.clear()

    with gzip.open(path, "rb") as f:
        for line in f:
            if line.startswith(b"-----"):
                break
        for line in f:
            line = line.rstrip(b"\r\n")
            if not line.strip():
                # A blank line after data closes the current section.
                if in_section:
                    flush()
                    section += 1
                    in_section = False
                continue
            if len(line) < 160:
                raise ValueError(f"short MPCORB record: {line[:60]!r}")
            in_section = True
            pending.append(line)
            if len(pending) >= chunk:
                flush()
                if not quiet:
                    print(f"    parsed {sum(len(s) for s in sections):,} lines", flush=True)
        flush()

    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    out["section"] = np.concatenate(sections)

    # Unpack numbers and epochs. There are only a few hundred distinct epochs,
    # so decode each once rather than 1.5 million times.
    numbered = out["section"] == 0
    number = np.zeros(len(numbered), dtype=np.int64)
    number[numbered] = [unpack_number(p) for p in
                        np.char.decode(out["packed"][numbered], "ascii").tolist()]
    out["number"] = number
    uniq, inverse = np.unique(out["epoch_packed"], return_inverse=True)
    jd = np.array([packed_date_to_jd(u.decode()) for u in uniq])
    out["epoch_jd"] = jd[inverse]
    return out


def load_mpcorb(*, refresh: bool = False, quiet: bool = False) -> dict[str, np.ndarray]:
    """Parsed MPCORB, from the ``.npz`` cache when it is newer than the file."""
    if (not refresh and MPCORB_PARSED.exists()
            and MPCORB_PARSED.stat().st_mtime > MPCORB_PATH.stat().st_mtime):
        with np.load(MPCORB_PARSED) as z:
            return {k: z[k] for k in z.files}
    if not MPCORB_PATH.exists():
        raise FileNotFoundError(f"{MPCORB_PATH} missing; run fetch_data.py first")
    if not quiet:
        print("  parsing MPCORB.DAT (once; cached afterwards)")
    data = parse_mpcorb(quiet=quiet)
    np.savez(MPCORB_PARSED, **data)
    return data


def orbit_type(flags: np.ndarray) -> np.ndarray:
    """The MPC's own orbit-type code (low 6 bits of the hex flags)."""
    return flags & 63
