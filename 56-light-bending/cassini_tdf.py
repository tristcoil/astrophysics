"""Reader for DSN Archival Tracking Data Files (ATDF/TDF, SFOC-NAV-2-25 format).

The ODFs of the Cassini SCE1 archive were made from these files, but they
leave out part of the closed-loop Doppler: the TDF of the X-band uplinks
("xmmm") also holds the records of DSS 25's Ka-band receiver channel. This
module turns the TDF's raw Doppler counts into the same 1-s observables as
the ODF, so the two can be checked against each other and the extra data
used.

Each 288-byte Tracking Data Logical Record (format 8; record type 91 =
high-rate data) carries an integer-second UTC time tag, the station, the
downlink and uplink bands, the ground mode (1/2/3-way), a good/bad flag, a
receiver lock flag, the Doppler bias B (kHz), the Doppler reference frequency
(sky level, the same number the ODF calls the reference frequency) and ten
cumulative Doppler counts at 0.1-s steps. The 1-s Doppler observable is

    F(t + 0.5 s) = [N(t + 1) - N(t)] / 1 s - B,

which reproduces the ODF observables exactly (checked: identical, to the
last digit, for all 74,000 X/X points of 2002 day 157 that both files
carry). The counters wrap around every ~10^8 cycles in a way that does not
follow one fixed modulus, so the ~1% of points where N jumps back are
dropped rather than unwrapped.

Field positions are from the PDS label of each TDF (CORS_0021-0028).
Caution, found on these files: the TDF labels every record of DSS 25's Ka
receiver channel with an X-band uplink, but the channel was tracking the
Ka/Ka link: the records are identical to the ODFs' Ka/Ka data and their
frequency is (14/15) x the Ka uplink, not (3344/749) x the X uplink (the two
differ by 1.6e-4; checked on six days across the experiment, including
June 14 and 17, whose Ka/Ka is in no ODF). ``cassini_doppler.load_observables``
therefore takes them as Ka/Ka; the one stretch where the channel held the
X-up/Ka-down carrier (June 8, 21:46-22:53, residual -1.57e-4) is screened out
with the other incoherent stretches (``cassini_doppler.MAX_PASS_OFFSET``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

RECORD = 288
#: |F| above this (Hz) marks a counter wrap or reset: drop the point
MAX_ABS_OBSERVABLE = 5e6
_DAYS_1950_2002 = int((np.datetime64("2002-01-01") - np.datetime64("1950-01-01")).astype("timedelta64[D]").astype(int))


def _uint(raw: np.ndarray, byte0: int, nbytes: int) -> np.ndarray:
    v = np.zeros(len(raw), dtype=np.uint64)
    for k in range(nbytes):
        v = (v << np.uint64(8)) | raw[:, byte0 + k].astype(np.uint64)
    return v


def _bits(raw: np.ndarray, byte0: int, start: int, n: int) -> np.ndarray:
    """Bit field: bytes counted from byte0 (0-based), bits from 1 (MSB first)."""
    b = byte0 + (start - 1) // 8
    off = (start - 1) % 8
    nb = (off + n + 7) // 8
    v = _uint(raw, b, nb)
    return ((v >> np.uint64(nb * 8 - off - n)) & np.uint64((1 << n) - 1)).astype(np.int64)


def read_tdf(path: str | Path, stations: tuple[int, ...] | None = None) -> pd.DataFrame:
    """High-rate Doppler records of a TDF as a table (one row per record)."""
    raw = np.fromfile(path, dtype=np.uint8)
    raw = raw[: len(raw) // RECORD * RECORD].reshape(-1, RECORD)
    fmt = _uint(raw, 0, 4).astype(np.int64)
    rtype = _uint(raw, 5, 4).astype(np.int64)
    sdt = _bits(raw, 18, 19, 6)
    keep = (fmt == 8) & (rtype == 91) & (sdt == 1)
    raw = raw[keep]
    sta = _bits(raw, 18, 1, 10)
    if stations is not None:
        m = np.isin(sta, stations)
        raw, sta = raw[m], sta[m]
    year = _bits(raw, 9, 1, 12) + 1900
    doy = _bits(raw, 9, 13, 16)
    hh, mm, ss = _bits(raw, 9, 29, 8), _bits(raw, 9, 37, 8), _bits(raw, 9, 45, 8)
    if len(raw) and not np.all(year == 2002):
        raise ValueError(f"{path}: records outside 2002")
    t1950 = ((_DAYS_1950_2002 + doy - 1) * 86400 + hh * 3600 + mm * 60 + ss).astype(np.float64)
    bias = _bits(raw, 27, 2, 18)
    bias = np.where(bias >= 1 << 17, bias - (1 << 18), bias)
    counts = np.column_stack(
        [_bits(raw, 36, 1, 24) * 1e6 + _bits(raw, 36, 25, 24) * 10.0 + _bits(raw, 36, 49, 24) * 1e-6]
        + [_bits(raw, 55, 245 + 72 * k, 24) * 1e6 + _bits(raw, 55, 269 + 72 * k, 24) * 10.0
           + _bits(raw, 55, 293 + 72 * k, 24) * 1e-6 for k in range(9)])
    ref = _bits(raw, 55, 149, 32) * 1e3 + _bits(raw, 55, 181, 32) * 1e-6
    df = pd.DataFrame({
        "t1950": t1950, "station": sta, "dl_band": _bits(raw, 18, 11, 8), "ul_band": raw[:, 180].astype(np.int64),
        "channel": _bits(raw, 18, 25, 4), "mode": _bits(raw, 18, 29, 4), "bad": _bits(raw, 27, 1, 1),
        "bias_hz": bias * 1e3, "freq_level": _bits(raw, 27, 21, 1), "out_of_lock": _bits(raw, 27, 23, 1),
        "rx_type": _bits(raw, 30, 1, 6), "sample_cs": _uint(raw, 32, 4).astype(np.int64),
        "ref_freq_hz": ref, "count0": counts[:, 0], "count9": counts[:, 9],
    })
    return df


def observables(tdf: pd.DataFrame) -> pd.DataFrame:
    """1-s Doppler observables (Hz, ODF convention) from consecutive records of
    each station/channel: F = N(t+1) - N(t) - B, time-tagged t + 0.5 s. Points
    whose two records are not one second apart, or either is flagged bad, out
    of lock, or with a different reference frequency or bias, are dropped."""
    out = []
    for (sta, ch, dl), g in tdf.groupby(["station", "channel", "dl_band"]):
        g = g.drop_duplicates("t1950").sort_values("t1950")
        t = g["t1950"].to_numpy()
        n = g["count0"].to_numpy()
        ok = (np.diff(t) == 1.0)
        for c in ("bad", "out_of_lock"):
            v = g[c].to_numpy()
            ok &= (v[:-1] == 0) & (v[1:] == 0)
        for c in ("ref_freq_hz", "bias_hz", "mode", "ul_band"):
            v = g[c].to_numpy()
            ok &= v[:-1] == v[1:]
        first = g.iloc[:-1][ok]
        obs = (n[1:] - n[:-1])[ok] - first["bias_hz"].to_numpy()
        sane = np.abs(obs) < MAX_ABS_OBSERVABLE
        first, obs, tt = first[sane], obs[sane], t[:-1][ok][sane]
        out.append(pd.DataFrame({
            "t1950": tt + 0.5, "station": sta, "channel": ch, "dl_band": dl,
            "ul_band_tdf": first["ul_band"].to_numpy(), "mode": first["mode"].to_numpy(),
            "ref_freq_hz": first["ref_freq_hz"].to_numpy(), "observable": obs}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


if __name__ == "__main__":
    import sys

    for p in sys.argv[1:]:
        d = read_tdf(p)
        o = observables(d)
        print(Path(p).name, len(d), "records,", len(o), "1-s observables")
        print(o.groupby(["station", "channel", "dl_band", "ul_band_tdf", "mode"]).size().to_string())
