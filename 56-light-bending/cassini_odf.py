"""Reader for DSN Orbit Data Files (ODF, TRK-2-18 format) from the PDS.

An ODF is a sequence of 36-byte big-endian records in groups, each starting
with a header record whose first word is the group's primary key and whose
fourth word is the header's own record number: file label (101), identifiers
(107), the orbit data (109: one record per observable), ramp groups (2030,
one per station: the transmitter frequency history), clock offsets (2040),
the data summary (105) and the end-of-file record (-1). Field positions follow
the PDS labels that accompany each file (CORS_0021-0028, "TRK-2-18, version
of 1 August 1996"). Times are UTC seconds since 1950-01-01T00:00, counted as
86,400 per day (no leap seconds), so they convert to calendar UTC directly.

Doppler data types: 11 one-way, 12 two-way, 13 three-way (Hz). Band IDs:
1 = S, 2 = X, 3 = Ka.

The two-way Doppler observable (checked on these files, see
``cassini_doppler.receive_frequency``) is

    F = K x f_ref - f_R,

with f_ref the record's reference frequency (items 18-19; for ramped data a
constant "receiver frequency" at the S-band-equivalent level, 2.117 GHz),
f_R the received sky frequency averaged over the count interval (item 21,
1 s here) centred on the time tag, and K = 880/221 for an X-band downlink,
3344/221 for Ka band. When f_ref changes in mid-pass, K f_ref - F stays
continuous, which is how the relation was confirmed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

EPOCH_1950 = np.datetime64("1950-01-01T00:00:00", "ns")
#: primary keys of group headers (TRK-2-18): label, identifier, data summary,
#: orbit data, ramps, clock offsets, (spare) and end of file
GROUP_KEYS = {101, 105, 107, 109, 2030, 2040, 2050, -1}


def utc_from_1950(t1950) -> np.ndarray:
    """ODF seconds since 1950 (UTC, 86400 s per day) -> numpy datetime64[ns]."""
    t = np.asarray(t1950, dtype=float)
    whole = np.floor(t)
    return (EPOCH_1950 + whole.astype(np.int64).astype("timedelta64[s]")
            + np.round((t - whole) * 1e9).astype(np.int64).astype("timedelta64[ns]"))


def _groups(raw: np.ndarray) -> list[tuple[int, int, int]]:
    """(record index, primary key, secondary key) of every group header."""
    keys = raw[:, 0].view(">i4")
    idx = np.flatnonzero(np.isin(keys, list(GROUP_KEYS)))
    # a header's word 4 is the "group start packet number", i.e. its own record
    # index; a data record that happens to start with a key value fails this
    return [(int(i), int(keys[i]), int(raw[i, 1])) for i in idx if raw[i, 3] == i or keys[i] == -1]


def _field(words: np.ndarray, start: int, n: int) -> np.ndarray:
    """Big-endian bit field (1-based start bit, n bits <= 32) from a stack of
    32-bit words; a field may straddle two adjacent words."""
    w0 = (start - 1) // 32
    off = (start - 1) % 32
    if off + n <= 32:
        v = words[:, w0].astype(np.uint64) >> np.uint64(32 - off - n)
    else:
        v = ((words[:, w0].astype(np.uint64) << np.uint64(32)) | words[:, w0 + 1].astype(np.uint64))
        v = v >> np.uint64(64 - off - n)
    return (v & np.uint64((1 << n) - 1)).astype(np.int64)


def read_odf(path: str | Path) -> dict:
    """{'doppler': orbit-data records (all data types), 'ramps': ramp tables}."""
    raw = np.fromfile(path, dtype=">u4").reshape(-1, 9)
    groups = _groups(raw)
    starts = [g[0] for g in groups] + [len(raw)]
    out = {"doppler": None, "ramps": []}
    for (i, key, sec), nxt in zip(groups, starts[1:]):
        recs = raw[i + 1:nxt]
        if key == 109:
            out["doppler"] = _orbit_data(recs)
        elif key == 2030 and len(recs):
            out["ramps"].append(_ramps(recs, sec))
    out["ramps"] = pd.concat(out["ramps"], ignore_index=True) if out["ramps"] else pd.DataFrame()
    return out


def _orbit_data(d: np.ndarray) -> pd.DataFrame:
    t_int = d[:, 0].astype(np.int64)
    t_ms = (d[:, 1] >> 22) & 0x3FF
    dl_delay_ns = d[:, 1] & 0x3FFFFF
    obs = d[:, 2].view(">i4").astype(np.float64) + d[:, 3].view(">i4").astype(np.float64) * 1e-9
    w = d[:, 4:7]
    f = {name: _field(w, s, n) for name, (s, n) in {
        "format": (1, 3), "rcv": (4, 7), "xmt": (11, 7), "network": (18, 2), "dtype": (20, 6),
        "dl_band": (26, 2), "ul_band": (28, 2), "ex_band": (30, 2), "bad": (32, 1),
        "item15": (33, 7), "sc_id": (40, 10), "rx_ex_indep": (50, 1), "ref_hi": (51, 22),
        "ref_lo": (73, 24)}.items()}
    w2 = d[:, 7:9]
    g = {name: _field(w2, s, n) for name, (s, n) in {
        "item20": (1, 20), "item21": (21, 22), "ul_delay_ns": (43, 22)}.items()}
    t1950 = t_int + t_ms * 1e-3
    df = pd.DataFrame({
        "utc": utc_from_1950(t1950), "t1950": t1950, "observable": obs,
        "ref_freq_hz": (f["ref_hi"].astype(np.float64) * (1 << 24) + f["ref_lo"]) * 1e-3,
        "tc_s": g["item21"] / 100.0, "dl_delay_ns": dl_delay_ns.astype(np.int64),
        "ul_delay_ns": g["ul_delay_ns"],
        **{k: v for k, v in f.items() if k not in ("ref_hi", "ref_lo")},
    })
    return df


def _ramps(d: np.ndarray, station: int) -> pd.DataFrame:
    t0 = d[:, 0].astype(np.int64) + d[:, 1].astype(np.float64) * 1e-9
    rate = d[:, 2].view(">i4").astype(np.float64) + d[:, 3].view(">i4").astype(np.float64) * 1e-9
    ghz = (d[:, 4] >> 10) & 0x3FFFFF
    sta = d[:, 4] & 0x3FF
    f0 = ghz.astype(np.float64) * 1e9 + d[:, 5].astype(np.float64) + d[:, 6].astype(np.float64) * 1e-9
    t1 = d[:, 7].astype(np.int64) + d[:, 8].astype(np.float64) * 1e-9
    return pd.DataFrame({"group_station": station, "station": sta, "t_start1950": t0, "t_end1950": t1,
                         "f_start_hz": f0, "rate_hz_s": rate, "sky_level": ghz > 0})


if __name__ == "__main__":
    import sys

    for p in sys.argv[1:]:
        r = read_odf(p)
        d = r["doppler"]
        print(Path(p).name, len(d), "records", d["utc"].min(), "->", d["utc"].max())
        print(d.groupby(["dtype", "rcv", "xmt", "ul_band", "dl_band"]).size().to_string())
        print("ramps:", len(r["ramps"]))
        if len(r["ramps"]):
            print(r["ramps"].head(8).to_string())
