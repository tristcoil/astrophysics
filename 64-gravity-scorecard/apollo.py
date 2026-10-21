"""Read APOLLO's lunar laser ranging normal points.

Each 79-character line is one normal point: the UTC time the laser fired, the
measured round-trip time to a reflector on the Moon and back (to 0.1 ps), its
uncertainty, and the weather at the telescope. The layout is documented on
``tmurphy.physics.ucsd.edu/apollo/norm_pts.html``; our column slices reproduce
the worked example there field for field (see ``test_physics.py``).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

#: (name, start, stop) of each fixed-width field; scale converts to SI.
FIELDS = [
    ("year", 2, 6), ("month", 6, 8), ("day", 8, 10), ("hour", 10, 12), ("minute", 12, 14),
    ("sec_1e7", 14, 23), ("rt_1e13", 23, 37), ("reflector", 37, 38), ("station", 38, 43),
    ("photons", 43, 46), ("sigma_01ps", 46, 52), ("snr_x10", 52, 55), ("grade", 55, 56),
    ("p_001mbar", 56, 62), ("t_01c", 62, 66), ("rh_pct", 66, 68), ("wavelength_A", 68, 73),
    ("duration_s", 74, 78),
]


def parse_line(line: str) -> dict:
    line = line.rstrip("\n")
    if len(line) != 79:
        raise ValueError(f"unexpected normal-point length {len(line)}: {line!r}")
    raw = {name: line[a:b] for name, a, b in FIELDS}
    out = {k: int(v) for k, v in raw.items() if k != "grade"}
    out["grade"] = raw["grade"].strip()
    return out


def load(groups: tuple[str, ...] | None = None) -> pd.DataFrame:
    """All normal points as a DataFrame, one row each, in SI units.

    Columns: utc (string ISO, fire time), group (a-f), reflector (name),
    rt_s (round trip, s; exact integer in 1e-13 s kept as rt_int),
    sigma_s (APOLLO's statistical uncertainty, s), pressure_hpa, temp_c, rh_pct,
    wavelength_um, photons, snr, grade, duration_s.
    """
    rows = []
    for group, name in C.APOLLO_FILES.items():
        if groups and group not in groups:
            continue
        for line in (C.APOLLO_DIR / name).read_text().splitlines():
            if not line.strip():
                continue
            r = parse_line(line)
            r["group"] = group
            rows.append(r)
    df = pd.DataFrame(rows)
    sec = df["sec_1e7"] // 10_000_000
    frac = (df["sec_1e7"] % 10_000_000) / 1e7
    df["utc_sec_int"] = sec
    df["utc_sec_frac"] = frac
    df["reflector_code"] = df["reflector"]
    df["reflector"] = df["reflector_code"].map(C.APOLLO_REFLECTOR_CODE)
    df["rt_int"] = df["rt_1e13"]
    df["rt_s"] = df["rt_1e13"] * 1e-13
    df["sigma_s"] = df["sigma_01ps"] * 1e-13
    df["pressure_hpa"] = df["p_001mbar"] / 100.0
    df["temp_c"] = df["t_01c"] / 10.0
    df["wavelength_um"] = df["wavelength_A"] / 1e4
    df["snr"] = df["snr_x10"] / 10.0
    df = df.drop(columns=["sec_1e7", "rt_1e13", "sigma_01ps", "p_001mbar", "t_01c", "wavelength_A", "snr_x10"])
    df = df.sort_values(["year", "month", "day", "hour", "minute", "utc_sec_int", "utc_sec_frac"]).reset_index(drop=True)
    if not (df["station"] == C.APOLLO_STATION_ID).all():
        raise ValueError("normal points from a station other than APOLLO")
    return df


def scaled_sigma(df: pd.DataFrame, how: str = "rss") -> np.ndarray:
    """APOLLO's recommended error inflation per group ("factor" or "rss"), in seconds."""
    out = np.empty(len(df))
    for g, (factor, rss) in C.APOLLO_ERROR_SCALE.items():
        m = (df["group"] == g).to_numpy()
        s = df.loc[m, "sigma_s"].to_numpy()
        out[m] = s * factor if how == "factor" else np.hypot(s, rss)
    return out


if __name__ == "__main__":
    d = load()
    print(len(d), "normal points;", d.groupby("group").size().to_dict())
    print(d.groupby("reflector").size().to_dict())
    print("round trip %.4f-%.4f s; median sigma %.1f ps (%.1f mm one-way)" % (
        d.rt_s.min(), d.rt_s.max(), 1e12 * d.sigma_s.median(), 1e3 * d.sigma_s.median() * C.C_LIGHT / 2))
