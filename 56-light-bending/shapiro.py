"""The Shapiro delay to Cassini during its 2002 solar conjunction.

In June 2002 Cassini, cruising to Saturn 7.4 au from the Sun, passed behind
the Sun as seen from the Earth. Radio signals between the Deep Space Network
and the spacecraft then passed within 1.7 solar radii of the Sun's centre
(1.67 along JPL's geometry; Bertotti et al. 2003 quote 1.6) and were delayed
by the curvature of space. For a signal between points at
distances r1 and r2 from the Sun, a straight-line distance R apart, the
one-way delay is (Shapiro 1964; Will 2014, eq. 64)

    dt = (1 + gamma) GM/c^3 ln[(r1 + r2 + R) / (r1 + r2 - R)]
       ~ (1 + gamma) GM/c^3 ln(4 r1 r2 / b^2)    near conjunction,

with b the impact parameter (closest approach of the ray to the Sun). The
round trip doubles it. Cassini's radio system could not see the delay itself
(it measured no ranges at the time) but the *rate of change* of the delay
shifts the frequency of the two-way radio link by y = d(dt)/dt (Bertotti,
Iess & Tortora 2003, eq. 2, which lacks a factor 2, as its authors later
noted): up to +-4.2e-10 about 13 hours either side of closest approach, and
zero at closest approach itself, where y changes sign. That frequency shift
is what Cassini measured.

The geometry comes from JPL Horizons (heliocentric ICRF vectors of Cassini
and the Earth every 10 minutes, fetched by fetch_data.py). Light time is
taken into account: the signal received at the Earth at time t3 left Cassini
at t2 = t3 - R/c and was sent from the Earth at t1 = t2 - R'/c.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C

GM_C3 = C.GM_SUN / C.C_LIGHT**3          # 4.925 microseconds
KM = 1e3

#: The first Cassini Solar Conjunction Experiment (SCE1): data volumes
#: CORS_0021-0028 at the PDS Atmospheres Node cover 2002-06-06 to 2002-07-05
#: (DATASET.CAT STOP_TIME; the stop below keeps all of 5 July).
SCE1_START, SCE1_STOP = "2002-06-06", "2002-07-06"

#: Bertotti, Iess & Tortora 2003 (Nature 425, 374): gamma - 1 = (2.1 +- 2.3)e-5
#: from 18 passes of Doppler data, about 2002 June 9-17 and June 22 - July 4
#: (read from their Fig. 3; closest approach 1.6 solar radii on June 21).


def load_vectors() -> pd.DataFrame:
    df = pd.read_csv(C.HORIZONS_DIR / "cassini_earth_2002.csv")
    df["t"] = (df["jd_tdb"] - df["jd_tdb"].iloc[0]) * C.DAY   # s since start (TDB)
    return df


def _interp(df: pd.DataFrame, body: str, t: np.ndarray) -> np.ndarray:
    return np.stack([np.interp(t, df["t"], df[f"{body}_{k}"]) for k in "xyz"], axis=-1) * KM


def one_way_delay(r1v: np.ndarray, r2v: np.ndarray, gamma: float = 1.0) -> np.ndarray:
    """Shapiro delay (s) between two heliocentric points (m)."""
    r1, r2 = np.linalg.norm(r1v, axis=-1), np.linalg.norm(r2v, axis=-1)
    R = np.linalg.norm(r2v - r1v, axis=-1)
    return (1 + gamma) * GM_C3 * np.log((r1 + r2 + R) / (r1 + r2 - R))


def impact_parameter(r1v: np.ndarray, r2v: np.ndarray) -> np.ndarray:
    """Distance (m) of the Sun's centre from the segment r1 -> r2."""
    d = r2v - r1v
    s = np.clip(-np.sum(r1v * d, axis=-1) / np.sum(d * d, axis=-1), 0.0, 1.0)
    return np.linalg.norm(r1v + s[..., None] * d, axis=-1)


def conjunction(step_s: float = 600.0, gamma: float = 1.0) -> pd.DataFrame:
    """Round-trip Shapiro delay and the two-way fractional frequency shift
    it causes, for signals received at the Earth (geocentre) at time t3."""
    df = load_vectors()
    t3 = np.arange(df["t"].iloc[0] + 3 * 3600, df["t"].iloc[-1] - 3 * 3600, step_s)
    e3 = _interp(df, "earth", t3)
    # light-time iteration for the downlink (Cassini at t2) and uplink (Earth at t1)
    t2 = t3.copy()
    for _ in range(4):
        c2 = _interp(df, "cassini", t2)
        t2 = t3 - np.linalg.norm(e3 - c2, axis=-1) / C.C_LIGHT
    t1 = t2.copy()
    for _ in range(4):
        e1 = _interp(df, "earth", t1)
        t1 = t2 - np.linalg.norm(c2 - e1, axis=-1) / C.C_LIGHT
    down = one_way_delay(c2, e3, gamma)
    up = one_way_delay(e1, c2, gamma)
    b_down = impact_parameter(c2, e3)
    b_up = impact_parameter(e1, c2)
    total = up + down
    y = np.gradient(total, t3)          # d(delay)/dt: fractional frequency shift of the 2-way link
    jd = df["jd_tdb"].iloc[0] + t3 / C.DAY
    # Bertotti et al.'s approximation y = 4 (1 + gamma) GM/c^3 (1/b) db/dt (their eq. 2 with the factor 2
    # restored; sign: delay grows as b shrinks)
    b_mean = 0.5 * (b_up + b_down)
    y_approx = -4 * (1 + gamma) * GM_C3 * np.gradient(b_mean, t3) / b_mean
    return pd.DataFrame({
        "jd_tdb": jd, "date": pd.to_datetime(jd - 2440587.5, unit="D"),
        "b_up_rsun": b_up / C.R_SUN, "b_down_rsun": b_down / C.R_SUN, "b_rsun": b_mean / C.R_SUN,
        "delay_up_us": up * 1e6, "delay_down_us": down * 1e6, "delay_2way_us": total * 1e6,
        "y_2way": y, "y_bertotti_approx": y_approx,
        "rtlt_h": (t3 - t1) / 3600.0,
        "r_cassini_au": np.linalg.norm(c2, axis=-1) / C.AU, "r_earth_au": np.linalg.norm(e3, axis=-1) / C.AU,
    })


def summary(df: pd.DataFrame | None = None) -> dict:
    df = conjunction() if df is None else df
    i = int(np.argmin(df["b_rsun"]))
    sce = df[(df["date"] >= SCE1_START) & (df["date"] <= SCE1_STOP)]
    return {
        "closest_date": str(df["date"].iloc[i])[:16], "b_min_rsun": float(df["b_rsun"].iloc[i]),
        "delay_2way_max_us": float(df["delay_2way_us"].max()),
        "delay_2way_at_sce_edges_us": [float(sce["delay_2way_us"].iloc[0]), float(sce["delay_2way_us"].iloc[-1])],
        "y_max": float(np.abs(df["y_2way"]).max()),
        "y_vs_approx_max_rel": float(np.max(np.abs(df["y_2way"] - df["y_bertotti_approx"])) / np.abs(df["y_2way"]).max()),
        "rtlt_h": float(df["rtlt_h"].iloc[i]), "r_cassini_au": float(df["r_cassini_au"].iloc[i]),
        "b_range_sce_rsun": [float(sce["b_rsun"].min()), float(sce["b_rsun"].max())],
    }


if __name__ == "__main__":
    s = summary()
    for k, v in s.items():
        print(f"{k:28s} {v}")
