"""The GPS clock budget: why satellite clocks are tuned slow before launch.

A clock in a circular orbit of radius r, compared with a clock on the geoid, runs

    fast by (W0 - GM/r)/c^2      because it sits higher in the Earth's potential, and
    slow by GM/(2 r c^2)         because it moves (v^2 = GM/r).

For GPS (r = 26,560 km) that is +45.7 and -7.2 microseconds a day: +38.6 net.
The two cancel at r = 3GM/(2 W0) = 1.5 Earth radii, about 3,200 km up; the ISS,
far below that, loses time. GPS clocks are therefore built to tick at
10.22999999543 MHz instead of 10.23 MHz (a fractional offset of -4.4647e-10).

The real-data check: for every GPS satellite in ESA's products on the first day
of each month (2016-2026), compute what general relativity says its clock gains
from its actual orbit (Earth monopole + J2, speed), and compare with the drift
the published clock actually shows against the ground reference. If the factory
offset is right, the published drift is near zero, and the drift plus the
offset reproduces the predicted gain.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
import gnss

C2 = C.C_LIGHT**2
US_PER_DAY = C.SECONDS_PER_DAY * 1e6

#: Reference orbits (circular-orbit radius in km) for the budget figure.
REFERENCE_ORBITS = {
    "ISS": 6_371 + 420,
    "GLONASS": 25_510,
    "GPS": 26_560,
    "BeiDou MEO": 27_906,
    "Galileo": 29_600,
    "geostationary": 42_164,
}


def circular_rates(r_m: np.ndarray) -> dict:
    """Rate of an orbiting clock relative to the geoid (microseconds per day)."""
    gravity = (C.W0_GEOID - C.GM_EARTH / r_m) / C2 * US_PER_DAY
    speed = -0.5 * C.GM_EARTH / r_m / C2 * US_PER_DAY
    return {"gravity": gravity, "speed": speed, "total": gravity + speed}


def crossover_radius() -> float:
    return 1.5 * C.GM_EARTH / C.W0_GEOID


def satellite_days() -> pd.DataFrame:
    rows = []
    for path in sorted(C.SP3_FULL_DIR.glob("*_ORB.SP3.gz")):
        day = gnss.read_sp3(path)
        if len(day.gps_s) != 289:
            continue
        rot = gnss.terrestrial_to_celestial(day.gps_s)
        sun_moon = gnss.sun_moon_gcrs(day.gps_s)
        t = day.gps_s - day.gps_s[0]
        for prn in sorted(p for p in day.pos if p.startswith("G")):
            od = gnss.orbit_day(day, prn, rot, sun_moon)
            if od is None:
                continue
            clk = od.clk
            ok = np.isfinite(clk)
            if ok.sum() < 250:
                continue
            el = gnss.kepler_elements(od.r, od.v)
            T = od.terms
            gain = (C.W0_GEOID / C2 - np.mean(T["rate_mono"] + T["rate_j2"] + T["rate_kin"])) * US_PER_DAY
            drift = np.polyfit(t[ok], clk[ok], 1)[0] * US_PER_DAY
            a, e = float(np.mean(el["a"])), float(np.mean(el["e"]))
            rows.append({"date": day.date.isoformat(), "prn": prn, "a_km": a / 1e3, "e": e,
                         "predicted_gain_us_day": gain, "published_drift_us_day": drift,
                         "periodic_amplitude_ns": 2 * np.sqrt(C.GM_EARTH * a) * e / C2 * 1e9})
    df = pd.DataFrame(rows)
    df["factory_offset_us_day"] = C.GPS_FACTORY_OFFSET * US_PER_DAY
    df["implied_gain_us_day"] = df["published_drift_us_day"] - df["factory_offset_us_day"]
    return df


def summary(df: pd.DataFrame) -> dict:
    r_gps = REFERENCE_ORBITS["GPS"] * 1e3
    circ = circular_rates(np.array([r_gps]))
    robust = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
    return {
        "gps_circular": {k: float(v[0]) for k, v in circ.items()},
        "crossover_radius_km": crossover_radius() / 1e3,
        "crossover_altitude_km": crossover_radius() / 1e3 - 6_371.0,
        "iss_total_us_day": float(circular_rates(np.array([REFERENCE_ORBITS["ISS"] * 1e3]))["total"][0]),
        "factory_offset_us_day": C.GPS_FACTORY_OFFSET * US_PER_DAY,
        "n_satellite_days": int(len(df)),
        "n_satellites": int(df["prn"].nunique()),
        "predicted_gain_median": float(df["predicted_gain_us_day"].median()),
        "predicted_gain_range": [float(df["predicted_gain_us_day"].min()), float(df["predicted_gain_us_day"].max())],
        "published_drift_median": float(df["published_drift_us_day"].median()),
        "published_drift_robust_sd": float(robust(df["published_drift_us_day"])),
        "published_drift_abs_p95": float(np.percentile(np.abs(df["published_drift_us_day"]), 95)),
        "implied_gain_median": float(df["implied_gain_us_day"].median()),
        "range_error_if_ignored_km_per_day": C.C_LIGHT * 38.6e-6 / 1e3,
        "periodic_amplitude_ns_max": float(df["periodic_amplitude_ns"].max()),
        "periodic_amplitude_ns_median": float(df["periodic_amplitude_ns"].median()),
    }
