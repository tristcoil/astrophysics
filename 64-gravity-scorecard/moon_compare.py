"""Where is the Moon? Along-track offsets of four lunar ephemerides from DE440.

For each ephemeris X the geocentric Moon is compared with DE440's every 15 days
(2000-2022) and the difference is expressed as the time shift tau that best maps
one onto the other, X(t) ~ DE440(t + tau): tau = (X - DE440) . v / |v|^2 with v
DE440's velocity. Positive tau means X's Moon is ahead in its orbit. What is left
after the shift (a few cm) is printed as a check that the difference really is
along the track.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import ephem

CACHE = C.DERIVED_DIR / "moon_compare.json"
HORIZONS = C.RAW_DIR / "horizons" / "moon_de441_2000_2022.txt"
DE430_ASCII = C.RAW_DIR / "jpl_ascii" / "ascp1950.430"


def _tau(day, sec, pos_km, e440):
    p, v = e440.posvel(ephem.MOON, day, sec, center=ephem.EARTH)
    p, v = p.T, v.T
    d = (pos_km - p) * 1e3
    vm = v * 1e3
    tau = np.einsum("ij,ij->i", d, vm) / np.einsum("ij,ij->i", vm, vm)
    rest = np.linalg.norm(d - vm * tau[:, None], axis=1)
    return tau, rest


def compute(force: bool = False) -> dict:
    if CACHE.exists() and not force:
        return json.loads(CACHE.read_text())
    e440 = ephem.default()
    jd = np.arange(2451544.5, 2459945.5, 15.0)
    day = np.floor(jd - 2451545.0).astype(np.int64)
    sec = (jd - 2451545.0 - day) * C.DAY
    year = 2000.0 + (jd - 2451544.5) / 365.25
    out = {}
    # DE421 (NAIF SPK)
    e421 = ephem.Ephemeris(C.NAIF_DIR / "de421.bsp", C.NAIF_DIR / "moon_pa_de421_1900-2050.bpc")
    p, _ = e421.posvel(ephem.MOON, day, sec, center=ephem.EARTH)
    tau, rest = _tau(day, sec, p.T, e440)
    out["DE421 (JPL 2008)"] = {"year": year.tolist(), "tau_ms": (tau * 1e3).tolist(), "rest_max_m": float(rest.max())}
    # DE430 (JPL ASCII, same record layout as DE440)
    recs = ephem.load_de_ascii(DE430_ASCII)
    de430 = ephem.DEAscii(recs)
    p, _ = de430("moon", day, sec)
    tau, rest = _tau(day, sec, p.T, e440)
    out["DE430 (JPL 2013)"] = {"year": year.tolist(), "tau_ms": (tau * 1e3).tolist(), "rest_max_m": float(rest.max())}
    # DE441 via Horizons (TDB, geocentric, ICRF)
    rows = []
    take = False
    for line in HORIZONS.read_text().splitlines():
        if line.startswith("$$SOE"):
            take = True
            continue
        if line.startswith("$$EOE"):
            break
        if take:
            f = [x.strip() for x in line.split(",")]
            rows.append([float(f[0])] + [float(x) for x in f[2:5]])
    h = np.array(rows)
    hd = np.floor(h[:, 0] - 2451545.0).astype(np.int64)
    hs = (h[:, 0] - 2451545.0 - hd) * C.DAY
    tau, rest = _tau(hd, hs, h[:, 1:4], e440)
    out["DE441 (JPL 2020, via Horizons)"] = {"year": (2000.0 + (h[:, 0] - 2451544.5) / 365.25).tolist(),
                                            "tau_ms": (tau * 1e3).tolist(), "rest_max_m": float(rest.max())}
    # EPM2021 (IAA RAS), 2005-2021 slice
    epm = ephem.epm2021().epm
    sel = (jd > epm.recs[0, 0] + 1) & (jd < epm.recs[-1, 1] - 1)
    p, _ = epm("moon", day[sel], sec[sel])
    tau, rest = _tau(day[sel], sec[sel], p.T, e440)
    out["EPM2021 (IAA RAS)"] = {"year": year[sel].tolist(), "tau_ms": (tau * 1e3).tolist(), "rest_max_m": float(rest.max())}
    CACHE.write_text(json.dumps(out))
    return out


def curves() -> dict:
    c = compute()
    return {k: (np.array(v["year"]), np.array(v["tau_ms"])) for k, v in c.items()}


if __name__ == "__main__":
    c = compute(force=True)
    for k, v in c.items():
        tau = np.array(v["tau_ms"])
        yr = np.array(v["year"])
        pick = [np.argmin(np.abs(yr - y)) for y in (2006, 2010, 2014, 2018)]
        print(f"{k:32s} tau(2006,2010,2014,2018) = {np.round(tau[pick], 3)} ms; max remainder {v['rest_max_m']:.3f} m")
