"""Turn the daily SP3 files into one array file per satellite.

For every satellite and day: the published clock, the IGS convention term
-2 r.v/c^2, the four parts of the general-relativistic proper time (Earth's
monopole, its flattening J2, lunisolar tides, speed), and the orbit geometry
(eccentric anomaly, argument of latitude, the Sun's elevation above the orbit
plane beta, the angle from the Sun mu, Earth-shadow flag, radius).

    python prepare.py          # ~5-10 min for ten years, cached in data/cache/
"""

from __future__ import annotations

import multiprocessing as mp
import time
from pathlib import Path

import numpy as np

import config as C
import gnss

N_EPOCHS = 289          # 00:00 to 24:00 inclusive, every 5 min
FIELDS_F8 = ("clk", "kc", "mono", "kin")
FIELDS_F4 = ("j2", "tide", "E", "M", "u", "mu", "beta", "r", "a", "e")


def cache_path(prn: str) -> Path:
    return C.CACHE_DIR / f"satdays_{prn}.npz"


def process_day(path: Path) -> dict:
    day = gnss.read_sp3(path)
    if len(day.gps_s) != N_EPOCHS:
        return {"path": str(path), "skip": f"{len(day.gps_s)} epochs"}
    rot = gnss.terrestrial_to_celestial(day.gps_s)
    sun, moon = gnss.sun_moon_gcrs(day.gps_s)
    out = {"path": str(path), "date": day.date.isoformat(), "series": day.series,
           "gps0": float(day.gps_s[0]), "sats": {}}
    for prn in C.KEEP_SATS:
        od = gnss.orbit_day(day, prn, rot, (sun, moon))
        if od is None:
            continue
        el = gnss.kepler_elements(od.r, od.v)
        geo = gnss.sun_geometry(od.r, od.v, sun)
        rec = {
            "clk": od.clk, "kc": od.terms["kepler_convention"],
            "mono": od.terms["mono"], "kin": od.terms["kin"],
            "j2": od.terms["j2"], "tide": od.terms["tide"],
            "E": el["E"], "M": el["M"], "u": el["u"], "mu": geo["mu"], "beta": geo["beta"],
            "r": el["r"], "a": el["a"], "e": el["e"], "shadow": geo["shadow"],
            "flags": "".join(sorted({f.strip() for f in day.flags.get(prn, []) if f.strip()})),
        }
        out["sats"][prn] = rec
    return out


def build(workers: int = 4) -> None:
    files = gnss.sp3_files()
    nf = len(files)
    print(f"{nf} daily SP3 files")
    t0 = time.time()
    # preallocate one block per satellite; rows of days without that satellite stay masked
    store = {}
    for prn in C.KEEP_SATS:
        s = {f: np.full((nf, N_EPOCHS), np.nan) for f in FIELDS_F8}
        s.update({f: np.full((nf, N_EPOCHS), np.nan, dtype=np.float32) for f in FIELDS_F4})
        s["shadow"] = np.zeros((nf, N_EPOCHS), bool)
        s["flags"] = np.full(nf, "", dtype=object)
        s["have"] = np.zeros(nf, bool)
        store[prn] = s
    dates = np.full(nf, "", dtype=object)
    series = np.full(nf, "", dtype=object)
    gps0 = np.full(nf, np.nan)
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        for k, res in enumerate(pool.imap(process_day, files, chunksize=8)):
            if "skip" in res:
                print("  skipped", res["path"], res["skip"])
                continue
            dates[k], series[k], gps0[k] = res["date"], res["series"], res["gps0"]
            for prn, rec in res["sats"].items():
                s = store[prn]
                for f in FIELDS_F8 + FIELDS_F4 + ("shadow",):
                    s[f][k] = rec[f]
                s["flags"][k] = rec["flags"]
                s["have"][k] = True
            if (k + 1) % 500 == 0:
                print(f"  {k + 1}/{nf} days, {time.time() - t0:.0f} s", flush=True)
    for prn, s in store.items():
        keep = s.pop("have")
        if not keep.any():
            continue
        arrays = {f: v[keep] for f, v in s.items()}
        arrays["flags"] = arrays["flags"].astype(str)
        arrays.update(date=dates[keep].astype(str), series=series[keep].astype(str), gps0=gps0[keep])
        np.savez_compressed(cache_path(prn), **arrays)
        print(f"  {prn}: {keep.sum()} days -> {cache_path(prn).name}")
    print(f"done in {time.time() - t0:.0f} s")


def load(prn: str) -> dict:
    with np.load(cache_path(prn)) as z:
        return {k: z[k] for k in z.files}


if __name__ == "__main__":
    build()
