"""Snapshot of the newest Type Ia supernovae: who classified what in the last
few weeks (Transient Name Server), and their ZTF photometry (ALeRCE broker).

    python fetch_live.py              # snapshot named after today's UTC date
    python fetch_live.py --days 60    # look further back than the default 45 days

Writes data/live/<YYYY-MM-DD>/:
    tns.csv                 TNS public search result (SN Ia and 91T-like)
    objects.json            TNS name -> ZTF object id, coordinates, redshifts
    alerce/<oid>_forced.json      ZTF forced photometry (from alert packets)
    alerce/<oid>_detections.json  ZTF alert detections (for reference)

A snapshot is never overwritten; the analysis uses the newest one. Re-running
this script on a later day adds a new snapshot: the Hubble diagram's "live"
points then move on to whatever was classified most recently.

Etiquette: TNS asks automated users of its API to register a bot; this script
makes a single request to the public search page per run (the same CSV a
person gets from the website's export button). ALeRCE requests are paced.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import time
from datetime import datetime, timedelta, timezone

import numpy as np
import pandas as pd
import requests

import config
from remotezip import USER_AGENT

#: TNS object types: 3 = "SN Ia", 104 = "SN Ia-91T-like". Peculiar classes
#: (91bg-like 103, Iax 105, Ia-CSM 106, Ia-pec 100, super-Chandra 102,
#: Ca-rich 118) are left out, as in the DR2 and Pantheon+ cosmology samples.
TNS_TYPES = (3, 104)

#: ZTF observes north of declination ~ -30 deg.
MIN_DEC = -31.0
ALERCE_PAUSE = 0.5  # seconds between ALeRCE requests


def _get(url, params=None, retries=4, timeout=60, **kw) -> requests.Response:
    delay = 3.0
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout, **kw)
            if r.status_code in (429, 500, 502, 503, 504):
                raise requests.HTTPError(f"HTTP {r.status_code}")
            return r
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as err:
            if attempt == retries - 1:
                raise
            print(f"    retry after {err}")
            time.sleep(delay)
            delay *= 2
    raise AssertionError("unreachable")


def fetch_tns(start: str, end: str) -> pd.DataFrame:
    """Classified SNe Ia discovered between two dates (inclusive), all pages."""
    frames, page = [], 0
    while True:
        params = [("discovery_date_start", start), ("discovery_date_end", end),
                  ("classified_sne", "1"), ("format", "csv"), ("num_page", "500"),
                  ("page", str(page))] + [("objtype[]", str(t)) for t in TNS_TYPES]
        # TNS serves this public search to browsers; a browser-like UA is needed.
        r = _get(config.TNS_SEARCH, params=params, headers={
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64; rv:130.0) Gecko/20100101 Firefox/130.0"})
        text = r.text
        if not text.lstrip().startswith('"ID"'):
            raise RuntimeError("TNS did not return CSV (rate-limited or changed?): "
                               + text[:200].replace("\n", " "))
        df = pd.read_csv(io.StringIO(text))
        frames.append(df)
        if len(df) < 500:
            break
        page += 1
        time.sleep(5)
    return pd.concat(frames, ignore_index=True)


def sexagesimal_to_deg(ra: str, dec: str) -> tuple[float, float]:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    return float(c.ra.deg), float(c.dec.deg)


def find_ztf_object(ra: float, dec: float, internal_name: str) -> dict | None:
    """The ZTF object at this position (nearest within 2 arcsec)."""
    if isinstance(internal_name, str) and internal_name.startswith("ZTF"):
        r = _get(f"{config.ALERCE_API}/objects/{internal_name}")
        if r.status_code == 200:
            o = r.json()
            return {"oid": o["oid"], "sep_arcsec": 0.0, "ndet": o.get("ndet"),
                    "firstmjd": o.get("firstmjd"), "lastmjd": o.get("lastmjd")}
    r = _get(f"{config.ALERCE_API}/objects/",
             params={"ra": ra, "dec": dec, "radius": 2, "page_size": 10})
    if r.status_code != 200:
        return None
    items = r.json().get("items") or []
    if not items:
        return None
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    target = SkyCoord(ra * u.deg, dec * u.deg)
    best = min(items, key=lambda o: target.separation(
        SkyCoord(o["meanra"] * u.deg, o["meandec"] * u.deg)).arcsec)
    sep = target.separation(SkyCoord(best["meanra"] * u.deg, best["meandec"] * u.deg)).arcsec
    return {"oid": best["oid"], "sep_arcsec": round(float(sep), 3), "ndet": best.get("ndet"),
            "firstmjd": best.get("firstmjd"), "lastmjd": best.get("lastmjd")}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--days", type=int, default=config.LIVE_DISCOVERY_WINDOW_DAYS)
    args = ap.parse_args()

    now = datetime.now(timezone.utc)
    snap = config.LIVE_DIR / now.strftime("%Y-%m-%d")
    (snap / "alerce").mkdir(parents=True, exist_ok=True)
    start = (now - timedelta(days=args.days)).strftime("%Y-%m-%d")
    end = now.strftime("%Y-%m-%d")

    tns_path = snap / "tns.csv"
    if tns_path.exists():
        tns = pd.read_csv(tns_path)
        print(f"TNS: cached snapshot {snap.name} ({len(tns)} objects)")
    else:
        print(f"TNS: SNe Ia discovered {start} .. {end} and already classified ...")
        tns = fetch_tns(start, end)
        tns.to_csv(tns_path, index=False)
        print(f"  {len(tns)} objects")

    objects_path = snap / "objects.json"
    objects = json.loads(objects_path.read_text()) if objects_path.exists() else {}
    n_new = 0
    for _, row in tns.iterrows():
        name = row["Name"]
        if name in objects and objects[name].get("done"):
            continue
        ra, dec = sexagesimal_to_deg(row["RA"], row["DEC"])
        entry = {
            "tns_name": name, "ra": ra, "dec": dec, "type": row["Obj. Type"],
            "redshift": None if pd.isna(row["Redshift"]) else float(row["Redshift"]),
            "host_redshift": None if pd.isna(row["Host Redshift"]) else float(row["Host Redshift"]),
            "host_name": None if pd.isna(row["Host Name"]) else str(row["Host Name"]),
            "discovery_date": row["Discovery Date (UT)"],
            "classifying_groups": None if pd.isna(row["Classifying Group/s"])
            else str(row["Classifying Group/s"]),
            "classification_bibcodes": None if pd.isna(row["Classification Bibcodes"])
            else str(row["Classification Bibcodes"]),
            "discovery_sources": str(row["Discovery Data Source/s"]),
            "internal_name": str(row["Disc. Internal Name"]),
        }
        if dec < MIN_DEC:
            entry.update(done=True, oid=None, note="south of ZTF's sky")
            objects[name] = entry
            continue
        ztf = find_ztf_object(ra, dec, entry["internal_name"])
        time.sleep(ALERCE_PAUSE)
        if ztf is None:
            entry.update(done=True, oid=None, note="no ZTF object within 2 arcsec")
            objects[name] = entry
            continue
        entry.update(ztf)
        oid = ztf["oid"]
        for kind, url in (("forced", config.ALERCE_FORCED.format(oid=oid)),
                          ("detections", f"{config.ALERCE_API}/objects/{oid}/lightcurve")):
            path = snap / "alerce" / f"{oid}_{kind}.json"
            if not path.exists():
                r = _get(url)
                path.write_text(r.text if r.status_code == 200 else "null")
                time.sleep(ALERCE_PAUSE)
        entry["done"] = True
        objects[name] = entry
        n_new += 1
        if n_new % 10 == 0:
            objects_path.write_text(json.dumps(objects, indent=1))
            print(f"  {len(objects)}/{len(tns)} objects matched", flush=True)
    meta = {"snapshot_utc": now.isoformat(timespec="seconds"), "discovered_from": start,
            "discovered_to": end, "tns_types": list(TNS_TYPES), "n_tns": int(len(tns))}
    objects_path.write_text(json.dumps(objects, indent=1))
    (snap / "snapshot.json").write_text(json.dumps(meta, indent=1))
    n_ztf = sum(1 for o in objects.values() if o.get("oid"))
    print(f"snapshot {snap.name}: {len(tns)} TNS SNe Ia, {n_ztf} with ZTF photometry")
    return 0


if __name__ == "__main__":
    sys.exit(main())
