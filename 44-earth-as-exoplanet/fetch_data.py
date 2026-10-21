"""Download and cache everything the analysis needs.

Public sources, no account needed (the calibrated level-1B files at NASA's
ASDC need an Earthdata Login, so this project does without them):

1. **EPIC API metadata** (epic.gsfc.nasa.gov): for every natural-colour image,
   its time, the latitude and longitude at the centre of the disc, and the
   positions of DSCOVR, the Moon and the Sun relative to Earth (J2000, km).
   One small JSON file per day, ~3,600 days since June 2015.
2. **Natural-colour images** from the same archive: 1080-pixel JPEGs for two
   whole years and three special days, and 120-pixel thumbnails of all
   ~55,000 images.
3. **Gallery sets** of special events (lunar transits and an occultation),
   full-size PNGs from the EPIC team's galleries.
4. **NOAA's land-sea mask** (0.25 degrees), to weigh the disc by surface type.
5. **JPL Horizons** vectors of DSCOVR, the Moon and the Sun (``ephemeris.py``),
   and **NASA's CMR catalogue**: how many granules of the login-gated EPIC
   collections exist during each predicted Moon event (metadata only).
6. **ASDC documents**: the L1B format book and the table of calibration
   factors that turns EPIC's counts per second into reflectance.

Every file lands under ``data/raw/``; a manifest (``data/manifest.json``)
records where each came from, when, its size and SHA-256 (image sets as one
entry each, with a hash over the per-file hashes). Nothing already cached is
fetched again, so after one run the project works offline.

    python fetch_data.py                  # fill the cache (~2.7 GB, a few hours)
    python fetch_data.py --metadata-only  # just the API metadata and documents
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import requests

import epicdata
from config import (
    API_DIR, API_NATURAL_DIR, ASDC_DOCS, CMR_GRANULES, DOCS, DOCS_DIR, EPIC_API, L1B_COLLECTION_V3,
    EPIC_POLITE_DELAY_S, EPIC_SITE, EPIC_WORKERS, GALLERIES, GALLERY_DIR,
    JPG_WINDOWS, LANDMASK_PATH, LANDMASK_URL, MANIFEST_PATH, RAW_DIR, SPECIAL_DAYS,
)

USER_AGENT = "astrophysics-atlas-A4/1.0 (hobbyist research; Earth as an exoplanet)"
TIMEOUT_S = 120
RETRIES = 6

_manifest_lock = threading.Lock()
_local = threading.local()


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}, "collections": {}}


def _save_manifest(manifest: dict) -> None:
    tmp = MANIFEST_PATH.with_name(MANIFEST_PATH.name + ".part")
    tmp.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    tmp.replace(MANIFEST_PATH)


def _entry(path: Path, url: str, source: str, **extra) -> dict:
    data = path.read_bytes()
    return {"url": url, "source": source, "retrieved_utc": _now_utc(), "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(), **extra}


def _record(path: Path, url: str, source: str, **extra) -> None:
    """Add one downloaded file to the provenance manifest."""
    with _manifest_lock:
        manifest = _load_manifest()
        manifest.setdefault("files", {})[str(path.relative_to(RAW_DIR.parent))] = \
            _entry(path, url, source, **extra)
        _save_manifest(manifest)


def _session() -> requests.Session:
    """One HTTP session per thread (requests sessions are not thread-safe)."""
    s = getattr(_local, "session", None)
    if s is None:
        s = requests.Session()
        s.headers["User-Agent"] = USER_AGENT
        _local.session = s
    return s


def http_get(url: str, *, params: dict | None = None, stream_to: Path | None = None) -> bytes:
    """GET with retries. DNS lookups on this machine fail intermittently, so
    connection errors are retried quickly; HTTP errors other than 429/5xx are not."""
    last = None
    for attempt in range(RETRIES):
        try:
            r = _session().get(url, params=params, timeout=TIMEOUT_S, stream=stream_to is not None)
            if r.status_code == 429 or r.status_code >= 500:
                last = f"HTTP {r.status_code}"
                time.sleep(3.0 * (attempt + 1))
                continue
            r.raise_for_status()
            if stream_to is None:
                return r.content
            tmp = stream_to.with_name(stream_to.name + ".part")
            with open(tmp, "wb") as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
            tmp.replace(stream_to)
            return b""
        except (requests.ConnectionError, requests.Timeout) as err:
            last = err
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"giving up on {url}: {last}")


def _write_atomic(path: Path, data: bytes) -> None:
    """Write via a temporary name so an interrupted download never leaves a
    truncated file where the analysis expects a complete one."""
    tmp = path.with_name(path.name + ".part")
    tmp.write_bytes(data)
    tmp.replace(path)


# ---------------------------------------------------------------------------
# 1. EPIC API metadata
# ---------------------------------------------------------------------------

DATES_PATH = API_DIR / "natural_all.json"


def fetch_dates(*, refresh: bool = False) -> list[str]:
    """Every date with natural-colour images, oldest first."""
    if refresh or not DATES_PATH.exists():
        url = f"{EPIC_API}/natural/all"
        data = http_get(url)
        json.loads(data)  # must parse
        _write_atomic(DATES_PATH, data)
        _record(DATES_PATH, url, "EPIC API: list of dates with natural-colour images")
    return sorted(d["date"] for d in json.loads(DATES_PATH.read_text()))


def day_path(date: str) -> Path:
    return API_NATURAL_DIR / f"{date}.json"


def _fetch_day(date: str, refresh: bool) -> tuple[str, bool]:
    path = day_path(date)
    if path.exists() and not refresh:
        return date, False
    data = http_get(f"{EPIC_API}/natural/date/{date}")
    images = json.loads(data)
    if not isinstance(images, list):
        raise RuntimeError(f"unexpected API answer for {date}: {data[:200]!r}")
    _write_atomic(path, data)
    time.sleep(EPIC_POLITE_DELAY_S)
    return date, True


def fetch_metadata(*, refresh: bool = False, quiet: bool = False) -> list[str]:
    """Per-day image metadata for every date the API lists. Returns the dates."""
    dates = fetch_dates(refresh=refresh)
    todo = [d for d in dates if refresh or not day_path(d).exists()]
    if todo and not quiet:
        print(f"EPIC metadata: {len(todo)} of {len(dates)} days to fetch")
    done = 0
    with ThreadPoolExecutor(EPIC_WORKERS) as pool:
        futures = [pool.submit(_fetch_day, d, refresh) for d in todo]
        for fut in as_completed(futures):
            fut.result()
            done += 1
            if not quiet and done % 200 == 0:
                print(f"  {done}/{len(todo)} days", flush=True)
    if todo:
        # One manifest entry for the whole collection: per-file hashes in order.
        h = hashlib.sha256()
        total = 0
        for d in dates:
            b = day_path(d).read_bytes()
            total += len(b)
            h.update(hashlib.sha256(b).digest())
        with _manifest_lock:
            manifest = _load_manifest()
            manifest.setdefault("collections", {})["api/natural"] = {
                "url_pattern": f"{EPIC_API}/natural/date/YYYY-MM-DD",
                "source": "EPIC API: metadata of every natural-colour image",
                "retrieved_utc": _now_utc(),
                "n_files": len(dates), "bytes": total,
                "first": dates[0], "last": dates[-1],
                "sha256_of_file_hashes": h.hexdigest(),
            }
            _save_manifest(manifest)
    return dates


# ---------------------------------------------------------------------------
# 2. Natural-colour images
# ---------------------------------------------------------------------------

def _fetch_image(kind: str, image: str) -> int:
    path = epicdata.image_path(kind, image)
    if path.exists():
        return 0
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        data = http_get(epicdata.image_url(kind, image))
    except requests.HTTPError as err:
        if err.response is not None and err.response.status_code == 404:
            return -1  # listed by the API but not on the server: recorded as missing
        raise
    _write_atomic(path, data)
    time.sleep(EPIC_POLITE_DELAY_S)
    return len(data)


def fetch_images(kind: str, images: list[str], label: str, *, quiet: bool = False) -> dict:
    """Download one size of a list of images; record the set in the manifest."""
    todo = [im for im in images if not epicdata.image_path(kind, im).exists()]
    missing = []
    if todo and not quiet:
        print(f"{label}: {len(todo)} of {len(images)} {kind} images to fetch", flush=True)
    done = 0
    with ThreadPoolExecutor(EPIC_WORKERS) as pool:
        futures = {pool.submit(_fetch_image, kind, im): im for im in todo}
        for fut in as_completed(futures):
            if fut.result() < 0:
                missing.append(futures[fut])
            done += 1
            if not quiet and done % 500 == 0:
                print(f"  {done}/{len(todo)}", flush=True)
    present = [im for im in images if epicdata.image_path(kind, im).exists()]
    h = hashlib.sha256()
    total = 0
    for im in present:
        b = epicdata.image_path(kind, im).read_bytes()
        total += len(b)
        h.update(hashlib.sha256(b).digest())
    with _manifest_lock:
        manifest = _load_manifest()
        old = manifest.get("collections", {}).get(f"natural/{kind}/{label}", {})
        manifest.setdefault("collections", {})[f"natural/{kind}/{label}"] = {
            "url_pattern": epicdata.image_url(kind, "epic_1b_YYYYMMDDhhmmss"),
            "source": f"EPIC natural-colour images ({kind}), {label}",
            "retrieved_utc": _now_utc() if todo else old.get("retrieved_utc", _now_utc()),
            "n_listed": len(images), "n_files": len(present),
            "missing_on_server": sorted(set(old.get("missing_on_server", [])) | set(missing)),
            "bytes": total, "sha256_of_file_hashes": h.hexdigest(),
        }
        _save_manifest(manifest)
    return {"n": len(present), "bytes": total, "missing": missing}


def fetch_all_images(meta, *, quiet: bool = False) -> None:
    for label, (start, end) in JPG_WINDOWS.items():
        fetch_images("jpg", epicdata.between(meta, start, end)["image"].tolist(), label, quiet=quiet)
    for label, day in SPECIAL_DAYS.items():
        fetch_images("jpg", epicdata.between(meta, day, day)["image"].tolist(), label, quiet=quiet)
    fetch_images("thumbs", meta["image"].tolist(), "all", quiet=quiet)


# ---------------------------------------------------------------------------
# 3. Gallery image sets (special events)
# ---------------------------------------------------------------------------

def gallery_files(slug: str, page: str) -> list[str]:
    """Full-size image URLs listed on one EPIC gallery page (cached page)."""
    html_path = GALLERY_DIR / slug / "page.html"
    html_path.parent.mkdir(parents=True, exist_ok=True)
    if not html_path.exists():
        url = f"{EPIC_SITE}/{page}"
        _write_atomic(html_path, http_get(url))
        _record(html_path, url, f"EPIC gallery page: {slug}")
    text = html_path.read_text(errors="replace")
    links = re.findall(r"""(?:href|src)=['"]([^'"]+)['"]""", text)
    full = sorted({u for u in links if "epic-galleries" in u and "/full/" in u})
    return [u if u.startswith("http") else f"{EPIC_SITE}/{u.lstrip('/')}" for u in full]


def fetch_galleries(*, quiet: bool = False) -> None:
    for slug, page in GALLERIES.items():
        for url in gallery_files(slug, page):
            path = GALLERY_DIR / slug / url.rsplit("/", 1)[1]
            if path.exists():
                continue
            if not quiet:
                print(f"gallery {slug}: {path.name}", flush=True)
            http_get(url, stream_to=path)
            _record(path, url, f"EPIC gallery: {slug}")
            time.sleep(EPIC_POLITE_DELAY_S)


# ---------------------------------------------------------------------------
# 4. A land-sea mask (for "what fraction of the dot's light comes from land")
# ---------------------------------------------------------------------------

def fetch_landmask(*, refresh: bool = False) -> Path:
    path = LANDMASK_PATH
    if refresh or not path.exists():
        http_get(LANDMASK_URL, stream_to=path)
        _record(path, LANDMASK_URL, "NOAA OISST v2 high-resolution land-sea mask (0.25 deg)")
    return path


# ---------------------------------------------------------------------------
# 5. Which Moon events EPIC imaged (NASA's public CMR catalogue)
# ---------------------------------------------------------------------------

CMR_COVERAGE_PATH = RAW_DIR / "cmr" / "moon_event_coverage.json"
CMR_COLLECTIONS = {
    "L1B_v3": L1B_COLLECTION_V3,
    "L1A_lunar": "C4025357001-LARC_CLOUD",  # DSCOVR_EPIC_L1A_LUNAR V1
    "L1A_misc": "C3866447285-LARC_CLOUD",   # DSCOVR_EPIC_L1A_MISC V1
}
N_CONTROL_WINDOWS = 40


def _cmr_count(collection: str, t0: str, t1: str) -> int:
    data = json.loads(http_get(CMR_GRANULES, params={
        "collection_concept_id": collection, "temporal": f"{t0},{t1}", "page_size": 200}))
    return len(data["feed"]["entry"])


def fetch_cmr_event_coverage(events, *, refresh: bool = False) -> dict:
    """Granule counts of the login-gated EPIC collections during each predicted
    transit/occultation, and in random 4-hour control windows away from events.

    Only catalogue metadata (public) is read; no data file is downloaded."""
    import numpy as np
    import pandas as pd

    if CMR_COVERAGE_PATH.exists() and not refresh:
        return json.loads(CMR_COVERAGE_PATH.read_text())
    CMR_COVERAGE_PATH.parent.mkdir(parents=True, exist_ok=True)
    fmt = "%Y-%m-%dT%H:%M:%SZ"
    rows = []
    for _, e in events[events["kind"] != "near miss"].iterrows():
        t0, t1 = e["overlap_start"].strftime(fmt), e["overlap_end"].strftime(fmt)
        rows.append({"kind": e["kind"], "start": t0, "end": t1,
                     **{k: _cmr_count(c, t0, t1) for k, c in CMR_COLLECTIONS.items()}})
    rng = np.random.default_rng(20150716)
    lo, hi = pd.Timestamp("2015-07-01"), pd.Timestamp("2026-09-01")
    controls = []
    while len(controls) < N_CONTROL_WINDOWS:
        t = lo + (hi - lo) * float(rng.random())
        near = ((events["fov_start"] - pd.Timedelta("1D") < t) & (events["fov_end"] + pd.Timedelta("1D") > t)).any()
        if near:
            continue
        t0, t1 = t.strftime(fmt), (t + pd.Timedelta("4h")).strftime(fmt)
        controls.append({"start": t0, "end": t1, "L1A_misc": _cmr_count(CMR_COLLECTIONS["L1A_misc"], t0, t1)})
    out = {"retrieved_utc": _now_utc(), "collections": CMR_COLLECTIONS, "events": rows, "controls": controls}
    _write_atomic(CMR_COVERAGE_PATH, json.dumps(out, indent=1).encode())
    _record(CMR_COVERAGE_PATH, CMR_GRANULES, "NASA CMR granule counts for EPIC collections during Moon events")
    return out


# ---------------------------------------------------------------------------
# 6. ASDC documents
# ---------------------------------------------------------------------------

def fetch_docs(*, refresh: bool = False) -> None:
    for name in DOCS.values():
        path = DOCS_DIR / name
        if path.exists() and not refresh:
            continue
        url = ASDC_DOCS + name
        http_get(url, stream_to=path)
        _record(path, url, "ASDC DSCOVR EPIC documentation")


# ---------------------------------------------------------------------------

def ensure_data(*, refresh: bool = False, quiet: bool = False, images: bool = True) -> None:
    fetch_docs(refresh=refresh)
    fetch_landmask()
    fetch_metadata(refresh=refresh, quiet=quiet)
    if images:
        meta = epicdata.load_metadata()
        fetch_galleries(quiet=quiet)
        fetch_all_images(meta, quiet=quiet)
    # Moon events from JPL Horizons (cached), then the catalogue check.
    import ephemeris
    from config import CACHE_DIR
    import pandas as pd
    ev_path = CACHE_DIR / "moon_events.parquet"
    if not ev_path.exists():
        ephemeris.mission_events().to_parquet(ev_path)
    fetch_cmr_event_coverage(pd.read_parquet(ev_path))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download the metadata")
    parser.add_argument("--metadata-only", action="store_true", help="skip the images")
    args = parser.parse_args(argv)
    ensure_data(refresh=args.refresh, images=not args.metadata_only)
    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
