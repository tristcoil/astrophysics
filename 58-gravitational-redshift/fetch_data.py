"""Download and cache every input this project uses.

    python fetch_data.py              # everything
    python fetch_data.py sp3 slr      # only some sources

1. **ESA/ESOC daily SP3 files** (``navigation-office.esa.int``, open access),
   2016-08-01 to 2026-10-06: precise orbits *and* satellite clocks every 5 min.
   Each ~0.6-1 MB file is streamed, and only the header, the epoch lines and the
   records of E14, E18 and seven nominal Galileo control satellites are kept
   (~15x smaller; the dropped satellites can be re-fetched any time). The full
   file is kept for the first day of every month (GPS satellites for the GPS
   clock budget). About 3.4 GB are streamed, ~0.3 GB stay on disk.
2. **ESA's weekly SLR residual summaries**: per-pass mean laser-ranging
   residuals of the same orbits, our independent check on orbit errors.
3. **IERS Earth orientation parameters** (finals2000A.all) for the
   terrestrial-to-celestial rotation.
4. **Literature** PDFs from arXiv whose numbers we check against.
5. **Solar atlas and laboratory iron wavelengths** for the solar redshift
   (see ``solar.py`` for what they are).

Files land in ``data/``; nothing is re-downloaded once present, so after one
run the analysis works offline. ``data/manifest.json`` records each file's
URL, size, SHA-256 and download date; ``data/raw/esa_sp3/log.jsonl`` records
the same for every daily SP3 file, including the original (unfiltered) file's
SHA-256.
"""

from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import gzip
import hashlib
import json
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

UA = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64) astro-atlas hobby project"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

class NotFound(Exception):
    pass


def fetch(url: str, *, timeout: float = 120.0, retries: int = 5, ua: str = UA) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404, 410):
                raise NotFound(f"HTTP {exc.code} for {url}") from exc
            last = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
        time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def save(dest: Path, payload: bytes) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(dest)
    return dest


_manifest_lock = threading.Lock()


def record(dest: Path, url: str, payload: bytes, note: str = "") -> None:
    with _manifest_lock:
        man = json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {"files": {}}
        man["files"][str(dest.relative_to(C.PROJECT_DIR))] = {
            "url": url, "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
            "downloaded": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
            **({"note": note} if note else {}),
        }
        C.MANIFEST_PATH.write_text(json.dumps(man, indent=1, sort_keys=True))


def fetch_file(url: str, dest: Path, *, note: str = "", ua: str = UA) -> Path:
    if dest.exists():
        return dest
    payload = fetch(url, ua=ua)
    save(dest, payload)
    record(dest, url, payload, note)
    print(f"  {dest.relative_to(C.PROJECT_DIR)}  ({len(payload) / 1e6:.2f} MB)")
    return dest


# ---------------------------------------------------------------------------
# 1. ESA daily SP3 (orbits + 5-min clocks)
# ---------------------------------------------------------------------------

def gps_week(day: dt.date) -> int:
    return (day - C.GPS_EPOCH).days // 7


def sp3_name(series: str, day: dt.date) -> str:
    return f"{series}_{day.year}{day.timetuple().tm_yday:03d}0000_01D_05M_ORB.SP3.gz"


def sp3_url(series: str, day: dt.date) -> str:
    return C.ESA_PRODUCTS_URL.format(week=gps_week(day)) + sp3_name(series, day)


def filtered_path(series: str, day: dt.date) -> Path:
    return C.SP3_DIR / str(day.year) / f"{series}_{day.year}{day.timetuple().tm_yday:03d}.sp3.gz"


def filter_sp3(text: str, keep: tuple[str, ...]) -> tuple[str, dict]:
    """Keep the header, epoch lines and the P/V records of the chosen satellites."""
    out, info = [], {"epochs": 0, "present": []}
    in_header = True
    seen = set()
    for line in text.splitlines():
        if line.startswith("*"):
            in_header = False
            info["epochs"] += 1
            out.append(line)
        elif in_header or line.startswith("EOF"):
            out.append(line)
        elif line[:1] in "PV" and line[1:4] in keep:
            out.append(line)
            seen.add(line[1:4])
    info["present"] = sorted(seen)
    return "\n".join(out) + "\n", info


_log_lock = threading.Lock()
SP3_LOG = C.SP3_DIR / "log.jsonl"


def _log(entry: dict) -> None:
    with _log_lock:
        with SP3_LOG.open("a") as fh:
            fh.write(json.dumps(entry) + "\n")


def fetch_sp3_day(day: dt.date) -> dict:
    for series in C.ESA_SERIES:
        if filtered_path(series, day).exists():
            return {"date": day.isoformat(), "series": series, "status": "cached"}
    for series in C.ESA_SERIES:
        url = sp3_url(series, day)
        try:
            payload = fetch(url)
        except NotFound:
            continue
        text = gzip.decompress(payload).decode("ascii", errors="replace")
        filtered, info = filter_sp3(text, C.KEEP_SATS)
        dest = filtered_path(series, day)
        save(dest, gzip.compress(filtered.encode("ascii"), compresslevel=6))
        if day.day == 1:   # keep the whole file once a month (all GPS satellites etc.)
            save(C.SP3_FULL_DIR / sp3_name(series, day), payload)
        entry = {"date": day.isoformat(), "series": series, "url": url,
                 "bytes_original": len(payload),
                 "sha256_original": hashlib.sha256(payload).hexdigest(),
                 "epochs": info["epochs"], "present": info["present"],
                 "downloaded": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")}
        _log(entry)
        return {**entry, "status": "new"}
    entry = {"date": day.isoformat(), "series": None, "status": "missing"}
    _log(entry)
    return entry


def fetch_sp3(workers: int = 6) -> None:
    days = [C.SP3_START + dt.timedelta(days=k)
            for k in range((C.SP3_END - C.SP3_START).days + 1)]
    known_missing = set()
    if SP3_LOG.exists():
        for line in SP3_LOG.read_text().splitlines():
            e = json.loads(line)
            if e.get("status") == "missing" or e.get("series") is None:
                known_missing.add(e["date"])
    todo = [d for d in days
            if not any(filtered_path(s, d).exists() for s in C.ESA_SERIES)
            and d.isoformat() not in known_missing]
    print(f"ESA SP3: {len(days)} days in range, {len(todo)} to fetch")
    t0, done, new_bytes = time.time(), 0, 0
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch_sp3_day, d): d for d in todo}
        for fut in cf.as_completed(futures):
            done += 1
            try:
                res = fut.result()
                new_bytes += res.get("bytes_original", 0)
            except Exception as exc:  # network trouble: leave the day for the next run
                print(f"  {futures[fut]}: {exc}")
            if done % 100 == 0 or done == len(todo):
                rate = new_bytes / max(time.time() - t0, 1) / 1e6
                print(f"  {done}/{len(todo)} days, {new_bytes / 1e9:.2f} GB streamed, {rate:.1f} MB/s",
                      flush=True)


# ---------------------------------------------------------------------------
# 2. ESA weekly SLR residual summaries
# ---------------------------------------------------------------------------

def fetch_slr() -> None:
    first, last = gps_week(C.SP3_START), gps_week(C.SP3_END)
    n_new = 0
    for week in range(first, last + 1):
        start = C.GPS_EPOCH + dt.timedelta(weeks=week)
        for series in C.ESA_SERIES:
            name = f"{series}_{start.year}{start.timetuple().tm_yday:03d}0000_07D_00U_SLR.SUM.gz"
            dest = C.SLR_DIR / name
            if dest.exists():
                break
            try:
                fetch_file(C.ESA_PRODUCTS_URL.format(week=week) + name, dest,
                           note="ESA/ESOC weekly SLR one-way residual summary")
                n_new += 1
                break
            except NotFound:
                continue
    print(f"ESA SLR summaries: {n_new} new, {len(list(C.SLR_DIR.glob('*.gz')))} on disk")


# ---------------------------------------------------------------------------
# 3. IERS Earth orientation parameters
# ---------------------------------------------------------------------------

def fetch_eop() -> None:
    fetch_file("https://datacenter.iers.org/data/9/finals2000A.all",
               C.EOP_DIR / "finals2000A.all",
               note="IERS Rapid Service/Prediction Center, Bulletin A/B combined, IAU2000")


# ---------------------------------------------------------------------------
# 4. Literature (arXiv PDFs we read numbers from)
# ---------------------------------------------------------------------------

ARXIV = {
    "delva2018.pdf": "1812.03711",        # Galileo redshift test, PRL 121, 231101
    "herrmann2018.pdf": "1812.09161",     # Galileo redshift test, PRL 121, 231102
    "joyce2018.pdf": "1809.01240",        # Sirius B's gravitational redshift, MNRAS 481, 2361
    "gonzalez_hernandez2020.pdf": "2009.10558",  # solar redshift with HARPS + laser comb
    "reiners2016.pdf": "1511.03014",      # IAG solar flux atlas, A&A 587, A65
}


def fetch_literature() -> None:
    for name, arxiv_id in ARXIV.items():
        try:
            fetch_file(f"https://arxiv.org/pdf/{arxiv_id}", C.LIT_DIR / name,
                       note=f"arXiv:{arxiv_id}", ua=UA_BROWSER)
        except Exception as exc:  # literature is for reference only
            print(f"  {name}: {exc}")


# ---------------------------------------------------------------------------

SOURCES = {"sp3": fetch_sp3, "slr": fetch_slr, "eop": fetch_eop, "literature": fetch_literature}


def main(argv: list[str]) -> None:
    wanted = argv or list(SOURCES)
    try:
        import solar  # solar atlas + line list live with their analysis code
        SOURCES["solar"] = solar.fetch
        if not argv:
            wanted.append("solar")
    except ImportError:
        pass
    for key in wanted:
        print(f"== {key}")
        SOURCES[key]()


if __name__ == "__main__":
    main(sys.argv[1:])
