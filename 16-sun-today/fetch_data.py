"""Download and cache everything the project needs for one epoch.

What is fetched, and from where (see ``config.py`` for the reasons behind each choice):

* **The snapshot**: all ten AIA channels plus an HMI magnetogram at the epoch, from
  the JSOC's near-real-time synoptic archive (Stanford). ~11 MB.
* **The week**: an AIA 193 A image and an HMI magnetogram every 3 hours for the 7
  days before the epoch, for the rotation measurement and the time-lapse. ~120 MB.
* **An HMI continuum image** (4096x4096), the most recent definitive one, through
  the Virtual Solar Observatory, for limb darkening. ~15 MB.
* **GOES X-rays**: NOAA SWPC's 7-day feeds (1-minute fluxes, flare list, active
  region list), NOAA NCEI's science-quality daily files for the same week, and the
  CHIANTI response table that turns the two X-ray channels into a temperature.
* **The flare of the week**: AIA cut-outs of the flare site in 8 channels every 3
  minutes. Each full frame is downloaded, cut, and thrown away. ~10 MB kept.

Everything lands in ``data/<epoch>/`` with a ``manifest.json`` recording the source
URL, size, checksum and retrieval time of every file. Nothing is downloaded twice,
so once the cache is complete the analysis runs offline.

    python fetch_data.py                  # the pinned epoch in config.py
    python fetch_data.py --epoch now      # the Sun as it is right now
    python fetch_data.py --epoch 2026-10-05T12:00
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import tempfile
import threading
import time
import warnings
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests
from astropy.io import fits

import config as C

USER_AGENT = "astrophysics-atlas/F2 (hobby project; every file is cached and fetched once)"
#: Parallel downloads from one archive. Two keeps the load polite (and JSOC's web
#: server stalls connections when asked for more at once).
WORKERS = 2
#: (connect, read) timeouts in seconds; a stalled transfer is abandoned and retried.
TIMEOUT = (15, 45)
RETRIES = 4


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 16), b""):
            h.update(block)
    return h.hexdigest()


def _iso(t: datetime) -> str:
    return t.astimezone(timezone.utc).isoformat(timespec="seconds")


# --- Frame quality ------------------------------------------------------------------------


def aia_frame_problem(header) -> str | None:
    """Why an AIA frame should not be used, or None if it is fine.

    Near-real-time frames carry quality bit 30 ("quicklook"), which is expected;
    any other bit flags a real problem. The eclipse and pointing flags catch the
    frames SDO takes while passing through Earth's shadow or off-point.
    """
    q = int(header.get("QUALITY", 0))
    if q & ~C.AIA_QUALITY_ALLOWED_MASK:
        return f"QUALITY={q:#x}"
    if str(header.get("ACS_ECLP", "NO")).strip() == "YES":
        return "eclipse"
    if str(header.get("ACS_MODE", "SCIENCE")).strip() != "SCIENCE":
        return f"ACS_MODE={header.get('ACS_MODE')}"
    if float(header.get("PERCENTD", 100.0) or 0.0) < 99.0:
        return f"PERCENTD={header.get('PERCENTD')}"
    if not float(header.get("EXPTIME", 0.0) or 0.0) > 0:
        return "no exposure time"
    return None


def hmi_frame_problem(header) -> str | None:
    """HMI sets the top bit of QUALITY when an observable could not be made."""
    q = int(header.get("QUALITY", 0))
    if q & 0x80000000:
        return f"QUALITY={q:#x}"
    return None


def _fits_header(path: Path):
    with fits.open(path) as hdul:
        return hdul[-1].header.copy()


# --- The fetcher --------------------------------------------------------------------------


class Fetcher:
    """Finds, downloads, checks and records files for one epoch."""

    def __init__(self, paths: C.EpochPaths, *, quiet: bool = False):
        self.p = paths.ensure()
        self.quiet = quiet
        self.session = requests.Session()
        self.session.headers["User-Agent"] = USER_AGENT
        self._lock = threading.Lock()
        self._listings: dict[str, list[str]] = {}
        if self.p.manifest.exists():
            self.manifest = json.loads(self.p.manifest.read_text())
        else:
            self.manifest = {"epoch": _iso(self.p.epoch), "created_utc": _utcnow(), "files": {}}

    # -- bookkeeping ------------------------------------------------------------------------

    def log(self, msg: str) -> None:
        if not self.quiet:
            print(msg, flush=True)

    def save_manifest(self) -> None:
        with self._lock:
            tmp = self.p.manifest.with_suffix(".json.part")
            tmp.write_text(json.dumps(self.manifest, indent=1, sort_keys=True) + "\n")
            tmp.replace(self.p.manifest)

    def record(self, path: Path, url: str, **extra) -> None:
        rel = str(path.relative_to(self.p.data))
        entry = {"url": url, "bytes": path.stat().st_size, "sha256": _sha256(path),
                 "retrieved_utc": _utcnow(), **extra}
        with self._lock:
            self.manifest["files"][rel] = entry
        self.save_manifest()  # after every file, so an interrupted run loses nothing

    def have(self, rel: str) -> bool:
        return rel in self.manifest["files"] and (self.p.data / rel).exists()

    # -- HTTP ------------------------------------------------------------------------------

    def get(self, url: str, *, timeout=TIMEOUT) -> requests.Response | None:
        """GET a small resource with retries. None for 404 (it does not exist)."""
        for attempt in range(RETRIES):
            try:
                r = self.session.get(url, timeout=timeout)
                if r.status_code == 404:
                    return None
                r.raise_for_status()
                return r
            except requests.RequestException as exc:
                if attempt == RETRIES - 1:
                    raise RuntimeError(f"giving up on {url}: {exc}") from exc
                time.sleep(3 * (attempt + 1))
        return None

    def download(self, url: str, dest: Path) -> bool:
        """Stream ``url`` to ``dest`` through a temporary file; False if it is not there.

        The whole transfer is retried, so a connection that stalls mid-file costs
        one timeout, not the run; the partial file never takes the final name.
        """
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".part")
        for attempt in range(RETRIES):
            try:
                with self.session.get(url, stream=True, timeout=TIMEOUT) as r:
                    if r.status_code == 404:
                        return False
                    r.raise_for_status()
                    with open(tmp, "wb") as fh:
                        for chunk in r.iter_content(chunk_size=1 << 16):
                            fh.write(chunk)
                tmp.replace(dest)
                return True
            except requests.RequestException as exc:
                if attempt == RETRIES - 1:
                    raise RuntimeError(f"giving up on {url}: {exc}") from exc
                time.sleep(3 * (attempt + 1))
        return False

    def listing(self, url: str) -> list[str]:
        """The file names in an Apache directory index (cached per run)."""
        with self._lock:
            if url in self._listings:
                return self._listings[url]
        r = self.get(url)
        names = [] if r is None else re.findall(r'href="([^"?/][^"]*)"', r.text)
        with self._lock:
            self._listings[url] = names
        return names

    # -- AIA near-real-time synoptic frames ------------------------------------------------

    def aia_candidates(self, t0: datetime, t1: datetime, wavelength: int) -> list[tuple]:
        """Every NRT synoptic frame of one channel between t0 and t1, as (time, url)."""
        out = []
        hour = t0.replace(minute=0, second=0, microsecond=0)
        pat = re.compile(rf"AIA(\d{{8}})_(\d{{6}})_{wavelength:04d}\.fits$")
        while hour <= t1:
            base = f"{C.AIA_NRT_ROOT}/{hour:%Y/%m/%d}/H{hour:%H}00/"
            for name in self.listing(base):
                m = pat.match(name)
                if m:
                    t = datetime.strptime(m.group(1) + m.group(2), "%Y%m%d%H%M%S").replace(
                        tzinfo=timezone.utc)
                    if t0 <= t <= t1:
                        out.append((t, base + name))
            hour += timedelta(hours=1)
        return sorted(out)

    def fetch_aia_near(self, target: datetime, wavelength: int, dest_dir: Path,
                       tolerance_min: float) -> str | None:
        """Download the good-quality frame nearest ``target``; return its relative path."""
        tol = timedelta(minutes=tolerance_min)
        cands = self.aia_candidates(target - tol, target + tol, wavelength)
        cands.sort(key=lambda c: abs(c[0] - target))
        for t, url in cands:
            dest = dest_dir / url.rsplit("/", 1)[1]
            rel = str(dest.relative_to(self.p.data))
            if self.have(rel):
                return rel
            if not self.download(url, dest):
                continue
            problem = aia_frame_problem(_fits_header(dest))
            if problem:
                self.log(f"  skipping {dest.name}: {problem}")
                dest.unlink()
                continue
            self.record(dest, url, t_obs=_iso(t))
            return rel
        return None

    # -- HMI magnetograms ------------------------------------------------------------------

    def hmi_candidates(self, t0: datetime, t1: datetime) -> list[tuple]:
        """Hourly HMI magnetograms between t0 and t1, as (time, url, definitive?).

        File times are TAI (37 s ahead of UTC), irrelevant at hourly spacing.
        """
        out = []
        pat = re.compile(r"hmi\.M_720s(_nrt)?\.(\d{8})_(\d{6})_TAI\.fits$")
        day = t0.replace(hour=0, minute=0, second=0, microsecond=0)
        while day <= t1:
            base = f"{C.HMI_FITS_ROOT}/{day:%Y/%m/%d}/"
            for name in self.listing(base):
                m = pat.match(name)
                if m:
                    t = datetime.strptime(m.group(2) + m.group(3), "%Y%m%d%H%M%S").replace(
                        tzinfo=timezone.utc)
                    if t0 <= t <= t1:
                        out.append((t, base + name, m.group(1) is None))
            day += timedelta(days=1)
        return sorted(out)

    def fetch_hmi_near(self, target: datetime, dest_dir: Path, tolerance_min: float) -> str | None:
        tol = timedelta(minutes=tolerance_min)
        cands = self.hmi_candidates(target - tol, target + tol)
        # nearest first; at equal time prefer the definitive file over the NRT one
        cands.sort(key=lambda c: (abs(c[0] - target), not c[2]))
        for t, url, definitive in cands:
            dest = dest_dir / url.rsplit("/", 1)[1]
            rel = str(dest.relative_to(self.p.data))
            if self.have(rel):
                return rel
            if not self.download(url, dest):
                continue
            problem = hmi_frame_problem(_fits_header(dest))
            if problem:
                self.log(f"  skipping {dest.name}: {problem}")
                dest.unlink()
                continue
            self.record(dest, url, t_obs=_iso(t), definitive=definitive)
            return rel
        return None

    # -- the pieces ------------------------------------------------------------------------

    def snapshot(self) -> None:
        """All ten AIA channels and a magnetogram at the epoch."""
        got = self.manifest.setdefault("snapshot", {})
        todo = [c.wavelength for c in C.CHANNELS if str(c.wavelength) not in got]
        if todo or "hmi" not in got:
            self.log(f"snapshot at {self.p.epoch:%Y-%m-%d %H:%M} UTC")
        for wl in todo:
            rel = self.fetch_aia_near(self.p.epoch, wl, self.p.snapshot, C.SNAPSHOT_TOLERANCE_MIN)
            if rel is None:
                raise RuntimeError(f"no good AIA {wl} frame within {C.SNAPSHOT_TOLERANCE_MIN} "
                                   f"min of {self.p.epoch}")
            got[str(wl)] = rel
            self.log(f"  AIA {wl:>4}: {Path(rel).name}")
        if "hmi" not in got:
            rel = self.fetch_hmi_near(self.p.epoch, self.p.snapshot, C.HMI_TOLERANCE_MIN)
            if rel is None:
                raise RuntimeError(f"no HMI magnetogram within {C.HMI_TOLERANCE_MIN} min")
            got["hmi"] = rel
            self.log(f"  HMI mag : {Path(rel).name}")
        self.save_manifest()

    def week(self) -> None:
        """An AIA image and a magnetogram every WEEK_CADENCE_HOURS for the week."""
        n = int(C.WEEK_DAYS * 24 / C.WEEK_CADENCE_HOURS)
        targets = [self.p.week_start + timedelta(hours=C.WEEK_CADENCE_HOURS * k) for k in range(n + 1)]
        wk = self.manifest.setdefault("week", {"aia": {}, "hmi": {}})
        todo_aia = [t for t in targets if _iso(t) not in wk["aia"]]
        todo_hmi = [t for t in targets if _iso(t) not in wk["hmi"]]
        if not todo_aia and not todo_hmi:
            return
        self.log(f"week: {len(targets)} times from {self.p.week_start:%Y-%m-%d %H:%M}, "
                 f"{len(todo_aia)} AIA {C.WEEK_AIA_CHANNEL} + {len(todo_hmi)} HMI to fetch")

        def one_aia(t):
            # 3-minute cadence normally; widen the search if a gap sits on the target
            try:
                for tol in (4, 30):
                    rel = self.fetch_aia_near(t, C.WEEK_AIA_CHANNEL, self.p.week_aia, tol)
                    if rel:
                        return t, rel
            except RuntimeError as exc:
                self.log(f"  {exc} (will retry on the next run)")
                return t, "retry"
            return t, None

        def one_hmi(t):
            try:
                return t, self.fetch_hmi_near(t, self.p.week_hmi, C.HMI_TOLERANCE_MIN)
            except RuntimeError as exc:
                self.log(f"  {exc} (will retry on the next run)")
                return t, "retry"

        with ThreadPoolExecutor(WORKERS) as pool:
            for t, rel in pool.map(one_aia, todo_aia):
                if rel != "retry":
                    wk["aia"][_iso(t)] = rel
            for t, rel in pool.map(one_hmi, todo_hmi):
                if rel != "retry":
                    wk["hmi"][_iso(t)] = rel
        self.save_manifest()
        missing = [k for k, v in {**wk["aia"], **wk["hmi"]}.items() if v is None]
        if missing:
            self.log(f"  {len(missing)} week slots have no usable frame (gaps are skipped)")
        unfinished = (len(targets) - len(wk["aia"])) + (len(targets) - len(wk["hmi"]))
        if unfinished:
            raise RuntimeError(f"{unfinished} week frames failed to download; run again to resume")

    def continuum(self) -> None:
        """The most recent definitive HMI continuum image, via the VSO."""
        if self.manifest.get("continuum") and self.have(self.manifest["continuum"]):
            return
        from sunpy.net import Fido
        from sunpy.net import attrs as a

        self.p.continuum.mkdir(parents=True, exist_ok=True)
        for back in range(C.CONTINUUM_SEARCH_BACK_DAYS + 1):
            day = (self.p.epoch - timedelta(days=back)).replace(hour=12, minute=0)
            if day > self.p.epoch:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = Fido.search(a.Time(day, day + timedelta(minutes=2)), a.Instrument.hmi,
                                  a.Physobs.intensity)
            if len(res) == 0 or len(res[0]) == 0:
                continue
            self.log(f"HMI continuum: fetching {day:%Y-%m-%d %H:%M} via the VSO")
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    files = Fido.fetch(res[0][:1], path=str(self.p.continuum / "{file}"),
                                       progress=False)
            except Exception as exc:  # the VSO hands off to mirrors that sometimes stall
                self.log(f"  failed ({type(exc).__name__}); trying the day before")
                continue
            if not files:
                self.log("  no file returned; trying the day before")
                continue
            path = Path(files[0])
            problem = hmi_frame_problem(_fits_header(path))
            if problem:
                self.log(f"  skipping {path.name}: {problem}")
                path.unlink()
                continue
            url = "VSO/JSOC: hmi.Ic_45s " + str(res[0][0]["Start Time"])
            self.record(path, url, t_obs=_iso(day))
            self.manifest["continuum"] = str(path.relative_to(self.p.data))
            self.save_manifest()
            return
        raise RuntimeError("no HMI continuum image found in the VSO for the last "
                           f"{C.CONTINUUM_SEARCH_BACK_DAYS} days")

    def goes(self) -> None:
        """SWPC's 7-day X-ray feeds, NCEI's science files, NOAA's regions, the response table."""
        g = self.manifest.setdefault("goes", {})
        self.p.goes.mkdir(parents=True, exist_ok=True)
        for key, url in (("xrays", C.SWPC_XRAYS_URL), ("flares", C.SWPC_FLARES_URL),
                         ("regions", C.SWPC_REGIONS_URL)):
            dest = self.p.goes / Path(url).name
            rel = str(dest.relative_to(self.p.data))
            if not self.have(rel):
                self.log(f"GOES/SWPC: {Path(url).name}")
                if not self.download(url, dest):
                    raise RuntimeError(f"SWPC feed not available: {url}")
                self.record(dest, url)
            g[key] = rel
        g["swpc_satellite"] = _swpc_satellite(self.p.data / g["xrays"])

        if "response" not in g or not self.have(g["response"]):
            dest = self.p.goes / "goes_chianti_response_latest.fits"
            for url in C.GOES_RESPONSE_URLS:
                if self.download(url, dest):
                    self.record(dest, url)
                    g["response"] = str(dest.relative_to(self.p.data))
                    break
            else:
                raise RuntimeError("could not download the GOES CHIANTI response table")

        # NCEI science files (fluxes and flare list), one per day, where archived yet
        sat = g["swpc_satellite"]
        for key_name, root, product in (("ncei", C.NCEI_XRS_ROOT, "avg1m"),
                                        ("ncei_flares", C.NCEI_FLSUM_ROOT, "flsum")):
            got = g.setdefault(key_name, {})
            day = self.p.week_start.replace(hour=0, minute=0, second=0, microsecond=0)
            while day <= self.p.epoch:
                key = f"{day:%Y%m%d}"
                if key not in got or (got[key] and not self.have(got[key])):
                    base = f"{root.format(sat=sat)}/{day:%Y/%m}/"
                    names = [n for n in self.listing(base)
                             if re.match(rf"sci_xrsf-l2-{product}_g{sat}_d{key}_v[\d-]+\.nc$", n)]
                    got[key] = None  # not archived yet (NCEI runs about a day behind)
                    if names:
                        name = sorted(names)[-1]
                        dest = self.p.goes / "ncei" / name
                        if self.download(base + name, dest):
                            self.record(dest, base + name)
                            got[key] = str(dest.relative_to(self.p.data))
                            self.log(f"GOES/NCEI: {name}")
                day += timedelta(days=1)
        self.save_manifest()

    def flare(self) -> None:
        """Pick the flare of the week, find it on the disc, and cut it out of every frame."""
        import goes
        import sdo

        fl = self.manifest.setdefault("flare", {})
        if "selection" not in fl:
            xrs = goes.load_xrs(self.p)
            events = goes.flare_list(self.p, xrs)
            sel = goes.select_flare(events, xrs, self.p)
            if sel is None:
                fl["selection"] = None
                self.save_manifest()
                self.log("no flare of class C or above this week")
                return
            fl["selection"] = sel
            self.save_manifest()
        sel = fl["selection"]
        if sel is None:
            return
        begin = datetime.fromisoformat(sel["begin"])
        peak = datetime.fromisoformat(sel["peak"])
        end = datetime.fromisoformat(sel["end"])

        # 1. Where on the Sun? Compare a hot-channel frame at the X-ray peak with one
        #    from before the flare: the flare is where the difference is largest.
        ctx = self.p.flare / "context"
        if "location" not in fl:
            self.log(f"flare {sel['class']} peaking {peak:%Y-%m-%d %H:%M}: locating it")
            pre_t = begin - timedelta(minutes=C.FLARE_BEFORE_MIN)
            ctx_files = {}
            for label, t in (("pre", pre_t), ("peak", peak)):
                for wl in (131, 94, 171):
                    rel = self.fetch_aia_near(t, wl, ctx, 3)
                    if rel is None:
                        rel = self.fetch_aia_near(t, wl, ctx, 10)
                    ctx_files[f"{label}_{wl}"] = rel
            fl["context"] = ctx_files
            loc = sdo.locate_flare(self.p.data / ctx_files["pre_131"],
                                   self.p.data / ctx_files["peak_131"])
            fl["location"] = loc
            self.save_manifest()
        loc = fl["location"]

        # 2. Cut-outs of the flare site in every channel, every 3 minutes.
        t0 = begin - timedelta(minutes=C.FLARE_BEFORE_MIN)
        t1 = end + timedelta(minutes=C.FLARE_AFTER_MIN)
        frames = fl.setdefault("frames", {})
        window = [_iso(t0), _iso(t1)]
        complete = (fl.get("window") == window
                    and all(frames.get(str(wl)) for wl in C.FLARE_CHANNELS)
                    and all(self.have(r) for rels in frames.values() for r in rels))
        if complete:  # nothing to list or fetch: this is what keeps re-runs offline
            return
        jobs = []
        for wl in C.FLARE_CHANNELS:
            have = set(frames.get(str(wl), []))
            for t, url in self.aia_candidates(t0, t1, wl):
                dest = self.p.flare / f"{wl}" / url.rsplit("/", 1)[1].replace(".fits", "_cut.fits")
                rel = str(dest.relative_to(self.p.data))
                if rel not in have or not dest.exists():
                    jobs.append((wl, t, url, dest))
        if jobs:
            self.log(f"flare cut-outs: {len(jobs)} frames in {len(C.FLARE_CHANNELS)} channels "
                     f"({t0:%H:%M}-{t1:%H:%M}); full frames are discarded after cutting")

        def one(job):
            """Fetch, check, cut and store one frame: (channel, path or None, retry?)."""
            wl, t, url, dest = job
            with tempfile.TemporaryDirectory(dir=self.p.data) as tmpd:
                full = Path(tmpd) / "full.fits"
                try:
                    if not self.download(url, full):
                        return wl, None, True
                except RuntimeError as exc:
                    self.log(f"  {exc} (will retry on the next run)")
                    return wl, None, True
                if aia_frame_problem(_fits_header(full)):
                    return wl, None, False  # a bad frame stays bad: not worth retrying
                dest.parent.mkdir(parents=True, exist_ok=True)
                sdo.save_cutout(full, dest, loc["hpc_x_arcsec"], loc["hpc_y_arcsec"],
                                C.FLARE_CUTOUT_HALF_WIDTH_ARCSEC)
            self.record(dest, url, t_obs=_iso(t), cutout_of=url.rsplit("/", 1)[1])
            return wl, str(dest.relative_to(self.p.data)), False

        retry = False
        with ThreadPoolExecutor(WORKERS) as pool:
            for wl, rel, again in pool.map(one, jobs):
                retry |= again
                if rel:
                    frames.setdefault(str(wl), []).append(rel)
        for wl in frames:
            frames[wl] = sorted(set(frames[wl]))
        if not retry:
            fl["window"] = window
        self.save_manifest()


def _swpc_satellite(xrays_json: Path) -> int:
    rows = json.loads(xrays_json.read_text())
    sats = {int(r["satellite"]) for r in rows}
    return max(sats, key=lambda s: sum(int(r["satellite"]) == s for r in rows))


def fetch_all(epoch: datetime | None = None, *, quiet: bool = False,
              skip_flare: bool = False) -> C.EpochPaths:
    """Make sure every input for ``epoch`` is cached; return its paths."""
    paths = C.paths_for(epoch or C.parse_epoch(None))
    f = Fetcher(paths, quiet=quiet)
    f.snapshot()
    f.goes()
    f.week()
    if not skip_flare:
        f.flare()
    try:
        f.continuum()
    except RuntimeError as exc:  # optional: only the limb-darkening measurement needs it
        print(f"WARNING: {exc}; limb darkening will be skipped")
    return paths


def load_manifest(paths: C.EpochPaths) -> dict:
    return json.loads(paths.manifest.read_text())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--epoch", default=None,
                        help="'now', or an ISO time (UTC); default: the pinned epoch")
    parser.add_argument("--skip-flare", action="store_true", help="do not fetch flare cut-outs")
    args = parser.parse_args(argv)
    paths = fetch_all(C.parse_epoch(args.epoch), skip_flare=args.skip_flare)
    m = load_manifest(paths)
    total = sum(v["bytes"] for v in m["files"].values())
    print(f"cache complete: {len(m['files'])} files, {total / 1e6:.0f} MB in {paths.data}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
