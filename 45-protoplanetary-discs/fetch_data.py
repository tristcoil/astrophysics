"""Download and cache everything this project uses (about 2.5 GB, nearly all of it ALMA images).

1. **DSHARP continuum images** (20 FITS files, 35-62 MB each, 0.8 GB): the "fiducial"
   1.25 mm images of the 20 discs of the ALMA Large Program DSHARP (Andrews et al. 2018),
   from its official data release at the ALMA regional centres.
2. **DSHARP radial profiles** (20 small text files): the deprojected, azimuthally averaged
   intensity profiles of DSHARP II (Huang et al. 2018), which our own profiles are checked
   against.
3. **The HD 163296 CO cube** (1.7 GB): 12CO J=2-1 channel maps, 105 channels of 0.32 km/s,
   used to weigh the star from its disc's rotation and to look for a planet's "kink".
   Skip with ``--skip-co``.
4. **HL Tau** (16 MB): the reference images from ALMA's 2014 Long Baseline Campaign,
   the first image of a ringed disc (ALMA Partnership 2015).
5. **Gaia DR3** astrometry of the 20 host stars (one ADQL query; parallaxes, RUWE,
   Bailer-Jones distances) for distances in au.
6. **NASA Exoplanet Archive** (pscomppars, ~2 MB): the known exoplanets, so the planets
   inferred from the gaps can be placed among them (as in A2).
7. **arXiv LaTeX sources** of DSHARP I, II and VII (~20 MB): their tables are parsed at
   run time, so the published numbers we compare against are never hand-typed.

Big files are fetched in parallel byte ranges with resume (single connections to the ALMA
servers stall now and then); every file's size is checked against the server's, and
``data/manifest.json`` records each URL, size, retrieval time and SHA-256.
Files already present and complete are not downloaded again (``--refresh`` forces it).

    python fetch_data.py               # everything
    python fetch_data.py --skip-co     # without the 1.7 GB CO cube
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import csv
import hashlib
import io
import json
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from config import (
    ARXIV_SOURCES,
    DSHARP_MIRRORS,
    EXOPLANET_PATH,
    EXOPLANET_TAP,
    GAIA_PATH,
    GAIA_TAP_FALLBACK,
    GAIA_TAP_SYNC,
    HLTAU,
    HLTAU_DIR,
    HLTAU_URL,
    IMAGE_DIR,
    MANIFEST_PATH,
    PAPER_DIR,
    PROFILE_DIR,
    RAW_DIR,
    TARGETS,
    ensure_dirs,
)

USER_AGENT = "Mozilla/5.0 (astrophysics-project-atlas/A5; hobby research)"
CHUNK = 1 << 20
SEGMENT_BYTES = 16 << 20  # one byte-range segment per this many bytes (up to the per-file limit)
READ_TIMEOUT_S = 60  # a stalled connection is dropped after this and resumed


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------
def _request(url: str, *, headers: dict | None = None, data: bytes | None = None) -> urllib.request.Request:
    h = {"User-Agent": USER_AGENT}
    if headers:
        h.update(headers)
    return urllib.request.Request(url, data=data, headers=h)


def remote_size(url: str) -> int:
    with urllib.request.urlopen(_request(url, headers={"Range": "bytes=0-0"}), timeout=READ_TIMEOUT_S) as r:
        rng = r.headers.get("Content-Range")
        if rng and "/" in rng:
            return int(rng.rsplit("/", 1)[1])
        return int(r.headers["Content-Length"])


def _fetch_range(url: str, part: Path, start: int, end: int, attempts: int = 8) -> None:
    """Fill ``part`` with bytes [start, end] of ``url``, resuming after stalls."""
    want = end - start + 1
    for attempt in range(1, attempts + 1):
        have = part.stat().st_size if part.exists() else 0
        if have == want:
            return
        if have > want:
            part.unlink()
            have = 0
        try:
            req = _request(url, headers={"Range": f"bytes={start + have}-{end}"})
            with urllib.request.urlopen(req, timeout=READ_TIMEOUT_S) as r, open(part, "ab") as f:
                if r.status != 206:
                    raise RuntimeError(f"server ignored the byte range (HTTP {r.status})")
                while chunk := r.read(CHUNK):
                    f.write(chunk)
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError, RuntimeError) as err:
            if attempt == attempts:
                raise RuntimeError(f"{url} bytes {start}-{end}: {err}") from err
            time.sleep(min(30.0, 2.0 * attempt))
    if part.stat().st_size != want:
        raise RuntimeError(f"{url}: segment {start}-{end} incomplete")


def download(urls: list[str], dest: Path, *, segments: int = 6) -> int:
    """Download the first working URL in ``urls`` to ``dest``; return its size in bytes."""
    last: Exception | None = None
    for url in urls:
        try:
            size = remote_size(url)
        except Exception as err:  # noqa: BLE001 - try the next mirror
            last = err
            continue
        tmp_dir = dest.parent / (dest.name + ".parts")
        tmp_dir.mkdir(parents=True, exist_ok=True)
        n = max(1, min(segments, -(-size // SEGMENT_BYTES)))
        bounds = [(i * size // n, (i + 1) * size // n - 1) for i in range(n)]
        try:
            with cf.ThreadPoolExecutor(max_workers=n) as pool:
                futures = [pool.submit(_fetch_range, url, tmp_dir / f"{i:03d}", a, b)
                           for i, (a, b) in enumerate(bounds)]
                for fut in futures:
                    fut.result()
        except Exception as err:  # noqa: BLE001 - keep the parts, try the next mirror
            last = err
            continue
        with open(dest, "wb") as out:
            for i in range(n):
                with open(tmp_dir / f"{i:03d}", "rb") as f:
                    shutil.copyfileobj(f, out, CHUNK)
        shutil.rmtree(tmp_dir)
        if dest.stat().st_size != size:
            raise RuntimeError(f"{dest.name}: got {dest.stat().st_size} bytes, server says {size}")
        return size
    raise RuntimeError(f"could not download {dest.name}: {last}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


class Manifest:
    def __init__(self) -> None:
        self.data = json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {"files": {}}

    def has(self, rel: str, path: Path) -> bool:
        entry = self.data["files"].get(rel)
        return bool(entry and path.exists() and path.stat().st_size == entry["bytes"])

    def add(self, rel: str, path: Path, url: str, description: str, **extra) -> None:
        self.data["files"][rel] = {
            "url": url,
            "description": description,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
            **extra,
        }
        self.save()

    def save(self) -> None:
        MANIFEST_PATH.write_text(json.dumps(self.data, indent=1))


def _rel(path: Path) -> str:
    return str(path.relative_to(RAW_DIR))


def fetch_file(manifest: Manifest, urls: list[str], dest: Path, description: str, refresh: bool,
               segments: int = 6) -> None:
    rel = _rel(dest)
    if not refresh and manifest.has(rel, dest):
        print(f"  cached   {rel}")
        return
    t0 = time.time()
    if not refresh and dest.exists():
        # adopt a complete file left by an earlier, interrupted run (size checked against the server)
        try:
            if dest.stat().st_size == remote_size(urls[0]):
                manifest.add(rel, dest, urls[0], description)
                print(f"  verified {rel}")
                return
        except Exception:  # noqa: BLE001
            pass
    size = download(urls, dest, segments=segments)
    manifest.add(rel, dest, urls[0], description)
    dt = time.time() - t0
    print(f"  fetched  {rel}  {size / 1e6:.1f} MB in {dt:.0f} s ({size / 1e6 / max(dt, 1e-3):.1f} MB/s)")


# ---------------------------------------------------------------------------
# The individual data sets
# ---------------------------------------------------------------------------
def fetch_dsharp(manifest: Manifest, refresh: bool, skip_co: bool) -> None:
    print("DSHARP continuum images and radial profiles (Andrews et al. 2018; Huang et al. 2018)")
    jobs = []
    for t in TARGETS:
        name = f"{t.key}_continuum.fits"
        jobs.append(([f"{m}/images/{name}" for m in DSHARP_MIRRORS], IMAGE_DIR / name,
                     f"DSHARP fiducial 1.25 mm continuum image of {t.name}"))
    # four images at a time, each in a few byte ranges
    with cf.ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(fetch_file, manifest, urls, dest, desc, refresh, 3) for urls, dest, desc in jobs]
        for fut in futures:
            fut.result()
    for t in TARGETS:
        names = [f"{t.key}.profile.txt"] if t.key != "AS205" else ["AS205_N.profile.txt", "AS205_S.profile.txt"]
        for name in names:
            fetch_file(manifest, [f"{m}/profiles/{name}" for m in DSHARP_MIRRORS], PROFILE_DIR / name,
                       f"DSHARP II radial intensity profile ({t.name})", refresh, 1)
    fetch_file(manifest, [f"{m}/profiles/README_profiles" for m in DSHARP_MIRRORS], PROFILE_DIR / "README_profiles",
               "README of the DSHARP radial profiles", refresh, 1)
    if skip_co:
        print("  skipped  HD163296_CO.fits (--skip-co)")
    else:
        print("HD 163296 12CO J=2-1 cube (DSHARP; Isella et al. 2018)")
        fetch_file(manifest, [f"{m}/images/HD163296_CO.fits" for m in DSHARP_MIRRORS], IMAGE_DIR / "HD163296_CO.fits",
                   "DSHARP fiducial 12CO J=2-1 image cube of HD 163296 (continuum subtracted)", refresh, 8)


def _safe_extract(tgz: Path, target: Path) -> list[str]:
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tgz, "r:*") as tar:
        members = [m for m in tar.getmembers() if m.isfile()]
        for m in members:
            if m.name.startswith("/") or ".." in m.name.split("/"):
                raise ValueError(f"refusing to extract suspicious path {m.name}")
        tar.extractall(target, members=members, filter="data")
    return sorted(m.name for m in members)


def fetch_hltau(manifest: Manifest, refresh: bool) -> None:
    print("HL Tau, ALMA 2014 Long Baseline Campaign reference images (ALMA Partnership 2015)")
    dest = HLTAU_DIR / Path(HLTAU_URL).name
    fetch_file(manifest, [HLTAU_URL, HLTAU_URL.replace("almascience.eso.org", "almascience.nrao.edu")], dest,
               "HL Tau Band 6+7 reference images (ADS/JAO.ALMA#2011.0.00015.SV)", refresh, 1)
    names = _safe_extract(dest, HLTAU_DIR)
    print(f"           {len(names)} files: " + ", ".join(Path(n).name for n in names))


def fetch_papers(manifest: Manifest, refresh: bool) -> None:
    print("arXiv LaTeX sources (tables parsed by literature.py)")
    for key, (arxiv_id, desc) in ARXIV_SOURCES.items():
        base = arxiv_id.split("v")[0]
        dest = PAPER_DIR / f"{base}.tar.gz"
        fetch_file(manifest, [f"https://arxiv.org/e-print/{arxiv_id}"], dest, desc, refresh, 1)
        out = PAPER_DIR / base
        if refresh or not out.exists():
            try:
                _safe_extract(dest, out)
            except tarfile.ReadError:  # single-file submissions are plain gzip
                import gzip

                out.mkdir(parents=True, exist_ok=True)
                (out / "main.tex").write_bytes(gzip.decompress(dest.read_bytes()))


def _tap_sync(url: str, query: str, fmt: str = "csv", timeout: float = 120.0) -> bytes:
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt, "QUERY": query}).encode()
    with urllib.request.urlopen(_request(url, data=data), timeout=timeout) as r:
        return r.read()


GAIA_COLUMNS = ("g.source_id, g.ra, g.dec, g.parallax, g.parallax_error, g.ruwe, g.phot_g_mean_mag, "
                "g.bp_rp, g.nu_eff_used_in_astrometry, g.pseudocolour, g.ecl_lat, g.astrometric_params_solved, "
                "g.pmra, g.pmdec")


def fetch_gaia(manifest: Manifest, refresh: bool) -> None:
    print("Gaia DR3 astrometry of the host stars")
    rel = _rel(GAIA_PATH)
    if not refresh and manifest.has(rel, GAIA_PATH):
        print(f"  cached   {rel}")
        return
    ids = [t.gaia_dr3 for t in TARGETS] + [HLTAU.gaia_dr3]
    id_list = ", ".join(str(i) for i in ids)
    query_esa = (f"SELECT {GAIA_COLUMNS}, d.r_med_geo, d.r_lo_geo, d.r_hi_geo "
                 "FROM gaiadr3.gaia_source AS g "
                 "LEFT OUTER JOIN external.gaiaedr3_distance AS d ON d.source_id = g.source_id "
                 f"WHERE g.source_id IN ({id_list})")
    try:
        raw = _tap_sync(GAIA_TAP_SYNC, query_esa)
        service = GAIA_TAP_SYNC
    except Exception as err:  # noqa: BLE001
        print(f"  ESA archive failed ({err}); trying ARI Heidelberg")
        query_ari = (f"SELECT {GAIA_COLUMNS}, d.r_med_geo, d.r_lo_geo, d.r_hi_geo "
                     "FROM gaiadr3.gaia_source AS g LEFT OUTER JOIN gaiaedr3_complements.distances AS d "
                     f"ON d.source_id = g.source_id WHERE g.source_id IN ({id_list})")
        raw = _tap_sync(GAIA_TAP_FALLBACK + "/sync", query_ari)
        service = GAIA_TAP_FALLBACK
    rows = list(csv.reader(io.StringIO(raw.decode())))
    header = [h.lower() for h in rows[0]]  # the archive's column-name case varies between responses
    by_id = {int(r[header.index("source_id")]): r for r in rows[1:]}
    missing = [i for i in ids if i not in by_id]
    if missing:
        raise RuntimeError(f"Gaia DR3 sources not returned: {missing}")
    key_of = {t.gaia_dr3: (t.key, t.tmass) for t in TARGETS + [HLTAU]}
    with open(GAIA_PATH, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["key", "tmass"] + header)
        for i in ids:
            w.writerow(list(key_of[i]) + by_id[i])
    manifest.add(rel, GAIA_PATH, service, "Gaia DR3 gaia_source rows (+ Bailer-Jones 2021 distances) of the hosts",
                 query=query_esa)
    print(f"  fetched  {rel}  ({len(ids)} stars)")


EXOPLANET_QUERY = (
    "select pl_name,hostname,discoverymethod,disc_year,pl_controv_flag,pl_orbsmax,pl_orbsmaxlim,pl_orbper,"
    "pl_bmassj,pl_bmasse,pl_bmasselim,pl_bmasse_reflink,pl_bmassprov,st_mass,st_age,sy_dist from pscomppars"
)


def fetch_exoplanets(manifest: Manifest, refresh: bool) -> None:
    print("NASA Exoplanet Archive, Planetary Systems Composite Parameters")
    rel = _rel(EXOPLANET_PATH)
    if not refresh and manifest.has(rel, EXOPLANET_PATH):
        print(f"  cached   {rel}")
        return
    url = EXOPLANET_TAP + "?" + urllib.parse.urlencode({"query": EXOPLANET_QUERY, "format": "csv"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(_request(url), timeout=180) as r:
                payload = r.read()
            break
        except (urllib.error.URLError, TimeoutError) as err:
            if attempt == 3:
                raise
            print(f"  retry ({err})")
            time.sleep(5)
    EXOPLANET_PATH.write_bytes(payload)
    n = payload.count(b"\n") - 1
    manifest.add(rel, EXOPLANET_PATH, url, "NASA Exoplanet Archive pscomppars (doi:10.26133/NEA13)",
                 query=EXOPLANET_QUERY, rows=n)
    print(f"  fetched  {rel}  ({n} planets)")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--refresh", action="store_true", help="download again even if cached")
    parser.add_argument("--skip-co", action="store_true", help="skip the 1.7 GB HD 163296 CO cube")
    args = parser.parse_args(argv)
    ensure_dirs()
    manifest = Manifest()
    t0 = time.time()
    fetch_papers(manifest, args.refresh)
    fetch_gaia(manifest, args.refresh)
    fetch_exoplanets(manifest, args.refresh)
    fetch_hltau(manifest, args.refresh)
    fetch_dsharp(manifest, args.refresh, args.skip_co)
    total = sum(e["bytes"] for e in manifest.data["files"].values())
    print(f"done in {time.time() - t0:.0f} s; {len(manifest.data['files'])} files, {total / 1e9:.2f} GB; "
          f"manifest: {MANIFEST_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
