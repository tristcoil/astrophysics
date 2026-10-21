"""Download and cache every input this project uses (about 1.6 GB, 1.1 GB of it MIST).

    python fetch_data.py              # the four clusters of the spec
    python fetch_data.py --all        # also every other HUGS cluster (+3 GB),
                                      # for the age-metallicity extension

Files land in ``data/raw/`` and are never re-downloaded once complete, so after
one run the whole analysis works offline. ``data/manifest.json`` records each
file's source, size, SHA-256 and download date. Large files are streamed to
disk and resumed with HTTP range requests if a transfer breaks. See the README
for citations.

What comes from where:

* MAST: the HUGS astro-photometric catalogues (method 1), five HST bands with
  proper-motion membership probabilities.
* mist.science: MIST v2.5 isochrones (HST ACS/WFC) and bolometric-correction
  tables (ACS/WFC and WFC3) with extinction.
* BaSTI-IAC (INAF Teramo): pre-computed HST ACS isochrones, alpha-enhanced,
  scaled-solar and helium-rich, one tarball per metallicity.
* Dartmouth (DSED): the 2012 HST ACS/WFC isochrone grid.
* McMaster: the Harris catalogue (2010 edition).
* arXiv: LaTeX tables of Baumgardt & Vasiliev (2021) distances, VandenBerg et
  al. (2013) and Valcin et al. (2025) ages, Milone et al. (2017) chromosome-map
  results and Massari et al. (2019) progenitor assignments.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import io
import json
import re
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
#: arXiv serves e-prints to browser-like agents more reliably.
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"
CHUNK = 1 << 20


def _open(url: str, *, ua: str, start: int = 0, timeout: float = 120.0):
    headers = {"User-Agent": ua}
    if start:
        headers["Range"] = f"bytes={start}-"
    return urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=timeout)


def _content_length(url: str, ua: str) -> int | None:
    req = urllib.request.Request(url, headers={"User-Agent": ua}, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            n = resp.headers.get("Content-Length")
            return int(n) if n else None
    except (urllib.error.URLError, TimeoutError):
        return None


def stream_download(url: str, dest: Path, *, ua: str = UA_TOOL, retries: int = 6) -> Path:
    """Stream ``url`` to ``dest``, resuming a partial ``.part`` file if one exists."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    total = _content_length(url, ua)
    last = None
    for attempt in range(retries):
        have = part.stat().st_size if part.exists() else 0
        if total is not None and have >= total:
            break
        try:
            with _open(url, ua=ua, start=have) as resp:
                if have and resp.status != 206:  # server ignored the range: start over
                    have = 0
                with open(part, "ab" if have else "wb") as fh:
                    t0, got = time.time(), 0
                    while block := resp.read(CHUNK):
                        fh.write(block)
                        got += len(block)
                        if total and got % (64 * CHUNK) < CHUNK:
                            done = (have + got) / total
                            rate = got / max(time.time() - t0, 1e-3) / 1e6
                            print(f"    {100 * done:5.1f}%  {rate:5.2f} MB/s", flush=True)
            if total is None or part.stat().st_size >= total:
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (403, 404, 410):
                raise RuntimeError(f"could not fetch {url}: {exc}") from exc
            time.sleep(3.0 * (attempt + 1))
    else:
        raise RuntimeError(f"could not fetch {url}: {last}")
    if total is not None and part.stat().st_size != total:
        raise RuntimeError(f"{url}: got {part.stat().st_size} of {total} bytes")
    part.replace(dest)
    return dest


def _fetch_bytes(url: str, *, ua: str = UA_TOOL, retries: int = 4) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            with _open(url, ua=ua) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (403, 404, 410):
                break
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(8 * CHUNK):
            h.update(block)
    return h.hexdigest()


def _record(manifest: dict, dest: Path, url: str, note: str) -> None:
    key = str(dest.relative_to(C.DATA_DIR))
    old = manifest["files"].get(key, {})
    size = dest.stat().st_size
    # Hashing 1 GB takes a while; reuse the stored hash if size and mtime match.
    mtime = int(dest.stat().st_mtime)
    if old.get("bytes") == size and old.get("mtime") == mtime and old.get("sha256"):
        digest = old["sha256"]
    else:
        digest = _sha256(dest)
    manifest["files"][key] = {
        "url": url,
        "note": note,
        "bytes": size,
        "mtime": mtime,
        "sha256": digest,
        "fetched": old.get("fetched") if old.get("sha256") == digest
        else dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def download(url: str, dest: Path, manifest: dict, note: str, **kw) -> Path:
    if not dest.exists():
        print(f"  downloading {url}")
        stream_download(url, dest, **kw)
    _record(manifest, dest, url, note)
    return dest


# ---------------------------------------------------------------------------
# Individual sources
# ---------------------------------------------------------------------------


def hugs_path(hugs_id: str, method: int = C.HUGS_METHOD) -> Path:
    return C.HUGS_DIR / f"hlsp_hugs_hst_wfc3-uvis-acs-wfc_{hugs_id}_multi_v1_catalog-meth{method}.txt"


def fetch_hugs(manifest: dict, every_cluster: bool) -> None:
    print("HUGS astro-photometric catalogues (Nardiello et al. 2018), photometric method 1")
    ids = [C.CLUSTERS[name]["hugs"] for name in C.CORE_CLUSTERS]
    if every_cluster:
        index = _fetch_bytes("https://archive.stsci.edu/hlsps/hugs/").decode()
        ids += sorted(set(re.findall(r'href="((?:ngc|ic|pal|e)\d+)/"', index)) - set(ids))
        ids = [i for i in ids if i != "ngc6791"]  # an old open cluster, not a globular
    for hid in ids:
        url = C.HUGS_URL.format(hugs=hid, method=C.HUGS_METHOD)
        download(url, hugs_path(hid), manifest, f"HUGS {hid} catalogue, method {C.HUGS_METHOD}")


def fetch_mist(manifest: dict) -> None:
    print("MIST v2.5 isochrones and bolometric corrections (Dotter et al. 2026)")
    download(C.MIST_ISO_URL, C.MIST_DIR / "HST_ACS_WFC.txz", manifest,
             "MIST v2.5 isochrones, HST ACS/WFC, all [Fe/H], [alpha/Fe] and both rotation rates")
    for name, url in C.MIST_BC_URLS.items():
        download(url, C.MIST_DIR / f"BC_{name}.txz", manifest,
                 f"MIST v2.5 bolometric-correction tables, {name}, A_V = 0-6 at R_V = 3.1")


def fetch_basti(manifest: dict) -> None:
    print("BaSTI-IAC pre-computed HST ACS isochrones (Hidalgo et al. 2018; Pietrinferni et al. 2021)")
    for label, grid in C.BASTI_GRIDS.items():
        for feh, node in C.BASTI_FEH_NODES.items():
            if label == "P04Y275" and feh > -0.55:
                continue  # the helium-rich grid stops at [Fe/H] = -0.6
            dest = C.BASTI_DIR / grid / f"{node}.tar.gz"
            if not dest.exists():
                option = _fetch_bytes(C.BASTI_OPTION.format(grid=grid, node=node)).decode()
                href = re.search(r'href="([^"]+\.tar\.gz)"', option)
                if not href:
                    raise RuntimeError(f"BaSTI: no ACS isochrones listed for {grid} {node}")
                url = C.BASTI_SITE + href.group(1)
            else:
                url = manifest["files"].get(str(dest.relative_to(C.DATA_DIR)), {}).get("url", "")
            download(url, dest, manifest, f"BaSTI-IAC {grid} [Fe/H] = {feh:+.2f}, HST ACS, ages 0.02-18.4 Gyr")


def fetch_dsed(manifest: dict) -> None:
    print("Dartmouth Stellar Evolution Database (Dotter et al. 2008), 2012 HST ACS/WFC grid")
    download(C.DSED_URL, C.DSED_DIR / "HST_ACSWF.tgz", manifest,
             "DSED isochrones, HST ACS/WFC (Vega), [Fe/H] -4.0 ... +0.5, [alpha/Fe] -0.2 ... +0.8")


def fetch_harris(manifest: dict) -> None:
    print("Harris (1996) catalogue of Milky Way globular clusters, 2010 edition")
    download(C.HARRIS_URL, C.LIT_DIR / "mwgc.dat", manifest, "Harris catalogue, 2010 edition")


def fetch_arxiv_tables(manifest: dict) -> None:
    print("Published tables from arXiv LaTeX sources")
    cache: dict[str, bytes] = {}
    for local, (arxiv_id, member) in C.ARXIV_SOURCES.items():
        dest = C.LIT_DIR / local
        url = C.ARXIV_EPRINT.format(id=arxiv_id)
        if not dest.exists():
            if arxiv_id not in cache:
                print(f"  downloading {url}")
                cache[arxiv_id] = _fetch_bytes(url, ua=UA_BROWSER)
            with tarfile.open(fileobj=io.BytesIO(cache[arxiv_id]), mode="r:*") as tar:
                payload = tar.extractfile(tar.getmember(member)).read()
            dest.write_bytes(payload)
        _record(manifest, dest, url, f"{member} from the arXiv:{arxiv_id} source")


def fetch_all(every_cluster: bool = False) -> dict:
    manifest = {"files": {}}
    if C.MANIFEST_PATH.exists():
        manifest = json.loads(C.MANIFEST_PATH.read_text())
    try:
        fetch_harris(manifest)
        fetch_arxiv_tables(manifest)
        fetch_hugs(manifest, every_cluster)
        fetch_dsed(manifest)
        fetch_basti(manifest)
        fetch_mist(manifest)
    finally:
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    total = sum(f["bytes"] for f in manifest["files"].values())
    print(f"{len(manifest['files'])} files cached, {total / 1e9:.2f} GB")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--all", action="store_true",
                        help="also fetch every other HUGS globular cluster (+3 GB)")
    args = parser.parse_args()
    fetch_all(every_cluster=args.all)
    sys.exit(0)
