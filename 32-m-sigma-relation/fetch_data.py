"""Download and cache every input this project uses (about 3 MB in all).

    python fetch_data.py

Files land in ``data/raw/`` and are never re-downloaded once present, so after
one run the whole analysis works offline. ``data/manifest.json`` records each
file's source, size, SHA-256 and download date. See the README for citations.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import http.cookiejar
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import config as C

#: A plain tool User-Agent works for GitHub, arXiv, CDS and Oxford. The AGN
#: database's PHP session is happier with a browser-like one.
UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"

VIZIER_TSV = "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source={src}&-out.all&-out.max=unlimited"
ARXIV_EPRINT = "https://arxiv.org/e-print/{id}"


def _fetch(url: str, *, data: bytes | None = None, opener=None, ua: str = UA_TOOL,
           timeout: float = 120.0, retries: int = 4) -> bytes:
    """GET (or POST, with ``data``) ``url`` and return the body, retrying transient failures."""
    opener = opener or urllib.request.build_opener()
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": ua})
            with opener.open(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (403, 404, 410):
                break
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def _save(dest: Path, payload: bytes) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(dest)
    return dest


def _record(manifest: dict, dest: Path, url: str, note: str) -> None:
    blob = dest.read_bytes()
    key = str(dest.relative_to(C.DATA_DIR))
    old = manifest["files"].get(key, {})
    manifest["files"][key] = {
        "url": url,
        "note": note,
        "bytes": len(blob),
        "sha256": hashlib.sha256(blob).hexdigest(),
        "fetched": old.get("fetched") if old.get("sha256") == hashlib.sha256(blob).hexdigest()
        else dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def download(url: str, dest: Path, manifest: dict, note: str, **kw) -> Path:
    """Fetch ``url`` to ``dest`` unless it is already cached."""
    if not dest.exists():
        print(f"  downloading {url}")
        _save(dest, _fetch(url, **kw))
    _record(manifest, dest, url, note)
    return dest


# ---------------------------------------------------------------------------
# Individual sources
# ---------------------------------------------------------------------------


def fetch_vdb16(manifest: dict) -> None:
    print("van den Bosch (2016) black hole mass compilation")
    for name in C.VDB16_FILES:
        download(C.VDB16_REPO + name, C.RAW_DIR / "vdb16" / name, manifest,
                 "van den Bosch 2016 compilation (author's GitHub repository)")
    for label, src in C.VDB16_VIZIER.items():
        download(VIZIER_TSV.format(src=src), C.RAW_DIR / "vdb16" / f"vizier_{label}.tsv", manifest,
                 f"van den Bosch 2016 as published (VizieR {src})")


def fetch_arxiv(manifest: dict) -> None:
    """Keep only the table files from each arXiv source bundle."""
    for arxiv_id, spec in C.ARXIV.items():
        label = spec["label"]
        outdir = C.RAW_DIR / label
        wanted = [outdir / name for name in spec["keep"]]
        print(f"arXiv:{arxiv_id} ({label}) source tables")
        url = ARXIV_EPRINT.format(id=arxiv_id)
        if not all(p.exists() for p in wanted):
            print(f"  downloading {url}")
            blob = _fetch(url)
            with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
                for name in spec["keep"]:
                    member = tar.getmember(name)
                    _save(outdir / name, tar.extractfile(member).read())
        for p in wanted:
            _record(manifest, p, url, f"{label}: {p.name} from the arXiv:{arxiv_id} source bundle")


def fetch_agn_database(manifest: dict) -> None:
    """The AGN Black Hole Mass Database's CSV export, for two choices of <f>.

    The export is a two-step PHP session: post the virial factor, post "Create
    File", then download.php serves the table for that session.
    """
    print("AGN Black Hole Mass Database (Bentz & Katz 2015)")
    for label, fvalue in C.AGN_DB_F.items():
        dest = C.RAW_DIR / "agn" / f"agnmass_{label}.csv"
        url = C.AGN_DB_URL + "download.php"
        if not dest.exists():
            print(f"  exporting masses with <f> = {fvalue}")
            jar = http.cookiejar.CookieJar()
            opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))
            _fetch(C.AGN_DB_URL, opener=opener, ua=UA_BROWSER)
            post = urllib.parse.urlencode({"otherf": fvalue, "submit": "Recalculate"}).encode()
            _fetch(C.AGN_DB_URL, data=post, opener=opener, ua=UA_BROWSER)
            post = urllib.parse.urlencode({"filetype": "csv", "submit": "Create File"}).encode()
            _fetch(C.AGN_DB_URL, data=post, opener=opener, ua=UA_BROWSER)
            body = _fetch(url, opener=opener, ua=UA_BROWSER)
            text = body.decode("utf-8", "replace")
            if f"<f> value applied = {fvalue}" not in text:
                raise RuntimeError(f"AGN database did not apply <f> = {fvalue}:\n{text[:300]}")
            _save(dest, body)
        _record(manifest, dest, url, f"AGN Black Hole Mass Database export with <f> = {fvalue}")


def fetch_atlas3d(manifest: dict) -> None:
    print("ATLAS3D volume-limited early-type galaxy sample")
    for label, name in C.ATLAS3D_FILES.items():
        download(C.ATLAS3D_URL + name, C.RAW_DIR / "atlas3d" / name, manifest,
                 f"ATLAS3D {label} table (project website)")


def fetch_all() -> dict:
    manifest = {"files": {}}
    if C.MANIFEST_PATH.exists():
        manifest = json.loads(C.MANIFEST_PATH.read_text())
    try:
        fetch_vdb16(manifest)
        fetch_arxiv(manifest)
        fetch_agn_database(manifest)
        fetch_atlas3d(manifest)
    finally:
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    total = sum(f["bytes"] for f in manifest["files"].values())
    print(f"{len(manifest['files'])} files cached, {total / 1e6:.2f} MB")
    return manifest


if __name__ == "__main__":
    fetch_all()
    sys.exit(0)
