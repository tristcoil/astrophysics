"""Download and cache everything the analysis needs (except the live snapshot,
which fetch_live.py takes).

    python fetch_data.py            # ~330 MB on first run, then a no-op

What it fetches:
  * ZTF SN Ia DR2: 6 tables and 3,592 light curves, read by byte range out of the
    release's 1.37 GB zip (~130 MB transferred). Light curves are stored gzipped.
  * Pantheon+: the distance table and covariance D1 fitted (copied from D1's
    cache when its checksum matches D1's manifest, otherwise downloaded), the
    SALT2 model the Pantheon+ team trained (D1's light-curve model), and 40
    of their Foundation light curves for a convention check.
  * The SFD Milky Way dust maps (2 x 67 MB).
  * sncosmo's SALT3 model and the ZTF g, r, i passbands, into data/sncosmo.

Every file's source, size, SHA-256 and retrieval time go into data/manifest.json.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import requests

import config
from remotezip import USER_AGENT, RemoteZip


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def load_manifest() -> dict:
    if config.MANIFEST_PATH.exists():
        return json.loads(config.MANIFEST_PATH.read_text())
    return {}


def save_manifest(manifest: dict) -> None:
    tmp = config.MANIFEST_PATH.with_suffix(".tmp")
    tmp.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    tmp.replace(config.MANIFEST_PATH)


def download(url: str, dest: Path, retries: int = 5) -> None:
    """Stream a URL to disk atomically, with retries."""
    delay = 2.0
    for attempt in range(retries):
        try:
            with requests.get(url, stream=True, timeout=120,
                              headers={"User-Agent": USER_AGENT,
                                       "Accept-Encoding": "identity"}) as r:
                r.raise_for_status()
                tmp = dest.with_suffix(dest.suffix + ".part")
                with open(tmp, "wb") as f:
                    for block in r.iter_content(1 << 20):
                        f.write(block)
                expected = r.headers.get("Content-Length")
                if expected is not None and tmp.stat().st_size != int(expected):
                    raise requests.ConnectionError("truncated download")
                tmp.replace(dest)
                return
        except (requests.ConnectionError, requests.Timeout) as err:
            if attempt == retries - 1:
                raise
            print(f"    retrying after {err!r}")
            time.sleep(delay)
            delay *= 2


# ---------------------------------------------------------------------------
# ZTF SN Ia DR2
# ---------------------------------------------------------------------------

def fetch_ztf_dr2(manifest: dict) -> None:
    tables_done = all((config.ZTF_DIR / t).exists() for t in config.ZTF_DR2_TABLES)
    entry = manifest.get("ztf_dr2", {})
    n_expected = entry.get("n_lightcurves")
    n_have = len(list(config.ZTF_LC_DIR.glob("*_lc.csv.gz")))
    if tables_done and n_expected and n_have == n_expected:
        print(f"ZTF DR2: cached ({n_have} light curves)")
        return

    print("ZTF DR2: reading the release's table of contents ...")
    rz = RemoteZip.open(config.ZTF_DR2_ZIP)
    archive = {
        "url": config.ZTF_DR2_ZIP,
        "filename": rz.remote.headers.get("Content-Disposition", ""),
        "bytes": rz.remote.size,
        "etag": rz.remote.headers.get("ETag"),
        "last_modified": rz.remote.headers.get("Last-Modified"),
        "n_members": len(rz.infos),
    }
    print(f"  archive {rz.remote.size / 1e9:.2f} GB, {len(rz.infos)} members, "
          f"last modified {archive['last_modified']}")

    # Tables: small, fetched one by one.
    tables = {}
    for rel in config.ZTF_DR2_TABLES:
        dest = config.ZTF_DIR / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        name = config.ZTF_DR2_PREFIX + rel
        info = rz.infos[name]
        if not dest.exists() or dest.stat().st_size != info.file_size:
            for _, raw in rz.fetch([name]):
                tmp = dest.with_suffix(dest.suffix + ".part")
                tmp.write_bytes(raw)
                tmp.replace(dest)
        tables[rel] = {"bytes": dest.stat().st_size, "sha256": sha256(dest),
                       "crc32_in_zip": f"{info.CRC:08x}", "retrieved_utc": _now()}
        print(f"  {rel}: {dest.stat().st_size / 1e6:.2f} MB")

    # Light curves: fetched in contiguous runs, stored gzipped, resumable.
    lc_names = rz.names(prefix=config.ZTF_DR2_PREFIX + "lightcurves/", suffix="_lc.csv")
    todo = [n for n in lc_names
            if not (config.ZTF_LC_DIR / (Path(n).name + ".gz")).exists()]
    print(f"  light curves: {len(lc_names)} in the release, {len(todo)} to fetch "
          f"({sum(rz.infos[n].compress_size for n in todo) / 1e6:.0f} MB compressed)")
    t_start = time.time()

    def progress(done, total):
        rate = done / max(time.time() - t_start, 1e-6)
        print(f"    {done}/{total} light curves ({rate:.0f}/s)", flush=True)

    for name, raw in rz.fetch(todo, chunk_bytes=8 << 20, progress=progress):
        dest = config.ZTF_LC_DIR / (Path(name).name + ".gz")
        tmp = dest.with_suffix(".part")
        with gzip.open(tmp, "wb", compresslevel=6) as f:
            f.write(raw)
        tmp.replace(dest)

    crcs = "".join(f"{Path(n).name}:{rz.infos[n].CRC:08x}\n" for n in lc_names)
    manifest["ztf_dr2"] = {
        "archive": archive,
        "tables": tables,
        "n_lightcurves": len(lc_names),
        "lightcurves_raw_bytes": sum(rz.infos[n].file_size for n in lc_names),
        # One hash over the archive's own CRC-32 of every light curve: changes
        # if any light curve in the release changes.
        "lightcurves_crc_digest": hashlib.sha256(crcs.encode()).hexdigest(),
        "http_requests": rz.remote.n_requests,
        "retrieved_utc": _now(),
    }
    save_manifest(manifest)
    print(f"  done: {rz.remote.n_requests} HTTP requests")


# ---------------------------------------------------------------------------
# Pantheon+ (D1's data)
# ---------------------------------------------------------------------------

def fetch_pantheon(manifest: dict) -> None:
    d1_manifest = {}
    d1_path = config.D1_DIR / "data" / "manifest.json"
    if d1_path.exists():
        d1_manifest = json.loads(d1_path.read_text())
    entries = manifest.setdefault("pantheon", {})
    for name, url in config.PANTHEON_FILES.items():
        dest = config.PANTHEON_DIR / name
        if dest.exists() and name in entries:
            print(f"Pantheon+ {name}: cached")
            continue
        d1_copy = config.D1_DIR / "data" / name
        d1_sha = d1_manifest.get(name, {}).get("sha256")
        how = "downloaded"
        if d1_copy.exists() and d1_sha and sha256(d1_copy) == d1_sha:
            shutil.copyfile(d1_copy, dest)
            how = "copied from D1's verified cache"
        else:
            download(url, dest)
        entries[name] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha256(dest),
                         "retrieved_utc": _now(), "how": how}
        print(f"Pantheon+ {name}: {how} ({dest.stat().st_size / 1e6:.1f} MB)")
    save_manifest(manifest)


def fetch_salt2_b21(manifest: dict) -> None:
    """D1's light-curve model (SALT2 retrained for Pantheon+), gunzipped into
    the layout sncosmo's SALT2Source reads."""
    entries = manifest.setdefault("salt2_b21", {})
    config.SALT2_B21_DIR.mkdir(parents=True, exist_ok=True)
    for name in config.SALT2_B21_FILES:
        dest = config.SALT2_B21_DIR / name.removesuffix(".gz")
        if dest.exists() and dest.name in entries:
            continue
        url = config.SALT2_B21_RAW + name
        tmp = config.SALT2_B21_DIR / (name + ".download")
        download(url, tmp)
        if name.endswith(".gz"):
            with gzip.open(tmp, "rb") as f_in, open(dest, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            tmp.unlink()
        else:
            tmp.replace(dest)
        entries[dest.name] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha256(dest),
                              "retrieved_utc": _now()}
    save_manifest(manifest)
    print(f"SALT2-B21 model: {len(entries)} files in "
          f"{config.SALT2_B21_DIR.relative_to(config.PROJECT_DIR)}")


def fetch_foundation_lcs(manifest: dict) -> None:
    """Pantheon+'s own Foundation light curves for the first N Foundation
    supernovae (alphabetical) that are in the Pantheon+ table."""
    import pandas as pd

    entries = manifest.setdefault("foundation_lcs", {})
    config.FOUNDATION_LC_DIR.mkdir(parents=True, exist_ok=True)
    if len(entries) >= config.N_FOUNDATION_CHECK:
        print(f"Foundation light curves: cached ({len(entries)})")
        return
    tree = requests.get("https://api.github.com/repos/PantheonPlusSH0ES/DataRelease/git/trees/"
                        "main?recursive=1", timeout=120, headers={"User-Agent": USER_AGENT}).json()
    files = sorted(x["path"].rsplit("/", 1)[1] for x in tree["tree"]
                   if "/Foundation_DJ17/" in x["path"] and x["path"].endswith(".txt"))
    table = pd.read_csv(config.PANTHEON_DIR / "Pantheon+SH0ES.dat", sep=r"\s+")
    found_ids = set(table.loc[table["IDSURVEY"] == 150, "CID"].astype(str))
    for name in files:
        if len(entries) >= config.N_FOUNDATION_CHECK:
            break
        dest = config.FOUNDATION_LC_DIR / name
        if name in entries:
            continue
        download(config.FOUNDATION_LC_RAW + name, dest)
        text = dest.read_text()
        snid = next((l.split()[1] for l in text.splitlines() if l.startswith("SNID:")), "")
        if snid not in found_ids:
            dest.unlink()
            continue
        entries[name] = {"url": config.FOUNDATION_LC_RAW + name, "snid": snid,
                         "bytes": dest.stat().st_size, "sha256": sha256(dest),
                         "retrieved_utc": _now()}
    save_manifest(manifest)
    print(f"Foundation light curves: {len(entries)}")


# ---------------------------------------------------------------------------
# SFD dust maps
# ---------------------------------------------------------------------------

def fetch_sfd(manifest: dict) -> None:
    entries = manifest.setdefault("sfd", {})
    for name in config.SFD_FILES:
        dest = config.DUST_DIR / name
        if dest.exists() and name in entries:
            print(f"SFD {name}: cached")
            continue
        url = config.SFD_RAW + name
        print(f"SFD {name}: downloading ...")
        download(url, dest)
        entries[name] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha256(dest),
                         "retrieved_utc": _now()}
    save_manifest(manifest)


# ---------------------------------------------------------------------------
# sncosmo model and passbands
# ---------------------------------------------------------------------------

def fetch_sncosmo(manifest: dict) -> None:
    import sncosmo

    sncosmo.conf.data_dir = str(config.SNCOSMO_DIR)
    source = sncosmo.get_source(config.SALT_SOURCE, version=config.SALT_VERSION)
    for band in (*config.ZTF_BANDS, "bessellb", "ps1::g", "ps1::r", "ps1::i", "ps1::z"):
        sncosmo.get_bandpass(band)
    files = sorted(p for p in config.SNCOSMO_DIR.rglob("*") if p.is_file())
    manifest["sncosmo"] = {
        "version": sncosmo.__version__,
        "source": f"{source.name} v{source.version}",
        "files": {str(p.relative_to(config.SNCOSMO_DIR)):
                  {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in files},
        "retrieved_utc": _now(),
    }
    save_manifest(manifest)
    print(f"sncosmo {sncosmo.__version__}: {source.name} v{source.version}, "
          f"{len(files)} files in {config.SNCOSMO_DIR.relative_to(config.PROJECT_DIR)}")


def main() -> int:
    manifest = load_manifest()
    fetch_pantheon(manifest)
    fetch_salt2_b21(manifest)
    fetch_foundation_lcs(manifest)
    fetch_sncosmo(manifest)
    fetch_ztf_dr2(manifest)
    fetch_sfd(manifest)
    total = sum(p.stat().st_size for p in config.DATA_DIR.rglob("*") if p.is_file())
    print(f"data/ now holds {total / 1e6:.0f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
