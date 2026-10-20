"""Download and cache everything the analysis needs.

Eight public sources, all free:

1. **BlackCAT** as published (Corral-Santana et al. 2016), from the CDS
   (VizieR J/A+A/587/A61): the 2016 table of dynamically confirmed black holes.
2. **BlackCAT online** (astro.puc.cl/BlackCAT), the living version: one page per
   black hole transient, each value linked to its paper on ADS.
3. **CatNoir** (Rocha et al. 2026): an independent living catalogue of 42
   Galactic black holes, served as JSON, used as a cross-check.
4. **Paulo Freire's neutron star masses** from radio pulsar timing (MPIfR).
5. **GWOSC**: the cumulative list of gravitational-wave events (GWTC-1 to
   GWTC-5.0) with source-frame masses, distances and signal-to-noise ratios.
6. **MAXI/GSC** daily light curves in 2-4, 4-10 and 10-20 keV (RIKEN/JAXA)
   for MAXI J1820+070, GX 339-4 and the Crab.
7. **Swift/BAT** hard X-ray transient monitor, 15-50 keV, for the same three.
8. **Gaia DR3** orbital solutions (``nss_two_body_orbit``) of Gaia BH1 and BH2.

Every file lands under ``data/raw/``; a manifest (``data/manifest.json``)
records where each came from, when, its size and SHA-256. Nothing already
cached is fetched again, so after one run the project works offline.

    python fetch_data.py            # fill the cache
    python fetch_data.py --refresh  # re-download everything
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from config import (
    BAT_URL, BLACKCAT_ONLINE, BLACKCAT_POLITE_DELAY_S, BLACKCAT_VIZIER,
    CATNOIR_URL, CDS_FTP, FREIRE_URL, GAIA_BH_SOURCE_IDS, GAIA_TAP_URLS,
    GWOSC_URL, MANIFEST_PATH, MAXI_URL, RAW_DIR, XRAY_SOURCES,
)

USER_AGENT = "astrophysics-atlas-I4/1.0 (hobbyist research; X-ray binaries)"
TIMEOUT_S = 180
RETRIES = 5

CDS_DIR = RAW_DIR / "blackcat_vizier"
BLACKCAT_DIR = RAW_DIR / "blackcat_online"
XRAY_DIR = RAW_DIR / "xray"
GAIA_DIR = RAW_DIR / "gaia"
for _d in (CDS_DIR, BLACKCAT_DIR, XRAY_DIR, GAIA_DIR):
    _d.mkdir(parents=True, exist_ok=True)

CATNOIR_PATH = RAW_DIR / "catnoir_BH_Catalog.json"
FREIRE_PATH = RAW_DIR / "freire_NS_masses.html"
GWOSC_PATH = RAW_DIR / "gwosc_GWTC.json"
GAIA_NSS_PATH = GAIA_DIR / "nss_two_body_orbit_gaia_bhs.csv"
BLACKCAT_INDEX_PATH = BLACKCAT_DIR / "transients.html"


def maxi_path(source: str) -> Path:
    return XRAY_DIR / f"maxi_{XRAY_SOURCES[source]['maxi']}_1day.dat"


def bat_path(source: str) -> Path:
    return XRAY_DIR / f"bat_{XRAY_SOURCES[source]['bat'].split('/')[-1]}.lc.txt"


def blackcat_page_path(object_id: int) -> Path:
    return BLACKCAT_DIR / f"info_{object_id:03d}.html"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(path: Path, url: str, source: str, **extra) -> None:
    """Add one downloaded file to the provenance manifest."""
    manifest = _load_manifest()
    data = path.read_bytes()
    manifest["files"][str(path.relative_to(RAW_DIR.parent))] = {
        "url": url,
        "source": source,
        "retrieved_utc": _now_utc(),
        "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(),
        **extra,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def http_get(url: str, *, data: dict | None = None, timeout: float = TIMEOUT_S) -> bytes:
    """GET (or POST form ``data``) with retries. DNS lookups on this machine
    fail intermittently, so a ``URLError`` is retried quickly; HTTP errors
    other than 5xx and 408 are not retried."""
    body = urllib.parse.urlencode(data).encode() if data else None
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code < 500 and err.code != 408:
                raise
            last = err
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last = err
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"giving up on {url}: {last}")


def _fetch(path: Path, url: str, source: str, *, refresh: bool, quiet: bool,
           check=None, **extra) -> Path:
    """Download ``url`` to ``path`` unless cached; ``check(bytes)`` may raise."""
    if path.exists() and path.stat().st_size > 0 and not refresh:
        return path
    raw = http_get(url)
    if check is not None:
        check(raw)
    path.write_bytes(raw)
    _record(path, url, source, **extra)
    if not quiet:
        print(f"  fetched {path.relative_to(RAW_DIR.parent)} ({len(raw) / 1e3:.0f} kB)")
    return path


def _must_contain(token: bytes):
    def check(raw: bytes) -> None:
        if token not in raw:
            raise RuntimeError(f"unexpected response (no {token!r})")
    return check


# ---------------------------------------------------------------------------
# Individual sources
# ---------------------------------------------------------------------------

def fetch_blackcat_vizier(*, refresh=False, quiet=False) -> None:
    cat = BLACKCAT_VIZIER["catalog"]
    for name in BLACKCAT_VIZIER["files"]:
        _fetch(CDS_DIR / name, f"{CDS_FTP}{cat}/{name}", BLACKCAT_VIZIER["citation"],
               refresh=refresh, quiet=quiet)


def blackcat_object_ids() -> list[int]:
    """IDs of every object listed on the online catalogue's index page."""
    html = BLACKCAT_INDEX_PATH.read_text(encoding="utf-8", errors="replace")
    return sorted({int(m) for m in re.findall(r'href="info\.php\?id=(\d+)"', html)})


def fetch_blackcat_online(*, refresh=False, quiet=False) -> None:
    _fetch(BLACKCAT_INDEX_PATH, BLACKCAT_ONLINE + "transients.php",
           "BlackCAT online (Corral-Santana et al.), index of transients",
           refresh=refresh, quiet=quiet, check=_must_contain(b"info.php?id="))
    ids = blackcat_object_ids()
    todo = [i for i in ids if refresh or not blackcat_page_path(i).exists()]
    if todo and not quiet:
        print(f"  BlackCAT online: {len(todo)} of {len(ids)} object pages to fetch")
    for i in todo:
        _fetch(blackcat_page_path(i), f"{BLACKCAT_ONLINE}info.php?id={i}",
               "BlackCAT online (Corral-Santana et al.), object page",
               refresh=True, quiet=True, check=_must_contain(b"ID#"))
        time.sleep(BLACKCAT_POLITE_DELAY_S)


def fetch_gaia_orbits(*, refresh=False, quiet=False) -> None:
    """Gaia DR3 two-body orbits of Gaia BH1 and BH2 (a two-row query)."""
    if GAIA_NSS_PATH.exists() and not refresh:
        return
    ids = ", ".join(str(s) for s in GAIA_BH_SOURCE_IDS.values())
    query = f"SELECT * FROM gaiadr3.nss_two_body_orbit WHERE source_id IN ({ids})"
    params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query}
    last = None
    for url in GAIA_TAP_URLS:
        try:
            raw = http_get(url, data=params, timeout=300)
        except Exception as err:  # try the next mirror
            last = err
            continue
        text = raw.decode("utf-8", errors="replace")
        n_rows = len([ln for ln in text.splitlines()[1:] if ln.strip()])
        if "source_id" not in text.splitlines()[0].lower() or n_rows != len(GAIA_BH_SOURCE_IDS):
            last = RuntimeError(f"{url}: unexpected answer ({n_rows} rows)")
            continue
        GAIA_NSS_PATH.write_bytes(raw)
        _record(GAIA_NSS_PATH, url, "Gaia Collaboration 2023 (DR3), nss_two_body_orbit",
                adql=query)
        if not quiet:
            print(f"  fetched Gaia DR3 orbits of {n_rows} sources from {url.split('/')[2]}")
        return
    raise RuntimeError(f"no Gaia TAP service answered: {last}")


def ensure_data(*, refresh: bool = False, quiet: bool = False) -> None:
    """Fetch whatever is missing (everything, with ``refresh``)."""
    if not quiet:
        print("Checking the data cache ...")
    fetch_blackcat_vizier(refresh=refresh, quiet=quiet)
    fetch_blackcat_online(refresh=refresh, quiet=quiet)
    _fetch(CATNOIR_PATH, CATNOIR_URL, "CatNoir v1.0 (Rocha et al. 2026, arXiv:2609.14154)",
           refresh=refresh, quiet=quiet, check=_must_contain(b'"m_bh"'))
    _fetch(FREIRE_PATH, FREIRE_URL, "P. Freire, neutron star masses from pulsar timing",
           refresh=refresh, quiet=quiet, check=_must_contain(b"individual NS mass"))
    _fetch(GWOSC_PATH, GWOSC_URL, "GWOSC event API, cumulative GWTC",
           refresh=refresh, quiet=quiet, check=_must_contain(b'"events"'))
    for source, names in XRAY_SOURCES.items():
        _fetch(maxi_path(source), MAXI_URL.format(name=names["maxi"]),
               f"MAXI/GSC daily light curve of {source} (RIKEN, JAXA, MAXI team)",
               refresh=refresh, quiet=quiet)
        _fetch(bat_path(source), BAT_URL.format(path=names["bat"]),
               f"Swift/BAT transient monitor (Krimm et al. 2013), {source}",
               refresh=refresh, quiet=quiet, check=_must_contain(b"Swift-BAT"))
    fetch_gaia_orbits(refresh=refresh, quiet=quiet)
    if not quiet:
        print("Data cache complete.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    args = parser.parse_args(argv)
    ensure_data(refresh=args.refresh)
    manifest = _load_manifest()["files"]
    total = sum(v["bytes"] for v in manifest.values())
    print(f"{len(manifest)} files, {total / 1e6:.1f} MB, listed in {MANIFEST_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
