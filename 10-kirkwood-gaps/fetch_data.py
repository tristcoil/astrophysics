"""Download and cache everything the analysis needs.

Three sources, all free and public:

1. **The Minor Planet Center's MPCORB.DAT** -- osculating orbital elements for
   every minor planet with an orbit, around 1.5 million of them. This is the
   primary data set: every gap, family and snapshot starts here.

2. **AstDyS synthetic proper elements and families** (Knezevic & Milani;
   Milani et al. 2014). Proper elements are orbital elements with the
   planets' periodic tugs averaged out, which is what makes collisional
   families stand out. Used for the families figure and as a sharper second
   look at the gaps.

3. **JPL Horizons** -- where the planets are at the snapshot instant, for the
   top-down map and the N-body simulation's starting conditions.

Everything is cached under ``data/`` with a provenance record (URL, retrieval
time, size, SHA-256, the server's Last-Modified date). Nothing already present
is downloaded again unless ``--refresh`` is given.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import requests

from config import (
    ASTDYS_BASE,
    ASTDYS_FILES,
    DATA_DIR,
    HORIZONS_BODIES,
    HORIZONS_PATH,
    MPCORB_PATH,
    MPCORB_URL,
    PROVENANCE_PATH,
)

USER_AGENT = "astrophysics-atlas/E1 (hobby research; python-requests)"


def _load_provenance() -> dict:
    if PROVENANCE_PATH.exists():
        return json.loads(PROVENANCE_PATH.read_text())
    return {}


def _save_provenance(prov: dict) -> None:
    PROVENANCE_PATH.write_text(json.dumps(prov, indent=2, sort_keys=True))


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, *, refresh: bool = False, quiet: bool = False) -> dict:
    """Stream ``url`` to ``dest`` and return its provenance record.

    Writes to a ``.part`` file first and renames on success, so an interrupted
    download never leaves a truncated file that looks complete. Text files are
    requested with gzip transfer encoding (``requests`` asks for it and
    decodes it transparently), which roughly quarters what the archive sends.
    """
    prov = _load_provenance()
    key = dest.name
    if dest.exists() and not refresh and key in prov:
        if not quiet:
            print(f"  cached: {dest.name} ({dest.stat().st_size / 1e6:.1f} MB, "
                  f"retrieved {prov[key]['retrieved_utc']})")
        return prov[key]

    if not quiet:
        print(f"  downloading {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=120,
                      headers={"User-Agent": USER_AGENT}) as r:
        r.raise_for_status()
        # Guard against HTML error or bot-challenge pages served with 200 OK.
        ctype = r.headers.get("Content-Type", "")
        if "html" in ctype and not url.endswith(".txt"):
            raise RuntimeError(f"{url} returned HTML ({ctype}), not data")
        n = 0
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
                n += len(chunk)
                if not quiet and n % (20 << 20) < (1 << 20):
                    print(f"    {n / 1e6:6.0f} MB", flush=True)
        headers = {k: r.headers.get(k) for k in ("Last-Modified", "ETag", "Content-Type")}
    tmp.replace(dest)

    record = {
        "url": url,
        "path": str(dest.relative_to(DATA_DIR.parent)),
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "bytes": dest.stat().st_size,
        "sha256": _sha256(dest),
        "last_modified": headers["Last-Modified"],
        "etag": headers["ETag"],
    }
    prov[key] = record
    _save_provenance(prov)
    if not quiet:
        print(f"    done: {record['bytes'] / 1e6:.1f} MB")
    return record


def mpcorb_header(path: Path = MPCORB_PATH) -> list[str]:
    """The MPC's own header lines (everything above the row of dashes)."""
    lines = []
    with gzip.open(path, "rt", encoding="latin-1") as f:
        for line in f:
            if line.startswith("-----"):
                break
            lines.append(line.rstrip("\n"))
    return lines


def fetch_mpcorb(*, refresh: bool = False, quiet: bool = False) -> dict:
    """The MPC orbit file, ~94 MB gzipped."""
    return download(MPCORB_URL, MPCORB_PATH, refresh=refresh, quiet=quiet)


def fetch_astdys(*, refresh: bool = False, quiet: bool = False) -> dict:
    """AstDyS proper elements and family classification, ~135 MB."""
    out = {}
    for name in ASTDYS_FILES:
        out[name] = download(ASTDYS_BASE + name, DATA_DIR / name,
                             refresh=refresh, quiet=quiet)
    return out


def snapshot_jd() -> float:
    """The snapshot instant: 0h TT on the day MPCORB was retrieved, as a JD.

    Tied to the retrieval date rather than to "now", so re-running the
    analysis later from the same cache reproduces the same figures.
    """
    from astropy.time import Time

    prov = _load_provenance()
    if MPCORB_PATH.name not in prov:
        raise FileNotFoundError("MPCORB not fetched yet; run fetch_data.py first")
    day = prov[MPCORB_PATH.name]["retrieved_utc"][:10]
    return float(Time(day + "T00:00:00", scale="tt").jd)


def fetch_horizons(*, refresh: bool = False, quiet: bool = False) -> dict:
    """Heliocentric ecliptic (J2000) states and osculating elements of the
    planets, Mercury to Saturn, at the snapshot instant."""
    if HORIZONS_PATH.exists() and not refresh:
        data = json.loads(HORIZONS_PATH.read_text())
        if not quiet:
            print(f"  cached: {HORIZONS_PATH.name} (epoch JD {data['epoch_jd_tdb']})")
        return data

    from astroquery.jplhorizons import Horizons

    jd = snapshot_jd()  # TT; TDB differs by < 2 ms, irrelevant here
    if not quiet:
        print(f"  querying JPL Horizons for {len(HORIZONS_BODIES)} planets at JD {jd}")
    bodies = {}
    for name, hid in HORIZONS_BODIES.items():
        obj = Horizons(id=hid, location="500@10", epochs=[jd])
        vec = obj.vectors(refplane="ecliptic")
        obj = Horizons(id=hid, location="500@10", epochs=[jd])
        el = obj.elements(refplane="ecliptic")
        bodies[name] = {
            "horizons_id": hid,
            "targetname": str(vec["targetname"][0]),
            # position (AU) and velocity (AU/day), heliocentric, ecliptic J2000
            "xyz_au": [float(vec[c][0]) for c in ("x", "y", "z")],
            "vxyz_au_per_day": [float(vec[c][0]) for c in ("vx", "vy", "vz")],
            "a_au": float(el["a"][0]),
            "e": float(el["e"][0]),
            "incl_deg": float(el["incl"][0]),
            "node_deg": float(el["Omega"][0]),
            "peri_deg": float(el["w"][0]),
            "mean_anomaly_deg": float(el["M"][0]),
            "period_days": float(el["P"][0]),
        }
    data = {
        "epoch_jd_tdb": jd,
        "center": "500@10 (Sun centre)",
        "refplane": "ecliptic and mean equinox of J2000",
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": "JPL Horizons via astroquery.jplhorizons",
        "bodies": bodies,
    }
    HORIZONS_PATH.write_text(json.dumps(data, indent=2))
    return data


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Download and cache E1's data.")
    parser.add_argument("--refresh", action="store_true",
                        help="download again even if cached (the MPC file changes daily)")
    parser.add_argument("--skip-astdys", action="store_true",
                        help="only the MPC file and Horizons")
    args = parser.parse_args(argv)

    print("Minor Planet Center orbit file")
    fetch_mpcorb(refresh=args.refresh)
    for line in mpcorb_header()[:3]:
        if line.strip():
            print(f"    | {line.strip()}")

    if not args.skip_astdys:
        print("AstDyS proper elements and families")
        fetch_astdys(refresh=args.refresh)

    print("JPL Horizons planet positions")
    fetch_horizons(refresh=args.refresh)

    total = sum(p.stat().st_size for p in DATA_DIR.iterdir() if p.is_file())
    print(f"\nCache: {DATA_DIR}  ({total / 1e6:.0f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
