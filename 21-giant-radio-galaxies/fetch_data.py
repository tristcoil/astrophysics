"""Download and cache everything the project needs.

About 0.4 GB in total, all public:

1. **Catalogues** (CDS): the 2,060 giants Oei et al. (2023) found by eye in
   LoTSS DR2, the 11,585 giants of Mostert et al. (2024), the Cosmic Web
   densities around giants and other radio galaxies of Oei et al. (2024), and
   Andernach et al.'s (2025) re-examination of the 142 giants larger than
   3 Mpc. Plus, from VizieR's TAP service, the resolved sources of the LoTSS DR2
   value-added catalogue (Hardcastle et al. 2023): sizes, radio luminosities and
   host galaxies of ordinary radio galaxies, our comparison sample.
2. **Images** of the giants we look at ourselves: LoTSS DR2 at 144 MHz (6" and
   20"), NVSS at 1.4 GHz, and Legacy Surveys colour images of the fields and of
   each host galaxy.
3. **SDSS** (SkyServer, DR17): the spectrum, if any, of every host galaxy (for
   velocity dispersions, hence black hole masses), and the Main Galaxy Sample
   around the radio galaxies (for their environments).

Everything lands in ``data/`` with ``manifest.json`` recording each file's
source URL, retrieval time and SHA-256 hash. Files already present are never
re-downloaded; ``--refresh`` forces it. Requests go one at a time, with
pauses, and back off when a server says it is busy.

    python fetch_data.py                   # everything
    python fetch_data.py --only catalogs   # catalogues | images | sdss
"""

from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import http.client
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

import config
from config import DATA_DIR, MANIFEST_PATH, SDSS_DIR, USER_AGENT, Source

GROUPS = ("catalogs", "images", "sdss")


class BadPayload(ValueError):
    """The server answered, but not with the file we asked for."""


def _check_payload(source: Source, head: bytes) -> None:
    """Cheap check that the bytes are the expected kind of file -- not an HTML
    error page or a service's error message served with status 200."""
    text = head[:4000].decode("latin-1", errors="replace")
    if source.kind == "fits" and not head.startswith(b"SIMPLE"):
        raise BadPayload(f"{source.filename}: not a FITS file ({' '.join(text[:300].split())!r})")
    if source.kind == "jpeg" and not head.startswith(b"\xff\xd8"):
        raise BadPayload(f"{source.filename}: not a JPEG ({text[:120]!r})")
    if source.kind == "gz" and not head.startswith(b"\x1f\x8b"):
        raise BadPayload(f"{source.filename}: not gzip ({text[:120]!r})")
    if source.kind in ("csv", "text"):
        stripped = text.lstrip()
        if not stripped:
            raise BadPayload(f"{source.filename}: empty response")
        if stripped.startswith("<") or "QUERY_STATUS" in text or stripped.startswith("ERROR"):
            raise BadPayload(f"{source.filename}: service returned an error ({stripped[:300]!r})")
    if source.kind == "skyserver" and not text.startswith("#Table1"):
        raise BadPayload(f"{source.filename}: SkyServer error ({text[:300]!r})")


def _download(source: Source, dest, *, attempts: int = 6) -> tuple[int, str]:
    """Stream one URL to ``dest``; return (bytes, sha256). Retries with backoff.

    HTTP 429 and 5xx answers are retried after a growing pause; other 4xx
    answers mean the request itself is wrong and fail at once. Dropped
    connections, failed DNS lookups and transfers cut off half-way are retried
    quickly.
    """
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            sha = hashlib.sha256()
            n = 0
            with urllib.request.urlopen(request, timeout=source.timeout) as response, open(dest, "wb") as out:
                first = True
                while chunk := response.read(1 << 20):
                    if first:
                        _check_payload(source, chunk)
                        first = False
                    out.write(chunk)
                    sha.update(chunk)
                    n += len(chunk)
            if n == 0:
                raise BadPayload(f"{source.filename}: empty response")
            return n, sha.hexdigest()
        except urllib.error.HTTPError as err:
            last_error = err
            if err.code != 429 and err.code < 500:
                break
            wait = float(err.headers.get("Retry-After") or 10 * attempt)
        except (urllib.error.URLError, http.client.HTTPException, OSError, BadPayload) as err:
            last_error = err
            if isinstance(err, BadPayload) and attempt >= 2:
                break  # the same wrong answer twice: the request is the problem
            wait = 3.0 * attempt
        if attempt < attempts:
            time.sleep(wait)
    raise RuntimeError(f"could not download {source.filename} from {source.url[:160]}...: {last_error}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(filename: str, entry: dict) -> None:
    """Add one manifest entry under a lock (re-read, update, write atomically)."""
    lock = MANIFEST_PATH.with_suffix(".lock")
    with open(lock, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        manifest = _load_manifest()
        manifest["files"][filename] = entry
        tmp = MANIFEST_PATH.with_suffix(f".tmp{os.getpid()}")
        tmp.write_text(json.dumps(manifest, indent=1))
        tmp.replace(MANIFEST_PATH)
        fcntl.flock(fh, fcntl.LOCK_UN)


def _sha256(path) -> str:
    sha = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            sha.update(chunk)
    return sha.hexdigest()


def fetch_sources(sources: list[Source], *, refresh: bool = False, quiet: bool = False,
                  keep_going: bool = True) -> list[str]:
    """Make sure every source file is cached. Returns the list of failures.

    Writes are atomic (download to a temporary name unique to this process,
    then rename), so an interrupted download never leaves a truncated file that
    later looks cached.
    """
    failures = []
    for source in sources:
        path = DATA_DIR / source.filename
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not refresh:
            continue
        if source.pause:
            time.sleep(source.pause)
        if not quiet:
            print(f"  fetching   {source.filename} ...", end="", flush=True)
        tmp = path.with_name(f"{path.name}.part{os.getpid()}")
        try:
            nbytes, sha = _download(source, tmp)
        except RuntimeError as err:
            tmp.unlink(missing_ok=True)
            failures.append(str(err))
            if not quiet:
                print(" FAILED")
            if not keep_going:
                break
            continue
        tmp.replace(path)
        _record(source.filename, {
            "url": source.url,
            "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "bytes": nbytes,
            "sha256": sha,
            "what": source.what,
            "reference": source.reference,
            "used_as": "analysis input" if source.input else "comparison / validation only",
        })
        if not quiet:
            print(f" {nbytes / 1e6:.2f} MB" if nbytes > 1e6 else f" {nbytes / 1e3:.0f} kB")
    return failures


# --------------------------------------------------------------------------
# SDSS queries built from the catalogues
# --------------------------------------------------------------------------

def skyserver_source(filename: str, sql: str, what: str, *, timeout: float = 600.0) -> Source:
    url = f"{config.SKYSERVER_URL}?" + urllib.parse.urlencode({"cmd": " ".join(sql.split()), "format": "csv"})
    return Source(filename, url, what, "SDSS DR17 (Abdurro'uf et al. 2022, ApJS 259, 35) via SkyServer",
                  kind="skyserver", pause=0.5, timeout=timeout)


def spectra_sql(rows: list[tuple[int, float, float]]) -> str:
    """Every SDSS spectrum within 3" of each position, with its velocity dispersion."""
    values = ",".join(f"({i},{ra:.6f},{dec:.6f})" for i, ra, dec in rows)
    return (f"SELECT v.id, s.specObjID, s.ra, s.dec, s.z, s.zErr, s.zWarning, s.class, s.subClass, "
            f"s.velDisp, s.velDispErr, s.snMedian, s.instrument, n.distance "
            f"FROM (VALUES {values}) AS v(id, ra, dec) "
            f"CROSS APPLY dbo.fGetNearbySpecObjEq(v.ra, v.dec, {config.SPEC_MATCH_ARCMIN}) AS n "
            f"JOIN SpecObj AS s ON s.specObjID = n.specObjID")


def spectra_sources(positions: list[tuple[float, float]], tag: str) -> list[Source]:
    """Batches of spectrum look-ups. Each batch file is named after a hash of
    its positions, so a changed host list fetches only the batches it needs."""
    out = []
    step = config.SKYSERVER_BATCH
    for start in range(0, len(positions), step):
        rows = [(start + i, ra, dec) for i, (ra, dec) in enumerate(positions[start:start + step])]
        sql = spectra_sql(rows)
        h = hashlib.sha1(sql.encode()).hexdigest()[:10]
        out.append(skyserver_source(f"sdss/spec_{tag}/batch{start:06d}_{h}.csv", sql,
                                    f"SDSS spectra within 3\" of {len(rows)} {tag} host galaxies"))
    return out


MGS_SQL = """
SELECT specObjID, ra, dec, z, petroMag_r, extinction_r, dered_g, dered_r
FROM SpecPhoto
WHERE (legacy_target1 & {bits}) != 0
  AND class = 'GALAXY' AND sciencePrimary = 1 AND mode = 1 AND zWarning = 0
  AND z BETWEEN {zmin} AND {zmax}
  AND petroMag_r - extinction_r < {rlim}
  AND ra >= {ra_lo} AND ra < {ra_hi} AND dec > {dec_min}
"""

#: The LoTSS DR2 footprint lies north of declination ~+15 deg; take a margin.
MGS_DEC_MIN = 10.0
MGS_RA_STEP = 15.0


def mgs_sources() -> list[Source]:
    out = []
    edges = [i * MGS_RA_STEP for i in range(int(360 / MGS_RA_STEP) + 1)]
    for lo, hi in zip(edges[:-1], edges[1:]):
        sql = MGS_SQL.format(bits=config.MGS_TARGET_BITS, zmin=config.MGS_Z_RANGE[0],
                             zmax=config.MGS_Z_RANGE[1], rlim=config.MGS_R_LIMIT,
                             ra_lo=lo, ra_hi=hi, dec_min=MGS_DEC_MIN)
        out.append(skyserver_source(f"sdss/mgs/mgs_ra{lo:05.1f}-{hi:05.1f}.csv", sql,
                                    f"SDSS Main Galaxy Sample, RA {lo:g}-{hi:g} deg, dec > {MGS_DEC_MIN:g}"))
    return out


def host_positions() -> dict[str, list[tuple[float, float]]]:
    """Positions whose SDSS spectra we want: every giant's host, and the hosts
    of radio-loud, resolved non-giants with an SDSS redshift (the comparison)."""
    import catalogs

    giants = catalogs.all_giants()
    pos = {"giants": list(zip(giants["ra"].round(6), giants["dec"].round(6)))}
    rg = catalogs.radio_galaxies()
    # Velocity dispersions are measurable in SDSS spectra only out to z ~ 0.7.
    rg = rg[rg["has_sdss_z"] & (rg["L144"] >= config.RADIO_LOUD_W_HZ) & ~rg["is_giant_candidate"]
            & (rg["z"] < 0.7)]
    pos["radiogalaxies"] = list(zip(rg["host_ra"].round(6), rg["host_dec"].round(6)))
    return pos


def sdss_sources() -> list[Source]:
    out = mgs_sources()
    for tag, positions in host_positions().items():
        out += spectra_sources(positions, tag)
    return out


# --------------------------------------------------------------------------

def sources_for(group: str) -> list[Source]:
    if group == "catalogs":
        return config.CATALOG_SOURCES + config.vac_sources()
    if group == "images":
        return config.image_sources()
    if group == "sdss":
        return sdss_sources()
    raise KeyError(group)


def missing(group: str) -> list[str]:
    try:
        sources = sources_for(group)
    except FileNotFoundError:
        return ["(catalogues needed first)"]
    return [s.filename for s in sources if not (DATA_DIR / s.filename).exists()]


def fetch_all(*, refresh: bool = False, only: str | None = None, quiet: bool = False) -> None:
    failures = []
    for group in GROUPS:
        if only and group != only:
            continue
        sources = sources_for(group)
        todo = [s for s in sources if refresh or not (DATA_DIR / s.filename).exists()]
        if not quiet:
            print(f"{group}: {len(sources)} files, {len(todo)} to fetch")
        failures += fetch_sources(sources, refresh=refresh, quiet=quiet)
    if failures:
        raise RuntimeError("some downloads failed:\n  " + "\n  ".join(failures))


def verify_cache() -> list[str]:
    """Names of cached files whose bytes no longer match the recorded hash."""
    manifest = _load_manifest()
    return [name for name, rec in manifest.get("files", {}).items()
            if (DATA_DIR / name).exists() and _sha256(DATA_DIR / name) != rec["sha256"]]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download files even if cached")
    parser.add_argument("--only", choices=GROUPS, help="fetch one group of files")
    parser.add_argument("--verify", action="store_true", help="check cached files against their hashes")
    args = parser.parse_args(argv)
    if args.verify:
        changed = verify_cache()
        print(f"{len(changed)} cached file(s) differ from their recorded hashes" + (f": {changed}" if changed else ""))
        return 1 if changed else 0
    try:
        fetch_all(refresh=args.refresh, only=args.only)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1
    total = sum(p.stat().st_size for p in DATA_DIR.rglob("*") if p.is_file())
    print(f"\nDone. {total / 1e6:.0f} MB in {DATA_DIR.name}/; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
