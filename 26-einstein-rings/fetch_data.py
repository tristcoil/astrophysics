"""Download and cache everything the project needs.

About 145 MB to download (87 MB kept once the Hubble zips are unpacked), all public:

1. **The SLACS catalogues** (CDS): the 131 candidates Hubble observed, with
   redshifts, velocity dispersions and photometry (Bolton et al. 2008), the
   lens models of the 63 grade-A lenses modelled there, and stellar masses for
   all 85 grade-A lenses (Auger et al. 2009).
2. **SDSS DR17 spectra** of all 131 candidates: the spectra in which the
   lenses were found, each holding two galaxies at two redshifts.
3. **Hubble images** of the 85 grade-A lenses: cutouts of MAST's Hubble
   Advanced Products mosaics, one per filter.
4. **Two published per-lens tables** that CDS does not hold (Auger et al.
   2010; Treu et al. 2010), from the papers' arXiv LaTeX sources. The analysis
   never reads them; they grade our numbers at the end.

Everything lands in ``data/`` with ``manifest.json`` recording each file's
source URL, retrieval time and SHA-256 hash. Files already present are never
re-downloaded; ``--refresh`` forces it.

    python fetch_data.py                 # everything
    python fetch_data.py --only tables   # tables | literature | sdss | hst
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import time
import urllib.error
import urllib.request
import zipfile
from datetime import datetime, timezone

import config
from config import DATA_DIR, MANIFEST_PATH, USER_AGENT, Source

GROUPS = ("tables", "literature", "sdss", "hst")


class BadPayload(ValueError):
    """The server answered, but not with the file we asked for."""


def _check_payload(source: Source, payload: bytes) -> None:
    """Cheap check that the bytes are the expected kind of file -- not an HTML
    error page served with status 200."""
    head = payload[:4000]
    text = head.decode("latin-1", errors="replace").lstrip()
    if not payload:
        raise BadPayload(f"{source.filename}: empty response")
    if source.kind == "cds" and (text.startswith("<") or "404 Not Found" in text):
        raise BadPayload(f"{source.filename}: not a CDS table ({text[:200]!r})")
    if source.kind == "fits" and not head.startswith(b"SIMPLE"):
        raise BadPayload(f"{source.filename}: not a FITS file ({text[:200]!r})")
    if source.kind == "gzip" and not head.startswith(b"\x1f\x8b"):
        raise BadPayload(f"{source.filename}: not gzip ({text[:200]!r})")
    if source.kind == "zip" and not head.startswith(b"PK"):
        raise BadPayload(f"{source.filename}: not a zip archive ({text[:200]!r})")


def _download(source: Source, *, attempts: int = 5, timeout: float = 180.0) -> bytes:
    """Fetch one URL, retrying with backoff; raise if every attempt fails.

    HTTP 429 and 5xx answers are retried after a growing pause; other 4xx
    answers mean the request itself is wrong and fail at once.
    """
    request = urllib.request.Request(source.url, headers={"User-Agent": USER_AGENT})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = response.read()
            _check_payload(source, payload)
            return payload
        except urllib.error.HTTPError as err:
            last_error = err
            if err.code != 429 and err.code < 500:
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError, BadPayload) as err:
            last_error = err
        if attempt < attempts:
            time.sleep(3.0 * attempt)
    raise RuntimeError(f"could not download {source.url}: {last_error}")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(manifest: dict, filename: str, source: Source, payload: bytes, **extra) -> None:
    manifest["files"][filename] = {
        "url": source.url,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "what": source.what,
        "reference": source.reference,
        "used_as": "analysis input" if source.input else "published comparison only",
        **extra,
    }


def _save_manifest(manifest: dict) -> None:
    tmp = MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(manifest, indent=2))
    tmp.replace(MANIFEST_PATH)


def _write_atomic(path, payload: bytes) -> None:
    """Download to a temporary name, then rename, so an interrupted download
    never leaves a truncated file that later looks cached."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    tmp.write_bytes(payload)
    tmp.replace(path)


def fetch_sources(sources: list[Source], *, refresh: bool = False, quiet: bool = False,
                  pause: float = 0.0) -> dict:
    """Make sure every plain-file source is cached; return the manifest."""
    manifest = _load_manifest()
    for source in sources:
        path = DATA_DIR / source.filename
        if path.exists() and not refresh:
            continue
        if not quiet:
            print(f"  fetching   {source.filename} ...", end="", flush=True)
        payload = _download(source)
        _write_atomic(path, payload)
        _record(manifest, source.filename, source, payload)
        _save_manifest(manifest)
        if not quiet:
            print(f" {len(payload) / 1024:.0f} kB")
        time.sleep(pause)
    return manifest


# --------------------------------------------------------------------------
# SDSS spectra: one per SLACS candidate
# --------------------------------------------------------------------------

def sdss_sources() -> list[Source]:
    import catalog  # needs the tables, so imported only once they are cached

    sources = []
    for row in catalog.read_candidates().itertuples():
        url = config.SDSS_SPECTRUM_URL.format(plate=row.plate, mjd=row.mjd, fiber=row.fiber)
        sources.append(Source(
            f"sdss/{row.name}.fits", url,
            f"SDSS spectrum of SLACS candidate {row.name} (plate {row.plate}, MJD {row.mjd}, "
            f"fibre {row.fiber}), with the pipeline's best-fit galaxy model",
            config.SDSS_REFERENCE, kind="fits"))
    return sources


# --------------------------------------------------------------------------
# Hubble cutouts: one zip per grade-A lens, unpacked to one FITS per filter
# --------------------------------------------------------------------------

def _hst_dir(name: str):
    return config.HST_DIR / name


def hst_done_marker(name: str):
    return _hst_dir(name) / "_fetched.json"


def _filter_key(member: str) -> str | None:
    """'hst_cutout_skycell-...-ra146d7361-decn10d1146_acs_wfc_f814w.fits' -> 'acs_wfc_f814w'.

    The IR 'coarse' mosaics (0.12"/pixel) duplicate the fine ones and are skipped."""
    stem = member.rsplit("/", 1)[-1].removesuffix(".fits")
    parts = stem.split("_")
    if parts[-1] == "coarse":
        return None
    try:
        start = next(i for i, p in enumerate(parts) if p in ("acs", "wfc3", "wfpc2"))
    except StopIteration:
        return None
    return "_".join(parts[start:])


def _usable_fraction(fits_bytes: bytes) -> float:
    """Share of the cutout's pixels with real data (non-zero weight)."""
    import numpy as np
    from astropy.io import fits

    with fits.open(io.BytesIO(fits_bytes)) as hdul:
        wht = hdul["WHT"].data if "WHT" in hdul else None
        sci = hdul["SCI"].data
        good = np.isfinite(sci)
        if wht is not None:
            good &= np.isfinite(wht) & (wht > 0)
        return float(good.mean())


def fetch_hst(*, refresh: bool = False, quiet: bool = False, names: list[str] | None = None) -> dict:
    """Cut every grade-A lens out of the Hubble mosaics and keep one FITS file
    per instrument and filter. A lens near a mosaic-tile boundary can come back
    with two cutouts of the same filter; we keep the one with more real pixels."""
    import catalog

    manifest = _load_manifest()
    lenses = catalog.read_grade_a()
    if names:
        lenses = lenses[lenses["name"].isin(names)]
    for row in lenses.itertuples():
        marker = hst_done_marker(row.name)
        if marker.exists() and not refresh:
            continue
        url = config.HAPCUT_URL.format(ra=row.ra, dec=row.dec, n=config.HST_CUTOUT_PIXELS)
        source = Source(f"hst/{row.name}/cutouts.zip", url,
                        f"Hubble mosaic cutouts around {row.name}", config.HAPCUT_REFERENCE,
                        kind="zip")
        if not quiet:
            print(f"  fetching   HST cutouts of {row.name} ...", end="", flush=True)
        try:
            payload = _download(source, timeout=300.0)
        except RuntimeError as err:
            print(f" FAILED ({err})")
            continue
        kept: dict[str, tuple[float, str, bytes]] = {}
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            for member in archive.namelist():
                key = _filter_key(member)
                if key is None or not member.endswith(".fits"):
                    continue
                data = archive.read(member)
                frac = _usable_fraction(data)
                if key not in kept or frac > kept[key][0]:
                    kept[key] = (frac, member, data)
        out_dir = _hst_dir(row.name)
        out_dir.mkdir(parents=True, exist_ok=True)
        for key, (frac, member, data) in sorted(kept.items()):
            fname = f"hst/{row.name}/{key}.fits"
            _write_atomic(DATA_DIR / fname, data)
            _record(manifest, fname, source, data, archive_member=member,
                    usable_pixel_fraction=round(frac, 4))
        marker.write_text(json.dumps({"filters": sorted(kept), "zip_bytes": len(payload),
                                      "zip_sha256": hashlib.sha256(payload).hexdigest()}))
        _save_manifest(manifest)
        if not quiet:
            print(f" {len(payload) / 1024:.0f} kB, filters: {', '.join(sorted(kept)) or 'none'}")
        time.sleep(1.0)
    return manifest


# --------------------------------------------------------------------------

def sources_for(group: str) -> list[Source]:
    if group == "tables":
        return config.TABLE_SOURCES
    if group == "literature":
        return config.LITERATURE_SOURCES
    if group == "sdss":
        return sdss_sources()
    raise ValueError(group)


def missing(group: str) -> list[str]:
    if group == "hst":
        import catalog
        return [n for n in catalog.read_grade_a()["name"] if not hst_done_marker(n).exists()]
    return [s.filename for s in sources_for(group) if not (DATA_DIR / s.filename).exists()]


def fetch_all(*, refresh: bool = False, only: str | None = None, quiet: bool = False) -> None:
    for group in GROUPS:
        if only and group != only:
            continue
        if not quiet:
            print(f"[{group}]")
        if group == "hst":
            fetch_hst(refresh=refresh, quiet=quiet)
        else:
            fetch_sources(sources_for(group), refresh=refresh, quiet=quiet,
                          pause=0.3 if group == "sdss" else 0.0)


def verify_cache() -> list[str]:
    """Names of cached files whose bytes no longer match the recorded hash."""
    manifest = _load_manifest()
    changed = []
    for name, record in manifest.get("files", {}).items():
        path = DATA_DIR / name
        if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            changed.append(name)
    return changed


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download every file even if it is cached")
    parser.add_argument("--only", choices=GROUPS, help="fetch just one group of files")
    args = parser.parse_args(argv)

    try:
        fetch_all(refresh=args.refresh, only=args.only)
    except RuntimeError as err:
        print(f"\n{err}", file=sys.stderr)
        return 1

    changed = verify_cache()
    if changed:
        print(f"\nWARNING: cached files differ from their recorded hashes: {changed}")
    gaps = {g: missing(g) for g in GROUPS if not args.only or g == args.only}
    for group, names in gaps.items():
        if names:
            print(f"  still missing in {group}: {len(names)} ({', '.join(names[:6])}"
                  f"{' ...' if len(names) > 6 else ''})")
    print(f"\nDone. Data in {DATA_DIR.name}/; the analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
