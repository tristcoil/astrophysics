"""Download and cache every input this project uses (about 50 MB kept).

    python fetch_data.py

Files land in ``data/raw/`` and are never re-downloaded once present, so after
one run the whole analysis works offline. ``data/manifest.json`` records each
file's source, size, SHA-256 and download date. See the README for citations.

What comes from where:

* eROSITA-DE DR1: the eRASS1 cosmology sample (5,259 clusters), the primary
  cluster catalogue (12,247, with exposure, contamination probability and mass
  PDFs), the optical follow-up catalogue (richness, Galactic N_H, the optical
  depth limit ZVLIM), and the 0.2-2.3 keV exposure map (a HiPS all-sky file).
* CDS: Legacy Surveys DR10 g/r/z coverage (MOCs) and the HI4PI N_H map.
* IRSA: Planck's PSZ2 Sunyaev-Zel'dovich cluster catalogue.
* arXiv: three published figures we read numbers from (see selection.py and
  footprint.py): the selection-function curves of Ghirardini et al. (2024,
  vector PDF), and the cleaned footprint and its exposure distribution from
  Bulbul et al. (2024).
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
#: The eROSITA web server serves its tarballs to a browser-like agent.
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"


def _fetch(url: str, *, ua: str = UA_TOOL, timeout: float = 300.0, retries: int = 4) -> bytes:
    """GET ``url`` and return the body, retrying transient failures."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
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
    digest = hashlib.sha256(blob).hexdigest()
    key = str(dest.relative_to(C.DATA_DIR))
    old = manifest["files"].get(key, {})
    manifest["files"][key] = {
        "url": url,
        "note": note,
        "bytes": len(blob),
        "sha256": digest,
        "fetched": old.get("fetched") if old.get("sha256") == digest
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


def fetch_erosita_catalogues(manifest: dict) -> None:
    print("eROSITA-DE DR1 cluster catalogues (Bulbul et al. 2024; Kluge et al. 2024)")
    for label, (tgz_path, member) in C.EROSITA_FILES.items():
        dest = C.RAW_DIR / member
        url = C.EROSITA_CATALOGUES + tgz_path
        if not dest.exists():
            print(f"  downloading {url}")
            blob = _fetch(url, ua=UA_BROWSER)
            with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
                _save(dest, tar.extractfile(tar.getmember(member)).read())
        _record(manifest, dest, url, f"eRASS1 {label} catalogue (extracted from the DR1 tarball)")


def fetch_exposure_map(manifest: dict) -> None:
    print("eRASS1 0.2-2.3 keV exposure map (HiPS, order 3 all-sky file) and its sky coverage")
    for name, local, note in (
        ("Norder3/Allsky.fits", "erass1_exposure_024_Norder3_Allsky.fits",
         "eRASS1 band-024 vignetted exposure, 768 HiPS tiles of 64x64 (HEALPix order 9)"),
        ("Moc.fits", "erass1_exposure_024_Moc.fits", "sky coverage of the eROSITA-DE half of eRASS1"),
        ("properties", "erass1_exposure_024_properties.txt", "HiPS metadata"),
    ):
        download(C.EXPOSURE_HIPS + name, C.RAW_DIR / local, manifest, note, ua=UA_BROWSER)


def fetch_legacy_survey_footprint(manifest: dict) -> None:
    print("Legacy Surveys DR10 g, r, z coverage (CDS MocServer)")
    for band in "grz":
        download(C.LS_MOC.format(band=band), C.RAW_DIR / f"ls_dr10_{band}_moc8.fits", manifest,
                 f"Legacy Surveys DR10 {band}-band HiPS coverage, MOC at order 8")


def fetch_hi4pi(manifest: dict) -> None:
    print("HI4PI Galactic neutral-hydrogen column density (CDS hips2fits)")
    download(C.HI4PI_NHI, C.RAW_DIR / "hi4pi_nhi_car_0p25deg.fits", manifest,
             "HI4PI N_HI (HI4PI Collaboration 2016) resampled to a 0.25 deg equatorial CAR grid")


def fetch_psz2(manifest: dict) -> None:
    print("Planck PSZ2 union catalogue (Planck Collaboration 2016, XXVII)")
    dest = C.RAW_DIR / "HFI_PCCS_SZ-union_R2.08.fits"
    if not dest.exists():
        print(f"  downloading {C.PSZ2_URL}")
        _save(dest, gzip.decompress(_fetch(C.PSZ2_URL)))
    _record(manifest, dest, C.PSZ2_URL, "PSZ2 union catalogue, 1,653 SZ detections (gunzipped)")


def fetch_published_figures(manifest: dict) -> None:
    print("Published figures read by this project")
    dest = C.RAW_DIR / "ghirardini2024_fig4_selection.pdf"
    url = C.ARXIV_EPRINT.format(id=C.GHIRARDINI_ARXIV)
    if not dest.exists():
        print(f"  downloading {url} (keeping only {C.GHIRARDINI_FIG4})")
        with tarfile.open(fileobj=io.BytesIO(_fetch(url)), mode="r:gz") as tar:
            _save(dest, tar.extractfile(tar.getmember(C.GHIRARDINI_FIG4)).read())
    _record(manifest, dest, url, "Ghirardini et al. 2024 Fig. 4: selection-function curves (vector PDF)")
    for url, local, note in (
        (C.BULBUL_FIG2A, "bulbul2024_fig2a_footprint.png", "Bulbul et al. 2024 Fig. 2 left: cleaned footprint"),
        (C.BULBUL_FIG2B, "bulbul2024_fig2b_area_exposure.png", "Bulbul et al. 2024 Fig. 2 right: area vs exposure"),
        (C.CLERC_FIG, "clerc2024_selection_200s_2000s.png", "Clerc et al. 2024: selection curves at 200 s and 2000 s"),
    ):
        download(url, C.RAW_DIR / local, manifest, note)


def fetch_all() -> dict:
    manifest = {"files": {}}
    if C.MANIFEST_PATH.exists():
        manifest = json.loads(C.MANIFEST_PATH.read_text())
    try:
        fetch_erosita_catalogues(manifest)
        fetch_exposure_map(manifest)
        fetch_legacy_survey_footprint(manifest)
        fetch_hi4pi(manifest)
        fetch_psz2(manifest)
        fetch_published_figures(manifest)
    finally:
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    total = sum(f["bytes"] for f in manifest["files"].values())
    print(f"{len(manifest['files'])} files cached, {total / 1e6:.1f} MB")
    return manifest


if __name__ == "__main__":
    fetch_all()
    sys.exit(0)
