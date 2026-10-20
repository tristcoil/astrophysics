"""Download and cache every input this project uses.

    python fetch_data.py              # everything
    python fetch_data.py apollo naif  # only some groups

1. **APOLLO lunar laser ranging normal points** (UC San Diego, ~0.3 MB): every
   current release, 3,846 laser round trips to the Moon from 2006 to 2020, plus
   the release notes (error scalings, detector changes).
2. **JPL DE440 from NAIF** (~45 MB): positions of the Sun, Earth and Moon
   (``de440s.bsp``, 1849-2150), the Moon's integrated orientation
   (``moon_pa_de440_200625.bpc``), its frame kernel and GM values.
3. **IERS Earth orientation** (EOP 20 C04, ~5 MB): daily pole position,
   UT1-UTC and celestial pole offsets.
4. **Documents** (~10 MB): the JPL lunar laser range model memo (Williams &
   Boggs 2020), the DE440 paper (reflector coordinates), IERS Conventions
   chapters 7-9 and their reference Fortran (whose coefficients we use), and
   the TeX source of Will's 2014 Living Review (its tables are our check).
5. **Triple-pulsar arrival times** (Zenodo 13899771, 0.5 MB): the 12,474
   Nancay times of PSR J0337+1715 published by Voisin et al. (2025), for the
   white-noise sensitivity of the strong-equivalence-principle test.
6. **Cross-checks** (~115 MB): other lunar ephemerides (EPM2021 slice, DE421,
   DE430, DE441 via Horizons), JPL's ASCII DE440, one ILRS CRD file, and GPS
   velocities near Apache Point.

The literature numbers themselves were transcribed into ``literature/*.csv``
and checked against the primary papers (see ``literature/checks/``). Those
papers were downloaded into ``data/literature/papers/`` while checking and are
not needed to run anything, except the Voisin et al. arrival times (group 5).

Files land in ``data/``; nothing is re-downloaded once present, so after one
run the analysis works offline. ``data/manifest.json`` records each file's URL,
size, SHA-256 and download date.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

UA = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64) astro-atlas hobby project"


class NotFound(Exception):
    pass


def fetch(url: str, *, timeout: float = 300.0, retries: int = 5, ua: str = UA) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (403, 404, 410):
                raise NotFound(f"HTTP {exc.code} for {url}") from exc
            last = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            last = exc
        time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def record(dest: Path, url: str, payload: bytes, note: str = "") -> None:
    man = json.loads(C.MANIFEST_PATH.read_text()) if C.MANIFEST_PATH.exists() else {"files": {}}
    man["files"][str(dest.relative_to(C.PROJECT_DIR))] = {
        "url": url, "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "downloaded": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        **({"note": note} if note else {}),
    }
    C.MANIFEST_PATH.write_text(json.dumps(man, indent=1, sort_keys=True))


def in_manifest(dest: Path) -> bool:
    if not C.MANIFEST_PATH.exists():
        return False
    return str(dest.relative_to(C.PROJECT_DIR)) in json.loads(C.MANIFEST_PATH.read_text())["files"]


def fetch_file(url: str, dest: Path, *, note: str = "", ua: str = UA) -> Path:
    if dest.exists():
        if not in_manifest(dest):   # fetched by hand while exploring: record it as it is
            record(dest, url, dest.read_bytes(), note or "recorded from an earlier download")
        return dest
    payload = fetch(url, ua=ua)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(dest)
    record(dest, url, payload, note)
    print(f"  {dest.relative_to(C.PROJECT_DIR)}  ({len(payload) / 1e6:.2f} MB)")
    return dest


def fetch_apollo() -> None:
    print("APOLLO normal points")
    for name in list(C.APOLLO_FILES.values()) + list(C.APOLLO_NOTES):
        fetch_file(C.APOLLO_URL.format(name=name), C.APOLLO_DIR / name)
    fetch_file(C.APOLLO_PAGE, C.APOLLO_DIR / "norm_pts.html", note="format description and release list")


def fetch_naif() -> None:
    print("JPL DE440 kernels (NAIF)")
    for name, path in C.NAIF_FILES.items():
        fetch_file(C.NAIF_URL.format(path=path), C.NAIF_DIR / name)


def fetch_iers() -> None:
    print("IERS EOP 20 C04")
    fetch_file(C.IERS_EOP_URL, C.IERS_DIR / C.IERS_EOP_FILE)


def fetch_documents() -> None:
    print("Documents (models, tables)")
    for rel, url in C.DOCUMENTS.items():
        ua = UA_BROWSER if "arxiv.org" in url else UA
        fetch_file(url, C.LIT_DATA_DIR / rel, ua=ua)
    # Will (2014): unpack the TeX source once (tables are read from article.tex)
    src = C.LIT_DATA_DIR / "will2014" / "src.tar.gz"
    out = src.parent / "x"
    if src.exists() and not (out / "article.tex").exists():
        out.mkdir(exist_ok=True)
        with tarfile.open(src) as tar:
            tar.extractall(out, filter="data")


def fetch_pulsar() -> None:
    """Voisin et al. (2025): the Data part of their Zenodo release (doi:10.5281/zenodo.13899771)."""
    print("Triple-pulsar arrival times (Voisin et al. 2025, Zenodo)")
    base = C.LIT_DATA_DIR / "papers" / "voisin-2025-zenodo"
    for name in ("Data.tar.bz2", "Readme.md"):
        fetch_file(f"https://zenodo.org/api/records/13899771/files/{name}/content", base / name,
                   note="Voisin et al. 2025, A&A 693, A143")
    if not (base / "x" / "Data").is_dir():
        (base / "x").mkdir(exist_ok=True)
        with tarfile.open(base / "Data.tar.bz2") as tar:
            tar.extractall(base / "x", filter="data")


def fetch_ftp_range(host: str, path: str, start: int, length: int, dest: Path, note: str) -> Path:
    """Fetch bytes [start, start+length) of a file on an anonymous FTP server."""
    if dest.exists():
        if not in_manifest(dest):
            record(dest, f"ftp://{host}{path}", dest.read_bytes(), note)
        return dest
    import ftplib

    class Enough(Exception):
        pass

    buf = bytearray()

    def take(chunk):
        buf.extend(chunk)
        if len(buf) >= length:
            raise Enough

    for attempt in range(4):
        buf.clear()
        ftp = ftplib.FTP(host, timeout=120)
        try:
            ftp.login()
            ftp.voidcmd("TYPE I")
            try:
                ftp.retrbinary(f"RETR {path}", take, rest=start)
            except Enough:
                pass
            break
        except (ftplib.all_errors) as exc:     # noqa: PERF203 - retry on network errors
            if attempt == 3:
                raise RuntimeError(f"could not fetch ftp://{host}{path}: {exc}") from exc
            time.sleep(5 * (attempt + 1))
        finally:
            try:
                ftp.close()
            except Exception:
                pass
    payload = bytes(buf[:length])
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(payload)
    record(dest, f"ftp://{host}{path} bytes {start}-{start + length - 1}", payload, note)
    print(f"  {dest.relative_to(C.PROJECT_DIR)}  ({len(payload) / 1e6:.2f} MB)")
    return dest


def fetch_crosscheck() -> None:
    """Independent lunar ephemerides for the along-track comparison (README, lunar ranging)."""
    print("Cross-check lunar ephemerides: EPM2021 slice, DE421, DE430, DE441 (Horizons)")
    fetch_ftp_range("ftp.iaaras.ru", "/pub/epm/EPM2021/DE/epm2021_de_ascii.txt", 67_000_000, 6_500_001,
                    C.RAW_DIR / "epm2021" / "epm2021_de_ascii_2005-2021.part",
                    "EPM2021 (IAA RAS), JPL-format ASCII, 2005-2021 records only (byte range)")
    fetch_file(C.NAIF_URL.format(path="spk/planets/a_old_versions/de421.bsp"), C.NAIF_DIR / "de421.bsp")
    fetch_file(C.NAIF_URL.format(path="pck/moon_pa_de421_1900-2050.bpc"), C.NAIF_DIR / "moon_pa_de421_1900-2050.bpc")
    fetch_file("https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de430/ascp1950.430",
               C.RAW_DIR / "jpl_ascii" / "ascp1950.430", note="DE430 ASCII, 1950-2050")
    # JPL's own ASCII DE440 (test_physics checks NAIF's SPK against it) and its header
    fetch_file("https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/ascp01950.440",
               C.RAW_DIR / "jpl_ascii" / "ascp01950.440", note="DE440 ASCII, 1950-2050")
    fetch_file("https://ssd.jpl.nasa.gov/ftp/eph/planets/ascii/de440/header.440", C.RAW_DIR / "jpl_ascii" / "header.440")
    # one ILRS CRD file with APOLLO normal points: its headers state the time-tag convention (ground transmit, UTC)
    fetch_file("ftp://edc.dgfi.tum.de/pub/slr/data/npt_crd/apollo15/2013/apollo15_201303.npt",
               C.RAW_DIR / "edc_crd" / "apollo15_201303.npt", note="EDC (DGFI-TUM) mirror of the ILRS archive")
    fetch_file("https://ssd.jpl.nasa.gov/ftp/eph/planets/ioms/de421_moon_coord_iom.pdf",
               C.LIT_DATA_DIR / "jpl" / "de421_moon_coord_iom.pdf", note="DE421 reflector coordinates (Table 3)")
    # DE441 geocentric Moon from JPL Horizons, every 15 days 2000-2022 (TDB)
    url = ("https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND='301'&OBJ_DATA='NO'&MAKE_EPHEM='YES'"
           "&EPHEM_TYPE='VECTORS'&CENTER='500@399'&START_TIME='2000-01-01'&STOP_TIME='2022-12-31'&STEP_SIZE='15%20d'"
           "&TIME_TYPE='TDB'&REF_PLANE='FRAME'&REF_SYSTEM='ICRF'&VEC_TABLE='2'&OUT_UNITS='KM-S'&CSV_FORMAT='YES'"
           "&VEC_LABELS='NO'")
    fetch_file(url, C.RAW_DIR / "horizons" / "moon_de441_2000_2022.txt", note="Horizons (DE441) geocentric Moon")
    # GPS velocities (Nevada Geodetic Laboratory MIDAS, IGS14): station P027 is 2.8 km from APOLLO
    fetch_file("https://geodesy.unr.edu/velocities/midas.IGS14.txt", C.RAW_DIR / "gps" / "midas.IGS14.txt",
               note="NGL MIDAS GPS velocities (Blewitt et al. 2016), IGS14")


GROUPS = {"apollo": fetch_apollo, "naif": fetch_naif, "iers": fetch_iers, "documents": fetch_documents,
          "pulsar": fetch_pulsar, "crosscheck": fetch_crosscheck}


def main(argv: list[str]) -> None:
    wanted = argv or list(GROUPS)
    for key in wanted:
        if key not in GROUPS:
            raise SystemExit(f"unknown group {key!r}; choose from {', '.join(GROUPS)}")
        GROUPS[key]()
    print("done; manifest at", C.MANIFEST_PATH.relative_to(C.PROJECT_DIR))


if __name__ == "__main__":
    main(sys.argv[1:])
