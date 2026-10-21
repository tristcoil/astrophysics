"""Download and cache every input this project uses.

    python fetch_data.py

1. **Dyson, Eddington & Davidson 1920** (Phil. Trans. A 220, 291): the 45 page
   images of the Royal Society's scan on the Internet Archive. The measured
   star displacements were transcribed by hand from these pages into
   ``transcribed/`` (committed); the images are kept so anyone can check.
2. **NASA ADS scans**: Campbell & Trumpler 1923 and 1928 (Lick Observatory
   Bulletins 346 and 397, the 1922 Wallal eclipse) and Harvey 1979 (the
   Observatory 99, 195: the Sobral plates remeasured).
3. **Gaia DR3** (VizieR I/355/gaiadr3): every star brighter than G = 10 within
   2.5 degrees of the 1919 field, which includes all 13 stars of Table I.
4. **JPL Horizons**: the Sun and the Moon seen from Sobral and Principe around
   totality in 1919; Cassini and the Earth around the Sun during the 2002
   solar conjunction.
5. **A sky image** of the 1919 field (DSS2 colour via CDS hips2fits).

Files land in ``data/raw/``; nothing is re-downloaded once present, so after
one run the analysis works offline. ``data/manifest.json`` records each file's
URL, size, SHA-256 and download date.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import io
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import config as C

UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _fetch(url: str, *, timeout: float = 300.0, retries: int = 4, ua: str = UA_TOOL) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (400, 403, 404, 410):
                raise RuntimeError(f"HTTP {exc.code} for {url}") from exc
            time.sleep(5.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def _save(dest: Path, payload: bytes) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(dest)
    return dest


def _load_manifest() -> dict:
    if C.MANIFEST_PATH.exists():
        return json.loads(C.MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(manifest: dict, dest: Path, url: str, note: str) -> None:
    h = hashlib.sha256(dest.read_bytes()).hexdigest()
    key = str(dest.relative_to(C.DATA_DIR))
    old = manifest["files"].get(key, {})
    manifest["files"][key] = {
        "url": url, "note": note, "bytes": dest.stat().st_size, "sha256": h,
        "fetched": old["fetched"] if old.get("sha256") == h
        else dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def download(url: str, dest: Path, manifest: dict, note: str, **kw) -> Path:
    if not dest.exists():
        print(f"  downloading {url}")
        _save(dest, _fetch(url, **kw))
    _record(manifest, dest, url, note)
    return dest


# ---------------------------------------------------------------------------
# 1-2. The papers
# ---------------------------------------------------------------------------

def fetch_dyson_pages(manifest: dict) -> None:
    """The 45 page images (2560 x 3300 px) of the 1920 paper."""
    for leaf in range(C.DYSON_N_LEAVES):
        url = C.DYSON_PAGE_URL.format(leaf=leaf)
        download(url, C.PAGES_DIR / f"leaf{leaf:03d}.jpg", manifest,
                 f"Dyson, Eddington & Davidson 1920, printed page {291 + leaf}")


def fetch_ads_scans(manifest: dict) -> None:
    for key, (bibcode, folder) in C.ADS_SCANS.items():
        download(C.ADS_PDF_URL.format(bibcode=bibcode), folder / f"{bibcode}.pdf", manifest,
                 f"{key} (ADS scan, image-only PDF)")


# ---------------------------------------------------------------------------
# 3. Gaia DR3 around the 1919 field
# ---------------------------------------------------------------------------

GAIA_COLUMNS = ["Source", "RA_ICRS", "DE_ICRS", "e_RA_ICRS", "e_DE_ICRS", "Plx", "e_Plx",
                "pmRA", "e_pmRA", "pmDE", "e_pmDE", "RV", "e_RV", "Gmag", "BPmag", "RPmag", "RUWE"]


def fetch_gaia_field(manifest: dict) -> Path:
    """Gaia DR3 stars with G < 10 within 2.5 deg of the 1919 field (epoch 2016.0)."""
    dest = C.GAIA_DIR / "gaia_dr3_hyades_field.csv"
    url = f"vizier:{C.VIZIER_GAIA} cone (66.0, +21.9) r=2.5 deg, Gmag<10"
    if not dest.exists():
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        from astroquery.vizier import Vizier

        print(f"  querying {url}")
        v = Vizier(columns=GAIA_COLUMNS, row_limit=-1)
        tab = v.query_region(SkyCoord(66.0 * u.deg, 21.9 * u.deg), radius=2.5 * u.deg,
                             catalog=C.VIZIER_GAIA, column_filters={"Gmag": "<10"})[0]
        df = tab.to_pandas()
        missing = [g for g, *_ in C.TABLE1_GAIA.values() if g not in set(df["Source"])]
        if missing:
            raise RuntimeError(f"Table I stars missing from the Gaia query: {missing}")
        _save(dest, df.to_csv(index=False).encode())
    _record(manifest, dest, url, "Gaia DR3 (ESA/Gaia/DPAC) via VizieR, epoch J2016.0")
    return dest


# ---------------------------------------------------------------------------
# 4. JPL Horizons
# ---------------------------------------------------------------------------

#: Observing sites, as (longitude E deg, latitude deg, altitude km) for
#: Horizons' topocentric 'coord' site. Sobral: the expedition's station by the
#: Sobral racecourse (Dyson et al. p.297); Principe: Roca Sundy, "latitude
#: 1deg 40' N., longitude 29m. 32s. E." at 500 feet (pp.313, 315); Oxford: the
#: University Observatory. Site errors move the Sun by < 0.01" per km, but the
#: Moon by ~0.6" per km (~1 s per km in the contact times).
SITES = {
    "sobral": (-40.35, -3.69, 0.07),
    "principe": (7.3833, 1.6667, 0.152),
    "oxford": (-1.25, 51.76, 0.07),
}


def _horizons(params: dict) -> str:
    q = {"format": "text", **{k: f"'{v}'" if isinstance(v, str) else v for k, v in params.items()}}
    return _fetch(C.HORIZONS_API + "?" + urllib.parse.urlencode(q), timeout=600).decode()


def _soe_block(text: str) -> list[str]:
    a, b = text.index("$$SOE"), text.index("$$EOE")
    return [ln for ln in text[a + 5:b].strip().splitlines() if ln.strip()]


def fetch_sun_1919(manifest: dict) -> Path:
    """Astrometric ICRF RA/Dec and distance of the Sun at 1-minute steps
    around totality, seen from Sobral and Principe (topocentric)."""
    dest = C.HORIZONS_DIR / "sun_1919_eclipse.csv"
    url = C.HORIZONS_API + " (Sun, OBSERVER, 1919-05-29 11:30-14:45 UT, sites Sobral/Principe)"
    if not dest.exists():
        rows = []
        for site in ("sobral", "principe"):
            lon, lat, alt = SITES[site]
            text = _horizons({
                "COMMAND": "10", "EPHEM_TYPE": "OBSERVER", "CENTER": "coord@399",
                "COORD_TYPE": "GEODETIC", "SITE_COORD": f"{lon},{lat},{alt}",
                "START_TIME": "1919-05-29 11:30", "STOP_TIME": "1919-05-29 14:45",
                "STEP_SIZE": "1m", "QUANTITIES": "1,4,13,20", "ANG_FORMAT": "DEG",
                "CSV_FORMAT": "YES", "EXTRA_PREC": "YES", "TIME_TYPE": "UT",
            })
            for ln in _soe_block(text):
                f = [s.strip() for s in ln.split(",")]
                # date, solar-presence, lunar-presence, RA, DEC, azimuth, elevation,
                # angular diameter, delta, deldot
                rows.append({"site": site, "utc": f[0], "ra_deg": float(f[3]), "dec_deg": float(f[4]),
                             "az_deg": float(f[5]), "el_deg": float(f[6]),
                             "ang_diam_arcsec": float(f[7]), "delta_au": float(f[8])})
        _save(dest, pd.DataFrame(rows).to_csv(index=False).encode())
    _record(manifest, dest, url, "JPL Horizons, Sun (10), astrometric RA/Dec ICRF, topocentric")
    return dest


#: Totality windows (UT, with a few minutes' margin) for the Moon-Sun track.
TOTALITY_WINDOWS = {"sobral": ("1919-05-29 11:55", "1919-05-29 12:06"),
                    "principe": ("1919-05-29 14:10", "1919-05-29 14:21")}


def fetch_moon_sun_1919(manifest: dict) -> Path:
    """The Moon and the Sun seen from each station through totality, every 10 s
    (66 equal steps of Horizons' 11-minute window): contact times, and the
    Moon's offset that dates Plate 1 of the 1920 paper (plate1.py)."""
    dest = C.HORIZONS_DIR / "moon_sun_sobral_1919.csv"
    url = C.HORIZONS_API + " (Moon 301 and Sun 10, OBSERVER, totality at Sobral and Principe, 10 s steps)"
    if not dest.exists():
        frames = []
        for site, (start, stop) in TOTALITY_WINDOWS.items():
            lon, lat, alt = SITES[site]
            cols = {}
            for body, name in (("301", "moon"), ("10", "sun")):
                text = _horizons({
                    "COMMAND": body, "EPHEM_TYPE": "OBSERVER", "CENTER": "coord@399",
                    "COORD_TYPE": "GEODETIC", "SITE_COORD": f"{lon},{lat},{alt}",
                    "START_TIME": start, "STOP_TIME": stop, "STEP_SIZE": "66",
                    "QUANTITIES": "1,13", "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES",
                    "EXTRA_PREC": "YES", "TIME_TYPE": "UT",
                })
                f = [[s.strip() for s in ln.split(",")] for ln in _soe_block(text)]
                # date, solar-presence, lunar-presence, RA, DEC, angular diameter
                cols["utc"] = [r[0] for r in f]
                cols[f"{name}_ra"] = [float(r[3]) for r in f]
                cols[f"{name}_dec"] = [float(r[4]) for r in f]
                cols[f"{name}_diam_arcsec"] = [float(r[5]) for r in f]
            frames.append(pd.DataFrame({"site": site, **cols}))
        _save(dest, pd.concat(frames).to_csv(index=False, float_format="%.9f").encode())
    _record(manifest, dest, url, "JPL Horizons, Moon (301) and Sun (10), astrometric RA/Dec ICRF, topocentric")
    return dest


def fetch_cassini_2002(manifest: dict) -> Path:
    """Heliocentric ICRF state vectors of Cassini (-82) and the Earth (399),
    every 10 minutes from 2002-05-20 to 2002-07-25 (the SCE1 conjunction)."""
    dest = C.HORIZONS_DIR / "cassini_earth_2002.csv"
    url = C.HORIZONS_API + " (VECTORS, -82 and 399 @ Sun, 2002-05-20..07-25, 10 min)"
    if not dest.exists():
        frames = []
        for body, name in ((C.CASSINI_ID, "cassini"), ("399", "earth")):
            text = _horizons({
                "COMMAND": body, "EPHEM_TYPE": "VECTORS", "CENTER": "500@10",
                "START_TIME": "2002-05-20", "STOP_TIME": "2002-07-25", "STEP_SIZE": "10m",
                "VEC_TABLE": "2", "REF_PLANE": "FRAME", "OUT_UNITS": "KM-S",
                "CSV_FORMAT": "YES", "VEC_LABELS": "NO", "TIME_TYPE": "TDB",
            })
            arr = []
            for ln in _soe_block(text):
                f = [s.strip() for s in ln.split(",")]
                arr.append([float(f[0])] + [float(v) for v in f[2:8]])
            arr = np.array(arr)
            frames.append(pd.DataFrame({
                "jd_tdb": arr[:, 0], **{f"{name}_{k}": arr[:, i + 1]
                                        for i, k in enumerate(("x", "y", "z", "vx", "vy", "vz"))}}))
        df = frames[0].merge(frames[1], on="jd_tdb")
        _save(dest, df.to_csv(index=False, float_format="%.12g").encode())
    _record(manifest, dest, url, "JPL Horizons, heliocentric ICRF vectors (km, km/s), TDB")
    return dest


# ---------------------------------------------------------------------------
# 5. A real image of the 1919 field
# ---------------------------------------------------------------------------

def fetch_sky_image(manifest: dict) -> Path:
    """DSS2 colour, 4.2 x 3.6 deg (TAN, north up), centred near the 1919 Sun."""
    dest = C.SKY_DIR / "hyades_dss2_color.jpg"
    params = {"hips": "CDS/P/DSS2/color", "width": 2100, "height": 1800, "projection": "TAN",
              "fov": 4.2, "ra": 66.05, "dec": 21.85, "format": "jpg"}
    url = C.HIPS2FITS_URL + "?" + urllib.parse.urlencode(params)
    download(url, dest, manifest, "DSS2 colour (CDS hips2fits); STScI/ESO/AAO/Caltech plates", timeout=600)
    return dest


# ---------------------------------------------------------------------------
# 6. (optional) Cassini's raw Doppler data from the 2002 conjunction
# ---------------------------------------------------------------------------

PDS_CORS = "https://pds-atmospheres.nmsu.edu/PDS/data/cors_{vol:04d}/"
CASSINI_DIR = C.RAW_DIR / "cassini_sce1"


def _listing(url: str) -> list[str]:
    import re

    html = _fetch(url, timeout=120).decode("latin-1")
    return sorted(set(re.findall(r'href="([^"?/][^"?]*)"', html)))


def fetch_cassini_sce1(manifest: dict) -> None:
    """Orbit Data Files (DSN Doppler/range observables, ~3.5 MB a day) for
    every day of the first Cassini Solar Conjunction Experiment, with their
    PDS labels, and the troposphere/ionosphere/Earth-orientation calibration
    files (PDS volumes CORS_0021-0028, PDS Atmospheres Node). The gamma fit
    starts from the navigation team's SCE1 short-arc SPK (fetch_cassini_sce1_extras);
    Horizons' reconstruction (fetch_cassini_ssb_vectors) serves the inventory-wide
    model and the horizons_* variants."""
    for vol in range(21, 29):
        base = PDS_CORS.format(vol=vol)
        for day in [d.rstrip("/") for d in _listing(base) if d.startswith("sce1_") and d[5:8].isdigit()]:
            try:
                files = _listing(base + day + "/odf/")
            except RuntimeError:      # e.g. CORS_0028 has no odf/ folder for day 186
                print(f"  (no ODF folder for {day} on CORS_{vol:04d})")
                continue
            for f in files:
                if f.endswith((".odf", ".lbl")):
                    download(base + day + "/odf/" + f, CASSINI_DIR / "odf" / f, manifest,
                             f"Cassini SCE1 ODF ({day}), PDS CORS_{vol:04d}")
        for sub in ("tro", "ion", "eop"):
            try:
                files = _listing(base + "sce1_ancillary/" + sub + "/")
            except RuntimeError:
                continue
            for f in files:
                if f.endswith((f".{sub}", ".lbl")) and not (CASSINI_DIR / sub / f).exists():
                    download(base + "sce1_ancillary/" + sub + "/" + f, CASSINI_DIR / sub / f, manifest,
                             f"Cassini SCE1 {sub} calibration, PDS CORS_{vol:04d}")


def fetch_cassini_ssb_vectors(manifest: dict) -> Path:
    """Barycentric (SSB) ICRF state vectors of the Earth's centre and of Cassini
    every 2 minutes, 2002-06-04 to 2002-07-09 (TDB), for the light-time
    solution of the Doppler model."""
    dest = C.HORIZONS_DIR / "cassini_earth_ssb_2002.csv"
    url = C.HORIZONS_API + " (VECTORS, -82 and 399 @ SSB, 2002-06-04..07-09, 2 min)"
    if not dest.exists():
        frames = []
        for body, name in ((C.CASSINI_ID, "cassini"), ("399", "earth")):
            text = _horizons({
                "COMMAND": body, "EPHEM_TYPE": "VECTORS", "CENTER": "500@0",
                "START_TIME": "2002-06-04", "STOP_TIME": "2002-07-09", "STEP_SIZE": "2m",
                "VEC_TABLE": "2", "REF_PLANE": "FRAME", "OUT_UNITS": "KM-S",
                "CSV_FORMAT": "YES", "VEC_LABELS": "NO", "TIME_TYPE": "TDB",
            })
            arr = np.array([[float(f[0])] + [float(v) for v in f[2:8]]
                            for f in ([s.strip() for s in ln.split(",")] for ln in _soe_block(text))])
            frames.append(pd.DataFrame({"jd_tdb": arr[:, 0], **{f"{name}_{k}": arr[:, i + 1]
                                        for i, k in enumerate(("x", "y", "z", "vx", "vy", "vz"))}}))
        df = frames[0].merge(frames[1], on="jd_tdb")
        _save(dest, df.to_csv(index=False, float_format="%.15g").encode())
    _record(manifest, dest, url, "JPL Horizons, SSB-centred ICRF vectors (km, km/s), TDB")
    return dest


def fetch_cassini_sun_vectors(manifest: dict) -> Path:
    """The Sun's barycentric (SSB) ICRF state vectors on the same 2-minute TDB
    grid. The Shapiro delay is computed about the Sun's centre, which the
    planets (Jupiter, Saturn) pull ~1.2 solar radii away from the SSB (1.15-1.17
    in June 2002). Nearly all of that lay along the line of sight, but the ~0.02
    across it still moves the impact parameter by 0.5% on the pass nearest the Sun."""
    dest = C.HORIZONS_DIR / "sun_ssb_2002.csv"
    url = C.HORIZONS_API + " (VECTORS, 10 @ SSB, 2002-06-04..07-09, 2 min)"
    if not dest.exists():
        text = _horizons({
            "COMMAND": "10", "EPHEM_TYPE": "VECTORS", "CENTER": "500@0",
            "START_TIME": "2002-06-04", "STOP_TIME": "2002-07-09", "STEP_SIZE": "2m",
            "VEC_TABLE": "2", "REF_PLANE": "FRAME", "OUT_UNITS": "KM-S",
            "CSV_FORMAT": "YES", "VEC_LABELS": "NO", "TIME_TYPE": "TDB",
        })
        arr = np.array([[float(f[0])] + [float(v) for v in f[2:8]]
                        for f in ([s.strip() for s in ln.split(",")] for ln in _soe_block(text))])
        df = pd.DataFrame({"jd_tdb": arr[:, 0], **{f"sun_{k}": arr[:, i + 1]
                           for i, k in enumerate(("x", "y", "z", "vx", "vy", "vz"))}})
        _save(dest, df.to_csv(index=False, float_format="%.15g").encode())
    _record(manifest, dest, url, "JPL Horizons, the Sun (10) about the SSB, ICRF vectors (km, km/s), TDB")
    return dest


#: Extra SCE1 products beyond the ODFs: the archive's documents (the AWVR
#: path-delay SIS DORS-002, the TRK-2-21/2-23 formats), the navigation team's
#: short-arc SCE1 trajectory (SPK "SM727D_SHORTARC_SCE1.BSP", on CORS_0021
#: only), the Archival Tracking Data Files of the X-band uplinks (TRK-2-25,
#: ~45 MB a day), which also hold DSS 25's Ka-band channel (labelled with an X
#: uplink but carrying the Ka/Ka link, including the passes of June 14 and 17
#: that the ODFs leave out; cassini_tdf.py), and the Advanced Water Vapour
#: Radiometer path delays at DSS 25 (pd1/pd2, the tropospheric calibration
#: Anderson et al. 2004 used).
SCE1_DOCUMENTS = ("casrssis.txt", "docinfo.txt", "dors_002_020501.txt", "rol_021101.txt",
                  "rsrlog_021101.txt", "spk_mm_sis_000525.txt", "trk_2_21_950831.txt",
                  "trk_2_23_000531.txt")


def fetch_cassini_sce1_extras(manifest: dict) -> None:
    """Documents, the SCE1 short-arc SPK, the X-uplink TDFs and the AWVR
    path-delay files (~1.4 GB, almost all of it the TDFs)."""
    base = PDS_CORS.format(vol=21)
    for f in SCE1_DOCUMENTS:
        download(base + "document/" + f, CASSINI_DIR / "document" / f, manifest, "Cassini SCE1 archive document")
    for f in ("c32easc2002_160_2002_186.spk", "c32easc2002_160_2002_186.lbl"):
        download(base + "sce1_ancillary/spk/" + f, CASSINI_DIR / "spk" / f, manifest,
                 "Cassini navigation short-arc SCE1 trajectory (SM727D_SHORTARC_SCE1.BSP)")
    for vol in range(21, 29):
        base = PDS_CORS.format(vol=vol)
        for day in [d.rstrip("/") for d in _listing(base) if d.startswith("sce1_") and d[5:8].isdigit()]:
            for sub, keep in (("tdf", "xmmmv0"), ("pd1", ""), ("pd2", "")):
                try:
                    files = _listing(base + day + f"/{sub}/")
                except RuntimeError:
                    continue
                for f in files:
                    if keep in f and f.endswith((f".{sub}", ".lbl")):
                        download(base + day + f"/{sub}/" + f, CASSINI_DIR / sub / f, manifest,
                                 f"Cassini SCE1 {sub.upper()} ({day}), PDS CORS_{vol:04d}", timeout=3600)


# ---------------------------------------------------------------------------

def fetch_all(cassini: bool = False) -> dict:
    manifest = _load_manifest()
    print("1. Dyson, Eddington & Davidson 1920 page images")
    fetch_dyson_pages(manifest)
    print("2. ADS scans (Wallal 1922; Harvey 1979)")
    fetch_ads_scans(manifest)
    print("3. Gaia DR3 around the 1919 field")
    fetch_gaia_field(manifest)
    print("4. JPL Horizons (Sun 1919; Cassini 2002)")
    fetch_sun_1919(manifest)
    fetch_moon_sun_1919(manifest)
    fetch_cassini_2002(manifest)
    print("5. Sky image of the 1919 field")
    fetch_sky_image(manifest)
    if cassini:
        print("6. Cassini SCE1 raw Doppler data (~1.5 GB)")
        fetch_cassini_sce1(manifest)
        fetch_cassini_ssb_vectors(manifest)
        fetch_cassini_sun_vectors(manifest)
        fetch_cassini_sce1_extras(manifest)
    C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    total = sum(f["bytes"] for f in manifest["files"].values())
    print(f"   {len(manifest['files'])} files, {total / 1e6:.1f} MB")
    return manifest


if __name__ == "__main__":
    import sys

    fetch_all(cassini="--cassini" in sys.argv)
