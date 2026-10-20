"""Download and cache every input this project uses (~0.75 GB in all).

    python fetch_data.py            # everything that is missing
    python fetch_data.py aavso      # one source only (names below)

Sources, by step name:

1. ``ogle``    OGLE-III/IV photometry of V1309 Sco 2001-2010, published with
               Tylenda et al. (2011) at astrouw.edu.pl/ogle/ogle3/V1309_SCO/.
2. ``kepler``  Kepler long-cadence light curves of KIC 9832227, Q0-Q17 (MAST).
3. ``tess``    TESS light curves of KIC 9832227 (SPOC 2-min; QLP and TESS-SPOC
               full-frame-image products), 2019-2024, via lightkurve/MAST.
   ``archival`` Older photometry of KIC 9832227: NSVS 1999-2000 (VizieR
               II/287), SuperWASP 2004-2010 (DR1 and VeSPA), ASAS-SN 2014-2018.
4. ``aavso``   Every AAVSO International Database observation of T CrB,
               1866 to today, through the public VSX API (no account needed).
5. ``ztf``     The ZTF catalogue of periodic variables (Chen, Wang & Deng
               2020, VizieR J/ApJS/249/18), and Chen et al. (2018)'s W UMa
               calibrators (J/ApJ/859/140).
6. ``xmatch``  CDS xMatch of the ZTF eclipsing binaries with Gaia DR3 (I/355).
7. ``allwise`` AllWISE W1/W2 for the ZTF binaries with good parallaxes (CDS xMatch);
   ``chen``    Chen et al.'s 183 calibrators against Gaia DR3 and AllWISE (CDS xMatch);
   ``veb``     Gaia DR3 eclipsing-binary candidates with P <= 2 d (VizieR I/358/veb).
8. ``lit``     arXiv sources of the papers whose tables we parse (literature.py).

Files land in ``data/raw/`` and are never re-downloaded once present, so after
one run the analysis works offline. ``data/manifest.json`` records each file's
source, size, SHA-256 and download date.
"""

from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

import numpy as np

import config as C


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _fetch(url: str, *, data: bytes | None = None, headers: dict | None = None,
           timeout: float = 180.0, retries: int = 4) -> bytes:
    """GET (or POST, with ``data``) ``url`` and return the body, retrying transient failures."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": C.USER_AGENT, **(headers or {})})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (400, 403, 404, 410):
                body = exc.read()[:500]
                raise RuntimeError(f"HTTP {exc.code} for {url}: {body!r}") from exc
            time.sleep(5.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def _check_not_html(payload: bytes, what: str) -> None:
    head = payload[:400].lstrip().lower()
    if head.startswith(b"<!doctype html") or head.startswith(b"<html"):
        raise RuntimeError(f"{what}: got an HTML page instead of data: {payload[:200]!r}")


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
    key = str(dest.relative_to(C.DATA_DIR))
    if key in manifest["files"]:
        return
    blob = dest.read_bytes()
    manifest["files"][key] = {
        "url": url, "note": note, "bytes": len(blob),
        "sha256": hashlib.sha256(blob).hexdigest(),
        "downloaded": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
    }


def _write_manifest(manifest: dict) -> None:
    """Merge our entries into the manifest on disk (several fetches may run at once)."""
    import fcntl
    with open(C.MANIFEST_PATH.with_suffix(".lock"), "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        on_disk = _load_manifest()
        on_disk["files"].update(manifest["files"])
        manifest["files"] = on_disk["files"]
        tmp = C.MANIFEST_PATH.with_suffix(".json.part")
        tmp.write_text(json.dumps(on_disk, indent=1, sort_keys=True))
        tmp.replace(C.MANIFEST_PATH)


def download(manifest: dict, url: str, dest: Path, note: str, *, check_html: bool = True,
             timeout: float = 180.0) -> Path:
    """Fetch ``url`` into ``dest`` unless it is already there, and record it."""
    if not dest.exists():
        payload = _fetch(url, timeout=timeout)
        if check_html:
            _check_not_html(payload, url)
        _save(dest, payload)
    _record(manifest, dest, url, note)
    return dest


# ---------------------------------------------------------------------------
# 1. OGLE photometry of V1309 Sco
# ---------------------------------------------------------------------------

def fetch_ogle(manifest: dict) -> None:
    print("OGLE: V1309 Sco photometry ...")
    download(manifest, C.OGLE_V1309_URL + "README", C.OGLE_DIR / "README",
             "OGLE V1309 Sco data README (Tylenda et al. 2011 supplementary material)")
    for key, name in C.V1309SCO["ogle_files"].items():
        download(manifest, C.OGLE_V1309_URL + name, C.OGLE_DIR / name,
                 f"OGLE V1309 Sco photometry ({key}): HJD, mag, error")


# ---------------------------------------------------------------------------
# 2. Kepler long-cadence light curves of KIC 9832227
# ---------------------------------------------------------------------------

def fetch_kepler(manifest: dict) -> None:
    kic = f"{C.KIC9832227['kic']:09d}"
    base = f"{C.KEPLER_LC_URL}{kic[:4]}/{kic}/"
    print(f"Kepler: listing {base} ...")
    listing = _fetch(base).decode("utf-8", "replace")
    import re
    names = sorted(set(re.findall(r'href="(kplr\d+-\d+_llc\.fits)"', listing)))
    if not names:
        raise RuntimeError(f"no long-cadence files listed at {base}")
    for n in names:
        download(manifest, base + n, C.KEPLER_DIR / n, "Kepler long-cadence light curve (PDCSAP/SAP)",
                 check_html=False)
    print(f"  {len(names)} quarters cached")


# ---------------------------------------------------------------------------
# 3. TESS light curves of KIC 9832227
# ---------------------------------------------------------------------------

#: Products kept, in order of preference for each sector: SPOC 2-min target
#: pixels photometry first, then the full-frame-image pipelines.
TESS_AUTHORS = ("SPOC", "TESS-SPOC", "QLP")


def fetch_tess(manifest: dict) -> None:
    import warnings
    warnings.filterwarnings("ignore")
    import lightkurve as lk
    print("TESS: searching MAST for TIC", C.KIC9832227["tic"], "...")
    res = lk.search_lightcurve(f"TIC {C.KIC9832227['tic']}", mission="TESS")
    tab = res.table
    for i in range(len(tab)):
        author = str(tab["author"][i])
        if author not in TESS_AUTHORS:
            continue
        uri = str(tab["dataURI"][i])
        name = uri.split("/")[-1]
        dest = C.TESS_DIR / name
        url = "https://mast.stsci.edu/api/v0.1/Download/file?uri=" + urllib.parse.quote(uri, safe="")
        download(manifest, url, dest, f"TESS {author} light curve, {tab['mission'][i]}, "
                 f"exptime {float(tab['exptime'][i]):.0f} s", check_html=False, timeout=300)
    print(f"  {len(list(C.TESS_DIR.glob('*.fits')))} TESS files cached")


# ---------------------------------------------------------------------------
# 3b. Older photometry of KIC 9832227: NSVS (1999-2000), SuperWASP (2004-2010),
#     ASAS-SN (2014-2018)
# ---------------------------------------------------------------------------

#: NSVS (Wozniak et al. 2004) on VizieR II/287: the star appears under two IDs
#: because it straddles two fields. The per-field files are sorted by object
#: then frame, so each is streamed only until the object has been passed.
NSVS_IDS = {5597755: "042d", 5620022: "043a"}
NSVS_FTP = "https://cdsarc.cds.unistra.fr/ftp/II/287/"
WASP_DR1_URL = ("https://exoplanetarchive.ipac.caltech.edu/data/ETSS/SuperWASP/FITS/DR1/tile288132/"
                "1SWASP%20J192916.01%2B463719.9.fits")
WASP_VESPA_URL = ("https://www.superwasp.org/media/sources/1SWASPJ192916.01%2B463719.9/"
                  "v0.9_1SWASPJ192916.01463719.9.fits")
ASASSN_VAR_URL = "https://asas-sn.osu.edu/variables/cadbbf9d-c501-5888-9792-bb1ea9e76030.json"


def _stream_nsvs(object_id: int, field: str) -> list[str]:
    import zlib
    url = f"{NSVS_FTP}id/skydotID_{field}.dat.gz"
    req = urllib.request.Request(url, headers={"User-Agent": C.USER_AGENT})
    dec = zlib.decompressobj(16 + zlib.MAX_WBITS)
    rows, buf = [], b""
    with urllib.request.urlopen(req, timeout=600) as resp:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            buf += dec.decompress(chunk)
            lines = buf.split(b"\n")
            buf = lines.pop()
            done = False
            for ln in lines:
                parts = ln.split()
                if not parts:
                    continue
                oid = int(parts[0])
                if oid == object_id:
                    rows.append(ln.decode())
                elif oid > object_id:
                    done = True
                    break
            if done:
                break
    if not rows:
        raise RuntimeError(f"NSVS object {object_id} not found in {url}")
    return rows


def fetch_kic_archival(manifest: dict) -> None:
    print("NSVS (VizieR II/287), SuperWASP and ASAS-SN photometry of KIC 9832227 ...")
    d = C.KEPLER_DIR.parent / "archival"
    d.mkdir(exist_ok=True)
    download(manifest, NSVS_FTP + "ReadMe", d / "II_287_ReadMe.txt", "VizieR II/287 (NSVS) ReadMe")
    download(manifest, NSVS_FTP + "frames.dat.gz", d / "II_287_frames.dat.gz",
             "NSVS frame table: the 'JD' column holds MJD (JD - 2400000.5) despite its label", check_html=False)
    for oid, field in NSVS_IDS.items():
        dest = d / f"nsvs_{oid}.dat"
        if not dest.exists():
            rows = _stream_nsvs(oid, field)
            _save(dest, ("\n".join(rows) + "\n").encode())
            print(f"  NSVS {oid}: {len(rows)} measurements")
        _record(manifest, dest, f"{NSVS_FTP}id/skydotID_{field}.dat.gz",
                f"NSVS photometry of object {oid} (rows of field file {field})")
    download(manifest, WASP_DR1_URL, d / "1SWASP_J192916.01+463719.9_DR1.fits",
             "SuperWASP DR1 light curve (Butters et al. 2010), NASA Exoplanet Archive", check_html=False)
    download(manifest, WASP_VESPA_URL, d / "1SWASP_J192916.01+463719.9_vespa_v0.9.fits",
             "SuperWASP light curve 2004-2010 from VeSPA (McMaster et al. 2021)", check_html=False)
    download(manifest, ASASSN_VAR_URL, d / "asassn_variables_J192915.86+463719.9.json",
             "ASAS-SN Variable Stars Database light curve (V band, 2014-2018)")
    # AAVSO observers' photometry 2015-2026 (VSX name ASAS J192916+4637.3), one request.
    dest = d / "aavso_kic9832227.json.gz"
    if not dest.exists():
        q = urllib.parse.urlencode({"view": "api.object", "ident": "KIC 9832227", "data": "500000",
                                    "format": "json"})
        payload = _fetch(f"{C.AAVSO_VSX_API}?{q}", timeout=600)
        _check_not_html(payload, "AAVSO VSX API")
        obs = json.loads(payload)["VSXObject"]["Data"]["Observation"]
        _save(dest, gzip.compress(json.dumps(obs).encode()))
        print(f"  AAVSO: {len(obs)} observations of KIC 9832227")
    _record(manifest, dest, C.AAVSO_VSX_API, "AAVSO International Database observations of KIC 9832227")


# ---------------------------------------------------------------------------
# 4. AAVSO observations of T CrB (VSX API)
# ---------------------------------------------------------------------------

class _WindowTooBig(RuntimeError):
    pass


def _aavso_window(jd0: float, jd1: float, attempts: int = 5) -> list[dict]:
    q = urllib.parse.urlencode({"view": "api.object", "ident": C.TCRB["name"], "data": str(10 * C.AAVSO_MAX_PER_REQUEST),
                                "fromjd": f"{jd0:.4f}", "tojd": f"{jd1:.4f}", "format": "json"})
    for attempt in range(attempts):
        payload = _fetch(f"{C.AAVSO_VSX_API}?{q}", timeout=600)
        _check_not_html(payload, "AAVSO VSX API")
        if b"Maximum execution time" in payload[:400]:
            # the server's PHP time limit (30 s): the window holds too many observations
            raise _WindowTooBig(payload[:200])
        try:
            obj = json.loads(payload)["VSXObject"]
            break
        except (json.JSONDecodeError, KeyError):
            # an empty or truncated body now and then: wait and ask again
            print(f"    JD {jd0:.0f}-{jd1:.0f}: unreadable response ({len(payload)} bytes), retrying")
            time.sleep(20.0 * (attempt + 1))
    else:
        raise RuntimeError(f"AAVSO VSX API: no readable response for JD {jd0}-{jd1}")
    data = obj.get("Data") or {}
    obs = data.get("Observation") or []
    if isinstance(obs, dict):  # a single observation comes back as an object
        obs = [obs]
    return obs


def fetch_aavso(manifest: dict) -> None:
    """Every observation of T CrB, saved window by window as gzipped JSON.

    The VSX API returns all observations between ``fromjd`` and ``tojd``. A
    window that returns more than ``AAVSO_MAX_PER_REQUEST`` rows is split in
    two and re-fetched, so no single response gets huge. Files are named by
    their JD range; the analysis concatenates them (overlaps are removed by
    observation Id).
    """
    print("AAVSO: T CrB observations from the VSX API ...")
    today_jd = 2440587.5 + time.time() / C.DAY_S
    todo = []
    for start, stop, step in C.AAVSO_WINDOWS:
        jd = start
        while jd < stop:
            todo.append((jd, min(jd + step, stop)))
            jd += step
    n_total = 0
    while todo:
        jd0, jd1 = todo.pop(0)
        dest = C.AAVSO_DIR / f"tcrb_{jd0:.0f}_{jd1:.0f}.json.gz"
        # The last window is re-fetched on every run until it lies in the past,
        # so a re-run picks up the newest observations.
        if dest.exists() and jd1 < today_jd - 2:
            _record(manifest, dest, C.AAVSO_VSX_API, f"AAVSO observations of T CrB, JD {jd0:.0f}-{jd1:.0f}")
            continue
        try:
            obs = _aavso_window(jd0, jd1)
        except _WindowTooBig:
            mid = 0.5 * (jd0 + jd1)
            todo[:0] = [(jd0, mid), (mid, jd1)]
            print(f"  JD {jd0:.0f}-{jd1:.0f}: server time limit, splitting")
            time.sleep(C.AAVSO_PAUSE_S)
            continue
        if len(obs) > C.AAVSO_MAX_PER_REQUEST and jd1 - jd0 > 2:
            mid = 0.5 * (jd0 + jd1)
            todo[:0] = [(jd0, mid), (mid, jd1)]
            print(f"  JD {jd0:.0f}-{jd1:.0f}: {len(obs)} rows, splitting")
            time.sleep(C.AAVSO_PAUSE_S)
            continue
        if dest.exists():
            dest.unlink()
            manifest["files"].pop(str(dest.relative_to(C.DATA_DIR)), None)
        _save(dest, gzip.compress(json.dumps(obs).encode()))
        _record(manifest, dest, C.AAVSO_VSX_API, f"AAVSO observations of T CrB, JD {jd0:.0f}-{jd1:.0f}")
        n_total += len(obs)
        print(f"  JD {jd0:.0f}-{jd1:.0f}: {len(obs)} observations")
        _write_manifest(manifest)
        time.sleep(C.AAVSO_PAUSE_S)
    print(f"  {n_total} new observations fetched")


# ---------------------------------------------------------------------------
# 5. ZTF periodic variables and Chen et al. (2018) calibrators (VizieR)
# ---------------------------------------------------------------------------

ZTF_COLUMNS = ("ID", "RAJ2000", "DEJ2000", "Per", "R21", "phi21", "T0", "gmag", "rmag",
               "Ng", "Nr", "gAmp", "rAmp", "log(FAP-g)", "log(FAP-r)", "R2-g", "R2-r",
               "Type", "Dmin-g", "Dmin-r")


def _vizier_tsv(source: str, columns, constraints: dict | None = None) -> bytes:
    params = {"-source": source, "-out": ",".join(columns), "-out.max": "unlimited",
              "-oc.form": "dec", "-out.meta": "hu"}
    params.update(constraints or {})
    payload = _fetch(f"{C.VIZIER_ASU}?{urllib.parse.urlencode(params)}", timeout=900)
    _check_not_html(payload, f"VizieR {source}")
    return payload


def fetch_ztf(manifest: dict) -> None:
    print("VizieR: ZTF periodic variables (Chen et al. 2020) ...")
    step = 15  # degrees of RA per request (~33,000 rows each)
    for ra0 in range(0, 360, step):
        dest = C.ZTF_DIR / f"ztf_periodic_ra{ra0:03d}.tsv"
        if not dest.exists():
            body = _vizier_tsv("J/ApJS/249/18/table2", ZTF_COLUMNS, {"RAJ2000": f"{ra0}..{ra0 + step}"})
            # VizieR range constraints are inclusive at both ends: drop the upper edge here and keep it
            # in the next slice, so no star appears twice.
            keep = []
            for ln in body.decode().splitlines():
                if ln.startswith("ZTFJ") and ra0 + step < 360 and float(ln.split("\t")[1]) >= ra0 + step:
                    continue
                keep.append(ln)
            _save(dest, ("\n".join(keep) + "\n").encode())
            print(f"  RA {ra0}-{ra0 + step}: {sum(1 for k in keep if k.startswith('ZTF'))} rows")
        _record(manifest, dest, C.VIZIER_ASU, f"VizieR J/ApJS/249/18/table2, {ra0} <= RA < {ra0 + step}")
    for table, cols in (("J/ApJ/859/140/table1", ("RAJ2000", "DEJ2000", "Period", "DM", "e_DM", "Amp", "Av",
                                                  "VMagMax", "GMagMax", "B-V", "W1MagMax", "W2MagMax",
                                                  "VMag", "Gmag", "W1Mag", "W2Mag", "KsMag")),
                        ("J/ApJ/859/140/table4", ("RAJ2000", "DEJ2000", "Sel", "Period", "AW1", "DM-W1",
                                                  "DMAvg", "e_DMAvg", "Amp", "Gmag", "W1mag"))):
        dest = C.ZTF_DIR / (table.replace("/", "_") + ".tsv")
        if not dest.exists():
            _save(dest, _vizier_tsv(table, cols))
        _record(manifest, dest, C.VIZIER_ASU, f"VizieR {table} (Chen et al. 2018)")


# ---------------------------------------------------------------------------
# 6. CDS xMatch: ZTF eclipsing binaries -> Gaia DR3
# ---------------------------------------------------------------------------

#: Gaia columns kept from the xMatch result (VizieR names). The rest of the
#: ~140 columns are dropped chunk by chunk to keep the cache small.
XMATCH_KEEP = ("angDist", "zid", "Source", "RAdeg", "DEdeg", "Plx", "e_Plx", "pmRA", "pmDE", "RUWE",
               "Gmag", "e_Gmag", "BPmag", "RPmag", "BP-RP", "E(BP/RP)", "nueff", "pscol", "Solved",
               "VarFlag", "Teff", "AG", "Dist", "NBPblend", "NRPblend", "IPDfmp", "Dup")


def _xmatch(csv_bytes: bytes) -> bytes:
    boundary = uuid.uuid4().hex
    fields = {"request": "xmatch", "distMaxArcsec": str(C.XMATCH_RADIUS_ARCSEC),
              "RESPONSEFORMAT": "csv", "colRA1": "ra", "colDec1": "dec", "cat2": C.VIZIER_GAIA_DR3}
    out = io.BytesIO()
    for k, v in fields.items():
        out.write(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode())
    out.write(f"--{boundary}\r\nContent-Disposition: form-data; name=\"cat1\"; filename=\"ztf.csv\"\r\n"
              f"Content-Type: text/csv\r\n\r\n".encode())
    out.write(csv_bytes)
    out.write(f"\r\n--{boundary}--\r\n".encode())
    return _fetch(C.CDS_XMATCH, data=out.getvalue(), timeout=900,
                  headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})


def fetch_xmatch(manifest: dict) -> None:
    import pandas as pd
    import catalogs
    ztf = catalogs.read_ztf(types=C.ZTF_BINARY_TYPES)
    print(f"CDS xMatch: {len(ztf)} ZTF eclipsing binaries vs Gaia DR3 ...")
    chunk = 10000  # ~80 s per request; 40,000-row requests hit the gateway's timeout (HTTP 504)
    for k in range(0, len(ztf), chunk):
        dest = C.GAIA_DIR / f"xmatch_ztf_gaiadr3_{k // chunk:03d}.csv.gz"
        if not dest.exists():
            sub = ztf.iloc[k:k + chunk]
            csv = "zid,ra,dec\n" + "\n".join(f"{z},{r:.6f},{d:.6f}" for z, r, d in
                                              zip(sub["zid"], sub["ra"], sub["dec"]))
            body = _xmatch((csv + "\n").encode())
            _check_not_html(body, "CDS xMatch")
            tab = pd.read_csv(io.BytesIO(body), low_memory=False)
            missing = [c for c in XMATCH_KEEP if c not in tab.columns]
            if missing:
                raise RuntimeError(f"xMatch result lacks columns {missing}")
            buf = io.BytesIO()
            tab[list(XMATCH_KEEP)].to_csv(buf, index=False)
            _save(dest, gzip.compress(buf.getvalue()))
            print(f"  chunk {k // chunk}: {len(sub)} stars -> {len(tab)} Gaia matches")
            time.sleep(2.0)
        _record(manifest, dest, C.CDS_XMATCH, f"CDS xMatch ZTF EW/EA rows {k}-{k + chunk - 1} with Gaia DR3 "
                f"(I/355), radius {C.XMATCH_RADIUS_ARCSEC}\"; columns trimmed to {len(XMATCH_KEEP)}")


# ---------------------------------------------------------------------------
# 7. Gaia DR3 and AllWISE through CDS (the ESA archive was congested in October 2026:
#    a 1,000-row upload join did not finish in two minutes)
# ---------------------------------------------------------------------------

#: Columns of the Gaia DR3 eclipsing-binary table on VizieR (I/358/veb).
VEB_COLUMNS = ("Source", "Rank", "Freq", "e_Freq", "magModRef", "sigPhaseGauss1", "DepthGauss1",
               "sigPhaseGauss2", "DepthGauss2", "AmpCHP", "ModelType", "rchi2", "DepthE1", "DepthE2",
               "RA_ICRS", "DE_ICRS")


def fetch_gaia_eb_vizier(manifest: dict) -> None:
    """Gaia DR3 eclipsing-binary candidates with orbital periods <= 2 d
    (frequency >= 0.5 per day), from VizieR I/358/veb in 24 RA slices."""
    print("VizieR: Gaia DR3 eclipsing-binary candidates (I/358/veb), P <= 2 d ...")
    step = 15
    total = 0
    for ra0 in range(0, 360, step):
        dest = C.GAIA_DIR / f"veb_ra{ra0:03d}.tsv"
        if not dest.exists():
            body = _vizier_tsv("I/358/veb", VEB_COLUMNS, {"RA_ICRS": f"{ra0}..{ra0 + step}", "Freq": ">=0.5"})
            keep = []
            for ln in body.decode().splitlines():
                parts = ln.split("\t")
                # (short Gaia IDs are right-aligned with leading blanks)
                if ln.strip() and ln.lstrip()[0].isdigit() and len(parts) > 14 and ra0 + step < 360:
                    try:
                        if float(parts[14]) >= ra0 + step:
                            continue
                    except ValueError:
                        pass
                keep.append(ln)
            _save(dest, ("\n".join(keep) + "\n").encode())
        n = sum(1 for ln in dest.read_text().splitlines() if ln.strip() and ln.lstrip()[0].isdigit())
        total += n
        _record(manifest, dest, C.VIZIER_ASU, f"VizieR I/358/veb (Gaia DR3 EB candidates), Freq >= 0.5/d, "
                f"{ra0} <= RA < {ra0 + step}")
    print(f"  {total} candidates")


def fetch_allwise_xmatch(manifest: dict) -> None:
    """AllWISE W1/W2 for the ZTF binaries with good Gaia parallaxes
    (parallax/error > 5), by CDS xMatch on Gaia positions (2" radius)."""
    import pandas as pd
    import catalogs
    xm = catalogs.read_xmatch()
    xm = xm.sort_values(["zid", "sep_arcsec"]).groupby("zid").head(1)
    good = xm[(xm["parallax"] / xm["parallax_error"] > 5)]
    print(f"CDS xMatch: {len(good)} ZTF binaries with parallax/error > 5 vs AllWISE ...")
    chunk = 10000
    for k in range(0, len(good), chunk):
        dest = C.GAIA_DIR / f"xmatch_allwise_{k // chunk:02d}.csv.gz"
        if not dest.exists():
            sub = good.iloc[k:k + chunk]
            csv = "zid,ra,dec\n" + "\n".join(f"{z},{r:.7f},{d:.7f}" for z, r, d in
                                               zip(sub["zid"], sub["ra_gaia"], sub["dec_gaia"]))
            body = _xmatch_cat((csv + "\n").encode(), "vizier:II/328/allwise", 2.0)
            _check_not_html(body, "CDS xMatch AllWISE")
            tab = pd.read_csv(io.BytesIO(body), low_memory=False)
            cols = [c for c in ("angDist", "zid", "AllWISE", "W1mag", "e_W1mag", "W2mag", "e_W2mag", "ccf", "ex",
                                "qph", "nW1", "mW1") if c in tab.columns]
            buf = io.BytesIO()
            tab[cols].to_csv(buf, index=False)
            _save(dest, gzip.compress(buf.getvalue()))
            print(f"  chunk {k // chunk}: {len(sub)} stars -> {len(tab)} AllWISE matches")
            time.sleep(2.0)
        _record(manifest, dest, C.CDS_XMATCH, f"CDS xMatch of ZTF binaries (Gaia positions) with AllWISE II/328, 2\"")


def fetch_chen_xmatch(manifest: dict) -> None:
    """Chen et al. (2018)'s 183 TGAS calibrators against Gaia DR3 and AllWISE."""
    import catalogs
    t1, _ = catalogs.read_chen2018()
    csv = "cid,ra,dec\n" + "\n".join(f"{i},{r:.6f},{d:.6f}" for i, (r, d) in
                                     enumerate(zip(t1["RAJ2000"], t1["DEJ2000"])))
    for cat, rad, name in (("vizier:I/355/gaiadr3", 3.0, "gaiadr3"), ("vizier:II/328/allwise", 3.0, "allwise")):
        dest = C.ZTF_DIR / f"chen2018_table1_x_{name}.csv"
        if not dest.exists():
            body = _xmatch_cat((csv + "\n").replace("cid,", "zid,").encode(), cat, rad)
            _check_not_html(body, "CDS xMatch")
            _save(dest, body)
        _record(manifest, dest, C.CDS_XMATCH, f"CDS xMatch of Chen et al. (2018) Table 1 with {cat}, {rad}\"")


def _xmatch_cat(csv_bytes: bytes, cat2: str, radius: float) -> bytes:
    boundary = uuid.uuid4().hex
    fields = {"request": "xmatch", "distMaxArcsec": str(radius), "RESPONSEFORMAT": "csv",
              "colRA1": "ra", "colDec1": "dec", "cat2": cat2}
    out = io.BytesIO()
    for k, v in fields.items():
        out.write(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode())
    out.write(f"--{boundary}\r\nContent-Disposition: form-data; name=\"cat1\"; filename=\"t.csv\"\r\n"
              f"Content-Type: text/csv\r\n\r\n".encode())
    out.write(csv_bytes)
    out.write(f"\r\n--{boundary}--\r\n".encode())
    return _fetch(C.CDS_XMATCH, data=out.getvalue(), timeout=900,
                  headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})


# ---------------------------------------------------------------------------
# 8. Literature: arXiv sources whose tables and figures we parse
# ---------------------------------------------------------------------------

#: arXiv identifier -> what we use from it.
ARXIV_SOURCES = {
    "1012.0163": "Tylenda et al. 2011, V1309 Sco: period law and the figures' period measurements",
    "1704.05502": "Molnar et al. 2017, KIC 9832227 merger prediction",
    "1809.02771": "Socia et al. 2018, KIC 9832227 eclipse times (Table 1)",
    "1909.00255": "Kovacs, Hartman & Bakos 2019, KIC 9832227 O-C compilation (Table 1)",
    "1804.04801": "Chen et al. 2018, contact-binary period-luminosity relations",
    "2005.08662": "Chen, Wang & Deng 2020, ZTF periodic variables",
    "0708.3020": "Rucinski 2007, the short-period end of contact binaries",
    "2303.04933": "Schaefer 2023, T CrB light curves 1842-2022 and the prediction",
    "2504.20592": "Merc et al. 2025, is T CrB late?",
    "2607.15245": "Pei et al. 2026, when will T CrB erupt? (review of published forecasts)",
    "2605.20991": "Pei et al. 2026, orbital-phase locking of T CrB's active phases",
}


def fetch_literature(manifest: dict) -> None:
    print("arXiv: source files of the papers we parse ...")
    for aid, note in ARXIV_SOURCES.items():
        dest = C.LIT_DIR / f"arXiv_{aid}.tar.gz"
        if not dest.exists():
            payload = _fetch(C.ARXIV_EPRINT + aid, timeout=300)
            _save(dest, payload)
            time.sleep(3.0)  # arXiv asks for at most one request every few seconds
        _record(manifest, dest, C.ARXIV_EPRINT + aid, note)


# ---------------------------------------------------------------------------

STEPS = {
    "ogle": fetch_ogle,
    "kepler": fetch_kepler,
    "tess": fetch_tess,
    "archival": fetch_kic_archival,
    "lit": fetch_literature,
    "ztf": fetch_ztf,
    "xmatch": fetch_xmatch,
    "allwise": fetch_allwise_xmatch,
    "chen": fetch_chen_xmatch,
    "veb": fetch_gaia_eb_vizier,
    "aavso": fetch_aavso,
}


def fetch_all(only: list[str] | None = None) -> dict:
    manifest = _load_manifest()
    try:
        for name, step in STEPS.items():
            if only and name not in only:
                continue
            step(manifest)
            _write_manifest(manifest)
    finally:
        _write_manifest(manifest)
    return manifest


if __name__ == "__main__":
    m = fetch_all(sys.argv[1:] or None)
    total = sum(f["bytes"] for f in m["files"].values())
    print(f"{len(m['files'])} files cached, {total / 1e6:.1f} MB")
