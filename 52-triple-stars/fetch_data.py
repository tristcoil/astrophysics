"""Download and cache every input this project uses (~0.3 GB, 10-20 minutes).

    python fetch_data.py            # everything
    python fetch_data.py dasch      # one group: msc, vizier, dasch, tess, kepler, gaia

Nothing is fetched twice: files already in ``data/raw`` are kept, so after one
run the whole analysis works offline. Each download is logged with its URL,
size, SHA-256 and date in ``data/manifest.json``.

Sources (citations in config.py and the README):

* Tokovinin's Multiple Star Catalog, release of 19 June 2026 (tarball).
* CDS/VizieR: Tokovinin 2017 (orbit alignment) and Tokovinin et al. 2006
  (tertiary companions to close spectroscopic binaries).
* DASCH DR7 (Harvard plate photometry, 1885-1992) through its public web API:
  SS Lac, HS Hya, V907 Sco and a few comparison stars around each.
* TESS: full-frame-image cutouts from MAST TESScut (TIC 290061484, V907 Sco,
  SS Lac) and SPOC 2-minute light curves (HS Hya).
* Kepler long-cadence light curves of HD 181068 and KOI-126 (MAST).
* Gaia DR3 (ESA archive): Proxima Centauri.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import config as C

#: Plain tool User-Agent (CDS refuses some browser UAs); DASCH's own client
#: identifies itself as "daschlab", which we mirror for its API.
UA = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _manifest() -> dict:
    if C.MANIFEST.exists():
        return json.loads(C.MANIFEST.read_text())
    return {}


def _log(key: str, url: str, path: Path, note: str = "") -> None:
    m = _manifest()
    m[key] = {
        "url": url,
        "file": str(path.relative_to(C.HERE)),
        "bytes": path.stat().st_size if path.is_file() else None,
        "sha256": _sha256(path) if path.is_file() else None,
        "retrieved_utc": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"),
        "note": note,
    }
    C.MANIFEST.write_text(json.dumps(m, indent=1, sort_keys=True))


def _get(url: str, *, retries: int = 5, data: bytes | None = None, timeout: int = 300,
         headers: dict | None = None) -> bytes:
    hdrs = {"User-Agent": UA}
    hdrs.update(headers or {})
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=hdrs)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            print(f"  attempt {attempt}/{retries} failed for {url[:120]}: {exc}", file=sys.stderr)
            time.sleep(4 * attempt)
    raise RuntimeError(f"could not download {url}")


def _slug(name: str) -> str:
    return name.lower().replace(" ", "_").replace("-", "_")


# ---------------------------------------------------------------------------
# Multiple Star Catalog and VizieR tables
# ---------------------------------------------------------------------------

def fetch_msc() -> None:
    dest_dir = C.RAW / "msc"
    tar_path = dest_dir / f"newmsc-{C.MSC_VERSION}.tar.gz"
    if (dest_dir / "sys.tsv").exists():
        return
    dest_dir.mkdir(parents=True, exist_ok=True)
    print("MSC: Tokovinin's Multiple Star Catalog", C.MSC_VERSION)
    tar_path.write_bytes(_get(C.MSC_URL))
    _log("msc:tarball", C.MSC_URL, tar_path, "Multiple Star Catalog, release 2026-06-19")
    with tarfile.open(tar_path) as tf:
        for member in tf.getmembers():
            name = Path(member.name).name
            if name in ("comp.tsv", "sys.tsv", "orb.tsv", "notes.tsv", "Readme") and member.isfile():
                (dest_dir / name).write_bytes(tf.extractfile(member).read())


def fetch_vizier() -> None:
    for cat, files in C.VIZIER_TABLES.items():
        folder = C.RAW / "vizier" / cat.replace("/", "_").replace("+", "")
        folder.mkdir(parents=True, exist_ok=True)
        for fname in files:
            dest = folder / fname
            if dest.exists():
                continue
            url = f"{C.CDS_FTP}/{cat}/{fname}"
            print("VizieR:", cat, fname)
            try:
                dest.write_bytes(_get(url, retries=2))
            except RuntimeError:
                # large tables are only served gzipped
                import gzip
                url += ".gz"
                dest.write_bytes(gzip.decompress(_get(url)))
            _log(f"vizier:{cat}/{fname}", url, dest)


# ---------------------------------------------------------------------------
# DASCH (Harvard plates)
# ---------------------------------------------------------------------------

def _dasch(endpoint: str, payload: dict) -> list:
    blob = _get(C.DASCH_API + endpoint, data=json.dumps(payload).encode(),
                headers={"User-Agent": "daschlab", "Accept": "application/json",
                         "Content-Type": "application/json"})
    out = json.loads(blob)
    if not isinstance(out, list):
        raise RuntimeError(f"DASCH {endpoint} failed: {str(out)[:300]}")
    return out


def _querycat(ra: float, dec: float, radius_arcsec: float) -> list[dict]:
    rows = _dasch("/dasch/dr7/querycat", {"refcat": C.DASCH_REFCAT, "ra_deg": ra, "dec_deg": dec,
                                          "radius_arcsec": radius_arcsec})
    hdr = rows[0].split(",")
    out = [dict(zip(hdr, line.split(","))) for line in rows[1:]]
    for d in out:
        d["dist_arcsec"] = math.hypot(float(d["dra_asec"]), float(d["ddec_asec"]))
    return out


#: comparison stars per target: similar brightness (+-0.6 mag), well observed
N_DASCH_COMPARISON = 3


def fetch_dasch() -> None:
    for name in ("SS Lac", "HS Hya", "V907 Sco"):
        sysd = C.SYSTEMS[name]
        folder = C.RAW / "dasch" / _slug(name)
        folder.mkdir(parents=True, exist_ok=True)
        cat_path = folder / "querycat.json"
        if not cat_path.exists():
            print("DASCH: catalogue around", name)
            cat = _querycat(sysd["ra"], sysd["dec"], 900.0)
            cat_path.write_text(json.dumps(cat))
            _log(f"dasch:{_slug(name)}:querycat", C.DASCH_API + "/dasch/dr7/querycat", cat_path,
                 f"APASS refcat within 900 arcsec of {name}")
        cat = json.loads(cat_path.read_text())
        target = min(cat, key=lambda d: d["dist_arcsec"])
        if target["dist_arcsec"] > 3.0:
            raise RuntimeError(f"no DASCH reference source within 3 arcsec of {name}")
        nm_t = int(target["num_matches"] or 0)
        mag_t = float(target["stdmag"])
        comps = [d for d in cat if d is not target and d["num_matches"] and float(d["stdmag"]) < 90
                 and abs(float(d["stdmag"]) - mag_t) < 0.6 and int(d["num_matches"]) > 0.7 * nm_t
                 and d["dist_arcsec"] > 60]
        comps.sort(key=lambda d: d["dist_arcsec"])
        picks = [("target", target)] + [(f"comp{k + 1}", d) for k, d in enumerate(comps[:N_DASCH_COMPARISON])]
        for label, d in picks:
            dest = folder / f"lc_{label}.json"
            if dest.exists():
                continue
            print(f"DASCH: light curve of {name} {label} ({d['ref_text']}, {d['stdmag']} mag, "
                  f"{d['num_matches']} matches)")
            lc = _dasch("/dasch/dr7/lightcurve", {"refcat": C.DASCH_REFCAT,
                                                   "gsc_bin_index": int(d["gsc_bin_index"]),
                                                   "ref_number": int(d["ref_number"])})
            dest.write_text(json.dumps({"source": d, "rows": lc}))
            _log(f"dasch:{_slug(name)}:{label}", C.DASCH_API + "/dasch/dr7/lightcurve", dest,
                 f"{d['ref_text']} ({d['stdmag']} mag), {len(lc) - 1} rows")


# ---------------------------------------------------------------------------
# TESS and Kepler (MAST)
# ---------------------------------------------------------------------------

def _tesscut(name: str, ra: float, dec: float, sector: int, folder: Path) -> None:
    done = list(folder.glob(f"tess-s{sector:04d}-*_astrocut.fits"))
    if done:
        return
    url = (f"{C.MAST_TESSCUT}?ra={ra}&dec={dec}&y={C.TESSCUT_SIZE}&x={C.TESSCUT_SIZE}"
           f"&sector={sector}")
    print(f"TESScut: {name} sector {sector}")
    blob = _get(url, timeout=900)
    with zipfile.ZipFile(io.BytesIO(blob)) as zf:
        names = [n for n in zf.namelist() if n.endswith(".fits")]
        if not names:
            raise RuntimeError(f"TESScut returned no FITS for {name} sector {sector}")
        for n in names:
            dest = folder / Path(n).name
            dest.write_bytes(zf.read(n))
            _log(f"tesscut:{_slug(name)}:s{sector}:{Path(n).name}", url, dest,
                 f"{C.TESSCUT_SIZE}x{C.TESSCUT_SIZE} px FFI cutout")


def _mast_products(target: str, mission: str, author: str, exptime: float) -> list[tuple[str, str]]:
    """(label, data URI) of the light-curve files lightkurve finds for a target."""
    import warnings
    warnings.filterwarnings("ignore")
    import lightkurve as lk
    res = lk.search_lightcurve(target, mission=mission, author=author, exptime=exptime)
    out = []
    for row in res.table:
        out.append((str(row["mission"]), str(row["dataURI"])))
    return out


def fetch_tess() -> None:
    for name, plan in C.TESS_PLAN.items():
        sysd = C.SYSTEMS[name]
        folder = C.RAW / "tess" / _slug(name)
        folder.mkdir(parents=True, exist_ok=True)
        if plan["kind"] == "tesscut":
            for s in plan["sectors"]:
                _tesscut(name, sysd["ra"], sysd["dec"], s, folder)
        else:
            have = {p.name for p in folder.glob("*.fits")}
            want = [s for s in plan["sectors"]
                    if not any(f"-s{s:04d}-" in h for h in have)]
            if not want:
                continue
            prods = _mast_products(f"TIC {sysd['tic']}", "TESS", "SPOC", 120)
            for label, uri in prods:
                fname = uri.split("/")[-1]
                sec = int(fname.split("-s")[1][:4])
                if sec not in want or (folder / fname).exists():
                    continue
                print(f"TESS SPOC: {name} {label}")
                dest = folder / fname
                url = C.MAST_DOWNLOAD + uri
                dest.write_bytes(_get(url))
                _log(f"tess:{_slug(name)}:s{sec}", url, dest, "SPOC 2-min light curve")


def fetch_kepler() -> None:
    for name, plan in C.KEPLER_PLAN.items():
        folder = C.RAW / "kepler" / _slug(name)
        folder.mkdir(parents=True, exist_ok=True)
        if len(list(folder.glob("*.fits"))) >= len(plan["quarters"]):
            continue
        prods = _mast_products(f"KIC {plan['kic']}", "Kepler", "Kepler", 1800)
        for label, uri in prods:
            q = int(label.split()[-1])
            fname = uri.split("/")[-1]
            if q not in plan["quarters"] or (folder / fname).exists():
                continue
            print(f"Kepler: {name} {label}")
            dest = folder / fname
            url = C.MAST_DOWNLOAD + uri
            dest.write_bytes(_get(url))
            _log(f"kepler:{_slug(name)}:q{q}", url, dest, "Kepler long-cadence light curve")


# ---------------------------------------------------------------------------
# Gaia DR3: Proxima
# ---------------------------------------------------------------------------

GAIA_COLUMNS = ("source_id, ref_epoch, ra, ra_error, dec, dec_error, parallax, parallax_error, "
                "pmra, pmra_error, pmdec, pmdec_error, ra_dec_corr, ra_parallax_corr, "
                "dec_parallax_corr, parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr, "
                "ra_pmra_corr, ra_pmdec_corr, dec_pmra_corr, dec_pmdec_corr, ruwe, "
                "phot_g_mean_mag, radial_velocity, radial_velocity_error")


def fetch_gaia() -> None:
    dest = C.RAW / "gaia_dr3_proxima.csv"
    if dest.exists():
        return
    adql = f"SELECT {GAIA_COLUMNS} FROM gaiadr3.gaia_source WHERE source_id = {C.PROXIMA_GAIA_DR3}"
    query = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                    "QUERY": adql}).encode()
    print("Gaia DR3: Proxima Centauri")
    for url in (C.GAIA_TAP, C.GAIA_TAP_MIRROR):
        try:
            text = _get(url, data=query, retries=3).decode()
        except RuntimeError:
            continue
        if text.count("\n") >= 2:
            dest.write_text(text)
            _log("gaia:dr3:proxima", url, dest, adql)
            return
    raise RuntimeError("Gaia DR3 query for Proxima failed on both servers")


def fetch_gaia_fields() -> None:
    """Gaia DR3 stars around each TESS cutout target (G < 19, 4 arcmin): the
    neighbours whose light dilutes the eclipses in our TESS apertures."""
    for name, plan in C.TESS_PLAN.items():
        if plan["kind"] != "tesscut":
            continue
        dest = C.RAW / "tess" / _slug(name) / "gaia_dr3_field.csv"
        if dest.exists():
            continue
        s = C.SYSTEMS[name]
        adql = ("SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax FROM gaiadr3.gaia_source "
                f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {s['ra']!r}, {s['dec']!r}, "
                f"{4 / 60!r})) AND phot_g_mean_mag < 19")
        query = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                        "QUERY": adql}).encode()
        print("Gaia DR3: field around", name)
        for url in (C.GAIA_TAP, C.GAIA_TAP_MIRROR):
            try:
                text = _get(url, data=query, retries=3).decode()
            except RuntimeError:
                continue
            if text.count("\n") >= 2:
                dest.write_text(text)
                _log(f"gaia:dr3:field:{_slug(name)}", url, dest, adql)
                break
        else:
            raise RuntimeError(f"Gaia field query failed for {name}")


# ---------------------------------------------------------------------------
# Hipparcos epoch photometry (1989-1993) and ASAS-3 (2000-2009)
# ---------------------------------------------------------------------------

HIP_EPOPHOT = f"{C.CDS_FTP}/I/239/epophot"
#: Hipparcos numbers and ASAS-3 designations of the switching systems
HIP_IDS = {"SS Lac": 108981, "HS Hya": 50966}      # V907 Sco (HD 163302) is not in Hipparcos
ASAS_IDS = {"HS Hya": "102437-1905.5", "V907 Sco": "175656-3445.0"}
ASAS_URL = "https://www.astrouw.edu.pl/cgi-asas/asas_cgi_get_data?{},asas3"


def fetch_hipparcos_epoch() -> None:
    """One star at a time from the block-gzipped hep.gz, using its index of
    byte offsets (the CDS 'hipepo' recipe) and an HTTP range request."""
    import gzip
    import zlib
    folder = C.RAW / "hipparcos"
    folder.mkdir(parents=True, exist_ok=True)
    idx_path = folder / "hep.gz.idx"
    if not idx_path.exists():
        idx_path.write_bytes(_get(f"{HIP_EPOPHOT}/hep.gz.idx"))
        _log("hipparcos:hep.idx", f"{HIP_EPOPHOT}/hep.gz.idx", idx_path, "byte offsets of each star in hep.gz")
    index = {}
    for line in idx_path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            index[k.strip()] = int(v.strip())
    for name, hip in HIP_IDS.items():
        dest = folder / f"hip{hip:06d}_epoch.txt"
        if dest.exists():
            continue
        off = index[f"{hip:06d}"]
        url = f"{HIP_EPOPHOT}/hep.gz"
        print(f"Hipparcos epoch photometry: {name} (HIP {hip}) at byte {off}")
        blob = _get(url, headers={"Range": f"bytes={off - 1}-{off - 1 + 200000}"})
        d = zlib.decompressobj(16 + zlib.MAX_WBITS)
        text = d.decompress(blob).decode("latin-1")
        # the star's records end at the first blank line
        block = text.split("\n\n")[0]
        if not block.lstrip().startswith(f"{hip}") and f"{hip}" not in block[:200]:
            raise RuntimeError(f"unexpected Hipparcos block for HIP {hip}: {block[:120]!r}")
        dest.write_text(block + "\n")
        _log(f"hipparcos:epoch:{hip}", url, dest, f"HIP {hip} ({name}); bytes {off - 1}-{off + 199999}")


def fetch_asas() -> None:
    folder = C.RAW / "asas"
    folder.mkdir(parents=True, exist_ok=True)
    for name, asas in ASAS_IDS.items():
        dest = folder / f"{_slug(name)}_asas3.txt"
        if dest.exists():
            continue
        url = ASAS_URL.format(asas)
        print(f"ASAS-3: {name} ({asas})")
        text = _get(url).decode("latin-1")
        if "LIGHT CURVE BEGINS" not in text:
            raise RuntimeError(f"unexpected ASAS answer for {name}: {text[:200]!r}")
        dest.write_text(text)
        _log(f"asas3:{_slug(name)}", url, dest, f"ASAS-3 V-band light curve {asas}")


GROUPS = {"msc": fetch_msc, "vizier": fetch_vizier, "dasch": fetch_dasch, "tess": fetch_tess,
          "kepler": fetch_kepler, "gaia": fetch_gaia, "gaia_fields": fetch_gaia_fields,
          "hipparcos": fetch_hipparcos_epoch, "asas": fetch_asas}


def ensure_data(groups=None) -> None:
    for g in (groups or GROUPS):
        GROUPS[g]()


if __name__ == "__main__":
    ensure_data(sys.argv[1:] or None)
    m = _manifest()
    total = sum((v.get("bytes") or 0) for v in m.values())
    print(f"done: {len(m)} items, {total / 1e6:.1f} MB in data/manifest.json")
