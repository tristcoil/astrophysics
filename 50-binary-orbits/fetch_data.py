"""Download and cache every input this project uses (~0.6 GB, 570 MB of it Gaia's scanning law).

    python fetch_data.py

1. **ORB6**, the Sixth Catalog of Orbits of Visual Binary Stars (GSU/USNO),
   and its format description.
2. **VizieR catalogues** from the CDS: Bond et al. 2017's 2,367 measures of
   Sirius B (J/ApJ/840/70); Eker et al. 2018's eclipsing-binary stars
   (J/MNRAS/479/5491); Hipparcos, new reduction (I/311).
3. **arXiv TeX sources** of the papers whose tables we read or check against
   (only the .tex files are kept), and Barstow et al. 2005 (PDF only).
4. **Mamajek's dwarf-star table** and the **Montreal white dwarf models**.
5. **Gaia DR3**: Sirius B's full astrometric record; CDS XMatch of every
   grade 1-2 ORB6 orbit (and the gallery systems) against Gaia DR3 within 60";
   the extra quality columns, Hipparcos cross-match and non-single-star
   solutions for the matched sources, by source_id; the DR3 commanded
   scanning law and DPAC's DR3 astrometric data gaps.

Files land in ``data/raw/``; nothing is re-downloaded once present, so after
one run the analysis works offline. ``data/manifest.json`` records each file's
URL, size, SHA-256 and download date.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import io
import json
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

import numpy as np
import pandas as pd

import config as C

UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64) astro-atlas hobby project"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _fetch(url: str, *, data: bytes | None = None, headers: dict | None = None,
           timeout: float = 300.0, retries: int = 4, ua: str = UA_TOOL) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": ua, **(headers or {})})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (400, 403, 404, 410):
                raise RuntimeError(f"HTTP {exc.code} for {url}: {exc.read()[:500]!r}") from exc
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
    h = hashlib.sha256()
    with open(dest, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    sha = h.hexdigest()
    key = str(dest.relative_to(C.DATA_DIR))
    old = manifest["files"].get(key, {})
    manifest["files"][key] = {
        "url": url, "note": note, "bytes": dest.stat().st_size, "sha256": sha,
        "fetched": old["fetched"] if old.get("sha256") == sha
        else dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def download(url: str, dest: Path, manifest: dict, note: str, **kw) -> Path:
    if not dest.exists():
        print(f"  downloading {url}")
        _save(dest, _fetch(url, **kw))
    _record(manifest, dest, url, note)
    return dest


def _stream(url: str, dest: Path, manifest: dict, note: str) -> Path:
    """Large file: stream to disk in blocks (resumable by re-running)."""
    if not dest.exists():
        print(f"  streaming {url} (large)")
        tmp = dest.with_suffix(dest.suffix + ".part")
        done = tmp.stat().st_size if tmp.exists() else 0
        req = urllib.request.Request(url, headers={"User-Agent": UA_TOOL, "Range": f"bytes={done}-"})
        with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "ab" if done else "wb") as fh:
            for block in iter(lambda: resp.read(1 << 22), b""):
                fh.write(block)
        tmp.replace(dest)
    _record(manifest, dest, url, note)
    return dest


# ---------------------------------------------------------------------------
# 1-4. Catalogues, papers, models
# ---------------------------------------------------------------------------

def fetch_orb6(manifest):
    print("ORB6 ...")
    download(C.ORB6_URL, C.ORB6_DIR / "orb6orbits.txt", manifest,
             "Sixth Catalog of Orbits of Visual Binary Stars (GSU/USNO), snapshot")
    download(C.ORB6_FORMAT_URL, C.ORB6_DIR / "orb6format.txt", manifest, "ORB6 format description")


def fetch_cds(manifest):
    print("VizieR catalogues (CDS) ...")
    for key, (cat, files) in C.CDS_CATALOGUES.items():
        base = C.HIP_DIR if key == "hip2" else C.LIT_DIR / key
        for f in files:
            download(f"{C.CDS_FTP}{cat}/{f}", base / f, manifest, f"VizieR {cat}: {f}")


def fetch_arxiv(manifest):
    print("arXiv TeX sources ...")
    for name, aid in C.ARXIV_SOURCES.items():
        d = C.LIT_DIR / name
        if d.is_dir() and any(d.glob("*.tex")):
            for tex in sorted(d.glob("*.tex")):
                _record(manifest, tex, f"https://arxiv.org/e-print/{aid}", f"arXiv:{aid} TeX source")
            continue
        blob = _fetch(f"https://arxiv.org/e-print/{aid}")
        d.mkdir(parents=True, exist_ok=True)
        try:
            with tarfile.open(fileobj=io.BytesIO(blob)) as tar:
                for m in tar.getmembers():
                    if m.isfile() and m.name.endswith(".tex") and "/" not in m.name.strip("./"):
                        (d / Path(m.name).name).write_bytes(tar.extractfile(m).read())
        except tarfile.ReadError:
            import gzip
            (d / "main.tex").write_bytes(gzip.decompress(blob))
        for tex in sorted(d.glob("*.tex")):
            _record(manifest, tex, f"https://arxiv.org/e-print/{aid}", f"arXiv:{aid} TeX source")
        time.sleep(3)
    # Barstow et al. 2005 has no TeX on arXiv: the e-print is the PDF itself.
    pdf = C.LIT_DIR / "barstow2005_src" / "barstow2005.pdf"
    download("https://arxiv.org/e-print/astro-ph/0506600", pdf, manifest,
             "Barstow et al. 2005, MNRAS 362, 1134 (PDF; Table 2 read with pdftotext)")


def fetch_models(manifest):
    print("Mamajek table and Montreal white dwarf models ...")
    download(C.MAMAJEK_URL, C.LIT_DIR / "mamajek" / "EEM_dwarf_UBVIJHK_colors_Teff.txt", manifest,
             "Pecaut & Mamajek 2013 dwarf sequence, version 2022.04.16", ua=UA_BROWSER)
    download(C.MONTREAL_BASE + C.MONTREAL_TABLE_DA, C.MONTREAL_DIR / "Table_DA", manifest,
             "Montreal synthetic photometry, pure-H white dwarf atmospheres", ua=UA_BROWSER)
    for m in C.MONTREAL_MASSES:
        name = f"seq_{m:03d}_thick.txt"
        download(f"{C.MONTREAL_BASE}CoolingModels/{name}", C.MONTREAL_DIR / name, manifest,
                 f"Bedard et al. 2020 cooling sequence, {m / 100:.2f} Msun, thick H", ua=UA_BROWSER)


# ---------------------------------------------------------------------------
# 5. Gaia DR3
# ---------------------------------------------------------------------------

GAIA_DETAIL_COLS = (
    "source_id", "ra", "dec", "ref_epoch", "parallax", "parallax_error", "pmra", "pmra_error",
    "pmdec", "pmdec_error", "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",
    "ruwe", "astrometric_params_solved", "astrometric_excess_noise", "astrometric_excess_noise_sig",
    "astrometric_n_good_obs_al", "astrometric_n_bad_obs_al", "astrometric_chi2_al", "astrometric_gof_al",
    "visibility_periods_used", "matched_transits", "ipd_frac_multi_peak", "ipd_gof_harmonic_amplitude",
    "ipd_frac_odd_win", "duplicated_source", "non_single_star", "nu_eff_used_in_astrometry",
    "pseudocolour", "ecl_lat", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
    "phot_bp_rp_excess_factor", "phot_variable_flag", "radial_velocity", "radial_velocity_error",
)


def gaia_sync(adql: str, fmt: str = "csv", service: str = "esa") -> bytes:
    """One synchronous ADQL query, ESA archive first, ARI Heidelberg mirror as fallback."""
    urls = {"esa": C.GAIA_TAP_ESA + "/sync", "ari": C.GAIA_TAP_ARI + "/sync"}
    order = [service] + [s for s in urls if s != service]
    last = None
    for s in order:
        data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt,
                                       "QUERY": adql, "MAXREC": 3_000_000}).encode()
        try:
            body = _fetch(urls[s], data=data, timeout=600, retries=2)
        except RuntimeError as exc:
            last = exc
            continue
        if b"OVERFLOW" in body[:4000] or b'value="ERROR"' in body[:4000]:
            last = RuntimeError(body[:1000])
            continue
        return body
    raise RuntimeError(f"Gaia query failed everywhere: {last}")


def _id_query(table_cols: str, table: str, ids, key="source_id", extra="") -> pd.DataFrame:
    """Query a Gaia table for a list of ids in chunks of 500 (no table upload needed)."""
    out = []
    ids = [int(i) for i in ids]
    for k in range(0, len(ids), 500):
        chunk = ",".join(str(i) for i in ids[k:k + 500])
        adql = f"SELECT {table_cols} FROM {table} WHERE {key} IN ({chunk}) {extra}"
        out.append(pd.read_csv(io.BytesIO(gaia_sync(adql))))
    df = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    df.columns = [c.lower() for c in df.columns]   # the archive's column case varies
    return df


def fetch_gaia_sirius(manifest):
    print("Gaia DR3: Sirius B ...")
    dest = C.GAIA_DIR / "sirius_b_dr3.csv"
    if not dest.exists():
        df = _id_query(", ".join(GAIA_DETAIL_COLS), "gaiadr3.gaia_source", [C.SIRIUS_B_GAIA_DR3])
        nb = pd.read_csv(io.BytesIO(gaia_sync(
            "SELECT n.dr2_source_id, n.angular_distance, n.magnitude_difference, g.parallax AS dr2_parallax, "
            "g.parallax_error AS dr2_parallax_error, g.pmra AS dr2_pmra, g.pmdec AS dr2_pmdec, "
            "g.pmra_error AS dr2_pmra_error, g.pmdec_error AS dr2_pmdec_error, g.ref_epoch AS dr2_ref_epoch "
            f"FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr2.gaia_source AS g ON g.source_id = n.dr2_source_id "
            f"WHERE n.dr3_source_id = {C.SIRIUS_B_GAIA_DR3}")))
        nb.columns = [c.lower() for c in nb.columns]
        if len(nb):
            best = nb.sort_values("angular_distance").iloc[0]
            for k, v in best.items():
                df[k] = v
        df.to_csv(dest, index=False)
    _record(manifest, dest, C.GAIA_TAP_ESA, "Gaia DR3 gaia_source row of Sirius B, plus its DR2 counterpart")


def fetch_scanlaw(manifest):
    print("Gaia DR3 commanded scanning law and data gaps ...")
    _stream(C.GAIA_SCANLAW_URL, C.GAIA_DIR / "CommandedScanLaw_001.csv.gz", manifest,
            "Gaia (E)DR3 commanded scanning law, 10-s sampling (ESA)")
    download(C.GAIA_GAPS_URL, C.GAIA_DIR / "dr3_astrometry_gaps.csv", manifest,
             "DPAC (E)DR3 astrometric data gaps (OBMT revolutions)")


def orb6_targets():
    """The ORB6 rows we need Gaia for: every grade 1-2 orbit plus the gallery systems."""
    import orb6
    o = orb6.read_orb6()
    keep = o["grade"].isin(C.MLR_GRADES)
    for wds, disc, _ in C.GALLERY:
        keep |= (o["wds"] == wds) & (o["disc"].str.replace(" ", "") == disc.replace(" ", ""))
    t = o[keep].copy()
    t["oid"] = t.index
    return t


def fetch_gaia_orb6(manifest):
    print("Gaia DR3 x ORB6 (CDS XMatch, then ESA by source_id) ...")
    t = orb6_targets()
    dest = C.GAIA_DIR / "orb6_xmatch_gaiadr3.csv"
    if not dest.exists():
        src = t[["oid", "ra_deg", "dec_deg"]].rename(columns={"ra_deg": "ra", "dec_deg": "dec"})
        boundary = uuid.uuid4().hex
        csv = src.to_csv(index=False).encode()
        fields = {"request": "xmatch", "distMaxArcsec": str(C.GAIA_SEARCH_RADIUS_ARCSEC),
                  "RESPONSEFORMAT": "csv", "colRA1": "ra", "colDec1": "dec", "cat2": "vizier:I/355/gaiadr3"}
        body = b""
        for k, v in fields.items():
            body += f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode()
        body += (f"--{boundary}\r\nContent-Disposition: form-data; name=\"cat1\"; filename=\"orb6.csv\"\r\n"
                 f"Content-Type: text/csv\r\n\r\n").encode() + csv + f"\r\n--{boundary}--\r\n".encode()
        blob = _fetch(C.CDS_XMATCH, data=body, headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
                      timeout=900)
        _save(dest, blob)
    _record(manifest, dest, C.CDS_XMATCH, "CDS XMatch: ORB6 grade 1-2 + gallery positions vs Gaia DR3 (60\")")
    xm = pd.read_csv(dest)

    # detail columns (ESA names) for the brighter candidates only (G < 17)
    det = C.GAIA_DIR / "orb6_gaia_details.csv"
    if not det.exists():
        ids = xm.loc[xm["Gmag"] < 17.0, "Source"].dropna().astype("int64").unique()
        _id_query(", ".join(GAIA_DETAIL_COLS), "gaiadr3.gaia_source", ids).to_csv(det, index=False)
    _record(manifest, det, C.GAIA_TAP_ESA, "Gaia DR3 gaia_source quality columns for matched candidates")

    nss = C.GAIA_DIR / "orb6_gaia_nss.csv"
    if not nss.exists():
        ids = xm.loc[xm["Gmag"] < 17.0, "Source"].dropna().astype("int64").unique()
        _id_query("source_id, nss_solution_type, parallax, parallax_error, period, period_error, "
                  "eccentricity, goodness_of_fit, significance", "gaiadr3.nss_two_body_orbit", ids
                  ).to_csv(nss, index=False)
    _record(manifest, nss, C.GAIA_TAP_ESA, "Gaia DR3 nss_two_body_orbit rows for matched candidates")

    hbn = C.GAIA_DIR / "orb6_hip_best_neighbour.csv"
    if not hbn.exists():
        hips = t.loc[t["hip"] > 0, "hip"].astype(int).unique()
        _id_query("source_id, original_ext_source_id, angular_distance, number_of_neighbours",
                  "gaiadr3.hipparcos2_best_neighbour", hips, key="original_ext_source_id").to_csv(hbn, index=False)
    _record(manifest, hbn, C.GAIA_TAP_ESA, "Gaia DR3 hipparcos2_best_neighbour for ORB6 HIP stars")

    hip1 = C.GAIA_DIR / "orb6_hipparcos1.csv"
    if not hip1.exists():
        hips = t.loc[t["hip"] > 0, "hip"].astype(int).unique()
        _id_query("hip, vmag, b_v, v_i, sptype, ccdm, n_ccdm, nsys, ncomp, multflag, hpmag, theta, rho, dhp",
                  "public.hipparcos", hips, key="hip").to_csv(hip1, index=False)
    _record(manifest, hip1, C.GAIA_TAP_ESA, "Hipparcos (1997) main catalogue rows: V, B-V, spectral type")


#: Gallery stars whose Gaia data the positional cross-match cannot reach: 61 Cyg
#: (5.2"/yr: 83" from its J2000 position by 2016) and xi UMa (two-parameter
#: solutions only in DR3, so its DR2 parallax is used).
GALLERY_DR3_IDS = (1872046609345556480, 1872046574983497216,      # 61 Cyg A, B
                   756853643638639104, 756853643637996160)        # xi UMa A, B


def fetch_gaia_gallery(manifest):
    print("Gaia DR3/DR2: gallery stars ...")
    dest = C.GAIA_DIR / "gallery_gaia.csv"
    if not dest.exists():
        ids = ",".join(str(i) for i in GALLERY_DR3_IDS)
        df = pd.read_csv(io.BytesIO(gaia_sync(
            "SELECT g.source_id, g.phot_g_mean_mag, g.parallax, g.parallax_error, g.ruwe, "
            "g.astrometric_params_solved, n.dr2_source_id, d.parallax AS dr2_parallax, "
            "d.parallax_error AS dr2_parallax_error, d.phot_g_mean_mag AS dr2_g "
            "FROM gaiadr3.gaia_source AS g "
            "LEFT OUTER JOIN gaiadr3.dr2_neighbourhood AS n ON n.dr3_source_id = g.source_id "
            "LEFT OUTER JOIN gaiadr2.gaia_source AS d ON d.source_id = n.dr2_source_id "
            f"WHERE g.source_id IN ({ids})")))
        df.columns = [c.lower() for c in df.columns]
        df.to_csv(dest, index=False)
    _record(manifest, dest, C.GAIA_TAP_ESA, "Gaia DR3 (+DR2 neighbours) rows of 61 Cyg A, B and xi UMa A, B")


def fetch_vizier_rows(manifest):
    """Small VizieR lookups: the Hipparcos DMSA/O orbit of Sirius (adopted from GG78)."""
    dest = C.HIP_DIR / "sirius_hip_dm_o.tsv"
    url = f"{C.VIZIER_ASU}?-source=I/239/hip_dm_o&HIP={C.SIRIUS_HIP}&-out.all"
    download(url, dest, manifest, "Hipparcos 1997 DMSA/O row for Sirius (photocentre orbit adopted from GG78)")


def fetch_all():
    manifest = _load_manifest()
    try:
        fetch_orb6(manifest)
        fetch_cds(manifest)
        fetch_arxiv(manifest)
        fetch_models(manifest)
        fetch_vizier_rows(manifest)
        fetch_gaia_sirius(manifest)
        fetch_scanlaw(manifest)
        fetch_gaia_orb6(manifest)
        fetch_gaia_gallery(manifest)
    finally:
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    total = sum(v["bytes"] for v in manifest["files"].values())
    print(f"All inputs cached: {len(manifest['files'])} files, {total / 1e6:.0f} MB")


if __name__ == "__main__":
    fetch_all()
