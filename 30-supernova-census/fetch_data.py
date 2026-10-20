"""Download and cache everything the analysis needs.

Everything is free and public; nothing needs an account.

1. **The ZTF Bright Transient Survey catalogue** (Fremling et al. 2020; Perley
   et al. 2020), from its Sample Explorer: every candidate ever saved, with
   classification, redshift, peak magnitude, timescales, host photometry and
   the survey's own quality/purity flags. Fetched twice: everything, and the
   subset inside SDSS imaging.
2. **ZTF light curves from the ALeRCE broker** (Förster et al. 2021) for every
   Type II supernova (the plateau analysis) and the brightest examples of each
   class (the gallery and the explorer).
3. **SDSS DR17 galaxies around each supernova** (SkyServer): photometry and
   shapes of every galaxy within 30 kpc, plus spectra and the MPA-JHU stellar
   masses and star-formation rates where they exist.
4. **The SDSS galaxy population** (MPA-JHU, z < 0.1): the background against
   which the hosts are judged.
5. **Cosmicflows-4** (Tully et al. 2023, CDS): redshift-independent distances
   for the nearest hosts, where a redshift says little about distance.
6. **The SFD dust map** (NASA LAMBDA): how much volume Galactic dust hides.
7. **arXiv LaTeX sources** of Smartt (2015) and Davies & Beasor (2018), whose
   progenitor tables we parse ourselves.
8. **ZTF's public exposure history** (IRSA), optional: if the service answers,
   rates.py measures the survey's active sky area itself.

Everything is cached under ``data/`` as served; re-running fetches only what
is missing.

    python fetch_data.py            # everything
    python fetch_data.py --only bts alerce
    python fetch_data.py --refresh-bts   # take a new snapshot of the catalogue
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import tarfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

import config as C

USER_AGENT = "astrophysics-atlas-L2/1.0 (hobby research; python-requests)"
#: (connect, read) timeouts: a short connect timeout lets requests fall back to
#: IPv4 quickly when a host publishes an IPv6 address this machine can't reach.
TIMEOUT = (15.0, 300.0)

_session_local = threading.local()


def _session() -> requests.Session:
    s = getattr(_session_local, "s", None)
    if s is None:
        s = requests.Session()
        s.headers["User-Agent"] = USER_AGENT
        _session_local.s = s
    return s


def http_get(url: str, *, params=None, tries: int = 5, stream: bool = False,
             ok_404: bool = False) -> requests.Response | None:
    """GET with retries and exponential backoff on network errors and 5xx."""
    delay = 3.0
    for attempt in range(1, tries + 1):
        try:
            r = _session().get(url, params=params, timeout=TIMEOUT, stream=stream)
            if r.status_code == 404 and ok_404:
                return None
            if r.status_code >= 500 or r.status_code == 429:
                raise requests.HTTPError(f"HTTP {r.status_code}", response=r)
            r.raise_for_status()
            return r
        except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
            if attempt == tries:
                raise
            print(f"    retry {attempt}/{tries - 1} after {exc.__class__.__name__}: {exc}")
            time.sleep(delay)
            delay *= 2
    return None


def _write_meta(path, **meta) -> None:
    meta["fetched_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    path.with_suffix(path.suffix + ".meta.json").write_text(json.dumps(meta, indent=1))


# ---------------------------------------------------------------------------
# 1. The BTS catalogue
# ---------------------------------------------------------------------------

BTS_ALL = C.BTS_DIR / "bts_all.csv"
BTS_SDSS = C.BTS_DIR / "bts_sdss_field.csv"


def fetch_bts(refresh: bool = False) -> None:
    for path, query, label in ((BTS_ALL, C.BTS_QUERY_ALL, "every candidate"),
                               (BTS_SDSS, C.BTS_QUERY_SDSS, "candidates in SDSS fields")):
        if path.exists() and not refresh:
            print(f"  [cached] {path.name}")
            continue
        print(f"  BTS explorer: {label} ...")
        r = http_get(C.BTS_EXPLORER_URL, params=query)
        text = r.text
        # The first line of the export is blank; the second is the header.
        n = sum(1 for line in text.splitlines() if line.startswith("ZTF"))
        if n < 1000:
            raise RuntimeError(f"BTS export looks wrong ({n} rows):\n{text[:500]}")
        path.write_text(text)
        _write_meta(path, url=r.url, rows=n)
        print(f"    {n:,} rows -> {path.relative_to(C.PROJECT_DIR)}")


def read_bts_csv(path=BTS_ALL) -> pd.DataFrame:
    """The explorer CSV as strings (limits like '>21.0' are kept for bts.py)."""
    return pd.read_csv(path, skiprows=1, dtype=str, keep_default_na=False)


# ---------------------------------------------------------------------------
# 2. ALeRCE light curves
# ---------------------------------------------------------------------------

#: TNS labels that are not extragalactic transients (no gallery for these).
NOT_TRANSIENT = {"-", "CV", "CV?", "AGN", "AGN?", "star", "varstar", "Other", "other", "bogus"}


def select_lightcurve_targets(bts: pd.DataFrame) -> list[str]:
    """Every Type II supernova, plus the brightest well-observed of each class.

    The gallery picks are rule-based, not hand-picked: for each TNS class, the
    ``LC_GALLERY_PER_CLASS`` brightest events at peak among those that pass
    the survey's quality cut and have a redshift.
    """
    ids = set(bts.loc[bts["type"].isin(C.LC_TYPE_II_CLASSES), "ZTFID"])
    peak = pd.to_numeric(bts["peakmag"], errors="coerce")
    good = (bts["qualcode"] == "111111") & (bts["redshift"] != "-") & peak.notna()
    for cls, grp in bts[good].groupby("type"):
        if cls in NOT_TRANSIENT:
            continue
        top = grp.assign(_m=peak[grp.index]).nsmallest(C.LC_GALLERY_PER_CLASS, "_m")
        ids.update(top["ZTFID"])
    return sorted(ids)


def _fetch_one_lc(oid: str) -> tuple[str, str]:
    path = C.ALERCE_DIR / f"{oid}.json"
    if path.exists():
        return oid, "cached"
    r = http_get(C.ALERCE_LC_URL.format(oid=oid), ok_404=True)
    if r is None:
        path.write_text(json.dumps({"oid": oid, "missing": True}))
        return oid, "missing"
    data = r.json()
    # Keep only the columns the analysis uses: the files stay ~4x smaller.
    keep_det = ("mjd", "fid", "magpsf", "sigmapsf", "isdiffpos", "diffmaglim", "rb", "drb",
                "dubious", "distnr")
    keep_non = ("mjd", "fid", "diffmaglim")
    slim = {
        "oid": oid,
        "detections": [{k: d.get(k) for k in keep_det} for d in data.get("detections", [])],
        "non_detections": [{k: d.get(k) for k in keep_non} for d in data.get("non_detections", [])],
    }
    path.write_text(json.dumps(slim, separators=(",", ":")))
    return oid, "ok"


def fetch_alerce(bts: pd.DataFrame) -> None:
    ids = select_lightcurve_targets(bts)
    todo = [i for i in ids if not (C.ALERCE_DIR / f"{i}.json").exists()]
    print(f"  ALeRCE: {len(ids):,} light curves wanted, {len(todo):,} to fetch")
    if not todo:
        return
    t0, done, failed = time.time(), 0, []
    with ThreadPoolExecutor(max_workers=C.ALERCE_WORKERS) as pool:
        futures = {pool.submit(_fetch_one_lc, oid): oid for oid in todo}
        for fut in as_completed(futures):
            oid = futures[fut]
            try:
                fut.result()
            except Exception as exc:  # keep going; a re-run picks these up
                failed.append(oid)
                print(f"    {oid}: {exc}")
            done += 1
            if done % 100 == 0:
                rate = done / (time.time() - t0)
                print(f"    {done:,}/{len(todo):,} ({rate:.1f}/s, "
                      f"~{(len(todo) - done) / rate / 60:.0f} min left)")
    print(f"    done in {(time.time() - t0) / 60:.1f} min; {len(failed)} failed")


# ---------------------------------------------------------------------------
# 3. SDSS galaxies around each supernova
# ---------------------------------------------------------------------------

def _skyserver(sql: str) -> pd.DataFrame:
    r = http_get(C.SKYSERVER_URL, params={"cmd": sql, "format": "csv"})
    text = r.text
    if text.startswith("#Table1"):
        text = text.split("\n", 1)[1]
    if not text.strip():
        return pd.DataFrame()
    # 64-bit SDSS IDs lose digits as floats, which pandas uses for any column
    # with a blank (galaxies without spectra): keep the IDs as text.
    return pd.read_csv(io.StringIO(text), dtype={"objID": str, "specObjID": str})


def host_search_radius_arcmin(z: np.ndarray) -> np.ndarray:
    """min(90", 30 kpc at the SN redshift), floored at 20" (config)."""
    from astropy.cosmology import FlatLambdaCDM
    cosmo = FlatLambdaCDM(H0=C.H0, Om0=C.OMEGA_M)
    z = np.clip(np.asarray(z, float), 1e-4, None)
    kpc_per_arcsec = cosmo.kpc_proper_per_arcmin(z).value / 60.0
    r = np.minimum(C.HOST_SEARCH_MAX_ARCSEC, C.HOST_SEARCH_MAX_KPC / kpc_per_arcsec)
    return np.maximum(r, C.HOST_SEARCH_MIN_ARCSEC) / 60.0


HOST_SQL = """
SELECT v.id AS sn_index, n.objID, n.distance*60.0 AS sep_arcsec,
  p.ra, p.dec, p.type, p.petroMag_r, p.extinction_r, p.modelMag_u, p.modelMag_g,
  p.modelMag_r, p.modelMag_i, p.modelMag_z, p.petroR50_r, p.petroR90_r,
  p.expRad_r, p.expAB_r, p.expPhi_r, p.deVRad_r, p.deVAB_r, p.deVPhi_r, p.fracDeV_r,
  s.specObjID, s.z AS spec_z, s.zWarning, s.class AS spec_class,
  e.lgm_tot_p16, e.lgm_tot_p50, e.lgm_tot_p84, e.sfr_tot_p16, e.sfr_tot_p50, e.sfr_tot_p84,
  e.specsfr_tot_p50, e.bptclass
FROM (VALUES {values}) AS v(id, ra, dec, r)
CROSS APPLY dbo.fGetNearbyObjEq(v.ra, v.dec, v.r) AS n
JOIN PhotoPrimary p ON p.objID = n.objID
LEFT JOIN SpecObj s ON s.bestObjID = p.objID AND s.sciencePrimary = 1
LEFT JOIN galSpecExtra e ON e.specObjID = s.specObjID
WHERE p.type = 3 AND p.modelMag_r < {rmax}
"""


def fetch_sdss_hosts(bts: pd.DataFrame) -> None:
    """Galaxies near every classified transient that lies in SDSS imaging."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    in_sdss = set(read_bts_csv(BTS_SDSS)["ZTFID"])
    z = pd.to_numeric(bts["redshift"], errors="coerce")
    want = bts["type"].map(C.CLASS_MAP).notna() & z.notna() & bts["ZTFID"].isin(in_sdss)
    tab = bts[want].reset_index(drop=True)
    c = SkyCoord(tab["RA"].values, tab["Dec"].values, unit=(u.hourangle, u.deg))
    radius = host_search_radius_arcmin(z[want].values)
    pos = pd.DataFrame({"ZTFID": tab["ZTFID"], "ra": c.ra.deg, "dec": c.dec.deg,
                        "r_arcmin": radius})
    pos.to_csv(C.SDSS_HOST_DIR / "positions.csv", index=False)
    batches = [pos.iloc[i:i + C.SKYSERVER_BATCH] for i in range(0, len(pos), C.SKYSERVER_BATCH)]
    print(f"  SDSS hosts: {len(pos):,} transients in SDSS fields, {len(batches)} batches")
    for k, b in enumerate(batches):
        path = C.SDSS_HOST_DIR / f"batch_{k:04d}.csv"
        if path.exists():
            continue
        values = ",".join(f"({i},{ra:.6f},{de:.6f},{r:.4f})"
                          for i, ra, de, r in zip(b.index, b.ra, b.dec, b.r_arcmin))
        df = _skyserver(HOST_SQL.format(values=values, rmax=C.HOST_MAX_RMAG))
        df.to_csv(path, index=False)
        if k % 20 == 0:
            print(f"    batch {k + 1}/{len(batches)}: {len(df):,} galaxies")
        time.sleep(0.5)


# ---------------------------------------------------------------------------
# 4. The SDSS galaxy population (MPA-JHU), z < 0.1
# ---------------------------------------------------------------------------

BACKGROUND_SQL = """
SELECT i.specObjID, i.ra, i.dec, i.z, i.primtarget, i.targettype,
  e.lgm_tot_p50, e.sfr_tot_p50, e.specsfr_tot_p50, e.bptclass
FROM galSpecInfo i JOIN galSpecExtra e ON e.specObjID = i.specObjID
WHERE i.specObjID >= {lo} AND i.specObjID < {hi}
  AND i.z BETWEEN 0.005 AND {zmax} AND i.reliable = 1
"""


def fetch_sdss_background() -> None:
    plates = range(C.GALSPEC_PLATE_MIN, C.GALSPEC_PLATE_MAX + 1, C.GALSPEC_PLATE_CHUNK)
    print(f"  SDSS background (MPA-JHU, z < {C.HOST_Z_MAX}): {len(plates)} plate ranges")
    for p0 in plates:
        p1 = min(p0 + C.GALSPEC_PLATE_CHUNK, C.GALSPEC_PLATE_MAX + 1)
        path = C.SDSS_BACKGROUND_DIR / f"plates_{p0:04d}_{p1 - 1:04d}.csv"
        if path.exists():
            continue
        # DR8+ specObjID = plate * 2^50 + ..., so a plate range is an ID range.
        df = _skyserver(BACKGROUND_SQL.format(lo=p0 << 50, hi=p1 << 50, zmax=C.HOST_Z_MAX))
        df.to_csv(path, index=False)
        print(f"    plates {p0}-{p1 - 1}: {len(df):,} galaxies")
        time.sleep(0.5)


# ---------------------------------------------------------------------------
# 5-7. Cosmicflows-4, the dust map, paper sources
# ---------------------------------------------------------------------------

def _download(url: str, path, *, headers=None, label: str = "") -> None:
    if path.exists():
        print(f"  [cached] {path.name}")
        return
    print(f"  {label or path.name} ...")
    r = http_get(url, stream=True)
    tmp = path.with_suffix(path.suffix + ".part")
    with open(tmp, "wb") as fh:
        for chunk in r.raw.stream(1 << 20, decode_content=False):
            fh.write(chunk)
    tmp.rename(path)
    _write_meta(path, url=url, bytes=path.stat().st_size)
    print(f"    {path.stat().st_size / 1e6:.1f} MB")


def fetch_cf4() -> None:
    for name in C.CF4_FILES:
        _download(C.CF4_BASE_URL + name, C.CF4_DIR / name, label=f"Cosmicflows-4 {name}")


def fetch_dust() -> None:
    _download(C.SFD_URL, C.DUST_DIR / "lambda_sfd_ebv.fits", label="SFD E(B-V) HEALPix map")


def fetch_papers() -> None:
    for key, arxiv_id in C.ARXIV_SOURCES.items():
        d = C.PAPERS_DIR / key
        d.mkdir(exist_ok=True)
        tarball = d / f"{arxiv_id}.tar.gz"
        _download(f"https://arxiv.org/src/{arxiv_id}", tarball, label=f"arXiv:{arxiv_id} ({key})")
        if not any(d.glob("*.tex")):
            with tarfile.open(tarball) as tf:
                members = [m for m in tf.getmembers() if m.name.endswith(".tex")]
                tf.extractall(d, members=members, filter="data")


# ---------------------------------------------------------------------------
# 8. ZTF exposure history (optional)
# ---------------------------------------------------------------------------

#: Four central readout quadrants: every exposure of a field appears in at
#: least one of them, so (field, filter, time) can be deduplicated locally.
EXPOSURE_SQL = """
SELECT field, fid, obsjd, maglimit, rcid FROM ztf.ztf_current_meta_sci
WHERE obsjd >= {lo} AND obsjd < {hi} AND programid = 1 AND rcid IN (21, 26, 37, 42)
"""


def fetch_ztf_exposures() -> bool:
    """Public-survey exposures month by month. Returns False if IRSA is down."""
    _download(C.ZTF_FIELDS_URL, C.ZTF_EXPOSURE_DIR / "ZTF_Fields.txt", label="ZTF field grid")
    months = pd.date_range(C.STAT_PEAK_START, C.STAT_PEAK_END, freq="MS")
    # Coverage cuts look up to 16.5 days before and 28.5 days after a peak.
    months = pd.date_range(months[0] - pd.DateOffset(months=1),
                           months[-1] + pd.DateOffset(months=2), freq="MS")
    from astropy.time import Time
    print(f"  ZTF exposure history (IRSA TAP): {len(months)} months")
    for m in months:
        path = C.ZTF_EXPOSURE_DIR / f"exposures_{m:%Y_%m}.csv"
        if path.exists():
            continue
        lo = Time(m.to_pydatetime()).jd
        hi = Time((m + pd.DateOffset(months=1)).to_pydatetime()).jd
        sql = " ".join(EXPOSURE_SQL.format(lo=lo, hi=hi).split())
        try:
            r = http_get(C.IRSA_TAP_URL, params={"QUERY": sql, "FORMAT": "csv"}, tries=3)
        except requests.RequestException as exc:
            print(f"    IRSA TAP unavailable ({exc}); rates will use Perley et al.'s coverage")
            return False
        path.write_text(r.text)
        print(f"    {m:%Y-%m}: {r.text.count(chr(10)) - 1:,} quadrant-exposures")
    return True


# ---------------------------------------------------------------------------

STEPS = ("bts", "alerce", "hosts", "background", "cf4", "dust", "papers", "exposures")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--only", nargs="+", choices=STEPS, help="run only these steps")
    ap.add_argument("--refresh-bts", action="store_true",
                    help="take a fresh snapshot of the BTS catalogue (changes every few hours)")
    args = ap.parse_args(argv)
    steps = args.only or STEPS

    if "bts" in steps or not BTS_ALL.exists():
        print("1. BTS catalogue")
        fetch_bts(refresh=args.refresh_bts)
    bts = read_bts_csv()
    if "alerce" in steps:
        print("2. ALeRCE light curves")
        fetch_alerce(bts)
    if "hosts" in steps:
        print("3. SDSS galaxies around the supernovae")
        fetch_sdss_hosts(bts)
    if "background" in steps:
        print("4. SDSS galaxy population")
        fetch_sdss_background()
    if "cf4" in steps:
        print("5. Cosmicflows-4")
        fetch_cf4()
    if "dust" in steps:
        print("6. Dust map")
        fetch_dust()
    if "papers" in steps:
        print("7. Paper sources")
        fetch_papers()
    if "exposures" in steps:
        print("8. ZTF exposure history (optional)")
        fetch_ztf_exposures()
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
