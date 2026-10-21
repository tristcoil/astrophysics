"""Download and cache everything the analysis needs.

All data are free and public:

1. **Catalogues from CDS** (Strasbourg): Yu et al. (2018) global oscillation
   parameters, masses and radii of 16,094 Kepler red giants; APOKASC-3
   (Pinsonneault et al. 2025), seismic and APOGEE spectroscopic properties of
   15,808 evolved stars; Vrard et al. (2016) gravity-mode period spacings of
   6,111 red giants; the Kepler LEGACY sample of main-sequence stars (Lund et
   al. 2017), for 16 Cyg A and B.
2. **Kepler light curves** from the Mikulski Archive for Space Telescopes
   (MAST): four years of 30-minute photometry for each red giant in the
   sample, and one year of 1-minute photometry of 16 Cyg A and B.
3. **The Sun as a star**: 25 years of 1-minute photometry of the whole Sun by
   the VIRGO sun photometers on ESA/NASA's SOHO, from the ESA SOHO Science
   Archive (green channel, 500 nm).
4. **Gaia DR3** astrometry and photometry of the sample from the ESA Gaia
   archive, and B1's 100-parsec Gaia sample (hard-linked from
   ``02-hr-diagram`` when present, so the HR diagram is literally B1's).

Everything is cached under ``data/``; nothing that is present is fetched
again, so once the cache is warm the analysis runs offline.

    python fetch_data.py            # catalogues, sample, light curves, Gaia
"""

from __future__ import annotations

import json
import os
import re
import shutil
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import config as C

# ---------------------------------------------------------------------------
# Download helper
# ---------------------------------------------------------------------------


class DownloadError(RuntimeError):
    pass


def download(url: str, dest: Path, *, retries: int = 5, timeout: float = 120.0,
             expect_binary: bool = False, quiet: bool = True) -> Path:
    """Fetch ``url`` to ``dest`` atomically (``.part`` file, then rename).

    Retries on network errors with a growing pause. If the server sends an
    HTML page where a data file was expected (CDS's bot challenge comes back
    with HTTP 200), the file is rejected rather than cached.
    """
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": C.USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r, open(part, "wb") as f:
                shutil.copyfileobj(r, f, length=1 << 20)
            head = part.read_bytes()[:512].lstrip().lower()
            if head.startswith(b"<!doctype html") or head.startswith(b"<html"):
                if not dest.suffix.lower() in (".html", ".htm"):
                    raise DownloadError(f"got an HTML page instead of data from {url}")
            if expect_binary and head[:6] != b"simple":
                raise DownloadError(f"not a FITS file: {url}")
            part.replace(dest)
            if not quiet:
                print(f"    {dest.name}: {dest.stat().st_size / 1e6:.2f} MB")
            return dest
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError, DownloadError) as e:
            last = e
            if isinstance(e, urllib.error.HTTPError) and e.code == 404:
                break
            time.sleep(min(60, 3 * 2 ** attempt))
    if part.exists():
        part.unlink()
    raise DownloadError(f"failed to download {url}: {last}")


def _fetch_text(url: str, *, timeout: float = 60.0, retries: int = 5) -> str:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": C.USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as e:
            last = e
            if isinstance(e, urllib.error.HTTPError) and e.code == 404:
                break
            time.sleep(min(60, 3 * 2 ** attempt))
    raise DownloadError(f"failed to fetch {url}: {last}")


# ---------------------------------------------------------------------------
# 1. CDS catalogues
# ---------------------------------------------------------------------------

#: One folder per catalogue, files under their original names: astropy's CDS
#: reader finds a table in the ReadMe by its file name.
CATALOGUES = {
    "yu2018": ("J/ApJS/236/42", ["table1.dat", "table2.dat"]),
    "apokasc3": ("J/ApJS/276/69", ["table4.dat", "table5.dat", "table9.dat"]),
    "vrard2016": ("J/A+A/588/A87", ["table2.dat"]),
    "lund2017": ("J/ApJ/835/172", ["table1.dat", "table6.dat"]),
}
CDS_FTP = "https://cdsarc.cds.unistra.fr/ftp/"


def fetch_catalogues(*, quiet: bool = False) -> None:
    """CDS tables, stored uncompressed under their original names.

    CDS serves large tables only gzipped (``table4.dat.gz``) and small ones
    plain, so both spellings are tried.
    """
    import gzip

    for name, (cat, files) in CATALOGUES.items():
        folder = C.CATALOG_DIR / name
        for fname in ["ReadMe"] + files:
            dest = folder / fname
            if dest.exists():
                continue
            if not quiet:
                print(f"  {cat}/{fname}")
            try:
                download(CDS_FTP + cat + "/" + fname, dest, quiet=quiet, retries=2)
            except DownloadError:
                gz = download(CDS_FTP + cat + "/" + fname + ".gz",
                              folder / (fname + ".gz"), quiet=quiet)
                part = dest.with_name(dest.name + ".part")
                with gzip.open(gz, "rb") as fi, open(part, "wb") as fo:
                    shutil.copyfileobj(fi, fo)
                part.replace(dest)
                gz.unlink()


# ---------------------------------------------------------------------------
# 2. Kepler light curves
# ---------------------------------------------------------------------------

MAST_KEPLER = "https://archive.stsci.edu/missions/kepler/lightcurves/"


def kepler_dir(kic: int) -> Path:
    return C.KEPLER_DIR / f"{int(kic):09d}"


def _kepler_listing(kic: int) -> list[str]:
    """Long-cadence file names in a star's MAST directory (cached as JSON)."""
    manifest = kepler_dir(kic) / "manifest.json"
    if manifest.exists():
        return json.loads(manifest.read_text())["files"]
    k9 = f"{int(kic):09d}"
    html = _fetch_text(f"{MAST_KEPLER}{k9[:4]}/{k9}/")
    files = sorted(set(re.findall(rf"kplr{k9}-\d+_llc\.fits", html)))
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps({
        "kic": int(kic), "files": files,
        "source": f"{MAST_KEPLER}{k9[:4]}/{k9}/",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=1))
    return files


def fetch_kepler_star(kic: int) -> list[Path]:
    """Every long-cadence quarter of one star (Q0-Q17 where observed)."""
    k9 = f"{int(kic):09d}"
    out = []
    for fname in _kepler_listing(kic):
        out.append(download(f"{MAST_KEPLER}{k9[:4]}/{k9}/{fname}", kepler_dir(kic) / fname,
                            expect_binary=True))
    return out


def fetch_kepler_sample(kics, *, workers: int = 8, quiet: bool = False) -> dict:
    """Download the long-cadence light curves of many stars in parallel."""
    kics = [int(k) for k in kics]
    done, failed = {}, {}
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(fetch_kepler_star, k): k for k in kics}
        for i, fut in enumerate(as_completed(futures), 1):
            k = futures[fut]
            try:
                done[k] = len(fut.result())
            except Exception as e:  # keep going; report at the end
                failed[k] = str(e)
            if not quiet and (i % 10 == 0 or i == len(kics)):
                print(f"    {i}/{len(kics)} stars ({time.time() - t0:.0f} s)")
    if failed and not quiet:
        print(f"  {len(failed)} stars failed: {failed}")
    return {"files": done, "failed": failed}


#: 16 Cyg A and B are too bright (Kp ~ 6) for standard apertures; Kepler
#: observed them with custom pixel masks, archived under their own "KIC-like"
#: target numbers, one per month. One year of short cadence (Quarters 7-10,
#: the data set of Metcalfe et al. 2012) resolves their modes amply.
SIXTEEN_CYG = {
    "16 Cyg A": {"kic": 12069424, "quarters": (7, 8, 9, 10)},
    "16 Cyg B": {"kic": 12069449, "quarters": (7, 8, 9, 10)},
}


def _sc_manifest(name: str, quarters) -> list[str]:
    """MAST URLs of a custom-aperture star's short-cadence files (cached)."""
    slug = name.replace(" ", "_").lower()
    manifest = C.KEPLER_DIR / f"{slug}_sc_manifest.json"
    if manifest.exists():
        return json.loads(manifest.read_text())["urls"]
    import warnings

    import lightkurve as lk

    warnings.filterwarnings("ignore")
    res = lk.search_lightcurve(name, mission="Kepler", cadence="short", author="Kepler")
    urls = []
    for row in res.table:
        q = int(re.search(r"Quarter (\d+)", row["mission"]).group(1))
        if q in quarters:
            uri = str(row["dataURI"])  # mast:KEPLER/url/missions/kepler/...
            urls.append("https://archive.stsci.edu/" + uri.split("/url/", 1)[1])
    urls = sorted(set(urls))
    manifest.write_text(json.dumps({
        "target": name, "quarters": list(quarters), "urls": urls,
        "source": "MAST, via lightkurve.search_lightcurve(mission='Kepler', cadence='short')",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=1))
    return urls


def fetch_16cyg(*, quiet: bool = False) -> dict:
    out = {}
    for name, info in SIXTEEN_CYG.items():
        urls = _sc_manifest(name, info["quarters"])
        folder = C.KEPLER_DIR / name.replace(" ", "_").lower()
        paths = []
        for u in urls:
            paths.append(download(u, folder / u.rsplit("/", 1)[1], expect_binary=True))
        out[name] = paths
        if not quiet:
            print(f"  {name}: {len(paths)} monthly short-cadence files")
    return out


# ---------------------------------------------------------------------------
# 3. The Sun: SOHO/VIRGO sun photometers
# ---------------------------------------------------------------------------

#: Mission-long level-2 time series of the VIRGO SPM green channel
#: (500 nm, 5 nm wide; 1-minute samples, 1996-2021), as linked from ESA's
#: SOHO "mission long files" page (cosmos.esa.int/web/soho/mission-long-files).
VIRGO_GREEN_URL = ("https://ssa.esac.esa.int/ssa-sl-tap/data?retrieval_type=OBSERVATION"
                   "&QUERY=WHERE+observation_oid+=194187590")
VIRGO_GREEN_PATH = C.VIRGO_DIR / "VIRGO-SPM-GREEN-L2-MISSIONLONG.fits"


def fetch_virgo(*, quiet: bool = False) -> Path:
    if not VIRGO_GREEN_PATH.exists() and not quiet:
        print("  VIRGO SPM green channel, 1996-2021 (229 MB)...")
    return download(VIRGO_GREEN_URL, VIRGO_GREEN_PATH, timeout=600, quiet=quiet)


# ---------------------------------------------------------------------------
# 4. Gaia
# ---------------------------------------------------------------------------

GAIA_SAMPLE_PATH = C.GAIA_DIR / "gaia_dr3_sample.fits"
B1_SOURCE = C.PROJECT_DIR.parent / "02-hr-diagram" / "data" / "gaia_dr3_parallax_gt_10mas.fits"
B1_LOCAL = C.GAIA_DIR / "gaia_dr3_parallax_gt_10mas.fits"

GAIA_COLUMNS = [
    "source_id", "ra", "dec", "l", "b", "parallax", "parallax_error", "ruwe",
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "nu_eff_used_in_astrometry", "pseudocolour", "ecl_lat", "astrometric_params_solved",
    "ag_gspphot", "ebpminrp_gspphot", "teff_gspphot",
]


def fetch_gaia_sample(source_ids, *, quiet: bool = False) -> Path:
    """Gaia DR3 rows for the sample, by source_id (one synchronous query)."""
    from astropy.table import Table

    ids = sorted({int(s) for s in source_ids if s and int(s) > 0})
    if GAIA_SAMPLE_PATH.exists():
        have = set(Table.read(GAIA_SAMPLE_PATH)["source_id"].astype(int))
        if set(ids) <= have:
            return GAIA_SAMPLE_PATH
    from astroquery.gaia import Gaia

    query = (f"SELECT {', '.join(GAIA_COLUMNS)} FROM gaiadr3.gaia_source "
             f"WHERE source_id IN ({', '.join(str(i) for i in ids)})")
    if not quiet:
        print(f"  Gaia DR3: {len(ids)} sources")
    last = None
    for attempt in range(4):
        try:
            tab = Gaia.launch_job(query, dump_to_file=False).get_results()
            break
        except Exception as e:  # the archive is sometimes busy
            last = e
            time.sleep(20 * (attempt + 1))
    else:
        raise DownloadError(f"Gaia query failed: {last}")
    tab.meta = {"QUERY": "gaiadr3.gaia_source WHERE source_id IN (<sample>)",
                "SOURCE": "ESA Gaia archive (gea.esac.esa.int), Gaia DR3",
                "RETRIEVD": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    for c in tab.colnames:
        tab.rename_column(c, c.lower())
    tab.write(GAIA_SAMPLE_PATH, overwrite=True)
    return GAIA_SAMPLE_PATH


EXTINCTION_PATH = C.GAIA_DIR / "irsa_dust_sfd.csv"
IRSA_DUST = "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr={ra:.6f}+{dec:.6f}+equ+j2000&regSize=2.0"


def fetch_extinction(gaia_path: Path = GAIA_SAMPLE_PATH, *, quiet: bool = False) -> Path:
    """Total line-of-sight reddening E(B-V) at each star, from IRSA's dust service.

    The Schlafly & Finkbeiner (2011) recalibration of the Schlegel, Finkbeiner
    & Davis (1998) dust map, at each Gaia position (cached as CSV). The giants
    are 0.4-4 kpc away at Galactic latitudes 6-21 deg, mostly beyond the dust
    layer, so the full column is a fair estimate (and an overestimate for the
    nearest few).
    """
    import pandas as pd
    from astropy.table import Table

    if EXTINCTION_PATH.exists():
        return EXTINCTION_PATH
    g = Table.read(gaia_path)
    rows = []
    for sid, ra, dec in zip(g["source_id"], g["ra"], g["dec"]):
        xml = _fetch_text(IRSA_DUST.format(ra=float(ra), dec=float(dec)))
        sf = re.search(r"<refPixelValueSandF>\s*([-\d.]+)", xml)
        sfd = re.search(r"<refPixelValueSFD>\s*([-\d.]+)", xml)
        rows.append({"source_id": int(sid), "ebv_sandf": float(sf.group(1)),
                     "ebv_sfd": float(sfd.group(1))})
    pd.DataFrame(rows).to_csv(EXTINCTION_PATH, index=False)
    if not quiet:
        print(f"  IRSA dust: E(B-V) for {len(rows)} stars")
    return EXTINCTION_PATH


SILSO_LOCAL = C.DATA_DIR / "silso" / "SN_y_tot_V2.0.csv"
SILSO_F1 = C.PROJECT_DIR.parent / "15-solar-cycles" / "data" / "silso" / "SN_y_tot_V2.0.csv"
SILSO_URL = "https://www.sidc.be/SILSO/DATA/SN_y_tot_V2.0.csv"


def fetch_sunspots(*, quiet: bool = False) -> Path:
    """Yearly international sunspot number (WDC-SILSO, Royal Observatory of
    Belgium), hard-linked from F1 (``15-solar-cycles``) when present."""
    if SILSO_LOCAL.exists():
        return SILSO_LOCAL
    SILSO_LOCAL.parent.mkdir(parents=True, exist_ok=True)
    if SILSO_F1.exists():
        try:
            os.link(SILSO_F1, SILSO_LOCAL)
        except OSError:
            shutil.copy2(SILSO_F1, SILSO_LOCAL)
        return SILSO_LOCAL
    return download(SILSO_URL, SILSO_LOCAL, quiet=quiet)


def link_b1_sample(*, quiet: bool = False) -> Path | None:
    """B1's Gaia DR3 sample (every star with parallax > 10 mas).

    Hard-linked from the B1 project when it is present (no copy, no
    download); otherwise re-queried with B1's exact query.
    """
    if B1_LOCAL.exists():
        return B1_LOCAL
    if B1_SOURCE.exists():
        try:
            os.link(B1_SOURCE, B1_LOCAL)
        except OSError:
            shutil.copy2(B1_SOURCE, B1_LOCAL)
        if not quiet:
            print("  linked B1's Gaia sample (02-hr-diagram)")
        return B1_LOCAL
    if not quiet:
        print("  B1's Gaia sample not found: querying parallax > 10 mas (~575k rows, slow)")
    import pyvo

    cols = ("source_id, ra, dec, parallax, parallax_error, parallax_over_error, ruwe, "
            "phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, bp_rp, "
            "phot_bp_rp_excess_factor")
    query = f"SELECT {cols} FROM gaiadr3.gaia_source WHERE parallax > 10"
    job = pyvo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap").submit_job(
        query, maxrec=5_000_000)
    job.run()
    job.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=7200)
    tab = job.fetch_result().to_table()
    tab.write(B1_LOCAL, overwrite=True)
    return B1_LOCAL


# ---------------------------------------------------------------------------


def fetch_all(*, quiet: bool = False) -> None:
    import sample as S

    if not quiet:
        print("Catalogues (CDS)")
    fetch_catalogues(quiet=quiet)
    smp = S.load_or_build_sample(quiet=quiet)
    if not quiet:
        print(f"Kepler long cadence: {len(smp)} red giants (MAST)")
    res = fetch_kepler_sample(smp["kic"], quiet=quiet)
    if res["failed"]:
        raise DownloadError(f"{len(res['failed'])} light curves failed: {res['failed']}")
    if not quiet:
        print("Kepler short cadence: 16 Cyg A and B (MAST)")
    fetch_16cyg(quiet=quiet)
    if not quiet:
        print("The Sun (SOHO/VIRGO)")
    fetch_virgo(quiet=quiet)
    if not quiet:
        print("Gaia DR3")
    fetch_gaia_sample(smp["gaia_dr3"], quiet=quiet)
    fetch_extinction(quiet=quiet)
    link_b1_sample(quiet=quiet)
    fetch_sunspots(quiet=quiet)


if __name__ == "__main__":
    fetch_all(quiet="-q" in sys.argv)
