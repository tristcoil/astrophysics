"""Download and cache every input this project uses.

    python fetch_data.py            # everything (the ZTF light curves take ~20 min)
    python fetch_data.py --no-ztf   # skip the ~700 ZTF light curves

Files land in ``data/`` and are never re-downloaded once present, so after one
run the whole analysis works offline. Each source is listed with what it is
and who made it; see the README for citations.
"""

from __future__ import annotations

import argparse
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

#: CDS answers browser-like User-Agents with an anti-bot challenge page, while
#: the AAS file server wants one. So: a plain tool UA by default, a browser UA
#: only for the hosts that need it.
UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"
UA_BROWSER = "Mozilla/5.0 (X11; Linux x86_64; rv:128.0) Gecko/20100101 Firefox/128.0"
BROWSER_HOSTS = ("content.cld.iop.org", "iopscience.iop.org")

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

#: Open Supernova Catalog (Guillochon et al. 2017, ApJ 835, 64): one JSON file
#: per supernova with every published photometric point the catalog collected.
OSC_FILES = {
    "SN2011fe": "https://raw.githubusercontent.com/astrocatalogs/sne-2010-2014/master/SN2011fe.json",
    "SN2014J": "https://raw.githubusercontent.com/astrocatalogs/sne-2010-2014/master/SN2014J.json",
    "SN1987A": "https://raw.githubusercontent.com/astrocatalogs/sne-pre-1990/master/SN1987A.json",
}

#: AAS journals' "data behind the figure" files for SN 2023ixf, which is newer
#: than the Open Supernova Catalog.
IOP_DBF = {
    # Hsu et al. 2025, ApJ 990, 148: Las Cumbres Observatory UBVgriz, days 1-363.
    "hsu2025_fig1.tar.gz": "https://content.cld.iop.org/journals/0004-637X/990/2/148/revision1/apjadf222f1.tar.gz",
    # Singh et al. 2024, ApJ 975, 132: Swift UV, optical and JHKs, days 1-150.
    "singh2024_fig2.tar.gz": "https://content.cld.iop.org/journals/0004-637X/975/1/132/revision2/apjad7955f2.tar.gz",
}

#: arXiv source bundles: tables and vector figures we read numbers from.
ARXIV = {
    # Matheson et al. 2012, ApJ 754, 19: JHKs light curves of SN 2011fe (Table 1).
    "1205.3828": "matheson2012",
    # Seitenzahl et al. 2014, ApJ 792, 10: Fig. 2 holds the late UVOIR light
    # curve of SN 1987A (Suntzeff et al. 1992, 1997) as vector graphics.
    "1408.5986": "seitenzahl2014",
    # Diehl et al. 2015, A&A 574, A72: INTEGRAL/SPI 56Co line fluxes of SN
    # 2014J in four epochs (table) with error bars (vector figures).
    "1409.5477": "diehl2015",
}

#: VizieR / CDS catalogues.
VIZIER = {
    # Amanullah et al. 2014, ApJ 788, L21: SN 2014J from HST (UV), NOT and the
    # Mount Abu infrared telescope (JHK).
    "J/ApJ/788/L21": ["ReadMe", "table1.dat"],
    # Shappee et al. 2017, ApJ 841, 48: late-time LBT and HST photometry of
    # SN 2011fe, including near-infrared.
    "J/ApJ/841/48": ["ReadMe", "table1.dat"],
    # Pereira et al. 2013, A&A 554, A27: Nearby Supernova Factory
    # spectrophotometry of SN 2011fe (3300-9700 A), 32 epochs.
    "J/A+A/554/A27": ["ReadMe", "list.dat"],
}

#: The published UVOIR bolometric light curve of SN 1987A (Suntzeff et al.
#: 1991, AJ 102, 1118, Tables 5-6; Bouchet et al. 1991, A&A 245, 490), as
#: transcribed by David Jeffery (UNLV) from the papers' tables.
SN1987A_UVOIR_URL = ("https://www.physics.unlv.edu/~jeffery/astro/supernovae/spectra/"
                     "d1980/sn1987a_II/lc_bol2.dat")

#: Filter profiles from the SVO Filter Profile Service (Rodrigo et al. 2012,
#: 2020): pivot wavelengths and Vega zero points for turning magnitudes into
#: fluxes, and transmission curves for synthetic photometry of spectra.
SVO_FILTERS = [
    "Swift/UVOT.UVW2", "Swift/UVOT.UVM2", "Swift/UVOT.UVW1", "Swift/UVOT.U", "Swift/UVOT.B",
    "Swift/UVOT.V",
    "Generic/Bessell.U", "Generic/Bessell.B", "Generic/Bessell.V", "Generic/Bessell.R",
    "Generic/Bessell.I",
    "SLOAN/SDSS.u", "SLOAN/SDSS.g", "SLOAN/SDSS.r", "SLOAN/SDSS.i", "SLOAN/SDSS.z",
    "2MASS/2MASS.J", "2MASS/2MASS.H", "2MASS/2MASS.Ks",
    "Palomar/ZTF.g", "Palomar/ZTF.r",
]

#: ZTF Bright Transient Survey (Fremling et al. 2020; Perley et al. 2020):
#: every transient brighter than 18.5 mag, spectroscopically classified.
BTS_URL = "https://sites.astro.caltech.edu/ztf/bts/explorer.php?format=csv"

#: ALeRCE broker (Forster et al. 2021): ZTF alert light curves per object.
ALERCE_URL = "https://api.alerce.online/ztf/v1/objects/{oid}/lightcurve"

#: ZTF names of the nearby supernovae (for their own ZTF light curves).
ZTF_NAMES = {"SN2023ixf": "ZTF23aaklqou"}

#: Coarse cut on the BTS catalog for the Arnett's-rule sample; the finer
#: quality cuts happen in ztf.py once the light curves are in hand.
ZTF_TYPES = ("SN Ia", "SN Ia-91T", "SN Ia-91bg")
ZTF_Z_RANGE = (0.02, 0.04)

# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------


def download(url: str, dest: Path, *, timeout: float = 120.0, retries: int = 4,
             quiet: bool = False) -> Path:
    """Fetch ``url`` to ``dest`` (atomically), unless it already exists.

    None of our inputs is a web page, so an HTML response (an error or a
    bot-check page served with status 200) is treated as a failure.
    """
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    ua = UA_BROWSER if any(h in url for h in BROWSER_HOSTS) else UA_TOOL
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": ua})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(part, "wb") as fh:
                while chunk := resp.read(1 << 16):
                    fh.write(chunk)
            head = part.read_bytes()[:512].lstrip().lower()
            if head.startswith(b"<!doctype html") or head.startswith(b"<html"):
                raise urllib.error.URLError(f"got an HTML page instead of data from {url}")
            part.replace(dest)
            if not quiet:
                print(f"  fetched {dest.relative_to(C.DATA_DIR)} ({dest.stat().st_size / 1e6:.2f} MB)")
            return dest
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            time.sleep(2.0 * (attempt + 1))
    part.unlink(missing_ok=True)
    raise RuntimeError(f"could not download {url}: {last}")


def fetch_osc() -> None:
    print("Open Supernova Catalog")
    for name, url in OSC_FILES.items():
        download(url, C.DATA_DIR / "osc" / f"{name}.json")


def fetch_iop() -> None:
    print("AAS data-behind-the-figure files (SN 2023ixf)")
    for fname, url in IOP_DBF.items():
        dest = download(url, C.DATA_DIR / "sn2023ixf" / fname)
        with tarfile.open(dest) as tar:
            for member in tar.getmembers():
                if member.isfile() and member.name.endswith(".txt"):
                    out = dest.parent / (fname.split("_")[0] + "_" + Path(member.name).name)
                    if not out.exists():
                        out.write_bytes(tar.extractfile(member).read())


def fetch_arxiv() -> None:
    print("arXiv source bundles")
    for arxiv_id, label in ARXIV.items():
        download(f"https://arxiv.org/e-print/{arxiv_id}", C.DATA_DIR / "arxiv" / f"{label}_{arxiv_id}.tar.gz")


def fetch_vizier() -> None:
    print("VizieR catalogues")
    for cat, files in VIZIER.items():
        folder = C.DATA_DIR / "vizier" / cat.replace("/", "_")
        for f in files:
            download(f"https://cdsarc.cds.unistra.fr/ftp/{cat}/{f}", folder / f)
    # Pereira et al. (2013): the individual spectra listed in list.dat.
    folder = C.DATA_DIR / "vizier" / "J_A+A_554_A27"
    for line in (folder / "list.dat").read_text().splitlines():
        fname = line.split()[0] if line.strip() else ""
        if fname.startswith("11fe") and fname.endswith(".fit"):
            download(f"https://cdsarc.cds.unistra.fr/ftp/J/A+A/554/A27/sp/{fname}",
                     folder / "sp" / fname, quiet=True)
    print("  Pereira et al. 2013 spectra present:", len(list((folder / "sp").glob("*.fit"))))


def fetch_sn1987a() -> None:
    print("SN 1987A published UVOIR light curve")
    download(SN1987A_UVOIR_URL, C.DATA_DIR / "sn1987a" / "jeffery_lc_bol2.dat")


def fetch_filters() -> None:
    print("SVO filter profiles")
    for fid in SVO_FILTERS:
        download(f"http://svo2.cab.inta-csic.es/theory/fps/fps.php?ID={fid}",
                 C.DATA_DIR / "filters" / (fid.replace("/", "_") + ".xml"))


def fetch_bts() -> Path:
    print("ZTF Bright Transient Survey catalogue")
    return download(BTS_URL, C.DATA_DIR / "ztf" / "bts.csv")


def bts_candidates(path: Path) -> list[dict]:
    """BTS rows passing the coarse cut (type, redshift, measured rise and fade)."""
    import csv

    rows = list(csv.DictReader(io.StringIO(path.read_text().lstrip())))
    out = []
    for r in rows:
        if r["type"] not in ZTF_TYPES:
            continue
        try:
            z = float(r["redshift"])
            float(r["rise"]), float(r["fade"])
        except ValueError:
            continue
        if ZTF_Z_RANGE[0] <= z <= ZTF_Z_RANGE[1]:
            out.append(r)
    return out


def fetch_alerce(oid: str) -> Path:
    return download(ALERCE_URL.format(oid=oid), C.DATA_DIR / "ztf" / "lightcurves" / f"{oid}.json",
                    timeout=60.0, quiet=True)


def fetch_ztf(limit: int | None = None) -> None:
    bts = fetch_bts()
    cands = bts_candidates(bts)
    if limit:
        cands = cands[:limit]
    oids = [r["ZTFID"] for r in cands] + list(ZTF_NAMES.values())
    todo = [o for o in oids if not (C.DATA_DIR / "ztf" / "lightcurves" / f"{o}.json").exists()]
    print(f"ALeRCE light curves: {len(oids)} objects, {len(todo)} to download")
    failed = []
    for i, oid in enumerate(todo, 1):
        try:
            fetch_alerce(oid)
        except RuntimeError as exc:
            failed.append(oid)
            print(f"  {oid}: {exc}")
        if i % 50 == 0:
            print(f"  {i}/{len(todo)}")
        time.sleep(0.2)  # be gentle with the public broker
    if failed:
        print(f"  {len(failed)} light curves failed; re-run to retry")


def fetch_all(ztf: bool = True) -> None:
    fetch_osc()
    fetch_iop()
    fetch_arxiv()
    fetch_vizier()
    fetch_sn1987a()
    fetch_filters()
    if ztf:
        fetch_ztf()
    manifest = sorted(str(p.relative_to(C.DATA_DIR)) for p in C.DATA_DIR.rglob("*")
                      if p.is_file() and "cache" not in p.parts and not p.name.endswith(".part"))
    (C.DATA_DIR / "manifest.json").write_text(json.dumps({"files": len(manifest)}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-ztf", action="store_true", help="skip the ZTF light curves")
    args = ap.parse_args()
    fetch_all(ztf=not args.no_ztf)
    sys.exit(0)
