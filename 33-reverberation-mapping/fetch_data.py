"""Download and cache every input this project uses.

    python fetch_data.py              # everything (~15 MB, ~10 minutes)
    python fetch_data.py --no-sdssrm  # skip the ~1,100 SDSS-RM light-curve files

Files land in ``data/`` and are never re-downloaded once present, so after one
run the whole analysis works offline. Each source is listed with what it is
and who made it; see the README for citations.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import config as C

#: CDS answers browser-like User-Agents with an anti-bot challenge page (with
#: HTTP status 200), so every request uses a plain tool User-Agent.
UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

#: VizieR / CDS copies of the journals' machine-readable tables.
VIZIER = {
    # Pei et al. 2017, ApJ 837, 131 (AGN STORM V): 2014 ground-based spectra of
    # NGC 5548 from five telescopes: 5100 A continuum, H-beta, He II 4686.
    "J/ApJ/837/131": ["ReadMe", "table2.dat"],
    # Fausnaugh et al. 2016, ApJ 821, 56 (AGN STORM III): 2014 optical
    # continuum light curves in nine bands (table3) and the HST/COS continuum
    # at 1158, 1367, 1479 and 1746 A (table4).
    "J/ApJ/821/56": ["ReadMe", "table3.dat", "table4.dat"],
    # De Rosa et al. 2015, ApJ 806, 128 (AGN STORM I): HST/COS 1367 A continuum
    # and the Ly-alpha, N V, Si IV, C IV and He II 1640 line light curves.
    "J/ApJ/806/128": ["ReadMe", "table2.dat"],
    # Peterson et al. 2013, ApJ 779, 109: NGC 5548 from 1988 to 2008 (AGN Watch
    # and later campaigns), host-subtracted 5100 A AGN continuum (table5) and
    # broad H-beta (table6), recalibrated for [O III] variability.
    "J/ApJ/779/109": ["ReadMe", "table5.dat", "table6.dat"],
}

#: The AGN Black Hole Mass Database (Bentz & Katz 2015, PASP 127, 67): every
#: published reverberation lag, line width and host-corrected luminosity, one
#: page per AGN. The index lists the pages.
AGNMASS_URL = "http://www.astro.gsu.edu/AGNmass/"

#: SDSS Reverberation Mapping final release (Shen et al. 2024, ApJS 272, 26):
#: a summary table of 849 quasars and per-quasar light curves.
SDSSRM_URL = "https://ariel.astro.illinois.edu/sdssrm/final_result/"
SDSSRM_LINES = {"HA": "ha", "HB": "hb", "MG2": "mg2", "C4": "c4"}

#: van den Bosch (2016, ApJ 831, 134) compilation of 294 dynamically measured
#: black hole masses with stellar velocity dispersions (the M-sigma sample).
VDB16_URL = ("https://raw.githubusercontent.com/remco-space/Black-Hole-Mass-compilation/"
             "master/BHcompilation.fits")

#: arXiv source bundles we read tables from.
ARXIV = {
    # Bentz et al. 2013, ApJ 767, 149: the 71 H-beta lags and host-corrected
    # luminosities behind the R-L relation (their Table 14), to validate our
    # regression against theirs.
    "1303.1742": "bentz2013",
    # Woo et al. 2015, ApJ 801, 38: virial products and stellar velocity
    # dispersions of 29 reverberation-mapped AGN (their Table A1).
    "1412.7225": "woo2015",
    # Shen et al. 2024, ApJS 272, 26: virial products and dynamical-model black
    # hole masses of 30 AGN (their Table 3), the virial factor calibration.
    "2305.01014": "shen2024",
    # Goad et al. 2016, ApJ 824, 11 (AGN STORM IV): the BLR holiday; Figure 1
    # (vector PostScript) colour-codes the epochs of the anomaly.
    "1603.08741": "goad2016",
}

# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------


def download(url: str, dest: Path, *, timeout: float = 120.0, retries: int = 4,
             quiet: bool = False, allow_html: bool = False) -> Path:
    """Fetch ``url`` to ``dest`` (atomically), unless it already exists.

    Apart from the AGN mass database, none of our inputs is a web page, so an
    HTML response (an error or a bot-check page served with status 200) is
    treated as a failure.
    """
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_suffix(dest.suffix + ".part")
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA_TOOL})
            with urllib.request.urlopen(req, timeout=timeout) as resp, open(part, "wb") as fh:
                while chunk := resp.read(1 << 16):
                    fh.write(chunk)
            head = part.read_bytes()[:512].lstrip().lower()
            if not allow_html and (head.startswith(b"<!doctype html") or head.startswith(b"<html")):
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


def fetch_vizier() -> None:
    print("VizieR catalogues (AGN STORM 2014, AGN Watch 1988-2008)")
    for cat, files in VIZIER.items():
        folder = C.DATA_DIR / "vizier" / cat.replace("/", "_")
        for f in files:
            download(f"https://cdsarc.cds.unistra.fr/ftp/{cat}/{f}", folder / f)


def fetch_agnmass() -> None:
    """The database index, then one details page per AGN (polite: one at a time)."""
    print("AGN Black Hole Mass Database")
    folder = C.DATA_DIR / "agnmass"
    index = download(AGNMASS_URL, folder / "index.html", allow_html=True)
    html = index.read_text(encoding="latin-1")
    pages = sorted(set(int(v) for v in re.findall(r"details\.php\?varname=(\d+)", html)))
    todo = [v for v in pages if not (folder / f"details_{v:03d}.html").exists()]
    print(f"  {len(pages)} AGN pages, {len(todo)} to download")
    for v in todo:
        download(f"{AGNMASS_URL}details.php?varname={v}", folder / f"details_{v:03d}.html",
                 allow_html=True, quiet=True)
        time.sleep(0.3)


def sdssrm_targets(summary: Path) -> dict[int, list[str]]:
    """RMIDs whose light curves we need, with the line files for each.

    Every quasar with a detected lag in any line, plus every quasar whose
    H-beta lag was attempted (detected or not), for the false-positive test.
    """
    import numpy as np
    from astropy.io import fits

    d = fits.getdata(summary, 1)
    out: dict[int, list[str]] = {}
    for key, short in SDSSRM_LINES.items():
        sel = d[f"{key}_LAG_DET"] > 0
        if key == "HB":
            sel |= d["HB_LAG_DONE"] == 1
        for rmid in d["RMID"][sel]:
            out.setdefault(int(rmid), []).append(short)
    return {k: sorted(set(v)) for k, v in sorted(out.items())} if out else {}


def fetch_sdssrm(light_curves: bool = True) -> None:
    print("SDSS-RM final release")
    folder = C.DATA_DIR / "sdssrm"
    summary = download(SDSSRM_URL + "summary.fits", folder / "summary.fits")
    download(SDSSRM_URL + "update_notes", folder / "update_notes.txt")
    if not light_curves:
        return
    targets = sdssrm_targets(summary)
    files: list[tuple[int, str]] = []
    for rmid, lines in targets.items():
        files.append((rmid, "cont"))
        files.extend((rmid, ln) for ln in lines)
    todo = [(r, f) for r, f in files if not (folder / f"rm{r:03d}" / f"{f}_lc.txt").exists()]
    print(f"  {len(targets)} quasars, {len(files)} light curves, {len(todo)} to download")
    failed = []
    for i, (rmid, f) in enumerate(todo, 1):
        try:
            download(f"{SDSSRM_URL}rm{rmid:03d}/{f}_lc.txt", folder / f"rm{rmid:03d}" / f"{f}_lc.txt",
                     quiet=True, timeout=60.0)
        except RuntimeError as exc:
            failed.append((rmid, f))
            print(f"  rm{rmid:03d} {f}: {exc}")
        if i % 100 == 0:
            print(f"  {i}/{len(todo)}")
        time.sleep(0.1)
    if failed:
        print(f"  {len(failed)} light curves failed; re-run to retry")


def fetch_arxiv() -> None:
    print("arXiv source bundles")
    for arxiv_id, label in ARXIV.items():
        download(f"https://arxiv.org/e-print/{arxiv_id}", C.DATA_DIR / "arxiv" / f"{label}_{arxiv_id}.tar.gz")
        time.sleep(1.0)


def fetch_vdb16() -> None:
    print("van den Bosch (2016) black hole mass compilation")
    download(VDB16_URL, C.DATA_DIR / "vdb16" / "BHcompilation.fits")


def fetch_all(sdssrm_light_curves: bool = True) -> None:
    fetch_vizier()
    fetch_arxiv()
    fetch_vdb16()
    fetch_agnmass()
    fetch_sdssrm(sdssrm_light_curves)
    manifest = sorted(str(p.relative_to(C.DATA_DIR)) for p in C.DATA_DIR.rglob("*")
                      if p.is_file() and "cache" not in p.parts and not p.name.endswith(".part"))
    (C.DATA_DIR / "manifest.json").write_text(json.dumps({"files": len(manifest),
                                                          "retrieved": time.strftime("%Y-%m-%d")}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-sdssrm", action="store_true", help="skip the SDSS-RM light curves")
    args = ap.parse_args()
    fetch_all(sdssrm_light_curves=not args.no_sdssrm)
    sys.exit(0)
