"""Download and cache every input this project uses (about 20 MB, a minute).

    python fetch_data.py

Nothing is fetched twice: files already in ``data/raw`` are kept, so after one
run the whole analysis works offline. Each download is logged with its URL,
size, SHA-256 and date in ``data/manifest.json``.

Sources (citations in config.py and the README):

* arXiv e-prints (LaTeX sources) of the astrometry papers. Only the .tex files
  holding the tables are kept; the tables are parsed by astrometry.py.
* orbitize!'s example data on GitHub: beta Pic b astrometry (Nielsen+2020).
* The Hipparcos-Gaia Catalog of Accelerations, EDR3 edition (Brandt 2021),
  as bundled with orbitize! (or from CDS if the bundled copy is missing).
* Gaia DR3 (ESA archive TAP): HR 8799 and beta Pic.
* whereistheplanet's HR 8799 orbit posterior (an outside check only).
"""

from __future__ import annotations

import hashlib
import io
import json
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import config as C

#: arXiv and CDS both answer plain tool User-Agents; arXiv asks for a contact.
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


def _get(url: str, *, retries: int = 4, data: bytes | None = None, timeout: int = 180) -> bytes:
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, OSError, TimeoutError) as exc:
            print(f"  attempt {attempt}/{retries} failed for {url}: {exc}", file=sys.stderr)
            time.sleep(3 * attempt)
    raise SystemExit(f"could not download {url}")


def fetch_arxiv() -> None:
    print("arXiv LaTeX sources")
    for src in C.ARXIV.values():
        folder = src.folder
        if all((folder / f).exists() for f in src.files if "*" not in f) and \
                (not any("*" in f for f in src.files) or any(folder.glob("*.tex"))):
            continue
        blob = _get(src.url)
        folder.mkdir(parents=True, exist_ok=True)
        tar_path = folder / "eprint.tar.gz"
        tar_path.write_bytes(blob)
        try:
            with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tar:
                names = tar.getnames()
                for want in src.files:
                    if "*" in want:               # every file matching the pattern
                        for n in names:
                            if Path(n).match(want) and tar.getmember(n).isfile():
                                with tar.extractfile(n) as fin, open(folder / Path(n).name, "wb") as fout:
                                    shutil.copyfileobj(fin, fout)
                        continue
                    hits = [n for n in names if Path(n).name == want]
                    if not hits:
                        raise SystemExit(f"{want} not in the e-print of {src.arxiv}: {names[:20]}")
                    member = tar.getmember(hits[0])
                    with tar.extractfile(member) as fin, open(folder / want, "wb") as fout:
                        shutil.copyfileobj(fin, fout)
        except tarfile.ReadError:
            # single-file submissions are a gzipped .tex, not a tarball
            import gzip
            text = gzip.decompress(blob)
            (folder / src.files[0]).write_bytes(text)
        _log(f"arxiv:{src.arxiv}", src.url, tar_path, src.citation)
        tar_path.unlink()          # keep only the table files (the PDFs are bulky)
        print(f"  {src.key:16s} {src.arxiv}: {', '.join(src.files)}")
        time.sleep(3)              # be polite to arXiv


def fetch_betapic() -> None:
    dest = C.RAW / C.ORBITIZE_BETAPIC
    if dest.exists():
        return
    print("orbitize! example data: beta Pic b astrometry")
    dest.write_bytes(_get(C.ORBITIZE_BETAPIC_URL))
    _log("orbitize:betaPic.csv", C.ORBITIZE_BETAPIC_URL, dest,
         "Nielsen et al. 2020 Table 1 + Snellen et al. 2014 RV, as shipped with orbitize!")


def fetch_hgca() -> None:
    """HR 8799's row of the HGCA (EDR3 edition) from CDS, kept as JSON.

    The copy of HGCA_vEDR3.fits bundled with orbitize! has the ra-dec
    correlation coefficients of the three proper motions set to zero for this
    star; the CDS catalogue gives 0.1087 (Gaia), 0.1291 (Hipparcos-Gaia) and
    0.07 (Hipparcos), as Brandt et al. (2021) used. We take the CDS row.
    """
    dest = C.DERIVED / "hgca_hr8799.json"
    if dest.exists():
        return
    print("Hipparcos-Gaia Catalog of Accelerations (EDR3): HR 8799 row from CDS")
    query = urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
        "QUERY": f'SELECT * FROM "J/ApJS/254/42/catalog" WHERE HIP={C.HIP_ID}'}).encode()
    text = _get(C.VIZIER_TAP, data=query).decode()
    raw = C.RAW / "hgca_hr8799_cds.csv"
    raw.write_text(text)
    import csv
    row = next(csv.DictReader(io.StringIO(text)))
    names = {  # CDS label -> HGCA FITS column name
        "HIP": "hip_id", "Gaia": "Gaia", "RA_ICRS": "RA_ICRS", "DE_ICRS": "DE_ICRS", "plx": "plx", "e_plx": "e_plx",
        "pmRA": "pmra_gaia", "pmDE": "pmdec_gaia", "e_pmRA": "pmra_gaia_error", "e_pmDE": "pmdec_gaia_error",
        "pmRApmDEcor": "pmra_pmdec_gaia", "pmRAhg": "pmra_hg", "pmDEhg": "pmdec_hg", "e_pmRAhg": "pmra_hg_error",
        "e_pmDEhg": "pmdec_hg_error", "pmRApmDEhg": "pmra_pmdec_hg", "pmRAhip": "pmra_hip", "pmDEhip": "pmdec_hip",
        "e_pmRAhip": "pmra_hip_error", "e_pmDEhip": "pmdec_hip_error", "pmRApmDEhip": "pmra_pmdec_hip",
        "EpochRAgaia": "epoch_ra_gaia", "EpochDEgaia": "epoch_dec_gaia", "EpochRAhip": "epoch_ra_hip",
        "EpochDEhip": "epoch_dec_hip", "chi2": "chi2"}
    out = {fits_name: float(row[cds]) for cds, fits_name in names.items()}
    out["_source"] = "CDS J/ApJS/254/42 (Brandt 2021), via TAPVizieR"
    dest.write_text(json.dumps(out, indent=1))
    _log("hgca:HR8799", C.VIZIER_TAP, raw, out["_source"])


GAIA_COLUMNS = ("source_id, ra, dec, parallax, parallax_error, pmra, pmra_error, pmdec, "
                "pmdec_error, pmra_pmdec_corr, ruwe, astrometric_excess_noise, "
                "astrometric_excess_noise_sig, astrometric_params_solved, phot_g_mean_mag, "
                "bp_rp, radial_velocity, radial_velocity_error, ipd_gof_harmonic_amplitude, "
                "ipd_frac_multi_peak, visibility_periods_used")


def _gaia_sync(adql: str) -> str:
    query = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv",
                                    "QUERY": adql}).encode()
    return _get(C.GAIA_TAP, data=query).decode()


def fetch_gaia() -> None:
    dest = C.RAW / "gaia_dr3_hr8799_betapic.csv"
    if dest.exists():
        return
    print("Gaia DR3: HR 8799 and beta Pic")
    # beta Pic by a cone search (it has no DR3 RV and its source_id is not in
    # our sources); the brightest source within 5 arcsec is the star.
    adql = (f"SELECT {GAIA_COLUMNS} FROM gaiadr3.gaia_source WHERE "
            f"source_id = {C.GAIA_DR3_SOURCE_ID} OR "
            "(1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 86.82119875, -51.06651143, 0.0014))"
            " AND phot_g_mean_mag < 6)")
    text = _gaia_sync(adql)
    if text.count("\n") < 3:
        raise SystemExit(f"unexpected Gaia answer:\n{text[:500]}")
    dest.write_text(text)
    _log("gaia:dr3", C.GAIA_TAP, dest, adql)


def fetch_gaia_bright() -> None:
    """Comparison sample for HR 8799's Gaia RUWE: stars of similar brightness
    and colour (G 5.6-6.2, BP-RP 0.1-0.7), where Gaia's calibration is hardest."""
    dest = C.RAW / "gaia_dr3_bright_comparison.csv"
    if dest.exists():
        return
    print("Gaia DR3: bright comparison stars for RUWE")
    adql = ("SELECT source_id, phot_g_mean_mag, bp_rp, ruwe, astrometric_excess_noise, "
            "astrometric_excess_noise_sig, astrometric_params_solved, parallax FROM gaiadr3.gaia_source "
            "WHERE phot_g_mean_mag BETWEEN 5.6 AND 6.2 AND bp_rp BETWEEN 0.1 AND 0.7")
    dest.write_text(_gaia_sync(adql))
    _log("gaia:bright_comparison", C.GAIA_TAP, dest, adql)


def fetch_witp() -> None:
    dest = C.RAW / "witp_post_hr8799.hdf5"
    if dest.exists():
        return
    print("whereistheplanet: current HR 8799 posterior (outside check)")
    blob = _get(C.WITP_URL)
    if blob.startswith(b"version https://git-lfs"):
        raise SystemExit("got a Git LFS pointer instead of the file")
    dest.write_bytes(blob)
    _log("witp:post_hr8799", C.WITP_URL, dest,
         "orbitize! posterior with unpublished GRAVITY data; not used as data")


def ensure_data() -> None:
    fetch_arxiv()
    fetch_betapic()
    fetch_hgca()
    fetch_gaia()
    fetch_gaia_bright()
    fetch_witp()


if __name__ == "__main__":
    ensure_data()
    m = _manifest()
    total = sum((v.get("bytes") or 0) for v in m.values())
    print(f"done: {len(m)} items in data/manifest.json")
