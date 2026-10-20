"""Download and cache everything the analysis needs.

    python fetch_data.py                 # everything (~0.25 GB, ~5-10 minutes)
    python fetch_data.py --only gaia     # one piece: gaia, dr16q, literature

1. **Gaia DR3's quasars** -- the 1,614,173 sources of the Gaia-CRF3, the
   celestial reference frame built from compact extragalactic sources
   (``gaiadr3.agn_cross_id``), with their astrometry from ``gaia_source`` and,
   where Gaia classified them itself, its spectro-photometric redshift from
   ``qso_candidates``. Fetched as 48+ small synchronous queries, one HEALPix
   sky pixel each, and checked against server-side counts: all sources, and
   the five- and six-parameter solutions separately.
2. **SDSS DR16Q redshifts** (Lyke et al. 2020): positions and spectroscopic
   redshifts of 750,414 quasars, to test that the drift is the same at every
   distance. Re-used from project I3 (34-quasars) when its cache exists,
   otherwise fetched from VizieR VII/289 in RA strips.
3. **The frame rotator's sources** (``gaiadr3.frame_rotator_source``): the
   428,034 quasars whose proper motions were used to set Gaia's reference
   frame spinning at zero during the astrometric solution, plus 1,215 that
   were considered and rejected -- to test the l = 1 rotation of the fit.
4. **Quasar look-alikes**: a random 0.5% of all Gaia DR3 sources with
   five-parameter solutions whose parallax and proper motion are both
   consistent with zero within 5 sigma, G > 16, that are *not* in the
   Gaia-CRF3 -- mostly distant stars. They show why the sample must be pure:
   Klioner et al. (2021) found such sources glide ~70 times faster than the
   quasars, in another direction. Ten asynchronous queries (~30 minutes).
5. **Quaia** (Storey-Fisher et al. 2024): the Gaia-unWISE quasar catalogue,
   1.3 million quasars to G = 20.5 with redshifts estimated from Gaia's
   low-resolution spectra and unWISE colours, calibrated on SDSS -- the
   redshifts used by the 2025-26 papers that report a redshift-dependent
   drift. One FITS file from Zenodo (171 MB), checked against its MD5.
6. **Literature**: the Gaia EDR3 acceleration paper (Klioner et al. 2021),
   kept for provenance; the numbers we grade against are in config.REFERENCE.

Every file is recorded in data/manifest.json with its size, SHA-256, row count
and retrieval date.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from astropy.table import Table

import gaia_archive as ga
from config import (
    AGN_COLUMNS, CHUNK_DIR, DR16Q_VIZIER, GAIA_COLUMNS, I3_DR16Q_PARQUET, LIT_DIR, MANIFEST_PATH,
    N_CRF3, N_CRF3_5P, N_CRF3_6P, N_DR16Q, QSOC_COLUMNS, RAW_DIR, VIZIER_ASU_TSV,
)

GAIA_PATH = RAW_DIR / "gaia_crf3.fits"
DR16Q_PATH = RAW_DIR / "dr16q_radecz.parquet"
DR16Q_STRIP_DIR = RAW_DIR / "dr16q_vizier"
#: The archive answers a level-1 sky pixel (~34,000 quasars) too slowly for
#: its synchronous limit, so queries start at level 2 (~8,400), four at a time.
START_LEVEL = 2
WORKERS = 4
#: Non-browser user agents are refused by some CDS mirrors; arXiv accepts this.
USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) astrophysics-atlas-H4"

_FROM = ("FROM gaiadr3.agn_cross_id AS a "
         "JOIN gaiadr3.gaia_source AS g ON g.source_id = a.source_id ")
_WHERE = "WHERE a.source_id >= {lo} AND a.source_id < {hi}"

GAIA_TEMPLATE = (
    "SELECT " + ", ".join([f"g.{c}" for c in GAIA_COLUMNS]
                         + [f"a.{c}" for c in AGN_COLUMNS]
                         + [f"q.{c}" for c in QSOC_COLUMNS]) + " "
    + _FROM
    + "LEFT OUTER JOIN gaiadr3.qso_candidates AS q ON q.source_id = a.source_id "
    + _WHERE
)
#: One server-side count per sky pixel, split by solution type (31 = five
#: parameters, 95 = six parameters: position, parallax, proper motion and, for
#: 95, an astrometric estimate of the colour).
COUNT_TEMPLATE = ("SELECT g.astrometric_params_solved, COUNT(*) AS n " + _FROM + _WHERE
                  + " GROUP BY g.astrometric_params_solved")
EXPECTED_COUNTS = {"all": N_CRF3, "five_parameter": N_CRF3_5P, "six_parameter": N_CRF3_6P}

FRAME_ROTATOR_PATH = RAW_DIR / "frame_rotator_source.fits"
FRAME_ROTATOR_ADQL = ("SELECT source_id, considered_for_reference_frame_spin, used_for_reference_frame_spin, "
                      "considered_for_reference_frame_orientation, used_for_reference_frame_orientation "
                      "FROM gaiadr3.frame_rotator_source")
N_FRAME_ROTATOR_SPIN_USED = 428_034

LOOKALIKE_PATH = RAW_DIR / "lookalikes_halfpct.fits"
#: random_index runs over 0 .. 1,811,709,770: the first 0.5% of it, fetched in
#: 10 slices (anonymous queries at ESA are cut off after 60 s, async included).
LOOKALIKE_N_RANDOM = 9_058_549
LOOKALIKE_SLICES = 10
LOOKALIKE_ADQL = (
    "SELECT g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.pmra, g.pmra_error, "
    "g.pmdec, g.pmdec_error, g.pmra_pmdec_corr, g.phot_g_mean_mag, g.nu_eff_used_in_astrometry, g.ruwe "
    "FROM gaiadr3.gaia_source AS g "
    "LEFT OUTER JOIN gaiadr3.agn_cross_id AS a ON a.source_id = g.source_id "
    "WHERE g.random_index >= {lo} AND g.random_index < {hi} AND g.astrometric_params_solved = 31 "
    "AND g.phot_g_mean_mag > 16 AND a.source_id IS NULL "
    "AND ABS((g.parallax + 0.017) / g.parallax_error) < 5 "
    "AND (POWER(g.pmra / g.pmra_error, 2) + POWER(g.pmdec / g.pmdec_error, 2) "
    "- 2 * g.pmra_pmdec_corr * (g.pmra / g.pmra_error) * (g.pmdec / g.pmdec_error)) "
    "/ (1 - POWER(g.pmra_pmdec_corr, 2)) < 25")

QUAIA_RECORD = "https://zenodo.org/api/records/10403370"
QUAIA_FILE = "quaia_G20.5.fits"
QUAIA_PATH = RAW_DIR / QUAIA_FILE

LITERATURE = {
    "klioner2021.pdf": "https://arxiv.org/pdf/2012.02036",
}


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def http_get(url: str, *, params: dict | None = None, timeout: float = 300,
             attempts: int = 4) -> bytes:
    if params:
        url = url + "?" + urllib.parse.urlencode(params)
    last = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001 -- network trouble of any kind
            last = exc
            time.sleep(4 * attempt)
    raise RuntimeError(f"GET {url[:120]} failed: {last}")


# --------------------------------------------------------------------------
# 1. Gaia DR3 quasars
# --------------------------------------------------------------------------
def fetch_gaia() -> dict:
    counts_path = RAW_DIR / "gaia_crf3_server_counts.json"
    if GAIA_PATH.exists() and counts_path.exists():
        print(f"  {GAIA_PATH.name}: cached")
        return json.loads(GAIA_PATH.with_suffix(".json").read_text())

    print("  Gaia-CRF3 sources, one sky pixel at a time ...", flush=True)
    t0 = time.time()
    table, log = ga.fetch_all_sky("crf3", GAIA_TEMPLATE, start_level=START_LEVEL,
                                  workers=WORKERS)
    print(f"  {len(table):,} rows in {time.time() - t0:.0f} s; verifying with server counts ...",
          flush=True)

    by_type = ga.count_all_sky(COUNT_TEMPLATE, start_level=START_LEVEL, workers=WORKERS)
    counts = {"all": sum(by_type.values()), "five_parameter": by_type.get(31, 0),
              "six_parameter": by_type.get(95, 0)}
    for key in counts:
        print(f"    server COUNT(*) {key:>15}: {counts[key]:,} (published {EXPECTED_COUNTS[key]:,})",
              flush=True)
    local = {
        "all": len(table),
        "five_parameter": int(np.sum(np.asarray(table["astrometric_params_solved"]) == 31)),
        "six_parameter": int(np.sum(np.asarray(table["astrometric_params_solved"]) == 95)),
    }
    for key in counts:
        if not (counts[key] == local[key] == EXPECTED_COUNTS[key]):
            raise RuntimeError(f"{key}: server {counts[key]}, downloaded {local[key]}, "
                               f"published {EXPECTED_COUNTS[key]}")
    ids = np.asarray(table["source_id"])
    if len(np.unique(ids)) != len(ids):
        raise RuntimeError("duplicate source_ids in the stacked download")
    counts_path.write_text(json.dumps({"server": counts, "downloaded": local,
                                       "checked_utc": ga.now_utc()}, indent=2) + "\n")

    table.sort("source_id")
    ga.write_with_provenance(table, GAIA_PATH, {
        "source": "ESA Gaia archive, Gaia DR3 (astrometry identical to EDR3)",
        "endpoint": "https://gea.esac.esa.int/tap-server/tap/sync",
        "adql_template": GAIA_TEMPLATE,
        "chunks": log,
        "server_counts": counts,
        "citation": "Gaia Collaboration, Vallenari et al. 2023, A&A 674, A1; "
                    "Gaia Collaboration, Klioner et al. 2022, A&A 667, A148 (Gaia-CRF3)",
    })
    return json.loads(GAIA_PATH.with_suffix(".json").read_text())


# --------------------------------------------------------------------------
# 2. SDSS DR16Q redshifts
# --------------------------------------------------------------------------
def parse_vizier_tsv(text: str) -> pd.DataFrame:
    """Parse a VizieR ASU tab-separated answer (header, units row, dashes row, data).

    The units row can be blank, so the parser anchors on the row of dashes
    (dropping blank lines and skipping two rows would silently lose a row).
    """
    if "QUERY_STATUS=OVERFLOW" in text or "#Truncated" in text:
        raise RuntimeError("VizieR truncated the answer; use smaller sky strips")
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        return pd.DataFrame()
    header = lines[i].split("\t")
    j = i + 1
    while j < len(lines) and not re.fullmatch(r"[\t ]*-[-\t ]*", lines[j]):
        j += 1
    if j >= len(lines):
        raise RuntimeError("VizieR answer has no header separator")
    body = "\n".join(ln for ln in lines[j + 1:] if ln.strip())
    df = pd.read_csv(io.StringIO(body), sep="\t", header=None, names=header, dtype=str,
                     keep_default_na=False, quoting=3)
    return df.apply(lambda col: col.str.strip())


def fetch_dr16q(strips: int = 36) -> dict:
    if DR16Q_PATH.exists():
        print(f"  {DR16Q_PATH.name}: cached")
        return json.loads(DR16Q_PATH.with_suffix(".json").read_text())
    if I3_DR16Q_PARQUET.exists():
        print(f"  DR16Q from project I3's cache ({I3_DR16Q_PARQUET.relative_to(I3_DR16Q_PARQUET.parents[3])})")
        df = pd.read_parquet(I3_DR16Q_PARQUET, columns=["recno", "RAJ2000", "DEJ2000", "z"])
        origin = f"project I3 cache {I3_DR16Q_PARQUET.name} (VizieR {DR16Q_VIZIER})"
    else:
        print(f"  DR16Q from VizieR {DR16Q_VIZIER} in {strips} RA strips ...", flush=True)
        DR16Q_STRIP_DIR.mkdir(parents=True, exist_ok=True)
        edges = np.linspace(0.0, 360.0, strips + 1)
        frames = []
        for k in range(strips):
            path = DR16Q_STRIP_DIR / f"dr16q_ra{k:02d}.tsv"
            if not path.exists():
                text = http_get(VIZIER_ASU_TSV, params={
                    "-source": DR16Q_VIZIER, "-out.max": "unlimited",
                    "-out": "recno,RAJ2000,DEJ2000,z",
                    "RAJ2000": f"{edges[k]:.6f}..{edges[k + 1]:.6f}"}).decode("utf-8")
                if not text.lstrip().startswith("#"):
                    raise RuntimeError(f"unexpected VizieR answer: {text[:200]!r}")
                path.write_text(text, encoding="utf-8")
            frames.append(parse_vizier_tsv(path.read_text(encoding="utf-8")))
        df = pd.concat(frames, ignore_index=True).apply(pd.to_numeric)
        df = df.drop_duplicates("recno").reset_index(drop=True)
        origin = f"VizieR {DR16Q_VIZIER}, {strips} RA strips"
    if len(df) != N_DR16Q:
        raise RuntimeError(f"DR16Q: {len(df)} rows, expected {N_DR16Q}")
    df = df.rename(columns={"RAJ2000": "ra", "DEJ2000": "dec"}).sort_values("recno")
    df.to_parquet(DR16Q_PATH, index=False)
    prov = {"rows": len(df), "origin": origin, "written_utc": ga.now_utc(),
            "citation": "Lyke et al. 2020, ApJS 250, 8 (SDSS DR16Q)"}
    DR16Q_PATH.with_suffix(".json").write_text(json.dumps(prov, indent=2) + "\n")
    return prov


# --------------------------------------------------------------------------
# 3. The frame rotator's sources
# --------------------------------------------------------------------------
def fetch_frame_rotator() -> dict:
    side = FRAME_ROTATOR_PATH.with_suffix(".json")
    if FRAME_ROTATOR_PATH.exists() and side.exists():
        print(f"  {FRAME_ROTATOR_PATH.name}: cached")
        return json.loads(side.read_text())
    table = ga.sync_query(FRAME_ROTATOR_ADQL, timeout=600)
    n_used = int(np.sum(np.asarray(table["used_for_reference_frame_spin"]).astype(bool)))
    if n_used != N_FRAME_ROTATOR_SPIN_USED:
        raise RuntimeError(f"frame rotator: {n_used} spin sources, expected {N_FRAME_ROTATOR_SPIN_USED}")
    ga.write_with_provenance(table, FRAME_ROTATOR_PATH, {
        "source": "ESA Gaia archive, gaiadr3.frame_rotator_source", "adql": FRAME_ROTATOR_ADQL,
        "used_for_spin": n_used,
        "citation": "Lindegren et al. 2021, A&A 649, A2 (the EDR3 astrometric solution)"})
    print(f"  {len(table):,} rows, {n_used:,} used for the spin")
    return json.loads(side.read_text())


# --------------------------------------------------------------------------
# 4. Quasar look-alikes (mostly stars)
# --------------------------------------------------------------------------
def _lookalike_slice(lo: int, hi: int, depth: int = 0) -> list:
    path = CHUNK_DIR / "lookalikes" / f"ri_{lo:010d}_{hi:010d}.fits"
    if path.exists():
        return [path]
    try:
        table = ga.async_query(LOOKALIKE_ADQL.format(lo=lo, hi=hi), poll_s=5, quiet=True)
    except ga.TransientError as exc:
        if depth >= 4:
            raise
        mid = (lo + hi) // 2
        print(f"    slice {lo:,}-{hi:,}: {str(exc)[:60]} -> halving", flush=True)
        return _lookalike_slice(lo, mid, depth + 1) + _lookalike_slice(mid, hi, depth + 1)
    path.parent.mkdir(parents=True, exist_ok=True)
    table.write(path.with_suffix(".part"), format="fits", overwrite=True)
    path.with_suffix(".part").rename(path)
    print(f"    random_index {lo:>10,}-{hi:>10,}: {len(table):>7,} rows", flush=True)
    return [path]


def fetch_lookalikes() -> dict:
    side = LOOKALIKE_PATH.with_suffix(".json")
    if LOOKALIKE_PATH.exists() and side.exists():
        print(f"  {LOOKALIKE_PATH.name}: cached")
        return json.loads(side.read_text())
    edges = np.linspace(0, LOOKALIKE_N_RANDOM, LOOKALIKE_SLICES + 1).astype(np.int64)
    paths = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        paths += _lookalike_slice(int(lo), int(hi))
    from astropy.table import vstack
    table = vstack([Table.read(p) for p in paths], metadata_conflicts="silent")
    if len(np.unique(np.asarray(table["source_id"]))) != len(table):
        raise RuntimeError("duplicate look-alikes")
    ga.write_with_provenance(table, LOOKALIKE_PATH, {
        "source": "ESA Gaia archive (async TAP), Gaia DR3 gaia_source, first 0.5% of random_index",
        "adql_template": LOOKALIKE_ADQL, "slices": [p.name for p in paths]})
    for p in paths:
        p.unlink()
    return json.loads(side.read_text())


# --------------------------------------------------------------------------
# 5. Quaia
# --------------------------------------------------------------------------
def fetch_quaia() -> dict:
    side = QUAIA_PATH.with_suffix(".json")
    if QUAIA_PATH.exists() and side.exists():
        print(f"  {QUAIA_FILE}: cached")
        return json.loads(side.read_text())
    record = json.loads(http_get(QUAIA_RECORD))
    entry = next(f for f in record["files"] if f["key"] == QUAIA_FILE)
    print(f"  {QUAIA_FILE} from Zenodo ({entry['size'] / 1e6:.0f} MB) ...", flush=True)
    body = http_get(entry["links"]["self"], timeout=1800)
    md5 = hashlib.md5(body).hexdigest()
    if f"md5:{md5}" != entry["checksum"] or len(body) != entry["size"]:
        raise RuntimeError(f"{QUAIA_FILE}: checksum/size mismatch ({md5}, {len(body)})")
    tmp = QUAIA_PATH.with_suffix(".part")
    tmp.write_bytes(body)
    tmp.rename(QUAIA_PATH)
    from astropy.io import fits
    rows = int(fits.getheader(QUAIA_PATH, 1)["NAXIS2"])
    prov = {"rows": rows, "origin": f"{QUAIA_RECORD} (version {record['metadata'].get('version')})",
            "md5": md5, "written_utc": ga.now_utc(),
            "citation": "Storey-Fisher et al. 2024, ApJ 964, 69 (Quaia)"}
    side.write_text(json.dumps(prov, indent=2) + "\n")
    return prov


# --------------------------------------------------------------------------
# 6. Literature
# --------------------------------------------------------------------------
def fetch_literature() -> None:
    for name, url in LITERATURE.items():
        path = LIT_DIR / "klioner2021" / name
        if path.exists():
            print(f"  {name}: cached")
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        body = http_get(url)
        if not body.startswith(b"%PDF"):
            raise RuntimeError(f"{url} did not return a PDF")
        path.write_bytes(body)
        print(f"  {name}: {len(body) / 1e6:.1f} MB")


def write_manifest() -> None:
    entries = []
    for path in sorted([GAIA_PATH, DR16Q_PATH, FRAME_ROTATOR_PATH, LOOKALIKE_PATH, QUAIA_PATH,
                        *LIT_DIR.rglob("*.pdf")]):
        if not path.exists():
            continue
        side = path.with_suffix(".json")
        meta = json.loads(side.read_text()) if side.exists() and path.suffix != ".pdf" else {}
        entries.append({
            "file": str(path.relative_to(MANIFEST_PATH.parent)),
            "bytes": path.stat().st_size, "sha256": sha256(path),
            "rows": meta.get("rows"), "written_utc": meta.get("written_utc"),
        })
    MANIFEST_PATH.write_text(json.dumps({"files": entries, "manifest_utc": ga.now_utc()},
                                        indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--only", choices=["gaia", "dr16q", "frame_rotator", "lookalikes", "quaia",
                                           "literature"])
    args = parser.parse_args()
    steps = {"gaia": fetch_gaia, "dr16q": fetch_dr16q, "frame_rotator": fetch_frame_rotator,
             "lookalikes": fetch_lookalikes, "quaia": fetch_quaia, "literature": fetch_literature}
    for name, step in steps.items():
        if args.only and args.only != name:
            continue
        print(f"== {name}", flush=True)
        step()
    write_manifest()
    print("done.")


if __name__ == "__main__":
    main()
