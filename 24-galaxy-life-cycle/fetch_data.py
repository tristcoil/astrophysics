"""Download and cache everything the analysis needs.

Four things are fetched, all free and public:

1. **SDSS Main Galaxy Sample with MPA-JHU measurements** (SkyServer, DR17).
   For each of the ~650,000 galaxies the survey targeted and got a clean
   redshift for: position, redshift, photometry, and the MPA-JHU group's
   stellar mass, star-formation rate (SFR), emission-line fluxes and the
   4000 A break (Kauffmann et al. 2003; Brinchmann et al. 2004; Salim et al.
   2007; Tremonti et al. 2004).

2. **Every Main Galaxy Sample target**, observed or not (SkyServer table
   ``sdssTiledTargetAll``). The survey could not put a fibre on every galaxy
   it selected; comparing targets with spectra measures how complete the
   sample is, where, and how much sky it covers.

3. **GSWLC-2** (Salim et al. 2016, 2018): independent masses and SFRs for the
   same galaxies from fitting their ultraviolet-to-infrared light.

4. **Madau & Dickinson (2014)**: the LaTeX source of their review, whose
   Table 1 (cosmic star-formation-rate density, 0 < z < 8) and Table 2
   (stellar mass density) we parse ourselves.

Everything is cached under ``data/`` exactly as served (SkyServer CSV per RA
strip, the GSWLC .dat.gz files, the arXiv tarball), so an interrupted
download resumes where it stopped and nothing on disk is fetched twice.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
import tarfile
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import requests

from config import (
    DATA_RELEASE,
    FETCH_Z_MAX,
    FETCH_Z_MIN,
    GALAXY_CATALOG,
    GALAXY_RAW_DIR,
    GALSPEC_RAW_DIR,
    GSWLC_BASE_URL,
    GSWLC_DIR,
    GSWLC_FILES,
    MD14_ARXIV_SRC,
    MD14_DIR,
    MGS_R_LIMIT,
    MGS_TARGET_BITS,
    RA_CHUNK_DEG,
    SKYSERVER_URL,
    TARGET_CATALOG,
    TARGET_RAW_DIR,
    TARGET_SPEC_RAW_DIR,
    TILING_BASE_URL,
    TILING_DIR,
    TILING_FILES,
)

USER_AGENT = "astrophysics-atlas-J3/1.0 (hobby research; python-requests)"

#: (connect, read) timeouts in seconds. The connect timeout is short on purpose:
#: some hosts publish IPv6 addresses this machine cannot reach, and requests
#: tries each address in turn, so a long connect timeout would stall for
#: minutes before falling back to IPv4.
CONNECT_TIMEOUT = 15.0

# ---------------------------------------------------------------------------
# 1. Galaxies with MPA-JHU measurements
# ---------------------------------------------------------------------------

#: The Main Galaxy Sample selection is D2's: targeted as a main-sample galaxy,
#: spectroscopically a galaxy, clean redshift, best spectrum of the object,
#: dust-corrected Petrosian r < 17.77.
GALAXY_QUERY = """
SELECT s.specObjID, s.objID, s.plate, s.mjd, s.fiberID,
       s.ra, s.dec, s.z, s.zErr,
       s.petroMag_r, s.extinction_r,
       s.dered_u, s.dered_g, s.dered_r, s.dered_i, s.dered_z,
       s.modelMagErr_g, s.modelMagErr_r
FROM SpecPhoto s
WHERE (s.legacy_target1 & {bits}) != 0
  AND s.class = 'GALAXY'
  AND s.sciencePrimary = 1
  AND s.mode = 1
  AND s.zWarning = 0
  AND s.z BETWEEN {zmin} AND {zmax}
  AND s.petroMag_r - s.extinction_r < {rlim}
  AND s.ra >= {ra_lo} AND s.ra < {ra_hi}
"""

#: The MPA-JHU tables, fetched whole and joined to the galaxies locally.
#: Joining them on the server is slow (a 2-degree strip took > 2 minutes), but
#: each table alone, read by ranges of specObjID, takes seconds: a DR8+
#: specObjID is the plate number times 2^50 plus fibre and MJD bits, so a
#: plate range is a contiguous index range. The tables hold one row per
#: spectrum of the original survey (1,477,411 rows, plates 266-2974).
#:
#: MPA-JHU columns (all medians of a likelihood distribution unless noted):
#:   lgm_tot_p50       log stellar mass, whole galaxy (M_sun, Kroupa IMF)
#:   sfr_tot_p50       log SFR, whole galaxy (M_sun/yr), fibre + out-of-fibre
#:   sfr_fib_p50       log SFR inside the 3" fibre only
#:   bptclass          their emission-line class (1 SF, 2 low-S/N SF,
#:                     3 composite, 4 AGN, 5 low-S/N LINER, -1 unclassifiable)
#:   d4000_n           narrow 4000 A break: big = old stars
#:   reliable          0 flags a spectrum MPA-JHU consider unreliable
GALSPEC_TABLES = {
    "galSpecExtra": ["lgm_tot_p16", "lgm_tot_p50", "lgm_tot_p84", "sfr_tot_p2p5",
                     "sfr_tot_p16", "sfr_tot_p50", "sfr_tot_p84", "sfr_tot_p97p5",
                     "specsfr_tot_p50", "sfr_fib_p50", "lgm_fib_p50", "bptclass"],
    "galSpecInfo": ["reliable"],
    "galSpecIndx": ["d4000_n", "d4000_n_err", "lick_hd_a"],
    "galSpecLine": ["h_alpha_flux", "h_alpha_flux_err", "h_alpha_eqw", "h_beta_flux",
                    "h_beta_flux_err", "oiii_5007_flux", "oiii_5007_flux_err",
                    "nii_6584_flux", "nii_6584_flux_err"],
}
GALSPEC_QUERY = """
SELECT specObjID, {cols} FROM {table}
WHERE specObjID >= {lo} AND specObjID < {hi}
"""
#: Plate chunks: 100 plates (~55,000 spectra) per query.
GALSPEC_PLATE_EDGES = list(range(200, 3101, 100))

GALAXY_COLUMNS = {
    "specObjID": np.int64, "objID": np.int64, "plate": np.int32, "mjd": np.int32,
    "fiberID": np.int32, "ra": np.float64, "dec": np.float64, "z": np.float64,
    "zErr": np.float32, "petroMag_r": np.float32, "extinction_r": np.float32,
    "dered_u": np.float32, "dered_g": np.float32, "dered_r": np.float32,
    "dered_i": np.float32, "dered_z": np.float32,
    "modelMagErr_g": np.float32, "modelMagErr_r": np.float32,
    "lgm_tot_p16": np.float32, "lgm_tot_p50": np.float32, "lgm_tot_p84": np.float32,
    "sfr_tot_p2p5": np.float32, "sfr_tot_p16": np.float32, "sfr_tot_p50": np.float32,
    "sfr_tot_p84": np.float32, "sfr_tot_p97p5": np.float32,
    "specsfr_tot_p50": np.float32, "sfr_fib_p50": np.float32, "lgm_fib_p50": np.float32,
    "bptclass": np.int16, "reliable": np.int16,
    "d4000_n": np.float32, "d4000_n_err": np.float32, "lick_hd_a": np.float32,
    "h_alpha_flux": np.float32, "h_alpha_flux_err": np.float32, "h_alpha_eqw": np.float32,
    "h_beta_flux": np.float32, "h_beta_flux_err": np.float32,
    "oiii_5007_flux": np.float32, "oiii_5007_flux_err": np.float32,
    "nii_6584_flux": np.float32, "nii_6584_flux_err": np.float32,
}

# ---------------------------------------------------------------------------
# 2. Every Main Galaxy Sample target
# ---------------------------------------------------------------------------

#: ``sdssTiledTargetAll`` lists every object SDSS-I/II's tiling software was
#: asked to put a fibre on, whether or not it got one. ``specObjID`` is 0 when
#: it did not. The photometry is joined from ``PhotoObjAll`` (DR17's version of
#: the imaging) so that targets and galaxies are cut on the same magnitudes.
TARGET_QUERY = """
SELECT t.targetID, t.photoObjID, t.specObjID, t.primTarget, t.tiMaskAll,
       t.ra, t.dec, p.petroMag_r, p.extinction_r
FROM sdssTiledTargetAll t
LEFT OUTER JOIN PhotoObjAll p ON p.objID = t.photoObjID
WHERE (t.primTarget & {bits}) != 0
  AND t.ra >= {ra_lo} AND t.ra < {ra_hi}
"""

#: For targets that were observed but are not clean galaxies in our sample:
#: what their spectrum turned out to be. Stars, quasars and galaxies outside
#: our redshift range legitimately have no place in the sample; failed
#: redshifts are genuine incompleteness. (Joining on the primary key is fast.)
TARGET_SPEC_QUERY = """
SELECT t.targetID, s.specObjID, s.class, s.zWarning, s.z
FROM sdssTiledTargetAll t JOIN SpecObjAll s ON s.specObjID = t.specObjID
WHERE (t.primTarget & {bits}) != 0
  AND (s.class != 'GALAXY' OR s.zWarning != 0 OR s.z < {zmin} OR s.z > {zmax})
  AND t.ra >= {ra_lo} AND t.ra < {ra_hi}
"""
#: spec_flag values stored per target
SPEC_FLAG = {"clean": 0, "star": 1, "qso": 2, "bad_z": 3, "z_low": 4, "z_high": 5}

TARGET_COLUMNS = {
    "targetID": np.int64, "photoObjID": np.int64, "specObjID": np.int64,
    "primTarget": np.int64, "tiMaskAll": np.int32, "ra": np.float64, "dec": np.float64,
    "petroMag_r": np.float32, "extinction_r": np.float32,
}


def skyserver_sql(sql: str, *, retries: int = 5, timeout: float = 900.0) -> str:
    """Run one SQL query on SkyServer and return the CSV text.

    SkyServer reports SQL errors inside a 200 response, so the body is checked,
    not just the status code. Transient failures are retried with backoff.
    """
    last_error = None
    for attempt in range(retries):
        try:
            r = requests.get(SKYSERVER_URL,
                             params={"cmd": " ".join(sql.split()), "format": "csv"},
                             headers={"User-Agent": USER_AGENT},
                             timeout=(CONNECT_TIMEOUT, timeout))
            if r.status_code >= 400:
                # SkyServer puts the SQL error in a JSON body: show it.
                raise RuntimeError(f"HTTP {r.status_code}: {r.text[:300]}")
            text = r.text
            if not text.startswith("#Table1"):
                raise RuntimeError(f"SkyServer returned an error: {text[:300]!r}")
            return text
        except (requests.RequestException, RuntimeError) as exc:
            last_error = exc
            if "Invalid column" in str(exc) or "Incorrect syntax" in str(exc):
                break  # a bug in the query, not a transient failure
            wait = 10 * 2**attempt
            print(f"    attempt {attempt + 1} failed ({exc}); retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"SkyServer query failed after {retries} attempts: {last_error}")


def _chunks() -> list[tuple[float, float]]:
    edges = np.arange(0.0, 360.0 + RA_CHUNK_DEG / 2, RA_CHUNK_DEG)
    return [(float(lo), float(hi)) for lo, hi in zip(edges[:-1], edges[1:])]


def _fetch_strips(kind: str, template: str, raw_dir, fmt: dict, *, refresh: bool,
                  quiet: bool) -> None:
    chunks = _chunks()
    if not quiet:
        print(f"SDSS {DATA_RELEASE} {kind}: {len(chunks)} strips of {RA_CHUNK_DEG:g} deg in RA")
    for i, (lo, hi) in enumerate(chunks, 1):
        path = raw_dir / f"ra{lo:05.1f}-{hi:05.1f}.csv"
        if path.exists() and not refresh:
            if not quiet:
                print(f"  [{i:2d}/{len(chunks)}] RA {lo:5.1f}-{hi:5.1f}: cached")
            continue
        t0 = time.time()
        text = skyserver_sql(template.format(ra_lo=lo, ra_hi=hi, **fmt))
        # Temporary name first, so a crash never leaves a truncated strip that
        # a later run would mistake for a complete one.
        tmp = path.with_suffix(".part")
        tmp.write_text(text)
        tmp.replace(path)
        if not quiet:
            print(f"  [{i:2d}/{len(chunks)}] RA {lo:5.1f}-{hi:5.1f}: "
                  f"{text.count(chr(10)) - 2:6,d} rows in {time.time() - t0:5.1f}s", flush=True)
        time.sleep(1.0)  # be polite to a free public service


def fetch_galaxies(*, refresh: bool = False, quiet: bool = False) -> None:
    _fetch_strips("Main Galaxy Sample", GALAXY_QUERY, GALAXY_RAW_DIR,
                  {"bits": MGS_TARGET_BITS, "zmin": FETCH_Z_MIN, "zmax": FETCH_Z_MAX,
                   "rlim": MGS_R_LIMIT}, refresh=refresh, quiet=quiet)


def _galspec_path(table: str, lo: int, hi: int):
    return GALSPEC_RAW_DIR / f"{table}_plates{lo:04d}-{hi:04d}.csv"


def fetch_galspec(*, refresh: bool = False, quiet: bool = False) -> None:
    """The four MPA-JHU tables, in plate chunks (see GALSPEC_TABLES)."""
    chunks = list(zip(GALSPEC_PLATE_EDGES[:-1], GALSPEC_PLATE_EDGES[1:]))
    for table, cols in GALSPEC_TABLES.items():
        if not quiet:
            print(f"MPA-JHU {table}: {len(chunks)} plate chunks")
        for lo, hi in chunks:
            path = _galspec_path(table, lo, hi)
            if path.exists() and not refresh:
                continue
            t0 = time.time()
            sql = GALSPEC_QUERY.format(cols=", ".join(cols), table=table,
                                       lo=lo * 2**50, hi=hi * 2**50)
            text = skyserver_sql(sql)
            tmp = path.with_suffix(".part")
            tmp.write_text(text)
            tmp.replace(path)
            if not quiet:
                print(f"  plates {lo}-{hi}: {text.count(chr(10)) - 2:6,d} rows in "
                      f"{time.time() - t0:4.1f}s", flush=True)
            time.sleep(0.5)


def _merge_galspec() -> pd.DataFrame:
    """All four tables joined on specObjID (one row per spectrum)."""
    merged = None
    chunks = list(zip(GALSPEC_PLATE_EDGES[:-1], GALSPEC_PLATE_EDGES[1:]))
    for table in GALSPEC_TABLES:
        frames = [pd.read_csv(_galspec_path(table, lo, hi), skiprows=1) for lo, hi in chunks]
        df = pd.concat([f for f in frames if len(f)], ignore_index=True)
        df = df.drop_duplicates("specObjID")
        merged = df if merged is None else merged.merge(df, on="specObjID", how="outer")
    return merged


def fetch_targets(*, refresh: bool = False, quiet: bool = False) -> None:
    _fetch_strips("Main Galaxy Sample targets", TARGET_QUERY, TARGET_RAW_DIR,
                  {"bits": MGS_TARGET_BITS}, refresh=refresh, quiet=quiet)
    _fetch_strips("spectra of targets that are not clean galaxies", TARGET_SPEC_QUERY,
                  TARGET_SPEC_RAW_DIR, {"bits": MGS_TARGET_BITS, "zmin": FETCH_Z_MIN,
                                        "zmax": FETCH_Z_MAX}, refresh=refresh, quiet=quiet)


def _target_spec_flags() -> pd.DataFrame:
    frames = [pd.read_csv(TARGET_SPEC_RAW_DIR / f"ra{lo:05.1f}-{hi:05.1f}.csv", skiprows=1)
              for lo, hi in _chunks()]
    df = pd.concat([f for f in frames if len(f)], ignore_index=True)
    cls = df["class"].str.strip().str.upper()
    flag = np.full(len(df), SPEC_FLAG["bad_z"], dtype=np.int16)
    gal_ok = (cls == "GALAXY") & (df["zWarning"] == 0)
    flag[gal_ok & (df["z"] < FETCH_Z_MIN)] = SPEC_FLAG["z_low"]
    flag[gal_ok & (df["z"] > FETCH_Z_MAX)] = SPEC_FLAG["z_high"]
    flag[cls == "STAR"] = SPEC_FLAG["star"]
    flag[cls == "QSO"] = SPEC_FLAG["qso"]
    out = pd.DataFrame({"targetID": df["targetID"].to_numpy(np.int64), "spec_flag": flag})
    return out.drop_duplicates("targetID")


def _merge_strips(raw_dir, columns: dict, key: str, extra: pd.DataFrame | None = None
                  ) -> tuple[dict, dict]:
    frames = []
    for lo, hi in _chunks():
        path = raw_dir / f"ra{lo:05.1f}-{hi:05.1f}.csv"
        if not path.exists():
            raise FileNotFoundError(f"{path} missing -- run `python fetch_data.py`")
        frame = pd.read_csv(path, skiprows=1)
        if len(frame):  # strips across the Milky Way's disc can be empty
            frames.append(frame)
    df = pd.concat(frames, ignore_index=True)
    n_raw = len(df)
    df = df.drop_duplicates(key)
    if extra is not None:  # LEFT join: every galaxy kept, matched or not
        df = df.merge(extra, on=key, how="left")
    # A LEFT join with no match gives an empty field; store it as the sentinel.
    df = df.fillna(-9999)
    missing = set(columns) - set(df.columns)
    if missing:
        raise RuntimeError(f"columns missing from SkyServer output: {sorted(missing)}")
    arrays = {c: df[c].to_numpy(dtype=t) for c, t in columns.items()}
    newest = max(p.stat().st_mtime for p in raw_dir.glob("ra*.csv"))
    meta = {"n_rows_served": n_raw, "n_duplicates_dropped": n_raw - len(df),
            "n_rows": len(df),
            "retrieved_utc": datetime.fromtimestamp(newest, timezone.utc)
            .isoformat(timespec="seconds")}
    return arrays, meta


def build_catalogs(*, quiet: bool = False) -> dict:
    """Merge the cached strips into two typed catalogues, with provenance."""
    galspec = _merge_galspec()
    gal, gmeta = _merge_strips(GALAXY_RAW_DIR, GALAXY_COLUMNS, "specObjID", extra=galspec)
    gmeta.update(source=f"SDSS {DATA_RELEASE} SkyServer: SpecPhoto, left-joined locally to "
                 "galSpecExtra/Info/Indx/Line", url=SKYSERVER_URL,
                 query=" ".join(GALAXY_QUERY.split()),
                 galspec_query=" ".join(GALSPEC_QUERY.split()),
                 n_galspec_rows=int(len(galspec)),
                 n_without_galspec=int(np.sum(gal["reliable"] < -9000)))
    np.savez_compressed(GALAXY_CATALOG, meta=json.dumps(gmeta), **gal)
    tgt, tmeta = _merge_strips(TARGET_RAW_DIR, {**TARGET_COLUMNS, "spec_flag": np.int16},
                               "targetID", extra=_target_spec_flags())
    tgt["spec_flag"] = np.where(tgt["spec_flag"] < 0, 0, tgt["spec_flag"]).astype(np.int16)
    tmeta.update(source=f"SDSS {DATA_RELEASE} SkyServer: sdssTiledTargetAll + PhotoObjAll",
                 url=SKYSERVER_URL, query=" ".join(TARGET_QUERY.split()))
    np.savez_compressed(TARGET_CATALOG, meta=json.dumps(tmeta), **tgt)
    if not quiet:
        print(f"  galaxies: {gmeta['n_rows_served']:,} rows -> {gmeta['n_rows']:,}")
        print(f"  targets:  {tmeta['n_rows_served']:,} rows -> {tmeta['n_rows']:,}")
    return {"galaxies": gmeta, "targets": tmeta}


def load_npz(path) -> tuple[dict, dict]:
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    with np.load(path) as f:
        meta = json.loads(str(f["meta"]))
        arrays = {k: f[k] for k in f.files if k != "meta"}
    return arrays, meta


def load_galaxies():
    return load_npz(GALAXY_CATALOG)


def load_targets():
    return load_npz(TARGET_CATALOG)


# ---------------------------------------------------------------------------
# 3. GSWLC-2
# ---------------------------------------------------------------------------

#: Column layout of every GSWLC-2 file (the authors' table2.pdf).
#: Missing values are -99.
GSWLC_NAMES = [
    "objid", "glxid", "plate", "mjd", "fiberid", "ra", "dec", "z", "chi2_r",
    "logmstar", "logmstar_err", "logsfr", "logsfr_err", "afuv", "afuv_err",
    "ab", "ab_err", "av", "av_err", "flag_sed", "uv_survey", "flag_uv",
    "flag_midir", "flag_mgs",
]


def _download(url: str, path, *, quiet: bool, refresh: bool) -> None:
    if path.exists() and not refresh:
        if not quiet:
            print(f"  {path.name}: cached ({path.stat().st_size / 1e6:.1f} MB)")
        return
    t0 = time.time()
    for attempt in range(5):
        try:
            with requests.get(url, stream=True, timeout=(CONNECT_TIMEOUT, 300),
                              headers={"User-Agent": USER_AGENT}) as r:
                r.raise_for_status()
                tmp = path.with_suffix(path.suffix + ".part")
                # Save the bytes exactly as served. The GSWLC server labels its
                # .gz files "Content-Encoding: x-gzip", and requests would
                # otherwise unzip them on the fly into a file still named .gz.
                with open(tmp, "wb") as f:
                    for block in r.raw.stream(1 << 20, decode_content=False):
                        f.write(block)
                tmp.replace(path)
            break
        except requests.RequestException as exc:
            wait = 10 * 2**attempt
            print(f"    attempt {attempt + 1} failed ({exc}); retrying in {wait}s")
            time.sleep(wait)
    else:
        raise RuntimeError(f"download failed: {url}")
    if not quiet:
        print(f"  {path.name}: {path.stat().st_size / 1e6:.1f} MB in {time.time() - t0:.0f}s")


def fetch_tiling(*, refresh: bool = False, quiet: bool = False) -> None:
    if not quiet:
        print("SDSS legacy tiling geometry (SAS tiling/final)")
    for name in TILING_FILES:
        _download(TILING_BASE_URL + name, TILING_DIR / name, quiet=quiet, refresh=refresh)


def fetch_gswlc(*, refresh: bool = False, quiet: bool = False) -> None:
    if not quiet:
        print("GSWLC-2 (Salim et al. 2016, 2018)")
    for name in GSWLC_FILES:
        _download(GSWLC_BASE_URL + name, GSWLC_DIR / name, quiet=quiet, refresh=refresh)


def load_gswlc(name: str = "GSWLC-X2.dat.gz") -> pd.DataFrame:
    path = GSWLC_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    return pd.read_csv(path, sep=r"\s+", header=None, names=GSWLC_NAMES,
                       compression="gzip")


# ---------------------------------------------------------------------------
# 4. Madau & Dickinson (2014): Tables 1 and 2 from the LaTeX source
# ---------------------------------------------------------------------------

MD14_TEX = MD14_DIR / "paper.tex"


def fetch_md14(*, refresh: bool = False, quiet: bool = False) -> None:
    """Keep only paper.tex from the arXiv tarball (the figures are not needed)."""
    if MD14_TEX.exists() and not refresh:
        if not quiet:
            print("Madau & Dickinson (2014) LaTeX source: cached")
        return
    r = requests.get(MD14_ARXIV_SRC, timeout=(CONNECT_TIMEOUT, 300),
                     headers={"User-Agent": USER_AGENT})
    r.raise_for_status()
    with tarfile.open(fileobj=io.BytesIO(r.content), mode="r:gz") as tar:
        member = tar.getmember("paper.tex")
        MD14_TEX.write_bytes(tar.extractfile(member).read())
    if not quiet:
        print(f"Madau & Dickinson (2014) LaTeX source: {MD14_TEX.stat().st_size / 1e3:.0f} kB "
              f"from {MD14_ARXIV_SRC}")




def _table_body(tex: str, label: str) -> list[str]:
    start = tex.index(label)
    body = tex[start:tex.index(r"\end{tabular*}", start)]
    return body.splitlines()


def parse_md14_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Table 1 (SFR density) and Table 2 (stellar mass density), Salpeter IMF.

    Each row: reference, z_lo, z_hi, log value, +err, -err (and for Table 1 the
    dust attenuation A_FUV, NaN for the infrared rows). Continuation rows with
    an empty reference column inherit the reference above them. Table 2 runs
    over two LaTeX tables ("Table 2 (cont.)"); both are read.
    """
    tex = MD14_TEX.read_text(encoding="latin-1")  # the source is not UTF-8

    def parse(lines: list[str], with_afuv: bool, band_split: bool) -> list[dict]:
        rows, ref, band = [], None, "UV"
        for line in lines:
            if band_split and line.strip().startswith(r"\hline") and rows:
                band = "IR"  # the double rule separating UV from IR rows
            if "&" not in line or "multicolumn" in line:
                continue
            # Split on unescaped "&" only: "Robotham \& Driver" is one name.
            cells = [c.strip() for c in re.split(r"(?<!\\)&", line.split(r"\\")[0])]
            name = cells[0]
            if name:
                ref = (name.replace(r"\etal", "et al.").replace(r"\&", "&")
                       .replace(r"{\'{a}}", "á").replace(r"\'{a}", "á")
                       .replace(r"\'{e}", "é").replace("$^{b}$", "").replace("$^c$", "")
                       .replace(r"\,", "").strip())
            zfield = cells[1]
            if "-" in zfield:
                z_lo, z_hi = (float(v) for v in zfield.split("-"))
            else:
                z_lo = z_hi = float(zfield)
            valcell = cells[3] if with_afuv else cells[2]
            m = re.match(r"\$?(-?[0-9.]+)(?:\^\{\+?([0-9.]+)\}_\{-?([0-9.]+)\})?\$?", valcell)
            if m is None:
                raise ValueError(f"cannot parse value cell {valcell!r}")
            val = float(m.group(1))
            up = float(m.group(2)) if m.group(2) else np.nan
            lo = float(m.group(3)) if m.group(3) else np.nan
            row = {"reference": ref, "z_lo": z_lo, "z_hi": z_hi, "log_value": val,
                   "err_plus": up, "err_minus": lo}
            if with_afuv:
                a = cells[2]
                row["a_fuv"] = float(a) if re.fullmatch(r"[0-9.]+", a) else np.nan
                row["band"] = band
            rows.append(row)
        return rows

    t1_lines = _table_body(tex, r"\label{tab:SFRD}")
    # The UV/IR split is the "\hline\hline" pair after the first data rows.
    t1 = parse(t1_lines, with_afuv=True, band_split=True)
    t2a = parse(_table_body(tex, r"\label{tab:SMD}"), with_afuv=False, band_split=False)
    cont = tex.index(r"Table 2: (cont.)")
    t2b_lines = tex[cont:tex.index(r"\end{tabular*}", cont)].splitlines()
    t2b = parse(t2b_lines, with_afuv=False, band_split=False)
    sfrd = pd.DataFrame(t1)
    smd = pd.DataFrame(t2a + t2b)
    return sfrd, smd


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true",
                        help="re-download even what is already cached")
    args = parser.parse_args(argv)

    fetch_md14(refresh=args.refresh)
    sfrd, smd = parse_md14_tables()
    print(f"  Table 1: {len(sfrd)} SFR-density points; Table 2: {len(smd)} mass-density points")
    print()
    fetch_gswlc(refresh=args.refresh)
    print()
    fetch_tiling(refresh=args.refresh)
    print()
    fetch_targets(refresh=args.refresh)
    print()
    fetch_galaxies(refresh=args.refresh)
    print()
    fetch_galspec(refresh=args.refresh)
    print("\nMerging strips")
    build_catalogs()
    print("\nDone. The analysis can now run offline.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
