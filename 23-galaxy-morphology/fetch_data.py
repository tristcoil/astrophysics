"""Download and cache everything the analysis needs.

All of it is free and public:

1. **Galaxy Zoo 2 votes** (Hart et al. 2016 debiasing): for 239,695 galaxies,
   how many volunteers chose each answer to each question, and the vote
   fractions corrected for the "distant galaxies look smoother" bias.
2. **SDSS DR17 data for every GZ2 galaxy**, joined on SkyServer: redshift,
   photometry, size and shape, plus the MPA-JHU stellar masses and star
   formation rates, and the original Willett et al. (2013) debiased fractions
   for comparison.
3. **The SDSS Main Galaxy Sample**, fetched with exactly project D2's query:
   the galaxies whose crowding defines each GZ2 galaxy's environment.
4. **Expert classifications** by Nair & Abraham (2010), from VizieR.
5. **The K-correction authors' code**, used only by the tests.

Images (for the tuning-fork poster, the explorer page and the neural network)
and the Galaxy Zoo DESI predictions are fetched later, by :func:`fetch_cutouts`
and :func:`fetch_gzdesi_subset`, once the analysis has decided which galaxies
it needs.

Everything is cached under ``data/`` as served, and every download writes to a
temporary name first, so an interrupted run resumes and never leaves a
truncated file behind. Once the cache is warm the analysis runs offline.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import io
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from config import (
    CATALOG_PATH, CUTOUT_DIR, DATA_DIR, DATA_RELEASE, GZ2_RA_CHUNK_DEG, GZ2_SDSS_DIR, GZ_DIR,
    GZDESI_DIR, GZDESI_URL, HART16_URL, IMGCUTOUT_URL, KCOR_CODE_URL, MGS_DIR, MGS_FETCH_Z,
    MGS_PATH, MGS_R_LIMIT, MGS_RA_CHUNK_DEG, MGS_TARGET_BITS, NA10_CATALOG, SKYSERVER_URL,
    VIZIER_DIR,
)

USER_AGENT = "astrophysics-atlas/J2 (hobbyist research script; python-requests)"
HART16_PATH = GZ_DIR / "gz2_hart16.csv.gz"
NA10_PATH = VIZIER_DIR / "nair_abraham_2010.csv"
KCOR_PATH = DATA_DIR / "calc_kcor.py"

#: The 37 answers of the GZ2 decision tree, as named in the catalogue columns:
#: task, answer -> short name used in the code.
ANSWERS = {
    "t01_smooth_or_features_a01_smooth": "smooth",
    "t01_smooth_or_features_a02_features_or_disk": "features",
    "t01_smooth_or_features_a03_star_or_artifact": "artifact",
    "t02_edgeon_a04_yes": "edgeon",
    "t02_edgeon_a05_no": "notedgeon",
    "t03_bar_a06_bar": "bar",
    "t03_bar_a07_no_bar": "nobar",
    "t04_spiral_a08_spiral": "spiral",
    "t04_spiral_a09_no_spiral": "nospiral",
    "t05_bulge_prominence_a10_no_bulge": "bulge_none",
    "t05_bulge_prominence_a11_just_noticeable": "bulge_noticeable",
    "t05_bulge_prominence_a12_obvious": "bulge_obvious",
    "t05_bulge_prominence_a13_dominant": "bulge_dominant",
    "t06_odd_a14_yes": "odd",
    "t06_odd_a15_no": "notodd",
    "t07_rounded_a16_completely_round": "round",
    "t07_rounded_a17_in_between": "inbetween",
    "t07_rounded_a18_cigar_shaped": "cigar",
    "t08_odd_feature_a19_ring": "ring",
    "t08_odd_feature_a20_lens_or_arc": "lens",
    "t08_odd_feature_a21_disturbed": "disturbed",
    "t08_odd_feature_a22_irregular": "irregular",
    "t08_odd_feature_a23_other": "other",
    "t08_odd_feature_a24_merger": "merger",
    "t08_odd_feature_a38_dust_lane": "dustlane",
    "t09_bulge_shape_a25_rounded": "edgebulge_round",
    "t09_bulge_shape_a26_boxy": "edgebulge_boxy",
    "t09_bulge_shape_a27_no_bulge": "edgebulge_none",
    "t10_arms_winding_a28_tight": "arms_tight",
    "t10_arms_winding_a29_medium": "arms_medium",
    "t10_arms_winding_a30_loose": "arms_loose",
    "t11_arms_number_a31_1": "arms_1",
    "t11_arms_number_a32_2": "arms_2",
    "t11_arms_number_a33_3": "arms_3",
    "t11_arms_number_a34_4": "arms_4",
    "t11_arms_number_a36_more_than_4": "arms_5plus",
    "t11_arms_number_a37_cant_tell": "arms_cant_tell",
}

#: Answers whose original Willett et al. (2013) debiased fraction is also
#: fetched, to compare the two debiasing methods.
W13_ANSWERS = ["smooth", "features", "edgeon", "notedgeon", "bar", "nobar", "spiral",
               "nospiral", "arms_1", "arms_2", "arms_3", "arms_4", "arms_5plus",
               "arms_cant_tell"]


# ---------------------------------------------------------------------------
# Generic helpers
# ---------------------------------------------------------------------------

def download(url: str, path: Path, *, retries: int = 5, quiet: bool = False) -> Path:
    """Stream ``url`` to ``path``, resuming a partial download if one exists."""
    if path.exists():
        return path
    part = path.with_suffix(path.suffix + ".part")
    for attempt in range(retries):
        try:
            have = part.stat().st_size if part.exists() else 0
            headers = {"User-Agent": USER_AGENT}
            if have:
                headers["Range"] = f"bytes={have}-"
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                if r.status_code == 416:  # already complete
                    break
                r.raise_for_status()
                mode = "ab" if (have and r.status_code == 206) else "wb"
                total = int(r.headers.get("Content-Length", 0)) + (have if mode == "ab" else 0)
                t0, done = time.time(), (have if mode == "ab" else 0)
                with open(part, mode) as fh:
                    for block in r.iter_content(1 << 20):
                        fh.write(block)
                        done += len(block)
                if total and done < total:
                    raise IOError(f"short read: {done} of {total} bytes")
                if not quiet:
                    print(f"    {path.name}: {done / 1e6:.1f} MB in {time.time() - t0:.0f}s")
            part.replace(path)
            return path
        except (requests.RequestException, IOError) as exc:
            wait = 5 * 2**attempt
            print(f"    attempt {attempt + 1} for {path.name} failed ({exc}); retry in {wait}s")
            time.sleep(wait)
    if part.exists() and not path.exists():
        part.replace(path)
    if not path.exists():
        raise RuntimeError(f"could not download {url}")
    return path


def skyserver_sql(sql: str, *, retries: int = 5, timeout: float = 600.0) -> str:
    """Run one SQL query on SkyServer and return the CSV text.

    SkyServer reports SQL errors inside a 200 response, so the body is checked,
    not just the status code. Transient failures are retried with backoff.
    """
    last = None
    for attempt in range(retries):
        try:
            r = requests.get(SKYSERVER_URL, params={"cmd": " ".join(sql.split()), "format": "csv"},
                             headers={"User-Agent": USER_AGENT}, timeout=timeout)
            r.raise_for_status()
            if not r.text.startswith("#Table1"):
                raise RuntimeError(f"SkyServer returned an error: {r.text[:300]!r}")
            return r.text
        except (requests.RequestException, RuntimeError) as exc:
            last = exc
            wait = 5 * 2**attempt
            print(f"    attempt {attempt + 1} failed ({exc}); retrying in {wait}s")
            time.sleep(wait)
    raise RuntimeError(f"SkyServer query failed after {retries} attempts: {last}")


def _parse_skyserver(text: str) -> pd.DataFrame:
    """SkyServer CSV: a ``#Table1`` line, then a normal header and rows."""
    return pd.read_csv(io.StringIO(text), skiprows=1)


def _strips(width: float) -> list[tuple[float, float]]:
    edges = np.arange(0.0, 360.0 + width / 2, width)
    return [(float(lo), float(hi)) for lo, hi in zip(edges[:-1], edges[1:])]


def _fetch_strips(template: str, out_dir: Path, prefix: str, width: float, fmt: dict,
                  quiet: bool) -> None:
    strips = _strips(width)
    for i, (lo, hi) in enumerate(strips, 1):
        path = out_dir / f"{prefix}_ra{lo:05.1f}-{hi:05.1f}.csv"
        if path.exists():
            continue
        t0 = time.time()
        text = skyserver_sql(template.format(ra_lo=lo, ra_hi=hi, **fmt))
        tmp = path.with_suffix(".part")
        tmp.write_text(text)
        tmp.replace(path)
        if not quiet:
            print(f"  [{i:2d}/{len(strips)}] RA {lo:5.1f}-{hi:5.1f}: "
                  f"{text.count(chr(10)) - 2:6,d} rows in {time.time() - t0:4.1f}s")
        time.sleep(1.0)  # be polite to a free public service


# ---------------------------------------------------------------------------
# 1-2. Galaxy Zoo 2 and its SDSS data
# ---------------------------------------------------------------------------

#: GZ2 joined with SDSS DR17 on the server, one right-ascension strip at a time.
#:
#: * ``zoo2MainSpecz`` is Willett et al.'s Table 5 (243,500 galaxies), keyed by
#:   the DR7 object ID (to match the Hart et al. table), the DR8+ object ID
#:   (photometry) and the spectrum ID (redshift, stellar mass).
#: * Positions and photometry come from ``PhotoObjAll``: ``petroMag_r`` is the
#:   flux GZ2 was selected on; ``dered_*`` are model magnitudes with Milky Way
#:   dust removed (the SDSS recommendation for colours); ``petroR50/R90`` are
#:   the radii holding 50% and 90% of the light; ``expAB_r`` is the axis ratio
#:   b/a of an exponential-disk fit (how inclined a disk is); ``fracDeV_r`` the
#:   fraction of light in a de Vaucouleurs (bulge-like) profile.
#: * ``galSpecExtra`` holds the MPA-JHU stellar masses (fits of stellar
#:   population models to the ugriz photometry; Kauffmann et al. 2003, Salim
#:   et al. 2007) and star formation rates (Brinchmann et al. 2004).
GZ2_QUERY = """
SELECT zz.dr7objid, zz.specobjid, zz.dr8objid, zz.sample,
       p.ra, p.dec, s.z, s.zErr, s.zWarning, s.class,
       p.petroMag_r, p.extinction_r, p.petroR50_r, p.petroR90_r,
       p.expAB_r, p.deVAB_r, p.fracDeV_r,
       p.dered_u, p.dered_g, p.dered_r, p.dered_i, p.dered_z,
       p.modelMagErr_u, p.modelMagErr_g, p.modelMagErr_r, p.modelMagErr_i,
       g.lgm_tot_p16, g.lgm_tot_p50, g.lgm_tot_p84,
       g.sfr_tot_p50, g.specsfr_tot_p50, g.bptclass,
       zz.t01_smooth_or_features_a01_smooth_weighted_fraction AS w13wf_smooth,
       {w13_columns}
FROM zoo2MainSpecz AS zz
JOIN SpecObjAll AS s ON s.specObjID = zz.specobjid
JOIN PhotoObjAll AS p ON p.objID = zz.dr8objid
LEFT OUTER JOIN galSpecExtra AS g ON g.specObjID = zz.specobjid
WHERE zz.ra >= {ra_lo} AND zz.ra < {ra_hi}
"""


#: 3,752 rows of ``zoo2MainSpecz`` carry no DR8 object ID (dr8objid = 0, a gap
#: in the original 2013 cross-match, mostly along the celestial equator), so the
#: ID join above loses them. They are recovered by position instead: the
#: nearest DR17 primary photometric object within 3", and its best spectrum.
#: The analysis keeps only matches closer than 1".
GZ2_NODR8_QUERY = """
SELECT zz.dr7objid, n.distance * 60.0 AS sep_arcsec, s.specObjID AS specobjid,
       p.objID AS dr8objid, zz.sample,
       p.ra, p.dec, s.z, s.zErr, s.zWarning, s.class,
       p.petroMag_r, p.extinction_r, p.petroR50_r, p.petroR90_r,
       p.expAB_r, p.deVAB_r, p.fracDeV_r,
       p.dered_u, p.dered_g, p.dered_r, p.dered_i, p.dered_z,
       p.modelMagErr_u, p.modelMagErr_g, p.modelMagErr_r, p.modelMagErr_i,
       g.lgm_tot_p16, g.lgm_tot_p50, g.lgm_tot_p84,
       g.sfr_tot_p50, g.specsfr_tot_p50, g.bptclass,
       zz.t01_smooth_or_features_a01_smooth_weighted_fraction AS w13wf_smooth,
       {w13_columns}
FROM zoo2MainSpecz AS zz
CROSS APPLY dbo.fGetNearestObjEq(zz.ra, zz.dec, 0.05) AS n
JOIN PhotoObjAll AS p ON p.objID = n.objID
JOIN SpecObj AS s ON s.bestObjID = p.objID
LEFT OUTER JOIN galSpecExtra AS g ON g.specObjID = s.specObjID
WHERE zz.dr8objid = 0 AND zz.ra >= {ra_lo} AND zz.ra < {ra_hi}
"""

#: Positional recoveries farther than this are rejected (the DR8 deblender
#: sometimes split or merged objects differently from DR7's).
NODR8_MAX_SEP_ARCSEC = 1.0


def _w13_columns() -> str:
    full = {short: long for long, short in ANSWERS.items()}
    return ",\n       ".join(f"zz.{full[a]}_debiased AS w13_{a}" for a in W13_ANSWERS)


def fetch_galaxy_zoo(*, quiet: bool = False) -> None:
    if not quiet:
        print("Galaxy Zoo 2, Hart et al. (2016) debiasing")
    download(HART16_URL, HART16_PATH, quiet=quiet)
    if not quiet:
        print(f"GZ2 x SDSS {DATA_RELEASE} on SkyServer, {GZ2_RA_CHUNK_DEG:g}-degree strips")
    _fetch_strips(GZ2_QUERY, GZ2_SDSS_DIR, "gz2sdss", GZ2_RA_CHUNK_DEG,
                  {"w13_columns": _w13_columns()}, quiet)
    if not quiet:
        print("GZ2 galaxies without a DR8 ID: positional match on SkyServer")
    _fetch_strips(GZ2_NODR8_QUERY, GZ2_SDSS_DIR, "gz2nodr8", GZ2_RA_CHUNK_DEG,
                  {"w13_columns": _w13_columns()}, quiet)


# ---------------------------------------------------------------------------
# 3. The density tracers: project D2's Main Galaxy Sample query, unchanged
# ---------------------------------------------------------------------------

MGS_QUERY = """
SELECT specObjID, ra, dec, z, zErr,
       petroMag_r, extinction_r,
       dered_u, dered_g, dered_r,
       modelMagErr_u, modelMagErr_g, modelMagErr_r
FROM SpecPhoto
WHERE (legacy_target1 & {bits}) != 0
  AND class = 'GALAXY'
  AND sciencePrimary = 1
  AND mode = 1
  AND zWarning = 0
  AND z BETWEEN {zmin} AND {zmax}
  AND petroMag_r - extinction_r < {rlim}
  AND ra >= {ra_lo} AND ra < {ra_hi}
"""


def fetch_mgs(*, quiet: bool = False) -> None:
    if not quiet:
        print(f"SDSS {DATA_RELEASE} Main Galaxy Sample (density tracers, as in project D2)")
    _fetch_strips(MGS_QUERY, MGS_DIR, "mgs", MGS_RA_CHUNK_DEG,
                  {"bits": MGS_TARGET_BITS, "zmin": MGS_FETCH_Z[0], "zmax": MGS_FETCH_Z[1],
                   "rlim": MGS_R_LIMIT}, quiet)


# ---------------------------------------------------------------------------
# 4-5. Experts, and the K-correction code
# ---------------------------------------------------------------------------

NA10_COLUMNS = ["SDSS", "_RA", "_DE", "zs", "TT", "Bar", "Ring", "Pair", "Int", "log(M)",
                "g'-r'", "b/a", "Rp90", "DR7"]


def fetch_nair_abraham(*, quiet: bool = False) -> None:
    if NA10_PATH.exists():
        return
    from astroquery.vizier import Vizier

    if not quiet:
        print("Nair & Abraham (2010) expert classifications, VizieR " + NA10_CATALOG)
    v = Vizier(columns=NA10_COLUMNS, row_limit=-1)
    table = v.get_catalogs(NA10_CATALOG)[0]
    df = table.to_pandas()
    tmp = NA10_PATH.with_suffix(".part")
    df.to_csv(tmp, index=False)
    tmp.replace(NA10_PATH)
    if not quiet:
        print(f"    {len(df):,d} galaxies")


def fetch_kcor_code(*, quiet: bool = False) -> None:
    try:
        download(KCOR_CODE_URL, KCOR_PATH, quiet=quiet)
    except RuntimeError as exc:  # only the tests use it; never fatal
        print(f"    (K-correction reference code unavailable: {exc})")


# ---------------------------------------------------------------------------
# Building the catalogues
# ---------------------------------------------------------------------------

def _read_strips(folder: Path, prefix: str) -> pd.DataFrame:
    frames = []
    for path in sorted(folder.glob(f"{prefix}_ra*.csv")):
        frame = _parse_skyserver(path.read_text())
        if len(frame):
            frames.append(frame)
    if not frames:
        raise FileNotFoundError(f"no {prefix} strips in {folder} -- run `python fetch_data.py`")
    return pd.concat(frames, ignore_index=True)


def build_catalog(*, quiet: bool = False) -> dict:
    """Merge the Hart et al. votes with the SDSS data into one typed catalogue."""
    expected = 2 * len(_strips(GZ2_RA_CHUNK_DEG))
    have = len(list(GZ2_SDSS_DIR.glob("gz2*_ra*.csv")))
    if have < expected:
        raise FileNotFoundError(f"{have} of {expected} GZ2 strips cached -- run fetch_data.py")

    usecols = ["dr7objid", "ra", "dec", "sample", "gz2_class", "total_classifications",
               "total_votes"]
    for long in ANSWERS:
        usecols += [f"{long}_count", f"{long}_weighted_fraction", f"{long}_debiased"]
    hart = pd.read_csv(HART16_PATH, usecols=usecols)
    n_hart = len(hart)

    sdss = _read_strips(GZ2_SDSS_DIR, "gz2sdss")
    n_sdss_raw = len(sdss)
    nodr8 = _read_strips(GZ2_SDSS_DIR, "gz2nodr8")
    n_nodr8 = len(nodr8)
    nodr8 = nodr8[(nodr8["sep_arcsec"] < NODR8_MAX_SEP_ARCSEC)
                  & ~nodr8["dr7objid"].isin(sdss["dr7objid"])]
    nodr8 = nodr8.sort_values("sep_arcsec").drop_duplicates("dr7objid")
    sdss = pd.concat([sdss, nodr8.drop(columns="sep_arcsec")], ignore_index=True)
    sdss = sdss.drop_duplicates("dr7objid")
    sdss = sdss.rename(columns={"ra": "ra_sdss", "dec": "dec_sdss", "sample": "sample_w13"})
    df = hart.merge(sdss, on="dr7objid", how="left", validate="one_to_one")
    matched = df["specobjid"].notna().to_numpy()

    arrays = {
        "dr7objid": df["dr7objid"].to_numpy(np.int64),
        "ra": df["ra"].to_numpy(np.float64),
        "dec": df["dec"].to_numpy(np.float64),
        "sample": np.array(df["sample"].astype(str).tolist(), dtype="U16"),
        "gz2_class": np.array(df["gz2_class"].astype(str).tolist(), dtype="U16"),
        "n_class": df["total_classifications"].to_numpy(np.int32),
        "n_votes": df["total_votes"].to_numpy(np.int32),
        "in_sdss": matched,
    }
    for long, short in ANSWERS.items():
        arrays[f"n_{short}"] = df[f"{long}_count"].fillna(0).to_numpy(np.int32)
        arrays[f"wf_{short}"] = df[f"{long}_weighted_fraction"].to_numpy(np.float32)
        arrays[f"p_{short}"] = df[f"{long}_debiased"].to_numpy(np.float32)
    float_cols = ["ra_sdss", "dec_sdss", "z", "zErr", "petroMag_r", "extinction_r", "petroR50_r",
                  "petroR90_r", "expAB_r", "deVAB_r", "fracDeV_r", "dered_u", "dered_g",
                  "dered_r", "dered_i", "dered_z", "modelMagErr_u", "modelMagErr_g",
                  "modelMagErr_r", "modelMagErr_i", "lgm_tot_p16", "lgm_tot_p50", "lgm_tot_p84",
                  "sfr_tot_p50", "specsfr_tot_p50", "bptclass", "w13wf_smooth"]
    float_cols += [f"w13_{a}" for a in W13_ANSWERS]
    for col in float_cols:
        a = df[col].to_numpy(np.float64)
        a[a <= -999] = np.nan  # SkyServer's -9999 / -999 "no measurement"
        arrays[col] = a
    arrays["zWarning"] = df["zWarning"].fillna(-1).to_numpy(np.int32)
    arrays["specobjid"] = df["specobjid"].fillna(-1).to_numpy(np.int64)
    arrays["dr8objid"] = df["dr8objid"].fillna(-1).to_numpy(np.int64)

    meta = {
        "sources": {"votes": HART16_URL, "sdss": SKYSERVER_URL + " (zoo2MainSpecz join)"},
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_hart16": int(n_hart),
        "n_sdss_rows": int(n_sdss_raw),
        "n_sdss_unique": int(len(sdss)),
        "n_nodr8_served": int(n_nodr8),
        "n_nodr8_recovered": int(len(nodr8)),
        "n_matched": int(matched.sum()),
    }
    np.savez_compressed(CATALOG_PATH, meta=json.dumps(meta), **arrays)
    if not quiet:
        print(f"  catalogue: {n_hart:,d} GZ2 galaxies, {matched.sum():,d} matched to "
              f"SDSS {DATA_RELEASE} ({len(sdss):,d} unique rows served)")
    return meta


MGS_COLUMNS = {"specObjID": np.int64, "ra": np.float64, "dec": np.float64, "z": np.float64,
               "zErr": np.float32, "petroMag_r": np.float32, "extinction_r": np.float32,
               "dered_u": np.float32, "dered_g": np.float32, "dered_r": np.float32,
               "modelMagErr_u": np.float32, "modelMagErr_g": np.float32,
               "modelMagErr_r": np.float32}


def build_mgs(*, quiet: bool = False) -> dict:
    """The Main Galaxy Sample, cleaned exactly as project D2 cleaned it."""
    df = _read_strips(MGS_DIR, "mgs")
    n_raw = len(df)
    df = df.drop_duplicates("specObjID")
    mags = ["petroMag_r", "dered_u", "dered_g", "dered_r"]
    bad = (df[mags] < -1000).any(axis=1) | df[mags].isna().any(axis=1)
    df = df[~bad].sort_values(["ra", "dec"]).reset_index(drop=True)
    arrays = {c: df[c].to_numpy(dtype=t) for c, t in MGS_COLUMNS.items()}
    meta = {"source": SKYSERVER_URL + " (project D2's query)", "n_rows_served": int(n_raw),
            "n_missing_photometry_dropped": int(bad.sum()), "n_galaxies": int(len(df)),
            "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    np.savez_compressed(MGS_PATH, meta=json.dumps(meta), **arrays)
    if not quiet:
        print(f"  tracers: {len(df):,d} Main Galaxy Sample galaxies "
              f"({int(bad.sum())} without photometry dropped)")
    return meta


def load_catalog() -> tuple[dict, dict]:
    with np.load(CATALOG_PATH, allow_pickle=False) as f:
        arrays = {k: f[k] for k in f.files if k != "meta"}
        meta = json.loads(str(f["meta"]))
    return arrays, meta


def load_mgs() -> tuple[dict, dict]:
    with np.load(MGS_PATH, allow_pickle=False) as f:
        arrays = {k: f[k] for k in f.files if k != "meta"}
        meta = json.loads(str(f["meta"]))
    return arrays, meta


def load_nair_abraham() -> pd.DataFrame:
    return pd.read_csv(NA10_PATH)


# ---------------------------------------------------------------------------
# Images: the very pictures the volunteers classified
# ---------------------------------------------------------------------------
#
# The Galaxy Zoo team archived the GZ2 subject images on Zenodo (record
# 3565489): 243,434 JPEGs, 424 x 424 pixels, in one 3.35 GB zip, plus a table
# mapping each DR7 object ID to its file. A zip file keeps an index (the
# "central directory") at its end listing where every member starts, and
# Zenodo serves byte ranges, so we read the 39 MB index once and then fetch
# each wanted image on its own -- a few hundred for the pictures, a few
# thousand for the neural network -- without downloading the rest.

GZ2_IMAGES_URL = "https://zenodo.org/api/records/3565489/files/images_gz2.zip/content"
GZ2_MAPPING_URL = "https://zenodo.org/api/records/3565489/files/gz2_filename_mapping.csv/content"
GZ2_MAPPING_PATH = GZ_DIR / "gz2_filename_mapping.csv"
GZ2_IMAGE_INDEX_PATH = GZ_DIR / "gz2_images_zip_index.npz"
GZ2_IMAGE_DIR = CUTOUT_DIR / "gz2"


def _range_get(url: str, lo: int, hi: int, *, retries: int = 12, session=None) -> bytes:
    """One byte range of ``url``, waiting out the server's rate limit.

    Zenodo allows 133 requests a minute and says so in its headers; when the
    limit is hit it answers 429 with a Retry-After, which is honoured.
    """
    get = (session or requests).get
    for attempt in range(retries):
        try:
            r = get(url, headers={"User-Agent": USER_AGENT, "Range": f"bytes={lo}-{hi}"},
                    timeout=120)
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 30)) + 1.0)
                continue
            if r.status_code != 206:
                raise IOError(f"HTTP {r.status_code} (expected 206 Partial Content)")
            if r.headers.get("X-RateLimit-Remaining") == "0":
                time.sleep(float(r.headers.get("Retry-After", 30)) + 1.0)
            return r.content
        except (requests.RequestException, IOError):
            time.sleep(min(60, 2 * 2**attempt))
    raise RuntimeError(f"range request {lo}-{hi} failed")


def _content_length(url: str, retries: int = 12) -> int:
    """Size of a remote file, waiting out rate limits (a 429's body is not the file)."""
    for attempt in range(retries):
        try:
            r = requests.head(url, headers={"User-Agent": USER_AGENT}, allow_redirects=True,
                              timeout=60)
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 30)) + 1.0)
                continue
            r.raise_for_status()
            return int(r.headers["Content-Length"])
        except requests.RequestException:
            time.sleep(min(60, 2 * 2**attempt))
    raise RuntimeError(f"could not get the size of {url}")


def gz2_image_index(*, quiet: bool = False) -> dict:
    """Byte offsets of every image inside the remote zip, read once and cached."""
    if GZ2_IMAGE_INDEX_PATH.exists():
        with np.load(GZ2_IMAGE_INDEX_PATH) as f:
            return {k: f[k] for k in f.files}
    import struct

    size = _content_length(GZ2_IMAGES_URL)
    tail = _range_get(GZ2_IMAGES_URL, size - 70_000, size - 1)
    i = tail.rfind(b"PK\x05\x06")
    # End-of-central-directory record. The entry count overflows its 16-bit
    # field in this file, so the directory is walked to its end instead.
    _, _, _, _, _, cd_size, cd_off, _ = struct.unpack("<IHHHHIIH", tail[i:i + 22])
    if not quiet:
        print(f"  GZ2 image zip: reading its {cd_size / 1e6:.0f} MB index")
    cd = _range_get(GZ2_IMAGES_URL, cd_off, cd_off + cd_size - 1)
    index = parse_central_directory(cd)
    np.savez(GZ2_IMAGE_INDEX_PATH, **index)
    if not quiet:
        print(f"    {len(asset):,d} images indexed")
    return index


def parse_central_directory(cd: bytes, prefix: str = "images/") -> dict:
    """Walk a zip's central directory: where each ``<prefix><n>.jpg`` member starts.

    Each entry is a fixed 46-byte record (signature PK\\1\\2, compression
    method, sizes, name length, extra and comment lengths, and the offset of
    the member's local header) followed by its name, extra field and comment.
    Other members (directories, macOS resource forks) are skipped.
    """
    import struct

    asset, offset, csize, method = [], [], [], []
    p = 0
    while p + 46 <= len(cd) and cd[p:p + 4] == b"PK\x01\x02":
        fields = struct.unpack("<IHHHHHHIIIHHHHHII", cd[p:p + 46])
        meth, cs, nlen, elen, clen, lho = (fields[4], fields[8], fields[10], fields[11],
                                           fields[12], fields[16])
        name = cd[p + 46:p + 46 + nlen].decode("utf-8", "replace")
        if name.startswith(prefix) and name.endswith(".jpg"):
            asset.append(int(name[len(prefix):-4]))
            offset.append(lho)
            csize.append(cs)
            method.append(meth)
        p += 46 + nlen + elen + clen
    return {"asset_id": np.array(asset, np.int64), "offset": np.array(offset, np.int64),
            "csize": np.array(csize, np.int64), "method": np.array(method, np.int16)}


def extract_member(blob: bytes, csize: int, method: int) -> bytes:
    """A zip member's bytes from a range starting at its local header.

    The local header is 30 bytes plus a name and an extra field whose lengths
    it states; the compressed data follow. Method 8 is DEFLATE, 0 is stored.
    """
    import struct
    import zlib

    _, _, _, _, _, _, _, _, _, nlen, elen = struct.unpack("<IHHHHHIIIHH", blob[:30])
    data = blob[30 + nlen + elen:30 + nlen + elen + csize]
    if method == 8:
        data = zlib.decompress(data, -15)
    elif method != 0:
        raise ValueError(f"unsupported zip compression method {method}")
    return data


def gz2_asset_ids(objids) -> np.ndarray:
    """Zenodo asset ID (file name) of each galaxy's normal-depth image, or -1."""
    download(GZ2_MAPPING_URL, GZ2_MAPPING_PATH, quiet=True)
    m = pd.read_csv(GZ2_MAPPING_PATH)
    m = m[m["sample"].isin(["original", "extra", "stripe82"])].drop_duplicates("objid")
    lookup = pd.Series(m["asset_id"].to_numpy(), index=m["objid"].to_numpy())
    return lookup.reindex(np.asarray(objids, np.int64)).fillna(-1).to_numpy(np.int64)


def gz2_image_path(objid: int) -> Path:
    return GZ2_IMAGE_DIR / f"{int(objid)}.jpg"


def _fetch_zip_member(offset: int, csize: int, method: int, path: Path, session) -> bool:
    import struct
    import zlib

    try:
        blob = _range_get(GZ2_IMAGES_URL, offset, offset + 30 + 512 + csize, session=session)
        data = extract_member(blob, csize, method)
        if not data.startswith(b"\xff\xd8"):
            return False
        tmp = path.with_suffix(".part")
        tmp.write_bytes(data)
        tmp.replace(path)
        return True
    except (RuntimeError, OSError, zlib.error, struct.error):
        return False


def fetch_gz2_images(objids, *, ra=None, dec=None, r90=None, workers: int = 8,
                     quiet: bool = False) -> dict:
    """The GZ2 subject image of each galaxy, cached as data/cutouts/gz2/<objid>.jpg.

    Galaxies whose image is not in the archive fall back to SkyServer's
    ImgCutout with GZ2's own recipe (needs ``ra``, ``dec``, ``r90``).
    """
    import threading

    GZ2_IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    objids = np.asarray(objids, np.int64)
    todo = np.array([not gz2_image_path(o).exists() for o in objids], bool)
    if not todo.any():
        return {"cached": int(len(objids)), "fetched": 0, "fallback": 0, "failed": 0}
    index = gz2_image_index(quiet=quiet)
    where = pd.Series(np.arange(len(index["asset_id"])), index=index["asset_id"])
    assets = gz2_asset_ids(objids[todo])
    pos = where.reindex(assets).to_numpy()
    jobs, fallback = [], []
    for o, a, k in zip(objids[todo], assets, pos):
        if a >= 0 and np.isfinite(k):
            k = int(k)
            jobs.append((int(index["offset"][k]), int(index["csize"][k]),
                         int(index["method"][k]), gz2_image_path(o)))
        else:
            fallback.append(o)
    if not quiet:
        print(f"  GZ2 images: {len(jobs):,d} to fetch from the Zenodo archive, "
              f"{len(fallback)} from SkyServer, {(~todo).sum():,d} cached")
    local = threading.local()

    def work(job):
        if not hasattr(local, "s"):
            local.s = requests.Session()
        return _fetch_zip_member(*job, local.s)

    ok, t0 = 0, time.time()
    with cf.ThreadPoolExecutor(workers) as pool:
        for i, good in enumerate(pool.map(work, jobs), 1):
            ok += bool(good)
            if not quiet and (i % 1000 == 0 or i == len(jobs)):
                print(f"    {i:,d}/{len(jobs):,d} ({i / (time.time() - t0):.1f}/s)")
    n_fb = 0
    if fallback and ra is not None:
        sel = np.isin(objids, fallback)
        res = fetch_cutouts("gz2", objids[sel], np.asarray(ra)[sel], np.asarray(dec)[sel],
                            0.02 * np.asarray(r90)[sel], 424, workers=4, quiet=quiet)
        n_fb = res["fetched"]
    return {"cached": int((~todo).sum()), "fetched": ok, "fallback": n_fb,
            "failed": int(len(jobs) - ok + len(fallback) - n_fb)}


def verify_image_mapping(objids, ra, dec, r90, *, n: int = 10, seed: int = 1) -> dict:
    """Do the archived images show the galaxies they are filed under?

    The archive's README asks users to check. For ``n`` galaxies with an
    archived image, a fresh cutout is drawn from SkyServer with GZ2's own
    recipe (424 pixels at 0.02 x petroR90 arcsec per pixel) and correlated
    with the archived one in grey, band-passed to structure 2-8 pixels across
    (every GZ2 image is a bright blob in the centre, so whole images of
    *different* galaxies correlate too; what is left -- field stars, arms,
    knots -- matches only for the same galaxy). The best correlation over a
    shift of up to +/-4 pixels and a zoom of 0.85-1.15 is taken: the 2009
    images were resampled with DR7's astrometry and scaled with DR7's
    petroR90, which DR8 onwards re-measured. As a control, each archived image
    is also matched, the same way, against the cutout of a different galaxy.
    """
    from PIL import Image
    from scipy.ndimage import gaussian_filter

    out_path = DATA_DIR / "image_mapping_check.json"
    if out_path.exists():
        return json.loads(out_path.read_text())
    objids = np.asarray(objids)
    have = np.array([gz2_image_path(o).exists() for o in objids])
    rng = np.random.default_rng(seed)
    pick = rng.choice(np.flatnonzero(have & np.isfinite(r90)), n, replace=False)
    fetch_cutouts("check", objids[pick], np.asarray(ra)[pick], np.asarray(dec)[pick],
                  0.02 * np.asarray(r90)[pick], 424, workers=2, quiet=True)

    def grey(path):
        a = np.asarray(Image.open(path).convert("L"), float)
        a = gaussian_filter(a, 1.5) - gaussian_filter(a, 8.0)
        return (a - a.mean()) / (a.std() + 1e-9)

    def zoomed(b, f):
        from scipy.ndimage import zoom

        z = zoom(b, f, order=1)
        n, m = b.shape[0], z.shape[0]
        if m >= n:
            k = (m - n) // 2
            return z[k:k + n, k:k + n]
        out = np.zeros_like(b)
        k = (n - m) // 2
        out[k:k + m, k:k + m] = z
        return out

    def best_corr(a, b, max_shift: int = 4) -> float:
        best = -1.0
        inner = slice(40, -40)
        for f in np.linspace(0.85, 1.15, 13):
            bz = zoomed(b, f)
            for dy in range(-max_shift, max_shift + 1):
                for dx in range(-max_shift, max_shift + 1):
                    s = np.roll(np.roll(bz, dy, 0), dx, 1)
                    aa, ss = a[inner, inner], s[inner, inner]
                    r = float(np.mean((aa - aa.mean()) * (ss - ss.mean())) / (aa.std() * ss.std() + 1e-9))
                    best = max(best, r)
        return best

    same, other = [], []
    for k, i in enumerate(pick):
        a = grey(gz2_image_path(objids[i]))
        b = grey(cutout_path("check", objids[i]))
        c = grey(cutout_path("check", objids[pick[(k + 1) % n]]))
        same.append(best_corr(a, b))
        other.append(best_corr(a, c))
    res = {"n": int(n), "dr7objid": [int(objids[i]) for i in pick],
           "corr_same_galaxy": same, "corr_other_galaxy": other}
    out_path.write_text(json.dumps(res, indent=1))
    return res


def cutout_path(kind: str, objid: int) -> Path:
    return CUTOUT_DIR / kind / f"{int(objid)}.jpg"


def _fetch_one_cutout(url: str, params: dict, path: Path, retries: int = 6) -> bool:
    for attempt in range(retries):
        try:
            r = requests.get(url, params=params, headers={"User-Agent": USER_AGENT}, timeout=90)
            r.raise_for_status()
            if not r.content.startswith(b"\xff\xd8"):
                raise RuntimeError("not a JPEG")
            tmp = path.with_suffix(".part")
            tmp.write_bytes(r.content)
            tmp.replace(path)
            return True
        except (requests.RequestException, RuntimeError):
            time.sleep(min(60, 2 * 2**attempt))
    return False


def fetch_cutouts(kind: str, objids, ra, dec, scale_arcsec, pixels: int, *,
                  workers: int = 4, quiet: bool = False) -> dict:
    """SDSS gri JPEG cutouts from SkyServer's ImgCutout (DR17), cached per galaxy.

    Used as a fallback, and to spot-check that the archived images belong to
    the galaxies they are filed under. ``scale_arcsec`` is arcsec per pixel.
    """
    out = CUTOUT_DIR / kind
    out.mkdir(parents=True, exist_ok=True)
    jobs = []
    for oid, a, d, s in zip(objids, ra, dec, np.broadcast_to(scale_arcsec, len(objids))):
        path = cutout_path(kind, oid)
        if not path.exists():
            params = {"ra": f"{a:.6f}", "dec": f"{d:.6f}", "scale": f"{s:.5f}",
                      "width": pixels, "height": pixels}
            jobs.append((params, path))
    n_cached = len(objids) - len(jobs)
    ok = 0
    with cf.ThreadPoolExecutor(workers) as pool:
        futures = [pool.submit(_fetch_one_cutout, IMGCUTOUT_URL, p, path) for p, path in jobs]
        for fut in cf.as_completed(futures):
            ok += bool(fut.result())
    return {"cached": n_cached, "fetched": ok, "failed": len(jobs) - ok}


# ---------------------------------------------------------------------------
# Galaxy Zoo DESI: Zoobot's published predictions
# ---------------------------------------------------------------------------

GZDESI_SUBSET_PATH = GZDESI_DIR / "gzdesi_gz2_matches.parquet"

#: The Zoobot predictions compared with: the first questions of the GZ DESI
#: decision tree, which mirror GZ2's (smooth or featured; edge-on; spiral arms;
#: bar -- split into strong and weak; arm count).
GZDESI_COLUMNS = [
    "smooth-or-featured_smooth_fraction", "smooth-or-featured_featured-or-disk_fraction",
    "smooth-or-featured_artifact_fraction", "disk-edge-on_yes_fraction",
    "has-spiral-arms_yes_fraction", "bar_strong_fraction", "bar_weak_fraction", "bar_no_fraction",
    "spiral-arm-count_1_fraction", "spiral-arm-count_2_fraction", "spiral-arm-count_3_fraction",
    "spiral-arm-count_4_fraction", "spiral-arm-count_more-than-4_fraction",
    "spiral-arm-count_cant-tell_fraction",
]


class RemoteFile(io.RawIOBase):
    """A read-only, seekable view of a remote file, one HTTP range per read.

    Enough of a file for pyarrow, which reads a parquet file's index at the
    end and then only the column chunks it is asked for. Built on
    :func:`_range_get`, so it waits out rate limits instead of failing.
    """

    def __init__(self, url: str):
        super().__init__()
        self.url, self.size, self.pos = url, _content_length(url), 0
        self.session = requests.Session()
        self.bytes_read = 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=io.SEEK_SET):
        self.pos = {io.SEEK_SET: 0, io.SEEK_CUR: self.pos, io.SEEK_END: self.size}[whence] + offset
        return self.pos

    def readinto(self, b):
        n = min(len(b), self.size - self.pos)
        if n <= 0:
            return 0
        data = _range_get(self.url, self.pos, self.pos + n - 1, session=self.session)
        b[:len(data)] = data
        self.pos += len(data)
        self.bytes_read += len(data)
        return len(data)


def fetch_gzdesi_subset(ra, dec, objids, columns: list[str] = GZDESI_COLUMNS, *,
                        max_sep_arcsec: float = 1.5, quiet: bool = False) -> Path:
    """Read only ``columns`` of the Galaxy Zoo DESI catalogue, by HTTP range.

    Parquet stores each column in its own chunk and keeps an index at the end
    of the file, so pyarrow can fetch the coordinates and a few morphology
    columns -- about 40% of the 659 MB -- without the rest. Rows are matched
    to our galaxies by position and only the matches are kept.
    """
    if GZDESI_SUBSET_PATH.exists():
        return GZDESI_SUBSET_PATH
    import pyarrow.parquet as pq
    from scipy.spatial import cKDTree

    if not quiet:
        print("Galaxy Zoo DESI: reading Zoobot predictions by byte range")
    t0 = time.time()
    fh = RemoteFile(GZDESI_URL)
    pf = pq.ParquetFile(fh)
    missing = [c for c in ["ra", "dec"] + columns if c not in pf.schema_arrow.names]
    if missing:
        raise KeyError(f"columns not in the GZ DESI catalogue: {missing}")
    desi = pf.read(columns=["ra", "dec"] + columns).to_pandas()
    if not quiet:
        print(f"    {len(desi):,d} rows, {fh.bytes_read / 1e6:.0f} MB read in {time.time() - t0:.0f}s")

    def unit(r, d):
        r, d = np.radians(np.asarray(r, float)), np.radians(np.asarray(d, float))
        return np.column_stack([np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)])

    tree = cKDTree(unit(desi["ra"], desi["dec"]))
    chord = 2 * np.sin(np.radians(max_sep_arcsec / 3600) / 2)
    dist, idx = tree.query(unit(ra, dec), distance_upper_bound=chord)
    hit = np.isfinite(dist)
    sub = desi.iloc[idx[hit]].reset_index(drop=True)
    sub.insert(0, "dr7objid", np.asarray(objids)[hit])
    sub["sep_arcsec"] = np.degrees(2 * np.arcsin(dist[hit] / 2)) * 3600
    tmp = GZDESI_SUBSET_PATH.with_suffix(".part")
    sub.to_parquet(tmp)
    tmp.replace(GZDESI_SUBSET_PATH)
    if not quiet:
        print(f"    {hit.sum():,d} of {len(hit):,d} GZ2 galaxies matched within {max_sep_arcsec}\"")
    return GZDESI_SUBSET_PATH


# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    fetch_galaxy_zoo(quiet=args.quiet)
    fetch_mgs(quiet=args.quiet)
    fetch_nair_abraham(quiet=args.quiet)
    fetch_kcor_code(quiet=args.quiet)
    build_catalog(quiet=args.quiet)
    build_mgs(quiet=args.quiet)


if __name__ == "__main__":
    main()
