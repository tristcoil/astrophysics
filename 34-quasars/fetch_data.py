"""Download and cache everything the analysis needs.

All sources are free and public:

1. **Wu & Shen (2022)**: continuum and emission-line measurements, black hole
   masses and Eddington ratios for all 750,414 quasars of SDSS DR16Q
   (1.39 GB gzipped FITS from the authors' site). Only the ~60 columns used are
   kept, read in one streaming pass (see ``fitsstream.py``).
2. **DR16Q itself** (Lyke et al. 2020) from VizieR VII/289: redshift sources and
   the target-selection bits that say which survey programme chose each quasar.
3. **The SDSS DR7 uniform quasar sample**: Shen et al. (2011) properties and the
   DR7 quasar catalogue (Schneider et al. 2010), both from VizieR.
4. **Richards et al. (2006)**: selection function, K-corrections, statistical
   sample and binned luminosity function; **Shen & Kelly (2012)**: binned DR7
   luminosity and mass functions. VizieR.
5. **Quasars beyond z = 5.7**: Farina et al. (2022), Shen et al. (2019),
   Yang et al. (2021), Banados et al. (2016). VizieR.
6. **Vanden Berk et al. (2001)**: the composite quasar spectrum (arXiv source).
7. **SDSS filter curves** and a patch of the **stellar locus** (SkyServer).
8. **~4,600 SDSS spectra** chosen evenly in redshift, for our own composite.

Everything is cached under ``data/`` as served (VizieR TSV files per sky strip,
the FITS catalogue, spectra), so an interrupted run resumes where it stopped.

    python fetch_data.py                # everything
    python fetch_data.py --no-spectra   # catalogues only
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import io
import re
import sys
import tarfile
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

import fitsstream
from config import (
    BANDS, BROWSER_AGENT, COMPOSITE_MAX_BAL_PROB, COMPOSITE_MIN_SN,
    COMPOSITE_N_ZBINS, COMPOSITE_PER_BIN, COMPOSITE_SEED, COMPOSITE_Z_RANGE, DERIVED,
    DR16Q_COLUMNS, DR16Q_DIR, DR16Q_N_EXPECTED, DR16Q_SOURCE, DR16Q_TABLE, DR7_DIR, DR7_TABLE,
    DR7Q_COLUMNS, DR7Q_SOURCE, FILTER_DIR, HIGHZ_TABLES, LIT_DIR, RICHARDS06_TABLES,
    SDSS_FILTER_URL, SHEN11_COLUMNS, SHEN11_N_EXPECTED, SHEN11_SOURCE, SHENKELLY12_TABLES,
    SKYSERVER_URL, SPEC_DIR, SPEC_LIST, SPEC_URL, STAR_QUERY, STARS_CSV, USER_AGENT,
    VB01_ARXIV_SRC, VB01_TAR, VIZIER_ASU, WUSHEN_COLUMNS_EXT1, WUSHEN_COLUMNS_EXT2, WUSHEN_GZ,
    WUSHEN_SIZE, WUSHEN_TABLE, WUSHEN_URL,
)

#: (connect, read) timeouts. The connect timeout is short so a host whose IPv6
#: address is unreachable from here falls back to IPv4 quickly.
TIMEOUT = (15.0, 600.0)
_local = threading.local()


def session() -> requests.Session:
    """One keep-alive session per thread (spectra are fetched in parallel)."""
    s = getattr(_local, "session", None)
    if s is None:
        s = requests.Session()
        s.headers["User-Agent"] = USER_AGENT
        _local.session = s
    return s


def get(url: str, *, params=None, headers=None, stream=False, retries: int = 6,
        timeout=TIMEOUT) -> requests.Response:
    """GET with retries on connection errors and on 429/5xx answers."""
    last = None
    for attempt in range(retries):
        try:
            r = session().get(url, params=params, headers=headers, stream=stream, timeout=timeout)
            if r.status_code in (200, 206):
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"HTTP {r.status_code}"
                wait = float(r.headers.get("Retry-After", 0) or 0) or 5.0 * (attempt + 1)
                time.sleep(min(wait, 120))
                continue
            r.raise_for_status()
        except (requests.ConnectionError, requests.Timeout,
                requests.exceptions.ChunkedEncodingError) as exc:
            last = exc
            time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"GET {url} failed after {retries} attempts: {last}")


def download(url: str, dest: Path, *, expected_size: int | None = None, headers=None,
             label: str | None = None) -> Path:
    """Resumable download: bytes go to ``dest.part`` and are renamed when complete."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and (expected_size is None or dest.stat().st_size == expected_size):
        return dest
    part = dest.with_name(dest.name + ".part")
    t0, last_report = time.time(), 0.0
    for attempt in range(20):
        pos = part.stat().st_size if part.exists() else 0
        if expected_size is not None and pos == expected_size:
            break
        hdr = dict(headers or {})
        if pos:
            hdr["Range"] = f"bytes={pos}-"
        try:
            r = get(url, headers=hdr, stream=True)
            if pos and r.status_code != 206:
                pos = 0  # server ignored the range: start over
            total = expected_size or (pos + int(r.headers.get("Content-Length", 0)))
            with open(part, "ab" if pos else "wb") as fh:
                for block in r.iter_content(1 << 20):
                    fh.write(block)
                    pos += len(block)
                    if label and time.time() - last_report > 20:
                        last_report = time.time()
                        rate = (pos - 0) / max(time.time() - t0, 1e-3) / 1e6
                        print(f"    {label}: {pos / 1e6:,.0f} / {total / 1e6:,.0f} MB "
                              f"({rate:.1f} MB/s)", flush=True)
            if expected_size is None or pos >= expected_size:
                break
        except (requests.ConnectionError, requests.Timeout,
                requests.exceptions.ChunkedEncodingError, RuntimeError) as exc:
            print(f"    {label or url}: interrupted at {pos / 1e6:.0f} MB ({exc}); resuming",
                  flush=True)
            time.sleep(5)
    size = part.stat().st_size
    if expected_size is not None and size != expected_size:
        raise RuntimeError(f"{dest.name}: got {size} bytes, expected {expected_size}")
    part.rename(dest)
    return dest


# ---------------------------------------------------------------------------
# VizieR
# ---------------------------------------------------------------------------

def parse_vizier_tsv(text: str) -> pd.DataFrame:
    """Parse a VizieR ASU tab-separated answer (header, units row, dashes row, data)."""
    if "QUERY_STATUS=OVERFLOW" in text or "#Truncated" in text:
        raise RuntimeError("VizieR truncated the answer; use smaller sky strips")
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    i = 0
    while i < len(lines) and not lines[i].strip():
        i += 1
    if i >= len(lines):
        return pd.DataFrame()
    header = lines[i].split("\t")
    # The units row can be blank (dimensionless tables), so anchor on the row of
    # dashes that always separates the column header from the data.
    j = i + 1
    while j < len(lines) and not re.fullmatch(r"[\t ]*-[-\t ]*", lines[j]):
        j += 1
    if j >= len(lines):
        raise RuntimeError("VizieR answer has no header separator")
    body = "\n".join(ln for ln in lines[j + 1:] if ln.strip())
    df = pd.read_csv(io.StringIO(body), sep="\t", header=None, names=header, dtype=str,
                     keep_default_na=False, quoting=3)
    return df.apply(lambda col: col.str.strip())


def vizier_query(source: str, columns: list[str] | None, constraints: dict | None = None) -> str:
    params = {"-source": source, "-out.max": "unlimited"}
    if columns:
        params["-out"] = ",".join(columns)
    else:
        params["-out.all"] = ""
    params.update(constraints or {})
    r = get(VIZIER_ASU, params=params, headers={"User-Agent": BROWSER_AGENT})
    r.encoding = "utf-8"
    return r.text


def fetch_vizier(source: str, dest_dir: Path, *, columns: list[str] | None = None,
                 ra_strips: int = 1, name: str | None = None) -> pd.DataFrame:
    """Fetch a VizieR table, optionally in RA strips, caching each answer as served."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    stem = name or source.replace("/", "_")
    frames = []
    edges = np.linspace(0.0, 360.0, ra_strips + 1)
    for k in range(ra_strips):
        path = dest_dir / (f"{stem}.tsv" if ra_strips == 1 else f"{stem}_ra{k:02d}.tsv")
        if not path.exists():
            cons = {}
            if ra_strips > 1:
                ra_col = "RAJ2000" if columns is None or "RAJ2000" in columns else "_RA"
                cons[ra_col] = f"{edges[k]:.6f}..{edges[k + 1]:.6f}"
            text = vizier_query(source, columns, cons)
            if not text.lstrip().startswith("#"):
                raise RuntimeError(f"unexpected VizieR answer for {source}: {text[:200]!r}")
            tmp = path.with_suffix(".tmp")
            tmp.write_text(text, encoding="utf-8")
            tmp.rename(path)
        frames.append(parse_vizier_tsv(path.read_text(encoding="utf-8")))
    df = pd.concat(frames, ignore_index=True)
    if "recno" in df.columns:  # strips share their boundary values: drop repeats
        df = df.drop_duplicates("recno").reset_index(drop=True)
    return df


def _numeric(df: pd.DataFrame, skip: tuple[str, ...] = ()) -> pd.DataFrame:
    """Convert every column that parses as numbers; leave text columns alone."""
    out = df.copy()
    for c in out.columns:
        if c in skip:
            continue
        conv = pd.to_numeric(out[c].replace("", np.nan), errors="coerce")
        if conv.notna().sum() >= 0.5 * (out[c] != "").sum() and (out[c] != "").any():
            out[c] = conv
    return out


def _flags(series: pd.Series) -> np.ndarray:
    """Target-bit columns hold up to 64-bit masks: parse exactly, as unsigned.

    A few entries are -1, DR16Q's marker for "this survey's flags do not
    apply" (e.g. EBOSS_TARGET0 outside SEQUELS); they carry no bits, like 0.
    """
    vals = [int(v) if v not in ("", None) else 0 for v in series]
    return np.array([v if v >= 0 else 0 for v in vals], dtype=np.uint64)


def fetch_dr16q() -> pd.DataFrame:
    if DR16Q_TABLE.exists():
        return pd.read_parquet(DR16Q_TABLE)
    print("  DR16Q columns from VizieR VII/289 (36 RA strips)...", flush=True)
    raw = fetch_vizier(DR16Q_SOURCE, DR16Q_DIR, columns=DR16Q_COLUMNS, ra_strips=36,
                       name="dr16q")
    flag_cols = ["Btarg1", "eBOSST0", "eBOSST1", "eBOSST2", "Atarg1", "Atarg2"]
    df = _numeric(raw, skip=("SDSS", "r_z", *flag_cols))
    for c in flag_cols:
        df[c] = _flags(raw[c])
    if len(df) != DR16Q_N_EXPECTED:
        raise RuntimeError(f"DR16Q: {len(df)} rows, expected {DR16Q_N_EXPECTED}")
    DERIVED.mkdir(parents=True, exist_ok=True)
    df.to_parquet(DR16Q_TABLE)
    print(f"    {len(df):,} quasars", flush=True)
    return df


def fetch_dr7() -> pd.DataFrame:
    """Shen et al. 2011 joined to the DR7 quasar catalogue (same objects, same order)."""
    if DR7_TABLE.exists():
        return pd.read_parquet(DR7_TABLE)
    print("  DR7 quasar properties (Shen+2011) and DR7Q from VizieR...", flush=True)
    s11 = _numeric(fetch_vizier(SHEN11_SOURCE, DR7_DIR, columns=SHEN11_COLUMNS, ra_strips=4,
                                name="shen11"), skip=("SDSS",))
    q7 = _numeric(fetch_vizier(DR7Q_SOURCE, DR7_DIR, columns=DR7Q_COLUMNS, ra_strips=4,
                               name="dr7q"), skip=("SDSS",))
    for name, df, n in (("Shen+2011", s11, SHEN11_N_EXPECTED), ("DR7Q", q7, SHEN11_N_EXPECTED)):
        if len(df) != n:
            raise RuntimeError(f"{name}: {len(df)} rows, expected {n}")
        if df["SDSS"].duplicated().any():
            raise RuntimeError(f"{name}: duplicated SDSS names")
    q7 = q7.drop(columns=["recno"]).rename(columns={c: f"dr7q_{c}" for c in q7.columns
                                                    if c not in ("SDSS",)})
    df = s11.merge(q7, on="SDSS", how="left", validate="one_to_one")
    if df["dr7q_z"].isna().any():
        raise RuntimeError("Shen+2011 objects missing from DR7Q")
    DERIVED.mkdir(parents=True, exist_ok=True)
    df.to_parquet(DR7_TABLE)
    print(f"    {len(df):,} DR7 quasars, {int((df['Uni'] == 1).sum()):,} with uniform flag 1",
          flush=True)
    return df


def fetch_small_vizier_tables() -> None:
    for group in (RICHARDS06_TABLES, SHENKELLY12_TABLES, HIGHZ_TABLES):
        for key, source in group.items():
            fetch_vizier(source, LIT_DIR, columns=None, name=source.replace("/", "_"))
    print("  literature tables cached (Richards+06, Shen & Kelly 12, z > 5.7 samples)", flush=True)


def load_vizier_table(source: str) -> pd.DataFrame:
    """A cached small VizieR table, numbers converted."""
    df = fetch_vizier(source, LIT_DIR, columns=None, name=source.replace("/", "_"))
    return _numeric(df)


# ---------------------------------------------------------------------------
# Wu & Shen 2022
# ---------------------------------------------------------------------------

def fetch_wushen() -> pd.DataFrame:
    if WUSHEN_TABLE.exists():
        return pd.read_parquet(WUSHEN_TABLE)
    print(f"  Wu & Shen 2022 catalogue ({WUSHEN_SIZE / 1e9:.2f} GB)...", flush=True)
    download(WUSHEN_URL, WUSHEN_GZ, expected_size=WUSHEN_SIZE, label="dr16q_prop")
    print("    extracting columns in one streaming pass (~5 min)...", flush=True)
    t0 = time.time()

    def progress(hdu, done, total):
        if done % 200000 < 20000 or done == total:
            print(f"      HDU {hdu}: {done:,}/{total:,} rows ({time.time() - t0:.0f} s)",
                  flush=True)

    cols = fitsstream.read_columns(WUSHEN_GZ, {1: WUSHEN_COLUMNS_EXT1, 2: WUSHEN_COLUMNS_EXT2},
                                   progress=progress)
    flat = {}
    for hdu in (1, 2):
        for name, arr in cols[hdu].items():
            if arr.ndim == 1:
                flat[name] = arr
            else:
                for j in range(arr.shape[1]):
                    flat[f"{name}_{j}"] = arr[:, j]
    df = pd.DataFrame(flat)
    DERIVED.mkdir(parents=True, exist_ok=True)
    df.to_parquet(WUSHEN_TABLE)
    print(f"    {len(df):,} rows x {df.shape[1]} columns -> {WUSHEN_TABLE.name}", flush=True)
    return df


# ---------------------------------------------------------------------------
# Small files: Vanden Berk composite, filters, stars
# ---------------------------------------------------------------------------

def fetch_vandenberk() -> Path:
    download(VB01_ARXIV_SRC, VB01_TAR)
    with tarfile.open(VB01_TAR) as tar:
        names = tar.getnames()
        for want in ("VandenBerk.tab1.dat", "VandenBerk.tab2.tex"):
            if want not in names:
                raise RuntimeError(f"{want} missing from the Vanden Berk arXiv source")
    return VB01_TAR


def read_vandenberk() -> tuple[pd.DataFrame, str]:
    """The median composite (wavelength, flux, 68% semi-interquantile range) and line table."""
    fetch_vandenberk()
    with tarfile.open(VB01_TAR) as tar:
        dat = tar.extractfile("VandenBerk.tab1.dat").read().decode("latin-1")
        lines = tar.extractfile("VandenBerk.tab2.tex").read().decode("latin-1")
    rows = []
    for ln in dat.splitlines():
        parts = ln.split()
        if len(parts) == 3:
            try:
                rows.append([float(p) for p in parts])
            except ValueError:
                continue
    comp = pd.DataFrame(rows, columns=["wave", "flux", "err"])
    return comp, lines


def fetch_filters() -> None:
    for band in BANDS:
        download(SDSS_FILTER_URL.format(band=band), FILTER_DIR / f"{band}.dat")


def fetch_stars() -> pd.DataFrame:
    if not STARS_CSV.exists():
        STARS_CSV.parent.mkdir(parents=True, exist_ok=True)
        r = get(SKYSERVER_URL, params={"cmd": " ".join(STAR_QUERY.split()), "format": "csv"})
        text = r.text
        if not text.startswith("#Table1"):
            raise RuntimeError(f"SkyServer error: {text[:300]!r}")
        STARS_CSV.write_text(text)
    return pd.read_csv(STARS_CSV, skiprows=1)


# ---------------------------------------------------------------------------
# Spectra for the composite
# ---------------------------------------------------------------------------

def select_composite_sample(ws: pd.DataFrame) -> pd.DataFrame:
    """Choose quasars evenly in log(1+z) among clean BOSS spectra.

    Clean means: a BOSS/eBOSS spectrum (3600-10400 A coverage), a systemic
    redshift, median S/N >= 4 per pixel, and a low probability of broad
    absorption lines (their troughs would dent the stack). Within each redshift
    bin the choice is random with a fixed seed, so it is reproducible and not
    biased toward bright objects beyond the S/N cut.
    """
    z = ws["Z_SYS"].to_numpy()
    good = ((ws["IF_BOSS_SDSS"] == "BOSS") & (z > COMPOSITE_Z_RANGE[0]) & (z < COMPOSITE_Z_RANGE[1])
            & (ws["SN_MEDIAN_ALL"] >= COMPOSITE_MIN_SN) & (ws["BAL_PROB"] < COMPOSITE_MAX_BAL_PROB)
            & (ws["Z_SYS_ERR"] > 0) & (ws["Z_SYS_ERR"] < 0.01))
    cand = ws.loc[good, ["SDSS_NAME", "PLATE", "MJD", "FIBERID", "Z_SYS", "Z_SYS_ERR", "EBV",
                         "SN_MEDIAN_ALL", "M_I"]].copy()
    lo, hi = np.log1p(COMPOSITE_Z_RANGE[0]), np.log1p(COMPOSITE_Z_RANGE[1])
    edges = np.linspace(lo, hi, COMPOSITE_N_ZBINS + 1)
    cand["zbin"] = np.clip(np.digitize(np.log1p(cand["Z_SYS"]), edges) - 1, 0,
                           COMPOSITE_N_ZBINS - 1)
    rng = np.random.default_rng(COMPOSITE_SEED)
    picks = []
    for b, grp in cand.groupby("zbin"):
        n = min(COMPOSITE_PER_BIN, len(grp))
        picks.append(grp.iloc[np.sort(rng.choice(len(grp), n, replace=False))])
    sample = pd.concat(picks).sort_values("Z_SYS").reset_index(drop=True)
    sample["n_candidates_in_bin"] = sample["zbin"].map(cand["zbin"].value_counts())
    return sample


def spec_path(plate: int, mjd: int, fiber: int) -> Path:
    return SPEC_DIR / f"{plate:04d}" / f"spec-{plate:04d}-{mjd:05d}-{fiber:04d}.fits"


def _fetch_one_spectrum(row) -> tuple[str, bool]:
    plate, mjd, fiber = int(row.PLATE), int(row.MJD), int(row.FIBERID)
    path = spec_path(plate, mjd, fiber)
    if path.exists() and path.stat().st_size > 20000:
        return path.name, True
    url = SPEC_URL.format(plate=plate, mjd=mjd, fiber=fiber)
    try:
        download(url, path)
        with open(path, "rb") as fh:
            ok = fh.read(6) == b"SIMPLE"
        if not ok:
            path.unlink()
        return path.name, ok
    except Exception as exc:  # report, keep going; a rerun retries the failures
        print(f"    {path.name}: {exc}", flush=True)
        return path.name, False


def fetch_spectra(ws: pd.DataFrame, workers: int = 12) -> pd.DataFrame:
    if SPEC_LIST.exists():
        sample = pd.read_parquet(SPEC_LIST)
    else:
        sample = select_composite_sample(ws)
        DERIVED.mkdir(parents=True, exist_ok=True)
        sample.to_parquet(SPEC_LIST)
    todo = [r for r in sample.itertuples()
            if not (spec_path(int(r.PLATE), int(r.MJD), int(r.FIBERID)).exists())]
    print(f"  composite sample: {len(sample):,} spectra in {sample['zbin'].nunique()} redshift "
          f"bins; {len(todo):,} to download (~0.2 MB each)", flush=True)
    t0, done, failed = time.time(), 0, 0
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        for name, ok in pool.map(_fetch_one_spectrum, todo):
            done += 1
            failed += (not ok)
            if done % 200 == 0 or done == len(todo):
                rate = done / max(time.time() - t0, 1e-3)
                print(f"    {done:,}/{len(todo):,} spectra ({rate:.1f}/s, {failed} failed)",
                      flush=True)
    return sample


# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-spectra", action="store_true", help="skip the composite's spectra")
    args = ap.parse_args(argv)
    t0 = time.time()
    print("Fetching I3 data", flush=True)
    fetch_small_vizier_tables()
    fetch_vandenberk()
    fetch_filters()
    fetch_stars()
    print("  Vanden Berk composite, SDSS filters and stellar-locus stars cached", flush=True)
    fetch_dr16q()
    fetch_dr7()
    ws = fetch_wushen()
    if not args.no_spectra:
        fetch_spectra(ws)
    print(f"Done in {(time.time() - t0) / 60:.1f} min", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
