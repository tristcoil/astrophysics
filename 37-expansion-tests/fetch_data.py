"""Download and cache everything the three tests need, and load it back.

1. **DES-SN5YR** (GitHub ``des-science/DES-SN5YR``, tag 1.3): photometry of
   every DES transient, the SALT3 fits of the 1,635 cosmology supernovae, and
   the Hubble-diagram metadata (peak dates, Type Ia probabilities). ~80 MB.
2. **White et al. (2024)'s light-curve widths** (GitHub
   ``ryanwhite1/DES-Time-Dilation``), used only as an answer key.
3. **Paper sources from arXiv**: the TeX tables of Lubin & Sandage's Tolman
   test (Papers I, III, IV), and the papers whose CMB-temperature values make
   up ``literature/tcmb_measurements.csv`` (checked by test_tcmb_table.py).

Everything lands in ``data/`` with a manifest recording each file's URL, size,
SHA-256 and retrieval time. Nothing already cached is downloaded again, so once
the cache is warm the project runs offline.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import pickle
import shutil
import sys
import tarfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from config import (ARXIV_SOURCES, DATA_DIR, DES_DIR, DES_FILES, DES_RAW, LAMDA_CO_URL,
                    MANIFEST_PATH, PAPER_DIR, WHITE_DIR, WHITE_FILES, WHITE_RAW)

#: arXiv asks automated clients to identify themselves and to pause between
#: requests. (A browser-like agent gets bot checks on some archives.)
USER_AGENT = "astro-hobby-project/1.0 (37-expansion-tests; python-urllib)"
ARXIV_PAUSE_S = 3.0


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {}


def _download(url: str, dest: Path, *, manifest: dict, refresh: bool = False,
              retries: int = 4) -> bool:
    """Fetch one URL into ``dest`` unless cached; record provenance. True if fetched."""
    key = str(dest.relative_to(DATA_DIR))
    if dest.exists() and not refresh:
        print(f"  cached   {key}")
        return False
    print(f"  fetching {key}  <-  {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=180) as response, open(tmp, "wb") as out:
                shutil.copyfileobj(response, out)
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if attempt == retries - 1:
                raise
            wait = 5 * (attempt + 1)
            print(f"           {err}; retrying in {wait} s")
            time.sleep(wait)
    head = tmp.read_bytes()[:200].lower()
    if b"<!doctype html" in head or b"<html" in head:
        tmp.unlink()
        raise RuntimeError(f"{url} returned an HTML page instead of data")
    tmp.replace(dest)  # atomic: a half-written file is never mistaken for a cached one
    manifest[key] = {
        "url": url,
        "bytes": dest.stat().st_size,
        "sha256": _sha256(dest),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    print(f"           {manifest[key]['bytes'] / 1e6:.2f} MB")
    return True


def _unpack_arxiv(src: Path, outdir: Path) -> None:
    """arXiv e-prints are a tarball, a single gzipped file, or a bare PDF."""
    outdir.mkdir(exist_ok=True)
    raw = src.read_bytes()
    if raw[:4] == b"%PDF":
        (outdir / "paper.pdf").write_bytes(raw)
        return
    data = gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw
    try:
        with tarfile.open(fileobj=io.BytesIO(data)) as tar:
            for member in tar.getmembers():
                # Only plain files, and never outside the target directory.
                name = Path(member.name)
                if not member.isfile() or name.is_absolute() or ".." in name.parts:
                    continue
                target = outdir / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(tar.extractfile(member).read())
    except tarfile.ReadError:
        (outdir / ("paper.pdf" if data[:4] == b"%PDF" else "main.tex")).write_bytes(data)


def fetch_all(*, refresh: bool = False) -> dict:
    manifest = _load_manifest()

    print(f"DES-SN5YR (tag of the release used by White et al. 2024)")
    for local, repo_path in DES_FILES.items():
        _download(DES_RAW + repo_path, DES_DIR / local, manifest=manifest, refresh=refresh)

    print("White et al. (2024) widths (answer key only)")
    for local, repo_path in WHITE_FILES.items():
        _download(WHITE_RAW + repo_path, WHITE_DIR / local, manifest=manifest, refresh=refresh)

    print("Paper sources from arXiv")
    for key, arxiv_id in ARXIV_SOURCES.items():
        dest = PAPER_DIR / f"{key}.eprint"
        if _download(f"https://arxiv.org/e-print/{arxiv_id}", dest, manifest=manifest,
                     refresh=refresh):
            time.sleep(ARXIV_PAUSE_S)
        if not (PAPER_DIR / key).exists() or refresh:
            _unpack_arxiv(dest, PAPER_DIR / key)

    print("Molecular data (LAMDA)")
    _download(LAMDA_CO_URL, PAPER_DIR / "lamda_co.dat", manifest=manifest, refresh=refresh)

    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


# ---------------------------------------------------------------------------
# Loaders -- read only from the cache, never from the network
# ---------------------------------------------------------------------------


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


def _decode(col):
    """FITS columns are big-endian bytes; pandas wants native numbers and str."""
    if col.dtype.kind in "SU":
        return np.char.strip(col.astype(str))
    return col.astype(col.dtype.newbyteorder("="))


def load_des_head():
    """One row per DES transient (19,706), as a pandas DataFrame."""
    import pandas as pd
    from astropy.io import fits

    with fits.open(_require(DES_DIR / "DES-SN5YR_DES_HEAD.FITS.gz")) as hdul:
        data = hdul[1].data
        cols = {name: _decode(np.asarray(data[name])) for name in data.columns.names
                if np.asarray(data[name]).ndim == 1}
    head = pd.DataFrame(cols)
    head["SNID"] = head["SNID"].astype(str).str.strip().astype(np.int64)
    return head


def load_des_phot():
    """Every DES flux measurement, as a DataFrame (pointers from the HEAD file).

    Rows with MJD = -777 are the end-of-light-curve markers SNANA writes; they
    are kept so that the HEAD file's 1-based pointers stay valid.
    """
    import pandas as pd
    from astropy.io import fits

    with fits.open(_require(DES_DIR / "DES-SN5YR_DES_PHOT.FITS.gz")) as hdul:
        data = hdul[1].data
        keep = ("MJD", "BAND", "FIELD", "PHOTFLAG", "PHOTPROB", "FLUXCAL", "FLUXCALERR")
        cols = {name: _decode(np.asarray(data[name])) for name in keep}
    phot = pd.DataFrame(cols)
    phot["BAND"] = phot["BAND"].str.replace("DES-", "", regex=False)
    return phot


def load_lcplot():
    """The SALT3 fits of the 1,635 cosmology supernovae: data and 2-day model.

    ``DATA_MODEL`` is 1 for data used in the fit, -1 for data excluded from it
    (outside the fit's phase window), and 0 for the model curve. ``PHASE`` is
    MJD minus the fitted peak date, in observer-frame days. Parsed once and
    cached as parquet in data/derived/.
    """
    import pandas as pd

    src = _require(DES_DIR / "DES5YR_SALT3_LCFIT.LCPLOT.gz")
    cache = DATA_DIR / "derived" / "lcplot.parquet"
    if cache.exists() and cache.stat().st_mtime >= src.stat().st_mtime:
        return pd.read_parquet(cache)
    lc = pd.read_csv(src, sep=r"\s+", comment="#", compression="gzip")
    lc = lc.rename(columns={"FLUXCALERR": "FLUXCAL_ERR"})
    lc["BAND"] = lc["BAND"].astype(str).str.replace("DES-", "", regex=False)
    lc.to_parquet(cache)
    return lc


def load_metadata():
    """SALT3 fit results and classification for the Hubble-diagram sample.

    The file is SNANA's whitespace format: a ``VARNAMES:`` line and one
    ``SN:`` row per supernova (DES and low-redshift samples together).
    """
    import pandas as pd

    path = _require(DES_DIR / "DES-SN5YR_HD+MetaData.csv")
    names, rows = None, []
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            names = line.split()[1:]
        elif line.startswith("SN:"):
            rows.append(line.split()[1:])
    if names is None:
        # Plain CSV variant.
        return pd.read_csv(path, comment="#")
    meta = pd.DataFrame(rows, columns=names[: len(rows[0])])
    for col in meta.columns:
        converted = pd.to_numeric(meta[col], errors="coerce")
        if converted.notna().mean() > 0.95:
            meta[col] = converted
    return meta


class _PandasOnlyUnpickler(pickle.Unpickler):
    """Unpickle pandas/numpy containers and nothing else.

    A pickle can run arbitrary code when loaded. White et al.'s widths are
    pickled DataFrames, so only the classes needed to rebuild a DataFrame are
    allowed; anything else raises.
    """

    ALLOWED_PREFIXES = ("numpy", "pandas")
    ALLOWED_BUILTINS = {("builtins", "slice"), ("builtins", "object"),
                        ("builtins", "list"), ("builtins", "dict"), ("builtins", "tuple"),
                        ("copyreg", "_reconstructor"), ("_codecs", "encode")}

    def find_class(self, module, name):
        if (module, name) in self.ALLOWED_BUILTINS or module.split(".")[0] in self.ALLOWED_PREFIXES:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f"refusing to load {module}.{name}")


def load_white_pickle(name: str):
    """One of White et al.'s pickled results, loaded with the restricted unpickler."""
    with open(_require(WHITE_DIR / name), "rb") as f:
        return _PandasOnlyUnpickler(f).load()


def load_white_widths(band: str, *, nonderedshifted: bool = False):
    """White et al.'s per-supernova widths in one band: CID, z, Width, Width_err."""
    suffix = "_non_deredshifted" if nonderedshifted else ""
    return load_white_pickle(f"{band}BAND_pickled_data{suffix}")


def paper_file(key: str, name: str) -> Path:
    """Path of a file inside an unpacked arXiv source (``key`` from ARXIV_SOURCES)."""
    return _require(PAPER_DIR / key / name)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    args = parser.parse_args(argv)
    manifest = fetch_all(refresh=args.refresh)
    total = sum(v["bytes"] for v in manifest.values())
    print(f"\n{len(manifest)} files, {total / 1e6:.0f} MB in {DATA_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
