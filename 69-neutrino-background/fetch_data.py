"""Download and cache everything the analysis needs, and load it back.

What this project reads, all public:

1. **Planck 2018 temperature spectrum** (D3's data): the binned TT band
   powers with their errors, and Planck's best-fit model.
2. **Planck 2018 plik_lite**: Planck's foreground-marginalised high-ℓ
   likelihood (TT, TE, EE band powers and their covariance), as distributed
   for the Cobaya sampler.
3. **CMB lensing**: Planck 2018's lensing band powers (the "consext8" set
   Planck's parameter paper used), and ACT DR6 + Planck PR4 lensing (the set
   DESI DR2 used). ACT's release is a 361 MB archive holding 2.9 GB of
   correction matrices; we stream through it once and keep only their
   projections onto ACT's bins (~2 MB, see :func:`reduce_act_lensing`).
4. **BAO**: DESI DR2 (D5's data), and the 2018 compilation Planck used
   (6dFGS, SDSS MGS, BOSS DR12).
5. **Supernovae**: DES-SN5YR and Pantheon+ (D5's data).
6. **Answer keys** (never inputs to a fit): DESI DR2's neutrino-mass runs
   (best fits and chains), Planck 2018's base_nnu/base_mnu chains (single files
   read out of the 11.8 GB grid archive by HTTP range requests), and ACT DR6's
   P-ACT-LB N_eff chain.
7. **NuFIT 6.0**: the global oscillation fit's Δχ² tables.

Files that an earlier atlas project already downloaded are **hard-linked** from
its cache instead of fetched again, and the link is recorded in the manifest.
Everything lands in ``data/raw/`` with each file's URL, size, SHA-256 and
retrieval time; nothing cached is fetched twice, so after one run the whole
project works offline.

    python fetch_data.py              # everything (~0.65 GB)
    python fetch_data.py --no-chains  # skip the MCMC chains (answer keys, ~0.3 GB)
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import re
import shutil
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import config as cfg

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

PLANCK_DIR = cfg.RAW_DIR / "planck"
BAO_DIR = cfg.RAW_DIR / "bao"
SN_DIR = cfg.RAW_DIR / "sn"
DESI_DIR = cfg.RAW_DIR / "desi_results"

PLA_COSMO = "https://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID={}"
PLANCK_FILES = {
    "COM_PowerSpect_CMB-TT-binned_R3.01.txt": 7143,
    "COM_PowerSpect_CMB-TT-full_R3.01.txt": 170547,
    "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt": 205647,
    "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum_R3.01.txt": 8842,
}

BAO_RAW = "https://raw.githubusercontent.com/CobayaSampler/bao_data/master/"
BAO_FILES = {
    "desi_dr2_mean.txt": "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt",
    "desi_dr2_cov.txt": "desi_bao_dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt",
}

SN_RAW = "https://raw.githubusercontent.com/CobayaSampler/sn_data/master/"
SN_FILES = {
    "PantheonPlus": {"Pantheon+SH0ES.dat": "PantheonPlus/Pantheon+SH0ES.dat",
                     "Pantheon+SH0ES_STAT+SYS.cov": "PantheonPlus/Pantheon+SH0ES_STAT+SYS.cov"},
    "DESY5": {"DES-SN5YR_HD.csv": "DESY5/DES-SN5YR_HD.csv",
              "covsys_000.txt": "DESY5/covsys_000.txt"},
}

DESI_RESULTS_URL = "https://data.desi.lbl.gov/public/papers/y3/bao-cosmo-params/"
DESI_CHECKSUMS = "dr2_vac_dr2_bao-cosmo-params_v1.0.sha256sum"
_LOWL = "planck2018-lowl-TT-clik_planck2018-lowl-EE-clik"
CMB_VARIANTS = {
    "CamSpec": f"{_LOWL}_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing",
    "plik": f"{_LOWL}_planck2018-highl-plik-TTTEEE_planck-act-dr6-lensing",
    "HiLLiPoP": ("planck2018-lowl-TT-clik_planck2020-lollipop-lowlE_planck2020-hillipop-TTTEEE_"
                 "planck-act-dr6-lensing"),
}
#: DESI runs whose full MCMC chains we fetch (about 80 MB each).
DESI_CHAIN_RUNS = [
    f"base_mnu/desi-bao-all_{CMB_VARIANTS['CamSpec']}",
    f"base_mnu/desi-bao-all_{CMB_VARIANTS['plik']}",
    f"base_mnu_w_wa/desi-bao-all_desy5sn_{CMB_VARIANTS['CamSpec']}",
]
#: Model folders whose small summary files (margestats, yaml, best fits) we fetch.
DESI_SUMMARY_MODELS = ("base_mnu", "base_mnu059", "base_mnu100", "base_mnu_w", "base_mnu_w_wa")

GITHUB_RELEASE = "https://github.com/CobayaSampler/planck_native_data/releases/download/v1/"
PLIK_LITE_ZIP = "plik_lite_2018_AL.zip"
PLIK_LITE_DIR = cfg.RAW_DIR / "plik_lite"
PLIK_LITE_FILES = ("cl_cmb_plik_v22.dat", "c_matrix_plik_v22.dat", "blmin.dat", "blmax.dat",
                   "bweight.dat", "plik_lite_v22.dataset")

#: Planck 2018 low-ℓ likelihoods in native Python form (Cobaya's translation of
#: Planck's clik files): Commander TT, ℓ = 2-29 (a Gaussianised Blackwell-Rao
#: likelihood) and SimAll EE, ℓ = 2-29 (a table of log-likelihoods).
LOWL_ZIPS = {"planck_2018_lowT.zip": (5_523_441, "lowl_tt"),
             "planck_2018_lowE.zip": (841_455, "lowl_ee")}
LOWL_DIR = cfg.RAW_DIR / "planck_lowl"

PLANCK_LENS_RAW = ("https://raw.githubusercontent.com/CobayaSampler/planck_supp_data_and_covmats/"
                   "v2.1/lensing/2018/")
PLANCK_LENS_ROOT = "smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8"
PLANCK_LENS_DIR = cfg.RAW_DIR / "planck_lensing_2018"
PLANCK_LENS_FILES = ([f"{PLANCK_LENS_ROOT}{x}" for x in
                      (".dataset", "_bandpowers.dat", "_cov.dat",
                       "_lensing_fiducial_correction.dat")]
                     + [f"{PLANCK_LENS_ROOT}_window/window{i}.dat" for i in range(1, 10)]
                     + [f"{PLANCK_LENS_ROOT}_lens_delta_window/window{i}.dat" for i in range(1, 10)])

ACT_LENS_URL = ("https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/"
                "ACT_dr6_likelihood_v1.2.tgz")
ACT_LENS_BYTES = 361_306_879
ACT_LENS_DIR = cfg.RAW_DIR / "act_dr6_lensing"
ACT_LENS_REDUCED = cfg.DERIVED_DIR / "act_lensing_reduced.npz"

BAO2018_FILES = {"sdss_DR12Consensus_bao.dat": "sdss_DR12Consensus_bao.dat",
                 "BAO_consensus_covtot_dM_Hz.txt": "BAO_consensus_covtot_dM_Hz.txt",
                 "sdss_MGS_prob.txt": "sdss_MGS_prob.txt"}

NUFIT_URL = "http://www.nu-fit.org/sites/default/files/"
NUFIT_DIR = cfg.RAW_DIR / "nufit"
#: "TByes" = IC24 with Super-Kamiokande atmospheric data; "TBoff" = IC19 without.
NUFIT_FILES = {"v60.release-TByes-NO.txt.xz": 4_222_132, "v60.release-TByes-IO.txt.xz": 4_293_576,
               "v60.release-TBoff-NO.txt.xz": 4_378_772, "v60.release-TBoff-IO.txt.xz": 4_427_612}

PLANCK_GRID_URL = ("https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/"
                   "cosmoparams/COM_CosmoParams_fullGrid_R3.01.zip")
PLANCK_GRID_DIR = cfg.RAW_DIR / "planck_grid"
#: Planck runs whose summary statistics we read (folder/run name in the grid).
PLANCK_GRID_RUNS = {
    "nnu_TT": "base_nnu/plikHM_TT_lowl_lowE/dist/base_nnu_plikHM_TT_lowl_lowE",
    "nnu_TTTEEE": "base_nnu/plikHM_TTTEEE_lowl_lowE/dist/base_nnu_plikHM_TTTEEE_lowl_lowE",
    "nnu_TTTEEE_lens_BAO": ("base_nnu/plikHM_TTTEEE_lowl_lowE_BAO/dist/"
                            "base_nnu_plikHM_TTTEEE_lowl_lowE_BAO_post_lensing"),
    "mnu_TT": "base_mnu/plikHM_TT_lowl_lowE/dist/base_mnu_plikHM_TT_lowl_lowE",
    "mnu_TTTEEE": "base_mnu/plikHM_TTTEEE_lowl_lowE/dist/base_mnu_plikHM_TTTEEE_lowl_lowE",
    "mnu_TTTEEE_lens": ("base_mnu/plikHM_TTTEEE_lowl_lowE_lensing/dist/"
                        "base_mnu_plikHM_TTTEEE_lowl_lowE_lensing"),
    "mnu_TTTEEE_lens_BAO": ("base_mnu/plikHM_TTTEEE_lowl_lowE_lensing_BAO/dist/"
                            "base_mnu_plikHM_TTTEEE_lowl_lowE_lensing_BAO"),
}
#: ... and whose chains we read (the "lensing + BAO" validation runs).
PLANCK_GRID_CHAINS = {
    "nnu_TTTEEE_lens_BAO": ("base_nnu/plikHM_TTTEEE_lowl_lowE_BAO/"
                            "base_nnu_plikHM_TTTEEE_lowl_lowE_BAO_post_lensing"),
    "mnu_TTTEEE_lens_BAO": ("base_mnu/plikHM_TTTEEE_lowl_lowE_lensing_BAO/"
                            "base_mnu_plikHM_TTTEEE_lowl_lowE_lensing_BAO"),
}

#: Escudero, Jackson, Laine & Sandner (2026; arXiv:2511.04747): electron-mass dependent
#: coefficients of the neutrino energy-transfer rates, from their public code (GPL-3.0).
NUDEC_URL = ("https://raw.githubusercontent.com/MiguelEA/nudec_BSM/v2/BasicModules_data/"
             "rate_coefficients_neutrinos.dat")
NUDEC_SHA256 = "dbe3d64936990ddd5cfd5730f19aa298159581833a1175cc8f6badf73d155851"

#: ... and whose best fits (with each likelihood's chi^2) test our likelihoods.
PLANCK_GRID_MINIMA = {
    "mnu_TTTEEE_lens_BAO": ("base_mnu/plikHM_TTTEEE_lowl_lowE_lensing_BAO/"
                            "base_mnu_plikHM_TTTEEE_lowl_lowE_lensing_BAO"),
    "mnu_TTTEEE_lens": ("base_mnu/plikHM_TTTEEE_lowl_lowE_lensing/"
                        "base_mnu_plikHM_TTTEEE_lowl_lowE_lensing"),
}

ACT_CHAIN_URL = "https://lambda.gsfc.nasa.gov/data/act/chains/neff/p-actlite-l-b_neff_camb.tar.gz"
ACT_CHAIN_DIR = cfg.RAW_DIR / "act_chains"

#: Sibling projects' copies: {our path: their path}. Same bytes, no download.
SIBLING_COPIES = {}
for _name in PLANCK_FILES:
    SIBLING_COPIES[PLANCK_DIR / _name] = cfg.D3_DIR / "data" / "planck" / _name
for _local in BAO_FILES:
    SIBLING_COPIES[BAO_DIR / _local] = cfg.D5_DIR / "data" / "bao" / _local
for _sample, _files in SN_FILES.items():
    for _local in _files:
        SIBLING_COPIES[SN_DIR / _sample / _local] = cfg.D5_DIR / "data" / "sn" / _sample / _local
SIBLING_COPIES[DESI_DIR / DESI_CHECKSUMS] = cfg.D5_DIR / "data" / "desi_results" / DESI_CHECKSUMS
SIBLING_COPIES[BAO_DIR / "sdss_MGS_prob.txt"] = cfg.D5_DIR / "data" / "bao" / "sdss_mgs_prob.txt"


# ---------------------------------------------------------------------------
# Manifest and downloads
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _load_manifest() -> dict:
    return json.loads(cfg.MANIFEST_PATH.read_text()) if cfg.MANIFEST_PATH.exists() else {}


def _save_manifest(manifest: dict) -> None:
    tmp = cfg.MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    tmp.replace(cfg.MANIFEST_PATH)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _get(url: str, retries: int = 5, timeout: int = 600, headers: dict | None = None) -> bytes:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT,
                                                       **(headers or {})})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if attempt == retries - 1:
                raise
            print(f"    retrying after {err!r}")
            time.sleep(2 + 3 * attempt)
    raise RuntimeError("unreachable")


def _sibling_sha256(sibling: Path) -> str | None:
    """The checksum a sibling project recorded for its copy, if it kept a manifest."""
    for parent in sibling.parents:
        man = parent / "manifest.json"
        if man.exists():
            try:
                entries = json.loads(man.read_text())
            except json.JSONDecodeError:
                return None
            for key in (str(sibling.relative_to(parent)),):
                entry = entries.get(key)
                if isinstance(entry, dict):
                    return entry.get("sha256")
            return None
        if parent.name == "data":
            break
    return None


def _link_from_sibling(dest: Path, manifest: dict, url: str) -> bool:
    src = SIBLING_COPIES.get(dest)
    if src is None or not src.exists():
        return False
    sha = _sha256(src)
    recorded = _sibling_sha256(src)
    if recorded is not None and recorded != sha:
        print(f"    {src} does not match its own manifest; downloading instead")
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(src, dest)
    except OSError:
        shutil.copy2(src, dest)
    key = str(dest.relative_to(cfg.DATA_DIR))
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha,
                     "linked_from": str(src.relative_to(cfg.ROOT_DIR)), "retrieved_utc": _now()}
    print(f"  linked   {key}  (from {src.relative_to(cfg.ROOT_DIR).parts[0]}, no download)")
    return True


def download(url: str, dest: Path, manifest: dict, *, refresh: bool = False,
             expected_sha256: str | None = None, expected_bytes: int | None = None) -> Path:
    """Fetch one file unless cached (or available from a sibling project)."""
    key = str(dest.relative_to(cfg.DATA_DIR))
    if dest.exists() and not refresh:
        if key not in manifest:
            manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                             "retrieved_utc": datetime.fromtimestamp(
                                 dest.stat().st_mtime, timezone.utc).isoformat(timespec="seconds")}
        return dest
    if not refresh and _link_from_sibling(dest, manifest, url):
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    blob = _get(url)
    if expected_bytes is not None and len(blob) != expected_bytes:
        raise IOError(f"{dest.name}: got {len(blob)} bytes, expected {expected_bytes}")
    tmp.write_bytes(blob)
    sha = _sha256(tmp)
    if expected_sha256 and sha != expected_sha256:
        tmp.unlink()
        raise IOError(f"{dest.name}: checksum mismatch (got {sha[:12]}, expected "
                      f"{expected_sha256[:12]})")
    tmp.replace(dest)
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha,
                     "retrieved_utc": _now()}
    print(f"  fetched  {key}  ({dest.stat().st_size / 1e6:.2f} MB)")
    return dest


def download_stream(url: str, dest: Path, manifest: dict, *, expected_bytes: int | None = None,
                    refresh: bool = False, chunk: int = 1 << 22) -> Path:
    """Fetch a large file to disk in chunks, resuming a partial download by Range."""
    key = str(dest.relative_to(cfg.DATA_DIR))
    if dest.exists() and not refresh:
        if key not in manifest:
            manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                             "retrieved_utc": _now()}
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(8):
        have = tmp.stat().st_size if tmp.exists() else 0
        if expected_bytes and have >= expected_bytes:
            break
        headers = {"User-Agent": cfg.USER_AGENT}
        if have:
            headers["Range"] = f"bytes={have}-"
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "ab" if have else "wb") as f:
                if have and r.status != 206:  # server ignored the range: start over
                    f.truncate(0)
                t0, got = time.monotonic(), 0
                while True:
                    block = r.read(chunk)
                    if not block:
                        break
                    f.write(block)
                    got += len(block)
                    if got % (64 * chunk) < chunk:
                        print(f"    {dest.name}: {(have + got) / 1e6:.0f} MB "
                              f"({got / 1e6 / max(time.monotonic() - t0, 1e-3):.1f} MB/s)")
            if expected_bytes is None or tmp.stat().st_size == expected_bytes:
                break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            print(f"    retrying after {err!r}")
            time.sleep(3 + 5 * attempt)
    if expected_bytes is not None and tmp.stat().st_size != expected_bytes:
        raise IOError(f"{dest.name}: {tmp.stat().st_size} bytes, expected {expected_bytes}")
    tmp.replace(dest)
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                     "retrieved_utc": _now()}
    print(f"  fetched  {key}  ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def fetch_plik_lite(manifest: dict, refresh: bool = False) -> None:
    import zipfile

    if all((PLIK_LITE_DIR / f).exists() for f in PLIK_LITE_FILES) and not refresh:
        return
    zpath = download(GITHUB_RELEASE + PLIK_LITE_ZIP, PLIK_LITE_DIR / PLIK_LITE_ZIP, manifest,
                     refresh=refresh, expected_bytes=1_546_224)
    with zipfile.ZipFile(zpath) as z:
        for info in z.infolist():
            name = Path(info.filename).name
            if name in PLIK_LITE_FILES:
                data = z.read(info)
                (PLIK_LITE_DIR / name).write_bytes(data)
                manifest[str((PLIK_LITE_DIR / name).relative_to(cfg.DATA_DIR))] = {
                    "member_of": PLIK_LITE_ZIP, "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(), "retrieved_utc": _now()}
    missing = [f for f in PLIK_LITE_FILES if not (PLIK_LITE_DIR / f).exists()]
    if missing:
        raise IOError(f"plik_lite archive lacks {missing}")


def fetch_lowl(manifest: dict, refresh: bool = False) -> None:
    import zipfile

    for name, (size, sub) in LOWL_ZIPS.items():
        zpath = download(GITHUB_RELEASE + name, LOWL_DIR / name, manifest, refresh=refresh,
                         expected_bytes=size)
        dest = LOWL_DIR / sub
        if dest.exists() and any(dest.iterdir()) and not refresh:
            continue
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zpath) as z:
            for info in z.infolist():
                if info.is_dir() or not info.filename.endswith(".txt"):
                    continue
                data = z.read(info)
                out = dest / Path(info.filename).name
                out.write_bytes(data)
                manifest[str(out.relative_to(cfg.DATA_DIR))] = {
                    "member_of": name, "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(), "retrieved_utc": _now()}


class _RangeFile(io.RawIOBase):
    """A remote file read by HTTP Range requests (for single members of big ZIPs)."""

    def __init__(self, url: str):
        self.url = url
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": cfg.USER_AGENT})
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])
        self.pos = 0

    def seekable(self):
        return True

    def readable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def readinto(self, b):
        n = len(b)
        if n == 0 or self.pos >= self.size:
            return 0
        end = min(self.pos + n, self.size) - 1
        data = None
        for attempt in range(5):
            try:
                req = urllib.request.Request(self.url, headers={"User-Agent": cfg.USER_AGENT,
                                                                "Range": f"bytes={self.pos}-{end}"})
                with urllib.request.urlopen(req, timeout=300) as r:
                    if r.status != 206:
                        raise IOError("server ignored the Range header")
                    data = r.read()
                break
            except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
                if attempt == 4:
                    raise
                print(f"    retrying range read after {err!r}")
                time.sleep(2 + 3 * attempt)
        b[: len(data)] = data
        self.pos += len(data)
        return len(data)


def fetch_planck_grid(manifest: dict, *, chains: bool, refresh: bool = False) -> None:
    """Single files out of Planck's 11.8 GB parameter-grid ZIP (IRSA serves byte ranges)."""
    import zipfile

    wanted = []
    for run in PLANCK_GRID_RUNS.values():
        wanted += [run + ".margestats"]
    for run in PLANCK_GRID_MINIMA.values():
        wanted += [run + ".minimum"]
    if chains:
        for run in PLANCK_GRID_CHAINS.values():
            wanted += [f"{run}_{i}.txt" for i in (1, 2, 3, 4)] + [run + ".paramnames"]
    todo = [w for w in wanted if refresh or not (PLANCK_GRID_DIR / w).exists()]
    if not todo:
        return
    print(f"  reading {len(todo)} files out of the Planck grid archive by HTTP range")
    zf = zipfile.ZipFile(io.BufferedReader(_RangeFile(PLANCK_GRID_URL), buffer_size=1 << 20))
    names = set(zf.namelist())
    for w in todo:
        member = next((n for n in names if n.endswith(w)), None)
        if member is None:
            raise KeyError(f"{w} not in the Planck grid archive")
        data = zf.read(member)
        dest = PLANCK_GRID_DIR / w
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
        manifest[str(dest.relative_to(cfg.DATA_DIR))] = {
            "url": PLANCK_GRID_URL, "member": member, "bytes": len(data),
            "sha256": hashlib.sha256(data).hexdigest(), "retrieved_utc": _now()}
        print(f"  fetched  {dest.relative_to(cfg.DATA_DIR)}  ({len(data) / 1e6:.2f} MB)")


def fetch_act_chain(manifest: dict, refresh: bool = False) -> None:
    import tarfile

    tgz = download(ACT_CHAIN_URL, ACT_CHAIN_DIR / Path(ACT_CHAIN_URL).name, manifest,
                   refresh=refresh, expected_bytes=9_854_053)
    if any(ACT_CHAIN_DIR.glob("*.1.txt")) and not refresh:
        return
    with tarfile.open(tgz, "r:gz") as tar:
        for m in tar.getmembers():
            name = Path(m.name).name
            if m.isfile() and re.search(r"\.(\d+\.txt|updated\.yaml|minimum\.txt|minimum)$", name):
                (ACT_CHAIN_DIR / name).write_bytes(tar.extractfile(m).read())


# ---------------------------------------------------------------------------
# ACT DR6 lensing: reduce 2.9 GB of correction matrices to their binned form
# ---------------------------------------------------------------------------

ACT_SMALL = ("clkk_bandpowers_act.txt", "clkk_bandpowers_planck.txt", "binning_matrix_act.txt",
             "binning_matrix_planck.txt", "covmat_actplanck.txt", "covmat_act.txt",
             "cosmo2017_10K_acc3_lensedCls.dat", "cosmo2017_10K_acc3_lenspotentialCls.dat",
             "n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt",
             "PLANCK_n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax3000.txt")
ACT_BIG = {"norm_correction_matrix_Lmin0_Lmax4000.npy": ("norm", ""),
           "P18_norm_correction_matrix_Lmin0_Lmax3000.npy": ("norm", "_planck")}
for _s in ("KK", "TT", "EE", "BB", "TE"):
    ACT_BIG[f"N1der_{_s}_lmin600_lmax3000_full.txt"] = ("n1_" + _s.lower(), "")
    ACT_BIG[f"N1_planck_der_{_s}_lmin100_lmax2048.txt"] = ("n1_" + _s.lower(), "_planck")
#: Theory vectors are standardised to L = 0 .. TRIM_LMAX + 1 (ACT's convention).
ACT_TRIM_LMAX = 2998


def _standardize(ls, cls, nlen=ACT_TRIM_LMAX + 2, dims="y"):
    """ACT's ``standardize``: place arrays that start at ls[0] onto L = 0..nlen-1."""
    c0 = int(ls[0])
    cend = nlen - c0
    if dims == "y":
        out = np.zeros(nlen)
        out[c0:] = cls[:cend]
    elif dims == "yy":
        out = np.zeros((nlen, nlen))
        out[c0:, c0:] = cls[:cend, :cend]
    elif dims == "xy":
        out = np.zeros((cls.shape[0], nlen))
        out[:, c0:] = cls[:, :cend]
    elif dims == "xyy":
        out = np.zeros((cls.shape[0], nlen, nlen))
        out[:, c0:, c0:] = cls[:, :cend, :cend]
    else:
        raise ValueError(dims)
    return out


def _act_members(tgz: Path, wanted: set):
    """Yield (name, file object) for wanted members, streaming through the archive once."""
    import tarfile

    with tarfile.open(tgz, mode="r|gz") as tar:
        for m in tar:
            name = Path(m.name).name
            if m.isfile() and name in wanted:
                yield name, tar.extractfile(m)


def _npy_blocks(fobj):
    """Yield arr[0], arr[1], ... of a C-ordered .npy array read from a stream.

    numpy's own reader wants a seekable real file; a member of a streamed tar
    archive is neither, and the full array here is 0.5 GB.
    """
    version = np.lib.format.read_magic(fobj)
    if version == (1, 0):
        shape, fortran, dtype = np.lib.format.read_array_header_1_0(fobj)
    else:
        shape, fortran, dtype = np.lib.format.read_array_header_2_0(fobj)
    if fortran:
        raise ValueError("Fortran-ordered array not supported")
    block_shape = shape[1:]
    nbytes = int(np.prod(block_shape)) * dtype.itemsize
    for _ in range(shape[0]):
        buf = bytearray()
        while len(buf) < nbytes:
            chunk = fobj.read(min(1 << 24, nbytes - len(buf)))
            if not chunk:
                raise IOError("truncated .npy stream")
            buf.extend(chunk)
        yield np.frombuffer(bytes(buf), dtype=dtype).reshape(block_shape)


def reduce_act_lensing(force: bool = False) -> Path:
    """Project ACT DR6 lensing's linear corrections onto the bins we use.

    ACT's likelihood (act_dr6_lenslike v1.2, variant ``actplanck_baseline``)
    corrects the theory C_L^κκ for the dependence of the lensing
    reconstruction on the CMB spectra and on C^κκ itself (normalisation and
    N1 bias), each correction linear in the theory spectra:

        C_corr = C_κκ + Σ_X [−2 (dA_L/dC_X · ΔC_X) / A_L^fid] C_κκ^fid
                 + dN1/dC_κκ · (C_κκ − C_κκ^fid) + Σ_X dN1/dC_X · ΔC_X,

    with X ∈ {TT, EE, BB, TE} and ΔC = C − C^fid. Binned, this is
    K · C_κκ + Σ_X M_X · C_X + const, with K and M_X small (bins × 3000)
    matrices, which is all we keep. Done separately for the ACT (10 bins,
    L = 40–763) and Planck PR4 (9 bins, L = 8–400) parts.
    """
    if ACT_LENS_REDUCED.exists() and not force:
        return ACT_LENS_REDUCED
    tgz = _require(ACT_LENS_DIR / Path(ACT_LENS_URL).name)
    small_dir = ACT_LENS_DIR / "small"
    small_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    if not all((small_dir / f).exists() for f in ACT_SMALL):
        print("  ACT lensing: pass 1, small files")
        for name, fobj in _act_members(tgz, set(ACT_SMALL)):
            (small_dir / name).write_bytes(fobj.read())
    missing = [f for f in ACT_SMALL if not (small_dir / f).exists()]
    if missing:
        raise IOError(f"ACT archive lacks {missing}")
    L = np.arange(ACT_TRIM_LMAX + 2)
    f_ls, f_tt, f_ee, f_bb, f_te = np.loadtxt(small_dir / "cosmo2017_10K_acc3_lensedCls.dat",
                                              unpack=True)
    fac = 2 * np.pi / (f_ls * (f_ls + 1.0))
    fid = {x: _standardize(f_ls, c * fac) for x, c in (("tt", f_tt), ("ee", f_ee), ("bb", f_bb),
                                                       ("te", f_te))}
    fd_ls, f_dd = np.loadtxt(small_dir / "cosmo2017_10K_acc3_lenspotentialCls.dat", unpack=True,
                             usecols=[0, 5])
    fid["kk"] = _standardize(fd_ls, f_dd * 2 * np.pi / 4.0)
    binmat = np.loadtxt(small_dir / "binning_matrix_act.txt")
    B = {"": _standardize(np.arange(binmat.shape[1]), binmat[2:-6], dims="xy"),
         "_planck": _standardize(np.arange(3000),
                                 np.loadtxt(small_dir / "binning_matrix_planck.txt"), dims="xy")}
    fAL = {}
    for suff, fname in (("", "n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt"),
                        ("_planck", "PLANCK_n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax3000.txt")):
        ls, a = np.loadtxt(small_dir / fname)
        fAL[suff] = _standardize(ls, a)
    # factor -2 f_kk / fAL (L >= 2; f_kk vanishes below L = 2 anyway)
    nfac = {}
    for suff in ("", "_planck"):
        g = np.zeros_like(L, dtype=float)
        ok = L >= 2
        g[ok] = -2.0 * fid["kk"][ok] / fAL[suff][ok]
        nfac[suff] = g
    tag = {"": "act", "_planck": "planck"}
    out = {f"{tag[suff]}_K": B[suff].copy() for suff in ("", "_planck")}  # identity part of K
    for suff in ("", "_planck"):
        for x in ("tt", "ee", "bb", "te"):
            out[f"{tag[suff]}_M_{x}"] = np.zeros_like(B[suff])
    print("  ACT lensing: pass 2, correction matrices (streaming ~2.9 GB)")
    seen = set()
    for name, fobj in _act_members(tgz, set(ACT_BIG)):
        kind, suff = ACT_BIG[name]
        b, t = B[suff], tag[suff]
        if kind == "norm":  # (4, n, n): TT, EE, BB, TE, read one n x n block at a time
            for i, block in enumerate(_npy_blocks(fobj)):
                x = ("tt", "ee", "bb", "te")[i]
                blk = _standardize(np.arange(block.shape[0]), block, dims="yy")
                out[f"{t}_M_{x}"] += (b * nfac[suff][None, :]) @ blk
                del block, blk
        else:
            import pandas as pd

            # pandas wants a seekable stream; a streamed tar member is not one
            mat = pd.read_csv(io.BytesIO(fobj.read()), sep=r"\s+", header=None,
                              engine="c").to_numpy(float)
            mat = _standardize(np.arange(mat.shape[0]), mat, dims="yy")
            spec = kind[3:]
            if spec == "kk":
                out[f"{t}_K"] += b @ mat
            else:
                out[f"{t}_M_{spec}"] += b @ mat
            del mat
        seen.add(name)
        print(f"    {name} ({time.monotonic() - t0:.0f} s)")
    if seen != set(ACT_BIG):
        raise IOError(f"ACT archive lacks {sorted(set(ACT_BIG) - seen)}")
    # constant: -B dN1_kk f_kk - sum_X M_X f_X  ==  (B - K) f_kk - sum_X M_X f_X
    for suff in ("", "_planck"):
        t = tag[suff]
        const = (B[suff] - out[f"{t}_K"]) @ fid["kk"]
        for x in ("tt", "ee", "bb", "te"):
            const -= out[f"{t}_M_{x}"] @ fid[x]
        out[f"{t}_const"] = const
        out[f"{t}_B"] = B[suff]
    data_act = np.loadtxt(small_dir / "clkk_bandpowers_act.txt")
    out["data"] = np.concatenate([data_act[2:-6],
                                  np.loadtxt(small_dir / "clkk_bandpowers_planck.txt")])
    fcov = np.loadtxt(small_dir / "covmat_actplanck.txt")
    n_act = data_act.size
    idx = [i for i in range(fcov.shape[0]) if not (i < 2 or n_act - 6 <= i < n_act)]
    out["cov"] = fcov[np.ix_(idx, idx)]
    for x, arr in fid.items():
        out[f"fid_{x}"] = arr
    np.savez(ACT_LENS_REDUCED, **out)
    print(f"  ACT lensing reduced in {time.monotonic() - t0:.0f} s -> "
          f"{ACT_LENS_REDUCED.relative_to(cfg.DATA_DIR)}")
    return ACT_LENS_REDUCED


# ---------------------------------------------------------------------------
# DESI's neutrino runs
# ---------------------------------------------------------------------------


def desi_checksums(manifest: dict, refresh: bool = False) -> dict[str, str]:
    path = download(DESI_RESULTS_URL + DESI_CHECKSUMS, DESI_DIR / DESI_CHECKSUMS, manifest,
                    refresh=refresh)
    sums = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            sums[parts[1].removeprefix("./")] = parts[0]
    return sums


def desi_wanted_files(sums: dict[str, str], *, chains: bool) -> list[str]:
    wanted = []
    for rel in sums:
        parts = rel.split("/")
        if len(parts) < 4 or parts[1] not in DESI_SUMMARY_MODELS:
            continue
        run = "/".join(parts[1:3])
        name = parts[-1]
        if parts[0] == "iminuit" and name == "bestfit.minimum.txt":
            wanted.append(rel)
        if parts[0] == "cobaya":
            if name in ("chain.updated.yaml", "chain.margestats"):
                wanted.append(rel)
            elif chains and run in DESI_CHAIN_RUNS and re.fullmatch(r"chain\.\d+\.txt", name):
                wanted.append(rel)
    return sorted(wanted)


# ---------------------------------------------------------------------------
# Fetch everything
# ---------------------------------------------------------------------------


def fetch_all(*, refresh: bool = False, chains: bool = True, act_lensing: bool = True) -> dict:
    manifest = _load_manifest()
    try:
        print("Planck 2018 temperature spectrum and best fit (D3's data)")
        for name, size in PLANCK_FILES.items():
            download(PLA_COSMO.format(name), PLANCK_DIR / name, manifest, refresh=refresh,
                     expected_bytes=size)
        print("Planck 2018 plik_lite (TT, TE, EE)")
        fetch_plik_lite(manifest, refresh)
        print("Planck 2018 low-ℓ likelihoods (TT and EE, ℓ = 2-29)")
        fetch_lowl(manifest, refresh)
        print("Planck 2018 lensing")
        for name in PLANCK_LENS_FILES:
            download(PLANCK_LENS_RAW + name, PLANCK_LENS_DIR / name, manifest, refresh=refresh)
        print("BAO: DESI DR2 (D5's data) and the 2018 compilation Planck used")
        for local, remote in BAO_FILES.items():
            download(BAO_RAW + remote, BAO_DIR / local, manifest, refresh=refresh)
        for local, remote in BAO2018_FILES.items():
            download(BAO_RAW + remote, BAO_DIR / local, manifest, refresh=refresh)
        print("Supernovae: Pantheon+ and DES-SN5YR (D5's data)")
        for sample, files in SN_FILES.items():
            for local, remote in files.items():
                download(SN_RAW + remote, SN_DIR / sample / local, manifest, refresh=refresh)
        print("Neutrino decoupling: rate coefficients (Escudero et al. 2026)")
        download(NUDEC_URL, cfg.RAW_DIR / "nudec" / "rate_coefficients_neutrinos.dat", manifest,
                 refresh=refresh, expected_sha256=NUDEC_SHA256)
        print("NuFIT 6.0 oscillation fits")
        for name, size in NUFIT_FILES.items():
            download(NUFIT_URL + name, NUFIT_DIR / name, manifest, refresh=refresh,
                     expected_bytes=size)
        print("DESI DR2's neutrino-mass runs (answer key)" + ("" if chains else ", no chains"))
        sums = desi_checksums(manifest, refresh=refresh)
        for rel in desi_wanted_files(sums, chains=chains):
            download(DESI_RESULTS_URL + rel, DESI_DIR / rel, manifest, refresh=refresh,
                     expected_sha256=sums[rel])
        print("Planck 2018 parameter grid: summaries" + (" and chains" if chains else ""))
        fetch_planck_grid(manifest, chains=chains, refresh=refresh)
        if chains:
            print("ACT DR6 P-ACT-LB N_eff chain (answer key)")
            fetch_act_chain(manifest, refresh)
        if act_lensing:
            print("ACT DR6 + Planck PR4 lensing likelihood data (361 MB archive)")
            download_stream(ACT_LENS_URL, ACT_LENS_DIR / Path(ACT_LENS_URL).name, manifest,
                            expected_bytes=ACT_LENS_BYTES, refresh=refresh)
    finally:
        _save_manifest(manifest)
    if act_lensing:
        reduce_act_lensing()
    return manifest


# ---------------------------------------------------------------------------
# Loaders -- read only from the cache, never from the network
# ---------------------------------------------------------------------------


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


def load_planck_binned_tt() -> np.ndarray:
    """Columns: l_eff, D_l, -dD_l, +dD_l, best fit (µK²); 83 bands, l = 30-2508."""
    return np.loadtxt(_require(PLANCK_DIR / "COM_PowerSpect_CMB-TT-binned_R3.01.txt"))


def load_planck_full_tt() -> np.ndarray:
    """Columns: l, D_l, -dD_l, +dD_l (µK²); every multipole l = 2-2508."""
    return np.loadtxt(_require(PLANCK_DIR / "COM_PowerSpect_CMB-TT-full_R3.01.txt"))


def load_planck_bestfit_theory() -> dict:
    """Planck 2018 baseline best-fit model: l, TT, TE, EE, BB, PP (D_l in µK²)."""
    a = np.loadtxt(_require(PLANCK_DIR / "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-"
                                         "lensing-minimum-theory_R3.01.txt"))
    return {"ell": a[:, 0], "tt": a[:, 1], "te": a[:, 2], "ee": a[:, 3], "bb": a[:, 4],
            "pp": a[:, 5]}


def load_planck_bestfit_params() -> dict:
    """Parameters of that best fit (name -> value), from the PLA .minimum file."""
    out = {}
    for line in _require(PLANCK_DIR / "COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-"
                                      "lensing-minimum_R3.01.txt").read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                out[parts[2]] = float(parts[1])
            except ValueError:
                pass
    return out


def load_plik_lite() -> dict:
    """Planck 2018 plik_lite v22: 613 binned C_ℓ (µK²; TT 215, TE 199, EE 199 bins).

    Bin edges and weights are stored three times (once per spectrum, identical);
    we keep the first copy. The covariance is a Fortran unformatted record of
    613 × 613 float64 with only the lower triangle filled.
    """
    d = PLIK_LITE_DIR
    raw = _require(d / "c_matrix_plik_v22.dat").read_bytes()
    n = 613
    if len(raw) != 8 * n * n + 8:
        raise ValueError("unexpected plik_lite covariance size")
    cov = np.frombuffer(raw[4:-4], "<f8").reshape(n, n)
    cov = np.tril(cov) + np.tril(cov, -1).T
    data = np.loadtxt(_require(d / "cl_cmb_plik_v22.dat"))
    blmin = np.loadtxt(d / "blmin.dat").astype(int)[:215]
    blmax = np.loadtxt(d / "blmax.dat").astype(int)[:215]
    bw = np.loadtxt(d / "bweight.dat")[: 2509 - 30]
    if not np.allclose(np.sqrt(np.diag(cov)), data[:, 2], rtol=1e-6):
        raise ValueError("plik_lite covariance diagonal does not match the data errors")
    return {"nbins": {"TT": 215, "TE": 199, "EE": 199}, "plmin": 30, "blmin": blmin,
            "blmax": blmax, "bweight": bw, "cl": data[:, 1], "err": data[:, 2],
            "ell_centre": data[:, 0], "cov": cov}


def load_lowl_tt() -> dict:
    """Planck 2018 Commander TT, ℓ = 2-29 (Gaussianised Blackwell-Rao): the mean
    and covariance of the Gaussianised variables x_ℓ, the x(C_ℓ) transform as
    tables (cl2x_1: D_ℓ grid, cl2x_2: x values; 1000 rows x one column per ℓ), and
    the reference spectrum used to normalise the likelihood."""
    d = LOWL_DIR / "lowl_tt"
    return {"cov": np.loadtxt(_require(d / "cov.txt")), "mu": np.loadtxt(d / "mu.txt"),
            "mu_sigma": np.loadtxt(d / "mu_sigma.txt"), "cl_grid": np.loadtxt(d / "cl2x_1.txt"),
            "x_grid": np.loadtxt(d / "cl2x_2.txt")}


def load_lowl_ee() -> np.ndarray:
    """Planck 2018 SimAll EE, ℓ = 2-29: log-likelihood table, rows D_ℓ^EE = i × 1e-4 µK²."""
    return np.loadtxt(_require(LOWL_DIR / "lowl_ee" / "prob_table.txt"))


def load_planck_lensing() -> dict:
    """Planck 2018 lensing (consext8): 9 band powers in [L(L+1)]² C^φφ/2π, L = 8-400.

    Each band is Σ_L W_b(L) PP(L) plus a linear correction for the CMB spectra
    (the "lens delta windows", acting on D_ℓ^TT, EE, TE in µK² and on PP for
    ℓ = 2-2500), minus the correction evaluated on Planck's fiducial model.
    """
    d = PLANCK_LENS_DIR
    root = PLANCK_LENS_ROOT
    bp = np.loadtxt(_require(d / f"{root}_bandpowers.dat"))
    cov = np.loadtxt(d / f"{root}_cov.dat")
    fid = np.loadtxt(d / f"{root}_lensing_fiducial_correction.dat")[:, 1]
    lmax = 2500
    W = np.zeros((9, lmax + 1))
    M = {x: np.zeros((9, lmax + 1)) for x in ("tt", "ee", "te", "pp")}
    for b in range(9):
        w = np.loadtxt(d / f"{root}_window/window{b + 1}.dat")
        W[b, w[:, 0].astype(int)] = w[:, 1]
        m = np.loadtxt(d / f"{root}_lens_delta_window/window{b + 1}.dat")
        m = m[(m[:, 0] >= 2) & (m[:, 0] <= lmax)]
        ls = m[:, 0].astype(int)
        for k, x in enumerate(("tt", "ee", "te", "pp")):
            M[x][b, ls] = m[:, 1 + k]
    return {"data": bp[:, 4], "cov": cov, "W": W, "M": M, "fid_correction": fid,
            "lmin": bp[:, 1], "lmax": bp[:, 2], "lmax_theory": lmax}


def load_act_lensing() -> dict:
    """ACT DR6 + Planck PR4 lensing, reduced (see :func:`reduce_act_lensing`)."""
    with np.load(_require(ACT_LENS_REDUCED)) as z:
        return {k: z[k] for k in z.files}


def load_bao2018() -> dict:
    """The BAO data of Planck 2018 VI: 6dFGS, SDSS DR7 MGS and BOSS DR12 (BAO only).

    As Cobaya's likelihoods use them: 6dF r_s/D_V(0.106) = 0.336 ± 0.015 with
    r_s rescaled by 1.027369826; MGS as a tabulated χ² of α = D_V/(r_d × 4.2972)
    on [0.8005, 1.1985]; DR12 consensus D_M (147.78/r_d) and H (r_d/147.78).
    """
    dr12 = []
    for line in _require(BAO_DIR / "sdss_DR12Consensus_bao.dat").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            z, v, kind = line.split()
            dr12.append((float(z), float(v), kind.removeprefix("bao_")))
    cov = np.loadtxt(BAO_DIR / "BAO_consensus_covtot_dM_Hz.txt")
    mgs = np.loadtxt(_require(BAO_DIR / "sdss_MGS_prob.txt"))
    return {"sixdf": {"z": 0.106, "value": 0.336, "err": 0.015, "rs_rescale": 1.027369826},
            "mgs": {"z": 0.15, "chi2": mgs, "alpha_bounds": (0.8005, 1.1985),
                    "rs_rescale": 4.29720761315},
            "dr12": {"z": np.array([r[0] for r in dr12]), "value": np.array([r[1] for r in dr12]),
                     "kind": np.array([r[2] for r in dr12]), "cov": cov, "rs_fid": 147.78}}


def load_nufit(analysis: str = "TByes") -> dict:
    """NuFIT 6.0 one-dimensional Δχ² profiles for both orderings.

    Returns {ordering: {param: (x, dchi2)}} plus each ordering's offset from
    the global minimum (NuFIT's files measure Δχ² from the best ordering).
    Parameters: T12, T13 (sin²θ), DMS (log10 Δm²21/eV²), DMA (Δm²3ℓ in
    10⁻³ eV²: Δm²31 > 0 for NO, Δm²32 < 0 for IO).
    """
    import lzma

    cache = cfg.DERIVED_DIR / f"nufit_{analysis}_1d.npz"
    if cache.exists():
        with np.load(cache) as z:
            flat = {k: z[k] for k in z.files}
        out = {}
        for k, v in flat.items():
            o, prm, col = k.split("__")
            out.setdefault(o, {}).setdefault(prm, {})[col] = v
        return {o: {prm: (d["x"], d["y"]) for prm, d in v.items()} for o, v in out.items()}
    out = {}
    for order in ("NO", "IO"):
        path = _require(NUFIT_DIR / f"v60.release-{analysis}-{order}.txt.xz")
        blocks, cur = {}, None
        with lzma.open(path, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    cur = None
                    if "projection" in line:
                        tag = line.split(":")[0].strip("# ").strip().split()[0]
                        if "/" not in tag:
                            cur = tag
                            blocks[cur] = []
                    continue
                if cur and line.strip():
                    blocks[cur].append([float(v) for v in line.split()[:2]])
        out[order] = {k: (np.array(v)[:, 0], np.array(v)[:, 1]) for k, v in blocks.items()}
    flat = {f"{o}__{prm}__{c}": arr for o, d in out.items() for prm, (x, y) in d.items()
            for c, arr in (("x", x), ("y", y))}
    np.savez(cache, **flat)
    return out


def load_planck_margestats(key: str) -> dict:
    """GetDist marginalised statistics of a Planck 2018 grid run."""
    path = _require(PLANCK_GRID_DIR / (PLANCK_GRID_RUNS[key] + ".margestats"))
    return _parse_margestats(path)


def _parse_margestats(path: Path) -> dict:
    lines = path.read_text().splitlines()
    head = next(i for i, ln in enumerate(lines) if ln.startswith("parameter"))
    cols = lines[head].split()
    out = {}
    for ln in lines[head + 1:]:
        parts = ln.split()
        if len(parts) < 4:
            continue
        row = {}
        for c, v in zip(cols[1:], parts[1:]):
            try:
                row[c] = float(v)
            except ValueError:
                row[c] = v
        out[parts[0]] = row
    return out


def load_planck_minimum(key: str) -> dict:
    """A Planck grid best fit: parameters and each likelihood's chi^2 (CosmoMC .minimum)."""
    path = _require(PLANCK_GRID_DIR / (PLANCK_GRID_MINIMA[key] + ".minimum"))
    params, chi2 = {}, {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                params[parts[2]] = float(parts[1])
            except ValueError:
                pass
        elif len(parts) >= 3 and parts[1].startswith("chi2_"):
            try:
                chi2[parts[1]] = float(parts[0])
            except ValueError:
                pass
        elif line.strip().startswith("-log(Like)") and "=" in line:
            params["-logLike"] = float(line.split("=")[1])
    return {"params": params, "chi2": chi2, "text": path.read_text()}


def load_planck_chain(key: str) -> dict[str, np.ndarray]:
    """A Planck 2018 grid chain (CosmoMC format; burn-in already removed)."""
    run = PLANCK_GRID_DIR / PLANCK_GRID_CHAINS[key]
    names = ["weight", "minuslogpost"] + [ln.split()[0].rstrip("*") for ln in
                                          _require(Path(str(run) + ".paramnames")).read_text()
                                          .splitlines() if ln.strip()]
    parts = [np.loadtxt(Path(f"{run}_{i}.txt")) for i in (1, 2, 3, 4)]
    a = np.vstack(parts)
    return {n: a[:, j] for j, n in enumerate(names[: a.shape[1]])}


def load_act_chain(burn: float = 0.3) -> dict[str, np.ndarray]:
    """ACT DR6 P-ACT-LB N_eff chains (Cobaya format), burn-in removed."""
    parts = []
    for path in sorted(ACT_CHAIN_DIR.glob("*.[0-9].txt")):
        t = _read_cobaya_table(path)
        n = len(t["weight"])
        parts.append({k: v[int(burn * n):] for k, v in t.items()})
    if not parts:
        raise FileNotFoundError("ACT chain not cached; run fetch_data.py")
    keys = set.intersection(*(set(p) for p in parts))
    return {k: np.concatenate([p[k] for p in parts]) for k in keys}


def load_bao() -> dict:
    z, value, kind = [], [], []
    for line in _require(BAO_DIR / "desi_dr2_mean.txt").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            a, b, c = line.split()
            z.append(float(a)), value.append(float(b)), kind.append(c)
    cov = np.loadtxt(_require(BAO_DIR / "desi_dr2_cov.txt"))
    if cov.shape != (len(z), len(z)):
        raise ValueError("BAO covariance does not match the data vector")
    return {"z": np.array(z), "value": np.array(value), "kind": np.array(kind), "cov": cov}


def _read_header_table(path: Path, sep: str | None = None) -> dict[str, np.ndarray]:
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    names = lines[0].lstrip("#").replace(",", " ").split() if sep == "," else \
        lines[0].lstrip("#").split()
    rows = [(ln.split(",") if sep == "," else ln.split()) for ln in lines[1:]
            if not ln.startswith("#")]
    out = {}
    for i, name in enumerate(names):
        col = [r[i] for r in rows]
        try:
            out[name] = np.array(col, dtype=float)
        except ValueError:
            out[name] = np.array(col)
    return out


def _read_flat_cov(path: Path, n: int) -> np.ndarray:
    """A covariance stored as its size followed by n² numbers (Pantheon+, DESY5)."""
    cache = cfg.DERIVED_DIR / (path.parent.name + "_" + path.stem + ".npy")
    if cache.exists():
        return np.load(cache)
    import pandas as pd

    flat = pd.read_csv(path, header=None, sep=r"\s+", engine="c").to_numpy(float).ravel()
    flat = flat[~np.isnan(flat)]
    if int(flat[0]) != n or flat.size != 1 + n * n:
        raise ValueError(f"{path.name}: covariance does not match the table ({flat.size})")
    cov = flat[1:].reshape(n, n)
    np.save(cache, cov)
    return cov


def load_sn(name: str) -> dict:
    """A supernova sample prepared as DESI's (Cobaya's) likelihood reads it.

    * **Pantheon+**: m_b_corr with the STAT+SYS covariance, z_HD > 0.01 only.
    * **DES-SN5YR**: MU with the systematic covariance plus MUERR_FINAL² on the
      diagonal; no redshift cut.
    """
    folder = SN_DIR / name
    if name == "PantheonPlus":
        t = _read_header_table(_require(folder / "Pantheon+SH0ES.dat"))
        cov = _read_flat_cov(_require(folder / "Pantheon+SH0ES_STAT+SYS.cov"), t["zHD"].size)
        keep = t["zHD"] > 0.01
        return {"label": "Pantheon+", "z_hd": t["zHD"][keep], "z_hel": t["zHEL"][keep],
                "mu": t["m_b_corr"][keep], "cov": cov[np.ix_(keep, keep)]}
    if name == "DESY5":
        t = _read_header_table(_require(folder / "DES-SN5YR_HD.csv"), sep=",")
        n = t["zHD"].size
        cov = _read_flat_cov(_require(folder / "covsys_000.txt"), n).copy()
        cov[np.diag_indices(n)] += t["MUERR_FINAL"] ** 2
        return {"label": "DES-SN5YR", "z_hd": t["zHD"], "z_hel": t["zHEL"], "mu": t["MU"],
                "cov": cov}
    raise KeyError(name)


def _read_cobaya_table(path: Path) -> dict[str, np.ndarray]:
    with open(path) as f:
        names = f.readline().lstrip("#").split()
    a = np.atleast_2d(np.loadtxt(path))
    return {n: a[:, j] for j, n in enumerate(names)}


def load_desi_chain(run: str, *, burn: float = 0.3) -> dict[str, np.ndarray]:
    """One of DESI's runs (four MCMC chains, burn-in removed, weighted)."""
    folder = DESI_DIR / "cobaya" / run
    cache = folder / f"chains_burn{burn:.2f}.npz"
    if cache.exists():
        with np.load(cache) as z:
            return {k: z[k] for k in z.files}
    parts = []
    for path in sorted(folder.glob("chain.[0-9]*.txt")):
        t = _read_cobaya_table(path)
        n = len(t["weight"])
        parts.append({k: v[int(burn * n):] for k, v in t.items()})
    if not parts:
        raise FileNotFoundError(f"no chains cached in {folder}; run fetch_data.py")
    keys = set.intersection(*(set(p) for p in parts))
    chain = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    np.savez(cache, **chain)
    return chain


def _parse_bound(text: str) -> dict:
    """One cell of a getdist LaTeX table: '< 0.064', '3.05\\pm 0.01', '0.97^{+0.01}_{-0.02}'."""
    t = text.replace(" ", "")
    m = re.fullmatch(r"([<>])([-+0-9.eE]+)", t)
    if m:
        return {"upper" if m.group(1) == "<" else "lower": float(m.group(2))}
    m = re.fullmatch(r"([-+0-9.eE]+)\\pm([-+0-9.eE]+)", t)
    if m:
        return {"mean": float(m.group(1)), "plus": float(m.group(2)), "minus": float(m.group(2))}
    m = re.fullmatch(r"([-+0-9.eE]+)\^\{\+?([-+0-9.eE]+)\}_\{-?([-+0-9.eE]+)\}", t)
    if m:
        return {"mean": float(m.group(1)), "plus": float(m.group(2)), "minus": float(m.group(3))}
    return {"text": text.strip()}


def load_desi_margestats(run: str) -> dict[str, dict]:
    """DESI's getdist summary of a run: name -> {"68": bound, "95": bound}.

    DESI ships a LaTeX-style table with 40-character columns (parameter, 68%,
    95%), e.g. ``mnu  < 0.0289  < 0.0642``.
    """
    path = _require(DESI_DIR / "cobaya" / run / "chain.margestats")
    out = {}
    for ln in path.read_text().splitlines()[1:]:
        if not ln.strip():
            continue
        name = ln[:40].strip()
        out[name] = {"68": _parse_bound(ln[40:80]), "95": _parse_bound(ln[80:120])}
    return out


def load_desi_bestfit(run: str) -> dict[str, float]:
    """DESI's MAP point for a run (the best of its minimiser runs)."""
    best = None
    for path in sorted((DESI_DIR / "iminuit" / run).glob("**/bestfit.minimum.txt")):
        row = {k: float(v[0]) for k, v in _read_cobaya_table(path).items()}
        if best is None or row["minuslogpost"] < best["minuslogpost"]:
            best = row
    if best is None:
        raise FileNotFoundError(f"no DESI best fit cached for {run}")
    return best


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refresh", action="store_true", help="download again even if cached")
    ap.add_argument("--no-chains", action="store_true", help="skip the MCMC chains (answer keys)")
    ap.add_argument("--no-act-lensing", action="store_true",
                    help="skip ACT's 361 MB lensing archive (its rungs are then skipped)")
    args = ap.parse_args(argv)
    fetch_all(refresh=args.refresh, chains=not args.no_chains,
              act_lensing=not args.no_act_lensing)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
