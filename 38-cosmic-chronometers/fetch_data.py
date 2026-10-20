"""Download and cache everything the analysis needs, and load it back.

Sources (all public, no account needed):

1. **Moresco's cosmic-chronometer primer** on GitLab (``mmoresco/CCcovariance``):
   his own 15 H(z) measurements to full precision, and the redshift-dependent
   systematic errors of Moresco et al. (2020) from which the covariance matrix
   is built. The primer's fit notebook is kept as an answer key.
2. **Three arXiv LaTeX sources**, for the tables they contain: the 32-point
   compilation of Moresco et al. (2022, Living Reviews in Relativity), the
   30-point table behind Moresco et al.'s (2016) transition redshift, and the
   D4000-redshift medians and D4000-age calibration of Moresco et al. (2012).
   Only the one ``.tex`` file needed is kept from each tarball.
3. **DESI DR2 BAO** distances (DESI Collaboration 2025), as packaged for the
   Cobaya sampler on GitHub.

It also reads, without copying, the cached outputs of earlier atlas projects:
D1's supernova fits, D3's CMB fit, and J3's SDSS galaxy table.

Everything downloaded is recorded in ``data/manifest.json`` (URL, bytes, SHA-256,
time). Nothing already cached is fetched again, so once the cache is warm the
whole project runs offline.

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import tarfile
import gzip
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import config as cfg

USER_AGENT = "Mozilla/5.0 (X11; Linux x86_64) astrophysics-atlas/D8"


# ---------------------------------------------------------------------------
# Downloading
# ---------------------------------------------------------------------------


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _load_manifest() -> dict:
    if cfg.MANIFEST_PATH.exists():
        return json.loads(cfg.MANIFEST_PATH.read_text())
    return {}


def _get(url: str, *, timeout: float = 180.0, retries: int = 3) -> bytes:
    last = None
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as response:
                return response.read()
        except OSError as err:  # URLError, timeouts, resets
            last = err
    raise RuntimeError(f"could not fetch {url}: {last}")


def _record(manifest: dict, key: str, url: str, data: bytes, **extra) -> None:
    manifest[key] = {
        "url": url,
        "bytes": len(data),
        "sha256": _sha256(data),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        **extra,
    }


def _fetch_file(url: str, dest: Path, manifest: dict, *, refresh: bool) -> None:
    key = str(dest.relative_to(cfg.DATA_DIR))
    if dest.exists() and not refresh:
        print(f"  cached   {key}")
        return
    print(f"  fetching {key}  <-  {url}")
    data = _get(url)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(data)
    tmp.replace(dest)  # atomic: a half-written file is never mistaken for a cached one
    _record(manifest, key, url, data)


def _fetch_arxiv_member(arxiv_id: str, member: str, dest: Path, manifest: dict,
                        *, refresh: bool) -> None:
    """Download an arXiv source tarball and keep one file from it."""
    key = str(dest.relative_to(cfg.DATA_DIR))
    if dest.exists() and not refresh:
        print(f"  cached   {key}")
        return
    url = f"https://arxiv.org/e-print/{arxiv_id}"
    print(f"  fetching {key}  <-  {url} ({member})")
    blob = _get(url)
    try:
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tar:
            text = tar.extractfile(member).read()
    except tarfile.ReadError:  # a single gzipped file rather than a tarball
        text = gzip.decompress(blob)
    dest.write_bytes(text)
    _record(manifest, key, url, text, tarball_bytes=len(blob),
            tarball_sha256=_sha256(blob), member=member)


def fetch_all(*, refresh: bool = False) -> dict:
    manifest = _load_manifest()
    print("Moresco's cosmic-chronometer primer (GitLab)")
    for local, path in cfg.CC_REPO_FILES.items():
        _fetch_file(cfg.CC_REPO_RAW + path, cfg.DATA_DIR / local, manifest, refresh=refresh)
    print("Tables from arXiv LaTeX sources")
    for local, (arxiv_id, member) in cfg.ARXIV_SOURCES.items():
        _fetch_arxiv_member(arxiv_id, member, cfg.SOURCE_DIR / local, manifest,
                            refresh=refresh)
    print("DESI DR2 BAO (Cobaya packaging)")
    for local, remote in cfg.DESI_FILES.items():
        _fetch_file(cfg.DESI_RAW + remote, cfg.DATA_DIR / local, manifest, refresh=refresh)
    cfg.MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))

    print("Earlier atlas projects (read in place, not copied)")
    for label, path in (("D1 supernova fits", cfg.D1_RESULTS), ("D3 CMB fit", cfg.D3_FIT),
                        ("J3 SDSS galaxies", cfg.J3_GALAXIES)):
        state = "found" if path.exists() else "MISSING (fallback values will be used)"
        print(f"  {label:20s} {state}: {path.relative_to(cfg.ATLAS_DIR)}")
    return manifest


# ---------------------------------------------------------------------------
# Loaders -- read only from the cache, never from the network
# ---------------------------------------------------------------------------


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


def load_moresco_table(models: str = "BC03") -> dict[str, np.ndarray]:
    """Moresco's own 15 measurements (GitLab), for BC03 or M11 models.

    Columns: z, H, err (the diagonal error: statistical plus metallicity, and for
    BC03 the split into those two parts), ref.
    """
    path = _require(cfg.DATA_DIR / f"HzTable_MM_{models}.dat")
    rows = [line.split(",") for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    out = {"z": np.array([float(r[0]) for r in rows]),
           "H": np.array([float(r[1]) for r in rows]),
           "err": np.array([float(r[2]) for r in rows]),
           "ref": np.array([r[-1].strip() for r in rows])}
    if models == "BC03":
        out["stat"] = np.array([float(r[3]) for r in rows])
        out["met"] = np.array([float(r[4]) for r in rows])
    if out["z"].size != 15:
        raise ValueError(f"expected 15 Moresco measurements, found {out['z'].size}")
    return out


def load_systematics() -> dict[str, np.ndarray]:
    """Moresco et al. (2020) systematic errors on H(z), per cent, vs redshift."""
    path = _require(cfg.DATA_DIR / "data_MM20.dat")
    z, imf, stlib, mod, mod_ooo = np.loadtxt(path, comments="#", unpack=True)
    return {"z": z, "imf": imf, "stlib": stlib, "mod": mod, "mod_ooo": mod_ooo}


def _table_block(tex: str, label: str) -> str:
    """The LaTeX of the table environment that carries ``\\label{label}``."""
    idx = tex.index(f"\\label{{{label}}}")
    start = tex.rfind("\\begin{table", 0, idx)
    end = tex.find("\\end{table", idx)
    return tex[start:end]


_NUM = r"([0-9]+(?:\.[0-9]*)?)"


def load_lrr_compilation() -> dict[str, np.ndarray]:
    """The 32-point table of Moresco et al. (2022), parsed from its LaTeX.

    Returns z, H, err (the diagonal error as tabulated), method (D/F/L) and the
    citation key of the original measurement.
    """
    tex = _require(cfg.SOURCE_DIR / "moresco2022_lrr.tex").read_text(errors="replace")
    block = _table_block(tex, "tab:CC1")
    pat = re.compile(rf"^\s*{_NUM}\s*&\s*{_NUM}\s*&\s*{_NUM}\s*&\s*([DFL])\s*&"
                     r"\s*\\cite\{(\w+)\}", re.M)
    rows = pat.findall(block)
    if len(rows) != 32:
        raise ValueError(f"expected 32 rows in the compilation, parsed {len(rows)}")
    return {"z": np.array([float(r[0]) for r in rows]),
            "H": np.array([float(r[1]) for r in rows]),
            "err": np.array([float(r[2]) for r in rows]),
            "method": np.array([r[3] for r in rows]),
            "ref": np.array([r[4] for r in rows])}


def load_m16_table() -> dict[str, np.ndarray]:
    """The 30-point table of Moresco et al. (2016), parsed from its LaTeX."""
    tex = _require(cfg.SOURCE_DIR / "moresco2016_boss.tex").read_text(errors="replace")
    block = _table_block(tex, "tab:Hzall")
    pat = re.compile(rf"^\s*{_NUM}\s*&\s*{_NUM}\s*&\s*{_NUM}\s*&\s*(.+?)\\\\", re.M)
    rows = pat.findall(block)
    if len(rows) != 30:
        raise ValueError(f"expected 30 rows in the 2016 table, parsed {len(rows)}")
    ref = [re.sub(r"\\cite\{(\w+)\}", r"\1", r[3]).strip() for r in rows]
    return {"z": np.array([float(r[0]) for r in rows]),
            "H": np.array([float(r[1]) for r in rows]),
            "err": np.array([float(r[2]) for r in rows]),
            "ref": np.array(ref)}


def load_m12_d4000() -> dict:
    """Moresco et al. (2012): D4000_n medians per redshift bin, and A(Z).

    Returns {"bins": {"low mass": rows, "high mass": rows}, "A": {model: rows}},
    where a bin row is (z, D4000, sigma, Z/Zsun, sigma_Z, log M, N) and an A row
    is (Z/Zsun, A_high, err_high, A_low, err_low) in Gyr^-1.
    """
    tex = _require(cfg.SOURCE_DIR / "moresco2012_d4000.tex").read_text(errors="replace")
    block = _table_block(tex, "tab:D4000med")
    bins: dict[str, list] = {}
    current = None
    row_pat = re.compile(rf"^\s*{_NUM}\s*&\s*{_NUM}\s*&\s*{_NUM}\s*&\s*{_NUM}\s*&\s*{_NUM}"
                         rf"\s*&\s*\$?{_NUM}\$?\s*&\s*([0-9]+)", re.M)
    for line in block.splitlines():
        if "low mass" in line:
            current = "low mass"
            bins[current] = []
        elif "high mass" in line:
            current = "high mass"
            bins[current] = []
        elif current:
            m = row_pat.match(line)
            if m:
                bins[current].append(tuple(float(v) for v in m.groups()))
    a_block = _table_block(tex, "tab:Zslope")
    a_pat = re.compile(r"(BC03|MaStro)\s*&\s*\$A\(Z/Z_\{\\odot\}=([0-9.]+)\)\$\s*&\s*"
                       rf"\$\s*{_NUM}\\pm\s*{_NUM}\$\s*&\s*\${_NUM}\\pm\s*{_NUM}\$")
    a_rows: dict[str, list] = {}
    for model, *vals in a_pat.findall(a_block):
        a_rows.setdefault(model, []).append(tuple(float(v) for v in vals))
    if len(bins.get("low mass", [])) != 13 or len(bins.get("high mass", [])) != 14:
        raise ValueError("unexpected shape of the D4000 table: "
                         f"{ {k: len(v) for k, v in bins.items()} }")
    if any(len(v) != 3 for v in a_rows.values()) or len(a_rows) != 2:
        raise ValueError(f"unexpected A(Z) table: {a_rows}")
    return {"bins": bins, "A": a_rows}


def load_desi_dr2() -> dict:
    """DESI DR2 BAO: the 13-element data vector and its covariance."""
    mean_path = _require(cfg.DATA_DIR / "desi_dr2_bao_mean.txt")
    rows = [line.split() for line in mean_path.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    cov = np.loadtxt(_require(cfg.DATA_DIR / "desi_dr2_bao_cov.txt"))
    out = {"z": np.array([float(r[0]) for r in rows]),
           "value": np.array([float(r[1]) for r in rows]),
           "quantity": np.array([r[2] for r in rows]),
           "cov": cov}
    if cov.shape != (out["z"].size,) * 2:
        raise ValueError("DESI covariance does not match the data vector")
    return out


def load_d1() -> dict:
    """D1's flat-LCDM supernova fit: H0, Omega_m, transition redshift, ...

    Values are (mean, std). Falls back to the numbers in D1's README when its
    cache is absent (data/ folders are not committed).
    """
    out = dict(cfg.D1_FALLBACK)
    out["source"] = "D1 README (fallback)"
    if cfg.D1_RESULTS.exists():
        r = json.loads(cfg.D1_RESULTS.read_text())
        cal = r["calibrated"]
        out.update({
            "H0": (cal["H0"]["mean"], cal["H0"]["std"]),
            "Om_calibrated": (cal["Om"]["mean"], cal["Om"]["std"]),
            "Om_flat": (r["flat_lcdm"]["statsys"]["Om"]["mean"],
                        r["flat_lcdm"]["statsys"]["Om"]["std"]),
            "z_transition": (r["flat_lcdm"]["z_transition"]["mean"],
                             r["flat_lcdm"]["z_transition"]["std"]),
            "q0_flat": (r["flat_lcdm"]["q0"]["mean"], r["flat_lcdm"]["q0"]["std"]),
            "w_flat": (r["flat_wcdm"]["w"]["mean"], r["flat_wcdm"]["w"]["std"]),
            "age_gyr": (cal["age_gyr"]["mean"], cal["age_gyr"]["std"]),
            "source": str(cfg.D1_RESULTS.relative_to(cfg.ATLAS_DIR)),
        })
    return out


def load_d3() -> dict:
    """D3's CAMB fit to the Planck TT spectrum: H0, Omega_m, r_drag, age."""
    out = dict(cfg.D3_FALLBACK)
    out["source"] = "D3 README (fallback)"
    if cfg.D3_FIT.exists():
        d = json.loads(cfg.D3_FIT.read_text())["derived"]
        out.update({
            "H0": (d["H0"]["value"], d["H0"]["error"]),
            "Om": (d["omegam"]["value"], d["omegam"]["error"]),
            "rdrag": (d["rdrag"]["value"], d["rdrag"]["error"]),
            "age_gyr": (d["age"]["value"], d["age"]["error"]),
            "source": str(cfg.D3_FIT.relative_to(cfg.ATLAS_DIR)),
        })
    return out


def load_j3_galaxies() -> dict[str, np.ndarray] | None:
    """J3's SDSS Main Galaxy Sample with MPA-JHU measurements, if cached."""
    if not cfg.J3_GALAXIES.exists():
        return None
    keep = ("z", "lgm_tot_p50", "d4000_n", "d4000_n_err", "h_alpha_eqw",
            "specsfr_tot_p50", "dered_u", "dered_r", "reliable")
    with np.load(cfg.J3_GALAXIES) as d:
        return {k: d[k] for k in keep}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    args = parser.parse_args()
    fetch_all(refresh=args.refresh)
    total = sum(v["bytes"] for v in _load_manifest().values())
    print(f"\nCache holds {total / 1e6:.2f} MB of downloaded files.")
