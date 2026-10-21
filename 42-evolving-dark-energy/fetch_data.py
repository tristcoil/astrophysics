"""Download and cache everything the analysis needs, and load it back.

Four groups of files, all public:

1. **DESI DR2 BAO** -- 13 numbers and their covariance: the sound-horizon ruler
   measured at seven redshifts. Plus DESI DR1 and SDSS for comparison.
2. **Four supernova samples** -- Pantheon+, Union3, DES-SN5YR and its 2025
   recalibration DES-Dovekie, each with its full covariance (or its inverse).
3. **Planck's early-universe chains** (Lemos & Lewis 2023), from which the
   compressed CMB prior on (theta_*, omega_b, omega_bc) is computed.
4. **DESI's own published results** -- best fits and MCMC chains -- used only
   as an answer key, never as an input.

Everything lands in ``data/`` with a manifest recording each file's URL, size,
SHA-256 and retrieval time. Nothing already cached is downloaded again, so
after one run the whole project works offline.

    python fetch_data.py            # everything (~0.55 GB, mostly DESI's chains)
    python fetch_data.py --no-chains   # skip DESI's chains (~0.15 GB)
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from config import (BAO_FILES, BAO_REPO_RAW, DATA_DIR, DESI_BESTFIT_MODELS, DESI_CHAINS,
                    LSS_FILES, LSS_URL,
                    DESI_CHECKSUMS, DESI_COMBOS, DESI_RESULTS_URL, MANIFEST_PATH,
                    PLANCK_EARLY_FILES, PLANCK_EARLY_RAW, PLANCK_EARLY_ROOT, SN_FILES,
                    SN_REPO_RAW)

BAO_DIR = DATA_DIR / "bao"
SN_DIR = DATA_DIR / "sn"
PLANCK_EARLY_DIR = DATA_DIR / "planck_early_lcdm"
DESI_DIR = DATA_DIR / "desi_results"
LSS_DIR = DATA_DIR / "desi_dr1_lss"

USER_AGENT = "astrophysics-atlas-D5/1.0 (educational; python-urllib)"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _load_manifest() -> dict:
    return json.loads(MANIFEST_PATH.read_text()) if MANIFEST_PATH.exists() else {}


def _save_manifest(manifest: dict) -> None:
    """Merge into the manifest on disk (another fetch may have run meanwhile)."""
    merged = _load_manifest()
    merged.update(manifest)
    tmp = MANIFEST_PATH.with_suffix(".json.part")
    tmp.write_text(json.dumps(merged, indent=2, sort_keys=True))
    tmp.replace(MANIFEST_PATH)


def download(url: str, dest: Path, manifest: dict, *, refresh: bool = False,
             expected_sha256: str | None = None, retries: int = 4) -> Path:
    """Fetch one file unless cached; record provenance; verify a checksum if given."""
    key = str(dest.relative_to(DATA_DIR))
    if dest.exists() and not refresh:
        if key not in manifest:  # cached by an earlier step: record it, dated by the file
            mtime = datetime.fromtimestamp(dest.stat().st_mtime, timezone.utc)
            manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": _sha256(dest),
                             "retrieved_utc": mtime.isoformat(timespec="seconds")}
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as out:
                shutil.copyfileobj(r, out, 1 << 20)
            break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            if attempt == retries - 1:
                raise
            print(f"    retrying {dest.name} after {err!r}")
            time.sleep(2 + 3 * attempt)
    sha = _sha256(tmp)
    if expected_sha256 and sha != expected_sha256:
        tmp.unlink()
        raise IOError(f"{dest.name}: checksum mismatch (got {sha[:12]}..., "
                      f"expected {expected_sha256[:12]}...)")
    tmp.replace(dest)  # atomic: a half-written file is never mistaken for a cached one
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha,
                     "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    print(f"  fetched {key}  ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


# ---------------------------------------------------------------------------
# DESI's published results: which files to fetch
# ---------------------------------------------------------------------------


def desi_checksums(manifest: dict, refresh: bool = False) -> dict[str, str]:
    """DESI's own SHA-256 list for its results release: {relative path: sha256}."""
    path = download(DESI_RESULTS_URL + DESI_CHECKSUMS, DESI_DIR / DESI_CHECKSUMS, manifest,
                    refresh=refresh)
    sums = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 2:
            sums[parts[1].removeprefix("./")] = parts[0]
    return sums


def desi_wanted_files(sums: dict[str, str], *, chains: bool) -> list[str]:
    """Best-fit tables for our model/dataset combinations, plus chains if asked."""
    wanted = []
    for rel in sums:
        parts = rel.split("/")
        if parts[0] == "iminuit" and parts[1] in DESI_BESTFIT_MODELS and \
                parts[2] in DESI_COMBOS and parts[-1] == "bestfit.minimum.txt":
            wanted.append(rel)
        if parts[0] == "cobaya" and "/".join(parts[1:3]) in DESI_CHAINS:
            name = parts[-1]
            if name in ("chain.updated.yaml", "chain.input.yaml") or (
                    chains and name.startswith("chain.") and name.endswith(".txt")
                    and name.split(".")[1].isdigit()):
                wanted.append(rel)
    return sorted(wanted)


# ---------------------------------------------------------------------------
# Fetch everything
# ---------------------------------------------------------------------------


def fetch_all(*, refresh: bool = False, chains: bool = True) -> dict:
    manifest = _load_manifest()
    try:
        print("DESI DR2 BAO (+ DR1, SDSS)")
        for local, repo_path in BAO_FILES.items():
            download(BAO_REPO_RAW + repo_path, BAO_DIR / local, manifest, refresh=refresh)

        print("Type Ia supernova samples")
        for folder, files in SN_FILES.items():
            for local, repo_path in files.items():
                url = SN_REPO_RAW + urllib.parse.quote(repo_path)  # '+' in Pantheon+ names
                download(url, SN_DIR / folder / local, manifest, refresh=refresh)

        print("Planck PR4 early-universe chains (Lemos & Lewis 2023)")
        for name in PLANCK_EARLY_FILES:
            download(PLANCK_EARLY_RAW + name, PLANCK_EARLY_DIR / name, manifest,
                     refresh=refresh)

        print("DESI's published best fits" + (" and chains" if chains else ""))
        sums = desi_checksums(manifest, refresh=refresh)
        for rel in desi_wanted_files(sums, chains=chains):
            download(DESI_RESULTS_URL + rel, DESI_DIR / rel, manifest, refresh=refresh,
                     expected_sha256=sums[rel])
    finally:
        _save_manifest(manifest)
    return manifest


# ---------------------------------------------------------------------------
# Loaders -- read only from the cache, never from the network
# ---------------------------------------------------------------------------


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


def load_bao(release: str = "dr2") -> dict:
    """A DESI BAO data vector: redshifts, values, observable names, covariance.

    Observables are ``DV_over_rs``, ``DM_over_rs`` and ``DH_over_rs``: the
    volume-averaged, transverse and line-of-sight distances in units of the
    sound horizon at the drag epoch.
    """
    mean = _require(BAO_DIR / f"desi_{release}_mean.txt")
    z, value, kind = [], [], []
    for line in mean.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            a, b, c = line.split()
            z.append(float(a)), value.append(float(b)), kind.append(c)
    cov = np.loadtxt(_require(BAO_DIR / f"desi_{release}_cov.txt"))
    if cov.shape != (len(z), len(z)):
        raise ValueError(f"BAO covariance shape {cov.shape} does not match {len(z)} values")
    return {"z": np.array(z), "value": np.array(value), "kind": np.array(kind), "cov": cov}


def load_sdss_bao() -> list[dict]:
    """SDSS BOSS DR12 / eBOSS DR16 consensus BAO points, for display only.

    Each file lists (z, value, observable) with D_M and D_H in units of r_d
    (DR16) or D_M * r_d,fid / r_d and H * r_d / r_d,fid (DR12). Returned as-is;
    the figure module converts them.
    """
    out = []
    for stem, cov in (("sdss_dr12_lrg_mean.dat", "sdss_dr12_lrg_cov.txt"),
                      ("sdss_dr16_lrg_mean.dat", "sdss_dr16_lrg_cov.txt"),
                      ("sdss_dr16_qso_mean.txt", "sdss_dr16_qso_cov.txt")):
        rows = [ln.split() for ln in _require(BAO_DIR / stem).read_text().splitlines()
                if ln.strip() and not ln.startswith("#")]
        out.append({"name": stem.split("_mean")[0],
                    "z": np.array([float(r[0]) for r in rows]),
                    "value": np.array([float(r[1]) for r in rows]),
                    "kind": np.array([r[2] for r in rows]),
                    "cov": np.loadtxt(_require(BAO_DIR / cov))})
    return out


def _read_header_table(path: Path, sep: str | None = None) -> dict[str, np.ndarray]:
    """Whitespace- or comma-separated table with a (possibly '#'-prefixed) header."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    header = lines[0].lstrip("#").strip()
    names = [c.strip() for c in (header.split(sep) if sep else header.split())]
    rows = [[v.strip() for v in (ln.split(sep) if sep else ln.split())] for ln in lines[1:]
            if not ln.startswith("#")]
    table = {}
    for i, name in enumerate(names):
        # Union3's header names 19 columns but its rows hold 18 values; the
        # missing trailing column (biascor) is unused, as in Cobaya's reader.
        col = [r[i] if i < len(r) else "nan" for r in rows]
        try:
            table[name] = np.array(col, dtype=float)
        except ValueError:
            table[name] = np.array(col)
    return table


def _read_flat_cov(path: Path, n: int) -> np.ndarray:
    """A covariance stored as N then N*N values (one per line), Cobaya's format."""
    import pandas as pd

    values = pd.read_csv(path, header=None, sep=r"\s+").to_numpy(float).ravel()
    if values[0] == n and values.size == 1 + n * n:
        values = values[1:]
    if values.size != n * n:
        raise ValueError(f"{path.name}: {values.size} values for a {n} x {n} matrix")
    cov = values.reshape(n, n)
    return 0.5 * (cov + cov.T)


def load_sn(name: str) -> dict:
    """One supernova sample, prepared exactly as DESI's (Cobaya's) likelihood reads it.

    Returns redshifts (``z_hd``: cosmological, CMB frame and corrected for
    peculiar velocities; ``z_hel``: heliocentric), standardised magnitudes or
    distance moduli ``mu`` (each up to an unknown constant offset, which the fit
    removes), and either ``cov`` or ``invcov``:

    * **Pantheon+**: m_b_corr with the STAT+SYS covariance, z_HD > 0.01 only.
    * **Union3**: 22 binned distance moduli from the UNITY1.5 Bayesian model,
      with their covariance; heliocentric = CMB-frame redshift.
    * **DES-SN5YR**: MU with the systematic covariance plus MUERR_FINAL^2 on the
      diagonal; no redshift cut.
    * **DES-Dovekie**: MU with the inverse of the total covariance, as released.
    """
    folder = SN_DIR / name
    if name == "PantheonPlus":
        t = _read_header_table(_require(folder / "Pantheon+SH0ES.dat"))
        cov = _read_flat_cov(_require(folder / "Pantheon+SH0ES_STAT+SYS.cov"), t["zHD"].size)
        keep = t["zHD"] > 0.01
        return {"name": name, "z_hd": t["zHD"][keep], "z_hel": t["zHEL"][keep],
                "mu": t["m_b_corr"][keep], "cov": cov[np.ix_(keep, keep)], "n_all": keep.size}
    if name == "Union3":
        t = _read_header_table(_require(folder / "lcparam_full.txt"))
        cov = _read_flat_cov(_require(folder / "mag_covmat.txt"), t["zcmb"].size)
        return {"name": name, "z_hd": t["zcmb"], "z_hel": t["zcmb"], "mu": t["mb"],
                "cov": cov, "n_all": t["zcmb"].size}
    if name == "DESY5":
        t = _read_header_table(_require(folder / "DES-SN5YR_HD.csv"), sep=",")
        n = t["zHD"].size
        cov = _read_flat_cov(_require(folder / "covsys_000.txt"), n)
        cov[np.diag_indices(n)] += t["MUERR_FINAL"] ** 2
        return {"name": name, "z_hd": t["zHD"], "z_hel": t["zHEL"], "mu": t["MU"],
                "cov": cov, "n_all": n, "idsurvey": t.get("IDSURVEY")}
    if name == "DES-Dovekie":
        t = _read_header_table(_require(folder / "DES-Dovekie_HD.csv"), sep=",")
        cols = {k.lower(): v for k, v in t.items()}
        d = np.load(_require(folder / "covtot_inv_000.npz"))
        n = int(d[d.files[0]][0])
        inv = np.zeros((n, n))
        inv[np.triu_indices(n)] = d[d.files[1]]
        inv = inv + np.triu(inv, 1).T
        if n != cols["zhd"].size:
            raise ValueError(f"Dovekie inverse covariance is {n} x {n} for {cols['zhd'].size} SNe")
        if np.any(cols["zhd"] <= 0.01):  # Cobaya would mask these; none exist
            raise ValueError("Dovekie sample has z <= 0.01 supernovae")
        return {"name": name, "z_hd": cols["zhd"], "z_hel": cols["zhel"], "mu": cols["mu"],
                "invcov": inv, "n_all": n, "idsurvey": cols.get("idsurvey")}
    raise KeyError(name)


SN_SAMPLES = ("PantheonPlus", "Union3", "DESY5", "DES-Dovekie")
SN_LABELS = {"PantheonPlus": "Pantheon+", "Union3": "Union3", "DESY5": "DESY5",
             "DES-Dovekie": "DES-Dovekie"}


def load_planck_early_chains(burn: float) -> dict[str, np.ndarray]:
    """Lemos & Lewis's four chains, burn-in removed, concatenated; with weights."""
    parts = []
    names = None
    for i in (1, 2, 3, 4):
        path = _require(PLANCK_EARLY_DIR / f"{PLANCK_EARLY_ROOT}.{i}.txt")
        with open(path) as f:
            names = f.readline().lstrip("#").split()
        a = np.loadtxt(path)
        parts.append(a[int(burn * len(a)):])
    a = np.vstack(parts)
    return {n: a[:, j] for j, n in enumerate(names)}


def _read_cobaya_table(path: Path) -> dict[str, np.ndarray]:
    with open(path) as f:
        names = f.readline().lstrip("#").split()
    a = np.atleast_2d(np.loadtxt(path))
    return {n: a[:, j] for j, n in enumerate(names)}


def load_desi_bestfits() -> dict[tuple[str, str], dict]:
    """DESI's MAP points: {(model, combo label): best of all minimiser runs}.

    Several folders hold more than one minimiser run (different starting
    covariances, reruns); the run with the lowest -log posterior is the MAP.
    """
    out: dict[tuple[str, str], dict] = {}
    root = DESI_DIR / "iminuit"
    for path in sorted(root.glob("*/*/**/bestfit.minimum.txt")):
        rel = path.relative_to(root).parts
        model, combo = rel[0], rel[1]
        if combo not in DESI_COMBOS:
            continue
        row = {k: float(v[0]) for k, v in _read_cobaya_table(path).items()}
        row["_file"] = str(path.relative_to(DATA_DIR))
        key = (model, DESI_COMBOS[combo])
        if key not in out or row["minuslogpost"] < out[key]["minuslogpost"]:
            out[key] = row
    return out


def load_desi_chain(model: str, combo: str, *, burn: float = 0.3) -> dict[str, np.ndarray]:
    """One of DESI's published chains (4 MCMC chains, burn-in removed, weighted).

    Converted once to a compact .npz cache so later loads take milliseconds.
    """
    folder = DESI_DIR / "cobaya" / model / combo
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


def fetch_lss(*, refresh: bool = False) -> dict:
    """DESI DR1 LRG clustering catalogues and randoms (~3.2 GB), for our own BAO."""
    manifest = _load_manifest()
    print("DESI DR1 LRG clustering catalogues (our own BAO measurement)")
    try:
        for name in LSS_FILES:
            download(LSS_URL + name, LSS_DIR / name, manifest, refresh=refresh)
    finally:
        _save_manifest(manifest)
    return manifest


def load_lss(cap: str, kind: str, index: int | None = None) -> dict[str, np.ndarray]:
    """Positions, redshifts and weights from one LRG catalogue (data or random).

    The clustering weight of each object is WEIGHT (completeness, imaging
    systematics, redshift failures) times WEIGHT_FKP, the Feldman-Kaiser-Peacock
    weight 1/(1 + n P0) that balances shot noise against cosmic variance.
    """
    from astropy.io import fits

    name = (f"LRG_{cap}_clustering.dat.fits" if kind == "data"
            else f"LRG_{cap}_{index}_clustering.ran.fits")
    with fits.open(_require(LSS_DIR / name), memmap=True) as hdul:
        d = hdul[1].data
        return {"ra": np.asarray(d["RA"], float), "dec": np.asarray(d["DEC"], float),
                "z": np.asarray(d["Z"], float),
                "w": np.asarray(d["WEIGHT"], float) * np.asarray(d["WEIGHT_FKP"], float)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    parser.add_argument("--no-chains", action="store_true",
                        help="skip DESI's MCMC chains (answer key for posteriors only)")
    parser.add_argument("--lss", action="store_true",
                        help="also fetch DESI DR1 LRG catalogues (~3.2 GB) for our own BAO")
    args = parser.parse_args(argv)
    if args.lss:
        fetch_lss(refresh=args.refresh)
    manifest = fetch_all(refresh=args.refresh, chains=not args.no_chains)
    total = sum(v["bytes"] for v in manifest.values())
    print(f"\n{len(manifest)} files, {total / 1e6:.0f} MB cached in {DATA_DIR}")
    for name in SN_SAMPLES:
        sn = load_sn(name)
        print(f"  {name:<12} {sn['z_hd'].size:5d} SNe, z = {sn['z_hd'].min():.4f}-"
              f"{sn['z_hd'].max():.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
