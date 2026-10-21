"""Download and cache everything the analysis needs, and load it back.

What this project reads, all public:

1. **Pantheon+** supernovae (D1's data): 1,701 light curves + full covariance.
2. **DESI DR2 BAO** (D5's data): 13 distance ratios at seven redshifts.
3. **Planck PR4 early-universe chains** (Lemos & Lewis 2023), from which the
   compressed CMB prior is rebuilt.
4. **DESI's released curvature runs** ("base_omegak"): MCMC chains for BAO
   alone, the full CMB alone and the two together, and best fits including
   BAO + Pantheon+. An answer key only: never an input to a fit.
5. Small literature tables: the 32 cosmic-chronometer H(z) points (Moresco
   et al. 2022, as in D8) and the most distant quasars (Yang et al. 2021, I3).

Files that an earlier atlas project already downloaded are **hard-linked** from
its cache instead of fetched again (same bytes, no extra disk), and the link is
recorded in the manifest. Everything lands in ``data/`` with each file's URL,
size, SHA-256 and retrieval time; nothing cached is fetched twice, so after one
run the whole project works offline.

    python fetch_data.py              # everything (~0.15 GB new, mostly DESI's chains)
    python fetch_data.py --no-chains  # skip DESI's chains (the answer key)
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import tarfile
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import config as cfg


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


def _get(url: str, retries: int = 4, timeout: int = 300) -> bytes:
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": cfg.USER_AGENT})
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
            entry = entries.get(str(sibling.relative_to(parent)))
            return entry.get("sha256") if isinstance(entry, dict) else None
        if parent.name == "data":
            break
    return None


def _link_from_sibling(dest: Path, manifest: dict, url: str) -> bool:
    """Hard-link an identical file from an earlier project's cache, if it has one."""
    src = cfg.SIBLING_COPIES.get(dest)
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
    except OSError:  # another filesystem: copy instead
        shutil.copy2(src, dest)
    key = str(dest.relative_to(cfg.DATA_DIR))
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha,
                     "linked_from": str(src.relative_to(cfg.ROOT_DIR)), "retrieved_utc": _now()}
    print(f"  linked   {key}  (from {src.relative_to(cfg.ROOT_DIR).parts[0]}, no download)")
    return True


def download(url: str, dest: Path, manifest: dict, *, refresh: bool = False,
             expected_sha256: str | None = None) -> Path:
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
    tmp.write_bytes(blob)
    sha = _sha256(tmp)
    if expected_sha256 and sha != expected_sha256:
        tmp.unlink()
        raise IOError(f"{dest.name}: checksum mismatch (got {sha[:12]}, expected "
                      f"{expected_sha256[:12]})")
    tmp.replace(dest)  # atomic: a half-written file is never mistaken for a cached one
    manifest[key] = {"url": url, "bytes": dest.stat().st_size, "sha256": sha,
                     "retrieved_utc": _now()}
    print(f"  fetched  {key}  ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def _fetch_arxiv_member(url: str, member: str, dest: Path, manifest: dict,
                        *, refresh: bool) -> None:
    """One file out of an arXiv source tarball (or the single gzipped file)."""
    key = str(dest.relative_to(cfg.DATA_DIR))
    if dest.exists() and not refresh:
        return
    if not refresh and _link_from_sibling(dest, manifest, url):
        return
    blob = _get(url)
    try:
        with tarfile.open(fileobj=io.BytesIO(blob), mode="r:*") as tar:
            text = tar.extractfile(member).read()
    except tarfile.ReadError:
        text = gzip.decompress(blob)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(text)
    manifest[key] = {"url": url, "member": member, "bytes": len(text),
                     "sha256": hashlib.sha256(text).hexdigest(),
                     "tarball_sha256": hashlib.sha256(blob).hexdigest(), "retrieved_utc": _now()}
    print(f"  fetched  {key}  (from the arXiv source of {url.rsplit('/', 1)[-1]})")


# ---------------------------------------------------------------------------
# DESI's curvature runs
# ---------------------------------------------------------------------------


def desi_checksums(manifest: dict, refresh: bool = False) -> dict[str, str]:
    """DESI's own SHA-256 list for the results release: {relative path: sha256}."""
    path = download(cfg.DESI_RESULTS_URL + cfg.DESI_CHECKSUMS,
                    cfg.DESI_DIR / cfg.DESI_CHECKSUMS, manifest, refresh=refresh)
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
        if len(parts) < 4 or parts[1] != "base_omegak":
            continue
        if parts[0] == "iminuit" and parts[2] in cfg.DESI_OMK_BESTFITS \
                and parts[-1] == "bestfit.minimum.txt":
            wanted.append(rel)
        if parts[0] == "cobaya" and parts[2] in cfg.DESI_OMK_CHAINS:
            name = parts[-1]
            if name in ("chain.updated.yaml", "chain.margestats") or (
                    chains and re.fullmatch(r"chain\.\d+\.txt", name)):
                wanted.append(rel)
    return sorted(wanted)


# ---------------------------------------------------------------------------
# Fetch everything
# ---------------------------------------------------------------------------


def fetch_all(*, refresh: bool = False, chains: bool = True) -> dict:
    manifest = _load_manifest()
    try:
        print("Pantheon+ supernovae (D1's data)")
        for local, remote in cfg.PANTHEON_FILES.items():
            download(cfg.PANTHEON_RAW + remote, cfg.PANTHEON_DIR / local, manifest,
                     refresh=refresh)
        print("DESI DR2 BAO (D5's data)")
        for local, remote in cfg.BAO_FILES.items():
            download(cfg.BAO_RAW + remote, cfg.BAO_DIR / local, manifest, refresh=refresh)
        print("Planck PR4 early-universe chains (Lemos & Lewis 2023)")
        for i in (1, 2, 3, 4):
            name = f"{cfg.PLANCK_EARLY_ROOT}.{i}.txt"
            download(cfg.PLANCK_EARLY_RAW + name, cfg.PLANCK_EARLY_DIR / name, manifest,
                     refresh=refresh)
        print("Literature tables (cosmic chronometers; the most distant quasars)")
        _fetch_arxiv_member(cfg.MORESCO_LRR_EPRINT, "main.tex",
                            cfg.LITERATURE_DIR / "moresco2022_lrr.tex", manifest, refresh=refresh)
        download(cfg.YANG21_VIZIER, cfg.LITERATURE_DIR / "J_ApJ_923_262_table2.tsv", manifest,
                 refresh=refresh)
        print("DESI's curvature runs (answer key)" + ("" if chains else ", best fits only"))
        sums = desi_checksums(manifest, refresh=refresh)
        for rel in desi_wanted_files(sums, chains=chains):
            download(cfg.DESI_RESULTS_URL + rel, cfg.DESI_DIR / rel, manifest, refresh=refresh,
                     expected_sha256=sums[rel])
    finally:
        _save_manifest(manifest)
    print("Earlier atlas projects (results read in place)")
    for label, path in (("D1 supernova fit", cfg.D1_RESULTS), ("D3 CAMB fit", cfg.D3_FIT),
                        ("B6 oldest stars", cfg.B6_RESULTS), ("D5 CMB priors", cfg.D5_RESULTS)):
        state = "found" if path.exists() else "MISSING (README fallback values will be used)"
        print(f"  {label:18s} {state}")
    return manifest


# ---------------------------------------------------------------------------
# Loaders -- read only from the cache, never from the network
# ---------------------------------------------------------------------------


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


def _read_header_table(path: Path) -> dict[str, np.ndarray]:
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    names = lines[0].lstrip("#").split()
    rows = [ln.split() for ln in lines[1:] if not ln.startswith("#")]
    out = {}
    for i, name in enumerate(names):
        col = [r[i] for r in rows]
        try:
            out[name] = np.array(col, dtype=float)
        except ValueError:
            out[name] = np.array(col)
    return out


def load_pantheon() -> dict:
    """Pantheon+ as the published cosmology likelihood uses it: z_HD > 0.01,
    standardised magnitudes m_b_corr with the STAT+SYS covariance."""
    t = _read_header_table(_require(cfg.PANTHEON_DIR / "Pantheon+SH0ES.dat"))
    n = t["zHD"].size
    cache = cfg.DERIVED_DIR / "pantheon_cov.npy"
    if cache.exists():
        cov = np.load(cache)
    else:
        import pandas as pd

        flat = pd.read_csv(_require(cfg.PANTHEON_DIR / "Pantheon+SH0ES_STAT+SYS.cov"),
                           header=None).to_numpy(float).ravel()
        if int(flat[0]) != n or flat.size != 1 + n * n:
            raise ValueError("Pantheon+ covariance does not match the table")
        cov = flat[1:].reshape(n, n)
        np.save(cache, cov)
    keep = t["zHD"] > 0.01
    return {"z_hd": t["zHD"][keep], "z_hel": t["zHEL"][keep], "m": t["m_b_corr"][keep],
            "cov": cov[np.ix_(keep, keep)], "n_all": n}


def load_bao() -> dict:
    """DESI DR2: redshifts, values, observable kinds (D_V, D_M, D_H over r_d), covariance."""
    z, value, kind = [], [], []
    for line in _require(cfg.BAO_DIR / "desi_dr2_mean.txt").read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            a, b, c = line.split()
            z.append(float(a)), value.append(float(b)), kind.append(c)
    cov = np.loadtxt(_require(cfg.BAO_DIR / "desi_dr2_cov.txt"))
    if cov.shape != (len(z), len(z)):
        raise ValueError("BAO covariance does not match the data vector")
    return {"z": np.array(z), "value": np.array(value), "kind": np.array(kind), "cov": cov}


def load_planck_early_chains(burn: float = cfg.PLANCK_EARLY_BURN) -> dict[str, np.ndarray]:
    """Lemos & Lewis's four chains, burn-in removed, concatenated, with weights."""
    parts, names = [], None
    for i in (1, 2, 3, 4):
        path = _require(cfg.PLANCK_EARLY_DIR / f"{cfg.PLANCK_EARLY_ROOT}.{i}.txt")
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


def load_desi_chain(combo: str, *, burn: float = cfg.DESI_CHAIN_BURN) -> dict[str, np.ndarray]:
    """One of DESI's base_omegak chains (four MCMC runs, burn-in removed, weighted)."""
    folder = cfg.DESI_DIR / "cobaya" / "base_omegak" / combo
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


def load_desi_bestfit(combo: str) -> dict[str, float]:
    """DESI's MAP point for a base_omegak run (the best of its minimiser runs)."""
    best = None
    for path in sorted((cfg.DESI_DIR / "iminuit" / "base_omegak" / combo).glob(
            "**/bestfit.minimum.txt")):
        row = {k: float(v[0]) for k, v in _read_cobaya_table(path).items()}
        if best is None or row["minuslogpost"] < best["minuslogpost"]:
            best = row
    if best is None:
        raise FileNotFoundError(f"no DESI best fit cached for {combo}")
    return best


def _table_block(tex: str, label: str) -> str:
    idx = tex.index(f"\\label{{{label}}}")
    return tex[tex.rfind("\\begin{table", 0, idx):tex.find("\\end{table", idx)]


def load_chronometers() -> dict[str, np.ndarray]:
    """The 32 cosmic-chronometer H(z) points of Moresco et al. (2022), Table 1.

    Parsed from the review's LaTeX with the same rule as project D8. The error
    is the diagonal one tabulated (statistical + systematic); D8 shows that a
    proper fit needs the full covariance, so here the points are only drawn.
    """
    tex = _require(cfg.LITERATURE_DIR / "moresco2022_lrr.tex").read_text(errors="replace")
    num = r"([0-9]+(?:\.[0-9]*)?)"
    pat = re.compile(rf"^\s*{num}\s*&\s*{num}\s*&\s*{num}\s*&\s*([DFL])\s*&", re.M)
    rows = pat.findall(_table_block(tex, "tab:CC1"))
    if len(rows) != 32:
        raise ValueError(f"expected 32 chronometer rows, parsed {len(rows)}")
    return {"z": np.array([float(r[0]) for r in rows]), "H": np.array([float(r[1]) for r in rows]),
            "err": np.array([float(r[2]) for r in rows])}


def load_distant_quasars() -> dict[str, np.ndarray]:
    """[C II] redshifts of the z > 6.3 quasars of Yang et al. (2021), via I3's table."""
    path = _require(cfg.LITERATURE_DIR / "J_ApJ_923_262_table2.tsv")
    text = path.read_text(encoding="utf-8")
    if not text.lstrip().startswith("#"):
        raise ValueError(f"{path.name} is not a VizieR table (a bot-check page?)")
    lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    i = next(k for k, ln in enumerate(lines) if ln.strip())
    header = lines[i].split("\t")
    j_name, j_z = header.index("SName"), header.index("z[CII]")
    # The units row can be blank, so anchor on the row of dashes (as I3 does).
    j = next(k for k in range(i + 1, len(lines)) if re.fullmatch(r"[\t ]*-[-\t ]*", lines[k]))
    names, zs = [], []
    for ln in lines[j + 1:]:
        parts = ln.split("\t")
        try:
            zs.append(float(parts[j_z]))
            names.append(parts[j_name].strip())
        except (ValueError, IndexError):
            continue
    return {"name": np.array(names), "z": np.array(zs)}


# --- Results of sibling projects -------------------------------------------


def _read_json(path: Path) -> dict | None:
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def load_b6_age_floor() -> dict:
    """B6's lower limit on the age of the universe from the oldest globular clusters."""
    r = _read_json(cfg.B6_RESULTS)
    if r is not None and "cosmology" in r:
        t = r["cosmology"]["t_universe_min"]
        return {"t_min_gyr": t["median"], "t_min_err": t["std"], "t_min_p5": t["p5"],
                "source": "49-oldest-stars/data/results.json"}
    return {**cfg.B6_FALLBACK, "source": "B6 README (fallback)"}


def load_d1_h0() -> dict:
    """D1's H0 from the distance ladder with no cosmological model (SH0ES's 277
    Hubble-flow supernovae with a kinematic expansion): usable for any universe."""
    r = _read_json(cfg.D1_RESULTS)
    if r is not None and "cross_check_sh0es_hf" in r:
        c = r["cross_check_sh0es_hf"]
        return {"H0": c["H0"], "err": c["H0_err"], "source": "05-supernova-cosmology/data/results.json"}
    H0, err = cfg.D1_FALLBACK["H0_kinematic"]
    return {"H0": H0, "err": err, "source": "D1 README (fallback)"}


def load_d3_fit() -> dict | None:
    """D3's LCDM fit to the Planck TT spectrum (parameters, errors, correlations)."""
    return _read_json(cfg.D3_FIT)


def load_d5_cmb_priors() -> dict | None:
    """D5's compressed CMB priors (for cross-checking our own rebuild)."""
    r = _read_json(cfg.D5_RESULTS)
    return None if r is None else r.get("cmb_priors")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    parser.add_argument("--no-chains", action="store_true",
                        help="skip DESI's MCMC chains (the answer key)")
    args = parser.parse_args(argv)
    manifest = fetch_all(refresh=args.refresh, chains=not args.no_chains)
    total = sum(v["bytes"] for v in manifest.values())
    new = sum(v["bytes"] for v in manifest.values() if "linked_from" not in v)
    print(f"\n{len(manifest)} files, {total / 1e6:.0f} MB cached "
          f"({new / 1e6:.0f} MB downloaded, the rest linked from earlier projects)")
    sn = load_pantheon()
    print(f"  Pantheon+: {sn['m'].size} supernovae at z > 0.01")
    q = load_distant_quasars()
    print(f"  most distant quasar: {q['name'][np.argmax(q['z'])]} at z = {q['z'].max():.4f}")
    cc = load_chronometers()
    print(f"  cosmic chronometers: {cc['z'].size} points, z = {cc['z'].min()}-{cc['z'].max()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
