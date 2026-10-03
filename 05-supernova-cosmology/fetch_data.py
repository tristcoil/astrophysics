"""Download and cache everything the analysis needs, and load it back.

All data come from the public Pantheon+ / SH0ES data release on GitHub
(https://github.com/PantheonPlusSH0ES/DataRelease):

1. **The distance table** -- 1701 Type Ia supernova light curves with their
   redshifts and standardised peak magnitudes. For the 77 light curves of
   supernovae that went off in galaxies whose distances the SH0ES team measured
   with Cepheid variable stars, it also gives that Cepheid distance.

2. **Two covariance matrices** of those magnitudes, 1701 x 1701: statistical
   errors only, and statistical plus systematic. Supernovae share calibration,
   dust and selection systematics, so their errors are correlated, and the data
   release is explicit that cosmology must be fitted with the full matrix.

3. **The Pantheon+ team's own MCMC chains**, used purely as an answer key: our
   likelihood is written independently, and these check it.

Everything is cached under ``data/`` with a manifest recording where each file
came from, when, and its SHA-256. Nothing already present is downloaded again, so
once the cache is warm the whole project runs offline.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from config import CHAIN_DIR, CHAIN_FILES, DATA_DIR, DATA_FILES, MANIFEST_PATH, REPO_RAW

#: Number of light curves in the release; the covariance files must match it.
N_LIGHT_CURVES = 1701


def _url(repo_path: str) -> str:
    # '+' is legal in a URL path but some servers read it as a space; encode it.
    return REPO_RAW + urllib.parse.quote(repo_path)


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


def _download(repo_path: str, dest: Path, *, refresh: bool, manifest: dict) -> None:
    """Fetch one file unless it is already cached, and record its provenance."""
    key = str(dest.relative_to(DATA_DIR))
    if dest.exists() and not refresh:
        print(f"  cached   {key}")
        return

    url = _url(repo_path)
    print(f"  fetching {key}  <-  {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as out:
        shutil.copyfileobj(response, out)
    tmp.replace(dest)  # atomic: a half-written file is never mistaken for a cached one

    manifest[key] = {
        "url": url,
        "bytes": dest.stat().st_size,
        "sha256": _sha256(dest),
        "retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    print(f"           {manifest[key]['bytes'] / 1e6:.1f} MB")


def fetch_all(*, refresh: bool = False, chains: bool = True) -> dict:
    """Populate ``data/``. Returns the manifest."""
    manifest = _load_manifest()
    print("Pantheon+ distances and covariances")
    for local, repo_path in DATA_FILES.items():
        _download(repo_path, DATA_DIR / local, refresh=refresh, manifest=manifest)
    if chains:
        print("Pantheon+ published MCMC chains (answer key only)")
        for name, repo_path in CHAIN_FILES.items():
            _download(repo_path, CHAIN_DIR / f"{name}.txt.gz", refresh=refresh,
                      manifest=manifest)
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


# ---------------------------------------------------------------------------
# Loaders -- read only from the cache, never from the network
# ---------------------------------------------------------------------------


def _require(path: Path) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


def load_table():
    """The distance table as a pandas DataFrame, one row per light curve."""
    import pandas as pd

    table = pd.read_csv(_require(DATA_DIR / "Pantheon+SH0ES.dat"), sep=r"\s+")
    if len(table) != N_LIGHT_CURVES:
        raise ValueError(f"expected {N_LIGHT_CURVES} light curves, found {len(table)}")
    return table


def load_covariance(kind: str = "statsys") -> np.ndarray:
    """The 1701 x 1701 magnitude covariance, in mag^2.

    ``kind`` is ``"statsys"`` (statistical + systematic; the one to fit with) or
    ``"stat"`` (statistical only). The text file -- a count N followed by N*N
    values, one per line, row by row -- is parsed once and cached as ``.npy``.
    """
    name = {"statsys": "Pantheon+SH0ES_STAT+SYS.cov",
            "stat": "Pantheon+SH0ES_STATONLY.cov"}[kind]
    text_path = _require(DATA_DIR / name)
    npy_path = text_path.with_suffix(".npy")
    if npy_path.exists() and npy_path.stat().st_mtime >= text_path.stat().st_mtime:
        return np.load(npy_path)

    import pandas as pd

    values = pd.read_csv(text_path, header=None).iloc[:, 0].to_numpy(dtype=float)
    n = int(values[0])
    if n != N_LIGHT_CURVES or values.size != 1 + n * n:
        raise ValueError(f"{name}: header says N={n}, file holds {values.size - 1} values")
    cov = values[1:].reshape(n, n)
    asym = np.max(np.abs(cov - cov.T))
    if asym > 1e-6 * np.max(np.abs(cov)):
        raise ValueError(f"{name} is not symmetric (max |C - C^T| = {asym:.3g})")
    cov = 0.5 * (cov + cov.T)  # remove round-off asymmetry from the text format
    np.save(npy_path, cov)
    return cov


def load_chain(name: str, *, burn_fraction: float = 0.3) -> dict[str, np.ndarray]:
    """One of the published chains, as {column: samples}, burn-in removed.

    The chains are emcee runs (16 walkers). The first ``burn_fraction`` of the
    samples is discarded while the walkers are still converging. Columns are
    renamed from CosmoSIS's ``section--name`` form to short names, and derived
    columns are added: ``OL`` = 1 - Om - Ok for the non-flat chain, and ``H0`` in
    km/s/Mpc.
    """
    path = _require(CHAIN_DIR / f"{name}.txt.gz")
    with gzip.open(path, "rt") as f:
        header = f.readline().lstrip("#").strip().split("\t")
    samples = np.loadtxt(path)
    samples = samples[int(burn_fraction * len(samples)):]

    short = {"cosmological_parameters--omega_m": "Om",
             "cosmological_parameters--h0": "h",
             "cosmological_parameters--omega_k": "Ok",
             "cosmological_parameters--w": "w",
             "supernova_params--m": "M",
             "prior": "log_prior", "post": "log_post"}
    chain = {short.get(col, col): samples[:, i] for i, col in enumerate(header)}
    chain["H0"] = 100.0 * chain["h"]
    if "Ok" in chain:
        chain["OL"] = 1.0 - chain["Om"] - chain["Ok"]
    return chain


def chain_prior_ranges(name: str) -> dict[str, tuple[float, float]]:
    """Parameter ranges (flat priors) recorded in a chain's header."""
    path = _require(CHAIN_DIR / f"{name}.txt.gz")
    ranges = {}
    in_values = False
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                break
            body = line.lstrip("#").strip()
            if body == "START_OF_VALUES_INI":
                in_values = True
            elif body == "END_OF_VALUES_INI":
                break
            elif in_values and "=" in body:
                key, _, rhs = body.partition("=")
                parts = rhs.split()
                if len(parts) == 3:  # "low start high" marks a varied parameter
                    ranges[key.strip()] = (float(parts[0]), float(parts[2]))
    return ranges


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--refresh", action="store_true", help="re-download everything")
    parser.add_argument("--no-chains", action="store_true",
                        help="skip the published chains (validation only)")
    args = parser.parse_args(argv)

    manifest = fetch_all(refresh=args.refresh, chains=not args.no_chains)
    total = sum(v["bytes"] for v in manifest.values())
    print(f"\n{len(manifest)} files, {total / 1e6:.0f} MB in {DATA_DIR}")

    # Parse the covariances now so the first analysis run is quick.
    for kind in ("statsys", "stat"):
        cov = load_covariance(kind)
        print(f"covariance {kind:<8} {cov.shape[0]} x {cov.shape[1]}, "
              f"median diagonal error {np.sqrt(np.median(np.diag(cov))):.3f} mag")
    return 0


if __name__ == "__main__":
    sys.exit(main())
