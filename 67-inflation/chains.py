"""Load the published MCMC chains (Planck 2018 CosmoMC, ACT DR6 and SPT-3G Cobaya)
into one format, and turn weighted samples into constraints and contours.

Every chain becomes a dict of numpy arrays with standard names -- ``weight``,
``ns``, ``r`` (tensor-to-scalar ratio at 0.05 Mpc^-1, where present), ``tau``,
``H0``, ``ombh2``, ``omch2``, ``logA`` -- cached as ``data/derived/chains/*.npz``
so the text files are parsed once.
"""

from __future__ import annotations

import io
import json
import math
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.special import erfinv

import config as cfg

CACHE_DIR = cfg.DERIVED_DIR / "chains"
#: Fraction of each Cobaya chain dropped as burn-in (the runs were made with
#: burn_in: 0; Planck's CosmoMC chains have burn-in removed already).
COBAYA_BURN = 0.3

#: Column names in each code's output -> our standard names.
RENAME = {
    "n_s": "ns", "ns": "ns", "r": "r", "tau": "tau", "tau_reio": "tau", "H0": "H0",
    "H0*": "H0", "omegabh2": "ombh2", "ombh2": "ombh2", "omega_b": "ombh2",
    "omegach2": "omch2", "omch2": "omch2", "omega_cdm": "omch2", "logA": "logA",
    "r02*": "r002", "chi2__bicep_keck_2018": "chi2_bk18", "chi2__BKPLANCK": "chi2_bk15",
    "BBdust": "BBdust", "BBsync": "BBsync", "BBalphadust": "BBalphadust",
    "BBalphasync": "BBalphasync", "BBbetadust": "BBbetadust", "BBbetasync": "BBbetasync",
    "BBdustsynccorr": "BBdustsynccorr", "omegamh2*": "ommh2", "omegamh2": "ommh2",
    "cosmomc_theta": "theta", "theta": "theta100",
}


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------


def _read_table(handle, names: list[str] | None = None) -> pd.DataFrame:
    return pd.read_csv(handle, sep=r"\s+", header=None, names=names, comment="#",
                       engine="c", dtype=np.float64)


def _standardise(df: pd.DataFrame, weight_col: str = "weight") -> dict[str, np.ndarray]:
    out = {"weight": df[weight_col].to_numpy(float)}
    for col in df.columns:
        std = RENAME.get(col)
        if std and std not in out:
            out[std] = df[col].to_numpy(float)
    if "minuslogpost" in df.columns:
        out["minuslogpost"] = df["minuslogpost"].to_numpy(float)
    return out


def _concat(parts: list[dict]) -> dict[str, np.ndarray]:
    keys = set(parts[0]).intersection(*parts[1:]) if len(parts) > 1 else set(parts[0])
    return {k: np.concatenate([p[k] for p in parts]) for k in sorted(keys)}


def read_cosmomc(folder: Path) -> dict[str, np.ndarray]:
    """A Planck-grid run: chain files root_1.txt ... with names in root.paramnames."""
    root = folder.name
    names = ["weight", "minuslogpost"] + [
        ln.split()[0] for ln in (folder / f"{root}.paramnames").read_text().splitlines()
        if ln.strip()]
    burn_removed = "burn_removed=T" in (folder / f"{root}.properties.ini").read_text()
    if not burn_removed:
        raise ValueError(f"{root}: burn-in not removed; handle it explicitly")
    parts = []
    for path in sorted(folder.glob(f"{root}_[0-9]*.txt")):
        parts.append(_standardise(_read_table(path, names)))
    return _concat(parts)


def _cobaya_part(handle, burn: float) -> dict[str, np.ndarray]:
    text = handle.read().decode() if hasattr(handle, "read") else Path(handle).read_text()
    header = text.split("\n", 1)[0].lstrip("#").split()
    df = _read_table(io.StringIO(text), header)
    df = df.iloc[int(burn * len(df)):]
    return _standardise(df)


def read_cobaya_dir(folder: Path, burn: float = COBAYA_BURN) -> dict[str, np.ndarray]:
    """An SPT-3G run: CLASS.1.txt ... (Cobaya format, burn-in not removed)."""
    parts = [_cobaya_part(p, burn) for p in sorted(folder.glob("*.[0-9].txt"))]
    if not parts:
        raise FileNotFoundError(folder)
    return _concat(parts)


def read_cobaya_tar(path: Path, burn: float = COBAYA_BURN) -> dict[str, np.ndarray]:
    """An ACT DR6 run: name.tar.gz holding name/name.1.txt ... (Cobaya format)."""
    parts = []
    with tarfile.open(path, "r:gz") as tar:
        for member in sorted(tar.getmembers(), key=lambda m: m.name):
            stem = member.name.rsplit("/", 1)[-1]
            bits = stem.split(".")
            if len(bits) >= 3 and bits[-1] == "txt" and bits[-2].isdigit():
                parts.append(_cobaya_part(tar.extractfile(member), burn))
    if not parts:
        raise FileNotFoundError(f"no chain files in {path}")
    return _concat(parts)


def cobaya_config(path: Path) -> str:
    """The updated.yaml of an ACT run (provenance: which likelihoods were used)."""
    with tarfile.open(path, "r:gz") as tar:
        for member in tar.getmembers():
            if member.name.endswith(".updated.yaml") and "minimize" not in member.name:
                return tar.extractfile(member).read().decode()
    return ""


# ---------------------------------------------------------------------------
# The catalogue of chains used in this project
# ---------------------------------------------------------------------------


def catalogue() -> dict[str, dict]:
    """key -> {label, source, loader args}."""
    cat = {}
    for key, (folder, label) in cfg.PLANCK_CHAINS.items():
        cat[key] = {"label": label, "kind": "cosmomc",
                    "path": cfg.PLANCK_DIR / folder.replace("/", "_"),
                    "source": f"Planck 2018 parameter grid: {folder}"}
    for name, label in {**cfg.ACT_R_CHAINS, **cfg.ACT_LCDM_CHAINS}.items():
        cat["act_" + name] = {"label": label, "kind": "cobaya_tar",
                              "path": cfg.ACT_DIR / f"{name}.tar.gz",
                              "source": f"ACT DR6.02 chains on LAMBDA: {name}"}
    for key, (_archive, folder, label) in cfg.SPT_CHAINS.items():
        cat[key] = {"label": label, "kind": "cobaya_dir",
                    "path": cfg.SPT_DIR / folder.split("/")[-1],
                    "source": f"SPT-3G D1 chains on LAMBDA: {folder}"}
    return cat


def load(key: str, *, force: bool = False) -> dict[str, np.ndarray]:
    """Load one chain by catalogue key (parsed once, then from the npz cache)."""
    entry = catalogue()[key]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache = CACHE_DIR / f"{key}.npz"
    if cache.exists() and not force:
        with np.load(cache) as z:
            return {k: z[k] for k in z.files}
    if entry["kind"] == "cosmomc":
        data = read_cosmomc(entry["path"])
    elif entry["kind"] == "cobaya_tar":
        data = read_cobaya_tar(entry["path"])
    else:
        data = read_cobaya_dir(entry["path"])
    keep = data["weight"] > 0
    data = {k: v[keep] for k, v in data.items()}
    np.savez_compressed(cache, **data)
    return data


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------


def n_eff(w: np.ndarray) -> float:
    """Kish effective sample size of weighted samples."""
    return float(w.sum() ** 2 / np.sum(w ** 2))


def mean_std(x: np.ndarray, w: np.ndarray) -> tuple[float, float]:
    m = np.average(x, weights=w)
    return float(m), float(math.sqrt(np.average((x - m) ** 2, weights=w)))


def quantile(x: np.ndarray, w: np.ndarray, q) -> np.ndarray:
    """Weighted quantiles (linear interpolation of the cumulative weight)."""
    order = np.argsort(x)
    xs, ws = x[order], w[order]
    c = (np.cumsum(ws) - 0.5 * ws) / ws.sum()
    return np.interp(q, c, xs)


def upper_limit(x: np.ndarray, w: np.ndarray, level: float = 0.95, *,
                lower_bound: float = 0.0, bins: int = 400) -> float:
    """One-tail upper limit for a parameter with a hard lower prior bound (r >= 0).

    A smoothed marginal density with a reflection at the bound (as GetDist's
    boundary correction does), integrated to ``level``.
    """
    hi = quantile(x, w, 0.9995) * 1.3
    edges = np.linspace(lower_bound, hi, bins + 1)
    h, _ = np.histogram(x, edges, weights=w)
    width = edges[1] - edges[0]
    sig = max(1.0, 0.25 * mean_std(x, w)[1] / width * n_eff(w) ** -0.2 * 4)
    hh = np.concatenate([h[::-1], h])
    sm = gaussian_filter(hh, sig, mode="constant")[bins:]
    cdf = np.cumsum(sm) / sm.sum()
    return float(np.interp(level, cdf, edges[1:]))


def summary(x: np.ndarray, w: np.ndarray) -> dict:
    m, s = mean_std(x, w)
    q = quantile(x, w, [0.02275, 0.15865, 0.5, 0.84135, 0.97725])
    return {"mean": m, "std": s, "median": float(q[2]), "lo68": float(q[1]),
            "hi68": float(q[3]), "lo95": float(q[0]), "hi95": float(q[4])}


# ---------------------------------------------------------------------------
# Two-dimensional densities and credible regions
# ---------------------------------------------------------------------------


def density_2d(x: np.ndarray, y: np.ndarray, w: np.ndarray, xr: tuple, yr: tuple, *,
               n: int = 300, y_lower_bound: float | None = 0.0,
               smooth: float = 1.0) -> dict:
    """Smoothed 2D posterior density on a grid, with a reflection boundary at
    y = y_lower_bound (r >= 0), and the density levels enclosing 68% and 95%.

    The Gaussian kernel follows Scott's rule, scaled per axis by each
    parameter's standard deviation; ``smooth`` multiplies the kernel width.
    """
    xe = np.linspace(*xr, n + 1)
    ye = np.linspace(*yr, n + 1)
    dx, dy = xe[1] - xe[0], ye[1] - ye[0]
    ne = n_eff(w)
    bw = smooth * ne ** (-1.0 / 6.0)
    sx = max(0.6, bw * mean_std(x, w)[1] / dx)
    sy = max(0.6, bw * mean_std(y, w)[1] / dy)
    if y_lower_bound is not None:
        # mirror the samples about the bound so the kernel does not leak density
        # into the excluded region (equivalent to a reflection boundary)
        x2 = np.concatenate([x, x])
        y2 = np.concatenate([y, 2 * y_lower_bound - y])
        w2 = np.concatenate([w, w])
    else:
        x2, y2, w2 = x, y, w
    h, _, _ = np.histogram2d(x2, y2, [xe, ye], weights=w2)
    p = gaussian_filter(h, [sx, sy], mode="constant")
    if y_lower_bound is not None:
        p[:, 0.5 * (ye[:-1] + ye[1:]) < y_lower_bound] = 0.0
    p /= p.sum()
    levels = hpd_levels(p, (0.6827, 0.9545))
    xc, yc = 0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:])
    inside = float(h.sum() / w2.sum())
    return {"x": xc, "y": yc, "p": p, "levels": levels, "kernel_bins": (sx, sy),
            "mass_in_grid": inside}


def hpd_levels(p: np.ndarray, masses) -> list[float]:
    """Density thresholds whose super-level sets enclose the given masses."""
    flat = np.sort(p.ravel())[::-1]
    c = np.cumsum(flat)
    c /= c[-1]
    return [float(flat[np.searchsorted(c, m)]) for m in masses]


def credibility(d: dict, x0: float, y0: float) -> float:
    """Posterior mass inside the highest-density contour through (x0, y0).

    0.68 means the point sits on the 68% contour; values near 1 mean it lies far
    outside. Points beyond the grid return 1.
    """
    x, y, p = d["x"], d["y"], d["p"]
    if not (x[0] <= x0 <= x[-1] and y[0] <= y0 <= y[-1]):
        return 1.0
    from scipy.interpolate import RegularGridInterpolator
    val = float(RegularGridInterpolator((x, y), p)((x0, y0)))
    return float(p[p > val].sum())


#: Beyond this credibility the density at a point rests on a handful of samples
#: smoothed by a narrow kernel (it drops to zero a few sigma out): report
#: "beyond 99.9%" (3.3 sigma) rather than a number.
CRED_RESOLVED = 0.999


def sigma_equivalent(cred: float) -> float:
    """Two-sided Gaussian sigma with the same probability (0.6827 -> 1, 0.9545 -> 2),
    the convention of the 68% / 95% contours. Capped where samples run out."""
    cred = min(cred, CRED_RESOLVED)
    return float(math.sqrt(2) * erfinv(cred))


def save_json(path: Path, obj) -> None:
    def conv(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    path.write_text(json.dumps(obj, indent=1, default=conv))
