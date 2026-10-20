"""One fit over many days: alpha plus Sun-locked nuisance signals shared across days.

A violation of the redshift would be locked to the perigee (the clock runs slow
where the potential is deepest). Orbit errors from solar radiation pressure and
thermal effects on the clock are locked to the Sun instead: they depend on where
the satellite is relative to the Sun (orbit angle mu) and on the Sun's elevation
above the orbit plane (beta). Because the Sun's direction moves around the orbit
once a year while the perigee stays put, ten years of data can tell the two apart.

Model, for every day d (per-arc clock offset/rate profiled out, GLS noise model
of that day as in clockfit.py):

    y_d(t) = offsets + rates + alpha T_d(t)
             + sum_{k, j} [c_kj cos(k mu) + s_kj sin(k mu)] phi_j(day d)

where phi_j are slowly varying weights: indicator functions of time windows,
optionally times powers of beta. The per-day normal equations are reduced
(Schur complement over the local parameters) and summed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.linalg import cholesky_banded, solve_banded

import clockfit as F
import config as C
import prepare


@dataclass(frozen=True)
class Nuisance:
    name: str
    harmonics: tuple = ()          # Sun-locked harmonics of mu (1 = once per orbit)
    window_days: float | None = None   # length of the time windows the coefficients live in
    beta_degree: int = 0           # polynomial degree in beta within each window
    perigee_quadrature: bool = False   # also fit a cos-E-like perigee-locked term (diagnostic)


NUISANCE_MODELS = {
    "none": Nuisance("none"),
    "sun1_1yr": Nuisance("sun1_1yr", (1,), 365.25, 0),
    "sun1_1yr_beta2": Nuisance("sun1_1yr_beta2", (1,), 365.25, 2),
    "sun12_1yr_beta2": Nuisance("sun12_1yr_beta2", (1, 2), 365.25, 2),
    "sun1_6mo_beta1": Nuisance("sun1_6mo_beta1", (1,), 182.6, 1),
    "sun12_6mo_beta2": Nuisance("sun12_6mo_beta2", (1, 2), 182.6, 2),
}


def _whiten_factory(q, sw, i0, i1):
    s2 = max(sw, F.SW_FLOOR) ** 2
    diag = q * (F.T_SEC[i1] - F.T_SEC[i0]) + 2 * s2
    shares = (i0[1:] == i1[:-1])
    off = np.where(shares, -s2, 0.0)
    ab = np.zeros((2, len(diag)))
    ab[0] = diag
    ab[1, :-1] = off
    L = cholesky_banded(ab, lower=True)
    return lambda M: solve_banded((1, 0), L, M)


def day_normals(y, labels, use, q, sw, local_cols, global_cols):
    """Reduced normal equations of one day for the global columns (differenced GLS)."""
    idx = np.flatnonzero(use)
    lab = labels[idx]
    same = lab[1:] == lab[:-1]
    i0, i1 = idx[:-1][same], idx[1:][same]
    whiten = _whiten_factory(q, sw, i0, i1)
    D = whiten((y[i1] - y[i0])[:, None])[:, 0]
    XL = whiten(local_cols[i1] - local_cols[i0])
    XG = whiten(global_cols[i1] - global_cols[i0])
    NLL = XL.T @ XL
    NLG = XL.T @ XG
    NGG = XG.T @ XG
    bL = XL.T @ D
    bG = XG.T @ D
    NLL_inv = np.linalg.pinv(NLL)
    return NGG - NLG.T @ NLL_inv @ NLG, bG - NLG.T @ NLL_inv @ bL, D, XL, XG


def global_columns(d: dict, i: int, alpha_tpl: np.ndarray, model: Nuisance, windows: np.ndarray,
                   n_windows: int) -> tuple[np.ndarray, list]:
    cols, names = [alpha_tpl], ["alpha"]
    mu = d["mu"][i].astype(float)
    beta = float(np.nanmean(d["beta"][i]))
    w = windows[i]
    for k in model.harmonics:
        for trig, fn in (("cos", np.cos), ("sin", np.sin)):
            base = fn(k * mu)
            for j in range(n_windows):
                for p in range(model.beta_degree + 1):
                    names.append(f"{trig}{k}_w{j}_b{p}")
                    cols.append(base * (np.radians(beta) ** p) if j == w else np.zeros_like(base))
    if model.perigee_quadrature:
        E = d["E"][i].astype(float)
        e = float(np.nanmean(d["e"][i]))
        a = float(np.nanmean(d["a"][i]))
        amp = np.sqrt(C.GM_EARTH * a) * e / C.C_LIGHT**2
        cols.append(-amp * np.cos(E))          # same size as T_alpha, 90 degrees out of phase
        names.append("perigee_quad")
    return np.vstack(cols).T, names


def _solve(N: np.ndarray, b: np.ndarray):
    live = np.diag(N) > 0
    cov = np.zeros_like(N)
    cov[np.ix_(live, live)] = np.linalg.pinv(N[np.ix_(live, live)])
    return cov @ b, cov, live


def perigee_angle(d: dict, i: int) -> float:
    """Argument of perigee (rad) of day i: u - true anomaly, circular mean over the day."""
    E = d["E"][i].astype(float)
    e = float(np.nanmean(d["e"][i]))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    return float(np.angle(np.nanmean(np.exp(1j * (d["u"][i].astype(float) - nu)))))


def fake_template(d: dict, i: int, omega_fn, amplitude: float) -> np.ndarray:
    """A perigee-locked-looking template with a made-up perigee angle omega_fn(day index).

    Same shape and size as the real one (-A sin E ~ -A sin(u - omega)), but with a
    perigee that does not exist: it carries no gravity signal, yet it is exposed to
    the same Sun-locked systematics. Used for null tests.
    """
    u = d["u"][i].astype(float)
    return -amplitude * np.sin(u - omega_fn(i))


def fit(prn: str, days: pd.DataFrame, model: Nuisance, *, template: str = "alpha",
        poly: int = 1, d: dict | None = None, n_boot: int = 500,
        block_days: int = C.BOOTSTRAP_BLOCK_DAYS, seed: int = C.RANDOM_SEED,
        template_fn=None) -> dict:
    """Global fit over the given days (rows of a clockfit.run_fits table with ok=True).

    Returns alpha with its noise-model error and a block-bootstrap error
    (30-day blocks of days resampled with replacement).
    """
    d = d if d is not None else prepare.load(prn)
    date_index = {s: k for k, s in enumerate(d["date"])}
    rows = np.array([date_index[s] for s in days["date"]])
    t_days = pd.to_datetime(days["date"]).to_numpy().astype("datetime64[D]").astype(np.int64)
    t_rel = t_days - t_days.min()
    win = (t_rel // model.window_days).astype(int) if model.window_days else np.zeros(len(rows), int)
    n_windows = int(win.max()) + 1
    windows = np.full(len(d["date"]), -1)
    windows[rows] = win
    block = (t_rel // block_days).astype(int)
    blocks = {}
    names = None
    n_used = 0
    for i, blk in zip(rows, block):
        y, tpl = F.series_for(d, i, "full")
        if template_fn is not None:
            tpl = {template: template_fn(d, i)}
        f = F.fit_day(y, {template: tpl[template]}, poly=poly, keep_resid=True)
        if not f.ok:
            continue
        labels = F.find_arcs(y, np.isfinite(y) & np.isfinite(tpl[template]))
        use = f.mask
        local = []
        for a in np.unique(labels[use]):
            m = (labels == a).astype(float)
            for k in range(1, poly + 1):
                local.append(m * F.TC**k)
        local = np.vstack(local).T
        G, names = global_columns(d, i, tpl[template], model, windows, n_windows)
        N, b, *_ = day_normals(y, labels, use, f.q, f.sw, local, G)
        if blk not in blocks:
            blocks[blk] = [np.zeros_like(N), np.zeros_like(b)]
        blocks[blk][0] += N
        blocks[blk][1] += b
        n_used += 1
    keys = sorted(blocks)
    Nb = np.stack([blocks[k][0] for k in keys])
    bb = np.stack([blocks[k][1] for k in keys])
    coef, cov, live = _solve(Nb.sum(0), bb.sum(0))
    k = names.index("alpha")
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    for r in range(n_boot):
        pick = rng.integers(0, len(keys), len(keys))
        cb, _, _ = _solve(Nb[pick].sum(0), bb[pick].sum(0))
        boot[r] = cb[k]
    sd = np.sqrt(np.clip(np.diag(cov), 1e-300, None))
    corr = cov[k] / (sd * sd[k])
    out = {"prn": prn, "model": model.name, "template": template, "n_days": n_used,
           "alpha": float(coef[k]), "err_model": float(np.sqrt(cov[k, k])),
           "err_boot": float(np.std(boot)), "n_params": int(live.sum()),
           "max_corr_with_nuisance": float(np.max(np.abs(np.delete(corr, k)))) if live.sum() > 1 else 0.0,
           "coef": dict(zip(names, coef.tolist())), "names": names}
    if "perigee_quad" in names:
        kq = names.index("perigee_quad")
        out.update(perigee_quad=float(coef[kq]), perigee_quad_err=float(np.sqrt(cov[kq, kq])))
    return out
