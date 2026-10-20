"""Fit the gravitational-redshift test to the eccentric Galileo clocks, one day at a time.

For each satellite-day we form

    y = clk + (-2 r.v/c^2) - tau_GR

the published clock, with the IGS convention term added back (giving the
physical clock offset), minus everything general relativity predicts
(Earth monopole + flattening + tides + speed). If GR is right, y is clock
noise plus whatever orbit errors leaked into the clock estimates.

Following Delva et al. (2018), a violation of the redshift is parametrised as

    d(tau)/dt = 1 - (1 + alpha) U/c^2 - v^2/(2c^2),

so it adds alpha times the *gravitational* part of the proper time,
T_alpha = -int (U - <U>)/c^2 dt (190 ns amplitude for E14/E18), to y. Each day
gets its own clock offset and rate (Delva's "daily linear fit"; one pair per arc
if the clock jumps). Herrmann et al. (2018) instead scaled the whole periodic
term (gravity + speed); that parameter is very nearly alpha/2 and is also reported.

Noise: the clock estimates behave like a random walk (white frequency noise of
the maser plus estimation noise) on top of white phase noise. Both levels are
estimated per day from the residuals, and give each day's alpha its exact
variance p^T C p under that model (p: the alpha row of the least-squares
solution). Days are combined with inverse-variance weights; a block bootstrap
over 30-day blocks checks the error bar against slow, correlated effects.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

import config as C
import prepare

DT = C.EPOCH_STEP_S
N = prepare.N_EPOCHS
T_SEC = np.arange(N) * DT
TC = (T_SEC - T_SEC.mean()) / T_SEC.mean()        # day time scaled to [-1, 1]


@dataclass
class DayFit:
    ok: bool
    reason: str = ""
    coef: dict | None = None        # template name -> coefficient
    sigma: dict | None = None       # template name -> 1-sigma (noise model)
    gls: dict | None = None         # GLS estimates under the same noise model
    gls_sigma: dict | None = None
    n: int = 0
    arcs: int = 0
    rms: float = np.nan             # s, residual rms
    q: float = np.nan               # s^2/s, random-walk diffusion
    sw: float = np.nan              # s, white phase noise
    resid: np.ndarray | None = None
    mask: np.ndarray | None = None


def robust_sigma(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) else np.nan


def find_arcs(y: np.ndarray, valid: np.ndarray) -> np.ndarray:
    """Label continuous arcs: split where the clock steps by more than C.JUMP_PS between epochs.

    A lone outlier gives two big opposite differences and is not a step; it is
    left for the outlier rejection. Returns an arc label per epoch (-1 = unusable).
    """
    labels = np.full(N, -1)
    idx = np.flatnonzero(valid)
    if len(idx) < 3:
        return labels
    yy = y[idx]
    # remove a straight line so the clock rate does not count as steps
    p = np.polyfit(idx, yy, 1)
    d = np.diff(yy - np.polyval(p, idx))
    s = robust_sigma(d)
    thr = max(C.JUMP_PS * 1e-12, 8 * s)
    big = np.abs(d) > thr
    steps = []
    for k in np.flatnonzero(big):
        lone_outlier = (k + 1 < len(d) and big[k + 1] and np.sign(d[k]) != np.sign(d[k + 1])
                        and abs(d[k] + d[k + 1]) < thr)
        prev_outlier = (k > 0 and big[k - 1] and np.sign(d[k]) != np.sign(d[k - 1])
                        and abs(d[k] + d[k - 1]) < thr)
        if not (lone_outlier or prev_outlier):
            steps.append(k)
    arc, start = 0, 0
    for k in steps + [len(idx) - 1]:
        seg = idx[start:k + 1]
        if len(seg) >= 48:                       # at least 4 h; shorter pieces are dropped
            labels[seg] = arc
            arc += 1
        start = k + 1
    return labels


def design(templates: dict, labels: np.ndarray, use: np.ndarray, poly: int = 1):
    """Columns: per-arc offset (and rate, ...), then the templates."""
    cols, names = [], []
    for a in np.unique(labels[use]):
        m = (labels == a).astype(float)
        for k in range(poly + 1):
            cols.append(m * TC**k)
            names.append(f"arc{a}_p{k}")
    for name, tpl in templates.items():
        cols.append(tpl)
        names.append(name)
    return np.vstack(cols).T[use], names


def noise_levels(res: np.ndarray, use: np.ndarray) -> tuple[float, float]:
    """Random-walk diffusion q (s^2/s) and white phase noise sw (s) from first differences.

    For x = random walk + white noise: var(dx) = q*dt + 2 sw^2 and lag-1
    autocovariance of dx = -sw^2. Robust (MAD-based) variance; gaps are skipped.
    """
    idx = np.flatnonzero(use)
    ok = np.diff(idx) == 1
    d = np.diff(res[idx])[ok]
    if len(d) < 50:
        return np.nan, np.nan
    var = robust_sigma(d) ** 2
    pairs = ok[:-1] & ok[1:]
    d0 = np.diff(res[idx])[:-1][pairs]
    d1 = np.diff(res[idx])[1:][pairs]
    # robust lag-1 autocovariance via the variance of sums and differences
    acov = 0.25 * (robust_sigma(d0 + d1) ** 2 - robust_sigma(d0 - d1) ** 2)
    sw2 = max(-acov, 0.0)
    q = max(var - 2 * sw2, 0.05 * var) / DT
    return q, np.sqrt(sw2)


SW_FLOOR = 1e-12      # s; keeps the noise model positive definite when no white noise is seen


def ols_variance(p: np.ndarray, t: np.ndarray, q: float, sw: float) -> float:
    """p^T C p for C = q * min(s_i, s_j) + sw^2 I (random walk + white phase noise), in O(n).

    min(s_i, s_j) = sum_k ds_k [k <= i][k <= j], so p^T K p = sum_k ds_k (sum_{i >= k} p_i)^2.
    The random walk's starting value does not matter: p sums to zero over every arc.
    """
    ds = np.diff(t, prepend=t[0] - DT)
    tail = np.cumsum(p[::-1])[::-1]
    return float(q * np.sum(ds * tail**2) + max(sw, SW_FLOOR) ** 2 * np.sum(p**2))


def gls_differenced(y: np.ndarray, X: np.ndarray, labels: np.ndarray, use: np.ndarray,
                    q: float, sw: float, n_offsets: int):
    """GLS with free per-arc offsets, done exactly on within-arc first differences.

    Differencing turns the random walk into white noise and the white phase noise
    into a tridiagonal (MA(1)) covariance, so the solve is O(n).
    Returns (coefficients, covariance) of the non-offset columns of X.
    """
    from scipy.linalg import solveh_banded
    idx = np.flatnonzero(use)
    lab = labels[idx]
    same = lab[1:] == lab[:-1]
    i0, i1 = idx[:-1][same], idx[1:][same]
    keep_cols = [c for c in range(X.shape[1]) if c >= n_offsets]
    Xu = np.full((N, X.shape[1]), np.nan)
    Xu[use] = X
    D = y[i1] - y[i0]
    XD = Xu[i1][:, keep_cols] - Xu[i0][:, keep_cols]
    s2 = max(sw, SW_FLOOR) ** 2
    diag = q * (T_SEC[i1] - T_SEC[i0]) + 2 * s2
    # consecutive differences share an epoch when the second starts where the first ended
    shares = (i0[1:] == i1[:-1])
    off = np.where(shares, -s2, 0.0)
    ab = np.zeros((2, len(D)))
    ab[0] = diag
    ab[1, :-1] = off
    rhs = np.column_stack([XD, D])
    sol = solveh_banded(ab, rhs, lower=True)
    A = XD.T @ sol[:, :-1]
    b = XD.T @ sol[:, -1]
    cov = np.linalg.inv(A)
    return cov @ b, cov


def fit_day(y: np.ndarray, templates: dict, *, poly: int = 1, keep_resid: bool = False) -> DayFit:
    valid = np.isfinite(y)
    for tpl in templates.values():
        valid &= np.isfinite(tpl)
    if valid.sum() < C.MIN_EPOCHS_PER_DAY:
        return DayFit(False, "too few epochs")
    labels = find_arcs(y, valid)
    use = labels >= 0
    if use.sum() < C.MIN_EPOCHS_PER_DAY:
        return DayFit(False, "too few epochs after arc split")
    # iterative robust fit (outliers beyond C.OUTLIER_SIGMA x MAD)
    for _ in range(5):
        X, names = design(templates, labels, use, poly)
        coef, *_ = np.linalg.lstsq(X, y[use], rcond=None)
        res = np.full(N, np.nan)
        res[use] = y[use] - X @ coef
        s = robust_sigma(res[use])
        bad = use & (np.abs(res) > C.OUTLIER_SIGMA * s)
        if not bad.any():
            break
        use &= ~bad
        if use.sum() < C.MIN_EPOCHS_PER_DAY:
            return DayFit(False, "too many outliers")
    X, names = design(templates, labels, use, poly)
    q, sw = noise_levels(res, use)
    if not np.isfinite(q):
        return DayFit(False, "noise estimate failed")
    tnames = list(templates)
    k0 = len(names) - len(tnames)
    pinv = np.linalg.pinv(X)
    coef_all = pinv @ y[use]
    t_use = T_SEC[use]
    sig = {n: np.sqrt(ols_variance(pinv[k0 + i], t_use, q, sw)) for i, n in enumerate(tnames)}
    # in the differenced GLS only the per-arc rate (and curvature) columns survive with the templates
    n_offsets = 0
    Xfull_names = names
    offset_cols = [k for k, nm in enumerate(Xfull_names) if nm.endswith("_p0")]
    order = offset_cols + [k for k in range(len(names)) if k not in offset_cols]
    Xo = X[:, order]
    n_offsets = len(offset_cols)
    coef_g, cov_g = gls_differenced(y, Xo, labels, use, q, sw, n_offsets)
    gnames = [Xfull_names[k] for k in order[n_offsets:]]
    gidx = {nm: gnames.index(nm) for nm in tnames}
    out = DayFit(True, n=int(use.sum()), arcs=int(len(np.unique(labels[use]))),
                 rms=float(np.std(res[use])), q=float(q), sw=float(sw),
                 coef={n: float(coef_all[k0 + i]) for i, n in enumerate(tnames)},
                 sigma={n: float(sig[n]) for n in tnames},
                 gls={n: float(coef_g[gidx[n]]) for n in tnames},
                 gls_sigma={n: float(np.sqrt(cov_g[gidx[n], gidx[n]])) for n in tnames})
    if keep_resid:
        out.resid, out.mask = res, use
    return out


# ---------------------------------------------------------------------------
# Building y and the templates from the cached arrays
# ---------------------------------------------------------------------------

def series_for(d: dict, i: int, model: str = "full") -> tuple[np.ndarray, dict]:
    """y and templates for day i of satellite arrays d.

    model "full": y = clk + kc - (mono + j2 + tide + kin); templates alpha (gravity) and
    alpha_tot (gravity + speed) are fitted in separate runs (see run_fits).
    model "kepler": y = clk + kc - (mono + kin): leaves the flattening and tides in, so
    their coefficients (expected 1) can be fitted: does the clock feel the Earth's J2?
    """
    clk, kc = d["clk"][i], d["kc"][i]
    mono, kin = d["mono"][i], d["kin"][i]
    j2, tide = d["j2"][i].astype(float), d["tide"][i].astype(float)
    y_phys = clk + kc
    if model == "full":
        return y_phys - (mono + j2 + tide + kin), {
            "alpha": mono + j2 + tide, "alpha_tot": mono + j2 + tide + kin}
    if model == "kepler":
        return y_phys - (mono + kin), {"j2": j2, "tide": tide, "alpha": mono + j2 + tide}
    if model == "published":
        # the clock as published (Keplerian term removed): what is left of GR in it is the part the
        # IGS convention misses -- J2 bending the orbit, J2's own potential, and the tides
        return clk, {"beyond": mono + kin - kc + j2 + tide}
    raise ValueError(model)


FIT_SETS = {
    # name: (model, templates kept, polynomial degree per arc)
    "alpha": ("full", ("alpha",), 1),
    "alpha_tot": ("full", ("alpha_tot",), 1),
    "alpha_quad": ("full", ("alpha",), 2),
    "j2_tide": ("kepler", ("j2", "tide", "alpha"), 1),
    "beyond": ("published", ("beyond",), 1),
}


def run_fits(prn: str, fit_set: str = "alpha", keep_resid: bool = False) -> tuple[pd.DataFrame, dict]:
    d = prepare.load(prn)
    model, keep, poly = FIT_SETS[fit_set]
    rows, resid = [], {}
    for i, date in enumerate(d["date"]):
        y, tpl = series_for(d, i, model)
        tpl = {k: tpl[k] for k in keep}
        f = fit_day(y, tpl, poly=poly, keep_resid=keep_resid)
        row = {"prn": prn, "date": date, "series": d["series"][i], "ok": f.ok, "reason": f.reason,
               "n": f.n, "arcs": f.arcs, "rms_ps": f.rms * 1e12, "q": f.q, "sw_ps": f.sw * 1e12,
               "beta_deg": float(np.degrees(np.nanmean(d["beta"][i]))),
               "shadow_epochs": int(d["shadow"][i].sum()),
               "e": float(np.nanmean(d["e"][i])), "a_km": float(np.nanmean(d["a"][i])) / 1e3}
        if f.ok:
            for k in keep:
                row[k] = f.coef[k]
                row[k + "_sig"] = f.sigma[k]
                row[k + "_gls"] = f.gls[k]
                row[k + "_gls_sig"] = f.gls_sigma[k]
            if keep_resid:
                resid[date] = (f.resid.astype(np.float32), f.mask)
        rows.append(row)
    return pd.DataFrame(rows), resid


# ---------------------------------------------------------------------------
# Combining days
# ---------------------------------------------------------------------------

def weighted_mean(x: np.ndarray, s: np.ndarray) -> tuple[float, float, float]:
    w = 1.0 / s**2
    m = np.sum(w * x) / np.sum(w)
    err = 1.0 / np.sqrt(np.sum(w))
    chi2 = float(np.sum(w * (x - m) ** 2) / max(len(x) - 1, 1))
    return float(m), float(err), chi2


def block_bootstrap(dates: np.ndarray, x: np.ndarray, s: np.ndarray, *, block_days: int = 30,
                    n_boot: int = 2000, seed: int = C.RANDOM_SEED) -> float:
    """Spread of the weighted mean when whole calendar blocks of days are resampled."""
    rng = np.random.default_rng(seed)
    t = pd.to_datetime(pd.Series(dates)).values.astype("datetime64[D]").astype(np.int64)
    blocks = (t - t.min()) // block_days
    ub = np.unique(blocks)
    groups = [np.flatnonzero(blocks == b) for b in ub]
    w = 1.0 / s**2
    sw_b = np.array([w[g].sum() for g in groups])
    swx_b = np.array([(w[g] * x[g]).sum() for g in groups])
    means = np.empty(n_boot)
    for k in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        means[k] = swx_b[pick].sum() / sw_b[pick].sum()
    return float(np.std(means))


def combine(df: pd.DataFrame, col: str = "alpha", *, use_gls: bool = False) -> dict:
    sel = df[df["ok"] & df[col].notna()]
    xcol, scol = (col + "_gls", col + "_gls_sig") if use_gls else (col, col + "_sig")
    x, s = sel[xcol].to_numpy(), sel[scol].to_numpy()
    m, err, chi2 = weighted_mean(x, s)
    boot = block_bootstrap(sel["date"].to_numpy(), x, s, block_days=C.BOOTSTRAP_BLOCK_DAYS)
    return {"value": m, "err_model": err, "chi2_nu": chi2, "err_model_scaled": err * np.sqrt(max(chi2, 1)),
            "err_boot": boot, "n_days": int(len(sel)),
            "first": str(sel["date"].min()), "last": str(sel["date"].max())}
