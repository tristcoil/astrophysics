"""Fitting the PPN parameter gamma to Cassini's 2002 conjunction Doppler.

The residual fractional frequency of each normal point (``cassini_doppler``;
model with the Shapiro rate switched off) is fitted as

    y = lambda y_GR + trajectory correction + biases + plasma template + noise,

where y_GR is general relativity's Shapiro rate (the 1PN term scales with
lambda = (1 + gamma)/2, so gamma - 1 = 2 (lambda - 1)).

**Trajectory.** The reconstructed trajectories (Horizons; the navigation
team's SM727D) were themselves fitted to Doppler data with gamma = 1, and
differ from each other by 12-15 mm/s along the line of sight. The fit
therefore corrects Cassini's state: a position and velocity offset at t_ref
plus a constant acceleration (RTG heat, solar pressure), 9 parameters,
through the exact partials of the two-way Doppler (the line-of-sight
velocity, and the position through the direction of the line of sight,
which the Earth's motion turns). Over 30 days at 7.4 au the gravity-gradient
terms of such a correction are < 1e-3 of it. The Shapiro signal is sharp
(~1/(t - t_c) around conjunction) and these are smooth; their correlation
with lambda is reported. Variants add per-track biases (and drifts), which
leave only the within-pass shape of the signal.

The base is the navigation team's SCE1 short-arc trajectory (June 9 - July
5): with the same 9-parameter correction it fits the data better than
Horizons' reconstruction (chi^2 lower by 40 for 884 points, same weights),
needs only a small correction, and gives the same lambda whatever the noise
model. Horizons' trajectory differs from it by a function that the
correction absorbs to 3e-15 rms, but that remainder, weighted by the
data-driven noise model, still moves lambda by ~-1.5e-4 (variants
``horizons_*``); half that difference is quoted as a systematic. The SPK's
own fit assumed GR, but its dynamical model (state, RTG, solar pressure) is
smooth and the correction spans it, so a deviation from GR would still
show in lambda.

The least-squares problems are solved by SVD of the whitened design with
normalised columns (``_solve``): the columns span 20 orders of magnitude and
some trajectory directions are nearly degenerate, and inverting the normal
matrix instead lost information and moved lambda by up to 5e-4.

The fitted acceleration is held physical by a Gaussian prior of 1e-8 m/s^2
per component (the RTG recoil is ~3e-9): left free, it runs to ~3e-8
m/s^2, almost all of it across the line of sight (2e-9 along it), nine times
the RTG recoil (lambda moves by only -0.13e-4).

**Noise.** Each pass gets its own excess (non-white) variance on top of the
white error of its normal points, with AR(1) correlation between the bins of
the pass; and each track a random offset common to its points, whose
variance is a floor plus a fixed fraction of the track's excess variance
(day-to-day errors: the troposphere's daily mean, slow plasma, unmodelled
loading). All are re-estimated from the fit's residuals until stable: the
scintillation of the solar wind grows steeply towards the Sun, differs from
day to day and is red, as is the troposphere's. Because these estimates
come from few tracks, the quoted statistical uncertainty is the largest of
the formal error, a delete-one-track jackknife and a track bootstrap.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import cassini_doppler as D
import cassini_plasma as P

#: trajectory-correction scales: position km, velocity mm/s, acceleration 1e-9 m/s^2
SCALE_R, SCALE_V, SCALE_A = 1e3, 1e-3, 1e-9
#: reference epoch of the trajectory correction: 2002-06-21 12:00 TDB (s since T0)
T_REF = 17.5 * 86400.0
#: the constant excess noise of noise="fixed" (typical of 600-s points far from the Sun)
FIXED_FLOOR = 1.5e-14


# ---------------------------------------------------------------------------
# Data sets
# ---------------------------------------------------------------------------

def rebin(df: pd.DataFrame, bin_s: float) -> pd.DataFrame:
    """Average normal points into longer bins within each pass (weights:
    the number of 1-s points). The white error combines in quadrature."""
    if bin_s <= 60:
        return df.copy()
    d = df.copy()
    d["bbin"] = np.floor(d["t1950"] / bin_s).astype(np.int64)
    num = [c for c in d.columns if c not in ("bbin", "tbin", "station", "ul_band", "dl_band", "link", "pass_id",
                                             "track", "n", "y_err_white", "y_sd", "utc", "source",
                                             "uplink_estimated") and pd.api.types.is_numeric_dtype(d[c])]
    w = d["n"].to_numpy(float)
    keys = ["station", "ul_band", "dl_band", "link", "pass_id", "track", "bbin"]
    for c in num:
        d[c] = d[c] * w
    d["_w"] = w
    d["_e2"] = (d["y_err_white"] * w) ** 2
    g = d.groupby(keys)
    out = g[num + ["_w", "_e2"]].sum()
    for c in num:
        out[c] = out[c] / out["_w"]
    out["y_err_white"] = np.sqrt(out["_e2"]) / out["_w"]
    out["n"] = out["_w"]
    out["n_np"] = g.size()
    out["uplink_estimated"] = g["uplink_estimated"].any()
    out = out.drop(columns=["_w", "_e2"]).reset_index()
    out = out[out["n"] >= 0.5 * bin_s].rename(columns={"bbin": "tbin"})
    out["utc"] = D.O.utc_from_1950(out["t1950"].to_numpy())
    return out.sort_values("t1950").reset_index(drop=True)


def dataset(np_: pd.DataFrame, kind: str = "combo", stations=(25,)) -> pd.DataFrame:
    """'combo' (X/X + Ka/Ka plasma-free, DSS 25), 'kaka', 'xx' (one link at
    the given stations) or 'xx+combo' (the combination where it exists and
    X/X elsewhere, i.e. on the days the Ka uplink was off)."""
    if kind == "combo":
        return np_[(np_["link"] == P.COMBO) & np_["station"].isin(stations)].reset_index(drop=True)
    if kind == "combo_np":
        return P.two_link(np_).reset_index(drop=True)
    if kind == "kaka":
        return np_[(np_["link"] == "Ka/Ka") & np_["station"].isin(stations)].reset_index(drop=True)
    if kind == "xx":
        return np_[(np_["link"] == "X/X") & np_["station"].isin(stations)].reset_index(drop=True)
    if kind == "xx+combo":
        c = np_[np_["link"] == P.COMBO]
        x = np_[(np_["link"] == "X/X") & np_["station"].isin(stations)]
        x = x[~((x["station"] == 25) & x["tbin"].isin(c["tbin"]))]
        return pd.concat([c, x], ignore_index=True).sort_values("t1950").reset_index(drop=True)
    raise ValueError(kind)


# ---------------------------------------------------------------------------
# The linear model
# ---------------------------------------------------------------------------

def design(d: pd.DataFrame, trajectory: str = "state+acc", track_bias: bool = False,
           track_drift: bool = False, station_bias: bool = True, plasma: str | None = None,
           plasma_index: float = 1.3) -> tuple[np.ndarray, list[str]]:
    """Columns: lambda first, then the nuisance parameters."""
    cols, names = [d["y_gr"].to_numpy()], ["lambda"]
    tau = d["t2_tdb"].to_numpy() - T_REF
    dr = d[["dy_dr_x", "dy_dr_y", "dy_dr_z"]].to_numpy()
    dv = d[["dy_dv_x", "dy_dv_y", "dy_dv_z"]].to_numpy()
    if trajectory in ("state", "state+acc"):
        for k, ax in enumerate("xyz"):
            cols.append(dr[:, k] * SCALE_R)
            names.append(f"dr0_{ax} [km]")
        for k, ax in enumerate("xyz"):
            cols.append((dv[:, k] + dr[:, k] * tau) * SCALE_V)
            names.append(f"dv0_{ax} [mm/s]")
    if trajectory == "state+acc":
        for k, ax in enumerate("xyz"):
            cols.append((dv[:, k] * tau + dr[:, k] * tau**2 / 2) * SCALE_A)
            names.append(f"da_{ax} [1e-9 m/s2]")
    if trajectory.startswith("poly"):
        deg = int(trajectory[4:])
        x = tau / (16 * 86400.0)
        for k in range(deg + 1):
            cols.append(x**k * 1e-10)
            names.append(f"poly{k} [1e-10]")
    if station_bias:
        for s in sorted(d["station"].unique())[1:] if len(d["station"].unique()) > 1 else []:
            cols.append((d["station"] == s).to_numpy(float) * 1e-11)
            names.append(f"bias_DSS{s} [1e-11]")
    tracks = sorted(d["track"].unique())
    est = set(d.loc[d["uplink_estimated"].astype(bool), "pass_id"])
    if track_bias:
        for tr in tracks:
            m = (d["track"] == tr).to_numpy()
            cols.append(m * 1e-11)
            names.append(f"bias_{tr} [1e-11]")
            if track_drift:
                tt = (d["t_tdb"].to_numpy() - d.loc[m, "t_tdb"].mean()) / 86400.0
                cols.append(m * tt * 1e-11)
                names.append(f"drift_{tr} [1e-11/day]")
    else:
        for pid in sorted(est):
            cols.append((d["pass_id"] == pid).to_numpy(float) * 1e-11)
            names.append(f"bias_pass{pid}_est_uplink [1e-11]")
    if plasma:
        tmpl = P.plasma_template(d, plasma_index)
        if plasma == "one":
            cols.append(tmpl * 1e-9)
            names.append("plasma")
        elif plasma == "ingress/egress":
            pre = (d["t_tdb"].to_numpy() < 17.5 * 86400.0)
            cols.append(tmpl * pre * 1e-9)
            names.append("plasma_ingress")
            cols.append(tmpl * ~pre * 1e-9)
            names.append("plasma_egress")
    return np.column_stack(cols), names


def _blocks(d: pd.DataFrame) -> list[np.ndarray]:
    return [np.flatnonzero((d["track"] == tr).to_numpy()) for tr in pd.unique(d["track"])]


def _track_cov(idx, e_w, pid, t, ex, rho, s_t2, bin_s):
    """Noise covariance of one track: white errors, each pass's excess noise
    with AR(1) correlation in time, and a common random offset."""
    c = np.diag(e_w[idx] ** 2 * (1 + 1e-9)) + s_t2
    for p in np.unique(pid[idx]):
        m = np.flatnonzero(pid[idx] == p)
        if ex.get(p, 0.0) <= 0:
            continue
        dt = np.abs(t[idx][m][:, None] - t[idx][m][None, :]) / bin_s
        c[np.ix_(m, m)] += ex[p] ** 2 * rho.get(p, 0.0) ** dt
    return c


def _prior_rows(names: list[str], priors: dict | None) -> tuple[np.ndarray, np.ndarray]:
    """Gaussian priors (zero mean) on nuisance parameters, as extra rows of
    the whitened system. ``priors`` maps a name prefix to sigma in the
    parameter's own units, e.g. {"da_": 10.0} = 1e-8 m/s^2."""
    rows = []
    for k, n in enumerate(names):
        for pre, sig in (priors or {}).items():
            if n.startswith(pre) and sig and np.isfinite(sig):
                r = np.zeros(len(names))
                r[k] = 1.0 / sig
                rows.append(r)
    if not rows:
        return np.zeros((0, len(names))), np.zeros(0)
    return np.array(rows), np.zeros(len(rows))


def _solve(Aw: np.ndarray, yw: np.ndarray, rcond: float = 1e-12) -> tuple[np.ndarray, np.ndarray]:
    """Least-squares solution and covariance of a whitened system, by SVD of
    the column-normalised matrix (the raw columns span 20 orders of magnitude
    and some trajectory directions are nearly degenerate; inverting the
    normal matrix instead squares the condition number and silently drops
    information)."""
    scale = np.linalg.norm(Aw, axis=0)
    scale[scale == 0] = 1.0
    U, sv, Vt = np.linalg.svd(Aw / scale, full_matrices=False)
    keep = sv > sv[0] * rcond
    V = Vt[keep].T
    beta = V @ ((U[:, keep].T @ yw) / sv[keep]) / scale
    cov = (V / sv[keep] ** 2) @ V.T / np.outer(scale, scale)
    return beta, cov


def fit(d: pd.DataFrame, n_iter: int = 40, noise: str = "ar1+offset", priors: dict | None = None,
        **design_kw) -> dict:
    """Generalized least squares with a noise model estimated from the data.

    noise='ar1+offset' (default): per pass, an excess noise with AR(1)
    correlation between bins (scintillation and troposphere are red), on top
    of the white errors of the normal points; and per track a random offset
    common to all its points (day-to-day errors: the troposphere's daily
    mean, slow plasma, unmodelled loading) whose variance is
    s_k = s0^2 + (kappa * excess_k)^2, i.e. a floor plus a part that grows
    with the track's own noise (towards the Sun). noise='ar1': no offsets.
    noise='gls': one offset variance for all tracks. noise='white':
    independent errors with a per-pass excess. The hyperparameters are
    re-estimated from the residuals (moments) until stable."""
    from scipy.optimize import nnls

    A, names = design(d, **design_kw)
    P_rows, P_y = _prior_rows(names, priors)
    y = d["y_obs"].to_numpy()
    e_w = d["y_err_white"].to_numpy()
    pid = d["pass_id"].to_numpy()
    t = d["t_tdb"].to_numpy()
    bin_s = float(np.median(np.diff(np.sort(t))[np.diff(np.sort(t)) > 0])) if len(t) > 1 else 60.0
    bin_s = max(60.0, round(bin_s / 60.0) * 60.0)
    blocks = _blocks(d)
    ex = {p: 0.0 for p in np.unique(pid)}
    rho = {p: 0.0 for p in np.unique(pid)}
    s_k = np.zeros(len(blocks))
    hyper = (0.0, 0.0)
    use_corr = noise in ("gls", "ar1", "ar1+offset")
    use_offset = noise in ("gls", "ar1+offset")
    if noise == "fixed":                       # white error + a constant 1.5e-14, nothing estimated
        ex = {p: FIXED_FLOOR for p in ex}
        n_iter = 0
    for it in range(n_iter):
        Aw, yw = [], []
        for k, idx in enumerate(blocks):
            c = _track_cov(idx, e_w, pid, t, ex, rho if use_corr else {}, s_k[k] if use_offset else 0.0, bin_s)
            L = np.linalg.cholesky(c)
            Aw.append(np.linalg.solve(L, A[idx]))
            yw.append(np.linalg.solve(L, y[idx]))
        Aw, yw = np.vstack(Aw + [P_rows]), np.concatenate(yw + [P_y])
        beta, _ = _solve(Aw, yw)
        r = y - A @ beta
        new_ex, new_rho, z, v, e2, use = {}, {}, [], [], [], []
        for k, idx in enumerate(blocks):
            c0 = _track_cov(idx, e_w, pid, t, ex, rho if use_corr else {}, 0.0, bin_s)
            ci1 = np.linalg.solve(c0, np.ones(len(idx)))
            vk = 1.0 / np.sum(ci1)
            mk = vk * float(ci1 @ r[idx])
            z.append(mk**2 - vk)
            v.append(vk)
            e2.append(float(np.mean([ex[p] ** 2 for p in pid[idx]])))
            use.append(len(idx) >= 5)
            rr = r[idx] - (mk if use_offset else 0.0)
            for p in np.unique(pid[idx]):
                q = np.flatnonzero(pid[idx] == p)
                new_ex[p] = float(np.sqrt(max(np.mean(rr[q] ** 2) - np.mean(e_w[idx][q] ** 2), 0.0)))
                if use_corr and len(q) > 6:
                    tq = t[idx][q]
                    adj = np.flatnonzero(np.abs(np.diff(tq) - bin_s) < 0.5 * bin_s)
                    if len(adj) > 4:
                        a, b = rr[q][adj], rr[q][adj + 1]
                        den = np.sqrt(np.sum(a * a) * np.sum(b * b))
                        if den > 0 and np.isfinite(den):
                            new_rho[p] = float(np.clip(np.sum(a * b) / den, 0.0, 0.95))
        z, v, e2, use = np.array(z), np.array(v), np.array(e2), np.array(use)
        new_s = np.zeros(len(blocks))
        if use_offset and use.sum() >= 3 and it > 0:
            wz = 1.0 / (v + s_k)[use]
            if noise == "gls" or not np.any(e2[use] > 0):
                s0 = max(float(np.sum(wz**2 * z[use]) / np.sum(wz**2)), 0.0)
                hyper = (s0, 0.0)
            else:
                X = np.column_stack([np.ones(use.sum()), e2[use]]) * wz[:, None]
                scale = np.max(np.abs(X), axis=0)
                sol, _ = nnls(X / scale, z[use] * wz)
                sol = sol / scale
                hyper = (float(sol[0]), float(sol[1]))
            new_s = hyper[0] + hyper[1] * e2
        conv = (all(abs(new_ex[p] - ex[p]) <= 1e-3 * max(ex[p], 1e-18) for p in ex)
                and np.allclose(new_s, s_k, rtol=1e-3, atol=1e-34)
                and all(abs(new_rho.get(p, 0) - rho.get(p, 0)) < 1e-3 for p in ex))
        ex, s_k = new_ex, new_s
        rho = {p: new_rho.get(p, 0.0) for p in ex} if use_corr else rho
        if conv and it > 2:
            break
    Aw, yw, var_diag = [], [], np.zeros(len(y))
    for k, idx in enumerate(blocks):
        c = _track_cov(idx, e_w, pid, t, ex, rho if use_corr else {}, s_k[k] if use_offset else 0.0, bin_s)
        var_diag[idx] = np.diag(c)
        L = np.linalg.cholesky(c)
        Aw.append(np.linalg.solve(L, A[idx]))
        yw.append(np.linalg.solve(L, y[idx]))
    Aw, yw = np.vstack(Aw + [P_rows]), np.concatenate(yw + [P_y])
    beta, cov = _solve(Aw, yw)
    r = y - A @ beta
    chi2 = float(np.sum((yw[:len(y)] - Aw[:len(y)] @ beta) ** 2))
    dof = len(y) - np.linalg.matrix_rank(A)
    lam, lam_err = float(beta[0]), float(np.sqrt(cov[0, 0]))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov[0] / np.sqrt(np.abs(np.diag(cov)) * cov[0, 0])
    corr = np.nan_to_num(corr)
    tracks = [str(d["track"].iloc[idx[0]]) for idx in blocks]
    return {"names": names, "beta": beta, "cov": cov, "lambda": lam, "lambda_err": lam_err,
            "gamma_minus_1": 2 * (lam - 1), "gamma_err": 2 * lam_err, "chi2": chi2, "dof": int(dof),
            "chi2_dof": chi2 / max(dof, 1), "resid": r, "weights": 1.0 / var_diag, "A": A,
            "excess": pd.Series(ex), "rho": pd.Series(rho),
            "track_offset_sigma": pd.Series(np.sqrt(s_k), index=tracks),
            "offset_floor": float(np.sqrt(hyper[0])), "offset_kappa": float(np.sqrt(hyper[1])),
            "bin_s": bin_s, "noise": noise,
            "n_points": int(len(y)), "n_passes": int(d["pass_id"].nunique()), "n_tracks": int(d["track"].nunique()),
            "max_abs_corr_lambda": float(np.max(np.abs(corr[1:]))) if len(corr) > 1 else 0.0,
            "corr_lambda": dict(zip(names[1:], np.round(corr[1:], 3).tolist()))}


def jackknife(d: pd.DataFrame, unit: str = "track", **kw) -> dict:
    """Delete-one-group jackknife of lambda (groups: track or day)."""
    groups = pd.unique(d[unit])
    lams = []
    for g in groups:
        sub = d[d[unit] != g]
        if sub["track"].nunique() < 3:
            continue
        lams.append(fit(sub, **kw)["lambda"])
    lams = np.array(lams)
    n = len(lams)
    err = float(np.sqrt((n - 1) / n * np.sum((lams - lams.mean()) ** 2)))
    return {"lambda_jk_err": err, "gamma_jk_err": 2 * err, "n_groups": n,
            "per_group": dict(zip([str(g) for g in groups[:n]], lams.tolist()))}


def measured_series(d: pd.DataFrame, f: dict) -> pd.DataFrame:
    """The data with the fitted nuisance terms removed: the frequency shift
    Cassini measured, beside GR's prediction."""
    A = f["A"]
    nuis = A[:, 1:] @ f["beta"][1:]
    out = pd.DataFrame({
        "utc": d["utc"].to_numpy(), "t_tdb": d["t_tdb"].to_numpy(), "b_rsun": d["b_down_rsun"].to_numpy(),
        "b_up_rsun": d["b_up_rsun"].to_numpy(), "station": d["station"].to_numpy(), "link": d["link"].to_numpy(),
        "track": d["track"].to_numpy(), "y_measured": d["y_obs"].to_numpy() - nuis,
        "sigma": 1 / np.sqrt(f["weights"]), "y_gr": d["y_gr"].to_numpy()})
    out["residual"] = out["y_measured"] - f["lambda"] * out["y_gr"]
    return out


def bootstrap(d: pd.DataFrame, n: int = 200, seed: int = 1, **kw) -> dict:
    """Bootstrap over whole tracks (resampled with replacement)."""
    rng = np.random.default_rng(seed)
    tracks = pd.unique(d["track"])
    lams = []
    for _ in range(n):
        pick = rng.choice(tracks, size=len(tracks), replace=True)
        parts = []
        for k, tr in enumerate(pick):
            g = d[d["track"] == tr].copy()
            g["track"] = f"{tr}#{k}"
            g["pass_id"] = g["pass_id"].astype(np.int64) * 1000 + k
            parts.append(g)
        b = pd.concat(parts, ignore_index=True)
        try:
            lams.append(fit(b, **kw)["lambda"])
        except np.linalg.LinAlgError:
            continue
    lams = np.array(lams)
    return {"lambda_bs_err": float(np.std(lams, ddof=1)), "gamma_bs_err": float(2 * np.std(lams, ddof=1)),
            "n_resamples": int(len(lams)), "lambda_bs_median": float(np.median(lams)),
            "lambda_bs_16_84": [float(np.percentile(lams, 16)), float(np.percentile(lams, 84))]}


# ---------------------------------------------------------------------------
# The analysis: main result and variants
# ---------------------------------------------------------------------------

#: the main analysis choices
MAIN = dict(config="nav", el_min=15.0, kind="combo", bin_s=600.0,
            fit=dict(trajectory="state+acc", noise="ar1+offset", priors={"da_": 10.0}))

#: robustness variants: (name, description, changes to MAIN; 'fit' entries
#: update the fit options, 'filter' selects normal points)
VARIANTS = [
    ("kaka_only", "Ka/Ka link alone (no plasma calibration)", dict(kind="kaka")),
    ("kaka_only_plasma_template", "Ka/Ka alone + fitted coronal template (b^-1.3)",
     dict(kind="kaka", fit=dict(plasma="one"))),
    ("xx_only", "X/X link alone at DSS 25 (no plasma calibration)", dict(kind="xx")),
    ("xx_only_plasma_template", "X/X alone + fitted coronal template", dict(kind="xx", fit=dict(plasma="one"))),
    ("combo_plasma_template", "combination + template for its plasma residue", dict(fit=dict(plasma="one"))),
    ("combo_plus_xx_no_ka_template", "combination + DSS 25 X/X where Ka was off (June 18-20), template",
     dict(kind="xx+combo", fit=dict(plasma="one"))),
    ("combo_plus_xx_other_stations", "combination + X/X of DSS 34/45/54/65/15 with b > 15 Rsun (station biases)",
     dict(kind="combo+others")),
    ("bin_300s", "300-s bins", dict(bin_s=300.0)),
    ("bin_1800s", "1800-s bins", dict(bin_s=1800.0)),
    ("noise_white", "independent errors (per-pass excess, no AR(1), no track offsets)",
     dict(fit=dict(noise="white"))),
    ("noise_ar1_no_offsets", "AR(1) per pass, no random track offsets", dict(fit=dict(noise="ar1"))),
    ("noise_common_offset", "AR(1) + random track offsets of one common variance", dict(fit=dict(noise="gls"))),
    ("trajectory_no_acceleration", "trajectory: state only", dict(fit=dict(trajectory="state"))),
    ("trajectory_acceleration_free", "trajectory: state + free acceleration", dict(fit=dict(priors={}))),
    ("trajectory_poly3", "line-of-sight velocity as a cubic in time", dict(fit=dict(trajectory="poly3", priors={}))),
    ("trajectory_poly5", "line-of-sight velocity as a quintic in time", dict(fit=dict(trajectory="poly5", priors={}))),
    ("track_offsets_2e-14", "per-track offset parameters with a 2e-14 Gaussian prior (instead of random offsets)",
     dict(fit=dict(track_bias=True, noise="ar1", priors={"da_": 10.0, "bias_25": 2e-3}))),
    ("track_offsets_free", "free per-track offsets (within-pass shape only)",
     dict(fit=dict(track_bias=True, noise="ar1"))),
    ("b_gt_5", "impact parameter > 5 Rsun", dict(filter="b_down_rsun > 5")),
    ("b_gt_8", "impact parameter > 8 Rsun", dict(filter="b_down_rsun > 8")),
    ("b_gt_12", "impact parameter > 12 Rsun", dict(filter="b_down_rsun > 12")),
    ("b_lt_25", "impact parameter < 25 Rsun", dict(filter="b_down_rsun < 25")),
    ("without_june_21_22", "the two closest tracks removed",
     dict(filter="~track.str.startswith('25_2002-06-21') & ~track.str.startswith('25_2002-06-22')")),
    ("ingress_only", "before conjunction only (June 9-17)", dict(filter="t_tdb < 1512000")),
    ("egress_only", "after conjunction only (June 21 - July 4)", dict(filter="t_tdb >= 1512000")),
    ("elevation_10", "elevation cut 10 deg (AWVR data only above 15)", dict(el_min=10.0)),
    ("elevation_20", "elevation cut 20 deg", dict(el_min=20.0)),
    ("elevation_25", "elevation cut 25 deg", dict(el_min=25.0)),
    ("tropo_dsn", "DSN troposphere calibration instead of the AWVRs", dict(config="nav_tropo_dsn")),
    ("awvr2", "water vapour radiometer unit 2 instead of unit 1", dict(config="nav_awvr2")),
    ("no_ionosphere", "no Earth-ionosphere calibration", dict(config="nav_no_iono")),
    ("no_solid_tides", "station not displaced by the solid-Earth tides (a check: the tides are real)",
     dict(config="nav_no_tides")),
    ("no_2pn", "the 2PN light-time terms left out of the model (a check)", dict(drop_2pn=True)),
    ("horizons_trajectory", "Horizons' reconstructed trajectory instead of the navigation SPK (June 6 - July 5)",
     dict(config="main")),
    ("horizons_trajectory_span", "Horizons' trajectory, same span as the main fit (June 9 - July 5)",
     dict(config="main", filter="t_tdb >= 489100")),
    ("horizons_simple_weights", "Horizons' trajectory, main span, b > 5 Rsun, fixed weights (white + 1.5e-14), "
     "no prior", dict(config="main", filter="b_down_rsun > 5 & t_tdb >= 489100", fit=dict(noise="fixed", priors={}))),
    ("nav_simple_weights", "b > 5 Rsun, fixed weights (white + 1.5e-14), no prior",
     dict(filter="b_down_rsun > 5", fit=dict(noise="fixed", priors={}))),
]

#: variants that are alternative reasonable analysis choices (their spread
#: about the main result is quoted as the systematic uncertainty)
SYSTEMATIC_SET = ["bin_300s", "bin_1800s", "noise_white", "noise_ar1_no_offsets", "track_offsets_2e-14",
                  "trajectory_no_acceleration", "trajectory_acceleration_free", "combo_plasma_template",
                  "elevation_20", "tropo_dsn", "awvr2", "b_gt_8"]
#: the alternative base trajectory, quoted as its own systematic (half the difference)
TRAJECTORY_BASE_VARIANT = "horizons_trajectory_span"


def _merge(base: dict, ch: dict) -> dict:
    out = {k: (dict(v) if isinstance(v, dict) else v) for k, v in base.items()}
    for k, v in ch.items():
        if k == "fit":
            out["fit"] = {**out["fit"], **v}
        else:
            out[k] = v
    return out


def build_dataset(spec: dict, cache: dict, log=print) -> pd.DataFrame:
    key = (spec["config"], spec["el_min"])
    if key not in cache:
        import gc

        cache[key] = D.normal_points(spec["config"], el_min=spec["el_min"], log=log)
        gc.collect()              # the 1-s tables of a model variant are ~1 GB
    np_ = cache[key]
    if spec["kind"] == "combo+others":
        c = dataset(np_, "combo")
        x = np_[(np_["link"] == "X/X") & (np_["station"] != 25) & (np_["b_down_rsun"] > 15)]
        d = pd.concat([rebin(c, spec["bin_s"]), rebin(x.reset_index(drop=True), spec["bin_s"])], ignore_index=True)
        d = d.sort_values("t1950").reset_index(drop=True)
    else:
        d = rebin(dataset(np_, spec["kind"]), spec["bin_s"])
    if spec.get("filter"):
        d = d.query(spec["filter"], engine="python").reset_index(drop=True)
    if spec.get("drop_2pn"):
        # the model without the 2PN terms predicts a frequency higher by y_2pn
        d = d.assign(y_obs=d["y_obs"] + d["y_2pn"])
    return d


def summarize(name: str, desc: str, spec: dict, d: pd.DataFrame, f: dict, jk: dict, bs: dict | None = None) -> dict:
    stat = max(f["lambda_err"], jk["lambda_jk_err"], (bs or {}).get("lambda_bs_err", 0.0))
    out = {"name": name, "description": desc, "dataset": spec["kind"], "model_config": spec["config"],
           "elevation_min_deg": spec["el_min"], "bin_s": spec["bin_s"], "filter": spec.get("filter"),
           "fit_options": {k: v for k, v in spec["fit"].items()},
           "n_points": f["n_points"], "n_passes": f["n_passes"], "n_tracks": f["n_tracks"],
           "chi2": round(f["chi2"], 2), "dof": f["dof"], "chi2_dof": round(f["chi2_dof"], 3),
           "lambda": f["lambda"], "lambda_minus_1": f["lambda"] - 1, "lambda_err_formal": f["lambda_err"],
           "lambda_err_jackknife": jk["lambda_jk_err"],
           "gamma_minus_1": 2 * (f["lambda"] - 1), "gamma_err_formal": 2 * f["lambda_err"],
           "gamma_err_jackknife": 2 * jk["lambda_jk_err"], "gamma_err_stat": 2 * stat,
           "max_abs_correlation_lambda_nuisance": round(f["max_abs_corr_lambda"], 3),
           "b_range_rsun": [float(d["b_down_rsun"].min()), float(d["b_down_rsun"].max())],
           "date_range_utc": [str(d["utc"].min())[:16], str(d["utc"].max())[:16]]}
    if bs:
        out.update({"lambda_err_bootstrap": bs["lambda_bs_err"], "gamma_err_bootstrap": bs["gamma_bs_err"],
                    "bootstrap_resamples": bs["n_resamples"]})
    if "plasma" in f["names"]:
        out["plasma_template_amplitude"] = float(f["beta"][f["names"].index("plasma")])
    return out


def run(rebuild: bool = False, variants: bool = True, n_boot: int = 200, log=print) -> dict:
    """Main fit (with jackknife and track bootstrap), every variant, and the
    outputs: data/cache/cassini_fit.json and cassini_measured_series.parquet."""
    import json

    cache: dict = {}
    if rebuild:
        D.normal_points("main", refresh=True, log=log)
    spec = _merge(MAIN, {})
    d = build_dataset(spec, cache, log)
    f = fit(d, **spec["fit"])
    jk = jackknife(d, **spec["fit"])
    bs = bootstrap(d, n=n_boot, **spec["fit"]) if n_boot else None
    main = summarize("main", "X/X + Ka/Ka plasma-free combination at DSS 25, navigation SPK trajectory + "
                     "state/acceleration correction, 600-s bins, AR(1) noise per pass + random track offsets",
                     spec, d, f, jk, bs)
    main["nuisance"] = {n: float(b) for n, b in zip(f["names"][1:], f["beta"][1:])}
    main["noise"] = {"excess_per_pass": {str(k): float(v) for k, v in f["excess"].items()},
                     "ar1_rho_per_pass": {str(k): float(v) for k, v in f["rho"].items()},
                     "track_offset_sigma": {str(k): float(v) for k, v in f["track_offset_sigma"].items()},
                     "offset_floor": f["offset_floor"], "offset_kappa": f["offset_kappa"]}
    main["jackknife_per_track_lambda"] = jk["per_group"]
    log(f"main: lambda - 1 = {main['lambda_minus_1']:+.2e} +- {main['lambda_err_formal']:.1e} (formal) "
        f"{main['lambda_err_jackknife']:.1e} (jk) {main.get('lambda_err_bootstrap', float('nan')):.1e} (bs)")
    results = {"main": main, "variants": []}
    series = [measured_series(d, f).assign(dataset="main: X/X+Ka/Ka combination")]
    if variants:
        for name, desc, ch in VARIANTS:
            sp = _merge(MAIN, ch)
            try:
                dv = build_dataset(sp, cache, log)
                fv = fit(dv, **sp["fit"])
                jv = jackknife(dv, **sp["fit"])
            except Exception as exc:                       # record, do not stop the run
                results["variants"].append({"name": name, "description": desc, "error": repr(exc)})
                log(f"  {name}: failed ({exc!r})")
                continue
            sv = summarize(name, desc, sp, dv, fv, jv)
            results["variants"].append(sv)
            log(f"  {name:32s} lambda - 1 = {sv['lambda_minus_1']:+.2e} +- {sv['lambda_err_formal']:.1e} "
                f"(jk {sv['lambda_err_jackknife']:.1e})  n={sv['n_points']}")
            if name in ("kaka_only_plasma_template", "xx_only_plasma_template", "kaka_only", "xx_only",
                        "horizons_trajectory"):
                series.append(measured_series(dv, fv).assign(dataset=f"{name}"))
        sys_names = [v for v in results["variants"] if v.get("name") in SYSTEMATIC_SET and "lambda" in v]
        dev = np.array([v["lambda"] - main["lambda"] for v in sys_names])
        syst = float(np.sqrt(np.mean(dev**2))) if len(dev) else float("nan")
        stat = max(main["lambda_err_formal"], main["lambda_err_jackknife"], main.get("lambda_err_bootstrap", 0.0))
        tb = next((v for v in results["variants"] if v.get("name") == TRAJECTORY_BASE_VARIANT and "lambda" in v), None)
        traj = 0.5 * abs(tb["lambda"] - main["lambda"]) if tb else float("nan")
        total = float(np.sqrt(stat**2 + syst**2 + traj**2))
        results["summary"] = {
            "lambda_minus_1": main["lambda_minus_1"], "lambda_err_stat": stat, "lambda_err_syst": syst,
            "lambda_err_trajectory_base": traj, "lambda_err_total": total,
            "gamma_minus_1": main["gamma_minus_1"], "gamma_err_stat": 2 * stat, "gamma_err_syst": 2 * syst,
            "gamma_err_trajectory_base": 2 * traj, "gamma_err_total": 2 * total,
            "systematic_from": [v["name"] for v in sys_names],
            "trajectory_base_from": TRAJECTORY_BASE_VARIANT,
            "notes": ("stat = max(formal, jackknife, bootstrap); syst = rms deviation of the analysis-choice "
                      "variants listed; trajectory_base = half the difference from the Horizons-based fit "
                      "over the same span"),
            "bertotti_2003": {"gamma_minus_1": 2.1e-5, "gamma_err": 2.3e-5},
        }
        log(f"summary: gamma - 1 = {main['gamma_minus_1']:+.2e} +- {2 * stat:.1e} (stat) +- {2 * syst:.1e} (syst) "
            f"+- {2 * traj:.1e} (trajectory base) = +- {2 * total:.1e}")
    out = pd.concat(series, ignore_index=True)
    out = out.rename(columns={"utc": "date_utc"})
    out.to_parquet(D.C.CACHE_DIR / "cassini_measured_series.parquet")
    (D.C.CACHE_DIR / "cassini_fit.json").write_text(json.dumps(results, indent=1, default=float))
    return results
