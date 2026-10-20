"""Measure frame dragging from the weekly node jumps of LAGEOS and LARES 2.

Each weekly jump J_s(k) (satellite s, between week k-1 and week k) is the part
of the orbital plane's motion over that interval that our model, which leaves
frame dragging out, does not explain. We model the jumps of both satellites
together:

    J_s(k) = mu * LT_s * dt_k                       frame dragging (common scale mu)
           + beta * Omega_J2,s * dt_k               an error in the Earth's flattening J2
           + sum_f [a_sf (cos th_sf(t_k) - cos th_sf(t_k-1)) + b_sf (sin ...)]
                                                    optional periodic signals per satellite
           + noise

LT_s = 2 G J / (c^2 a^3 (1 - e^2)^1.5) is GR's Lense-Thirring node rate (30.7
mas/yr for both satellites). Omega_J2,s is the bulge-driven node rate, +125 deg/yr
for LAGEOS and -125 deg/yr for LARES 2: an error in J2 pushes the two nodes in
opposite directions and frame dragging pushes them the same way, which is what
lets the pair separate the two. The periodic terms are the arguments at which
errors in the tide model and in the radiation-pressure model would appear in a
node: the satellite's node angle (K1), twice it (K2), its angle from the Sun
(radiation pressure, S1/P1) and twice that (S2).

Noise: if each weekly fit's node is off by an independent error, consecutive
jumps are anti-correlated (an "MA(1)" pattern); if instead our force model makes
an error each week, the jumps are independent. We fit a mixture of the two,
with variances estimated by maximum likelihood for each satellite.
"""

from __future__ import annotations

import erfa
import numpy as np
from scipy.optimize import minimize

import config as C
import nodes

MAS_PER_RAD = 180 / np.pi * 3600e3


def lt_rate_mas_yr(a_m: float, e: float) -> float:
    """GR's Lense-Thirring node rate for an orbit around the Earth (mas/yr)."""
    rate = 2 * C.GM_EARTH * C.EARTH_SPIN_PER_MASS / (C.C_LIGHT ** 2 * a_m ** 3 * (1 - e * e) ** 1.5)
    return rate * C.JULIAN_YEAR_S * MAS_PER_RAD


def j2_node_rate_mas_yr(a_m: float, e: float, inc: float, j2: float = 1.0826e-3) -> float:
    n = np.sqrt(C.GM_EARTH / a_m ** 3)
    rate = -1.5 * n * j2 * (C.R_EARTH / a_m) ** 2 * np.cos(inc) / (1 - e * e) ** 2
    return rate * C.JULIAN_YEAR_S * MAS_PER_RAD


def sun_longitude(mjd) -> np.ndarray:
    """Mean longitude of the Sun (rad), adequate for the arguments of periodic terms."""
    return np.radians(280.460 + 0.9856474 * (np.asarray(mjd) - 51_544.5))


def arc_angles(sat: str, ac: str, tag: str = "") -> dict:
    """Node angle (true equator, CIRS) and mean elements at each arc start."""
    out = {"mjd": [], "node": [], "a": [], "e": [], "i": []}
    for a in nodes.load_arcs(sat, ac, tag):
        t = float(a["mjd0"])
        y = a["y_week"]
        x, yy = erfa.xy06(2_400_000.5, t)
        s = erfa.s06(2_400_000.5, t, x, yy)
        Q = erfa.c2ixys(x, yy, s)            # CIRS <- GCRS
        rc = Q @ a["y0"][:3]
        vc = Q @ a["y0"][3:]
        h = np.cross(rc, vc)
        out["mjd"].append(t)
        out["node"].append(np.arctan2(h[0], -h[1]))
        from frames import kepler_elements
        el = kepler_elements(y[:, :3], y[:, 3:])
        out["a"].append(np.mean(el["a"]))
        out["e"].append(np.mean(el["e"]))
        out["i"].append(np.arccos(h[2] / np.linalg.norm(h)))
    return {k: np.array(v) for k, v in out.items()}


def series(ac: str = C.PRIMARY_AC, tag: str = "", sats=("lageos1", "lares2")) -> dict:
    """Jumps of both satellites on their common week boundaries, in mas."""
    js = {s: nodes.jumps(s, ac, tag, max_gap_days=14.01) for s in sats}
    ang = {s: arc_angles(s, ac, tag) for s in sats}
    common = js[sats[0]]["mjd"]
    for s in sats[1:]:
        common = np.intersect1d(np.round(common, 6), np.round(js[s]["mjd"], 6))
    out = {"mjd": common, "sats": sats}
    for s in sats:
        idx = np.searchsorted(np.round(js[s]["mjd"], 6), common)
        gap = js[s]["gap"][idx]
        out[s] = {"jump": js[s]["dnode"][idx] * MAS_PER_RAD, "jump_inc": js[s]["dinc"][idx] * MAS_PER_RAD,
                  "gap_days": gap, "mjd_prev": common - gap,
                  "a": float(np.median(ang[s]["a"])), "e": float(np.median(ang[s]["e"])),
                  "i": float(np.median(ang[s]["i"])),
                  "node_at": lambda m, s=s: np.interp(m, ang[s]["mjd"], np.unwrap(ang[s]["node"]))}
    return out


NUISANCE_SETS = {
    "none": (),
    "K1": ("K1",),
    "K1+S2": ("K1", "S2"),
    "tides": ("K1", "K2", "S2"),
    "tides+srp": ("K1", "K2", "S2", "S1"),
    "full": ("K1", "K2", "S2", "S1", "P1", "SA"),
}


def _argument(kind: str, node, lam):
    return {"K1": node, "K2": 2 * node, "S2": 2 * node - 2 * lam, "S1": node - lam,
            "P1": node - 2 * lam, "SA": lam}[kind]


def design(ser: dict, nuisance=(), fit_j2: bool = True):
    """Design matrix for the stacked jumps [sat0..., sat1...] (mas)."""
    sats = ser["sats"]
    n = ser["mjd"].size
    cols, names = [], []
    rows = []
    for s in sats:
        d = ser[s]
        dt_yr = d["gap_days"] / 365.25
        rows.append({"lt": lt_rate_mas_yr(d["a"], d["e"]) * dt_yr,
                     "j2": j2_node_rate_mas_yr(d["a"], d["e"], d["i"]) * dt_yr * 1e-6})
    names.append("mu_LT")
    cols.append(np.concatenate([r["lt"] for r in rows]))
    if fit_j2:
        names.append("dJ2/J2 [1e-6]")
        cols.append(np.concatenate([r["j2"] for r in rows]))
    for si, s in enumerate(sats):
        d = ser[s]
        t1, t0 = ser["mjd"], d["mjd_prev"]
        for kind in nuisance:
            th1 = _argument(kind, d["node_at"](t1), sun_longitude(t1))
            th0 = _argument(kind, d["node_at"](t0), sun_longitude(t0))
            for fn, lab in ((np.cos, "cos"), (np.sin, "sin")):
                col = np.zeros(n * len(sats))
                col[si * n:(si + 1) * n] = fn(th1) - fn(th0)
                cols.append(col)
                names.append(f"{s}:{kind}:{lab}")
    return np.array(cols).T, names


def noise_cov(n: int, s_white: float, s_rw: float, s_ar: float = 0.0, rho: float = 0.0) -> np.ndarray:
    """Covariance of n consecutive jumps.

    s_white: independent error of each week's fitted node (enters two consecutive
             jumps with opposite signs: an MA(1) pattern);
    s_rw:    independent error of each week's model (white jumps);
    s_ar, rho: slowly varying unmodelled forces (AR(1) jumps, correlation rho per week).
    """
    Cm = np.eye(n) * (2 * s_white ** 2 + s_rw ** 2)
    i = np.arange(n - 1)
    Cm[i, i + 1] = Cm[i + 1, i] = -s_white ** 2
    if s_ar > 0:
        k = np.arange(n)
        Cm = Cm + s_ar ** 2 * rho ** np.abs(k[:, None] - k[None, :])
    return Cm


def fit(ser: dict, nuisance=(), fit_j2: bool = True, noise: str = "ar1",
        channels: str = "sumdiff") -> dict:
    """Generalised least squares for mu (and nuisances), noise variances by REML.

    channels="sumdiff" rotates the two satellites' jumps into their sum and
    difference, (J1 + J2)/sqrt2 and (J1 - J2)/sqrt2, and gives each channel its own
    noise model: an error in J2 is common to both satellites with opposite signs,
    so it is large in the difference and nearly absent in the sum, where frame
    dragging lives. channels="sats" keeps one noise model per satellite.
    The covariance is block-diagonal (one block per channel); whitening is done
    block by block with Cholesky factors.
    """
    from scipy.linalg import solve_triangular

    sats = ser["sats"]
    n = ser["mjd"].size
    y = np.concatenate([ser[s]["jump"] for s in sats])
    X, names = design(ser, nuisance, fit_j2)
    y_orig, X_orig = y.copy(), X.copy()
    if channels == "sumdiff" and len(sats) == 2:
        y = np.concatenate([y[:n] + y[n:], y[:n] - y[n:]]) / np.sqrt(2)
        X = np.concatenate([X[:n] + X[n:], X[:n] - X[n:]]) / np.sqrt(2)
        chan_names = ("sum", "difference")
    else:
        chan_names = sats

    npar = 4 if noise == "ar1" else 2

    def cov_k(p, k):
        q = p[npar * k: npar * (k + 1)]
        if noise == "ar1":
            return noise_cov(n, np.exp(q[0]), np.exp(q[1]), np.exp(q[2]), 1 / (1 + np.exp(-q[3])))
        return noise_cov(n, np.exp(q[0]), np.exp(q[1]))

    def whiten(p):
        Xw = np.empty_like(X)
        yw = np.empty_like(y)
        logdet = 0.0
        for k in range(len(sats)):
            sl = slice(k * n, (k + 1) * n)
            Lk = np.linalg.cholesky(cov_k(p, k))
            Xw[sl] = solve_triangular(Lk, X[sl], lower=True)
            yw[sl] = solve_triangular(Lk, y[sl], lower=True)
            logdet += 2 * np.sum(np.log(np.diag(Lk)))
        return Xw, yw, logdet

    def gls(p):
        Xw, yw, logdet = whiten(p)
        beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        r = yw - Xw @ beta
        XtX = Xw.T @ Xw
        return beta, XtX, r, logdet

    def nll(p):
        beta, XtX, r, logdet = gls(p)
        return 0.5 * (r @ r + logdet + np.linalg.slogdet(XtX)[1])     # REML

    lo, hi = np.log(1e-4), np.log(30.0)
    if noise == "ar1":
        starts = [[np.log(0.3), np.log(0.5), np.log(1.0), 1.0],
                  [np.log(0.05), np.log(1.0), np.log(1.5), 2.0],
                  [np.log(0.5), np.log(0.05), np.log(1.5), 3.0]]
        bounds = [(lo, hi), (lo, hi), (lo, hi), (-4.0, 6.0)] * len(sats)
    else:
        starts = [[np.log(0.3), np.log(0.3)], [np.log(1.0), np.log(0.05)], [np.log(0.05), np.log(1.0)]]
        bounds = [(lo, hi)] * (2 * len(sats))
    best = None
    for p0 in starts:
        res = minimize(nll, np.array(p0 * len(sats)), method="L-BFGS-B", bounds=bounds)
        if best is None or res.fun < best.fun:
            best = res
    p_best = best.x
    beta, XtX, r, _ = gls(p_best)
    cov = np.linalg.inv(XtX)
    err = np.sqrt(np.diag(cov))
    resid = y_orig - X_orig @ beta
    out = {"names": names, "beta": beta, "err": err, "cov": cov,
           "mu": float(beta[0]), "mu_err": float(err[0]),
           "noise": {s: ({"node_error_mas": float(np.exp(p_best[npar * k])),
                          "weekly_model_error_mas": float(np.exp(p_best[npar * k + 1]))}
                         | ({"slow_error_mas": float(np.exp(p_best[npar * k + 2])),
                             "rho": float(1 / (1 + np.exp(-p_best[npar * k + 3])))}
                            if noise == "ar1" else {}))
                     for k, s in enumerate(chan_names)},
           "reml": float(best.fun), "noise_model": noise, "channels": channels,
           "resid": resid, "y": y_orig, "X": X_orig, "n": n, "chi2_white": float(r @ r),
           "dof": int(y.size - X.shape[1])}
    if fit_j2:
        out["dj2_rel"] = float(beta[1]) * 1e-6
        out["dj2_rel_err"] = float(err[1]) * 1e-6
    return out


def cumulative(ser: dict, result: dict | None = None, keep=("mu_LT",)) -> dict:
    """Cumulative node residuals per satellite and their sum (mas).

    With a fit result, every fitted term except those in ``keep`` is subtracted
    from the jumps first (e.g. the J2 error and periodic terms), so what remains
    is frame dragging plus noise.
    """
    sats = ser["sats"]
    n = ser["mjd"].size
    y = np.concatenate([ser[s]["jump"] for s in sats])
    if result is not None:
        drop = np.array([nm not in keep for nm in result["names"]])
        y = y - result["X"][:, drop] @ result["beta"][drop]
    out = {s: np.cumsum(y[k * n:(k + 1) * n]) for k, s in enumerate(sats)}
    out["sum"] = sum(out[s] for s in sats)
    out["mjd"] = ser["mjd"]
    return out


def j2_from_difference(ser: dict, smooth_days: float = 35.0) -> dict:
    """The Earth's J2 error, week by week, from LAGEOS minus LARES 2.

    Frame dragging moves both nodes alike and cancels in the difference; a J2
    error moves them oppositely. beta_k = D_k / ((rate_J2,1 - rate_J2,2) dt) is the
    fractional J2 error of our model over week k. A Gaussian smoothing over about a
    month suppresses the weekly noise. Returns MJD, beta (raw, smoothed) and the
    matching correction to the normalised C20 (dCbar20 = -J2 beta / sqrt 5).
    """
    s1, s2 = ser["sats"]
    d1, d2 = ser[s1], ser[s2]
    rate = (j2_node_rate_mas_yr(d1["a"], d1["e"], d1["i"])
            - j2_node_rate_mas_yr(d2["a"], d2["e"], d2["i"]))           # mas/yr per unit beta
    lt = (lt_rate_mas_yr(d1["a"], d1["e"]) - lt_rate_mas_yr(d2["a"], d2["e"]))
    dt_yr = d1["gap_days"] / 365.25
    D = d1["jump"] - d2["jump"] - lt * dt_yr          # remove GR's (tiny) frame-dragging difference
    beta = D / (rate * dt_yr)
    t_mid = ser["mjd"] - d1["gap_days"] / 2
    w = np.exp(-0.5 * ((t_mid[:, None] - t_mid[None, :]) / (smooth_days / 2.355)) ** 2)
    beta_s = (w @ beta) / w.sum(axis=1)
    j2 = 1.0826e-3
    return {"mjd": t_mid, "beta": beta, "beta_smooth": beta_s,
            "dcbar20": -j2 * beta_s / np.sqrt(5.0), "smooth_days": smooth_days}


#: The adopted fit ("baseline") has a periodic term at every argument where errors of the
#: tide and radiation-pressure models must appear: K1, K2, S2 and P1 (tides), S1 (radiation
#: pressure, node minus Sun) and the annual term. With them the summed residuals are white
#: (no slowly varying component left); without P1 and S1 (significant at ~6 and ~3 sigma)
#: those signals end up in the noise model and inflate the error. The other sets are
#: kept as systematics.
VARIANTS = {
    # name: (nuisance set, noise model, channels)
    "baseline": ("full", "ar1", "sumdiff"),
    "none": ("none", "ar1", "sumdiff"),
    "K1": ("K1", "ar1", "sumdiff"),
    "K1+S2": ("K1+S2", "ar1", "sumdiff"),
    "tides": ("tides", "ar1", "sumdiff"),
    "tides+srp": ("tides+srp", "ar1", "sumdiff"),
    "baseline, white jumps": ("full", "ml", "sumdiff"),
    "baseline, per satellite": ("full", "ar1", "sats"),
}


def _run_variant(args):
    ac, tag, name, (nuis, noise, channels) = args
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:
        pass
    ser = series(ac, tag)
    f = fit(ser, NUISANCE_SETS[nuis], noise=noise, channels=channels)
    return name, {"mu": f["mu"], "mu_err": f["mu_err"], "dj2_rel": f.get("dj2_rel"),
                  "dj2_rel_err": f.get("dj2_rel_err"), "reml": f["reml"], "noise": f["noise"],
                  "params": {n: [float(b), float(e)] for n, b, e in zip(f["names"], f["beta"], f["err"])},
                  "nuisance": nuis, "noise_model": noise, "channels": channels, "dof": f["dof"]}


def analyse(ac: str = C.PRIMARY_AC, tag: str = C.BASE_TAG, variants: dict | None = None,
            force: bool = False, workers: int = 8) -> dict:
    """Fit every variant for one analysis centre; cached in data/cache/framedrag_<ac><tag>.json."""
    import json
    import multiprocessing as mp
    path = C.CACHE_DIR / f"framedrag_{ac}{tag}.json"
    variants = variants or VARIANTS
    old = json.loads(path.read_text()) if path.exists() and not force else {"variants": {}}
    todo = [(ac, tag, k, v) for k, v in variants.items() if k not in old["variants"]]
    ser = series(ac, tag)
    gap = ser[ser["sats"][0]]["gap_days"]
    S = sum(ser[s]["jump"] for s in ser["sats"])
    out = {"ac": ac, "tag": tag, "n_weeks": int(ser["mjd"].size),
           "first": float(ser["mjd"][0]), "last": float(ser["mjd"][-1]),
           "lt_rate": {s: lt_rate_mas_yr(ser[s]["a"], ser[s]["e"]) for s in ser["sats"]},
           "raw_rate": {s: float(ser[s]["jump"].sum() / gap.sum() * 365.25) for s in ser["sats"]},
           "raw_sum_rate": float(S.sum() / gap.sum() * 365.25),
           "elements": {s: {"a_km": ser[s]["a"] / 1e3, "e": ser[s]["e"], "i_deg": float(np.degrees(ser[s]["i"]))}
                        for s in ser["sats"]},
           "variants": dict(old["variants"])}
    if todo:
        ctx = mp.get_context("fork")
        with ctx.Pool(min(workers, len(todo))) as pool:
            for name, res in pool.imap_unordered(_run_variant, todo):
                out["variants"][name] = res
                print(f"    {ac}{tag} {name:22s} mu = {res['mu']:.3f} +- {res['mu_err']:.3f}")
    path.write_text(json.dumps(out, indent=1, default=float))
    return out


def curves(ac: str = C.PRIMARY_AC, tag: str = C.BASE_TAG, variant: str = "baseline") -> dict:
    """Cumulative node residuals (mas) for the figures, from a cached analyse() result.

    'raw': the summed jumps as measured; 'clean': with the fitted J2 offset and
    periodic terms removed (frame dragging and noise remain); 'model': mu * LT * t.
    """
    import json
    res = json.loads((C.CACHE_DIR / f"framedrag_{ac}{tag}.json").read_text())
    v = res["variants"][variant]
    ser = series(ac, tag)
    X, names = design(ser, NUISANCE_SETS[v["nuisance"]], True)
    beta = np.array([v["params"][nm][0] for nm in names])
    sats = ser["sats"]
    n = ser["mjd"].size
    y = np.concatenate([ser[s]["jump"] for s in sats])
    nuis = np.array([nm != "mu_LT" for nm in names])
    clean = y - X[:, nuis] @ beta[nuis]
    gap_yr = ser[sats[0]]["gap_days"] / 365.25
    t_yr = np.cumsum(gap_yr)
    out = {"mjd": ser["mjd"], "t_yr": t_yr, "mu": v["mu"], "mu_err": v["mu_err"],
           "lt": {s: lt_rate_mas_yr(ser[s]["a"], ser[s]["e"]) for s in sats}}
    for k, s in enumerate(sats):
        sl = slice(k * n, (k + 1) * n)
        out[s] = {"raw": np.cumsum(y[sl]), "clean": np.cumsum(clean[sl]),
                  "nuisance": np.cumsum((X[:, nuis] @ beta[nuis])[sl])}
    out["sum"] = {kk: out[sats[0]][kk] + out[sats[1]][kk] for kk in ("raw", "clean", "nuisance")}
    out["lt_sum"] = sum(out["lt"].values())
    return out


def j2_comparison(ac: str = C.PRIMARY_AC, tag: str = "") -> dict:
    """The pair's J2 (as normalised C20) against CSR's monthly SLR C20 and GOCO06s.

    Uses the pass-1 fits (tag "") where our model's C20 is GOCO06s as published, so
    C20(pair) = C20(GOCO06s, t) - J2 * beta / sqrt 5.
    """
    import gravity
    from frames import mjd_to_year
    ser = series(ac, tag)
    j = j2_from_difference(ser)
    goco = np.array([gravity.load(C.GRAVITY_MODEL, 2, mjd=m).cbar[2, 0] for m in j["mjd"]])
    pair = goco - 1.0826e-3 * j["beta"] / np.sqrt(5.0)
    pair_s = goco - 1.0826e-3 * j["beta_smooth"] / np.sqrt(5.0)
    csr = np.loadtxt(C.GRAVITY_DIR / "C20_RL06.txt", comments="#", usecols=(0, 1, 3))
    yr = mjd_to_year(j["mjd"])
    sel = (csr[:, 0] >= yr.min() - 0.05) & (csr[:, 0] <= yr.max() + 0.05)
    csr_on = np.interp(yr, csr[:, 0], csr[:, 1])
    # monthly means of the pair, for a like-for-like comparison with CSR's months
    months = csr[sel, 0]
    pm = np.array([np.mean(pair[np.abs(yr - m) < 0.042]) if np.any(np.abs(yr - m) < 0.042) else np.nan
                   for m in months])
    ok = np.isfinite(pm)
    d = pm[ok] - csr[sel, 1][ok]
    return {"year": yr, "pair": pair, "pair_smooth": pair_s, "goco": goco,
            "csr_year": csr[sel, 0], "csr": csr[sel, 1], "csr_err": csr[sel, 2] * 1e-10,
            "pair_monthly": pm, "offset": float(np.mean(d)), "offset_err": float(np.std(d) / np.sqrt(d.size)),
            "scatter": float(np.std(d)),
            "corr": float(np.corrcoef(pm[ok] - pm[ok].mean(), csr[sel, 1][ok] - csr[sel, 1][ok].mean())[0, 1]),
            "goco_minus_csr": float(np.mean(goco - csr_on))}


def span_test(ac: str = C.PRIMARY_AC, tag: str = C.BASE_TAG, days: float = 1050.0) -> dict:
    """The baseline fit on the first and the last `days` of data (one K1 period each)."""
    ser = series(ac, tag)
    out = {}
    for name, sel in (("first", ser["mjd"] <= ser["mjd"][0] + days),
                      ("last", ser["mjd"] >= ser["mjd"][-1] - days)):
        sub = {"mjd": ser["mjd"][sel], "sats": ser["sats"]}
        for s in ser["sats"]:
            d = dict(ser[s])
            for k in ("jump", "jump_inc", "gap_days", "mjd_prev"):
                d[k] = ser[s][k][sel]
            sub[s] = d
        f = fit(sub, NUISANCE_SETS["full"], noise="ar1")
        out[name] = {"mu": f["mu"], "mu_err": f["mu_err"], "n_weeks": int(sel.sum()),
                     "from": float(sub["mjd"][0]), "to": float(sub["mjd"][-1])}
    return out


def inclination_rates(ac: str = C.PRIMARY_AC, tag: str = C.BASE_TAG) -> dict:
    """Mean inclination jump rate per satellite (mas/yr): frame dragging does not change it."""
    ser = series(ac, tag)
    gap = ser[ser["sats"][0]]["gap_days"]
    return {s: float(ser[s]["jump_inc"].sum() / gap.sum() * 365.25) for s in ser["sats"]}


def inclination_fit(ac: str = C.PRIMARY_AC, tag: str = C.BASE_TAG,
                    nuisance=("K1", "K2", "S2", "S1", "P1", "SA")) -> dict:
    """Constant inclination drift per satellite after the same periodic terms (white jumps).

    Frame dragging by the Earth's spin does not change inclinations, so what remains is a
    measure of forces our model misses (thermal thrust, radiation pressure, tides).
    """
    ser = series(ac, tag)
    out = {}
    for s in ser["sats"]:
        d = ser[s]
        dt_yr = d["gap_days"] / 365.25
        cols = [dt_yr]
        t1, t0 = ser["mjd"], d["mjd_prev"]
        for kind in nuisance:
            th1 = _argument(kind, d["node_at"](t1), sun_longitude(t1))
            th0 = _argument(kind, d["node_at"](t0), sun_longitude(t0))
            cols += [np.cos(th1) - np.cos(th0), np.sin(th1) - np.sin(th0)]
        X = np.array(cols).T
        y = d["jump_inc"]
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ beta
        cov = np.linalg.inv(X.T @ X) * np.var(r, ddof=X.shape[1])
        rho = np.sum(r[1:] * r[:-1]) / np.sum(r * r)
        infl = np.sqrt((1 + rho) / (1 - rho)) if rho > 0 else 1.0     # AR(1) error inflation
        out[s] = {"rate_mas_yr": float(beta[0]), "err": float(np.sqrt(cov[0, 0]) * infl),
                  "raw_rate": float(y.sum() / d["gap_days"].sum() * 365.25), "lag1": float(rho)}
    return out


def series_combined(acs, tag: str = "_c20", min_centres: int = 4) -> dict:
    """Jumps averaged over analysis centres at each common week boundary.

    Every centre fits the same laser ranges with different software and models; averaging
    their weekly jumps reduces what is specific to one centre (frame realisation, models)
    but not the noise common to the data. Elements and node angles are taken from the first.
    """
    sers = [series(ac, tag) for ac in acs]
    allm = np.unique(np.round(np.concatenate([s["mjd"] for s in sers]), 6))
    base = sers[0]
    out = {"mjd": [], "sats": base["sats"]}
    rows = {s: {"jump": [], "jump_inc": [], "gap_days": [], "mjd_prev": [], "n": []} for s in base["sats"]}
    for m in allm:
        vals = {s: [] for s in base["sats"]}
        incs = {s: [] for s in base["sats"]}
        gaps = []
        for sr in sers:
            k = np.where(np.isclose(sr["mjd"], m, atol=1e-5))[0]
            if k.size == 0:
                continue
            k = k[0]
            g = sr[base["sats"][0]]["gap_days"][k]
            if abs(g - 7.0) > 0.01:          # keep plain one-week jumps only
                continue
            gaps.append(g)
            for s in base["sats"]:
                vals[s].append(sr[s]["jump"][k])
                incs[s].append(sr[s]["jump_inc"][k])
        if len(gaps) < min_centres:
            continue
        out["mjd"].append(m)
        for s in base["sats"]:
            rows[s]["jump"].append(np.mean(vals[s]))
            rows[s]["jump_inc"].append(np.mean(incs[s]))
            rows[s]["gap_days"].append(7.0)
            rows[s]["mjd_prev"].append(m - 7.0)
            rows[s]["n"].append(len(gaps))
    out["mjd"] = np.array(out["mjd"])
    for s in base["sats"]:
        d = {k: np.array(v) for k, v in rows[s].items()}
        for k in ("a", "e", "i", "node_at"):
            d[k] = base[s][k]
        out[s] = d
    return out


def de_sitter_node_rate(ser: dict) -> dict:
    """GR's de Sitter node rate of each satellite over the data (mas/yr), from its node geometry.

    The orbit's angular momentum h turns as dh/dt = w x h, w = 19.2 mas/yr along the ecliptic
    pole (GCRS axes); for each week the node rate follows from h's geometry.
    """
    import literature as L
    eps = np.radians(23.4392911)
    w = L.DE_SITTER_RATE * np.array([0.0, -np.sin(eps), np.cos(eps)])
    out = {}
    for s in ser["sats"]:
        d = ser[s]
        node = d["node_at"](ser["mjd"])
        inc = d["i"]
        h = np.stack([np.sin(inc) * np.sin(node), -np.sin(inc) * np.cos(node),
                      np.full_like(node, np.cos(inc))], 1)
        hd = np.cross(w, h)
        rate = (-hd[:, 0] * h[:, 1] + hd[:, 1] * h[:, 0]) / (h[:, 0] ** 2 + h[:, 1] ** 2)
        out[s] = float(np.mean(rate))
    out["sum"] = sum(out[s] for s in ser["sats"])
    return out
