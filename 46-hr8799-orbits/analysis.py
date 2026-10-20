"""Numbers for the README: reproduction of Wang et al. (2018), model comparison,
fit diagnostics, the prediction test, stability, masses.

Each function returns plain dicts (saved to data/results.json by run.py).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
import diagnostics as D
import fits
import kepler as K
import predict
import published as P

N_GOF = 1000          # posterior draws for chi2_nu and BIC (as Wang et al.)


def q(x, pct=(16, 50, 84)) -> dict:
    lo, med, hi = np.percentile(np.asarray(x, float), pct)
    return {"median": float(med), "plus": float(hi - med), "minus": float(med - lo)}


def angle(x):
    """Wang et al.'s convention: omega and Omega folded into [0, 180)."""
    return K.wrap180(x)


PARAMS = [f"{q_}_{p}" for p in C.PLANETS for q_ in ("a", "tau", "omega", "Omega", "i", "e")] + ["plx", "mstar"]


def posterior_summary(fit: fits.Fit) -> dict:
    t = fit.table()
    out = {}
    for name in PARAMS + ["P_b/P_c", "P_c/P_d", "P_d/P_e"] + [f"P_{p}" for p in C.PLANETS]:
        x = t[name].to_numpy()
        if name.startswith(("omega", "Omega")):
            x = angle(x)
        out[name] = q(x)
    return out


def goodness_of_fit(fit: fits.Fit, n: int = N_GOF) -> dict:
    """chi2_nu and BIC on n random posterior draws (Wang et al. Sec. 3.5)."""
    rng = np.random.default_rng(C.RANDOM_SEED)
    s = fit.equal[rng.choice(len(fit.equal), min(n, len(fit.equal)), replace=False)]
    chi2 = np.array([fit.model.chi2(th) for th in s])
    n_data = 2 * len(fit.model.data.pidx)
    k = fit.model.ndim
    return {"chi2": chi2, "chi2nu": chi2 / (n_data - k), "bic": chi2 + k * np.log(n_data),
            "n_data": n_data, "k": k}


def w18_reproduction(prefix: str = "W18") -> dict:
    """Our four models against Wang et al.'s Table 2."""
    W = P.wang2018_models()
    res = {}
    gof = {}
    for kind in C.MODELS:
        f = fits.load(f"{prefix}/{kind}")
        res[kind] = {"posterior": posterior_summary(f), "logz": f.logz, "logzerr": float(f.raw["logzerr"]),
                     "ess": float(f.raw["ess"]), "max_lnL": float(f.raw["logl"].max())}
        gof[kind] = goodness_of_fit(f)
    bic0 = np.median(gof["unconstrained"]["bic"])
    for kind in C.MODELS:
        g = gof[kind]
        res[kind]["chi2nu"] = q(g["chi2nu"])
        res[kind]["dBIC"] = q(g["bic"] - bic0)
        res[kind]["ln_bayes_vs_unconstrained"] = res[kind]["logz"] - res["unconstrained"]["logz"]
        res[kind]["wang"] = W[kind]
    return res


def pulls(ours: dict, theirs: dict) -> dict:
    """(ours - Wang) / combined sigma, for every parameter both report."""
    out = {}
    for name, w in theirs.items():
        if name not in ours or not isinstance(w, tuple):
            continue
        o = ours[name]
        s1 = 0.5 * (o["plus"] + o["minus"])
        s2 = 0.5 * (w[1] + w[2])
        out[name] = (o["median"] - w[0]) / np.hypot(s1, s2)
    return out


def consistency(prefix: str = "W18") -> dict:
    """Nesting check, repeat-run agreement and box margins for one data set."""
    keys = [f"{prefix}/{k}" for k in C.MODELS]
    loaded = {k: fits.load(k) for k in keys if fits.cache_path(k).exists()}
    report = D.check(loaded)
    out = {}
    for k, r in report.items():
        out[k] = {"posterior_max_lnL": r["posterior_max_lnL"], "polished_max_lnL": r["polished_max_lnL"],
                  "shortfall": r["shortfall"], "ok": bool(r["ok"]), "best_from": r["source_of_best"],
                  "box_margin_min_sigma": float(min(D.box_margins(loaded[k]).values()))}
    # repeat runs
    for k in keys:
        k2 = k + "#2"
        if fits.cache_path(k2).exists() and k in loaded:
            f1, f2 = loaded[k], fits.load(k2)
            t1, t2 = f1.table(), f2.table()
            shifts = {}
            for name in PARAMS:
                x1, x2 = t1[name].to_numpy(), t2[name].to_numpy()
                if name.startswith(("omega", "Omega")):
                    x1, x2 = angle(x1), angle(x2)
                shifts[name] = float((np.median(x1) - np.median(x2)) / np.hypot(np.std(x1), np.std(x2)))
            out[k]["repeat_run"] = {"logz": [f1.logz, f2.logz],
                                    "max_abs_median_shift_sigma": float(np.max(np.abs(list(shifts.values())))),
                                    "worst_param": max(shifts, key=lambda n: abs(shifts[n]))}
    return out


# ---------------------------------------------------------------------------
# Stability (Wang et al. 2018, Sec. 4)
# ---------------------------------------------------------------------------

def _binom(k: int, n: int) -> dict:
    """Rate per million with a 68% interval (Jeffreys) and a 95% upper limit."""
    from scipy.stats import beta
    lo, hi = beta.ppf([0.16, 0.84], k + 0.5, n - k + 0.5)
    up95 = beta.ppf(0.95, k + 0.5, n - k + 0.5)
    return {"k": int(k), "n": int(n), "per_million": 1e6 * k / n, "lo": 1e6 * lo, "hi": 1e6 * hi,
            "upper95": 1e6 * up95}


def stability_summary(key: str, tag: str = "") -> dict:
    import stability as St
    df = St.read_results(St.results_path(key, tag))
    st = df.outcome.eq("stable")
    out = {"draws": int(len(df)), "outcomes": df.outcome.value_counts().to_dict(),
           "stable": _binom(int(st.sum()), len(df)),
           "survival_yr_percentiles": {str(p): float(np.percentile(df.t_end_yr, p)) for p in (10, 50, 90, 99)},
           "cpu_s_mean": float(df.cpu_s.mean())}
    bins = np.arange(4, 12)
    df["bin"] = pd.cut(df.m_cde, bins)
    by = df.groupby("bin", observed=True).outcome.apply(lambda s: (s == "stable").mean() * 1e6)
    out["per_million_by_mcde"] = {str(k): float(v) for k, v in by.items()}
    return out


def stable_posterior(key: str, tag: str = "") -> dict:
    """Elements of the stable orbits (at MJD 56609), as Wang et al.'s Table 3."""
    import stability as St
    f = fits.load(key)
    t = f.table()
    stab = St.stable_systems(key, tag)
    sub = t.iloc[stab["sample"].to_numpy()].reset_index(drop=True)
    out = {}
    for name in PARAMS + ["P_b/P_c", "P_c/P_d", "P_d/P_e"]:
        x = sub[name].to_numpy()
        if name.startswith(("omega", "Omega")):
            x = angle(x)
        out[name] = q(x)
    dw = np.mod(sub["omega_e"].to_numpy() - sub["omega_d"].to_numpy(), 360.0)
    dw = np.where(dw > 180, dw - 360, dw)
    out["omega_e-omega_d"] = q(np.abs(dw))
    m_cde, m_b = stab.m_cde.to_numpy(), stab.m_b.to_numpy()
    out["masses"] = {
        "n_stable": int(len(stab)),
        "frac_mcde_lt6": float(np.mean(m_cde < 6)), "frac_mcde_lt5": float(np.mean(m_cde < 5)),
        "mcde_999": float(np.percentile(m_cde, 99.9)), "mb_999": float(np.percentile(m_b, 99.9)),
        "mcde_max": float(m_cde.max()), "mb_max": float(m_b.max()),
    }
    # hot-start luminosity prior (Wang et al.: Baraffe et al. 2003 tracks at 42 +- 5 Myr)
    from scipy.stats import norm
    w = (norm.pdf(m_b, *C.HOTSTART_MASS["b"]) * norm.pdf(m_cde, *C.HOTSTART_MASS["cde"])
         / (1.0 / (m_cde - C.MASS_B_MIN)))       # divide out the M_b ~ U(3, M_cde) prior density
    w /= w.sum()

    def wq(x):
        o = np.argsort(x)
        c = np.cumsum(w[o])
        lo, med, hi = np.interp([0.16, 0.5, 0.84], c, x[o])
        return {"median": float(med), "plus": float(hi - med), "minus": float(med - lo)}

    out["masses"]["luminosity_weighted_b"] = wq(m_b)
    out["masses"]["luminosity_weighted_cde"] = wq(m_cde)
    out["masses"]["luminosity_weighted_ess"] = float(1.0 / np.sum(w ** 2))
    return out


# ---------------------------------------------------------------------------
# How Keplerian are the orbits over the data span?
# ---------------------------------------------------------------------------

def keplerian_vs_nbody(key: str = "ALL/near1248_coplanar", n: int = 100, m_b: float = 6.0,
                       m_cde: float = 8.7) -> dict:
    """The fits treat each planet as a fixed Keplerian (Jacobi) ellipse. Here the
    planets also pull on each other: start a REBOUND run (IAS15) from each orbit's
    state on MJD NBODY_EPOCH_MJD with hot-start masses, run it forwards and (with
    reversed velocities) backwards to every observing date, and measure how far the
    N-body positions drift from the Keplerian ones with the same masses."""
    import stability as St
    f = fits.load(key)
    rng = np.random.default_rng(C.RANDOM_SEED)
    s = f.equal[rng.choice(len(f.equal), n, replace=False)]
    phys = f.model.physical(s)
    data = __import__("astrometry").dataset("ALL")
    epochs = np.sort(data.mjd.unique())
    t0 = C.NBODY_EPOCH_MJD
    mj = {"b": m_b, "c": m_cde, "d": m_cde, "e": m_cde}
    msun = {p: mj[p] * C.MJUP_MSUN for p in C.INNER_TO_OUTER}
    dev = {p: np.zeros((n, len(epochs))) for p in C.INNER_TO_OUTER}
    for k in range(n):
        el = {p: tuple(float(np.asarray(v)[k]) for v in phys["elements"][p]) for p in C.INNER_TO_OUTER}
        mstar, plx = float(phys["mstar"][k]), float(phys["plx"][k])
        for direction, sel in ((1, epochs >= t0), (-1, epochs < t0)):
            sim = St.build_sim(el, mstar, mj, reverse=direction < 0)
            sim.integrator = "ias15"
            order = np.argsort(np.abs(epochs[sel] - t0))
            idx = np.where(sel)[0][order]
            for j in idx:
                sim.integrate(abs(epochs[j] - t0) / 365.25)
                star = sim.particles[0]
                kep = K.system_offsets(epochs[j], el, plx, mstar, msun)
                for q, p in enumerate(C.INNER_TO_OUTER, start=1):
                    pq = sim.particles[q]
                    dev[p][k, j] = np.hypot((pq.x - star.x) * plx - kep[p][0], (pq.y - star.y) * plx - kep[p][1])
    out = {"masses_mjup": mj, "n_orbits": n, "epochs": [float(epochs.min()), float(epochs.max())]}
    for p in C.PLANETS:
        worst = dev[p].max(axis=1)                 # worst epoch for each orbit
        out[p] = {"median_of_max_mas": float(np.median(worst)), "p90_of_max_mas": float(np.percentile(worst, 90)),
                  "median_rms_mas": float(np.median(np.sqrt((dev[p] ** 2).mean(axis=1))))}
    return out


def goodness_with_nuisances(fit: fits.Fit, n: int = 1000) -> dict:
    """chi2 of the all-data fits counting what the model knows about the instruments:
    each instrument's fitted jitter added to the errors and its best-fitting
    north/scale offsets removed from the residuals (plain chi2 ignores both)."""
    m = fit.model
    d = m.data
    rng = np.random.default_rng(C.RANDOM_SEED)
    s = fit.equal[rng.choice(len(fit.equal), min(n, len(fit.equal)), replace=False)]
    out = []
    for th in s:
        phys = m.physical(th)
        _, offs = m.loglike_phys(phys, return_offsets=True)
        r, A, _ = m.residuals(phys)
        for src, (sol, _) in offs.items():
            rows = d.sidx == m._off_index[src]
            r[rows] -= np.einsum("nij,j->ni", A[rows], sol)
        s0, s1 = d.err[:, 0].copy(), d.err[:, 1].copy()
        if len(d.jitter_sources):
            j = np.where(m._row_jit_idx >= 0, phys["jitter"][np.maximum(m._row_jit_idx, 0)], 0.0)
            s0, s1 = np.sqrt(s0 ** 2 + j ** 2), np.sqrt(s1 ** 2 + j ** 2)
        c01 = d.corr * s0 * s1
        det = s0 ** 2 * s1 ** 2 - c01 ** 2
        out.append(((s1 ** 2 * r[:, 0] ** 2 + s0 ** 2 * r[:, 1] ** 2 - 2 * c01 * r[:, 0] * r[:, 1]) / det).sum())
    out = np.array(out)
    n_data = 2 * len(d.pidx)
    k = m.ndim + 2 * len(d.offset_sources)          # sampled parameters + the profiled offsets
    return {"chi2": q(out), "n_data": n_data, "k": k, "chi2nu": q(out / (n_data - k))}
