"""Our own lags for SDSS-RM quasars, from the raw light curves.

For every quasar whose lag Shen et al. (2024) detected in H-beta, Mg II, C IV
or H-alpha, and for every quasar whose H-beta lag they attempted, we run our
ICCF with FR/RSS errors and the alias weighting of Grier et al. (2017, 2019)
over the same lag search window. Settings follow pyPETAL, the pipeline behind
the release: a 2-day lag grid, threshold 0.8, r_max >= 0.2 per realization.
Differences from the release: 500 instead of 3000 realizations per lag, and
no 5-sigma outlier clipping of the continuum against a damped-random-walk fit.

A lag counts as detected (our ICCF-only version of Shen et al.'s rule) when
it is positive at >= 2 sigma, the ICCF reaches r > 0.4 within +/-1 sigma of
it, and the alias weighting rejects less than half of the realizations.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import config as C
import iccf
import lightcurves as LC
import rl

SHORT = {"HB": "hb", "MG2": "mg2", "C4": "c4", "HA": "ha"}
CACHE = C.CACHE_DIR / "sdssrm_iccf.parquet"


def _init_worker():
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(1)
    except ImportError:
        pass
    os.environ["OMP_NUM_THREADS"] = "1"


def measure(rmid: int, line: str, search: tuple, n_mc: int = C.SDSSRM_N_MC, seed: int = 0) -> dict:
    cont = LC.sdssrm(rmid, "cont")
    ln = LC.sdssrm(rmid, SHORT[line])
    base = dict(RMID=rmid, line=line)
    if cont is None or ln is None or len(ln) < 10:
        return dict(base, ok=False)
    lo, hi = search
    lags = iccf.lag_grid(np.floor(lo / C.SDSSRM_LAG_STEP) * C.SDSSRM_LAG_STEP,
                         np.ceil(hi / C.SDSSRM_LAG_STEP) * C.SDSSRM_LAG_STEP, C.SDSSRM_LAG_STEP)
    res = iccf.fr_rss(cont.t, cont.y, cont.e, ln.t, ln.y, ln.e, lags, n_mc=n_mc, seed=seed)
    w = iccf.weighted_lag(res, cont.t, cont.y, ln.t, k=C.ALIAS_K, gap=C.ALIAS_GAP_DAYS, smooth=C.ALIAS_SMOOTH_DAYS)
    unw = res.cent
    out = dict(base, ok=True, n_cont=len(cont), n_line=len(ln), search_lo=lo, search_hi=hi,
               lag=w["lag"], lag_lo=w["err"][0], lag_hi=w["err"][1], frac_rej=w["frac_rej"],
               rmax_near=w.get("rmax_near", np.nan), r_max=res.observed.r_max, success=res.success,
               lag_unweighted=unw[0], lag_unweighted_lo=unw[1], lag_unweighted_hi=unw[2],
               bound_lo=w["bounds"][0], bound_hi=w["bounds"][1])
    out["detected"] = bool(np.isfinite(out["lag"]) and out["lag"] > 0 and out["lag"] >= 2 * out["lag_lo"]
                           and out["rmax_near"] > 0.4 and out["frac_rej"] < 0.5)
    return out


def _job(args):
    return measure(*args)


def targets() -> list[tuple]:
    """(RMID, line, search window) for every lag we measure."""
    t = rl.sdssrm_table()
    jobs = []
    for line in ("HB", "MG2", "C4", "HA"):
        sel = t[f"{line}_lag"] > 0
        if line == "HB":
            sel |= t["HB_done"] == 1
        for _, r in t[sel].iterrows():
            lo, hi = r.search_lo, r.search_hi
            pub = r[f"{line}_lag"]
            # Shen et al. widened the window for a few lines whose lag hit 90% of it.
            if np.isfinite(pub) and pub > 0.9 * hi:
                lo, hi = -2500.0, 4000.0
            jobs.append((int(r.RMID), line, (float(lo), float(hi))))
    return jobs


def run_all(workers: int = 6, n_mc: int = C.SDSSRM_N_MC, force: bool = False) -> pd.DataFrame:
    """Measure every target (cached: ~20-40 minutes on 6 cores the first time)."""
    if CACHE.exists() and not force:
        return pd.read_parquet(CACHE)
    jobs = [(rmid, line, search, n_mc, 1000 + i) for i, (rmid, line, search) in enumerate(targets())]
    print(f"  SDSS-RM: measuring {len(jobs)} lags with our ICCF on {workers} cores ...")
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as ex:
        rows = list(ex.map(_job, jobs, chunksize=2))
    df = pd.DataFrame(rows)
    df.to_parquet(CACHE)
    return df


def compare(ours: pd.DataFrame) -> dict:
    """Our ICCF lags against Shen et al.'s ICCF and PyROA lags (observed frame)."""
    t = rl.sdssrm_table().set_index("RMID")
    out = {}
    for line in ("HB", "MG2", "C4", "HA"):
        o = ours[(ours.line == line) & ours.ok].set_index("RMID")
        pub_det = t.loc[o.index, f"{line}_lag"] > 0
        both = o[pub_det.values]
        iccf_pub = t.loc[both.index, f"{line}_iccf"]
        pyroa = t.loc[both.index, f"{line}_lag"]
        ratio = (both.lag / iccf_pub).replace([np.inf, -np.inf], np.nan)
        err = np.sqrt((0.5 * (both.lag_lo + both.lag_hi)) ** 2
                      + (0.5 * (t.loc[both.index, f"{line}_iccf_lo"] + t.loc[both.index, f"{line}_iccf_hi"])) ** 2)
        within = np.abs(both.lag - iccf_pub) < err
        out[line] = dict(
            n_measured=int(len(o)), n_pub_det=int(pub_det.sum()),
            n_ours_det=int(o.detected.sum()), n_both_det=int((o.detected & pub_det.values).sum()),
            median_ratio_iccf=float(np.nanmedian(ratio)), frac_within_1sig_iccf=float(np.mean(within)),
            median_ratio_pyroa=float(np.nanmedian(both.lag / pyroa)),
            median_abs_log_diff=float(np.nanmedian(np.abs(np.log10(both.lag.clip(lower=0.1) / iccf_pub.clip(lower=0.1))))),
        )
    return out


def false_positive_test(ours: pd.DataFrame) -> dict:
    """Shen et al.'s asymmetry test on every attempted H-beta lag: with a
    symmetric search window, noise gives as many negative as positive lags, so
    the negative/positive ratio at r_max > 0.4 estimates the false-positive rate."""
    h = ours[(ours.line == "HB") & ours.ok & np.isfinite(ours.lag)]
    hi_r = h.rmax_near > 0.4
    pos = (h.lag > 0) & hi_r
    neg = (h.lag < 0) & hi_r
    lo_pos = (h.lag > 0) & ~hi_r
    lo_neg = (h.lag < 0) & ~hi_r
    return dict(n=int(len(h)), pos_hi=int(pos.sum()), neg_hi=int(neg.sum()), pos_lo=int(lo_pos.sum()),
                neg_lo=int(lo_neg.sum()),
                fpr_hi=float(neg.sum() / max(pos.sum(), 1)), fpr_lo=float(lo_neg.sum() / max(lo_pos.sum(), 1)))
