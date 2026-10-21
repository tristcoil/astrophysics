"""Stage 1: measure every time delay three ways and cache the results.

    python measure_delays.py            # all six lenses (~1 h on 4 cores)
    python measure_delays.py HE0435     # one lens
    python measure_delays.py --workers 2 --mocks 100

For each lens:

1. **Blind start.** The high-passed ICCF over +-200 days (+-60 for PG 1115's
   single season) gives every delay a starting value without looking at the
   published numbers. An image that correlates poorly with the reference
   image (r_max < 0.5) is located through the image it correlates with best.
2. **Joint spline** (``delays.py``), refined from the blind start, for a grid
   of knot settings, each with 100 simulated light curves for the errors.
3. **ICCF** with flux-randomization / random-subset-selection errors
   (``ccf.py``), in a +-40-day window around the blind start.
4. **PyCS3** (Millon et al. 2020c), COSMOGRAIL's own curve-shifting toolbox:
   its free-knot spline and regression-difference estimators.

Each task writes ``results/delays/<lens>_<task>.json``; finished tasks are
skipped on re-runs, so the stage can be interrupted and resumed.
"""

from __future__ import annotations

import argparse
import json
import logging
import multiprocessing as mp
import os
import sys
import time
import warnings

import numpy as np

import config as C

OUT = C.RESULTS / "delays"

#: Reference image (the first to arrive, or the brightest) per lens.
REF = {"RXJ1131": "A", "HE0435": "A", "PG1115": "A", "WFI2033": "B", "J1206": "A", "DES0408": "A"}

#: Joint-spline settings. The intrinsic knot step is chosen per lens from the
#: data (``knot_selection``: lowest BIC on this ladder at the shorter ML knot
#: step), then the best step and its two neighbours on the ladder are each run
#: with two microlensing knot steps: six settings per lens.
LADDER_MULTI = [10, 15, 20, 25, 35, 45, 55]
LADDER_SINGLE = [5, 8, 10, 15, 25]
ML_MULTI = [150, 300]
ML_SINGLE = [100, 200]
SINGLE_SEASON = {"PG1115", "DES0408"}
#: Mock true delays are drawn within +-SPREAD days of the measurement (and the
#: blind re-fit searches +-SPREAD+3 days); wider where the errors are large.
SPREAD = {"DES0408": 12.0, "J1206": 8.0}

#: PyCS3 settings: spline (knotstep, ML knotstep) and regdiff GP kernel scale.
PYCS_SPL = [(25, 150), (35, 300)]
PYCS_SPL_SINGLE = [(15, 40), (25, 55)]  # one 6-9 month season: PyCS3 needs >= 3 ML knot steps inside it


def _limit_threads():
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(1)
    except ImportError:
        pass
    warnings.filterwarnings("ignore")
    logging.disable(logging.WARNING)


def _save(name: str, obj: dict):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(obj, indent=1))


def _load(name: str):
    p = OUT / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


# ---------------------------------------------------------------------------
# Tasks
# ---------------------------------------------------------------------------


def blind_start(key: str) -> dict:
    """Starting delays from the wide ICCF search (cached)."""
    import ccf
    import lightcurves as L

    cached = _load(f"{key}_start")
    if cached:
        return cached
    lc = L.read_lens(key)
    ref = REF[key]
    win = 60.0 if key == "PG1115" else 200.0
    pairs = {}
    for a in lc.images:
        for b in lc.images:
            if a < b:
                r = ccf.ccf_delay(lc, a, b, -win, win, step=0.5, n_mc=0, min_overlap=0.5)
                pairs[f"{a}{b}"] = dict(lag=r.tau_cent if r.ok else r.tau_peak, rmax=r.r_max,
                                        lags=r.lags.tolist(), r=r.r.tolist())
    start, via = {}, {}
    for im in lc.images:
        if im == ref:
            continue
        k = f"{min(ref, im)}{max(ref, im)}"
        sgn = 1.0 if ref < im else -1.0
        start[im], via[im] = sgn * pairs[k]["lag"], ref
    # weakly correlated images: go through the best-correlated other image
    for im in list(start):
        k = f"{min(ref, im)}{max(ref, im)}"
        if pairs[k]["rmax"] >= 0.5:
            continue
        best = None
        for other in start:
            if other == im:
                continue
            ko = f"{min(other, im)}{max(other, im)}"
            ko_ref = f"{min(ref, other)}{max(ref, other)}"
            if pairs[ko_ref]["rmax"] < 0.5:
                continue
            if best is None or pairs[ko]["rmax"] > best[1]:
                sgn = 1.0 if other < im else -1.0
                best = (start[other] + sgn * pairs[ko]["lag"], pairs[ko]["rmax"], other)
        if best is not None and best[1] > pairs[k]["rmax"]:
            start[im], via[im] = best[0], best[2]
    out = dict(lens=key, ref=ref, start=start, via=via, pairs=pairs)
    _save(f"{key}_start", out)
    return out



def knot_selection(key: str) -> dict:
    """Chi^2 and BIC of the joint spline at its optimum for each knot step on the ladder.

    BIC = chi^2 + k ln N with k the number of intrinsic-curve coefficients
    (span / knot step) and N the number of data points; errors are rescaled
    once at a 20-day (10-day for one season) fit so all chi^2 are comparable.
    """
    import delays as D
    import lightcurves as L

    cached = _load(f"{key}_knots")
    if cached:
        return cached
    single = key in SINGLE_SEASON
    ladder = LADDER_SINGLE if single else LADDER_MULTI
    ml = (ML_SINGLE if single else ML_MULTI)[0]
    st = blind_start(key)
    lc = L.read_lens(key)
    ref = st["ref"]
    s0 = D.SplineSettings(knot_int=10.0 if single else 20.0, knot_ml=float(ml))
    t0 = D.JointSpline(lc, ref, s0).optimize(st["start"], width=12.0 if single else 20.0)
    scale = D.rescale_errors(lc, ref, t0, s0)
    n = sum(len(v) for v in lc.mag.values())
    span = lc.t.max() - lc.t.min()
    rows = []
    for k in ladder:
        m = D.JointSpline(lc, ref, D.SplineSettings(knot_int=float(k), knot_ml=float(ml)), err_scale=scale)
        best = m.optimize(t0, width=8.0)
        f = m.fit(best)
        rows.append(dict(knot_int=k, chi2=f.chi2_data, n_coef=span / k, bic=f.chi2_data + span / k * np.log(n), taus=best))
    ib = int(np.argmin([r["bic"] for r in rows]))
    chosen = [ladder[i] for i in range(max(ib - 1, 0), min(ib + 2, len(ladder)))]
    out = dict(lens=key, knot_ml=ml, n=n, rows=rows, best=ladder[ib], chosen=chosen)
    _save(f"{key}_knots", out)
    return out


def settings_for(key: str) -> list[tuple[int, int]]:
    ks = knot_selection(key)
    ml = ML_SINGLE if key in SINGLE_SEASON else ML_MULTI
    return [(k, m) for m in ml for k in ks["chosen"]]


def task_knots(key):
    _limit_threads()
    knot_selection(key)
    return f"{key}_knots"


def task_spline(args):
    key, knot_int, knot_ml, n_mock = args
    _limit_threads()
    import delays as D
    import lightcurves as L

    name = f"{key}_spline_{knot_int}_{knot_ml}"
    if _load(name):
        return name
    t0 = time.time()
    st = blind_start(key)
    lc = L.read_lens(key)
    s = D.SplineSettings(knot_int=float(knot_int), knot_ml=float(knot_ml))
    m = D.measure(lc, st["ref"], st["start"], s, n_mock=n_mock, width=12.0 if key in SINGLE_SEASON else 20.0,
                  spread=SPREAD.get(key, 5.0),
                  seed=knot_int * 1000 + knot_ml)
    _save(name, dict(lens=key, method="joint spline", knot_int=knot_int, knot_ml=knot_ml, ref=m.ref,
                     others=m.others, taus=m.taus, err=m.err, rand=m.rand, bias=m.bias, cov=m.cov.tolist(),
                     chi2=m.chi2, ndata=m.ndata, n_mock=n_mock,
                     err_scale={f"{a}|{b}": v for (a, b), v in m.err_scale.items()},
                     mock_dev={k: v.tolist() for k, v in (m.mock_dev or {}).items()},
                     seconds=time.time() - t0))
    return name


def task_ccf(args):
    key, n_mc = args
    _limit_threads()
    import ccf
    import lightcurves as L

    name = f"{key}_ccf"
    if _load(name):
        return name
    st = blind_start(key)
    lc = L.read_lens(key)
    ref = st["ref"]
    res = {}
    for im, g in st["start"].items():
        r = ccf.ccf_delay(lc, ref, im, g - 40.0, g + 40.0, step=0.25, n_mc=n_mc, min_overlap=0.5)
        res[im] = dict(peak=r.tau_peak, cent=r.tau_cent, rmax=r.r_max, med=r.cent[0], minus=r.cent[1],
                       plus=r.cent[2], success=r.success, lags=r.lags.tolist(), r=r.r.tolist())
    _save(name, dict(lens=key, method="ICCF", ref=ref, delays=res))
    return name


def task_pycs(args):
    key, tech, p1, p2 = args
    _limit_threads()
    import lightcurves as L
    import pycs3.gen.lc_func as lcf
    import pycs3.gen.splml as splml
    import pycs3.regdiff.multiopt as rmo
    import pycs3.spl.topopt as topopt

    name = f"{key}_pycs_{tech}_{p1}_{p2}"
    if _load(name):
        return name
    t0 = time.time()
    st = blind_start(key)
    lc = L.read_lens(key)
    ref = st["ref"]
    # Telescopes sit on different zero points for some lenses (WFI 2033's four
    # instruments); PyCS3 wants one curve per image, so each telescope is put
    # on the first one's zero point using the joint spline's fitted offsets.
    offs = _load(f"{key}_offsets") or {}
    lcs = []
    for im in lc.images:
        mag = lc.mag[im].copy()
        for tel, o in offs.get(im, {}).items():
            mag[lc.tel == tel] -= o
        lo = lcf.factory(lc.t.copy(), mag, lc.err[im].copy(), telescopename="merged", object=im)
        lo.shifttime(-st["start"].get(im, 0.0))  # PyCS shifts a lagging curve back in time
        lcs.append(lo)
    if tech == "spl":
        for lo in lcs:
            if lo.object != ref:
                splml.addtolc(lo, knotstep=float(p2))
        topopt.opt_fine(lcs, nit=5, knotstep=float(p1), verbose=False)
    else:
        rmo.opt_ts(lcs, pd=2, covkernel="matern", pow=1.5, amp=0.5, scale=float(p1), errscale=1.0, verbose=False)
    ts = {lo.object: lo.timeshift for lo in lcs}
    taus = {im: float(ts[ref] - ts[im]) for im in lc.images if im != ref}
    _save(name, dict(lens=key, method=f"PyCS3 {tech}", p1=p1, p2=p2, ref=ref, taus=taus, seconds=time.time() - t0))
    return name


def telescope_offsets(key: str):
    """Per-image telescope offsets from a joint-spline fit at the blind start (for PyCS3)."""
    import delays as D
    import lightcurves as L

    if _load(f"{key}_offsets") is not None:
        return
    st = blind_start(key)
    lc = L.read_lens(key)
    m = D.JointSpline(lc, st["ref"], D.SplineSettings())
    taus = m.optimize(st["start"], width=12.0 if key in SINGLE_SEASON else 20.0)
    fit = m.fit(taus)
    out = {}
    for im in lc.images:
        o = m.offsets(fit.coef, im)
        out[im] = {tel: float(o[lc.tel == tel][0]) for tel in lc.telescopes[1:]}
    _save(f"{key}_offsets", out)


# ---------------------------------------------------------------------------
# Combination across settings
# ---------------------------------------------------------------------------


def combine(key: str) -> dict:
    """Joint-spline result marginalized over the settings grid.

    Central value: mean over settings. Covariance: mean of the per-setting
    mock covariances (random + bias) plus the scatter between settings
    (the systematic from the estimator's choices).
    """
    runs = [_load(f"{key}_spline_{a}_{b}") for a, b in settings_for(key)]
    runs = [r for r in runs if r]
    others = runs[0]["others"]
    T = np.array([[r["taus"][im] for im in others] for r in runs])
    covs = []
    for r in runs:
        M = np.array([r["mock_dev"][im] for im in others])
        b = M.mean(axis=1)
        covs.append(np.atleast_2d(np.cov(M)) + np.outer(b, b))
    cov_stat = np.mean(covs, axis=0)
    cov_sys = np.atleast_2d(np.cov(T.T)) if len(runs) > 1 else np.zeros_like(cov_stat)
    cov = cov_stat + cov_sys
    mean = T.mean(axis=0)
    out = dict(lens=key, ref=runs[0]["ref"], others=others, taus=dict(zip(others, mean.tolist())),
               err=dict(zip(others, np.sqrt(np.diag(cov)).tolist())),
               err_stat=dict(zip(others, np.sqrt(np.diag(cov_stat)).tolist())),
               err_sys=dict(zip(others, np.sqrt(np.diag(cov_sys)).tolist())),
               cov=cov.tolist(), n_settings=len(runs),
               settings=[dict(knot_int=r["knot_int"], knot_ml=r["knot_ml"], taus=r["taus"], err=r["err"]) for r in runs])
    _save(f"{key}_combined", out)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("lenses", nargs="*", default=C.MEASURED)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--mocks", type=int, default=100)
    ap.add_argument("--ccf-mc", type=int, default=500)
    ap.add_argument("--no-pycs", action="store_true")
    args = ap.parse_args()
    _limit_threads()
    for key in args.lenses:
        st = blind_start(key)
        print(f"{key}: blind start {', '.join(f'{im} {v:+.1f} (via {st['via'][im]})' for im, v in st['start'].items())}", flush=True)
        if not args.no_pycs:
            telescope_offsets(key)
    fn = {"spline": task_spline, "ccf": task_ccf, "pycs": task_pycs, "knots": task_knots}
    ctx = mp.get_context("fork")
    t0 = time.time()
    # Phase 1: choose knot steps (BIC) for every lens.
    with ctx.Pool(args.workers, maxtasksperchild=1) as pool:
        for name in pool.imap_unordered(task_knots, args.lenses):
            print(f"  done {name}  ({(time.time() - t0) / 60:.1f} min)", flush=True)
    for key in args.lenses:
        ks = knot_selection(key)
        print(f"{key}: BIC-best intrinsic knot step {ks['best']} d -> settings {settings_for(key)}", flush=True)
    # Phase 2: everything else.
    tasks = []
    for key in args.lenses:
        tasks += [("spline", (key, a, b, args.mocks)) for a, b in settings_for(key)]
        tasks.append(("ccf", (key, args.ccf_mc)))
        if not args.no_pycs:
            spl = PYCS_SPL_SINGLE if key in SINGLE_SEASON else PYCS_SPL
            tasks += [("pycs", (key, "spl", a, b)) for a, b in spl]
            tasks.append(("pycs", (key, "regdiff", 200, 0)))
    with ctx.Pool(args.workers, maxtasksperchild=1) as pool:
        futs = [pool.apply_async(fn[kind], (a,)) for kind, a in tasks]
        for f in futs:
            try:
                print(f"  done {f.get()}  ({(time.time() - t0) / 60:.1f} min)", flush=True)
            except Exception as exc:  # keep going; a failed task is reported and can be re-run
                print(f"  FAILED: {exc!r}", flush=True)
    for key in args.lenses:
        c = combine(key)
        print(f"{key} (vs {c['ref']}): " + ", ".join(f"{im} {c['taus'][im]:+.2f} +- {c['err'][im]:.2f}" for im in c["others"]))


if __name__ == "__main__":
    sys.exit(main())
