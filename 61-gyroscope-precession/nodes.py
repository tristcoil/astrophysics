"""Weekly orbit fits and the chain of node "jumps" that measures frame dragging.

For every week k of an analysis centre's orbits of a satellite we

1. fit our orbit model (frame dragging OFF, ``mu_LT = 0``) to the 7-day orbit
   (``odfit.fit_arc``);
2. propagate the fitted orbit half a day past the end of the week, and half a
   day back before its start.

Where week k's forward propagation and week k+1's backward propagation overlap
(one day centred on the boundary between the weeks), both describe the same
satellite at the same times. The difference of their orbital planes is the
"node jump": how far the real orbital plane has turned in one week beyond what
our model predicts. Each fit lands on the true node at the middle of its week,
so the jump equals (true node rate - model node rate) x 7 days, whatever the
evaluation time. Summed over the years, the jumps give the node *residual*: the
part of the orbital plane's motion our Newtonian (plus de Sitter) model leaves
unexplained. Frame dragging should make it grow by ~31 mas per year for each
satellite.

Results are cached per arc in ``data/cache/arcs/<sat>/<ac>/<week>.npz`` so a
run can be interrupted and resumed.
"""

from __future__ import annotations

import datetime as dt
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np

import config as C
import frames as F
import odfit

WINDOW_HALF_S = 0.5 * C.SECONDS_PER_DAY      # half-width of the overlap window
OUT_STEP_S = 120.0                            # sampling of the stored windows
ARC_DIR = C.CACHE_DIR / "arcs"


def week_date(week: str) -> dt.date:
    return dt.date(2000 + int(week[:2]), int(week[2:4]), int(week[4:]))


MIN_SP3_BYTES = 20_000      # a week of 2-min records is ~0.24 MB; a few EDC files hold only a header


def arc_files(sat: str, ac: str) -> list[Path]:
    """Weekly SP3 files within the analysis span, one (latest version) per week.

    Truncated files (e.g. jcet.orb.lageos1.250405.v90, 98 bytes: a header line only) are
    skipped, so the week before propagates across the gap instead of breaking the chain.
    """
    out = []
    for p in F.orbit_files(sat, ac):
        d = week_date(p.name.split(".")[3])
        if C.ANALYSIS_START <= d <= C.ANALYSIS_END + dt.timedelta(days=7) \
                and p.stat().st_size >= MIN_SP3_BYTES:
            out.append(p)
    return out


def sp3_start_mjd_tt(path: Path) -> float:
    """First epoch of an SP3 file (from its header line), MJD TT.

    The EDC folder name is the *last* day of the week, so read the header rather
    than guess (e.g. folder 230121 holds 2023-01-15 00:00 to 2023-01-22 00:00).
    """
    import gzip
    with gzip.open(path, "rt") as fh:
        first = fh.readline()
    y, mo, d, h, mi = (int(v) for v in first[3:31].split()[:5])
    sec = float(first[3:31].split()[5])
    mjd_utc = F.date_to_mjd(dt.datetime(y, mo, d, h, mi)) + sec / C.SECONDS_PER_DAY
    return float(F.utc_to_tt_mjd(mjd_utc))


def cache_path(sat: str, ac: str, week: str, tag: str = "") -> Path:
    return ARC_DIR / f"{sat}{tag}" / ac / f"{week}.npz"


def process_arc(path: Path, next_start_mjd: float | None = None, model_kw: dict | None = None,
                tag: str = "", force: bool = False, fit_kw: dict | None = None) -> Path:
    """Fit one week and store the fit plus the two overlap windows."""
    orbit = F.read_sp3(path)
    out = cache_path(orbit.sat, orbit.ac, orbit.week, tag)
    if out.exists() and not force:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    mjd0 = float(orbit.mjd_tt[0])
    span_days = (orbit.mjd_tt[-1] - mjd0)
    fwd_days = span_days if next_start_mjd is None else max(span_days, next_start_mjd - mjd0)
    kw = dict(model_kw or {})
    kw.setdefault("t_before", WINDOW_HALF_S + 3600.0)
    t0 = time.time()
    fit, model = odfit.fit_arc(orbit, extra_days=fwd_days - span_days + 0.6, model_kw=kw,
                               **(fit_kw or {}))
    stride = int(round(OUT_STEP_S / C.INTEGRATION_STEP_S))
    # forward through the week and half a day past the next week's start
    t_fwd_end = fwd_days * C.SECONDS_PER_DAY + WINDOW_HALF_S
    Tf, Yf = model.run(fit.y0, fit.prm, 0.0, t_fwd_end, stride=stride)
    # backward half a day from the start
    Tb, Yb = model.run(fit.y0, fit.prm, 0.0, -WINDOW_HALF_S, stride=stride)
    # the forward window around the next start, and the backward window around this start
    lo = fwd_days * C.SECONDS_PER_DAY - WINDOW_HALF_S
    sel_f = Tf >= lo - 1e-6
    Tb, Yb = Tb[::-1], Yb[::-1]
    sel_b_extra = (Tf >= 0) & (Tf <= WINDOW_HALF_S + 1e-6)
    tb = np.concatenate([Tb[:-1], Tf[sel_b_extra]])
    yb = np.concatenate([Yb[:-1, :6], Yf[sel_b_extra, :6]])
    # weekly mean osculating elements (for the bulk drift plots)
    inweek = (Tf >= 0) & (Tf <= span_days * C.SECONDS_PER_DAY + 1e-6)
    np.savez_compressed(
        out, sat=orbit.sat, ac=orbit.ac, week=orbit.week, mjd0=mjd0, span_days=span_days,
        fwd_days=fwd_days, y0=fit.y0, prm=fit.prm, cov=fit.cov, rms_rtn=fit.rms_rtn,
        n_obs=fit.n_obs, n_iter=fit.n_iter, converged=fit.converged,
        t_fwd=Tf[sel_f], y_fwd=Yf[sel_f, :6], t_bwd=tb, y_bwd=yb,
        t_week=Tf[inweek][::5], y_week=Yf[inweek][::5, :6],
        model_kw=str(kw), runtime=time.time() - t0)
    return out


def _worker(job):
    path, next_start, model_kw, tag, force, fit_kw = job
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:
        pass
    try:
        return str(process_arc(path, next_start, model_kw, tag, force, fit_kw)), None
    except Exception as exc:                    # keep going; report at the end
        return str(path), f"{type(exc).__name__}: {exc}"


def run_all(sats=("lageos1", "lares2"), acs=(C.PRIMARY_AC,), model_kw=None, tag="",
            workers: int = 8, force: bool = False, fit_kw: dict | None = None) -> list:
    jobs = []
    for sat in sats:
        for ac in acs:
            files = arc_files(sat, ac)
            weeks = [p.name.split(".")[3] for p in files]
            mjd_starts = [sp3_start_mjd_tt(p) for p in files]
            for i, p in enumerate(files):
                nxt = mjd_starts[i + 1] if i + 1 < len(files) else None
                if not force and cache_path(sat, ac, weeks[i], tag).exists():
                    continue
                jobs.append((p, nxt, model_kw, tag, force, fit_kw))
    print(f"  {len(jobs)} arcs to fit with {workers} workers")
    if not jobs:
        return []
    errors = []
    ctx = mp.get_context("fork")
    t0 = time.time()
    with ctx.Pool(workers) as pool:
        for i, (path, err) in enumerate(pool.imap_unordered(_worker, jobs), 1):
            if err:
                errors.append((path, err))
                print(f"    ERROR {Path(path).name}: {err}")
            if i % 25 == 0 or i == len(jobs):
                print(f"    {i}/{len(jobs)} arcs, {time.time() - t0:.0f} s")
    return errors


# ---------------------------------------------------------------------------
# Jumps and residuals
# ---------------------------------------------------------------------------

def load_arcs(sat: str, ac: str, tag: str = "") -> list[dict]:
    arcs = []
    for p in sorted((ARC_DIR / f"{sat}{tag}" / ac).glob("*.npz")):
        with np.load(p, allow_pickle=False) as z:
            arcs.append({k: z[k] for k in z.files})
    arcs.sort(key=lambda a: float(a["mjd0"]))
    return arcs


def _wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi


def jumps(sat: str, ac: str, tag: str = "", max_gap_days: float = 7.01) -> dict:
    """Node and inclination jumps between consecutive weekly fits.

    Returns arrays over boundaries: mjd (boundary, TT), dnode, dinc (rad), the
    scatter of the node difference within the window (rad), and the gap length.
    """
    arcs = load_arcs(sat, ac, tag)
    rows = []
    for a, b in zip(arcs[:-1], arcs[1:]):
        gap = float(b["mjd0"]) - float(a["mjd0"])
        if gap > float(a["fwd_days"]) + 1e-3 or gap > max_gap_days:
            continue
        ta = float(a["mjd0"]) * C.SECONDS_PER_DAY + a["t_fwd"]
        tb = float(b["mjd0"]) * C.SECONDS_PER_DAY + b["t_bwd"]
        common, ia, ib = np.intersect1d(np.round(ta, 3), np.round(tb, 3), return_indices=True)
        if common.size < 100:
            continue
        ea = F.kepler_elements(a["y_fwd"][ia, :3], a["y_fwd"][ia, 3:])
        eb = F.kepler_elements(b["y_bwd"][ib, :3], b["y_bwd"][ib, 3:])
        dn = _wrap(eb["node"] - ea["node"])
        di = eb["i"] - ea["i"]
        rows.append((float(b["mjd0"]), np.mean(dn), np.mean(di), np.std(dn), gap,
                     float(np.mean(a["rms_rtn"][2])), float(np.mean(b["rms_rtn"][2]))))
    r = np.array(rows)
    return {"mjd": r[:, 0], "dnode": r[:, 1], "dinc": r[:, 2], "dnode_scatter": r[:, 3],
            "gap": r[:, 4], "rms_n_a": r[:, 5], "rms_n_b": r[:, 6]}


def weekly_elements(sat: str, ac: str, tag: str = "") -> dict:
    """Mean osculating elements per week (GCRS and true-equator node) for the drift plots."""
    import erfa
    out = {k: [] for k in ("mjd", "node_gcrs", "inc_gcrs", "node_cip", "inc_cip", "a", "e")}
    for a in load_arcs(sat, ac, tag):
        t = float(a["mjd0"]) + a["t_week"] / C.SECONDS_PER_DAY
        y = a["y_week"]
        el = F.kepler_elements(y[:, :3], y[:, 3:])
        # rotate into the CIRS (pole = CIP, origin = CIO): r_cirs = Q^T r_gcrs
        x, yy = erfa.xy06(2_400_000.5, t)
        s = erfa.s06(2_400_000.5, t, x, yy)
        Q = np.transpose(erfa.c2ixys(x, yy, s), (0, 2, 1))
        rc = np.einsum("nji,nj->ni", Q, y[:, :3])
        vc = np.einsum("nji,nj->ni", Q, y[:, 3:])
        ec = F.kepler_elements(rc, vc)
        out["mjd"].append(np.mean(t))
        out["node_gcrs"].append(np.angle(np.mean(np.exp(1j * el["node"]))))
        out["inc_gcrs"].append(np.mean(el["i"]))
        out["node_cip"].append(np.angle(np.mean(np.exp(1j * ec["node"]))))
        out["inc_cip"].append(np.mean(ec["i"]))
        out["a"].append(np.mean(el["a"]))
        out["e"].append(np.mean(el["e"]))
    return {k: np.array(v) for k, v in out.items()}


if __name__ == "__main__":
    import sys
    acs = tuple(sys.argv[1:]) or (C.PRIMARY_AC,)
    run_all(acs=acs, workers=int(os.environ.get("N4_WORKERS", "8")))
