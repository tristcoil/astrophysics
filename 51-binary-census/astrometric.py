"""Would Gaia DR3 have found this binary astrometrically? (gaiamock)

gaiamock (El-Badry et al. 2024) predicts the along-scan positions Gaia would
have measured for a binary at a given place on the sky -- with the real DR3
scan times and angles there, the published per-transit uncertainties for its
magnitude, and the light centre's orbit -- and then runs the cascade DR3's
non-single-star pipeline ran (Halbwachs et al. 2023):

1. a 5-parameter single-star fit; if its RUWE < 1.4 the source stays "single";
2. otherwise 9- and 7-parameter "acceleration" models (curved proper motion,
   for orbits longer than the 34-month baseline), accepted if significant
   and well fitted;
3. otherwise a full 12-parameter orbit, which DR3 published only if it passed
   its acceptance cuts (below).

We call it for every mock binary and record which solution DR3 would have
published. The binary's noise realisation uses numpy's global random state
(gaiamock's convention), so each task seeds it explicitly.
"""

from __future__ import annotations

import multiprocessing as mp
import os
import sys
import time

import numpy as np

import config as C

_GM = None
_CF = None


def _load():
    global _GM, _CF
    if _GM is None:
        sys.path.insert(0, str(C.GAIAMOCK_DIR))
        import gaiamock  # noqa: E402  (external code, pinned commit)
        _GM = gaiamock
        _CF = gaiamock.read_in_C_functions()
    return _GM, _CF


def dr3_orbit_accepted(period, a0_over_err, plx_over_err, ecc_err, f2) -> bool:
    """DR3's acceptance cuts for astrometric orbits (Halbwachs et al. 2023), as
    implemented in gaiamock: goodness of fit F2 < 25, significance
    a0/sigma_a0 > max(5, 158/sqrt(P)), parallax/sigma > 20000/P, and
    sigma_e < 0.079 ln P - 0.244 (P in days)."""
    return bool((f2 < 25) and (a0_over_err > 5) and (a0_over_err > 158 / np.sqrt(period))
                and (plx_over_err > 20000 / period) and (ecc_err < 0.079 * np.log(period) - 0.244))


def run_one(task) -> tuple:
    """One mock binary through the cascade.

    Returns (outcome, ruwe, x, y): outcome 0 too few visibility periods,
    1 single (RUWE < 1.4), 7/9 acceleration solution (x = significance s,
    y = F2), 12 accepted orbit / -12 orbit fitted but rejected by the DR3 cuts
    (x = a0 in mas, y = fitted period).
    """
    gm, cf = _load()
    (seed, ra, dec, plx, pmra, pmdec, m1, m2, period, t_peri, ecc, node, inc_deg, omega, gmag, fratio) = task
    np.random.seed(seed)
    try:
        res = gm.run_full_astrometric_cascade(
            ra=ra, dec=dec, parallax=plx, pmra=pmra, pmdec=pmdec, m1=m1, m2=m2, period=period,
            Tp=t_peri, ecc=ecc, omega=node, inc_deg=inc_deg, w=omega, phot_g_mean_mag=gmag, f=fratio,
            data_release="dr3", c_funcs=cf)
    except Exception:  # a failed fit counts as no published orbit
        return (-99, np.nan, np.nan, np.nan)
    first = res[0]
    if first == 0 and all(r == 0 for r in res):
        return (0, np.nan, np.nan, np.nan)
    if first == -1:
        return (1, float(res[1]), np.nan, np.nan)
    # accelerations: keep the significance s and goodness of fit F2, because
    # DR3's *final* selection (Halbwachs et al. 2023, Table 1) kept only
    # s > 20 (and F2 < 22 for the 7-parameter model); main processing used s > 12
    if first == -9:
        return (9, float(res[12]), float(res[1]), float(res[13]))
    if first == -7:
        return (7, float(res[8]), float(res[1]), float(res[9]))
    (plx_f, sig_plx, _a, _sa, _b, _sb, _f, _sf, _g, _sg, p_fit, _sp, _phi, _sphi, _e, sig_e, _inc,
     a0, sig_a0, _nvis, _nobs, f2, ruwe) = res
    ok = dr3_orbit_accepted(p_fit, a0 / sig_a0, plx_f / sig_plx, sig_e, f2)
    return (12 if ok else -12, float(ruwe), float(a0), float(p_fit))


def _init_worker():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:
        pass
    _load()


def run_many(tasks: list, workers: int = 6, chunksize: int = 8, progress_every: int = 2000) -> np.ndarray:
    """Run tasks in a fork pool (Python 3.14 defaults to forkserver)."""
    out = np.zeros((len(tasks), 4))
    ctx = mp.get_context("fork")
    t0 = time.time()
    with ctx.Pool(workers, initializer=_init_worker) as pool:
        for i, r in enumerate(pool.imap(run_one, tasks, chunksize=chunksize)):
            out[i] = r
            if progress_every and (i + 1) % progress_every == 0:
                rate = (i + 1) / (time.time() - t0)
                print(f"    gaiamock {i + 1:,}/{len(tasks):,} ({rate:.1f}/s)", flush=True)
    return out


os.environ.setdefault("OMP_NUM_THREADS", "1")


def dr3_acceleration_published(code, s, f2) -> np.ndarray:
    """DR3 final selection for acceleration solutions: s > 20, and F2 < 22
    (7-parameter) or F2 < 25 (9-parameter)."""
    code, s, f2 = np.asarray(code), np.asarray(s, dtype=float), np.asarray(f2, dtype=float)
    return ((code == 7) & (s > 20) & (f2 < 22)) | ((code == 9) & (s > 20) & (f2 < 25))
