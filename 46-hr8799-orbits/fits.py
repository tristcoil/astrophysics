"""Nested-sampling orbit fits (dynesty), cached in data/cache/.

    python fits.py                      # every fit in FITS (skips cached ones)
    python fits.py W18/coplanar         # one fit

Each fit is a data set + orbit model + options. Nested sampling gives both
posterior samples and the Bayesian evidence Z, so models can be compared by
Bayes factors - the comparison Wang et al. (2018) wanted but could not afford
with MCMC ("the Bayes factor is computationally difficult to calculate").

Fits:

* ``W18/<model>``: Wang et al.'s data (Keck 2004-2014 + GPI 2013-2016 with its
  original calibration) and priors - the reproduction, and the training set of
  the prediction test (nothing after 2016 Sep 19 is used).
* ``W18R/<model>``: the same, but with GPI's revised calibration (De Rosa et
  al. 2020): does the 0.3 deg north-angle correction change the predictions?
* ``ALL/<model>``: all precise data 1998-2025, Gaia DR3 parallax, a broad
  stellar-mass prior, free planet masses (Jacobi geometry), per-instrument
  calibration offsets and jitter.
* ``ALLGPIW18/coplanar``: ALL, but with GPI's *original* calibration and a wide
  prior on its north angle, so the planets themselves measure GPI's north error.
"""

from __future__ import annotations

import multiprocessing as mp
import sys
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

import astrometry as A
import config as C
import models as M


@dataclass(frozen=True)
class FitSpec:
    dataset: str
    kind: str
    options: dict = field(default_factory=dict)
    nlive: int = C.NLIVE

    def sources(self):
        if self.dataset in ("W18",):
            return C.W18_SOURCES
        if self.dataset == "W18R":
            return ("NIRC2-K16", "GPI-D20")
        if self.dataset == "ALL":
            return C.ALL_SOURCES
        if self.dataset == "ALLGPIW18":
            return tuple(s if s != "GPI-D20" else "GPI-W18" for s in C.ALL_SOURCES)
        raise KeyError(self.dataset)


ALL_OPTIONS = dict(parallax=C.DR3_PARALLAX, mstar=C.ALL_MSTAR, fit_masses=True)

FITS = {}
for kind in C.MODELS:
    FITS[f"W18/{kind}"] = FitSpec("W18", kind)
    FITS[f"W18R/{kind}"] = FitSpec("W18R", kind)
for kind in C.MODELS:
    FITS[f"ALL/{kind}"] = FitSpec("ALL", kind, ALL_OPTIONS)
FITS["ALLGPIW18/coplanar"] = FitSpec("ALLGPIW18", "coplanar", ALL_OPTIONS)

#: the GPI north angle is left nearly free in the ALLGPIW18 fit
CALIBRATION_OVERRIDE = {"ALLGPIW18": {"GPI-W18": (1.0, 0.007 / 14.166, "wide: measured by the fit")}}


def _spec(key: str) -> FitSpec:
    """'W18/coplanar#2' -> the spec of W18/coplanar (run 2 differs only in its seed)."""
    return FITS[key.split("#")[0]]


def build_model(key: str) -> M.OrbitModel:
    spec = _spec(key)
    df = A.load(spec.sources())
    calibrate = spec.dataset.startswith("ALL")
    saved = dict(C.CALIBRATION)
    C.CALIBRATION.update(CALIBRATION_OVERRIDE.get(spec.dataset, {}))
    try:
        data = M.Data.from_frame(df, calibrate=calibrate)
        model = M.OrbitModel(spec.kind, data, **spec.options)
    finally:
        C.CALIBRATION.clear()
        C.CALIBRATION.update(saved)
    return model


def cache_path(key: str):
    return C.CACHE / f"nested_{key.replace('/', '_').replace('#', '_run')}.npz"


# The model lives in a module global so forked workers inherit it.
_MODEL: M.OrbitModel | None = None


def _loglike(theta):
    return _MODEL.loglike(theta)


def _ptform(u):
    return _MODEL.prior_transform(u)


def run(key: str, nproc: int = 1, verbose: bool = True) -> dict:
    import dynesty
    path = cache_path(key)
    if path.exists():
        return load_raw(key)
    global _MODEL
    _MODEL = build_model(key)
    spec = _spec(key)
    # a stable per-fit seed (str hash() changes between Python runs)
    rng = np.random.default_rng([C.RANDOM_SEED, sum(map(ord, key)) * 7919])
    ckpt = C.CACHE / f"dynesty_{key.replace('/', '_').replace('#', '_run')}.save"
    t0 = time.time()
    # The likelihood costs ~20-60 us, less than dynesty's own overhead per call,
    # so a process pool is slower than one process per fit: we run fits side by
    # side instead (nproc = 1 is the default).
    pool = mp.get_context("fork").Pool(nproc, initializer=_init_worker) if nproc > 1 else None
    try:
        if ckpt.exists():
            sampler = dynesty.NestedSampler.restore(str(ckpt), pool=pool)
            print(f"[{key}] resuming from checkpoint")
            sampler.run_nested(resume=True, dlogz=C.DLOGZ, checkpoint_file=str(ckpt),
                               checkpoint_every=300, print_progress=verbose, print_func=_quiet_print)
        else:
            kw = dict(pool=pool, queue_size=nproc) if pool else {}
            sampler = dynesty.NestedSampler(_loglike, _ptform, _MODEL.ndim, nlive=spec.nlive,
                                            bound="multi", sample="rslice",
                                            periodic=_MODEL.periodic or None, rstate=rng, **kw)
            sampler.run_nested(dlogz=C.DLOGZ, checkpoint_file=str(ckpt), checkpoint_every=300,
                               print_progress=verbose, print_func=_quiet_print)
    finally:
        if pool:
            pool.close()
    res = sampler.results
    weights = res.importance_weights()
    equal = res.samples_equal(rstate=np.random.default_rng(C.RANDOM_SEED))
    out = dict(samples=res.samples, logwt=res.logwt, logl=res.logl, weights=weights,
               equal=equal, logz=res.logz[-1], logzerr=res.logzerr[-1],
               ncall=int(np.sum(res.ncall)), runtime_s=time.time() - t0,
               labels=np.array(_MODEL.labels), ess=float(1.0 / np.sum(weights ** 2)))
    np.savez_compressed(path, **out)
    ckpt.unlink(missing_ok=True)
    print(f"[{key}] logZ = {out['logz']:.2f} +- {out['logzerr']:.2f}, ESS {out['ess']:.0f}, "
          f"{out['ncall']:.3g} calls, {out['runtime_s'] / 60:.1f} min")
    return out


def _quiet_print(itresult, niter, ncall, add_live_it=None, dlogz=None, **kw):
    """dynesty progress, one line per 2000 iterations (instead of every one)."""
    if niter % 2000 == 0 or add_live_it == 1:
        dl = getattr(itresult, "delta_logz", float("nan"))
        print(f"  it {niter:7d}  calls {ncall:11d}  logL* {itresult.loglstar:14.2f}  "
              f"logZ {itresult.logz:12.2f}  dlogz {dl:.3g}", flush=True)


def _init_worker():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def load_raw(key: str) -> dict:
    with np.load(cache_path(key), allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


@dataclass
class Fit:
    key: str
    model: M.OrbitModel
    raw: dict

    @property
    def equal(self) -> np.ndarray:
        return self.raw["equal"]

    @property
    def logz(self) -> float:
        return float(self.raw["logz"])

    def table(self, n: int | None = None, seed: int = C.RANDOM_SEED) -> pd.DataFrame:
        s = self.equal
        if n is not None and n < len(s):
            s = s[np.random.default_rng(seed).choice(len(s), n, replace=False)]
        return M.derived(self.model, s)


def load(key: str) -> Fit:
    return Fit(key, build_model(key), load_raw(key) if cache_path(key).exists() else run(key))


if __name__ == "__main__":
    keys = sys.argv[1:] or list(FITS)
    for k in keys:
        run(k)
