"""Fit geometric models to closure quantities by nested sampling.

The likelihood is Gaussian in the closure phases and log closure amplitudes,
with their full covariance (M87 Paper VI eqs. 13 and 16). Closure phase
residuals are wrapped into (-180, 180] degrees. Nested sampling (dynesty, the
code Paper VI used for its closure-only fits) returns both the posterior and the
Bayesian evidence Z, which ranks models of different complexity.
"""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix

from uvdata import Prepared

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")


class ClosureLikelihood:
    """-0.5 chi^2 of a model against one observation's closure phases and log closure amplitudes."""

    def __init__(self, prep: Prepared, model, *, use=("cphase", "lcamp")):
        self.model = model
        self.u, self.v = prep.avg.u, prep.avg.v
        self.sets = []
        if "cphase" in use and len(prep.cphase):
            self.sets.append(("cphase", csr_matrix(prep.cphase.design), prep.cphase.values, prep.cphase.whitener()))
        if "lcamp" in use and len(prep.lcamp):
            self.sets.append(("lcamp", csr_matrix(prep.lcamp.design), prep.lcamp.values, prep.lcamp.whitener()))
        self.n_data = sum(len(s[2]) for s in self.sets)

    def chi2_parts(self, theta) -> dict:
        V = self.model.vis(theta, self.u, self.v)
        out = {}
        for kind, D, obs, W in self.sets:
            if kind == "cphase":
                r = np.angle(np.exp(1j * (obs - D @ np.angle(V))))
            else:
                r = obs - D @ np.log(np.abs(V) + 1e-30)
            w = W @ r
            out[kind] = float(w @ w)
        return out

    def __call__(self, theta) -> float:
        V = self.model.vis(theta, self.u, self.v)
        chi2 = 0.0
        for kind, D, obs, W in self.sets:
            if kind == "cphase":
                r = np.angle(np.exp(1j * (obs - D @ np.angle(V))))
            else:
                r = obs - D @ np.log(np.abs(V) + 1e-30)
            w = W @ r
            chi2 += w @ w
        if not np.isfinite(chi2):
            return -1e100
        return -0.5 * chi2


class MultiLikelihood:
    """Sum of closure likelihoods over several observations (e.g. two nights of Sgr A*)."""

    def __init__(self, preps, model, *, use=("cphase", "lcamp")):
        self.parts = [ClosureLikelihood(p, model, use=use) for p in preps]
        self.model = model
        self.n_data = sum(p.n_data for p in self.parts)

    def chi2_parts(self, theta) -> dict:
        out: dict = {}
        for part in self.parts:
            for k, v in part.chi2_parts(theta).items():
                out[k] = out.get(k, 0.0) + v
        return out

    def __call__(self, theta) -> float:
        return float(sum(part(theta) for part in self.parts))


@dataclass
class FitResult:
    label: str
    model: str
    names: list
    samples: np.ndarray  # equal-weight posterior samples
    logz: float
    logz_err: float
    max_logl: float
    best: np.ndarray
    chi2: dict
    n_data: dict
    ncall: int
    seconds: float

    def summary(self) -> dict:
        out = {
            "label": self.label, "model": self.model, "logz": self.logz, "logz_err": self.logz_err,
            "max_logl": self.max_logl, "chi2": self.chi2, "n_data": self.n_data,
            "reduced_chi2": {k: self.chi2[k] / max(self.n_data[k], 1) for k in self.chi2},
            "ncall": self.ncall, "seconds": self.seconds, "best": dict(zip(self.names, map(float, self.best))),
        }
        return out


def _quiet_threads():
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(1)
    except Exception:  # pragma: no cover
        pass


def run_nested(prep, model, *, nlive: int = 500, dlogz: float = 0.1, seed: int = 1,
               njobs: int = 4, use=("cphase", "lcamp"), sample: str = "rwalk") -> FitResult:
    """Nested sampling of `model` against one observation, or a list of them (closure quantities)."""
    import dynesty
    import dynesty.pool as dypool

    preps = prep if isinstance(prep, (list, tuple)) else [prep]
    if isinstance(model, AmplitudeModel):
        like = AmplitudeLikelihood(preps, model)
    elif isinstance(model, NoisyModel):
        like = NoiseModelLikelihood(preps, model)
    elif len(preps) > 1:
        like = MultiLikelihood(preps, model, use=use)
    else:
        like = ClosureLikelihood(preps[0], model, use=use)
    periodic = [i for i, p in enumerate(model.params) if p.kind == "angle"] or None
    rstate = np.random.default_rng(seed)
    t0 = time.time()
    if njobs > 1:
        with dypool.Pool(njobs, like, model.prior_transform) as pool:
            sampler = dynesty.NestedSampler(pool.loglike, pool.prior_transform, model.ndim, nlive=nlive,
                                            bound="multi", sample=sample, periodic=periodic, rstate=rstate,
                                            pool=pool, queue_size=njobs)
            sampler.run_nested(dlogz=dlogz, print_progress=False)
    else:
        sampler = dynesty.NestedSampler(like, model.prior_transform, model.ndim, nlive=nlive, bound="multi",
                                        sample=sample, periodic=periodic, rstate=rstate)
        sampler.run_nested(dlogz=dlogz, print_progress=False)
    res = sampler.results
    samples = res.samples_equal(rstate=np.random.default_rng(seed + 1))
    ibest = int(np.argmax(res.logl))
    best = res.samples[ibest]
    chi2 = like.chi2_parts(best)
    n_data = {"cphase": sum(len(p.cphase) for p in preps), "lcamp": sum(len(p.lcamp) for p in preps),
              "amp": getattr(like, "n_data", 0)}
    n_data = {k: n_data[k] for k in chi2}
    label = " + ".join(p.label for p in preps)
    return FitResult(label, model.name, model.names, samples, float(res.logz[-1]), float(res.logzerr[-1]),
                     float(res.logl[ibest]), best, chi2, n_data, int(np.sum(res.ncall)), time.time() - t0)


# ---------------------------------------------------------------------------
# Posterior summaries
# ---------------------------------------------------------------------------
def summarize(x: np.ndarray) -> dict:
    lo, med, hi = np.percentile(x, [15.865, 50, 84.135])
    return {"median": float(med), "lo68": float(lo), "hi68": float(hi), "mean": float(np.mean(x)),
            "sd": float(np.std(x))}


def summarize_angle_deg(x_deg: np.ndarray) -> dict:
    """Circular summary: centre on the circular mean, then take percentiles."""
    mean = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(x_deg))))) % 360
    unwrapped = (x_deg - mean + 180.0) % 360.0 - 180.0 + mean
    s = summarize(unwrapped)
    s["circular_mean"] = float(mean)
    return s


def derived_samples(model, samples: np.ndarray) -> dict:
    """Apply model.derived to every posterior sample."""
    rows = [model.derived(t) for t in samples]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}


def save_fit(path, fit: FitResult, model) -> None:
    der = derived_samples(model, fit.samples) if hasattr(model, "derived") else {}
    np.savez_compressed(path.with_suffix(".npz"), samples=fit.samples, names=np.array(fit.names),
                        **{f"derived_{k}": v for k, v in der.items()})
    path.with_suffix(".json").write_text(json.dumps(fit.summary(), indent=1))


def load_fit(path):
    if not path.with_suffix(".npz").exists() or not path.with_suffix(".json").exists():
        return None
    z = np.load(path.with_suffix(".npz"))
    summ = json.loads(path.with_suffix(".json").read_text())
    der = {k[len("derived_"):]: z[k] for k in z.files if k.startswith("derived_")}
    return {"samples": z["samples"], "names": list(z["names"]), "derived": der, "summary": summ}


# ---------------------------------------------------------------------------
# Sgr A*: noise budget fitted together with the source (Sgr A* Paper IV, sec. 7.2.1)
# ---------------------------------------------------------------------------
class NoisyModel:
    """A source model plus Paper IV's parameterised noise budget.

    sigma^2 = sigma_th^2 + (f |V|)^2 + sigma_ref^2 + sigma_var(|u|)^2, with sigma_var
    the broken power law of uvdata.variability_sigma (b = 2.4, c = 2, u0 = 1 Glambda
    fixed). f, sigma_ref and the variability amplitude a4 are free, with log-uniform
    priors spanning Paper IV's ranges. Scattering is applied to the model (the data
    are not deblurred), as in Paper IV's full-track fits.
    """

    def __init__(self, source, *, b: float = 2.4, c: float = 2.0, u0: float = 1.0e9, f_max: float = 0.1):
        from models import Param

        self.source, self.b, self.c, self.u0 = source, b, c, u0
        self.noise_params = [Param("log10_f", -3.0, float(np.log10(f_max)), label="log10 fractional systematic"),
                             Param("log10_sref", -4.0, -1.5, label="log10 refractive floor (of total flux)"),
                             Param("log10_a4", -2.5, -1.0, label="log10 variability at 4 Glambda")]
        self.params = list(source.params) + self.noise_params
        self.names = [p.name for p in self.params]
        self.ndim = len(self.params)
        self.name = source.name + " + noise model"

    def prior_transform(self, cube):
        k = self.source.ndim
        return np.concatenate([self.source.prior_transform(cube[:k]),
                               [p.from_unit(c) for p, c in zip(self.noise_params, cube[k:])]])

    def split(self, theta):
        k = self.source.ndim
        return theta[:k], 10.0 ** np.asarray(theta[k:])

    def derived(self, theta):
        th, (f, sref, a4) = self.split(theta)
        out = self.source.derived(th)
        out.update({"f": f, "sigma_ref": sref, "a4": a4})
        return out


class _BlockSet:
    """Closures of one kind from one observation, grouped by scan for batched linear algebra."""

    def __init__(self, cs, amp_obs):
        scans = np.unique(cs.scan)
        rows = [np.flatnonzero(cs.scan == s) for s in scans]
        cols = [np.flatnonzero(np.any(cs.design[r] != 0, axis=0)) for r in rows]
        M = max(len(r) for r in rows)
        K = max(len(c) for c in cols)
        nb = len(rows)
        self.D = np.zeros((nb, M, K))
        self.idx = np.zeros((nb, K), dtype=int)
        self.colmask = np.zeros((nb, K))
        self.obs = np.zeros((nb, M))
        self.rowmask = np.zeros((nb, M), dtype=bool)
        for b, (r, c) in enumerate(zip(rows, cols)):
            self.D[b, :len(r), :len(c)] = cs.design[np.ix_(r, c)]
            self.idx[b, :len(c)] = c
            self.colmask[b, :len(c)] = 1.0
            self.obs[b, :len(r)] = cs.values[r]
            self.rowmask[b, :len(r)] = True
        # padding rows of each block get a unit diagonal: no effect on chi^2 or the determinant
        self.pad_eye = np.eye(M)[None] * (~self.rowmask)[:, :, None]
        # outer products D_ik D_jk, so each covariance is one batched matrix product
        self.P = (self.D[:, :, None, :] * self.D[:, None, :, :]).reshape(nb, M * M, K)
        self.M = M
        self.kind = cs.kind
        self.n = len(cs.values)
        self.inv_amp2 = 1.0 / np.where(amp_obs > 0, amp_obs, 1.0) ** 2

    def loglike(self, model_vals_per_vis, var_vis):
        """Gaussian log likelihood; covariance D diag(var / amp^2) D^T within each scan."""
        w = (var_vis * self.inv_amp2)[self.idx] * self.colmask
        S = np.matmul(self.P, w[:, :, None]).reshape(-1, self.M, self.M) + self.pad_eye
        r = self.obs - np.matmul(self.D, (model_vals_per_vis[self.idx] * self.colmask)[:, :, None])[..., 0]
        if self.kind == "cphase":
            r = np.angle(np.exp(1j * r))
        r = np.where(self.rowmask, r, 0.0)
        L = np.linalg.cholesky(S)
        z = np.linalg.solve(L, r[..., None])[..., 0]
        chi2 = float(np.sum(z**2))
        logdet = 2.0 * np.sum(np.log(np.diagonal(L, axis1=1, axis2=2)))
        return -0.5 * (chi2 + logdet + self.n * np.log(2 * np.pi)), chi2


class NoiseModelLikelihood:
    """Closure-quantity likelihood for Sgr A*, with the noise budget as free parameters.

    Expects observations prepared by uvdata.prepare_sgra_raw: scan-averaged, light-curve
    normalised, NOT deblurred, errors thermal only.
    """

    def __init__(self, preps, model: NoisyModel):
        from uvdata import debiased_amplitude, scattering_kernel

        self.model = model
        self.parts = []
        for p in preps:
            d = p.avg
            sets = []
            if len(p.cphase):
                sets.append(_BlockSet(p.cphase, np.abs(d.vis)))
            if len(p.lcamp):
                sets.append(_BlockSet(p.lcamp, debiased_amplitude(np.abs(d.vis), d.sigma)))
            self.parts.append((d.u, d.v, d.uvdist, np.abs(d.vis), d.sigma**2, scattering_kernel(d), sets))
        self.n_data = sum(s.n for part in self.parts for s in part[-1])

    def _evaluate(self, theta):
        from uvdata import variability_sigma

        th, (f, sref, a4) = self.model.split(theta)
        total, chi = 0.0, {"cphase": 0.0, "lcamp": 0.0}
        for u, v, uvd, absv, th2, kern, sets in self.parts:
            V = kern * self.model.source.vis(th, u, v)
            var = th2 + (f * absv) ** 2 + sref**2 + variability_sigma(uvd, a4=a4, b=self.model.b, u0=self.model.u0,
                                                                     c=self.model.c) ** 2
            for s in sets:
                vals = np.angle(V) if s.kind == "cphase" else np.log(np.abs(V) + 1e-30)
                ll, c2 = s.loglike(vals, var)
                total += ll
                chi[s.kind] += c2
        return total, chi

    def __call__(self, theta) -> float:
        ll, _ = self._evaluate(theta)
        return float(ll) if np.isfinite(ll) else -1e100

    def chi2_parts(self, theta) -> dict:
        return self._evaluate(theta)[1]


class AmplitudeLikelihood:
    """Sgr A*: fit the visibility *amplitudes* beyond `min_uv` (Paper IV sec. 4.2, Fig. 7).

    Closure quantities alone leave Sgr A*'s size degenerate with a broad halo (a small,
    lopsided crescent fits them as well), so the ring is also measured where it leaves
    its clearest mark: the nulls of the amplitudes. The light-curve-normalised data
    have calibrated amplitudes; residual station-gain errors (~5-10%) are absorbed in
    the fitted fractional term f, and the model's flux scale is free (it takes up any
    emission resolved out beyond `min_uv`). Baselines below 2 Glambda are left out:
    there the normalisation and the LMT's size-based calibration shape the amplitudes.
    """

    def __init__(self, preps, model: NoisyModel, *, min_uv: float = 2.0e9, min_snr: float = 1.0):
        from uvdata import debiased_amplitude, scattering_kernel

        self.model = model
        self.parts = []
        skip = set(getattr(model, "exclude_stations", ()))
        for p in preps:
            d = p.avg
            amp = debiased_amplitude(np.abs(d.vis), d.sigma)
            use = (d.uvdist > min_uv) & (np.abs(d.vis) / d.sigma >= min_snr)
            use &= np.array([a not in skip and b not in skip for a, b in zip(d.t1, d.t2)])
            kern = scattering_kernel(d) if getattr(model, "scattered", True) else np.ones(len(d))
            self.parts.append((d.u[use], d.v[use], d.uvdist[use], amp[use], np.abs(d.vis[use]), d.sigma[use] ** 2,
                               kern[use]))
        self.n_data = sum(len(p[3]) for p in self.parts)

    def _evaluate(self, theta):
        from uvdata import variability_sigma

        th, (f, sref, a4) = self.model.split(theta)
        ll, chi = 0.0, 0.0
        for u, v, uvd, amp, absv, th2, kern in self.parts:
            m = np.abs(kern * self.model.source.vis(th, u, v))
            var = th2 + (f * absv) ** 2 + sref**2 + variability_sigma(uvd, a4=a4, b=self.model.b, u0=self.model.u0,
                                                                      c=self.model.c) ** 2
            r2 = (amp - m) ** 2 / var
            chi += float(np.sum(r2))
            ll += -0.5 * float(np.sum(r2 + np.log(2 * np.pi * var)))
        return ll, chi

    def __call__(self, theta) -> float:
        ll, _ = self._evaluate(theta)
        return ll if np.isfinite(ll) else -1e100

    def chi2_parts(self, theta) -> dict:
        return {"amp": self._evaluate(theta)[1]}


class AmplitudeModel(NoisyModel):
    """A NoisyModel fitted with AmplitudeLikelihood.

    `scattered` applies the Galactic Centre scattering kernel (Sgr A* only);
    `exclude_stations` drops baselines to stations with unreliable amplitude gains.
    """

    def __init__(self, source, *, scattered: bool = True, exclude_stations=(), **kw):
        super().__init__(source, **kw)
        self.scattered = scattered
        self.exclude_stations = tuple(exclude_stations)
