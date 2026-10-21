"""Straight lines through scattered data with errors in both coordinates.

The model for every black hole scaling relation is

    y = a + b x + (Gaussian intrinsic scatter of width eps, in y)

with y = log M and x = log(sigma / 200 km/s) (or log(M_bulge / 1e11 Msun)).
Each point also has measurement errors ey and ex, so its total spread about
the line is  s^2 = ey^2 + b^2 ex^2 + eps^2.

Two estimators, because the literature uses both:

* ``fit_fitexy``: least squares with that s^2 in the denominator, eps tuned
  until chi^2 per degree of freedom is 1. This is Tremaine et al. (2002)'s
  method, used by Kormendy & Ho (2013), and MPFITEXY (Williams et al. 2010),
  used by McConnell & Ma (2013). Errors on a and b come from the curvature of
  chi^2 at fixed eps, as MPFIT reports them.

* ``fit_bayes``: the full Gaussian likelihood (including the normalisation
  term that makes eps a real parameter rather than a tuning knob), with upper
  limits treated as censored points: P(true log M below the limit). Sampled
  with emcee; flat priors on a, b and on eps > 0. Two versions:

  - ``model="structural"`` (our default) follows Kelly (2007), whose
    LINMIX_ERR van den Bosch (2016) used: the *true* x values are drawn from
    a mixture of Gaussians fitted alongside the line, so a noisy sigma is
    shrunk towards where sigmas usually are before it is used. Then (x, y)
    of each galaxy is a two-dimensional Gaussian and the likelihood is exact.
  - ``model="conditional"`` treats each measured x as given and adds b^2 ex^2
    to the variance, the same approximation FITEXY makes. It is slightly
    biased when the x errors are not small next to the spread of x.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy import optimize, special

import config as C

# emcee logs a note whenever a chain is shorter than 50 autocorrelation times;
# we report the autocorrelation time with every fit instead.
logging.getLogger("emcee").setLevel(logging.ERROR)


@dataclass
class LineFit:
    a: float
    b: float
    eps: float
    ea: float
    eb: float
    eeps: float | None
    n: int
    n_upper: int
    method: str
    pivot: float
    extra: dict = field(default_factory=dict)

    def predict(self, x):
        return self.a + self.b * np.asarray(x)

    def summary(self) -> dict:
        d = {"a": self.a, "b": self.b, "eps": self.eps, "ea": self.ea, "eb": self.eb,
             "eeps": self.eeps, "n": self.n, "n_upper": self.n_upper, "method": self.method,
             "pivot": self.pivot}
        d.update({k: v for k, v in self.extra.items() if np.isscalar(v)})
        return d


def _arrays(x, ex, y, ey, upper=None, w=None):
    x, ex, y, ey = (np.asarray(v, float) for v in (x, ex, y, ey))
    upper = np.zeros_like(x, bool) if upper is None else np.asarray(upper, bool)
    w = np.ones_like(x) if w is None else np.asarray(w, float)
    ex = np.nan_to_num(ex, nan=0.0)
    ey = np.nan_to_num(ey, nan=0.0)
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], ex[ok], y[ok], ey[ok], upper[ok], w[ok]


# ---------------------------------------------------------------------------
# Least squares with eps tuned to chi^2/dof = 1 (Tremaine et al. 2002; MPFITEXY)
# ---------------------------------------------------------------------------

def _residuals(p, x, ex, y, ey, eps, sw):
    a, b = p
    return sw * (y - a - b * x) / np.sqrt(ey**2 + b**2 * ex**2 + eps**2)


def _lsq(x, ex, y, ey, eps, w, p0):
    sw = np.sqrt(w)
    res = optimize.least_squares(_residuals, p0, args=(x, ex, y, ey, eps, sw), method="lm",
                                 xtol=1e-12, ftol=1e-12)
    return res, float(np.sum(res.fun**2))


def fit_fitexy(x, ex, y, ey, w=None) -> LineFit:
    """MPFITEXY-style fit: minimise chi^2, then tune eps so chi^2 = sum(w) - 2."""
    x, ex, y, ey, _, w = _arrays(x, ex, y, ey, None, w)
    dof = float(np.sum(w)) - 2.0
    p0 = np.polyfit(x, y, 1)[::-1]
    _, chi0 = _lsq(x, ex, y, ey, 0.0, w, p0)
    if chi0 <= dof:
        eps = 0.0
    else:
        eps = optimize.brentq(lambda e: _lsq(x, ex, y, ey, e, w, p0)[1] - dof, 0.0, 5.0, xtol=1e-6)
    res, chi2 = _lsq(x, ex, y, ey, eps, w, p0)
    jac = res.jac
    cov = np.linalg.inv(jac.T @ jac)
    return LineFit(a=res.x[0], b=res.x[1], eps=eps, ea=np.sqrt(cov[0, 0]), eb=np.sqrt(cov[1, 1]),
                   eeps=None, n=len(x), n_upper=0, method="fitexy", pivot=np.nan,
                   extra={"chi2": chi2, "dof": dof, "rho_ab": cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])})


# ---------------------------------------------------------------------------
# Bayesian fit with censored upper limits
# ---------------------------------------------------------------------------

_LN2PI = np.log(2.0 * np.pi)


def _unpack(thetas, n_gauss):
    """Split a (walkers, ndim) array into line and x-distribution parameters."""
    a, b, eps = thetas[:, 0:1], thetas[:, 1:2], thetas[:, 2:3]
    mus = thetas[:, 3:3 + n_gauss]
    ltau = thetas[:, 3 + n_gauss:3 + 2 * n_gauss]
    logits = np.concatenate([np.zeros((len(thetas), 1)), thetas[:, 3 + 2 * n_gauss:]], axis=1)
    return a, b, eps, mus, ltau, logits


def log_posterior(thetas, x, ex, y, ey, upper, w, model="structural", n_gauss=1, xlim=(-10.0, 10.0)):
    """Log posterior for a batch of parameter vectors, shape (walkers, ndim) -> (walkers,).

    Flat priors: -10 < a < 20, -20 < b < 30, 0 < eps < 3; for the structural
    model the x-distribution's means lie within the data's range +/- 1 dex
    (and in increasing order, so components cannot swap labels), its widths
    between 0.001 and 5 dex, and the mixture logits within +/- 10.
    """
    thetas = np.atleast_2d(thetas)
    a, b, eps, mus, ltau, logits = _unpack(thetas, n_gauss if model == "structural" else 0)
    ok = ((eps[:, 0] > 0) & (eps[:, 0] < 3) & (b[:, 0] > -20) & (b[:, 0] < 30)
          & (a[:, 0] > -10) & (a[:, 0] < 20))
    if model == "conditional":
        s2 = ey**2 + (b * ex) ** 2 + eps**2
        z = (y - a - b * x) / np.sqrt(s2)
        lp = np.where(upper, special.log_ndtr(z), -0.5 * (z**2 + _LN2PI + np.log(s2)))
        out = np.sum(w * lp, axis=1)
    else:
        ok &= (np.all((mus > xlim[0]) & (mus < xlim[1]), axis=1)
               & np.all((ltau > np.log(1e-3)) & (ltau < np.log(5.0)), axis=1)
               & np.all(np.abs(logits) <= 10, axis=1))
        if n_gauss > 1:
            ok &= np.all(np.diff(mus, axis=1) > 0, axis=1)
        log_pi = logits - special.logsumexp(logits, axis=1, keepdims=True)
        terms = []
        for k in range(n_gauss):
            t2 = np.exp(2.0 * ltau[:, k:k + 1])
            vx = t2 + ex**2
            cxy = b * t2
            dx = x - mus[:, k:k + 1]
            lpx = -0.5 * (dx**2 / vx + _LN2PI + np.log(vx))
            my = a + b * mus[:, k:k + 1] + cxy / vx * dx
            sy2 = b * b * t2 + eps**2 + ey**2 - cxy**2 / vx
            z = (y - my) / np.sqrt(sy2)
            lpy = np.where(upper, special.log_ndtr(z), -0.5 * (z**2 + _LN2PI + np.log(sy2)))
            terms.append(log_pi[:, k:k + 1] + lpx + lpy)
        out = np.sum(w * special.logsumexp(np.stack(terms), axis=0), axis=1)
    return np.where(ok & np.isfinite(out), out, -np.inf)


def fit_bayes(x, ex, y, ey, upper=None, w=None, *, model="structural", n_gauss=1, walkers=None,
              steps=None, burn=None, seed=None, keep_chain=False) -> LineFit:
    """Posterior medians and 68% intervals for (a, b, eps); upper limits censored.

    The structural model defaults to a single Gaussian for the distribution of
    true x: on van den Bosch's sample, two or three Gaussians (as LINMIX_ERR
    allows) move the slope by less than 0.01 and only slow the sampler down.
    """
    import emcee

    x, ex, y, ey, upper, w = _arrays(x, ex, y, ey, upper, w)
    steps = steps or C.MCMC_STEPS
    burn = burn or C.MCMC_BURN
    rng = np.random.default_rng(C.RNG_SEED if seed is None else seed)
    if model == "conditional":
        n_gauss = 0
    ndim = 3 + (3 * n_gauss - 1 if n_gauss else 0)
    walkers = max(walkers or C.MCMC_WALKERS, 4 * ndim)
    span = (x.min() - 1.0, x.max() + 1.0)
    args = (x, ex, y, ey, upper, w, model, n_gauss, span)

    # Start at the maximum of the posterior; the x-distribution starts as
    # equal Gaussians placed at the data's quantiles.
    start = fit_fitexy(x[~upper], ex[~upper], y[~upper], ey[~upper], w[~upper])
    p_init = [start.a, start.b, max(start.eps, 0.05)]
    if n_gauss:
        q = np.quantile(x, (np.arange(n_gauss) + 0.5) / n_gauss)
        p_init += list(q) + [np.log(max(np.std(x) / n_gauss, 0.02))] * n_gauss + [0.0] * (n_gauss - 1)
    nll = lambda t: -float(log_posterior(t, *args)[0])  # noqa: E731
    best = np.asarray(p_init, float)
    for _ in range(3):  # Nelder-Mead restarts settle the x-distribution parameters
        best = optimize.minimize(nll, best, method="Nelder-Mead",
                                 options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 40000}).x
    p0 = best + 1e-3 * rng.standard_normal((walkers, ndim))
    p0[:, 2] = np.abs(p0[:, 2])
    if n_gauss > 1:
        p0[:, 3:3 + n_gauss] = np.sort(p0[:, 3:3 + n_gauss], axis=1)
    sampler = emcee.EnsembleSampler(walkers, ndim, log_posterior, args=args, vectorize=True)
    # emcee copies numpy's global random state when the sampler is built; give
    # it its own seeded state instead, so every run is exactly reproducible.
    sampler.random_state = np.random.RandomState(int(rng.integers(2**31 - 1))).get_state()
    sampler.run_mcmc(p0, steps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)[:, :3]
    lo, med, hi = np.percentile(chain, [16, 50, 84], axis=0)
    try:
        tau = float(np.max(sampler.get_autocorr_time(discard=burn, quiet=True)[:3]))
    except Exception:  # noqa: BLE001
        tau = np.nan
    extra = {"a_lo": lo[0], "a_hi": hi[0], "b_lo": lo[1], "b_hi": hi[1], "eps_lo": lo[2],
             "eps_hi": hi[2], "rho_ab": float(np.corrcoef(chain[:, 0], chain[:, 1])[0, 1]),
             "ml_a": best[0], "ml_b": best[1], "ml_eps": best[2], "autocorr_steps": tau,
             "n_independent": float(len(chain) / tau) if np.isfinite(tau) and tau > 0 else np.nan,
             "acceptance": float(np.mean(sampler.acceptance_fraction)), "n_gauss": n_gauss}
    if keep_chain:
        extra["chain"] = chain[:: max(1, len(chain) // 20000)]
    return LineFit(a=med[0], b=med[1], eps=med[2], ea=0.5 * (hi[0] - lo[0]), eb=0.5 * (hi[1] - lo[1]),
                   eeps=0.5 * (hi[2] - lo[2]), n=len(x), n_upper=int(upper.sum()),
                   method=f"bayes-{model}", pivot=np.nan, extra=extra)


# ---------------------------------------------------------------------------
# Helpers for the M-sigma and M-M_bulge relations
# ---------------------------------------------------------------------------

def msigma_xy(df, *, sigma_col="logsig", esigma_col="elogsig", equal_logm_err=None):
    """x = log(sigma/200), y = log M from a catalogue table."""
    x = df[sigma_col].to_numpy(float) - np.log10(C.SIGMA_PIVOT)
    ex = df[esigma_col].to_numpy(float)
    y = df["logM"].to_numpy(float)
    ey = df["elogM"].to_numpy(float) if equal_logm_err is None else np.full(len(df), equal_logm_err)
    upper = df["upper"].to_numpy(bool) if "upper" in df else np.zeros(len(df), bool)
    w = df["weight"].to_numpy(float) if "weight" in df else np.ones(len(df))
    return x, ex, y, ey, upper, w


def fit_msigma(df, *, method="bayes", equal_logm_err=None, **kw) -> LineFit:
    x, ex, y, ey, upper, w = msigma_xy(df, equal_logm_err=equal_logm_err)
    if method == "fitexy":
        keep = ~upper
        fit = fit_fitexy(x[keep], ex[keep], y[keep], ey[keep], w[keep])
    else:
        fit = fit_bayes(x, ex, y, ey, upper, w, **kw)
    fit.pivot = C.SIGMA_PIVOT
    return fit


def predict_logm(fit: LineFit, sigma_kms):
    return fit.a + fit.b * (np.log10(sigma_kms) - np.log10(C.SIGMA_PIVOT))
