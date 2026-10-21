"""Straight-line fits with errors in both coordinates and intrinsic scatter.

The model is y = a + b (x - x0) + e, with e drawn from a Gaussian of width
sigma_int (the "intrinsic scatter": real object-to-object differences plus any
unrecognized measurement error). Each point's measured x and y are scattered
around the truth by their own errors. Integrating over the unknown true x
gives the familiar likelihood

    ln L = -1/2 sum[ (y - a - b(x - x0))^2 / s^2 + ln(2 pi s^2) ],
    s^2 = sigma_y^2 + b^2 sigma_x^2 + sigma_int^2,

which is exact when the true x values are spread much more widely than their
errors (true for every relation here: x spans 2-4 dex, errors ~0.1 dex). Kelly
(2007)'s LINMIX_ERR, used by Bentz et al. (2013) and Shen et al. (2024),
reaches the same answer in that limit by modelling the x distribution as a
Gaussian mixture; ``test_physics.py`` checks that the two agree on Bentz et
al.'s own table.

Posteriors are sampled with emcee: flat priors on a and b, and a flat prior on
sigma_int >= 0.
"""

from __future__ import annotations

from dataclasses import dataclass

import emcee
import numpy as np
from scipy import optimize


@dataclass
class LineFit:
    a: tuple  # (median, minus, plus)
    b: tuple
    sig: tuple
    x0: float
    n: int
    samples: np.ndarray
    ml: tuple  # maximum-likelihood (a, b, sig)

    def predict(self, x):
        return self.ml[0] + self.ml[1] * (np.asarray(x) - self.x0)

    def summary(self, nd=3) -> str:
        f = lambda t: f"{t[0]:.{nd}f} -{t[1]:.{nd}f} +{t[2]:.{nd}f}"  # noqa: E731
        return f"a = {f(self.a)}, b = {f(self.b)}, sigma_int = {f(self.sig)} (N = {self.n})"


def _lnlike(p, x, sx, y, sy, x0, fix_b):
    if fix_b is None:
        a, b, s = p
    else:
        (a, s), b = p, fix_b
    if s < 0:
        return -np.inf
    var = sy**2 + (b * sx) ** 2 + s**2
    r = y - a - b * (x - x0)
    return -0.5 * np.sum(r * r / var + np.log(2 * np.pi * var))


def fit_line(x, sx, y, sy, x0: float = 0.0, fix_b: float | None = None, nwalkers: int = 32,
             nsteps: int = 3000, seed: int = 7) -> LineFit:
    """Fit y = a + b(x - x0) with intrinsic scatter; ``fix_b`` holds the slope fixed."""
    x, sx, y, sy = (np.asarray(v, float) for v in (x, sx, y, sy))
    nll = lambda p: -_lnlike(p, x, sx, y, sy, x0, fix_b)  # noqa: E731
    if fix_b is None:
        b0, a0 = np.polyfit(x - x0, y, 1)
        p0 = np.array([a0, b0, max(np.std(y - a0 - b0 * (x - x0)), 0.05)])
    else:
        a0 = np.mean(y - fix_b * (x - x0))
        p0 = np.array([a0, max(np.std(y - a0 - fix_b * (x - x0)), 0.05)])
    ml = optimize.minimize(nll, p0, method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-9, maxiter=20000)).x
    ml[-1] = abs(ml[-1])
    rng = np.random.default_rng(seed)
    ndim = len(p0)
    start = ml + 1e-3 * rng.standard_normal((nwalkers, ndim))
    start[:, -1] = np.abs(start[:, -1])
    sampler = emcee.EnsembleSampler(nwalkers, ndim, _lnlike, args=(x, sx, y, sy, x0, fix_b))
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(start, nsteps, progress=False)
    ch = sampler.get_chain(discard=nsteps // 3, thin=5, flat=True)
    if fix_b is not None:
        ch = np.column_stack([ch[:, 0], np.full(len(ch), fix_b), ch[:, 1]])
        ml = np.array([ml[0], fix_b, ml[1]])

    def q(v):
        lo, md, hi = np.percentile(v, [15.865525, 50, 84.134475])
        return float(md), float(md - lo), float(hi - md)

    return LineFit(q(ch[:, 0]), q(ch[:, 1]), q(ch[:, 2]), x0, len(x), ch, tuple(float(v) for v in ml))


def linmix_fit(x, sx, y, sy, x0: float = 0.0, K: int = 2, seed: int = 2, maxiter: int = 10000) -> LineFit:
    """Kelly (2007)'s LINMIX_ERR (the Python port by J. Meyers), as used by
    Bentz et al. 2013 and Shen et al. 2024; same summary format as ``fit_line``.

    LINMIX models the true x values as a mixture of K Gaussians instead of
    leaving them free, which steepens the slope slightly when the x errors are
    not negligible (by 0.011 on Bentz et al.'s table). ``sigma_int`` is the
    square root of LINMIX's intrinsic variance.
    """
    import linmix

    x, sx, y, sy = (np.asarray(v, float) for v in (x, sx, y, sy))
    lm = linmix.LinMix(x - x0, y, sx, sy, K=K, seed=seed, parallelize=False)
    # The port ignores ``seed`` when parallelize=False: each chain draws its own
    # seed from the operating system. Re-seed and re-initialize every chain so
    # that a fit is exactly reproducible.
    for i, chain in enumerate(lm._chains):
        chain.rng = np.random.RandomState(seed + i)
        chain.initial_guess()
    lm.run_mcmc(silent=True, maxiter=maxiter)
    ch = lm.chain
    sam = np.column_stack([ch["alpha"], ch["beta"], np.sqrt(ch["sigsqr"])])

    def q(v):
        lo, md, hi = np.percentile(v, [15.865525, 50, 84.134475])
        return float(md), float(md - lo), float(hi - md)

    med = tuple(float(np.median(sam[:, i])) for i in range(3))
    return LineFit(q(sam[:, 0]), q(sam[:, 1]), q(sam[:, 2]), x0, len(x), sam, med)


def log_err(value, minus, plus):
    """Symmetric error in log10 from an asymmetric linear error (mean of both sides)."""
    value, minus, plus = (np.asarray(v, float) for v in (value, minus, plus))
    return 0.5 * (np.log10(value + plus) - np.log10(np.clip(value - minus, 1e-3 * value, None)))


def gaussian_mean_scatter(v, e, nwalkers: int = 24, nsteps: int = 3000, seed: int = 3):
    """Mean and intrinsic scatter of a Gaussian population measured with errors e.

    Used for <log f>: ln L = -1/2 sum[(v - mu)^2 / (e^2 + s^2) + ln(e^2 + s^2)].
    Returns ((mu, -, +), (s, -, +), maximum-likelihood (mu, s)).
    """
    v, e = np.asarray(v, float), np.asarray(e, float)

    def lnl(p):
        mu, s = p
        if s < 0:
            return -np.inf
        var = e**2 + s**2
        return -0.5 * np.sum((v - mu) ** 2 / var + np.log(2 * np.pi * var))

    ml = optimize.minimize(lambda p: -lnl(p), [np.mean(v), np.std(v)], method="Nelder-Mead").x
    ml[1] = abs(ml[1])
    rng = np.random.default_rng(seed)
    start = ml + 1e-3 * rng.standard_normal((nwalkers, 2))
    start[:, 1] = np.abs(start[:, 1])
    s = emcee.EnsembleSampler(nwalkers, 2, lnl)
    s.random_state = np.random.RandomState(seed).get_state()
    s.run_mcmc(start, nsteps, progress=False)
    ch = s.get_chain(discard=nsteps // 3, thin=5, flat=True)

    def q(x):
        lo, md, hi = np.percentile(x, [15.865525, 50, 84.134475])
        return float(md), float(md - lo), float(hi - md)

    return q(ch[:, 0]), q(ch[:, 1]), (float(ml[0]), float(ml[1]))
