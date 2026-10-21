"""Gaussian-process reconstruction of H(z), with the full data covariance.

A Gaussian process (GP) is a way of fitting a smooth curve through data without
choosing a formula. Instead of "H(z) is a parabola" it says only "H(z) is a
smooth function whose values at nearby redshifts are correlated": the
*kernel* k(z, z') sets how strongly, through an amplitude sigma_f and a
correlation length ell. Conditioning that prior on the data gives a posterior
distribution over curves. Every draw is a plausible expansion history, and
the spread of the draws is the uncertainty, including the uncertainty of
derived quantities such as H0 (the curve at z = 0), the deceleration parameter
(its slope), and the transition redshift (the minimum of H / (1 + z)).

Nothing here assumes a cosmological model. The kernel and the prior mean are
still choices, and the results depend on them somewhat; ``analysis.py`` tests
two kernels and marginalises over the hyperparameters instead of fixing them.

The data enter with their full covariance C (not just error bars):

    mean(z*) = m(z*) + K(z*, z) [K(z, z) + C]^-1 (H - m(z)),
    cov(z*)  = K(z*, z*) - K(z*, z) [K(z, z) + C]^-1 K(z, z*).
"""

from __future__ import annotations

import emcee
import numpy as np
from scipy import linalg, optimize

import config as cfg


def kernel(x1, x2, ln_sf, ln_ell, kind="se"):
    """Covariance between f(x1) and f(x2).

    ``se``: squared exponential, sigma_f^2 exp(-r^2 / (2 ell^2)).
    ``m92``: Matern nu = 9/2 (four times differentiable; less smooth than SE),
    in the form used by GaPP (Seikel, Clarkson & Smith 2012).
    """
    sf2 = np.exp(2.0 * ln_sf)
    ell = np.exp(ln_ell)
    r = np.abs(np.subtract.outer(np.asarray(x1, float), np.asarray(x2, float)))
    if kind == "se":
        return sf2 * np.exp(-0.5 * (r / ell) ** 2)
    if kind == "m92":
        u = 3.0 * r / ell
        return sf2 * np.exp(-u) * (1 + u + 27 * r**2 / (7 * ell**2)
                                   + 18 * r**3 / (7 * ell**3) + 27 * r**4 / (35 * ell**4))
    raise ValueError(f"unknown kernel {kind!r}")


def zero_mean(z):
    """The default prior mean: zero (the data alone set the level)."""
    return np.zeros_like(np.asarray(z, dtype=float))


class GaussianProcess:
    """GP regression of y(x) with data covariance ``cov`` and prior mean ``mean``.

    ``mean_sd=(c, l)`` instead marginalises a constant (prior sd c) and a
    straight-line slope (prior sd l) in the mean: equivalent to adding
    c^2 + l^2 x x' to the kernel, so their uncertainty is propagated rather
    than fixed at a fitted value.
    """

    def __init__(self, x, y, cov, *, kind="se", mean=None, mean_sd=(0.0, 0.0)):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.cov = np.asarray(cov, dtype=float)
        self.kind = kind
        self.mean = zero_mean if mean is None else mean
        self.mean_sd = mean_sd
        self.resid = self.y - self.mean(self.x)

    def _k(self, x1, x2, ln_sf, ln_ell):
        c, l = self.mean_sd
        k = kernel(x1, x2, ln_sf, ln_ell, self.kind)
        if c or l:
            k = k + c**2 + l**2 * np.multiply.outer(np.asarray(x1, float), np.asarray(x2, float))
        return k

    # -- hyperparameters ---------------------------------------------------

    def log_marginal(self, theta) -> float:
        """ln p(y | theta): how well the GP prior with these hyperparameters
        explains the data, automatically penalising needless flexibility."""
        ln_sf, ln_ell = theta
        K = self._k(self.x, self.x, ln_sf, ln_ell) + self.cov
        try:
            L = linalg.cholesky(K, lower=True)
        except linalg.LinAlgError:
            return -np.inf
        alpha = linalg.cho_solve((L, True), self.resid)
        return float(-0.5 * self.resid @ alpha - np.sum(np.log(np.diag(L)))
                     - 0.5 * self.x.size * np.log(2 * np.pi))

    def log_posterior(self, theta) -> float:
        (a, b), (c, d) = cfg.GP_LN_SIGMA_F, cfg.GP_LN_LENGTH
        if not (a < theta[0] < b and c < theta[1] < d):
            return -np.inf
        return self.log_marginal(theta)

    def optimize(self) -> np.ndarray:
        """Maximum of the marginal likelihood (type-II maximum likelihood)."""
        # scan the prior box, then polish the best grid point
        (a, b), (c, d) = cfg.GP_LN_SIGMA_F, cfg.GP_LN_LENGTH
        grid = [(s, l) for s in np.linspace(a + 0.05, b - 0.05, 25)
                for l in np.linspace(c + 0.05, d - 0.05, 25)]
        start = max(grid, key=self.log_posterior)
        res = optimize.minimize(lambda t: -self.log_posterior(t), start, method="Nelder-Mead",
                                options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 4000})
        return res.x

    def sample_hyperparameters(self, n_keep, seed, *, n_walkers=24, n_steps=3000,
                               burn=1000) -> np.ndarray:
        """Draws of (ln sigma_f, ln ell) from their posterior (flat priors on
        the logarithms inside the box in config), by emcee."""
        start = self.optimize()
        rng = np.random.default_rng(seed)
        p0 = start + 1e-3 * rng.standard_normal((n_walkers, 2))
        sampler = emcee.EnsembleSampler(n_walkers, 2, self.log_posterior)
        # emcee copies numpy's global RNG at construction; seed it explicitly
        sampler.random_state = np.random.RandomState(seed).get_state()
        sampler.run_mcmc(p0, n_steps)
        chain = sampler.get_chain(discard=burn, flat=True)
        pick = rng.choice(chain.shape[0], size=n_keep, replace=chain.shape[0] < n_keep)
        return chain[pick]

    # -- prediction --------------------------------------------------------

    def predict(self, xs, theta):
        """Posterior mean and covariance at ``xs``."""
        ln_sf, ln_ell = theta
        K = self._k(self.x, self.x, ln_sf, ln_ell) + self.cov
        Ks = self._k(xs, self.x, ln_sf, ln_ell)
        Kss = self._k(xs, xs, ln_sf, ln_ell)
        L = linalg.cholesky(K, lower=True)
        alpha = linalg.cho_solve((L, True), self.resid)
        V = linalg.solve_triangular(L, Ks.T, lower=True)
        return self.mean(xs) + Ks @ alpha, Kss - V.T @ V

    def sample(self, xs, theta, n, rng) -> np.ndarray:
        """``n`` posterior draws of the curve at ``xs`` (shape n x len(xs))."""
        mu, cov = self.predict(xs, theta)
        vals, vecs = linalg.eigh(0.5 * (cov + cov.T))
        root = vecs * np.sqrt(np.clip(vals, 0.0, None))
        return mu + rng.standard_normal((n, xs.size)) @ root.T


def derived_from_curves(zs, curves, *, z_lo=0.0, z_hi=cfg.Z_MAX_DATA) -> dict:
    """H0, q(z), and the transition redshift from sampled H(z) curves.

    The transition redshift is where adot = H / (1 + z) is smallest. A draw whose
    minimum sits on an edge of [z_lo, z_hi] has no transition inside the data:
    at z_lo it is decelerating today; at z_hi it is accelerating throughout.
    The minimum is refined with a parabola through the three nearest points.
    """
    zs = np.asarray(zs)
    adot = curves / (1.0 + zs)
    dH = np.gradient(curves, zs, axis=1, edge_order=2)
    q = (1.0 + zs) * dH / curves - 1.0
    inside = (zs >= z_lo) & (zs <= z_hi)
    idx_in = np.nonzero(inside)[0]
    k = idx_in[np.argmin(adot[:, inside], axis=1)]
    interior = (k > idx_in[0]) & (k < idx_in[-1])
    zt = np.full(curves.shape[0], np.nan)
    rows = np.nonzero(interior)[0]
    km = k[rows]
    y0, y1, y2 = adot[rows, km - 1], adot[rows, km], adot[rows, km + 1]
    step = zs[1] - zs[0]
    denom = y0 - 2 * y1 + y2
    shift = np.where(denom > 0, 0.5 * (y0 - y2) / np.where(denom > 0, denom, 1), 0.0)
    zt[rows] = zs[km] + np.clip(shift, -1, 1) * step
    return {
        "H0": curves[:, 0] if zs[0] == 0 else np.interp(0.0, zs, curves),
        "q0": q[:, 0],
        "q": q,
        "adot": adot,
        "zt": zt,
        "frac_decelerating_today": float(np.mean(k == idx_in[0])),
        "frac_no_turn_in_data": float(np.mean(k == idx_in[-1])),
        "frac_transition": float(np.mean(interior)),
    }
