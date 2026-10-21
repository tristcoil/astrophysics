"""A small Gaussian process for interpolating light curves.

Supernova light curves change fast early (days) and slowly late (months). A GP
in the logarithm of time since explosion sees both at a similar scale, so one
length scale serves the whole curve. The kernel is Matern-3/2 (smooth but not
infinitely so); an extra "jitter" term absorbs the few-hundredths of a
magnitude by which different telescopes disagree. Hyperparameters are set by
maximising the marginal likelihood. Used only to interpolate, never to
extrapolate.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg, optimize


def _matern32(x1, x2, amp, ell):
    r = np.abs(x1[:, None] - x2[None, :]) / ell
    s = np.sqrt(3.0) * r
    return amp**2 * (1.0 + s) * np.exp(-s)


class LightCurveGP:
    """GP fit of magnitude vs ln(phase) for one band."""

    def __init__(self, phase, mag, err, ell_bounds=(0.02, 1.5), jitter_bounds=(0.003, 0.3)):
        phase = np.asarray(phase, float)
        keep = phase > 0
        self.x = np.log(phase[keep])
        self.y = np.asarray(mag, float)[keep]
        self.e = np.asarray(err, float)[keep]
        self.t_min, self.t_max = phase[keep].min(), phase[keep].max()
        w = 1.0 / self.e**2
        self.mean = float(np.sum(w * self.y) / np.sum(w))
        self.bounds = [(np.log(0.01), np.log(10.0)), tuple(np.log(ell_bounds)), tuple(np.log(jitter_bounds))]
        self._fit()

    def _nll(self, p):
        amp, ell, jit = np.exp(p)
        k = _matern32(self.x, self.x, amp, ell)
        k[np.diag_indices_from(k)] += self.e**2 + jit**2
        try:
            cf = linalg.cho_factor(k, lower=True, check_finite=False)
        except linalg.LinAlgError:
            return 1e25
        r = self.y - self.mean
        alpha = linalg.cho_solve(cf, r, check_finite=False)
        return 0.5 * r @ alpha + np.sum(np.log(np.diag(cf[0])))

    def _fit(self):
        p0 = np.log([max(np.std(self.y), 0.05), 0.15, 0.03])
        best = None
        for ell0 in (0.05, 0.15, 0.4):
            p0[1] = np.log(ell0)
            res = optimize.minimize(self._nll, p0, method="L-BFGS-B", bounds=self.bounds)
            if best is None or res.fun < best.fun:
                best = res
        self.amp, self.ell, self.jitter = np.exp(best.x)
        k = _matern32(self.x, self.x, self.amp, self.ell)
        k[np.diag_indices_from(k)] += self.e**2 + self.jitter**2
        self._cf = linalg.cho_factor(k, lower=True)
        self._alpha = linalg.cho_solve(self._cf, self.y - self.mean)

    def predict(self, phase):
        """Mean magnitude and its 1-sigma uncertainty at the given phases (days)."""
        xs = np.log(np.atleast_1d(np.asarray(phase, float)))
        ks = _matern32(xs, self.x, self.amp, self.ell)
        mu = self.mean + ks @ self._alpha
        v = linalg.cho_solve(self._cf, ks.T)
        var = self.amp**2 - np.sum(ks * v.T, axis=1)
        return mu, np.sqrt(np.clip(var, 0.0, None))

    def covered(self, phase, frac: float = 0.10, min_days: float = 1.5):
        """True where ``phase`` lies within the data and near enough to a point.

        Near enough means within max(min_days, frac * phase) of an observation,
        so a ten-day gap at day 300 is bridged but one at day 10 is not.
        """
        t = np.atleast_1d(np.asarray(phase, float))
        tol = np.maximum(min_days, frac * t)
        obs = np.exp(self.x)
        nearest = np.min(np.abs(t[:, None] - obs[None, :]), axis=1)
        inside = (t >= self.t_min - 0.5 * tol) & (t <= self.t_max + 0.5 * tol)
        return inside & (nearest <= tol)
