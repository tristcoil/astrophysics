"""Simulated light curves (LABELLED SIMULATION, used only to test the code).

AGN continuum variability is well described by a damped random walk (DRW;
Kelly et al. 2009; MacLeod et al. 2010): fluctuations with a characteristic
timescale tau_d and long-term rms sigma. An emission line that reverberates
is the continuum convolved with a transfer function; here a top-hat of given
mean lag and width. ``test_physics.py`` uses these to check that our ICCF
recovers a known lag with honest error bars. No simulated curve ever enters
the analysis of real data.
"""

from __future__ import annotations

import numpy as np


def drw(t, tau_d: float = 100.0, sigma: float = 0.2, mean: float = 1.0, seed: int = 0):
    """Exact DRW (Ornstein-Uhlenbeck) realization at the sorted times ``t``."""
    rng = np.random.default_rng(seed)
    t = np.asarray(t, float)
    x = np.empty_like(t)
    x[0] = rng.normal(0.0, sigma)
    for i in range(1, len(t)):
        a = np.exp(-(t[i] - t[i - 1]) / tau_d)
        x[i] = a * x[i - 1] + rng.normal(0.0, sigma * np.sqrt(1 - a * a))
    return mean + x


def reverberate(t_fine, cont_fine, lag: float, width: float, scale: float = 1.0):
    """Line = continuum convolved with a top-hat of mean ``lag`` and full ``width``
    (both in days), on the same uniform fine grid ``t_fine``."""
    dt = t_fine[1] - t_fine[0]
    n0 = int(round((lag - width / 2) / dt))
    n1 = int(round((lag + width / 2) / dt))
    k = np.zeros(n1 + 1)
    k[max(n0, 0):n1 + 1] = 1.0
    k /= k.sum()
    line = np.convolve(cont_fine - cont_fine.mean(), k, mode="full")[:len(cont_fine)]
    return cont_fine.mean() + scale * line


def campaign(lag: float = 10.0, width: float = 4.0, days: float = 180.0, cadence: float = 1.0,
             noise: float = 0.01, tau_d: float = 80.0, sigma: float = 0.15, seed: int = 0, gaps: bool = True):
    """A reverberation-mapping campaign: continuum and line sampled at the same
    nightly-ish epochs (random 30% of nights lost to weather when ``gaps``),
    each with fractional Gaussian noise ``noise``. The first 100 days of the
    simulated DRW are a burn-in before observations start."""
    rng = np.random.default_rng(seed + 12345)
    t_fine = np.arange(-100.0, days + 1.0, 0.1)
    c_fine = drw(t_fine, tau_d, sigma, 1.0, seed)
    l_fine = reverberate(t_fine, c_fine, lag, width)
    t_obs = np.arange(0.0, days, cadence) + rng.uniform(-0.2, 0.2, int(np.ceil(days / cadence)))[: len(np.arange(0.0, days, cadence))]
    if gaps:
        t_obs = t_obs[rng.random(len(t_obs)) > 0.3]
    t_obs = np.sort(t_obs)
    c = np.interp(t_obs, t_fine, c_fine)
    ln = np.interp(t_obs, t_fine, l_fine)
    ec = noise * np.abs(c)
    el = noise * np.abs(ln)
    c = c + rng.normal(0, ec)
    ln = ln + rng.normal(0, el)
    return t_obs, c, ec, t_obs.copy(), ln, el
