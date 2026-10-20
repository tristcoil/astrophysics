"""Mock companions: the population we push through each detection model.

To ask "what fraction of companions at period P would Gaia have found around
*these* stars?", we give real sample stars (their sky positions, distances,
magnitudes, number of Gaia visits) imaginary companions and simulate what
Gaia would have seen. The companions' mass ratios and eccentricities follow
Moe & Di Stefano (2017)'s fits for the primary's mass and the orbit's period,
so the completeness of a period bin is an average over realistic companions;
we vary that assumption as a systematic (flat mass ratios, as Raghavan et al.
found, and thermal eccentricities).

Orbits are oriented at random (cos i uniform), seen at a random moment.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import binaries as B
import literature as L
import stellar

Q_MIN = 0.1  # companions are counted down to q = 0.1 (stars, not brown dwarfs, for Sun-like primaries)
TWIN_Q = 0.95


@dataclass
class MockSet:
    """A batch of mock binaries, one row per system (all arrays of equal length)."""
    star_index: np.ndarray   # row of the primary in its sample table
    logp: np.ndarray
    period: np.ndarray       # days
    m1: np.ndarray
    q: np.ndarray
    ecc: np.ndarray
    inc: np.ndarray          # radians, 0 = face-on
    omega: np.ndarray        # argument of periastron (radians)
    node: np.ndarray         # longitude of the ascending node (radians)
    t_peri: np.ndarray       # days relative to the Gaia DR3 reference epoch J2016.0

    @property
    def m2(self):
        return self.q * self.m1

    def __len__(self):
        return len(self.logp)


def sample_q(rng, logp, m1, *, model: str = "mds") -> np.ndarray:
    """Mass ratios q in [0.1, 1].

    ``mds``: Moe & Di Stefano's broken power law (slope gamma_smallq below
    q = 0.3, gamma_largeq above) plus their excess of twins at q > 0.95.
    ``flat``: uniform in q (Raghavan et al. 2010's description for 0.2-0.95).
    """
    n = len(logp)
    if model == "flat":
        return rng.uniform(Q_MIN, 1.0, n)
    out = np.empty(n)
    for i in range(n):
        gs = float(L.mds_gamma_smallq(np.array([logp[i]]), m1[i])[0])
        gl = float(L.mds_gamma_largeq(np.array([logp[i]]), m1[i])[0])
        ft = float(L.mds_ftwin(np.array([logp[i]]), m1[i])[0])
        out[i] = _draw_q(rng, gs, gl, ft)
    return out


def _powerlaw_mass(lo, hi, g):
    return np.log(hi / lo) if abs(g + 1) < 1e-9 else (hi ** (g + 1) - lo ** (g + 1)) / (g + 1)


def _powerlaw_draw(rng, lo, hi, g):
    x = rng.uniform()
    if abs(g + 1) < 1e-9:
        return lo * (hi / lo) ** x
    return (lo ** (g + 1) + x * (hi ** (g + 1) - lo ** (g + 1))) ** (1 / (g + 1))


def _draw_q(rng, g_small, g_large, f_twin):
    # Continuous broken power law on [0.1, 1] (MDS Eq. 2) plus an excess of
    # twins spread uniformly over [0.95, 1]. F_twin is the excess as a
    # fraction of *all* binaries with q > 0.3 (El-Badry et al. 2019's wording
    # of the same definition), so w_twin / (w_large + w_twin) = F_twin.
    w_large = _powerlaw_mass(0.3, 1.0, g_large)
    c_small = 0.3 ** (g_large - g_small)  # continuity at q = 0.3
    w_small = c_small * _powerlaw_mass(Q_MIN, 0.3, g_small)
    w_twin = f_twin / (1.0 - f_twin) * w_large
    u = rng.uniform() * (w_small + w_large + w_twin)
    if u < w_small:
        return _powerlaw_draw(rng, Q_MIN, 0.3, g_small)
    if u < w_small + w_large:
        return _powerlaw_draw(rng, 0.3, 1.0, g_large)
    return rng.uniform(TWIN_Q, 1.0)


def sample_ecc(rng, logp, m1, *, model: str = "mds") -> np.ndarray:
    """Eccentricities: p(e) ~ e^eta below e_max(P) (MDS Eqs. 3, 17-18);
    circular below P = 2 d. ``thermal``: p(e) = 2e (eta = 1)."""
    period = 10 ** np.asarray(logp)
    emax = B.e_max_envelope(period)
    eta = np.full(len(logp), 1.0) if model == "thermal" else L.mds_eta(logp, m1)
    eta = np.clip(eta, -0.95, 5)
    x = rng.uniform(size=len(logp))
    e = emax * x ** (1 / (eta + 1))
    return np.where(period > 2.0, np.minimum(e, 0.97), 0.0)


def make_mocks(rng, n_per_bin: int, logp_edges, primaries_m1, *, q_model="mds", e_model="mds",
               star_weights=None) -> MockSet:
    """``n_per_bin`` mock companions in each log P bin, around primaries drawn
    (with replacement) from the sample; logP uniform within each bin."""
    logp_edges = np.asarray(logp_edges)
    nb = len(logp_edges) - 1
    n = n_per_bin * nb
    logp = np.concatenate([rng.uniform(logp_edges[k], logp_edges[k + 1], n_per_bin) for k in range(nb)])
    p = None if star_weights is None else star_weights / star_weights.sum()
    star = rng.choice(len(primaries_m1), size=n, replace=True, p=p)
    m1 = np.asarray(primaries_m1)[star]
    q = sample_q(rng, logp, m1, model=q_model)
    ecc = sample_ecc(rng, logp, m1, model=e_model)
    inc, omega, node = B.random_orientation(rng, n)
    period = 10 ** logp
    t_peri = rng.uniform(0, 1, n) * period
    return MockSet(star, logp, period, m1, q, ecc, inc, omega, node, t_peri)


def flux_ratios(mock: MockSet, grvs_coeffs) -> tuple[np.ndarray, np.ndarray]:
    """G-band and RVS-band flux ratios F2/F1 (companions below Mamajek's table,
    i.e. brown dwarfs, count as dark)."""
    return stellar.flux_ratio_g(mock.m1, mock.m2), stellar.flux_ratio_rvs(mock.m1, mock.m2, grvs_coeffs)
