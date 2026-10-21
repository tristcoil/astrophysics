"""Cluster membership: a mixture model that knows about measurement errors.

The classical way to decide which stars belong to a cluster (Vasilevskis,
Klemola & Preston 1958; Sanders 1971) treats the stars' motions as two
overlapping populations -- a tight clump, the cluster, on a broad background,
the field -- and asks of each star which population it more plausibly came
from. The answer is a membership *probability*, not a yes/no.

We do the same, with one refinement that matters for Gaia: every star carries
its own measurement uncertainty, and a faint star's is ten times a bright
star's. "Extreme deconvolution" (Bovy, Hogg & Roweis 2011, Ann. Appl. Stat. 5,
1657) fits each population's *intrinsic* distribution while treating every
star as a blurred measurement of its true value. Two things follow that a plain
mixture model cannot deliver:

* membership probabilities that are fair to faint stars, whose large error
  bars would otherwise push them out of a narrow clump;
* the cluster's intrinsic spread in motion -- a physical property, the
  internal velocity dispersion -- instead of that spread plus noise.

The model is a sum of Gaussians: two for the cluster and several for the
field, whose distribution of motions is lumpy and long-tailed. For star i and
component j, the chance of seeing measurement x_i is a Gaussian whose width is
the component's own spread V_j *plus* the star's error ellipse S_i.

Why two for the cluster: real clusters are not single Gaussians. They have a
dense core and a looser halo (or "corona") of members that share the cluster's
motion but spread several times further in space (Meingast et al. 2021, A&A
645, A84). Given only one cluster Gaussian, the model invents a "field"
component at the cluster's own motion and distance to hold the halo, and each
member's probability is then split between the two. Giving the cluster a core
and a halo, and counting membership as belonging to either, keeps the physics
and the bookkeeping aligned.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

LOG_2PI = np.log(2.0 * np.pi)


@dataclass
class Mixture:
    """A fitted (or initial) Gaussian mixture.

    ``cluster`` lists the components that describe the cluster; the first is
    its core, whose mean defines the cluster's motion and distance.
    """

    weights: np.ndarray  # (K,)
    means: np.ndarray    # (K, D)
    covs: np.ndarray     # (K, D, D) -- intrinsic, i.e. deconvolved
    cluster: tuple = (0,)
    log_likelihood: float = float("nan")
    n_iter: int = 0
    converged: bool = False

    @property
    def core(self) -> int:
        return self.cluster[0]

    @property
    def dim(self) -> int:
        return self.means.shape[1]

    def copy(self) -> "Mixture":
        return Mixture(self.weights.copy(), self.means.copy(), self.covs.copy(),
                       tuple(self.cluster), self.log_likelihood, self.n_iter, self.converged)


def _inverse_and_logdet(t: np.ndarray):
    """Batched inverse and log-determinant of symmetric positive-definite matrices.

    Closed forms for 2x2 and 3x3 -- the only sizes used here -- are several
    times faster than a general LAPACK call over ~10^5 tiny matrices.
    """
    d = t.shape[-1]
    if d == 2:
        a, b, c = t[:, 0, 0], t[:, 0, 1], t[:, 1, 1]
        det = a * c - b * b
        inv = np.empty_like(t)
        inv[:, 0, 0] = c / det
        inv[:, 1, 1] = a / det
        inv[:, 0, 1] = inv[:, 1, 0] = -b / det
        return inv, np.log(det)
    if d == 3:
        a, b, c = t[:, 0, 0], t[:, 0, 1], t[:, 0, 2]
        e, f, i = t[:, 1, 1], t[:, 1, 2], t[:, 2, 2]
        c00 = e * i - f * f
        c01 = c * f - b * i
        c02 = b * f - c * e
        det = a * c00 + b * c01 + c * c02
        inv = np.empty_like(t)
        inv[:, 0, 0] = c00 / det
        inv[:, 0, 1] = inv[:, 1, 0] = c01 / det
        inv[:, 0, 2] = inv[:, 2, 0] = c02 / det
        inv[:, 1, 1] = (a * i - c * c) / det
        inv[:, 1, 2] = inv[:, 2, 1] = (b * c - a * f) / det
        inv[:, 2, 2] = (a * e - b * b) / det
        return inv, np.log(det)
    sign, logdet = np.linalg.slogdet(t)
    return np.linalg.inv(t), logdet


def _component_terms(x, s, mean, cov):
    """Per-star log-density under one component, plus the pieces the M-step needs."""
    t_inv, logdet = _inverse_and_logdet(cov[None, :, :] + s)
    r = x - mean
    t_inv_r = np.einsum("nij,nj->ni", t_inv, r)
    maha = np.einsum("ni,ni->n", r, t_inv_r)
    logpdf = -0.5 * (x.shape[1] * LOG_2PI + logdet + maha)
    return logpdf, t_inv, t_inv_r


def _logsumexp(a, axis):
    m = np.max(a, axis=axis, keepdims=True)
    return (m + np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True))).squeeze(axis)


def responsibilities(x, s, mix: Mixture):
    """Probability that each star belongs to each component, and the log-likelihood."""
    k = len(mix.weights)
    logp = np.empty((len(x), k))
    for j in range(k):
        logp[:, j] = np.log(mix.weights[j]) + _component_terms(x, s, mix.means[j], mix.covs[j])[0]
    norm = _logsumexp(logp, axis=1)
    return np.exp(logp - norm[:, None]), float(norm.sum())


def fit(x, s, init: Mixture, *, max_iter: int = 400, tol: float = 1e-8,
        regularisation: float = 1e-4) -> Mixture:
    """Extreme-deconvolution EM (Bovy et al. 2011, their Eqs. 17-21).

    E-step: for each star and component, the responsibility q_ij, and the
    expected *true* value b_ij and its covariance B_ij given the measurement:

        T_ij = V_j + S_i
        b_ij = m_j + V_j T_ij^-1 (x_i - m_j)
        B_ij = V_j - V_j T_ij^-1 V_j

    M-step: weights, means and covariances from the q-weighted b and B. The
    covariance update includes B, which is exactly the part a naive mixture
    leaves out -- and why this recovers intrinsic widths.

    Only the q-weighted *sum* of the B_ij is ever needed, and because V_j is
    the same for every star it factors out:

        sum_i q_ij B_ij = Q_j V_j - V_j (sum_i q_ij T_ij^-1) V_j

    so no per-star matrix product is formed -- the difference between minutes
    and seconds on a 10^5-star field.

    Stops when the mean log-likelihood per star changes by less than ``tol``.
    """
    mix = init.copy()
    n, d = x.shape
    k = len(mix.weights)
    eye = np.eye(d)
    previous = -np.inf

    for iteration in range(1, max_iter + 1):
        logp = np.empty((n, k))
        cache = []
        for j in range(k):
            logpdf, t_inv, t_inv_r = _component_terms(x, s, mix.means[j], mix.covs[j])
            logp[:, j] = np.log(mix.weights[j]) + logpdf
            cache.append((t_inv, t_inv_r))
        norm = _logsumexp(logp, axis=1)
        q = np.exp(logp - norm[:, None])
        total = float(norm.sum())

        for j in range(k):
            t_inv, t_inv_r = cache[j]
            v = mix.covs[j]
            qj = q[:, j]
            qsum = qj.sum()
            if qsum < 1e-8:  # a component that has lost every star
                continue
            b = mix.means[j] + t_inv_r @ v  # V is symmetric, so (V T^-1 r)^T = r^T T^-1 V
            mean = (qj @ b) / qsum
            delta = b - mean
            summed_t_inv = np.tensordot(qj, t_inv, axes=1)
            cov = ((qj[:, None] * delta).T @ delta
                   + qsum * v - v @ summed_t_inv @ v) / qsum
            mix.means[j] = mean
            mix.covs[j] = 0.5 * (cov + cov.T) + regularisation * eye
            mix.weights[j] = qsum / n

        mix.weights /= mix.weights.sum()
        mix.log_likelihood = total
        mix.n_iter = iteration
        if abs(total - previous) / n < tol:
            mix.converged = True
            break
        previous = total
    return mix


def membership_probability(x, s, mix: Mixture, *, core_only: bool = False) -> np.ndarray:
    """P(cluster | measurement) for every star: core plus halo, or core alone."""
    q, _ = responsibilities(x, s, mix)
    return q[:, [mix.core] if core_only else list(mix.cluster)].sum(axis=1)


def posterior_mean(x, s, mix: Mixture, component: int | None = None) -> np.ndarray:
    """Best estimate of each star's *true* value, given that it is a member.

    For a cluster member this combines the star's own measurement with the
    knowledge that it shares the cluster's motion and distance, each weighted
    by its precision. A bright star with a superb parallax keeps essentially
    its own value; a faint star with a noisy one leans on the cluster.

    With a core and a halo, the estimate from each is weighted by how likely
    the star is to belong to that part.
    """
    parts = [component] if component is not None else list(mix.cluster)
    q, _ = responsibilities(x, s, mix)
    w = q[:, parts]
    w = w / np.maximum(w.sum(axis=1, keepdims=True), 1e-300)
    out = np.zeros_like(x, dtype=float)
    for k, j in enumerate(parts):
        _, _, t_inv_r = _component_terms(x, s, mix.means[j], mix.covs[j])
        out += w[:, [k]] * (mix.means[j] + t_inv_r @ mix.covs[j])
    return out


def mean_uncertainty(s, mix: Mixture, weights) -> np.ndarray:
    """Approximate covariance of the cluster core's mean, from the summed precision.

    Each member constrains the mean with precision (V + S_i)^-1 -- the cluster's
    own spread plus its measurement error -- weighted by how surely it is a
    member. Adding those up and inverting gives the uncertainty of the mean.
    """
    t_inv, _ = _inverse_and_logdet(mix.covs[mix.core][None, :, :] + s)
    return np.linalg.inv(np.einsum("n,nij->ij", weights, t_inv))


def cluster_components(mix: Mixture, *, depth_axis: int = 2, max_depth_fraction: float = 0.15,
                       max_offset: float = 1.0) -> tuple:
    """Which components describe the cluster: the core, plus any halo-like ones.

    Labels given at the start of a fit are not reliable -- a component seeded
    as the halo can wander off, while one seeded on the field settles onto
    the halo. So components are classified afterwards, by what they are:

    * **centred** on the core: its mean within ``max_offset`` of its own
      standard deviations of the core's, in motion *and* parallax;
    * **compact in depth**: its parallax spread below ``max_depth_fraction``
      of the cluster's parallax -- i.e. +/-15% in distance, a generous
      corona. Broad field components in the same patch of motion space span
      25-50%, because they mix stars from ~70 pc to kiloparsecs.

    Needs the parallax dimension; a motion-only model keeps just its core.
    """
    core = mix.core
    if mix.dim <= depth_axis:
        return (core,)
    plx_core = mix.means[core][depth_axis]
    members = [core]
    for j in range(len(mix.weights)):
        if j == core:
            continue
        sd = np.sqrt(np.diag(mix.covs[j]))
        centred = np.all(np.abs(mix.means[j] - mix.means[core]) < max_offset * sd)
        compact = sd[depth_axis] < max_depth_fraction * abs(plx_core)
        if centred and compact:
            members.append(j)
    return tuple(members)


def initialise(x, seed_mask, n_field: int, *, halo: bool = True,
               random_state: int = 0) -> Mixture:
    """Starting point: the cluster from the seed stars, the field from the rest.

    The seed is the knot found in proper-motion space -- the stars within a few
    mas/yr of it. The core starts as their robust spread; the halo, if used,
    starts at the same centre and three times wider. The field components come
    from an ordinary Gaussian mixture fitted to every other star. Extreme
    deconvolution then refines everything jointly, including which stars
    really belong where.
    """
    from sklearn.mixture import GaussianMixture

    seed = x[seed_mask]
    centre = np.median(seed, axis=0)
    spread = np.maximum(1.4826 * np.median(np.abs(seed - centre), axis=0), 1e-2)

    gmm = GaussianMixture(n_components=n_field, covariance_type="full",
                          random_state=random_state, max_iter=300)
    gmm.fit(x[~seed_mask])

    w_cluster = seed_mask.mean()
    if halo:
        c_weights = [0.7 * w_cluster, 0.3 * w_cluster]
        c_means = [centre, centre]
        c_covs = [np.diag(spread ** 2), np.diag((3 * spread) ** 2)]
    else:
        c_weights, c_means, c_covs = [w_cluster], [centre], [np.diag(spread ** 2)]
    n_c = len(c_weights)
    weights = np.concatenate([c_weights, gmm.weights_ * (1.0 - w_cluster)])
    means = np.vstack([np.array(c_means), gmm.means_])
    covs = np.concatenate([np.array(c_covs), gmm.covariances_])
    return Mixture(weights=weights, means=means, covs=covs, cluster=tuple(range(n_c)))
