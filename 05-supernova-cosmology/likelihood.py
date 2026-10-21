"""Scoring a universe against the supernovae: chi-squared with the full covariance.

The data are standardised peak magnitudes m_i. A universe predicts

    m_i = M + mu_i,    mu_i = 5 log10(d_L(z_i) / Mpc) + 25

where M is the absolute magnitude every standardised Type Ia supernova shares.
The fit quality is chi^2 = r^T C^-1 r with r the vector of data-minus-model and C
the 1701 x 1701 covariance published with the data. C is not diagonal: supernovae
observed with the same telescope share its calibration errors, supernovae at
similar redshift share selection effects, and so on. Treating the errors as
independent would be wrong, and the data release says so in capitals.

**The degeneracy that shapes everything.** d_L is proportional to c/H0, so

    m_i = [M + 5 log10(c/H0) + 25] + 5 log10(d_L,i in units of c/H0)

and M and H0 enter only through the bracket. Supernovae alone measure the
*shape* of the Hubble diagram -- how distance grows with redshift, which is what
tells matter from dark energy -- but not its absolute scale. Getting H0 needs
the supernova luminosity M pinned down independently, which is what the Cepheid
distances to the calibrator galaxies do (``fit_calibrated``).

Both kinds of nuisance parameter enter linearly, so they are solved exactly by
generalised least squares at each point of a grid of universes, rather than
sampled. With flat priors, solving for them this way ("profiling") and
integrating over them ("marginalising") give the same answer here, because the
Gaussian's width in those directions is identical for every universe.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cholesky, solve_triangular

import cosmology as cosmo
from config import C_KM_S, Z_MIN_COSMOLOGY


@dataclass
class Sample:
    """A subset of light curves, prepared for fast repeated chi^2 evaluation."""

    name: str
    rows: np.ndarray  # indices into the full 1701-row table
    z_hd: np.ndarray
    z_hel: np.ndarray
    m: np.ndarray  # standardised peak magnitude m_b_corr
    cov: np.ndarray
    is_calibrator: np.ndarray
    ceph_mu: np.ndarray  # Cepheid distance modulus (calibrators only; NaN otherwise)
    whiten: np.ndarray  # W = L^-1 with C = L L^T, so chi^2 = |W r|^2
    nodes: np.ndarray  # distance-integration nodes (Hubble-flow redshifts)
    index: np.ndarray  # node index of each Hubble-flow supernova

    @property
    def n(self) -> int:
        return self.m.size

    @property
    def hubble_flow(self) -> np.ndarray:
        return ~self.is_calibrator


def build_sample(table, cov, *, name: str, calibrators: bool = False,
                 diagonal_only: bool = False, rows=None) -> Sample:
    """Select light curves exactly as the published likelihoods do.

    Supernova-only fits keep z_HD > 0.01. Calibrated fits also keep every
    supernova in a Cepheid host, at any redshift, and for those the model
    distance is the Cepheid one rather than the cosmological one.
    ``rows`` overrides the selection (used for the SH0ES Hubble-flow cross-check).
    """
    z_hd = table["zHD"].to_numpy()
    is_cal = table["IS_CALIBRATOR"].to_numpy().astype(bool)
    if rows is None:
        keep = z_hd > Z_MIN_COSMOLOGY
        if calibrators:
            keep |= is_cal
        rows = np.flatnonzero(keep)
    rows = np.asarray(rows)

    sub = cov[np.ix_(rows, rows)]
    if diagonal_only:
        sub = np.diag(np.diag(sub))
    lower = cholesky(sub, lower=True)
    whiten = solve_triangular(lower, np.eye(rows.size), lower=True)

    cal = is_cal[rows] if calibrators else np.zeros(rows.size, bool)
    ceph = np.where(cal, table["CEPH_DIST"].to_numpy()[rows], np.nan)
    zh = z_hd[rows]
    nodes, index = cosmo.integration_nodes(zh[~cal])
    return Sample(name=name, rows=rows, z_hd=zh, z_hel=table["zHEL"].to_numpy()[rows],
                  m=table["m_b_corr"].to_numpy()[rows], cov=sub, is_calibrator=cal,
                  ceph_mu=ceph, whiten=whiten, nodes=nodes, index=index)


def shape_magnitudes(sample: Sample, Om, OL, w=-1.0) -> np.ndarray:
    """5 log10(d_L / (c/H0)) for each Hubble-flow supernova; shape (n_models, n_hf)."""
    hf = sample.hubble_flow
    dl = cosmo.luminosity_distance_shape(sample.z_hd[hf], sample.z_hel[hf], Om, OL, w,
                                         nodes=sample.nodes, index=sample.index)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 5.0 * np.log10(dl)


# ---------------------------------------------------------------------------
# Supernovae alone: profile over the one degenerate offset
# ---------------------------------------------------------------------------


def chi2_profile_offset(sample: Sample, shape: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """chi^2 minimised over the offset A = M + 5 log10(c/H0) + 25, per model.

    For residuals d = m - shape - A, the best A and the chi^2 there are closed
    form: A = (1^T C^-1 d)/(1^T C^-1 1), chi^2 = d^T C^-1 d - (1^T C^-1 d)^2/(1^T C^-1 1).
    Computed in whitened space so one matrix product handles a whole batch.
    Invalid models (NaN shapes) get chi^2 = inf.
    """
    if sample.is_calibrator.any():
        raise ValueError("supernova-only chi^2 needs a sample without calibrators")
    shape = np.atleast_2d(shape)
    bad = ~np.all(np.isfinite(shape), axis=1)
    safe = np.where(bad[:, None], 0.0, shape)
    u = sample.whiten @ np.ones(sample.n)
    r = (sample.whiten @ sample.m)[:, None] - sample.whiten @ safe.T  # (n, k)
    uu = u @ u
    ur = u @ r
    chi2 = np.einsum("ij,ij->j", r, r) - ur**2 / uu
    offset = ur / uu
    chi2[bad] = np.inf
    offset[bad] = np.nan
    return chi2, offset


def scan(sample: Sample, Om, OL, w=-1.0, *, chunk: int = 2000):
    """chi^2 and best offset for many universes (flattened parameter arrays)."""
    Om, OL, w = np.broadcast_arrays(np.asarray(Om, float), np.asarray(OL, float),
                                    np.asarray(w, float))
    shape_in = Om.shape
    Om, OL, w = Om.ravel(), OL.ravel(), w.ravel()
    chi2 = np.empty(Om.size)
    offset = np.empty(Om.size)
    for start in range(0, Om.size, chunk):
        sl = slice(start, start + chunk)
        s = shape_magnitudes(sample, Om[sl], OL[sl], w[sl])
        chi2[sl], offset[sl] = chi2_profile_offset(sample, s)
    return chi2.reshape(shape_in), offset.reshape(shape_in)


# ---------------------------------------------------------------------------
# With the Cepheid calibrators: M and H0 separately
# ---------------------------------------------------------------------------


@dataclass
class CalibratedFit:
    Om: np.ndarray  # grid
    chi2: np.ndarray  # chi^2 at the best (M, 5 log10 H0) for each Om
    M: np.ndarray  # best-fit supernova absolute magnitude at each Om
    h5: np.ndarray  # best-fit 5 log10(H0 / km s^-1 Mpc^-1) at each Om
    cov_beta: np.ndarray  # 2x2 covariance of (M, h5); the same for every Om
    n: int


def fit_calibrated(sample: Sample, Om_grid, OL_grid=None, *, dl_shape_fn=None) -> CalibratedFit:
    """Fit (M, H0) jointly with the calibrators, on a grid of Omega_m (flat LCDM).

    Calibrator rows:   m_i - mu_Cepheid,i                   = M
    Hubble-flow rows:  m_i - 5 log10(c d_L,i[c/H0]) - 25    = M - 5 log10 H0

    Linear in beta = (M, 5 log10 H0), so beta is solved exactly by generalised
    least squares. ``dl_shape_fn`` replaces the cosmological distance with any
    other function of the sample returning d_L in c/H0 units (used for the
    kinematic cross-check); ``Om_grid`` is then a dummy one-element grid.
    """
    Om_grid = np.atleast_1d(np.asarray(Om_grid, float))
    OL_grid = 1.0 - Om_grid if OL_grid is None else np.atleast_1d(OL_grid)
    cal, hf = sample.is_calibrator, sample.hubble_flow
    X = np.zeros((sample.n, 2))
    X[:, 0] = 1.0
    X[hf, 1] = -1.0
    Xw = sample.whiten @ X
    F = Xw.T @ Xw
    cov_beta = np.linalg.inv(F)

    if dl_shape_fn is None:
        shape = shape_magnitudes(sample, Om_grid, OL_grid)  # (k, n_hf)
    else:
        shape = np.atleast_2d(5.0 * np.log10(dl_shape_fn(sample)))
    y = np.repeat(sample.m[None, :], shape.shape[0], axis=0)
    y[:, cal] -= sample.ceph_mu[cal]
    y[:, hf] -= shape + 5.0 * np.log10(C_KM_S) + 25.0
    yw = sample.whiten @ y.T  # (n, k)
    beta = cov_beta @ (Xw.T @ yw)  # (2, k)
    resid = yw - Xw @ beta
    chi2 = np.einsum("ij,ij->j", resid, resid)
    return CalibratedFit(Om=Om_grid, chi2=chi2, M=beta[0], h5=beta[1], cov_beta=cov_beta,
                         n=sample.n)


def kinematic_dl(sample: Sample, q0: float = -0.55, j0: float = 1.0) -> np.ndarray:
    """Model-independent d_L (c/H0 units) from a Taylor expansion in z.

    d_L = z [1 + (1 - q0) z / 2 - (1 - q0 - 3 q0^2 + j0) z^2 / 6], as used by SH0ES
    (Riess et al. 2016, 2022) with q0 = -0.55, j0 = 1. Only valid at low z, which
    is the point: it measures H0 without assuming any particular cosmology.
    """
    hf = sample.hubble_flow
    z = sample.z_hd[hf]
    series = z * (1.0 + 0.5 * (1.0 - q0) * z - (1.0 - q0 - 3.0 * q0**2 + j0) * z**2 / 6.0)
    return series * (1.0 + sample.z_hel[hf]) / (1.0 + z)


def h0_posterior(fit: CalibratedFit, h0_grid) -> tuple[np.ndarray, np.ndarray]:
    """Posterior density of H0 on a grid, marginalised over Omega_m and M.

    At fixed Omega_m, 5 log10 H0 is Gaussian around the GLS solution with variance
    cov_beta[1, 1]; integrating over Omega_m mixes those Gaussians with weights
    exp(-chi^2/2). Flat prior in H0, as in the published chains.
    Returns (density, weights over the Omega_m grid).
    """
    wts = np.exp(-0.5 * (fit.chi2 - fit.chi2.min()))
    wts /= wts.sum()
    h5_of_grid = 5.0 * np.log10(np.asarray(h0_grid, float))
    sig = np.sqrt(fit.cov_beta[1, 1])
    dens = (wts[:, None] * np.exp(-0.5 * ((h5_of_grid[None, :] - fit.h5[:, None]) / sig) ** 2)).sum(0)
    dens /= np.trapezoid(dens, h0_grid)
    return dens, wts


# ---------------------------------------------------------------------------
# Turning chi^2 grids into probabilities
# ---------------------------------------------------------------------------

#: Probability contained within 1, 2 and 3 sigma of a Gaussian.
SIGMA_MASS = (0.682689, 0.954500, 0.997300)


def grid_posterior(chi2: np.ndarray, prior_mask=None) -> np.ndarray:
    """Normalised posterior (sums to 1) from chi^2 on a uniform grid, flat prior."""
    p = np.exp(-0.5 * (chi2 - np.min(chi2[np.isfinite(chi2)])))
    p[~np.isfinite(p)] = 0.0
    if prior_mask is not None:
        p = np.where(prior_mask, p, 0.0)
    return p / p.sum()


def credible_levels(p: np.ndarray, masses=SIGMA_MASS) -> list[float]:
    """Density thresholds whose superlevel sets hold the given probability mass.

    These are highest-posterior-density regions: the smallest region containing
    68.3% (95.4%, 99.7%) of the probability -- what "1, 2, 3 sigma contours" mean.
    """
    flat = np.sort(p.ravel())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return [float(flat[min(np.searchsorted(cum, m), flat.size - 1)]) for m in masses]


def summarize_1d(x: np.ndarray, p: np.ndarray) -> dict:
    """Mean, standard deviation and central 68% interval of a 1D grid density."""
    p = np.asarray(p, float) / np.sum(p)
    mean = float(np.sum(x * p))
    std = float(np.sqrt(np.sum((x - mean) ** 2 * p)))
    cdf = np.cumsum(p)
    lo, med, hi = np.interp([0.158655, 0.5, 0.841345], cdf, x)
    return {"mean": mean, "std": std, "median": float(med),
            "minus": float(med - lo), "plus": float(hi - med)}


def summarize_samples(x: np.ndarray) -> dict:
    lo, med, hi = np.percentile(x, [15.8655, 50, 84.1345])
    return {"mean": float(np.mean(x)), "std": float(np.std(x)), "median": float(med),
            "minus": float(med - lo), "plus": float(hi - med)}


def sigma_from_delta_chi2(dchi2: float) -> float:
    """Equivalent Gaussian significance of a one-parameter nested-model comparison."""
    return float(np.sqrt(max(dchi2, 0.0)))


# ---------------------------------------------------------------------------
# Binning -- for display only; every fit uses the unbinned data
# ---------------------------------------------------------------------------


def binned_residuals(sample: Sample, resid: np.ndarray, edges) -> dict:
    """Compress residuals into redshift bins with the full covariance.

    Fits a constant per bin by generalised least squares -- the same estimator as
    the cosmology fits -- so each binned point is an honest summary of its
    supernovae, correlations included. Returns bin redshifts, values, and the
    bins' covariance. Empty bins are dropped.

    The redshift of each bin is an inverse-variance weighted mean. It is *not*
    the same GLS projection applied to z: with correlated errors some GLS weights
    are negative, which can put a bin's "redshift" outside the bin itself.
    """
    z = sample.z_hd
    which = np.digitize(z, edges) - 1
    used = [k for k in range(len(edges) - 1) if np.any(which == k)]
    B = np.zeros((sample.n, len(used)))
    for j, k in enumerate(used):
        B[which == k, j] = 1.0
    Bw = sample.whiten @ B
    cov_b = np.linalg.inv(Bw.T @ Bw)
    proj = cov_b @ Bw.T @ sample.whiten
    ivar = 1.0 / np.diag(sample.cov)
    return {"z": (B.T @ (ivar * z)) / (B.T @ ivar), "value": proj @ resid, "cov": cov_b,
            "err": np.sqrt(np.diag(cov_b)), "count": B.sum(0).astype(int)}
