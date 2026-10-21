"""Turning SALT parameters into distances: the Tripp relation.

    mu = mB + alpha * x1 - beta * c - M          (Tripp 1998)

Slower supernovae (larger x1) are intrinsically brighter, so they are made
fainter by alpha * x1; redder ones (larger c) look dimmer, so they are
brightened by beta * c. M is the absolute magnitude of a standard supernova.

"D1's alpha, beta and M". D1 fitted Pantheon+'s standardised magnitudes
m_b_corr, which the Pantheon+ team built from raw SALT2 parameters. The
recipe is recoverable exactly from D1's own table: a least-squares fit gives

    m_b_corr = mB + 0.14753 x1 - 3.09222 c - biasCor_m_b - 0.0524

with residuals of 4e-5 mag (rounding), i.e. alpha = 0.148 and beta = 3.09,
as in Brout et al. (2022, Table 2), plus a simulation-based bias correction
(biasCor_m_b) and a constant. D1's M then turns m_b_corr into mu. At low
redshift that bias correction is mostly not about selection: in Pantheon+'s
"BS21" dust model it carries a colour term and a host-mass step, so applied
to raw SALT parameters D1's recipe behaves like beta ~ 2.7 plus a ~0.07 mag
mass step. We test both readings: the nominal Tripp relation, and the full
recipe with an emulated bias correction.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import optimize

NMAD_SCALE = 1.4826  # median absolute deviation -> sigma for a Gaussian


def nmad(x) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(NMAD_SCALE * np.median(np.abs(x - np.median(x))))


def chauvenet_mask(x, n_sigma: float | None = None) -> np.ndarray:
    """Keep points within n_sigma (robust) of the median.

    Chauvenet's criterion: reject a point if fewer than half a point is
    expected that far out among N Gaussian draws. For N = 2,600 that is ~4.1
    sigma, the cut Rigault et al. (2025) applied. Sigma is the nMAD, so the
    cut is not dragged around by the outliers it is looking for.
    """
    from scipy import stats

    x = np.asarray(x, float)
    if n_sigma is None:
        n_sigma = float(stats.norm.isf(0.25 / np.isfinite(x).sum()))
    s = nmad(x)
    return np.abs(x - np.nanmedian(x)) < n_sigma * s


# ---------------------------------------------------------------------------
# D1's recipe, recovered from the Pantheon+ table
# ---------------------------------------------------------------------------

@dataclass
class D1Recipe:
    alpha: float
    beta: float
    const: float               # constant inside m_b_corr (absorbed by M)
    regression_rms: float      # how exactly the recipe reproduces m_b_corr
    bias_coef: dict = field(default_factory=dict)  # emulated low-z biasCor
    bias_rms: float = np.nan
    bias_n: int = 0

    def bias_correction(self, x1, c, host_mass, z) -> np.ndarray:
        """Emulated Pantheon+ bias correction for a low-z supernova.

        A linear fit to Pantheon+'s biasCor_m_b at 0.01 < z < 0.1 in colour,
        stretch, a host-mass step at 10^10 Msun, and redshift. Unknown host
        mass counts as half way (the step's midpoint).
        """
        b = self.bias_coef
        hm = np.where(np.isfinite(host_mass), (np.asarray(host_mass) > 10).astype(float), 0.5)
        return (b["const"] + b["c"] * np.asarray(c) + b["x1"] * np.asarray(x1)
                + b["high_mass"] * hm + b["z"] * np.asarray(z))


def recover_d1_recipe(pantheon: pd.DataFrame, z_range=(0.01, 0.1)) -> D1Recipe:
    """Read alpha, beta and the bias-correction behaviour out of Pantheon+."""
    p = pantheon
    y = p["m_b_corr"] - p["mB"] + p["biasCor_m_b"]
    A = np.column_stack([p["x1"], p["c"], np.ones(len(p))])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    rms = float(np.std(y - A @ coef))
    recipe = D1Recipe(alpha=float(coef[0]), beta=float(-coef[1]), const=float(coef[2]),
                      regression_rms=rms)

    lo = p[(p["zHD"] > z_range[0]) & (p["zHD"] < z_range[1])]
    hm = (lo["HOST_LOGMASS"] > 10).astype(float)
    B = np.column_stack([lo["c"], lo["x1"], hm, lo["zHD"], np.ones(len(lo))])
    bc, *_ = np.linalg.lstsq(B, lo["biasCor_m_b"], rcond=None)
    recipe.bias_coef = dict(zip(["c", "x1", "high_mass", "z", "const"], map(float, bc)))
    recipe.bias_rms = float(np.std(lo["biasCor_m_b"] - B @ bc))
    recipe.bias_n = int(len(lo))
    return recipe


# ---------------------------------------------------------------------------
# Distance moduli and their errors
# ---------------------------------------------------------------------------

def tripp(fits: pd.DataFrame, alpha: float, beta: float, M: float,
          offset=0.0) -> np.ndarray:
    """mu = mB + alpha x1 - beta c - M - offset (offset: any extra per-SN term)."""
    return fits["mB"] + alpha * fits["x1"] - beta * fits["c"] - M - offset


def tripp_error(fits: pd.DataFrame, alpha: float, beta: float) -> np.ndarray:
    """Measurement error on mu from the SALT fit's covariance.

    var = var(mB) + a^2 var(x1) + b^2 var(c) + 2a cov(mB,x1) - 2b cov(mB,c)
          - 2ab cov(x1,c)
    """
    v = (fits["mB_err"] ** 2 + alpha ** 2 * fits["x1_err"] ** 2 + beta ** 2 * fits["c_err"] ** 2
         + 2 * alpha * fits["cov_mB_x1"] - 2 * beta * fits["cov_mB_c"]
         - 2 * alpha * beta * fits["cov_x1_c"])
    return np.sqrt(np.clip(v, 0, None))


# ---------------------------------------------------------------------------
# Fitting the standardisation on ZTF itself
# ---------------------------------------------------------------------------

@dataclass
class TrippFit:
    params: dict
    errors: dict
    sigma_int: float
    chi2: float
    n: int
    residuals: np.ndarray
    names: list
    n_clipped: int = 0


def fit_tripp(fits: pd.DataFrame, mu_model: np.ndarray, extra_var: np.ndarray,
              host_mass=None, mass_split: float = 10.0, n_iter: int = 10,
              clip: float | None = None) -> TrippFit:
    """Fit alpha, beta (and a host-mass step gamma) with an intrinsic scatter.

    Same scheme as Ginolin et al. (2025a): at a fixed intrinsic scatter
    sigma_int, minimise chi2 = sum (mB + a x1 - b c - M [- gamma step] - mu)^2 /
    (sigma_mu(a, b)^2 + sigma_int^2 + extra), then reset sigma_int so that
    chi2 per degree of freedom is 1, and repeat. M is a free offset (it
    absorbs the blinded ZTF zero point and H0). Errors come from the inverse
    Hessian of chi2 at the final sigma_int. With `clip`, supernovae further
    than clip x their total error from the model are dropped at each
    iteration (Hubble residuals have heavier tails than a Gaussian).
    """
    use_step = host_mass is not None
    names = ["alpha", "beta", "M"] + (["gamma"] if use_step else [])
    if use_step:
        hm = np.asarray(host_mass, float)
        step = np.where(hm > mass_split, 0.5, -0.5)  # +gamma/2 high mass, -gamma/2 low
    y_obs = np.asarray(fits["mB"], float) - np.asarray(mu_model, float)

    def resid(p):
        a, b, M = p[:3]
        r = y_obs + a * fits["x1"].to_numpy() - b * fits["c"].to_numpy() - M
        if use_step:
            r = r + p[3] * step
        return r

    extra_var = np.asarray(extra_var, float)
    keep = np.ones(len(y_obs), bool)

    def pulls(p, s_int):
        return resid(p) / np.sqrt(tripp_error(fits, p[0], p[1]) ** 2 + s_int ** 2 + extra_var)

    def chi2(p, s_int):
        return float(np.sum(pulls(p, s_int)[keep] ** 2))

    p = np.array([0.15, 3.0, float(np.median(y_obs))] + ([0.1] if use_step else []))
    s_int = 0.12
    for _ in range(n_iter):
        dof = int(keep.sum()) - len(p)
        p = optimize.minimize(chi2, p, args=(s_int,), method="Nelder-Mead",
                              options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 20000}).x
        p = optimize.minimize(chi2, p, args=(s_int,), method="BFGS").x
        # sigma_int that makes chi2/dof = 1
        f = lambda s: chi2(p, s) / dof - 1
        if f(0.0) <= 0:
            s_int = 0.0
        else:
            s_int = optimize.brentq(f, 0.0, 2.0)
        if clip is not None:
            keep = np.abs(np.asarray(pulls(p, s_int))) < clip
    H = _hessian(lambda q: chi2(q, s_int), p)
    cov = 2 * np.linalg.inv(H)
    errs = np.sqrt(np.diag(cov))
    return TrippFit(params=dict(zip(names, map(float, p))), errors=dict(zip(names, map(float, errs))),
                    sigma_int=float(s_int), chi2=chi2(p, s_int), n=int(keep.sum()),
                    residuals=resid(p), names=names, n_clipped=int((~keep).sum()))


def _hessian(f, x, rel: float = 1e-4) -> np.ndarray:
    x = np.asarray(x, float)
    n = len(x)
    h = rel * np.maximum(np.abs(x), 1e-2)
    H = np.empty((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ei[i] = h[i]
            ej = np.zeros(n); ej[j] = h[j]
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej)
                                 + f(x - ei - ej)) / (4 * h[i] * h[j])
    return H
