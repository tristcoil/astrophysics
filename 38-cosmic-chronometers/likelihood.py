"""The chronometer data set, its covariance, and the likelihood fits.

**Data.** The 32 measurements of Moresco et al. (2022, Table 1). For the 15 that
come from Moresco's own papers, the table's rounded values are replaced by the
precise ones from his GitLab primer (matched by redshift and H; e.g. the table's
z = 0.445 is the measurement at z = 0.4497).

**Covariance.** Following Moresco et al. (2020) and the primer notebook,

    C = diag(sigma_i^2) + sum_terms (H_i f(z_i)) (H_j f(z_j)),

where sigma_i is the tabulated error (statistical, metallicity, star-formation
history and young-star contamination; the 2022 review says these are already
in it) and f(z) is the fractional systematic error from the choice of initial
mass function and of stellar-population-synthesis model. Those two are
treated as fully correlated between all redshifts, because the same model
would shift every measurement the same way.

**Fits.** In every model H(z) = H0 E(z; theta), so the chi^2 is an exact
quadratic in H0: chi^2 = a H0^2 - 2 b H0 + c, with a = E.C^-1.E, b = E.C^-1.H,
c = H.C^-1.H. H0 can therefore be integrated out analytically -- over a flat
prior box, or against a Gaussian prior -- leaving a grid over theta alone.
Posterior samples (for derived quantities such as the transition redshift)
draw theta from the grid and then H0 from its exact conditional distribution.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import linalg, special, stats

import config as cfg
import cosmology as cosmo
import fetch_data


# ---------------------------------------------------------------------------
# The data set
# ---------------------------------------------------------------------------


def compile_chronometers() -> dict[str, np.ndarray]:
    """The 32 chronometer H(z) values, sorted by redshift.

    Keys: z, H, err (diagonal error used), method (D/F/L), ref, is_moresco,
    and the 2022 table's own z_table, H_table, err_table for comparison.
    """
    lrr = fetch_data.load_lrr_compilation()
    mm = fetch_data.load_moresco_table("BC03")
    z, H, err = lrr["z"].copy(), lrr["H"].copy(), lrr["err"].copy()
    is_mm = np.zeros(z.size, dtype=bool)
    used = []
    for i in np.nonzero(lrr["method"] == "D")[0]:
        hit = np.nonzero((np.abs(mm["z"] - z[i]) < cfg.MATCH_DZ)
                         & (np.abs(mm["H"] - H[i]) < cfg.MATCH_DH))[0]
        if hit.size != 1:
            raise ValueError(f"table row z={z[i]}, H={H[i]} matches {hit.size} GitLab rows")
        j = hit[0]
        used.append(j)
        z[i], H[i], err[i] = mm["z"][j], mm["H"][j], mm["err"][j]
        is_mm[i] = True
    if sorted(used) != list(range(mm["z"].size)):
        raise ValueError("Moresco's 15 measurements were not matched one-to-one")
    order = np.argsort(z, kind="stable")
    out = {"z": z, "H": H, "err": err, "method": lrr["method"], "ref": lrr["ref"],
           "is_moresco": is_mm, "z_table": lrr["z"], "H_table": lrr["H"],
           "err_table": lrr["err"]}
    return {k: v[order] for k, v in out.items()}


def chronometer_covariance(z, H, err, *, apply_to=None, terms=cfg.SYSTEMATIC_TERMS,
                           syst=None) -> np.ndarray:
    """Moresco et al. (2020) covariance: diagonal errors plus fully correlated
    systematic terms (per cent of H, interpolated in redshift and held constant
    beyond the tabulated range 0.075-1.475, as in the primer notebook).

    The per-cent errors are turned into km/s/Mpc by multiplying with ``H``. The
    published recipe passes the measured values; ``smooth_scaled_covariance``
    passes a smooth reconstruction instead (see there for why it matters).
    ``apply_to`` (boolean mask) limits the systematic terms to some points.
    """
    syst = fetch_data.load_systematics() if syst is None else syst
    z, H, err = (np.asarray(v, dtype=float) for v in (z, H, err))
    cov = np.diag(err**2)
    mask = np.ones(z.size) if apply_to is None else np.asarray(apply_to, dtype=float)
    for term in terms:
        f = np.interp(z, syst["z"], syst[term]) / 100.0 * mask
        cov = cov + np.outer(H * f, H * f)
    return cov


def smooth_scaled_covariance(z, H, err, *, terms=cfg.SYSTEMATIC_TERMS, syst=None,
                             n_iter=6) -> tuple[np.ndarray, np.ndarray]:
    """The same covariance, but with the shared systematic scaled by a smooth
    reconstruction of H(z) instead of by each measured value.

    Why: a point that scattered high would otherwise get a larger correlated
    error than one that scattered low, so a fit leans on the low points and is
    pulled below the truth -- "Peelle's pertinent puzzle" of nuclear-data
    fitting (d'Agostini 1994). For n points sharing a fractional error f the
    pull is a factor of about 1 / (1 + n f^2): ~7% for 32 points at 5%.
    The remedy is to scale by the best estimate of the true values. Here that
    is the Gaussian-process mean (no cosmological model), iterated until it
    stops changing. Returns (covariance, the smooth H used for scaling).
    """
    import gp  # local import: gp imports nothing from here, but keep the layers apart

    syst = fetch_data.load_systematics() if syst is None else syst
    scale = np.asarray(H, dtype=float)
    for _ in range(n_iter):
        cov = chronometer_covariance(z, scale, err, terms=terms, syst=syst)
        g = gp.GaussianProcess(z, H, cov)
        new, _ = g.predict(np.asarray(z, float), g.optimize())
        if np.max(np.abs(new - scale) / scale) < 1e-4:
            scale = new
            break
        scale = new
    return chronometer_covariance(z, scale, err, terms=terms, syst=syst), scale


# ---------------------------------------------------------------------------
# Posterior summaries
# ---------------------------------------------------------------------------


def summarize(samples, weights=None) -> dict:
    """Mean, std, median, 68% and 95% intervals of a 1-D sample."""
    x = np.asarray(samples, dtype=float)
    keep = np.isfinite(x)
    x = x[keep]
    if weights is not None:
        w = np.asarray(weights, dtype=float)[keep]
        order = np.argsort(x)
        x, w = x[order], w[order]
        cdf = np.cumsum(w) / w.sum()
        q = lambda p: float(np.interp(p, cdf, x))
        mean = float(np.average(x, weights=w))
        std = float(np.sqrt(np.average((x - mean) ** 2, weights=w)))
    else:
        q = lambda p: float(np.quantile(x, p))
        mean, std = float(np.mean(x)), float(np.std(x))
    med = q(0.5)
    return {"mean": mean, "std": std, "median": med,
            "minus": med - q(0.15865), "plus": q(0.84135) - med,
            "lo95": q(0.02275), "hi95": q(0.97725), "n": int(x.size)}


def grid_summary(axis, density) -> dict:
    """Summary of a 1-D posterior tabulated on a grid."""
    axis = np.asarray(axis, dtype=float)
    p = np.clip(np.asarray(density, dtype=float), 0, None)
    # cumulative probability with the trapezoid rule
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(axis))])
    cdf /= cdf[-1]
    q = lambda prob: float(np.interp(prob, cdf, axis))
    norm = np.trapezoid(p, axis)
    mean = float(np.trapezoid(axis * p, axis) / norm)
    std = float(np.sqrt(np.trapezoid((axis - mean) ** 2 * p, axis) / norm))
    med = q(0.5)
    return {"mean": mean, "std": std, "median": med, "minus": med - q(0.15865),
            "plus": q(0.84135) - med, "lo95": q(0.02275), "hi95": q(0.97725),
            "mode": float(axis[np.argmax(p)])}


# ---------------------------------------------------------------------------
# Grid fits with H0 integrated out analytically
# ---------------------------------------------------------------------------


@dataclass
class GridFit:
    """Posterior on a grid of theta, with H0 marginalised analytically."""

    names: tuple                    # names of the theta axes
    axes: tuple                     # 1-D arrays, one per theta axis
    ln_like: np.ndarray             # ln of the H0-marginal likelihood on the grid
    mu: np.ndarray                  # conditional mean of H0 at each theta
    sigma: np.ndarray               # conditional std of H0 at each theta
    chi2_profile: np.ndarray        # min over H0 (inside the prior box) of chi^2
    h0_box: tuple                   # (lo, hi) of the flat H0 prior, or (-inf, inf)
    n_data: int
    extra: dict = field(default_factory=dict)

    @property
    def weights(self) -> np.ndarray:
        w = np.exp(self.ln_like - np.nanmax(self.ln_like))
        return np.where(np.isfinite(w), w, 0.0)

    def best(self) -> dict:
        i = np.unravel_index(np.nanargmin(self.chi2_profile), self.chi2_profile.shape)
        out = {n: float(ax[k]) for n, ax, k in zip(self.names, self.axes, i)}
        out["H0"] = float(np.clip(self.mu[i], *self.h0_box))
        out["chi2"] = float(self.chi2_profile[i])
        out["dof"] = self.n_data - len(self.names) - 1
        return out

    def marginal(self, name) -> tuple[np.ndarray, np.ndarray]:
        k = self.names.index(name)
        other = tuple(j for j in range(len(self.names)) if j != k)
        return self.axes[k], self.weights.sum(axis=other) if other else self.weights

    def h0_density(self, h0_grid) -> np.ndarray:
        """Exact marginal posterior of H0: a weighted sum of truncated Gaussians."""
        w = self.weights.ravel()
        keep = w > 1e-10 * w.max()
        mu, sig, w = self.mu.ravel()[keep], self.sigma.ravel()[keep], w[keep]
        lo, hi = self.h0_box
        h = np.asarray(h0_grid, dtype=float)
        mass = stats.norm.cdf(hi, mu, sig) - stats.norm.cdf(lo, mu, sig)
        coef = w / (mass * sig * np.sqrt(2 * np.pi))
        total = np.zeros(h.size)
        for start in range(0, mu.size, 2000):  # chunks keep memory small
            sl = slice(start, start + 2000)
            u = (h[:, None] - mu[sl]) / sig[sl]
            total += np.exp(-0.5 * u**2) @ coef[sl]
        return np.where((h >= lo) & (h <= hi), total, 0.0) / w.sum()

    def sample(self, n, rng) -> dict[str, np.ndarray]:
        """Posterior draws of theta (jittered within grid cells) and H0."""
        w = self.weights.ravel()
        idx = rng.choice(w.size, size=n, p=w / w.sum())
        out = {}
        for name, ax, k in zip(self.names, self.axes,
                               np.unravel_index(idx, self.weights.shape)):
            step = np.median(np.diff(ax)) if ax.size > 1 else 0.0
            out[name] = ax[k] + rng.uniform(-0.5, 0.5, n) * step
        mu, sig = self.mu.ravel()[idx], self.sigma.ravel()[idx]
        lo, hi = self.h0_box
        a, b = (lo - mu) / sig, (hi - mu) / sig
        out["H0"] = stats.truncnorm.rvs(a, b, loc=mu, scale=sig, random_state=rng)
        return out


def _quadratic_terms(E, H, Cinv):
    """a = E.Ci.E, b = E.Ci.H for E with shape (..., n); c = H.Ci.H."""
    CiH = Cinv @ H
    a = np.einsum("...i,ij,...j->...", E, Cinv, E)
    b = E @ CiH
    c = float(H @ CiH)
    return a, b, c


def fit_grid(E, names, axes, H, cov, *, h0_box=cfg.H0_PRIOR, h0_prior=None,
             valid=None) -> GridFit:
    """Fit H = H0 * E on a grid of theta.

    ``E`` has shape grid_shape + (n,). ``h0_box`` is a flat prior on H0;
    ``h0_prior=(mean, sigma)`` adds a Gaussian prior instead (the box is then
    removed unless given explicitly). ``valid`` masks impossible universes.
    """
    Cinv = linalg.inv(cov)
    a, b, c = _quadratic_terms(np.nan_to_num(E), np.asarray(H, float), Cinv)
    if h0_prior is not None:
        m0, s0 = h0_prior
        a, b, c = a + 1 / s0**2, b + m0 / s0**2, c + m0**2 / s0**2
    mu, sigma = b / a, 1.0 / np.sqrt(a)
    chi2_min = c - b**2 / a
    lo, hi = h0_box
    log_mass = np.log(np.clip(stats.norm.cdf(hi, mu, sigma) - stats.norm.cdf(lo, mu, sigma),
                              1e-300, None))
    ln_like = -0.5 * chi2_min + np.log(sigma) + log_mass
    h0_star = np.clip(mu, lo, hi)
    chi2_prof = a * h0_star**2 - 2 * b * h0_star + c
    if valid is None:
        valid = np.all(np.isfinite(E), axis=-1)
    ln_like = np.where(valid, ln_like, -np.inf)
    chi2_prof = np.where(valid, chi2_prof, np.nan)
    return GridFit(tuple(names), tuple(axes), ln_like, mu, sigma, chi2_prof, (lo, hi),
                   n_data=np.size(H))


def fit_flat_lcdm(z, H, cov, *, om_axis=None, **kw) -> GridFit:
    om = np.linspace(*cfg.OM_PRIOR, 1961) if om_axis is None else om_axis
    E = np.sqrt(om[:, None] * (1 + z) ** 3 + (1 - om[:, None]))
    return fit_grid(E, ("Om",), (om,), H, cov, **kw)


def _lcdm_grid_E(z, om, ol2d):
    """E(z) for matter + Lambda universes on a 2-D grid (ol2d has shape
    (len(om), m)), and a mask of universes whose H^2 stays positive from today
    to the highest data redshift. Built one Omega_m row at a time."""
    zz = np.linspace(0, z.max(), 200)
    E = np.empty(ol2d.shape + (z.size,))
    valid = np.empty(ol2d.shape, dtype=bool)
    for i, o in enumerate(om):
        valid[i] = np.all(cosmo.E2(zz[None, :], o, ol2d[i][:, None]) > 0, axis=-1)
        e2 = cosmo.E2(z[None, :], o, ol2d[i][:, None])
        E[i] = np.sqrt(np.where(e2 > 0, e2, np.nan))
    return E, valid


def fit_lcdm(z, H, cov, *, om_axis=None, ol_axis=None, **kw) -> GridFit:
    """Matter + Lambda with free curvature. Universes whose H^2 turns negative
    anywhere between today and the highest data redshift are excluded."""
    om = np.linspace(0.0, 1.0, 401) if om_axis is None else om_axis
    ol = np.linspace(*cfg.OL_RANGE, 526) if ol_axis is None else ol_axis
    E, valid = _lcdm_grid_E(np.asarray(z, float), om, np.broadcast_to(ol, (om.size, ol.size)))
    return fit_grid(E, ("Om", "OL"), (om, ol), H, cov, valid=valid, **kw)


def fit_lcdm_zt(z, H, cov, *, om_axis=None, zt_axis=None, **kw) -> GridFit:
    """Moresco et al. (2016)'s parametrisation of the same model: flat priors on
    Omega_m in [0, 1] and on the transition redshift z_t in [0, 1.5], with
    Omega_L = Omega_m (1 + z_t)^3 / 2 and curvature free."""
    om = np.linspace(0.0, 1.0, 401)[1:] if om_axis is None else om_axis
    zt = np.linspace(0.0, 1.5, 601) if zt_axis is None else zt_axis
    E, valid = _lcdm_grid_E(np.asarray(z, float), om, om[:, None] * (1 + zt[None, :]) ** 3 / 2.0)
    return fit_grid(E, ("Om", "zt"), (om, zt), H, cov, valid=valid, **kw)


def flat_lcdm_plane(z, H, cov, h0_axis, om_axis) -> np.ndarray:
    """chi^2 on an (H0, Omega_m) grid, for contour plots. Shape (n_om, n_h0)."""
    Cinv = linalg.inv(cov)
    E = np.sqrt(om_axis[:, None] * (1 + z) ** 3 + (1 - om_axis[:, None]))
    a, b, c = _quadratic_terms(E, np.asarray(H, float), Cinv)
    h = np.asarray(h0_axis)[None, :]
    return a[:, None] * h**2 - 2 * b[:, None] * h + c


def chi2(resid, cov) -> float:
    resid = np.asarray(resid, dtype=float)
    return float(resid @ linalg.solve(cov, resid, assume_a="pos"))


# ---------------------------------------------------------------------------
# A model-independent transition: two straight lines in adot(z)
# ---------------------------------------------------------------------------


def adot_data(z, H, cov):
    """adot = H / (1 + z) and its covariance."""
    J = 1.0 / (1.0 + np.asarray(z))
    return np.asarray(H) * J, cov * np.outer(J, J)


def piecewise_transition(z, y, cov, *, zt_axis=None, anchor=None) -> dict:
    """Fit adot(z) with two straight lines that meet at z_t (Moresco et al. 2016),
    and a single straight line (no transition), by generalised least squares.

    ``anchor=(H0, sigma)`` adds adot(0) = H0 as an independent data point.
    The amplitude and both slopes are solved exactly at each z_t; the z_t
    posterior marginalises over them with flat priors (the log-determinant
    term), and the profile chi^2 is kept as well.
    """
    zt = np.linspace(*cfg.ZT_GRID) if zt_axis is None else zt_axis
    z, y = np.asarray(z, float), np.asarray(y, float)
    if anchor is not None:
        z = np.concatenate([[0.0], z])
        y = np.concatenate([[anchor[0]], y])
        cov = linalg.block_diag([[anchor[1] ** 2]], cov)
    Cinv = linalg.inv(cov)
    n = z.size
    chi2s, logdet, slopes = np.empty(zt.size), np.empty(zt.size), np.empty((zt.size, 3))
    for k, t in enumerate(zt):
        X = np.column_stack([np.ones(n), np.minimum(z - t, 0.0), np.maximum(z - t, 0.0)])
        F = X.T @ Cinv @ X
        sign, ld = np.linalg.slogdet(F)
        if sign <= 0:  # a segment with fewer than two points: slope undetermined
            chi2s[k], logdet[k], slopes[k] = np.nan, np.nan, np.nan
            continue
        beta = np.linalg.solve(F, X.T @ Cinv @ y)
        r = y - X @ beta
        chi2s[k], logdet[k], slopes[k] = r @ Cinv @ r, ld, beta
    ok = np.isfinite(chi2s)
    ln_post = np.where(ok, -0.5 * chi2s - 0.5 * logdet, -np.inf)
    post = np.exp(ln_post - np.max(ln_post))
    X1 = np.column_stack([np.ones(n), z])
    F1 = X1.T @ Cinv @ X1
    b1 = np.linalg.solve(F1, X1.T @ Cinv @ y)
    r1 = y - X1 @ b1
    chi2_line = float(r1 @ Cinv @ r1)
    k_best = int(np.nanargmin(chi2s))
    chi2_pw = float(chi2s[k_best])

    def ic(chi2_min, k):
        aic = chi2_min + 2 * k
        return {"chi2": chi2_min, "k": k, "AICc": aic + 2 * k * (k + 1) / (n - k - 1),
                "BIC": chi2_min + k * np.log(n)}

    prof = np.exp(-0.5 * (chi2s - chi2_pw))
    return {
        "zt_axis": zt, "chi2": chi2s, "posterior": post,
        "summary_marginal": grid_summary(zt, post),
        "summary_profile": grid_summary(zt, np.nan_to_num(prof)),
        "best": {"zt": float(zt[k_best]), "amplitude": float(slopes[k_best, 0]),
                 "slope_low_z": float(slopes[k_best, 1]),
                 "slope_high_z": float(slopes[k_best, 2]), "chi2": chi2_pw},
        "line": {"intercept": float(b1[0]), "slope": float(b1[1]), "chi2": chi2_line},
        "piecewise_ic": ic(chi2_pw, 4), "line_ic": ic(chi2_line, 2),
        "delta_chi2_line": chi2_line - chi2_pw, "n": n,
    }
