"""How long do giants grow? The length distribution and its selection effects.

The model (Oei et al. 2023, A&A 672, A163), in words:

1. **Intrinsic lengths are Pareto distributed** above some minimum: the chance
   of a radio galaxy being longer than l falls as l^(xi + 1). The tail index
   xi is the one number of physical interest.
2. **We see them projected.** A jet pair at angle Theta to our line of sight
   has projected length l_p = l sin(Theta), with random orientations. Oei et
   al. proved that projected lengths are then Pareto distributed with the
   *same* xi (above the cut-off), so xi can be read from projected lengths.
3. **Surveys are not perfect.** Two selection effects bend the observed
   distribution: (a) the search only inspected sources larger than about 5'
   on the sky (a "fuzzy" threshold, ramping linearly from phi_min to
   phi_max), which removes small, distant giants; (b) lobes fade as they grow
   (surface brightness ~ l^-2 for self-similar growth) and with redshift
   ((1+z)^-4 for spectral index -1), so long or distant giants sink below the
   noise. Surface brightness at the reference length varies from giant to
   giant (lognormal, median b_ref, width sigma_ref).

The likelihood of the observed projected lengths is the Pareto density times
the completeness C(l_p) -- the fraction of giants of that length, anywhere in
the volume out to z_max, that the survey would have caught -- renormalised.
We evaluate it on a regular grid over (xi, phi half-width, b_ref, sigma_ref),
with Oei et al.'s priors, and marginalise. ``test_physics.py`` checks the
machinery on simulated catalogues drawn from the same model.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import special, stats

import config
import cosmology

LM = config.LENGTH_MODEL


# --------------------------------------------------------------------------
# Pareto distributions and projection
# --------------------------------------------------------------------------

def pareto_mle(lp, l_min: float = config.GIANT_THRESHOLD_MPC) -> tuple[float, float]:
    """Maximum-likelihood tail index of lengths above ``l_min`` ignoring all
    selection effects (Oei et al. 2023, Eq. 17), with its standard error."""
    lp = np.asarray(lp, float)
    lp = lp[lp > l_min]
    n = lp.size
    s = np.log(lp / l_min).sum()
    xi = -n / s - 1
    return float(xi), float(n / s / np.sqrt(n))


def survival(lp, xi: float, l_min: float = config.GIANT_THRESHOLD_MPC) -> np.ndarray:
    """P(L_p > l_p) for a Pareto distribution with tail index xi."""
    return (np.asarray(lp, float) / l_min) ** (xi + 1)


def empirical_survival(lp) -> tuple[np.ndarray, np.ndarray]:
    """The empirical survival function: fraction of the sample longer than each value."""
    x = np.sort(np.asarray(lp, float))
    n = x.size
    return x, 1.0 - np.arange(n) / n


def _int_x_pow_over_sqrt(a: float) -> float:
    """integral_1^inf x^a / sqrt(x^2 - 1) dx  (for a < 0)."""
    return float(np.sqrt(np.pi) / 2 * special.gamma(-a / 2) / special.gamma((1 - a) / 2))


def deprojection(lp: float, xi: float) -> dict:
    """What a projected length implies for the true length, given xi.

    For random orientations, a giant of true length l is seen with projected
    length l_p < l with density l_p / (l sqrt(l^2 - l_p^2)). Combined with a
    Pareto prior on l, the posterior of x = l / l_p is proportional to
    x^(xi - 1) / sqrt(x^2 - 1) for x > 1 (Oei et al. 2023, Sect. 2.4).
    """
    mean_x = _int_x_pow_over_sqrt(xi) / _int_x_pow_over_sqrt(xi - 1)
    x = 1 + np.geomspace(1e-7, 1e3, 200_000)
    pdf = x ** (xi - 1) / np.sqrt(x**2 - 1)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(x))])
    # add the analytic piece below x[0]: integral_1^x0 ~ sqrt(2 (x0 - 1))
    head = np.sqrt(2 * (x[0] - 1))
    cdf = (cdf + head) / (cdf[-1] + head + _tail_beyond(x[-1], xi))
    q = np.interp([0.16, 0.5, 0.84], cdf, x)
    return {"mean": lp * mean_x, "median": lp * q[1], "p16": lp * q[0], "p84": lp * q[2],
            "mean_factor": mean_x}


def _tail_beyond(x0: float, xi: float) -> float:
    # x^(xi-1)/sqrt(x^2-1) ~ x^(xi-2) for large x
    return x0 ** (xi - 1) / (1 - xi)


# --------------------------------------------------------------------------
# Selection effects and completeness
# --------------------------------------------------------------------------

@dataclass
class SelectionGrid:
    """Everything about the survey volume that does not depend on the
    model parameters, precomputed once."""

    lp: np.ndarray  # projected-length grid (Mpc), log spaced
    z: np.ndarray  # redshift grid
    w: np.ndarray  # comoving-volume weights r^2 / E dz, normalised
    phi: np.ndarray  # (n_lp, n_z) angular length, arcmin
    log_s0: np.ndarray  # (n_lp, n_z) ln( b_th (l_p/l_ref)^-zeta (1+z)^(3-alpha) )
    z_max: float


def selection_grid(z_max: float, cosmo: config.Cosmology = config.PLANCK18, *,
                   lp_max: float = 60.0, n_lp: int = 320, n_z: int = 400) -> SelectionGrid:
    lp = np.geomspace(config.GIANT_THRESHOLD_MPC, lp_max, n_lp)
    z = np.linspace(z_max / n_z / 2, z_max, n_z)
    dz = z[1] - z[0]
    w = cosmology.comoving_volume_element(z, cosmo) * dz
    w = w / w.sum()
    phi = cosmology.angle_arcmin(lp[:, None], z[None, :], cosmo)
    zeta, alpha = LM["zeta"], LM["lobe_spectral_index"]
    log_s0 = (np.log(LM["b_threshold_jy_deg2"]) - zeta * np.log(lp[:, None] / LM["l_ref_mpc"])
              + (3 - alpha) * np.log1p(z[None, :]))
    return SelectionGrid(lp, z, w, phi, log_s0, z_max)


def angular_detection(phi, half_width: float, centre: float = LM["phi_threshold_arcmin"]) -> np.ndarray:
    """Fuzzy angular-length threshold: probability 0 below phi_min, 1 above
    phi_max, linear in between (Oei et al. 2023, Eq. 22)."""
    if half_width <= 0:
        return (np.asarray(phi) >= centre).astype(float)
    return np.clip((np.asarray(phi) - (centre - half_width)) / (2 * half_width), 0.0, 1.0)


def sb_detection_table(sigma_ref: float, zeta: float = LM["zeta"]) -> tuple[np.ndarray, np.ndarray]:
    """g(ln s_min) = P(detected) for a giant whose minimal detectable
    surface-brightness factor is s_min (Oei et al. 2023, Eq. A.83):

        g = integral_{s_min}^inf sqrt(1 - (s_min / s)^(-2/zeta)) f_S(s) ds,

    with S lognormal (median 1, dispersion sigma_ref). The square root is the
    probability, over random orientations, that a giant with surface-brightness
    factor s is seen bright enough: sin(Theta) must exceed (s_min/s)^(-1/zeta).
    Returned as a table in ln s_min for fast interpolation.
    """
    sig = max(sigma_ref, 0.05)
    # Far below the lognormal's bulk every giant is bright enough and only the
    # orientation factor matters; it reaches 1 only once s_min/s < ~e^-15.
    ln_smin = np.linspace(-(15 + 10 * sig), 10 * sig + 1, 900)
    # integrate in t = ln s - ln s_min >= 0 (the square root ~ sqrt(t) near 0)
    t = np.concatenate([[0.0], np.geomspace(1e-6, 1.0, 120),
                        np.linspace(1.0 + 1e-3, 15 + 25 * sig, 1500)])
    ln_s = ln_smin[:, None] + t[None, :]
    root = np.sqrt(np.clip(1 - np.exp(t * (2 / zeta)), 0, 1))  # (s_min/s)^(-2/zeta) = exp(-t * (-2/zeta))
    # f_S(s) ds = normal pdf in ln s
    dens = stats.norm.pdf(ln_s, scale=sigma_ref)
    integrand = root[None, :] * dens
    g = np.trapezoid(integrand, t, axis=1)
    return ln_smin, np.clip(g, 0.0, 1.0)


def completeness(grid: SelectionGrid, half_width: float, b_ref: float, sigma_ref: float) -> np.ndarray:
    """C(l_p): fraction of giants of each projected length, distributed with
    constant comoving density out to z_max, that the survey detects."""
    p_al = angular_detection(grid.phi, half_width)
    ln_smin_tab, g = sb_detection_table(sigma_ref)
    p_sb = np.interp(grid.log_s0 - np.log(b_ref), ln_smin_tab, g, left=1.0, right=0.0)
    return (p_al * p_sb) @ grid.w


# --------------------------------------------------------------------------
# The posterior on a grid
# --------------------------------------------------------------------------

@dataclass
class Posterior:
    xi: np.ndarray
    half_width: np.ndarray
    b_ref: np.ndarray
    sigma_ref: np.ndarray
    log_post: np.ndarray  # (xi, half_width, b_ref, sigma_ref)
    n: int
    z_max: float

    def marginal(self, axis_name: str) -> tuple[np.ndarray, np.ndarray]:
        names = ["xi", "half_width", "b_ref", "sigma_ref"]
        axis = names.index(axis_name)
        p = np.exp(self.log_post - self.log_post.max())
        other = tuple(i for i in range(4) if i != axis)
        m = p.sum(axis=other)
        x = getattr(self, axis_name)
        m = m / np.trapezoid(m, x)
        return x, m

    def joint(self, a: str, b: str) -> np.ndarray:
        names = ["xi", "half_width", "b_ref", "sigma_ref"]
        ia, ib = names.index(a), names.index(b)
        p = np.exp(self.log_post - self.log_post.max())
        other = tuple(i for i in range(4) if i not in (ia, ib))
        j = p.sum(axis=other)
        return j if ia < ib else j.T

    def summary(self, axis_name: str) -> dict:
        x, m = self.marginal(axis_name)
        cdf = np.concatenate([[0], np.cumsum(0.5 * (m[1:] + m[:-1]) * np.diff(x))])
        cdf /= cdf[-1]
        mean = np.trapezoid(x * m, x)
        sd = np.sqrt(np.trapezoid((x - mean) ** 2 * m, x))
        q = np.interp([0.025, 0.16, 0.5, 0.84, 0.975], cdf, x)
        return {"mean": float(mean), "sd": float(sd), "median": float(q[2]),
                "p16": float(q[1]), "p84": float(q[3]), "p2.5": float(q[0]), "p97.5": float(q[4]),
                "mode": float(x[np.argmax(m)])}

    def map_point(self) -> dict:
        i = np.unravel_index(np.argmax(self.log_post), self.log_post.shape)
        return {"xi": float(self.xi[i[0]]), "half_width": float(self.half_width[i[1]]),
                "b_ref": float(self.b_ref[i[2]]), "sigma_ref": float(self.sigma_ref[i[3]])}


def default_axes() -> dict:
    return {"xi": np.round(np.arange(-6.0, -1.5 + 1e-9, 0.025), 4),
            "half_width": np.round(np.arange(0.0, 4.0 + 1e-9, 0.1), 4),
            "b_ref": np.arange(25.0, 3500.0 + 1e-9, 25.0),
            "sigma_ref": np.round(np.arange(0.1, 3.0 + 1e-9, 0.05), 4)}


def log_prior(axes: dict) -> np.ndarray:
    """Oei et al.'s priors: flat in xi and the half-width; Gaussian in b_ref
    and sigma_ref (68% intervals [250, 1750] Jy/deg^2 and [0.65, 1.95])."""
    mb, sb = LM["b_ref_prior_jy_deg2"]
    ms, ss = LM["sigma_ref_prior"]
    lb = stats.norm.logpdf(axes["b_ref"], mb, sb)
    ls = stats.norm.logpdf(axes["sigma_ref"], ms, ss)
    return lb[None, None, :, None] + ls[None, None, None, :]


def grid_posterior(lp_data, z_max: float, cosmo: config.Cosmology = config.PLANCK18,
                   axes: dict | None = None, *, progress: bool = False) -> Posterior:
    """Posterior over (xi, half_width, b_ref, sigma_ref) given projected lengths."""
    axes = axes or default_axes()
    lp_data = np.asarray(lp_data, float)
    lp_data = lp_data[lp_data > config.GIANT_THRESHOLD_MPC]
    n = lp_data.size
    grid = selection_grid(z_max, cosmo)
    lp = grid.lp
    # interpolation of ln C at the data, in ln l_p
    x = np.log(lp)
    xd = np.log(lp_data)
    j = np.clip(np.searchsorted(x, xd) - 1, 0, x.size - 2)
    t = (xd - x[j]) / (x[j + 1] - x[j])
    sum_log_lp = np.log(lp_data).sum()
    # normalisation integral Z(xi) = int C(l) l^xi dl on the log grid (trapezoid in ln l)
    dlnl = np.gradient(x)
    pow_tab = np.exp(np.outer(x, axes["xi"] + 1)) * dlnl[:, None]  # l^xi * l dln l
    p_al_all = np.stack([angular_detection(grid.phi, h) for h in axes["half_width"]])  # (h, lp, z)
    xi = axes["xi"]
    logL = np.empty((xi.size, axes["half_width"].size, axes["b_ref"].size, axes["sigma_ref"].size))
    for k, sig in enumerate(axes["sigma_ref"]):
        if progress and k % 10 == 0:
            print(f"    sigma_ref {sig:.2f} ({k + 1}/{axes['sigma_ref'].size})", flush=True)
        ln_smin_tab, g = sb_detection_table(sig)
        for ib, b in enumerate(axes["b_ref"]):
            p_sb = np.interp(grid.log_s0 - np.log(b), ln_smin_tab, g, left=1.0, right=0.0)  # (lp, z)
            C = np.einsum("hlz,lz,z->hl", p_al_all, p_sb, grid.w, optimize=True)  # (h, lp)
            lnC = np.log(np.maximum(C, 1e-300))
            lnC_data = ((1 - t) * lnC[:, j] + t * lnC[:, j + 1]).sum(axis=1)  # (h,)
            Z = C @ pow_tab  # (h, xi)
            logL[:, :, ib, k] = (xi[:, None] * sum_log_lp + lnC_data[None, :]
                                 - n * np.log(np.maximum(Z.T, 1e-300)))
    log_post = logL + log_prior(axes)
    return Posterior(xi, axes["half_width"], axes["b_ref"], axes["sigma_ref"], log_post, n, z_max)


def predicted_survival(grid: SelectionGrid, C: np.ndarray, xi: float) -> np.ndarray:
    """Model survival function of observed projected lengths on ``grid.lp``."""
    lnl = np.log(grid.lp)
    dens = C * grid.lp ** (xi + 1)  # per unit ln l
    cum = np.concatenate([[0], np.cumsum(0.5 * (dens[1:] + dens[:-1]) * np.diff(lnl))])
    return 1 - cum / cum[-1]


def posterior_predictive(post: Posterior, cosmo: config.Cosmology = config.PLANCK18,
                         n_draws: int = 200, seed: int = 1) -> dict:
    """Survival curves for parameter sets drawn from the posterior."""
    rng = np.random.default_rng(seed)
    p = np.exp(post.log_post - post.log_post.max()).ravel()
    p /= p.sum()
    idx = rng.choice(p.size, size=n_draws, p=p)
    grid = selection_grid(post.z_max, cosmo)
    curves, comps = [], []
    for flat in idx:
        i = np.unravel_index(flat, post.log_post.shape)
        C = completeness(grid, post.half_width[i[1]], post.b_ref[i[2]], post.sigma_ref[i[3]])
        curves.append(predicted_survival(grid, C, post.xi[i[0]]))
        comps.append(C)
    return {"lp": grid.lp, "survival": np.array(curves), "completeness": np.array(comps)}


# --------------------------------------------------------------------------
# Simulated catalogues (for testing the machinery)
# --------------------------------------------------------------------------

def simulate(n_target: int, xi: float, half_width: float, b_ref: float, sigma_ref: float,
             z_max: float, cosmo: config.Cosmology = config.PLANCK18, l_min: float = 0.3,
             seed: int = 0) -> dict:
    """Draw giants from the full forward model: Pareto intrinsic lengths above
    ``l_min``, isotropic orientations, constant comoving density out to z_max,
    lognormal surface brightness; keep what the survey would detect."""
    rng = np.random.default_rng(seed)
    zg = np.linspace(1e-4, z_max, 4000)
    cdf_z = np.cumsum(cosmology.comoving_volume_element(zg, cosmo))
    cdf_z /= cdf_z[-1]
    kept_lp, kept_z = [], []
    zeta, alpha = LM["zeta"], LM["lobe_spectral_index"]
    while sum(len(k) for k in kept_lp) < n_target:
        m = 200_000
        L = l_min * rng.random(m) ** (1 / (xi + 1))  # Pareto: P(L > l) = (l/l_min)^(xi+1)
        cos_t = rng.random(m)
        lp = L * np.sqrt(1 - cos_t**2)
        z = np.interp(rng.random(m), cdf_z, zg)
        s = np.exp(rng.normal(0, sigma_ref, m))
        b = b_ref * s * (L / LM["l_ref_mpc"]) ** zeta / (1 + z) ** (3 - alpha)
        phi = cosmology.angle_arcmin(lp, z, cosmo)
        p_al = angular_detection(phi, half_width)
        keep = (lp > config.GIANT_THRESHOLD_MPC) & (b > LM["b_threshold_jy_deg2"]) & (rng.random(m) < p_al)
        kept_lp.append(lp[keep])
        kept_z.append(z[keep])
    lp = np.concatenate(kept_lp)[:n_target]
    z = np.concatenate(kept_z)[:n_target]
    return {"lp": lp, "z": z}


# --------------------------------------------------------------------------
# The analysis
# --------------------------------------------------------------------------

def oei_sample(z_max: float) -> np.ndarray:
    """Oei et al.'s analysis sample: their giants with a definite host (not a
    lower limit) below z_max; 1,473 at z < 0.5 and 811 at z < 0.25."""
    import catalogs

    o = catalogs.oei2023()
    keep = (o["z"] < z_max) & ~o["host_is_nearest_candidate"]
    return o.loc[keep, "lp_mpc"].to_numpy(float)


def cosmology_effect() -> dict:
    """How much the adopted cosmology moves every length: our three cosmologies
    applied to the 11,585 giants' angles and redshifts."""
    import catalogs

    g = catalogs.all_giants()
    ok = np.isfinite(g["z"]) & (g["z"] > 0)
    base = cosmology.projected_length_mpc(g.loc[ok, "las_arcmin"], g.loc[ok, "z"], config.PLANCK18)
    out = {"n": int(ok.sum()),
           "catalogue_vs_planck_median_ratio": float(np.median(base / g.loc[ok, "lp_mpc"])),
           "catalogue_vs_planck_p1_p99": [float(v) for v in np.percentile(base / g.loc[ok, "lp_mpc"], [1, 99])]}
    for c in config.COSMOLOGIES:
        lp = cosmology.projected_length_mpc(g.loc[ok, "las_arcmin"], g.loc[ok, "z"], c)
        out[c.name] = {"median_ratio_to_planck": float(np.median(lp / base)),
                       "n_giants": int((lp >= config.GIANT_THRESHOLD_MPC).sum()),
                       "n_above_5mpc": int((lp >= 5).sum())}
    return out
