"""Weigh the diffuse baryons: fit the Macquart relation to localized FRBs.

For each burst with a host redshift z we know DM_EG, the dispersion measure left
after removing the Milky Way. It is the sum of two unknowns:

    DM_EG = DM_cosmic + DM_host / (1 + z)

The model (Macquart et al. 2020, Methods) gives both a probability distribution:

* DM_cosmic follows p(Delta) around the mean relation <DM_cosmic>(z), scaled by
  **f_d**, the fraction of all baryons that are in diffuse ionized gas, with a
  spread sigma = F / sqrt(z) (see ``frb_physics``).
* DM_host is lognormal with median exp(mu) and log-width sigma_host, free.

On top of Macquart et al.'s model we allow for one more thing they fixed: the
Milky Way term we subtracted is itself uncertain on each sightline (NE2001 errors,
direction-dependent halo gas). We treat it as Gaussian scatter of width
``DM_MW_SCATTER`` about the subtracted value. Without that term the fit, having
no other way to explain a nearby burst whose DM_EG is off by a few tens of
pc cm^-3, inflates the *cosmic* scatter instead -- see the README.

The likelihood of one burst is the convolution of the three, and the total is
the product over bursts. Cosmology is held at Planck 2018 (H0, Omega_b, Omega_m),
so f_d is the only thing the mean relation can adjust. Because the mean DM
scales as f_d * Omega_b * H0, the same fit can be read the other way round:
assume f_d from the census of stars and cold gas and solve for Omega_b h70.

The posterior is computed exhaustively on a 4-D grid with flat priors over the
ranges in ``config.FIT_GRID`` -- no sampler, no random numbers, fully
reproducible.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.cosmology import Planck18

import frb_physics as fp
from config import (
    DM_MW_SCATTER,
    F_DIFFUSE_FIDUCIAL,
    FIT_GRID,
    FIT_N_CONV,
    FIT_N_MW_NODES,
    MACQUART_ALPHA,
    MACQUART_BETA,
)

PARAM_NAMES = ("f_d", "log10_F", "mu", "sigma_host")


def grid_axes(spec: dict | None = None) -> dict[str, np.ndarray]:
    spec = spec or FIT_GRID
    return {k: np.linspace(*spec[k]) for k in PARAM_NAMES}


def omega_b_h70_planck(cosmo=Planck18) -> float:
    """Planck 2018's baryon density in the units FRB papers quote, Omega_b * H0/70."""
    return float(cosmo.Ob0 * cosmo.H0.value / 70.0)


def _likelihood_one(d: float, zi: float, m1: float, axes: dict[str, np.ndarray], n_conv: int,
                    alpha: float, beta: float) -> np.ndarray:
    """p(DM_EG = d | z, params) on the grid, shape (n_f_d * n_F, n_mu * n_sigma_host).

    The host's observed-frame contribution x runs over a midpoint grid on
    (0, d); the cosmic part is d - x. The cosmic PDF depends on (f_d, F) and the
    host PDF on (mu, sigma_host), so the convolution at every grid point at once
    is a single matrix product.
    """
    f_d, mu, sh = axes["f_d"], axes["mu"], axes["sigma_host"]
    F = 10.0 ** axes["log10_F"]
    k = np.arange(n_conv) + 0.5
    x = d * k / n_conv
    # Cosmic part, shape (n_f_d, n_F, n_conv)
    mean = f_d * m1
    delta = (d - x)[None, None, :] / mean[:, None, None]
    sig = (F / np.sqrt(zi))[None, :, None]
    pc = np.exp(fp.macquart_log_pdf_delta(delta, sig, alpha, beta)) / mean[:, None, None]
    # Host part, shape (n_mu, n_sigma_host, n_conv); its rest-frame DM is x (1+z)
    host = x * (1 + zi)
    lh = np.log(host)[None, None, :]
    px = ((1 + zi) * np.exp(-((lh - mu[:, None, None]) ** 2) / (2 * sh[None, :, None] ** 2))
          / (host[None, None, :] * sh[None, :, None] * np.sqrt(2 * np.pi)))
    return (pc.reshape(-1, n_conv) @ px.reshape(-1, n_conv).T) * (d / n_conv)


def log_likelihood_grid(dm_eg, z, axes: dict[str, np.ndarray], *, sigma_mw: float = DM_MW_SCATTER,
                        n_mw: int = FIT_N_MW_NODES, n_conv: int = FIT_N_CONV,
                        alpha: float = MACQUART_ALPHA, beta: float = MACQUART_BETA) -> np.ndarray:
    """Summed log-likelihood over bursts, on the full parameter grid.

    Returns an array of shape (n_f_d, n_F, n_mu, n_sigma_host).

    The Milky Way scatter is integrated with Gauss-Hermite quadrature: the
    likelihood of the observed DM_EG is the model density averaged over
    ``n_mw`` true values DM_EG + sqrt(2) sigma_mw x_k with weights w_k/sqrt(pi).
    Nodes that would make the true value negative carry zero density.
    """
    dm_eg = np.asarray(dm_eg, dtype=float)
    z = np.asarray(z, dtype=float)
    shape = tuple(len(axes[n]) for n in PARAM_NAMES)
    mean_unit = fp.mean_dm_cosmic(z, f_d=1.0)
    if sigma_mw > 0:
        nodes, weights = np.polynomial.hermite.hermgauss(n_mw)
        offsets, weights = np.sqrt(2.0) * sigma_mw * nodes, weights / np.sqrt(np.pi)
    else:
        offsets, weights = np.array([0.0]), np.array([1.0])

    total = np.zeros((shape[0] * shape[1], shape[2] * shape[3]))
    for d, zi, m1 in zip(dm_eg, z, mean_unit):
        like = np.zeros_like(total)
        for off, w in zip(offsets, weights):
            if d + off > 0:
                like += w * _likelihood_one(d + off, zi, m1, axes, n_conv, alpha, beta)
        total += np.log(np.maximum(like, 1e-300))
    return total.reshape(shape)


def _quantiles(x: np.ndarray, pdf: np.ndarray, qs=(0.025, 0.16, 0.5, 0.84, 0.975)) -> dict:
    if len(x) == 1:  # parameter held fixed
        return {f"q{int(round(q * 1000)):03d}": float(x[0]) for q in qs}
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    # Centre each cell: the cumulative mass reaches cdf[i] at the cell's upper edge.
    edges = np.concatenate([[x[0] - (x[1] - x[0]) / 2], (x[:-1] + x[1:]) / 2, [x[-1] + (x[-1] - x[-2]) / 2]])
    cdf_e = np.concatenate([[0.0], cdf])
    return {f"q{int(round(q * 1000)):03d}": float(np.interp(q, cdf_e, edges)) for q in qs}


@dataclass
class FitResult:
    label: str
    names: list[str]
    z: np.ndarray
    dm_eg: np.ndarray
    axes: dict[str, np.ndarray]
    log_post: np.ndarray  # normalised: logsumexp over the grid = 0
    marginals: dict[str, np.ndarray] = field(default_factory=dict)
    summary: dict[str, dict] = field(default_factory=dict)
    map_params: dict[str, float] = field(default_factory=dict)
    edge_warnings: list[str] = field(default_factory=list)
    sigma_mw: float = DM_MW_SCATTER

    def to_json(self) -> dict:
        ob = omega_b_h70_planck()
        s = self.summary
        fd = s["f_d"]
        return {
            "label": self.label,
            "n_frb": len(self.names),
            "sigma_mw": self.sigma_mw,
            "z_range": [float(self.z.min()), float(self.z.max())],
            "f_d": fd,
            "omega_b_h70_if_fd_fiducial": {k: v * ob / F_DIFFUSE_FIDUCIAL for k, v in fd.items()},
            "log10_F": s["log10_F"],
            "host_median_dm": {k: float(np.exp(v)) for k, v in s["mu"].items()},
            "sigma_host": s["sigma_host"],
            "map": self.map_params,
            "prior_edge_warnings": self.edge_warnings,
        }


def fit(names, z, dm_eg, *, label: str = "baseline", axes: dict | None = None,
        sigma_mw: float = DM_MW_SCATTER, n_conv: int = FIT_N_CONV) -> FitResult:
    axes = axes or grid_axes()
    logl = log_likelihood_grid(dm_eg, z, axes, sigma_mw=sigma_mw, n_conv=n_conv)
    log_post = logl - logl.max()
    log_post -= np.log(np.exp(log_post).sum())
    post = np.exp(log_post)

    res = FitResult(label=label, names=list(names), z=np.asarray(z), dm_eg=np.asarray(dm_eg),
                    axes=axes, log_post=log_post, sigma_mw=sigma_mw)
    for i, name in enumerate(PARAM_NAMES):
        other = tuple(j for j in range(4) if j != i)
        marg = post.sum(axis=other)
        res.marginals[name] = marg
        res.summary[name] = _quantiles(axes[name], marg)
        # A posterior still high at a grid edge means the prior range, not the
        # data, is setting the answer there.
        peak = marg.max()
        for side, val in ((("lower", marg[0]), ("upper", marg[-1])) if len(marg) > 1 else ()):
            if val > 0.1 * peak:
                res.edge_warnings.append(f"{name}: posterior at {side} grid edge is {val / peak:.2f} of peak")
    idx = np.unravel_index(np.argmax(log_post), log_post.shape)
    res.map_params = {n: float(axes[n][j]) for n, j in zip(PARAM_NAMES, idx)}
    return res


def posterior_draws(res: FitResult, n: int = 400, seed: int = 20260929) -> list[dict]:
    """Parameter sets drawn from the grid posterior (fixed seed: reproducible)."""
    rng = np.random.default_rng(seed)
    p = np.exp(res.log_post).ravel()
    flat = rng.choice(p.size, size=n, p=p / p.sum())
    out = []
    for f in flat:
        idx = np.unravel_index(f, res.log_post.shape)
        out.append({k: float(res.axes[k][j]) for k, j in zip(PARAM_NAMES, idx)})
    return out


def predictive_band(res: FitResult, z_grid, *, qs=(0.05, 0.16, 0.5, 0.84, 0.95), n_draws: int = 300,
                    dm_min: float = -300.0, dm_max: float = 4000.0, n_dm: int = 4300) -> dict[str, np.ndarray]:
    """Posterior-predictive quantiles of the *measured* DM_EG versus redshift.

    For each z the model density -- cosmic convolved with host, then with the
    Milky Way scatter -- is averaged over parameter sets drawn from the
    posterior, so the band carries both the physical scatter between sightlines
    and our uncertainty in the parameters. Because the Milky Way term is
    uncertain, a measured DM_EG can in principle be slightly negative; the grid
    extends below zero to allow it. Convolutions are done by FFT.
    """
    draws = posterior_draws(res, n_draws)
    step = (dm_max - dm_min) / n_dm
    dm = dm_min + (np.arange(n_dm) + 0.5) * step
    pos = dm > 0
    n_pos = int(pos.sum())
    dm_pos = dm[pos]
    kern_x = (np.arange(n_dm) - n_dm // 2) * step
    kern = np.exp(-0.5 * (kern_x / res.sigma_mw) ** 2) if res.sigma_mw > 0 else (kern_x == 0).astype(float)
    kern /= kern.sum()
    out = {f"q{int(round(q * 100)):02d}": np.zeros(len(z_grid)) for q in qs}
    out["mean_cosmic_fit"] = np.zeros(len(z_grid))
    for j, zj in enumerate(z_grid):
        acc = np.zeros(n_pos)
        m1 = fp.mean_dm_cosmic(zj, f_d=1.0)
        for p in draws:
            mean = p["f_d"] * m1
            pc = fp.macquart_pdf_delta(dm_pos / mean, 10 ** p["log10_F"] / np.sqrt(zj)) / mean
            px = (1 + zj) * fp.lognormal_pdf(dm_pos * (1 + zj), p["mu"], p["sigma_host"])
            acc += np.fft.irfft(np.fft.rfft(pc, 2 * n_pos) * np.fft.rfft(px, 2 * n_pos))[:n_pos] * step
        full = np.zeros(n_dm)
        full[pos] = acc / len(draws)
        # Milky Way scatter: a kernel centred on the grid keeps the alignment.
        smeared = _circular_to_linear(full, kern, n_dm)
        cdf = np.cumsum(smeared)
        cdf /= cdf[-1]
        for q in qs:
            out[f"q{int(round(q * 100)):02d}"][j] = np.interp(q, cdf, dm)
        out["mean_cosmic_fit"][j] = res.summary["f_d"]["q500"] * m1
    out["z"] = np.asarray(z_grid)
    return out


def _circular_to_linear(signal: np.ndarray, kernel_centred: np.ndarray, n: int) -> np.ndarray:
    """Linear convolution of ``signal`` with a kernel centred at index n//2,
    trimmed back to the signal's own grid ('same' mode)."""
    full = np.fft.irfft(np.fft.rfft(signal, 2 * n) * np.fft.rfft(kernel_centred, 2 * n))
    return full[n // 2: n // 2 + n]
