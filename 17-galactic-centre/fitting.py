"""Orbit fitting: likelihood, best fit, and MCMC for S2.

A *model* is a list of free parameters plus fixed values for the rest. The
full parameter set is

* the black hole: ``gm6`` (GM in 10^6 GM_sun), ``r0`` (distance, kpc);
* S2's orbit: ``period`` (yr), ``t0`` (the 2018 pericentre passage, yr),
  ``e``, ``inc``, ``omega``, ``Omega`` (deg);
* relativity switches: ``f`` (redshift, 0 = Newton, 1 = GR) and ``f_sp``
  (Schwarzschild precession, same convention);
* for **each team's reference frame** ``k``: the black hole's position
  ``x0_k``, ``y0_k`` (mas) and drift ``vx0_k``, ``vy0_k`` (mas/yr) relative
  to that frame's origin. The frames were built independently (UCLA from
  radio masers, Sakai et al. 2019; MPE from Plewa et al. 2015), so they get
  independent offsets -- the "fit an offset rather than mixing them naively"
  of the project brief;
* ``vz0``: a constant line-of-sight velocity offset (Sgr A*'s own motion is
  < 1 km/s; this absorbs LSR and line-shape systematics);
* ``nirc2_offset``: the extra correction Do et al. found for Keck's NIRC2
  spectrograph (2002-2005), added to those velocities; ``vlt_rv_offset`` is
  the same idea for the VLT's SINFONI velocities (used only as a systematics
  test);
* noise terms: an error-scale factor ``lnscale_*`` (Gillessen et al.'s
  approach: rescale until the fit is statistically consistent), or the
  correlated-noise model of Do et al. (2019): nearby positions on the sky
  share part of their error, because the same unseen neighbouring stars
  confuse them. ``corr_length_*`` (mas) and ``corr_mix_*`` (0-1).

With W walkers every function works on a ``(W, D)`` array of parameter
vectors at once, so emcee can evaluate a whole ensemble in one numpy call.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import optimize

import orbit_physics as op
from s2data import Astrometry, RadialVelocities

# --------------------------------------------------------------------------
# Parameter bookkeeping
# --------------------------------------------------------------------------

#: Uniform prior ranges, in the spirit of Do et al. (2019): wide boxes, flat
#: inside. Angles are not wrapped: the posteriors are a fraction of a degree wide.
BOUNDS = {
    "gm6": (1.0, 10.0),
    "r0": (4.0, 14.0),
    "period": (15.0, 17.0),
    "t0": (2017.8, 2019.0),
    "e": (0.70, 0.99),
    "inc": (90.0, 180.0),
    "omega": (0.0, 180.0),
    "Omega": (135.0, 315.0),
    "f": (-10.0, 10.0),
    "f_sp": (-50.0, 50.0),
    "vz0": (-300.0, 300.0),
    "nirc2_offset": (-500.0, 500.0),
    "vlt_rv_offset": (-300.0, 300.0),
}

#: Instrument-specific velocity corrections: parameter name -> instrument label.
RV_OFFSETS = {"nirc2_offset": "NIRC2", "vlt_rv_offset": "VLT"}
FRAME_BOUNDS = {"x0": (-30.0, 30.0), "y0": (-30.0, 30.0), "vx0": (-3.0, 3.0), "vy0": (-3.0, 3.0)}
NOISE_BOUNDS = {"lnscale": (-3.0, 3.0), "corr_length": (0.5, 200.0), "corr_mix": (0.0, 0.999)}


def param_bounds(name: str) -> tuple[float, float]:
    if name in BOUNDS:
        return BOUNDS[name]
    base = name.split("_")[0]
    if base in FRAME_BOUNDS:
        return FRAME_BOUNDS[base]
    for key, b in NOISE_BOUNDS.items():
        if name.startswith(key):
            return b
    raise KeyError(name)


# --------------------------------------------------------------------------
# Data + model description
# --------------------------------------------------------------------------


@dataclass
class Frame:
    """One team's astrometry, with its own reference frame and noise model."""

    key: str  # "keck" or "vlt": suffix of this frame's parameter names
    data: Astrometry
    t_ref: float  # reference epoch of x0/y0 (2000.0 for Keck as Do et al., 2009.0 for VLT)
    #: "fixed" (published errors), "scaled" (one free factor on all errors) or
    #: "correlated" (Do et al. 2019's covariance model)
    noise: str = "fixed"
    #: pairwise sky distances between the measured positions (for "correlated")
    dist: np.ndarray = field(default=None, repr=False)

    def __post_init__(self):
        dx = self.data.x[:, None] - self.data.x[None, :]
        dy = self.data.y[:, None] - self.data.y[None, :]
        self.dist = np.hypot(dx, dy)


@dataclass
class Model:
    """Which parameters are free, which are fixed, and which effects are on."""

    name: str
    frames: list[Frame]
    rv: RadialVelocities | None
    free: list[str]
    fixed: dict
    roemer: bool = True
    rv_noise: str = "fixed"  # or "scaled" (free lnscale_rv)
    description: str = ""
    #: Optional Gaussian priors {name: (mean, sigma)} on top of the flat boxes, e.g.
    #: the coordinate-system priors Gillessen et al. (2017) adopt from Plewa et al. (2015).
    gaussian_priors: dict = field(default_factory=dict)

    def __post_init__(self):
        self.lo = np.array([param_bounds(n)[0] for n in self.free])
        self.hi = np.array([param_bounds(n)[1] for n in self.free])
        if self.rv is not None:
            self.is_nirc2 = self.rv.instrument == "NIRC2"
            self.offset_masks = {n: self.rv.instrument == inst for n, inst in RV_OFFSETS.items()}
        names = set(self.free) | set(self.fixed)
        for needed in ("gm6", "r0", "period", "t0", "e", "inc", "omega", "Omega"):
            if needed not in names:
                raise ValueError(f"{self.name}: {needed} is neither free nor fixed")

    @property
    def ndim(self) -> int:
        return len(self.free)

    @property
    def n_data(self) -> int:
        n = sum(2 * f.data.n for f in self.frames)
        return n + (self.rv.n if self.rv is not None else 0)

    def unpack(self, theta) -> dict:
        """(W, D) array -> dict of (W, 1) arrays (fixed values broadcast)."""
        theta = np.atleast_2d(theta)
        p = {n: theta[:, k:k + 1] for k, n in enumerate(self.free)}
        for n, v in self.fixed.items():
            p.setdefault(n, np.full((theta.shape[0], 1), float(v)))
        return p

    def pack(self, values: dict) -> np.ndarray:
        return np.array([values[n] for n in self.free], dtype=float)

    # ---- the likelihood -------------------------------------------------

    def model_observables(self, p: dict) -> dict:
        """Predicted positions for each frame and velocities, for (W, 1) params."""
        out = {}
        t_all = np.concatenate([f.data.t for f in self.frames]) if self.frames else None
        pred = op.predict(p, t_all, self.rv.t if self.rv is not None else None, roemer=self.roemer)
        start = 0
        for fr in self.frames:
            n = fr.data.n
            dt = fr.data.t - fr.t_ref
            xs = pred["x"][:, start:start + n]
            ys = pred["y"][:, start:start + n]
            out[f"x_{fr.key}"] = xs + p[f"x0_{fr.key}"] + p[f"vx0_{fr.key}"] * dt
            out[f"y_{fr.key}"] = ys + p[f"y0_{fr.key}"] + p[f"vy0_{fr.key}"] * dt
            start += n
        if self.rv is not None:
            rv = pred["rv"] + p["vz0"]
            for name, mask in self.offset_masks.items():
                if name in p:
                    rv = rv - p[name] * mask
            out["rv"] = rv
        return out

    def log_likelihood(self, theta) -> np.ndarray:
        p = self.unpack(theta)
        m = self.model_observables(p)
        logl = np.zeros(np.atleast_2d(theta).shape[0])
        for fr in self.frames:
            rx = (fr.data.x - m[f"x_{fr.key}"]) / fr.data.x_err
            ry = (fr.data.y - m[f"y_{fr.key}"]) / fr.data.y_err
            n = fr.data.n
            if fr.noise == "fixed":
                logl += -0.5 * np.sum(rx**2 + ry**2, axis=1)
            elif fr.noise == "scaled":
                s2 = np.exp(2.0 * p[f"lnscale_astro_{fr.key}"][:, 0])
                logl += -0.5 * np.sum(rx**2 + ry**2, axis=1) / s2 - n * np.log(s2)
            elif fr.noise == "correlated":
                logl += _correlated_loglike(rx, ry, fr.dist, p[f"corr_length_{fr.key}"][:, 0],
                                            p[f"corr_mix_{fr.key}"][:, 0])
            else:
                raise ValueError(fr.noise)
            logl += -np.sum(np.log(fr.data.x_err) + np.log(fr.data.y_err)) - n * np.log(2 * np.pi)
        if self.rv is not None:
            rr = (self.rv.v - m["rv"]) / self.rv.v_err
            if self.rv_noise == "scaled":
                s2 = np.exp(2.0 * p["lnscale_rv"][:, 0])
                logl += -0.5 * np.sum(rr**2, axis=1) / s2 - 0.5 * self.rv.n * np.log(s2)
            else:
                logl += -0.5 * np.sum(rr**2, axis=1)
            logl += -np.sum(np.log(self.rv.v_err)) - 0.5 * self.rv.n * np.log(2 * np.pi)
        return logl

    def log_prob(self, theta) -> np.ndarray:
        theta = np.atleast_2d(theta)
        inside = np.all((theta > self.lo) & (theta < self.hi), axis=1)
        out = np.full(theta.shape[0], -np.inf)
        if np.any(inside):
            with np.errstate(all="ignore"):
                ll = self.log_likelihood(theta[inside])
            for name, (mu, sig) in self.gaussian_priors.items():
                ll = ll - 0.5 * ((theta[inside, self.free.index(name)] - mu) / sig) ** 2
            out[inside] = np.where(np.isfinite(ll), ll, -np.inf)
        return out

    def chi2(self, theta) -> dict:
        """Plain chi-square pieces at one parameter vector (published errors, unscaled)."""
        p = self.unpack(theta)
        m = self.model_observables(p)
        out = {}
        for fr in self.frames:
            rx = (fr.data.x - m[f"x_{fr.key}"][0]) / fr.data.x_err
            ry = (fr.data.y - m[f"y_{fr.key}"][0]) / fr.data.y_err
            out[f"astro_{fr.key}"] = float(np.sum(rx**2 + ry**2))
            out[f"n_astro_{fr.key}"] = 2 * fr.data.n
        if self.rv is not None:
            rr = (self.rv.v - m["rv"][0]) / self.rv.v_err
            out["rv"] = float(np.sum(rr**2))
            out["n_rv"] = self.rv.n
        return out


def _correlated_loglike(rx, ry, dist, length, mix) -> np.ndarray:
    """Do et al. (2019) eq. S42: rho_ij = (1 - p) delta_ij + p exp(-d_ij / Lambda).

    ``rx``, ``ry`` are residuals already divided by each point's error, so the
    covariance of the normalised residuals is just rho. Same rho for x and y.
    Returns log L without the constant and the sum of log sigma.
    """
    n = dist.shape[0]
    rho = (1.0 - mix)[:, None, None] * np.eye(n) + mix[:, None, None] * np.exp(
        -dist[None, :, :] / length[:, None, None])
    chol = np.linalg.cholesky(rho)
    rhs = np.stack([rx, ry], axis=-1)  # (W, n, 2)
    z = np.linalg.solve(chol, rhs)
    logdet = 2.0 * np.sum(np.log(np.diagonal(chol, axis1=1, axis2=2)), axis=1)
    return -0.5 * np.sum(z**2, axis=(1, 2)) - logdet  # logdet counted once each for x and y


# --------------------------------------------------------------------------
# Finding the best fit, then sampling around it
# --------------------------------------------------------------------------


def best_fit(model: Model, start: dict, *, n_restarts: int = 3, seed: int = 1) -> dict:
    """Maximise the posterior (flat priors: the likelihood) from ``start``.

    Nelder-Mead first (robust to the orbit's strong correlations), then
    Powell, repeated from the previous optimum until the value stops improving.
    """
    x = model.pack(start)
    scale = np.maximum(np.abs(x) * 1e-3, 1e-3)

    def neg(z):
        v = model.log_prob((x0 + z * scale)[None, :])[0]
        return 1e30 if not np.isfinite(v) else -v

    best = -np.inf
    x0 = x.copy()
    for _ in range(n_restarts):
        res = optimize.minimize(neg, np.zeros_like(x0), method="Nelder-Mead",
                                options={"maxiter": 40_000, "maxfev": 40_000, "xatol": 1e-6,
                                         "fatol": 1e-7, "adaptive": True})
        x0 = x0 + res.x * scale
        res = optimize.minimize(neg, np.zeros_like(x0), method="Powell",
                                options={"maxiter": 40_000, "xtol": 1e-7, "ftol": 1e-10})
        x0 = x0 + res.x * scale
        val = model.log_prob(x0[None, :])[0]
        if val - best < 1e-3:
            best = max(best, val)
            break
        best = val
    return {"theta": x0, "log_prob": float(best), "values": dict(zip(model.free, x0))}


def numerical_hessian_cov(model: Model, theta: np.ndarray) -> np.ndarray:
    """Covariance from the curvature of log L at the optimum (the Laplace
    approximation), used only to start the walkers in a sensible ball."""
    d = theta.size
    h = np.maximum(np.abs(theta) * 1e-5, 1e-6)
    f0 = model.log_prob(theta[None, :])[0]
    hess = np.zeros((d, d))
    for i in range(d):
        for j in range(i, d):
            pts = []
            for si, sj in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
                t = theta.copy()
                t[i] += si * h[i]
                t[j] += sj * h[j]
                pts.append(t)
            if i == j:
                tp, tm = theta.copy(), theta.copy()
                tp[i] += h[i]
                tm[i] -= h[i]
                fp, fm = model.log_prob(np.vstack([tp, tm]))
                hess[i, i] = (fp - 2 * f0 + fm) / h[i] ** 2
            else:
                f = model.log_prob(np.vstack(pts))
                hess[i, j] = hess[j, i] = (f[0] - f[1] - f[2] + f[3]) / (4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(-hess)
        np.linalg.cholesky(cov)
        return cov
    except np.linalg.LinAlgError:
        return np.diag(np.maximum(-1.0 / np.minimum(np.diag(hess), -1e-12), 1e-12))


def run_mcmc(model: Model, start_theta: np.ndarray, *, n_walkers: int | None = None,
             n_burn: int = 4_000, n_steps: int = 12_000, thin: int = 10, seed: int = 42,
             progress: bool = False) -> dict:
    """Affine-invariant ensemble MCMC (emcee), vectorised over walkers.

    Walkers start in a tight Gaussian ball drawn from the Laplace covariance at
    the best fit; the burn-in is discarded. Returns the flattened, thinned
    chain, the log-probabilities, the integrated autocorrelation times and the
    acceptance fraction, so convergence can be checked and reported.
    """
    import emcee

    rng = np.random.default_rng(seed)
    d = model.ndim
    n_walkers = n_walkers or max(4 * d, 48)
    cov = numerical_hessian_cov(model, start_theta)
    p0 = rng.multivariate_normal(start_theta, 0.1 * cov, size=n_walkers)
    p0 = np.clip(p0, model.lo + 1e-9 * (model.hi - model.lo), model.hi - 1e-9 * (model.hi - model.lo))
    bad = ~np.isfinite(model.log_prob(p0))
    while np.any(bad):
        p0[bad] = rng.multivariate_normal(start_theta, 0.01 * cov, size=int(bad.sum()))
        bad = ~np.isfinite(model.log_prob(p0))

    sampler = emcee.EnsembleSampler(n_walkers, d, model.log_prob, vectorize=True,
                                    moves=[(emcee.moves.DEMove(), 0.8),
                                           (emcee.moves.DESnookerMove(), 0.2)])
    sampler.random_state = np.random.RandomState(seed).get_state()
    # Burn-in is not stored at all; the production run keeps every thin-th step,
    # which keeps memory small (the machine is shared) without losing information:
    # consecutive steps are strongly correlated anyway.
    state = sampler.run_mcmc(p0, n_burn, progress=progress, store=False)
    sampler.reset()
    sampler.run_mcmc(state, n_steps // thin, progress=progress, thin_by=thin)
    chain = sampler.get_chain(flat=True)
    logp = sampler.get_log_prob(flat=True)
    try:
        tau = sampler.get_autocorr_time(quiet=True) * thin  # in steps
    except Exception:  # noqa: BLE001 -- emcee raises a plain Exception subclass
        tau = np.full(d, np.nan)
    return {
        "names": list(model.free),
        "chain": chain,
        "log_prob": logp,
        "tau": tau,
        "n_steps": n_steps,
        "n_walkers": n_walkers,
        "acceptance": float(np.mean(sampler.acceptance_fraction)),
    }


def summarize(samples: np.ndarray, names: list[str], weights=None) -> dict:
    """Median and 68% interval (and mean/sd) of every parameter."""
    out = {}
    for k, n in enumerate(names):
        x = samples[:, k]
        if weights is None:
            q = np.percentile(x, [2.5, 16, 50, 84, 97.5])
            mean, sd = float(np.mean(x)), float(np.std(x))
        else:
            o = np.argsort(x)
            cw = np.cumsum(weights[o])
            cw /= cw[-1]
            q = np.interp([0.025, 0.16, 0.5, 0.84, 0.975], cw, x[o])
            mean = float(np.sum(weights * x))
            sd = float(np.sqrt(np.sum(weights * (x - mean) ** 2)))
        out[n] = {"median": float(q[2]), "lo68": float(q[1]), "hi68": float(q[3]),
                  "lo95": float(q[0]), "hi95": float(q[4]), "mean": mean, "sd": sd}
    return out
