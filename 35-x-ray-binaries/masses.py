"""From a velocity curve to a black hole mass, and from masses to a population.

The physics, all in one place:

* the **mass function** ``f(M) = P K^3 (1 - e^2)^{3/2} / (2 pi G)``, the one
  number a velocity curve alone gives;
* the **mass** ``M1 = f (1 + q)^2 / sin^3 i`` once the mass ratio
  ``q = M2/M1`` and the inclination ``i`` are known (or for a massive
  companion, the root of ``M^3 sin^3 i / (M + M2)^2 = f``);
* Monte Carlo draws of (f, q, i) from their measurement distributions, which
  turn an orbit into a probability distribution for the mass ("marginalising
  over inclination");
* Eggleton's Roche-lobe radius and the highest inclination without eclipses;
* the hierarchical population models of Farr et al. (2011): power law,
  decaying exponential and Gaussian, with their priors, and the 1% quantile
  M_1% that defines "the lightest black hole".
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy.special import ndtr, ndtri

from config import CACHE_DIR, DAY_S, GM_SUN, POP_ALPHA_RANGE, POP_MASS_MAX

# ---------------------------------------------------------------------------
# Orbits
# ---------------------------------------------------------------------------


def mass_function(porb_d, k_kms, e=0.0):
    """f(M) in solar masses from the period (days), the companion's velocity
    semi-amplitude (km/s) and the eccentricity."""
    porb_d, k_kms, e = np.asarray(porb_d), np.asarray(k_kms), np.asarray(e)
    return porb_d * DAY_S * (k_kms * 1e3) ** 3 * (1 - e**2) ** 1.5 / (2 * np.pi * GM_SUN)


def m1_from_fm(fm, q, incl_deg):
    """Mass of the unseen object: M1 = f (1 + q)^2 / sin^3 i, with q = M2/M1."""
    return fm * (1 + q) ** 2 / np.sin(np.radians(incl_deg)) ** 3


def solve_m_bh(fm, m_companion, incl_deg, iterations: int = 80):
    """Unseen mass M from f = M^3 sin^3 i / (M + Mc)^2 when the companion mass
    Mc is known (massive companions, where q = Mc/M is not a separate
    measurement). Vectorised bisection; the root is unique for M > 0."""
    fm, mc, s3 = np.broadcast_arrays(np.asarray(fm, float), np.asarray(m_companion, float),
                                     np.sin(np.radians(incl_deg)) ** 3)
    lo = np.zeros_like(fm)
    hi = np.full_like(fm, 1e4)
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        positive = mid**3 * s3 - fm * (mid + mc) ** 2 > 0
        hi = np.where(positive, mid, hi)
        lo = np.where(positive, lo, mid)
    return 0.5 * (lo + hi)


def roche_lobe_fraction(q):
    """Eggleton (1983): volume-equivalent Roche-lobe radius over separation for
    the star of mass M2, with q = M2/M1."""
    q13 = np.cbrt(np.asarray(q, float))
    q23 = q13**2
    return 0.49 * q23 / (0.6 * q23 + np.log1p(q13))


def max_inclination_no_eclipse(q):
    """Highest inclination (degrees) at which a Roche-lobe-filling donor does not
    eclipse a point-like compact object: cos(i_max) = R_L/a."""
    return np.degrees(np.arccos(roche_lobe_fraction(q)))


def astrometric_mass_function(a0_mas, parallax_mas, period_d):
    """(a0/parallax)^3 / P^2 in AU, years and solar masses: equals
    M2^3/(M1 + M2)^2 for a dark companion (a0 is then the star's own orbit)."""
    return (a0_mas / parallax_mas) ** 3 / (period_d / 365.25) ** 2


def solve_m2_astrometric(f_ast, m1, iterations: int = 80):
    """Dark companion mass M2 from f_ast = M2^3/(M1 + M2)^2 (bisection)."""
    return solve_m_bh(f_ast, m1, 90.0, iterations)


def veiling_mass_bias(disc_fraction):
    """Mass overestimate if disc light is ignored in an ellipsoidal fit.

    The ellipsoidal amplitude scales roughly as sin^2 i. If a fraction ``d`` of
    the optical light comes from the disc but the model assumes none, the
    diluted amplitude is matched by sin^2 i_fit = (1 - d) sin^2 i_true, so the
    mass, proportional to 1/sin^3 i, comes out (1 - d)^(-3/2) times too large.
    A first-order illustration of Kreidberg et al.'s (2012) point, not a light
    curve model."""
    return (1.0 - np.asarray(disc_fraction)) ** -1.5


# ---------------------------------------------------------------------------
# Measurement distributions
# ---------------------------------------------------------------------------


def draw(spec, n: int, rng: np.random.Generator, *, lo: float = -np.inf, hi: float = np.inf):
    """Draw ``n`` values from a measurement distribution, truncated to (lo, hi].

    spec: ("N", mu, sigma) | ("SN", median, plus, minus) | ("U", a, b) |
          ("I", a_deg, b_deg) [cos i uniform] | ("fixed", value)
    """
    kind = spec[0]

    def one(m):
        if kind == "N":
            return rng.normal(spec[1], spec[2], m)
        if kind == "SN":
            u = rng.standard_normal(m)
            return spec[1] + np.where(u > 0, spec[2] * u, spec[3] * u)
        if kind == "U":
            return rng.uniform(spec[1], spec[2], m)
        if kind == "I":
            c = rng.uniform(np.cos(np.radians(spec[2])), np.cos(np.radians(spec[1])), m)
            return np.degrees(np.arccos(c))
        if kind == "fixed":
            return np.full(m, float(spec[1]))
        raise ValueError(f"unknown distribution {spec!r}")

    out = one(n)
    bad = ~((out > lo) & (out <= hi))
    for _ in range(200):
        if not bad.any():
            break
        out[bad] = one(int(bad.sum()))
        bad = ~((out > lo) & (out <= hi))
    if bad.any():
        raise ValueError(f"cannot draw {spec!r} inside ({lo}, {hi}]")
    return out


def spec_summary(spec) -> str:
    """Human-readable form of a distribution, for tables and the README."""
    kind = spec[0]
    if kind == "N":
        return f"N({spec[1]:g}, {spec[2]:g})"
    if kind == "SN":
        return f"{spec[1]:g} +{spec[2]:g} -{spec[3]:g}"
    if kind == "U":
        return f"U({spec[1]:g}, {spec[2]:g})"
    if kind == "I":
        return f"I({spec[1]:g}, {spec[2]:g})"
    return f"{spec[1]:g}"


def sample_lmxb_mass(fm_spec, q_spec, i_spec, n, rng):
    """Mass samples for a Roche-lobe-overflow system from (f, q, i) distributions
    (Farr et al. 2011, Sect. 2). ``fm_spec`` may also be an array of f draws."""
    fm = fm_spec if isinstance(fm_spec, np.ndarray) else draw(fm_spec, n, rng, lo=0.0)
    q = draw(q_spec, n, rng, lo=0.0)
    incl = draw(i_spec, n, rng, lo=0.0, hi=90.0)
    return m1_from_fm(fm, q, incl)


def summarise(samples, *, probs=(0.16, 0.5, 0.84)) -> dict:
    lo, med, hi = np.quantile(samples, probs)
    return {"median": float(med), "plus": float(hi - med), "minus": float(med - lo),
            "p05": float(np.quantile(samples, 0.05)), "p95": float(np.quantile(samples, 0.95))}


# ---------------------------------------------------------------------------
# Population models (Farr et al. 2011, Sect. 3.3)
# ---------------------------------------------------------------------------


def _triangle_transform(u):
    """Unit square -> uniform on x >= 0, y > 0, x + 2y <= 40 (the priors of the
    exponential and Gaussian models: most of the mass below 40 M_sun)."""
    x = POP_MASS_MAX * (1.0 - np.sqrt(1.0 - u[0]))
    return np.array([x, u[1] * (POP_MASS_MAX - x) / 2.0])


class PowerLaw:
    """p(M) = A M^alpha on [M_min, M_max]; uniform priors 0 <= M_min <= M_max <= 40,
    -15 <= alpha <= 13."""

    name = "power_law"
    labels = ("M_min", "M_max", "alpha")

    @staticmethod
    def log_prior(t):
        mmin, mmax, alpha = t[:, 0], t[:, 1], t[:, 2]
        ok = (mmin >= 0) & (mmin <= mmax) & (mmax <= POP_MASS_MAX) & \
             (alpha >= POP_ALPHA_RANGE[0]) & (alpha <= POP_ALPHA_RANGE[1])
        return np.where(ok, 0.0, -np.inf)

    @staticmethod
    def pdf(m, t):
        mmin, mmax, alpha = (t[:, k:k + 1] for k in range(3))
        a1 = 1.0 + alpha
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            norm = np.where(np.abs(a1) > 1e-8,
                            a1 / (mmax**a1 - mmin**a1),
                            1.0 / np.log(mmax / mmin))
            p = norm * m[None, :] ** alpha
        inside = (m[None, :] >= mmin) & (m[None, :] <= mmax)
        return np.where(inside & np.isfinite(p), p, 0.0)

    @staticmethod
    def cdf(m, t):
        mmin, mmax, alpha = (t[:, k:k + 1] for k in range(3))
        a1 = 1.0 + alpha
        mc = np.clip(m[None, :], mmin, mmax)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            f = np.where(np.abs(a1) > 1e-8,
                         (mc**a1 - mmin**a1) / (mmax**a1 - mmin**a1),
                         np.log(mc / mmin) / np.log(mmax / mmin))
        return f  # NaN where the model cannot be normalised (e.g. M_min = 0, alpha <= -1)

    @staticmethod
    def quantile(t, prob):
        mmin, mmax, alpha = t[:, 0], t[:, 1], t[:, 2]
        a1 = 1.0 + alpha
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            q = np.where(np.abs(a1) > 1e-8,
                         (mmin**a1 + prob * (mmax**a1 - mmin**a1)) ** (1.0 / a1),
                         mmin * (mmax / mmin) ** prob)
        return q

    @staticmethod
    def prior_transform(u):
        """Unit cube -> uniform on 0 <= M_min <= M_max <= 40, alpha in range."""
        mmax = POP_MASS_MAX * np.sqrt(u[0])
        return np.array([u[1] * mmax, mmax,
                         POP_ALPHA_RANGE[0] + (POP_ALPHA_RANGE[1] - POP_ALPHA_RANGE[0]) * u[2]])


class Exponential:
    """p(M) = exp(-(M - M_min)/M0)/M0 for M >= M_min (Fryer & Kalogera 2001);
    prior uniform with M_min + 2 M0 <= 40."""

    name = "exponential"
    labels = ("M_min", "M0")

    @staticmethod
    def log_prior(t):
        mmin, m0 = t[:, 0], t[:, 1]
        ok = (mmin >= 0) & (m0 > 0) & (mmin + 2 * m0 <= POP_MASS_MAX)
        return np.where(ok, 0.0, -np.inf)

    @staticmethod
    def pdf(m, t):
        mmin, m0 = t[:, 0:1], t[:, 1:2]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            p = np.exp(-(m[None, :] - mmin) / m0) / m0
        return np.where((m[None, :] >= mmin) & np.isfinite(p), p, 0.0)

    @staticmethod
    def cdf(m, t):
        mmin, m0 = t[:, 0:1], t[:, 1:2]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            return -np.expm1(-np.clip(m[None, :] - mmin, 0, None) / m0)

    @staticmethod
    def quantile(t, prob):
        return t[:, 0] - t[:, 1] * np.log1p(-prob)

    @staticmethod
    def prior_transform(u):
        return _triangle_transform(u)


class Gaussian:
    """p(M) = N(mu, sigma); prior 0 <= mu <= 40, sigma > 0, mu + 2 sigma <= 40."""

    name = "gaussian"
    labels = ("mu", "sigma")

    @staticmethod
    def log_prior(t):
        mu, sig = t[:, 0], t[:, 1]
        ok = (mu >= 0) & (mu <= POP_MASS_MAX) & (sig > 0) & (mu + 2 * sig <= POP_MASS_MAX)
        return np.where(ok, 0.0, -np.inf)

    @staticmethod
    def pdf(m, t):
        mu, sig = t[:, 0:1], t[:, 1:2]
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            p = np.exp(-0.5 * ((m[None, :] - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        return np.where(np.isfinite(p), p, 0.0)

    @staticmethod
    def cdf(m, t):
        mu, sig = t[:, 0:1], t[:, 1:2]
        with np.errstate(divide="ignore", invalid="ignore"):
            return ndtr((m[None, :] - mu) / sig)

    @staticmethod
    def quantile(t, prob):
        return t[:, 0] + t[:, 1] * ndtri(prob)

    @staticmethod
    def prior_transform(u):
        return _triangle_transform(u)


MODELS = {cls.name: cls for cls in (PowerLaw, Exponential, Gaussian)}


#: Log-spaced mass grid on which each system's samples are histogrammed for
#: the population likelihood: 0.25% wide bins from 0.3 to 300 M_sun.
GRID_EDGES = np.geomspace(0.3, 300.0, 2766)
GRID_CENTRES = np.sqrt(GRID_EDGES[:-1] * GRID_EDGES[1:])


class PopulationLikelihood:
    """log p(d | theta) = sum_i log( integral p_i(M) p(M | theta) dM ), Farr et
    al. (2011) Eq. (5): each system's Monte Carlo mass samples stand in for its
    measurement p_i(M). No selection effects are modelled (neither did Farr et
    al.).

    Farr et al. approximate the integral by averaging p(M | theta) over the
    samples (their Eq. 6). Here each system's samples are histogrammed once on
    a fine log-mass grid (0.25% bins) and the model is integrated exactly over
    every bin through its cumulative distribution:
    integral ~ sum_b [w_ib / width_b] [F(e_b+1) - F(e_b)]. With hundreds of
    thousands of samples this is the many-sample limit of the same estimator
    at the cost of a few thousand bins, and it stays finite when a model
    becomes narrower than a bin (an exponential or Gaussian collapsing onto
    one mass), where evaluating p(M | theta) at points would diverge."""

    def __init__(self, samples: dict[str, np.ndarray], model):
        self.model = model
        self.names = list(samples)
        counts = np.array([np.histogram(np.asarray(samples[k], float), GRID_EDGES)[0]
                           / len(samples[k]) for k in self.names])  # (systems, bins)
        self.density = counts / np.diff(GRID_EDGES)[None, :]

    def __call__(self, theta):
        theta = np.atleast_2d(theta)
        lp = self.model.log_prior(theta)
        out = np.full(len(theta), -np.inf)
        ok = np.isfinite(lp)
        if ok.any():
            prob = np.diff(self.model.cdf(GRID_EDGES, theta[ok]), axis=1)  # (walkers, bins)
            bad = ~np.isfinite(prob).all(axis=1)
            prob = np.where(np.isfinite(prob), prob, 0.0)
            per_system = prob @ self.density.T  # (walkers, systems)
            with np.errstate(divide="ignore"):
                ll = np.log(per_system).sum(axis=1) + lp[ok]
            out[ok] = np.where(bad, -np.inf, ll)
        return out


def _cache_key(samples: dict, model_name: str, extra: dict) -> str:
    h = hashlib.sha256()
    for k in sorted(samples):
        h.update(k.encode())
        h.update(np.ascontiguousarray(np.asarray(samples[k], np.float64)).tobytes())
    h.update(model_name.encode())
    h.update(json.dumps(extra, sort_keys=True).encode())
    return h.hexdigest()[:16]


def fit_population(samples: dict[str, np.ndarray], model_name: str, *, label: str,
                   nlive: int, seed: int, quiet: bool = False) -> dict:
    """Nested sampling (dynesty) for one population model; results cached.

    Nested sampling rather than MCMC because the power-law posterior has two
    branches (falling slopes with M_min pinned at the lightest systems, rising
    slopes with M_min free) whose relative weight an ensemble MCMC sampler
    gets wrong; it also returns the evidence used to compare models."""
    import dynesty
    from dynesty.utils import resample_equal

    model = MODELS[model_name]
    settings = {"nlive": nlive, "seed": seed, "sampler": "dynesty-rwalk", "dlogz": 0.01,
                "likelihood": "bin-integrated"}
    key = _cache_key(samples, model_name, settings)
    path = CACHE_DIR / f"pop_{label}_{model_name}_{key}.npz"
    if path.exists():
        z = np.load(path)
        post, logz, logz_err = z["post"], float(z["logz"]), float(z["logz_err"])
    else:
        like = PopulationLikelihood(samples, model)

        def loglike(t):
            v = like(t[None, :])[0]
            return v if np.isfinite(v) else -1e300

        ndim = len(model.labels)
        sampler = dynesty.NestedSampler(loglike, model.prior_transform, ndim, nlive=nlive,
                                        sample="rwalk", rstate=np.random.default_rng(seed))
        sampler.run_nested(dlogz=0.01, print_progress=False)
        res = sampler.results
        w = np.exp(res.logwt - res.logz[-1])
        post = resample_equal(res.samples, w / w.sum(), rstate=np.random.default_rng(seed + 1))
        logz, logz_err = float(res.logz[-1]), float(res.logzerr[-1])
        np.savez_compressed(path, post=post, logz=logz, logz_err=logz_err)
        if not quiet:
            print(f"    {label} {model_name}: ln Z = {logz:.2f}")
    m1 = model.quantile(post, 0.01)
    return {
        "model": model_name,
        "n_systems": len(samples),
        "systems": list(samples),
        "params": {lab: summarise(post[:, k]) for k, lab in enumerate(model.labels)},
        "m1pct_q10_q50_q90": [float(x) for x in np.quantile(m1, [0.1, 0.5, 0.9])],
        "ln_evidence": logz, "ln_evidence_err": logz_err,
        "_chain": post,
        "_m1pct": m1,
    }


def posterior_predictive(chain, model_name: str, grid, *, max_draws: int = 2000, seed: int = 0):
    """<p(M)> averaged over the posterior (Farr et al. Eq. 15), with 68% band."""
    model = MODELS[model_name]
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(chain), size=min(max_draws, len(chain)), replace=False)
    p = model.pdf(np.asarray(grid, float), chain[idx])
    return p.mean(axis=0), np.quantile(p, 0.16, axis=0), np.quantile(p, 0.84, axis=0)
