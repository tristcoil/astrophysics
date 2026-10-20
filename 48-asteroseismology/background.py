"""Granulation background + oscillation envelope: the fit that gives numax.

The power density of a red giant is modelled (after Kallinger et al. 2014,
A&A 570, A41) as

    P(nu) = T(nu) * [ W + eta^2(nu) * ( S(nu) + G1(nu) + G2(nu)
                                        + H exp(-(nu-numax)^2 / 2 sigma^2) ) ]

* G_i = xi a_i^2 / b_i / (1 + (nu/b_i)^4), xi = 2 sqrt(2) / pi: two
  granulation components ("super-Lorentzians"); a_i is the rms brightness
  fluctuation in ppm and b_i the characteristic frequency;
* S = (2/pi) a0^2 / b0 / (1 + (nu/b0)^2): a slow component for whatever
  varies on timescales of days to weeks besides granulation (starspots
  rotating, activity, residual instrumental drifts). Without it, a star with
  strong slow variability spends one granulation component on it and the
  envelope fit goes wrong;
* W: white (photon) noise;
* the Gaussian: the hump of oscillation power, centred on numax;
* eta^2 = sinc^2(pi nu / 2 nu_Nyq): the attenuation from integrating over each
  cadence (``spectrum.apodization``);
* T: the transfer function of our own detrending filter
  (``lightcurves.filter_transfer``), ~1 above 1 uHz.

Each component's normalisation makes a_i its rms amplitude: the integral of
the component over frequency is a_i^2.

The power in each frequency bin of a stochastic signal scatters like
chi-squared with two degrees of freedom around the model, so the fit
maximises the matching likelihood. To keep it fast the spectrum is first
averaged into narrow logarithmic bins; a mean of n independent chi^2_2 bins
follows a gamma distribution, so the log-likelihood is

    ln L = - sum_k n_k [ ln M_k + Pbar_k / M_k ]   (+ const).

The search for numax is blind: the fit is started from trial values of numax
between 3 and 270 uHz (other starting values follow from Kallinger et al.'s
scaling relations) and the best likelihood wins. Its significance is the
likelihood gain over the best fit without an envelope. Posteriors come from
emcee.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

import config as C
from lightcurves import filter_transfer
from spectrum import Spectrum, apodization, log_bin

XI = 2.0 * np.sqrt(2.0) / np.pi
PARAM_NAMES = ("W", "a0", "b0", "a1", "b1", "a2", "b2", "H", "numax", "sigma")
IDX = {n: i for i, n in enumerate(PARAM_NAMES)}


@dataclass
class BinnedSpectrum:
    nu: np.ndarray
    p: np.ndarray
    n: np.ndarray
    eta2: np.ndarray
    transfer: np.ndarray
    nyquist: float
    log_lo: np.ndarray = None  # ln of the parameter bounds, filled by bin_for_fit
    log_hi: np.ndarray = None


def _bounds(nu_lo: float, nu_hi: float) -> dict:
    return {
        "W": (1e-4, 1e7), "a0": (1e-2, 1e5), "b0": (1e-2, nu_hi),
        "a1": (1.0, 1e5), "b1": (0.2, 3 * nu_hi), "a2": (1.0, 1e5), "b2": (0.2, 3 * nu_hi),
        "H": (1e-4, 1e9), "numax": (nu_lo, nu_hi), "sigma": (0.02, 0.6 * nu_hi),
    }


def bin_for_fit(spec: Spectrum, *, trend_sigma_days: float, nu_min: float = C.FIT_NU_MIN_UHZ,
                nu_max: float | None = None, dlnnu: float = C.FIT_LOG_BIN) -> BinnedSpectrum:
    hi = spec.nyquist if nu_max is None else nu_max
    lo = max(nu_min, spec.nu[0])
    nu, p, n = log_bin(spec.nu, spec.psd, nu_min=lo, nu_max=hi, dlnnu=dlnnu)
    # with an oversampled grid, neighbouring bins are not independent
    n = np.maximum(n / spec.oversample, 1.0 / spec.oversample)
    bs = BinnedSpectrum(nu=nu, p=p, n=n, eta2=apodization(nu, spec.nyquist),
                        transfer=filter_transfer(nu, trend_sigma_days), nyquist=spec.nyquist)
    b = _bounds(nu.min(), nu.max())
    bs.log_lo = np.log([b[k][0] for k in PARAM_NAMES])
    bs.log_hi = np.log([b[k][1] for k in PARAM_NAMES])
    return bs


def harvey(nu, a, b):
    """Super-Lorentzian with exponent 4; integrates to a^2."""
    return XI * a ** 2 / b / (1.0 + (nu / b) ** 4)


def lorentzian(nu, a, b):
    """Lorentzian (exponent 2); integrates to a^2."""
    return (2.0 / np.pi) * a ** 2 / b / (1.0 + (nu / b) ** 2)


def components(theta, nu, eta2, transfer) -> dict:
    W, a0, b0, a1, b1, a2, b2, H, numax, sig = theta
    k = transfer * eta2
    return {"white": transfer * W, "slow": k * lorentzian(nu, a0, b0),
            "gran1": k * harvey(nu, a1, b1), "gran2": k * harvey(nu, a2, b2),
            "envelope": k * H * np.exp(-0.5 * ((nu - numax) / sig) ** 2)}


def model(theta, nu, eta2, transfer):
    return sum(components(theta, nu, eta2, transfer).values())


def background_only(theta, nu, eta2, transfer):
    c = components(theta, nu, eta2, transfer)
    return c["white"] + c["slow"] + c["gran1"] + c["gran2"]


def granulation_at(theta, nu):
    """Intrinsic granulation power density at nu (no apodization, no filter)."""
    return harvey(nu, theta[IDX["a1"]], theta[IDX["b1"]]) + harvey(nu, theta[IDX["a2"]], theta[IDX["b2"]])


def kallinger_guess(numax: float, bs: BinnedSpectrum) -> np.ndarray:
    """Starting values from empirical scalings with numax.

    Kallinger et al. 2014: b1 = 0.317 numax^0.970, b2 = 0.948 numax^0.992,
    a = 3382 numax^-0.609 ppm; Mosser et al. 2012: envelope FWHM =
    0.66 numax^0.88.
    """
    b1 = 0.317 * numax ** 0.970
    b2 = 0.948 * numax ** 0.992
    a = 3382.0 * numax ** -0.609
    sig = 0.66 * numax ** 0.88 / 2.3548
    top = bs.nu > 0.9 * bs.nu.max()
    W = max(np.median(bs.p[top]) * 0.8, 1e-3)
    b0 = max(0.05, min(1.0, b1 / 4.0))
    th = np.array([W, 0.3 * a, b0, a, b1, a, b2, 1.0, numax, sig])
    eta2_nm = float(np.interp(numax, bs.nu, bs.eta2))
    bg = float(background_only(th, np.array([numax]), np.array([eta2_nm]), np.array([1.0]))[0])
    near = np.abs(bs.nu - numax) < sig
    excess = np.median(bs.p[near]) - bg if near.any() else 0.0
    th[IDX["H"]] = max(excess / max(eta2_nm, 0.05), 0.05 * bg)
    lb = _bounds(bs.nu.min(), bs.nu.max())
    th = np.array([np.clip(v, lb[k][0] * 1.01, lb[k][1] * 0.99) for v, k in zip(th, PARAM_NAMES)])
    return th


def neg_lnlike(log_theta, bs: BinnedSpectrum, with_env: bool = True) -> float:
    th = np.exp(log_theta)
    if not with_env:
        th = th.copy()
        th[IDX["H"]] = 0.0
    m = model(th, bs.nu, bs.eta2, bs.transfer)
    if not np.all(np.isfinite(m)) or np.any(m <= 0):
        return 1e30
    return float(np.sum(bs.n * (np.log(m) + bs.p / m)))


def _in_bounds(log_theta, bs) -> bool:
    if np.any(log_theta < bs.log_lo) or np.any(log_theta > bs.log_hi):
        return False
    # b0 < b1 < b2; 0.02 numax < sigma < 0.5 numax
    return (log_theta[2] < log_theta[4] < log_theta[6]) and \
        (-3.912 < log_theta[9] - log_theta[8] < -0.693)


def _optimise(log0, bs, with_env=True, maxiter=6000):
    def f(x):
        if not _in_bounds(x, bs):
            return 1e30
        return neg_lnlike(x, bs, with_env)

    best = minimize(f, log0, method="Nelder-Mead",
                    options={"maxiter": maxiter, "maxfev": maxiter, "xatol": 1e-6, "fatol": 1e-6})
    # restart from the optimum: Nelder-Mead often stalls early in 10 dimensions
    for _ in range(2):
        again = minimize(f, best.x, method="Nelder-Mead",
                         options={"maxiter": maxiter, "maxfev": maxiter, "xatol": 1e-7, "fatol": 1e-7})
        if again.fun >= best.fun - 1e-3:
            break
        best = again
    return best


@dataclass
class EnvelopeFit:
    theta: np.ndarray  # linear parameters at maximum likelihood
    nll: float
    nll_no_env: float
    search: list = field(default_factory=list)  # (start numax, end numax, nll)
    samples: np.ndarray | None = None  # posterior samples, linear

    @property
    def delta_lnl(self) -> float:
        return self.nll_no_env - self.nll

    def p(self, name: str) -> float:
        return float(self.theta[IDX[name]])

    def summary(self) -> dict:
        out = {"ml": dict(zip(PARAM_NAMES, map(float, self.theta))),
               "delta_lnL_envelope": float(self.delta_lnl)}
        if self.samples is not None:
            q = np.percentile(self.samples, [16, 50, 84], axis=0)
            out["posterior"] = {n: [float(q[1, i]), float(q[1, i] - q[0, i]), float(q[2, i] - q[1, i])]
                                for i, n in enumerate(PARAM_NAMES)}
        return out


def fit_envelope(bs: BinnedSpectrum, *, search=C.NUMAX_SEARCH, numax_range=None) -> EnvelopeFit:
    """Blind multi-start maximum-likelihood fit of background + envelope."""
    lo, hi, n = search
    lo = max(lo, bs.nu.min() * 1.5)
    hi = min(hi, bs.nu.max() * 0.97)
    if numax_range is not None:
        lo, hi = max(lo, numax_range[0]), min(hi, numax_range[1])
    trials = np.exp(np.linspace(np.log(lo), np.log(hi), n))

    def f(x):
        if not _in_bounds(x, bs):
            return 1e30
        return neg_lnlike(x, bs)

    # stage 1: a quick fit from every starting value
    results = []
    for nm in trials:
        r = minimize(f, np.log(kallinger_guess(nm, bs)), method="Nelder-Mead",
                     options={"maxiter": 1500, "maxfev": 1500, "xatol": 1e-5, "fatol": 1e-4})
        results.append((nm, float(np.exp(r.x[IDX["numax"]])), float(r.fun), r.x))
    # stage 2: converge the best few distinct solutions properly
    ranked = sorted(results, key=lambda r: r[2])
    finalists = []
    for r in ranked:
        if all(abs(np.log(r[1] / q[1])) > 0.05 for q in finalists):
            finalists.append(r)
        if len(finalists) == 4:
            break
    polished = []
    for r in finalists:
        rr = _optimise(r[3], bs, maxiter=12000)
        polished.append((rr.x, float(rr.fun)))
    best_x, best_f = min(polished, key=lambda r: r[1])
    # background-only fit, started from the best fit's background
    r0 = _optimise(best_x, bs, with_env=False)
    return EnvelopeFit(theta=np.exp(best_x), nll=best_f, nll_no_env=float(r0.fun),
                       search=[(r[0], r[1], r[2]) for r in results])


def sample_posterior(fit: EnvelopeFit, bs: BinnedSpectrum, *, seed: int,
                     walkers: int = C.MCMC_WALKERS, steps: int = C.MCMC_STEPS,
                     burn: int = C.MCMC_BURN) -> np.ndarray:
    """emcee posterior (flat priors in ln of each parameter, within bounds)."""
    import emcee

    def lnprob(x):
        if not _in_bounds(x, bs):
            return -np.inf
        v = neg_lnlike(x, bs)
        return -v if v < 1e29 else -np.inf

    rng = np.random.default_rng(seed)
    x0 = np.log(fit.theta)
    p0 = x0 + 1e-4 * rng.standard_normal((walkers, len(x0)))
    for i in range(walkers):  # keep every walker inside the bounds
        while not np.isfinite(lnprob(p0[i])):
            p0[i] = x0 + 1e-5 * rng.standard_normal(len(x0))
    sampler = emcee.EnsembleSampler(walkers, len(x0), lnprob)
    # emcee copies numpy's *global* RNG state at construction; seed it here
    # explicitly or runs are not reproducible.
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(p0, steps, progress=False)
    chain = sampler.get_chain(discard=burn, flat=True)
    return np.exp(chain)
