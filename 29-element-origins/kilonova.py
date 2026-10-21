"""AT2017gfo: weigh the r-process elements made in GW170817's kilonova.

The model is the one Villar et al. (2017) fitted with MOSFiT, re-implemented
from their paper and MOSFiT's source:

* Heating (Korobkin et al. 2012): freshly made r-process nuclei decay and
  heat the ejecta at L_in = 4e18 M [0.5 - arctan((t - t0)/s)/pi]^1.3 erg/s
  (M in grams, t0 = 1.3 s, s = 0.11 s): ~constant for a second, then ~t^-1.3.
* Thermalisation (Barnes et al. 2016): only a fraction
  0.36 [exp(-a t) + ln(1 + 2 b t^d) / (2 b t^d)] of it heats the gas, falling
  as the ejecta thin out; a, b, d interpolated in their table of mass and
  velocity (as MOSFiT does).
* Diffusion (Arnett 1982): the heat leaks out over t_d = sqrt(2 kappa M /
  (beta v c)), beta = 13.7 (MOSFiT's value; the paper prints 13.4), times
  MOSFiT's gamma-ray trapping factor 1 - exp(-A/t^2), A = 3 kappa_g M /
  (4 pi v^2) with kappa_g = 100 cm^2/g as in its kilonova model.
* Photosphere: a blackbody expanding at v until it cools to a floor T_c, then
  receding at fixed temperature.
* Three independent components with fixed grey opacities: blue (lanthanide-
  poor, kappa = 0.5 cm^2/g), purple (3) and red (lanthanide-rich, 10).
  Free: mass, velocity and floor temperature of each, plus a scatter term
  added to every error bar -- ten parameters, MOSFiT's priors.

Fluxes are integrated over each filter's real transmission curve, with
Milky Way extinction (E(B-V) = 0.105, Cardelli et al. 1989 / O'Donnell
1994) and the redshift of NGC 4993; upper limits enter as one-sided
likelihoods. The posterior is sampled with emcee.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.special import log_ndtr

import kn_data
from config import KN_DISTANCE_MPC, KN_EBV_MW, KN_KAPPA, KN_REDSHIFT, KN_RV

C = 2.99792458e10
H = 6.62607015e-27
K_B = 1.380649e-16
SIGMA_SB = 5.670374419e-5
M_SUN = 1.98847e33
DAY = 86400.0
MPC = 3.0856775814913673e24
KM = 1e5

COMPONENTS = ("blue", "purple", "red")
PARAMS = ["m_blue", "v_blue", "T_blue", "m_purple", "v_purple", "T_purple",
          "m_red", "v_red", "T_red", "sigma"]
#: MOSFiT's priors for the three-component model (masses in Msun, speeds in
#: km/s, temperatures in K); "log" = uniform in the logarithm.
PRIORS = {
    "m_blue": (1e-3, 5e-2, "lin"), "v_blue": (6e4, 1.2e5, "lin"), "T_blue": (1e2, 1e4, "log"),
    "m_purple": (1e-3, 1e-1, "log"), "v_purple": (1e4, 1e5, "lin"), "T_purple": (1e2, 1e4, "log"),
    "m_red": (1e-3, 1e-1, "log"), "v_red": (1e4, 6e4, "lin"), "T_red": (1e2, 1e4, "log"),
    "sigma": (1e-3, 1e2, "log"),
}
UPPER_LIMIT_ERROR = 0.1      # MOSFiT's default error on an upper limit (mag)

# Barnes et al. (2016) thermalisation coefficients, as tabulated in MOSFiT.
_B_V = np.array([0.1, 0.2, 0.3, 0.4])
_B_M = np.array([1e-3, 5e-3, 1e-2, 5e-2, 1e-1])
_B_A = np.array([[2.01, 4.52, 8.16, 16.3], [0.81, 1.9, 3.2, 5.0], [0.56, 1.31, 2.19, 3.0],
                 [.27, .55, .95, 2.0], [0.20, 0.39, 0.65, 0.9]])
_B_B = np.array([[0.28, 0.62, 1.19, 2.4], [0.19, 0.28, 0.45, 0.65], [0.17, 0.21, 0.31, 0.45],
                 [0.10, 0.13, 0.15, 0.17], [0.06, 0.11, 0.12, 0.12]])
_B_D = np.array([[1.12, 1.39, 1.52, 1.65], [0.86, 1.21, 1.39, 1.5], [0.74, 1.13, 1.32, 1.4],
                 [0.6, 0.9, 1.13, 1.25], [0.63, 0.79, 1.04, 1.5]])
_THERM = [RegularGridInterpolator((_B_M, _B_V), t, bounds_error=False, fill_value=None)
          for t in (_B_A, _B_B, _B_D)]

#: Dense rest-frame time grid (days) for the light-curve integrals.
T_GRID = np.logspace(-4, np.log10(45.0), 700)
DIFF_BETA = 13.7
KAPPA_GAMMA = 100.0


def heating_rate(t_day, mass_msun):
    """Korobkin et al. (2012) radioactive power (erg/s), before thermalisation."""
    t = t_day * DAY
    return 4e18 * mass_msun * M_SUN * (0.5 - np.arctan((t - 1.3) / 0.11) / np.pi) ** 1.3


def thermalisation(t_day, mass_msun, v_kms):
    a, b, d = (float(f([[mass_msun, v_kms / (C / KM)]])[0]) for f in _THERM)
    bt = 2.0 * b * t_day ** d
    return 0.36 * (np.exp(-a * t_day) + np.log1p(bt) / bt)


def diffusion_time(kappa, mass_msun, v_kms):
    """Arnett diffusion time in days."""
    return np.sqrt(2.0 * kappa * mass_msun * M_SUN / (DIFF_BETA * C * v_kms * KM)) / DAY


def bolometric(t, mass_msun, v_kms, kappa):
    """Arnett (1982) output luminosity on the grid t (days):

        L(t) = (2/td^2) exp(-t^2/td^2) int_0^t Q(t') t' exp(t'^2/td^2) dt',

    Q = heating x thermalisation. With u = t'^2 the integrand becomes
    (Q/2) exp(-(u_t - u)/td^2) du; Q is smooth, so taking it linear in u over
    each grid step makes every step's integral exact, however long the step
    is compared with td (a plain trapezoid overshoots by ~15% late on)."""
    q = 0.5 * heating_rate(t, mass_msun) * thermalisation(t, mass_msun, v_kms)
    a = 1.0 / diffusion_time(kappa, mass_msun, v_kms) ** 2
    u = t ** 2
    du = np.diff(u)
    x = a * du
    decay = np.exp(-x)
    one_minus = -np.expm1(-x)
    # h(x) = 1 - e^-x (1 + x), by its series where it would cancel
    h = np.where(x < 1e-3, x ** 2 / 2 - x ** 3 / 3 + x ** 4 / 8, one_minus - x * decay)
    seg = q[1:] * one_minus / a + (q[:-1] - q[1:]) * h / (a ** 2 * du)
    out = np.zeros_like(t)
    acc = 0.0
    for k in range(1, t.size):                  # 700 steps
        acc = acc * decay[k - 1] + seg[k - 1]
        out[k] = acc
    trap_a = 3.0 * KAPPA_GAMMA * mass_msun * M_SUN / (4 * np.pi * (v_kms * KM) ** 2) / DAY ** 2
    return 2.0 * a * out * -np.expm1(-trap_a / t ** 2)


def photosphere(t, lum, v_kms, t_floor):
    """Radius (cm) and temperature (K) of the blackbody photosphere."""
    r_exp = v_kms * KM * t * DAY
    t_exp = (lum / (4 * np.pi * SIGMA_SB * r_exp ** 2)) ** 0.25
    hot = t_exp > t_floor
    radius = np.where(hot, r_exp, np.sqrt(lum / (4 * np.pi * SIGMA_SB * t_floor ** 4)))
    return radius, np.where(hot, t_exp, t_floor)


# ---------------------------------------------------------------------------
# Filters, extinction, blackbody band fluxes
# ---------------------------------------------------------------------------


def ccm_odonnell(lam_angstrom, r_v=KN_RV):
    """A_lambda / A_V of Cardelli, Clayton & Mathis (1989) with O'Donnell's
    (1994) optical coefficients; 0.3 <= x <= 10 per micron (zero beyond 3.3 um)."""
    x = 1e4 / np.asarray(lam_angstrom, float)
    a, b = np.zeros_like(x), np.zeros_like(x)
    ir = (x >= 0.3) & (x < 1.1)
    a[ir], b[ir] = 0.574 * x[ir] ** 1.61, -0.527 * x[ir] ** 1.61
    op = (x >= 1.1) & (x < 3.3)
    y = x[op] - 1.82
    a[op] = (1 + 0.104 * y - 0.609 * y ** 2 + 0.701 * y ** 3 + 1.137 * y ** 4 - 1.718 * y ** 5
             - 0.827 * y ** 6 + 1.647 * y ** 7 - 0.505 * y ** 8)
    b[op] = (1.952 * y + 2.908 * y ** 2 - 3.989 * y ** 3 - 7.985 * y ** 4 + 11.102 * y ** 5
             + 5.491 * y ** 6 - 10.805 * y ** 7 + 3.347 * y ** 8)
    uv = (x >= 3.3) & (x <= 8.0)
    xu = x[uv]
    fa = np.where(xu >= 5.9, -0.04473 * (xu - 5.9) ** 2 - 0.009779 * (xu - 5.9) ** 3, 0.0)
    fb = np.where(xu >= 5.9, 0.2130 * (xu - 5.9) ** 2 + 0.1207 * (xu - 5.9) ** 3, 0.0)
    a[uv] = 1.752 - 0.316 * xu - 0.104 / ((xu - 4.67) ** 2 + 0.341) + fa
    b[uv] = -3.090 + 1.825 * xu + 1.206 / ((xu - 4.62) ** 2 + 0.263) + fb
    fuv = x > 8.0
    xf = x[fuv] - 8.0
    a[fuv] = -1.073 - 0.628 * xf + 0.137 * xf ** 2 - 0.070 * xf ** 3
    b[fuv] = 13.670 + 4.257 * xf - 0.420 * xf ** 2 + 0.374 * xf ** 3
    return a + b / r_v


LOG_T_TABLE = np.linspace(1.5, 5.3, 600)


class Bands:
    """For each filter: the band-averaged f_nu of a unit-radius blackbody
    (pi B_nu, redshifted, extinguished) as a function of temperature, so a
    model magnitude is a table lookup."""

    def __init__(self, svo_ids, z=KN_REDSHIFT, ebv=KN_EBV_MW):
        self.ids = list(svo_ids)
        self.table = np.empty((len(self.ids), LOG_T_TABLE.size))
        self.lam_eff = np.empty(len(self.ids))
        T = 10 ** LOG_T_TABLE
        for i, sid in enumerate(self.ids):
            lam, tr = kn_data.read_filter(sid)
            nu_rest = C / (lam * 1e-8) * (1 + z)
            ext = 10 ** (-0.4 * KN_RV * ebv * ccm_odonnell(lam))
            x = np.clip(H * nu_rest[None, :] / (K_B * T[:, None]), 1e-8, 700)
            bnu = 2 * H * nu_rest ** 3 / C ** 2 / np.expm1(x)        # rest frame
            fnu = (1 + z) * np.pi * bnu * ext                            # observed, per cm^2 of R^2
            w = tr / lam                                                 # photon counting, AB
            self.table[i] = np.trapezoid(fnu * w, lam, axis=1) / np.trapezoid(w, lam)
            self.lam_eff[i] = np.trapezoid(tr * lam, lam) / np.trapezoid(tr, lam)
        self.log_table = np.log10(np.clip(self.table, 1e-300, None))

    def log_flux(self, idx, temp):
        """log10 of the table at each (filter index, temperature)."""
        lt = np.clip(np.log10(temp), LOG_T_TABLE[0], LOG_T_TABLE[-1])
        j = (lt - LOG_T_TABLE[0]) / (LOG_T_TABLE[1] - LOG_T_TABLE[0])
        j0 = np.clip(j.astype(int), 0, LOG_T_TABLE.size - 2)
        f = j - j0
        return self.log_table[idx, j0] * (1 - f) + self.log_table[idx, j0 + 1] * f


# ---------------------------------------------------------------------------
# The model and the likelihood
# ---------------------------------------------------------------------------


def component_curves(theta, t_grid=T_GRID):
    """Bolometric luminosity, photospheric radius and temperature of each
    component on the rest-frame time grid."""
    p = dict(zip(PARAMS, theta))
    out = {}
    for c in COMPONENTS:
        lum = bolometric(t_grid, p[f"m_{c}"], p[f"v_{c}"], KN_KAPPA[c])
        r, temp = photosphere(t_grid, np.clip(lum, 1e10, None), p[f"v_{c}"], p[f"T_{c}"])
        out[c] = dict(lum=lum, radius=r, temp=temp)
    return out


class Kilonova:
    """Villar et al. (2017)'s data with the three-component model."""

    def __init__(self, table=None, distance_mpc=KN_DISTANCE_MPC, z=KN_REDSHIFT):
        t = kn_data.read_villar_table() if table is None else table
        u = t["used"]
        self.t_obs = t["t"][u]
        self.t_rest = self.t_obs / (1 + z)
        self.mag = t["mag"][u]
        self.err = np.where(t["upper"][u], UPPER_LIMIT_ERROR, t["err"][u])
        self.upper = t["upper"][u]
        self.svo = t["svo"][u]
        self.band = t["band"][u]
        self.telescope = t["telescope"][u]
        self.bands = Bands(sorted(set(self.svo)), z=z)
        self.idx = np.array([self.bands.ids.index(s) for s in self.svo])
        self.log_d2 = np.log10((distance_mpc * MPC) ** 2)
        self.distance_mpc = distance_mpc

    def model_mag(self, theta, t_rest=None, idx=None):
        t_rest = self.t_rest if t_rest is None else t_rest
        idx = self.idx if idx is None else idx
        curves = component_curves(theta)
        flux = np.zeros(t_rest.size)
        for c in COMPONENTS:
            r = np.interp(t_rest, T_GRID, curves[c]["radius"])
            temp = np.interp(t_rest, T_GRID, curves[c]["temp"])
            flux += 10 ** (self.bands.log_flux(idx, temp) + 2 * np.log10(r) - self.log_d2)
        return -2.5 * np.log10(np.clip(flux, 1e-60, None)) - 48.6       # AB

    def log_prior(self, theta):
        lp = 0.0
        for name, val in zip(PARAMS, theta):
            lo, hi, kind = PRIORS[name]
            if not lo <= val <= hi:
                return -np.inf
            lp += -np.log(val) if kind == "log" else 0.0
        return lp

    def log_likelihood(self, theta):
        sigma = theta[-1]
        m = self.model_mag(theta)
        var = self.err ** 2 + sigma ** 2
        det = ~self.upper
        ll = -0.5 * np.sum((self.mag[det] - m[det]) ** 2 / var[det] + np.log(2 * np.pi * var[det]))
        # An upper limit: the source must be fainter than the limit.
        ll += np.sum(log_ndtr((m[~det] - self.mag[~det]) / np.sqrt(var[~det])))
        return ll

    def log_posterior(self, theta):
        lp = self.log_prior(theta)
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.log_likelihood(theta)


def find_map(kn: Kilonova, *, seed=0, n_starts=12):
    """Maximum a posteriori from several random starts across the prior
    (no published value is used as a starting point)."""
    from scipy.optimize import minimize

    rng = np.random.default_rng(seed)
    best = None
    for _ in range(n_starts):
        theta0 = []
        for name in PARAMS:
            lo, hi, kind = PRIORS[name]
            theta0.append(np.exp(rng.uniform(np.log(lo), np.log(hi))) if kind == "log"
                          else rng.uniform(lo, hi))
        theta0[-1] = 0.3
        x0 = np.log(theta0)
        r = minimize(lambda x: -kn.log_posterior(np.exp(x)), x0, method="Nelder-Mead",
                     options=dict(maxiter=6000, maxfev=6000, xatol=1e-4, fatol=1e-3))
        if best is None or r.fun < best.fun:
            best = r
    return np.exp(best.x), -best.fun


def sample(kn: Kilonova, theta_map, *, n_walkers=40, n_steps=4000, burn=2000, seed=1):
    """emcee in log-parameters around the MAP; returns the flattened chain
    (natural units) and the acceptance fraction."""
    import emcee

    rng = np.random.default_rng(seed)
    ndim = len(PARAMS)

    def lp(x):
        th = np.exp(x)
        v = kn.log_posterior(th)
        return v + np.sum(x) if np.isfinite(v) else -np.inf   # Jacobian of the log transform

    p0 = np.log(theta_map) + 1e-3 * rng.standard_normal((n_walkers, ndim))
    sampler = emcee.EnsembleSampler(n_walkers, ndim, lp)
    sampler.run_mcmc(p0, n_steps, progress=False)
    chain = np.exp(sampler.get_chain(discard=burn, flat=True))
    try:
        tau = sampler.get_autocorr_time(discard=burn, quiet=True)
    except Exception:
        tau = np.full(ndim, np.nan)
    return chain, float(np.mean(sampler.acceptance_fraction)), tau


def prior_transform(u):
    """Unit cube -> parameters under MOSFiT's priors (for nested sampling)."""
    out = np.empty_like(u)
    for i, name in enumerate(PARAMS):
        lo, hi, kind = PRIORS[name]
        out[i] = (np.exp(np.log(lo) + u[i] * np.log(hi / lo)) if kind == "log"
                  else lo + u[i] * (hi - lo))
    return out


def nested(kn: Kilonova, *, nlive=500, dlogz=0.1, seed=2):
    """Nested sampling (dynesty) of the full posterior. The purple and red
    components trade mass, so the posterior has separate modes; nested
    sampling finds them and weighs them by their evidence, which an MCMC
    started in one mode cannot. Returns equally weighted samples and log Z."""
    import dynesty

    rng = np.random.default_rng(seed)
    s = dynesty.NestedSampler(kn.log_likelihood, prior_transform, len(PARAMS), nlive=nlive,
                              bound="multi", sample="rwalk", rstate=rng)
    s.run_nested(dlogz=dlogz, print_progress=False)
    r = s.results
    samples = r.samples_equal(rstate=rng) if hasattr(r, "samples_equal") else \
        dynesty.utils.resample_equal(r.samples, np.exp(r.logwt - r.logz[-1]))
    return dict(samples=samples, logz=float(r.logz[-1]), logz_err=float(r.logzerr[-1]),
                n_calls=int(np.sum(r.ncall)), max_loglike=float(np.max(r.logl)),
                best=r.samples[int(np.argmax(r.logl))])
