"""The physics behind every FRB number in this project.

Three pieces:

1. **Dispersion.** Radio waves crossing ionized gas travel slower the lower their
   frequency. A burst emitted at one instant arrives spread out in time, low
   frequencies last, with a delay ``t = k_DM * DM / nu^2``. The dispersion
   measure ``DM = integral of n_e dl`` is the column density of free electrons
   between us and the source -- a direct count of all the electrons the signal
   passed on its way here.

2. **The Macquart relation.** Most of the universe's ordinary matter (baryons)
   sits in thin, hot, ionized gas between galaxies -- too diffuse and too hot to
   see in emission or absorption, which is why it was called "missing". Every
   one of its electrons disperses a passing FRB, so the mean DM grows with
   distance in a way fixed by the cosmic baryon density. Measuring DM against
   host-galaxy redshift weighs that gas directly.

3. **The scatter around it.** Sightlines differ: some cross galaxy halos, some
   thread voids. Macquart et al. (2020, Nature 581, 391) give a probability
   distribution for the cosmic DM along a random sightline, which together with
   a lognormal host-galaxy contribution turns the relation into a likelihood we
   can fit.

All DMs are in pc cm^-3, frequencies in MHz, times in seconds unless stated.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import Planck18
from scipy import integrate, interpolate, optimize

from config import F_DIFFUSE_FIDUCIAL, MACQUART_ALPHA, MACQUART_BETA, Y_HE

# ---------------------------------------------------------------------------
# 1. Dispersion
# ---------------------------------------------------------------------------

#: The dispersion constant e^2 / (2 pi m_e c), from first principles. The group
#: delay of a wave at frequency nu through a column DM of cold plasma is
#: K_DM * DM / nu^2. Equal to 4.148808 ms GHz^2 pc^-1 cm^3.
K_DM_EXACT = (
    const.e.gauss ** 2 / (2 * np.pi * const.m_e * const.c)
).to(u.s * u.MHz ** 2 * u.cm ** 3 / u.pc).value

#: The value pulsar and FRB astronomers actually use, 1/2.41e-4 s MHz^2 pc^-1
#: cm^3, a convention dating from Manchester & Taylor (1972). It is 0.014% larger
#: than the exact constant, so DMs quoted under it are 0.014% smaller than the
#: true electron column. CHIME's catalogue DMs use it, so every delay we compute
#: from a catalogue DM must use it too, or the two would disagree.
K_DM_CONVENTION = 1.0 / 2.41e-4


def dispersion_delay(dm, freq_mhz, ref_freq_mhz=np.inf, k=K_DM_CONVENTION):
    """Arrival delay (s) at ``freq_mhz`` relative to ``ref_freq_mhz``.

    Positive when ``freq_mhz`` is below the reference: low frequencies arrive
    later. With the default infinite reference this is the absolute delay
    relative to a signal that crossed no plasma at all.
    """
    freq = np.asarray(freq_mhz, dtype=float)
    return k * np.asarray(dm, dtype=float) * (freq ** -2.0 - np.asarray(ref_freq_mhz, float) ** -2.0)


# ---------------------------------------------------------------------------
# 2. The mean cosmic DM (the Macquart relation)
# ---------------------------------------------------------------------------


def electrons_per_baryon_mass(y_he: float = Y_HE) -> float:
    """Free electrons per proton mass of fully ionized primordial gas.

    Hydrogen (mass fraction 1-Y) gives one electron per m_p; helium (Y) gives two
    per 4 m_p. Hence (1 - Y) + Y/2 = 1 - Y/2 = 0.877 for Y = 0.2454.
    """
    return 1.0 - y_he / 2.0


def cosmic_dm_prefactor(cosmo=Planck18, y_he: float = Y_HE) -> float:
    """``3 c H0 Omega_b (1 - Y/2) / (8 pi G m_p)`` in pc cm^-3.

    This is c/H0 times the present-day mean electron density n_e,0 if *every*
    baryon were in diffuse ionized gas: the DM scale of the universe.
    ``mean_dm_cosmic`` multiplies it by f_d and a dimensionless integral.
    """
    h0 = cosmo.H0.to(1 / u.s)
    val = (
        3 * const.c * h0 * cosmo.Ob0 * electrons_per_baryon_mass(y_he)
        / (8 * np.pi * const.G * const.m_p)
    )
    return val.to(u.pc / u.cm ** 3).value


@lru_cache(maxsize=8)
def _cumulative_integral(om0: float, ode0: float, ok0: float, z_max: float = 6.0, n: int = 6001):
    """Tabulate I(z) = integral_0^z (1+z') / E(z') dz' for a (flat or curved)
    matter + Lambda universe, cached per cosmology. Radiation is negligible below
    z ~ 6 and is ignored."""
    z = np.linspace(0.0, z_max, n)
    e = np.sqrt(om0 * (1 + z) ** 3 + ok0 * (1 + z) ** 2 + ode0)
    integrand = (1 + z) / e
    cum = integrate.cumulative_trapezoid(integrand, z, initial=0.0)
    return z, cum


def dm_integral(z, cosmo=Planck18):
    """The dimensionless integral ``I(z) = int_0^z (1+z')/E(z') dz'``."""
    zz, cum = _cumulative_integral(float(cosmo.Om0), float(cosmo.Ode0), float(cosmo.Ok0))
    z = np.asarray(z, dtype=float)
    if np.any(z > zz[-1]):
        raise ValueError(f"redshift beyond tabulated range (max {zz[-1]})")
    return np.interp(z, zz, cum)


def mean_dm_cosmic(z, f_d: float = F_DIFFUSE_FIDUCIAL, cosmo=Planck18, y_he: float = Y_HE):
    """Mean DM contributed by diffuse ionized gas out to redshift z.

    ``<DM_cosmic>(z) = f_d * [3 c H0 Omega_b (1-Y/2) / 8 pi G m_p] * I(z)``.

    The factor (1+z) inside I(z) comes from two competing effects: the proper
    electron density grows as (1+z)^3 into the past, while each electron's
    contribution is weighted by 1/(1+z) for the redshifted frequency and by the
    proper path length dl = c dz / ((1+z) H(z)).

    ``f_d`` is the fraction of all baryons in diffuse ionized gas -- the
    intergalactic medium plus galaxy halos, but not stars or cold ISM.
    """
    return f_d * cosmic_dm_prefactor(cosmo, y_he) * dm_integral(z, cosmo)


def z_from_mean_dm(dm_cosmic, f_d: float = F_DIFFUSE_FIDUCIAL, cosmo=Planck18,
                   y_he: float = Y_HE):
    """Invert the mean relation: the redshift at which <DM_cosmic> = ``dm_cosmic``.

    Negative inputs return 0. This is the redshift an FRB would have if its
    excess DM were entirely intergalactic, i.e. an upper limit in the sense that
    any host-galaxy contribution would place it nearer.
    """
    zz, cum = _cumulative_integral(float(cosmo.Om0), float(cosmo.Ode0), float(cosmo.Ok0))
    dm_grid = f_d * cosmic_dm_prefactor(cosmo, y_he) * cum
    dm = np.clip(np.asarray(dm_cosmic, dtype=float), 0.0, None)
    return np.interp(dm, dm_grid, zz, right=np.nan)


# ---------------------------------------------------------------------------
# 3. The scatter: Macquart et al. (2020) p(DM_cosmic | z) and the host term
# ---------------------------------------------------------------------------
#
# Writing Delta = DM_cosmic / <DM_cosmic>(z), Macquart et al. propose
#
#     p(Delta) = A Delta^-beta exp[ -(Delta^-alpha - C0)^2 / (2 alpha^2 sigma^2) ],
#
# with sigma = F z^-1/2. F encodes how strongly galaxies' feedback has blown gas
# out of halos into a smoother intergalactic medium; the z^-1/2 is the usual
# square-root averaging over more structures along longer sightlines. C0 shifts
# the peak so that <Delta> = 1, and A normalises. At small sigma it tends to a
# Gaussian of width sigma; at large sigma it is skewed, with a long tail to high
# DM from sightlines that happen to pierce galaxy halos.
#
# The integrals are done in u = Delta^-alpha, where the exponential is a plain
# Gaussian in u and the power law becomes an integrable u^p singularity at 0:
#
#     int_0^inf Delta^(k-beta) e^{...} dDelta
#         = (1/alpha) int_0^inf u^((beta-k-1)/alpha - 1) e^{-(u-C0)^2/(2 s^2)} du,
#
# with s = alpha * sigma. For broad distributions (large sigma, i.e. low z or
# large F) C0 is large and negative and the Gaussian factor underflows at every
# u >= 0 -- e^-G with G = C0^2/(2 s^2) reaches hundreds. So whenever C0 < 0 we
# factor e^-G out analytically and carry it separately. What remains, the
# "reduced" exponent -(u^2 - 2 u C0)/(2 s^2), is never positive and never
# underflows where the integrand matters.


def _shift(c0, s):
    """G = C0^2 / (2 s^2) where C0 < 0, else 0: the factor carried analytically."""
    c0 = np.asarray(c0, dtype=float)
    return np.where(c0 < 0, c0 * c0 / (2 * s * s), 0.0)


def _log_kernel(u, c0, s):
    """-(u - C0)^2 / (2 s^2) + G, evaluated without cancellation."""
    c0 = np.asarray(c0, dtype=float)
    return np.where(c0 < 0, -(u * u - 2 * u * c0) / (2 * s * s), -((u - c0) ** 2) / (2 * s * s))


def _reduced_moment(k: int, c0: float, sigma: float, alpha: float, beta: float) -> float:
    """``e^G * int_0^inf Delta^(k - beta) exp(...) dDelta`` via the u substitution."""
    p = (beta - k - 1.0) / alpha - 1.0
    s = alpha * sigma

    def f(uu):
        return uu ** p * np.exp(_log_kernel(uu, c0, s))

    # Split the range where the integrand's mass sits so quad resolves both the
    # u^p singularity at 0 and the peak: at C0 for C0 > 0; for C0 < 0 the
    # integrand decays from u = 0 on the scale min(s, s^2/|C0|).
    if c0 > 0:
        pieces = [(0.0, c0), (c0, c0 + 12 * s), (c0 + 12 * s, np.inf)]
    else:
        scale = s if c0 == 0 else min(s, s * s / abs(c0))
        pieces = [(0.0, 30 * scale), (30 * scale, np.inf)]
    total = sum(integrate.quad(f, a, b, limit=200)[0] for a, b in pieces)
    return total / alpha


def macquart_c0(sigma: float, alpha: float = MACQUART_ALPHA, beta: float = MACQUART_BETA) -> float:
    """The C0 that makes the mean of p(Delta) exactly 1 for this sigma."""

    def mean_minus_one(c0):
        m0 = _reduced_moment(0, c0, sigma, alpha, beta)
        m1 = _reduced_moment(1, c0, sigma, alpha, beta)
        return m1 / m0 - 1.0

    # The mean falls monotonically as C0 grows (the peak sits near C0^(-1/alpha)).
    # For small sigma the root is C0 ~ 1; broad distributions need C0 far below
    # zero. Walk outward from 1 with a growing step until the root is bracketed.
    step = 0.5 * alpha * sigma
    hi = 1.0 + step
    while mean_minus_one(hi) > 0:
        hi += step
    lo = 1.0 - step
    while mean_minus_one(lo) < 0:
        step *= 1.5
        lo -= step
    return optimize.brentq(mean_minus_one, lo, hi, xtol=1e-12, rtol=1e-13, maxiter=300)


def macquart_log_norm(sigma: float, c0: float, alpha: float = MACQUART_ALPHA,
                      beta: float = MACQUART_BETA) -> float:
    """log of the reduced normalisation: p = exp(log_norm + log_kernel) Delta^-beta."""
    return -np.log(_reduced_moment(0, c0, sigma, alpha, beta))


@lru_cache(maxsize=4)
def _macquart_table(alpha: float, beta: float, s_min: float = 0.003, s_max: float = 40.0,
                    n: int = 140):
    """C0(sigma) and the log normalisation on a log grid, for fast interpolation."""
    sig = np.geomspace(s_min, s_max, n)
    c0 = np.array([macquart_c0(s, alpha, beta) for s in sig])
    log_a = np.array([macquart_log_norm(s, c, alpha, beta) for s, c in zip(sig, c0)])
    ls = np.log(sig)
    return (
        interpolate.CubicSpline(ls, c0),
        interpolate.CubicSpline(ls, log_a),
        (s_min, s_max),
    )


def macquart_params(sigma, alpha: float = MACQUART_ALPHA, beta: float = MACQUART_BETA):
    """(C0, log normalisation) for an array of sigma, from a cached exact table."""
    c0_spline, loga_spline, (s_min, s_max) = _macquart_table(alpha, beta)
    s = np.asarray(sigma, dtype=float)
    if np.any((s < s_min) | (s > s_max)):
        raise ValueError(f"sigma outside tabulated range [{s_min}, {s_max}]")
    ls = np.log(s)
    return c0_spline(ls), loga_spline(ls)


def macquart_log_pdf_delta(delta, sigma, alpha: float = MACQUART_ALPHA, beta: float = MACQUART_BETA):
    """log p(Delta | sigma); -inf for Delta <= 0. Broadcasts over both arguments."""
    delta = np.asarray(delta, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    c0, log_a = macquart_params(sigma, alpha, beta)
    d = np.where(delta > 0, delta, 1.0)
    u = d ** -alpha
    val = log_a + _log_kernel(u, c0, alpha * sigma) - beta * np.log(d)
    return np.where(delta > 0, val, -np.inf)


def macquart_pdf_delta(delta, sigma, alpha: float = MACQUART_ALPHA, beta: float = MACQUART_BETA):
    """p(Delta | sigma), normalised, with <Delta> = 1. Zero for Delta <= 0."""
    return np.exp(macquart_log_pdf_delta(delta, sigma, alpha, beta))


def lognormal_pdf(x, mu, sigma):
    """Lognormal density in x: median exp(mu), log-width sigma. Zero for x <= 0."""
    x = np.asarray(x, dtype=float)
    xs = np.where(x > 0, x, 1.0)
    val = np.exp(-((np.log(xs) - mu) ** 2) / (2 * sigma ** 2)) / (xs * sigma * np.sqrt(2 * np.pi))
    return np.where(x > 0, val, 0.0)


def pdf_dm_extragalactic(dm_eg, z, f_d, F, mu, sigma_host, *, cosmo=Planck18, n_conv: int = 400):
    """p(DM_EG | z) for one parameter set, by direct convolution.

    DM_EG = DM_cosmic + DM_host / (1+z): the host's rest-frame contribution is
    diluted by (1+z) because both the electron density and the path length are
    measured in the host's frame while the frequency is redshifted. Used for
    the figure's predictive band and to test the vectorised fit in analysis.
    """
    dm_eg = np.atleast_1d(np.asarray(dm_eg, dtype=float))
    mean = mean_dm_cosmic(z, f_d, cosmo)
    sigma = F / np.sqrt(z)
    out = np.zeros_like(dm_eg)
    for i, d in enumerate(dm_eg):
        if d <= 0:
            continue
        x = d * (np.arange(n_conv) + 0.5) / n_conv  # observed-frame host DM
        pc = macquart_pdf_delta((d - x) / mean, sigma) / mean
        px = (1 + z) * lognormal_pdf(x * (1 + z), mu, sigma_host)
        out[i] = np.sum(pc * px) * (d / n_conv)
    return out
