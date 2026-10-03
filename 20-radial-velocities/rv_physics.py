"""Radial-velocity physics: from a star's wobble to a planet's mass.

A star and its planet orbit their common centre of mass. Seen from Earth, the
star's velocity along the line of sight swings back and forth. For a single
planet on a Keplerian orbit that velocity is

    v(t) = gamma + K [cos(nu(t) + omega) + e cos(omega)]

where ``nu`` is the planet's true anomaly (its angle from periastron), ``omega``
the argument of periastron *of the star's orbit*, ``e`` the eccentricity, and
``K`` the semi-amplitude. ``nu(t)`` comes from Kepler's equation, which has no
closed-form solution, so it is solved numerically here.

K, the period P and the star's mass M give the planet's mass through the
*mass function*:

    (m sin i)^3 / (M + m)^2 = P K^3 (1 - e^2)^(3/2) / (2 pi G)

Only m sin i appears, because the velocity is measured along our line of sight
and the orbit's tilt ``i`` is unknown: every mass derived here is a minimum
mass. Constants are the IAU 2015 nominal values (exact by definition); masses
are carried as GM products, which are known far better than G or M alone.

Periodograms: ``gls`` is the generalised Lomb-Scargle periodogram
(Zechmeister & Kurster 2009) -- at each trial frequency a weighted least-
squares fit of a sinusoid -- extended so that data from several instruments
each keep their own zero point, fitted afresh at every frequency. With one
instrument it reduces exactly to the usual GLS (tested against astropy).
"""

from __future__ import annotations

import numpy as np
from scipy.special import gammaln

# --- Constants (IAU 2015 Resolution B3 nominal values; m^3 s^-2) -------------
GM_SUN = 1.3271244e20
GM_JUP = 1.2668653e17
GM_EARTH = 3.986004e14
AU_M = 1.495978707e11
DAY_S = 86400.0
#: A star returns to the same place in the sky every sidereal day, so nightly
#: observations are spaced by multiples of it, not of the 24-h solar day.
SIDEREAL_DAY_D = 0.99726957
#: ...and to the same place relative to the Sun every (sidereal) year.
YEAR_D = 365.25636
#: Mean density of Earth and Jupiter, g/cm^3, for density comparisons.
RHO_EARTH = 5.513


# ---------------------------------------------------------------------------
# Kepler's equation and the velocity curve
# ---------------------------------------------------------------------------

def solve_kepler(M, e, *, tol: float = 1e-12, max_iter: int = 30):
    """Eccentric anomaly E from mean anomaly M: solves ``M = E - e sin E``.

    Newton's method from Danby's starting point ``E0 = M + 0.85 e sign(sin M)``,
    which converges for every ``0 <= e < 1``. Fully vectorised: ``M`` and ``e``
    broadcast against each other.
    """
    M = np.asarray(M, dtype=float)
    e = np.asarray(e, dtype=float)
    M = np.mod(M, 2 * np.pi)
    E = M + 0.85 * e * np.sign(np.sin(M))
    for _ in range(max_iter):
        f = E - e * np.sin(E) - M
        fp = 1.0 - e * np.cos(E)
        dE = f / fp
        E = E - dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def true_anomaly(E, e):
    """True anomaly from eccentric anomaly (quadrant-safe half-angle form)."""
    return 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(E / 2.0),
                            np.sqrt(1.0 - e) * np.cos(E / 2.0))


def tc_to_tp(tc, P, e, omega):
    """Time of periastron from the time of inferior conjunction.

    Conjunction -- the moment the planet passes in front of the star, when a
    transit would happen -- is at true anomaly ``pi/2 - omega`` (with omega of
    the star's orbit). Parameterising orbits by ``tc`` instead of ``tp`` keeps
    the phase well defined as e -> 0, where periastron stops meaning anything.
    """
    f = np.pi / 2.0 - omega
    E = 2.0 * np.arctan(np.sqrt((1.0 - e) / (1.0 + e)) * np.tan(f / 2.0))
    M = E - e * np.sin(E)
    return tc - P / (2.0 * np.pi) * M


def rv_keplerian(t, P, K, e, omega, tc):
    """Stellar radial velocity (same units as K) of one Keplerian orbit.

    ``v = K [cos(nu + omega) + e cos(omega)]``, zero-mean over an orbit.
    Positive = moving away from us (redshift). All arguments broadcast.
    """
    tp = tc_to_tp(tc, P, e, omega)
    M = 2.0 * np.pi * (t - tp) / P
    E = solve_kepler(M, e)
    nu = true_anomaly(E, e)
    return K * (np.cos(nu + omega) + e * np.cos(omega))


def he_to_eomega(h, k):
    """(sqrt(e) cos w, sqrt(e) sin w) -> (e, w). Sampling in these variables
    gives uniform priors on e and w and avoids the e = 0 boundary pile-up."""
    e = h * h + k * k
    return e, np.arctan2(k, h)


# ---------------------------------------------------------------------------
# Masses and orbits
# ---------------------------------------------------------------------------

def mass_function_gm(P_days, K, e):
    """G times the mass function, ``P K^3 (1-e^2)^1.5 / (2 pi)``, in m^3 s^-2."""
    return P_days * DAY_S * np.asarray(K, float) ** 3 * (1.0 - np.asarray(e, float) ** 2) ** 1.5 / (2.0 * np.pi)


def msini_gm(P_days, K, e, gm_star, *, n_iter: int = 60):
    """G * (m sin i) solving the mass function exactly.

    ``(m sin i)^3 / (M + m)^2 = f`` has no closed form in m; iterate
    ``m <- [f (M + m)^2]^(1/3)`` from ``m = (f M^2)^(1/3)``. Converges
    geometrically because m << M. As is standard, sin i = 1 is assumed inside
    the (M + m) term -- a correction of order m/M, below 0.1% for every planet
    here.
    """
    f = mass_function_gm(P_days, K, e)
    gm_star = np.asarray(gm_star, float)
    m = np.cbrt(f * gm_star ** 2)
    for _ in range(n_iter):
        m = np.cbrt(f * (gm_star + m) ** 2)
    return m


def msini_mjup(P_days, K, e, mstar_msun):
    return msini_gm(P_days, K, e, np.asarray(mstar_msun) * GM_SUN) / GM_JUP


def msini_mearth(P_days, K, e, mstar_msun):
    return msini_gm(P_days, K, e, np.asarray(mstar_msun) * GM_SUN) / GM_EARTH


def semi_amplitude(P_days, m_gm, gm_star, e=0.0, sin_i=1.0):
    """The star's velocity semi-amplitude K (m/s) caused by a planet of mass
    ``m_gm`` (G times mass). The inverse of ``msini_gm``: used to check that
    Jupiter moves the Sun by 12.5 m/s and Earth by 0.09 m/s."""
    P = np.asarray(P_days, float) * DAY_S
    return ((2 * np.pi / P) ** (1 / 3) * m_gm * sin_i / (gm_star + m_gm) ** (2 / 3)
            / np.sqrt(1.0 - np.asarray(e, float) ** 2))


def semi_major_axis_au(P_days, gm_star, m_gm=0.0):
    """Kepler's third law: a^3 = G(M + m) P^2 / (4 pi^2)."""
    P = np.asarray(P_days, float) * DAY_S
    return np.cbrt((gm_star + m_gm) * P ** 2 / (4 * np.pi ** 2)) / AU_M


def orbital_speed(P_days, gm_star, m_gm=0.0):
    """Relative orbital speed of a circular orbit, ``2 pi a / P`` (m/s)."""
    P = np.asarray(P_days, float) * DAY_S
    return np.cbrt(2 * np.pi * (gm_star + m_gm) / P)


# ---------------------------------------------------------------------------
# Periodograms
# ---------------------------------------------------------------------------

def frequency_grid(t, period_min: float, period_max: float | None = None,
                   oversample: int = 10):
    """Evenly spaced frequencies (cycles/day) from 1/period_max to 1/period_min.

    Peaks have width ~1/T for a baseline T, so a step of 1/(oversample T)
    samples each peak about ``oversample`` times. ``period_max`` defaults to
    twice the baseline.
    """
    T = float(np.ptp(t))
    fmin = 1.0 / (period_max if period_max else 2.0 * T)
    fmax = 1.0 / period_min
    df = 1.0 / (oversample * T)
    return np.arange(fmin, fmax + df, df)


def _instrument_design(inst, n_inst):
    J = np.zeros((inst.size, n_inst))
    J[np.arange(inst.size), inst] = 1.0
    return J


def gls(t, y, err, freqs, inst=None, *, return_fit: bool = False, chunk_elems: int = 3_000_000):
    """Generalised Lomb-Scargle periodogram with one zero point per instrument.

    At each frequency ``f`` the model ``a cos(2 pi f t) + b sin(2 pi f t) +
    gamma_k`` (one ``gamma`` per instrument ``k``) is fitted by weighted least
    squares, weights ``1/err^2``. The power is the fraction of the weighted
    variance about the best constant-offsets model that the sinusoid removes::

        p(f) = (chi2_offsets - chi2(f)) / chi2_offsets,    0 <= p <= 1

    The offsets are eliminated analytically (for fixed a, b each gamma_k is a
    weighted mean), leaving a 2x2 system per frequency. With a single
    instrument this is exactly Zechmeister & Kurster's (2009) GLS.

    Returns the power, and with ``return_fit`` also the amplitude and the
    phase ``phi`` of the fitted sinusoid ``A cos(2 pi f t - phi)``.
    """
    t = np.asarray(t, float)
    y = np.asarray(y, float)
    w = 1.0 / np.asarray(err, float) ** 2
    if inst is None:
        inst = np.zeros(t.size, int)
    inst = np.asarray(inst, int)
    n_inst = int(inst.max()) + 1
    J = _instrument_design(inst, n_inst)
    W_k = w @ J
    if np.any(W_k <= 0):
        raise ValueError("an instrument index has no data")
    Y_k = (w * y) @ J
    chi2_0 = np.sum(w * y * y) - np.sum(Y_k ** 2 / W_k)
    wy = w * y
    WJ = w[:, None] * J
    t0 = t - t.mean()  # a centred time axis keeps the phases numerically clean

    freqs = np.asarray(freqs, float)
    power = np.empty(freqs.size)
    amp = np.empty(freqs.size) if return_fit else None
    phase = np.empty(freqs.size) if return_fit else None
    step = max(1, chunk_elems // t.size)
    for i0 in range(0, freqs.size, step):
        f = freqs[i0:i0 + step]
        z = np.exp(2j * np.pi * np.outer(f, t0))
        c, s = z.real, z.imag
        C_k = c @ WJ
        S_k = s @ WJ
        CC = (c * c) @ w - np.sum(C_k ** 2 / W_k, axis=1)
        SS = (s * s) @ w - np.sum(S_k ** 2 / W_k, axis=1)
        CS = (c * s) @ w - np.sum(C_k * S_k / W_k, axis=1)
        YC = c @ wy - np.sum(C_k * Y_k / W_k, axis=1)
        YS = s @ wy - np.sum(S_k * Y_k / W_k, axis=1)
        D = CC * SS - CS ** 2
        a = (SS * YC - CS * YS) / D
        b = (CC * YS - CS * YC) / D
        power[i0:i0 + step] = (a * YC + b * YS) / chi2_0
        if return_fit:
            amp[i0:i0 + step] = np.hypot(a, b)
            # a cos(x) + b sin(x) with x = 2 pi f (t - tbar)  ->  A cos(2 pi f t - phi)
            phase[i0:i0 + step] = np.arctan2(b, a) + 2 * np.pi * f * t.mean()
    if return_fit:
        return power, amp, np.mod(phase, 2 * np.pi)
    return power


def spectral_window(t, freqs, *, chunk_elems: int = 3_000_000):
    """Power of the sampling pattern alone, ``|sum_k exp(2 pi i f t_k)|^2 / N^2``.

    A perfect sinusoid sampled at times t produces a periodogram peak at its
    frequency *plus a copy of this window centred on it* -- that is what an
    alias is. Peaks here at 1 cycle per sidereal day and 1 per year mean every
    real signal has ghosts at f +- 1.0027/day and f +- 1/year.
    """
    t = np.asarray(t, float) - np.mean(t)
    freqs = np.asarray(freqs, float)
    out = np.empty(freqs.size)
    step = max(1, chunk_elems // t.size)
    for i0 in range(0, freqs.size, step):
        z = np.exp(2j * np.pi * np.outer(freqs[i0:i0 + step], t)).sum(axis=1)
        out[i0:i0 + step] = np.abs(z) ** 2 / t.size ** 2
    return out


def top_peaks(freqs, power, n: int = 5, min_separation: float | None = None):
    """The ``n`` highest local maxima, at least ``min_separation`` apart in
    frequency (default: 1.5x the grid step times 10). Returns indices."""
    freqs = np.asarray(freqs)
    power = np.asarray(power)
    is_max = np.r_[False, (power[1:-1] >= power[:-2]) & (power[1:-1] >= power[2:]), False]
    cand = np.flatnonzero(is_max)
    cand = cand[np.argsort(power[cand])[::-1]]
    if min_separation is None:
        min_separation = 15 * np.median(np.diff(freqs))
    chosen: list[int] = []
    for i in cand:
        if all(abs(freqs[i] - freqs[j]) > min_separation for j in chosen):
            chosen.append(int(i))
        if len(chosen) == n:
            break
    return np.array(chosen, int)


def refine_peak(t, y, err, f0, df, inst=None, *, n: int = 201):
    """Locate a peak's maximum to ~1% of the grid step with a fine local grid."""
    fine = np.linspace(f0 - df, f0 + df, n)
    p = gls(t, y, err, fine, inst)
    i = int(np.argmax(p))
    return fine[i], p[i]


def _weighted_time_variance(t, err):
    w = 1.0 / np.asarray(err, float) ** 2
    tm = np.average(t, weights=w)
    return np.average((t - tm) ** 2, weights=w)


def fap_baluev(power, t, err, fmax: float, n_offsets: int = 1, n_eff: int | None = None):
    """Analytic false-alarm probability of the highest peak (Baluev 2008).

    The chance that pure white noise, sampled like the data, would produce a
    peak at least this high *anywhere* up to frequency ``fmax``. Baluev's
    bound, ``1 - (1 - FAP_single) exp(-tau)``, accounts for the number of
    effectively independent frequencies searched via the effective baseline
    ``sqrt(4 pi var(t))``. For the "standard" normalisation (power in [0, 1])
    with ``n_offsets`` fitted zero points: the null model has dH = n_offsets
    parameters and the periodic model dK = n_offsets + 2. With one offset this
    is the formula astropy implements (tested against it).

    ``n_eff`` replaces the number of points N when points are not independent:
    forty exposures in one night share that night's noise and count once.
    Using N there makes the FAP absurdly small for night-clumped data.
    """
    z = np.asarray(power, float)
    N = len(t) if n_eff is None else int(n_eff)
    dH, dK = n_offsets, n_offsets + 2
    NH, NK = N - dH, N - dK
    fap_single = (1.0 - z) ** (0.5 * NK)
    T_eff = np.sqrt(4 * np.pi * _weighted_time_variance(t, err))
    W = fmax * T_eff
    log_gamma = 0.5 * np.log(2.0 / NH) + gammaln(NH / 2.0) - gammaln((NH - 1) / 2.0)
    tau = np.exp(log_gamma) * W * (1.0 - z) ** (0.5 * (NK - 1)) * np.sqrt(0.5 * NH * z)
    return np.clip(-np.expm1(-tau) + fap_single * np.exp(-tau), 0.0, 1.0)


def alias_frequencies(f_signal: float, f_sampling: float, orders=(1,)):
    """Frequencies where a signal at ``f_signal`` reappears through a sampling
    frequency ``f_sampling``: ``|f_signal +- k f_sampling|``."""
    out = []
    for k in orders:
        out += [abs(f_signal - k * f_sampling), f_signal + k * f_sampling]
    return np.array(out)
