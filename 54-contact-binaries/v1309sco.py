"""V1309 Sco: a contact binary watched until it merged (OGLE photometry 2001-2010).

The progenitor was observed by OGLE about 1,340 times before the 2008 eruption
(Tylenda et al. 2011). We redo the analysis from the published photometry:

1. Periods per observing season and per ~50-point subsample, with our own
   multi-harmonic periodogram (Fourier series plus a polynomial trend, the
   least-squares relative of Schwarzenberg-Czerny's AoV method that Tylenda
   et al. used).
2. Eclipse timing (O - C) per subsample from the phase of the second Fourier
   harmonic, which follows the eclipses and is insensitive to the changing
   heights of the two maxima (Molnar et al. 2017 used the same idea).
3. Tylenda et al.'s exponential period law P = P0 exp(c / (t - t0)) fitted to
   our periods and, independently, to our O - C phases.
4. The changing light curve, season by season, and the onset of the 2008
   brightening.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import least_squares

import config as C


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

#: The 2008 eruption was discovered on 2008 Sep 2.5 (JD 2454712; Nakano 2008).
JD_DISCOVERY = 2454712.0


def load_ogle() -> dict:
    """OGLE-III (2001-2009) and OGLE-IV (2010) I-band photometry: HJD, I, error."""
    out = {}
    for key in ("OIII_I", "OIV_I", "OIII_V", "OIV_V"):
        d = np.loadtxt(C.OGLE_DIR / C.V1309SCO["ogle_files"][key])
        t = d[:, 0]
        if t[0] < 1e6:
            t = t + 2450000.0
        out[key] = {"t": t, "mag": d[:, 1], "err": d[:, 2]}
    i = {k: np.concatenate([out["OIII_I"][k], out["OIV_I"][k]]) for k in ("t", "mag", "err")}
    o = np.argsort(i["t"])
    out["I"] = {k: v[o] for k, v in i.items()}
    return out


def jd_to_year(jd):
    return 2000.0 + (np.asarray(jd) - C.JD_J2000) / C.YEAR_D


def seasons(t, gap: float = 60.0):
    """Split times into observing seasons separated by gaps longer than ``gap`` days."""
    br = np.nonzero(np.diff(t) > gap)[0] + 1
    return np.split(np.arange(t.size), br)


def subsamples(idx, size: int = 50):
    """Consecutive chunks of ~``size`` points (the last chunk absorbs the remainder)."""
    n = len(idx)
    k = max(1, int(round(n / size)))
    edges = np.linspace(0, n, k + 1).round().astype(int)
    return [idx[edges[j]:edges[j + 1]] for j in range(k)]


# ---------------------------------------------------------------------------
# Periodogram: Fourier series + polynomial trend, linear least squares per frequency
# ---------------------------------------------------------------------------


def _design(t, freq, n_harm, deg, t_mid, t_scale):
    cols = [((t - t_mid) / t_scale) ** j for j in range(deg + 1)]
    for k in range(1, n_harm + 1):
        arg = 2 * np.pi * k * freq * t
        cols += [np.cos(arg), np.sin(arg)]
    return np.column_stack(cols)


def periodogram(t, mag, err, freqs, n_harm: int = 4, deg: int = 1):
    """chi^2 of the best Fourier-plus-trend model at each trial frequency.

    Returns (chi2, chi2_trend_only). A power spectrum is 1 - chi2/chi2_trend.
    """
    w = 1.0 / err
    t_mid, t_scale = 0.5 * (t.min() + t.max()), max(1.0, 0.5 * (t.max() - t.min()))
    chi2 = np.empty(freqs.size)
    for n, f in enumerate(freqs):
        A = _design(t, f, n_harm, deg, t_mid, t_scale) * w[:, None]
        coef, res, rank, _ = np.linalg.lstsq(A, mag * w, rcond=None)
        chi2[n] = res[0] if res.size else np.sum((A @ coef - mag * w) ** 2)
    A0 = _design(t, 0.0, 0, deg, t_mid, t_scale) * w[:, None]
    c0, r0, _, _ = np.linalg.lstsq(A0, mag * w, rcond=None)
    chi2_0 = r0[0] if r0.size else np.sum((A0 @ c0 - mag * w) ** 2)
    return chi2, chi2_0


@dataclass
class PeriodResult:
    t_mean: float
    t_min: float
    t_max: float
    n: int
    period: float
    period_err: float
    power: float
    half_period_power: float


def best_period(t, mag, err, p_lo=1.38, p_hi=1.47, n_harm=4, deg=1, oversample=40) -> PeriodResult:
    """Best period in [p_lo, p_hi] with a curvature (Delta chi^2 = 1) uncertainty.

    The chi^2 is rescaled so that the best fit has reduced chi^2 = 1: OGLE's
    nominal errors do not include the star's cycle-to-cycle variations.
    """
    span = t.max() - t.min()
    df = 1.0 / (span * oversample)
    freqs = np.arange(1.0 / p_hi, 1.0 / p_lo, df)
    chi2, chi2_0 = periodogram(t, mag, err, freqs, n_harm, deg)
    j = int(np.argmin(chi2))
    # refine on a fine grid around the minimum
    fine = np.linspace(freqs[max(j - 2, 0)], freqs[min(j + 2, freqs.size - 1)], 201)
    c2f, _ = periodogram(t, mag, err, fine, n_harm, deg)
    jf = int(np.argmin(c2f))
    dof = t.size - (deg + 1 + 2 * n_harm) - 1
    scale = c2f[jf] / dof
    # parabola through the minimum for the curvature
    lo, hi = max(jf - 15, 0), min(jf + 16, fine.size)
    pp = np.polyfit(fine[lo:hi] - fine[jf], c2f[lo:hi] / scale, 2)
    f_best = fine[jf] - pp[1] / (2 * pp[0]) if pp[0] > 0 else fine[jf]
    sig_f = 1.0 / np.sqrt(pp[0]) if pp[0] > 0 else np.nan
    # power of the half-period alias (a one-maximum light curve at P/2)
    c2h, _ = periodogram(t, mag, err, np.array([2 * f_best]), 2, deg)
    return PeriodResult(t_mean=float(np.mean(t)), t_min=float(t.min()), t_max=float(t.max()), n=int(t.size),
                        period=float(1 / f_best), period_err=float(sig_f / f_best**2),
                        power=float(1 - c2f[jf] / chi2_0), half_period_power=float(1 - c2h[0] / chi2_0))


# ---------------------------------------------------------------------------
# The exponential period law and its integral (cycle count)
# ---------------------------------------------------------------------------


def exp_period(t, p0, c, t0):
    """Tylenda et al. (2011): P(t) = P0 exp(c / (t - t0)), valid for t < t0."""
    return p0 * np.exp(c / (t - t0))


def cycles(t, p0, c, t0, t_ref):
    """Number of orbits between t_ref and t for the exponential law (numerical integral)."""
    t = np.atleast_1d(t)
    grid = np.linspace(min(t.min(), t_ref), max(t.max(), t_ref), 4001)
    n = cumulative_trapezoid(1.0 / exp_period(grid, p0, c, t0), grid, initial=0.0)
    n -= np.interp(t_ref, grid, n)
    return np.interp(t, grid, n)


def fit_period_law(t, p, perr, p0_guess=1.4456, c_guess=15.29, t0_guess=2455233.5):
    """Weighted least-squares fit of the exponential law to period measurements."""
    def resid(x):
        return (exp_period(t, *x) - p) / perr
    r = least_squares(resid, [p0_guess, c_guess, t0_guess],
                      bounds=([1.40, 0.01, t.max() + 10], [1.50, 500.0, t.max() + 20000]))
    J = r.jac
    cov = np.linalg.inv(J.T @ J) * max(1.0, 2 * r.cost / max(1, t.size - 3))
    return r.x, np.sqrt(np.diag(cov)), cov, 2 * r.cost


# ---------------------------------------------------------------------------
# Eclipse timing from the second harmonic
# ---------------------------------------------------------------------------


def harmonic_phases(t, mag, err, p_ref, t_ref, n_harm=4, deg=1):
    """Fourier fit at a fixed period; returns amplitudes, phases and their errors.

    Model: m = trend + sum_k A_k cos(2 pi k (t - t_ref)/p_ref - phi_k).
    The k = 2 term carries the two eclipses per orbit. In magnitudes the
    star is faintest when the cosine is +1, so the deeper eclipse falls at
    orbital phase phi_2 / (4 pi) (mod 1/2).
    """
    w = 1.0 / err
    t_mid, t_scale = 0.5 * (t.min() + t.max()), max(1.0, 0.5 * (t.max() - t.min()))
    A = _design(t - t_ref, 1.0 / p_ref, n_harm, deg, t_mid - t_ref, t_scale)
    coef, _, _, _ = np.linalg.lstsq(A * w[:, None], mag * w, rcond=None)
    res = mag - A @ coef
    dof = max(1, t.size - A.shape[1])
    s2 = np.sum((res * w) ** 2) / dof
    cov = np.linalg.inv((A * w[:, None]).T @ (A * w[:, None])) * s2
    out = {}
    for k in range(1, n_harm + 1):
        ia, ib = deg + 1 + 2 * (k - 1), deg + 2 + 2 * (k - 1)
        a, b = coef[ia], coef[ib]
        amp = np.hypot(a, b)
        phi = np.arctan2(b, a)
        # error propagation for phi = atan2(b, a)
        ga = -b / amp**2
        gb = a / amp**2
        sphi = np.sqrt(ga * ga * cov[ia, ia] + gb * gb * cov[ib, ib] + 2 * ga * gb * cov[ia, ib])
        out[k] = (amp, phi, sphi)
    out["rms"] = float(np.sqrt(np.mean(res**2)))
    return out


def oc_from_phase(phi2, sphi2, p_ref):
    """O - C (days) of the deeper eclipse from the second-harmonic phase, in (-P/4, P/4]."""
    oc = p_ref * phi2 / (4 * np.pi)
    return oc, p_ref * sphi2 / (4 * np.pi)


# ---------------------------------------------------------------------------
# Light-curve shape
# ---------------------------------------------------------------------------


def fold(t, p0, c, t0, t_ref):
    """Orbital phase with the exponential ephemeris (phase 0 at t_ref)."""
    return np.mod(cycles(t, p0, c, t0, t_ref), 1.0)


def binned_curve(phase, mag, nbins=25):
    edges = np.linspace(0, 1, nbins + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    med = np.full(nbins, np.nan)
    err = np.full(nbins, np.nan)
    for j in range(nbins):
        s = (phase >= edges[j]) & (phase < edges[j + 1])
        if s.sum() >= 3:
            med[j] = np.median(mag[s])
            err[j] = 1.2533 * np.std(mag[s]) / np.sqrt(s.sum())
    return mid, med, err


def maxima_difference(t, mag, err, p_ref, t_ref):
    """Height difference of the two maxima (mag) from a 4-harmonic fit: positive
    when the maximum at phase 0.25 is brighter than the one at 0.75."""
    w = 1.0 / err
    t_mid, t_scale = 0.5 * (t.min() + t.max()), max(1.0, 0.5 * (t.max() - t.min()))
    A = _design(t - t_ref, 1.0 / p_ref, 4, 1, t_mid - t_ref, t_scale)
    coef, _, _, _ = np.linalg.lstsq(A * w[:, None], mag * w, rcond=None)
    ph = np.linspace(0, 1, 1001)
    tt = t_ref + ph * p_ref
    model = _design(tt - t_ref, 1.0 / p_ref, 4, 1, t_mid - t_ref, t_scale)[:, 2:] @ coef[2:]
    m1 = model[(ph > 0.0) & (ph < 0.5)].min()  # brightest point in the first half
    m2 = model[(ph > 0.5) & (ph < 1.0)].min()
    amp = model.max() - model.min()
    return float(m2 - m1), float(amp), ph, model


# ---------------------------------------------------------------------------
# The 2008 eruption: onset of the exponential brightening
# ---------------------------------------------------------------------------


def eruption_onset(t, mag, t_lo=2454400.0, t_hi=JD_DISCOVERY):
    """Fit flux = F_q + F_1 exp((t - t_on)/tau) to the 2008 pre-discovery data.

    Returns t_on (when the rising component equals the quiescent flux),
    the e-folding time tau, and the model. Flux is in units of I = 0.
    """
    s = (t > t_lo) & (t < t_hi)
    tt, ff = t[s], 10 ** (-0.4 * mag[s])

    def resid(x):
        fq, t_on, tau = x
        model = fq * (1.0 + np.exp((tt - t_on) / tau))
        return (model - ff) / ff  # relative residuals (the rise spans 2.5 decades)

    fq0 = np.median(ff[tt < tt.min() + 60])
    r = least_squares(resid, [fq0, 2454550.0, 30.0], bounds=([0, 2454300, 1], [np.inf, 2454712, 500]))
    J = r.jac
    cov = np.linalg.inv(J.T @ J) * (2 * r.cost / max(1, tt.size - 3))
    return {"f_quiescent": r.x[0], "t_on": r.x[1], "tau": r.x[2], "err": np.sqrt(np.diag(cov)), "n": int(tt.size)}


# ---------------------------------------------------------------------------
# Bootstrap uncertainties and the O - C against the exponential ephemeris
# ---------------------------------------------------------------------------


def bootstrap_period(t, mag, err, n_boot=200, seed=C.SEED, **kw):
    """Period uncertainty from resampling the points of a subsample with
    replacement (captures cycle-to-cycle light-curve changes that the
    nominal OGLE errors do not)."""
    rng = np.random.default_rng(seed)
    best = best_period(t, mag, err, **kw)
    p = best.period
    span = t.max() - t.min()
    # search only around the best period: +-4 resolution elements
    dp = 4 * p * p / span
    vals = []
    for _ in range(n_boot):
        k = rng.integers(0, t.size, t.size)
        if np.unique(t[k]).size < 12:
            continue
        r = best_period(t[k], mag[k], err[k], p_lo=p - dp, p_hi=p + dp, oversample=20, **kw)
        vals.append(r.period)
    vals = np.array(vals)
    lo, hi = np.percentile(vals, [15.87, 84.13])
    return best, 0.5 * (hi - lo), vals


def oc_exponential(t, mag, err, groups, law, t_ref, n_harm=4, deg=1):
    """Eclipse-time offsets (days) of each group of points relative to the
    exponential ephemeris ``law`` = (P0, c, t0), with phase zero at t_ref.

    Each group is fitted with a Fourier series in the ephemeris phase; the
    phase of the k = 2 term gives the offset of the deeper eclipse.
    """
    p0, c, t0 = law
    out = []
    for g in groups:
        ph = cycles(t[g], p0, c, t0, t_ref)
        # fit in phase units: period 1, reference 0
        h = harmonic_phases(ph, mag[g], err[g], 1.0, 0.0, n_harm=n_harm, deg=deg)
        dphi = h[2][1] / (4 * np.pi)          # cycles, in (-1/4, 1/4]
        sdphi = h[2][2] / (4 * np.pi)
        pmid = exp_period(np.mean(t[g]), p0, c, t0)
        out.append({"t": float(np.mean(t[g])), "dphase": float(dphi), "dphase_err": float(sdphi),
                    "oc_days": float(dphi * pmid), "oc_err_days": float(sdphi * pmid),
                    "a1": float(h[1][0]), "a2": float(h[2][0]), "rms": h["rms"], "n": int(len(g))})
    return out


def fit_law_to_phases(t, dphase, dphase_err, law_ref, t_ref):
    """Refit (P0, c, t0) and a phase zero point to eclipse phases measured
    against a reference law: minimises the cycle-count residuals."""
    n_ref = cycles(t, *law_ref, t_ref)
    n_obs = n_ref + dphase

    def resid(x):
        p0, c, t0, z = x
        return (cycles(t, p0, c, t0, t_ref) + z - n_obs) / dphase_err

    x0 = [*law_ref, 0.0]
    r = least_squares(resid, x0, x_scale=[1e-3, 1.0, 50.0, 0.01],
                      bounds=([1.40, 0.01, t.max() + 10, -0.5], [1.50, 500.0, t.max() + 20000, 0.5]))
    J = r.jac
    chi2 = 2 * r.cost
    dof = max(1, t.size - 4)
    cov = np.linalg.inv(J.T @ J) * max(1.0, chi2 / dof)
    return r.x, np.sqrt(np.diag(cov)), chi2, dof
