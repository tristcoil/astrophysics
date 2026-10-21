"""The published spectra: FIRAS's blackbody and Planck's angular power spectrum.

Everything in this module runs on ~0.4 MB of downloads. It answers three
questions:

1. **Is the CMB a blackbody, and how hot is it?** Fit a Planck curve to the
   FIRAS monopole spectrum (``fit_firas``).
2. **Where are the acoustic peaks, and how high?** Measure them in Planck's
   power spectrum ourselves, with the same method applied to Planck's best-fit
   model as a control (``measure_peaks``).
3. **What angle do they imply?** The peaks form a harmonic series whose spacing
   is set by the angular size of the sound horizon (``acoustic_scale``).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import curve_fit, least_squares

import cmb_physics as cp
from config import (ENVELOPE_ORDERS, FIRST_ENVELOPE_PEAK, PEAK_GUESSES, PEAK_HALF_WINDOWS,
                    PEAK_MC_DRAWS, RANDOM_SEED, product)

# Peaks 1-4 are fitted in D_l directly; from FIRST_ENVELOPE_PEAK on, the damping
# tail is divided out first (see ``damping_envelope``). Planck's own peak table
# was made the same way: its peak 1-4 positions match direct fits while its 5-8
# positions match envelope-divided ones.

# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------


@dataclass
class Spectrum:
    """A TT power spectrum: D_l in uK^2 with (possibly asymmetric) 1-sigma errors."""

    ell: np.ndarray
    dl: np.ndarray
    err_lo: np.ndarray
    err_hi: np.ndarray

    @property
    def err(self) -> np.ndarray:
        """Symmetrised error. Exact for l >= 30, where Planck's errors are symmetric."""
        return 0.5 * (self.err_lo + self.err_hi)

    def select(self, lmin: float, lmax: float = np.inf) -> "Spectrum":
        m = (self.ell >= lmin) & (self.ell <= lmax)
        return Spectrum(self.ell[m], self.dl[m], self.err_lo[m], self.err_hi[m])


@dataclass
class Theory:
    """Planck's best-fit LCDM model: lensed D_l in uK^2 for l = 2..2508."""

    ell: np.ndarray
    tt: np.ndarray
    te: np.ndarray
    ee: np.ndarray


def load_tt_full() -> Spectrum:
    """Planck 2018 TT spectrum at every multipole from 2 to 2508.

    Two different pipelines stitched at l = 30. Below it, Commander (a
    pixel-based Bayesian method) -- its errors are asymmetric because at low l
    the distribution of a power estimate is a skewed chi-squared, not a
    Gaussian. Above it, the co-added high-l spectrum with Planck's foreground
    model subtracted.
    """
    a = np.loadtxt(product("tt_full").path)
    return Spectrum(a[:, 0].astype(int), a[:, 1], a[:, 2], a[:, 3])


def load_tt_binned() -> tuple[Spectrum, np.ndarray]:
    """Planck's binned spectrum (bins of 30 from l = 30) and its binned best fit."""
    a = np.loadtxt(product("tt_binned").path)
    return Spectrum(a[:, 0], a[:, 1], a[:, 2], a[:, 3]), a[:, 4]


def load_theory() -> Theory:
    a = np.loadtxt(product("theory").path)
    return Theory(a[:, 0].astype(int), a[:, 1], a[:, 2], a[:, 3])


def load_bestfit_parameters() -> dict[str, float]:
    """Parameters of Planck's best-fit model, keyed by their CosmoMC names.

    The file lists each parameter as ``index  value  name  latex``; derived
    quantities (H0, sound horizon, age, ...) follow the fitted ones.
    """
    params = {}
    for line in product("bestfit").path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                params[parts[2]] = float(parts[1])
            except ValueError:
                pass
    return params


@dataclass
class Firas:
    wavenumber: np.ndarray  # cm^-1
    ghz: np.ndarray
    intensity: np.ndarray  # MJy/sr
    residual: np.ndarray  # kJy/sr, relative to a 2.725 K blackbody
    sigma: np.ndarray  # kJy/sr
    galaxy: np.ndarray  # kJy/sr, the Galactic emission model removed


def load_firas() -> Firas:
    a = np.loadtxt(product("firas").path, comments="#")
    return Firas(a[:, 0], cp.wavenumber_to_ghz(a[:, 0]), a[:, 1], a[:, 2], a[:, 3], a[:, 4])


# --------------------------------------------------------------------------
# FIRAS: the blackbody
# --------------------------------------------------------------------------


def fit_firas(firas: Firas) -> dict:
    """Fit a blackbody, then a blackbody plus a Compton-y distortion.

    The intensity column is computed by LAMBDA as a 2.725 K blackbody plus
    FIRAS's measured residuals, rounded to 1 kJy/sr; we rebuild it from the
    unrounded residuals instead, which changes nothing but the last digit.

    **What this fit can and cannot tell us.** FIRAS measured the *difference*
    between the sky and an internal blackbody calibrator with exquisite
    precision, so the spectral *shape* is known to ~50 parts per million of
    the peak. The *absolute* temperature rests on the calibrator's
    thermometry, good to about 1 mK. The statistical error our fit returns
    (tens of microkelvin) is therefore real for the shape and meaningless for
    the absolute scale; the honest FIRAS-only number is 2.725 +/- 0.001 K
    (Fixsen & Mather 2002).
    """
    t_ref = 2.725
    sky = cp.planck_intensity(firas.ghz, t_ref) + firas.residual / 1e3
    sigma = firas.sigma / 1e3

    def resid_bb(p):
        return (cp.planck_intensity(firas.ghz, p[0]) - sky) / sigma

    def resid_y(p):
        t, y = p
        return (cp.planck_intensity(firas.ghz, t) + y * cp.y_distortion(firas.ghz, t) - sky) / sigma

    out = {}
    for name, fn, p0 in (("blackbody", resid_bb, [2.7]), ("with_y", resid_y, [2.7, 0.0])):
        sol = least_squares(fn, p0, x_scale="jac", xtol=1e-14, ftol=1e-14)
        cov = np.linalg.inv(sol.jac.T @ sol.jac)
        chi2 = float(np.sum(sol.fun**2))
        out[name] = {
            "params": sol.x.tolist(),
            "errors": np.sqrt(np.diag(cov)).tolist(),
            "chi2": chi2,
            "dof": len(sky) - len(p0),
        }

    t_best = out["blackbody"]["params"][0]
    model = cp.planck_intensity(firas.ghz, t_best)
    peak = float(cp.planck_intensity(cp.wien_peak_ghz(t_best), t_best))
    resid_kjy = (sky - model) * 1e3
    out.update({
        "temperature_k": t_best,
        "temperature_stat_err_k": out["blackbody"]["errors"][0],
        "y": out["with_y"]["params"][1],
        "y_err": out["with_y"]["errors"][1],
        "peak_ghz": cp.wien_peak_ghz(t_best),
        "peak_mjy_sr": peak,
        "rms_residual_kjy_sr": float(np.sqrt(np.mean(resid_kjy**2))),
        "max_abs_residual_fraction_of_peak": float(np.max(np.abs(resid_kjy)) / 1e3 / peak),
        "n_points": int(len(sky)),
        "sky_mjy_sr": sky.tolist(),
        "residual_kjy_sr": resid_kjy.tolist(),
    })
    return out


# --------------------------------------------------------------------------
# Acoustic peaks
# --------------------------------------------------------------------------


def _gauss_lin(x, amp, mu, s, b0, b1):
    """A Gaussian bump on a sloping straight line."""
    return amp * np.exp(-0.5 * ((x - mu) / s) ** 2) + b0 + b1 * (x - mu)


def _fit_extremum(x, y, w, guess: float, half_window: float):
    """Position and value of the maximum of a Gaussian-plus-line fit in a window.

    We report the maximum of the fitted curve rather than the Gaussian's centre:
    with a sloping baseline the two differ, and the maximum is the thing a
    reader of the plot means by "where the peak is".
    """
    m = (x >= guess - half_window) & (x <= guess + half_window)
    p0 = [y[m].max() - y[m].min(), guess, half_window / 1.5, y[m].min(), 0.0]
    p, _ = curve_fit(_gauss_lin, x[m], y[m], p0=p0, sigma=w[m], maxfev=20000)
    xx = np.linspace(guess - half_window, guess + half_window, 4001)
    yy = _gauss_lin(xx, *p)
    i = int(np.argmax(yy))
    return float(xx[i]), float(yy[i]), p


def damping_envelope(ell, dl, err, *, order: int = 5, lmin: float = 300.0):
    """A smooth curve through the damping tail, to divide out before finding peaks.

    Photons diffuse a finite distance before recombination, erasing structure
    on small scales -- "Silk damping". It makes the spectrum fall steeply
    beyond l ~ 1000, steeply enough that the higher oscillation maxima stop
    being local maxima of D_l at all: in Planck's best-fit model there is no
    local maximum anywhere near the 7th or 8th peak, just shoulders on a
    falling curve. Dividing by a featureless fit to the fall exposes the
    oscillation itself. This mirrors the Planck team's procedure ("we have to
    remove a featureless damping tail to fit the higher C_l^TT region",
    Planck 2015 results I, Sect. 8.2.3).

    The envelope is a polynomial of degree ``order`` in l fitted to log D_l --
    too stiff to follow oscillations with a period of ~300 in l over the
    fitted range, flexible enough to follow the decline.
    """
    ell = np.asarray(ell, float)
    m = ell >= lmin
    t = (ell[m] - 1400.0) / 1100.0
    design = np.vander(t, order + 1)
    wts = np.asarray(dl, float)[m] / np.asarray(err, float)[m]  # 1/sigma of log D
    coef = np.linalg.lstsq(design * wts[:, None],
                           np.log(np.clip(dl[m], 1.0, None)) * wts, rcond=None)[0]
    return lambda x: np.exp(np.polyval(coef, (np.asarray(x, float) - 1400.0) / 1100.0))


def _measure_one(ell, dl, err, n: int, *, envelope_order: int = 5):
    g, hw = PEAK_GUESSES[n - 1], PEAK_HALF_WINDOWS[n - 1]
    if n < FIRST_ENVELOPE_PEAK:
        pos, height, _ = _fit_extremum(ell, dl, err, g, hw)
        return pos, height
    env = damping_envelope(ell, dl, err, order=envelope_order)
    e = env(ell)
    pos, ratio, _ = _fit_extremum(ell, dl / e, err / e, g, hw)
    return pos, float(ratio * env(pos))


def measure_peaks(spec: Spectrum, *, n_peaks: int = 8, draws: int = PEAK_MC_DRAWS,
                  seed: int = RANDOM_SEED, errors: bool = True) -> list[dict]:
    """Measure the positions and heights of the TT acoustic peaks.

    Errors come from re-running the whole measurement on ``draws`` copies of
    the spectrum, each with Gaussian noise of Planck's quoted size added at
    every multipole. This treats neighbouring multipoles as independent, which
    they are not quite (masking correlates them at the few-percent level), so
    the errors are somewhat optimistic.

    For peaks measured through the damping envelope a systematic error is added:
    half the spread of positions obtained with envelopes of polynomial order 4,
    5 and 6. Where it matters (peaks 7 and 8) it is comparable to the
    statistical error.
    """
    spec = spec.select(30)
    ell, dl, err = spec.ell.astype(float), spec.dl, spec.err
    rng = np.random.default_rng(seed)
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for n in range(1, n_peaks + 1):
            pos, height = _measure_one(ell, dl, err, n)
            rec = {"n": n, "ell": pos, "dl": height,
                   "method": "direct" if n < FIRST_ENVELOPE_PEAK else "envelope"}
            if errors:
                mc = []
                for _ in range(draws):
                    noisy = dl + rng.normal(0.0, err)
                    try:
                        mc.append(_measure_one(ell, noisy, err, n))
                    except RuntimeError:
                        continue  # a failed fit on an extreme noise draw
                mc = np.array(mc)
                rec["ell_stat"] = float(np.std(mc[:, 0]))
                rec["dl_stat"] = float(np.std(mc[:, 1]))
                rec["n_mc"] = int(len(mc))
                sys_l = sys_d = 0.0
                if n >= FIRST_ENVELOPE_PEAK:
                    alt = np.array([_measure_one(ell, dl, err, n, envelope_order=o)
                                    for o in ENVELOPE_ORDERS])
                    sys_l = float(np.ptp(alt[:, 0]) / 2)
                    sys_d = float(np.ptp(alt[:, 1]) / 2)
                rec["ell_sys"], rec["dl_sys"] = sys_l, sys_d
                rec["ell_err"] = float(np.hypot(rec["ell_stat"], sys_l))
                rec["dl_err"] = float(np.hypot(rec["dl_stat"], sys_d))
            out.append(rec)
    return out


def measure_theory_peaks(theory: Theory, *, n_peaks: int = 8) -> list[dict]:
    """The same measurement applied to Planck's noise-free best-fit model.

    This is the control: any difference between it and the data is a
    difference between the sky and the model, not an artefact of the method.
    The model has no error bars, so the fit uses equal fractional weights.
    """
    m = theory.ell >= 30
    ell = theory.ell[m].astype(float)
    dl = theory.tt[m]
    spec = Spectrum(ell, dl, 0.01 * dl, 0.01 * dl)
    return measure_peaks(spec, n_peaks=n_peaks, errors=False)


def acoustic_scale(peaks: list[dict], *, use: int | None = None) -> dict:
    """Fit the peak positions as a harmonic series: ``l_n = l_A (n - phi)``.

    If the peaks are the overtones of one sound wave frozen at recombination,
    their positions are evenly spaced, with spacing ``l_A = pi / theta_s`` set
    by the angle ``theta_s`` that the sound horizon subtends on the sky.

    The offset ``phi`` is the other half of the story. Adiabatic initial
    conditions -- all components perturbed together, as inflation predicts --
    start each mode at maximum compression and give phi of about 0.25 (the
    exact value is nudged by the gravitational driving of the oscillations).
    The alternative, isocurvature initial conditions, starts the modes at zero
    and would give peaks in the ratio 1:3:5 (phi = 0.5 + shift). The measured
    offset settles which one nature chose.

    A caveat that the model-based fit (``cosmology_fit``) does not share:
    driving effects make phi drift slightly from peak to peak, so the spacing
    is not perfectly constant and l_A from this fit carries a ~1% systematic.
    """
    rows = [p for p in peaks if use is None or p["n"] <= use]
    n = np.array([p["n"] for p in rows], float)
    ell = np.array([p["ell"] for p in rows], float)
    sig = np.array([p.get("ell_err", 1.0) or 1.0 for p in rows], float)
    design = np.column_stack([n, -np.ones_like(n)])  # l = l_A n - (l_A phi)
    w = 1.0 / sig**2
    cov = np.linalg.inv(design.T @ (design * w[:, None]))
    slope, offset = cov @ design.T @ (w * ell)
    phi = offset / slope
    # Error on phi = offset/slope by linear propagation.
    j = np.array([-offset / slope**2, 1.0 / slope])
    phi_err = float(np.sqrt(j @ cov @ j))
    theta = np.pi / slope
    return {
        "n_peaks": int(len(rows)),
        "ell_A": float(slope),
        "ell_A_err": float(np.sqrt(cov[0, 0])),
        "phi": float(phi),
        "phi_err": phi_err,
        "theta_deg": float(np.degrees(theta)),
        "theta_deg_err": float(np.degrees(theta) * np.sqrt(cov[0, 0]) / slope),
        "chi2": float(np.sum(w * (ell - design @ np.array([slope, offset])) ** 2)),
        "ratios_to_first": (ell / ell[0]).tolist(),
    }


def binned_chi2(spec: Spectrum, model: np.ndarray) -> dict:
    """Chi-squared of a binned spectrum against a binned model, diagonal errors."""
    r = (spec.dl - model) / spec.err
    return {"chi2": float(np.sum(r**2)), "n": int(len(r)),
            "worst_bin_ell": float(spec.ell[np.argmax(np.abs(r))]),
            "worst_bin_sigma": float(r[np.argmax(np.abs(r))])}
