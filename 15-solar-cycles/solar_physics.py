"""The formulas this project rests on, kept free of any data so they can be tested.

* SILSO's 13-month smoothing, and finding cycle minima and maxima in it.
* The Hathaway, Wilson & Reichmann (1994) cycle shape, whose fitted start time t0
  is the clock Hathaway (2011) measured the drift of the sunspot zones against.
* The exponential "standard law" for that drift.
* Solar rotation: synodic vs sidereal periods, the Carrington frame, and the
  differential-rotation law Omega = A + B sin^2(latitude).
* Foreshortening: how a spot's measured area depends on where it sits on the disc.
* Small statistics helpers: correlations with bootstrap intervals and a straight-line
  fit with a prediction interval, used for every "rule of the solar cycle".

``test_physics.py`` pins each of these to a case with a known answer.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, signal, stats

# --- Time ---------------------------------------------------------------------------------

#: Earth's sidereal orbital period (days). The Sun's apparent rotation as seen from
#: the moving Earth (synodic) is slower than its true (sidereal) rotation by this.
EARTH_SIDEREAL_YEAR_DAYS = 365.25636
JD_UNIX_EPOCH = 2440587.5  # JD of 1970-01-01 00:00 UT


def decimal_year(dates) -> np.ndarray:
    """Calendar dates (anything numpy can make datetime64 from) -> decimal years.

    Uses the actual length of each year, so 1 July is 2000.497 in a leap year.
    """
    d = np.asarray(dates, dtype="datetime64[s]")
    years = d.astype("datetime64[Y]")
    start = years.astype("datetime64[s]")
    end = (years + 1).astype("datetime64[s]")
    frac = (d - start).astype(float) / (end - start).astype(float)
    return years.astype(int) + 1970 + frac


def julian_date(dates) -> np.ndarray:
    d = np.asarray(dates, dtype="datetime64[ms]")
    return d.astype("int64") / 86400000.0 + JD_UNIX_EPOCH


def datetime_from_jd(jd) -> np.ndarray:
    ms = np.round((np.asarray(jd, dtype=float) - JD_UNIX_EPOCH) * 86400000.0)
    return ms.astype("int64").astype("datetime64[ms]")


# --- The 13-month smoothed sunspot number ---------------------------------------------------

#: The "tapered boxcar": 13 months, the two end months at half weight, total 12.
SMOOTHING_WEIGHTS = np.r_[0.5, np.ones(11), 0.5] / 12.0


def smooth_13_month(monthly) -> np.ndarray:
    """SILSO's 13-month smoothed series from monthly means.

    Each value is centred on its month. The first and last six months, and any
    window that touches a missing (NaN) month, come out NaN -- as SILSO leaves them.
    """
    x = np.asarray(monthly, dtype=float)
    out = np.full_like(x, np.nan)
    if x.size < 13:
        return out
    windows = np.lib.stride_tricks.sliding_window_view(x, 13)
    out[6:-6] = windows @ SMOOTHING_WEIGHTS
    return out


@dataclass(frozen=True)
class Extrema:
    minima: np.ndarray  # indices into the series
    maxima: np.ndarray  # maxima[i] follows minima[i]; the last may be provisional


def find_cycle_extrema(smoothed, months_per_step: float = 1.0,
                       min_separation_yr: float = 7.0,
                       min_prominence: float = 10.0) -> Extrema:
    """Cycle minima and maxima of a smoothed sunspot series.

    A minimum is the lowest point of a trough that is at least ``min_prominence``
    deep and ``min_separation_yr`` from any deeper one. A maximum is the highest
    point between two consecutive minima (or after the last one). Ties -- a run of
    identical lowest values, as in the Dalton minimum where the smoothed number
    sits at 0.0 -- resolve to the middle of the run.
    """
    s = np.asarray(smoothed, dtype=float)
    valid = np.isfinite(s)
    filled = np.where(valid, s, np.nanmax(s))  # gaps can never be minima
    distance = max(1, int(round(min_separation_yr * 12 / months_per_step)))
    # For a flat-bottomed trough find_peaks reports the middle sample (rounded down).
    troughs, _ = signal.find_peaks(-filled, distance=distance, prominence=min_prominence)
    minima = np.asarray(troughs, dtype=int)
    maxima = []
    for k, i in enumerate(minima):
        j = minima[k + 1] if k + 1 < len(minima) else len(s)
        seg = np.where(valid[i:j], s[i:j], -np.inf)
        maxima.append(int(i + np.argmax(seg)))
    return Extrema(minima=minima, maxima=np.array(maxima, dtype=int))


# --- The shape of a cycle (Hathaway, Wilson & Reichmann 1994; Hathaway 2011) ----------------


def shape_width_months(rmax_v2, b0: float = 27.12, b1: float = 25.15,
                       v1_over_v2: float = 0.6):
    """Hathaway's width-amplitude relation, b = b0 + b1 (100 / Rmax)^(1/4) months.

    It was calibrated on version-1 sunspot numbers, so a version-2 amplitude is
    first put on that scale.
    """
    rmax_v1 = np.asarray(rmax_v2, dtype=float) * v1_over_v2
    return b0 + b1 * (100.0 / rmax_v1) ** 0.25


def cycle_shape(t_months, t0_months: float, rmax: float, b_months: float | None = None,
                c: float = 0.71):
    """F(t) = Rmax * 2 x^3 / (exp(x^2) - c), x = (t - t0) / b; zero before t0.

    With c = 0.71 the peak is 1.007 Rmax, reached at x = 1.09, i.e. ~1.09 b after t0.
    """
    b = shape_width_months(rmax) if b_months is None else b_months
    x = (np.asarray(t_months, dtype=float) - t0_months) / b
    xp = np.clip(x, 0.0, None)
    with np.errstate(over="ignore"):
        f = rmax * 2.0 * xp**3 / (np.exp(xp**2) - c)
    return np.where(x > 0, f, 0.0)


def fit_cycle_shape(t_months, sn, t0_guess: float, rmax_guess: float,
                    t0_bounds: tuple[float, float]) -> dict:
    """Least-squares fit of ``cycle_shape`` (free t0 and Rmax, width from the
    relation) to monthly sunspot numbers. Returns the parameters and 1-sigma errors."""
    t = np.asarray(t_months, dtype=float)
    y = np.asarray(sn, dtype=float)

    def model(tt, t0, rmax):
        return cycle_shape(tt, t0, rmax)

    p, cov = optimize.curve_fit(
        model, t, y, p0=(t0_guess, rmax_guess),
        bounds=((t0_bounds[0], 5.0), (t0_bounds[1], 1000.0)), maxfev=20000)
    resid = y - model(t, *p)
    err = np.sqrt(np.diag(cov))
    return {"t0_months": float(p[0]), "rmax": float(p[1]),
            "t0_err_months": float(err[0]), "rmax_err": float(err[1]),
            "b_months": float(shape_width_months(p[1])),
            "rms_residual": float(np.sqrt(np.mean(resid**2)))}


# --- The equatorward drift of the sunspot zones ----------------------------------------------


def drift_law(t_months, lambda0_deg: float, tau_months: float):
    """Hathaway's standard law: lambda(t) = lambda0 exp(-t / tau), t from cycle start."""
    return lambda0_deg * np.exp(-np.asarray(t_months, dtype=float) / tau_months)


def fit_drift_law(t_months, lat_deg, weights=None) -> dict:
    """Fit lambda0 and tau to centroid latitudes (absolute values, degrees)."""
    t = np.asarray(t_months, dtype=float)
    y = np.asarray(lat_deg, dtype=float)
    sigma = None if weights is None else 1.0 / np.sqrt(np.asarray(weights, dtype=float))
    p, cov = optimize.curve_fit(drift_law, t, y, p0=(28.0, 90.0), sigma=sigma,
                                absolute_sigma=False, maxfev=20000)
    resid = y - drift_law(t, *p)
    err = np.sqrt(np.diag(cov))
    return {"lambda0_deg": float(p[0]), "tau_months": float(p[1]),
            "lambda0_err": float(err[0]), "tau_err": float(err[1]),
            "rms_deg": float(np.sqrt(np.mean(resid**2))), "n": int(t.size)}


# --- Solar rotation -----------------------------------------------------------------------------


def synodic_from_sidereal(p_sidereal_days, p_orbit_days: float = EARTH_SIDEREAL_YEAR_DAYS):
    """The Earth orbits in the same sense the Sun spins, so a feature takes longer
    to come back to the same place *as seen from Earth*: 1/P_syn = 1/P_sid - 1/P_orb."""
    return 1.0 / (1.0 / np.asarray(p_sidereal_days, dtype=float) - 1.0 / p_orbit_days)


def sidereal_from_synodic(p_synodic_days, p_orbit_days: float = EARTH_SIDEREAL_YEAR_DAYS):
    return 1.0 / (1.0 / np.asarray(p_synodic_days, dtype=float) + 1.0 / p_orbit_days)


def period_days(omega_deg_per_day):
    return 360.0 / np.asarray(omega_deg_per_day, dtype=float)


def sidereal_rate(carrington_drift_deg_per_day, carrington_rate: float = 360.0 / 25.38):
    """A feature's sidereal angular velocity from its drift in Carrington longitude.

    Carrington longitude increases toward the west limb -- the direction of
    rotation -- so a feature spinning faster than the Carrington frame drifts to
    higher longitude: Omega = Omega_Carrington + dL/dt.
    """
    return carrington_rate + np.asarray(carrington_drift_deg_per_day, dtype=float)


def rotation_law(lat_deg, A: float, B: float):
    """Omega(latitude) = A + B sin^2(latitude), degrees per day."""
    return A + B * np.sin(np.radians(np.asarray(lat_deg, dtype=float))) ** 2


def fit_rotation_law(lat_deg, omega, sigma=None) -> dict:
    """Weighted linear least squares of Omega against sin^2(latitude)."""
    x = np.sin(np.radians(np.asarray(lat_deg, dtype=float))) ** 2
    y = np.asarray(omega, dtype=float)
    w = np.ones_like(y) if sigma is None else 1.0 / np.asarray(sigma, dtype=float) ** 2
    X = np.column_stack([np.ones_like(x), x])
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    A, B = cov @ (X.T @ (w * y))
    resid = y - (A + B * x)
    dof = max(1, y.size - 2)
    chi2 = float(np.sum(w * resid**2))
    if sigma is None:  # errors from the scatter itself
        cov = cov * chi2 / dof
    err = np.sqrt(np.diag(cov))
    return {"A": float(A), "B": float(B), "A_err": float(err[0]), "B_err": float(err[1]),
            "chi2": chi2, "dof": int(dof)}


def unwrap_degrees(lon_deg) -> np.ndarray:
    """Remove 360-degree jumps from a sequence of longitudes."""
    return np.degrees(np.unwrap(np.radians(np.asarray(lon_deg, dtype=float))))


# --- Sunspot areas --------------------------------------------------------------------------------


def corrected_area(observed_millionths_disk, distance_from_centre):
    """Foreshortening correction: millionths of the visible disc -> millionths of
    the solar hemisphere.

    A spot at heliocentric angle theta (sin theta = distance from disc centre in
    units of the radius) is seen at cos(theta) of its true area, and a hemisphere
    has twice the area of the disc: A_hem = A_disc / (2 cos theta).
    """
    r = np.asarray(distance_from_centre, dtype=float)
    cos_theta = np.sqrt(np.clip(1.0 - r**2, 1e-6, None))
    return np.asarray(observed_millionths_disk, dtype=float) / (2.0 * cos_theta)


def carrington_rotation_number(jd):
    """Fractional Carrington rotation number (the mean-rate estimate used by sunpy)."""
    return (np.asarray(jd, dtype=float) - 2398167.4) / 27.2753 + 1.0


# --- Statistics ------------------------------------------------------------------------------------


def correlation(x, y, n_boot: int = 10000, seed: int = 1) -> dict:
    """Pearson r with a two-sided p-value, Spearman rho, and a bootstrap 95% interval
    for r (resampling cycles, so it reflects how few cycles there are)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r, p = stats.pearsonr(x, y)
    rho, p_rho = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n_boot, x.size))
    xb, yb = x[idx], y[idx]
    xb = xb - xb.mean(axis=1, keepdims=True)
    yb = yb - yb.mean(axis=1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        rb = (xb * yb).sum(1) / np.sqrt((xb**2).sum(1) * (yb**2).sum(1))
    lo, hi = np.nanpercentile(rb, [2.5, 97.5])
    return {"r": float(r), "p": float(p), "rho": float(rho), "p_rho": float(p_rho),
            "r_ci95": [float(lo), float(hi)], "n": int(x.size)}


def line_fit(x, y) -> dict:
    """Ordinary least squares y = a + b x, with what a prediction interval needs."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    b, a = np.polyfit(x, y, 1)
    resid = y - (a + b * x)
    s = float(np.sqrt(np.sum(resid**2) / (n - 2)))
    return {"a": float(a), "b": float(b), "s": s, "n": int(n),
            "x_mean": float(x.mean()), "sxx": float(np.sum((x - x.mean()) ** 2))}


def predict(fit: dict, x0, level: float = 0.68) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Point prediction and the ``level`` prediction interval of a ``line_fit``
    for a *new* observation at x0 (scatter of the data plus uncertainty of the line)."""
    x0 = np.asarray(x0, dtype=float)
    yhat = fit["a"] + fit["b"] * x0
    se = fit["s"] * np.sqrt(1.0 + 1.0 / fit["n"] + (x0 - fit["x_mean"]) ** 2 / fit["sxx"])
    tq = stats.t.ppf(0.5 + level / 2.0, fit["n"] - 2)
    return yhat, yhat - tq * se, yhat + tq * se
