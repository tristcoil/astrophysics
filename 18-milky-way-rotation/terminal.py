"""Measure the terminal velocity: the edge of the hydrogen emission.

Looking into the inner Galaxy at longitude l, the gas moving fastest toward or
away from us sits at the tangent point, where the line of sight grazes the
smallest orbit it meets (radius R0 |sin l|). In a perfectly cold disc the
emission would stop dead at that velocity. Real gas has random motions of
several km/s, which smear the cliff into a slope, so "where the emission ends"
needs a definition.

We follow McClure-Griffiths & Dickey (2007, 2016):

1. Average the spectra over |b| <= 0.5 deg.
2. Find the most extreme velocity where the average still exceeds a threshold
   brightness (10 K; also 20 K as a check). This only locates the cliff.
3. Fit the cliff with a sum of two complementary error functions sharing one
   centre: the shape a step takes after smoothing by Gaussian random motions
   with two different widths. The shared centre, where the emission has fallen
   to half its level just inside the edge, is the terminal velocity.

Unlike a pure threshold, step 3 does not depend on how bright the gas happens
to be at the edge; the threshold's influence is checked by re-running with
other thresholds, fitting windows, and latitude statistics.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.special import erfc

import config as cfg


def two_erfc(u, u0, a1, w1, a2, w2):
    """Emission beyond the edge (u increasing away from allowed velocities)."""
    return a1 * erfc((u - u0) / w1) + a2 * erfc((u - u0) / w2)


@dataclass
class EdgeFit:
    v_threshold: float  # threshold-crossing velocity, km/s (signed, LSRK)
    v_terminal: float   # fitted edge centre, km/s (signed, LSRK)
    v_terminal_err: float
    amplitudes: tuple[float, float]
    widths: tuple[float, float]  # erfc widths, km/s (sigma = width / sqrt 2)
    rms: float          # fit residual rms, K
    ok: bool


BAD = EdgeFit(np.nan, np.nan, np.nan, (np.nan, np.nan), (np.nan, np.nan), np.nan, False)


def threshold_crossing(u: np.ndarray, T: np.ndarray, threshold: float,
                       min_run: int = 3) -> float:
    """Largest u at which T falls through ``threshold``, requiring the
    emission to have stayed above it for ``min_run`` channels (so a single
    noise spike or small cloud beyond the edge cannot set it). Linear
    interpolation between channels. ``u`` must be increasing."""
    above = T >= threshold
    if not above.any():
        return np.nan
    # last index i with above[i-min_run+1 .. i] all True
    run = np.convolve(above.astype(int), np.ones(min_run, int), mode="full")[:above.size]
    good = np.flatnonzero(run >= min_run)
    if good.size == 0:
        return np.nan
    i = good[-1]
    if i + 1 >= u.size:
        return np.nan
    t0, t1 = T[i], T[i + 1]
    frac = (t0 - threshold) / (t0 - t1) if t0 != t1 else 0.0
    return float(u[i] + frac * (u[i + 1] - u[i]))


def fit_edge(v: np.ndarray, T: np.ndarray, sign: int, *,
             threshold: float = cfg.TERMINAL_THRESHOLD_K,
             w_in: float = cfg.TERMINAL_W_IN_KMS,
             w_out: float = cfg.TERMINAL_W_OUT_KMS) -> EdgeFit:
    """Terminal velocity from one latitude-averaged spectrum.

    ``sign`` is +1 in the first quadrant (edge at positive velocity) and -1 in
    the fourth. The fit uses velocities from ``w_in`` km/s inside the threshold
    crossing to ``w_out`` km/s beyond it.
    """
    u = sign * v
    order = np.argsort(u)
    u, T = u[order], np.asarray(T, dtype=float)[order]
    finite = np.isfinite(T)
    u, T = u[finite], T[finite]
    u_thr = threshold_crossing(u, T, threshold)
    if not np.isfinite(u_thr):
        return BAD
    sel = (u >= u_thr - w_in) & (u <= u_thr + w_out)
    if sel.sum() < 12:
        return BAD
    uu, tt = u[sel], T[sel]
    plateau = max(np.max(tt[uu <= u_thr]), threshold)
    p0 = [u_thr + 3.0, 0.4 * plateau, 6.0, 0.1 * plateau, 18.0]
    lo = [u_thr - w_in, 0.0, 0.5, 0.0, 3.0]
    hi = [u_thr + 25.0, 5 * plateau, 20.0, 5 * plateau, 80.0]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", OptimizeWarning)
            popt, pcov = curve_fit(two_erfc, uu, tt, p0=p0, bounds=(lo, hi), maxfev=20000)
    except (RuntimeError, ValueError):
        return BAD
    u0, a1, w1, a2, w2 = popt
    if w2 < w1:  # keep component 1 the narrow one
        a1, w1, a2, w2 = a2, w2, a1, w1
    resid = tt - two_erfc(uu, *popt)
    err = float(np.sqrt(pcov[0, 0])) if np.all(np.isfinite(pcov)) else np.nan
    ok = bool(lo[0] + 0.5 < u0 < hi[0] - 0.5)
    return EdgeFit(sign * u_thr, sign * u0, err, (a1, a2), (w1, w2),
                   float(np.sqrt(np.mean(resid ** 2))), ok)


def quadrant_sign(l_deg) -> np.ndarray:
    return np.where(np.sin(np.radians(l_deg)) >= 0, 1, -1)


def measure(lv: dict, l_range: tuple[float, float], **kwargs) -> dict:
    """Fit every longitude pixel in ``l_range`` of a latitude-averaged l-v
    diagram (from ``survey.lv_diagram``). Returns arrays keyed by field."""
    l_all, v, T = lv["l"], lv["v"], lv["T"]
    idx = np.flatnonzero((l_all >= l_range[0]) & (l_all <= l_range[1]))
    fits = [fit_edge(v, T[i], int(quadrant_sign(l_all[i])), **kwargs) for i in idx]
    return {
        "l": l_all[idx],
        "v_threshold": np.array([f.v_threshold for f in fits]),
        "v_terminal": np.array([f.v_terminal for f in fits]),
        "v_terminal_err": np.array([f.v_terminal_err for f in fits]),
        "a1": np.array([f.amplitudes[0] for f in fits]),
        "a2": np.array([f.amplitudes[1] for f in fits]),
        "w1": np.array([f.widths[0] for f in fits]),
        "w2": np.array([f.widths[1] for f in fits]),
        "rms": np.array([f.rms for f in fits]),
        "ok": np.array([f.ok for f in fits]),
    }


def half_power_edge(v: np.ndarray, T: np.ndarray, sign: int, *,
                    search: float = 40.0, threshold: float = 10.0) -> float:
    """An independent, model-free edge: where the emission falls to half of
    the highest peak within ``search`` km/s inside the threshold crossing
    (the definition of Sinha 1978). Used only as a cross-check."""
    u = sign * v
    order = np.argsort(u)
    u, T = u[order], np.asarray(T, float)[order]
    u_thr = threshold_crossing(u, T, threshold)
    if not np.isfinite(u_thr):
        return np.nan
    win = (u >= u_thr - search) & (u <= u_thr)
    i_pk = np.flatnonzero(win)[np.argmax(T[win])]
    half = 0.5 * T[i_pk]
    beyond = np.flatnonzero((np.arange(u.size) > i_pk) & (T < half))
    if beyond.size == 0:
        return np.nan
    j = beyond[0]
    frac = (T[j - 1] - half) / (T[j - 1] - T[j])
    return float(sign * (u[j - 1] + frac * (u[j] - u[j - 1])))
