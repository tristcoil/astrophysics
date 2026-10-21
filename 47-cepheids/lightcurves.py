"""Our own light-curve analysis: blind period search, Fourier fits, intensity means.

A Cepheid's light curve repeats exactly (to ~10^-6 in period over OGLE's
decade), so a truncated Fourier series in phase,

    m(t) = A0 + sum_k [a_k cos(2 pi k t / P) + b_k sin(2 pi k t / P)],

describes it to the photometric noise. From the fit we take:

* the **period**, found blind with a Lomb-Scargle periodogram and then refined
  by minimising the Fourier fit's chi^2;
* the **intensity-mean magnitude**: the magnitude of the time-averaged *flux*
  (not the average magnitude), which is what the Leavitt law is built on;
* amplitude, the time of maximum light, and the Fourier shape parameters
  R21 = A2/A1, phi21 = phi2 - 2 phi1 (and R31, phi31), in OGLE's convention
  m = A0 + sum A_k cos(k omega (t - T0) + phi_k).
"""

from __future__ import annotations

import numpy as np

#: Fewer points than this and we do not attempt a fit (a handful of OGLE stars
#: have only a dozen epochs).
MIN_POINTS = 30

#: Frequency grid of the blind search: periods from 0.2 to 250 days. OGLE's
#: ~3,660-day baseline resolves 1/3660 d^-1; we oversample that tenfold.
F_MIN, F_MAX = 1 / 250.0, 1 / 0.2
OVERSAMPLE = 10


def design(t, period, order):
    w = 2 * np.pi * t / period
    cols = [np.ones_like(t)]
    for k in range(1, order + 1):
        cols += [np.cos(k * w), np.sin(k * w)]
    return np.column_stack(cols)


def fourier_lsq(t, y, e, period, order):
    """Weighted linear least squares; returns (coefficients, covariance, chi2)."""
    A = design(t, period, order)
    wgt = 1.0 / e
    coef, *_ = np.linalg.lstsq(A * wgt[:, None], y * wgt, rcond=None)
    r = (y - A @ coef) * wgt
    cov = np.linalg.pinv((A * wgt[:, None]).T @ (A * wgt[:, None]))
    return coef, cov, float(r @ r)


def model(coef, phase):
    """Magnitude at phase(s) in [0, 1) for coefficients from ``fourier_lsq`` (t0 = 0)."""
    order = (len(coef) - 1) // 2
    out = np.full_like(np.asarray(phase, float), coef[0])
    for k in range(1, order + 1):
        out = out + coef[2 * k - 1] * np.cos(2 * np.pi * k * phase) + coef[2 * k] * np.sin(2 * np.pi * k * phase)
    return out


def shape_parameters(coef):
    """Amplitudes A_k and phases phi_k for m = A0 + sum A_k cos(k w t + phi_k)."""
    order = (len(coef) - 1) // 2
    a = coef[1::2][:order]
    b = coef[2::2][:order]
    amp = np.hypot(a, b)
    # a cos x + b sin x = A cos(x + phi) with A cos phi = a, -A sin phi = b
    phi = np.arctan2(-b, a)
    return amp, phi


def fourier_shape(coef):
    """(R21, phi21, R31, phi31) with phases in [0, 2 pi), OGLE's cosine convention."""
    amp, phi = shape_parameters(coef)
    out = [np.nan] * 4
    if len(amp) >= 2 and amp[0] > 0:
        out[0] = amp[1] / amp[0]
        out[1] = (phi[1] - 2 * phi[0]) % (2 * np.pi)
    if len(amp) >= 3 and amp[0] > 0:
        out[2] = amp[2] / amp[0]
        out[3] = (phi[2] - 3 * phi[0]) % (2 * np.pi)
    return tuple(out)


def intensity_mean(coef, n=2000):
    """Magnitude of the phase-averaged flux of the Fourier model."""
    ph = (np.arange(n) + 0.5) / n
    m = model(coef, ph)
    return float(-2.5 * np.log10(np.mean(10 ** (-0.4 * m))))


def extrema(coef, n=4000):
    """(brightest magnitude, faintest magnitude, phase of maximum light)."""
    ph = (np.arange(n) + 0.5) / n
    m = model(coef, ph)
    k = int(np.argmin(m))
    return float(m.min()), float(m.max()), float(ph[k])


def best_order(t, y, e, period, max_order=8):
    """Fourier order by the Bayesian information criterion (at least 1)."""
    n = len(t)
    best = (np.inf, 1)
    for order in range(1, max_order + 1):
        if 2 * order + 1 > n / 4:
            break
        _, _, chi2 = fourier_lsq(t, y, e, period, order)
        bic = chi2 + (2 * order + 1) * np.log(n)
        if bic < best[0]:
            best = (bic, order)
    return best[1]


def refine_period(t, y, e, period, order, span=None):
    """Minimise chi^2 over period near ``period``: golden-section search within +-2 grid steps."""
    from scipy.optimize import minimize_scalar
    baseline = t.max() - t.min()
    dp = (span or 2.0) * period ** 2 / (baseline * OVERSAMPLE)
    res = minimize_scalar(lambda p: fourier_lsq(t, y, e, p, order)[2],
                          bounds=(period - dp, period + dp), method="bounded",
                          options={"xatol": period * 1e-9})
    return float(res.x)


def period_error(t, y, e, period, order):
    """1-sigma period error from the chi^2 curvature, scaled by the reduced chi^2."""
    h = period * 2e-6
    c = [fourier_lsq(t, y, e, period + d, order)[2] for d in (-h, 0.0, h)]
    curv = (c[0] - 2 * c[1] + c[2]) / h ** 2
    dof = max(len(t) - (2 * order + 2), 1)
    scale = max(c[1] / dof, 1.0)
    return float(np.sqrt(2.0 * scale / curv)) if curv > 0 else np.nan


def clip(t, y, e, period, order, nsigma=4.0, iters=5):
    """Mask of points kept after iterative clipping on Fourier-fit residuals."""
    keep = np.ones(len(t), bool)
    for _ in range(iters):
        coef, _, _ = fourier_lsq(t[keep], y[keep], e[keep], period, order)
        r = y - model(coef, (t / period) % 1.0)
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        new = np.abs(r) < nsigma * max(s, 1e-3) + 2 * e
        if np.array_equal(new, keep):
            break
        keep = new
    return keep


def blind_peaks(t, y, e, n_peaks=3):
    """Top Lomb-Scargle peaks (periods) between 0.2 and 250 days, distinct by >1%."""
    from astropy.timeseries import LombScargle
    baseline = t.max() - t.min()
    df = 1.0 / (baseline * OVERSAMPLE)
    freq = np.arange(F_MIN, F_MAX, df)
    power = LombScargle(t, y, e).power(freq, method="fast", assume_regular_frequency=True)
    order = np.argsort(power)[::-1]
    peaks = []
    for k in order:
        f = freq[k]
        if all(abs(f - g) / g > 0.01 for g in peaks):
            peaks.append(f)
        if len(peaks) == n_peaks:
            break
    return [1.0 / f for f in peaks]


def analyse(t, y, e, *, period_hint=None, blind=True, fix_period=None, max_order=8):
    """Full treatment of one light curve. Returns a dict (NaNs where impossible).

    The blind search's best peak is refined with a 2nd-order Fourier fit, then
    the order is chosen by BIC, outliers are clipped, and the period refined
    again. With ``blind=False`` the search starts from ``period_hint`` instead.
    With ``fix_period`` the period is held fixed (used for the sparser V-band
    light curves, which take the I-band period).
    """
    out = {"n": len(t), "P_blind": np.nan, "P": np.nan, "eP": np.nan, "order": 0, "mean": np.nan,
           "amp": np.nan, "max": np.nan, "min": np.nan, "T_max": np.nan, "R21": np.nan, "phi21": np.nan,
           "R31": np.nan, "phi31": np.nan, "rms": np.nan, "n_used": 0}
    if len(t) < MIN_POINTS:
        return out
    e = np.maximum(e, 0.003)  # OGLE's smallest quoted errors understate the scatter of bright stars
    if fix_period is not None:
        p = float(fix_period)
    else:
        if blind:
            p0 = blind_peaks(t, y, e, n_peaks=1)[0]
            out["P_blind"] = p0
        else:
            p0 = period_hint
        p = refine_period(t, y, e, p0, order=2)
    order = best_order(t, y, e, p, max_order)
    keep = clip(t, y, e, p, order)
    tk, yk, ek = t[keep], y[keep], e[keep]
    order = best_order(tk, yk, ek, p, max_order)
    if fix_period is None:
        p = refine_period(tk, yk, ek, p, order, span=0.5)
    coef, cov, chi2 = fourier_lsq(tk, yk, ek, p, order)
    resid = yk - model(coef, (tk / p) % 1.0)
    mx, mn, ph_max = extrema(coef)
    r21, phi21, r31, phi31 = fourier_shape(_shift_to_max(coef, ph_max))
    # epoch of maximum light nearest the middle of the data
    t_mid = np.median(tk)
    t_max = p * (np.floor(t_mid / p) + ph_max)
    out.update({"P": p, "eP": np.nan if fix_period is not None else period_error(tk, yk, ek, p, order),
                "order": order,
                "mean": intensity_mean(coef), "amp": mn - mx, "max": mx, "min": mn, "T_max": t_max,
                "R21": r21, "phi21": phi21, "R31": r31, "phi31": phi31,
                "rms": float(np.std(resid)), "n_used": int(keep.sum())})
    return out


def _shift_to_max(coef, ph_max):
    """Coefficients re-expressed with phase zero at maximum light (phi21 is shift-invariant anyway)."""
    order = (len(coef) - 1) // 2
    out = coef.copy()
    for k in range(1, order + 1):
        a, b = coef[2 * k - 1], coef[2 * k]
        d = 2 * np.pi * k * ph_max
        # cos(k w (t + t_max)) = cos(kwt) cos d - sin(kwt) sin d, etc.
        out[2 * k - 1] = a * np.cos(d) + b * np.sin(d)
        out[2 * k] = -a * np.sin(d) + b * np.cos(d)
    return out


def fold(t, period, t0):
    """Phase in [0, 1) with maximum light at phase 0 if t0 is a time of maximum."""
    return ((t - t0) / period) % 1.0


# ---------------------------------------------------------------------------
# Batch processing of every OGLE light curve (parallel, checkpointed)
# ---------------------------------------------------------------------------

FIELDS = ("n", "n_used", "P_blind", "P", "eP", "order", "mean", "amp", "max", "min", "T_max",
          "R21", "phi21", "R31", "phi31", "rms")


def _one_star(args):
    gal, sid = args
    import ogle
    row = {"id": sid, "gal": gal}
    lc = ogle.read_lightcurve(gal, sid, "I")
    p_i = np.nan
    if lc is not None:
        r = analyse(*lc)
        row.update({f"I_{k}": r[k] for k in FIELDS})
        p_i = r["P"]
    lv = ogle.read_lightcurve(gal, sid, "V")
    if lv is not None and np.isfinite(p_i):
        r = analyse(*lv, fix_period=p_i, max_order=6)
        row.update({f"V_{k}": r[k] for k in FIELDS if k not in ("P_blind", "P", "eP")})
    elif lv is not None:  # saturated in I (the brightest Cepheids): search V blind
        r = analyse(*lv, max_order=6)
        row.update({f"V_{k}": r[k] for k in FIELDS})
    return row


def _init_worker():
    from threadpoolctl import threadpool_limits
    threadpool_limits(1)


def run_all(workers: int = 6, chunksize: int = 20):
    """Analyse every OGLE Magellanic Cepheid light curve; resumable via a JSON-lines checkpoint."""
    import json
    import multiprocessing as mp
    import config as C
    import ogle
    out = C.CACHE_DIR / "lightcurve_fits.jsonl"
    done = set()
    if out.exists():
        for line in out.read_text().splitlines():
            if line.strip():
                done.add(json.loads(line)["id"])
    todo = [(gal, sid) for gal in C.GALAXIES for sid in ogle.catalog(gal)["id"] if sid not in done]
    if not todo:
        return out
    print(f"Fitting {len(todo)} light curves with {workers} workers ({len(done)} already done) ...",
          flush=True)
    ctx = mp.get_context("fork")  # Python 3.14 defaults to forkserver
    with ctx.Pool(workers, initializer=_init_worker) as pool, open(out, "a") as fh:
        for k, row in enumerate(pool.imap_unordered(_one_star, todo, chunksize=chunksize), 1):
            clean = {key: (None if isinstance(v, float) and not np.isfinite(v) else v) for key, v in row.items()}
            fh.write(json.dumps(clean) + "\n")
            if k % 500 == 0:
                fh.flush()
                print(f"  {k}/{len(todo)}", flush=True)
    return out


def load_fits() -> dict:
    """Light-curve fit results keyed by OGLE ID."""
    import json
    import config as C
    rows = {}
    for line in (C.CACHE_DIR / "lightcurve_fits.jsonl").read_text().splitlines():
        if line.strip():
            r = json.loads(line)
            rows[r["id"]] = {k: (np.nan if v is None else v) for k, v in r.items()}
    return rows


if __name__ == "__main__":
    import sys
    run_all(int(sys.argv[1]) if len(sys.argv) > 1 else 6)


def model_for_plot(t, y, e, period, order=None, max_order=8):
    """Fourier model of one light curve at a fixed period, phased to maximum light.

    Returns (data phase, data mag, model phase grid, model mag, intensity mean).
    """
    e = np.maximum(e, 0.003)
    order = order or best_order(t, y, e, period, max_order)
    keep = clip(t, y, e, period, order)
    coef, _, _ = fourier_lsq(t[keep], y[keep], e[keep], period, order)
    _, _, ph_max = extrema(coef)
    grid = np.linspace(0, 1, 400)
    mod = model(coef, (grid + ph_max) % 1.0)
    phase = ((t[keep] / period) - ph_max) % 1.0
    return phase, y[keep], grid, mod, intensity_mean(coef)
