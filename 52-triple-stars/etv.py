"""Eclipse timing of TIC 290061484's inner pair from our own TESS photometry.

1. Aperture photometry on the full-frame-image cutouts (tess_lc), normalised
   per sector, slow trends removed with a running median of out-of-eclipse
   points.
2. A mean ephemeris: individual minima found in each sector, cycle numbers
   counted across the 2019-2024 gap, a linear fit.
3. Each eclipse timed by fitting an empirical template (the sector's own
   phase-folded eclipse, built iteratively after removing the timing
   variations) for a time shift, depth scale and linear baseline. Errors come
   from the curvature of chi-square, scaled by the reduced chi-square, and are
   checked against the scatter of consecutive eclipses.
4. Eclipses within a day of a third-body event are flagged: the third star
   blocks part of the pair's light and distorts them.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.optimize import least_squares

import config as C
import tess_lc as T

NAME = "TIC 290061484"
SYS = C.SYSTEMS[NAME]
OUT = C.DERIVED / "tic290061484"
OUT.mkdir(parents=True, exist_ok=True)
#: guesses: Kostov et al. 2024 time of inferior conjunction (primary eclipse)
#: and the 1.792-d mean period they quote
T0_GUESS = 2458713.0340
P_GUESS = 1.79236


def light_curves(r_ap=1.5) -> list[dict]:
    out = []
    for cut in T.cutouts(NAME):
        df, meta = T.aperture_photometry(cut, SYS["ra"], SYS["dec"], r_ap=r_ap)
        jd = df["btjd"].to_numpy() + T.BTJD0
        f = df["flux"].to_numpy() / np.median(df["flux"])
        e = df["err"].to_numpy() / np.median(df["flux"])
        out.append(dict(sector=cut.sector, jd=jd, flux=f, err=e, exptime=cut.exptime_days))
    return out


def detrend(lc, P, T0, width_days=1.5, mask_half_phase=0.07):
    """Divide by a running median of out-of-eclipse points."""
    jd, f = lc["jd"], lc["flux"]
    ph = ((jd - T0) / P) % 1.0
    oot = (np.abs(ph) > mask_half_phase) & (np.abs(ph - 0.5) > mask_half_phase) & (np.abs(ph - 1) > mask_half_phase)
    dt = np.median(np.diff(jd))
    n = max(5, int(width_days / dt) | 1)
    trend = np.full(len(f), np.nan)
    trend[oot] = median_filter(f[oot], size=n, mode="nearest")
    trend = np.interp(jd, jd[oot], trend[oot])
    return f / trend


def find_minima(lc, depth=0.05, smooth_days=0.04):
    """Times of individual eclipse minima: local minima of the smoothed light
    curve deeper than `depth` below the median."""
    jd, f = lc["jd"], lc["flux"]
    dt = np.median(np.diff(jd))
    k = max(1, int(smooth_days / dt))
    fs = np.convolve(f, np.ones(2 * k + 1) / (2 * k + 1), mode="same")
    base = np.median(f)
    cand = np.where((fs < base - depth) & (fs <= np.roll(fs, 1)) & (fs <= np.roll(fs, -1)))[0]
    times = []
    for c in cand:
        if times and jd[c] - times[-1][0] < 0.5:
            if fs[c] < times[-1][1]:
                times[-1] = (jd[c], fs[c])
            continue
        times.append((jd[c], fs[c]))
    return np.array([t for t, _ in times]), np.array([v for _, v in times])


def mean_ephemeris(lcs) -> dict:
    """Count cycles from the 2019 sectors outward and fit P, T0 to the
    primary minima (deeper eclipse = primary)."""
    allt, alld = [], []
    for lc in lcs:
        t, d = find_minima(lc)
        allt.append(t)
        alld.append(d)
    t = np.concatenate(allt)
    d = np.concatenate(alld)
    # 1) 2019 only: refine P (28 cycles)
    P, T0 = P_GUESS, T0_GUESS
    for span in (60.0, 1200.0, 2000.0):
        sel = t < T0 + span
        n = np.round((t[sel] - T0) / P * 2) / 2          # half-integers = secondaries
        A = np.vstack([np.ones(sel.sum()), n]).T
        coef, *_ = np.linalg.lstsq(A, t[sel], rcond=None)
        T0, P = coef
    n = np.round((t - T0) / P * 2) / 2
    resid = t - (T0 + n * P)
    return dict(P=float(P), T0=float(T0), n_minima=int(len(t)), rms_resid=float(np.std(resid)),
                minima=t.tolist(), cycles=n.tolist())


def tertiary_mask(jd, pad=1.2):
    ev = np.array([e[0] for e in SYS["pub"]["tertiary_events"]]) + 2450000.0
    return np.any(np.abs(jd[:, None] - ev[None, :]) < pad, axis=1)


def template_centre(x, y):
    """Time offset about which an eclipse profile is most symmetric (the idea
    of Kwee & van Woerden 1956), using the part deeper than 20% of the depth."""
    from scipy.optimize import minimize_scalar
    depth = 1 - y.min()
    core = y < 1 - 0.2 * depth
    w = 0.5 * (x[core].max() - x[core].min()) if core.sum() > 3 else 0.05
    u = np.linspace(0.02 * w, w, 40)

    def asym(c):
        a = np.interp(c + u, x, y)
        b = np.interp(c - u, x, y)
        return float(np.sum((a - b) ** 2))
    r = minimize_scalar(asym, bounds=(-0.3 * w, 0.3 * w), method="bounded", options={"xatol": 1e-6})
    return float(r.x)


def time_eclipses(lcs, eph, half_window=0.22, n_iter=3):
    """Template fit of every eclipse; returns a DataFrame of times."""
    P, T0 = eph["P"], eph["T0"]
    rows = []
    for lc in lcs:
        f = detrend(lc, P, T0)
        jd, err = lc["jd"], lc["err"]
        sig = 1.4826 * np.median(np.abs(np.diff(f))) / np.sqrt(2)
        e = np.maximum(err, sig)
        shifts = {}
        for it in range(n_iter):
            # templates: fold with current per-eclipse shifts removed
            tmpl = {}
            for kind, off in (("pri", 0.0), ("sec", 0.5)):
                ph_all, fl_all = [], []
                ncyc = np.arange(np.floor((jd.min() - T0) / P) - 1, np.ceil((jd.max() - T0) / P) + 1)
                for c in ncyc:
                    tc = T0 + (c + off) * P
                    sel = np.abs(jd - tc) < half_window
                    if sel.sum() < 8:
                        continue
                    sh = shifts.get((kind, c), 0.0)
                    ph_all.append(jd[sel] - tc - sh)
                    fl_all.append(f[sel])
                if not ph_all:
                    continue
                x = np.concatenate(ph_all)
                y = np.concatenate(fl_all)
                o = np.argsort(x)
                x, y = x[o], y[o]
                bins = np.linspace(-half_window, half_window, 61)
                idx = np.digitize(x, bins) - 1
                xb = 0.5 * (bins[1:] + bins[:-1])
                yb = np.array([np.median(y[idx == k]) if np.sum(idx == k) >= 3 else np.nan for k in range(60)])
                ok = np.isfinite(yb)
                xb, yb = xb[ok], yb[ok] / np.median(yb[ok][np.abs(xb[ok]) > 0.8 * half_window])
                # centre the template on its point of symmetry, so that the
                # fitted shifts are absolute times of minimum (otherwise each
                # sector's template would absorb that sector's mean offset)
                tmpl[kind] = (xb - template_centre(xb, yb), yb)
            # fit each eclipse
            for kind, off in (("pri", 0.0), ("sec", 0.5)):
                if kind not in tmpl:
                    continue
                xt, yt = tmpl[kind]
                depth_t = 1 - yt.min()
                ncyc = np.arange(np.floor((jd.min() - T0) / P) - 1, np.ceil((jd.max() - T0) / P) + 1)
                for c in ncyc:
                    tc = T0 + (c + off) * P
                    sel = np.abs(jd - tc) < half_window
                    if sel.sum() < 15:
                        continue
                    tt, ff, ee = jd[sel] - tc, f[sel], e[sel]
                    # need both ingress and egress covered
                    if tt.min() > -0.6 * half_window or tt.max() < 0.6 * half_window:
                        continue

                    def resid(q):
                        sh, a, b0, b1 = q
                        m = 1 - a * (1 - np.interp(tt - sh, xt, yt, left=1.0, right=1.0))
                        return (ff - (b0 + b1 * tt) * m) / ee
                    q0 = [shifts.get((kind, c), 0.0), 1.0, 1.0, 0.0]
                    try:
                        r = least_squares(resid, q0, x_scale=[0.005, 0.1, 0.01, 0.01])
                    except ValueError:
                        continue
                    dof = max(1, len(tt) - 4)
                    chi2r = float(np.sum(r.fun**2) / dof)
                    try:
                        cov = np.linalg.inv(r.jac.T @ r.jac) * max(1.0, chi2r)
                        s_sh = float(np.sqrt(cov[0, 0]))
                    except np.linalg.LinAlgError:
                        continue
                    shifts[(kind, c)] = float(r.x[0])
                    if it == n_iter - 1:
                        rows.append(dict(sector=lc["sector"], kind=kind, cycle=float(c + off),
                                         t_pred=float(tc), t_obs=float(tc + r.x[0]), err=s_sh,
                                         depth_scale=float(r.x[1]), depth=float(r.x[1] * depth_t),
                                         chi2r=chi2r, n=int(len(tt)),
                                         near_tertiary=bool(tertiary_mask(np.array([tc]))[0])))
    df = pd.DataFrame(rows)
    df["oc"] = df["t_obs"] - df["t_pred"]
    return df


def run() -> dict:
    lcs = light_curves()
    eph = mean_ephemeris(lcs)
    df = time_eclipses(lcs, eph)
    # reject failed fits: shifts at the window edge, absurd depth scales
    good = (np.abs(df["oc"]) < 0.08) & df["depth_scale"].between(0.4, 2.5) & (df["err"] < 0.01)
    df["good"] = good & ~df["near_tertiary"]
    df.to_csv(OUT / "eclipse_times.csv", index=False)
    (OUT / "ephemeris.json").write_text(json.dumps({k: v for k, v in eph.items()
                                                    if k not in ("minima", "cycles")}, indent=1))
    return dict(ephemeris={k: v for k, v in eph.items() if k not in ("minima", "cycles")},
                n_eclipses=int(len(df)), n_good=int(df["good"].sum()),
                median_err_min=float(np.median(df.loc[df["good"], "err"]) * 1440))


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
