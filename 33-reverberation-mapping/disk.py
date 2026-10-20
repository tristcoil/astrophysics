"""The accretion disk echoes too: continuum lags against wavelength (NGC 5548, 2014).

If the disk is heated from above by the far-UV/X-ray source near the black
hole, each ring re-emits at its own temperature: hot inner rings in the UV,
cooler outer rings in the optical. Light from farther out arrives later, and
for a standard thin disk (T ~ R^-3/4) the delay grows as lambda^(4/3).

We repeat Fausnaugh et al. (2016)'s measurement: every ground-based optical
band and the other HST continuum bands against the HST 1367 A light curve,
over the HST campaign only, each light curve detrended with a second-order
polynomial, only the 1367 A curve interpolated (0.1-day grid), lags in the
rest frame; then the fit tau = alpha [(lambda / 1367 A)^beta - 1]. The thin
disk's alpha follows their equation 8 with the black hole mass from this
project.
"""

from __future__ import annotations

import numpy as np
from scipy import optimize

import config as C
import iccf
import lightcurves as LC

#: Pivot wavelengths (A) used by Fausnaugh et al. 2016 (their Table 6).
PIVOT = {"1157.5": 1158, "1478.5": 1479, "1746": 1746, "u": 3472, "B": 4369, "g": 4776, "V": 5404, "r": 6176,
         "R": 6440, "i": 7648, "I": 8561, "z": 9157}
#: Fausnaugh et al. 2016, Table 6: ICCF centroid lags relative to 1367 A (rest frame, days).
FAUSNAUGH16 = {"1157.5": (-0.17, 0.16, 0.16), "1478.5": (0.15, 0.16, 0.18), "1746": (0.22, 0.19, 0.16),
               "u": (2.03, 0.39, 0.43), "B": (1.42, 0.33, 0.36), "g": (1.98, 0.29, 0.34), "V": (2.04, 0.20, 0.22),
               "r": (3.13, 0.46, 0.41), "R": (3.22, 0.29, 0.30), "i": (3.99, 0.29, 0.29), "I": (3.59, 0.54, 0.53),
               "z": (3.93, 0.40, 0.44)}
#: Their Table 7: alpha with beta fixed at 4/3 (all bands; without u and U).
FAUSNAUGH16_ALPHA43 = {"all": (0.43, 0.02), "no_u": (0.42, 0.02)}
FAUSNAUGH16_FREE = {"all": ((0.97, 0.24), (0.90, 0.12)), "no_u": ((0.79, 0.22), (0.99, 0.14))}
LAMBDA0 = 1367.0


def _detrend2(t, y):
    p = np.polyfit(t - t.mean(), y, 2)
    return y - np.polyval(p, t - t.mean()) + y.mean()


def _merge_same_time(lc: LC.LightCurve, tol: float = 0.02) -> LC.LightCurve:
    """Average points closer than ``tol`` days (several telescopes on one night)."""
    t, y, e = lc.t, lc.y, lc.e
    groups = np.concatenate([[0], np.cumsum(np.diff(t) > tol)])
    w = 1.0 / e**2
    sw = np.bincount(groups, w)
    return LC.LightCurve(np.bincount(groups, w * t) / sw, np.bincount(groups, w * y) / sw, 1 / np.sqrt(sw),
                         lc.label, lc.unit)


def disk_lags(n_mc: int = 500) -> dict:
    ref = LC.storm_hst_continuum(1367.0)
    t0, t1 = ref.t[0], ref.t[-1]
    bands = {"1157.5": LC.storm_hst_continuum(1157.5), "1478.5": LC.storm_hst_continuum(1478.5),
             "1746": LC.storm_hst_continuum(1746.0)}
    phot = LC.storm_photometry()
    for b in ("u", "B", "g", "V", "r", "R", "i", "I", "z"):
        bands[b] = _merge_same_time(phot[b])
    yr = _detrend2(ref.t, ref.y)
    lags = iccf.lag_grid(-6.0, 12.0, 0.1)
    out = {}
    for k, (name, lc) in enumerate(bands.items()):
        w = lc.window(t0, t1)
        yw = _detrend2(w.t, w.y)
        # Only the 1367 A curve is interpolated: the "interp1" direction.
        rng = np.random.default_rng(700 + k)
        r_obs, _ = iccf.ccf(ref.t, yr, w.t, yw, lags, mode="interp1")
        cents = []
        for _ in range(n_mc):
            a = iccf._rss(rng, ref.t, yr, ref.e)
            b = iccf._rss(rng, w.t, yw, w.e)
            ra, _ = iccf.ccf(a[0], rng.normal(a[1], a[2]), b[0], rng.normal(b[1], b[2]), lags, mode="interp1")
            pc = iccf.peak_and_centroid(lags, ra)
            if pc.ok:
                cents.append(pc.tau_cent)
        cents = np.array(cents) / (1 + C.Z_NGC5548)
        lo, med, hi = np.percentile(cents, [15.865525, 50, 84.134475])
        out[name] = dict(lam=PIVOT[name], tau=(float(med), float(med - lo), float(hi - med)), n=len(w),
                         r_max=float(r_obs.max()), published=FAUSNAUGH16[name])
    return out


def _fit(lam, tau, err, beta=None):
    x = np.asarray(lam, float) / LAMBDA0
    if beta is None:
        f = lambda xx, a, b: a * (xx**b - 1)  # noqa: E731
        p, cov = optimize.curve_fit(f, x, tau, p0=[0.5, 1.2], sigma=err, absolute_sigma=True)
        chi2 = float(np.sum(((tau - f(x, *p)) / err) ** 2))
        return dict(alpha=(float(p[0]), float(np.sqrt(cov[0, 0]))), beta=(float(p[1]), float(np.sqrt(cov[1, 1]))),
                    rho=float(cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])), chi2=chi2, dof=len(x) - 2)
    g = x**beta - 1
    w = 1 / err**2
    a = float(np.sum(w * g * tau) / np.sum(w * g * g))
    sa = float(1 / np.sqrt(np.sum(w * g * g)))
    chi2 = float(np.sum(((tau - a * g) / err) ** 2))
    return dict(alpha=(a, sa), beta=(beta, 0.0), chi2=chi2, dof=len(x) - 1)


def fits(lags: dict) -> dict:
    names = list(lags)
    lam = np.array([lags[n]["lam"] for n in names], float)
    tau = np.array([lags[n]["tau"][0] for n in names])
    err = np.array([0.5 * (lags[n]["tau"][1] + lags[n]["tau"][2]) for n in names])
    no_u = np.array([n != "u" for n in names])
    return {"all_free": _fit(lam, tau, err), "all_43": _fit(lam, tau, err, 4 / 3),
            "no_u_free": _fit(lam[no_u], tau[no_u], err[no_u]), "no_u_43": _fit(lam[no_u], tau[no_u], err[no_u], 4 / 3)}


def thin_disk_alpha(mbh_msun: float, mdot_edd: float = 0.1, eta: float = 0.1, kappa: float = 1.0,
                    X: float = 2.49, lam0: float = LAMBDA0) -> float:
    """Fausnaugh et al. 2016, eq. 8: the lag scale (days) of an irradiated thin disk."""
    h, k, sig = 6.62607015e-27, 1.380649e-16, 5.670374419e-5
    M = mbh_msun * C.MSUN_G
    ledd = C.L_EDD_PER_MSUN * mbh_msun
    term = (C.G_CGS * M / (8 * np.pi * sig)) * (ledd / (eta * C.C_CGS**2)) * (3 + kappa) * mdot_edd
    return float((X * k * lam0 * 1e-8 / (h * C.C_CGS)) ** (4 / 3) * term ** (1 / 3) / C.C_CGS / C.DAY_S)
