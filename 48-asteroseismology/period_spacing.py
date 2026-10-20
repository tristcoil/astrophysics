"""Gravity-mode period spacings: is the core burning helium or not?

A red giant's dipole (l = 1) modes are *mixed*: sound waves in the envelope
couple to gravity waves in the dense core. Pure gravity modes are evenly
spaced in **period** by Delta Pi_1, which depends on the core's size and
density. Hydrogen-shell-burning giants (red giant branch) have small, dense,
degenerate helium cores and Delta Pi_1 ~ 50-100 s; core-helium-burning
giants (red clump) have Delta Pi_1 ~ 250-350 s (Bedding et al. 2011, Nature
471, 608). On the surface the two kinds of star look the same.

Mixing with the pressure modes squeezes the periods near each pure p mode,
so the observed spacings are not constant. Mosser et al. (2015) and Vrard et
al. (2016, A&A 588, A87) undo this: with

    zeta(nu) = [1 + (q / N) / (q^2 cos^2 theta_p + sin^2 theta_p)]^-1,
    N = Dnu / (nu^2 Delta Pi_1),  theta_p = pi (nu - nu_p) / Dnu,

(q the coupling strength, nu_p the pure dipole p modes from the universal
pattern) the "stretched period" d tau = d(1/nu) / zeta makes the mixed modes
exactly evenly spaced in tau, by Delta Pi_1. A Fourier transform of the
spectrum as a function of tau then peaks at 1/Delta Pi_1. Because zeta barely
depends on the trial Delta Pi_1 and q, a few iterations converge.

Two hypotheses are tried for stars with Dnu < 9.5 uHz (red giant branch or
clump) and the self-consistent one with the higher peak is kept, as in Vrard
et al.; above 9.5 uHz only the red giant branch is possible.
"""

from __future__ import annotations

import numpy as np

from deltanu import universal_pattern


def zeta_function(nu, dnu, numax, eps, dpi1_s, q, *, d01=None, alpha=None) -> np.ndarray:
    """zeta(nu) from the asymptotic mixed-mode relation (all frequencies in uHz)."""
    up = universal_pattern(dnu)
    d01 = up["d01"] if d01 is None else d01
    alpha = up["alpha"] if alpha is None else alpha
    n_cont = nu / dnu - eps
    n_max = numax / dnu - eps
    # phase of the nearest pure dipole p mode: integer where nu = nu_p(n, l=1)
    phase = n_cont - 0.5 + d01 - 0.5 * alpha * (n_cont - n_max) ** 2
    theta = np.pi * phase
    N = (dnu * 1e-6) / ((nu * 1e-6) ** 2 * dpi1_s)
    return 1.0 / (1.0 + (q / N) / (q ** 2 * np.cos(theta) ** 2 + np.sin(theta) ** 2))


def stretched_period(nu, zeta) -> np.ndarray:
    """tau(nu) = integral of d(1/nu) / zeta, in seconds (increasing with period)."""
    inv = 1.0 / (nu * 1e-6)  # period, s
    # integrate along decreasing frequency so tau increases with period
    dP = -np.diff(inv)  # positive
    z_mid = 0.5 * (zeta[1:] + zeta[:-1])
    tau = np.concatenate([[0.0], np.cumsum(dP / z_mid)])
    return tau[-1] - tau  # same orientation as the period: larger at low frequency


def dipole_mask(nu, dnu, eps, numax, *, lo=0.06, hi=0.80) -> np.ndarray:
    """Keep the part of each order away from the radial and quadrupole modes.

    Vrard et al. 2016, their Eq. 9: reduced frequency in [0.06, 0.80].
    """
    up = universal_pattern(dnu)
    n_cont = nu / dnu - eps
    n_max = numax / dnu - eps
    x = (n_cont - 0.5 * up["alpha"] * (n_cont - n_max) ** 2) % 1.0
    return (x >= lo) & (x <= hi)


def tau_spectrum(tau, signal, dtau, periods_s):
    """Fourier power of signal(tau) at trial periods: |sum s exp(-2 pi i tau / P) dtau|^2.

    ``dtau`` is each frequency bin's own extent in tau (so masked-out stretches
    of spectrum, set to zero, leave no gaps or jumps in the integral).
    """
    f = 1.0 / np.asarray(periods_s)
    amp = np.exp(-2j * np.pi * np.outer(f, tau)) @ (signal * dtau)
    return np.abs(amp) ** 2


def _stretch(nu, dnu, numax, eps, dpi, q, up):
    z = zeta_function(nu, dnu, numax, eps, dpi, q, d01=up["d01"], alpha=up["alpha"])
    dnu_bin = nu[1] - nu[0]
    dtau = (dnu_bin * 1e-6) / ((nu * 1e-6) ** 2 * z)  # |d(1/nu)| / zeta per bin, s
    tau = np.cumsum(dtau[::-1])[::-1]  # increases towards low frequency, like the period
    return tau, dtau


#: Trial periods for the Fourier analysis of the stretched spectrum (s).
TRIAL_PERIODS = np.round(np.arange(30.0, 450.0001, 0.2), 4)
#: Short periods where no period-spacing signal is expected: the noise level.
NOISE_PERIODS = np.round(np.arange(10.0, 25.0001, 0.05), 4)
#: Detection threshold on the peak height over the mean Fourier power in the
#: hypothesis' period range. With ~100-300 independent trial periods per
#: range and exponentially distributed noise powers, a pure-noise maximum
#: exceeds 8 times the mean with a probability of ~3-10%, and 10 times with
#: ~0.5-1.5% (the threshold is discussed in the README).
DETECTION_THRESHOLD = 8.0


#: Period ranges allowed for each evolutionary hypothesis (s), after Vrard et
#: al. 2016's Fig. 1 and Mosser et al. 2014: hydrogen-shell burning (red
#: giant branch) 45-130 s; helium-core burning 130-420 s, which covers both
#: the red clump (~250-350 s) and the secondary clump (stars massive enough,
#: above ~2 Msun, to ignite helium without a degenerate core; ~150-250 s).
#: The two helium-burning groups are one hypothesis, so that their solutions
#: are compared on the same footing (separately, an alias of a red-clump
#: peak ~25-50 s below it wins the secondary-clump range). Each hypothesis
#: has several starting values and a typical coupling q.
HYPOTHESES = {"RGB": ((60.0, 75.0, 90.0), 0.12, (45.0, 130.0)),
              "helium burning": ((150.0, 190.0, 230.0, 270.0, 310.0, 350.0), 0.25, (130.0, 420.0))}
#: Below this Dnu (uHz) the measurement is not attempted: in luminous giants
#: the coupling is so weak that the dipole modes are almost pure p modes and
#: the gravity-dominated ones are invisible. Vrard et al. 2016's 6,111 stars
#: have Dnu > 3.3 uHz (99.5% of them).
MIN_DNU = 3.3
Q_GRID = (0.03, 0.05, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5)


def _iterate(nu_w, s_w, dnu, numax, eps, dpi, q, up, lo, hi, n_iter=10, half_width=0.2):
    """Local fixed-point iteration: re-stretch, re-measure within +-20%."""
    P = TRIAL_PERIODS
    for _ in range(n_iter):
        sel = (P >= max(lo, dpi * (1 - half_width))) & (P <= min(hi, dpi * (1 + half_width)))
        if sel.sum() < 3:
            return None
        tau, dtau = _stretch(nu_w, dnu, numax, eps, dpi, q, up)
        pw = tau_spectrum(tau, s_w, dtau, P[sel])
        j = int(np.argmax(pw))
        new = float(P[sel][j])
        at_edge = j == 0 or j == sel.sum() - 1
        converged = abs(new - dpi) < 0.05
        dpi = new
        if converged:
            return dpi if not at_edge else None
    return None


def measure_dpi1(nu, snr, dnu, numax, eps, sigma_env, *, n_orders: float = 3.5) -> dict:
    """Delta Pi_1 by the stretched-period method; returns the best hypothesis.

    ``snr`` is the background-divided spectrum. Frequencies within
    ``n_orders`` Dnu of numax are used; outside the dipole part of each order
    the signal is set to zero (Vrard et al.'s mask removes the radial and
    quadrupole modes).

    Each hypothesis (hydrogen-shell or helium-core burning) is iterated
    from a few starting values, each step searching within +-20% of the
    current value, as a fixed point: stretch with Delta Pi_1, measure the
    period of the resulting comb, re-stretch. As in Vrard et al. 2016 a
    hypothesis counts only if it converges to a genuine peak inside its own
    period range ("coherent"). The peak's height is measured against the mean
    Fourier power in that range: the stretched spectrum also carries
    long-period power from the order-to-order modulation of the mode heights,
    which a single global noise level would mistake for signal. Only the red
    giant branch is possible above Dnu = 9.5 uHz.
    """
    up = universal_pattern(dnu)
    if dnu < MIN_DNU:
        return {"detected": False, "attempted": False, "reason": f"Dnu < {MIN_DNU} uHz"}
    win = np.abs(nu - numax) < n_orders * dnu
    if win.sum() < 200:
        return {"detected": False, "attempted": False, "reason": "too little spectrum"}
    nu_w = nu[win]
    keep = dipole_mask(nu_w, dnu, eps, numax)
    s_w = np.where(keep, np.clip(snr[win] - 1.0, -1.0, 50.0), 0.0)
    names = ["RGB"] if dnu >= 9.5 else list(HYPOTHESES)
    P = TRIAL_PERIODS
    hyps = []
    for label in names:
        starts, q0, (lo, hi) = HYPOTHESES[label]
        in_range = (P >= lo) & (P <= hi)
        best = None
        for start in starts:
            dpi = _iterate(nu_w, s_w, dnu, numax, eps, start, q0, up, lo, hi)
            if dpi is None:
                continue
            # coupling: pick q maximising the peak at fixed Delta Pi_1, then re-converge
            near = np.abs(P - dpi) < 0.03 * dpi
            qs = []
            for qq in Q_GRID:
                tau, dtau = _stretch(nu_w, dnu, numax, eps, dpi, qq, up)
                qs.append(tau_spectrum(tau, s_w, dtau, P[near]).max())
            q = Q_GRID[int(np.argmax(qs))]
            dpi2 = _iterate(nu_w, s_w, dnu, numax, eps, dpi, q, up, lo, hi)
            dpi = dpi if dpi2 is None else dpi2
            tau, dtau = _stretch(nu_w, dnu, numax, eps, dpi, q, up)
            pw = tau_spectrum(tau, s_w, dtau, P)
            k = int(np.argmin(np.abs(P - dpi)))
            off_peak = in_range & (np.abs(P - dpi) > 0.05 * dpi)
            local = float(np.mean(pw[off_peak]))
            h = float(pw[k] / local)
            if best is None or h > best["height"]:
                noise = float(np.median(tau_spectrum(tau, s_w, dtau, NOISE_PERIODS)) / np.log(2.0))
                err = dpi ** 2 / float(np.ptp(tau)) / max(np.sqrt(h), 1.0)
                best = {"hypothesis": label, "dpi1": float(dpi), "q": float(q), "height": h,
                        "height_vs_white": float(pw[k] / noise), "coherent": True,
                        "dpi1_err": float(err), "power": (pw / local).astype(np.float32)}
        if best is None:
            best = {"hypothesis": label, "dpi1": np.nan, "q": np.nan, "height": 0.0,
                    "height_vs_white": 0.0, "coherent": False, "dpi1_err": np.nan,
                    "power": np.zeros_like(P, dtype=np.float32)}
        hyps.append(best)
    ok = [h for h in hyps if h["coherent"]] or hyps
    best = dict(max(ok, key=lambda h: h["height"]))
    best["detected"] = bool(best["coherent"] and best["height"] >= DETECTION_THRESHOLD)
    best["attempted"] = True
    best["all"] = [{k: v for k, v in h.items() if k != "power"} for h in hyps]
    return best


def measure_for_star(arrays: dict, star: dict) -> dict:
    """Delta Pi_1 from a star's cached spectrum (``stars.load_star``)."""
    nu = arrays["nu"]
    snr = arrays["psd"] / arrays["bg"]
    sig = star["envelope_fwhm"] / 2.3548
    res = measure_dpi1(nu, snr, star["dnu"], star["numax"], star["epsilon"], sig)
    res.pop("power", None)
    return res


def observed_period_spacing(nu, snr, dnu, numax, eps, *, n_orders: float = 3.5) -> dict:
    """Bedding et al. (2011)'s Delta P_obs: the same Fourier analysis, unstretched.

    The mean period spacing of the observed l = 1 modes, measured from the
    spectrum as a function of period 1/nu (zeta = 1). Smaller than Delta Pi_1
    because coupling squeezes the periods near the p modes.
    """
    win = np.abs(nu - numax) < n_orders * dnu
    nu_w = nu[win]
    keep = dipole_mask(nu_w, dnu, eps, numax)
    s_w = np.where(keep, np.clip(snr[win] - 1.0, -1.0, 50.0), 0.0)
    period = 1.0 / (nu_w * 1e-6)
    dper = (nu_w[1] - nu_w[0]) * 1e-6 / (nu_w * 1e-6) ** 2
    trial = np.round(np.arange(25.0, 250.0001, 0.2), 4)
    pw = tau_spectrum(period, s_w, dper, trial)
    k = int(np.argmax(pw))
    noise = np.median(tau_spectrum(period, s_w, dper, NOISE_PERIODS)) / np.log(2.0)
    return {"dp_obs": float(trial[k]), "height": float(pw[k] / noise), "trial": trial, "power": pw}
