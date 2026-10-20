"""Checks of the physics and the methods (synthetic inputs are for tests only).

    python -m pytest test_physics.py -q      (or: python test_physics.py)
"""

from __future__ import annotations

import numpy as np

import background as B
import config as C
import deltanu as D
import lightcurves as L
import period_spacing as P
import scaling as SC
import sonify as SO
import spectrum as S


def _gappy_times(n=20000, dt=C.KEPLER_LC_CADENCE_S, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(n) * dt / 86400.0
    keep = rng.random(n) > 0.1  # 10% of cadences missing at random
    keep[5000:5600] = False  # and one long gap
    return t[keep]


def test_parseval_white_noise():
    t = _gappy_times()
    y = np.random.default_rng(1).normal(0, 80.0, len(t))
    sp = S.power_density(t, y, cadence_s=C.KEPLER_LC_CADENCE_S)
    ratio = np.sum(sp.psd) * sp.df / np.var(y)
    assert abs(ratio - 1) < 0.02, ratio


def test_sine_power_is_half_amplitude_squared():
    """A sine of amplitude A carries power A^2/2 (its variance)."""
    t = np.arange(40000) * C.KEPLER_LC_CADENCE_S / 86400.0
    A, f0 = 50.0, 60.0  # ppm, uHz
    y = A * np.sin(2 * np.pi * f0 * 1e-6 * t * 86400.0)
    sp = S.power_density(t, y, cadence_s=C.KEPLER_LC_CADENCE_S, oversample=4)
    m = np.abs(sp.nu - f0) < 0.2
    assert abs(np.sum(sp.psd[m]) * sp.df / (A ** 2 / 2) - 1) < 0.03


def test_apodization_limits():
    assert S.apodization(0.0, 283.0) == 1.0
    assert abs(S.apodization(283.0, 283.0) - (2 / np.pi) ** 2) < 1e-12


def test_background_components_integrate_to_rms_squared():
    nu = np.linspace(1e-4, 5e4, 2_000_001)
    for f, a, b in ((B.harvey, 120.0, 30.0), (B.lorentzian, 40.0, 0.5)):
        integral = np.trapezoid(f(nu, a, b), nu)
        assert abs(integral / a ** 2 - 1) < 0.01, (f.__name__, integral / a ** 2)


def test_filter_transfer_matches_detrending():
    """Detrending a sine reduces its power by filter_transfer(nu)."""
    dt = C.KEPLER_LC_CADENCE_S
    n = 70000
    t = np.arange(n) * dt / 86400.0
    for f0 in (0.8, 1.5, 3.0):  # uHz
        y = 1.0 + 1e-3 * np.sin(2 * np.pi * f0 * 1e-6 * t * 86400.0)
        ts = L.prepare_segments([(t, y, np.arange(n))], cadence_s=dt, trend_sigma_days=5.0, clip_sigma=50)
        mid = slice(n // 4, 3 * n // 4)  # away from the edges
        amp = np.std(ts.flux[mid]) * np.sqrt(2) / 1e3
        assert abs(amp ** 2 / L.filter_transfer(f0, 5.0) - 1) < 0.03, (f0, amp ** 2, L.filter_transfer(f0, 5.0))


def test_scaling_relations_give_the_sun():
    kw = dict(numax_sun=3090.0, dnu_sun=135.1, teff_sun=5777.0)
    assert abs(SC.mass(3090.0, 135.1, 5777.0, **kw) - 1) < 1e-12
    assert abs(SC.radius(3090.0, 135.1, 5777.0, **kw) - 1) < 1e-12
    assert abs(SC.logg(3090.0, 5777.0, numax_sun=3090.0, teff_sun=5777.0) - SC.LOGG_SUN) < 1e-12


def test_apokasc3_masses_reproduced():
    """Their own inputs + our formula (their convention) give their masses."""
    import catalogs as K

    a = K.apokasc3()
    g = a[a["category"] == "Gold"].dropna(subset=["numax", "dnu", "f_dnu", "f_numax", "mass", "teff"])
    m = SC.mass(g["numax"], g["dnu"], g["teff"], numax_sun=C.NUMAX_SUN_APOKASC, dnu_sun=C.DNU_SUN_APOKASC,
                teff_sun=C.TEFF_SUN_NOMINAL, f_numax=g["f_numax"], f_dnu=g["f_dnu"])
    r = np.median(m / g["mass"])
    assert abs(r - 1) < 0.001, r


def _synthetic_mixed_spectrum(dnu=8.0, numax=90.0, dpi1=80.0, q=0.15, eps=0.15, df=0.008, seed=3):
    """Asymptotic dipole mixed modes (Mosser et al. 2015) + radial/quadrupole modes, as Lorentzians."""
    rng = np.random.default_rng(seed)
    up = D.universal_pattern(dnu)
    nu = np.arange(numax - 5 * dnu, numax + 5 * dnu, df)
    snr = np.ones_like(nu)
    # mixed modes: solve tan(pi (nu - nu_p)/dnu) = q tan(pi / (dpi1 nu))  on a fine grid
    fine = np.arange(nu[0], nu[-1], df / 20)
    n_cont = fine / dnu - eps
    theta_p = np.pi * (n_cont - 0.5 + up["d01"] - 0.5 * up["alpha"] * (n_cont - numax / dnu + eps) ** 2)
    theta_g = np.pi / (dpi1 * fine * 1e-6)
    g = np.tan(theta_p) - q * np.tan(theta_g)
    roots = fine[:-1][(np.sign(g[:-1]) != np.sign(g[1:])) & (np.abs(g[:-1]) < 5) & (np.abs(g[1:]) < 5)]
    z = P.zeta_function(roots, dnu, numax, eps, dpi1, q)
    env = np.exp(-0.5 * ((roots - numax) / (2.0 * dnu)) ** 2)
    for f, zz, e in zip(roots, z, env):
        # p-dominated mixed modes (small zeta) are broad and tall like radial
        # modes; gravity-dominated ones (zeta -> 1) are narrow and weak
        h = 60 * e * np.sqrt(1 - zz) + 2
        w = 0.05 * (1 - zz) + 0.003
        snr += h / (1 + ((nu - f) / w) ** 2)
    for n in range(int(nu[0] / dnu) - 1, int(nu[-1] / dnu) + 2):
        for d0l, vis in ((0.0, 1.0), (up["d02"], 0.6)):
            f = (n + eps - d0l) * dnu
            snr += 80 * vis * np.exp(-0.5 * ((f - numax) / (2.0 * dnu)) ** 2) / (1 + ((nu - f) / 0.05) ** 2)
    return nu, snr * rng.exponential(1.0, nu.size)


def test_stretched_period_recovers_dpi1():
    nu, snr = _synthetic_mixed_spectrum()
    res = P.measure_dpi1(nu, snr, 8.0, 90.0, 0.15, 12.0)
    assert res["hypothesis"] == "RGB"
    assert abs(res["dpi1"] - 80.0) < 1.0, res["dpi1"]


def test_stretched_period_helium_burning():
    """Clump-like (300 s) and secondary-clump-like (250 s) cores are recovered and classified."""
    for dpi1, q, dnu, numax, seed in ((300.0, 0.3, 4.2, 35.0, 7), (250.0, 0.25, 6.0, 60.0, 8)):
        eps = D.universal_pattern(dnu)["eps"] % 1.0
        nu, snr = _synthetic_mixed_spectrum(dnu=dnu, numax=numax, dpi1=dpi1, q=q, eps=eps, seed=seed)
        res = P.measure_dpi1(nu, snr, dnu, numax, eps, 0.2 * numax / 2.35)
        assert res["hypothesis"] == "helium burning" and res["detected"], res
        assert abs(res["dpi1"] - dpi1) < 2.0, (dpi1, res["dpi1"])


def test_dnu_autocorrelation_and_refinement():
    """On a clump-like spectrum (dense, strongly coupled dipole modes) the plain
    autocorrelation is pulled by a few percent, as for real red-clump stars;
    the universal-pattern refinement removes the bias."""
    nu, snr = _synthetic_mixed_spectrum(dpi1=300.0, q=0.3)
    psd = snr.copy()
    res = D.measure_dnu(nu, psd, np.ones_like(psd), 90.0, 12.0, n_perturb=10)
    assert abs(res["dnu"] / 8.0 - 1) < 0.03, res["dnu"]
    ref = D.refine_dnu(nu, psd, np.ones_like(psd), 90.0, 12.0, res["dnu"], n_perturb=5)
    assert abs(ref["dnu"] / 8.0 - 1) < 0.003, ref["dnu"]
    # a red-giant-branch-like spectrum: the autocorrelation alone is good to 1%
    nu, snr = _synthetic_mixed_spectrum(dpi1=75.0, q=0.12, seed=5)
    res = D.measure_dnu(nu, snr, np.ones_like(snr), 90.0, 12.0, n_perturb=10)
    assert abs(res["dnu"] / 8.0 - 1) < 0.01, res["dnu"]


def test_note_names():
    assert SO.note_name(440.0).startswith("A4 +0")
    assert SO.note_name(261.63).startswith("C4")


if __name__ == "__main__":
    import sys

    fails = 0
    for name, f in list(globals().items()):
        if name.startswith("test_") and callable(f):
            try:
                f()
                print(f"ok    {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL  {name}: {e}")
    sys.exit(1 if fails else 0)
