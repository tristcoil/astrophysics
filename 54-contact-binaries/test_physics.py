"""Checks of the physics and the methods against known answers.

    python test_physics.py          # or: python -m pytest test_physics.py -q

Each test names what it checks. They need the cached data only where noted.
"""

from __future__ import annotations

import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("NUMBA_NUM_THREADS", "4")

import numpy as np

import config as C


# ---------------------------------------------------------------------------
# Roche geometry
# ---------------------------------------------------------------------------


def test_lagrange_points_equal_masses():
    """q = 1: L1 midway; L2 and L3 symmetric at 1.69841 (known value); grad Omega = 0 there."""
    import roche as R
    lp = R.lagrange_points(1.0)
    assert abs(lp["L1"] - 0.5) < 1e-12
    assert abs(lp["L2"] - 1.69841) < 1e-5 and abs(lp["L3"] + 0.69841) < 1e-5
    for name in ("L1", "L2", "L3"):
        gx, gy, gz = R.grad_omega(lp[name], 0.0, 0.0, 1.0)
        assert abs(gx) < 1e-9
    # L4 is a stationary point too
    gx, gy, gz = R.grad_omega(0.5, math.sqrt(3) / 2, 0.0, 1.0)
    assert math.hypot(gx, gy) < 1e-9


def test_lobe_volumes_match_eggleton():
    """Our integrated Roche-lobe volumes agree with Eggleton's (1983) formula to 1%,
    the accuracy he quotes, over q = 0.05..1 for both stars."""
    import roche as R
    for q in (1.0, 0.5, 0.1, 0.05):
        for which, qq in (("primary", 1 / q), ("secondary", q)):
            r = R.lobe_volume_radius(q, which, 48, 48)
            assert abs(r / R.eggleton_radius(qq) - 1) < 0.01, (q, which, r)


def test_fillout_definition():
    """f = 0 at the inner critical potential, f = 1 at the outer one."""
    import roche as R
    lp = R.lagrange_points(0.3)
    assert abs(R.fillout_potential(0.3, 0.0) - lp["omega_L1"]) < 1e-12
    assert abs(R.fillout_potential(0.3, 1.0) - lp["omega_L2"]) < 1e-12
    assert abs(R.fillout_factor(0.3, R.fillout_potential(0.3, 0.4)) - 0.4) < 1e-12


def test_renderer_limb_darkened_sphere():
    """A tiny, nearly spherical star: the summed image equals pi R^2 I0 (1 - u/3)
    for linear limb darkening, to 0.1%."""
    import roche as R
    q, om = 1e-4, 20.0
    b = R.ContactBinary.__new__(R.ContactBinary)
    b.q, b.lp, b.om = q, R.lagrange_points(q), om
    b.f = R.fillout_factor(q, om)
    b.t1 = b.t2 = 5800.0
    b.lam_m, b.beta, b.ld, b.spot = 640e-9, 0.0, (0.6, 0.0), None
    b.g_ref1, b.g_ref2 = b._pole_gravity(0.0), b._pole_gravity(1.0)
    b.cx, b.rb = 0.0, 0.08
    img, hit, area = b.render(60, 0.1, npix=400)
    r = R._surface_radius(q, om, 0.0, 0.0, 1.0, 0.0, 0.5)
    planck = 1 / (math.exp(R._HC_K / (640e-9 * 5800)) - 1)
    assert abs(img.sum() * area / (math.pi * r * r * planck * (1 - 0.6 / 3)) - 1) < 1e-3


def test_equal_twins_have_equal_minima():
    """q = 1, equal temperatures: the two eclipses are identical (symmetry)."""
    import roche as R
    b = R.ContactBinary(1.0, 0.3, t1=5800)
    assert abs(b.flux(80, 0.0, npix=140) / b.flux(80, 0.5, npix=140) - 1) < 1e-6


def test_kic9832227_model_amplitude():
    """Molnar et al.'s (2017) parameters give Kepler's mean amplitude to 0.01 mag (needs the cache)."""
    res = _results()
    if res is None:
        return
    k = res["roche"]["kic_model"]
    assert abs(k["amp_model_mag"] - k["amp_kepler_mag"]) < 0.01


# ---------------------------------------------------------------------------
# Timing and period laws
# ---------------------------------------------------------------------------


def test_template_shift_recovery():
    """Inject a known 0.0123-cycle shift into a noisy synthetic contact-binary
    light curve; the template fit recovers it within its error."""
    import timing as T
    rng = np.random.default_rng(1)
    coef = np.zeros(1 + 2 * 6)
    coef[0], coef[3], coef[7] = 1.0, -0.08, -0.01  # cos 4 pi phi and cos 8 pi phi terms
    coef[1] = 0.01
    tpl = T.Template(coef)
    p, t0 = 0.4579515, 2455000.0
    t = np.sort(rng.uniform(t0, t0 + 60, 600))
    true = 0.0123
    flux = tpl((t - t0) / p - true) * 1.03 + 0.02 + rng.normal(0, 0.01, t.size)
    r = T.fit_shift(t, flux, np.full(t.size, 0.01), tpl, p, t0)
    assert abs(r["delta"] - true) < 3 * r["delta_err"] and r["delta_err"] < 0.002


def test_no_negative_scale_alias():
    """A pure second-harmonic light curve flipped and shifted by a quarter cycle
    looks identical; with the positive-scale rule the fit must not take that alias."""
    import timing as T
    coef = np.zeros(1 + 2 * 4)
    coef[0], coef[3] = 1.0, -0.1
    tpl = T.Template(coef)
    t = np.linspace(0, 10, 400)
    flux = tpl(t / 0.5 - 0.02)
    r = T.fit_shift(t, flux, np.full(t.size, 0.002), tpl, 0.5, 0.0)
    assert abs(r["delta"] - 0.02) < 1e-3


def test_exponential_law_cycle_count():
    """Integrated cycle count of P = P0 exp(c / (t - t0)) matches direct
    summation of periods far from t0."""
    import v1309sco as V
    p0, c, t0, tref = 1.4456, 15.29, 2455233.5, 2452500.0
    t = 2452500.0
    n = 0
    while t < 2454000.0:
        # midpoint step: the period halfway through the orbit
        t += V.exp_period(t + 0.5 * V.exp_period(t, p0, c, t0), p0, c, t0)
        n += 1
    assert abs(V.cycles(np.array([t]), p0, c, t0, tref)[0] - n) < 1e-3


def test_linearised_exponential_oc():
    """For tau << |t - t0| the exponential O - C equals a(t-T) + b + c ln(t0 - t);
    our profile fit relies on it."""
    import kic9832227 as K
    p0, tau, t0 = 0.4579615, 0.114, K.jd_of(2022.2)
    t = np.linspace(K.jd_of(1999.5), K.jd_of(2021.0), 50)
    exact = K.oc_exponential(t, p0, tau, t0, 0.0)
    A = np.column_stack([t - K.T_REF, np.ones_like(t), np.log((t0 - t) / (t0 - K.T_REF))])
    coef, *_ = np.linalg.lstsq(A, exact, rcond=None)
    assert np.max(np.abs(A @ coef - exact)) * 1440 < 0.05  # minutes


def test_nsvs_half_day():
    """The NSVS frame-table 'JD' column is MJD: the ISO date of frame 80850
    (1999-04-06T08:53:47.35 UT) converts to 51274.37069 only as JD - 2400000.5."""
    from astropy.time import Time
    t = Time("1999-04-06T08:53:47.35", scale="utc")
    assert abs(t.mjd - 51274.37069) < 1e-5
    assert abs((t.jd - 2400000.0) - 51274.37069 - 0.5) < 1e-5


def test_bjd_conversion_size():
    """HJD(UTC) -> BJD(TDB) for KIC 9832227 differs by the ~65 s TDB-UTC offset plus
    less than 4 s of helio/barycentre difference."""
    import timing as T
    hjd = np.array([2454300.0, 2455000.0])
    bjd = T.hjd_utc_to_bjd_tdb(hjd, C.KIC9832227["ra"], C.KIC9832227["dec"])
    d = (bjd - hjd) * 86400
    assert np.all((d > 60) & (d < 72)), d


# ---------------------------------------------------------------------------
# Recurrence statistics and catalogue conventions
# ---------------------------------------------------------------------------


def test_pei2026_probabilities():
    """Our normal-recurrence calculation reproduces Pei et al. (2026): 30.2% for
    the rest of 2026 and 56.9% for the twelve months after 2026 July 11."""
    import literature as L
    import tcrb as TC
    last = 2000.0 + (L.SCHAEFER2023["peak_jd"][1946] - C.JD_J2000) / C.YEAR_D
    jul11 = 2000.0 + (2461232.5 - C.JD_J2000) / C.YEAR_D
    iv = L.PEI2026["intervals_yr"]
    assert abs(TC.recurrence_hazard(iv, last, jul11, 2027.0 - jul11)["p"] - 0.302) < 0.003
    assert abs(TC.recurrence_hazard(iv, last, jul11, 1.0)["p"] - 0.569) < 0.003


def test_chen2018_table2_reproduced():
    """Chen et al.'s (2018) W1 relation from their own Table 1 by ordinary least
    squares: slope -5.83, intercept -0.05, scatter 0.17 (needs the cache)."""
    res = _results()
    if res is None:
        return
    r = res["pl"]["chen2018"]["reproduce_W1_mean"]
    assert abs(r["slope"] + 5.83) < 0.005 and abs(r["b0"] + 0.05) < 0.005 and abs(r["scatter"] - 0.17) < 0.005


def _results():
    import json
    if not C.RESULTS_PATH.exists():
        return None
    return json.loads(C.RESULTS_PATH.read_text())


if __name__ == "__main__":
    import sys
    import time
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for name, fn in tests:
        t0 = time.time()
        try:
            fn()
            print(f"PASS {name} ({time.time() - t0:.1f} s)")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL {name}: {exc}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
