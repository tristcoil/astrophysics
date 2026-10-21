"""Sanity checks of the physics and the published formulas we transcribed.

    python -m pytest test_physics.py -q      (or: python test_physics.py)
"""

from __future__ import annotations

import numpy as np

import binaries as B
import literature as L
import stellar


def test_kepler_earth_and_jupiter():
    # The Earth: 1 AU, 1 year around 1 Msun.
    assert abs(B.semi_major_axis_au(365.25, 1.0) - 1.0) < 1e-12
    # Jupiter: 5.20 AU, 11.86 yr.
    assert abs(B.semi_major_axis_au(11.862 * 365.25, 1.0 + 9.55e-4) - 5.204) < 0.005
    assert abs(B.period_days_from_a(1.0, 1.0) - 365.25) < 1e-9


def test_k1_sun_from_jupiter():
    # Jupiter makes the Sun wobble at 12.5 m/s (edge-on).
    k = B.k1_kms(11.862 * 365.25, 1.0, 9.55e-4, ecc=0.048)
    assert abs(k * 1000 - 12.5) < 0.2


def test_mass_function_round_trip():
    p, m1, m2, e = 100.0, 1.1, 0.4, 0.3
    k = B.k1_kms(p, m1, m2, e)
    fm = B.mass_function(p, k, e)
    assert abs(fm - m2 ** 3 / (m1 + m2) ** 2) < 1e-9


def test_twin_photocentre_does_not_move():
    assert B.photocentre_factor(1.0, 1.0) == 0.0
    # a dark companion: the light centre follows the primary, a0 = a q/(1+q)
    assert abs(B.photocentre_factor(0.5, 0.0) - 1 / 3) < 1e-12
    # and the blended radial velocity of twins does not vary either
    assert B.light_weighted_k(30.0, 1.0, 1.0) == 0.0


def test_e_max_envelope():
    assert B.e_max_envelope(2.0) == 0.0
    assert abs(B.e_max_envelope(16.0) - (1 - 8 ** (-2 / 3))) < 1e-12


def test_kepler_solver():
    rng = np.random.default_rng(1)
    m = rng.uniform(0, 2 * np.pi, 1000)
    e = rng.uniform(0, 0.97, 1000)
    big_e = B.solve_kepler(m, e)
    assert np.max(np.abs(big_e - e * np.sin(big_e) - m)) < 1e-10


def test_projected_separation_statistics():
    # For circular orbits seen from random directions, <s/a> = pi/4.
    rng = np.random.default_rng(2)
    s = B.projected_separation_factor(rng, np.zeros(400_000))
    assert abs(s.mean() - np.pi / 4) < 0.003


def test_radial_velocity_curve_amplitude():
    t = np.linspace(0, 10, 20001)
    v = B.radial_velocity(t, 10.0, 0.0, 0.5, 0.7, 20.0)
    assert abs((v.max() - v.min()) / 2 - 20.0) < 0.01


def test_meibom_mathieu_pcirc():
    # e(Pcirc) = 0.01 by construction
    pc = B.meibom_mathieu_pcirc(8.0)
    assert abs(B.meibom_mathieu_e(pc, 8.0) - 0.01) < 1e-12


def test_moe_di_stefano_integrals():
    # Sec. 9.4: f_mult;q>0.3 = 0.36 (1 Msun), 0.63 (7 Msun), 1.02 (28 Msun)
    x = np.linspace(0.2, 7.99999, 100001)
    for m, ref in ((1, 0.36), (7, 0.63), (28, 1.02)):
        assert abs(np.trapezoid(L.mds_flogp_q03(x, m), x) - ref) < 0.01
    # Table 13 twin fractions for solar primaries at log P = 1, 3, 5
    assert np.allclose(L.mds_ftwin(np.array([1.0, 3.0, 5.0]), 1.0), [0.30, 0.20, 0.10])


def test_mamajek_sun():
    # The Sun (G2V): ~1 Msun, M_G ~ 4.6, BP-RP ~ 0.82
    assert abs(stellar.mass_from_bprp(0.823) - 1.0) < 0.01
    assert abs(stellar.mg_from_mass(1.0) - 4.635) < 0.01


def test_rv_error_inversion():
    import spectro
    # Katz et al. 2023 Eq. 1 forwards, then inverted
    n, sigma = 20, 0.4
    eps = np.sqrt((np.sqrt(np.pi / (2 * n)) * sigma) ** 2 + spectro.ERROR_FLOOR ** 2)
    assert abs(spectro.epoch_sigma_from_error(eps, n) - sigma) < 1e-12


def test_f2_statistic_is_standard_normal_for_constant_stars():
    import spectro
    rng = np.random.default_rng(3)
    nu = 19
    f2 = spectro.f2_statistic(rng.chisquare(nu, 200_000), nu)
    assert abs(f2.mean()) < 0.02 and abs(f2.std() - 1) < 0.02


def test_phase_gap():
    import spectro
    t = np.arange(10) * 1.0
    assert abs(spectro.max_phase_gap(t, 10.0) - 0.1) < 1e-12
    assert spectro.max_phase_gap(np.array([0.0, 0.1]), 10.0) > 0.98


def test_sensitivity_curve():
    import wide
    assert abs(wide.sensitivity(2.0, 2.0, 10.0) - 0.5) < 1e-12
    assert wide.sensitivity(4.0, 2.0, 10.0) > 0.99 and wide.sensitivity(1.0, 2.0, 10.0) < 0.01


def test_meibom_mathieu_fit_recovers_circularization_period():
    import pe
    rng = np.random.default_rng(4)
    p = 10 ** rng.uniform(0, 2.5, 3000)
    e_true = B.meibom_mathieu_e(p, 8.0)
    e = np.clip(e_true * rng.uniform(0, 2, p.size), 0, 0.95)  # scatter about the curve, as real eccentricities
    r = pe.fit_meibom_mathieu(p, e, p_grid=np.arange(2, 20, 0.02))
    assert abs(r["p_circ"] - B.meibom_mathieu_pcirc(8.0)) < 0.6


def test_light_weighted_velocity_of_blended_twin_is_zero():
    import spectro
    rng = np.random.default_rng(5)
    t = np.sort(rng.uniform(0, 1000, 25))
    # twins with blended lines (K1 + K2 below the resolution): nothing to detect
    res = spectro.simulate(t, 0.3, 300.0, 0.0, 0.2, 1.0, 15.0, 1.0, 1.0, rng)
    assert res in ("", "variable")


if __name__ == "__main__":
    import sys
    failures = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok   {name}")
            except AssertionError as exc:
                failures += 1
                print(f"FAIL {name}: {exc}")
    sys.exit(1 if failures else 0)
