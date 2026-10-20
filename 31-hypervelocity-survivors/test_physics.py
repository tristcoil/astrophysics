"""Checks of the physics and the plumbing, against independent calculations.

    python -m pytest test_physics.py -q      (or: python test_physics.py)
"""

from __future__ import annotations

import warnings

import numpy as np

warnings.filterwarnings("ignore")

import binary
import kinematics as K
import orbits
import population as P
from config import R0_KPC, V_SUN_KMS, Z_SUN_KPC


def test_fast_transform_matches_astropy():
    rng = np.random.default_rng(1)
    n = 500
    for ra, dec in ((249.38, -74.34), (0.0, 0.0), (300.0, 60.0)):
        d = rng.uniform(0.05, 20, n)
        pa, pd, rv = rng.normal(0, 200, n), rng.normal(0, 200, n), rng.normal(0, 1000, n)
        x1, v1 = K.to_galactocentric(ra, dec, d, pa, pd, rv)
        x2, v2 = K.to_galactocentric_astropy(ra, dec, d, pa, pd, rv)
        assert np.allclose(x1, x2, atol=1e-10) and np.allclose(v1, v2, atol=1e-8)


def test_sun_in_frame():
    # a star at the Sun's place moving with the Sun
    x, v = K.to_galactocentric(10.0, 20.0, np.array([1e-9]), np.array([0.0]), np.array([0.0]),
                               np.array([0.0]))
    assert np.allclose(x[:, 0], [-np.sqrt(R0_KPC ** 2 - Z_SUN_KPC ** 2), 0, Z_SUN_KPC], atol=1e-6)
    assert np.allclose(v[:, 0], V_SUN_KMS, atol=1e-6)


def test_roche_lobe_known_values():
    # Eggleton (1983): R_L/a = 0.3789 for equal masses, 0.2054 at q = 0.1 (wrt the lighter star)
    assert abs(binary.roche_lobe_fraction(1.0) - 0.37892) < 1e-4
    assert abs(binary.roche_lobe_fraction(0.1) - 0.20644) < 2e-3


def test_contact_orbit_kepler():
    a, v, p = binary.contact_orbit(0.6, 1.0, 0.0125)
    # Kepler's third law and the donor's share of the relative velocity
    v_rel = 2 * np.pi * a / p
    assert abs(v - v_rel * 1.0 / 1.6) / v < 1e-9


def test_shen2018_runaway_speeds():
    # Shen et al. 2018, sec. 2.4: 0.6 on 1.0 Msun -> 1,800 km/s; 1.0 + 1.0 -> 2,200 km/s
    assert abs(binary.ejection_speed(0.6, 1.0, "eggleton", with_kick=False) - 1800) < 60
    assert abs(binary.ejection_speed(1.0, 1.0, "eggleton", with_kick=False) - 2200) < 90


def test_white_dwarf_radii_sane():
    # Sirius B-like: ~0.0081 Rsun at 1.02 Msun (hot); our cold/10 kK radii bracket ~0.008
    assert 0.0070 < binary.eggleton_radius(1.0) < 0.0085
    assert 0.0075 < binary.bedard_radius(1.0) < 0.0090
    # heavier = smaller
    m = np.linspace(0.3, 1.2, 10)
    assert np.all(np.diff(binary.bedard_radius(m)) < 0)


def test_posterior_recovers_precise_distance():
    rng = np.random.default_rng(3)
    row = dict(parallax=2.0, parallax_error=0.02, pmra=100.0, pmdec=-50.0, pmra_error=0.02,
               pmdec_error=0.02, parallax_pmra_corr=0.0, parallax_pmdec_corr=0.0, pmra_pmdec_corr=0.0,
               ra=150.0, dec=30.0, phot_g_mean_mag=15.0, bp_rp=0.3)
    s = K.posterior(row, 50.0, 1.0, 0.0, rng, 3000, zero_point=False)
    assert abs(np.median(s["d"]) - 0.5) < 0.01
    assert abs(np.median(s["v_tan"]) - K.K_KMS * np.hypot(100, 50) * 0.5) < 2


def test_posterior_uses_parallax_pm_correlation():
    """With a strong parallax-pmra correlation, the conditional draws must
    reproduce the brute-force joint posterior (importance sampling)."""
    rng = np.random.default_rng(4)
    row = dict(parallax=0.5, parallax_error=0.2, pmra=60.0, pmdec=10.0, pmra_error=0.3,
               pmdec_error=0.3, parallax_pmra_corr=0.8, parallax_pmdec_corr=0.0, pmra_pmdec_corr=0.0,
               ra=150.0, dec=30.0, phot_g_mean_mag=18.0, bp_rp=0.0)
    s = K.posterior(row, 0.0, 10.0, 0.0, rng, 20000, zero_point=False, v_max=1e9)
    # brute force: draw (d from prior on grid, pm from wide box), weight by the joint likelihood
    n = 400000
    d = np.exp(rng.uniform(np.log(0.3), np.log(40), n))           # flat in log d == 1/d prior
    pm = rng.uniform(55, 65, n)
    cov = np.array([[0.04, 0.8 * 0.2 * 0.3], [0.8 * 0.2 * 0.3, 0.09]])
    r = np.vstack([0.5 - 1 / d, 60.0 - pm])
    w = np.exp(-0.5 * np.einsum("in,ij,jn->n", r, np.linalg.inv(cov), r))
    ref_mean = np.sum(w * pm) / w.sum()
    assert abs(np.mean(s["pmra"]) - ref_mean) < 0.03


def test_launch_speed_zero_for_disc_star():
    # a star at the Sun's position moving with the local circular velocity
    x = np.array([[-R0_KPC], [0.0], [0.0]])
    v = K.rotation_velocity(x, orbits.vcirc)
    # plus 1 km/s upward, 1 milliparsec up: it crossed the plane ~1,000 yr ago
    # (1 km/s is ~1 pc per Myr; much longer and the straight line would drift
    # off the curved disc orbit, the straight-line method's known limitation)
    t, xc, vej = K.launch(x + np.array([[0], [0], [1e-6]]), v + np.array([[0], [0], [1.0]]),
                          orbits.vcirc)
    assert abs(t[0] - 1e-6 / K.KPC_PER_KMS_MYR) < 1e-12
    assert abs(vej[0] - 1.0) < 0.01


def test_rotation_direction_and_energy():
    # a circular orbit at the Sun keeps its radius and moves toward +y first
    x = np.array([[-R0_KPC], [0.0], [0.0]])
    v = K.rotation_velocity(x, orbits.vcirc)
    assert v[1, 0] > 200
    pos, vel = orbits.integrate(x, v, np.linspace(0, 50, 501))
    R = np.hypot(pos[0, 0], pos[1, 0])
    assert np.ptp(R) < 0.02
    assert pos[1, 0, 1] > 0
    # energy conservation for a fast, unbound star
    v2 = np.array([[300.0], [1500.0], [800.0]])
    pos, vel = orbits.integrate(x, v2, np.linspace(0, -5, 501))
    E = 0.5 * np.sum(vel[:, 0] ** 2, axis=0) + orbits.phi(np.hypot(pos[0, 0], pos[1, 0]), pos[2, 0])
    assert np.ptp(E) / abs(E[0]) < 1e-6


def test_potential_grid_matches_direct():
    rng = np.random.default_rng(7)
    R, z = rng.uniform(0.3, 40, 30), rng.uniform(-20, 20, 30)
    assert np.max(np.abs(orbits.phi_fast(R, z) / orbits.phi(R, z) - 1)) < 1e-3


def test_escape_speed_reasonable():
    v = orbits.escape_speed(np.array([[-R0_KPC], [0.0], [0.0]]))[0]
    assert 550 < v < 650


def test_poisson_interval():
    lo, hi = P.poisson_interval(1)
    assert abs(lo - 0.1727) < 1e-3 and abs(hi - 3.2995) < 1e-3
    lo, hi = P.poisson_interval(0)
    assert lo == 0 and abs(hi - 1.8410) < 1e-3


def test_dust_layer():
    sun = np.array([-R0_KPC, 0.0, 0.0])
    # 1 kpc straight up through an exponential layer of scale height h: A_V = k h (1 - e^(-1/h))
    from config import DUST_AV_PER_KPC, DUST_SCALE_HEIGHT_KPC as h
    av = P.dust_av(sun, np.array([[-R0_KPC], [0.0], [1.0]]), n_steps=4000)[0]
    assert abs(av - DUST_AV_PER_KPC * h * (1 - np.exp(-1 / h))) < 1e-3
    # 1 kpc in the midplane: k x 1 kpc
    av = P.dust_av(sun, np.array([[-R0_KPC + 1], [0.0], [0.0]]))[0]
    assert abs(av - DUST_AV_PER_KPC) < 1e-9


if __name__ == "__main__":
    import sys
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"ok    {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL  {name} {e}")
    sys.exit(1 if fails else 0)
