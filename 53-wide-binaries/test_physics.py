"""Checks of the project's building blocks against exact or published values.

    python -m pytest test_physics.py -q      (or: python test_physics.py)
"""

from __future__ import annotations

import numpy as np

import kinematics as KIN
import masses as ML
import qumond as Q


def test_galactic_field_matches_banik():
    """g_e = v_c^2/R0 = 1.785 a0 and g_Ne = 1.144 a0 (Banik et al. 2024, Sect. 3.1)."""
    gf = Q.galactic_field()
    assert abs(gf["ge_a0"] - 1.785) < 0.002
    assert abs(gf["gNe_a0"] - 1.144) < 0.001
    # the inversion is exact for the simple interpolating function
    assert abs(Q.nu_simple(gf["gNe_a0"]) * gf["gNe_a0"] - gf["ge_a0"]) < 1e-12


def test_qumond_limits():
    """Near the star: isolated MOND; far away: the analytic EFE-dominated
    field nu_e (1 + K_e sin^2(theta)/2), angle average nu_e (1 + K_e/3)."""
    gf = Q.galactic_field()
    f = Q.solve(gf["gNe_a0"], lmax=96, nr=4000)
    r = np.exp(f.lnr)
    i_near = np.argmin(abs(r - 0.03))
    assert abs(f.ratio_r[i_near, 180] - Q.nu_simple(1 / r[i_near] ** 2)) < 2e-3
    i_far = np.argmin(abs(r - 300))
    mu = f.mu
    pred = gf["nu_e"] * (1 + gf["K_e"] / 2 * (1 - mu**2))
    assert np.max(np.abs(f.ratio_r[i_far] - pred)) < 0.01
    assert abs(f.eta[i_far] - gf["eta_inf"]) < 2e-3
    # no net phantom force on the star itself
    assert abs(f.self_g) < 1e-5


def test_newtonian_library_orbit():
    """RK4 in the field tables with alpha = 0 conserves energy and reaches
    the exact peri/apocentre."""
    import orbits as OB
    tabs = OB.field_tables()
    out = np.empty((OB.N_STORE, 6))
    rmin, rmax, revs, steps = OB._orbit(0.5, 0.9, 0.5, 0.0, *tabs, 20, OB.DT_FAC, OB.N_STORE, 1, out)
    E = 0.5 * np.sum(out[:, 3:] ** 2, 1) - 1 / np.linalg.norm(out[:, :3], axis=1)
    assert np.max(np.abs(E + 1.0)) < 1e-6
    assert abs(rmin - 0.05) < 1e-4 and abs(rmax - 0.95) < 1e-4


def test_mass_luminosity_inverse():
    m = np.array([0.15, 0.3, 0.6, 1.0, 1.5, 2.5])
    back = ML.banik_mass(ML.banik_MG(m), flame=False)
    assert np.max(np.abs(back / m - 1)) < 1e-6
    # Banik et al.'s FLAME calibration leaves masses near 0.75 Msun untouched
    assert abs(ML.flame_adjust(np.array([0.75]))[0] - 0.75) < 1e-12


def test_perspective_cancels_for_comoving_pair():
    """Two stars moving with one 3D velocity have v_sky ~ 0 once the systemic
    RV is put along each line of sight (and a spurious ~V_r theta without)."""
    ra1, dec1 = np.array([30.0]), np.array([10.0])
    ra2, dec2 = np.array([30.3]), np.array([10.2])           # ~0.36 deg apart
    d = np.array([100.0])
    n0 = KIN.unit(ra1, dec1)[0][0]
    V = 40.0 * n0 + np.array([20.0, -35.0, 15.0])              # km/s, shared; ~40 km/s radial
    pm = []
    for ra, dec in ((ra1, dec1), (ra2, dec2)):
        n, ea, ed = KIN.unit(ra, dec)
        vt = V - np.sum(V * n, axis=1, keepdims=True) * n
        pm.append((np.sum(vt * ea, 1) / (KIN.K_PM * d / 1000), np.sum(vt * ed, 1) / (KIN.K_PM * d / 1000)))
    n1 = KIN.unit(ra1, dec1)[0]
    n2 = KIN.unit(ra2, dec2)[0]
    nm = (n1 + n2) / np.linalg.norm(n1 + n2)
    rv_sys = np.sum(V * nm, 1)
    v_ok, _ = KIN.sky_velocity(ra1, dec1, pm[0][0], pm[0][1], ra2, dec2, pm[1][0], pm[1][1], d, d, rv_sys)
    v_bad, _ = KIN.sky_velocity(ra1, dec1, pm[0][0], pm[0][1], ra2, dec2, pm[1][0], pm[1][1], d, d, rv_sys,
                                perspective=False)
    # residual is second order in the separation (~V theta^2); spurious term ~V_r theta
    assert np.linalg.norm(v_ok) < 2e-3, np.linalg.norm(v_ok)
    assert np.linalg.norm(v_bad) > 0.15, np.linalg.norm(v_bad)


def test_deprojection_prior_normalised():
    rng = np.random.default_rng(0)
    x = KIN.sample_deprojection(rng, 200000)
    assert x.min() >= 1.0
    # P(x) ~ x^-2.6 / sqrt(x^2 - 1): median ~1.13 (computed by quadrature)
    from scipy.integrate import quad
    f = lambda t: t**-2.6 / np.sqrt(t * t - 1)
    tot = quad(f, 1, np.inf)[0]
    med = np.median(x)
    assert abs(quad(f, 1, med)[0] / tot - 0.5) < 0.01


if __name__ == "__main__":
    import sys
    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except AssertionError as e:
                fails += 1
                print(f"FAIL {name}: {e}")
    sys.exit(1 if fails else 0)
