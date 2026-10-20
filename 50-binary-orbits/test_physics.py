"""Checks of the physics and the data plumbing.

    python test_physics.py          # (or python -m pytest test_physics.py, if pytest is installed)
"""

from __future__ import annotations

import numpy as np
import pandas as pd

try:
    import pytest
except ImportError:      # the shared venv has no pytest: a minimal stand-in
    class _Skip(Exception):
        pass

    class pytest:  # noqa: N801
        @staticmethod
        def approx(x, rel=None, abs=None):
            class _A:
                def __eq__(self, y):
                    tol = max(rel * builtins_abs(x) if rel else 0.0, abs if abs else 0.0) or 1e-12 * builtins_abs(x)
                    return builtins_abs(y - x) <= tol

                def __repr__(self):
                    return f"approx({x})"
            return _A()

        @staticmethod
        def skip(msg):
            raise _Skip(msg)

    import builtins
    builtins_abs = builtins.abs

import config as C
import literature as L
import orb6
import orbits


def test_kepler_solver():
    rng = np.random.default_rng(1)
    M = rng.uniform(-10, 10, 5000)
    for e in (0.0, 0.3, 0.59, 0.9, 0.99):
        E = orbits.solve_kepler(M, e)
        assert np.max(np.abs(E - e * np.sin(E) - np.mod(M, 2 * np.pi))) < 1e-10


def test_thiele_innes_round_trip():
    rng = np.random.default_rng(2)
    for _ in range(200):
        a, i, O, w = rng.uniform(0.1, 10), rng.uniform(1, 179), rng.uniform(0, 180), rng.uniform(0, 360)
        A, B, F, G = orbits.thiele_innes(a, i, O, w)
        a2, i2, O2, w2 = orbits.campbell_from_thiele_innes(A, B, F, G)
        assert a2 == pytest.approx(a, rel=1e-10)
        assert i2 == pytest.approx(i, abs=1e-8)
        assert O2 == pytest.approx(O, abs=1e-8)
        assert orbits.wrap_angle_deg(w2 - w) == pytest.approx(0, abs=1e-8)


def test_bond_elements_reproduce_hubble_2016():
    """Bond et al.'s Table 4 orbit predicts their last HST measure (2016.6361) to < 1 mas."""
    e = L.BOND17_ELEMENTS
    T = orbits.besselian_to_year(e["T_besselian"][0])
    P = e["P"][0] * orbits.TROPICAL_YEAR / C.JULIAN_YEAR      # their years are Besselian
    rho, pa = orbits.rho_theta(orbits.besselian_to_year(2016.6361), P, T, e["e"][0], e["a"][0], e["i"][0],
                               e["Omega"][0], e["omega"][0])
    assert abs(rho - 10.6983) < 0.001
    assert abs(orbits.wrap_angle_deg(pa - 74.645)) < 0.01


def test_kepler_third_law_units():
    # the Earth: 1 au at 1 pc subtends 1", one sidereal year (365.25636 d)
    assert orbits.total_mass(1.0, 1.0, 365.256363 / C.JULIAN_YEAR) == pytest.approx(1.0, rel=2e-5)
    # Bond et al. 2017: 3.081 Msun from their a, P and adopted parallax
    M = orbits.total_mass(7.4957, 0.3789, 50.1284 * orbits.TROPICAL_YEAR / C.JULIAN_YEAR)
    assert M == pytest.approx(3.081, abs=0.001)


def test_relative_velocity_matches_finite_difference():
    el = (50.13, 1994.57, 0.5914, 7.4957, 136.34, 45.4, 149.16)
    t = np.array([1991.25, 2016.0])
    h = 1e-4
    x1, y1 = orbits.relative_xy(t - h, *el)
    x2, y2 = orbits.relative_xy(t + h, *el)
    vx, vy = orbits.relative_velocity(t, *el)
    assert np.allclose(vx, (x2 - x1) / (2 * h), atol=1e-7)
    assert np.allclose(vy, (y2 - y1) / (2 * h), atol=1e-7)


def test_orb6_overflowing_columns():
    """Long periods in days and big errors spill out of their documented columns."""
    o = orb6.read_orb6()
    stf73 = orb6.find(o, "00550+2338", "STF  73AB")
    assert stf73["P_yr"] == pytest.approx(61183.0 / C.JULIAN_YEAR)      # not 1183 d
    aand = o[(o["wds"] == "00084+2905") & o["disc"].str.contains("MKT")].iloc[0]
    assert aand["a_arcsec"] == pytest.approx(0.02443)                  # not lost to the P error
    assert aand["a_err_arcsec"] == pytest.approx(0.00039)
    sirius = orb6.find(o, "06451-1643", "AGC   1AB")
    assert (sirius["P_yr"], sirius["a_arcsec"], sirius["e"], sirius["grade"]) == (50.1284, 7.4957, 0.59142, 2)


def test_chandrasekhar_limit():
    import whitedwarf as WD
    M, R = WD.chandrasekhar_star(60.0, 2.0)
    assert 1.44 < M < 1.46 and R < 0.002
    # the non-relativistic end follows R ~ M^(-1/3)
    (M1, R1), (M2, R2) = WD.chandrasekhar_star(0.02), WD.chandrasekhar_star(0.04)
    assert np.log(R2 / R1) / np.log(M2 / M1) == pytest.approx(-1 / 3, abs=0.01)


def test_gaia_simulation_unbiased_without_orbit():
    """With no orbital motion the five-parameter fit returns the true parallax exactly."""
    import gaia_sirius as GS
    path = C.CACHE_DIR / "sirius_b_transits_dr3.csv"
    if not path.exists():
        pytest.skip("transits not computed yet (run run.py)")
    tr = pd.read_csv(path)
    r = GS.simulate_fit(tr, 101.2866, -16.7209, (50.0, 1994.6, 0.59, 7.5, 136.3, 45.4, 149.2), 0.0, 379.0)
    assert abs(r["plx_bias"]) < 1e-6 and r["rms_res"] < 1e-6


def test_hipparcos_leak_vanishes_for_gg78_orbit():
    """If our orbit of A were GG78's, there would be nothing to leak into the proper motion."""
    import gaia_sirius as GS
    d = L.HIP1_DMSA_O_SIRIUS
    P = d["P_d"] / C.JULIAN_YEAR
    T = orbits.jd_to_year(d["T_jd2440000"] + 2440000.0)
    # relative orbit whose primary share equals GG78's photocentre orbit: a = a0 / f, omega_B = omega_A - 180
    f = 0.33
    theta = (P, T, d["e"], d["a0_mas"] / 1e3 / f, d["i"], d["Omega"], d["omega"] - 180.0)
    b = GS.hipparcos_barycentre(theta, d["a0_mas"] / 1e3)
    assert abs(b["leak_pmra"]) < 1e-6 and abs(b["leak_pmdec"]) < 1e-6


def test_block_marginalisation_matches_brute_force():
    """The Sherman-Morrison block likelihood equals a dense multivariate normal."""
    from scipy.stats import multivariate_normal
    import sirius
    rng = np.random.default_rng(3)
    n = 7
    r = rng.normal(0, 1, n)
    sig = rng.uniform(0.5, 2.0, n)
    tau = 0.8
    cov = np.diag(sig**2) + tau**2 * np.ones((n, n))
    expected = multivariate_normal(np.zeros(n), cov).logpdf(r)

    class Dummy(sirius.HierLikelihood):
        def __init__(self):
            self.gid = np.zeros(n, int)
            self.ngroups = 1
    got = Dummy()._block_loglike(r, sig, np.ones(n, bool), np.array([tau]))
    assert got == pytest.approx(expected, rel=1e-10)


def test_photometric_radius_recovers_model():
    """A Montreal model star observed at its own magnitude gives back its own radius."""
    import whitedwarf as WD
    t = WD.table_da()
    row = t[(t["Teff"] == 25000.0) & (t["logg"] == 8.5)].iloc[0]
    plx = 0.1                                    # 10 pc: apparent = absolute magnitude
    r = WD.photometric_radius(row["V"], 1e-6, plx, 1e-9, 25000.0, 1e-6, 8.5)
    assert np.median(r["R"]) == pytest.approx(row["R_model"], rel=1e-3)


def main() -> int:
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    failed = 0
    for name, f in tests:
        try:
            f()
            print(f"PASS  {name}")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL  {name}: {exc}")
        except Exception as exc:  # noqa: BLE001
            if type(exc).__name__ in ("_Skip", "Skipped"):
                print(f"SKIP  {name}: {exc}")
            else:
                failed += 1
                print(f"ERROR {name}: {type(exc).__name__}: {exc}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
