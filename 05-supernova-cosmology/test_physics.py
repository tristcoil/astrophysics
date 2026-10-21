"""Self-checks for the cosmology module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas, not the data.

Every result in this project is a comparison between measured distances and the
distances these functions predict. A slip here -- a missing (1 + z), a sin where
a sinh belongs -- would produce a perfectly plausible-looking wrong answer for
the contents of the universe, so each function is checked two independent ways:
against textbook closed forms, and against astropy's implementation.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy.integrate import quad

import cosmology as cosmo

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "",
          *, absolute: bool = False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if absolute or expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    kind = "abs" if absolute or expected == 0 else "rel"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<48} {got:>12.6g} vs {expected:>12.6g} "
          f"{unit:<6} ({kind} err {err:.1e}, tol {tol:.0e})")


def assert_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    z_test = np.array([0.01, 0.1, 0.5, 1.0, 2.26])

    print("Units")
    # c/H0 for H0 = 70 km/s/Mpc is the textbook 4283 Mpc, 1/H0 is 13.97 Gyr.
    check("c/H0 at H0 = 70", cosmo.hubble_distance_mpc(70.0), 4282.7494, 1e-6, "Mpc")
    check("1/H0 at H0 = 70", cosmo.hubble_time_gyr(70.0), 13.968, 1e-3, "Gyr")
    check("mu at 10 pc", float(cosmo.distance_modulus(1e-5)), 0.0, 1e-12, "mag")
    check("mu at 1 Mpc", float(cosmo.distance_modulus(1.0)), 25.0, 1e-12, "mag")

    print("\nNumerical distances vs closed forms (the three textbook universes)")
    nodes, idx = cosmo.integration_nodes(z_test)
    for label, (om, ol), exact in [
        ("empty / Milne (0, 0)", (0.0, 0.0), cosmo.dl_milne),
        ("Einstein-de Sitter (1, 0)", (1.0, 0.0), cosmo.dl_einstein_de_sitter),
        ("de Sitter (0, 1)", (0.0, 1.0), cosmo.dl_de_sitter),
    ]:
        got = cosmo.luminosity_distance_shape(z_test, z_test, om, ol, nodes=nodes, index=idx)[0]
        worst = np.max(np.abs(got / exact(z_test) - 1.0))
        check(f"{label}: worst d_L over z", worst, 0.0, 1e-6, absolute=True)

    print("\nAgainst astropy (radiation off), open, flat, closed and w != -1")
    from astropy.cosmology import LambdaCDM, wCDM

    for om, ol, w in [(0.3, 0.7, -1.0), (0.3, 0.0, -1.0), (0.1, 0.4, -1.0),
                      (1.2, 0.5, -1.0), (0.3, 0.7, -0.8), (0.25, 0.75, -1.4)]:
        if w == -1.0:
            ref = LambdaCDM(H0=70, Om0=om, Ode0=ol, Tcmb0=0)
        else:
            ref = wCDM(H0=70, Om0=om, Ode0=ol, w0=w, Tcmb0=0)
        want = ref.luminosity_distance(z_test).value
        got = cosmo.luminosity_distance_shape(z_test, z_test, om, ol, w)[0] \
            * cosmo.hubble_distance_mpc(70.0)
        dmu = np.max(np.abs(5 * np.log10(got / want)))
        check(f"Om={om}, OL={ol}, w={w}: worst |d mu|", dmu, 0.0, 2e-5, "mag", absolute=True)

    print("\nIntegration accuracy vs adaptive quadrature (closed universe, z = 2.26)")
    om, ol = 0.6, 0.9
    exact_dc, _ = quad(lambda z: 1 / np.sqrt(cosmo.e_squared(z, om, ol)), 0, 2.26,
                       epsabs=1e-13, epsrel=1e-13)
    nodes, idx = cosmo.integration_nodes(np.array([2.26]))
    got_dc = cosmo.comoving_distance_nodes(nodes, om, ol)[0, idx[0]]
    check("D_C(2.26) trapezoid vs quad", got_dc, exact_dc, 1e-6)

    print("\nLeaving out radiation is harmless at supernova redshifts")
    from astropy.cosmology import FlatLambdaCDM

    zz = np.linspace(0.01, 2.3, 200)
    with_rad = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.7255, Neff=3.046, m_nu=0.0)
    got = cosmo.distance_modulus_curve(zz, 0.3, 0.7, H0=70.0)
    worst = np.max(np.abs(got - with_rad.distmod(zz).value))
    check("worst |d mu| from ignoring radiation", worst, 0.0, 1e-3, "mag", absolute=True)

    print("\nDeceleration")
    check("q0 Einstein-de Sitter", float(cosmo.q0(1.0, 0.0)), 0.5, 1e-12)
    check("q0 empty universe", float(cosmo.q0(0.0, 0.0)), 0.0, 1e-12)
    check("q0 (0.3, 0.7)", float(cosmo.q0(0.3, 0.7)), -0.55, 1e-12)
    # A w = -1/3 fluid exerts no net force: it behaves like curvature.
    check("q0 w = -1/3 exerts no force", float(cosmo.q0(0.0, 1.0, -1.0 / 3.0)), 0.0, 1e-12)
    zt = float(cosmo.transition_redshift(0.3, 0.7))
    check("transition z for (0.3, 0.7)", zt, (2 * 0.7 / 0.3) ** (1 / 3) - 1, 1e-12)
    check("q(z_t) = 0 there", float(cosmo.deceleration_parameter(zt, 0.3, 0.7)), 0.0,
          1e-12, absolute=True)

    print("\nAges")
    # Einstein-de Sitter: t0 = (2/3)/H0. Empty universe: t0 = 1/H0 exactly.
    check("age, Einstein-de Sitter", cosmo.age_gyr(1.0, 0.0, H0=70), 2 / 3 * 13.968, 1e-4, "Gyr")
    check("age, empty universe", cosmo.age_gyr(0.0, 0.0, H0=70), cosmo.hubble_time_gyr(70), 1e-6, "Gyr")
    # Flat LambdaCDM has a closed form: t0 = 2/(3 H0 sqrt(OL)) asinh(sqrt(OL/Om)).
    om, ol = 0.3, 0.7
    exact = 2 / (3 * np.sqrt(ol)) * np.arcsinh(np.sqrt(ol / om)) * cosmo.hubble_time_gyr(70)
    check("age, flat (0.3, 0.7) vs closed form", cosmo.age_gyr(om, ol, H0=70), exact, 1e-5, "Gyr")
    ref = LambdaCDM(H0=70, Om0=0.4, Ode0=0.8, Tcmb0=0)
    check("age, closed (0.4, 0.8) vs astropy", cosmo.age_gyr(0.4, 0.8, H0=70),
          ref.age(0).value, 1e-5, "Gyr")
    check("lookback to z=1, (0.3, 0.7) vs astropy",
          cosmo.lookback_time_gyr(1.0, 0.3, 0.7, H0=70),
          LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0).lookback_time(1.0).value, 1e-5, "Gyr")

    print("\nThe 'no Big Bang' boundary: closed form vs direct analysis")
    om_grid = np.linspace(0.02, 1.5, 75)
    edge = cosmo.no_big_bang_boundary(om_grid)
    above = cosmo.bounces(om_grid, edge + 1e-3)
    below = cosmo.bounces(om_grid, edge - 1e-3)
    assert_true(f"just above the line bounces at all {om_grid.size} test points", bool(above.all()))
    assert_true("just below the line has a Big Bang at all of them", bool((~below).all()))
    check("boundary as Om -> 0 tends to 1", float(cosmo.no_big_bang_boundary(1e-6)), 1.0, 1e-3)
    check("boundary at Om = 1/2 is 2", float(cosmo.no_big_bang_boundary(0.5)), 2.0, 1e-9)
    # A bouncing universe cannot show us light from beyond its bounce: the
    # distance integral must refuse (NaN) rather than return a number.
    nodes, idx = cosmo.integration_nodes(np.array([2.0]))
    assert_true("distance to z=2 in a bouncing (0.05, 1.5) universe is NaN",
                bool(np.isnan(cosmo.comoving_distance_nodes(nodes, 0.05, 1.5)).all()))

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
