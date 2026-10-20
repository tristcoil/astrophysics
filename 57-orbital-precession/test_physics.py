"""Physics checks: every formula and convention the project relies on, tested
against an independent calculation or a published number.

    python test_physics.py
"""

from __future__ import annotations

import sys
import traceback

import numpy as np

import ephemeris as eph
import orbit_physics as op
from config import ARCSEC_PER_RAD, J2000_JD, RUNS_DIR


def test_elements_match_horizons():
    """Our osculating elements from DE440 = Horizons' printed elements (J2000)."""
    r, v = eph.heliocentric_state("mercury", J2000_JD)
    rot = eph.rot_icrf_to_ecliptic()
    el = op.elements(rot @ r, rot @ v, eph.gm("sun") + eph.gm("mercury"))
    # Horizons API, Mercury (199) wrt Sun body centre, ecliptic J2000, JD 2451545.0 (DE441).
    assert abs(el["e"][0] - 2.056302933705567e-01) < 1e-13
    assert abs(el["a"][0] - 3.870982121841597e-01) < 1e-13
    assert abs(np.degrees(el["inc"][0]) - 7.005014303275355) < 1e-9
    assert abs(np.degrees(el["Omega"][0]) - 48.33053855197922) < 1e-9
    assert abs(np.degrees(el["omega"][0]) - 29.12428165737838) < 1e-9


def test_in_plane_angle_is_frame_independent():
    """psi from a precessing, tilting orbit is the same in any rotated frame."""
    t = np.linspace(0, 1, 500)
    Om, inc, w = 0.3 - 0.8 * t, 0.12 + 0.05 * t, 1.0 + 1.3 * t
    e_vec, h_vec = [], []
    for O, i, ww in zip(Om, inc, w):
        P = np.array([np.cos(O) * np.cos(ww) - np.sin(O) * np.sin(ww) * np.cos(i),
                      np.sin(O) * np.cos(ww) + np.cos(O) * np.sin(ww) * np.cos(i), np.sin(ww) * np.sin(i)])
        H = np.array([np.sin(O) * np.sin(i), -np.cos(O) * np.sin(i), np.cos(i)])
        e_vec.append(P)
        h_vec.append(H)
    e_vec, h_vec = np.array(e_vec).T, np.array(h_vec).T
    psi1 = op.in_plane_angle(e_vec, h_vec)
    R = eph.rot_icrf_to_ecliptic() @ np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1.0]])
    psi2 = op.in_plane_angle(R @ e_vec, R @ h_vec)
    assert np.abs(psi1 - psi2).max() < 1e-12
    # ... and equals the integral of domega + cos(i) dOmega.
    ref = np.concatenate([[0], np.cumsum(np.diff(w) + np.cos(0.5 * (inc[1:] + inc[:-1])) * np.diff(Om))])
    assert np.abs(psi1 - ref).max() < 1e-5


def test_schwarzschild_mercury():
    """6 pi GM / (c^2 a (1-e^2)) per orbit = 42.98"/cy for Mercury (textbook value)."""
    r, v = eph.heliocentric_state("mercury", J2000_JD)
    mu = eph.gm("sun") + eph.gm("mercury")
    el = op.elements(r, v, mu)
    rate = op.schwarzschild_rate(mu, el["a"][0], el["e"][0], eph.c_au_per_day())
    assert abs(rate - 42.98) < 0.005


def test_laplace_coefficient():
    """b_{3/2}^{(1)} against its series: 3 alpha (1 + 15/8 alpha^2 + 175/64 alpha^4 + ...).

    The series is cut after alpha^6; the next term is ~3e-11 at alpha = 0.05."""
    a = 0.05
    series = 3 * a * (1 + 15 / 8 * a**2 + 175 / 64 * a**4 + 3675 / 1024 * a**6)
    assert abs(op.laplace_coefficient(1.5, 1, a) - series) < 1e-9


def test_gauss_equals_laplace_lagrange_for_circular_orbits():
    """Gauss's averaging -> the Laplace-Lagrange free precession when e, i -> 0."""
    mu, gm_out = 1.0, 1e-3
    r_in, v_in = np.array([0.4 * 0.999, 0, 0]), np.array([0, np.sqrt(mu / 0.4 * 1.001 / 0.999), 0])
    r_out, v_out = np.array([1.0, 0, 0]), np.array([0, np.sqrt(mu), 0])
    g = op.gauss_secular_rate(mu, r_in, v_in, gm_out, r_out, v_out, n=720)["rate"]
    el = op.elements(r_in, v_in, mu)
    ll = op.laplace_lagrange_rate(mu, 0.0, el["a"][0], el["e"][0], 0.0, gm_out, 1.0, 0.0, 0.0)["total"]
    assert abs(g / ll - 1) < 2e-5


def test_sterne_small_eccentricity():
    """T2 - T1 -> P/2 + (2P/pi) e cos(omega) for small e."""
    from binaries import sterne_spacing
    P, e, w = 10.0, 1e-4, 0.7
    assert abs(sterne_spacing(P, e, w) - (P / 2 + 2 * P / np.pi * e * np.cos(w))) < 1e-10


def test_baroch_gr_rates():
    """Our 6 pi x per orbit reproduces Baroch et al.'s Table 9 GR column (8 of 9 systems)."""
    import binaries
    res = binaries.analyse()
    bad = [n for n, s in res["systems"].items()
           if abs(s["gr_ours"] - s["wdot_gr"]) > 2e-7 and n != "KW Hya"]
    assert not bad, bad


def test_pulsar_masses_from_published_route():
    """omega-dot + gamma give Weisberg & Huang's masses (1.438 / 1.390 Msun)."""
    from scipy.optimize import fsolve
    import pulsars as P
    from config import PUBLISHED
    p = PUBLISHED["b1913"]
    pb = p["pb_days"] * 86400
    target = np.radians(p["omdot_deg_yr"]) / (365.25 * 86400)

    def eqs(m):
        return [P.omdot_gr(pb, p["e"], m[0] + m[1]) / target - 1,
                P.gamma_gr(pb, p["e"], m[0], m[1]) / p["gamma_s"] - 1]
    mp, mc = fsolve(eqs, [1.4, 1.4], xtol=1e-13)
    assert abs(mp - 1.438) < 0.001 and abs(mc - 1.390) < 0.001


def test_double_pulsar_masses():
    """Shapiro shape + orbital decay give Kramer et al.'s masses (1.338176 / 1.248842)."""
    import pulsars as P
    from config import PUBLISHED
    p = PUBLISHED["j0737"]
    ma, mb = P.masses_s_pbdot(p["pb_days"] * 86400, p["e"], p["x_a_s"], p["s"], p["pbdot_gw"])
    assert abs(ma - 1.338176) < 2e-4 and abs(mb - 1.248842) < 1e-4


def test_s2_apocentre_shift():
    """GRAVITY 2020: GR moves S2's 2026 apocentre by 0.78 mas, mostly in RA."""
    import s2
    sig = s2.precession_signal()
    assert abs(sig["shift_at_apocentre_mas"] - 0.78) < 0.01
    assert abs(sig["shift_at_apocentre_ra_dec_mas"][0]) > abs(sig["shift_at_apocentre_ra_dec_mas"][1])
    assert abs(sig["precession_per_orbit_arcmin"] - 12.1) < 0.1


def test_two_body_runs():
    """Cached N-body runs: lone Mercury does not precess; GR alone gives 42.98"/cy."""
    if not (RUNS_DIR / "two_body_gr.npz").exists():
        raise SkipTest("runs not computed")
    import budget
    zero = budget.run_series("two_body_newton")
    assert abs(op.fit_secular(zero["t_cy"], zero["psi"])["rate"]) < 1e-6
    for name in ("two_body_gr", "two_body_gr_full", "two_body_gr_potential"):
        s = budget.run_series(name)
        assert abs(op.fit_secular(s["t_cy"], s["psi"])["rate"] - 42.9807) < 1e-3, name


def test_park_recipe_reproduces_messenger():
    """Park et al.'s own fit on our six-millennium run (+ asteroid estimate) = 575.3100 +- 0.0015."""
    if not (RUNS_DIR / "millennia_full.npz").exists():
        raise SkipTest("six-millennium run not computed")
    import park_check
    r = park_check.analyse()
    assert abs(r["plus_asteroid_estimate"] - 575.3100) < 0.0015
    assert abs(r["with_periodic_terms"]["Q"] - park_check.PARK_Q) < 2e-4


class SkipTest(Exception):
    pass


def main() -> int:
    tests = [(k, v) for k, v in globals().items() if k.startswith("test_") and callable(v)]
    failed, skipped = [], []
    for name, fn in tests:
        label = name[5:].replace("_", " ")
        try:
            fn()
            print(f"  [ok  ] {label}")
        except SkipTest:
            skipped.append(label)
            print(f"  [skip] {label}")
        except Exception:  # noqa: BLE001
            failed.append(label)
            print(f"  [FAIL] {label}")
            traceback.print_exc(limit=2)
    print()
    if failed:
        print(f"{len(failed)} CHECK(S) FAILED: {', '.join(failed)}")
        return 1
    print(f"All {len(tests) - len(skipped)} physics checks passed"
          + (f" ({len(skipped)} skipped)." if skipped else "."))
    return 0


if __name__ == "__main__":
    sys.exit(main())
