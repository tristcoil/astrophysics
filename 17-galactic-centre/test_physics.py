"""Self-checks for the orbit model, the likelihood and the dark-cluster arithmetic.

Run with ``python test_physics.py`` (no test framework needed). The checks that
use data -- grading our conventions against Do et al.'s own posterior chains --
need ``python fetch_data.py`` to have run once, and are skipped otherwise.

Each check guards against a mistake that would still produce plausible-looking
plots: a flipped sign in the sky projection, a light-travel delay applied the
wrong way round, a precession that turns the ellipse by the wrong amount.
"""

from __future__ import annotations

import sys
import traceback

import numpy as np

import black_hole_case as bh
import fitting as ft
import orbit_physics as op


def approx(got, expected, *, rel=None, abs=None):  # noqa: A002 -- mirrors pytest.approx
    tol = max(rel * np.abs(expected) if rel else 0.0, abs or 0.0)
    if tol == 0.0:
        tol = 1e-12 * max(1.0, np.abs(expected))
    return np.abs(got - expected) <= tol


def parametrize(values):
    def deco(fn):
        fn.cases = values
        return fn
    return deco

S2 = dict(gm6=4.1, period=16.05, t0=2018.38, e=0.886, inc_deg=133.8, omega_deg=66.1, Omega_deg=227.5)


def _el(**kw):
    p = dict(S2) | kw
    return (p["gm6"], p["period"], p["t0"], p["e"], p["inc_deg"], p["omega_deg"], p["Omega_deg"])


def _have_data():
    from s2data import DO2019_DIR
    return (DO2019_DIR / "S0-2_rv.csv").exists()


def needs_data(fn):
    fn.needs_data = True
    return fn


# --------------------------------------------------------------------------
# Kepler
# --------------------------------------------------------------------------


def test_kepler_solver_accuracy():
    m = np.linspace(-30, 30, 20001)
    for e in (0.0, 0.3, 0.886, 0.97):
        big_e, _ = op.solve_kepler(m, e)
        assert np.max(np.abs(big_e - e * np.sin(big_e) - m)) < 1e-11


def test_keplers_third_law_earth():
    # The Earth's sidereal year is 365.25636 days: a = 1 AU to within 1e-5
    # (the residual is the Earth-Moon mass the Sun's GM does not include).
    a = op.semi_major_axis_au(1e-6, 365.256_363 / 365.25)
    assert approx(a, 1.0, abs=1e-5)
    assert approx(op.mass_from_orbit(1.0, 365.256_363 / 365.25), 1.0, rel=3e-5)


def test_pericentre_apocentre_and_closure():
    t = np.array([S2["t0"], S2["t0"] - S2["period"] / 2, S2["t0"] - S2["period"]])
    x, y, z, *_, r = op.orbit_state(t, *_el())
    a = op.semi_major_axis_au(S2["gm6"], S2["period"])
    assert approx(r[0], a * (1 - S2["e"]), rel=1e-12)
    assert approx(r[1], a * (1 + S2["e"]), rel=1e-12)
    # A Newtonian orbit closes exactly after one period.
    assert np.allclose([x[0], y[0], z[0]], [x[2], y[2], z[2]], atol=1e-8)


@parametrize([0.0, 1.0])
def test_velocity_is_derivative_of_position(f_sp):
    t = np.linspace(2010.0, 2019.0, 3001)
    x, y, z, vx, vy, vz, _ = op.orbit_state(t, *_el(), f_sp=f_sp)
    dt = 1e-6
    x2, y2, z2, *_ = op.orbit_state(t + dt, *_el(), f_sp=f_sp)
    x1, y1, z1, *_ = op.orbit_state(t - dt, *_el(), f_sp=f_sp)
    for num, ana in (((x2 - x1), vx), ((y2 - y1), vy), ((z2 - z1), vz)):
        num_kms = num / (2 * dt) * op.KMS_PER_AUYR
        assert np.max(np.abs(num_kms - ana)) < 0.05  # km/s, out of ~8000 at pericentre


def test_vis_viva_energy_conserved():
    t = np.linspace(2000.0, 2020.0, 5001)
    x, y, z, vx, vy, vz, r = op.orbit_state(t, *_el())
    gm = S2["gm6"] * 1e6 * op.GM_SUN_KM3S2
    energy = 0.5 * (vx**2 + vy**2 + vz**2) - gm / (r * op.AU_KM)
    a_km = op.semi_major_axis_au(S2["gm6"], S2["period"]) * op.AU_KM
    assert np.allclose(energy, -gm / (2 * a_km), rtol=1e-10)


def test_sky_motion_is_clockwise_for_s2():
    # i > 90 deg: position angle decreases with time, as S2's does on the sky.
    t = np.linspace(2008.0, 2012.0, 50)
    x, y, *_ = op.orbit_state(t, *_el())
    pa = np.unwrap(np.arctan2(x, y))
    assert np.all(np.diff(pa) < 0)


# --------------------------------------------------------------------------
# Relativity
# --------------------------------------------------------------------------


def test_precession_per_orbit_is_twelve_arcminutes():
    a = op.semi_major_axis_au(S2["gm6"], S2["period"])
    dw = np.degrees(op.precession_per_orbit(S2["gm6"], a, S2["e"])) * 60
    assert 11.0 < dw < 12.5


def test_closed_form_precession_turns_the_ellipse():
    """With f_sp = 1 the pericentre direction advances by 6 pi GM / (c^2 a (1 - e^2)) per orbit."""
    t = np.array([S2["t0"], S2["t0"] + S2["period"]])
    x, y, z, *_ = op.orbit_state(t, *_el(), f_sp=1.0)
    v1, v2 = np.array([x[0], y[0], z[0]]), np.array([x[1], y[1], z[1]])
    angle = np.arccos(v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2))
    a = op.semi_major_axis_au(S2["gm6"], S2["period"])
    assert approx(angle, op.precession_per_orbit(S2["gm6"], a, S2["e"]), rel=1e-6)


@needs_data
def test_closed_form_matches_numerical_1pn_orbit():
    """The fast model reproduces a brute-force 1PN integration below the data's precision."""
    from scipy import optimize

    import s2data as sd

    ch = sd.load_do2019_chain("GR")
    w = ch["weight"]
    m = {k: float(np.sum(w * ch[k])) for k in ch if k != "weight"}
    el = (m["gm6"], m["period"], m["t0"], m["e"], m["inc"], m["omega"], m["Omega"])
    t = np.linspace(1992.0, 2020.0, 1500)
    ref = op.integrate_1pn(t, *el)

    def resid(q):
        s = np.array(op.orbit_state(t, m["gm6"], *q, f_sp=1.0)[:6])
        return np.concatenate([((s[:3] - ref[:3]) / m["r0"]).ravel(), ((s[3:] - ref[3:]) / 10).ravel()])

    q0 = [m["period"], m["t0"], m["e"], m["inc"], m["omega"], m["Omega"]]
    sol = optimize.least_squares(resid, q0, x_scale=[1e-3, 1e-4, 1e-4, 1e-2, 1e-2, 1e-2])
    s = np.array(op.orbit_state(t, m["gm6"], *sol.x, f_sp=1.0)[:6])
    assert np.max(np.abs(s[:2] - ref[:2])) / m["r0"] < 0.05  # mas on the sky (data: >= 0.15)
    assert np.max(np.abs(s[5] - ref[5])) < 3.0  # km/s along the line of sight (data: >= 13)


def test_redshift_at_s2_pericentre_is_about_200_kms():
    rp = op.pericentre_au(S2["gm6"], S2["period"], S2["e"])
    vp = op.pericentre_speed_kms(S2["gm6"], S2["period"], S2["e"])
    z = op.relativistic_redshift_kms(S2["gm6"], rp, vp)
    assert 7000 < vp < 8000
    assert 190 < z < 215


def test_shadow_size_for_published_mass_and_distance():
    # GRAVITY 2022's M and R0 give theta_g = 5.13 uas and a 53 uas shadow.
    assert approx(op.gravitational_angle_uas(4.297, 8.277), 5.13, abs=0.01)
    assert approx(op.shadow_diameter_uas(4.297, 8.277), 53.3, abs=0.2)
    # M87*: 6.5e9 Msun at 16.8 Mpc -> ~40 uas shadow (EHT ring: 42 +- 3 uas).
    assert approx(op.shadow_diameter_uas(6500.0, 16_800.0), 39.7, abs=0.5)


# --------------------------------------------------------------------------
# Data conventions (graded against Do et al.'s own chains)
# --------------------------------------------------------------------------


@needs_data
def test_do2019_newtonian_solution_fits_their_data():
    """Their chain mean, through our model and loaders, fits their data -- which
    pins every sign convention (East, x0, NIRC2 offset, velocity direction)."""
    import s2data as sd

    ka = sd.load_keck_astrometry().used()
    kr = sd.load_keck_rv()
    ch = sd.load_do2019_chain("Newton")
    m = {k: float(np.sum(ch["weight"] * ch[k])) for k in ch if k != "weight"}
    p = {"gm6": m["gm6"], "r0": m["r0"], "period": m["period"], "t0": m["t0"], "e": m["e"], "inc": m["inc"],
         "omega": m["omega"], "Omega": m["Omega"]}
    pred = op.predict(p, ka.t, kr.t, roemer=False)
    x = pred["x"] + m["x0"] + m["vx0"] * (ka.t - 2000)
    y = pred["y"] + m["y0"] + m["vy0"] * (ka.t - 2000)
    chi_a = np.sum(((ka.x - x) / ka.x_err) ** 2 + ((ka.y - y) / ka.y_err) ** 2)
    rv = pred["rv"] + m["vz0"] - m["nirc2_offset"] * (kr.instrument == "NIRC2")
    chi_v = np.sum(((kr.v - rv) / kr.v_err) ** 2)
    assert chi_a < 130  # 90 values; their likelihood also has correlated noise
    assert chi_v < 100  # 116 velocities


@needs_data
def test_roemer_delay_sign():
    """Light from S2 behind the black hole arrives late: t_em = t_obs - Z/c with Z away
    from us. The opposite sign fits Do et al.'s GR solution far worse."""
    import s2data as sd

    kr = sd.load_keck_rv()
    ch = sd.load_do2019_chain("GR")
    m = {k: float(np.sum(ch["weight"] * ch[k])) for k in ch if k != "weight"}
    el = (m["gm6"], m["period"], m["t0"], m["e"], m["inc"], m["omega"], m["Omega"])
    chis = []
    for sign in (1.0, -1.0):
        s = op.integrate_1pn(kr.t, *el)
        s = op.integrate_1pn(kr.t - sign * s[2] / op.C_AUYR, *el)
        v = np.sqrt(s[3] ** 2 + s[4] ** 2 + s[5] ** 2)
        r = np.sqrt(s[0] ** 2 + s[1] ** 2 + s[2] ** 2)
        rv = s[5] + op.relativistic_redshift_kms(m["gm6"], r, v) + m["vz0"] - m["nirc2_offset"] * (kr.instrument == "NIRC2")
        chis.append(np.sum(((kr.v - rv) / kr.v_err) ** 2))
    assert chis[0] < 90 < chis[1]


# --------------------------------------------------------------------------
# Likelihood
# --------------------------------------------------------------------------


def test_correlated_likelihood_matches_scipy():
    from scipy.stats import multivariate_normal

    rng = np.random.default_rng(3)
    n = 12
    xy = rng.normal(0, 50, (n, 2))
    dist = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1])
    rx, ry = rng.normal(size=(1, n)), rng.normal(size=(1, n))
    length, mix = np.array([20.0]), np.array([0.6])
    ours = ft._correlated_loglike(rx, ry, dist, length, mix)[0]
    rho = (1 - mix[0]) * np.eye(n) + mix[0] * np.exp(-dist / length[0])
    ref = (multivariate_normal(np.zeros(n), rho).logpdf(rx[0]) + multivariate_normal(np.zeros(n), rho).logpdf(ry[0])
           + n * np.log(2 * np.pi))  # ours omits the 2 pi constant (added by the caller)
    assert approx(ours, ref, rel=1e-10)


# --------------------------------------------------------------------------
# The dark-cluster argument (Maoz 1998)
# --------------------------------------------------------------------------


def test_maoz_worked_example():
    # Maoz: "the maximum possible lifetime of a cluster with half mass of 1e7 Msun and
    # half-mass density of 1e9 Msun/pc^3 is ~1e10 years".
    t = bh.max_lifetime(1e7, 1e9)["lifetime_yr"]
    assert 5e9 < t < 3e10


def test_maoz_best_case_objects():
    assert bh.max_lifetime(1e8, 1e12)["class"] == "neutron stars"
    wd = bh.max_lifetime(1e6, 1e12)
    assert wd["class"] == "white dwarfs" and 0.4 < wd["m_star"] < 1.0


def test_maoz_primordial_black_hole_limits():
    # Maoz sect. 3: clusters survive 10 Gyr only for m_BH < 0.005 Msun (Galaxy, 1997)
    # and < 0.04 Msun (NGC 4258).
    assert approx(bh.primordial_bh_mass_limit(1.3e6, 6e11), 0.005, rel=0.35)
    assert approx(bh.primordial_bh_mass_limit(1.7e7, 1e12), 0.04, rel=0.2)


def test_galactic_centre_1997_was_already_short_lived():
    assert bh.max_lifetime(1.3e6, 6e11)["lifetime_yr"] < 1e9


def main() -> int:
    tests = [(n, f) for n, f in globals().items() if n.startswith("test_") and callable(f)]
    failed, skipped = [], []
    for name, fn in tests:
        if getattr(fn, "needs_data", False) and not _have_data():
            skipped.append(name)
            print(f"  [skip] {name} (run fetch_data.py first)")
            continue
        for case in getattr(fn, "cases", [None]):
            label = name if case is None else f"{name}[{case}]"
            try:
                fn() if case is None else fn(case)
                print(f"  [ok  ] {label}")
            except Exception:  # noqa: BLE001
                failed.append(label)
                print(f"  [FAIL] {label}")
                traceback.print_exc(limit=2)
    print()
    if failed:
        print(f"{len(failed)} CHECK(S) FAILED: {', '.join(failed)}")
        return 1
    print(f"All {len(tests) - len(skipped)} physics checks passed" + (f" ({len(skipped)} skipped)." if skipped else "."))
    return 0


if __name__ == "__main__":
    sys.exit(main())
