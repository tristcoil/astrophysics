"""Self-checks for the geometry, the velocity frames and the edge fit.

Run with ``python test_physics.py``. No test framework and no network. Most
checks need no survey data; the last group checks the cached HI4PI mosaic and
is skipped if the data have not been fetched.

Everything this project concludes passes through a handful of sign
conventions: which way the Galaxy turns, which way longitude runs, what the
survey's "LSR" velocity frame means, and how a terminal velocity becomes a
rotation speed. A sign slip in any of them produces plots that still look
reasonable, so each is pinned here to something known independently: a
numerical search along a sightline, a published number, or a round trip.
"""

from __future__ import annotations

import sys

import numpy as np

import config as cfg
import galaxy
import terminal

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (absolute) of ``expected``."""
    got = float(got)
    ok = abs(got - expected) <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<58} {got:>11.4f} vs {expected:>10.4f} "
          f"{unit:<6} (tol {tol:g})")


def check_true(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def test_frames() -> None:
    print("Velocity frames")
    std = galaxy.standard_solar_motion()
    # the standard solar motion in Galactic components, as tabulated e.g. by
    # Reid et al. (2009): (10.3, 15.3, 7.7) km/s
    check("standard solar motion U", std[0], 10.27, 0.02, "km/s")
    check("standard solar motion V", std[1], 15.32, 0.02, "km/s")
    check("standard solar motion W", std[2], 7.74, 0.02, "km/s")
    obs = galaxy.Observer.reid19()
    # the LSRK frame's orbital speed: the Sun's (Omega * R0) minus the standard V
    check("LSRK frame speed = 30.32 x 8.15 - 15.32", obs.v, 30.32 * 8.15 - std[1], 1e-6, "km/s")
    # textbook observer reproduces v = (Theta R0/R - Theta0) sin l
    tb = galaxy.Observer.textbook(8.5, 220.0)
    l, d = 40.0, 3.0
    R = galaxy.radius_from_ld(l, d, 8.5)
    v_classic = (220.0 * 8.5 / R - 220.0) * np.sin(np.radians(l))
    check("textbook frame = classic formula",
          galaxy.v_lsr(l, 0, d, lambda r: galaxy.theta_flat(r, 220.0), tb), v_classic, 1e-9, "km/s")


def test_rotation_curves() -> None:
    print("Rotation curves")
    # Reid et al. (2019) Table 3, fit A5: Theta0 = 236 km/s from their curve
    check("Reid 2019 universal curve at R0", galaxy.theta_reid19(cfg.R0_KPC), 236.0, 0.5, "km/s")
    # Clemens (1985) was built for Theta0 = 220 at R0 = 8.5; the polynomial
    # passes within ~1 km/s of it
    check("Clemens 1985 polynomial at 8.5 kpc", galaxy.theta_clemens85(np.array([8.5]))[0],
          220.0, 1.5, "km/s")
    check("MGD 2016 joint fit at R0 = 8.5 (0.171 + 0.889)", galaxy.theta_mgd(8.5), 220 * 1.06,
          1e-9, "km/s")
    check("mass inside 8.15 kpc at 236 km/s (1e11 Msun)",
          galaxy.enclosed_mass(8.15, 236.0) / 1e11, 1.0555, 0.001)


def test_terminal_geometry() -> None:
    print("Terminal velocity = the extreme velocity along the sightline")
    obs = galaxy.Observer.reid19()
    for l in (25.0, 50.0, 70.0, 290.0, 320.0, 340.0):
        d = np.linspace(1e-3, 2 * obs.r0 * abs(np.cos(np.radians(l))), 400001)
        v = galaxy.v_lsr(l, 0.0, d, galaxy.theta_reid19, obs)
        s = np.sign(np.sin(np.radians(l)))
        i = np.argmax(s * v)
        # the extreme is at the tangent point, d = R0 cos l ...
        check(f"l = {l:5.1f}: extreme sits at d = R0 cos l", d[i], obs.r0 * np.cos(np.radians(l)),
              2e-3, "kpc")
        # ... and inverting it recovers the rotation speed there
        check(f"l = {l:5.1f}: Theta from V_t", galaxy.theta_from_terminal(l, v[i], obs),
              galaxy.theta_reid19(galaxy.tangent_radius(l, obs.r0)), 1e-6, "km/s")


def test_outer_distance() -> None:
    print("Kinematic distances outside the solar circle")
    obs = galaxy.Observer.reid19()
    for l, d in ((60.0, 12.0), (120.0, 4.0), (150.0, 9.0), (230.0, 6.0), (300.0, 20.0)):
        v = galaxy.v_lsr(l, 0.0, d, galaxy.theta_reid19, obs)
        R, d_back = galaxy.outer_distance(l, 0.0, v, galaxy.theta_reid19, obs)
        check(f"l = {l:5.1f}, d = {d:4.1f}: round trip", d_back, d, 0.01, "kpc")


def test_arms() -> None:
    print("Spiral-arm model against Reid et al. (2019) Table 2 tangencies")
    obs = galaxy.Observer.reid19()
    published = {"Norma": 327.5, "Scutum-Centaurus": 306.1, "Sagittarius-Carina": 285.6}
    for arm in galaxy.REID19_ARMS:
        if arm.name not in published:
            continue
        beta = np.linspace(-80, 5, 20001)
        x, y = arm.xy(beta)
        l, _ = galaxy.ld_from_xy(x, y, obs.r0)
        fourth = (l > 270) & (l < 360)
        check(f"{arm.name}: fourth-quadrant tangent longitude", l[fourth].min(),
              published[arm.name], 1.0, "deg")


def test_edge_fit() -> None:
    print("Edge fit recovers a known edge (synthetic spectrum, for testing only)")
    rng = np.random.default_rng(0)
    v = np.arange(-50, 200, 1.288)
    for sign, v0 in ((1, 110.0), (-1, -80.0)):
        vv = v if sign > 0 else -v
        T = terminal.two_erfc(sign * vv, sign * v0, 35.0, 6.0, 8.0, 20.0)
        T = T + 0.05 * rng.standard_normal(v.size)
        f = terminal.fit_edge(vv, T, sign)
        check(f"edge at {v0:+.0f} km/s", f.v_terminal, v0, 0.1, "km/s")
    u = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    t = np.array([30.0, 25.0, 20.0, 5.0, 1.0])
    check("threshold crossing interpolates between channels",
          terminal.threshold_crossing(u, t, 10.0, min_run=2), 2.0 + 10 / 15, 1e-9)


def test_masers() -> None:
    print("Maser orbital speeds: a star on a circular orbit comes back circular")
    import astropy.units as u
    from astropy.coordinates import CartesianDifferential, Galactocentric, ICRS, SkyCoord

    import reference

    frame = Galactocentric(galcen_distance=cfg.R0_KPC * u.kpc, z_sun=5.5 * u.pc,
                           galcen_v_sun=CartesianDifferential(
                               [cfg.U_SUN_KMS, cfg.OMEGA_SUN * cfg.R0_KPC, cfg.W_SUN_KMS]
                               * u.km / u.s))
    # stars at R = 6 and 11 kpc; clockwise seen from +z means v_y < 0 at +x
    xs, ys = np.array([-6.0, 3.0]), np.array([0.0, -10.57])
    R = np.hypot(xs, ys)
    th = galaxy.theta_reid19(R)
    vx, vy = th * ys / R, -th * xs / R
    c = SkyCoord(x=xs * u.kpc, y=ys * u.kpc, z=np.zeros(2) * u.kpc,
                 v_x=vx * u.km / u.s, v_y=vy * u.km / u.s, v_z=np.zeros(2) * u.km / u.s,
                 frame=frame).transform_to(ICRS())
    std = galaxy.standard_solar_motion()
    gal = c.galactic
    lr, br = np.radians(gal.l.deg), np.radians(gal.b.deg)
    vlsr = c.radial_velocity.to_value(u.km / u.s) + (
        (std[0] * np.cos(lr) + std[1] * np.sin(lr)) * np.cos(br) + std[2] * np.sin(br))
    m = reference.Masers(
        name=np.array(["t1", "t2"]), arm=np.array(["x", "x"]), l=gal.l.deg, b=gal.b.deg,
        plx=1 / c.distance.to_value(u.kpc), e_plx=np.full(2, 1e-9),
        pm_e=c.pm_ra_cosdec.to_value(u.mas / u.yr), e_pm_e=np.full(2, 1e-9),
        pm_n=c.pm_dec.to_value(u.mas / u.yr), e_pm_n=np.full(2, 1e-9),
        vlsr=vlsr, e_vlsr=np.full(2, 1e-9), ra=c.ra.deg, dec=c.dec.deg)
    out = reference.maser_rotation(m, n_mc=5)
    for i in range(2):
        check(f"star at R = {R[i]:.1f} kpc: orbital speed", out["vrot"][i], th[i], 1e-3, "km/s")
        check(f"star at R = {R[i]:.1f} kpc: radial speed", out["vrad"][i], 0.0, 1e-3, "km/s")


def test_survey_mosaic() -> None:
    print("HI4PI mosaic (cached data)")
    from fetch_data import missing_tiles

    if missing_tiles():
        print("  [skip] survey cubes not fetched")
        return
    import survey

    seen = np.zeros(survey.N_LON, dtype=int)
    for t in survey.tiles():
        _, k = survey.own_columns(t)
        seen[k] += 1
    check_true("every 5' longitude column comes from exactly one cube", bool(np.all(seen == 1)))
    e = survey.open_tile("CAR_E04.fits")
    f = survey.open_tile("CAR_F04.fits")
    top_e = e.lat[survey.own_rows(e)].max()
    bottom_f = f.lat[survey.own_rows(f)].min()
    check("rows E and F meet without overlap (gap = 1 pixel)", (bottom_f - top_e) * 12, 1.0, 1e-6,
          "pix")


if __name__ == "__main__":
    for test in (test_frames, test_rotation_curves, test_terminal_geometry, test_outer_distance,
                 test_arms, test_edge_fit, test_masers, test_survey_mosaic):
        test()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
        sys.exit(1)
    print("All checks passed.")
