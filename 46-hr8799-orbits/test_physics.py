"""Checks of the orbit mathematics and data conventions.

    python -m pytest test_physics.py -q      (or: python test_physics.py)

* our Kepler solver and sky positions match orbitize! (the package Wang et
  al.'s conventions live in) to micro-arcseconds, including e = 0.95;
* the analytic velocities match numerical derivatives, and the planet RV sign
  matches orbitize!'s (positive = receding);
* Jacobi state vectors handed to REBOUND reproduce the fitted sky positions,
  and REBOUND recovers the Jacobi semi-major axes and eccentricities;
* the astrometry table keeps each paper's sign conventions;
* state vectors -> elements inverts orbit_xyz, and the phase-space Jacobians
  that carry Wang et al.'s priors into the sampled coordinates are right
  (checked against finite differences);
* the compiled likelihood equals the numpy reference, prior terms included,
  for every model, and the diagnostics' inverse map round-trips.
"""

from __future__ import annotations

import numpy as np

import config as C
import kepler as K


def _random_elements(rng, n):
    return dict(
        a=10 ** rng.uniform(0.5, 2.0, n), e=rng.uniform(0, 0.95, n),
        inc=np.arccos(rng.uniform(-1, 1, n)), aop=rng.uniform(0, 2 * np.pi, n),
        pan=rng.uniform(0, 2 * np.pi, n), tau=rng.uniform(0, 1, n),
        plx=rng.normal(24.5, 0.5, n), mass=rng.normal(1.5, 0.1, n))


def test_sky_offsets_match_orbitize():
    from orbitize import kepler as OK
    rng = np.random.default_rng(1)
    el = _random_elements(rng, 400)
    epochs = np.array([51116.0, 55000.0, 58358.24, 61000.0])
    ra_o, dec_o, vz_o = OK.calc_orbit(epochs, el["a"], el["e"], el["inc"], el["aop"], el["pan"],
                                      el["tau"], el["plx"], el["mass"], tau_ref_epoch=C.TAU_REF_MJD)
    ra, dec = K.sky_offsets(epochs[:, None], el["a"], el["e"], el["inc"], el["aop"], el["pan"],
                            el["tau"], el["plx"], el["mass"])
    assert np.max(np.abs(ra - ra_o)) < 1e-4 and np.max(np.abs(dec - dec_o)) < 1e-4   # mas
    # planet RV: our dz/dt (au/yr -> km/s) equals orbitize's vz for the planet
    (_, _, _), (_, _, vz) = K.orbit_xyz(epochs[:, None], el["a"], el["e"], el["inc"], el["aop"],
                                        el["pan"], el["tau"], el["mass"], velocity=True)
    kms = C.AU / C.YEAR / 1e3
    assert np.max(np.abs(vz * kms - vz_o)) < 1e-3 * np.max(np.abs(vz_o))


def test_velocity_is_derivative_of_position():
    rng = np.random.default_rng(2)
    el = _random_elements(rng, 50)
    t, h = 57000.0, 0.05   # days
    args = (el["a"], el["e"], el["inc"], el["aop"], el["pan"], el["tau"], el["mass"])
    (x, y, z), (vx, vy, vz) = K.orbit_xyz(t, *args, velocity=True)
    xp, yp, zp = K.orbit_xyz(t + h, *args)
    xm, ym, zm = K.orbit_xyz(t - h, *args)
    num = np.array([xp - xm, yp - ym, zp - zm]) / (2 * h / 365.25)
    ana = np.array([vx, vy, vz])
    speed = np.sqrt((ana ** 2).sum(axis=0))
    assert np.max(np.abs(num - ana) / speed) < 1e-7      # central differences: error ~ h^2


def _hr8799_like():
    deg = np.radians
    el = {"e": (16.2, 0.12, deg(26.8), deg(110), deg(67.9), 0.12),
          "d": (26.2, 0.13, deg(26.8), deg(17), deg(67.9), 0.84),
          "c": (43.1, 0.02, deg(26.8), deg(67), deg(67.9), 0.43),
          "b": (70.8, 0.02, deg(26.8), deg(87), deg(67.9), 0.46)}
    masses = {p: m * C.MJUP_MSUN for p, m in zip("edcb", (7.2, 7.2, 7.2, 5.8))}
    return el, masses


def test_state_vectors_reproduce_sky_positions_in_rebound():
    import rebound
    el, masses = _hr8799_like()
    t0, plx, mstar = C.NBODY_EPOCH_MJD, 24.3, 1.47
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    for m, x, y, z, vx, vy, vz in K.state_vectors(t0, el, mstar, masses):
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    star = sim.particles[0]
    sky = K.system_offsets(t0, el, plx, mstar, masses)
    for k, p in enumerate(C.INNER_TO_OUTER, start=1):
        q = sim.particles[k]
        assert abs((q.x - star.x) * plx - sky[p][0]) < 1e-9
        assert abs((q.y - star.y) * plx - sky[p][1]) < 1e-9
    # REBOUND's own Jacobi elements (primary = barycentre of the inner bodies)
    for k, p in enumerate(C.INNER_TO_OUTER, start=1):
        com = sim.com(last=k)
        o = sim.particles[k].orbit(primary=com)
        assert abs(o.a - el[p][0]) < 1e-9 * el[p][0] and abs(o.e - el[p][1]) < 1e-9


def test_test_particle_limit_follows_kepler():
    """Tiny masses: a 50-yr N-body run must track the Keplerian prediction."""
    import rebound
    el, _ = _hr8799_like()
    tiny = {p: 1e-12 for p in el}
    t0, mstar = C.NBODY_EPOCH_MJD, 1.47
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.G = C.G_AU_YR
    for m, x, y, z, vx, vy, vz in K.state_vectors(t0, el, mstar, tiny):
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.integrator = "ias15"
    sim.integrate(50.0)
    t1 = t0 + 50.0 * 365.25
    sky = K.system_offsets(t1, el, 1.0, mstar, tiny)
    star = sim.particles[0]
    for k, p in enumerate(C.INNER_TO_OUTER, start=1):
        q = sim.particles[k]
        assert abs((q.x - star.x) - sky[p][0]) < 1e-6 and abs((q.y - star.y) - sky[p][1]) < 1e-6


def test_astrometry_conventions():
    import astrometry as A
    df = A.load(include_excluded=True)
    k = df[(df.source == "NIRC2-K16") & (df.date == "2014-07-17")].set_index("planet")
    g = df[(df.source == "GPI-W18") & (df.date == "2014-09-12")].set_index("planet")
    # Keck (Delta x flipped) and GPI (sep, PA) agree on b's direction within 0.5 deg
    assert abs(k.loc["b", "pa"] - g.loc["b", "pa"]) < 0.5
    assert len(A.dataset("W18")) == 55


def test_elements_from_state_inverts_orbit_xyz():
    rng = np.random.default_rng(5)
    el = _random_elements(rng, 2000)
    (x, y, z), (vx, vy, vz) = K.orbit_xyz(C.STATE_REF_MJD, el["a"], el["e"], el["inc"], el["aop"],
                                          el["pan"], el["tau"], el["mass"], velocity=True)
    a, e, i, w, O, M = K.elements_from_state(np.stack([x, y, z], -1), np.stack([vx, vy, vz], -1),
                                             C.G_AU_YR * el["mass"])
    P = K.period_days(el["a"], el["mass"])
    M0 = np.mod(K.mean_anomaly(C.STATE_REF_MJD, P, el["tau"]), 2 * np.pi)
    dang = lambda u, v: np.abs(np.mod(u - v + np.pi, 2 * np.pi) - np.pi)
    assert np.max(np.abs(a / el["a"] - 1)) < 1e-10 and np.max(np.abs(e - el["e"])) < 1e-10
    for got, want in ((i, el["inc"]), (O, el["pan"]), (w, el["aop"]), (M, M0)):
        assert np.max(dang(got, want)) < 1e-9


def _state(el6, mu):
    a, e, i, O, w, Mn = el6
    mass = mu / C.G_AU_YR
    tau = np.mod(np.mod((0.0 - C.TAU_REF_MJD) / K.period_days(a, mass), 1.0) - Mn / (2 * np.pi), 1.0)
    (x, y, z), (vx, vy, vz) = K.orbit_xyz(0.0, a, e, i, w, O, tau, mass, velocity=True)
    return np.array([x, y, z, vx, vy, vz])


def test_phase_space_jacobians():
    """d^3x d^3v = mu^1.5 a^0.5 e sin(i)/2 da de di dO dw dM, and the planar
    (sky-projected) version mu e cos(i)^2 / (2 sqrt(1-e^2)) da de dM dw."""
    rng = np.random.default_rng(6)
    for _ in range(100):
        a = 10 ** rng.uniform(0.5, 2)
        e = rng.uniform(0.01, 0.9)
        i = np.arccos(rng.uniform(-0.95, 0.95))
        O, w, Mn = rng.uniform(0, 2 * np.pi, 3)
        mu = C.G_AU_YR * rng.uniform(1, 2)
        el6 = np.array([a, e, i, O, w, Mn])
        h = 1e-6 * np.array([a, 1, 1, 1, 1, 1])
        cols = [(_state(el6 + h[k] * np.eye(6)[k], mu) - _state(el6 - h[k] * np.eye(6)[k], mu)) / (2 * h[k])
                for k in range(6)]
        J = np.column_stack(cols)
        assert abs(abs(np.linalg.det(J)) / (0.5 * mu ** 1.5 * np.sqrt(a) * e * np.sin(i)) - 1) < 1e-5
        J4 = np.column_stack([cols[k][[0, 1, 3, 4]] for k in (0, 1, 5, 4)])
        assert abs(abs(np.linalg.det(J4)) / (mu * e / (2 * np.sqrt(1 - e * e)) * np.cos(i) ** 2) - 1) < 1e-4


def test_compiled_likelihood_matches_numpy():
    import astrometry as A
    import diagnostics as D
    import models as M
    rng = np.random.default_rng(7)
    for name, kw in (("W18", {}), ("ALL", dict(parallax=C.DR3_PARALLAX, mstar=C.ALL_MSTAR, fit_masses=True))):
        data = M.Data.from_frame(A.dataset(name), calibrate=(name == "ALL"))
        for kind in C.MODELS:
            m = M.OrbitModel(kind, data, **kw)
            n_valid = 0
            for _ in range(3000):
                th = m.prior_transform(rng.uniform(size=m.ndim))
                a, b = m.loglike(th), m.loglike_reference(th)
                assert (a > -1e11) == (b > -1e299)
                if a > -1e11:
                    n_valid += 1
                    assert abs(a - b) <= 1e-9 * max(1.0, abs(b))
                    th2 = D.theta_from_physical(m, D._single(m.physical(th[None, :]), 0))
                    if th2 is not None:
                        assert abs(m.loglike(th2) - a) < 1e-6 * max(1.0, abs(a))
            assert n_valid > 0 or kind.startswith("near1248")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok ", name)
