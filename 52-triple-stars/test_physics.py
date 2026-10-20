"""Checks of the machinery on cases with known answers.

    python test_physics.py
"""

from __future__ import annotations

import numpy as np

import config as C


def test_occultation_uniform_disc():
    """Occulted fraction of a uniform disc = analytic overlap area / pi."""
    import eclipse_model as E

    def overlap(z, p):
        if z >= 1 + p:
            return 0.0
        if z <= abs(1 - p):
            return min(p, 1) ** 2
        k0 = np.arccos((p * p + z * z - 1) / (2 * p * z))
        k1 = np.arccos((1 - p * p + z * z) / (2 * z))
        return (p * p * k0 + k1 - 0.5 * np.sqrt(max(0, 4 * z * z - (1 + z * z - p * p) ** 2))) / np.pi
    rng = np.random.default_rng(1)
    worst = max(abs(E.occulted_fraction(z, p, 0.0, 0.0) - overlap(z, p))
                for z, p in zip(rng.uniform(0, 2.2, 5000), rng.uniform(0.05, 1.5, 5000)))
    assert worst < 5e-5, worst
    return f"max |error| {worst:.1e}"


def test_occultation_limb_darkened():
    """Against a brute-force pixel integration of a quadratically darkened disc."""
    import eclipse_model as E
    x = np.linspace(-1, 1, 2501)
    X, Y = np.meshgrid(x, x)
    R2 = X**2 + Y**2
    mu = np.sqrt(np.clip(1 - R2, 0, 1))
    worst = 0.0
    for u1, u2 in ((0.45, 0.25), (0.2, 0.3)):
        I = (1 - u1 * (1 - mu) - u2 * (1 - mu) ** 2) * (R2 < 1)
        for z, p in ((0.0, 0.5), (0.3, 0.8), (0.95, 0.9), (1.3, 0.6), (0.2, 1.1)):
            brute = (I * (((X - z) ** 2 + Y**2) < p * p)).sum() / I.sum()
            worst = max(worst, abs(E.occulted_fraction(z, p, u1, u2) - brute))
    assert worst < 2e-4, worst
    return f"max |error| {worst:.1e}"


def test_precessing_orbit_geometry():
    """sky_positions: with no precession the orbit's normal makes the angle
    the inclination formula says, and positions stay on the unit circle."""
    import eclipse_model as E
    t = np.linspace(0, 10, 2001)
    eps, I, h0 = np.radians(17.0), np.radians(70.0), np.radians(35.0)
    xyz = E.sky_positions(t, 0.0, 0.3, 2 * np.pi / 3.0, h0, 0.0, 0.0, 0.0, 0.0, eps, I)
    r = np.linalg.norm(xyz, axis=1)
    n = np.cross(xyz[10], xyz[30])
    n /= np.linalg.norm(n)
    i_num = np.degrees(np.arccos(n[2]))
    i_formula = np.degrees(E.inclination_to_sky(0.0, 0.0, h0, 0.0, eps, I))
    assert np.allclose(r, 1, atol=1e-12) and abs(i_num - i_formula) < 1e-8, (i_num, i_formula)
    return f"i = {i_formula:.4f} deg both ways"


def test_conjunction_lambda():
    """The mean longitude chosen for a primary eclipse puts the secondary in
    front at the smallest sky separation."""
    import eclipse_model as E
    import switching as S
    for e, g in ((0.0, 0.0), (0.136, 2.0)):
        h, eps, I = 0.7, np.radians(15), np.radians(84)
        lam = S._conj_lambda(h, eps, I, g, e)
        n = 2 * np.pi / 14.4
        t = np.linspace(-0.5, 0.5, 20001)
        xyz = E.sky_positions(t, 0.0, lam, n, h, 0.0, g, 0.0, e, eps, I)
        rho = np.hypot(xyz[:, 0], xyz[:, 1])
        rho[xyz[:, 2] < 0] = np.inf
        tmin = t[np.argmin(rho)]
        assert abs(tmin) < 2e-3 * 14.4, (e, tmin)
    return "minimum within 0.2% of a period"


def test_secular_rate_vs_nbody():
    """The quadrupole nodal-precession rate against a direct REBOUND
    integration for an SS-Lac-like triple (expect agreement to ~10%,
    Soderhjelm 1982)."""
    import rebound
    import switching as S
    m1, m2, m3 = 2.6, 2.5, 0.8
    P_in, P_out, e_out = 14.4, 678.8, 0.159
    J = np.radians(27.0)
    eps_in, eps_out, nu, gam, a_in, a_out = S.secular(m1, m2, m3, P_in, P_out, 0.0, e_out, J)
    sim = rebound.Simulation()
    sim.G = C.G_AU_DAY
    sim.add(m=m1)
    sim.add(m=m2, P=P_in, e=0.0, inc=0.0)
    sim.add(m=m3, P=P_out, e=e_out, inc=J, Omega=0.0)
    sim.move_to_com()
    L = np.zeros(3)
    for p in sim.particles:
        L += p.m * np.cross([p.x, p.y, p.z], [p.vx, p.vy, p.vz])
    L /= np.linalg.norm(L)
    sim.integrator = "whfast"
    sim.dt = P_in / 40
    ts = np.linspace(0, 60 * 365.25, 61)
    ang = []
    e1 = np.cross(L, [1, 0, 0]); e1 /= np.linalg.norm(e1)
    e2 = np.cross(L, e1)
    for t in ts:
        sim.integrate(t)
        r = np.array([sim.particles[1].x - sim.particles[0].x, sim.particles[1].y - sim.particles[0].y,
                      sim.particles[1].z - sim.particles[0].z])
        v = np.array([sim.particles[1].vx - sim.particles[0].vx, sim.particles[1].vy - sim.particles[0].vy,
                      sim.particles[1].vz - sim.particles[0].vz])
        h = np.cross(r, v)
        ang.append(np.arctan2(h @ e2, h @ e1))
    rate = np.polyfit(ts, np.unwrap(ang), 1)[0]
    ratio = rate / nu
    assert 0.85 < ratio < 1.15, ratio
    return f"N-body / secular nodal rate = {ratio:.3f} (P_node {2 * np.pi / abs(rate) / 365.25:.0f} yr)"


def test_ma01():
    """Mardling & Aarseth (2001): equal masses, circular, coplanar ->
    a_out/a_in = 2.8 * 1.5^0.4 = 3.29; retrograde lowers it by 30%."""
    r = float(C.ma01_ratio(0.5, 0.0, 0.0))
    r_retro = float(C.ma01_ratio(0.5, 0.0, 180.0))
    assert abs(r - 2.8 * 1.5**0.4) < 1e-12 and abs(r_retro / r - 0.7) < 1e-12
    return f"{r:.3f} (prograde), {r_retro:.3f} (retrograde)"


def test_mutual_inclination():
    """Both node solutions against vector arithmetic."""
    import catalogs as K
    import nbody as NB
    rng = np.random.default_rng(2)
    for _ in range(200):
        i1, O1, i2, O2 = rng.uniform(0, 180), rng.uniform(0, 360), rng.uniform(0, 180), rng.uniform(0, 360)
        p1, p2 = K.mutual_inclinations(i1, O1, i2, O2)
        n1 = NB.unit_normal(np.radians(i1), np.radians(O1))
        n2 = NB.unit_normal(np.radians(i2), np.radians(O2))
        n2b = NB.unit_normal(np.radians(i2), np.radians(O2 + 180))
        a = np.degrees(np.arccos(np.clip(n1 @ n2, -1, 1)))
        b = np.degrees(np.arccos(np.clip(n1 @ n2b, -1, 1)))
        assert abs(p1 - a) < 1e-8 and abs(p2 - b) < 1e-8
    return "200 random pairs"


def test_proxima_state_roundtrip():
    """Elements from a state vector reproduce the state (Kepler, visual-binary
    convention)."""
    import proxima as X
    rng = np.random.default_rng(3)
    ra, dec = 219.9, -60.8
    u, e, n = X._basis(np.array([ra]), np.array([dec]))
    for _ in range(50):
        r = rng.normal(size=3) * 8000
        v = rng.normal(size=3) * 0.15
        el = X.elements(r[None, :], v[None, :], np.array([2.1]), (u, e, n))
        a, ecc, i, Om, w = (float(el[k][0]) for k in ("a", "e", "i", "Omega", "omega"))
        if a <= 0:
            continue
        rr = np.linalg.norm(r)
        # reconstruct the present position from the elements
        M = -2 * np.pi * float(el["t_to_peri_yr"][0]) / float(el["P_yr"][0])
        E_ = M
        for _ in range(60):
            E_ = E_ - (E_ - ecc * np.sin(E_) - M) / (1 - ecc * np.cos(E_))
        nu = 2 * np.arctan2(np.sqrt(1 + ecc) * np.sin(E_ / 2), np.sqrt(1 - ecc) * np.cos(E_ / 2))
        rad = a * (1 - ecc * np.cos(E_))
        uu = nu + np.radians(w)
        Omr, ir = np.radians(Om), np.radians(i)
        N = rad * (np.cos(uu) * np.cos(Omr) - np.sin(uu) * np.sin(Omr) * np.cos(ir))
        Ea = rad * (np.cos(uu) * np.sin(Omr) + np.sin(uu) * np.cos(Omr) * np.cos(ir))
        Z = rad * np.sin(uu) * np.sin(ir)
        assert np.allclose([N, Ea, Z], [r @ n[0], r @ e[0], r @ u[0]], atol=1e-6 * rr), (N, Ea, Z)
    return "50 random states"


def test_eclipse_finder_two_body():
    """TIC 290061484 code on an isolated binary: eclipses every period,
    exactly where Kepler puts them."""
    import rebound
    import tic290_nbody as N
    sim = rebound.Simulation()
    sim.G = C.G_AU_DAY
    sim.t = 2458711.0
    sim.add(m=6.85)
    sim.add(m=6.11, P=1.79, e=0.0, inc=np.radians(86.0), Omega=0.0, theta=np.pi / 2)  # Ab in front now
    sim.add(m=1e-12, P=1e5, e=0.0)   # a massless placeholder for "B"
    sim.add(m=1e-12, P=1e7, e=0.0)
    sim.move_to_com()
    tm, kind, _ = N.all_eclipses(sim, 2458711.0 + 20, per_orbit=120)
    pri = tm[kind > 0]
    d = np.diff(pri)
    assert np.allclose(d, 1.79, atol=2e-5) and abs(pri[0] - (2458711.0 + 1.79)) < 2e-4, (pri[:3], d)
    return f"spacing {d.mean():.6f} d"


def test_dasch_flags():
    """Our DASCH rejection masks are exactly daschlab 1.0's flags (1-based
    bits 7, 12, 13, 14, 16) plus 8, 10, 11."""
    std = sum(1 << (b - 1) for b in (7, 12, 13, 14, 16))
    extra = sum(1 << (b - 1) for b in (8, 10, 11))
    assert C.DASCH_BAD_AFLAGS == std and C.DASCH_EXTRA_BAD_AFLAGS == extra
    return "masks match"


def test_escaper_identification():
    """A tight pair and a third star leaving: the escaper is the third star,
    not a member of the pair (whatever their order in the simulation)."""
    import rebound
    import stability as St
    for order in ((0, 1, 2), (2, 0, 1), (1, 2, 0)):
        sim = rebound.Simulation()
        sim.G = 1.0
        stars = {0: dict(m=6.85), 1: dict(m=6.11), 2: dict(m=7.90)}
        pair = rebound.Simulation()
        pair.G = 1.0
        pair.add(m=6.85)
        pair.add(m=6.11, a=0.5, e=0.6, M=1.0)
        pair.move_to_com()
        ps = pair.particles
        stars[0].update(x=ps[0].x, y=ps[0].y, vx=ps[0].vx, vy=ps[0].vy)
        stars[1].update(x=ps[1].x, y=ps[1].y, vx=ps[1].vx, vy=ps[1].vy)
        stars[2].update(x=-80.0, y=30.0, vx=-1.2, vy=0.4)       # B, far and leaving
        for k in order:
            sim.add(**stars[k])
        found = St.find_escaper(sim, 2.5)
        assert found is not None and order[found[0]] == 2, (order, found)
        assert abs(found[2].a - 0.5) < 1e-6
    # nobody leaving: a bound hierarchical triple
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.add(m=6.85)
    sim.add(m=6.11, a=1.0)
    sim.add(m=7.90, a=6.0, e=0.2)
    assert St.find_escaper(sim, 6.0) is None
    return "escaper found in all orders; bound triple: none"


def test_tic290_timing_converged():
    """Model eclipse times do not depend on how finely the orbit is sampled
    (the minimum is refined; IAS15), for Kostov et al.'s solution in 2019."""
    import tic290_nbody as N
    obs = N.observed_times()
    obs = obs[obs["sector"].isin((15, 16))].reset_index(drop=True)
    t1 = N.model_times_cycles(N.build(None), obs, fine=120)
    t2 = N.model_times_cycles(N.build(None), obs, fine=360)
    d = np.nanmax(np.abs(t1 - t2)) * 86400
    assert np.isfinite(t1).all() and d < 3.0, d
    return f"max difference {d:.2f} s between 1/120 and 1/360-orbit sampling"


TESTS = [test_occultation_uniform_disc, test_occultation_limb_darkened, test_precessing_orbit_geometry,
         test_conjunction_lambda, test_secular_rate_vs_nbody, test_ma01, test_mutual_inclination,
         test_proxima_state_roundtrip, test_eclipse_finder_two_body, test_dasch_flags,
         test_escaper_identification, test_tic290_timing_converged]

if __name__ == "__main__":
    n_ok = 0
    for t in TESTS:
        try:
            msg = t()
            n_ok += 1
            print(f"PASS {t.__name__}: {msg}")
        except AssertionError as exc:
            print(f"FAIL {t.__name__}: {exc}")
    print(f"{n_ok}/{len(TESTS)} passed")
