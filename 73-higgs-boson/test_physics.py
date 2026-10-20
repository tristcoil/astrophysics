"""Self-tests of the physics and numerics. Run with:  python test_physics.py

Each test checks a piece of the machinery against something known exactly or
published: kinematics, the outreach code's integer quirk, the loop functions,
the matching and running against Buttazzo et al. (2013), the light-cone volume
against its closed form, the integrator against an adaptive one, and the fast
pseudo-experiment statistic against the full fit.
"""

from __future__ import annotations

import math
import sys

import numpy as np
from scipy.integrate import quad, solve_ivp

import events
import rge
import stability as st

TESTS = []


def test(f):
    TESTS.append(f)
    return f


@test
def invariant_mass_back_to_back():
    # two massless leptons, 45.6 GeV each, back to back in the transverse plane
    v1 = events.p4(np.array([45.6]), np.array([0.0]), np.array([0.0]), np.array([0.0]))
    v2 = events.p4(np.array([45.6]), np.array([0.0]), np.array([np.pi]), np.array([0.0]))
    assert abs(events.mass(v1 + v2)[0] - 91.2) < 1e-12
    # boost invariance: the same pair moved along the beam (both eta + 1.3)
    v1 = events.p4(np.array([45.6]), np.array([1.3]), np.array([0.0]), np.array([0.0]))
    v2 = events.p4(np.array([45.6]), np.array([1.3]), np.array([np.pi]), np.array([0.0]))
    assert abs(events.mass(v1 + v2)[0] - 91.2) < 1e-9
    # massless pair in general: m^2 = 2 pT1 pT2 (cosh(d eta) - cos(d phi))
    v1 = events.p4(np.array([30.0]), np.array([1.0]), np.array([0.4]), np.array([0.0]))
    v2 = events.p4(np.array([20.0]), np.array([-0.7]), np.array([2.9]), np.array([0.0]))
    m = math.sqrt(2 * 30 * 20 * (math.cosh(1.7) - math.cos(2.5)))
    assert abs(events.mass(v1 + v2)[0] - m) < 1e-9


@test
def delta_r_wraps_phi():
    assert abs(events.delta_r(0.0, 3.1, 0.0, -3.1) - (2 * np.pi - 6.2)) < 1e-12


@test
def integer_pairing_quirk():
    """skim.cxx keeps the best Z mass in an int: a pair at 91.9 GeV (stored as
    91, |91.2 - 91| = 0.2) beats a later pair at 90.8 GeV (0.4), although the
    latter is closer to 91.2 than 91.9 is."""
    def lepton_pair(m, eta0):
        # two massless leptons with invariant mass m, back to back in phi
        pt = m / 2
        return [(pt, eta0, 0.0), (pt, eta0, np.pi)]
    # leptons 0,1: pair A (91.9 GeV); leptons 2,3: pair B (90.8 GeV); charges + - + -
    leps = lepton_pair(91.9, 0.0) + lepton_pair(90.8, 2.0)
    pt = np.array([[l[0] for l in leps]])
    eta = np.array([[l[1] for l in leps]])
    phi = np.array([[l[2] for l in leps]])
    v = events.p4(pt, eta, phi, np.zeros_like(pt))
    charge = np.array([[1, -1, 1, -1]])
    i1, i2 = events.pair_closest_to_z(v, charge, emulate_int=True)
    j1, j2 = events.pair_closest_to_z(v, charge, emulate_int=False)
    # combinations (0,1), (0,2), (0,3), (1,2), (1,3), (2,3): masses of the cross pairs
    # are far from 91.2 here, so the choice is between (0,1) and (2,3)
    assert (i1[0], i2[0]) == (0, 1), (i1, i2)
    assert (j1[0], j2[0]) == (2, 3), (j1, j2)


@test
def passarino_veltman():
    mu = 173.34
    # A0 definition
    assert abs(rge.A0(80.0, mu) - 80.0**2 * (1 - math.log(80.0**2 / mu**2))) < 1e-9
    # B0 with a massless line: closed form against direct quadrature
    p, M = 173.34, 80.384
    f = lambda x: -math.log(abs((1 - x) * (M**2 - x * p**2)) / mu**2)
    num = quad(f, 0, 1, points=[M**2 / p**2], limit=400)[0]
    assert abs(rge.B0(p, 0.0, M, mu) - num) < 1e-9
    # continuous through the threshold p = 2m (a double root of the integrand); the
    # real part has a square-root cusp there, ~ pi sqrt(2 dp / p) below threshold
    at = rge.B0(140.0, 70.0, 70.0, mu)
    for dp in (-1e-9, 1e-9):
        assert math.isfinite(at) and abs(rge.B0(140.0 + dp, 70.0, 70.0, mu) - at) < 1e-4
    # closed form against direct quadrature away from thresholds
    for p_, m1, m2 in ((125.0, 80.0, 91.0), (173.0, 125.0, 173.0), (91.0, 30.0, 40.0)):
        a, b, c = p_**2, m1**2 - m2**2 - p_**2, m2**2
        d = b * b - 4 * a * c
        pts = [r for r in (((-b - d**0.5) / (2 * a), (-b + d**0.5) / (2 * a)) if d > 0 else ())
               if 0 < r < 1]
        num = quad(lambda x: -math.log(abs(x * m1**2 + (1 - x) * m2**2 - x * (1 - x) * p_**2)
                                       / mu**2), 0, 1, points=pts or None, limit=400)[0]
        assert abs(rge.B0(p_, m1, m2, mu) - num) < 1e-9
    # B0(0; m, m) = -ln(m^2/mu^2)
    assert abs(rge.B0(1e-9, 91.0, 91.0, mu) + math.log(91.0**2 / mu**2)) < 1e-6


def _buttazzo_inputs():
    R = rge.REF
    return rge.Inputs(Mh=R["Mh"], Mt=R["Mt"], MW=R["MW"], MZ=R["MZ"], alpha_s=R["alpha_s"],
                      GF=1 / (math.sqrt(2) * R["V"] ** 2))


@test
def matching_reproduces_buttazzo_table_3():
    """Table 3 of Buttazzo et al. (2013): couplings at mu = M_t, LO/NLO/NNLO."""
    inp = _buttazzo_inputs()
    paper = {"LO": (0.12917, 0.99561, 0.65294, 0.34972),
             "NLO": (0.12774, 0.95113, 0.64754, 0.35940),
             "NNLO": (0.12604, 0.94018, 0.64779, 0.35830)}
    for order, (lam, yt, g2, gy) in paper.items():
        c = rge.match(inp, order)
        if order == "NNLO":     # the table's NNLO y_t excludes the 3-loop QCD term
            a = c["g3"] ** 2 / (4 * math.pi**2)
            yt_ours = c["yt"] + 80.4045 * a**3 * math.sqrt(2) * inp.Mt / inp.V
        else:
            yt_ours = c["yt"]
        for ours, theirs, name in ((c["lam"], lam, "lambda"), (yt_ours, yt, "yt"),
                                   (c["g2"], g2, "g2"), (c["gY"], gy, "gY")):
            assert abs(ours - theirs) < 1.5e-5, (order, name, ours, theirs)
    # their eq. (3.3) with the three-loop QCD term
    assert abs(rge.match(inp, "NNLO")["yt"] - 0.93690) < 1e-5


@test
def gy_typo_reading():
    """The printed 36 M_W^2/(M_Z^2-M_W^2) in the one-loop gY threshold must read
    36 M_W^4/(...) to reproduce the paper's own NLO value 0.35940."""
    inp = _buttazzo_inputs()
    R = rge.REF
    tree = 2 * math.sqrt(R["MZ"] ** 2 - R["MW"] ** 2) / inp.V
    good = tree + rge.gY_1loop(R["Mh"], R["Mt"], R["MW"], R["MZ"], inp.V, R["Mt"], True)
    bad = tree + rge.gY_1loop(R["Mh"], R["Mt"], R["MW"], R["MZ"], inp.V, R["Mt"], False)
    assert abs(good - 0.35940) < 2e-5 and abs(bad - 0.35940) > 5e-4


@test
def one_loop_beta_coefficients():
    """Textbook one-loop SM coefficients: b = (41/10, -19/6, -7) for g1 (SU(5)
    normalised), g2, g3; d lambda/d ln mu^2 has 12 lambda^2 - 3 y_t^4."""
    X = np.array([0.2, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0])
    b = rge.betas(X, loops=1) / rge.K
    assert np.allclose(b[:3] / X[:3] ** 2, [41 / 10, -19 / 6, -7.0])
    X = np.array([0.0, 0.0, 0.0, 0.81, 0.0, 0.0, 0.1, 0.0])
    b = rge.betas(X, loops=1) / rge.K
    assert abs(b[6] - (0.1 * (12 * 0.1 + 6 * 0.81) - 3 * 0.81**2)) < 1e-12


@test
def running_reproduces_buttazzo_planck_values():
    """Their eq. (4.1): couplings at M_Pl = 1.22e19 GeV for their inputs."""
    tr = st.trajectory(_buttazzo_inputs(), mu_max=1.3e19)
    L = st.landmarks(tr)
    assert abs(L["lambda_MPl"] - (-0.0143)) < 1e-4, L["lambda_MPl"]
    assert abs(L["yt_MPl"] - 0.3825) < 1e-4
    assert abs(L["g2_MPl"] - 0.5055) < 1e-4
    assert abs(L["g3_MPl"] - 0.4873) < 1e-4
    assert abs(L["g1_MPl"] - 0.6154) < 2e-4
    # instability scales: log10 Lambda_V = 9.5, Lambda_lambda ~ 2 Lambda_V, Lambda_I ~ 13 Lambda_V
    assert abs(math.log10(L["Lambda_V"]) - 9.5) < 0.05
    assert 1.7 < L["Lambda_lambda"] / L["Lambda_V"] < 2.3
    assert 10 < L["Lambda_I"] / L["Lambda_V"] < 16


@test
def rk4_matches_adaptive_integrator():
    inp = _buttazzo_inputs()
    X0 = rge.state(rge.match(inp))
    mu, X = rge.run(X0, inp.Mt, 1e19, loops=3)
    sol = solve_ivp(lambda s, x: rge.betas(x, 3), (2 * math.log(inp.Mt), 2 * math.log(1e19)), X0,
                    method="DOP853", rtol=1e-12, atol=1e-14)
    assert np.max(np.abs(sol.y[:, -1] - X[-1])) < 1e-10


@test
def light_cone_volume_matter_era():
    """Flat matter-only universe: H0^4 x past-light-cone volume = 16 pi / 1485
    (Buttazzo et al. eq. 6.5)."""
    v = st.light_cone_volume(om_l=0.0, om_r=1e-12)
    assert abs(v - 16 * math.pi / 1485) < 2e-5, v


@test
def fast_toy_statistic_matches_full_fit():
    """The vectorised Newton q0 of the pseudo-experiments equals a direct
    maximisation of the same likelihood (one channel, fixed nuisances)."""
    import higgsfit as hf
    rng = np.random.default_rng(3)
    grid = hf.GRID
    fb = np.ones_like(grid) / (grid[-1] - grid[0])
    fs = np.exp(-0.5 * ((grid - 125) / 1.6) ** 2)
    fs /= np.trapezoid(fs, grid)
    m = np.concatenate([rng.uniform(grid[0], grid[-1], 12), rng.normal(125, 1.6, 6)])
    S, B = 5.0, 12.0
    r = (np.interp(m, grid, fs) * S)[None, :]
    b = np.interp(m, grid, fb) * B
    q_fast = hf._q0_newton([(S, r, b)])[0]
    mus = np.linspace(0, 6, 60001)
    ll = np.array([-mu * S + np.sum(np.log1p(mu * r[0] / b)) for mu in mus])
    assert abs(q_fast - 2 * ll.max()) < 1e-6


@test
def afb_quark_level_signs():
    import dimuon
    a = dimuon.afb_quark_level(np.array([70.0, 150.0]))
    assert a[0] < -0.3 and a[1] > 0.3
    m = np.linspace(80, 95, 15001)
    z = m[np.argmin(np.abs(dimuon.afb_quark_level(m)))]
    assert 88.5 < z < 90.0


@test
def collins_soper_axis():
    """A pair moving slowly along +z whose negative lepton goes along +z in the
    pair frame has cos(theta*) = +1."""
    m, pz = 91.0, 1.0
    E = math.hypot(m, pz)
    gamma, beta = E / m, pz / E
    # in the pair frame: lepton (+m/2 along z), antilepton (-m/2 along z); boost along z
    lep = np.array([[gamma * (m / 2 + beta * m / 2)], [0.0], [0.0], [gamma * (m / 2 + beta * m / 2)]])
    alep = np.array([[gamma * (m / 2 - beta * m / 2)], [0.0], [0.0], [gamma * (-m / 2 + beta * m / 2)]])
    c = events.collins_soper_cos(lep, alep)
    assert abs(c[0] - 1.0) < 1e-9, c


def main() -> int:
    failed = 0
    for t in TESTS:
        try:
            t()
            print(f"  ok    {t.__name__}")
        except AssertionError as err:
            failed += 1
            print(f"  FAIL  {t.__name__}: {err}")
    print(f"{len(TESTS) - failed}/{len(TESTS)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
