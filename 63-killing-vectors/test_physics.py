"""Self-checks for project N6:  python test_physics.py

Each test compares a piece of our code with an independent route or a textbook
number. They run in about two minutes (the Kerr algebra dominates).
"""

from __future__ import annotations

import random
import sys
import time
from fractions import Fraction

import numpy as np
import sympy as sp

import config as C
import geometry as geo
import jets
import killing as kl
import orbits as ob

R = sp.Rational
RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append((name, bool(ok)))
    print(f"[{'ok' if ok else 'FAIL'}] {name} {detail}")


def test_lie_vs_covariant():
    """Killing's equation two ways: Lie derivative and symmetrised covariant derivative."""
    S = geo.schwarzschild()
    f = [geo.r**2 * geo.x, geo.t * geo.x, geo.r * geo.phi, geo.x**3]
    d = (S.lie_derivative(f) - S.killing_operator(f)).applyfunc(lambda e: sp.cancel(sp.together(e)))
    check("Lie derivative = nabla_a xi_b + nabla_b xi_a", d == sp.zeros(4))


def test_taylor():
    e = (geo.r**2 + 3 * geo.x) / (geo.r - 2 * geo.x**2 + 5)
    ser = jets.taylor_series(e, (geo.r, geo.x), (R(1, 2), R(1, 3)), 4)
    y0, y1 = sp.symbols("y0 y1")
    ref = sp.series(sp.series(e.subs({geo.r: R(1, 2) + y0, geo.x: R(1, 3) + y1}), y0, 0, 5).removeO(), y1, 0, 5).removeO()
    refp = sp.Poly(sp.expand(ref), y0, y1).as_dict()
    ok = all(abs(Fraction(int(sp.Rational(refp.get(k_, 0)).p), int(sp.Rational(refp.get(k_, 0)).q)) - v_) == 0
             for k_, v_ in ser.items())
    check("exact Taylor series of a rational function", ok)


def test_jets():
    b = jets.killing_vector_bound(geo.minkowski(), (R(1, 3), R(1, 5), R(2, 7), R(1, 2)), 3)["bound"]
    check("flat spacetime: 10 Killing vectors allowed", b == 10, f"(bound {b})")
    b4 = jets.killing_vector_bound(geo.schwarzschild(), (0, R(7, 2), R(1, 3), 0), 4, {geo.M: 1})["bound"]
    check("Schwarzschild: bound 4 at order 4", b4 == 4, f"(bound {b4})")
    bF = jets.killing_vector_bound(geo.minkowski_spherical(), (0, 3, R(1, 3), 0), 4)["bound"]
    check("flat spacetime in spherical coordinates: still 10", bF == 10)


def test_explicit():
    sol = kl.solve_polynomial(kl.chart_minkowski(), 2)
    sp_ = kl.same_span(sol["basis"], kl.minkowski_named(), (geo.t, geo.x, geo.y, geo.z))
    check("Minkowski: polynomial search finds exactly the 10 Poincare generators", sp_["equal"])
    S = geo.schwarzschild_angles()
    v = kl.verify(S, {"d_t": [1, 0, 0, 0], **kl.rotations_spherical()})
    check("Schwarzschild: d_t and the three rotations satisfy L_xi g = 0", all(v.values()))
    check("Schwarzschild: r d_r is not a Killing vector", not kl.verify(S, {"rd": [0, geo.r, 0, 0]})["rd"])
    F = geo.flrw()
    a_t = sp.Function("a")(geo.t)
    psi = F.conformal_factor([a_t, 0, 0, 0])
    check("FLRW: a(t) d_t is a conformal Killing vector with factor a'", sp.simplify(psi - sp.diff(a_t, geo.t)) == 0)


def test_kerr_hidden_symmetry():
    K_ = geo.kerr()
    Y = kl.kerr_killing_yano()
    check("Kerr: Penrose-Floyd tensor is Killing-Yano", kl.killing_yano_check(K_, Y))
    Kc = kl.carter_tensor(K_, Y)
    check("Kerr: Carter's K = Y Y is a Killing tensor", kl.killing_tensor_check(K_, Kc))
    check("Kerr is Ricci-flat (vacuum)", all(sp.cancel(sp.together(K_.ricci[i, j])) == 0
                                             for i in range(4) for j in range(4)))


def test_special_radii():
    sch = ob.schwarzschild_special_radii()
    check("Schwarzschild ISCO = 6, photon sphere = 3, b = sqrt(27)",
          abs(sch["r_isco"] - 6) < 1e-6 and abs(sch["r_photon"] - 3) < 1e-6 and abs(sch["b_crit"] - 27**0.5) < 1e-9)
    check("Schwarzschild efficiency 1 - sqrt(8/9)", abs(sch["efficiency"] - C.PUBLISHED["efficiency_schwarzschild"]) < 1e-9)
    dev = 0
    for a in (0.1, 0.5, 0.9, 0.99, -0.5, -0.99):
        d = ob.isco_kerr_derived(a)
        dev = max(dev, abs(d["r_isco"] - float(ob.isco_bpt(a))), abs(d["E"] ** 2 - (1 - 2 / (3 * d["r_isco"]))),
                  abs((d["L"] - a * d["E"]) ** 2 - d["r_isco"] ** 2 / 3),
                  abs(ob.photon_orbit_derived(a)["r"] - float(ob.photon_orbit_bpt(a))))
    check("Kerr ISCO and photon orbits from R(r) = Bardeen-Press-Teukolsky", dev < 1e-9, f"(max dev {dev:.1e})")
    check("extremal Kerr efficiencies 42.26% / 3.775%",
          abs(ob.isco_kerr_derived(1.0)["efficiency"] - C.PUBLISHED["efficiency_extremal_prograde"]) < 1e-12 and
          abs(ob.isco_kerr_derived(-0.999999999)["efficiency"] - C.PUBLISHED["efficiency_extremal_retrograde"]) < 1e-5)
    E, L = ob.circular_from_R(10.0, 0.7)
    E0, L0 = ob.circular_kerr(10.0, 0.7)
    check("circular-orbit E, L from R = R' = 0 match BPT", abs(E - E0) < 1e-12 and abs(L - L0) < 1e-12)


def test_precession():
    for p, e in [(30.0, 0.4), (5541.0, 0.8859)]:
        E, L = ob.bound_orbit_EL(p, e)
        a_ = ob.precession_exact(p, e)
        b_ = ob.precession_by_quadrature(E, L)
        check(f"precession exact = quadrature (p = {p:.0f})", abs(a_ - b_) < 1e-11 * max(1, a_))
    check("precession -> 6 pi / p for wide orbits", abs(ob.precession_exact(1e6, 0.5) / (6 * np.pi / 1e6) - 1) < 1e-5)


def test_noether_quick():
    S = geo.schwarzschild()
    tau, X, U = ob.integrate_geodesic(S, {geo.M: 1}, np.array([0, 10.0, 0.0, 0]),
                                      np.array([1.2, 0.0, -0.02, 0.03]), 200.0, n_out=50)
    c = ob.conserved_along(S, {geo.M: 1}, X, U, {"T": [1, 0, 0, 0], "Phi": [0, 0, 0, 1]})
    check("Noether: E and L_z constant along an integrated geodesic",
          np.ptp(c["T"]) < 1e-9 and np.ptp(c["Phi"]) < 1e-9)


def test_shadow():
    import shadows as sh
    d0 = sh.mean_diameter(*sh.critical_curve(0.0, 45))["mean_diameter"]
    check("Schwarzschild shadow diameter 2 sqrt(27)", abs(d0 - 2 * 27**0.5) < 1e-6)
    d1 = sh.mean_diameter(*sh.critical_curve(0.998, 1))["mean_diameter"]
    d2 = sh.mean_diameter(*sh.critical_curve(0.998, 90))["mean_diameter"]
    check("near-extremal Kerr: face-on 9.66, edge-on 9.86 (EHT range 9.6-10.4)", 9.6 < d1 < 9.7 and 9.8 < d2 < 9.9,
          f"({d1:.3f}, {d2:.3f})")


def test_thermo():
    import thermo as th
    t = th.kerr_thermo(1.0, 0.0)
    check("Hawking T for 1 Msun = 6.1703e-8 K", abs(t["T_K"] / 6.1703e-8 - 1) < 2e-4)
    check("entropy of 1 Msun = 1.0489e77 k_B", abs(t["S_kB"] / 1.0489e77 - 1) < 2e-4)
    check("Page 1976 lifetime 2.16e66 yr", abs(th.evaporation_time_yr(1.0, 2.011e-4) / 2.16e66 - 1) < 0.01)
    k = th.kerr_thermo(1.0, 1.0 - 1e-15)
    check("extremal Kerr: T -> 0, M_irr -> M/sqrt(2)", k["T_K"] < 1e-12 and abs(k["M_irr_frac"] - 2**-0.5) < 1e-6)


def test_cosmology():
    import cosmology as cs
    import json
    f = cs.tcmb_fit()
    d7 = json.loads(C.D7_RESULTS.read_text())["cmb_temperature"]["baseline"] if C.D7_RESULTS.exists() else None
    if d7:
        check("CMB beta reproduces D7 exactly", abs(f["beta"] - d7["beta"]) < 1e-9 and abs(f["chi2_min"] - d7["chi2"]) < 1e-6)
    g = cs.geodesics_demo()
    check("photon a E constant through LCDM (conformal Killing vector)", g["photon"]["aE_spread"] < 1e-6)
    check("massive a^2 dx/dlambda constant (spatial Killing vector)", g["massive"]["a2dx_spread"] < 1e-9)


def test_s2_link():
    import h1link
    import s2
    m = h1link.model()
    ch = h1link.chain()
    theta = np.array([ch["_best"][n_] for n_ in m.free])
    c2 = m.chi2(theta)
    check("H1's model reproduces H1's best-fit chi^2 (Keck 94.8, VLT 289.3, RV 84.7)",
          abs(c2["astro_keck"] - 94.76) < 0.05 and abs(c2["astro_vlt"] - 289.31) < 0.05 and abs(c2["rv"] - 84.72) < 0.05)
    op = h1link.modules().op
    b = ch["_best"]
    X, Y, Z, VX, VY, VZ, rr = op.orbit_state(np.array([2010.0]), b["gm6"], b["period"], b["t0"], b["e"], b["inc"],
                                             b["omega"], b["Omega"], f_sp=1.0)
    h_state = ((X * VY - Y * VX) / op.KMS_PER_AUYR / b["r0"] ** 2)[0]
    gm = b["gm6"] * 1e6 * op.GM_SUN_AU3YR2
    a_ = op.semi_major_axis_au(b["gm6"], b["period"])
    check("sky angular-momentum sign convention", abs(h_state - s2._h_model(gm, a_, b["e"], b["inc"], b["r0"], b["gm6"])) < 1e-6)
    g = s2.newton_vs_gr_along_orbit()
    check("Killing E and L constant along S2's Schwarzschild orbit", g["E_killing_spread"] < 1e-12 and
          g["L_killing_spread"] < 1e-9)


def test_gw250114():
    import thermo as th
    gw = th.gw250114()
    p = gw["paper"]
    ok = (abs(gw["M_final"]["median"] - p["M_final"][0]) < 0.06 and abs(gw["chi_final"]["median"] - p["chi_final"][0]) < 0.006
          and abs(gw["m1"]["median"] - p["m1"][0]) < 0.06 and abs(gw["m2"]["median"] - p["m2"][0]) < 0.06)
    check("GW250114 samples reproduce the paper's medians", ok)


if __name__ == "__main__":
    t0 = time.time()
    for fn in (test_lie_vs_covariant, test_taylor, test_jets, test_explicit, test_special_radii, test_precession,
               test_noether_quick, test_shadow, test_thermo, test_cosmology, test_s2_link, test_gw250114,
               test_kerr_hidden_symmetry):
        fn()
    n_fail = sum(not ok for _, ok in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} passed in {time.time() - t0:.0f} s")
    sys.exit(1 if n_fail else 0)
