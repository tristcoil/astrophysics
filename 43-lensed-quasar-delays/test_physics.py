"""Checks of the physics and the numerics. Run: python test_physics.py (or pytest).

The synthetic light curves below are *tests of the code* (does the estimator
return a delay we planted?); no synthetic data enters any result.
"""

from __future__ import annotations

import numpy as np

import config as C


def test_distances_match_astropy():
    from astropy.cosmology import FlatLambdaCDM

    import h0 as H

    for zd, zs, om, hh in [(0.295, 0.654, 0.3, 70.0), (0.745, 1.789, 0.2, 65.0), (0.6304, 1.394, 0.45, 80.0)]:
        cos = FlatLambdaCDM(H0=hh, Om0=om)
        dd, ddt = H.distances(zd, zs, om, hh)
        Dd = cos.angular_diameter_distance(zd).value
        Ds = cos.angular_diameter_distance(zs).value
        Dds = cos.angular_diameter_distance(zd, zs).value
        assert abs(dd[0] / Dd - 1) < 1e-6
        assert abs(ddt[0] / ((1 + zd) * Dd * Ds / Dds) - 1) < 1e-6


def test_epl_identities():
    import lensmodel as LM

    m = LM.EPLShear(1.2, 0.75, 0.4, 2.15, 0.05, -0.03, 0.03, -0.02)
    h = 1e-5
    for x, y in [(0.7, -0.9), (-1.3, 0.2), (0.1, 1.6)]:
        ax, ay = m.deflection(x, y)
        px = (m.potential(x + h, y) - m.potential(x - h, y)) / (2 * h)
        py = (m.potential(x, y + h) - m.potential(x, y - h)) / (2 * h)
        assert abs(px - ax) < 1e-7 and abs(py - ay) < 1e-7  # alpha = grad psi
        div = ((m.deflection(x + h, y)[0] - m.deflection(x - h, y)[0])
               + (m.deflection(x, y + h)[1] - m.deflection(x, y - h)[1])) / (2 * h)
        # the shear is traceless, so div(alpha)/2 is the EPL's own convergence
        assert abs(div / 2 - m.convergence(x, y)) < 1e-5
    sis = LM.EPLShear(1.0, 1.0, 0.0, 2.0, 0.0, 0.0)
    assert abs(np.hypot(*sis.deflection(2.3, 0.0)) - 1.0) < 1e-12  # SIS: |alpha| = theta_E


def test_mass_sheet_keeps_images_scales_delays():
    import lensmodel as LM

    img = LM.he0435_positions()
    m, _ = LM.fit_positions(img, 2.0)
    d1 = LM.fermat_differences(m, img)
    for lam in (0.8, 1.2):
        mm = m.mst(lam)
        xs = np.array([v[0] for v in img.values()])
        ys = np.array([v[1] for v in img.values()])
        sx0, sy0 = m.source(xs, ys)
        sx, sy = mm.source(xs, ys)
        # the four images still map to one (rescaled) source position
        assert np.allclose(sx, lam * sx0, atol=1e-12) and np.allclose(sy, lam * sy0, atol=1e-12)
        dl = LM.fermat_differences(mm, img)
        for k in d1:
            assert abs(dl[k] / d1[k] - lam) < 1e-3 * max(1, abs(1 / d1[k]))


def test_witt_formula_isothermal():
    """Isothermal lens, any angular structure, no shear: tau_i - tau_j = (r_j^2 - r_i^2)/2 (Witt et al. 2000)."""
    import lensmodel as LM

    m = LM.EPLShear(1.1, 0.7, 0.3, 2.0, 0.0, 0.0)
    # images of a source: find stationary points by solving the lens equation on a grid + polish
    from scipy.optimize import root

    bx, by = 0.05, 0.08
    sols = []
    for r0 in np.linspace(0.5, 1.8, 8):
        for th in np.linspace(0, 2 * np.pi, 24, endpoint=False):
            s = root(lambda p: np.array(m.source(p[0], p[1])) - (bx, by), [r0 * np.cos(th), r0 * np.sin(th)])
            if s.success and np.hypot(*s.x) > 0.2 and not any(np.hypot(*(s.x - q)) < 1e-4 for q in sols):
                sols.append(s.x)
    assert len(sols) >= 4
    tau = [m.fermat(x, y, bx, by) for x, y in sols]
    r2 = [x * x + y * y for x, y in sols]
    for i in range(len(sols)):
        for j in range(len(sols)):
            assert abs((tau[i] - tau[j]) - 0.5 * (r2[j] - r2[i])) < 1e-6


def test_solver_matches_dense():
    """The banded Schur-complement solve equals a direct dense solve."""
    import delays as D
    import lightcurves as L

    lc = L.read_lens("HE0435")
    m = D.JointSpline(lc, "A", D.SplineSettings(knot_int=25, knot_ml=150), objective="chi2")
    taus = {"B": 9.0, "C": 0.8, "D": 13.8}
    tot, coef, _ = m.solve(taus)
    A = [m.design(im, 0.0 if im == "A" else taus[im]) for im in lc.images]
    N = m.P.copy()
    b = np.zeros(m.ncol)
    yWy = 0.0
    for im, Ai in zip(lc.images, A):
        Ad = Ai.toarray()
        N += Ad.T @ (Ad * m.w[im][:, None])
        b += Ad.T @ (m.w[im] * lc.mag[im])
        yWy += float(lc.mag[im] @ (m.w[im] * lc.mag[im]))
    c2 = np.linalg.solve(N, b)
    assert np.allclose(coef, c2, rtol=1e-7, atol=1e-6)  # normal matrix condition number ~3e10
    assert abs(tot - (yWy - b @ c2)) < 1e-6 * tot
    me = D.JointSpline(lc, "A", D.SplineSettings(knot_int=25, knot_ml=150), objective="evidence")
    tot_e = me.chi2(taus)
    sign, logdet = np.linalg.slogdet(N)
    assert sign > 0 and abs(tot_e - (tot + logdet)) < 1e-6 * abs(tot_e)


def test_joint_spline_recovers_planted_delays():
    """Plant delays in curves built from a real image's shape; the estimator must find them."""
    import delays as D
    import lightcurves as L

    lc = L.read_lens("HE0435")
    rng = np.random.default_rng(7)
    m0 = D.JointSpline(lc, "A", D.SplineSettings(knot_int=25, knot_ml=150))
    fit = m0.fit({"B": 9.0, "C": 0.8, "D": 13.8})
    true = {"B": 11.3, "C": -2.4, "D": 17.1}
    y = {}
    for im in lc.images:
        tau = 0.0 if im == "A" else true[im]
        y[im] = m0.intrinsic(fit.coef, lc.t - tau) + 0.02 * np.sin(lc.t / 900.0 + len(im)) + rng.normal(0, 0.01, lc.n)
    lc2 = L.LightCurves("test", "test", lc.images, lc.t, y, {k: np.full(lc.n, 0.01) for k in lc.images}, lc.tel)
    m = D.JointSpline(lc2, "A", D.SplineSettings(knot_int=25, knot_ml=150))
    est = m.optimize({"B": 9.0, "C": 0.8, "D": 13.8}, width=8)
    for k in true:
        assert abs(est[k] - true[k]) < 0.5, (k, est[k], true[k])


def test_iccf_finds_a_shift():
    import ccf

    t = np.sort(np.random.default_rng(1).uniform(0, 600, 300))
    f = lambda x: np.sin(x / 23.0) + 0.5 * np.sin(x / 7.0)
    r, _ = ccf.iccf(t, f(t), t, f(t - 12.0), np.arange(0, 25, 0.25))
    ok, rmax, tp, tc = ccf.peak_and_centroid(np.arange(0, 25, 0.25), r)
    assert ok and abs(tp - 12.0) <= 0.25 and abs(tc - 12.0) < 0.5


def test_delay_scale_and_swap():
    import h0 as H

    pub = {("A", "B"): 10.0, ("A", "C"): 20.0}
    ours = {k: 1.05 * v for k, v in pub.items()}
    sc = H.delay_scale(ours, np.diag([0.25, 0.25]), pub, {k: 0.5 for k in pub})
    assert abs(sc["r"] - 1.05) < 1e-12
    rng = np.random.default_rng(2)
    lik = H.LensLikelihood("x", 0.5, 1.5, "1d", ddt=rng.lognormal(np.log(3000), 0.06, 20000)).build()
    same = H.swap_delays(lik, 1.0, 0.03, 0.03)
    assert np.allclose(same.ddt, lik.ddt)
    up = H.swap_delays(lik, 1.05, 0.03, 0.03)
    assert abs(np.median(up.ddt) / np.median(lik.ddt) - 1.05) < 1e-9


def test_published_delay_conventions():
    """Sign conversions in config: arrival order must match the physics of each lens."""
    L = C.LENSES
    assert L["RXJ1131"]["model_delays"][("B", "D")][0] > 85  # D arrives last, ~3 months late
    assert L["HE0435"]["model_delays"][("A", "D")][0] > 0  # A first
    assert L["PG1115"]["model_delays"][("A", "C")][0] < 0  # C leads A
    assert L["WFI2033"]["model_delays"][("B", "C")][0] > 0  # B first
    for key, v in L.items():
        for pairs in (v.get("model_delays") or {}, v.get("cosmograil19") or {}):
            for (a, b), (val, em, ep) in pairs.items():
                assert em > 0 and ep > 0


if __name__ == "__main__":
    import sys

    fails = 0
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print("ok  ", name)
            except AssertionError as exc:
                fails += 1
                print("FAIL", name, exc)
    sys.exit(1 if fails else 0)
