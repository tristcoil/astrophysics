"""Self-checks: our Friedmann solver and likelihoods against closed forms, CAMB,
independent methods, and DESI's own published numbers.

    python test_physics.py

Each check prints PASS/FAIL with the measured deviation and its tolerance, and
the script exits non-zero if anything fails.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
from scipy import integrate

import background as bg
import config as cfg
import fetch_data as fd
import friedmann as fr
import likelihoods as lk

RESULTS = []


def check(name: str, value: float, tol: float) -> None:
    ok = bool(np.isfinite(value) and abs(value) <= tol)
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<66} {value: .2e}  (tol {tol:.0e})")


U = fr.Universe


def test_closed_forms() -> None:
    print("Ages and lifetimes against closed-form solutions (units of 1/H0)")
    check("Einstein-de Sitter: t0 = 2/3", U(1, 0)._t_from_zero(1.0) - 2 / 3, 1e-12)
    check("Milne (empty): t0 = 1", U(0, 0)._t_from_zero(1.0) - 1.0, 1e-12)
    check("radiation only: t0 = 1/2", U(0, 0, Or=1.0)._t_from_zero(1.0) - 0.5, 1e-12)
    check("de Sitter: no beginning (t0 infinite)", 0.0 if np.isinf(U(0, 1)._t_from_zero(1.0))
          else 1.0, 0.0)
    for om in (0.05, 0.3, 0.9):
        check(f"flat LCDM, Om = {om}", U(om, 1 - om)._t_from_zero(1.0) - fr.age_flat_lcdm(om), 1e-11)
    for om in (0.3, 0.7):
        check(f"open matter only, Om = {om}", U(om, 0)._t_from_zero(1.0) - fr.age_open_matter(om),
              1e-11)
    for om in (1.5, 3.0):
        u = U(om, 0)
        check(f"closed matter only, Om = {om}: age", u._t_from_zero(1.0) - fr.age_closed_matter(om),
              1e-11)
        check(f"closed matter only, Om = {om}: Big Bang to Big Crunch",
              u.fate()["lifetime_gyr"] / u.t_hubble_gyr - fr.lifetime_closed_matter(om), 1e-10)
    check("flat radiation + matter, t(a = 0.5)",
          U(0.7, 0, Or=0.3)._t_from_zero(0.5) - fr.t_radiation_matter(0.5, 0.3, 0.7), 1e-12)
    # pure phantom energy, flat: t_rip - t0 = 2 / (3 |1 + w|) exactly
    for w in (-1.2, -1.5, -2.0):
        u = U(0.0, 1.0, w=w)
        check(f"Big Rip with phantom energy only, w = {w}",
              u.time_to_rip_gyr() / u.t_hubble_gyr - 2 / (3 * abs(1 + w)), 1e-9)
    u = U(0.3, 0.7, w=-1.5, H0=70.0)
    check("Big Rip, w = -1.5, H0 = 70: CKW's '22 Gyr' (matter makes it 3% sooner)",
          u.time_to_rip_gyr() / 22.0 - 1.0, 0.05)


def test_boundaries() -> None:
    print("The geography of the (Omega_m, Omega_Lambda) plane")
    om = np.array([0.05, 0.1, 0.3, 0.5, 0.8, 1.5, 2.5])
    x = (1 - om) / om
    cpt = np.where(om < 0.5, 4 * om * np.cosh(np.arccosh(np.maximum(x, 1)) / 3) ** 3,
                   4 * om * np.cos(np.arccos(np.clip(x, -1, 1)) / 3) ** 3)
    check("no Big Bang edge vs Carroll, Press & Turner 1992 eq. 12",
          np.max(np.abs(fr.no_big_bang_boundary(om) - cpt)), 1e-6)
    om2 = np.array([1.1, 1.2, 2.0, 3.0, 5.0])
    rc = 4 * om2 * np.cos(np.arccos((1 - om2) / om2) / 3 + 4 * np.pi / 3) ** 3
    check("recollapse edge vs Carroll, Press & Turner 1992 eq. 11",
          np.max(np.abs(fr.recollapse_boundary(om2) - rc)), 1e-6)
    Or = fr.OR_GALLERY
    a = np.array([0.2, 0.5, 0.8, 1.5, 3.0])
    Om, OL = fr.double_root_curve(a, Or)
    for ai, o, l in zip(a, Om, OL):
        u = U(o, l, Or=Or)
        check(f"double root at a* = {ai}: P(a*) and P'(a*) vanish",
              max(abs(u.P(ai)), abs(u.dP(ai))), 1e-12)
    check("with radiation the edge meets Om = 0 at OL = (1 + sqrt(Or))^2",
          fr.no_big_bang_boundary([0.0], Or)[0] - (1 + np.sqrt(Or)) ** 2, 1e-9)
    check("smallest static size a* = (1 + Or^-1/2)^-1/2",
          fr.smallest_physical_static_a(Or) - (1 + Or**-0.5) ** -0.5, 1e-15)
    # the same edge for general w, by the generic root finder
    Om_w, OL_w = fr.double_root_curve(np.array([0.4]), Or, w=-0.8)
    u = U(float(Om_w[0]), float(OL_w[0]), Or=Or, w=-0.8)
    check("double root for w = -0.8: P(0.4) and P'(0.4) vanish",
          max(abs(u.P(0.4)), abs(u.dP(0.4))), 1e-12)
    # vectorised classification and bounce sizes vs one-by-one analysis
    rng = np.random.default_rng(3)
    om, ol = rng.uniform(0, 3, 400), rng.uniform(-1, 3, 400)
    codes = fr.classify_grid(om, ol, Or)
    a_min = fr.bounce_size_grid(om, ol, Or)
    wrong, da = 0, []
    for o, l, c, am in zip(om, ol, codes, a_min):
        f = U(o, l, Or=Or).fate()
        exp = 3 if "a_min" in f else (2 if f["future"] == "Big Crunch" else (0 if f["q0"] < 0 else 1))
        wrong += exp != c
        if "a_min" in f:
            da.append(abs(f["a_min"] - am))
    check("vectorised fate classification: mismatches out of 400", wrong, 0)
    check("vectorised bounce size vs polynomial roots", max(da), 1e-9)


def test_ode_and_grids() -> None:
    print("Independent methods: the acceleration equation, vectorised ages")
    for name, u in (("our universe", U(0.31, 0.69, Or=fr.OR_GALLERY)), ("closed", U(3.0, 0.0)),
                    ("bounce", U(0.3, 1.8)), ("negative Lambda", U(0.3, -0.5)),
                    ("phantom", U(0.3, 0.7, w=-1.5))):
        e = []
        for span in ((0, -3), (0, 3)):
            sol = fr.integrate_ode(u, span)
            for tg, a in u.history(-3 * u.t_hubble_gyr, 3 * u.t_hubble_gyr, n=400):
                th = tg / u.t_hubble_gyr
                msk = (th >= sol.t.min()) & (th <= sol.t.max()) & (a > 2e-3) & (a < 50)
                if msk.any():
                    e.append(np.max(np.abs(sol.sol(th[msk])[0] - a[msk]) / np.maximum(a[msk], 1)))
        check(f"a(t) from a'' = (aP' - 2P)/(2a^3) vs quadrature: {name}", max(e), 1e-7)
    rng = np.random.default_rng(4)
    om, ol = rng.uniform(0.05, 2.0, 60), rng.uniform(-0.5, 1.0, 60)
    t_grid = fr.age_grid(om, ol, fr.OR_GALLERY)
    err = [abs(t_grid[i] - U(om[i], ol[i], Or=fr.OR_GALLERY)._t_from_zero(1.0))
           for i in range(om.size)]
    check("vectorised Gauss-Legendre ages vs adaptive quadrature", max(err), 1e-8)


def test_background_vs_camb() -> None:
    import camb

    print("Physical-density background vs CAMB 2.0.4, curvature included")
    rec = bg.Recombination()
    zs = np.array([0.295, 0.51, 0.706, 0.934, 1.321, 1.484, 2.33])
    for H0, ob, oc, ok in ((67.36, 0.02237, 0.1200, 0.0), (55.0, 0.0222, 0.118, -0.2),
                           (75.0, 0.0225, 0.121, 0.03), (70.0, 0.0221, 0.125, 0.2)):
        res = camb.get_background(bg.camb_params(ob, oc, H0=H0, omk=ok))
        d = res.get_derived_params()
        m = bg.Model(h=H0 / 100, ombh2=ob, omch2=oc, omk=ok)
        tag = f"H0 = {H0}, Omega_k = {ok:+.2f}"
        check(f"theta_* ({tag})", rec.theta_star(m)[0] / (d["thetastar"] / 100) - 1, 2e-5)
        check(f"r_drag ({tag})", rec.r_drag(m)[0] / d["rdrag"] - 1, 2e-5)
        dm = bg.transverse_distance(bg.comoving_distance(zs, m), m)[0]
        check(f"D_M at the BAO redshifts ({tag})",
              np.max(np.abs(dm / (res.angular_diameter_distance(zs) * (1 + zs)) - 1)), 1e-6)
        check(f"H(z) at the BAO redshifts ({tag})",
              np.max(np.abs(bg.hubble(zs, m)[0] / res.hubble_parameter(zs) - 1)), 1e-6)
        t_camb = integrate.quad(lambda xx: 1.0 / res.hubble_parameter(np.expm1(xx)), 0,
                                np.log(1e8), epsrel=1e-12, epsabs=0, limit=500)[0] * \
            cfg.HUBBLE_TIME_GYR_X_H0
        check(f"age vs a tight integral of CAMB's H(z), Gyr ({tag})",
              bg.age_gyr(m)[0] - t_camb, 1e-5)
    m = bg.Model(h=np.array([0.68, 0.68]), ombh2=0.0224, omch2=0.12, omk=np.array([0.0, 1e-12]))
    sd = bg.SupernovaDistances(np.linspace(0.01, 2.3, 50))
    dc_gl = bg.comoving_distance(np.linspace(0.01, 2.3, 50), m)
    check("supernova distances (trapezoid) vs Gauss-Legendre",
          np.max(np.abs(sd.comoving(m) / dc_gl - 1)), 2e-6)
    m2 = bg.Model(h=0.68, ombh2=0.0224, omch2=0.12, omk=np.array([1e-10, 2e-9]))
    dc = np.array([[3000.0], [3000.0]])
    exact = bg.C_H100 / 0.68 / np.sqrt(m2.omk) * np.sinh(np.sqrt(m2.omk) * 3000.0 * 0.68 / bg.C_H100)
    check("D_M: small-curvature series vs sinh",
          np.max(np.abs(bg.transverse_distance(dc, m2)[:, 0] / exact - 1)), 1e-9)


def test_likelihoods() -> None:
    print("Likelihoods: DESI's own numbers, D5's prior, the Laplace approximation")
    bao = lk.BAOLikelihood(fd.load_bao())
    try:
        bf = fd.load_desi_bestfit("desi-bao-all")
        mb = bg.Model.from_omk(bf["omm"], bf["omk"], 0.6736, 0.02237)
        check("BAO chi^2 at DESI's best fit (LCDM + curvature, BAO alone)",
              bao.chi2_calibrated(mb, [bf["hrdrag"] / 0.6736])[0] - bf["chi2__BAO"], 1e-4)
        bf = fd.load_desi_bestfit("desi-bao-all_pantheonplus")
        mb = bg.Model.from_omk(bf["omm"], bf["omk"], 0.6736, 0.02237)
        sn = lk.SupernovaLikelihood(fd.load_pantheon())
        check("Pantheon+ chi^2 at DESI's BAO + Pantheon+ best fit",
              sn.chi2(mb)[0][0] - bf["chi2__SN"], 1e-3)
        check("BAO chi^2 at DESI's BAO + Pantheon+ best fit",
              bao.chi2_calibrated(mb, [bf["hrdrag"] / 0.6736])[0] - bf["chi2__BAO"], 1e-4)
    except FileNotFoundError:
        print("  (DESI best fits not cached: skipped)")
    cmb = lk.cmb_prior_from_chains(fd.load_planck_early_chains())
    d5 = fd.load_d5_cmb_priors()
    if d5 is not None:
        check("compressed CMB prior rebuilt = D5's consistent prior (means, relative)",
              np.max(np.abs(cmb.mean / np.array(d5["consistent"]["mean"]) - 1)), 1e-12)
    rec = bg.Recombination()
    for use_bao in (False, True):
        marg = lk.PhysicalDensityMarginal(cmb, rec, bao if use_bao else None)
        om = np.array([0.30, 0.31, 0.32, 0.29])
        ok = np.array([0.000, 0.002, -0.004, 0.006])
        res = marg(om, ok)
        worst = 0.0
        for i in range(om.size):
            x0 = np.array([res["ombh2"][i], res["omega_bc"][i]])
            s = np.array([res["sig_ombh2"][i], res["sig_omega_bc"][i]])
            g1 = x0[0] + s[0] * np.linspace(-7, 7, 141)
            g2 = x0[1] + s[1] * np.linspace(-7, 7, 141)
            G1, G2 = np.meshgrid(g1, g2, indexing="ij")
            X = np.column_stack([G1.ravel(), G2.ravel()])
            r, _ = marg.predictions(np.full(X.shape[0], om[i]), np.full(X.shape[0], ok[i]), X)
            chi2 = np.einsum("nk,nk->n", r, r)
            L = np.exp(-0.5 * (chi2 - chi2.min())).reshape(G1.shape)
            integ = np.trapezoid(np.trapezoid(L, g2, axis=1), g1)
            brute = chi2.min() - 2 * np.log(integ) + 2 * np.log(2 * np.pi)
            worst = max(worst, abs(res["m2lnl"][i] - brute))
        check(f"Laplace marginal vs brute-force integral ({'BAO + ' if use_bao else ''}CMB)",
              worst, 0.02)


def main() -> int:
    warnings.simplefilter("ignore", category=RuntimeWarning)
    test_closed_forms()
    test_boundaries()
    test_ode_and_grids()
    test_background_vs_camb()
    test_likelihoods()
    n_fail = RESULTS.count(False)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
