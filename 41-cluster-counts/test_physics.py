"""Check the physics, the selection model and the data readers against answers known independently.

    python test_physics.py

Each check prints PASS or FAIL with the numbers compared; the script exits
non-zero if any check fails. The comparisons with colossus and CAMB need those
packages (both in the shared venv); the data checks need fetch_data.py to have run.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import catalog as K  # noqa: E402
import config as C  # noqa: E402
import cosmology as CO  # noqa: E402
import footprint as F  # noqa: E402
import model as MD  # noqa: E402
import selection as S  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


# --------------------------------------------------------------------------


def test_cosmology():
    print("Cosmology against colossus")
    from colossus.cosmology import cosmology
    from colossus.lss import mass_function

    c = CO.fiducial()
    col = cosmology.setCosmology("t", flat=True, H0=c.h * 100, Om0=c.Om, Ob0=c.Ob, sigma8=c.sigma8, ns=c.ns,
                                 Tcmb0=C.T_CMB, Neff=C.N_EFF, relspecies=True)
    z = np.array([0.1, 0.3, 0.5, 0.8])
    g = c.growth(z) / np.array([col.growthFactor(zz) for zz in z])
    check("growth factor D(z)/D(0) matches colossus to 1e-4", np.all(np.abs(g - 1) < 1e-4), f"max dev {np.abs(g - 1).max():.1e}")
    d = c.comoving_distance(z) / (col.comovingDistance(np.zeros(4), z) / c.h)
    check("comoving distance matches colossus to 1e-4", np.all(np.abs(d - 1) < 1e-4), f"max dev {np.abs(d - 1).max():.1e}")
    M = np.geomspace(1e13, 3e15, 8)
    R = (3 * M / (4 * np.pi * c.rho_m0)) ** (1 / 3)
    s, _ = c.sigma_M(M, np.array([0.0]))
    r = s[0] / col.sigma(R * c.h, 0)
    check("sigma(M) (Eisenstein & Hu 1998) matches colossus to 0.2%", np.all(np.abs(r - 1) < 2e-3),
          f"max dev {np.abs(r - 1).max():.1e}")
    # Tinker 2008 with colossus's convention (linear in Delta, rounded) must agree closely
    orig = CO.Cosmology.__dict__["tinker08_params"]  # the staticmethod object itself

    def lin(delta_m, zz):
        dm = np.round(delta_m)
        D = 10**CO._T08_DELTA
        A0, a0, b0, c0 = (np.interp(dm, D, t) for t in (CO._T08_A, CO._T08_a, CO._T08_b, CO._T08_c))
        alpha = 10 ** (-((0.75 / np.log10(dm / 75.0)) ** 1.2))
        zp = 1 + np.asarray(zz, float)
        return A0 * zp**-0.14, a0 * zp**-0.06, b0 * zp**-alpha, c0

    try:
        CO.Cosmology.tinker08_params = staticmethod(lin)
        Mm = np.geomspace(5e13, 1.5e15, 6)
        worst = 0.0
        for zz in (0.1, 0.4, 0.7):
            ours = c.dndlnM(Mm, zz)[0]
            theirs = mass_function.massFunction(Mm * c.h, zz, mdef="500c", model="tinker08", q_out="dndlnM") * c.h**3
            worst = max(worst, np.abs(ours / theirs - 1).max())
        check("Tinker 2008 M500c mass function matches colossus (same table interpolation) to 1%", worst < 0.01,
              f"max dev {worst:.1e}")
    finally:
        CO.Cosmology.tinker08_params = orig
    zz = np.array([0.5, 1.0, 3.0])
    r = c.growth_ode(zz) / c.growth(zz)
    check("growth equation with radiation agrees with the integral formula at z <= 3 (to 0.3%)",
          np.all(np.abs(r - 1) < 3e-3), f"{np.round(r, 5)}")
    eds = CO.Cosmology(Om=1 - CO.OMEGA_R_H2 / 0.38**2, h=0.38, Ode=0.0)
    zz = np.array([1.0, 5.0, 20.0])
    r = eds.growth_ode(zz) * (1 + zz)
    check("matter-only universe: growth proportional to the scale factor (to 1% at z <= 20)", np.all(np.abs(r - 1) < 0.01),
          f"{np.round(r, 4)}")


def test_camb_link():
    print("sigma_8 through CAMB")
    import cmb_link

    s8 = cmb_link.planck_check()
    check("Planck 2018 best fit gives its published sigma_8 = 0.8120 (to 0.002)", abs(s8 - 0.8120) < 0.002, f"{s8:.4f}")


def test_survey():
    print("Survey footprint and exposure")
    fp = F.build()
    P = K.load_primary()
    T, _, _, _ = F.conditions_at(fp, P.ra, P.dec)
    r = np.nanmedian(np.abs(np.log(T / P.exp.values)))
    check("exposure map reproduces the catalogue's per-cluster exposures (median |ln ratio| < 2%)", r < 0.02, f"{r:.4f}")
    check("primary footprint area within 0.5% of the published 13,116 deg^2",
          abs(fp["area_primary"] / C.AREA_PRIMARY_DEG2 - 1) < 0.005, f"{fp['area_primary']:.0f}")
    check("cosmology footprint area within 2% of the published 12,791 deg^2",
          abs(fp["area_cosmo"] / C.AREA_COSMO_DEG2 - 1) < 0.02, f"{fp['area_cosmo']:.0f}")
    clean, _ = F.digitize_fig2b()
    e = np.sort(fp["exposure"][fp["mask_primary"] & np.isfinite(fp["exposure"])])[::-1]
    worst = 0.0
    for tq in (100, 200, 500, 1000):
        ours = (e >= tq).sum() * F.PIX_AREA
        j = np.argmin(np.abs(np.log(clean[:, 0] / tq)))
        worst = max(worst, abs(ours / clean[j, 1] - 1))
    check("area deeper than 100-1000 s matches Bulbul et al. Fig. 2 to 12%", worst < 0.12, f"worst {worst:.3f}")


def test_selection():
    print("Selection function")
    m = S.SelectionModel()
    v = S.validate(m)
    worst = max(v["fig4_max_abs_dP"].values())
    check("model reproduces the nine published curves to |dP| < 0.015", worst < 0.015, f"{worst:.4f}")
    rat = np.array([r["ratio"] for r in v["clerc"].values()])
    check("model matches Clerc et al.'s 200 s and 2000 s curves to 6%", np.all(np.abs(rat - 1) < 0.06),
          f"{rat.min():.3f}-{rat.max():.3f}")
    p = m(np.array([100.0]), 0.3, 150.0)[0]
    check("bright clusters are almost always detected (P > 0.99 at 100 cts/s)", p > 0.99, f"{p:.4f}")
    lo, hi = m(np.array([0.3]), 0.3, 150.0, nh=1e20)[0], m(np.array([0.3]), 0.3, 150.0, nh=1e21)[0]
    check("more Galactic absorption makes detection harder", lo > hi, f"P = {lo:.3f} (1e20) vs {hi:.3f} (1e21)")


def test_catalogue():
    print("Catalogues")
    df = K.load_cosmology_sample()
    check("5,259 published clusters, 14 with zero count rate dropped", df.attrs["n_published"] == 5259 and len(df) == 5245,
          f"{df.attrs['n_published']} -> {len(df)}")
    check("purity weights give the published ~95% purity", 0.93 < df.weight.mean() < 0.96, f"{df.weight.mean():.3f}")
    check("median count rate 0.33 cts/s and redshift 0.29 as published", abs(df.cr.median() - 0.33) < 0.01
          and abs(df.z.median() - 0.29) < 0.01, f"{df.cr.median():.3f}, {df.z.median():.3f}")
    fit = K.absorption_fit(df)
    a_ref = np.exp(-fit["k"] * 3.5**fit["p"])
    check("absorbed / unabsorbed count rate is ~0.86 at N_H = 3.5e20", 0.82 < a_ref < 0.90, f"{a_ref:.3f}")
    check("every cluster lies below its optical depth limit (ZVLIM)", np.all(df.z <= df.zvlim + 1e-6))


def test_counts():
    print("Counts model")
    cm = MD.CountsModel()
    df = K.load_cosmology_sample()
    data = MD.data_counts(df)
    N = cm.expected(CO.fiducial())
    check("published best fit predicts the observed (purity-weighted) number to 2%", abs(N.sum() / data.sum() - 1) < 0.02,
          f"{N.sum():.0f} vs {data.sum():.0f}")
    lo, hi = cm.expected(CO.fiducial(sigma8=0.80)).sum(), cm.expected(CO.fiducial(sigma8=0.96)).sum()
    check("more clumpy universe, more clusters", lo < N.sum() < hi, f"{lo:.0f} < {N.sum():.0f} < {hi:.0f}")
    # Injection and recovery: Poisson mock from a known cosmology, refit by maximum likelihood
    from scipy.optimize import minimize

    rng = np.random.default_rng(4)
    truth = (0.32, 0.83)
    mock = rng.poisson(cm.expected(CO.Cosmology(Om=truth[0], sigma8=truth[1])))
    f = lambda p: -MD.loglike(cm.expected(CO.Cosmology(Om=p[0], sigma8=p[1])), mock)  # noqa: E731
    sol = minimize(f, [0.29, 0.88], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 0.01})
    s8 = sol.x[1] * np.sqrt(sol.x[0] / 0.3)
    s8_true = truth[1] * np.sqrt(truth[0] / 0.3)
    check("a Poisson mock from (0.32, 0.83) is recovered (S_8 within 0.01)", abs(s8 - s8_true) < 0.01,
          f"Om {sol.x[0]:.3f}, sigma_8 {sol.x[1]:.3f}, S_8 {s8:.3f} (true {s8_true:.3f})")


if __name__ == "__main__":
    for t in (test_cosmology, test_camb_link, test_survey, test_selection, test_catalogue, test_counts):
        t()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    sys.exit(1 if n_fail else 0)
