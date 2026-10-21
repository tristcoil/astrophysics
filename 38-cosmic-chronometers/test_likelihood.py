"""Self-checks of the data assembly, the likelihood shortcuts and the Gaussian process.

    python test_likelihood.py
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
from scipy import linalg

import config as cfg
import cosmology as cosmo
import fetch_data
import gp
import likelihood as lk

warnings.filterwarnings("ignore", category=RuntimeWarning)
FAILED = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILED.append(name)


def main() -> None:
    print("The data set")
    cc = lk.compile_chronometers()
    check("32 measurements, 15 of them Moresco's", cc["z"].size == 32 and cc["is_moresco"].sum() == 15)
    check("sorted by redshift, 0.07 to 1.965",
          np.all(np.diff(cc["z"]) >= 0) and cc["z"][0] == 0.07 and cc["z"][-1] == 1.965)
    check("methods: 15 D4000, 16 full-spectrum, 1 Lick",
          [int(np.sum(cc["method"] == m)) for m in "DFL"] == [15, 16, 1])
    k = np.argmin(np.abs(cc["z"] - 0.4497))
    check("the table's z = 0.445 is the measurement at 0.4497",
          abs(cc["z"][k] - 0.4497) < 1e-9 and abs(cc["z_table"][k] - 0.445) < 1e-9)
    m16 = fetch_data.load_m16_table()
    check("2016 table: 30 rows", m16["z"].size == 30)
    m12 = fetch_data.load_m12_d4000()
    check("2012 D4000 table: first low-mass bin as printed",
          m12["bins"]["low mass"][0][:3] == (0.16, 1.961, 0.002))

    print("Covariance")
    syst = fetch_data.load_systematics()
    mm = fetch_data.load_moresco_table("BC03")
    det = np.linalg.det(lk.chronometer_covariance(mm["z"], mm["H"], mm["err"], syst=syst))
    check("15-point determinant = primer notebook's 4.560e+33", abs(det / 4.560e33 - 1) < 1e-3,
          f"{det:.4e}")
    cov = lk.chronometer_covariance(cc["z"], cc["H"], cc["err"], syst=syst)
    check("covariance is symmetric positive definite",
          np.allclose(cov, cov.T) and np.all(np.linalg.eigvalsh(cov) > 0))

    print("Analytic H0 marginalisation vs brute force")
    z, H = cc["z"], cc["H"]
    fit = lk.fit_flat_lcdm(z, H, cov, om_axis=np.linspace(0.05, 0.95, 181))
    h_ax = np.linspace(50, 100, 2001)
    om_ax = fit.axes[0]
    chi2 = lk.flat_lcdm_plane(z, H, cov, h_ax, om_ax)
    p = np.exp(-0.5 * (chi2 - chi2.min()))
    brute_h = lk.grid_summary(h_ax, p.sum(axis=0))
    brute_om = lk.grid_summary(om_ax, p.sum(axis=1))
    ours_h = lk.grid_summary(h_ax, fit.h0_density(h_ax))
    ours_om = lk.grid_summary(*fit.marginal("Om"))
    check("H0 median and width", abs(brute_h["median"] - ours_h["median"]) < 0.02
          and abs(brute_h["std"] / ours_h["std"] - 1) < 2e-3,
          f"{ours_h['median']:.3f} vs {brute_h['median']:.3f}")
    check("Om median", abs(brute_om["median"] - ours_om["median"]) < 1e-3)
    check("H0 density integrates to 1", abs(np.trapezoid(fit.h0_density(h_ax), h_ax) - 1) < 1e-3)
    rng = np.random.default_rng(1)
    s = fit.sample(100_000, rng)
    check("posterior samples reproduce the marginal",
          abs(np.median(s["H0"]) - ours_h["median"]) < 0.1 and abs(np.median(s["Om"]) - ours_om["median"]) < 3e-3)

    print("Two-line fit vs a generic least-squares solver")
    y, cy = lk.adot_data(z, H, cov)
    zt = 0.5
    res = lk.piecewise_transition(z, y, cy, zt_axis=np.array([zt]))
    X = np.column_stack([np.ones(z.size), np.minimum(z - zt, 0), np.maximum(z - zt, 0)])
    Lc = linalg.cholesky(cy, lower=True)
    beta, *_ = np.linalg.lstsq(linalg.solve_triangular(Lc, X, lower=True),
                               linalg.solve_triangular(Lc, y, lower=True), rcond=None)
    r = y - X @ beta
    check("same chi^2 at fixed z_t", abs(res["chi2"][0] - r @ linalg.solve(cy, r)) < 1e-8)

    print("Gaussian process")
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import RBF, ConstantKernel

        dcov = np.diag(cc["err"] ** 2)
        g = gp.GaussianProcess(z, H, dcov)
        theta = np.array([np.log(130.0), np.log(1.8)])
        xs = np.linspace(0, 2, 41)
        mu, cv = g.predict(xs, theta)
        sk = GaussianProcessRegressor(ConstantKernel(130.0**2, "fixed") * RBF(1.8, "fixed"),
                                      alpha=cc["err"] ** 2, optimizer=None, normalize_y=False)
        sk.fit(z[:, None], H)
        mu_sk, sd_sk = sk.predict(xs[:, None], return_std=True)
        check("posterior mean and sd match scikit-learn (diagonal errors)",
              np.allclose(mu, mu_sk, atol=1e-6) and np.allclose(np.sqrt(np.diag(cv)), sd_sk, atol=1e-6))
        check("log marginal likelihood matches scikit-learn",
              abs(g.log_marginal(theta) - sk.log_marginal_likelihood_value_) < 1e-6)
    except ImportError:
        print("  (scikit-learn not installed: skipped)")
    om = 0.31
    zs = np.linspace(0, 2.5, 251)
    curve = cosmo.hubble(zs, 68.0, om)
    d = gp.derived_from_curves(zs, curve[None, :])
    check("derived z_t of an exact LCDM curve", abs(d["zt"][0] - cosmo.z_transition_flat(om)) < 2e-3,
          f"{d['zt'][0]:.4f} vs {float(cosmo.z_transition_flat(om)):.4f}")
    check("derived q0 of an exact LCDM curve", abs(d["q0"][0] - (1.5 * om - 1)) < 2e-3)
    check("derived H0", abs(d["H0"][0] - 68.0) < 1e-9)

    print("Why the scaling of the shared systematic matters")
    y2 = np.array([1.0, 1.5])
    C2 = np.diag((0.1 * y2) ** 2) + 0.2**2 * np.outer(y2, y2)
    one = np.ones(2)
    est = (one @ linalg.solve(C2, y2)) / (one @ linalg.solve(C2, one))
    check("Peelle's puzzle: data-scaled GLS lands below both points", est < y2.min(), f"{est:.3f}")
    C3 = np.diag((0.1 * y2) ** 2) + 0.2**2 * np.outer(np.full(2, 1.25), np.full(2, 1.25))
    est3 = (one @ linalg.solve(C3, y2)) / (one @ linalg.solve(C3, one))
    check("scaled by a common estimate instead, it sits between them", y2.min() < est3 < y2.max(),
          f"{est3:.3f}")

    print("DESI conversion")
    desi = fetch_data.load_desi_dr2()
    i = np.nonzero((desi["quantity"] == "DH_over_rs") & np.isclose(desi["z"], 0.51))[0][0]
    check("H(0.51) = c / (D_H/r_d x r_d)", abs(cfg.C_KM_S / (desi["value"][i] * 147.09) - 93.22) < 0.05,
          f"{cfg.C_KM_S / (desi['value'][i] * 147.09):.2f} km/s/Mpc with r_d = 147.09")

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {FAILED}")
        sys.exit(1)
    print("all likelihood checks passed")


if __name__ == "__main__":
    main()
