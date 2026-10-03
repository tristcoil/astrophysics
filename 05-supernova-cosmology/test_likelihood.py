"""Self-checks for the statistics in ``likelihood.py``.

Run with ``python test_likelihood.py``. No network, no cached data needed.

The analysis leans on three shortcuts: whitening with a Cholesky factor instead
of inverting the covariance, solving for the nuisance parameters (M, H0) in
closed form instead of numerically, and reading credible regions off a grid.
Each is checked here against the slow, obvious way of doing the same thing.

The inputs are small synthetic vectors and covariances. They test algebra only:
no number in the project's results comes from them.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
from scipy.optimize import minimize, minimize_scalar

import cosmology as cosmo
import likelihood as lk

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float) -> None:
    err = abs(got - expected) if expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<52} {got:>13.7g} vs {expected:>13.7g} "
          f"(err {err:.1e}, tol {tol:.0e})")


def _toy_table(rng, n=60, n_cal=8):
    """A small Pantheon+-shaped table: n_cal calibrators at low z, the rest beyond 0.01."""
    z = np.sort(np.concatenate([rng.uniform(0.002, 0.009, n_cal), rng.uniform(0.015, 1.5, n - n_cal)]))
    om_true, h0_true, m_true = 0.3, 72.0, -19.25
    dl = cosmo.luminosity_distance_shape(z, z, om_true, 1 - om_true)[0] * cosmo.hubble_distance_mpc(h0_true)
    mu = cosmo.distance_modulus(dl)
    # Correlated errors: a shared calibration offset plus individual scatter.
    A = rng.normal(size=(n, n)) * 0.02
    cov = A @ A.T + np.diag(rng.uniform(0.1, 0.2, n) ** 2) + 0.01 ** 2
    m = mu + m_true + np.linalg.cholesky(cov) @ rng.normal(size=n)
    is_cal = np.zeros(n, int)
    is_cal[:n_cal] = 1
    ceph = np.where(is_cal == 1, mu + rng.normal(0, 0.05, n), -9.0)
    table = pd.DataFrame({"zHD": z, "zHEL": z, "m_b_corr": m, "IS_CALIBRATOR": is_cal,
                          "CEPH_DIST": ceph})
    return table, cov


def main() -> int:
    rng = np.random.default_rng(1998)
    table, cov = _toy_table(rng)

    print("Whitening: |W r|^2 equals r^T C^-1 r")
    sn = lk.build_sample(table, cov, name="toy")
    r = rng.normal(size=sn.n)
    check("chi^2 via Cholesky vs via solve", float(np.sum((sn.whiten @ r) ** 2)),
          float(r @ np.linalg.solve(sn.cov, r)), 1e-10)

    print("\nProfiling the degenerate offset in closed form")
    shape = lk.shape_magnitudes(sn, 0.3, 0.7)
    chi2, offset = lk.chi2_profile_offset(sn, shape)

    def chi2_at(a):
        d = sn.m - shape[0] - a
        return d @ np.linalg.solve(sn.cov, d)

    brute = minimize_scalar(chi2_at, bracket=(40.0, 45.0), tol=1e-12)
    check("best offset: closed form vs numerical minimum", float(offset[0]), float(brute.x), 1e-8)
    check("chi^2 at that offset", float(chi2[0]), float(brute.fun), 1e-9)

    print("\nSolving (M, H0) with the calibrators in closed form")
    cal = lk.build_sample(table, cov, name="toy+cal", calibrators=True)
    fit = lk.fit_calibrated(cal, np.array([0.3]))
    s_hf = lk.shape_magnitudes(cal, 0.3, 0.7)[0]

    def chi2_beta(beta):
        M, h5 = beta
        model = np.empty(cal.n)
        model[cal.is_calibrator] = cal.ceph_mu[cal.is_calibrator] + M
        model[cal.hubble_flow] = s_hf + 5 * np.log10(299_792.458) + 25 - h5 + M
        d = cal.m - model
        return d @ np.linalg.solve(cal.cov, d)

    brute = minimize(chi2_beta, x0=[-19.0, 9.0], method="Nelder-Mead",
                     options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 20000})
    check("M: closed form vs Nelder-Mead", float(fit.M[0]), float(brute.x[0]), 1e-7)
    check("H0: closed form vs Nelder-Mead", 10 ** (fit.h5[0] / 5), 10 ** (brute.x[1] / 5), 1e-7)
    check("chi^2 at the solution", float(fit.chi2[0]), float(brute.fun), 1e-8)

    print("\nH0 posterior is a normalised density")
    grid = np.linspace(60, 85, 5001)
    dens, _ = lk.h0_posterior(fit, grid)
    check("integral of p(H0) dH0", float(np.trapezoid(dens, grid)), 1.0, 1e-9)

    print("\nCredible regions read off a grid")
    x = np.linspace(-6, 6, 601)
    X, Y = np.meshgrid(x, x, indexing="ij")
    rho = 0.6  # a correlated 2D Gaussian: its HPD region is an ellipse of fixed chi^2
    chi2_2d = (X ** 2 - 2 * rho * X * Y + Y ** 2) / (1 - rho ** 2)
    p = lk.grid_posterior(chi2_2d)
    lev68, lev95, _ = lk.credible_levels(p)
    # For two parameters the 68.3% and 95.4% regions are chi^2 < 2.30 and < 6.18.
    check("68.3% region boundary at chi^2 = 2.30", float(-2 * np.log(lev68 / p.max())), 2.2957, 5e-3)
    check("95.4% region boundary at chi^2 = 6.18", float(-2 * np.log(lev95 / p.max())), 6.1801, 5e-3)
    s = lk.summarize_1d(x, p.sum(1))
    check("marginal standard deviation of a unit Gaussian", s["std"], 1.0, 1e-4)
    check("its central 68% half-width", 0.5 * (s["plus"] + s["minus"]), 1.0, 1e-3)

    print("\nBinning with a diagonal covariance reduces to weighted means")
    diag = np.diag(np.diag(cov))
    sd = lk.build_sample(table, diag, name="toy, diagonal")
    resid = rng.normal(size=sd.n)
    edges = np.array([0.01, 0.1, 0.5, 2.0])
    b = lk.binned_residuals(sd, resid, edges)
    w = 1 / np.diag(sd.cov)
    k = (sd.z_hd >= 0.1) & (sd.z_hd < 0.5)
    check("middle bin: GLS value vs weighted mean", float(b["value"][1]),
          float(np.sum(w[k] * resid[k]) / np.sum(w[k])), 1e-10)
    check("middle bin: error vs 1/sqrt(sum of weights)", float(b["err"][1]),
          float(1 / np.sqrt(np.sum(w[k]))), 1e-10)
    inside = bool(np.all((b["z"] >= edges[:-1]) & (b["z"] < edges[1:])))
    print(f"  [{'ok  ' if inside else 'FAIL'}] every bin's redshift lies inside its bin")
    if not inside:
        FAILURES.append("bin redshifts inside bins")

    print("\nThe kinematic expansion SH0ES uses matches LambdaCDM where it is used")
    z = np.linspace(0.023, 0.15, 50)
    toy = lk.Sample(name="k", rows=np.arange(z.size), z_hd=z, z_hel=z, m=z, cov=np.eye(z.size),
                    is_calibrator=np.zeros(z.size, bool), ceph_mu=z, whiten=np.eye(z.size),
                    nodes=z, index=np.arange(z.size))
    # q0 = -0.55 is flat LambdaCDM with Om = 0.3; j0 = 1 for any flat LambdaCDM.
    exact = cosmo.luminosity_distance_shape(z, z, 0.3, 0.7)[0]
    worst = float(np.max(np.abs(5 * np.log10(lk.kinematic_dl(toy) / exact))))
    check("worst |d mu| over 0.023 < z < 0.15 (mag)", worst, 0.0, 2e-4)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All likelihood checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
