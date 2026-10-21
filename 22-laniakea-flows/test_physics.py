"""Self-checks for the physics and the measuring machinery, against known answers.

Run with ``python test_physics.py``. No test framework and no network: each
check compares a formula or a code path with something known independently
(astropy's coordinates and cosmology, closed-form results, synthetic data with
a planted answer). The last group needs the cached catalogues.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np

warnings.filterwarnings("ignore")

import catalogs
import config
import cosmology as cg

FAILURES: list[str] = []


def check(label: str, got, expected, tol: float, unit: str = "", *, absolute: bool = False) -> None:
    """``got`` within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    got, expected = float(got), float(expected)
    err = abs(got - expected) if absolute or expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    print(f"  [{'ok' if ok else 'FAIL'}] {label}: {got:.6g} {unit} (expected {expected:.6g}, tol {tol:g}"
          f"{'' if absolute else ' frac'})")
    if not ok:
        FAILURES.append(label)


def test_coordinates():
    print("Coordinates and frames")
    l, b = catalogs.lonlat(catalogs.sg_to_gal([0, 0, 1]))
    check("supergalactic pole longitude", l, 47.37, 0.01, "deg", absolute=True)
    check("supergalactic pole latitude", b, 6.32, 0.01, "deg", absolute=True)
    v = np.random.default_rng(0).normal(size=(50, 3))
    check("gal -> SG -> gal round trip", np.max(np.abs(catalogs.sg_to_gal(catalogs.gal_to_sg(v)) - v)), 0, 1e-12,
          absolute=True)
    # Toward the CMB apex the CMB-frame velocity exceeds the heliocentric one by the Sun's 369.82 km/s.
    n = catalogs.unit_vector(config.SUN_CMB.l, config.SUN_CMB.b)[None, :]
    z = catalogs.helio_to_frame(np.array([0.01]), n, config.SUN_CMB)
    check("CMB-frame boost toward the apex (cz_cmb - cz_hel)", (z[0] - 0.01) * config.C_KMS / 1.01, 369.82, 1e-3, "km/s")
    zb = catalogs.frame_to_helio(z, n, config.SUN_CMB)
    check("boost round trip", zb[0], 0.01, 1e-12, absolute=True)
    lg = catalogs.local_group_motion()
    check("Local Group speed through the CMB (Planck: 620 +- 15)", lg["speed"], 620.0, 15.0, "km/s", absolute=True)
    check("... direction l (Planck: 271.9)", lg["l"], 271.9, 2.0, "deg", absolute=True)


def test_cosmology():
    print("Cosmology")
    from astropy.cosmology import FlatLambdaCDM

    ast = FlatLambdaCDM(H0=100, Om0=config.PLANCK18.Om0, Tcmb0=0)
    for z in (0.01, 0.05, 0.1):
        check(f"comoving distance at z = {z} vs astropy", cg.comoving_distance(z), ast.comoving_distance(z).value,
              1e-5, "h^-1 Mpc")
    check("z_from_comoving inverts comoving_distance", cg.z_from_comoving(cg.comoving_distance(0.07)), 0.07, 1e-6)
    check("growth rate f vs Omega_m^0.55", cg.growth_rate(0.3111), 0.3111**0.55, 0.01)
    check("growth rate f(Omega_m = 1) = 1", cg.growth_rate(1.0), 1.0, 1e-4)
    check("CF4 curvature factor f(0)", cg.cf4_curvature_factor(0.0), 1.0, 1e-12)
    dz = 1e-4
    check("... first-order slope (1 - q0)/2", (cg.cf4_curvature_factor(dz) - 1) / dz, (1 - config.CF4_Q0) / 2, 1e-3)
    check("Nusser & Tully K -> cz as z -> 0", cg.nt_K(1e-4, 1e-4) / (1e-4 * config.C_KMS), 1.0, 1e-3)


def test_estimators():
    print("Peculiar-velocity estimators")
    cz, H0 = 10000.0, 75.0
    for v in (50.0, 200.0):
        d = (cz - v) / H0
        f1 = 1.0  # isolate the logarithm from the curvature factor
        ratio = cg.vpec_wf(cz, d, H0, f=f1) / cg.vpec_direct(cz, d, H0, f=f1)
        check(f"Watkins-Feldman (ln) / direct for v = {v:g} km/s", ratio, 1.0, 0.02)
        ratio10 = cg.vpec_wf(cz, d, H0, f=f1, log=np.log10) / cg.vpec_direct(cz, d, H0, f=f1)
        check(f"... with log10 (CF4 Table 4) for v = {v:g}", ratio10, 1 / np.log(10), 0.02)
    # Gaussian errors in DM: the WF estimator stays unbiased, the direct one skews negative.
    rng = np.random.default_rng(1)
    dm_err = rng.normal(0, 0.4, 200_000)
    d = (cz / H0) * 10 ** (dm_err / 5)
    check("WF estimator unbiased under 0.4 mag errors (mean, km/s)", np.mean(cg.vpec_wf(cz, d, H0, f=1.0)), 0, 30,
          absolute=True)
    # The direct estimate inherits the lognormal mean of d: E[d] = d0 exp(sigma^2 / 2), sigma = 0.4 ln10 / 5.
    sig = 0.4 * np.log(10) / 5
    check("... while the direct estimator is biased low by cz (1 - exp(sigma^2/2))",
          np.mean((cz - H0 * d)), cz * (1 - np.exp(sig**2 / 2)), 0.05, "km/s")


def test_lcdm_velocities():
    print("Linear-theory velocity statistics")
    k, pk = cg.linear_power()
    check("sigma8 normalisation", cg.sigma_r(k, pk, 8.0), config.PLANCK18.sigma8, 1e-6)
    s, par, perp = cg.velocity_correlations()
    check("Psi_par(0) = Psi_perp(0)", par[0], perp[0], 1e-6)
    f = cg.growth_rate(config.PLANCK18.Om0)
    from scipy import integrate

    sig2 = (100 * f) ** 2 / (6 * np.pi**2) * integrate.trapezoid(pk * k, np.log(k))
    check("Psi(0) = (H0 f)^2 / 6 pi^2 int P dk", par[0], sig2, 0.01, "(km/s)^2")
    check("bulk flow of a vanishing sphere = sqrt(3) sigma_v", cg.bulk_flow_rms(0.01)[0], np.sqrt(3 * par[0]), 0.01, "km/s")
    check("velocities 300 h^-1 Mpc apart nearly uncorrelated: Psi_perp(300)/Psi_perp(0) < 0.1",
          float(np.interp(300, s, perp) / perp[0] < 0.1), 1.0, 0, absolute=True)


def test_bulk_fit():
    print("Bulk-flow fits on synthetic data")
    import bulkflow as bf
    import pandas as pd

    rng = np.random.default_rng(2)
    n = 6000
    rhat = rng.normal(size=(n, 3))
    rhat /= np.linalg.norm(rhat, axis=1)[:, None]
    r = rng.uniform(40, 300, n)
    V = np.array([-300.0, 120.0, 50.0])
    m = 80.0
    df = pd.DataFrame({"sgx": rhat[:, 0], "sgy": rhat[:, 1], "sgz": rhat[:, 2], "r": r,
                       "u": m + rhat @ V, "sig_u": np.full(n, 500.0)})
    sums = bf.shell_sums(df)
    M = np.zeros((sums["b"].size, sums["b"].size))
    fit = bf.combine(sums, M, bf._shells_for(40, 300))
    check("recovered monopole", fit["beta"][0], m, 1e-9)
    check("recovered dipole SGX", fit["beta"][1], V[0], 1e-9)
    # Shell sums combine to the same answer as one direct weighted least squares.
    df["u"] = df["u"] + rng.normal(0, 500, n)
    sums = bf.shell_sums(df)
    fit = bf.combine(sums, M, bf._shells_for(40, 300))
    X = np.column_stack([np.ones(n), rhat])
    direct = np.linalg.lstsq(X, df["u"].to_numpy(), rcond=None)[0]
    check("shell-sum combination equals direct least squares (SGY)", fit["beta"][2], direct[2], 1e-9)


def test_linear_theory_field():
    print("Linear-theory velocities and basins")
    import basins
    import density as dn

    _, cell = dn.grid_axes()
    n = config.RECON_NGRID
    c, _ = dn.grid_axes()
    X, Y, Z = np.meshgrid(c, c, c, indexing="ij")
    x0 = np.array([30.0, -20.0, 10.0])
    sig = 5.0
    blob = 2.0 * np.exp(-0.5 * ((X - x0[0]) ** 2 + (Y - x0[1]) ** 2 + (Z - x0[2]) ** 2) / sig**2)
    mass = blob.sum() * cell**3
    v = dn.velocity_field(blob, 0.5, cell)
    p = x0 + np.array([0.0, 40.0, 0.0])
    vp = dn.sample_field(v, p[None, :])[0]
    expected = 0.5 * 100 / (4 * np.pi) * mass / 40.0**2
    check("infall onto an overdense blob, 40 h^-1 Mpc away (inverse square)", -vp[1], expected, 0.02, "km/s")
    check("... and pointing straight at it (sideways component)", np.hypot(vp[0], vp[2]) / abs(vp[1]), 0, 0.01,
          absolute=True)
    # No periodic images: the far side of the box feels the blob as 1/r^2, not its image.
    q = x0 + np.array([-150.0, 0.0, 0.0])
    vq = dn.sample_field(v, q[None, :])[0]
    check("no periodic image at 150 h^-1 Mpc", vq[0], 0.5 * 100 / (4 * np.pi) * mass / 150.0**2, 0.03, "km/s")
    # Two wells: the watershed splits the box at the plane between them.
    x1 = np.array([-40.0, 0.0, 0.0])
    x2 = np.array([40.0, 0.0, 0.0])
    two = sum(np.exp(-0.5 * ((X - q_[0]) ** 2 + (Y - q_[1]) ** 2 + (Z - q_[2]) ** 2) / 8.0**2) for q_ in (x1, x2))
    psi = basins.potential(two, 0.5, cell)
    labels, minima = basins.watershed(psi)
    i1 = np.argmin(np.abs(c - (-40))), n // 2, n // 2
    i2 = np.argmin(np.abs(c - 40)), n // 2, n // 2
    check("two wells give two basins near the centre", float(labels[i1] != labels[i2]), 1.0, 0, absolute=True)
    left = labels[np.argmin(np.abs(c + 10)), n // 2, n // 2] == labels[i1]
    right = labels[np.argmin(np.abs(c - 10)), n // 2, n // 2] == labels[i2]
    check("the divide runs between them", float(left and right), 1.0, 0, absolute=True)


def test_luminosity_weights():
    print("Luminosity weights and bias")
    import density as dn

    lf = {"alpha": -0.95, "Mstar": -23.4}
    check("no light missed when the survey sees everything", dn.luminosity_weight(-5.0, lf), 1.0, 1e-6)
    check("bias of all light = 0.73 + 0.24 (alpha + 2)", dn.bias_factor(-5.0, lf), 0.73 + 0.24 * (lf["alpha"] + 2), 1e-6)
    from scipy import integrate

    M = np.linspace(-30, -5, 200_001)
    x = 10 ** (0.4 * (lf["Mstar"] - M))
    lphi = x * x ** (lf["alpha"] + 1) * np.exp(-x)
    tot = integrate.trapezoid(lphi, M)
    seen = integrate.trapezoid(np.where(M < -23.0, lphi, 0), M)
    check("weight = total light / light above the limit (numerical)", dn.luminosity_weight(-23.0, lf), tot / seen, 1e-3)


def test_data_checks():
    print("Checks on the cached catalogues")
    if not (config.CF4_DIR / "table4.dat").exists():
        print("  (catalogues not cached; skipped)")
        return
    import flows

    e = flows.estimator_check()
    check("CF4 Vpds reproduced (median |diff|, km/s)", e["Vpds_reproduced"]["median_abs_diff"], 0, 1.0, absolute=True)
    check("CF4 Vpwf reproduced with log10 (median |diff|)", e["Vpwf_reproduced"]["median_abs_diff"], 0, 1.0, absolute=True)
    g = catalogs.cf4_groups()
    off = (g["DMzp"] - g["DMav"])[g["DMav"] > 0]
    check("CF4 zero point is one offset (-0.147 mag)", np.median(off), -0.147, 1e-3, absolute=True)
    if (config.DATA_DIR / "howlett2022" / "SDSS_PV_public.dat").exists():
        import cf4fix

        ev = cf4fix.evidence()
        check("Howlett sizes are computed at the group redshift distance (MAD, dex)",
              ev["size_at_zcmb_group"]["mad"], 0, 0.001, absolute=True)
        check("CF4 SDSS moduli follow +5 eta (slope)", ev["cf4_vs_eta"]["slope"], 5.0, 0.01)


def main() -> int:
    for t in (test_coordinates, test_cosmology, test_estimators, test_lcdm_velocities, test_bulk_fit,
              test_linear_theory_field, test_luminosity_weights, test_data_checks):
        t()
    print(f"\n{'All checks passed' if not FAILURES else f'{len(FAILURES)} FAILED: ' + ', '.join(FAILURES)}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
