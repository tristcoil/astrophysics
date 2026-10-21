"""Checks of the physics and the mathematics the measurement rests on.

    python test_physics.py     # or: python -m pytest test_physics.py

1. Unit conversions: 1 uas/yr of drift = 4.60566e-11 m/s^2 = 1.45343 km/s/Myr
   (the conversion printed in Klioner et al. 2021, Sect. 2).
2. The vector spherical harmonics are orthonormal over the sphere.
3. Their l = 1 terms are a rigid rotation and a glide, with Mignard & Klioner's
   coefficient relations (Klioner et al. 2021, Eq. 6).
4. The least-squares fitter recovers an injected glide and rotation from
   noisy fake proper motions with the right error bars, and clipping removes
   injected outliers.
5. The ICRS -> Galactic matrix maps Sgr A* and the Galactic poles where they
   belong.
6. The aberration drift itself: differentiating the first-order aberration
   formula numerically gives the glide formula.
7. The pictures: the Hammer projection inverts exactly; sky averages and the
   ring profile return a noise-free glide exactly as theory says.
8. Conditioning proper motions on zero parallax is unbiased and shrinks the
   errors by sqrt(1 - rho^2).
"""

from __future__ import annotations

import math

import numpy as np

import config as C
import vsh


def _sphere_grid(n: int = 400):
    """Gauss-Legendre in sin(dec) x uniform in ra: exact for low-degree products."""
    x, w = np.polynomial.legendre.leggauss(n)
    ra = (np.arange(2 * n) + 0.5) * np.pi / n
    dec = np.arcsin(x)
    RA, DEC = np.meshgrid(ra, dec)
    W = np.outer(w, np.full(2 * n, np.pi / n))
    return RA.ravel(), DEC.ravel(), W.ravel()


def test_units() -> None:
    assert abs(C.UAS_YR_TO_M_S2 / 4.60566e-11 - 1) < 2e-6
    assert abs(C.UAS_YR_TO_KMS_MYR / 1.45343 - 1) < 1e-5
    # 5.05 uas/yr is 2.32e-10 m/s^2 (Klioner et al. 2021, Table 2).
    assert abs(5.05 * C.UAS_YR_TO_M_S2 / 1e-10 - 2.326) < 0.001


def test_vsh_orthonormal() -> None:
    ra, dec, w = _sphere_grid(60)
    lmax = 6
    b_ra, b_dec = vsh.basis(ra, dec, lmax)
    gram = (b_ra * w[:, None]).T @ b_ra + (b_dec * w[:, None]).T @ b_dec
    assert np.allclose(gram, np.eye(vsh.n_params(lmax)), atol=1e-10), np.abs(gram - np.eye(len(gram))).max()


def test_l1_is_rotation_and_glide() -> None:
    rng = np.random.default_rng(1)
    ra = rng.uniform(0, 2 * np.pi, 500)
    dec = np.arcsin(rng.uniform(-1, 1, 500))
    tl = vsh.terms(1)
    b_ra, b_dec = vsh.basis(ra, dec, 1)
    for kind, field in (("S", vsh.glide_field), ("T", vsh.rotation_field)):
        mvec = vsh.l1_vector_matrix(kind, tl)
        for vec in np.eye(3):
            mu_ra, mu_dec = field(ra, dec, vec)
            coef, *_ = np.linalg.lstsq(np.vstack([b_ra, b_dec]), np.concatenate([mu_ra, mu_dec]),
                                       rcond=None)
            assert np.allclose(mvec @ coef, vec, atol=1e-12)
            # and the field is purely of that kind
            other = vsh.l1_vector_matrix("T" if kind == "S" else "S", tl)
            assert np.allclose(other @ coef, 0, atol=1e-12)
            # Mignard & Klioner's coefficients (Klioner et al. 2021, Eq. 6)
            mk = {t.label: c / t.mk_factor for t, c in zip(tl, coef)}
            k = kind.lower()
            s8, s4 = math.sqrt(8 * math.pi / 3), math.sqrt(4 * math.pi / 3)
            assert math.isclose(mk[f"{k}10"], s8 * vec[2], abs_tol=1e-12)
            assert math.isclose(mk[f"{k}11R"], -s4 * vec[0], abs_tol=1e-12)
            assert math.isclose(mk[f"{k}11I"], s4 * vec[1], abs_tol=1e-12)


def test_glide_points_to_apex() -> None:
    """A glide towards +x moves sources 90 deg away towards +x, by |g|."""
    g = np.array([5.0, 0.0, 0.0])
    mu_ra, mu_dec = vsh.glide_field(np.array([np.pi / 2]), np.array([0.0]), g)
    # At (ra, dec) = (90, 0) the direction towards +x is -e_ra.
    assert np.allclose([mu_ra[0], mu_dec[0]], [-5.0, 0.0])
    mu_ra, mu_dec = vsh.glide_field(np.array([0.0]), np.array([0.0]), g)
    assert np.allclose([mu_ra[0], mu_dec[0]], [0.0, 0.0])  # at the apex itself: nothing


def test_aberration_drift_is_a_glide() -> None:
    """Differentiate first-order aberration u' = u + v/c - (v.u/c) u in time."""
    rng = np.random.default_rng(2)
    u = rng.normal(size=(200, 3))
    u /= np.linalg.norm(u, axis=1)[:, None]
    a = np.array([1e-3, 2e-3, -0.5e-3])  # arbitrary units per unit time
    dt = 1e-3

    def apparent(vel):
        up = u + vel - (u @ vel)[:, None] * u
        return up / np.linalg.norm(up, axis=1)[:, None]

    expected = a - (u @ a)[:, None] * u
    # Directions referred to an observer at rest: exact to first order.
    rate = (apparent(a * dt / 2) - apparent(-a * dt / 2)) / dt
    assert np.allclose(rate, expected, rtol=0, atol=1e-9 * np.abs(a).max())
    # Gaia moves at ~30 km/s around the barycentre: Eq. 2 then holds to
    # ~|v|/c = 1e-4 of the drift (Klioner et al. 2021, Sect. 2).
    v = np.array([20.0, -15.0, 15.0]) / 299792.458
    rate = (apparent(v + a * dt / 2) - apparent(v - a * dt / 2)) / dt
    err = np.abs(rate - expected).max() / np.abs(a).max()
    assert 1e-6 < err < 3e-4, err


def test_galactic_matrix() -> None:
    A = np.array(C.A_G)
    assert np.allclose(A @ A.T, np.eye(3), atol=1e-12)
    ra, dec = np.radians(C.SGRA_RA_DEG), np.radians(C.SGRA_DEC_DEG)
    u = np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])
    x, y, z = A @ u
    l, b = math.degrees(math.atan2(y, x)) % 360, math.degrees(math.asin(z))
    # Sgr A* sits 3.3' west and 2.7' south of (l, b) = (0, 0) (Klioner et al. 2021, Sect. 3.4)
    assert abs((l - 360) * 60 - (-3.3)) < 0.1 and abs(b * 60 - (-2.7)) < 0.1, (l, b)
    # The North Galactic Pole is at ICRS (192.85948, +27.12825)
    ngp = A.T @ np.array([0, 0, 1.0])
    assert abs(math.degrees(math.atan2(ngp[1], ngp[0])) % 360 - 192.85948) < 1e-5
    assert abs(math.degrees(math.asin(ngp[2])) - 27.12825) < 1e-5


def test_fit_recovers_injected_field() -> None:
    import fit
    rng = np.random.default_rng(3)
    n = 60_000
    ra = rng.uniform(0, 2 * np.pi, n)
    dec = np.arcsin(rng.uniform(-1, 1, n))
    g = np.array([1.0, -40.0, -25.0])  # uas/yr, exaggerated so the test is quick
    w = np.array([3.0, 0.0, -2.0])
    s_ra = rng.uniform(300, 800, n)
    s_dec = rng.uniform(300, 800, n)
    rho = rng.uniform(-0.5, 0.5, n)
    # correlated Gaussian noise with the stated covariance
    z1, z2 = rng.normal(size=(2, n))
    e_ra = s_ra * z1
    e_dec = s_dec * (rho * z1 + np.sqrt(1 - rho**2) * z2)
    g_ra, g_dec = vsh.glide_field(ra, dec, g)
    w_ra, w_dec = vsh.rotation_field(ra, dec, w)
    mu_ra, mu_dec = g_ra + w_ra + e_ra, g_dec + w_dec + e_dec
    # 0.5% gross outliers
    bad = rng.random(n) < 0.005
    mu_ra[bad] += rng.normal(0, 20 * s_ra[bad])
    sample = fit.Sample(ra=ra, dec=dec, pmra=mu_ra, pmdec=mu_dec, pmra_error=s_ra,
                        pmdec_error=s_dec, pmra_pmdec_corr=rho)
    res = fit.fit_vsh(sample, lmax=3, kappa=3.0)
    pull_g = (res.glide - g) / np.sqrt(np.diag(res.glide_cov))
    pull_w = (res.rotation - w) / np.sqrt(np.diag(res.rotation_cov))
    assert np.all(np.abs(pull_g) < 4) and np.all(np.abs(pull_w) < 4), (pull_g, pull_w)
    assert 0.9 < res.uwe < 1.1
    caught = res.outlier[bad].mean()
    assert caught > 0.9, caught
    assert res.outlier[~bad].mean() < 0.005


def _fake_quasars(n, rng, glide_icrs=None, sigma=300.0):
    import pandas as pd
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    pmra = np.zeros(n)
    pmdec = np.zeros(n)
    if glide_icrs is not None:
        gra, gdec = vsh.glide_field(np.radians(ra), np.radians(dec), np.asarray(glide_icrs, float))
        pmra, pmdec = gra / 1000.0, gdec / 1000.0  # mas/yr
    return pd.DataFrame({"ra": ra, "dec": dec, "pmra": pmra, "pmdec": pmdec,
                         "pmra_error": np.full(n, sigma / 1000), "pmdec_error": np.full(n, sigma / 1000),
                         "pmra_pmdec_corr": rng.uniform(-0.3, 0.3, n)})


def test_hammer_roundtrip() -> None:
    import skymap
    rng = np.random.default_rng(4)
    l = rng.uniform(0, 360, 2000)
    b = np.degrees(np.arcsin(rng.uniform(-1, 1, 2000)))
    l2, b2 = skymap.unproject(*skymap.project(l, b))
    assert np.allclose(((l2 - l + 180) % 360) - 180, 0, atol=1e-9) and np.allclose(b2, b, atol=1e-9)
    assert skymap.project(90, 0)[0] < 0 < skymap.project(270, 0)[0]  # l increases to the left


def test_sky_averages_and_profile_of_a_pure_glide() -> None:
    import analysis
    rng = np.random.default_rng(5)
    gc = analysis.gc_unit_icrs()
    g = 5.0 * gc  # a glide of 5 uas/yr straight at Sgr A*
    df = _fake_quasars(40_000, rng, g)
    mask = np.ones(len(df), bool)
    l = np.array([30.0, 120.0, 250.0])
    b = np.array([10.0, -40.0, 60.0])
    lr, br = np.radians(l), np.radians(b)
    centres = np.column_stack([np.cos(br) * np.cos(lr), np.cos(br) * np.sin(lr), np.sin(br)])
    cap = analysis.cap_averages(df, mask, centres, 20.0, model_glide_icrs=g)
    assert np.allclose(cap["east"], cap["model_east"], atol=1e-9)
    assert np.allclose(cap["north"], cap["model_north"], atol=1e-9)
    # the average over a 20-degree cap is close to the field at its centre
    A = np.array(C.A_G)
    for i in range(3):
        u = centres[i]
        gg = A @ g
        e_l = np.array([-np.sin(lr[i]), np.cos(lr[i]), 0])
        e_b = np.array([-np.sin(br[i]) * np.cos(lr[i]), -np.sin(br[i]) * np.sin(lr[i]), np.cos(br[i])])
        assert abs(cap["east"][i] - gg @ e_l) < 0.15 and abs(cap["north"][i] - gg @ e_b) < 0.15
    prof = analysis.drift_profile(df, mask, g, None, bin_deg=10.0)
    for r in prof["rings"]:
        assert abs(r["along"] - 5.0 * np.sin(np.radians(r["theta_mean"]))) < 0.05, r
        assert abs(r["across"]) < 0.01  # the drift varies inside a 10-degree ring


def test_conditioning_on_zero_parallax() -> None:
    import pandas as pd
    import analysis
    rng = np.random.default_rng(6)
    n = 200_000
    s_a, s_d, s_p = 0.40, 0.35, 0.30  # mas(/yr)
    r_ap, r_dp, r_ad = 0.4, -0.3, 0.1
    cov = np.array([[s_a**2, r_ad * s_a * s_d, r_ap * s_a * s_p],
                    [r_ad * s_a * s_d, s_d**2, r_dp * s_d * s_p],
                    [r_ap * s_a * s_p, r_dp * s_d * s_p, s_p**2]])
    e = rng.multivariate_normal(np.zeros(3), cov, n)
    zp = -0.017
    df = pd.DataFrame({"ra": np.zeros(n), "dec": np.zeros(n), "pmra": e[:, 0], "pmdec": e[:, 1],
                       "parallax": zp + e[:, 2], "pmra_error": s_a, "pmdec_error": s_d, "parallax_error": s_p,
                       "pmra_pmdec_corr": r_ad, "parallax_pmra_corr": r_ap, "parallax_pmdec_corr": r_dp})
    s = analysis.conditioned_sample(df, np.ones(n, bool), zp)
    # unbiased, with the conditional errors that are actually realised
    assert abs(np.mean(s.pmra)) < 4 * s.pmra_error[0] / np.sqrt(n)
    assert abs(np.std(s.pmra) / s.pmra_error[0] - 1) < 0.01
    assert abs(np.std(s.pmdec) / s.pmdec_error[0] - 1) < 0.01
    assert abs(s.pmra_error[0] / 1000 - s_a * np.sqrt(1 - r_ap**2)) < 1e-12
    assert abs(np.corrcoef(s.pmra, s.pmdec)[0, 1] - s.pmra_pmdec_corr[0]) < 0.01


if __name__ == "__main__":
    import sys
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  ok   {t.__name__}")
        except Exception as exc:  # noqa: BLE001
            failed += 1
            print(f"  FAIL {t.__name__}: {exc!r}")
    sys.exit(1 if failed else 0)
