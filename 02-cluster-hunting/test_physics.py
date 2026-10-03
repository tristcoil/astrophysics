"""Self-checks for the formulas and the membership model, against known answers.

Run with ``python test_physics.py``. No test framework and no network needed.

Every result in this project flows through these functions. A sign error in the
perspective geometry, or a mixture model that quietly under-deconvolves, would
produce figures that look entirely plausible and are wrong -- so each is checked
against a case whose answer is known independently.

The membership-model checks use simulated stars. That is the one place this
project uses synthetic data: to prove the code recovers a *known* answer. No
figure is ever made from it.
"""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd

import gaia_physics as gp
import membership as mb

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "",
          absolute: bool = False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if absolute else (
        abs(got - expected) / abs(expected) if expected else abs(got))
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    shown = f"{err:.3g} off" if absolute else f"{err * 100:.2f}% off"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<50} {got:>11.5g} vs {expected:>10.5g} "
          f"{unit:<7} ({shown})")


def claim(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    print("Units and distances")
    # 1 AU per year in km/s -- the standard 4.74047 of every astronomy textbook.
    check("K = 1 AU/yr in km/s", gp.K_KMS, 4.74047, 1e-6, "km/s")
    check("distance modulus at 10 pc", float(gp.distance_modulus(10.0)), 0.0, 1e-12, "mag",
          absolute=True)
    check("distance modulus at 100 pc", float(gp.distance_modulus(100.0)), 5.0, 1e-12, "mag")
    # The Sun (M_G = 4.67) seen from the Pleiades' ~135 pc would shine at G ~ 10.3.
    check("absolute mag round trip", float(gp.absolute_magnitude(10.322, 1000 / 135.0)),
          4.67, 1e-3, "mag")

    print("\nProper motion -> velocity")
    # Barnard's Star, Gaia DR3: total proper motion 10,393 mas/yr at parallax
    # 546.98 mas. Its tangential velocity is the textbook ~90 km/s.
    check("Barnard's Star tangential velocity",
          float(gp.tangential_velocity(10393.3, 546.98)), 90.07, 2e-3, "km/s")

    print("\nThree-dimensional motion")
    # Round trip: velocity from the six observables, then back again.
    v = gp.space_velocity(56.75, 24.12, 20.0, -45.5, 7.4, 5.7)
    pa, pd_, rv = gp.observed_motion(v, 56.75, 24.12, 7.4)
    check("round trip pmra", float(pa), 20.0, 1e-12, "mas/yr")
    check("round trip pmdec", float(pd_), -45.5, 1e-12, "mas/yr")
    check("round trip rv", float(rv), 5.7, 1e-12, "km/s")

    # Independent check against a published computation: the Gaia Collaboration
    # (2018, A&A 616, A10, Table A.3) Pleiades centre, motion and parallax,
    # turned into their listed ICRS velocity (U,V,W) = (-1.311, 21.39, -24.457).
    # They fitted all members jointly, so agreement at the ~0.5% level is right.
    v_gc = gp.space_velocity(56.44, 23.86, 19.997, -45.548, 7.364, 5.65)
    speed = float(np.linalg.norm([-1.311, 21.39, -24.457]))
    for comp, ours, theirs in zip("xyz", v_gc, (-1.311, 21.39, -24.457)):
        check(f"Pleiades velocity {comp} vs GC2018 (tol 1% of |v|)", float(ours), theirs,
              0.01 * speed, "km/s", absolute=True)
    ra_cp, dec_cp = gp.convergent_point(v_gc)
    check("Pleiades convergent point RA vs GC2018", ra_cp, 93.5183, 0.5, "deg", absolute=True)
    check("Pleiades convergent point Dec vs GC2018", dec_cp, -48.7831, 0.5, "deg", absolute=True)

    # Perspective: a group moving straight away from us at 10 km/s, at 100 pc.
    # A star 1 degree north of centre must show a proper motion pointing back
    # towards the centre (the group "contracts" as it recedes), of size
    # v_r * sin(1 deg) * parallax / K.
    v_away = gp.space_velocity(0.0, 0.0, 0.0, 0.0, 10.0, 10.0)
    dpa, dpd = gp.perspective_offset(0.0, 1.0, 0.0, 0.0, v_away, 10.0)
    expected = -10.0 * np.sin(np.radians(1.0)) * 10.0 / gp.K_KMS
    check("receding group contracts (pmdec, 1 deg N)", float(dpd), expected, 1e-9, "mas/yr")
    check("...with no east-west component", float(dpa), 0.0, 1e-12, "mas/yr", absolute=True)
    dpa0, dpd0 = gp.perspective_offset(0.0, 0.0, 0.0, 0.0, v, 7.4)
    claim("perspective offset vanishes at the centre", abs(dpa0) < 1e-12 and abs(dpd0) < 1e-12)

    print("\nGalactic orientation")
    # A velocity pointing at the Galactic centre is pure +U; at the north
    # Galactic pole, pure +W.
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    gc = SkyCoord(l=0 * u.deg, b=0 * u.deg, frame="galactic").icrs
    uvw = gp.galactic_velocity(10 * gp.sky_basis(gc.ra.deg, gc.dec.deg)[0])
    check("towards Galactic centre is +U", float(uvw[0]), 10.0, 1e-9, "km/s")
    ngp = SkyCoord(l=0 * u.deg, b=90 * u.deg, frame="galactic").icrs
    uvw = gp.galactic_velocity(10 * gp.sky_basis(ngp.ra.deg, ngp.dec.deg)[0])
    check("towards north Galactic pole is +W", float(uvw[2]), 10.0, 1e-9, "km/s")

    print("\nSky geometry")
    check("separation along a meridian", float(gp.angular_separation(10, 1, 10, 0)), 1.0,
          1e-12, "deg")
    check("separation across the pole", float(gp.angular_separation(0, 89, 180, 89)), 2.0,
          1e-9, "deg")
    check("position angle due north", gp.position_angle(10, 20, 10, 25), 0.0, 1e-9, "deg",
          absolute=True)
    check("position angle due east", gp.position_angle(10, 0, 11, 0), 90.0, 1e-9, "deg",
          absolute=True)
    xi, eta = gp.gnomonic(56.75, 25.12, 56.75, 24.12)
    check("gnomonic: 1 deg north", float(eta), np.degrees(np.tan(np.radians(1.0))), 1e-12, "deg")
    check("gnomonic: ...and no east offset", float(xi), 0.0, 1e-12, "deg", absolute=True)

    print("\nPhotometric excess (Riello et al. 2021)")
    # The three polynomials in their Table 2 must join up where the colour
    # ranges meet. A mistyped coefficient shows up as a jump here.
    for edge in (0.5, 4.0):
        jump = float(gp.expected_excess(edge - 1e-9) - gp.expected_excess(edge))
        check(f"C(BP-RP) continuous at {edge}", jump, 0.0, 1e-4, "", absolute=True)
    claim("clean photometry passes; a blend fails",
          bool(gp.clean_photometry(gp.expected_excess(1.0), 1.0, 15.0))
          and not bool(gp.clean_photometry(gp.expected_excess(1.0) + 0.2, 1.0, 15.0)))

    print("\nRadial-velocity magnitude trend (Katz et al. 2023)")
    # Zero at and brighter than G_RVS = 11 (the fitted parabola must meet zero
    # there), and "about 400 m/s at grvs_mag = 14" per the paper's text.
    check("RV correction at G_RVS = 11", float(gp.rv_magnitude_correction(11.0)), 0.0,
          1e-3, "km/s", absolute=True)
    check("RV correction at G_RVS = 14", float(gp.rv_magnitude_correction(14.0)), 0.40,
          0.05, "km/s", absolute=True)
    check("RV correction at G_RVS = 9", float(gp.rv_magnitude_correction(9.0)), 0.0,
          1e-12, "km/s", absolute=True)

    print("\nParallax zero point (Lindegren et al. 2021)")
    df = pd.DataFrame({
        "phot_g_mean_mag": [12.0, 4.0], "nu_eff_used_in_astrometry": [1.5, 1.5],
        "pseudocolour": [np.nan, np.nan], "ecl_lat": [0.0, 0.0],
        "astrometric_params_solved": [31, 31],
    })
    zp, extrap = gp.parallax_zero_point(df)
    claim(f"typical G=12 star: zero point {zp[0] * 1000:.1f} uas, in the expected -20..-40 range",
          -0.040 < zp[0] < -0.020)
    claim("G=4 star flagged as outside the calibration", bool(extrap[1]) and not bool(extrap[0]))

    print("\nExtreme deconvolution -- recovering a known answer")
    # One-dimensional analytic case: a single Gaussian prior N(7.4, 0.2^2) and
    # a measurement 8.0 +/- 0.4 give the textbook posterior mean
    # 7.4 + 0.04 / (0.04 + 0.16) * 0.6 = 7.52.
    one = mb.Mixture(weights=np.array([1.0]), means=np.array([[7.4]]),
                     covs=np.array([[[0.04]]]))
    post = mb.posterior_mean(np.array([[8.0]]), np.array([[[0.16]]]), one)
    check("posterior mean (analytic)", float(post[0, 0]), 7.52, 1e-12)

    rng = np.random.default_rng(42)

    def simulate(centre, sigma_true, n_cl, field_parts):
        """Cluster + field, every star blurred by its own error of 0.1-1.5."""
        truth = [rng.normal(centre, sigma_true, (n_cl, 2))]
        truth += [rng.normal(c, w, (n, 2)) for c, w, n in field_parts]
        truth = np.vstack(truth)
        errs = rng.uniform(0.1, 1.5, len(truth))
        x = truth + rng.normal(0, 1, truth.shape) * errs[:, None]
        s = np.zeros((len(x), 2, 2))
        s[:, 0, 0] = s[:, 1, 1] = errs ** 2
        return x, s, np.arange(len(x)) < n_cl

    # (a) A hard case: the cluster sits in the wing of a dense field, so the
    # two genuinely overlap. A plain fit sees the cluster ~2x too wide;
    # deconvolution must recover the true 0.5. Some contamination is
    # unavoidable here -- the test is that the model *knows how much*.
    x, s, is_cluster = simulate([5.0, 5.0], 0.5, 2000, [([0.0, 0.0], 8.0, 20000)])
    naive = float(np.sqrt(np.mean(np.var(x[is_cluster], axis=0))))
    init = mb.initialise(x, np.hypot(x[:, 0] - 5, x[:, 1] - 5) < 1.5, n_field=2, halo=False)
    fitted = mb.fit(x, s, init, max_iter=500, tol=1e-9)
    sigma_fit = float(np.sqrt(np.mean(np.diag(fitted.covs[fitted.core]))))
    print(f"         (a) overlapping populations; a plain fit sees width {naive:.3f}")
    check("deconvolved cluster width", sigma_fit, 0.5, 0.05)
    p = mb.membership_probability(x, s, fitted)
    contamination = float(np.mean(~is_cluster[p > 0.5]))
    estimate = float(np.sum(1 - p[p > 0.5]) / np.sum(p > 0.5))
    check(f"model's contamination estimate vs truth ({contamination:.1%})", estimate,
          contamination, 0.02, "", absolute=True)

    # (b) A Pleiades-like case: the cluster moves ~45 units away from a field
    # that is itself lumpy (a core plus a broad halo). Here membership should be
    # nearly perfect.
    x, s, is_cluster = simulate([20.0, -45.0], 1.0, 1300,
                                [([1.6, -3.0], 5.0, 40000), ([0.0, -5.0], 20.0, 4000)])
    init = mb.initialise(x, np.hypot(x[:, 0] - 20, x[:, 1] + 45) < 3.0, n_field=3, halo=False)
    fitted = mb.fit(x, s, init, max_iter=500, tol=1e-9)
    p = mb.membership_probability(x, s, fitted)
    recovered = float(np.mean(p[is_cluster] > 0.5))
    contamination = float(np.mean(~is_cluster[p > 0.5]))
    estimate = float(np.sum(1 - p[p > 0.5]) / np.sum(p > 0.5))
    print("         (b) well-separated cluster, lumpy field")
    claim(f"recovers {recovered:.1%} of true members (need > 97%)", recovered > 0.97)
    claim(f"contamination {contamination:.2%} of selected stars, model estimates "
          f"{estimate:.2%} (need < 2%)", contamination < 0.02)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
