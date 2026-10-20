"""Self-checks for every piece of the pipeline that a wrong sign or factor could break.

    python test_physics.py          # ~1 minute, needs data/raw from fetch_data.py

Each test prints PASS/FAIL with the number it checked. The reference
implementation for data handling is ehtim, the EHT's own library.
"""

from __future__ import annotations

import os
import sys
import warnings

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402

import models as M  # noqa: E402
import uvdata as U  # noqa: E402
from config import (C, CACHE_DIR, G, M_SUN, MPC, PAPER6_COUNTS, PC, RAD_PER_UAS, RAW_DIR,  # noqa: E402
                    SHADOW_FACTOR_SCHWARZSCHILD)

M87_DIR = RAW_DIR / "m87_2017" / "EHTC_FirstM87Results_Apr2019" / "uvfits"
APR11 = M87_DIR / "SR1_M87_2017_101_lo_hops_netcal_StokesI.uvfits"

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {detail}")


def _quiet_ehtim():
    import contextlib
    import io

    with contextlib.redirect_stdout(io.StringIO()):
        import ehtim as eh
    return eh


def _load_ehtim(path):
    import contextlib
    import io

    eh = _quiet_ehtim()
    with contextlib.redirect_stdout(io.StringIO()):
        obs = eh.obsdata.load_uvfits(str(path))
    return eh, obs


# ---------------------------------------------------------------------------
def test_reader_and_averaging():
    import contextlib
    import io

    d = U.read_uvfits(APR11)
    eh, obs = _load_ehtim(APR11)
    i = np.lexsort((obs.data["t2"], obs.data["t1"], obs.data["time"]))
    j = np.lexsort((d.t2, d.t1, d.time))
    od = obs.data[i]
    dv = np.max(np.abs(od["vis"] - d.vis[j]))
    check("uvfits reader = ehtim", len(d) == len(od) and dv < 1e-5 and np.all(od["t1"] == d.t1[j]),
          f"{len(d)} records, max |dV| = {dv:.1e} Jy")
    with contextlib.redirect_stdout(io.StringIO()):
        obs.add_scans()
        oa = obs.avg_coherent(0.0, scan_avg=True)
    a = U.scan_average(d)
    emap = {}
    for r in oa.data:
        t1, t2, vis = r["t1"], r["t2"], r["vis"]
        if t1 > t2:
            t1, t2, vis = t2, t1, np.conj(vis)
        emap.setdefault((t1, t2), []).append((r["time"], vis, r["sigma"]))
    worst = 0.0
    for k in range(len(a)):
        best = min(emap[(a.t1[k], a.t2[k])], key=lambda x: abs(x[0] - a.time[k]))
        worst = max(worst, abs(best[1] - a.vis[k]), abs(best[2] - a.sigma[k]))
    check("scan averaging = ehtim avg_coherent", len(a) == len(oa.data) and worst < 1e-6,
          f"{len(a)} scan-averaged visibilities, max difference {worst:.1e} Jy")


def test_counts_against_paper6():
    worst_vis, worst_cp, worst_la = 0, 0, 0
    for (day, band), (ncp, nla) in PAPER6_COUNTS.items():
        prep = U.prepare(M87_DIR / f"SR1_M87_2017_{day}_{band}_hops_netcal_StokesI.uvfits")
        worst_cp = max(worst_cp, abs(len(prep.cphase) - ncp))
        worst_la = max(worst_la, abs(len(prep.lcamp) - nla))
    # visibility counts from Paper VI Table 1 (with and without intra-site baselines)
    table1 = {("095", "hi"): (168, 152), ("096", "hi"): (284, 250), ("096", "lo"): (274, 242),
              ("100", "hi"): (96, 86), ("100", "lo"): (91, 82), ("101", "hi"): (223, 194), ("101", "lo"): (216, 189)}
    for (day, band), (n_all, n_inter) in table1.items():
        a = U.scan_average(U.read_uvfits(M87_DIR / f"SR1_M87_2017_{day}_{band}_hops_netcal_StokesI.uvfits"))
        worst_vis = max(worst_vis, abs(len(a) - n_all), abs(int((~a.intra_site()).sum()) - n_inter))
    check("visibility counts = Paper VI Table 1", worst_vis == 0, f"largest difference {worst_vis}")
    check("closure counts ~ Paper VI Table 1", worst_cp <= 5 and worst_la <= 5,
          f"largest difference: {worst_cp} closure phases, {worst_la} log closure amplitudes")


def test_closures_match_ehtim():
    import contextlib
    import io

    prep = U.prepare(APR11)
    eh, obs = _load_ehtim(APR11)
    with contextlib.redirect_stdout(io.StringIO()):
        obs.add_scans()
        obs = obs.avg_coherent(0.0, scan_avg=True).add_fractional_noise(0.01)
        cp = obs.c_phases(mode="all", count="max", vtype="vis", ang_unit="rad")
    emap = {}
    for r in cp:
        tri = (r["t1"], r["t2"], r["t3"])
        st = tuple(sorted(tri))
        cyc = [tri, (tri[1], tri[2], tri[0]), (tri[2], tri[0], tri[1])]
        emap.setdefault(st, []).append((r["time"], (1 if st in cyc else -1) * r["cphase"], r["sigmacp"]))
    c = prep.cphase
    dmax, rmax = 0.0, 0.0
    for k in range(len(c)):
        best = min(emap[tuple(sorted(c.stations[k]))], key=lambda x: abs(x[0] - c.time[k]))
        dmax = max(dmax, abs(np.angle(np.exp(1j * (best[1] - c.values[k])))))
        rmax = max(rmax, abs(best[2] / c.sigma[k] - 1))
    check("closure phases = ehtim c_phases", np.degrees(dmax) < 1e-3 and rmax < 1e-6,
          f"{len(c)} closure phases, max difference {np.degrees(dmax):.1e} deg, sigma ratio within {rmax:.0e}")


def test_gain_invariance():
    """Closure quantities must not change when every station gets an arbitrary complex gain."""
    prep = U.prepare(APR11)
    avg = prep.avg
    rng = np.random.default_rng(7)
    gains = {(s, st): rng.uniform(0.5, 1.5) * np.exp(1j * rng.uniform(-np.pi, np.pi))
             for s in np.unique(avg.scan) for st in avg.stations}
    g = np.array([gains[(s, a)] * np.conj(gains[(s, b)]) for s, a, b in zip(avg.scan, avg.t1, avg.t2)])
    corrupted = U.replace(avg, vis=avg.vis * g, sigma=avg.sigma * np.abs(g))
    cp2 = U.closure_phases(corrupted)
    la2 = U.log_closure_amplitudes(corrupted, exclude_intra_site=False)
    same_sets = len(cp2) == len(prep.cphase) and len(la2) == len(prep.lcamp)
    dcp = np.max(np.abs(np.angle(np.exp(1j * (cp2.values - prep.cphase.values))))) if same_sets else np.inf
    # debiasing breaks exact invariance at very low S/N; compare where it is negligible
    hi = prep.lcamp.sigma < 0.2
    dla = np.max(np.abs(la2.values - prep.lcamp.values)[hi]) if same_sets else np.inf
    check("closure quantities ignore station gains", dcp < 1e-9 and dla < 2e-3,
          f"max change {dcp:.1e} rad (phases), {dla:.1e} (log amplitudes, sigma < 0.2)")


def test_chi2_is_basis_independent():
    """With the full covariance, any minimal set of closures gives the same chi^2.

    (Exactly so while residuals are small; phase wrapping breaks it for a model far from the data.)
    """
    prep = U.prepare(APR11)
    rng = np.random.default_rng(1)
    V = prep.avg.vis * np.exp(1j * rng.normal(0, 0.05, len(prep.avg)))  # a 'model' near the data
    chi_a = prep.cphase.chi2(prep.cphase.design @ np.angle(V))
    # a different basis: reverse the S/N preference when choosing the minimal set
    orig = U._greedy_basis
    U._greedy_basis = lambda rows, sig, tol=1e-9: orig(rows, [-s for s in sig], tol)
    try:
        alt = U.closure_phases(prep.avg)
    finally:
        U._greedy_basis = orig
    chi_b = alt.chi2(alt.design @ np.angle(V))
    check("closure-phase chi^2 independent of the chosen minimal set", abs(chi_a - chi_b) / chi_a < 1e-6,
          f"chi^2 = {chi_a:.4f} vs {chi_b:.4f} with two different bases")


def test_models_against_brute_force():
    from ehtim.observing.obs_helpers import ftmatrix
    from ehtim.observing.pulses import deltaPulse2D

    npix, ps = 256, 0.6
    xl = np.arange(0, -npix, -1) * ps + (ps * npix) / 2.0 - ps / 2.0
    X, Y = np.meshgrid(xl, xl)
    rng = np.random.default_rng(1)
    rho, ang = rng.uniform(0.2e9, 9e9, 60), rng.uniform(0, 2 * np.pi, 60)
    u, v = rho * np.sin(ang), rho * np.cos(ang)
    A = ftmatrix(ps * RAD_PER_UAS, npix, npix, np.column_stack([u, v]), pulse=deltaPulse2D)

    def err(model, theta, image=None):
        img = image if image is not None else model.image(theta, X, Y)
        img = img / img.sum()
        Va = model.vis(theta, u, v) / model.vis(theta, np.array([1e-3]), np.array([0.0]))[0]
        return np.max(np.abs(A @ img.ravel() - Va))

    mg = M.MGRing(2, large=False, n_gauss=1)
    th = mg.prior_transform(np.full(mg.ndim, 0.5))
    p = dict(d=43, w=10, a1=0.3, p1=np.radians(-170), a2=0.1, p2=0.7, f_floor=0.2, w_floor=20, g0_f=0.3,
             g0_x=20, g0_y=-15, g0_w=15, g0_q=0.5, g0_pa=0.4)
    for k, val in p.items():
        th[mg._index[k]] = val
    e_ring = err(mg, th)
    cr = M.Crescent(large=False)
    e_cr = err(cr, np.array([24.0, 0.35, 0.3, np.radians(160), 2.0]))
    xs = M.XSRing(n_gauss=0, large=False)
    thx = np.array([25.0, 0.4, 0.3, np.radians(160), 0.0, 0.2, 0.3])
    npx2, ps2 = 400, 0.3
    xl2 = np.arange(0, -npx2, -1) * ps2 + (ps2 * npx2) / 2.0 - ps2 / 2.0
    X2, Y2 = np.meshgrid(xl2, xl2)
    A2 = ftmatrix(ps2 * RAD_PER_UAS, npx2, npx2, np.column_stack([u, v]), pulse=deltaPulse2D)
    img = xs.crescent_image(xs.unpack(thx), X2, Y2)
    e_xs = np.max(np.abs(A2 @ img.ravel() - xs.vis(thx, u, v)))
    check("mG-ring + Gaussian = brute-force DFT of its image (ehtim)", e_ring < 1e-6, f"max |dV| = {e_ring:.1e}")
    check("crescent = brute-force DFT", e_cr < 5e-3, f"max |dV| = {e_cr:.1e} (pixelised edges)")
    check("xs-ring (slash + floor + offset) = brute-force DFT", e_xs < 5e-3, f"max |dV| = {e_xs:.1e}")
    # brightness position angle convention: eta = -arg(beta_1) points at the brightest side
    m1 = M.MGRing(1, large=False, floor=False)
    t1 = np.array([43.0, 8.0, 0.3, np.radians(-200.0)])
    img = m1.image(t1, X, Y)
    iy, ix = np.unravel_index(np.argmax(img), img.shape)
    pa = np.degrees(np.arctan2(X[iy, ix], Y[iy, ix])) % 360
    check("ring brightness PA convention (east of north)", abs(pa - 200.0) < 3, f"brightest pixel at {pa:.1f} deg, eta = 200")


def test_bessel_and_ring_null():
    from scipy.special import jv, jn_zeros

    z = np.concatenate([np.logspace(-7, -1, 100), np.linspace(0.1, 20, 300)])
    B = M._bessel_stack(z, 4)
    worst = max(np.max(np.abs(B[k] - jv(k, z))) for k in range(5))
    check("fast Bessel recurrence = scipy jv", worst < 1e-8, f"max |error| {worst:.1e} for orders 0-4")
    d = 43.0
    first_zero = jn_zeros(0, 1)[0]
    rho0 = first_zero / (np.pi * d * RAD_PER_UAS)
    ring = M.SymmetricShape("ring", large=False)
    val = abs(ring.vis(np.array([d]), np.array([rho0]), np.array([0.0]))[0])
    check("thin ring's first null at 0.766/d", val < 1e-10 and abs(rho0 / 1e9 - 3.67) < 0.01,
          f"{rho0 / 1e9:.3f} Glambda for d = 43 uas")


def test_imager_gradient():
    from imaging import ClosureImager

    preps = [U.prepare(APR11, exclude_intra_site=True)]
    im = ClosureImager(preps, npix=24, fov_uas=120, flux=0.6, amplitudes=True)
    rng = np.random.default_rng(3)
    x = np.concatenate([rng.normal(size=24 * 24) * 0.1, rng.normal(size=im.n_gain) * 0.05])
    w = dict(cphase=1, lcamp=1, amp=1, mem=1, tsv=1, tv=0.3, centroid=0.1)
    _, g = im.objective(x, w)
    worst = 0.0
    for k in list(rng.integers(0, 24 * 24, 4)) + list(24 * 24 + rng.integers(0, im.n_gain, 4)):
        e = np.zeros_like(x)
        e[k] = 1e-6
        fd = (im.objective(x + e, w)[0] - im.objective(x - e, w)[0]) / 2e-6
        worst = max(worst, abs(fd - g[k]) / max(abs(fd), 1e-3))
    check("imager analytic gradient = finite differences", worst < 1e-4, f"worst relative error {worst:.1e}")


def test_units_and_shadow():
    theta = G * 6.5e9 * M_SUN / C**2 / (16.8 * MPC) / RAD_PER_UAS
    check("M87: 6.5e9 Msun at 16.8 Mpc -> theta_g = 3.8 uas", abs(theta - 3.8) < 0.05, f"{theta:.3f} uas")
    th_s = G * 3.9565e6 * M_SUN / C**2 / (7931.1 * PC) / RAD_PER_UAS
    sh = SHADOW_FACTOR_SCHWARZSCHILD * th_s
    check("H1's Sgr A* -> shadow 51.2 uas (as H1 reported)", abs(sh - 51.17) < 0.05, f"{sh:.2f} uas")


def test_scattering_kernel():
    """Short baselines: the kernel is the published Gaussian blur, 1.380 x 0.703 mas (lambda/cm)^2 at PA 81.9."""
    k = U.Johnson2018Kernel()
    lam = 2.99792458e10 / 227.07e9
    fwhm = {"major": 1.380e3 * lam**2, "minor": 0.703e3 * lam**2}  # uas
    rho = 0.5e9
    worst = 0.0
    for axis, pa in (("major", 81.9), ("minor", 81.9 + 90.0)):
        a = np.radians(pa)
        kv = k(np.array([rho * np.sin(a)]), np.array([rho * np.cos(a)]), lam)[0]
        gauss = np.exp(-((np.pi * fwhm[axis] * RAD_PER_UAS * rho) ** 2) / (4 * np.log(2)))
        worst = max(worst, abs((1 - kv) / (1 - gauss) - 1))
    long_k = k(np.array([8e9 * np.sin(np.radians(81.9))]), np.array([8e9 * np.cos(np.radians(81.9))]), lam)[0]
    long_g = np.exp(-((np.pi * fwhm["major"] * RAD_PER_UAS * 8e9) ** 2) / (4 * np.log(2)))
    check("scattering kernel = published Gaussian blur on short baselines", worst < 0.05,
          f"{fwhm['major']:.1f} x {fwhm['minor']:.1f} uas at 1.3 mm; worst deviation {100 * worst:.1f}% at 0.5 Glambda; "
          f"at 8 Glambda {long_k:.3f} vs Gaussian {long_g:.3f} (the 800-km inner scale)")
    # and against ehtim's own ScatteringModel, made to run on numpy >= 2.5 by handing its
    # integrator a scalar integrand (the only incompatibility)
    import contextlib
    import io

    import scipy.integrate as si

    with contextlib.redirect_stdout(io.StringIO()):
        import ehtim.scattering.stochastic_optics as so

    class _Integrate:
        def __getattr__(self, name):
            return getattr(si, name)

        @staticmethod
        def quad(f, *a, **kw):
            return si.quad(lambda x: float(np.squeeze(f(x))), *a, **kw)

    saved = so.integrate
    so.integrate = _Integrate()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sm = so.ScatteringModel()
        worst_e = 0.0
        for rho in (0.5e9, 2e9, 4e9, 6e9, 8e9):
            for pa in (81.9, 171.9, 30.0):
                a = np.radians(pa)
                e = float(np.squeeze(sm.Ensemble_Average_Kernel_Visibility(rho * np.sin(a), rho * np.cos(a), lam)))
                worst_e = max(worst_e, abs(k(np.array([rho * np.sin(a)]), np.array([rho * np.cos(a)]), lam)[0] - e))
    finally:
        so.integrate = saved
    check("scattering kernel = ehtim's ScatteringModel", worst_e < 1e-5, f"max |difference| {worst_e:.1e} up to 8 Glambda")


def test_noise_model_likelihood():
    """Sgr A*'s batched per-scan covariance must reproduce the direct full-covariance chi^2."""
    import fitting as F
    import analysis as A

    preps = U.prepare_sgra_raw(A.OBS["sgra_2017_lo"]["paths"][:1], "sgra test")
    model = F.NoisyModel(M.MGRing(2, large=False))
    theta = model.prior_transform(np.random.default_rng(4).uniform(size=model.ndim))
    theta[-3:] = [-7.0, -8.0, -8.0]  # switch the extra noise off
    like = F.NoiseModelLikelihood(preps, model)
    prep = preps[0]
    V = U.scattering_kernel(prep.avg) * model.source.vis(theta[: model.source.ndim], prep.avg.u, prep.avg.v)
    direct = prep.cphase.chi2(prep.cphase.design @ np.angle(V))
    batched = like.chi2_parts(theta)["cphase"]
    check("Sgr A* noise-model likelihood = full-covariance chi^2", abs(batched / direct - 1) < 1e-5,
          f"{batched:.2f} vs {direct:.2f}")


def test_ehtim_image_fits_our_closures():
    from astropy.io import fits
    from imaging import ClosureImager

    path = CACHE_DIR / "ehtim" / "m87_2017_101.fits"
    if not path.exists():
        print("  skip  EHT pipeline image not made yet (run.py makes it)")
        return
    img = fits.open(path)[0].data.astype(float)[::-1, :]  # FITS rows run south->north
    worst = 0.0
    for band in ("lo", "hi"):
        prep = U.prepare(M87_DIR / f"SR1_M87_2017_101_{band}_hops_netcal_StokesI.uvfits", exclude_intra_site=True)
        im = ClosureImager([prep], npix=img.shape[0], fov_uas=128.0, flux=img.sum())
        red, _ = im.data_terms(img.ravel())
        worst = max(worst, red["cphase"], red["lcamp"])
    check("EHT pipeline's April 11 image fits our closure data", worst < 1.2, f"worst reduced chi^2 {worst:.2f}")


def main() -> int:
    tests = [test_reader_and_averaging, test_counts_against_paper6, test_closures_match_ehtim, test_gain_invariance,
             test_chi2_is_basis_independent, test_models_against_brute_force, test_bessel_and_ring_null,
             test_imager_gradient, test_units_and_shadow, test_scattering_kernel, test_noise_model_likelihood,
             test_ehtim_image_fits_our_closures]
    for t in tests:
        print(t.__name__)
        t()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
