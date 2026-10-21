"""Check the code against answers known independently.

    python test_physics.py

Each check prints PASS or FAIL with the numbers compared; the script exits
non-zero if any check fails. Light curves called "simulated" come from
``simulate.py`` (labelled simulation) and are used here only.
"""

from __future__ import annotations

import sys

import numpy as np

import config as C
import iccf
import physics as P
import regression as R
import simulate as S

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


# ---------------------------------------------------------------------------


def loop_ccf(t1, y1, t2, y2, lags):
    """A literal, unvectorized transcription of the ICCF (White & Peterson 1994):
    for each lag, interpolate one curve at the other's shifted times, keep the
    pairs inside the interpolated curve's span, Pearson r from those pairs;
    both directions averaged."""
    out = []
    for tau in lags:
        rs = []
        for ta, ya, tb, yb, s in ((t1, y1, t2, y2, tau), (t2, y2, t1, y1, -tau)):
            xs, ys = [], []
            for ti, yi in zip(ta, ya):
                tt = ti + s
                if tb[0] <= tt <= tb[-1]:
                    xs.append(yi)
                    ys.append(np.interp(tt, tb, yb))
            if len(xs) < 2 or np.std(xs) == 0 or np.std(ys) == 0:
                rs.append(0.0)
            else:
                rs.append(np.corrcoef(xs, ys)[0, 1])
        out.append(0.5 * (rs[0] + rs[1]))
    return np.array(out)


def test_ccf():
    print("Interpolated cross-correlation function")
    rng = np.random.default_rng(1)
    t1 = np.sort(rng.uniform(0, 100, 60))
    t2 = np.sort(rng.uniform(0, 100, 45))
    y1 = rng.normal(size=60)
    y2 = rng.normal(size=45)
    lags = iccf.lag_grid(-20, 30, 0.5)
    r, _ = iccf.ccf(t1, y1, t2, y2, lags)
    r_loop = loop_ccf(t1, y1, t2, y2, lags)
    check("vectorized ICCF equals a literal loop on random data", np.allclose(r, r_loop, atol=1e-12),
          f"max diff {np.max(np.abs(r - r_loop)):.1e}")
    # An exact shifted copy of a smooth signal peaks at the shift.
    f = lambda t: np.sin(2 * np.pi * t / 37.0) + 0.5 * np.sin(2 * np.pi * t / 11.3)  # noqa: E731
    ta = np.sort(rng.uniform(0, 200, 150))
    tb = np.sort(rng.uniform(0, 200, 120))
    r, _ = iccf.ccf(ta, f(ta), tb, f(tb - 7.5), iccf.lag_grid(-20, 30, 0.25))
    pc = iccf.peak_and_centroid(iccf.lag_grid(-20, 30, 0.25), r)
    check("a copy shifted by 7.5 d peaks at 7.5 d", abs(pc.tau_peak - 7.5) <= 0.25 and r.max() > 0.95,
          f"peak {pc.tau_peak:.2f} d, r = {r.max():.4f}")
    rr, _ = iccf.ccf(ta, f(ta), ta, f(ta), np.array([0.0]))
    check("autocorrelation at zero lag is 1", abs(rr[0] - 1) < 1e-12, f"{rr[0]:.15f}")
    # Direction convention: interp2(1 vs 2, +tau) == interp1(2 vs 1, -tau)
    a, _ = iccf.ccf(t1, y1, t2, y2, lags, mode="interp2")
    b, _ = iccf.ccf(t2, y2, t1, y1, -lags, mode="interp1")
    check("swapping the light curves mirrors the lag axis", np.allclose(a, b, atol=1e-12))
    # Centroid of a symmetric peak is its centre
    lg = iccf.lag_grid(-10, 30, 0.25)
    rg = 0.9 * np.exp(-0.5 * ((lg - 8.25) / 3.0) ** 2)
    pc = iccf.peak_and_centroid(lg, rg)
    check("centroid of a symmetric peak on the grid is its centre", abs(pc.tau_cent - 8.25) < 1e-9,
          f"{pc.tau_cent:.6f} d")
    rg2 = 0.9 * np.exp(-0.5 * ((lg - 8.3) / 3.0) ** 2)
    pc2 = iccf.peak_and_centroid(lg, rg2)
    check("off the grid, the centroid is within half a grid step", abs(pc2.tau_cent - 8.3) < 0.125,
          f"{pc2.tau_cent:.3f} d for a peak at 8.30 d (0.25 d grid)")
    check("r > 0.8 r_max window of a Gaussian peak spans 2 x 0.668 sigma",
          abs((pc.window[1] - pc.window[0]) - 2 * 3.0 * np.sqrt(2 * np.log(1 / 0.8))) < 0.6,
          f"{pc.window[1] - pc.window[0]:.2f} d vs {2 * 3.0 * np.sqrt(2 * np.log(1.25)):.2f}")
    flat = iccf.peak_and_centroid(lg, np.full_like(lg, 0.1))
    check("a CCF below r = 0.2 is a failed realization", not flat.ok)


def test_fr_rss():
    print("Lag recovery and error bars on simulated campaigns (labelled simulation)")
    lags = iccf.lag_grid(-20, 40, 0.5)
    hits, meds = 0, []
    n = 30
    for k in range(n):
        tc, yc, ec, tl, yl, el = S.campaign(lag=10.0, width=4.0, days=150, seed=k)
        res = iccf.fr_rss(tc, yc, ec, tl, yl, el, lags, n_mc=200, seed=k)
        m, lo, hi = res.cent
        meds.append(m)
        hits += (m - lo) <= 10.0 <= (m + hi)
    meds = np.array(meds)
    check("median recovered lag within 1 d of the true 10 d (30 campaigns)", abs(np.median(meds) - 10.0) < 1.0,
          f"{np.median(meds):.2f} d")
    check("FR/RSS 68% intervals cover the truth 50-90% of the time", 0.5 <= hits / n <= 0.9, f"{hits}/{n}")
    # RSS keeps on average 63.2% of the points (1 - 1/e)
    rng = np.random.default_rng(3)
    t = np.arange(1000.0)
    frac = np.mean([len(iccf._rss(rng, t, t, np.ones_like(t))[0]) / 1000 for _ in range(200)])
    check("random subset selection keeps 1 - 1/e of the points", abs(frac - (1 - np.exp(-1))) < 0.005, f"{frac:.4f}")
    idx = np.random.default_rng(4).integers(0, 50, 50)
    _, counts = np.unique(idx, return_counts=True)
    _, _, e = iccf._rss(np.random.default_rng(4), np.arange(50.0), np.zeros(50), np.ones(50))
    check("a point drawn m times keeps one copy with error / sqrt(m)", np.allclose(e, 1 / np.sqrt(counts)),
          f"multiplicities up to {counts.max()}")


def test_alias():
    print("Alias weights for seasonal light curves")
    years = np.arange(6)
    t_cont = np.concatenate([365.25 * y + np.sort(np.random.default_rng(y).uniform(0, 180, 40)) for y in years])
    t_line = np.concatenate([365.25 * y + np.sort(np.random.default_rng(10 + y).uniform(10, 170, 20)) for y in years[1:5]])
    lags = iccf.lag_grid(-400, 400, 2.0)
    n = iccf.overlap_counts(t_cont, t_line, lags, 30.0)
    n0 = iccf.overlap_counts(t_cont, t_line, np.array([0.0]), 30.0)[0]
    check("N(0) counts every line epoch inside observing seasons", n0 == len(t_line), f"{n0}/{len(t_line)}")
    i180 = np.argmin(np.abs(lags - 182))
    i365 = np.argmin(np.abs(lags - 365.25 + 1))
    check("half a year of shift puts the line in the seasonal gaps; a full year does not",
          n[i180] < 0.1 * n0 and n[i365] > 0.6 * n0, f"N(182 d) = {n[i180]}, N(364 d) = {n[i365]}, N(0) = {n0}")


def test_physics():
    print("Physics and units")
    vp = P.virial_product(4.17, 4278.0)
    check("Pei et al.'s virial product: 4.17 d and 4278 km/s give 1.49e7 Msun", abs(vp / 1.49e7 - 1) < 0.01,
          f"{vp:.4e}")
    check("one light-day is 2.59e15 cm (173 AU)", abs(C.LIGHT_DAY_CM / 2.590e15 - 1) < 1e-3,
          f"{C.LIGHT_DAY_CM / C.AU_CM:.1f} AU")
    import literature as LIT

    b = LIT.bentz13_table()
    n = b[b.object == "NGC 5548"]
    pred = P.log_lambda_l5100(n.f_agn.values)
    check("Bentz et al.'s 16 NGC 5548 luminosities from their fluxes (75 Mpc)",
          np.max(np.abs(pred - n.logL.values)) <= 0.011, f"max diff {np.max(np.abs(pred - n.logL.values)):.3f} dex")
    a = P.galactic_extinction_5100(1.0) / 3.1
    check("CCM89: A(5100)/A(V) = 1.09 for R_V = 3.1", abs(a - 1.094) < 0.005, f"{a:.4f}")
    rs = P.schwarzschild_radius_cm(1.0) / 1e5
    check("Schwarzschild radius of the Sun is 2.95 km", abs(rs - 2.953) < 0.005, f"{rs:.3f} km")


def test_regression():
    print("Regression with errors and intrinsic scatter")
    rng = np.random.default_rng(5)
    x_true = rng.uniform(-2, 2, 300)
    sx = np.full(300, 0.05)
    sy = np.full(300, 0.08)
    y = 1.5 + 0.55 * x_true + rng.normal(0, 0.2, 300) + rng.normal(0, sy)
    x = x_true + rng.normal(0, sx)
    f = R.fit_line(x, sx, y, sy, nsteps=1500)
    ok = abs(f.a[0] - 1.5) < 3 * f.a[1] and abs(f.b[0] - 0.55) < 3 * f.b[1] and abs(f.sig[0] - 0.2) < 3 * f.sig[1]
    check("simulated line (a=1.5, b=0.55, scatter 0.2) recovered", ok, f.summary())
    v = rng.normal(0.6, 0.3, 200)
    e = np.full(200, 0.15)
    mu, s, _ = R.gaussian_mean_scatter(v + rng.normal(0, e), e, nsteps=1500)
    check("population mean and intrinsic scatter recovered (0.6, 0.3)",
          abs(mu[0] - 0.6) < 3 * mu[1] and abs(s[0] - 0.3) < 3 * s[1], f"{mu[0]:.3f}, {s[0]:.3f}")
    import literature as LIT

    b = LIT.bentz13_table()
    lm = R.linmix_fit(b.logL - 44, b.logL_err, np.log10(b.tau), R.log_err(b.tau, b.tau_lo, b.tau_hi))
    check("LINMIX on Bentz et al.'s table gives their K = 1.527, alpha = 0.533, var = 0.036",
          abs(lm.a[0] - 1.527) < 0.006 and abs(lm.b[0] - 0.533) < 0.006 and abs(lm.sig[0] ** 2 - 0.036) < 0.003,
          lm.summary())
    lm2 = R.linmix_fit(b.logL - 44, b.logL_err, np.log10(b.tau), R.log_err(b.tau, b.tau_lo, b.tau_hi))
    check("LINMIX with a fixed seed is exactly reproducible", np.array_equal(lm.samples, lm2.samples))


def test_parsers():
    print("Parsed tables")
    import agnmass as A
    import literature as LIT

    check("Bentz et al. 2013: 71 lags of 41 AGN", len(LIT.bentz13_table()) == 71
          and LIT.bentz13_table().object.nunique() == 41)
    check("Woo et al. 2015: 29 AGN with sigma_*", len(LIT.woo15_table()) == 29)
    check("Shen et al. 2024: 30 dynamical-model f values", len(LIT.shen24_f_table()) == 30)
    d = A.measurements()
    n = d[(d.object == "NGC5548") & (d.line == "Hb") & (d.jd_start > 56600)]
    check("AGN mass database: 86 AGN; NGC 5548 2014 H-beta lag 4.17 d", d.object.nunique() == 86
          and len(n) == 1 and abs(n.tau_cent.iloc[0] - 4.17) < 1e-9, f"{d.object.nunique()} AGN")


def test_disk():
    print("Thin accretion disk")
    import disk as D

    a = [D.thin_disk_alpha(5.2e7, md) for md in (0.1, 1.0, 10.0)]
    check("Fausnaugh et al.'s thin-disk lag scale: 0.14, 0.30, 0.65 d for 5.2e7 Msun",
          all(abs(x - y) / y < 0.03 for x, y in zip(a, (0.14, 0.30, 0.65))), ", ".join(f"{x:.3f}" for x in a))
    check("the disk lag scale grows as M^(2/3) at fixed Eddington ratio",
          abs(D.thin_disk_alpha(8 * 5.2e7) / D.thin_disk_alpha(5.2e7) - 4.0) < 1e-9)


def test_holiday_model():
    print("Line-response model of the BLR holiday")
    import ngc5548 as N

    rng = np.random.default_rng(8)
    x = np.log10(rng.uniform(20, 80, 100))
    y = 0.4 + 0.3 * x
    (a, eta), _ = N._wls(x, y, np.full(100, 0.01))
    check("weighted fit recovers a responsivity of 0.30 exactly", abs(eta - 0.3) < 1e-10 and abs(a - 0.4) < 1e-10,
          f"eta = {eta:.6f}")


def main() -> int:
    for t in (test_ccf, test_fr_rss, test_alias, test_physics, test_regression, test_parsers, test_disk,
              test_holiday_model):
        t()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail} of {len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
