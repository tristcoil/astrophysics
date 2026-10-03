"""Self-checks against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network needed.

Groups:

1. **Physics** -- Kepler's equation, the velocity curve's conventions, and the
   mass function against textbook cases (Jupiter moves the Sun at 12.5 m/s,
   Earth at 0.09 m/s, Earth's year gives 1 AU).
2. **Periodograms** -- our generalised Lomb-Scargle and false-alarm
   probability against astropy's independent implementations, and the
   multi-instrument version against explicit least squares.
3. **Likelihood** -- the fast block-diagonal (night-noise) likelihood and the
   whitened residuals against brute-force dense linear algebra.
4. **Fitting** -- a Keplerian injected into simulated noise is recovered.
   These simulated inputs are test fixtures; no figure uses them.
5. **Transit model** -- known limits of the limb-darkened transit.
6. **Cached data** -- sanity checks on the downloads, skipped if absent.
"""

from __future__ import annotations

import sys

import numpy as np

import rv_physics as phys

FAILURES: list[str] = []


def check(name: str, condition: bool, detail: str = "") -> None:
    status = "ok  " if condition else "FAIL"
    print(f"  [{status}] {name}" + (f"  ({detail})" if detail else ""))
    if not condition:
        FAILURES.append(name)


# ---------------------------------------------------------------------------
# 1. Physics
# ---------------------------------------------------------------------------

def test_kepler_solver():
    M = np.linspace(0, 2 * np.pi, 2001)
    worst = 0.0
    for e in (0.0, 0.1, 0.5, 0.8, 0.95):
        E = phys.solve_kepler(M, e)
        resid = np.mod(E - e * np.sin(E) - M + np.pi, 2 * np.pi) - np.pi  # angles, modulo 2 pi
        worst = max(worst, np.max(np.abs(resid)))
    check("Kepler's equation solved to < 1e-11 for e up to 0.95", worst < 1e-11, f"worst {worst:.1e}")


def test_circular_orbit_is_a_sinusoid():
    t = np.linspace(0, 20, 500)
    v = phys.rv_keplerian(t, 4.2, 55.0, 0.0, 0.0, 1.3)
    expect = -55.0 * np.sin(2 * np.pi * (t - 1.3) / 4.2)
    check("circular orbit = -K sin(2 pi (t - Tc)/P)", np.allclose(v, expect, atol=1e-9))
    before = phys.rv_keplerian(1.3 - 0.1, 4.2, 55.0, 0.0, 0.0, 1.3)
    check("star recedes (v > 0) just before the planet passes in front", before > 0)


def test_eccentric_orbit_properties():
    P, K, e, w, tc = 10.0, 30.0, 0.4, 1.1, 3.0
    t = np.linspace(0, P, 200001)[:-1]
    v = phys.rv_keplerian(t, P, K, e, w, tc)
    check("velocity averages to zero over one orbit", abs(v.mean()) < 1e-6 * K, f"mean {v.mean():.2e}")
    check("velocity range is 2K", abs((v.max() - v.min()) - 2 * K) < 1e-4)
    tp = phys.tc_to_tp(tc, P, e, w)
    E = phys.solve_kepler(2 * np.pi * (tc - tp) / P, e)
    nu = phys.true_anomaly(E, e)
    check("conjunction is at true anomaly pi/2 - omega",
          abs(np.mod(nu - (np.pi / 2 - w) + np.pi, 2 * np.pi) - np.pi) < 1e-9)


def test_mass_function():
    k_jup = phys.semi_amplitude(4332.59, phys.GM_JUP, phys.GM_SUN, e=0.0489)
    check("Jupiter moves the Sun at 12.5 m/s", abs(k_jup - 12.5) < 0.1, f"{k_jup:.2f} m/s")
    k_earth = phys.semi_amplitude(365.256, phys.GM_EARTH, phys.GM_SUN, e=0.0167)
    check("Earth moves the Sun at 0.09 m/s", abs(k_earth - 0.0895) < 0.001, f"{k_earth:.4f} m/s")
    m = phys.msini_gm(4332.59, k_jup, 0.0489, phys.GM_SUN)
    check("mass function inverts K exactly (Jupiter)", abs(m / phys.GM_JUP - 1) < 1e-9)
    # the usual m << M approximation vs the exact solution, for a hot Jupiter
    P, K, M = 4.2308, 55.9, 1.07 * phys.GM_SUN
    approx = K * (P * phys.DAY_S / (2 * np.pi)) ** (1 / 3) * M ** (2 / 3)
    exact = phys.msini_gm(P, K, 0.0, M)
    check("approximate and exact m sin i agree to 0.1% for a hot Jupiter",
          abs(approx / exact - 1) < 1e-3, f"{approx / exact - 1:+.1e}")


def test_kepler_third_law():
    a = phys.semi_major_axis_au(365.256363, phys.GM_SUN, phys.GM_EARTH)
    check("Earth's sidereal year gives 1 AU", abs(a - 1) < 1e-4, f"{a:.6f} AU")


def test_alias_arithmetic():
    f = 1 / 4.23079
    p_alias = 1 / abs(1 / phys.SIDEREAL_DAY_D - f)
    check("51 Peg b's sidereal-day alias is at 1.3048 d", abs(p_alias - 1.3048) < 5e-4, f"{p_alias:.4f} d")
    p_e = 1 / abs(1 / phys.SIDEREAL_DAY_D - 1 / 0.736547)
    check("55 Cnc e's daily alias is near 2.8 d", abs(p_e - 2.81) < 0.02, f"{p_e:.3f} d")


# ---------------------------------------------------------------------------
# 2. Periodograms
# ---------------------------------------------------------------------------

def _uneven_sample(rng, n=150, span=400.0):
    """Nights drawn at random, each observation near the star's meridian
    transit, which comes 3.9 minutes earlier every night: sidereal spacing."""
    nights = np.sort(rng.choice(np.arange(int(span)), n, replace=False)).astype(float)
    return nights * phys.SIDEREAL_DAY_D + rng.normal(0, 0.05, n)


def test_gls_matches_astropy():
    from astropy.timeseries import LombScargle
    rng = np.random.default_rng(1)
    t = _uneven_sample(rng)
    err = rng.uniform(1, 3, t.size)
    y = 5 * np.sin(2 * np.pi * t / 7.3) + rng.normal(0, err) + 20
    f = np.linspace(0.005, 1.5, 3000)
    ours = phys.gls(t, y, err, f)
    ref = LombScargle(t, y, err, fit_mean=True, center_data=True, normalization="standard").power(f, method="slow")
    check("GLS equals astropy's Lomb-Scargle (floating mean, weights)", np.max(np.abs(ours - ref)) < 1e-9,
          f"max diff {np.max(np.abs(ours - ref)):.1e}")
    check("GLS finds the injected 7.3-d period", abs(1 / f[np.argmax(ours)] - 7.3) < 0.05)


def test_multi_offset_gls():
    rng = np.random.default_rng(2)
    t = _uneven_sample(rng, 200, 600)
    inst = (t > 300).astype(int)
    err = np.full(t.size, 1.5)
    y = 4 * np.cos(2 * np.pi * t / 11.1 + 0.4) + np.where(inst == 1, 40.0, -15.0) + rng.normal(0, err)
    f = np.array([1 / 11.1, 1 / 5.0, 1 / 300.0])
    power = phys.gls(t, y, err, f, inst)
    # explicit weighted least squares with indicator columns
    w = 1 / err ** 2
    J = np.column_stack([(inst == 0), (inst == 1)]).astype(float)

    def chi2(X):
        beta = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
        return np.sum(w * (y - X @ beta) ** 2)
    c0 = chi2(J)
    expect = [(c0 - chi2(np.column_stack([np.cos(2 * np.pi * fi * t), np.sin(2 * np.pi * fi * t), J]))) / c0
              for fi in f]
    check("multi-instrument GLS equals explicit least squares", np.allclose(power, expect, atol=1e-10))
    single = phys.gls(t, y, err, f)
    check("ignoring a 55 m/s instrument offset swamps a 4 m/s signal",
          power[0] > 0.5 and single[0] < 0.1, f"multi {power[0]:.2f}, single {single[0]:.3f}")


def test_fap_matches_astropy():
    from astropy.timeseries import LombScargle
    rng = np.random.default_rng(3)
    t = _uneven_sample(rng, 120, 500)
    err = rng.uniform(1, 2, t.size)
    y = rng.normal(0, err)
    ls = LombScargle(t, y, err)
    # astropy snaps the maximum frequency onto its own grid; use its value
    fmax = ls.autofrequency(maximum_frequency=1.0, return_freq_limits=True)[1]
    for z in (0.08, 0.12, 0.2):
        ref = ls.false_alarm_probability(z, method="baluev", maximum_frequency=1.0)
        ours = phys.fap_baluev(z, t, err, fmax)
        check(f"Baluev FAP at power {z} equals astropy's", abs(ours - ref) < 1e-10 * max(1, ref) + 1e-15,
              f"{ours:.3e} vs {ref:.3e}")


def test_window():
    rng = np.random.default_rng(4)
    t = _uneven_sample(rng)
    w = phys.spectral_window(t, np.array([0.0, 1 / phys.SIDEREAL_DAY_D]))
    check("spectral window is 1 at zero frequency", abs(w[0] - 1) < 1e-12)
    check("nightly sampling puts a strong window peak at 1 sidereal day", w[1] > 0.5, f"{w[1]:.2f}")


# ---------------------------------------------------------------------------
# 3. Likelihood with night-to-night noise
# ---------------------------------------------------------------------------

def _toy_data(rng):
    import rv_fit as fit
    nights = np.sort(rng.choice(np.arange(300), 40, replace=False)).astype(float)
    t = np.concatenate([n + np.sort(rng.uniform(0, 0.25, rng.integers(1, 6))) for n in nights])
    t2 = np.sort(rng.uniform(0, 300, 30)) + 0.5
    tt = np.concatenate([t, t2])
    inst = np.r_[np.zeros(t.size, int), np.ones(t2.size, int)]
    err = rng.uniform(0.8, 2.0, tt.size)
    y = 20 * np.sin(2 * np.pi * tt / 9.7) + np.where(inst == 1, 5.0, 0.0) + rng.normal(0, 2.0, tt.size)
    return fit.RVData("toy", tt, y, err, inst, ["multi", "single"])


def test_block_likelihood():
    import rv_fit as fit
    rng = np.random.default_rng(5)
    d = _toy_data(rng)
    m = fit.KeplerModel(d, [fit.PlanetSpec(P=9.7, K=20, tc=2.0, h=0.1, k=-0.2)], trend=1)
    check("night term switched on for the multi-exposure instrument", m.night_insts == ["multi"])
    th = m.x0.copy()
    th[m.names.index("jit[multi]")] = 1.3
    th[m.names.index("jit[single]")] = 2.1
    th[m.names.index("night[multi]")] = 2.4
    r = d.y - m.keplerian(th)[0]
    C = np.diag(d.err ** 2 + m.jitter(th)[0, d.inst] ** 2)
    c = m._night_variance_per_group(th)[0]
    o, s = m._g_order, m._g_starts
    for g, (a, b) in enumerate(zip(s, np.r_[s[1:], o.size])):
        idx = o[a:b]
        C[np.ix_(idx, idx)] += c[g]
    Ci = np.linalg.inv(C)
    X = m.X
    A = X.T @ Ci @ X
    beta = np.linalg.solve(A, X.T @ Ci @ r)
    res = r - X @ beta
    ll = -0.5 * (res @ Ci @ res + np.linalg.slogdet(C)[1] + d.n * np.log(2 * np.pi))
    llm = ll + 0.5 * X.shape[1] * np.log(2 * np.pi) - 0.5 * np.linalg.slogdet(A)[1]
    check("profile likelihood equals dense linear algebra", abs(m.log_like(th, marginal=False)[0] - ll) < 1e-8)
    check("marginal likelihood equals dense linear algebra", abs(m.log_like(th)[0] - llm) < 1e-8)
    z = m.whitened_residuals(th)
    check("whitened residuals reproduce r^T C^-1 r", abs(z @ z - res @ Ci @ res) < 1e-8)


# ---------------------------------------------------------------------------
# 4. Fitting (simulated data: test fixture only)
# ---------------------------------------------------------------------------

def test_injection_recovery():
    import rv_fit as fit
    rng = np.random.default_rng(6)
    t = _uneven_sample(rng, 120, 800)
    err = np.full(t.size, 3.0)
    truth = dict(P=4.2308, K=55.0, e=0.15, w=0.7, tc=100.3)
    y = phys.rv_keplerian(t, truth["P"], truth["K"], truth["e"], truth["w"], truth["tc"]) + rng.normal(0, 4.0, t.size) - 12
    d = fit.RVData("sim", t, y, err, np.zeros(t.size, int), ["sim"])
    freqs = phys.frequency_grid(t, 0.5, None, 10)
    steps, res = fit.search_planets(d, freqs, n_max=1, fap_threshold=1e-3, quiet=True)
    check("simulated planet found first by the blind search", abs(steps[0].peak_period - truth["P"]) < 0.01,
          f"{steps[0].peak_period:.4f} d")
    q = res.planets()[0]
    spec = fit.PlanetSpec(P=q["P"], K=q["K"], tc=q["tc"], h=0.2, k=0.2, P_bounds=(4.1, 4.4))
    ecc = fit.fit_map(fit.KeplerModel(d, [spec]))
    p = ecc.planets()[0]
    err_, _ = fit.laplace_errors(ecc.model, ecc.theta)
    kerr = err_[ecc.model.names.index("K1")]
    check("K recovered within 3 sigma", abs(p["K"] - truth["K"]) < 3 * kerr, f"{p['K']:.2f} +- {kerr:.2f}")
    check("e recovered within 0.05", abs(p["e"] - truth["e"]) < 0.05, f"{p['e']:.3f}")
    jit = ecc.params()["jit[sim]"]
    check("jitter recovers the excess noise (sqrt(4^2 - 3^2) = 2.6)", abs(jit - np.sqrt(16 - 9)) < 1.0,
          f"{jit:.2f} m/s")


# ---------------------------------------------------------------------------
# 5. Transit model
# ---------------------------------------------------------------------------

def test_transit_model():
    import transit as tr
    p = 0.1
    d_uniform = 1 - tr.transit_flux(np.array([0.0]), p, 0.0, 0.0)[0]
    check("uniform disc: central depth = (Rp/R*)^2", abs(d_uniform - p * p) < 1e-10)
    u1, u2 = 0.41, 0.22
    small = 0.01
    d_ld = 1 - tr.transit_flux(np.array([0.0]), small, u1, u2)[0]
    check("limb darkened: small planet at centre blocks p^2 I(0)/<I>",
          abs(d_ld / (small ** 2 / (1 - u1 / 3 - u2 / 6)) - 1) < 1e-4)
    z = np.array([1.2, 0.5, 0.95])
    f = tr.transit_flux(z, 0.1, u1, u2)
    check("no dip once the planet is off the disc", f[0] == 1.0)
    check("dip shallower near the limb than at the centre", (1 - f[2]) < (1 - f[1]))
    a = tr.a_over_r_from_density(365.256, 1.0)
    check("Kepler's third law in stellar units: Earth sits at 215 solar radii", abs(a - 215.0) < 1.5, f"{a:.1f}")


# ---------------------------------------------------------------------------
# 6. Cached data
# ---------------------------------------------------------------------------

def test_cached_data():
    from config import KECK, LICK
    if not KECK.local("table1.dat.gz").exists():
        print("  [skip] cached data not present (run fetch_data.py)")
        return
    import rvdata
    d, log = rvdata.build("51peg")
    check("51 Peg: four spectrographs loaded", {n.split()[0] for n in d.inst_names} ==
          {"ELODIE", "Lick", "HIRES", "HARPS"}, ", ".join(d.inst_names))
    check("51 Peg: 24 years of data", d.baseline > 8500, f"{d.baseline:.0f} d")
    check("times are barycentric Julian dates in the right era", 2449000 < d.t.min() < d.t.max() < 2460000)
    lick = rvdata.load_lick(rvdata.STARS["51peg"], rvdata.CutLog())
    raw = rvdata.load_lick(rvdata.STARS["51peg"], rvdata.CutLog(), convert_time=False)
    shift = np.abs(lick["t"].values - raw["t"].values).max() * 1440
    check("Lick JD -> BJD correction within the +-8.3 min light-travel time plus ~1 min UTC->TDB",
          shift < 10.5, f"max {shift:.1f} min")
    import transit as tr
    u1, u2 = tr.claret_quadratic(5196.0, 4.45, 0.31)
    check("TESS limb darkening for a K0 dwarf in the expected range", 0.3 < u1 < 0.5 and 0.1 < u2 < 0.3,
          f"u1 {u1:.3f}, u2 {u2:.3f}")


def main() -> int:
    groups = [
        ("Physics", [test_kepler_solver, test_circular_orbit_is_a_sinusoid, test_eccentric_orbit_properties,
                     test_mass_function, test_kepler_third_law, test_alias_arithmetic]),
        ("Periodograms", [test_gls_matches_astropy, test_multi_offset_gls, test_fap_matches_astropy, test_window]),
        ("Likelihood", [test_block_likelihood]),
        ("Fitting (simulated fixtures)", [test_injection_recovery]),
        ("Transit model", [test_transit_model]),
        ("Cached data", [test_cached_data]),
    ]
    for title, tests in groups:
        print(title)
        for test in tests:
            test()
    print(f"\n{'All checks passed' if not FAILURES else f'{len(FAILURES)} FAILED: ' + ', '.join(FAILURES)}")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
