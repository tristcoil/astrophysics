"""Checks of the physics and of the timing pipeline against published numbers and
injected signals.  Run with ``python -m pytest -q test_physics.py`` or ``python test_physics.py``.

The injection test needs the cached ZTF photometry (it reuses the real exposure
times and errors); the others run on their own.
"""

from __future__ import annotations

import numpy as np

import config as C
import fates
import gw
import population as pop
import timing as T


def test_robson_j0651_example():
    # Robson, Cornish & Liu (2019), section 4: SDSS J0651 at ~1 kpc, 0.5 + 0.25 Msun,
    # f = 2.6 mHz, 4-year mission -> SNR 140.
    mc = gw.chirp_mass(0.5, 0.25)
    A = gw.amplitude(2.6e-3, mc, 1.0 * C.KPC)
    snr = gw.snr(2.6e-3, A, 4 * C.YEAR, t_obs_yr=4.0)
    assert abs(snr - 140) < 3, snr


def test_kupfer2018_amplitudes():
    # Kupfer et al. (2018) Tables 1-2: script-A in units of 1e-23
    cases = [("AM CVn", 1028.73, 0.68, 0.125, 299, 28.3), ("HM Cnc", 321.529, 0.55, 0.27, 5000, 6.4),
             ("SDSS J0651", 765.5, 0.247, 0.49, 933, 16.2), ("ES Cet", 620.21, 0.8, 0.161, 1584, 10.7)]
    for name, P, m1, m2, d, a_pub in cases:
        A = gw.amplitude(2 / P, gw.chirp_mass(m1, m2), d * C.PC) / 1e-23
        assert abs(A - a_pub) / a_pub < 0.01, (name, A, a_pub)


def test_peters_timescale_and_inversion():
    P, m1, m2 = 414.79, 0.61, 0.21
    mc = gw.chirp_mass(m1, m2)
    pdot = gw.pdot_gw(P, mc)
    # chirp mass round trip
    assert abs(gw.chirp_mass_from_pdot(P, pdot) - mc) < 1e-12
    # for pure GR decay the time to coalescence is (3/8) P / |P-dot|
    assert abs(gw.time_to_coalescence(P, m1, m2) / (3 / 8 * P / abs(pdot)) - 1) < 1e-9
    # f-dot and P-dot are the same physics: f = 2/P -> fdot = -2 Pdot / P^2
    assert abs(gw.fdot_gw(2 / P, mc) / (-2 * pdot / P**2) - 1) < 1e-9


def test_eclipse_geometry():
    # an occulter twice the size, centred: total eclipse; far apart: nothing hidden
    assert abs(T.blocked_fraction([0.0], 0.1, 0.2)[0] - 1) < 1e-6
    assert T.blocked_fraction([0.5], 0.1, 0.2)[0] == 0.0
    # uniform disc half covered by a huge occulter edge-on: ~0.5
    assert abs(T.blocked_fraction([100.0], 1.0, 100.0)[0] - 0.5) < 0.01


def test_marsh_boundaries_and_contact():
    # Marsh et al. (2004) fig. 1: at M1 = 1.0 the "always unstable" line is near M2 = 0.53-0.55
    lo, hi = fates.boundary_curves(np.array([1.0]))
    assert 0.50 < hi[0] < 0.58 and 0.15 < lo[0] < 0.25
    st = fates.stability(0.61, 0.21)
    assert st["label"] == "depends on tides"
    # Burdge et al. (2019): ~130,000 years until ZTF J1539 reaches P ~ 5 min
    t5 = gw.time_to_coalescence(414.79, 0.61, 0.21) - gw.time_to_coalescence(300.0, 0.61, 0.21)
    assert 1.2e5 < t5 / C.YEAR < 1.4e5


def test_maoz_distribution_integral():
    # Korol et al. (2022) eq. 9: the closed-form integral of Maoz et al.'s N(x)
    a = -1.3
    G = lambda x: x**a / (a + 5) * (x * (1 + x**4) * (1 + x**-4) ** ((a + 1) / 4) - x**5)
    tab = pop.SeparationTable(a)
    for x1, x2 in [(0.0125, 0.45), (0.1, 1.0), (1.0, 10.0)]:
        assert abs((tab.F(x2) - tab.F(x1)) / abs(G(x2) - G(x1)) - 1) < 2e-3


def test_kupfer_table_counts():
    import kupfer

    df = kupfer.parse()
    assert (df["status"] == "verification").sum() == 18
    assert (df["status"] == "detectable").sum() == 22
    assert np.allclose(df["f_mHz"], df["f_mHz_table"], atol=1.5e-3)
    v = df[df["status"] == "verification"]
    assert (v["type"] == "AMCVn").sum() == 10 and (v["type"] == "DWD").sum() == 8


def injection_recovery(n_real=4, seed=11):
    """Simulate ZTF data with a known ephemeris on the real exposure times and errors,
    then fit it exactly as the real data are fitted. Returns (true, recovered, errors)."""
    import pickle

    cache = C.DATA_DIR / "analysis_cache.pkl"
    res, art = pickle.loads(cache.read_bytes())
    fit = art["ztf_fit"]
    shapes = art["shapes"]
    rng = np.random.default_rng(seed)
    th_true = np.array(fit["theta"], float)
    th_true[2] *= 1.05  # a P-dot 5% different from the real one
    out = []
    for _ in range(n_real):
        sim_bands = []
        for b in fit["bands"]:
            B = T.smeared_basis(b, th_true)
            y = B @ shapes[b.name]["coef"] + b.G @ np.array(list(shapes[b.name]["offsets"].values()))
            y = y + rng.normal(0, b.sig)
            nb = T.BandData(b.name, b.t, b.exptime_s, y, b.sig, b.groups, b.group_names, b.n_harm, b.sines,
                            b.eclipses)
            sim_bands.append(nb)
        start = th_true + np.array([0.02, 2e-6, 0.0])
        th, cov, c2, dof = T.fit_ephemeris(sim_bands, start)
        out.append((-th[2] / th[1] ** 2, np.sqrt(cov[2, 2]) / th[1] ** 2, (th[0] - th_true[0]) / th[1] * 86400,
                    np.sqrt(cov[0, 0]) / th[1] * 86400))
    pdot_true = -th_true[2] / th_true[1] ** 2
    return pdot_true, np.array(out)


def test_injection():
    pdot_true, out = injection_recovery()
    pulls = (out[:, 0] - pdot_true) / out[:, 1]
    assert np.all(np.abs(pulls) < 4), pulls
    assert np.all(np.abs(out[:, 2] / out[:, 3]) < 4)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    p, o = injection_recovery(n_real=8)
    print(f"injection: true P-dot {p:.5e}; recovered {np.mean(o[:, 0]):.5e} +- {np.std(o[:, 0]):.1e} "
          f"(mean quoted error {np.mean(o[:, 1]):.1e}); pulls {np.round((o[:, 0] - p) / o[:, 1], 2)}; "
          f"phase-zero offsets {np.round(o[:, 2], 2)} s (errors ~{np.mean(o[:, 3]):.2f} s)")
