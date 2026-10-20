"""Checks on the physics and the parsers (python test_physics.py, or pytest).

Each test builds something whose answer is known and checks we get it back.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Kilonova model
# ---------------------------------------------------------------------------


def test_heating_power_law():
    """After the first seconds, r-process heating falls as t^-1.3."""
    import kilonova as K

    t = np.array([1.0, 10.0])
    q = K.heating_rate(t, 0.01)
    slope = np.log(q[1] / q[0]) / np.log(10.0)
    assert abs(slope + 1.3) < 1e-3, slope


def test_arnett_matches_quadrature():
    """The recursion equals the Arnett integral done directly."""
    from scipy.integrate import quad

    import kilonova as K

    m, v, kappa = 0.02, 8e4, 0.5
    td = K.diffusion_time(kappa, m, v)
    lum = K.bolometric(K.T_GRID, m, v, kappa)
    trap = lambda t: -np.expm1(-3 * K.KAPPA_GAMMA * m * K.M_SUN / (4 * np.pi * (v * K.KM) ** 2) / K.DAY ** 2 / t ** 2)
    for t in (0.5, 2.0, 8.0):
        f = lambda tp: (K.heating_rate(tp, m) * K.thermalisation(tp, m, v) * tp
                        * np.exp((tp ** 2 - t ** 2) / td ** 2))
        direct = 2 / td ** 2 * quad(f, 1e-4, t, limit=400, points=[1e-3, 1e-2, 0.1])[0] * trap(t)
        ours = np.interp(t, K.T_GRID, lum)
        assert abs(ours / direct - 1) < 0.01, (t, ours, direct)


def test_arnett_instantaneous_limit():
    """Long after the diffusion time, the light out equals the heat put in
    (Arnett's L -> Q for t >> t_d)."""
    import kilonova as K

    m, v, kappa = 0.002, 1.2e5, 0.01
    td = K.diffusion_time(kappa, m, v)
    t = K.T_GRID
    trap = -np.expm1(-3 * K.KAPPA_GAMMA * m * K.M_SUN / (4 * np.pi * (v * K.KM) ** 2) / K.DAY ** 2 / t ** 2)
    ratio = K.bolometric(t, m, v, kappa) / trap / (K.heating_rate(t, m) * K.thermalisation(t, m, v))
    for tt in (20 * td, 60 * td):
        assert abs(np.interp(tt, t, ratio) - 1) < 0.02, (tt, np.interp(tt, t, ratio))


def test_photosphere_stefan_boltzmann():
    import kilonova as K

    t = np.array([0.5, 2.0, 10.0])
    lum = np.array([1e42, 1e41, 1e40])
    r, temp = K.photosphere(t, lum, 3e4, 2500.0)
    assert np.allclose(4 * np.pi * r ** 2 * K.SIGMA_SB * temp ** 4, lum, rtol=1e-10)
    assert temp[-1] == 2500.0 and temp[0] > 2500.0


def test_ab_zero_point():
    """A flat 3631 Jy spectrum is magnitude 0 in any filter with our weights."""
    import kn_data

    lam, tr = kn_data.read_filter("SLOAN/SDSS.r")
    w = tr / lam
    fnu = np.full_like(lam, 3631e-23)
    mean = np.trapezoid(fnu * w, lam) / np.trapezoid(w, lam)
    assert abs(-2.5 * np.log10(mean) - 48.6) < 1e-6 + 0.002   # 48.60 vs 48.594 zero point


def test_extinction_law():
    import kilonova as K

    assert abs(K.ccm_odonnell(5495.0) - 1.0) < 0.01           # A_V / A_V at V
    assert 1.2 < K.ccm_odonnell(4400.0) < 1.4                  # A_B / A_V ~ 1.32
    assert K.ccm_odonnell(2000.0) > 2.5                        # strong in the UV


def test_villar_table():
    import kn_data

    t = kn_data.read_villar_table()
    assert t["used"].sum() == 646 and len(t["mjd"]) == 737
    assert np.all(t["svo"][t["used"]] != "")


# ---------------------------------------------------------------------------
# Knees and the two-process model on synthetic stars
# ---------------------------------------------------------------------------


def test_broken_line_recovers_knee():
    import knee

    rng = np.random.default_rng(3)
    x = rng.uniform(-2.2, -0.3, 1500)
    y = 0.32 - 0.25 * np.clip(x + 1.3, 0, None) + rng.normal(0, 0.06, x.size)
    f = knee.broken_line_fit(x, y, np.full(x.size, 0.03))
    assert abs(f["knee"] + 1.3) < 0.1 and abs(f["plateau"] - 0.32) < 0.02, f["knee"]


def test_two_sequence_recovers_knee():
    import knee

    rng = np.random.default_rng(4)
    n = 30000
    x = rng.uniform(-1.5, 0.4, n)
    w = 1 / (1 + np.exp((x + 0.2) / 0.15))
    high = rng.random(n) < w
    y = np.where(high, 0.31 - 0.35 * np.clip(x + 0.8, 0, None), 0.02 - 0.08 * x) + rng.normal(0, 0.04, n)
    f = knee.two_sequence_fit(x, y, np.full(n, 0.02), coarse=np.round(np.arange(-1.2, -0.39, 0.1), 2),
                              fine_step=0.01)
    assert abs(f["knee"] + 0.8) < 0.05 and abs(f["plateau"] - 0.31) < 0.02, (f["knee"], f["plateau"])


def test_two_process_recovers_fcc():
    """Stars built from known process vectors give back the prompt share."""
    import twoprocess as tp

    rng = np.random.default_rng(5)
    n = 60000
    plateau = 0.30
    q_cc_fe = 10 ** -plateau
    a_cc = 10 ** rng.uniform(-0.55, 0.55, n)             # [Mg/H]
    ratio = np.where(rng.random(n) < 0.3, rng.uniform(0, 0.25, n), rng.uniform(0.8, 1.3, n))
    fe = a_cc * (q_cc_fe + (1 - q_cc_fe) * ratio)
    true_fcc = 0.65                                    # element X: q_cc = 0.65, q_Ia = 0.35
    xx = a_cc * (0.65 + 0.35 * ratio)
    noise = lambda: rng.normal(0, 0.02, n)
    mg_h = np.log10(a_cc) + noise()
    fe_h = np.log10(fe) + noise()
    d = {"mg_fe": (mg_h - fe_h).astype(np.float32), "fe_h": fe_h, "mg_h": mg_h,
         "x_fe": (np.log10(xx) - fe_h + noise()).astype(np.float32)}
    low_ia = ratio < 0.5
    sel = np.ones(n, bool)
    out = tp.decompose(d, sel, low_ia, ["x"], plateau, n_boot=5)
    assert abs(out["fe"]["f_cc"] - q_cc_fe) < 0.03
    assert abs(out["x"]["f_cc"] - true_fcc) < 0.04, out["x"]["f_cc"]


# ---------------------------------------------------------------------------
# Literature inputs and the budget
# ---------------------------------------------------------------------------


def test_solar_europium():
    import budget as B

    sh = B.europium_shares()
    assert 3.0e-10 < sh["x_eu_sun"] < 4.2e-10          # Eu mass fraction in the Sun
    assert 0.9 < sh["eu_r_fraction"] < 0.98
    assert B.history_factor()[0] > 1.0                  # mergers commoner in the past


def test_asplund_and_jinabase():
    import pandas as pd

    import ncapture as nc
    from config import LIT_DIR

    a = pd.read_csv(LIT_DIR / "asplund2009_table1.csv").set_index("element")
    assert a.loc["Fe", "photosphere"] == 7.50 and a.loc["Eu", "meteorites"] == 0.51
    j = nc.load_jinabase()
    assert j["eu_upper"].sum() > 100 and np.all(np.isfinite(j["eu_fe"][j["eu_upper"]]))
    assert abs(nc.pure_r_ba_eu() + 0.80) < 0.05


if __name__ == "__main__":
    import sys
    import time

    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_") and callable(v)]
    failed = 0
    for name, fn in tests:
        t0 = time.time()
        try:
            fn()
            print(f"PASS {name} ({time.time() - t0:.1f} s)")
        except Exception as exc:   # report and keep going
            failed += 1
            print(f"FAIL {name}: {exc!r}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
