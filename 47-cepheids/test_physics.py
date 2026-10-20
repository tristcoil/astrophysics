"""Tests of the machinery on synthetic inputs with known answers.

    python test_physics.py        (or: python -m pytest -q test_physics.py)

Each check prints PASS or FAIL; the script exits non-zero if any fails. One test
borrows the time stamps of a cached OGLE light curve for realistic sampling
(random times if the cache is empty). Synthetic data appear only here, to check
the code; no figure uses them.
"""

from __future__ import annotations

import numpy as np
from scipy.special import i0

import config as C
import gaia_cal as G
import geometry as GEO
import leavitt
import lightcurves as LC
import plfit

RNG = np.random.default_rng(12345)


class approx:
    """Minimal stand-in for pytest.approx (absolute tolerance), so no extra package is needed."""

    def __init__(self, expected, abs):
        self.expected, self.tol = expected, abs

    def __eq__(self, other):
        return bool(np.all(np.abs(np.asarray(other, float) - np.asarray(self.expected, float)) <= self.tol))

    def __repr__(self):
        return f"{self.expected} ± {self.tol}"


def test_intensity_mean_of_a_sinusoid():
    # m = m0 + A cos(wt): <10^(-0.4 m)> = 10^(-0.4 m0) I0(0.4 ln10 A), a modified Bessel function
    m0, A = 15.0, 0.4
    coef = np.array([m0, A, 0.0])
    expect = m0 - 2.5 * np.log10(i0(0.4 * np.log(10) * A))
    assert LC.intensity_mean(coef) == approx(expect, abs=1e-6)
    assert LC.intensity_mean(coef) < m0  # the flux mean is brighter than the magnitude mean


def _ogle_times():
    path = C.OGLE_DIR / "lmc" / "phot" / "I" / "OGLE-LMC-CEP-0002.dat"
    if path.exists():
        return np.loadtxt(path)[:, 0]
    return np.sort(RNG.uniform(5260, 8920, 420))


def test_period_and_mean_recovered_with_ogle_sampling():
    t = _ogle_times()
    P, m0 = 7.3315, 14.2
    coef = np.array([m0, 0.25, -0.12, 0.05, 0.06, -0.02, 0.03])  # a 3-harmonic Cepheid-like shape
    y = LC.model(coef, (t / P) % 1.0) + RNG.normal(0, 0.008, len(t))
    e = np.full(len(t), 0.008)
    r = LC.analyse(t, y, e)
    assert abs(r["P"] - P) / P < 2e-5
    assert abs(r["P_blind"] - P) / P < 1e-3
    assert r["mean"] == approx(LC.intensity_mean(coef), abs=0.002)
    assert r["eP"] > 0


def test_fourier_shape_phase_convention_is_shift_invariant():
    coef = np.array([0.0, 0.3, 0.1, 0.08, -0.05, 0.02, 0.01])
    r1 = LC.fourier_shape(coef)
    r2 = LC.fourier_shape(LC._shift_to_max(coef, 0.37))
    assert r1 == approx(r2, abs=1e-9)


def test_clipped_line_ignores_outliers():
    x = RNG.uniform(0, 1.5, 2000)
    y = 15.9 - 3.31 * x + RNG.normal(0, 0.08, x.size)
    y[:60] -= RNG.uniform(0.5, 1.5, 60)  # blends: too bright
    f = plfit.fit_line(x, y, x0=0.0)
    assert f["b"] == approx(-3.31, abs=0.02)
    assert f["sigma"] == approx(0.08, abs=0.01)
    assert f["n_clipped"] >= 55


def test_wesenheit_is_reddening_free():
    V, I = 15.5, 14.7
    for E in (0.0, 0.1, 0.5):  # E(V - I); A_I = 1.55 E
        Vr, Ir = V + 2.55 * E, I + 1.55 * E
        assert Ir - C.R_WI * (Vr - Ir) == approx(I - C.R_WI * (V - I), abs=1e-12)


def test_disk_fit_recovers_an_injected_plane():
    ra0, dec0 = 80.05, -69.30
    n = 3000
    rho = np.radians(np.sqrt(RNG.uniform(0, 1, n)) * 5.0)
    pa = RNG.uniform(0, 2 * np.pi, n)
    # back to RA/Dec from (rho, PA) on the sphere
    d0, a0 = np.radians(dec0), np.radians(ra0)
    dec = np.arcsin(np.sin(d0) * np.cos(rho) + np.cos(d0) * np.sin(rho) * np.cos(pa))
    ra = a0 + np.arctan2(np.sin(pa) * np.sin(rho) * np.cos(d0), np.cos(rho) - np.sin(d0) * np.sin(dec))
    dmu = 5 * np.log10(GEO.disk_distance_ratio(rho, pa, 24.0, 151.0)) + RNG.normal(0, 0.07, n)
    f = GEO.fit_disk(np.degrees(ra), np.degrees(dec), dmu, ra0, dec0)
    assert f["incl"] == approx(24.0, abs=2.0)
    assert f["theta"] == approx(151.0, abs=4.0)
    assert 30 < f["near_side_pa"] < 90  # north-east


def test_distance_modulus_conversions():
    assert GEO.distance_from_mu(18.477) == approx(49.59, abs=0.01)
    assert GEO.mu_from_distance(62.44) == approx(18.977, abs=0.001)


def _synthetic_lmc(n=2500, R=1.8):
    """An LMC with a W_G law, intrinsic colours and small reddening, at mu = 18.48."""
    logp = RNG.uniform(0.0, 1.5, n)
    col0 = 0.88 + 0.27 * (logp - 1) + RNG.normal(0, 0.06, n)  # intrinsic BP - RP
    E = RNG.uniform(0.05, 0.2, n)
    Mw = -5.80 - 3.33 * (logp - 1) + RNG.normal(0, 0.07, n)  # absolute W at R
    G0 = Mw + R * col0
    g = G0 + 18.48 + R * E
    bp = g + 0.4 * (col0 + E)  # any split of the colour into BP, RP works
    rp = bp - (col0 + E)
    return logp, g, bp, rp


def test_joint_fit_recovers_the_lmc_distance():
    global RNG
    RNG = np.random.default_rng(2026)  # independent of the order the tests run in
    R_true, mu_true, gamma_true, delta_true, k_true = 1.8, 18.48, -0.3, 0.020, 1.5
    logp_l, g_l, bp_l, rp_l = _synthetic_lmc(R=R_true)
    shapes = {"F": G.LMCShape(logp_l, g_l, bp_l, rp_l)}
    n = 700
    logp = RNG.uniform(0.3, 1.4, n)
    feh = RNG.normal(0.05, 0.15, n)
    E = RNG.uniform(0.1, 1.5, n)
    col0 = 0.88 + 0.27 * (logp - 1) + RNG.normal(0, 0.06, n)
    d_kpc = RNG.uniform(0.3, 4.0, n)
    mu = 5 * np.log10(d_kpc * 1e3) - 5
    Mw = -5.80 - 3.33 * (logp - 1) + gamma_true * (feh + 0.409) + RNG.normal(0, 0.07, n)
    g = Mw + R_true * col0 + mu + R_true * E
    bp = g + 0.4 * (col0 + E)
    rp = bp - (col0 + E)
    e_plx = 0.012 + 0.01 * RNG.uniform(0, 1, n)
    zp = np.full(n, -0.030)
    plx_true = 1.0 / d_kpc
    plx = plx_true + zp + delta_true + RNG.normal(0, k_true * e_plx)
    mw = {"int_average_g": g, "int_average_bp": bp, "int_average_rp": rp, "logP": logp, "parallax": plx,
          "parallax_error": e_plx, "zp_l21": zp, "mode": np.array(["F"] * n)}
    jf = G.JointFit(mw, np.ones(n, bool), feh, shapes, -0.409,
                    free=("mu", "gamma", "delta", "R", "sigR", "k"), fixed={"lnsig": np.log(0.07)})
    jf.fit()
    p = jf.params(jf.theta)
    got = {k: round(float(v), 4) for k, v in p.items()}
    assert p["mu"] == approx(mu_true, abs=0.05), got
    assert p["delta"] == approx(delta_true, abs=0.006), got
    assert p["R"] == approx(R_true, abs=0.05), got
    assert p["k"] == approx(k_true, abs=0.25), got


def test_leavitt_table_is_self_consistent():
    t = leavitt.table()
    assert len(t["H"]) == 25 and len(set(t["H"])) == 25
    assert np.all(np.diff(t["P"]) > 0)  # "arranged in the order of their length"
    assert int(t["new"].sum()) == 9
    # her residuals are observed - line: subtracting them must leave two straight lines
    for key, res in (("max", "res_max"), ("min", "res_min")):
        line = t[key] - t[res]
        ok = ~((t["H"] == 1374) & (key == "min"))  # the probable sign misprint
        coef = np.polyfit(t["logP"][ok], line[ok], 1)
        scatter = np.std(line[ok] - np.polyval(coef, t["logP"][ok]))
        assert scatter < 0.06, (key, scatter)
        assert coef[0] == approx(-1 / leavitt.LEAVITT_DLOGP_PER_MAG, abs=0.1)


def test_alias_detector():
    t = {"H": np.array([1, 2]), "P": np.array([1.6637, 4.0]), "P_ogle": np.array([1.6560782, 4.0]),
         "e_P": np.array([1e-4, 1e-3])}
    found = leavitt.aliases(t)
    assert len(found) == 1 and found[0]["alias"] == "1 cycle/year"


def main() -> int:
    import inspect
    import sys
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and inspect.isfunction(f)]
    failed = 0
    for name, f in tests:
        try:
            f()
            print(f"PASS  {name}")
        except AssertionError as exc:
            failed += 1
            print(f"FAIL  {name}: {exc}")
    print(f"{len(tests) - failed}/{len(tests)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
