"""Self-checks of the physics and the machinery (offline, about a minute).

    python test_physics.py
"""

from __future__ import annotations

import io
import json
import math
import sys
import time
import traceback
import zipfile

import numpy as np

import config as cfg

RESULTS = []


def check(name):
    def deco(fn):
        RESULTS.append((name, fn))
        return fn
    return deco


@check("Mukhanov-Sasaki: power-law inflation reproduces the exact n_s and r")
def _():
    from scipy.integrate import solve_ivp

    import inflation as inf

    for lam in (0.1, 0.2):
        eps = lam ** 2 / 2  # constant: a ∝ t^p with p = 1/eps

        def rhs(N, y):
            v, dv, _ = inf.potential(inf.EXPONENTIAL, lam, 0.0, y[0])
            e = 0.5 * y[1] ** 2
            return [y[1], -(3 - e) * (y[1] + dv / v)]

        sol = solve_ivp(rhs, (0, 60), [0.0, lam], method="DOP853", rtol=1e-13, atol=1e-14,
                        dense_output=True)
        res = inf.spectra_at(inf.Background(inf.Model("exponential", lam), 60.0, sol), 30.0)
        assert abs(res["ns"] - (1 - 2 * eps / (1 - eps))) < 2e-7, res["ns"]
        assert abs(res["r"] - 16 * eps) < 2e-7, res["r"]
        assert abs(res["alpha_s"]) < 1e-8


@check("Mukhanov-Sasaki: converged in the starting depth (300 vs 1000 times inside the horizon)")
def _():
    import inflation as inf

    bg = inf.background(inf.Model("power", 2.0))
    lnk = bg.ln_ah(bg.n_end - 55)
    a = inf.mode(bg, lnk, k_start=300)
    b = inf.mode(bg, lnk, k_start=1000)
    assert abs(a[0] - b[0]) < 1e-5 and abs(a[1] - b[1]) < 1e-5, (a, b)


@check("Slow roll: phi^2 exact vs second-order (< 1e-5) and first-order formulas (~5e-4)")
def _():
    import inflation as inf

    m = inf.Model("power", 2.0)
    bg = inf.background(m)
    res = inf.spectra_at(bg, 60.0)
    so = inf.second_order(res["eps1"], res["eps2"], res["eps3"])
    sr = m.slow_roll(60.0)
    assert abs(res["ns"] - so["ns"]) < 1e-5, (res["ns"], so["ns"])
    assert abs(res["r"] - so["r"]) / res["r"] < 2e-4
    assert 2e-4 < sr["ns"] - res["ns"] < 1e-3
    # the textbook closed forms: n_s = 1 - 2/(N + 1/2), r = 8/(N + 1/2)
    assert abs(sr["ns"] - (1 - 2 / 60.5)) < 1e-9 and abs(sr["r"] - 8 / 60.5) < 1e-9


@check("Models: the atlas spec's slow-roll numbers (N = 60) hold to the quoted precision")
def _():
    s = json.loads((cfg.DERIVED_DIR / "models.json").read_text())["models"]
    t = {k: s[k]["table"]["N60"] for k in ("phi2", "phi4", "starobinsky")}
    assert round(t["phi2"]["ns"], 3) == 0.966 and round(t["phi2"]["r"], 2) == 0.13
    assert round(t["phi4"]["ns"], 2) == 0.95 and 0.26 < t["phi4"]["r"] < 0.27
    assert round(t["starobinsky"]["ns"], 3) == 0.968 and round(t["starobinsky"]["r"], 3) == 0.003


@check("Reheating: Starobinsky's instant-reheating N* matches Planck 2018 X eq. 47")
def _():
    s = json.loads((cfg.DERIVED_DIR / "models.json").read_text())["models"]["starobinsky"]
    n_ours = s["reheating"]["instant"]["n_star"]
    # Planck 2018 X, eq. 47 with w-term = 0 (instant): N* = 67 - ln(k*/a0H0)
    #   + (1/4) ln(V*^2 / (M^4 rho_end)) - (1/12) ln g_th
    h = json.loads((cfg.DERIVED_DIR / "horizons.json").read_text())["starobinsky_instant"]
    v_star = (h["V_star_quarter_gev"] / cfg.M_PL_GEV) ** 4
    import inflation as inf

    bg = inf.background(inf.Model("starobinsky"))
    res = inf.spectra_at(bg, n_ours)
    e = inf.energy_scale(res, cfg.LN10_10_AS, cfg.M_PL_GEV)
    rho_end = e["rho_end_gev4"] / cfg.M_PL_GEV ** 4
    k_over = cfg.K_PIVOT / (cfg.H0_FID / cfg.C_KM_S)
    n_planck = 67 - math.log(k_over) + 0.25 * math.log(v_star ** 2 / rho_end) \
        - math.log(106.75) / 12
    assert abs(n_ours - n_planck) < 0.25, (n_ours, n_planck)


@check("g*(T): 106.75 above the electroweak scale, 10.75 at 10 MeV, 3.36 / 3.91 today")
def _():
    import horizons as hz

    assert abs(hz.g_star(1e4)[0] - 106.75) < 0.01
    assert abs(hz.g_star(0.01)[0] - 10.75) < 0.05
    gr, gs = hz.g_star(2.35e-13)
    assert abs(gr - 3.363) < 0.002 and abs(gs - 3.909) < 0.002, (gr, gs)


@check("Horizon: the particle horizon at last scattering subtends 1.16 degrees")
def _():
    s = json.loads((cfg.DERIVED_DIR / "horizons.json").read_text())
    h = s["horizon"]
    assert 1.1 < h["theta_H_deg"] < 1.2 and 270 < h["eta_star_mpc"] < 290
    assert 60 < s["history"]["n_min_horizon"] < 63


@check("BK18: our likelihood reproduces Cobaya's chi^2 sample by sample")
def _():
    s = json.loads((cfg.DERIVED_DIR / "bk18_validation.json").read_text())
    assert s["n"] >= 100 and s["rms_diff"] < 0.02 and s["max_abs_diff"] < 0.06, s


@check("BK18: our r posterior gives r < 0.036 (95%) and BK18's curve shape")
def _():
    s = json.loads((cfg.DERIVED_DIR / "bk18.json").read_text())
    assert abs(s["r95"] - 0.036) < 0.0015, s["r95"]
    assert s["max_abs_diff_vs_published"] < 0.05
    assert abs(s["L0_over_Lpeak"] - s["L0_over_Lpeak_published"]) < 0.04


@check("BK18 likelihood: the noise-free fiducial model is a good fit to the data (chi^2 ~ N)")
def _():
    import bk18

    like = bk18.BK18()
    chi2 = like.chi2_from_model(like.fiducial)
    n = like.nbins * like.ncl
    assert abs(chi2 - n) < 6 * math.sqrt(2 * n), (chi2, n)


@check("Chains: Planck 2018 n_s matches Planck's own GetDist summary (margestats)")
def _():
    import chains as ch

    d = ch.load("planck18_lcdm")
    mean, std = ch.mean_std(d["ns"], d["weight"])
    root = "base_plikHM_TTTEEE_lowl_lowE_lensing"
    for line in (cfg.PLANCK_DIR / root / "dist" / f"{root}.margestats").read_text().splitlines():
        if line.startswith("ns "):
            m_pub, s_pub = map(float, line.split()[1:3])
    assert abs(mean - m_pub) < 1e-4 and abs(std - s_pub) < 1e-4, (mean, std, m_pub, s_pub)
    assert abs(mean - 0.9649) < 1e-4 and abs(std - 0.0042) < 1e-4


@check("Remote zip: a stored zip inside a zip is read in place through the byte window")
def _():
    import fetch_data as fd

    inner = io.BytesIO()
    with zipfile.ZipFile(inner, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("LCDM/run/chain.1.txt", "# weight ns\n1 0.965\n")
    outer = io.BytesIO()
    with zipfile.ZipFile(outer, "w", zipfile.ZIP_STORED) as z:
        z.writestr("Compressed_Chains/LCDM.zip", inner.getvalue())
    outer.seek(0)
    zo = zipfile.ZipFile(outer)
    info = zo.getinfo("Compressed_Chains/LCDM.zip")
    outer.seek(info.header_offset)
    h = outer.read(30)
    import struct
    start = info.header_offset + 30 + sum(struct.unpack("<HH", h[26:30]))
    zi = zipfile.ZipFile(fd._SubFile(outer, start, info.compress_size))
    assert zi.read("LCDM/run/chain.1.txt").decode().endswith("0.965\n")


@check("Credible regions: a point on a Gaussian's 1-sigma ellipse reads as 68% (2D)")
def _():
    import chains as ch

    rng = np.random.default_rng(3)
    x = rng.normal(0, 1, 400_000)
    y = rng.normal(5, 1, 400_000)
    d = ch.density_2d(x, y, np.ones_like(x), (-5, 5), (0, 10), n=300, y_lower_bound=None)
    # the 68.27% HPD region of a 2D Gaussian has radius sqrt(-2 ln(1 - 0.6827)) = 1.515
    cred = ch.credibility(d, 1.515, 5.0)
    assert abs(cred - 0.6827) < 0.02, cred
    assert abs(ch.sigma_equivalent(0.6827) - 1.0) < 1e-3


@check("Reflection boundary: a half-Gaussian's 95% upper limit is 1.96 sigma")
def _():
    import chains as ch

    rng = np.random.default_rng(5)
    r = np.abs(rng.normal(0, 0.02, 400_000))
    lim = ch.upper_limit(r, np.ones_like(r))
    assert abs(lim - 1.96 * 0.02) < 0.0015, lim


def main() -> int:
    t0 = time.monotonic()
    failed = 0
    for name, fn in RESULTS:
        t = time.monotonic()
        try:
            fn()
            print(f"  PASS  {name}  ({time.monotonic() - t:.1f} s)")
        except Exception:  # noqa: BLE001
            failed += 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print(f"{len(RESULTS) - failed}/{len(RESULTS)} passed in {time.monotonic() - t0:.0f} s")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
