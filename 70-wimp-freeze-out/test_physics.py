"""Self-checks: closed forms, conservation laws, and independent codes and published numbers.

    python test_physics.py

Each check prints PASS/FAIL with the measured deviation and its tolerance; the script exits
non-zero if anything fails. Checks that need a cached file (PRyMordial runs, results.json) are
skipped, not failed, when it is missing.
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402
from scipy import integrate, special  # noqa: E402

import config as C  # noqa: E402
import freezeout as fo  # noqa: E402
import thermo  # noqa: E402

RESULTS = []


def check(name: str, value: float, tol: float) -> None:
    ok = bool(np.isfinite(value) and abs(value) <= tol)
    RESULTS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<72} {value: .2e}  (tol {tol:.0e})")


def skip(name: str, why: str) -> None:
    print(f"  [skip] {name:<72} ({why})")


def test_equilibrium():
    print("Equilibrium yield and the Bessel function")
    for x in (5.0, 20.0):
        # n_eq = g/(2 pi^2) int p^2 exp(-E/T) dp in units T = 1, m = x
        n = integrate.quad(lambda p: p * p * math.exp(-math.sqrt(p * p + x * x)), 0, x + 60, epsabs=0, epsrel=1e-12)[0]
        n *= 2 / (2 * math.pi ** 2)
        s = 2 * math.pi ** 2 / 45 * 100.0
        check(f"Y_eq formula = direct phase-space integral (x = {x:g})", math.exp(fo.ln_yeq(x, 2.0, 100.0)) / (n / s) - 1, 1e-9)
    for x in (300.0, 999.0):
        check(f"ln K2 e^x: series agrees with SciPy below the switch (x = {x:g})",
              float(fo.ln_kve2(np.array([x]))[0]) - math.log(special.kve(2, x)), 1e-12)
    check("asymptotic series continuous at x = 1000", float(fo.ln_kve2(np.array([1000.0001]))[0] - math.log(special.kve(2, 1000.0))), 1e-7)


def test_plasma():
    print("The plasma: g*(T) tables and the clock")
    e = thermo.load("ss2020")
    check("Saikawa & Shirai g*s today = 3.931", e.g_s0 - 3.931, 1e-3)
    check("high-temperature limit g*rho(10 TeV) ~ 104 (interacting, below 106.75)", e.g_rho(1e4) - 104.2, 0.4)
    T, t = e.time_of(np.logspace(-6, 4, 3000))
    t1 = float(np.exp(np.interp(0.0, -np.log(T / 1e-3), np.log(t))))
    check("t(1 MeV) from S&S g* matches Escudero 2020 (0.7450 s)", t1 / 0.7450 - 1, 0.01)
    import neutrons as nu

    bg = nu.instant_background()
    check("instant decoupling: T_gamma/T_nu -> (11/4)^(1/3)", bg.extra["Tg_over_Tnu_final"] / (11 / 4) ** (1 / 3) - 1, 2e-4)


def test_solver():
    print("The Boltzmann solver")
    r = fo.Relic(100.0, a=2.11e-26)
    a, b = fo.solve(r, rtol=1e-8).omega_h2, fo.solve(r, rtol=1e-11).omega_h2
    check("Radau converged (rtol 1e-8 vs 1e-11)", a / b - 1, 1e-7)
    for m in (0.05, 100.0, 2e4):
        rad = fo.solve(fo.Relic(m, a=2.2e-26)).omega_h2
        imp = float(fo.omega_grid([m], 2.2e-26)[0])
        check(f"implicit-Euler grid stepper = Radau (m = {m:g} GeV)", imp / rad - 1, 2e-5)
    o1 = fo.solve(fo.Relic(100.0, a=2e-26)).omega_h2
    o2 = fo.solve(fo.Relic(100.0, a=4e-26)).omega_h2
    check("doubling <sigma v> roughly halves Omega (ratio - 0.5)", o2 / o1 - 0.5, 0.03)
    tab = np.loadtxt(C.RAW_DIR / "relic" / "ss2020_S_2.dat")
    k = int(np.argmin(np.abs(tab[:, 0] - 100.0)))
    sv = fo.thermal_cross_section(float(tab[k, 0]))["sigv"]
    check("thermal <sigma v> at 100 GeV = Saikawa & Shirai 2020 (rel.)", sv / (tab[k, 1] * 1e-26) - 1, 3e-3)
    pw = fo.solve(fo.Relic(100.0, b=1e-25)).omega_h2
    sw = fo.solve(fo.Relic(100.0, a=1e-25)).omega_h2
    # After freeze-out 1/Y grows by the integral of lambda: a/x_f for s-wave, 3b/x_f^2 for p-wave
    # (<sigma v> = 6b/x keeps falling), so equal coefficients leave ~x_f/3 ~ 7 times more behind.
    check("p-wave (b = a) leaves ~x_f/3 ~ 7x the s-wave relic (ratio - 7)", pw / sw - 7.0, 1.5)


def test_neutrons():
    print("Neutrons: weak rates and the network")
    import neutrons as nu

    check("decay phase space lambda0 = 1.636098 (PRyMordial's Fn_Born)", nu.lambda0() - 1.636098, 2e-6)
    for T in (3e-3, 1e-3, 0.4e-3):
        a, b = nu.weak_rates(T, T)
        check(f"detailed balance lambda_pn / lambda_np = e^(-Q/T) at T = {T * 1e3:g} MeV", b / a / math.exp(-C.Q_NP / T) - 1, 1e-8)
    a, _ = nu.weak_rates(2e-6, 1.4e-6)
    check("cold limit: lambda_np = 1/tau_n", a * C.TAU_N - 1, 1e-6)
    src = (C.RAW_DIR / "bbn" / "prym_init.py").read_text()
    for k, i, o in nu.REACTIONS[:5] + [nu.REACTIONS[7]]:
        m = re.search(rf"alpha_{k},beta_{k},gamma_{k} = ([^\n#]+)", src)
        ref = [float(eval(v)) for v in m.group(1).split(",")]  # noqa: S307 - our own downloaded constants
        ours = nu.reverse_coefficients(i, o)
        check(f"reverse-rate coefficient alpha for {k} = PRyMordial's", ours[0] / ref[0] - 1, 2e-4)
    cache = C.DERIVED_DIR / "prymordial_runs.json"
    res = C.RESULTS_PATH
    if cache.exists() and res.exists():
        r = json.loads(res.read_text())["neutrons"]
        if "born_small_instant_a" in r.get("prymordial", {}):
            p = r["prymordial"]["born_small_instant_a"]
            o = r["runs"]["escudero"]
            check("Y_p: our network = PRyMordial with the same physics (abs.)", o["Y_p"] - p["Y_p"], 3e-4)
            check("D/H: our network = PRyMordial with the same physics (rel.)", o["D_H"] / p["D_H"] - 1, 0.02)
    else:
        skip("network vs PRyMordial", "run analysis.py first")


def test_recombination_and_gamma():
    print("Recombination, dwarf galaxies, the Higgs portal")
    if C.RESULTS_PATH.exists():
        rc = json.loads(C.RESULTS_PATH.read_text())["recombination"]
        check("x_e(z = 200): Peebles atom = CAMB (RECFAST)", rc["x_e_200"] / rc["x_e_200_camb"] - 1, 0.03)
    import gamma_dwarfs as gd

    r = gd.reproduce_mcd_profile("Reticulum_2", "bb", masses_idx=(20,))
    check("McDaniel TS profile of Reticulum II rebuilt from its SED (near max)", r["rows"][0]["max_abs_dTS_near_max"], 0.05)
    import pdf_curves as pc

    ack = pc.ackermann2015(C.LIT_DIR / "1503.02641" / "p302_6yr_source_composite_limit_data_plot_bb.pdf")
    check("Ackermann 2015 figure: axis calibration residual (decades)", ack["calibration"]["max_tick_residual"], 1e-3)
    sdb = np.loadtxt(C.RAW_DIR / "relic" / "sdb12_m-sv-num.txt")
    t = ack["thermal_sdb12"]
    tm, ts = np.array(t["mass"]), np.array(t["sigmav"])
    sel = (tm > sdb[0, 0]) & (tm < sdb[-1, 0])
    ref = np.exp(np.interp(np.log(tm[sel]), np.log(sdb[:, 0]), np.log(sdb[:, 1]))) * 1e-26
    check("their thermal line read back = SDB12's table (max rel.)", float(np.max(np.abs(ts[sel] / ref - 1))), 0.02)
    import higgs_portal as hp

    check("SM Higgs width at 125.2 GeV from YR4 = 4.10 MeV", hp.gamma_h() / 4.10e-3 - 1, 0.01)
    m, lam = 1000.0, 0.3
    th = hp.thermal_average(m, lam, 3000.0)
    thr = hp.sigv_total(4 * m * m * (1 + 1e-8), m, lam)
    check("thermal average -> threshold value when cold (x = 3000)", th / thr - 1, 5e-3)
    sig = hp.sigma_si(1000.0, 0.31)
    check("sigma_SI(1 TeV, lambda = 0.31) = 8.7e-46 cm^2 (by hand)", sig / 8.72e-46 - 1, 5e-3)


def main() -> int:
    test_equilibrium()
    test_plasma()
    test_solver()
    test_neutrons()
    test_recombination_and_gamma()
    n_fail = RESULTS.count(False)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
