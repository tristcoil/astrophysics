"""Self-tests of the physics: run with  python test_physics.py  (no pytest needed).

Each test checks a piece against an exact result or a published number, so a change that breaks
the physics shows up here before it reaches a figure.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import config as C  # noqa: E402

PI2 = np.pi ** 2
FAILED = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILED.append(name)


def test_ideal_gas():
    import thermo

    for st in (+1, -1):
        f = thermo.functions(np.array([1e-7]), st)
        r0, p0, n0 = thermo.MASSLESS[st]
        check(f"massless limit, stat {st:+d}", abs(f["rho"][0] / r0 - 1) < 1e-9 and abs(f["p"][0] / p0 - 1) < 1e-9)
        x = np.array([0.01, 0.3, 1.0, 4.0, 20.0])
        r, p = thermo.bessel_series(x, st, 3000)
        g = thermo.functions(x, st)
        err = max(np.max(np.abs(g["rho"] / r - 1)), np.max(np.abs(g["p"] / p - 1)))
        check(f"tabulated integrals = Bessel series, stat {st:+d}", err < 1e-8, f"max rel. diff {err:.1e}")
    # thermodynamic consistency: s = dp/dT for a massive species
    sp = thermo.Species("x", "x", "test", 1.0, 2, -1)
    T = np.array([0.30, 0.30 * (1 + 1e-5)])
    d = thermo.gas(sp, T)
    dpdT = (d["p"][1] - d["p"][0]) / (T[1] - T[0])
    check("s = dp/dT for a massive fermion", abs(dpdT / d["s"][0] - 1) < 1e-4)


def test_standard_model_count():
    import particles
    import thermo

    sm = particles.standard_model()
    g = sum(sp.g_eff_massless for sp in sm)
    check("Standard Model count = 106.75", abs(g - 106.75) < 1e-12, f"{g}")
    T = np.array([1e6])
    g_rho = sum(thermo.gstar(sp, T)["rho"][0] for sp in sm)
    check("g*rho(1 PeV) of the free gas -> 106.75", abs(g_rho - 106.75) < 1e-3, f"{g_rho:.5f}")
    check("photons + e+- + 3 nu = 10.75", 2 + 7 / 8 * (4 + 6) == 10.75)


def test_hadron_gas():
    import particles
    import qcd

    hs = particles.hadrons()
    pids = {h.pdgid for h in hs}
    check("K_S and K_L excluded, K0 and anti-K0 kept", 130 not in pids and 310 not in pids and 311 in pids and -311 in pids)
    check("every hadron has its antiparticle or is its own",
          all((-h.pdgid in pids) or (h.pdgid in (111, 113, 221, 223, 331, 333, 443, 553)) or
              _self_conjugate(h.pdgid) for h in hs))
    # at very low T only pions matter: 3 bosons
    T = np.array([0.02])
    h = qcd.hrg(T, hs)
    pi = qcd.hrg(T, [x for x in hs if x.pdgid in (111, 211, -211)])
    check("hadron gas at 20 MeV = pion gas", abs(h["p"][0] / pi["p"][0] - 1) < 2e-3, f"{h['p'][0] / pi['p'][0]:.5f}")
    # HotQCD 2014 pressure at 130 MeV (0.439 +0.065 -0.044) against our hadron gas
    p130 = qcd.hrg(np.array([0.130]), hs)["p"][0]
    check("hadron gas p/T^4 at 130 MeV within HotQCD's band", 0.439 - 0.044 <= p130 <= 0.439 + 0.065, f"{p130:.3f}")


def _self_conjugate(pid):
    from particle import Particle

    return Particle.from_pdgid(pid).is_self_conjugate


def test_lattice_tables():
    import qcd

    hot, wb = qcd.hotqcd2014(), qcd.wb2014()
    check("HotQCD 2014 table parsed (130-400 MeV)", len(hot["T"]) == 55 and abs(hot["T"][0] - 0.130) < 1e-9
          and abs(hot["p"][0] - 0.439) < 1e-9 and abs(hot["p_lo"][0] - 0.044) < 1e-9)
    # thermodynamic identity within the published table: s/T^3 = (e + p)/T^4
    check("HotQCD table: s/T^3 = e/T^4 + p/T^4 to rounding", np.max(np.abs(hot["s"] - hot["rho"] - hot["p"])) < 0.02)
    check("WB 2014 table: s = e + p", np.max(np.abs(wb["s"] - wb["rho"] - wb["p"]) / wb["s"]) < 0.02)
    c = qcd.borsanyi2016_params()["2+1+1"]
    check("Borsanyi 2016 2+1+1 coefficients read from the TeX", c["h0"] == 0.353 and c["g2"] == 0.160)


def test_strong_coupling():
    import qcd

    a = float(qcd.alpha_s(C.M_Z, 5))
    check("alpha_s(M_Z) = 0.1180", abs(a - 0.1180) < 1e-9)
    a_tau = float(qcd.alpha_s(1.777, 3))
    check("alpha_s(m_tau) ~ 0.31-0.33 (PDG tau-decay value 0.314 +- 0.014)", 0.30 < a_tau < 0.34, f"{a_tau:.3f}")
    # the g^2 term of the pressure: p = p_SB (1 - (15/4)(1+5nf/12)/(1+21nf/32) alpha/pi + ...) at tiny alpha
    T = np.array([1e8])
    p = qcd.p_pqcd_massless(T, 6)[0]
    check("pQCD pressure below Stefan-Boltzmann at 10^8 GeV, within 5%",
          0.95 < p / qcd.p_sb(6) < 1.0, f"{p / qcd.p_sb(6):.4f}")


def test_qed_corrections():
    import neutrinos as N

    q = N.qed(np.array([1.0]))
    expect = -5 * N.E2 / 288 + N.E3 / (36 * np.sqrt(3) * np.pi)
    check("QED pressure: -5e^2/288 + e^3/(36 sqrt3 pi) for T >> m_e", abs(q["P"][0] / expect - 1) < 1e-6)
    check("QED energy density = 3 P for T >> m_e", abs(q["rho"][0] / (3 * q["P"][0]) - 1) < 1e-4)


def test_neutrino_decoupling():
    import neutrinos as N

    r = N.solve(collisions=False, qed_order=0, T_start=0.3)
    check("instantaneous decoupling: T_gamma/T_nu = (11/4)^(1/3)", abs(r["ratio_e"] / (11 / 4) ** (1 / 3) - 1) < 2e-6,
          f"{r['ratio_e']:.7f}")
    g = N.low_T_gstar(r)
    check("then g*rho = 3.3626 and g*s = 3.9091", abs(g["g_rho"] - 3.36264) < 2e-4 and abs(g["g_s"] - 3.90909) < 2e-4)
    r1 = N.solve(rates="fd", qed_order=3, separate=False)
    check("Escudero 2020 Table 1 row 'FD + NLO-QED': ratio 1.3954, N_eff 3.049",
          abs(r1["ratio_common"] - 1.3954) < 6e-5 and abs(r1["N_eff"] - 3.049) < 6e-4,
          f"{r1['ratio_common']:.5f}, {r1['N_eff']:.4f}")
    r2 = N.solve(rates="fd_me", qed_order=3, separate=True)
    g2 = N.low_T_gstar(r2)
    check("Escudero+2025 'no oscillations': N_eff 3.0443, g*rho 3.3828, g*s 3.9302",
          abs(r2["N_eff"] - 3.0443) < 6e-5 and abs(g2["g_rho"] - 3.3828) < 6e-5 and abs(g2["g_s"] - 3.9302) < 6e-5,
          f"{r2['N_eff']:.5f}, {g2['g_rho']:.5f}, {g2['g_s']:.5f}")
    check("within 0.0005 of Bennett et al. 2021 (3.0440)", abs(r2["N_eff"] - 3.0440) < 5e-4)


def test_clock():
    import cosmic_time as ct

    # constant g: t = 1/(2H) exactly; check the integrator on a pure power law
    u = np.linspace(-3, 3, 601)
    num = ct._cumint_from_top(u, np.ones_like(u), -2.0)
    exact = (np.exp(-2 * u) - np.exp(-2 * u[-1])) / 2
    check("exact-exponential integrator", np.max(np.abs(num[:-1] / exact[:-1] - 1)) < 1e-12)
    # the classic t = 2.42 s / sqrt(g) / T_MeV^2 for constant g
    g, T = 10.75, 1e-3
    H = np.sqrt(8 * np.pi / 3 * PI2 / 30 * g * T ** 4) / C.M_PL
    t = 1 / (2 * H) * C.HBAR_GEV_S
    check("t(1 MeV, g = 10.75) = 2.42 s / sqrt(10.75)", abs(t / (2.42 / np.sqrt(10.75)) - 1) < 2e-3, f"{t:.4f} s")


def test_results_file():
    import json

    if not C.RESULTS_PATH.exists():
        print("  skip results checks (run analysis.py first)")
        return
    r = json.loads(C.RESULTS_PATH.read_text())
    age = r["calendar"]["age_today_Gyr"]
    check("age today 13.787 Gyr (astropy Planck18: 13.7869)", abs(age - 13.7869) < 0.002, f"{age}")
    check("our g* inside Saikawa & Shirai's band on > 75% of the grid",
          r["comparisons"]["fraction_of_grid_within_SS_error"] > 0.75)
    mh = r["calendar"]["byrnes2018_MH_formula_check"]
    check("horizon mass matches Byrnes et al.'s formula to 0.5%",
          all(abs(a / b - 1) < 5e-3 for a, b in zip(mh["ours_Msun"], mh["formula_Msun"])))


if __name__ == "__main__":
    for fn in (test_ideal_gas, test_standard_model_count, test_hadron_gas, test_lattice_tables,
               test_strong_coupling, test_qed_corrections, test_neutrino_decoupling, test_clock,
               test_results_file):
        print(fn.__name__)
        fn()
    print(f"\n{'all passed' if not FAILED else str(len(FAILED)) + ' FAILED: ' + ', '.join(FAILED)}")
    sys.exit(1 if FAILED else 0)
