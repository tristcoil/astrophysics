"""Self-checks of the formulas against known answers.

    python test_physics.py          # all checks (~1-2 min; CAMB accuracy check included)
    python test_physics.py --quick  # skip the checks that run CAMB

No pytest needed (the shared venv has none). Each check prints PASS/FAIL.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import config as cfg  # noqa: E402
import masses as ms  # noqa: E402
import relic  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))


def close(a, b, rtol=0.0, atol=0.0):
    return abs(a - b) <= atol + rtol * abs(b)


# ---------------------------------------------------------------------------
# Relic neutrinos
# ---------------------------------------------------------------------------


def test_relic():
    print("Relic neutrinos")
    T = cfg.T_CMB_K
    n_g = relic.photon_number_density(T)
    check("photons per cm^3 at 2.7255 K = 410.7", close(n_g, 410.7, atol=0.1), f"{n_g:.2f}")
    t_nu = relic.neutrino_temperature_K(T)
    check("T_nu = (4/11)^(1/3) T_cmb = 1.9454 K", close(t_nu, 1.9454, atol=1e-4), f"{t_nu:.5f}")
    n_nu = relic.neutrino_number_density(t_nu)
    check("nu + nubar per flavour = 3/11 n_gamma = 112.0 per cm^3",
          close(n_nu, 3 / 11 * n_g, rtol=1e-12) and close(n_nu, 112.0, atol=0.05), f"{n_nu:.3f}")
    check("mean momentum = 3.1514 kT (Fermi-Dirac)", close(relic.MEAN_P_FD, 3.1514, atol=1e-4),
          f"{relic.MEAN_P_FD:.5f}")
    rc = relic.RHO_CRIT_H2
    rho_g = rc * relic.EV_J * 1e7 / (relic.C_CM_S**2)  # eV/cm^3 -> g/cm^3
    check("critical density = 1.87834e-29 h^2 g/cm^3 (PDG)", close(rho_g, 1.87834e-29, rtol=2e-5),
          f"{rho_g:.5e} g/cm^3 = {rc:.1f} h^2 eV/cm^3")
    f = relic.mass_to_density_factor(t_nu)
    check("Omega_nu h^2 = sum m / 94.1 eV (instantaneous decoupling)", close(f, 94.1, atol=0.1),
          f"{f:.2f}")
    # exact energy-density integral -> m n in the non-relativistic limit
    m = 1.0
    rho = relic.energy_density_massive(m, t_nu)
    check("rho of a 1 eV state -> m n (non-relativistic limit)",
          close(rho, m * n_nu, rtol=5e-6), f"ratio {rho / (m * n_nu):.7f}")
    # massless limit of the same integral
    rho0 = relic.energy_density_massive(1e-12, t_nu)
    check("rho of a massless state = 7/8 x 2 x pi^2/30 T^4", close(rho0, relic.massless_energy_density(t_nu), rtol=1e-9))
    # speed distribution normalises and matches the momentum quantiles
    v, pdf = relic.speed_distribution(0.05, t_nu, n=200001)
    norm = np.trapezoid(pdf, v)
    check("speed distribution integrates to 1 (m = 0.05 eV)", close(norm, 1.0, atol=2e-4), f"{norm:.6f}")
    # PTOLEMY: Long, Lunardini & Sabancilar 2014 find ~4 per year for 100 g (Dirac)
    rate = relic.tritium_capture_rate_per_year(cfg.TRITIUM_SIGMA_V_CM2, n_nu / 2, 100.0,
                                               majorana=False)
    check("PTOLEMY: ~4 captures/yr per 100 g of tritium (Dirac)", close(rate, 4.06, atol=0.1),
          f"{rate:.2f}")
    knr = relic.k_nonrelativistic(1.0, t_nu, 1.0)
    check("k_nr = 0.018 Omega_m^1/2 (m/eV)^1/2 h/Mpc", close(knr, 0.018, atol=0.0005),
          f"{knr:.4f}")


# ---------------------------------------------------------------------------
# Masses
# ---------------------------------------------------------------------------


def test_masses():
    print("Neutrino masses")
    o = cfg.NUFIT["NO"]
    i = cfg.NUFIT["IO"]
    s_no = ms.minimum_sum(o["dm21"][0], o["dm3l"][0], "NO")
    s_io = ms.minimum_sum(i["dm21"][0], i["dm3l"][0], "IO")
    check("minimum sum, normal ordering ~0.059 eV", 0.0585 < s_no < 0.0595, f"{s_no:.4f}")
    check("minimum sum, inverted ordering ~0.10 eV", 0.098 < s_io < 0.101, f"{s_io:.4f}")
    for order, p in (("NO", o), ("IO", i)):
        m1, m2, m3 = ms.masses(0.02, p["dm21"][0], p["dm3l"][0], order)
        d21 = m2**2 - m1**2
        d3 = m3**2 - (m1**2 if order == "NO" else m2**2)
        check(f"{order}: masses reproduce both splittings",
              close(d21, p["dm21"][0], rtol=1e-12) and close(d3, p["dm3l"][0], rtol=1e-12))
        total = ms.sum_masses(0.037, p["dm21"][0], p["dm3l"][0], order)
        back = ms.lightest_from_sum(total, p["dm21"][0], p["dm3l"][0], order)
        check(f"{order}: lightest mass from the sum round-trips", close(back, 0.037, atol=1e-12))
    # degenerate limit: m_beta -> m, sum -> 3 m
    mb = ms.m_beta(1.0, o["dm21"][0], o["dm3l"][0], "NO", o["s12sq"][0], o["s13sq"][0])
    check("m_beta -> m for quasi-degenerate masses", close(mb, 1.0, rtol=2e-3), f"{mb:.5f}")
    lo, hi = ms.m_betabeta_range(0.0, o["dm21"][0], o["dm3l"][0], "NO", o["s12sq"][0],
                                 o["s13sq"][0])
    check("NO with m1 = 0: m_bb between ~1 and ~4 meV", 0.0005 < lo < 0.003 and 0.003 < hi < 0.005,
          f"{lo * 1e3:.2f}-{hi * 1e3:.2f} meV")


# ---------------------------------------------------------------------------
# Decoupling
# ---------------------------------------------------------------------------


def test_decoupling():
    print("Neutrino decoupling")
    import decoupling as dc

    plasma = dc.Plasma()
    T = 1e3  # MeV: electrons ultra-relativistic (m_e/T clipped at 1e-3)
    p2 = plasma.qed(T, 2)[0]
    expect2 = -5.0 / 288.0 * dc.E2 * T**4
    check("O(e^2) QED pressure -> -(5/288) e^2 T^4 when m_e -> 0", close(p2, expect2, rtol=1e-5),
          f"ratio {p2 / expect2:.7f}")
    p3 = plasma.qed(T, 3)[0] - p2
    expect3 = dc.E2**1.5 * T**4 / (36 * np.sqrt(3) * np.pi)
    check("O(e^3) QED pressure -> e^3 T^4 / (36 sqrt3 pi) when m_e -> 0",
          close(p3, expect3, rtol=1e-5), f"ratio {p3 / expect3:.7f}")
    # the spline derivatives against finite differences, where the mass matters
    Tm = 0.4
    P, dP, d2P = plasma.qed(Tm, 3)
    h = 1e-4 * Tm
    fd1 = (plasma.qed(Tm + h, 3)[0] - plasma.qed(Tm - h, 3)[0]) / (2 * h)
    fd2 = (plasma.qed(Tm + h, 3)[0] - 2 * P + plasma.qed(Tm - h, 3)[0]) / h**2
    check("dP_int/dT and d2P_int/dT2 consistent with P_int at T = 0.4 MeV",
          close(dP, fd1, rtol=1e-6) and close(d2P, fd2, rtol=1e-4),
          f"{dP / fd1 - 1:.1e}, {d2P / fd2 - 1:.1e}")
    rho_e = plasma.rho_e(T)
    check("e+- energy density -> 7/8 x 4 x pi^2/30 T^4 when m_e -> 0",
          close(rho_e, 7 / 8 * 4 * np.pi**2 / 30 * T**4, rtol=1e-5))
    res = dc.solve(rates="off", qed_order=0, plasma=plasma)
    check("instantaneous decoupling: T_nu/T_gamma -> (4/11)^(1/3), N_eff -> 3",
          close(res["neff"], 3.0, atol=2e-5), f"N_eff = {res['neff']:.6f}")
    table = dc.load_rate_table()
    for label, rates, qed in cfg.DECOUPLING_RUNS:
        r = dc.solve(rates=rates, qed_order=qed, plasma=plasma, table=table)
        ref = cfg.ESCUDERO_2020_TABLE1[label]
        check(f"Escudero 2020 Table 1: {label} = {ref:.3f}", abs(r["neff"] - ref) <= 0.0006,
              f"ours {r['neff']:.5f}")
    for k, v in cfg.ESCUDERO_2020_A12.items():
        check(f"Escudero 2020 eq. A.12-13: {k} = {v}", close(r[k], v, atol=2e-5 if k != "neff"
                                                              else 1e-4), f"ours {r[k]:.5f}")


# ---------------------------------------------------------------------------
# Likelihood plumbing
# ---------------------------------------------------------------------------


def test_binning():
    print("Planck binning (D3's recipe)")
    import fetch_data as fd
    import likelihoods as lk

    tt = lk.BinnedTT(fd.load_planck_binned_tt())
    th = fd.load_planck_bestfit_theory()
    dl = np.zeros(tt.lmax + 1)
    dl[th["ell"].astype(int)] = th["tt"]
    model = tt.W @ dl
    rel = np.max(np.abs(model / tt.bestfit_binned - 1))
    check("binning Planck's best-fit theory reproduces Planck's binned best fit", rel < 1e-7,
          f"max rel diff {rel:.1e}")


def _theory_from_planck_file():
    import fetch_data as fd

    th = fd.load_planck_bestfit_theory()
    theory = {k: np.zeros(3001) for k in ("tt", "te", "ee", "bb", "pp")}
    ell = th["ell"].astype(int)
    m = ell <= 3000
    for k in theory:
        theory[k][ell[m]] = th[k][m]
    return theory


def test_likelihoods():
    print("Likelihoods against published chi-squared values")
    import fetch_data as fd
    import likelihoods as lk

    theory = _theory_from_planck_file()
    params = fd.load_planck_bestfit_params()
    lens = lk.PlanckLensing2018(fd.load_planck_lensing())
    r = lens.residual(theory, {"A_planck": params["calPlanck"]})
    check("Planck 2018 lensing at Planck's best fit: chi2 = 8.868", abs(r @ r - 8.868) < 0.02,
          f"ours {r @ r:.3f} (theory file rounded to 6 digits)")
    pl = lk.PlikLite(fd.load_plik_lite())
    r = pl.residual(theory, {"A_planck": 1.0})
    check("plik_lite TTTEEE chi2 at Planck's best fit ~ number of bins", 550 < r @ r < 640,
          f"{r @ r:.1f} for 613 bins")
    cut = lk.PlikLite(fd.load_plik_lite(), lmax_tt=999, lmax_pol=594)
    check("plik_lite with ACT's cut keeps 252 bins (TT < 1000, TE/EE < 600)", cut.n == 252,
          f"{cut.n}")
    a = params["calPlanck"]
    tt = lk.PlanckLowlTT(fd.load_lowl_tt())
    c = tt.chi2(theory, a)
    check("Planck low-l TT at Planck's best fit: chi2 = 23.2515", abs(c - 23.2515) < 0.002,
          f"ours {c:.4f}")
    r = tt.residual(theory, {"A_planck": a})
    check("low-l TT residuals reproduce its chi2", abs(r @ r + tt.chi2_constant - c) < 1e-8)
    ee = lk.PlanckLowlEE(fd.load_lowl_ee())
    c = ee.chi2(theory, a)
    check("Planck low-l EE (Planck's lookup) at Planck's best fit: chi2 = 396.0499",
          abs(c - 396.0499) < 0.002, f"ours {c:.4f}")
    cs = ee.chi2(theory, a, truncate=False)
    check("smooth low-l EE within 0.1 of the step-function lookup", abs(cs - c) < 0.1,
          f"{cs:.4f}")


def test_desi_bestfit(quick: bool):
    """DESI DR2's MAP for DESI BAO + Planck plik + ACT lensing, through our code."""
    if quick:
        return
    print("DESI DR2 best fit (base_mnu, plik + ACT lensing) through our likelihoods")
    import cmb_model as cm
    import fetch_data as fd
    import likelihoods as lk

    bf = fd.load_desi_bestfit("base_mnu/desi-bao-all_" + fd.CMB_VARIANTS["plik"])
    p = {"logA": bf["logA"], "ns": bf["ns"], "theta100": bf["theta_MC_100"],
         "ombh2": bf["ombh2"], "omch2": bf["omch2"], "tau": bf["tau"], "mnu": bf["mnu"]}
    th = cm.theory(p, "degenerate")
    d = cm.derived_dict(th)
    check("H0 from theta_MC agrees with DESI's (CAMB 2.0.4 vs 1.5.4)",
          close(d["H0"], bf["H0"], rtol=5e-5), f"{d['H0']:.4f} vs {bf['H0']:.4f}")
    check("r_drag agrees with DESI's", close(d["rdrag"], bf["rdrag"], rtol=2e-5),
          f"{d['rdrag']:.4f} vs {bf['rdrag']:.4f}")
    r = lk.BAO(fd.load_bao()).residual(th, p)
    check("DESI DR2 BAO chi2 at DESI's best fit", abs(r @ r - bf["chi2__BAO"]) < 0.02,
          f"ours {r @ r:.4f}, DESI {bf['chi2__BAO']:.4f}")
    c = lk.PlanckLowlTT(fd.load_lowl_tt()).chi2(th, bf["A_planck"])
    ref = bf["chi2__planck_2018_lowl.TT_clik"]
    check("Planck low-l TT chi2 at DESI's best fit", abs(c - ref) < 0.05,
          f"ours {c:.3f}, DESI {ref:.3f}")
    c = lk.PlanckLowlEE(fd.load_lowl_ee()).chi2(th, bf["A_planck"])
    ref = bf["chi2__planck_2018_lowl.EE_clik"]
    check("Planck low-l EE chi2 at DESI's best fit", abs(c - ref) < 0.3,
          f"ours {c:.3f}, DESI {ref:.3f}")
    if fd.ACT_LENS_REDUCED.exists():
        r = lk.ACTLensing(fd.load_act_lensing()).residual(th, p)
        ref = bf["chi2__act_dr6_lenslike_v1_2.ACTDR6LensLike"]
        check("ACT DR6 + Planck PR4 lensing chi2 at DESI's best fit", abs(r @ r - ref) < 0.1,
              f"ours {r @ r:.3f}, DESI {ref:.3f}")


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    quick = "--quick" in argv
    test_relic()
    test_masses()
    test_decoupling()
    test_binning()
    test_likelihoods()
    test_desi_bestfit(quick)
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
