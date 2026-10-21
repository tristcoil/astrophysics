"""Check the physics modules against answers known independently.

    python test_physics.py

Each check prints PASS or FAIL with the numbers compared; the script exits
non-zero if any check fails.
"""

from __future__ import annotations

import sys

import numpy as np
from scipy import integrate

import config as C
import nickel as N

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


def close(a, b, rtol=1e-6, atol=0.0):
    return bool(np.all(np.isclose(a, b, rtol=rtol, atol=atol)))


# --------------------------------------------------------------------------

def test_decay_chain():
    print("Radioactive decay chain 56Ni -> 56Co -> 56Fe")
    t = np.linspace(0, 800, 4001)
    n = N.nuclei(t, 1.0)
    check("Ni + Co + Fe is conserved", close(n["Ni"] + n["Co"] + n["Fe"], N.N_PER_MSUN, 1e-12))
    check("one solar mass of 56Ni is 2.14e55 nuclei", close(N.N_PER_MSUN, 2.1406e55, 1e-3),
          f"{N.N_PER_MSUN:.4e}")
    # Cobalt peaks when dN_Co/dt = 0: t = ln(lNi/lCo) / (lNi - lCo).
    t_peak = np.log(N.LAMBDA_NI / N.LAMBDA_CO) / (N.LAMBDA_NI - N.LAMBDA_CO)
    check("cobalt abundance peaks where the analytic formula says",
          abs(t[np.argmax(n["Co"])] - t_peak) < 0.25, f"{t_peak:.2f} d")
    check("mean lifetimes 8.76 d and 111.4 d", close([N.TAU_NI, N.TAU_CO], [8.764, 111.43], 1e-3),
          f"{N.TAU_NI:.3f}, {N.TAU_CO:.2f} d")


def test_power():
    print("Radioactive power")
    a, b = N.two_exponential_coefficients()
    check("A matches Nadyozhin (1994): 6.45e43 erg/s/Msun", close(a, 6.45e43, 0.01), f"{a:.4e}")
    check("B matches Nadyozhin (1994): 1.45e43 erg/s/Msun", close(b, 1.45e43, 0.01), f"{b:.4e}")
    t = np.linspace(0.0, 600.0, 301)
    q = N.power(t)["total"]
    q2 = a * np.exp(-t / N.TAU_NI) + b * np.exp(-t / N.TAU_CO)
    check("Bateman solution equals the two-exponential form", close(q, q2, 1e-10))

    # Energy released over all time = N0 (E_Ni + E_Co) per solar mass.
    tot, _ = integrate.quad(lambda x: N.power(x)["total"], 0, 5000, limit=400, points=[10, 50, 200])
    expect = N.N_PER_MSUN * (C.E_GAMMA_NI_MEV + C.E_GAMMA_CO_MEV + C.E_POSITRON_CO_MEV) * C.MEV_ERG
    check("time-integrated power = total decay energy (1.88e50 erg/Msun)",
          close(tot * C.DAY_S, expect, 1e-5), f"{tot * C.DAY_S:.4e} vs {expect:.4e}")

    check("cobalt's slope is 0.98 mag per 100 days",
          abs(N.CO_SLOPE_MAG_PER_100D - 0.974) < 0.002, f"{N.CO_SLOPE_MAG_PER_100D:.4f}")

    # Arnett's rule for a textbook Type Ia: 0.6 Msun and a 19-day rise gives
    # ~1.2e43 erg/s (Stritzinger & Leibundgut 2005: 2.0e43 erg/s per Msun at 19 d).
    q19 = N.power(19.0)["total"]
    check("power at day 19 is 2.0e43 erg/s per Msun", close(q19, 2.0e43, 0.03), f"{q19:.3e}")


def test_trapping():
    print("Gamma-ray trapping")
    check("t << T0: everything trapped", close(N.gamma_trapping(1.0, 40.0), 1.0, 1e-12))
    f = N.gamma_trapping(400.0, 40.0)
    check("t >> T0: deposited fraction -> (T0/t)^2", close(f, 0.01, 0.01), f"{f:.5f}")
    check("no T0 means complete trapping", close(N.gamma_trapping([5.0, 500.0], None), 1.0))
    late = N.deposition(2000.0, 30.0) / N.power(2000.0)["co_positron"]
    check("at very late times positrons dominate the heating", late < 1.1, f"ratio {late:.3f}")


def test_arnett():
    print("Arnett's diffusion model")
    t = np.linspace(0.0, 400.0, 8001)
    m, tau = 0.6, 15.0
    lum = N.arnett_light_curve(t, m, tau, None, dt=0.02)
    dep = m * N.deposition(t, None)
    i = np.argmax(lum)
    # In Arnett's model the light curve crosses the deposition curve exactly at
    # peak (dL/dt = 0 <=> L = D): "Arnett's rule".
    check("Arnett's rule holds at the model's own peak", close(lum[i], dep[i], 0.01),
          f"L={lum[i]:.4e}, D={dep[i]:.4e} at t={t[i]:.2f} d")
    check("late light curve follows the instantaneous deposition", close(lum[-1], dep[-1], 0.01),
          f"{lum[-1] / dep[-1]:.4f}")
    # dL/dt = (2t/tau^2)(D - L)  =>  Int t L dt = Int t D dt - (tau^2/2) L(T)
    lhs = np.trapezoid(t * lum, t)
    rhs = np.trapezoid(t * dep, t) - 0.5 * tau**2 * lum[-1]
    check("the Katz integral identity holds in the model", close(lhs, rhs, 2e-3), f"{lhs / rhs:.5f}")
    m_back = N.katz_ratio(t[1:], lum[1:], None)
    check("katz_ratio recovers the input nickel mass", close(m_back, m, 5e-3), f"{m_back:.4f}")
    m_rule = N.arnett_rule_mass(lum[i], t[i])
    check("arnett_rule_mass inverts the model's peak", close(m_rule, m, 0.01), f"{m_rule:.4f}")


def test_tail_fit():
    print("Tail fit recovers known inputs")
    t = np.linspace(60, 400, 40)
    lum = 0.55 * N.deposition(t, 38.0)
    r = N.fit_tail(t, lum)
    check("M_Ni recovered", close(r["m_ni"], 0.55, 1e-4), f"{r['m_ni']:.5f}")
    check("T0 recovered", close(r["t0"], 38.0, 1e-4), f"{r['t0']:.4f}")


def test_arnett_fit():
    print("Arnett fit with leakage recovers known inputs")
    t = np.concatenate([np.arange(2.0, 40.0, 1.0), np.arange(40.0, 150.0, 5.0)])
    lum = N.arnett_light_curve(t, 0.55, 14.0, 38.0, dt=0.05)
    r = N.fit_arnett(t, lum, 0.01 * lum, t_max_fit=150.0)
    check("M_Ni", close(r["m_ni"], 0.55, 2e-3), f"{r['m_ni']:.4f}")
    check("diffusion time", close(r["tau_m"], 14.0, 2e-3), f"{r['tau_m']:.3f}")
    check("T0", close(r["t0"], 38.0, 5e-3), f"{r['t0']:.2f}")


def test_extinction():
    import extinction as E

    print("Dust: Cardelli, Clayton & Mathis (1989)")
    check("A_lambda / A_V = 1 at 5495 A", close(E.ccm89(5495.0, 3.1), 1.0, 0.01), f"{E.ccm89(5495.0, 3.1)[0]:.4f}")
    for rv in (3.1, 1.4):
        ab, av = E.a_lambda(4405.0, 1.0, rv)[0], E.a_lambda(5495.0, 1.0, rv)[0]
        check(f"E(B-V) = A_B - A_V recovered for R_V = {rv}", close(ab - av, 1.0, 0.04), f"{ab - av:.3f}")
    check("low R_V means steeper extinction in the blue",
          E.ccm89(3600.0, 1.4)[0] > E.ccm89(3600.0, 3.1)[0])


def test_gamma_rays():
    import config as C
    import gammarays as G

    print("Gamma-ray lines")
    check("Klein-Nishina -> Thomson at low energy", close(G.klein_nishina(0.1), C.SIGMA_T, 1e-3))
    check("Klein-Nishina at 1 MeV is 0.3175 sigma_T", close(G.klein_nishina(1000.0) / C.SIGMA_T, 0.3175, 1e-3),
          f"{G.klein_nishina(1000.0) / C.SIGMA_T:.4f}")
    ej = G.Ejecta("uniform", m_ni_msun=0.6)
    p_early, p_late = ej.escape_probability([5.0, 2000.0], 847.0)
    check("opaque early, transparent late", p_early < 0.01 and p_late > 0.99, f"{p_early:.2e}, {p_late:.4f}")
    # Uniform sphere, central point source: escape = exp(-kappa rho R).
    sph = G.Ejecta("uniform", m_ni_msun=1e-4, m_core_msun=0.0)
    t = 60.0 * C.DAY_S
    kappa = G.line_opacity(847.0)
    tau = kappa * (sph.m / (4 / 3 * np.pi * (sph.v_scale * t) ** 3)) * sph.v_scale * t
    check("uniform sphere: central escape = exp(-tau)", close(sph.escape_probability(60.0, 847.0)[0], np.exp(-tau), 0.01),
          f"{sph.escape_probability(60.0, 847.0)[0]:.4f} vs {np.exp(-tau):.4f}")
    d = 3.0e25
    f = G.line_flux(80.0, 0.5, 847, d)[0]
    n_co = N.nuclei(80.0, 0.5)["Co"]
    check("transparent line flux = decays x yield / 4 pi d^2",
          close(f, N.LAMBDA_CO / C.DAY_S * n_co * 0.999399 / (4 * np.pi * d * d), 1e-9))


def test_photometry():
    import photometry as P

    print("Synthetic photometry")
    wave = np.linspace(3000.0, 10000.0, 7001)
    flat_fnu = 3631e-23 * P.C_ANG / wave**2  # an AB magnitude-zero source
    for band in ("ztf_g", "ztf_r", "B"):
        check(f"flat f_nu of 3631 Jy is AB 0 in {band}", abs(P.synthetic_mag(wave, flat_fnu, band, "AB")) < 1e-3)
    z = 0.04
    k = P.synthetic_mag(wave * (1 + z), flat_fnu / (1 + z), "ztf_g") - P.synthetic_mag(wave, flat_fnu, "ztf_g")
    check("K-correction of a flat f_nu source is -2.5 log10(1+z)", close(k, -2.5 * np.log10(1 + z), 1e-3),
          f"{k:.4f} vs {-2.5 * np.log10(1 + z):.4f}")
    check("mag_to_flam(0, AB) = 3631 Jy at the pivot",
          close(P.mag_to_flam(0.0, "ztf_g", "AB"), 3631e-23 * P.C_ANG / P.pivot("ztf_g") ** 2, 1e-12))


def test_bolometric_integration():
    import bolometric as B
    import config as C

    print("Bolometric integration of a blackbody")
    temp = 15000.0
    wave = np.array([2055.0, 2246.0, 2581.0, 3585.0, 4371.0, 5478.0, 6498.0, 8020.0, 12393.0, 16495.0, 21639.0])
    shape = B.planck_lam(wave, np.array([temp]))[:, 0]
    exact = C.SIGMA_SB * temp**4 / np.pi  # Int B_lambda over all wavelengths (per steradian)
    # planck_lam drops 2hc^2 and uses Angstrom; restore both for an absolute check.
    flam = shape * 2 * C.H_CGS * C.C_CGS**2 * 1e40 * 1e-8
    out = B.Bolometric._integrate(wave, flam[None, :], hot_uv=True)
    total = out["uvoi"][0] + out["nir"][0]
    check("UV-to-Ks trapezoid + tails recovers sigma T^4 within 5%", close(total, exact, 0.05),
          f"ratio {total / exact:.3f}")


def test_gp():
    from gp import LightCurveGP

    print("Light-curve GP")
    rng = np.random.default_rng(4)
    t = np.sort(rng.uniform(2, 200, 120))
    truth = lambda x: 15 + 2.5 * np.log10(1 + (x / 18) ** 2) + 0.01 * x
    y = truth(t) + rng.normal(0, 0.03, t.size)
    g = LightCurveGP(t, y, np.full(t.size, 0.03))
    tt = np.linspace(20, 190, 50)  # where the points are dense
    m, s = g.predict(tt)
    check("interpolation within 0.05 mag where sampled", np.max(np.abs(m - truth(tt))) < 0.05,
          f"max error {np.max(np.abs(m - truth(tt))):.3f} mag")
    tt = np.linspace(3, 199, 200)  # including the sparse first days
    m, s = g.predict(tt)
    check("errors stay within 3 sigma of the GP's own uncertainty", np.max(np.abs(m - truth(tt)) / s) < 3.0,
          f"max |error|/sigma = {np.max(np.abs(m - truth(tt)) / s):.2f}")
    check("coverage flags a gap and the edges", (not g.covered([400.0])[0]) and g.covered([100.0])[0])


def main() -> int:
    for fn in (test_decay_chain, test_power, test_trapping, test_arnett, test_tail_fit, test_arnett_fit,
               test_extinction, test_gamma_rays, test_photometry, test_bolometric_integration, test_gp):
        fn()
        print()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
