"""Self-checks for the physics and statistics, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas and the fitting code, not the data.

Every FRB number in this project passes through these functions. A slip in the
dispersion constant, a factor of (1+z) in the wrong place, or a probability
density that does not integrate to one would all produce plots that look
perfectly reasonable and are wrong -- hence checking each against something
independent: a closed-form limit, a hand calculation from raw constants, or a
published value.

The last two checks feed *simulated* numbers through the fitting code to confirm
it recovers the parameters that generated them. They are labelled as such;
nothing simulated reaches a figure.
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
from astropy.cosmology import FlatLambdaCDM, Planck18
from scipy import special

warnings.filterwarnings("ignore")

import analysis
import frb_physics as fp
import macquart_fit as mf
import waterfalls as wf

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "", absolute: bool = False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if absolute or expected == 0 else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    kind = "abs" if absolute or expected == 0 else "rel"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<52} {got:>12.6g} vs {expected:>11.6g} {unit:<9}"
          f"({kind} err {err:.1e}, tol {tol:.0e})")


def truth(label: str, condition: bool) -> None:
    if not condition:
        FAILURES.append(label)
    print(f"  [{'ok  ' if condition else 'FAIL'}] {label}")


def main() -> int:
    # ------------------------------------------------------------------ dispersion
    print("Dispersion constant")
    # e^2/(2 pi m_e c) = 4.148808 ms GHz^2 pc^-1 cm^3 is the textbook value
    # (Lorimer & Kramer, Handbook of Pulsar Astronomy, eq. 4.7).
    check("k_DM from constants (s MHz^2 pc^-1 cm^3)", fp.K_DM_EXACT, 4148.808, 2e-6)
    check("convention 1/2.41e-4 exceeds exact by", fp.K_DM_CONVENTION / fp.K_DM_EXACT - 1, 1.377e-4, 0.01)

    print("\nDispersion delays")
    # Hand calculation: DM 100 across CHIME's band, (1/2.41e-4) * 100 * (400^-2 - 800^-2) = 1.945 s.
    check("DM 100, 400 vs 800 MHz", float(fp.dispersion_delay(100, 400, 800)), 1.9450208, 1e-6, "s")
    # The Lorimer burst (FRB 010724, DM 375): Lorimer et al. 2007 show the sweep
    # across 1.23-1.52 GHz taking ~0.35 s.
    check("Lorimer burst sweep 1230-1520 MHz", float(fp.dispersion_delay(375, 1230, 1520)), 0.355, 0.02, "s")
    truth("low frequencies arrive later (positive delay)", fp.dispersion_delay(300, 400, 800) > 0)

    # ------------------------------------------------------------------ cosmology
    print("\nMean cosmic DM (the Macquart relation)")
    check("electrons per m_p, Y = 0.25", fp.electrons_per_baryon_mass(0.25), 0.875, 1e-12)
    # Prefactor 3 c H0 Omega_b (1 - Y/2) / (8 pi G m_p) from raw CGS constants.
    c, G, m_p, pc = 2.99792458e10, 6.67430e-8, 1.67262192e-24, 3.0856775814913673e18
    h0 = Planck18.H0.value * 1e5 / (1e6 * pc)
    by_hand = 3 * c * h0 * Planck18.Ob0 * (1 - 0.2454 / 2) / (8 * np.pi * G * m_p) / pc
    check("DM scale for Planck 2018, all baryons diffuse", fp.cosmic_dm_prefactor(), by_hand, 1e-6, "pc cm^-3")
    # Einstein-de Sitter: E(z) = (1+z)^1.5, so I(z) = 2 (sqrt(1+z) - 1) exactly.
    eds = FlatLambdaCDM(H0=70, Om0=1.0, Ob0=0.05)
    for zz in (0.5, 1.0, 3.0):
        check(f"I(z={zz}) matter-only closed form", float(fp.dm_integral(zz, eds)), 2 * (np.sqrt(1 + zz) - 1), 2e-6)
    # Near z = 0, <DM> -> f_d * prefactor * z.
    check("low-z slope, f_d = 0.844", float(fp.mean_dm_cosmic(1e-3)) / 1e-3,
          0.844 * fp.cosmic_dm_prefactor(), 2e-3, "pc cm^-3")
    for zz in (0.05, 0.5, 1.0, 2.5):
        check(f"z -> DM -> z round trip at z = {zz}", float(fp.z_from_mean_dm(fp.mean_dm_cosmic(zz))), zz, 1e-5)

    # ------------------------------------------------------------------ Macquart PDF
    print("\nMacquart et al. (2020) DM_cosmic distribution")
    d = np.geomspace(1e-7, 1e10, 3_000_001)
    for s in (0.02, 0.3, 2.0, 10.0, 30.0):
        p = fp.macquart_pdf_delta(d, s)
        check(f"sigma = {s}: normalisation", np.trapezoid(p, d), 1.0, 2e-4)
        check(f"sigma = {s}: mean Delta", np.trapezoid(p * d, d), 1.0, 2e-4)
    p = fp.macquart_pdf_delta(d, 0.02)
    check("small sigma -> Gaussian of width sigma", np.sqrt(np.trapezoid(p * (d - 1) ** 2, d)), 0.02, 0.02)
    # Large sigma: C0 -> -inf and p -> Delta^-3 exp(-lambda Delta^-3), i.e.
    # v = lambda Delta^-3 ~ Gamma(2/3). With <Delta> = 1 the median is
    # (Gamma(2/3)/Gamma(1/3)) * [median of Gamma(2/3)]^(-1/3).
    med_limit = special.gamma(2 / 3) / special.gamma(1 / 3) * special.gammaincinv(2 / 3, 0.5) ** (-1 / 3)
    p = fp.macquart_pdf_delta(d, 38.0)
    cdf = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(d))])
    check("large sigma: median -> analytic limit", float(np.interp(0.5, cdf, d)), med_limit, 2e-3)
    x = np.geomspace(1e-3, 1e6, 400_001)
    check("lognormal integrates to 1", np.trapezoid(fp.lognormal_pdf(x, np.log(120), 0.9), x), 1.0, 1e-5)

    print("\nLikelihood code")
    # The vectorised grid likelihood must equal the direct convolution.
    axes = {"f_d": np.array([0.7, 1.1]), "log10_F": np.array([np.log10(0.3)]),
            "mu": np.array([np.log(90.0)]), "sigma_host": np.array([0.8])}
    for z_, dm_ in ((0.3, 420.0), (0.05, 160.0)):
        grid = np.exp(mf.log_likelihood_grid([dm_], [z_], axes, sigma_mw=0.0))
        for i, f_d in enumerate(axes["f_d"]):
            direct = fp.pdf_dm_extragalactic(dm_, z_, f_d, 0.3, np.log(90.0), 0.8)[0]
            check(f"grid = direct, z={z_}, DM={dm_:.0f}, f_d={f_d}", float(grid[i, 0, 0, 0]), float(direct), 1e-9)
    dm_axis = np.linspace(1, 4000, 4000)
    total = np.trapezoid(fp.pdf_dm_extragalactic(dm_axis, 0.4, 0.84, 0.3, np.log(100), 0.9), dm_axis)
    check("p(DM_EG | z) integrates to 1", total, 1.0, 5e-3)
    # With the Milky Way term switched on, the density still integrates to one.
    axes1 = {"f_d": np.array([0.84]), "log10_F": np.array([np.log10(0.3)]),
             "mu": np.array([np.log(100.0)]), "sigma_host": np.array([0.9])}
    dm_axis2 = np.linspace(-150, 3000, 900)
    dens = [float(np.exp(mf.log_likelihood_grid([v], [0.4], axes1, sigma_mw=20.0))[0, 0, 0, 0]) if v > -80 else 0.0
            for v in dm_axis2]
    check("with Milky Way scatter: still integrates to 1", np.trapezoid(dens, dm_axis2), 1.0, 1e-2)

    # ------------------------------------------------------------------ sky helpers
    print("\nSky geometry and statistics")
    check("north Galactic pole has b = +90", float(analysis._galactic_latitude(192.85948, 27.12825)), 90.0, 1e-6,
          "deg", absolute=True)
    # The origin of Galactic coordinates (RA 17h45m37.2s, Dec -28d56m10s, J2000).
    check("Galactic coordinate origin has b = 0", float(analysis._galactic_latitude(266.40499, -28.93617)), 0.0, 0.01,
          "deg", absolute=True)
    even = np.arange(0, 360, 1.0)
    truth("Rayleigh: evenly spaced angles are not clustered (p > 0.9)", analysis._rayleigh(even, 1)["pvalue"] > 0.9)
    truth("Rayleigh: angles within 30 deg of each other are (p < 1e-10)",
          analysis._rayleigh(np.linspace(0, 30, 200), 1)["pvalue"] < 1e-10)

    # ------------------------------------------------------------------ simulated recoveries
    print("\nRecovery tests on SIMULATED numbers (code checks only; nothing here is data)")
    rng = np.random.default_rng(1)
    # Dispersion law: exact t = t0 + k DM f^-2 plus 0.05 ms noise must give alpha = -2.
    f = np.linspace(405, 795, 60)
    arr = {"freq_mhz": f, "t_arrival_s": 0.3 + fp.K_DM_CONVENTION * 400.0 * f ** -2 + rng.normal(0, 5e-5, f.size),
           "t_err_s": np.full(f.size, 5e-5), "k_dm": fp.K_DM_CONVENTION, "dm_chime": 400.0}
    fit = wf.fit_dispersion_law(arr)
    check("dispersion index recovered (3 sigma)", fit["alpha"], -2.0, 3 * fit["alpha_err"], absolute=True)
    check("DM recovered with index fixed", fit["dm_alpha_fixed"], 400.0, 5 * fit["dm_alpha_fixed_err"], absolute=True)

    # Macquart fit: draw 250 bursts from the model itself and refit f_d.
    true = {"f_d": 0.84, "F": 0.3, "mu": np.log(100.0), "sigma_host": 0.8}
    zs = rng.uniform(0.05, 1.0, 250)
    delta_grid = np.geomspace(1e-3, 50, 20001)
    dm_eg = []
    for z_ in zs:
        p = fp.macquart_pdf_delta(delta_grid, true["F"] / np.sqrt(z_))
        cdf = np.concatenate([[0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(delta_grid))])
        cosmic = np.interp(rng.uniform(0, cdf[-1]), cdf, delta_grid) * fp.mean_dm_cosmic(z_, true["f_d"])
        host = np.exp(rng.normal(true["mu"], true["sigma_host"])) / (1 + z_)
        dm_eg.append(cosmic + host)
    axes = mf.grid_axes({"f_d": (0.5, 1.2, 36), "log10_F": (-0.9, 0.0, 4), "mu": (3.9, 5.3, 8),
                         "sigma_host": (0.4, 1.3, 6)})
    res = mf.fit([f"sim{i}" for i in range(zs.size)], zs, np.array(dm_eg), axes=axes, sigma_mw=0.0, n_conv=250)
    q = res.summary["f_d"]
    half_width = 0.5 * (q["q840"] - q["q160"])
    print(f"         simulated f_d = {true['f_d']}; recovered {q['q500']:.3f} "
          f"(68% interval {q['q160']:.3f}-{q['q840']:.3f})")
    check("f_d recovered within 2 sigma", q["q500"], true["f_d"], 2 * half_width, absolute=True)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
