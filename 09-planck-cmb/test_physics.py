"""Self-checks for the physics module, against cases with known answers.

Run with ``python test_physics.py``. No test framework and no network: these
check the formulas, not the data.

The strongest checks are the last group. Planck publishes, alongside its
best-fit model, derived quantities computed by the Boltzmann code CAMB: the
redshift of last scattering, the sound horizon, the distance to last scattering.
Feeding Planck's best-fit densities through our own few-line formulas has to
land on the same numbers -- if it does, the physics in ``cmb_physics`` is the
same physics the professionals use, not an approximation of it.
"""

from __future__ import annotations

import sys

import numpy as np

import cmb_physics as cp

FAILURES: list[str] = []

#: Planck 2018 baseline best fit (TT,TE,EE+lowE+lensing), from
#: COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum_R3.01.txt.
#: Hard-coded so these tests need no downloaded data.
PLANCK_BF = {
    "ombh2": 0.0223828, "omch2": 0.1201075, "h": 0.6732117, "omnuh2": 0.0006451439,
    "zstar": 1089.914, "rstar": 144.3938, "thetastar": 1.041085e-2,
    "DMstar": 13869.55, "zeq": 3405.080, "zdrag": 1059.971, "rdrag": 147.0493,
}


def check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    """Assert ``got`` is within ``tol`` (fractional) of ``expected``."""
    rel = abs(got - expected) / abs(expected) if expected else abs(got)
    ok = rel <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<44} {got:>13.6g} vs {expected:>12.6g} "
          f"{unit:<7} ({rel * 100:.3f}% off, tol {tol * 100:g}%)")


def main() -> int:
    t0 = 2.7255

    print("Blackbody radiation")
    # 1 cm^-1 is 29.979 GHz, by definition of c.
    check("1 cm^-1 in GHz", float(cp.wavenumber_to_ghz(1.0)), 29.9792458, 1e-12, "GHz")
    # Wien's law in frequency form: 58.789 GHz per kelvin.
    check("Wien peak at 2.7255 K", cp.wien_peak_ghz(t0), 58.78925757 * t0, 1e-6, "GHz")
    # B_nu must be maximal exactly at the Wien peak.
    nu = np.linspace(100, 250, 150001)
    check("B_nu numerical maximum", float(nu[np.argmax(cp.planck_intensity(nu, t0))]),
          cp.wien_peak_ghz(t0), 1e-5, "GHz")
    # Normalisation, independently: integrated over all frequencies, a
    # blackbody surface radiates sigma T^4 (the Stefan-Boltzmann law).
    grid = np.linspace(1e-3, 5000, 2_000_001)
    total = np.pi * np.trapezoid(cp.planck_intensity(grid, t0), grid * 1e9) * cp.MJY_PER_SR
    import astropy.constants as const
    check("pi * integral of B_nu = sigma T^4", total, const.sigma_sb.si.value * t0**4, 1e-6,
          "W/m^2")
    # 411 photons in every cubic centimetre of the universe.
    check("photon number density", cp.photon_number_density(t0), 410.7, 1e-3, "cm^-3")
    # dB/dT must match a finite difference of B.
    nu0, dT = 150.0, 1e-4
    fd = (cp.planck_intensity(nu0, t0 + dT) - cp.planck_intensity(nu0, t0 - dT)) / (2 * dT)
    check("dB/dT vs finite difference", float(cp.planck_intensity_dT(nu0, t0)), float(fd),
          1e-7, "MJy/sr/K")
    # A y-distortion conserves photon number and changes sign at 217 GHz.
    nus = np.linspace(150, 300, 300001)
    zero = nus[np.argmin(np.abs(cp.y_distortion(nus, t0)))]
    check("y-distortion null frequency", float(zero), 217.5, 3e-3, "GHz")

    print("\nDensities")
    # The critical density is 1.878e-26 h^2 kg/m^3 (about 5 protons per m^3 for h=0.7).
    check("rho_crit / h^2", cp.RHO_CRIT_H2, 1.87847e-26, 1e-4, "kg/m^3")
    check("Omega_gamma h^2 at 2.7255 K", cp.omega_gamma_h2(t0), 2.4728e-5, 1e-3)
    check("T_nu / T_gamma", cp.NEUTRINO_TEMPERATURE_RATIO, 0.713765, 1e-5)
    # Matter-radiation equality for the Planck best fit. The massive neutrino is
    # still relativistic at z_eq, so it counts as radiation, not matter.
    zeq = cp.z_equality(PLANCK_BF["ombh2"] + PLANCK_BF["omch2"], t0)
    check("z_eq, Planck 2018 best fit", zeq, PLANCK_BF["zeq"], 2e-3)

    print("\nPower spectrum conventions")
    ell = np.arange(2, 3000)
    cl = 1.0 / ell**2
    check("D_l -> C_l -> D_l round trip", float(np.max(np.abs(
        cp.cl_from_dl(ell, cp.dl_from_cl(ell, cl)) / cl - 1))) + 1.0, 1.0, 1e-12)
    # Five modes at l = 2: the quadrupole is uncertain to 63% however good the
    # instrument.
    check("cosmic variance at l = 2", float(cp.cosmic_variance_fraction(2)), np.sqrt(0.4), 1e-12)
    check("cosmic variance at l = 1000, f_sky 0.5", float(cp.cosmic_variance_fraction(1000, 0.5)),
          np.sqrt(2 / (2001 * 0.5)), 1e-12)
    check("l = 200 is about a degree", float(cp.multipole_to_degrees(200)), 0.9, 1e-12, "deg")
    # For a flat D_l the map variance is D * sum (2l+1)/(2 l (l+1)), which
    # telescopes: sum (1/l + 1/(l+1)) / 2.
    lmax = 1000
    ls = np.arange(2, lmax + 1)
    expected = np.sqrt(np.sum(1.0 / ls + 1.0 / (ls + 1.0)) / 2.0 * 100.0)
    check("map rms from a flat D_l = 100", cp.map_rms_from_cl(ls, cp.cl_from_dl(ls, 100.0)),
          expected, 1e-10, "uK")
    # Binning a constant C_l returns that constant, at the weighted-mean l.
    ell_all = np.arange(2, 2509)
    leff, db = cp.planck_bin(ell_all, cp.dl_from_cl(ell_all, np.full(ell_all.size, 1e-3)))
    check("binning a constant C_l", float(cp.cl_from_dl(leff, db)[10]), 1e-3, 1e-12)
    check("first bin l_eff (l(l+1) weights)", float(leff[0]), 47.71122, 1e-6)
    check("number of Planck bins", float(len(leff)), 83.0, 0)

    print("\nBeams")
    # A 5 arcmin beam at l = 1000 keeps 82.6% of the amplitude.
    check("5' Gaussian beam at l = 1000", float(cp.gaussian_beam(1000, 5.0)), 0.8261882, 1e-6)
    import healpy as hp
    check("Gaussian beam matches healpy", float(cp.gaussian_beam(2500, 5.0)),
          float(hp.gauss_beam(np.radians(5 / 60), lmax=2500)[2500]), 1e-6)
    # A half-and-half mask: f_sky 0.5, w2 0.5.
    m = cp.mask_moments(np.r_[np.ones(500), np.zeros(500)])
    check("mask moments: f_sky", m["fsky"], 0.5, 1e-12)
    check("mask moments: f_sky,eff of a binary mask", m["fsky_eff"], 0.5, 1e-12)
    # The cut-sky correction, through its whole correlation-function round trip:
    # a uniform "mask" W = 0.8 multiplies every C_l by 0.64 and must be undone
    # exactly, even for a spectrum spanning four orders of magnitude.
    ell_m = np.arange(0, 2001)
    cl_true = np.where(ell_m >= 2, 1.0 / np.maximum(ell_m, 1) ** 2.2, 0.0)
    w_cl = np.zeros(8)
    w_cl[0] = 4 * np.pi * 0.64  # the power spectrum of a constant map W = 0.8
    back = cp.deconvolve_mask([0.64 * cl_true], w_cl, lmax=2000)[0]
    check("mask correction round trip (worst l)",
          float(np.max(np.abs(back[2:] / cl_true[2:] - 1))) + 1.0, 1.0, 1e-4)

    print("\nThe acoustic scale, against Planck's own derived parameters")
    bf = PLANCK_BF
    omh2 = bf["ombh2"] + bf["omch2"] + bf["omnuh2"]
    check("z_* (Hu & Sugiyama fit)", cp.z_star_hu_sugiyama(bf["ombh2"], omh2), bf["zstar"], 2e-3)
    cosmo = cp.make_cosmology(h=bf["h"], ombh2=bf["ombh2"], omch2=bf["omch2"])
    check("Omega_k of the constructed model", float(cosmo.Ok0) + 1.0, 1.0, 1e-12)
    acoustic = cp.angular_scale(cosmo, bf["zstar"], bf["ombh2"])
    check("sound horizon r_*", acoustic["r_s_mpc"], bf["rstar"], 2e-3, "Mpc")
    check("distance to last scattering D_M(z_*)", acoustic["d_m_mpc"], bf["DMstar"], 1e-3, "Mpc")
    check("100 theta_*", 100 * acoustic["theta_star"], 100 * bf["thetastar"], 2e-3)
    check("sound horizon at the drag epoch", cp.sound_horizon(cosmo, bf["zdrag"], bf["ombh2"]),
          bf["rdrag"], 2e-3, "Mpc")
    # Curvature moves the peaks: a closed universe makes the ruler look bigger.
    closed = cp.angular_scale(cp.make_cosmology(h=bf["h"], ombh2=bf["ombh2"], omch2=bf["omch2"],
                                                omk=-0.05), bf["zstar"], bf["ombh2"])
    assert closed["theta_star"] > acoustic["theta_star"], "curvature has the wrong sign"
    print(f"  [ok  ] closed universe (Omega_k = -0.05) enlarges theta_* by "
          f"{(closed['theta_star'] / acoustic['theta_star'] - 1) * 100:.1f}% (correct sign)")

    print()
    if FAILURES:
        print(f"{len(FAILURES)} CHECK(S) FAILED: {', '.join(FAILURES)}")
        return 1
    print("All physics checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
