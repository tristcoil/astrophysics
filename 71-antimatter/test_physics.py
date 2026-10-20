"""Self-checks: the physics pieces against exact results and published numbers.

    python test_physics.py
"""

from __future__ import annotations

import json
import math
import sys

import numpy as np

import baryons
import config as cfg
import literature as lit

CHECKS = []


def check(fn):
    CHECKS.append(fn)
    return fn


def close(a, b, rel, what):
    if not abs(a - b) <= rel * abs(b):
        raise AssertionError(f"{what}: {a} vs {b} (tolerance {rel:.1%})")


# -- baryons and photons ------------------------------------------------------

@check
def photon_density():
    close(baryons.n_gamma() * 1e-6, 410.7, 0.001, "CMB photons per cm^3")


@check
def eta_conversion():
    # Pitrou et al. 2021 eq. 1: Omega_b h^2 = 0.0224 <-> eta = 6.13197e-10 at Y_p = 0.2471
    eta = float(baryons.eta_from_ombh2(0.0224, yp=0.2471, t_cmb=2.7255))
    close(eta, 6.13197e-10, 2e-4, "eta for Omega_b h^2 = 0.0224 (Pitrou et al. 2021)")
    back = float(baryons.ombh2_from_eta(eta, yp=0.2471, t_cmb=2.7255))
    close(back, 0.0224, 1e-12, "eta <-> Omega_b h^2 round trip")


@check
def entropy_per_photon_today():
    # Kolb & Turner's s = 7.04 n_gamma uses g*s = 3.909 (instantaneous neutrino
    # decoupling); Saikawa & Shirai's 3.931 gives 7.08
    close(baryons.entropy_per_photon(3.909), 7.04, 0.001, "s / n_gamma today")


# -- deuterium ------------------------------------------------------------------

@check
def dh_averages():
    close(lit.dh_weighted_mean()["mean"], 2.508, 0.001, "PDG 2025/26 D/H mean")
    close(lit.dh_weighted_mean()["err_scaled"], 0.029, 0.03, "PDG 2025/26 D/H error")
    close(lit.dh_log_mean([r for r in lit.DH_SYSTEMS if r[5]])["mean"], 2.527, 0.001,
          "Cooke et al. 2018 D/H mean")


@check
def bbn_against_literature():
    r = json.loads(cfg.RESULTS_PATH.read_text())["bbn"]["rates"]
    # PRIMAT rates at Planck's eta: D/H = 2.442-2.444e-5 (Launders et al. 2026; LINX)
    close(r["primat"]["dh_at_planck_eta"], 2.443, 0.01, "PRyMordial PRIMAT D/H at Planck eta")
    # Pitrou et al. 2021: Omega_b h^2 = 0.02195 +- 0.00022 from Cooke's D/H
    close(r["primat"]["cooke2018_7"]["ombh2_median"], 0.02195, 0.005, "BBN Omega_b h^2 (PRIMAT, Cooke)")


# -- domains --------------------------------------------------------------------

@check
def vacuum_sink_flux():
    import domains as dm

    rho, p = 1.0, 0.6
    c = math.sqrt(5 / 3 * p / rho)
    mass, _, _ = dm._sink_flux(rho, 0.0, p)
    close(-mass, 0.75**4 * rho * c, 1e-12, "rarefaction into vacuum: (3/4)^4 rho c")


@check
def klein_nishina_values():
    import domains as dm

    close(float(dm.klein_nishina(1e-4)) / cfg.SIGMA_T, 1.0, 0.001, "Thomson limit")
    close(float(dm.klein_nishina(1.0)) * 1e28, 0.2112, 0.005, "sigma_KN(1 MeV) in barn")


@check
def einstein_de_sitter_age():
    import domains as dm

    u = dm.cdg_fiducial()
    exact = 2 / (3 * cfg.h0_si(75.0) * 20.0**1.5)
    close(float(u.age(20.0)), exact, 1e-5, "EdS age at y = 20")


@check
def electron_range_matches_cdg():
    import domains as dm

    eh = dm.ElectronHeating(dm.cdg_fiducial(), 1100.0)
    close(eh.range_m / cfg.MPC, 0.8e-6, 0.1, "electron range at y = 1100 (CDG eq. A4)")


@check
def phase_space_conserves_energy():
    import domains as dm

    rng = np.random.default_rng(1)
    e, w = dm._phase_space_energies([139.57, 139.57, 134.98, 139.57, 139.57], 1876.5, 2000, rng)
    close(float(np.max(np.abs(e.sum(axis=1) - 1876.5))) + 1.0, 1.0, 1e-9, "5-pion energy sum")


@check
def cdg_reproduced():
    r = json.loads(cfg.RESULTS_PATH.read_text())["domains"]["cdg_check"]
    y, j = np.array(r["y"]), np.array(r["J_ratio"])
    assert np.all(np.abs(np.log(j[y >= 30])) < math.log(1.35)), f"J(y) vs CDG: {j}"
    s = np.array(r["spec_ratio_their_J"])[:7]
    assert np.all(np.abs(s - 1) < 0.1), f"spectrum vs CDG below 10 MeV: {s}"


@check
def photon_spectrum_matches_backenstoss():
    r = json.loads(cfg.RESULTS_PATH.read_text())["domains"]["photons"]
    assert all(abs(x - 1) < 0.05 for x in r["backenstoss_ratio_30_500MeV"]), r


# -- cosmic rays ----------------------------------------------------------------

@check
def pbar_cross_section_korsmeier_fig4():
    import pbar_production as pp

    v = pp.dsigma_dT(20.0, 2.0, galaxy=False) * 1e-31
    close(v, 1.3e-33, 0.1, "prompt dsigma/dT, T_p = 20 GeV, T_pbar = 2 GeV (Korsmeier+18 Fig. 4)")
    close(pp.multiplicity_pp(158.0), 0.04, 0.05, "pbar multiplicity in pp at 158 GeV/c")


@check
def boron_cross_sections_rw18():
    import cosmic_rays as cr

    c = float(cr.sigma_rw18(("C", "11B"), 10.0) + cr.sigma_rw18(("C", "10B"), 10.0))
    close(c, 68.6, 0.01, "C -> B at 10 GeV/n (Reinert & Winkler 2018, Table 3)")


# -- Sakharov -------------------------------------------------------------------

@check
def jarlskog_invariant():
    import sakharov as sk

    c = lit.CKM["pdg2026"]
    J, _ = sk.jarlskog_from_wolfenstein(c["lambda"][0], c["A"][0], c["rhobar"][0], c["etabar"][0])
    close(J, 3.16e-5, 0.01, "J from the Wolfenstein parameters (PDG 2026)")


@check
def sphaleron_freezeout():
    r = json.loads(cfg.RESULTS_PATH.read_text())["sakharov"]
    assert abs(r["T_star_ours"] - 131.7) < 2.3, r["T_star_ours"]


def main() -> int:
    fails = 0
    for fn in CHECKS:
        try:
            fn()
            print(f"  ok    {fn.__name__}")
        except Exception as e:                     # noqa: BLE001
            fails += 1
            print(f"  FAIL  {fn.__name__}: {e}")
    print(f"{len(CHECKS) - fails}/{len(CHECKS)} checks passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
