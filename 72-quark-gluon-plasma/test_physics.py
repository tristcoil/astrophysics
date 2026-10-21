"""Self-checks: physics limits, bookkeeping, read-back figures, and every
published number the project claims to reproduce.

    python test_physics.py        (no pytest needed)
"""

from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
from scipy.special import kn, kve  # noqa: E402

import config as cfg  # noqa: E402
import datasets as D  # noqa: E402
import figparse  # noqa: E402
import hadrons  # noqa: E402
import hrg  # noqa: E402
import lattice  # noqa: E402
import thermal_fit as F  # noqa: E402
import universe as U  # noqa: E402

CHECKS = []


def check(fn):
    CHECKS.append(fn)
    return fn


def _series(m, T, mu, stat, K=80):
    k = np.arange(1, K + 1)
    s = np.where((stat == 1) & (k % 2 == 0), -1.0, 1.0)
    f = s * np.exp(k * (mu - m) / T)
    n = m * m * T / (2 * np.pi ** 2) * np.sum(f * kve(2, k * m / T) / k)
    e = m * m * T / (2 * np.pi ** 2) * np.sum(f * (3 * T * kve(2, k * m / T) / k ** 2
                                                    + m * kve(1, k * m / T) / k))
    return n, e


@check
def quadrature_matches_bessel_series():
    """64-point Laguerre momentum integrals = exact Bessel series (1e-8)."""
    gas = hrg.HRG(widths="zero")
    worst = 0.0
    for T in (0.10, 0.1565, 0.20):
        for pdg, mu in ((211, 0.0), (321, 0.05), (2212, 0.4), (1000020040, 0.0), (3334, 0.0)):
            i = gas.hl.index[pdg]
            mu_i = np.zeros(len(gas.hl))
            mu_i[i] = mu
            st = gas._state_mu(T, mu_i)
            n, e = _series(gas.hl.mass[i], T, mu, gas.hl.stat[i])
            g = gas.g[i] / cfg.HBARC ** 3
            worst = max(worst, abs(st.n[i] / (g * n) - 1), abs(st.e[i] / (g * e) - 1))
    assert worst < 1e-8, worst
    return f"worst relative difference {worst:.1e}"


@check
def boltzmann_limit():
    gas = hrg.HRG(widths="zero", stats="boltzmann")
    T = 0.150
    st = gas.state(T)
    i = gas.hl.index[211]
    m = gas.hl.mass[i]
    exact = m * m * T * kn(2, m / T) / (2 * np.pi ** 2) / cfg.HBARC ** 3
    assert abs(st.n[i] / exact - 1) < 1e-9
    return f"pi+ density {st.n[i]:.6f} fm^-3 = m^2 T K2/(2 pi^2)"


@check
def decay_chains_conserve_charges():
    for name in ("PDG2025", "PDG2014"):
        hl = hadrons.load(name)
        Fm = hadrons.feeddown_matrix(name, "strong")
        assert np.abs(Fm @ hl.B - hl.B).max() < 1e-12
        assert np.abs(Fm @ hl.Q - hl.Q).max() < 1e-12
        assert not hadrons.conserved_in_decays(name)
    return "baryon number and charge conserved along every chain (both lists)"


@check
def conservation_laws_solved():
    gas = hrg.model()
    st = gas.solve_mu(0.150, 0.300, cfg.Q_OVER_B["Au"])
    t = gas.totals(st)
    assert abs(t["nS"]) < 1e-12 and abs(t["nQ"] / t["nB"] - cfg.Q_OVER_B["Au"]) < 1e-10
    st0 = gas.solve_mu(0.1565, 0.0)
    assert abs(st0.muS) < 1e-15 and abs(st0.muQ) < 1e-15
    return f"mu_S = {st.muS * 1e3:.2f} MeV, mu_Q = {st.muQ * 1e3:.2f} MeV at mu_B = 300 MeV"


@check
def excluded_volume_is_a_common_factor():
    gas = hrg.model(stats="boltzmann")
    st = gas.state(0.1565)
    sup = gas.excluded_volume(st)["suppression"][st.n > 0]
    assert sup.max() / sup.min() - 1 < 1e-9
    return f"every density x {sup.mean():.4f} (r0 = 0.3 fm)"


@check
def pions_mostly_from_decays():
    r = json.loads((cfg.DERIVED_DIR / "main_fit.json").read_text())
    f = [y["decay_fraction"] for y in r["yields"] if y["key"] == "pi+"][0]
    assert 0.6 < f < 0.75
    return f"{100 * f:.0f}% of pi+ from decays (Andronic et al.: 70%)"


@check
def hrg_matches_lattice_below_tpc():
    gas = hrg.model()
    h = lattice.hotqcd_2014()
    worst = 0.0
    for T, p, e in zip(h["T"], h["p"], h["p_err"]):
        if T > 150:
            break
        tot = gas.totals(gas.state(T / 1e3))
        worst = max(worst, abs(tot["p"] * cfg.HBARC ** 3 / (T / 1e3) ** 4 - p) / e)
    assert worst < 1.0
    return f"p/T^4 within {worst:.2f} sigma of HotQCD 2014 for T <= 150 MeV"


@check
def hotqcd_ansatz_reproduces_table():
    h = lattice.hotqcd_2014()
    a = lattice.hotqcd_ansatz_thermo(h["T"])
    dev = max(np.abs(a["p"] / h["p"] - 1).max(), np.abs(a["e"] / h["e"] - 1).max())
    assert dev < 0.015
    return f"ansatz vs Table I: {100 * dev:.2f}% at worst"


@check
def compilation_matches_andronic_figure():
    fig = figparse.andronic_fig1()
    tf = D.thermal_fist_compilation()
    worst = 0.0
    for k, e in fig["entries"].items():
        y = D.alice_andronic_compilation().get(k)
        worst = max(worst, abs(e["data"] / y.value - 1))
    assert worst < 0.003 and fig["axis_fit_residual_decades"] < 1e-3
    return f"Fig. 1 data = Thermal-FIST compilation to {100 * worst:.2f}%"


@check
def star_figure_points_match_hepdata():
    r = json.loads((cfg.DERIVED_DIR / "phase_diagram.json").read_text())
    dT = max(abs(m["dT"]) for m in r["star_figure_vs_hepdata"])
    assert dT < 0.5 and len(r["star_figure_vs_hepdata"]) == 7
    return f"7 STAR points read from Andronic Fig. 4 within {dT:.2f} MeV of HEPData"


@check
def same_thermal_densities_as_andronic():
    r = json.loads((cfg.DERIVED_DIR / "andronic_comparison.json").read_text())
    ratio = r["variants"]["default"]["ratio"]
    pure = [ratio[k] for k in ("Omega-", "Omegabar+", "d", "dbar", "He3", "He4")]
    assert max(abs(x - 1) for x in pure) < 0.01
    return "Omega and nuclei (no feed-down) agree with their model to < 1%"


@check
def andronic_fit_reproduced():
    r = json.loads((cfg.DERIVED_DIR / "andronic_comparison.json").read_text())["fit_their_data"]
    A = cfg.ANDRONIC_2018
    zT = (r["T"] - A["T"]) / A["T_err"]
    zmu = (r["muB"] - A["muB"]) / A["muB_err"]
    zV = (r["V_EV"] - A["V"]) / A["V_err"]
    assert abs(zT) < 1 and abs(zmu) < 1 and abs(zV) < 1.5
    return (f"their data: T {r['T']:.1f} ({zT:+.2f} sigma), mu_B {r['muB']:.1f} ({zmu:+.2f}), "
            f"V {r['V_EV']:.0f} ({zV:+.2f})")


@check
def star_fits_reproduced():
    r = json.loads((cfg.DERIVED_DIR / "star.json").read_text())
    out = []
    for e in ("7.7", "39.0"):
        a, p = r["replica"][e], r["published"][e]
        assert abs(a["T"] - p["T"]) < p["T_err"] and abs(a["muB"] - p["muB"]) < p["muB_err"]
        assert max(abs(v - 1) for v in a["model_ratio_to_THERMUS"].values()) < 0.03
        out.append(f"{e} GeV: T {a['T']:.1f}/{p['T']:.1f}, mu_B {a['muB']:.0f}/{p['muB']:.0f}")
    return "; ".join(out)


@check
def photon_slope_reproduced():
    r = json.loads((cfg.DERIVED_DIR / "photons.json").read_text())
    p = cfg.ALICE_PHOTONS_2016
    assert abs(r["T_eff_MeV"] - p["T_eff_no_pqcd_subtraction"]) < 1.0
    assert abs(r["stat_MeV"] - p["stat_ns"]) < 1.0
    return f"T_eff {r['T_eff_MeV']:.1f} +- {r['stat_MeV']:.1f} (ALICE: 304 +- 11)"


@check
def universe_age_matches_p2():
    if U.p2_calendar() is None:
        return "P2 not run: skipped"
    a = U.age_at(156.5)[0]
    b = U.friedmann_age(np.array([156.5]))[0]
    assert abs(b / a - 1) < 0.02
    return f"{a * 1e6:.2f} us (P2) vs {b * 1e6:.2f} us (own integration)"


@check
def baryon_to_photon_ratio():
    eta = U.eta_from_omega_b()
    assert 6.0e-10 < eta < 6.25e-10
    return f"eta = {eta:.4e} from omega_b = {cfg.PLANCK_OMEGA_B[0]}"


@check
def charm_decays_sensible():
    pp = hadrons.charm_feeddown_per_pair()
    assert 0.3 < pp[321] < 0.7 and abs(pp[321] - pp[-321]) < 1e-12
    return f"per c-cbar pair: {pp[321]:.2f} K+, {pp[211]:.2f} pi+"


@check
def yields_compilation():
    ds = D.alice("0-10")
    tf = D.thermal_fist_compilation()
    for k, pdg in (("pi+", 211), ("K-", -321), ("Xi-", 3312), ("d", 1000010020),
                   ("H3L", 1010010030), ("He4bar", -1000020040)):
        assert abs(ds.get(k).value / tf[pdg][0] - 1) < 4e-3, k   # its file keeps 3 digits
    return "pi, K, Xi, d, hypertriton, anti-4He equal Thermal-FIST's file"


def main() -> int:
    fails = 0
    for fn in CHECKS:
        try:
            msg = fn()
            print(f"  PASS  {fn.__name__}: {msg}")
        except Exception:
            fails += 1
            print(f"  FAIL  {fn.__name__}")
            traceback.print_exc()
    print(f"{len(CHECKS) - fails}/{len(CHECKS)} checks passed")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
