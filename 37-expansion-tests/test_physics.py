"""Self-checks of the physics and the statistics against published numbers and brute force.

    python test_physics.py
"""

from __future__ import annotations

import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

import reference_method as rm
import tcmb
import template_fit as tf
import tolman
from config import TOLMAN_CLUSTERS

FAILURES = []


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        FAILURES.append(name)


def test_distances():
    print("Distances for the Tolman test (Lubin & Sandage 2001, Paper IV Tables 1 and 8)")
    published = {  # cluster: (q0=1/2, q0=0, q0=1, tired light) log A
        "Cl 1324+3011": (3.910, 3.992, 3.853, 4.215),
        "Cl 1604+4304": (3.924, 4.021, 3.860, 4.270),
        "Cl 1604+4321": (3.926, 4.025, 3.861, 4.280),
    }
    for cl, vals in published.items():
        z = TOLMAN_CLUSTERS[cl]["z"]
        for model, v in zip(("q0=0.5", "q0=0", "q0=1", "tired"), vals):
            ours = float(tolman.log_a(z, model))
            check(f"log A {cl} {model}", abs(ours - v) < 0.0015, f"{ours:.4f} vs {v}")
    # Paper I's distance factor f for the nearby cluster A76 (z = 0.0378): 3.013.
    f = float(tolman.log_a(0.0378, "q0=0.5"))
    check("Paper I f for A76", abs(f - 3.013) < 0.001, f"{f:.4f} vs 3.013")
    # q0 = 1/2 Mattig equals flat Omega_m = 1 in astropy.
    from astropy.cosmology import FlatLambdaCDM
    eds = FlatLambdaCDM(H0=50, Om0=1.0, Tcmb0=0)
    for z in (0.1, 0.9, 2.0):
        a = tolman.angular_diameter_distance_pc(z, "q0=0.5")
        b = eds.angular_diameter_distance(z).value * 1e6
        check(f"Mattig q0=1/2 vs astropy EdS at z={z}", abs(a / b - 1) < 1e-6, f"ratio {a / b:.8f}")
    # Tired light's linear radii are larger by the published 0.305 dex for Cl 1324+3011.
    d = float(tolman.log_a(0.7565, "tired") - tolman.log_a(0.7565, "q0=0.5"))
    check("tired-light radius offset, Cl 1324+3011", abs(d - 0.305) < 0.002, f"{d:.3f} vs 0.305")


def test_sz():
    print("Sunyaev-Zel'dovich spectrum")
    from scipy.optimize import brentq
    x0 = brentq(tcmb.sz_spectrum, 3.0, 4.5)
    check("null of x^4 e^x/(e^x-1)^2 (x coth(x/2) - 4)", abs(x0 - 3.8300) < 5e-4, f"x0 = {x0:.4f}")
    nu = tcmb.null_frequency_ghz(2.72548, 0.0)
    check("null frequency today", abs(nu - 217.5) < 0.1, f"{nu:.2f} GHz")
    for z in (0.3, 1.0):
        same = tcmb.null_frequency_ghz(2.72548 * (1 + z), z)
        tired = tcmb.null_frequency_ghz(2.72548, z)
        check(f"expanding: null independent of z ({z})", abs(same - nu) < 1e-9)
        check(f"tired light: null at 217.5/(1+z) ({z})", abs(tired - nu / (1 + z)) < 1e-9)


def test_co_equilibrium():
    print("CO statistical equilibrium (LAMDA)")
    mol = tcmb.read_lamda()
    check("CO J=1-0 frequency", abs(mol["nu_ghz"][0] - 115.2712) < 1e-3, f"{mol['nu_ghz'][0]} GHz")
    for t_cmb in (2.725, 7.4, 10.1):
        tex = tcmb.co_excitation_temperature(mol, t_cmb, 1e-6, 100.0)
        check(f"no collisions: T_ex = T_CMB ({t_cmb} K)", abs(tex - t_cmb) < 1e-3, f"{tex:.4f}")
    tex = tcmb.co_excitation_temperature(mol, 2.725, 1e8, 60.0, n_levels=25)
    check("very dense gas: T_ex -> T_kin (60 K)", abs(tex / 60.0 - 1) < 0.02, f"{tex:.2f} K")
    lo = tcmb.co_excitation_temperature(mol, 7.4, 30.0, 100.0)
    hi = tcmb.co_excitation_temperature(mol, 7.4, 300.0, 100.0)
    check("excess grows with density", 7.4 < lo < hi, f"{lo:.2f} < {hi:.2f} K")


def test_beta_fit():
    print("T(z) fits")
    z = np.array([0.0, 0.2, 0.5, 0.9, 1.7, 2.4, 3.0])
    for beta in (-0.03, 0.0, 0.05):
        df = pd.DataFrame({"z": z, "T": 2.72548 * (1 + z) ** (1 - beta),
                           "err_plus": 0.05 * (1 + z), "err_minus": 0.05 * (1 + z),
                           "err_sys": 0.0, "kind": "value"})
        df.loc[0, ["err_plus", "err_minus"]] = 0.00057
        r = tcmb.fit_beta(df)
        check(f"recover beta = {beta} from noiseless data", abs(r["beta"] - beta) < 2e-5,
              f"{r['beta']:.5f}")
    # Split-normal: a measurement above the model uses its lower error.
    d = pd.DataFrame({"z": [1.0], "T": [6.0], "err_plus": [1.0], "err_minus": [0.5],
                      "err_sys": [0.0], "kind": ["value"]})
    c = tcmb.chi2(d, 0.0)
    check("split-normal error choice", abs(c - ((5.45096 - 6.0) / 0.5) ** 2) < 1e-6, f"{c:.4f}")


def test_reference_method():
    print("Reference-curve method")
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(-20, 60, 3000))
    f = np.exp(-0.5 * (t / 12) ** 2) + rng.normal(0, 0.05, t.size)
    prof = rm.MedianProfile(t, f, tau=2.0, step=0.01)
    x = rng.uniform(-18, 58, 200)
    brute = np.array([np.median(f[(t >= xi - 2) & (t <= xi + 2)]) for xi in x])
    err = np.max(np.abs(prof.value(x) - brute))
    check("precomputed running median = brute force", err < 0.02, f"max diff {err:.4f}")
    flux = np.arange(30.0)
    d = rm.binned_dispersion(np.arange(30.0), flux, 3)
    check("binned dispersion (3 chunks of 10)", abs(d - np.std(np.arange(10.0))) < 1e-12, f"{d:.4f}")
    # Width of a stretched light curve against its own unstretched reference.
    tt = np.linspace(-15, 40, 12)
    for w_true in (1.0, 1.6):
        ff = np.exp(-0.5 * (tt / w_true / 12) ** 2)
        grid = rm.chi2_on_grid(prof, tt, ff, np.full(tt.size, 0.05), rm.W_GRID)
        w = rm.W_GRID[np.argmin(grid)]
        check(f"width recovered (true {w_true})", abs(w - w_true) < 0.06, f"{w:.3f}")


def test_template_fitter():
    print("Global template fitter on simulated light curves")
    rng = np.random.default_rng(5)
    rows = []
    bands = {"g": 473.0, "r": 642.0, "i": 784.0, "z": 926.0}
    for i in range(160):
        z = rng.uniform(0.1, 1.0)
        for band, lam in bands.items():
            lam_rest = lam / (1 + z)
            if lam_rest < 320 or lam_rest > 840:
                continue
            t = np.arange(-20, 80, 5.0) + rng.uniform(-2, 2)
            width = 14 * (1 + z) * (1 + 0.3 * (lam_rest - 500) / 500)   # true b = 1
            flux = 100 * np.exp(-0.5 * (t / width) ** 2)
            err = np.full(t.size, 2.0)
            rows.append(pd.DataFrame({"cid": i, "band": band, "z": z, "t": t,
                                      "flux": flux + rng.normal(0, 2.0, t.size), "err": err}))
    pts = pd.concat(rows, ignore_index=True)
    F = tf.TemplateFitter(pts)
    r = F.fit_power(0.5)
    check("recovers b = 1 from simulated data", abs(r["b"] - 1) < 0.02,
          f"b = {r['b']:.4f} ± {r['b_err_curvature']:.4f}")
    S0 = F.template_grid([0.0], [350, 500, 800])
    check("template normalised to 1 at peak at every wavelength", np.allclose(S0, 1.0, atol=1e-9))
    c0, c1 = F.chi2_power(0.0), F.chi2_power(1.0)
    check("b = 0 fits much worse than b = 1", c0 > c1 + 100, f"chi2 {c0:.0f} vs {c1:.0f}")


def main() -> int:
    test_distances()
    test_sz()
    test_co_equilibrium()
    test_beta_fit()
    test_reference_method()
    test_template_fitter()
    print(f"\n{len(FAILURES)} failures")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
