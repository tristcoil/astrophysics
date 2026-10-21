"""Physics checks: textbook universes, identities, the external codes, and the published numbers.

    python test_physics.py        # prints PASS/FAIL per check, exit code 1 on any failure

Runs offline once ``fetch_data.py`` has cached the Davis & Lineweaver source.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402

import config as C  # noqa: E402
from cosmology import Universe  # noqa: E402

warnings.filterwarnings("ignore")
FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))
    if not ok:
        FAILS.append(name)


def rel(a, b):
    return float(np.max(np.abs(np.asarray(a) / np.asarray(b) - 1)))


def main() -> int:
    z = np.array([0.05, 0.3, 1.0, 2.0, 5.0, 20.0])

    # --- 1. Textbook universes with closed-form answers ------------------------------
    eds = Universe(C.EDS)
    dh = eds.hubble_distance
    check("EdS age = 2/(3 H0)", abs(eds.age / (2 / 3 * dh) - 1) < 1e-9, f"{eds.age:.6f} Gyr")
    check("EdS particle horizon = 2c/H0 = 3 c t0", abs(eds.eta0 / (2 * dh) - 1) < 1e-9)
    check("EdS chi(z) = (2c/H0)(1 - 1/sqrt(1+z))", rel(eds.chi_of_z(z), 2 * dh * (1 - 1 / np.sqrt(1 + z))) < 1e-9)
    check("EdS: D_A peaks at z = 5/4", abs(eds.z_max_angular_diameter() - 1.25) < 1e-6,
          f"{eds.z_max_angular_diameter():.8f}")
    check("EdS has no event horizon", not np.isfinite(eds.eta_inf))

    milne = Universe(C.MILNE)
    check("Milne chi = (c/H0) ln(1+z)", rel(milne.chi_of_z(z), milne.hubble_distance * np.log1p(z)) < 1e-9)
    check("Milne D_L = (c/H0) z (1 + z/2)", rel(milne.dl_of_z(z), milne.hubble_distance * z * (1 + z / 2)) < 1e-9)
    check("Milne: D_A grows forever (no turnover)", not np.isfinite(milne.z_max_angular_diameter()))
    check("Milne: v_then = c ln(1+z) crosses c at z = e - 1",
          abs(milne.z_light_cone_meets_hubble_sphere() - (np.e - 1)) < 1e-7)

    ds = Universe(C.DE_SITTER)
    check("de Sitter chi = cz/H0", rel(ds.chi_of_z(z), ds.hubble_distance * z) < 1e-9)
    check("de Sitter: event horizon = Hubble sphere = c/H0",
          abs((ds.eta_inf - ds.eta0) / ds.hubble_distance - 1) < 1e-8)
    check("de Sitter: v_then = z/(1+z) < c", rel(ds.v_then(z), z / (1 + z)) < 1e-8)

    rad = Universe(C.PLANCK18.but(Om=0.0, m_nu=()))  # photons + 3.046 massless neutrinos + Lambda
    a = np.array([1e-8, 1e-6, 1e-5])
    # radiation + Lambda + curvature are negligible this early: t = a^2/(2 H0 sqrt(Or))
    check("radiation era: t = a^2 / (2 H0 sqrt(Omega_r))",
          rel(rad.t_of_a(a), a**2 / (2 * rad.H0 * np.sqrt(rad.Or))) < 1e-6)

    # --- 2. Identities in a flat universe ----------------------------------------------
    u = Universe(C.PLANCK18)
    zt = u.z_light_cone_meets_hubble_sphere()
    check("flat: light cone meets Hubble sphere where D_A peaks",
          abs(zt - u.z_max_angular_diameter()) < 1e-6, f"{zt:.6f}")
    t_turn = float(u.t_of_z(zt))
    tt = np.linspace(t_turn - 0.3, t_turn + 0.3, 6001)
    d_ph = u.a_of_t(tt) * (u.eta0 - u.eta_of_t(tt))
    check("photon proper distance is largest when it crosses the Hubble sphere",
          abs(tt[np.argmax(d_ph)] - t_turn) < 2e-4, f"{tt[np.argmax(d_ph)]:.4f} vs {t_turn:.4f} Gyr")
    # Hubble's law is exact: d(a chi)/dt = H a chi at any distance
    t0, chi, h = 5.0, 30.0, 1e-4
    ddt = (u.a_of_t(t0 + h) - u.a_of_t(t0 - h)) * chi / (2 * h)
    check("Hubble's law v = H D holds at 30 billion ly", abs(ddt / (u.H_of_t(t0) * u.a_of_t(t0) * chi) - 1) < 1e-6)
    zz = np.array([0.51, 1.484, 2.33])
    ratio = u.dm_of_z(zz) / ((1 + zz) * (1.0 / u.H_of_z(zz)))
    check("DESI's D_M/((1+z) D_H) is the recession speed at emission", rel(ratio, u.v_then(zz)) < 1e-10)
    check("event horizon today = eta_inf - eta0 = chi(z_EH)",
          abs(u.chi_of_z(u.z_event_horizon_now()) - (u.eta_inf - u.eta0)) < 1e-8)
    tg = np.linspace(4.0, 12.0, 40001)
    z_peak = 1 / u.a_of_t(tg[np.argmax(1 / (u.a_of_t(tg) * u.H_of_t(tg)))]) - 1
    check("comoving Hubble radius 1/(aH) peaks when the acceleration starts (q = 0)",
          abs(z_peak - u.z_acceleration_start()) < 1e-3,
          f"peak z = {z_peak:.4f}, q = 0 at z = {u.z_acceleration_start():.4f}")
    tl = u.t_last_seen(np.array([u.eta_inf - u.eta0]))[0]
    check("a galaxy on today's event horizon is last seen as it is today", abs(tl - u.age) < 1e-6)

    # --- 3. Against astropy and CAMB ------------------------------------------------------
    from astropy.cosmology import Planck18 as AP

    zc = np.array([0.1, 1.0, 3.0, 10.0, 1089.8])
    check("comoving distances = astropy Planck18 (to 3e-6)",
          rel(u.chi_of_z(zc) / C.GLY_PER_MPC, AP.comoving_distance(zc).value) < 3e-6)
    check("age = astropy (to 1e-6)", abs(u.age / AP.age(0).value - 1) < 1e-6)
    try:
        import camb

        from analysis import planck18_camb_params

        res = camb.get_background(planck18_camb_params())
        check("conformal time today = CAMB tau0 (to 2e-5)",
              abs(u.eta0 / (res.conformal_time(0) * C.GLY_PER_MPC) - 1) < 2e-5)
        check("age = CAMB (to 3e-5)", abs(u.age / res.get_derived_params()["age"] - 1) < 3e-5)
    except ImportError:
        print("SKIP  CAMB not installed")

    # --- 4. The numbers this project is checked against -------------------------------------
    s = u.summary()
    for k, v in C.SPEC_CHECKS.items():
        tol = 0.01 if k == "reachable_fraction" else (0.05 if "gly" in k else 0.005)
        check(f"spec: {k} = {v}", abs(s[k] - v) <= tol, f"{s[k]:.4f}")
    check("Planck 2018 age 13.787 +- 0.020", abs(u.age - 13.787) < 0.02, f"{u.age:.4f}")
    check("Planck 2018 z_eq 3387 +- 21", abs(u.z_equality() - 3387) < 21, f"{u.z_equality():.1f}")
    dl = Universe(C.DL2004)
    check("D&L: Hubble sphere at z ~ 1.46", abs(dl.z_hubble_sphere_now() - 1.46) < 0.005)
    check("D&L: event horizon at z ~ 1.8", abs(dl.z_event_horizon_now() - 1.8) < 0.05)
    check("D&L: CMB (z = 1100) recedes at 3.2c now and 58.1c then",
          abs(dl.v_now(1100) - 3.2) < 0.05 and abs(dl.v_then(1100) - 58.1) < 0.1,
          f"{float(dl.v_now(1100)):.3f}, {float(dl.v_then(1100)):.2f}")

    # --- 5. Reading Davis & Lineweaver's figures back ----------------------------------------
    src = C.LIT_DIR / "dl2004_src"
    if (src / "fig1a.eps").exists():
        import dl2004

        cur = dl2004.figure1_curves(src / "fig1a.eps")
        check("EPS axis calibration from labelled ticks (rms < 1e-4 page units)",
              cur["calibration"].rms_x < 1e-4 and cur["calibration"].rms_y < 1e-4)
        lc = np.concatenate(cur["light_cone"])
        check("their light cone closes at their age, 13.47 Gyr", abs(lc[:, 1].max() - dl.age) < 0.012,
              f"{lc[:, 1].max():.3f} vs {dl.age:.3f}")
        f5 = dl2004.figure5_points(src / "fig5.eps")
        check("their Fig. 5 holds the 42 Perlmutter et al. supernovae", len(f5["z"]) == 42)
    else:
        print("SKIP  D&L source not cached")

    print(f"\n{len(FAILS)} failure(s)" if FAILS else "\nall checks passed")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
