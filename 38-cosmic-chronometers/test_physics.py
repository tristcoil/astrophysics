"""Self-checks of the cosmology formulas against closed forms and astropy.

    python test_physics.py
"""

from __future__ import annotations

import sys
import warnings

import numpy as np

import config as cfg
import cosmology as cosmo

warnings.filterwarnings("ignore", category=RuntimeWarning)
FAILED = []


def check(name, ok, detail=""):
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILED.append(name)


def main() -> None:
    print("Expansion rate and deceleration")
    check("E(0) = 1 for a flat universe", abs(cosmo.E2(0.0, 0.3) - 1) < 1e-14)
    check("E(0) = 1 with curvature and CPL dark energy",
          abs(cosmo.E2(0.0, 0.25, 0.6, -0.8, -0.5, 1e-4) - 1) < 1e-14)
    for om in (0.2, 0.3, 0.45):
        check(f"flat LCDM q0 = 1.5 Om - 1 (Om = {om})",
              abs(cosmo.deceleration(0.0, om) - (1.5 * om - 1)) < 1e-12)
    check("matter only: q = 1/2 at every z",
          np.allclose(cosmo.deceleration(np.linspace(0, 3, 7), 1.0), 0.5))
    check("empty universe: q = 0 (curvature does not decelerate)",
          np.allclose(cosmo.deceleration(np.linspace(0, 3, 7), 0.0, 0.0), 0.0))
    check("de Sitter: q = -1", np.allclose(cosmo.deceleration(np.linspace(0, 3, 7), 0.0, 1.0), -1))
    for om in (0.25, 0.32):
        zt_closed = float(cosmo.z_transition_flat(om))
        zt_num = cosmo.z_transition(om)
        check(f"transition: closed form = root of q (Om = {om})", abs(zt_closed - zt_num) < 1e-8,
              f"{zt_closed:.6f} vs {zt_num:.6f}")
        zz = np.linspace(0, 3, 300001)
        adot = cosmo.hubble(zz, 70, om) / (1 + zz)
        check(f"minimum of adot = transition redshift (Om = {om})",
              abs(zz[np.argmin(adot)] - zt_closed) < 2e-5)
    check("no transition when Om > 2/3", np.isnan(cosmo.z_transition_flat(0.7)))
    w = -1.3
    a = np.array([0.5, 1.0, 3.0])
    check("constant w: dark energy scales as a^(-3(1+w))",
          np.allclose(cosmo.de_density(a, w, 0.0), a ** (-3 * (1 + w))))
    check("CPL density from integrating w(a) numerically",
          abs(cosmo.de_density(0.4, -0.8, -0.6)
              - np.exp(3 * np.trapezoid((1 + cosmo.w_of_a(np.linspace(0.4, 1, 200001), -0.8, -0.6))
                                        / np.linspace(0.4, 1, 200001), np.linspace(0.4, 1, 200001))))
          < 1e-8)
    big = cosmo.E2_of_a(1e30, 0.3, None, -1.5)
    check("far-future E^2 stays finite (no 1/a - 1 round trip)", np.isfinite(big) and big > 0)

    print("Ages")
    h0, om = 67.0, 0.32
    t_closed = (2 / 3) * cosmo.T_H / (h0 * np.sqrt(1 - om)) * np.arcsinh(np.sqrt((1 - om) / om))
    t_num = cosmo.age(h0, om, Or=0.0)
    check("flat LCDM age: integral = closed form", abs(t_num / t_closed - 1) < 1e-6,
          f"{t_num:.6f} vs {t_closed:.6f} Gyr")
    check("matter-only age = 2/(3 H0)",
          abs(cosmo.age(h0, 1.0, Or=0.0) / (2 / 3 * cosmo.T_H / h0) - 1) < 1e-6)
    try:
        from astropy.cosmology import FlatLambdaCDM

        ac = FlatLambdaCDM(H0=h0, Om0=om, Tcmb0=2.7255, Neff=3.046, m_nu=0.0)
        check("age with radiation vs astropy", abs(cosmo.age(h0, om) / ac.age(0).value - 1) < 2e-4,
              f"{cosmo.age(h0, om):.4f} vs {ac.age(0).value:.4f} Gyr")
        zz = np.array([0.3, 1.0, 2.5])
        check("lookback times vs astropy",
              np.allclose(cosmo.lookback_time(zz, h0, om), ac.lookback_time(zz).value, rtol=3e-4))
        check("H(z) vs astropy (radiation included)",
              np.allclose(cosmo.hubble(zz, h0, om, Or=float(cosmo.omega_r(h0))),
                          ac.H(zz).value, rtol=1e-4))
    except ImportError:
        print("  (astropy not installed: skipped)")

    print("The Big Rip")
    for w in (-1.2, -1.5, -2.0):
        exact = cosmo.big_rip_time(70.0, 0.0, w, Or=0.0)
        closed = (2 / 3) / abs(1 + w) * cosmo.T_H / 70.0
        check(f"pure phantom universe: t_rip = 2/(3|1+w| H0) (w = {w})", abs(exact / closed - 1) < 1e-5,
              f"{exact:.4f} vs {closed:.4f} Gyr")
    c03 = cfg.CALDWELL03
    approx = cosmo.big_rip_time_approx(c03["H0"], c03["Om"], c03["w"])
    check("Caldwell et al.'s approximation gives their 22 Gyr", abs(approx - 22.0) < 0.5,
          f"{approx:.2f} Gyr")
    check("lead-time factor for w = -3/2 is ~0.3 P (Caldwell et al.)",
          abs(cosmo.rip_lead_time(1.0, -1.5) - 0.2807) < 1e-3)
    check("no rip for w >= -1", np.isinf(cosmo.big_rip_time(70, 0.3, -1.0)))

    print("Acceleration that ends")
    check("LCDM accelerates forever", cosmo.acceleration_end(67, 0.31, -1.0, 0.0) is None)
    end = cosmo.acceleration_end(67.5, 0.311, -0.838, -0.62)
    q_at = cosmo.deceleration_of_a(end["a"], 0.311, None, -0.838, -0.62, float(cosmo.omega_r(67.5)))
    check("CPL with wa < 0: q returns to zero in the future", end is not None and abs(q_at) < 1e-8,
          f"a = {end['a']:.3f}, in {end['t_from_now_gyr']:.2f} Gyr")

    print()
    if FAILED:
        print(f"{len(FAILED)} check(s) FAILED: {FAILED}")
        sys.exit(1)
    print("all physics checks passed")


if __name__ == "__main__":
    main()
