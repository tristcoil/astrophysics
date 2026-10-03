"""Checks of the analysis algorithms on inputs built so the answer is known.

These are *algorithm* tests: the inputs are synthetic (random numbers with a
gap cut into them at a chosen place), never stand-ins for data. The question
each answers is simply "if the answer were X, would the code say X?"

Imported and run by ``test_physics.py``; also runnable on its own.
"""

from __future__ import annotations

import sys

import numpy as np

import analysis
import orbit_physics as op

FAILURES: list[str] = []


def _check(label: str, got: float, expected: float, tol: float, unit: str = "") -> None:
    ok = abs(got - expected) <= tol
    if not ok:
        FAILURES.append(label)
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<52} {got:>10.5f} vs {expected:>9.5f} "
          f"{unit:<4} (tol {tol:g})")


def _synthetic_belt(rng, n=600_000, lo=2.0, hi=3.0, gaps=(), slope=0.0):
    """Uniform (or linearly sloped) semi-major axes with gaps cut in.
    Each gap is (centre, half-width, fraction removed)."""
    a = rng.uniform(lo, hi, n)
    if slope:
        # Thin the sample so density rises linearly with a.
        keep = rng.uniform(0, 1, n) < (1 + slope * (a - lo)) / (1 + slope * (hi - lo))
        a = a[keep]
    for c, hw, frac in gaps:
        inside = np.abs(a - c) < hw
        drop = inside & (rng.uniform(0, 1, a.size) < frac)
        a = a[~drop]
    return np.sort(a)


def run() -> list[str]:
    rng = np.random.default_rng(42)
    print("\nAlgorithms, on planted gaps (synthetic inputs)")

    belt = _synthetic_belt(rng, gaps=[(2.40, 0.015, 1.0), (2.70, 0.004, 0.4)])
    found = analysis.detect_gaps(belt, 2.05, 2.95)
    centers = sorted(g["a_detect"] for g in found)
    _check("blind search finds exactly the 2 planted gaps", len(found), 2, 0)
    null = analysis.detect_gaps(np.sort(rng.uniform(2.0, 3.0, 600_000)), 2.05, 2.95)
    _check("...and nothing at all in pure noise", len(null), 0, 0)
    deep = analysis.refine_gap(belt, centers[0], 0.03)
    weak = analysis.refine_gap(belt, centers[-1], 0.008)
    _check("deep gap: centre", deep["center"], 2.40, 0.001, "AU")
    _check("deep gap: depth", deep["depth"], 1.0, 0.05)
    _check("deep gap: full width at half depth", deep["fwhm"], 0.030, 0.003, "AU")
    _check("weak gap: centre", weak["center"], 2.70, 0.001, "AU")
    _check("weak gap (40% removed): depth", weak["depth"], 0.40, 0.08)

    # A gap on a steep slope -- density doubling across 1 AU -- must not be
    # dragged sideways by the slope.
    sloped = _synthetic_belt(rng, gaps=[(2.50, 0.012, 1.0)], slope=1.0)
    g = analysis.refine_gap(sloped, 2.50, 0.024)
    _check("gap on a sloped background: centre", g["center"], 2.50, 0.001, "AU")

    # A lopsided gap: the left wall ramps down gently over 0.05 AU (like the
    # inner side of the real 3:1), the right wall is a cliff. The floor is
    # centred on 2.50; the reported centre should stay within the floor.
    lop = _synthetic_belt(rng, gaps=[(2.50, 0.01, 1.0)])
    ramp = (lop > 2.44) & (lop < 2.49)
    keep = ~ramp | (rng.uniform(0, 1, lop.size) < (2.49 - lop) / 0.05)
    g = analysis.refine_gap(np.sort(lop[keep]), 2.50, 0.02)
    _check("lopsided gap: centre stays inside the empty floor", g["center"], 2.50, 0.01, "AU")

    print("\nNaming gaps")
    ident = analysis.identify(op.resonance_a(3, 1) + 0.0004)
    _check("a gap at 2.5009 AU is named 3:1 (p)", ident["p"], 3, 0)
    _check("a gap at 2.5009 AU is named 3:1 (q)", ident["q"], 1, 0)
    ident = analysis.identify(op.resonance_a(7, 3))
    _check("a gap at the 7:3 location is named 7:3 (p)", ident["p"], 7, 0)

    print("\nChance-alignment test")
    # No point in the test range, 2.0-3.5 AU, is more than 0.223 AU from a
    # resonance with q <= 3 (the worst spot is 3.5 AU, between the 2:1 at
    # 3.277 and the 5:3 at 3.700), so "within 0.25 AU" is certain.
    far = analysis.chance_of_alignment([0.25, 0.25, 0.25, 0.25], trials=200_000,
                                       chunk=100_000)
    _check("gaps 0.25 AU from resonances: chance = 1", far["probability"], 1.0, 1e-9)
    near = analysis.chance_of_alignment([0.001, 0.001, 0.001, 0.001], trials=2_000_000,
                                        chunk=1_000_000)
    # Analytically: each random point lands within 1 mAU of one of the 9
    # resonances in range with probability ~9 * 0.002 / 1.5, so all four do
    # with probability ~(0.012)^4 = 2e-8 -- far too rare to see in 2e6 trials.
    _check("gaps 1 mAU from resonances: chance ~ 0", near["probability"], 0.0, 1e-6)

    try:
        import rebound  # noqa: F401
    except ImportError:
        print("  (rebound not installed: skipping the orbit-conversion check)")
    else:
        import simulation
        print("\nSimulation: positions and velocities -> orbital elements")
        sim = rebound.Simulation()
        sim.units = ("yr", "AU", "Msun")
        sim.add(m=1.0)
        elems = [(2.2, 0.05, 3.0), (2.8, 0.25, 12.0), (3.1, 0.12, 25.0)]
        for a, e, i in elems:
            sim.add(primary=sim.particles[0], a=a, e=e, inc=np.radians(i), Omega=0.7,
                    omega=1.9, M=4.0)
        xyz = np.array([p.xyz for p in sim.particles[1:]]) - np.array(sim.particles[0].xyz)
        v = np.array([p.vxyz for p in sim.particles[1:]]) - np.array(sim.particles[0].vxyz)
        a, e, inc = simulation.xyz_to_aei(xyz, v)
        worst = max(np.abs(a - [x[0] for x in elems]).max(),
                    np.abs(e - [x[1] for x in elems]).max(),
                    np.abs(inc - [x[2] for x in elems]).max() / 100)
        _check("round trip through REBOUND (worst error)", float(worst), 0.0, 1e-9)
    return FAILURES


if __name__ == "__main__":
    fails = run()
    print(f"\n{len(fails)} FAILED: {fails}" if fails else "\nall algorithm checks passed")
    sys.exit(1 if fails else 0)
