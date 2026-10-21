"""Neutrino masses from oscillations, and what each kind of experiment measures.

Oscillation experiments measure differences of squared masses, never the
masses themselves:

* Δm²₂₁ ≈ 7.5 × 10⁻⁵ eV² ("solar"), with m₂ > m₁ known from matter effects
  in the Sun;
* |Δm²₃ℓ| ≈ 2.5 × 10⁻³ eV² ("atmospheric"), whose *sign* is unknown: m₃ is
  either the heaviest (normal ordering, NO) or the lightest (inverted, IO).

Given the lightest mass m_l, everything follows:

* NO: m₁ = m_l, m₂ = √(m_l² + Δm²₂₁), m₃ = √(m_l² + Δm²₃₁)
* IO: m₃ = m_l, m₂ = √(m_l² + |Δm²₃₂|), m₁ = √(m₂² − Δm²₂₁)

(NuFIT quotes Δm²₃ℓ = Δm²₃₁ > 0 for NO and Δm²₃₂ < 0 for IO.) Three
experiments see three different combinations:

* cosmology: the sum Σm_ν = m₁ + m₂ + m₃ (gravity counts every mass);
* β decay (KATRIN): the effective electron-neutrino mass
  m_β = √(Σ |U_ei|² m_i²), with |U_e1|² = c₁₂²c₁₃², |U_e2|² = s₁₂²c₁₃², |U_e3|² = s₁₃²;
* neutrinoless double-β decay (if neutrinos are their own antiparticles):
  m_ββ = |Σ U_ei² m_i|, which depends on two unknown Majorana phases.
"""

from __future__ import annotations

import numpy as np


def masses(m_light, dm21: float, dm3l: float, ordering: str):
    """(m1, m2, m3) in eV for lightest mass ``m_light`` (array OK).

    ``dm3l`` is NuFIT's Δm²₃ℓ: Δm²₃₁ (> 0) for NO, Δm²₃₂ (< 0) for IO.
    """
    m = np.asarray(m_light, float)
    if ordering == "NO":
        if dm3l <= 0:
            raise ValueError("normal ordering needs dm3l = dm31 > 0")
        return m, np.sqrt(m**2 + dm21), np.sqrt(m**2 + dm3l)
    if ordering == "IO":
        if dm3l >= 0:
            raise ValueError("inverted ordering needs dm3l = dm32 < 0")
        m2 = np.sqrt(m**2 - dm3l)
        return np.sqrt(m2**2 - dm21), m2, m
    raise ValueError(ordering)


def sum_masses(m_light, dm21, dm3l, ordering):
    m1, m2, m3 = masses(m_light, dm21, dm3l, ordering)
    return m1 + m2 + m3


def minimum_sum(dm21, dm3l, ordering) -> float:
    return float(sum_masses(0.0, dm21, dm3l, ordering))


def lightest_from_sum(total, dm21, dm3l, ordering) -> float:
    """Invert Σm_ν(m_l) (monotonic) by bisection; NaN below the minimum."""
    lo, hi = 0.0, float(total)
    if total < minimum_sum(dm21, dm3l, ordering):
        return float("nan")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if sum_masses(mid, dm21, dm3l, ordering) < total:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def ue_squared(s12sq: float, s13sq: float):
    """|U_e1|², |U_e2|², |U_e3|² from sin²θ₁₂ and sin²θ₁₃."""
    c13sq = 1.0 - s13sq
    return (1.0 - s12sq) * c13sq, s12sq * c13sq, s13sq


def m_beta(m_light, dm21, dm3l, ordering, s12sq, s13sq):
    """Effective electron-neutrino mass probed by β-decay endpoints (KATRIN)."""
    m1, m2, m3 = masses(m_light, dm21, dm3l, ordering)
    u1, u2, u3 = ue_squared(s12sq, s13sq)
    return np.sqrt(u1 * m1**2 + u2 * m2**2 + u3 * m3**2)


def m_betabeta_range(m_light, dm21, dm3l, ordering, s12sq, s13sq):
    """Smallest and largest |Σ U_ei² m_i| over the two Majorana phases.

    With a_i = |U_ei|² m_i the extremes are: max = a1 + a2 + a3; min = the
    largest a minus the other two, or 0 if no single term dominates (three
    "sides" that can close a triangle can cancel completely).
    """
    m1, m2, m3 = masses(m_light, dm21, dm3l, ordering)
    u1, u2, u3 = ue_squared(s12sq, s13sq)
    a = np.stack([u1 * m1, u2 * m2, u3 * m3])
    total = a.sum(axis=0)
    big = a.max(axis=0)
    lo = np.maximum(2 * big - total, 0.0)
    return lo, total


def mass_band(m_grid, ranges: dict, ordering: str, quantity: str = "sum", n_draw: int = 4000,
              seed: int = 12021):
    """Spread of a mass combination from the oscillation parameters' 3σ ranges.

    ``ranges`` maps each parameter to (best, lo3, hi3). We take the extremes
    over the corners of the 3σ box (the quantities are monotonic in each
    parameter), which is what the usual "3σ band" plots show.
    """
    keys = ("dm21", "dm3l", "s12sq", "s13sq")
    corners = np.array(np.meshgrid(*[[ranges[k][1], ranges[k][2]] for k in keys],
                                   indexing="ij")).reshape(len(keys), -1).T
    vals = []
    for c in corners:
        d = dict(zip(keys, c))
        if quantity == "sum":
            vals.append(sum_masses(m_grid, d["dm21"], d["dm3l"], ordering))
        elif quantity == "mbeta":
            vals.append(m_beta(m_grid, d["dm21"], d["dm3l"], ordering, d["s12sq"], d["s13sq"]))
        else:
            raise ValueError(quantity)
    vals = np.array(vals)
    return vals.min(axis=0), vals.max(axis=0)
