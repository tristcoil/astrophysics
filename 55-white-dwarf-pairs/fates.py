"""What happens when a white-dwarf pair touches: contact, stability and outcome.

Gravitational waves shrink a detached pair until the less massive white dwarf (the
larger one: white dwarfs grow as they lose mass) fills its Roche lobe. Mass then
flows to the companion, and the outcome depends on the mass ratio q = M2/M1:

* stable transfer -> an AM CVn binary: the orbit *widens* again as the donor is
  whittled down (HM Cnc, AM CVn, ES Cet ... in Kupfer et al.'s list are these);
* unstable transfer -> a merger within a few orbits: one more massive white dwarf,
  an R CrB or extreme-helium star, or (for massive enough pairs) a thermonuclear
  explosion, one route to Type Ia supernovae (L5).

The boundaries follow Marsh, Nelemans & Steeghs (2004, their eqs. 10, 26, 27):

* guaranteed unstable:  q > 1 + (zeta_2 - zeta_rL)/2 - k r_1^2 (1 + q) lambda
  (the classic "q > 2/3", slightly modified by the accretor's moment of inertia);
* guaranteed stable:    q < 1 + (zeta_2 - zeta_rL)/2 - sqrt((1 + q) r_h)
  (r_h -> r_1 once the stream misses the accretor and a disc forms);
* in between, stability depends on how strongly tides return the accretor's spin
  to the orbit -- a "maybe" zone that ZTF J1539+5027 sits in.

Radii use Eggleton's zero-temperature mass-radius relation (as quoted by Verbunt &
Rappaport 1988), as Marsh et al. did; the Roche lobe is Eggleton's (1983) fit; the
stream's specific angular momentum radius r_h is Verbunt & Rappaport's eq. 13 and the
minimum stream radius r_min is Nelemans et al. (2001) eq. 6, both for q = M2/M1.
"""

from __future__ import annotations

import numpy as np

import config as C
import gw

M_CH = 1.44
M_P = 0.00057


def radius_eggleton(m):
    """Zero-temperature white-dwarf radius (Rsun) for mass m (Msun)."""
    m = np.asarray(m, float)
    x = m / M_CH
    with np.errstate(invalid="ignore"):
        return (0.0114 * np.sqrt(np.clip(x ** (-2 / 3) - x ** (2 / 3), 0, None))
                * (1 + 3.5 * (m / M_P) ** (-2 / 3) + (M_P / m)) ** (-2 / 3))


def roche_lobe_fraction(q):
    """Eggleton (1983): R_L / a for the star of mass M2 = q M1."""
    q = np.asarray(q, float)
    q23 = q ** (2 / 3)
    return 0.49 * q23 / (0.6 * q23 + np.log1p(q ** (1 / 3)))


def zeta_rl(q):
    """d log(R_L/a) / d log M2 for conservative transfer (Marsh et al. 2004, eq. 9)."""
    q = np.asarray(q, float)
    q13 = q ** (1 / 3)
    return (1 + q) / 3 * (2 * np.log1p(q13) - q13 / (1 + q13)) / (0.6 * q13**2 + np.log1p(q13))


def dlog(func, m, h=1e-4):
    return (np.log(func(m * (1 + h))) - np.log(func(m * (1 - h)))) / (np.log1p(h) - np.log1p(-h))


def k_inertia(m1):
    """Moment-of-inertia factor (Marsh et al. 2004 fit to Chandrasekhar models)."""
    return 0.1939 * (1.44885 - np.asarray(m1, float)) ** 0.1917


def r_h(q):
    """Verbunt & Rappaport (1988) eq. 13 in Nelemans et al.'s (2001) form, q = M2/M1."""
    lq = np.log10(q)
    return 0.0883 - 0.04858 * lq + 0.11489 * lq**2 + 0.020475 * lq**3


def r_min(q):
    """Minimum distance of the stream from the accretor (Nelemans et al. 2001 eq. 6), / a."""
    lq = np.log10(q)
    return 0.04948 - 0.03815 * lq + 0.04752 * lq**2 - 0.006973 * lq**3


def contact_separation(m1, m2, r2_rsun=None):
    """Separation (m) at which the donor (m2) fills its Roche lobe."""
    r2 = radius_eggleton(m2) if r2_rsun is None else r2_rsun
    return r2 * C.R_SUN / roche_lobe_fraction(np.asarray(m2) / np.asarray(m1))


def contact_period(m1, m2, r2_rsun=None):
    a = contact_separation(m1, m2, r2_rsun)
    return 2 * np.pi * np.sqrt(a**3 / (C.GM_SUN * (np.asarray(m1) + np.asarray(m2))))


def time_to_contact(P, m1, m2, r2_rsun=None):
    """GW-driven time (s) from period P until contact; 0 if already in contact."""
    a0 = gw.semi_major_axis(P, np.asarray(m1) + np.asarray(m2))
    ac = contact_separation(m1, m2, r2_rsun)
    return np.clip(gw.time_between_separations(a0, ac, m1, m2), 0, None)


def stability(m1, m2):
    """Classify mass transfer at contact (accretor m1 >= donor m2).

    Returns q, the two Marsh et al. limits, the accretion mode and a label in
    {"stable", "depends on tides", "unstable"}."""
    m1, m2 = float(max(m1, m2)), float(min(m1, m2))
    q = m2 / m1
    a_c = contact_separation(m1, m2)
    r1 = radius_eggleton(m1) * C.R_SUN / a_c
    zeta2 = float(dlog(radius_eggleton, m2))
    zrl = float(zeta_rl(q))
    lam = 1 + 2 * float(dlog(radius_eggleton, m1)) + float(dlog(k_inertia, m1))
    disc = r_min(q) > r1  # stream misses the accretor
    r_ang = r1 if disc else r_h(q)
    q_unstable = 1 + (zeta2 - zrl) / 2 - k_inertia(m1) * r1**2 * (1 + q) * lam
    q_stable = 1 + (zeta2 - zrl) / 2 - np.sqrt((1 + q) * r_ang)
    if q > q_unstable:
        label = "unstable"
    elif q < q_stable:
        label = "stable"
    else:
        label = "depends on tides"
    return {"q": q, "q_stable_below": float(q_stable), "q_unstable_above": float(q_unstable),
            "mode": "disc" if disc else "direct impact", "label": label}


def boundary_curves(m1_grid):
    """Donor masses on Marsh et al.'s two lines, for each accretor mass (figure)."""
    lo, hi = [], []
    for m1 in m1_grid:
        m2s = np.linspace(0.005, min(m1, 1.4) * 0.999, 1500)
        lab = [stability(m1, m2) for m2 in m2s]
        q = m2s / m1
        qs = np.array([d["q_stable_below"] for d in lab])
        qu = np.array([d["q_unstable_above"] for d in lab])
        i_lo = np.where(q < qs)[0]
        i_hi = np.where(q > qu)[0]
        lo.append(m2s[i_lo[-1]] if i_lo.size else np.nan)
        hi.append(m2s[i_hi[0]] if i_hi.size else np.nan)
    return np.array(lo), np.array(hi)
