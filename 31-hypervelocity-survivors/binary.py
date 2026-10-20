"""The orbital-speed argument: how fast is a white dwarf released when its
companion explodes?

Two white dwarfs spiral together by emitting gravitational waves until the
lighter one (the *donor*) is so close that the heavier one (the *accretor*)
starts pulling matter off it: the donor then exactly fills its Roche lobe,
the teardrop-shaped region inside which its own gravity wins. Roche-lobe
size depends only on the orbit's size a and the mass ratio q = M_don/M_acc
(Eggleton 1983):

    R_L / a = 0.49 q^(2/3) / (0.6 q^(2/3) + ln(1 + q^(1/3)))

so the donor's radius fixes the separation, and Kepler's third law then fixes
how fast the donor is moving around the pair's centre of mass:

    v_orb = M_acc * sqrt(G / ((M_acc + M_don) a)).

If the accretor explodes and is completely destroyed, nothing holds the donor
any more and it flies off in a straight line at v_orb (plus a small sideways
kick from the debris hitting it). White dwarfs are tiny, a few thousand km,
because degenerate matter gets *smaller* as you add mass, so these orbits are
minutes long and v_orb is 1,000-3,000 km/s: faster than any other star could
be launched.

Two mass-radius relations are used:

* Bedard et al. (2020) evolutionary models (Montreal cooling sequences, CO
  cores, thick hydrogen layers), at a chosen effective temperature -- what
  Bauer et al. (2021) and El-Badry et al. (2023) used;
* Eggleton's zero-temperature fit (quoted by Verbunt & Rappaport 1988), the
  coldest (smallest) a white dwarf of that mass can be, so the fastest launch.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from config import AUX_DIR

G_MSUN_KM3_S2 = 1.32712440018e11     # G * M_sun, km^3 / s^2
R_SUN_KM = 695700.0
#: Bedard et al. (2020) masses with published cooling sequences.
BEDARD_MASSES = np.round(np.arange(0.20, 1.3001, 0.05), 2)


# --------------------------------------------------------------------------
# White dwarf radii
# --------------------------------------------------------------------------

def eggleton_radius(m):
    """Zero-temperature white dwarf radius in R_sun (Eggleton's fit to the
    Hamada-Salpeter relation, as given by Verbunt & Rappaport 1988)."""
    m = np.asarray(m, float)
    mch, mp = 1.44, 0.00057
    x = m / mch
    return (0.0114 * np.sqrt(np.clip(x ** (-2 / 3) - x ** (2 / 3), 0, None))
            * (1 + 3.5 * (m / mp) ** (-2 / 3) + (mp / m)) ** (-2 / 3))


def bedard_path(mass: float, layer: str = "thick") -> Path:
    return AUX_DIR / "bedard2020" / f"seq_{int(round(mass * 100)):03d}_{layer}.txt"


def read_bedard(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Teff (K) and radius (cm) along one cooling sequence. Each model takes
    three lines; the first holds #, Teff, log g, R, age, L."""
    teff, rad = [], []
    lines = [ln for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("=")
             and not ln.lstrip().startswith(("#", "Log", "Lnu"))]
    for ln in lines:
        parts = ln.split()
        # model lines start with an integer model number
        if len(parts) == 6 and parts[0].isdigit():
            teff.append(float(parts[1]))
            rad.append(float(parts[3]))
    return np.array(teff), np.array(rad)


@lru_cache(maxsize=8)
def _bedard_mr(teff: float, layer: str) -> PchipInterpolator:
    masses, radii = [], []
    for m in BEDARD_MASSES:
        t, r = read_bedard(bedard_path(m, layer))
        order = np.argsort(t)
        if not (t.min() <= teff <= t.max()):
            continue
        masses.append(m)
        radii.append(np.interp(np.log(teff), np.log(t[order]), np.log(r[order])))
    return PchipInterpolator(np.array(masses), np.array(radii), extrapolate=False)


def bedard_radius(m, teff: float = 1.0e4, layer: str = "thick"):
    """Radius in R_sun of a CO white dwarf of mass m at temperature teff,
    interpolated in the Bedard et al. (2020) sequences (0.20-1.30 Msun)."""
    lnr_cm = _bedard_mr(float(teff), layer)(np.asarray(m, float))
    return np.exp(lnr_cm) / 1e5 / R_SUN_KM


def radius(m, model: str = "bedard", teff: float = 1.0e4):
    if model == "bedard":
        return bedard_radius(m, teff)
    if model == "eggleton":
        return eggleton_radius(m)
    raise ValueError(model)


# --------------------------------------------------------------------------
# Roche geometry and orbital speed
# --------------------------------------------------------------------------

def roche_lobe_fraction(q):
    """Eggleton (1983) volume-equivalent Roche-lobe radius over separation,
    for the star of mass ratio q = M_this / M_other."""
    q = np.asarray(q, float)
    q23 = q ** (2 / 3)
    return 0.49 * q23 / (0.6 * q23 + np.log1p(q ** (1 / 3)))


def contact_orbit(m_don, m_acc, r_don_rsun):
    """Separation (km), donor orbital speed (km/s) and period (s) when a donor
    of radius r_don exactly fills its Roche lobe."""
    m_don = np.asarray(m_don, float)
    m_acc = np.asarray(m_acc, float)
    a = np.asarray(r_don_rsun, float) * R_SUN_KM / roche_lobe_fraction(m_don / m_acc)
    v = m_acc * np.sqrt(G_MSUN_KM3_S2 / ((m_acc + m_don) * a))
    period = 2 * np.pi * np.sqrt(a ** 3 / (G_MSUN_KM3_S2 * (m_acc + m_don)))
    return a, v, period


def kick_speed(m_don, m_acc, r_don_rsun, eta: float = 0.5, v_ejecta: float = 1.0e4):
    """Sideways kick from supernova debris striking the donor (El-Badry et al.
    2023, eq. 7): a fraction eta of the momentum of the ejecta (mass M_acc at
    v_ejecta) that falls on the donor's cross-section pi R^2 at distance a."""
    a, _, _ = contact_orbit(m_don, m_acc, r_don_rsun)
    r = np.asarray(r_don_rsun, float) * R_SUN_KM
    return eta * (m_acc * v_ejecta / m_don) * (np.pi * r ** 2) / (4 * np.pi * a ** 2)


def ejection_speed(m_don, m_acc, model: str = "bedard", teff: float = 1.0e4,
                   with_kick: bool = True):
    """Launch speed of the surviving donor: orbital speed at contact, with the
    debris kick (perpendicular to it) added in quadrature."""
    r = radius(m_don, model, teff)
    _, v, _ = contact_orbit(m_don, m_acc, r)
    if with_kick:
        v = np.hypot(v, kick_speed(m_don, m_acc, r))
    return v


def min_donor_mass(v_ej: float, m_acc_max: float, model: str = "bedard", teff: float = 1.0e4,
                   with_kick: bool = True) -> float:
    """Smallest donor mass that can reach launch speed v_ej when the accretor
    is at most m_acc_max (and at least as heavy as the donor). The speed rises
    with both masses, so the lightest donor uses the heaviest accretor.
    Returns nan if even two m_acc_max stars are too slow."""
    lo = 0.20 if model == "bedard" else 0.05
    f = lambda m: ejection_speed(m, m_acc_max, model, teff, with_kick) - v_ej
    if f(lo) > 0:
        return lo
    if f(m_acc_max) < 0:
        return float("nan")
    return brentq(f, lo, m_acc_max, xtol=1e-4)


def max_speed(m_acc_max: float, model: str = "bedard", teff: float = 1.0e4,
              with_kick: bool = True) -> float:
    """Fastest launch with both stars at m_acc_max."""
    return float(ejection_speed(m_acc_max, m_acc_max, model, teff, with_kick))
