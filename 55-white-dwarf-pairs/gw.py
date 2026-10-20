"""Gravitational waves from circular binaries, and LISA's sensitivity to them.

Conventions (all SI unless stated; masses in solar masses where named ``m``):

* GW frequency f = 2 / P_orb (the quadrupole radiates at twice the orbital frequency).
* Chirp mass Mc = (m1 m2)^(3/5) / (m1 + m2)^(1/5).
* Strain amplitude (Kupfer et al. 2018, 2024, their eq. for script-A):
      A = 2 (G Mc)^(5/3) (pi f)^(2/3) / (c^4 d),
  with the two polarisations h+ = A (1 + cos^2 i) cos(Phi) and hx = 2 A cos(i) sin(Phi).
  (A face-on binary therefore has strain amplitude 2A, the "4 (G Mc)^(5/3)..." of
  textbooks such as Maggiore; edge-on it has A in h+ only.)
* LISA sensitivity: the analytic curve of Robson, Cornish & Liu (2019, CQG 36,
  105011), their eqs. (1)-(14): sky-, polarisation-averaged, both low-frequency TDI
  channels combined, plus their fit to the unresolved Galactic foreground.
* Sky- and polarisation-averaged SNR of a monochromatic source observed for T:
      rho^2 = 4 T A^2 [((1 + cos^2 i)/2)^2 + cos^2 i] / S_n(f),
  which averaged over inclination gives Robson et al.'s rho^2 = (16/5) T A^2 / S_n
  (their h_GB, eq. 27). The unit tests check both against published numbers.
"""

from __future__ import annotations

import numpy as np

import config as C

# ---------------------------------------------------------------------------
# Binary orbits and the quadrupole formula
# ---------------------------------------------------------------------------


def chirp_mass(m1, m2):
    """Chirp mass (same units as the inputs)."""
    m1, m2 = np.asarray(m1, float), np.asarray(m2, float)
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def semi_major_axis(P, mtot):
    """Kepler's third law: separation (m) for period P (s) and total mass (Msun)."""
    return (C.GM_SUN * np.asarray(mtot) * (np.asarray(P) / (2 * np.pi)) ** 2) ** (1 / 3)


def pdot_gw(P, mc):
    """Orbital-period derivative (s/s) from gravitational radiation alone (circular orbit).

    dP/dt = -(192 pi / 5) (2 pi G Mc / (c^3 P))^(5/3)    (Peters & Mathews 1963)
    """
    tc = C.T_SUN * np.asarray(mc, float)
    return -(192 * np.pi / 5) * (2 * np.pi * tc / np.asarray(P, float)) ** (5 / 3)


def chirp_mass_from_pdot(P, pdot):
    """Invert ``pdot_gw``: the chirp mass (Msun) implied if all the decay is gravitational."""
    x = (-np.asarray(pdot, float) / (192 * np.pi / 5)) ** 0.6  # = 2 pi G Mc / (c^3 P)
    return x * np.asarray(P, float) / (2 * np.pi) / C.T_SUN


def fdot_gw(f, mc):
    """GW-frequency derivative (Hz/s): (96/5) pi^(8/3) (G Mc/c^3)^(5/3) f^(11/3)."""
    tc = C.T_SUN * np.asarray(mc, float)
    return 96 / 5 * np.pi ** (8 / 3) * tc ** (5 / 3) * np.asarray(f, float) ** (11 / 3)


def amplitude(f, mc, d):
    """Strain amplitude script-A for GW frequency f (Hz), chirp mass (Msun), distance (m)."""
    tc = C.T_SUN * np.asarray(mc, float)
    return 2 * tc ** (5 / 3) * (np.pi * np.asarray(f, float)) ** (2 / 3) * C.C_LIGHT / np.asarray(d, float)


def inclination_factor(incl):
    """[((1 + cos^2 i)/2)^2 + cos^2 i]: 2 face-on, 1/4 edge-on, 4/5 averaged over the sky."""
    c2 = np.cos(np.asarray(incl, float)) ** 2
    return ((1 + c2) / 2) ** 2 + c2


def time_to_coalescence(P, m1, m2):
    """Time (s) for point masses on a circular orbit to spiral together from period P."""
    tc = C.T_SUN * chirp_mass(m1, m2)
    f = 2.0 / np.asarray(P, float)
    return 5 / 256 * tc ** (-5 / 3) * (np.pi * f) ** (-8 / 3)


def time_between_separations(a0, a1, m1, m2):
    """Time (s) for GW emission to shrink a circular orbit from a0 to a1 (m)."""
    beta = 64 / 5 * C.GM_SUN**3 * m1 * m2 * (m1 + m2) / C.C_LIGHT**5
    return (np.asarray(a0) ** 4 - np.asarray(a1) ** 4) / (4 * beta)


def shklovskii_pdot(P, mu_mas_yr, d):
    """Apparent P-dot from the transverse motion (Shklovskii 1970): P mu^2 d / c."""
    mu = np.asarray(mu_mas_yr) * (np.pi / 180 / 3600e3) / C.YEAR  # rad/s
    return np.asarray(P) * mu**2 * np.asarray(d) / C.C_LIGHT


# ---------------------------------------------------------------------------
# LISA: Robson, Cornish & Liu (2019)
# ---------------------------------------------------------------------------

LISA_ARM = 2.5e9  # m
F_STAR = 19.09e-3  # Hz, transfer frequency c / (2 pi L)

#: Galactic confusion-noise fit (their Table 1): alpha, beta, kappa, gamma, f_k for an
#: observation of T years. Amplitude fixed at 9e-45 (two-channel normalisation).
CONFUSION_FIT = {
    0.5: (0.133, 243.0, 482.0, 917.0, 0.00258),
    1.0: (0.171, 292.0, 1020.0, 1680.0, 0.00215),
    2.0: (0.165, 299.0, 611.0, 1340.0, 0.00173),
    4.0: (0.138, -221.0, 521.0, 1680.0, 0.00113),
}
CONFUSION_AMPLITUDE = 9e-45


def p_oms(f):
    """Single-link optical metrology noise (m^2/Hz)."""
    f = np.asarray(f, float)
    return (1.5e-11) ** 2 * (1 + (2e-3 / f) ** 4)


def p_acc(f):
    """Single test-mass acceleration noise ((m/s^2)^2/Hz)."""
    f = np.asarray(f, float)
    return (3e-15) ** 2 * (1 + (0.4e-3 / f) ** 2) * (1 + (f / 8e-3) ** 4)


def sn_instrument(f):
    """Sky-averaged strain sensitivity from instrument noise alone (1/Hz), eq. (13)."""
    f = np.asarray(f, float)
    return (10 / (3 * LISA_ARM**2)
            * (p_oms(f) + 2 * (1 + np.cos(f / F_STAR) ** 2) * p_acc(f) / (2 * np.pi * f) ** 4)
            * (1 + 0.6 * (f / F_STAR) ** 2))


def _confusion_at(f, key):
    alpha, beta, kappa, gamma, fk = CONFUSION_FIT[key]
    return (CONFUSION_AMPLITUDE * f ** (-7 / 3) * np.exp(-(f**alpha) + beta * f * np.sin(kappa * f))
            * (1 + np.tanh(gamma * (fk - f))))


def sn_confusion(f, t_obs_yr: float = 4.0):
    """Robson et al.'s fit to the unresolved Galactic foreground (1/Hz) after t_obs_yr.

    The fit exists for 0.5, 1, 2 and 4 yr. In between, log S_c is interpolated linearly
    in log T; shorter observations use the 0.5-yr curve (the foreground is only
    subtracted as sources get resolved, so it is higher early) and longer ones the 4-yr.
    """
    f = np.asarray(f, float)
    keys = sorted(CONFUSION_FIT)
    t = min(max(t_obs_yr, keys[0]), keys[-1])
    if t in CONFUSION_FIT:
        return _confusion_at(f, t)
    lo = max(k for k in keys if k < t)
    hi = min(k for k in keys if k > t)
    w = np.log(t / lo) / np.log(hi / lo)
    s_lo, s_hi = _confusion_at(f, lo), _confusion_at(f, hi)
    tiny = 1e-300
    return np.exp((1 - w) * np.log(s_lo + tiny) + w * np.log(s_hi + tiny))


def sn_total(f, t_obs_yr: float = 4.0, confusion: bool = True):
    s = sn_instrument(f)
    return s + sn_confusion(f, t_obs_yr) if confusion else s


def snr(f, amp, t_obs_s, incl=None, sn=None, t_obs_yr=None):
    """Sky- and polarisation-averaged SNR; inclination-averaged too if ``incl`` is None."""
    if sn is None:
        sn = sn_total(f, t_obs_yr if t_obs_yr is not None else t_obs_s / C.YEAR)
    fac = 4 / 5 if incl is None else inclination_factor(incl)
    return np.sqrt(4 * t_obs_s * np.asarray(amp) ** 2 * fac / sn)


def characteristic_strain_source(f, amp, t_obs_s, incl=None):
    """Height of a monochromatic source on a sqrt(f S_n) plot, chosen so SNR = h_src / h_noise."""
    fac = 4 / 5 if incl is None else inclination_factor(incl)
    return np.sqrt(4 * fac * np.asarray(f) * t_obs_s) * np.asarray(amp)


def characteristic_strain_noise(f, t_obs_yr: float = 4.0, confusion: bool = True):
    return np.sqrt(np.asarray(f) * sn_total(f, t_obs_yr, confusion))
