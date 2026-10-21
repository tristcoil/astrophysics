"""How cosmic rays make antiprotons: the production cross-section.

When a cosmic-ray proton with more than ~6 GeV of kinetic energy hits a
hydrogen atom in interstellar gas, it can make a proton-antiproton pair
(p + p -> p + p + p + pbar needs 6 m_p of kinetic energy in the lab). The
antiprotons AMS-02 sees are mostly such secondaries; the question is whether
anything is left over.

We use the published parametrisation of Korsmeier, Donato & di Mauro (2018,
PRD 97, 103019; "Param. II-B"), fitted to NA49, NA61, Dekkers and LHCb data:

* the Lorentz-invariant cross-section E d^3sigma/dp^3 for prompt pbar in
  p + p, as a function of the CM energy sqrt(s), x_R = E*/E*_max and p_T
  (their eq. 12, Winkler's 2017 form; Table III, Param. II);
* antineutrons (which decay to antiprotons) and antihyperons (Lambda-bar ->
  pbar pi+) added by the factor 2 + Delta_IS + 2 Delta_Lambda (eqs. 8-10);
* helium, as projectile or target, by the nuclear factor f^{A1 A2} with the
  projectile/target overlap functions of NA49 (eqs. 15-17, Tables VI, IX).

From these we build dsigma/dT_pbar for a cosmic ray of given kinetic energy
per nucleon, by integrating over the antiproton's angle in the lab, and the
source term q(T_pbar) for interstellar gas (their eqs. 1 and 3).
"""

from __future__ import annotations

import math

import numpy as np

MP = 0.938272          # GeV
MN_AMU = 0.931494      # GeV, per nucleon in 4He (m_He / 4)

#: Param. II (Winkler-type) best fit to p + p -> pbar + X, Table III of
#: Korsmeier et al. 2018. C4 was fixed to Winkler's value. (See the note on
#: the low-energy factor R in sigma_inv_pp: a typo in their eq. 14.)
PARAM_II = {"C1": 5.02e-2, "C2": 7.790, "C3": 1.649e-1, "C4": 3.800e-2, "C5": 4.74e-4,
            "C6": 3.70}
#: Total inelastic pp cross-section (mb), eq. 13: c1 + c2 ln sqrt(s) + c3 ln^2 sqrt(s).
SIGMA_IN = (30.9, -1.74, 0.71)
#: Antineutron and antihyperon factors (Winkler 2017, as adopted in eqs. 9-10).
ISOSPIN = {"c1": 0.114, "c2": 144.0**2, "c3": 0.51}
HYPERON = {"scale": 0.81, "c1": 0.31, "c2": 0.30, "c3": 146.0**2, "c4": 0.9}
#: Nuclear factor parameters D1, D2 for Param. II with NA49 pC + LHCb pHe
#: data ("II-B"), Table IX.
NUCLEAR = {"D1": 0.828, "D2": 0.145}
#: Projectile overlap function F_pro(x_F) (Table VI, from NA49 / Baatar et al. 2013);
#: F_tar(x_F) = F_pro(-x_F) = 1 - F_pro(x_F).
_XF = np.array([-0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0.0,
                0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25])
_FPRO = np.array([0.0, 0.0003, 0.0008, 0.0027, 0.010, 0.035, 0.110, 0.197, 0.295, 0.4, 0.5,
                  0.6, 0.705, 0.803, 0.890, 0.965, 0.990, 0.9973, 0.9992, 0.9997, 1.0])


#: Param. I (di Mauro et al. 2014 form, eq. 12; refit by Korsmeier et al.
#: 2018, Table III) with its own nuclear factor ("I-B", Table IX).
PARAM_I = {"C1": 3.50, "C2": 5.59, "C3": 4.00e-2, "C4": -0.251, "C5": 2.651, "C6": 3.78e-2,
           "C7": 4.3e-2, "C8": 2.695}
NUCLEAR_I = {"D1": 0.825, "D2": 0.167}
#: di Mauro et al. 2014 (PRD 90, 085017), Appendix B, Table V: PDG-form fits
#: sigma = Z + B ln^2(s/s_M) + Y1 (s_M/s)^eta1 - Y2 (s_M/s)^eta2 for the total
#: and elastic pp cross-sections; sigma_in = total - elastic (mb, GeV).
DIMAURO_PP = {"total": (2.06, 33.44, 13.53, 6.38, 0.324, 0.324),
              "elastic": (3.06, 144.98, 2.64, 137.27, 1.57, -4.65e-3)}
HBARC2_MB = 0.389379                    # (hbar c)^2 in mb GeV^2


def sigma_inelastic_pp_dimauro(sqrt_s):
    s = np.asarray(sqrt_s, float) ** 2
    out = 0.0
    for key, sign in (("total", 1.0), ("elastic", -1.0)):
        M, Z, Y1, Y2, e1, e2 = DIMAURO_PP[key]
        sm = (2 * MP + M) ** 2
        B = math.pi * HBARC2_MB / M**2
        out = out + sign * (Z + B * np.log(s / sm) ** 2 + Y1 * (sm / s) ** e1 - Y2 * (sm / s) ** e2)
    return out


def sigma_inv_pp_param1(sqrt_s, x_r, p_t, params=PARAM_I):
    """Prompt E d^3sigma/dp^3 (mb GeV^-2), di Mauro et al.'s form (Param. I)."""
    c = params
    val = (sigma_inelastic_pp_dimauro(sqrt_s) * np.clip(1 - x_r, 0, None) ** c["C1"]
           * np.exp(-c["C2"] * x_r)
           * (c["C3"] * sqrt_s ** c["C4"] * np.exp(-c["C5"] * p_t)
              + c["C6"] * sqrt_s ** c["C7"] * np.exp(-c["C8"] * p_t * p_t)))
    return np.where((x_r < 1) & (x_r > 0), val, 0.0)


#: Which parametrisation dsigma_dT uses: "II" (default) or "I".
PARAMETRISATION = {"name": "II"}


def f_pro(xf):
    return np.interp(xf, _XF, _FPRO, left=0.0, right=1.0)


def sigma_inelastic_pp(sqrt_s):
    c1, c2, c3 = SIGMA_IN
    ls = np.log(sqrt_s)
    return c1 + c2 * ls + c3 * ls * ls


def delta_isospin(s):
    c = ISOSPIN
    return c["c1"] / (1 + (s / c["c2"]) ** c["c3"])


def delta_hyperon(s):
    c = HYPERON
    return c["scale"] * (c["c1"] + c["c2"] / (1 + (c["c3"] / s) ** c["c4"]))


def sigma_inv_pp(sqrt_s, x_r, p_t, params=PARAM_II):
    """Prompt E d^3sigma/dp^3 for p + p -> pbar + X, in mb GeV^-2."""
    c = params
    s = sqrt_s**2
    e_max = (s - 8 * MP * MP) / (2 * sqrt_s)
    x_rmin = MP / e_max
    lnr = np.log(sqrt_s / (4 * MP))
    X = c["C4"] * lnr * lnr
    m_t = np.sqrt(p_t * p_t + MP * MP)
    low = sqrt_s < 10.0
    d = np.where(low, 10.0 - sqrt_s, 0.0)
    xr_in = np.clip(x_r, 0.0, 1.0)                     # outside (0, 1) the result is zeroed below
    # Korsmeier et al.'s eq. 14 prints (10 - sqrt s)^2 in the exponent, but
    # Winkler's original form (JCAP 02 (2017) 048, eq. 13) has the first power,
    # and only the first power reproduces Korsmeier et al.'s own Fig. 4
    # (T_p = 20 GeV: 1.30e-33 m^2/GeV at T_pbar = 2 GeV, 1.0e-35 at 10 GeV); the
    # square gives a non-monotonic pbar multiplicity near threshold and a
    # 3x bump in the Galactic source term at 10 GeV. We use the first power.
    R = np.where(low, (1 + c["C5"] * d**5) * np.exp(c["C6"] * d * (xr_in - x_rmin) ** 2), 1.0)
    pt_term = np.where(X > 1e-12, (1 + X * (m_t - MP)) ** (-1.0 / (c["C3"] * np.maximum(X, 1e-12))),
                       np.exp(-(m_t - MP) / c["C3"]))
    val = sigma_inelastic_pp(sqrt_s) * R * c["C1"] * np.clip(1 - x_r, 0, None) ** c["C2"] * pt_term
    return np.where((x_r < 1) & (x_r > 0), val, 0.0)


def nuclear_factor(a1, a2, xf, s):
    """f^{A1 A2}(x_F) (eq. 17): projectile A1, target A2 (A = 1 gives pp = 1)."""
    nuc = NUCLEAR if PARAMETRISATION["name"] == "II" else NUCLEAR_I
    d1, d2 = nuc["D1"], nuc["D2"]
    n1, n2 = a1 - (1 if a1 == 1 else a1 // 2), a2 - (1 if a2 == 1 else a2 // 2)
    dis = delta_isospin(s)
    fp = f_pro(xf)
    return (a1 * a2) ** d1 * (a1**d2 * (1 + n1 / a1 * dis) * fp
                              + a2**d2 * (1 + n2 / a2 * dis) * (1 - fp))


#: Lab angles for the angular integral: energetic antiprotons come out in a
#: cone of width ~ p_T / p (milliradians at 100 GeV), so the grid is
#: logarithmic down to 1e-7 rad and linear at wide angles.
_THETA = np.unique(np.concatenate([np.geomspace(1e-7, 0.3, 260), np.linspace(0.3, math.pi, 140)]))
_THETA = np.concatenate([[0.0], _THETA])


def dsigma_dT(t_proj, t_pbar, a1: int = 1, a2: int = 1, galaxy: bool = True):
    """Energy-differential cross-section dsigma/dT_pbar (mb / GeV) for a
    projectile of kinetic energy per nucleon t_proj (GeV) on a target at rest:
    p_pbar * 2 pi int dcos(theta) sigma_inv (Korsmeier et al. eq. 3), times
    the nuclear factor and, if ``galaxy``, the antineutron + antihyperon factor.
    Scalars t_proj, t_pbar."""
    e_p = t_proj + MP
    s = 2 * MP * MP + 2 * MP * e_p
    sqrt_s = math.sqrt(s)
    if sqrt_s <= 4 * MP:
        return 0.0
    gamma_cm = (e_p + MP) / sqrt_s
    bg_cm = math.sqrt(e_p * e_p - MP * MP) / sqrt_s
    beta_cm = bg_cm / gamma_cm
    e = t_pbar + MP
    p = math.sqrt(e * e - MP * MP)
    th = _THETA
    mu = np.cos(th)
    p_t = p * np.sin(th)
    e_star = gamma_cm * (e - beta_cm * p * mu)
    pl_star = gamma_cm * (p * mu - beta_cm * e)
    e_max = (s - 8 * MP * MP) / (2 * sqrt_s)
    x_r = e_star / e_max
    xf = 2 * pl_star / sqrt_s
    if PARAMETRISATION["name"] == "II":
        inv = sigma_inv_pp(sqrt_s, x_r, p_t)
    else:
        inv = sigma_inv_pp_param1(sqrt_s, x_r, p_t)
    if a1 > 1 or a2 > 1:
        inv = inv * nuclear_factor(a1, a2, xf, s)
    if galaxy:
        inv = inv * (2 + delta_isospin(s) + 2 * delta_hyperon(s))
    return 2 * math.pi * p * np.trapezoid(inv * np.sin(th), th)


def multiplicity_pp(p_lab: float, n_t: int = 400) -> float:
    """Mean prompt pbar multiplicity in p + p at beam momentum p_lab (GeV/c):
    int dT dsigma/dT / sigma_inel. A check against measured values."""
    e_p = math.hypot(p_lab, MP)
    t_proj = e_p - MP
    t = np.geomspace(1e-3, t_proj, n_t)
    ds = np.array([dsigma_dT(t_proj, ti, galaxy=False) for ti in t])
    sqrt_s = math.sqrt(2 * MP * MP + 2 * MP * e_p)
    return float(np.trapezoid(ds, t) / sigma_inelastic_pp(sqrt_s))


def source_term_kernel(t_pbar: np.ndarray, t_proj: np.ndarray, a1: int, a2: int) -> np.ndarray:
    """Matrix K[i, j] = dsigma/dT_pbar(t_proj[j] -> t_pbar[i]) in mb/GeV."""
    K = np.zeros((len(t_pbar), len(t_proj)))
    for j, tp in enumerate(t_proj):
        for i, tb in enumerate(t_pbar):
            if tb < tp:
                K[i, j] = dsigma_dT(tp, tb, a1, a2)
    return K
