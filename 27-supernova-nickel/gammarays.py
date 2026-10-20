"""Weighing SN 2014J's nickel with gamma rays.

Each 56Co decay emits an 847 keV photon (99.9% of decays) and, two times in
three, a 1238 keV photon. Count them at Earth and you count the cobalt, and so
the nickel it came from. Dust does not dim MeV gamma rays, so this route skips
the large and uncertain dust correction of the optical light.

What a gamma-ray telescope sees is reduced by the ejecta themselves: a line
photon that Compton-scatters even once loses energy and leaves the line. The
escape probability is computed here for a simple, standard description of a
Type Ia explosion (1.38 Msun of ejecta, 1.2e51 erg kinetic energy, exponential
density profile like the W7 model), with the nickel in a shell outside a core
of stable iron-group elements. A uniform-density sphere is the alternative used
to show how much the answer depends on that choice.
"""

from __future__ import annotations

import numpy as np

import config as C
import nickel as N

R_E2 = (C.SIGMA_T * 3.0 / (8.0 * np.pi))  # classical electron radius squared, cm^2

#: INTEGRAL/SPI line fluxes, 1e-4 photons/cm^2/s.
#: Churazov et al. 2014 (Nature 512, 406), Extended Data Table 1: averaged over
#: days 50-100 after explosion (median day 75), 1-sigma errors.
CHURAZOV2014 = {"t_range": (50.0, 100.0), 847: (2.34, 0.74), 1238: (2.78, 0.74)}

#: Diehl et al. 2015 (A&A 574, A72), Table 1: four epochs (days since
#: explosion). The paper quotes errors only for epoch 3 (the brightest), in
#: the text; for the other epochs we scale that error by sqrt(duration ratio),
#: since INTEGRAL observed almost continuously (an approximation, noted).
DIEHL2015 = {
    "epochs": [(16.3, 41.3), (41.3, 66.3), (66.3, 99.1), (134.8, 164.0)],
    847: [2.34, 2.74, 3.65, 1.90],
    1238: [0.91, 1.11, 2.27, 0.38],
    "err_epoch3": {847: 1.21, 1238: 0.69},
}

#: Fiducial Type Ia ejecta (Chandrasekhar-mass explosion such as W7 or the
#: delayed detonations Churazov et al. compared with).
M_EJ_MSUN = 1.38
E_K_ERG = 1.2e51
Y_E = 0.5  # electrons per nucleon (Fe-group and intermediate-mass elements)
M_STABLE_CORE_MSUN = 0.1  # stable iron-group core inside the 56Ni (W7-like)


def klein_nishina(e_kev):
    """Total Compton cross-section per electron (cm^2) at photon energy E."""
    k = np.asarray(e_kev, float) / C.ME_C2_KEV
    lg = np.log1p(2 * k)
    term = ((1 + k) / k**2) * (2 * (1 + k) / (1 + 2 * k) - lg / k) + lg / (2 * k) \
        - (1 + 3 * k) / (1 + 2 * k) ** 2
    return 2 * np.pi * R_E2 * term


def line_opacity(e_kev) -> float:
    """Mass opacity (cm^2/g) for removing a photon from the line: any Compton scatter."""
    return float(klein_nishina(e_kev) * Y_E / C.AMU_G)


class Ejecta:
    """Homologous ejecta, density rho(v, t) = rho_0 (t_ref/t)^3 g(v / v_scale)."""

    def __init__(self, profile: str = "exponential", m_ej_msun: float = M_EJ_MSUN,
                 e_k: float = E_K_ERG, m_core_msun: float = M_STABLE_CORE_MSUN,
                 m_ni_msun: float = 0.6):
        self.profile = profile
        m = m_ej_msun * C.MSUN_G
        self.m = m
        if profile == "exponential":
            self.v_scale = np.sqrt(e_k / (6.0 * m))  # E_k = 6 M v_e^2
            self.rho0_v3 = m / (8.0 * np.pi * self.v_scale**3)  # rho * t^3 at v = 0
        elif profile == "uniform":
            self.v_scale = np.sqrt(10.0 * e_k / (3.0 * m))  # edge velocity
            self.rho0_v3 = m / (4.0 / 3.0 * np.pi * self.v_scale**3)
        else:
            raise ValueError(profile)
        # Nickel occupies the mass shell [m_core, m_core + m_ni].
        self.v_in = self.velocity_at_mass(m_core_msun * C.MSUN_G)
        self.v_out = self.velocity_at_mass(min(m_core_msun + m_ni_msun, 0.999 * m_ej_msun) * C.MSUN_G)

    def g(self, x):
        if self.profile == "exponential":
            return np.exp(-x)
        return (x <= 1.0).astype(float)

    def enclosed_mass(self, v):
        x = np.asarray(v, float) / self.v_scale
        if self.profile == "exponential":
            return self.m * (1.0 - np.exp(-x) * (1.0 + x + 0.5 * x * x))
        return self.m * np.clip(x, 0, 1) ** 3

    def velocity_at_mass(self, m_target):
        lo, hi = 0.0, 50.0 * self.v_scale
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self.enclosed_mass(mid) < m_target else (lo, mid)
        return 0.5 * (lo + hi)

    def column_t2(self, n_v: int = 48, n_mu: int = 48, n_s: int = 600) -> tuple[np.ndarray, np.ndarray]:
        """Column density x t^2 (g cm^-2 s^2) along rays from the nickel shell.

        Returns the columns for emission points (v, mu) and their weights
        (nickel mass x solid angle), so that any escape-type average is a
        weighted mean of a function of the column.
        """
        v = np.linspace(self.v_in, self.v_out, n_v + 1)
        vc = 0.5 * (v[1:] + v[:-1])
        dm = np.diff(self.enclosed_mass(v))  # nickel mass per shell (uniform Ni fraction)
        mu, wmu = np.polynomial.legendre.leggauss(n_mu)
        s_max = (30.0 if self.profile == "exponential" else 2.0) * self.v_scale
        s = np.linspace(0.0, s_max, n_s)
        vv = np.sqrt(vc[:, None, None] ** 2 + s[None, None, :] ** 2
                     + 2 * vc[:, None, None] * s[None, None, :] * mu[None, :, None])
        rho_t3 = self.rho0_v3 * self.g(vv / self.v_scale)
        col = np.trapezoid(rho_t3, s, axis=2)  # Int rho t^3 d(v) = column * t^2
        weights = dm[:, None] * wmu[None, :] / 2.0
        return col, weights / weights.sum()

    def escape_probability(self, t_days, e_kev):
        """Probability that a line photon emitted at time t leaves without scattering."""
        col, w = self.column_t2()
        kappa = line_opacity(e_kev)
        t = np.atleast_1d(np.asarray(t_days, float)) * C.DAY_S
        tau = kappa * col[None, :, :] / t[:, None, None] ** 2
        return np.sum(w[None] * np.exp(-tau), axis=(1, 2))

    def t0_days(self, kappa_eff: float = 0.025) -> float:
        """Gamma-ray escape time T0 implied by this ejecta for an effective opacity.

        Late on, the deposited fraction is ~ kappa_eff <column> = (T0/t)^2.
        kappa_eff ~ 0.025-0.03 cm^2/g is the standard effective absorption
        opacity for 56Co gamma rays (Sutherland & Wheeler 1984; Swartz et al. 1995).
        """
        col, w = self.column_t2()
        return float(np.sqrt(kappa_eff * np.sum(w * col)) / C.DAY_S)


def line_flux(t_days, m_ni: float, line: int, dist_cm: float, ejecta: Ejecta | None = None):
    """Line photon flux at Earth (photons/cm^2/s) from M_Ni of nickel at time t.

    Without ``ejecta`` the ejecta are transparent (a lower limit on the
    nickel when inverted).
    """
    t = np.atleast_1d(np.asarray(t_days, float))
    n_co = N.nuclei(t, m_ni)["Co"]
    rate = N.LAMBDA_CO / C.DAY_S * n_co * C.CO56_LINES[line]["yield"]
    flux = rate / (4.0 * np.pi * dist_cm**2)
    if ejecta is not None:
        flux = flux * ejecta.escape_probability(t, C.CO56_LINES[line]["energy_kev"])
    return flux


def averaged_line_flux(t_range, m_ni: float, line: int, dist_cm: float, ejecta=None, n: int = 26):
    """Mean of ``line_flux`` over an observing window (uniform exposure)."""
    t = np.linspace(*t_range, n)
    return float(np.mean(line_flux(t, m_ni, line, dist_cm, ejecta)))


def _lines_to_mass(dist_cm: float, ejecta: Ejecta | None) -> dict:
    out = {}
    for line in (847, 1238):
        f_obs, f_err = (v * 1e-4 for v in CHURAZOV2014[line])
        unit = averaged_line_flux(CHURAZOV2014["t_range"], 1.0, line, dist_cm, ejecta)
        out[line] = (f_obs / unit, f_err / unit)
    w = np.array([1 / out[k][1] ** 2 for k in (847, 1238)])
    m = np.array([out[k][0] for k in (847, 1238)])
    out["combined"] = (float(np.sum(w * m) / w.sum()), float(1 / np.sqrt(w.sum())))
    return out


def nickel_from_lines(dist_cm: float, profile: str | None, **ejecta_kw) -> dict:
    """Nickel mass from the Churazov et al. (2014) line fluxes, each line and combined.

    The model flux is linear in M_Ni, so each line gives M = F_obs / F_model(1 Msun);
    the two lines are combined with inverse-variance weights. Because the
    nickel shell's extent depends on the nickel mass, the ejecta model is
    iterated until the two agree. ``profile=None`` means transparent ejecta.
    """
    if profile is None:
        return _lines_to_mass(dist_cm, None)
    m_ni = 0.6
    for _ in range(20):
        ej = Ejecta(profile, m_ni_msun=m_ni, **ejecta_kw)
        out = _lines_to_mass(dist_cm, ej)
        new = out["combined"][0]
        if abs(new - m_ni) < 1e-3:
            break
        m_ni = min(new, 1.2)
    out["escape_75d"] = {line: float(ej.escape_probability(75.0, C.CO56_LINES[line]["energy_kev"])[0])
                         for line in (847, 1238)}
    out["ejecta"] = ej
    return out


def diehl_points() -> dict:
    """Diehl et al. (2015) epoch fluxes with the scaled errors described above."""
    ep = np.array(DIEHL2015["epochs"])
    dur = ep[:, 1] - ep[:, 0]
    out = {"t_mid": ep.mean(axis=1), "t_lo": ep[:, 0], "t_hi": ep[:, 1]}
    for line in (847, 1238):
        err3 = DIEHL2015["err_epoch3"][line]
        out[line] = np.array(DIEHL2015[line]) * 1e-4
        out[f"{line}_err"] = err3 * np.sqrt(dur[2] / dur) * 1e-4
    return out
