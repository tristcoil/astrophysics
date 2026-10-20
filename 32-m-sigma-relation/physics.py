"""The physics behind the numbers: scales, feedback and energy.

Three ideas recur in the analysis:

1. **Sphere of influence.** A black hole of mass M dominates the motions of
   stars moving with velocity dispersion sigma only inside r_infl = G M / sigma^2.
   To measure M from stellar or gas motions, that radius must be resolved.

2. **Momentum-driven feedback (King 2003).** A black hole radiating at its
   Eddington limit pushes on the surrounding gas with force L_Edd / c. The
   push can sweep the galaxy's gas out of an isothermal halo (gas fraction
   f_g, velocity dispersion sigma) once it beats gravity. That happens at

       M_sigma = f_g kappa sigma^4 / (pi G^2),

   with kappa the electron-scattering opacity. The black hole then cuts off
   its own food supply: a prediction of both the slope (4) and, with no free
   parameter, the zero point.

3. **Energy budget.** Growing a black hole of mass M released roughly
   eps_r M c^2 of energy (eps_r ~ 0.1); binding a bulge of mass M_bulge takes
   about M_bulge sigma^2. Their ratio says whether the black hole *could*
   have shaped its galaxy.
"""

from __future__ import annotations

import numpy as np

import config as C


def r_influence_pc(m_bh, sigma_kms):
    """Sphere-of-influence radius G M / sigma^2 in parsecs."""
    return C.G_PC * np.asarray(m_bh, float) / np.asarray(sigma_kms, float) ** 2


def angular_size_arcsec(r_pc, dist_mpc):
    """Angle subtended by a physical radius r (pc) at distance D (Mpc)."""
    return np.asarray(r_pc, float) / (np.asarray(dist_mpc, float) * C.PC_PER_MPC) * C.ARCSEC_PER_RAD


def schwarzschild_radius_pc(m_bh):
    """Event-horizon radius 2 G M / c^2 of a non-spinning black hole, in parsecs."""
    return 2.0 * C.G_PC * np.asarray(m_bh, float) / C.C_KMS**2


def min_measurable_mass(sigma_kms, dist_mpc, resolution_arcsec=C.RESOLUTION_ARCSEC):
    """The smallest M whose sphere of influence spans ``resolution_arcsec`` at distance D."""
    r_pc = resolution_arcsec / C.ARCSEC_PER_RAD * np.asarray(dist_mpc, float) * C.PC_PER_MPC
    return r_pc * np.asarray(sigma_kms, float) ** 2 / C.G_PC


def kappa_es(hydrogen_x=C.HYDROGEN_X):
    """Electron-scattering opacity of fully ionised gas: (sigma_T / m_p) (1 + X) / 2 ~ 0.2 (1 + X) cm^2/g."""
    return C.SIGMA_T / C.M_P * (1.0 + hydrogen_x) / 2.0


def king_msigma(sigma_kms, f_gas=C.F_GAS_COSMIC, kappa=None):
    """King (2003) momentum-driven black hole mass, f_g kappa sigma^4 / (pi G^2), in Msun."""
    kappa = kappa_es() if kappa is None else kappa
    sig = np.asarray(sigma_kms, float) * C.KM_CM
    return f_gas * kappa * sig**4 / (np.pi * C.G_CGS**2) / C.MSUN_G


def eddington_luminosity(m_bh):
    """L_Edd = 4 pi G M m_p c / sigma_T (erg/s), hydrogen gas."""
    return 4.0 * np.pi * C.G_CGS * np.asarray(m_bh, float) * C.MSUN_G * C.M_P * C.C_CGS / C.SIGMA_T


def energy_ratio(m_bh, m_bulge, sigma_kms, efficiency=C.RADIATIVE_EFFICIENCY):
    """Energy radiated growing the black hole over the bulge's binding energy, ~ eps M c^2 / (M_b sigma^2)."""
    return efficiency * np.asarray(m_bh, float) * C.C_KMS**2 / (np.asarray(m_bulge, float)
                                                               * np.asarray(sigma_kms, float) ** 2)
