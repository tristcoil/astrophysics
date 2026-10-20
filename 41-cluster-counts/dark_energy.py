"""Dark energy's brake on structure growth.

Take two universes with the same physical matter and baryon densities
(Omega_m h^2, Omega_b h^2), the same primordial fluctuations (A_s, n_s) and
the same CMB temperature today. They are indistinguishable at z ~ 1100: the
CMB of each shows the same density ripples at recombination. One also contains
a cosmological constant (flat LCDM with Planck 2018's parameters); the other
does not (flat and matter-only, so h = sqrt(Omega_m h^2) ~ 0.38).

In the matter-only universe, density contrasts grow in proportion to the scale
factor forever. With Lambda, the expansion starts to accelerate at z ~ 0.6
and growth slows, then nearly stops. Same seeds, same matter: the only
difference is the late-time expansion, so the ratio of the two growth factors
since z = 1100 isolates what dark energy does. Because the most massive halos
are exponentially rare, a modest brake on growth removes most of them.

Masses here are M200m (200 times the mean matter density), so the two
universes, with equal mean density, share the same mass definition. Abundances
are per comoving Gpc^3.
"""

from __future__ import annotations

import numpy as np

import cosmology as CO
import literature as L

Z_CMB = 1100.0


def universes():
    p = L.PLANCK2018_BESTFIT
    omh2 = p["ombh2"] + p["omch2"] + 0.06 / 93.14  # + one 0.06 eV neutrino, as matter today
    h = 0.6732  # Planck 2018 best fit H0 = 67.32 (same point)
    lcdm = CO.Cosmology(Om=omh2 / h**2, sigma8=p["sigma8"], h=h, Ob=p["ombh2"] / h**2, ns=p["ns"])
    h0 = np.sqrt(omh2 + CO.OMEGA_R_H2)
    eds = CO.Cosmology(Om=omh2 / h0**2, sigma8=p["sigma8"], h=h0, Ob=p["ombh2"] / h0**2, ns=p["ns"], Ode=0.0)
    return lcdm, eds


def brake() -> dict:
    lcdm, eds = universes()
    z = np.array([0.0, 0.25, 0.5, 1.0, 2.0, 3.0, Z_CMB])
    D_l = lcdm.growth_ode(z)
    D_e = eds.growth_ode(z)
    # growth since recombination, ratio of the two universes at each z
    since_l = D_l / D_l[-1]
    since_e = D_e / D_e[-1]
    ratio = since_l / since_e
    # sigma_8 of the no-Lambda universe, in the same *physical* sphere
    # (8/h_LCDM Mpc): same seeds and transfer function, so it scales by the
    # growth ratio at z = 0.
    sigma8_noL_phys = lcdm.sigma8 / ratio[0]
    # Halo abundances: same P(k) shape in Mpc^-1 (EH98 depends on omega_m,
    # omega_b only), amplitudes tied at z = 1100.
    eds_same_seeds = eds.with_(sigma8=1.0)
    # normalise EdS so that sigma(R_phys = 8/h_L Mpc, z=0) = sigma8_noL_phys
    R8 = 8.0 / lcdm.h
    M8 = 4 / 3 * np.pi * eds.rho_m0 * R8**3
    s_unit, _ = eds_same_seeds.sigma_M(np.array([M8]), np.array([0.0]))
    eds_n = eds.with_(sigma8=float(sigma8_noL_phys / s_unit[0, 0]))
    M = np.geomspace(1e13, 5e15, 200)
    zz = np.array([0.0, 0.5, 1.0])
    out = {"z": z.tolist(), "D_lcdm": D_l.tolist(), "D_eds": D_e.tolist(), "ratio_since_cmb": ratio.tolist(),
           "sigma8_lcdm": lcdm.sigma8, "sigma8_noL_phys": float(sigma8_noL_phys),
           "h_lcdm": lcdm.h, "h_noL": eds.h, "Om_lcdm": lcdm.Om}
    nl = lcdm.dndlnM_mean(M, zz)
    ne = eds_n.dndlnM_mean(M, zz)
    dlnm = np.gradient(np.log(M))
    cum = lambda n: np.cumsum((n * dlnm)[:, ::-1], axis=1)[:, ::-1] * 1e9  # noqa: E731  per Gpc^3, N(>M)
    out.update({"M": M, "zz": zz, "n_gt_lcdm": cum(nl), "n_gt_noL": cum(ne)})
    for mcut in (1e14, 1e15):
        i = np.searchsorted(M, mcut)
        out[f"boost_gt_{mcut:.0e}"] = (out["n_gt_noL"][:, i] / out["n_gt_lcdm"][:, i]).tolist()
    # growth histories for the figure
    a = np.geomspace(1e-3, 1.0, 400)
    zc = 1 / a - 1
    out["a"] = a
    out["Da_lcdm"] = lcdm.growth_ode(zc) / lcdm.growth_ode(np.array([Z_CMB]))[0] * (1 / (1 + Z_CMB))
    out["Da_noL"] = eds.growth_ode(zc) / eds.growth_ode(np.array([Z_CMB]))[0] * (1 / (1 + Z_CMB))
    return out
