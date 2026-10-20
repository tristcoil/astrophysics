"""Antimatter in cosmic rays: are the antiprotons just secondaries?

AMS-02 on the International Space Station counts about one antiproton for
every 5,000 protons at high energy. Cosmic-ray protons hitting interstellar
gas make antiprotons, so most or all of them are expected. To see how many,
we build the simplest model that uses the *same* gas for everything, the
"leaky box": cosmic rays wander in the Galaxy's magnetic fields and leak out
after crossing a column of gas, the grammage X(R) (g/cm^2), that depends on
their rigidity R = pc/Ze.

1. **The grammage from boron.** Boron is (almost) absent from stars, so
   cosmic-ray boron is made by carbon, nitrogen and oxygen breaking up on
   interstellar gas. The measured boron-to-carbon ratio therefore *measures*
   X(R). We fit it to AMS-02's B/C (2011-2018).
2. **The antiprotons.** With the same X(R), the antiproton flux follows from
   the measured proton and helium fluxes and the production cross-section
   (``pbar_production.py``), minus antiprotons lost by annihilating on the
   gas.
3. **Solar modulation.** Inside the heliosphere, cosmic rays lose energy
   climbing out against the solar wind; the force-field approximation shifts
   every spectrum by the modulation potential phi (CRDB lists the value for
   each data-taking period).

Steady-state leaky box, kinetic energy per nucleon T conserved in
fragmentation (straight-ahead approximation):

    N_i [1/X + 1/lambda_i] = sum_j N_j / lambda_(j->i) (+ primary source)

with lambda = m_ISM / sigma the interaction length in g/cm^2 for gas of
helium-to-hydrogen number ratio 0.1.
"""

from __future__ import annotations

import math

import numpy as np

import fetch_data as fd

MP = 0.938272
MU = 0.931494
#: Interstellar gas: helium-to-hydrogen number ratio and mass per hydrogen atom (g).
HE_TO_H = 0.1
M_ISM_PER_H = 1.6726e-24 * (1.0 + 4.0 * HE_TO_H)

#: Species: (Z, A) for rigidity <-> kinetic energy per nucleon conversions.
SPECIES = {"p": (1, 1), "pbar": (1, 1), "He": (2, 4), "B": (5, 10.67), "C": (6, 12),
           "N": (7, 14), "O": (8, 16), "Ne": (10, 20), "Mg": (12, 24), "Si": (14, 28)}
#: Boron is about two thirds 11B and one third 10B (mean A = 10.67).


def r_to_tn(R, species: str):
    """Rigidity (GV) -> kinetic energy per nucleon (GeV/n)."""
    z, a = SPECIES[species]
    m = MP if a == 1 else MU
    pn = np.asarray(R) * z / a
    return np.sqrt(pn**2 + m * m) - m


def tn_to_r(T, species: str):
    z, a = SPECIES[species]
    m = MP if a == 1 else MU
    T = np.asarray(T)
    return np.sqrt(T * (T + 2 * m)) * a / z


def beta_of_tn(T, species: str):
    z, a = SPECIES[species]
    m = MP if a == 1 else MU
    T = np.asarray(T)
    return np.sqrt(T * (T + 2 * m)) / (T + m)


def dR_dTn(T, species: str):
    z, a = SPECIES[species]
    m = MP if a == 1 else MU
    T = np.asarray(T)
    return (a / z) * (T + m) / np.sqrt(T * (T + 2 * m))


class Spectrum:
    """A measured flux J(R) (per m^2 s sr GV), interpolated in log-log and
    extended beyond the data with power laws fitted to the last decade.
    ``per_tn`` returns the flux per GeV/n at kinetic energy per nucleon T."""

    def __init__(self, R, J, species: str, phi_mv: float, err=None):
        o = np.argsort(R)
        self.R, self.J = np.asarray(R)[o], np.asarray(J)[o]
        self.err = None if err is None else np.asarray(err)[o]
        self.species = species
        self.phi_gv = phi_mv / 1000.0
        lr, lj = np.log(self.R), np.log(self.J)
        hi = self.R > self.R[-1] / 10
        self.slope_hi = np.polyfit(lr[hi], lj[hi], 1)[0]
        lo = self.R < self.R[0] * 3
        self.slope_lo = np.polyfit(lr[lo], lj[lo], 1)[0] if lo.sum() > 1 else 0.0

    def __call__(self, R):
        R = np.asarray(R, float)
        lr = np.log(R)
        lr0, lr1 = math.log(self.R[0]), math.log(self.R[-1])
        lj = np.interp(lr, np.log(self.R), np.log(self.J))
        lj = np.where(lr > lr1, math.log(self.J[-1]) + self.slope_hi * (lr - lr1), lj)
        lj = np.where(lr < lr0, math.log(self.J[0]) + self.slope_lo * (lr - lr0), lj)
        return np.exp(lj)

    def per_tn(self, T):
        """Top-of-atmosphere flux per (m^2 s sr GeV/n)."""
        return self(tn_to_r(T, self.species)) * dR_dTn(T, self.species)

    def lis_per_tn(self, T):
        """Interstellar flux by the force-field approximation (Gleeson & Axford
        1968): J_LIS(T) = J_TOA(T - Phi) T(T+2m) / ((T-Phi)(T-Phi+2m)), with
        Phi = (Z/A) phi per nucleon."""
        z, a = SPECIES[self.species]
        m = MP if a == 1 else MU
        T = np.asarray(T, float)
        big_phi = z / a * self.phi_gv
        t_toa = np.maximum(T - big_phi, 1e-4)
        return self.per_tn(t_toa) * T * (T + 2 * m) / (t_toa * (t_toa + 2 * m))


def modulate_ratio(ratio_lis_fn, T, phi_gv: float):
    """For two species with the same Z/A (pbar and p), the force field shifts
    the ratio in energy only: (pbar/p)_TOA(T) = (pbar/p)_LIS(T + phi)."""
    return ratio_lis_fn(np.asarray(T) + phi_gv)


def load_spectrum(stem: str, subexp: str, species: str) -> Spectrum:
    d = fd.load_crdb(stem)
    d = d[(d.subexp == subexp) & ~d.is_ul]
    return Spectrum(d.e_mean.values, d.value.values, species, float(d.phi_mv.iloc[0]),
                    d.err.values)


def load_ratio(stem: str, subexp: str):
    d = fd.load_crdb(stem)
    d = d[(d.subexp == subexp) & ~d.is_ul].sort_values("e_mean")
    return d


# ---------------------------------------------------------------------------
# Nuclear cross-sections (mb)
# ---------------------------------------------------------------------------

#: Boron production on hydrogen, ghost nuclei included (11C -> 11B,
#: 10C -> 10B), Reinert & Winkler 2018 (JCAP 01, 055), Table 2:
#: sigma(T) = s0 G^2 (T-Eth)^2 / ((T^2-M^2)^2 + G^2 M^2)
#:            + s1 (1 - Eth/T)^xi (1 + D / (1 + (Th/T)^2)),  T in MeV/n, Th = 2 GeV.
#: (Eth, s0, M, G, s1, xi, D)
RW18_BORON = {
    ("C", "11B"): (17.3, 330.0, 18.9, 23.5, 60.2, 0.47, -0.16),
    ("C", "10B"): (21.4, 34.1, 46.9, 46.9, 18.4, 0.75, -0.03),
    ("O", "11B"): (23.6, 240.0, 45.9, 34.2, 31.4, 4.0, -0.36),
    ("O", "10B"): (26.9, 104.0, 55.0, 29.0, 10.9, 1.3, 0.10),
    ("N", "11B"): (3.1, 193.0, 10.6, 7.3, 31.6, 9.6, -0.24),
    ("N", "10B"): (12.5, 360.0, 16.9, 14.1, 11.6, 3.7, -0.02),
    ("Ne", "11B"): (23.6, 0.0, 1.0, 1.0, 23.3, 6.7, 0.17),
    ("Ne", "10B"): (23.6, 0.0, 1.0, 1.0, 13.1, 14.3, 0.13),
    ("Mg", "11B"): (21.0, 0.0, 1.0, 1.0, 17.6, 7.4, 0.29),
    ("Mg", "10B"): (22.5, 0.0, 1.0, 1.0, 8.7, 10.4, 0.24),
    ("Si", "11B"): (26.7, 0.0, 1.0, 1.0, 15.1, 7.1, 0.30),
    ("Si", "10B"): (29.0, 0.0, 1.0, 1.0, 7.5, 6.4, 0.12),
    ("11B", "10B"): (10.1, 0.0, 1.0, 1.0, 49.0, 2.7, -0.04),
}
#: Helium-target factor for production (Reinert & Winkler's refit of the
#: Ferrando et al. 1988 form): exp(mu |dZ - f delta|^0.57), mu = 0.29 - 0.056 T,
#: delta = 1.93 - 0.45 T (T in GeV/n, frozen above 1.6).
RW18_HE_F = {"C": -0.35, "N": -0.8, "O": 0.04, "Ne": 0.24, "Mg": 0.45, "Si": 0.65, "11B": 0.0}
CHARGE = {"B": 5, "C": 6, "N": 7, "O": 8, "Ne": 10, "Mg": 12, "Si": 14, "11B": 5, "10B": 5}
MASSNUM = {"10B": 10, "11B": 11, "C": 12, "N": 14, "O": 16, "Ne": 20, "Mg": 24, "Si": 28}
#: Boron isotopes in cosmic rays: about one third 10B (Genolini et al. 2018);
#: AMS assumes an 11B fraction of 0.7 when converting to rigidity.
F_11B = 0.67
MB_PER_GCM2 = 1e-27 / (1.6726e-24 * (1.0 + 4.0 * HE_TO_H))   # (cm^2 per mb) / (g per H atom)


def sigma_rw18(channel, t_gevn):
    eth, s0, M, G, s1, xi, D = RW18_BORON[channel]
    T = np.maximum(np.asarray(t_gevn, float) * 1e3, eth * 1.0001)
    res = s0 * G**2 * (T - eth) ** 2 / ((T**2 - M**2) ** 2 + G**2 * M**2) if s0 > 0 else 0.0
    return res + s1 * (1 - eth / T) ** xi * (1 + D / (1 + (2000.0 / T) ** 2))


def he_factor_production(parent, t_gevn):
    t = np.minimum(np.asarray(t_gevn, float), 1.6)
    mu, delta = 0.29 - 0.056 * t, 1.93 - 0.45 * t
    dz = CHARGE[parent] - 5
    return np.exp(mu * np.abs(dz - RW18_HE_F[parent] * delta) ** 0.57)


def sigma_inelastic_letaw(a, t_gevn):
    """Total inelastic cross-section on hydrogen (Letaw, Silberberg & Tsao
    1983, ApJS 51, 271, eqs. 2 and 5)."""
    e = np.asarray(t_gevn, float) * 1e3
    he = 45.0 * a**0.7 * (1 + 0.016 * math.sin(5.3 - 2.63 * math.log(a)))
    return he * (1 - 0.62 * np.exp(-e / 200.0) * np.sin(10.9 * e ** (-0.28)))


def he_factor_inelastic(a):
    """sigma(He target) / sigma(H target) = 2.10 A^0.055 (Ferrando et al. 1988)."""
    return 2.10 * a**0.055


def boron_production_per_g(parent, t_gevn):
    """kappa_(parent -> B) in cm^2/g of interstellar gas (He/H = 0.1)."""
    tot = 0.0
    for iso in ("11B", "10B"):
        s = sigma_rw18((parent, iso), t_gevn)
        tot = tot + s * (1 + HE_TO_H * he_factor_production(parent, t_gevn))
    return tot * MB_PER_GCM2


def boron_loss_per_g(t_gevn):
    """Destruction of elemental boron: inelastic minus 11B -> 10B (stays boron)."""
    l11 = (sigma_inelastic_letaw(11, t_gevn) * (1 + HE_TO_H * he_factor_inelastic(11))
           - sigma_rw18(("11B", "10B"), t_gevn) * (1 + HE_TO_H * 1.29))
    l10 = sigma_inelastic_letaw(10, t_gevn) * (1 + HE_TO_H * he_factor_inelastic(10))
    return (F_11B * l11 + (1 - F_11B) * l10) * MB_PER_GCM2


def pbar_annihilation_per_g(t_gev):
    """Antiproton annihilation on interstellar gas, cm^2/g: Tan & Ng (1983) as
    given by Moskalenko et al. 2002 (ApJ 565, 280, eq. A3) below 14.64 GeV,
    36 T^-0.5 mb above (the branches cross there); helium x 4^(2/3)."""
    t = np.asarray(t_gev, float)
    low = 661.0 * (1 + 0.0115 * t ** (-0.774) - 0.948 * t**0.0151)
    sig = np.where(t < 14.64, low, 36.0 * t ** (-0.5))
    return sig * (1 + HE_TO_H * 4 ** (2 / 3)) * MB_PER_GCM2


# ---------------------------------------------------------------------------
# The leaky box: grammage from boron, antiprotons from the same grammage
# ---------------------------------------------------------------------------

#: AMS-02 data sets used (CRDB sub-experiment names). The 2011-2018 Physics
#: Reports set measures every species over the same period (phi = 637 MV).
AMS_PERIOD = "AMS02 (2011/05-2018/05)"
PARENTS = {"C": ("c_flux", AMS_PERIOD), "N": ("n_flux", AMS_PERIOD),
           "O": ("o_flux", AMS_PERIOD), "Ne": ("ne_flux", AMS_PERIOD),
           "Mg": ("mg_flux", AMS_PERIOD), "Si": ("si_flux", AMS_PERIOD)}
#: Hardening of the escape time seen in primaries and secondaries above
#: ~200-300 GV (AMS-02, PRL 120, 021101; break from Genolini et al. 2017).
R_BREAK, S_BREAK = 312.0, 0.1


def grammage(R, x0, delta, d_delta, beta=None):
    """X(R) = X0 beta (R / 10 GV)^-delta [1 + (R/R_b)^(d_delta/s)]^s  (g/cm^2)."""
    R = np.asarray(R, float)
    b = 1.0 if beta is None else beta
    return x0 * b * (R / 10.0) ** (-delta) * (1 + (R / R_BREAK) ** (d_delta / S_BREAK)) ** S_BREAK


class LeakyBox:
    """Everything that depends on the measured spectra, set up once."""

    def __init__(self):
        self.parents = {p: load_spectrum(stem, sub, p) for p, (stem, sub) in PARENTS.items()}
        self.phi_gv = self.parents["C"].phi_gv
        self.bc = load_ratio("b_over_c", AMS_PERIOD)
        self.proton = load_spectrum("h_flux", AMS_PERIOD, "p")
        self.helium = load_spectrum("he_flux", AMS_PERIOD, "He")

    # -- boron -------------------------------------------------------------
    def boron_lis(self, T, params, xs_scale: float = 1.0):
        """Interstellar boron flux per GeV/n at kinetic energy per nucleon T."""
        x = grammage(tn_to_r(T, "B"), *params, beta=beta_of_tn(T, "B"))
        prod = sum(self.parents[p].lis_per_tn(T) * boron_production_per_g(p, T)
                   for p in self.parents)
        return xs_scale * prod / (1.0 / x + boron_loss_per_g(T))

    def bc_toa(self, R, params, phi_gv=None, xs_scale: float = 1.0):
        """B/C at rigidity R (GV) as AMS measures it (both at the same R)."""
        phi = self.phi_gv if phi_gv is None else phi_gv
        out = []
        for s in ("B", "C"):
            z, a = SPECIES[s]
            T = r_to_tn(R, s)
            big = z / a * phi
            t_is = T + big
            if s == "B":
                j_is = self.boron_lis(t_is, params, xs_scale)
            else:
                # carbon's interstellar spectrum, rebuilt from its own TOA data
                j_is = self.parents["C"].lis_per_tn(t_is)
            m = MU
            j_toa = j_is * T * (T + 2 * m) / (t_is * (t_is + 2 * m))
            out.append(j_toa / dR_dTn(T, s))
        return out[0] / out[1]

    def fit_bc(self, r_min: float = 4.0, xs_scale: float = 1.0):
        from scipy import optimize

        d = self.bc[self.bc.e_mean >= r_min]
        R, y, e = d.e_mean.values, d.value.values, d.err.values

        def resid(p):
            return (self.bc_toa(R, p, xs_scale=xs_scale) - y) / e

        sol = optimize.least_squares(resid, x0=[10.0, 0.5, 0.15], x_scale=[1.0, 0.1, 0.05])
        J = sol.jac
        cov = np.linalg.inv(J.T @ J)
        chi2 = float(np.sum(sol.fun**2))
        return {"params": sol.x, "cov": cov, "chi2": chi2, "ndof": len(R) - 3, "R": R, "data": y,
                "err": e, "r_min": r_min, "xs_scale": xs_scale}

    # -- antiprotons -------------------------------------------------------
    def pbar_source_matrices(self, t_pbar: np.ndarray, n_proj: int = 160):
        """Kernels dsigma/dT for p+H, p+He, He+H, He+He on a projectile grid
        (GeV/n), cached on disk (they take a minute)."""
        import pickle

        import config as cfg
        import pbar_production as pp

        path = cfg.DERIVED_DIR / "pbar_kernels.pkl"
        t_proj = np.geomspace(6.0, 2e5, n_proj)
        if path.exists():
            d = pickle.load(open(path, "rb"))
            if np.allclose(d["t_pbar"], t_pbar) and np.allclose(d["t_proj"], t_proj):
                return d
        K = {ch: pp.source_term_kernel(t_pbar, t_proj, *ch) for ch in
             ((1, 1), (1, 4), (4, 1), (4, 4))}
        d = {"t_pbar": t_pbar, "t_proj": t_proj, "K": K}
        pickle.dump(d, open(path, "wb"))
        return d

    def pbar_lis(self, t_pbar, params, xs_scale: float = 1.0, src_scale: float = 1.0):
        """Interstellar antiproton flux per GeV (leaky box, same grammage as boron):
        J = X / (1 + X kappa_ann) * sum_ij (f_j / m_ISM) int dT_i J_i dsigma_ij/dT."""
        k = self.pbar_source_matrices(t_pbar)
        tp = k["t_proj"]
        jp = self.proton.lis_per_tn(tp)
        jhe = self.helium.lis_per_tn(tp)
        q = np.zeros_like(t_pbar)
        for (a1, a2), K in k["K"].items():
            j = jp if a1 == 1 else jhe
            f_t = 1.0 if a2 == 1 else HE_TO_H
            q += f_t * np.trapezoid(K * j[None, :], tp, axis=1)
        q = q * MB_PER_GCM2 * src_scale                     # per gram of ISM
        x = grammage(tn_to_r(t_pbar, "pbar"), *params, beta=beta_of_tn(t_pbar, "pbar"))
        x = x / xs_scale     # a boron cross-section scale is absorbed in the fitted X
        return x * q / (1 + x * pbar_annihilation_per_g(t_pbar))

    def pbar_over_p(self, R, params, phi_gv, **kw):
        """Antiproton-to-proton ratio at rigidity R (GV) for modulation phi."""
        t_grid = np.geomspace(0.1, 2000.0, 70)
        jb = self.pbar_lis(t_grid, params, **kw)
        T = r_to_tn(R, "p") + phi_gv                      # force field: same Z/A
        jb_t = np.exp(np.interp(np.log(T), np.log(t_grid), np.log(jb)))
        return jb_t / self.proton.lis_per_tn(T)
