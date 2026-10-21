"""How many survivors should Gaia see, if every Type Ia supernova made one?

A forward model in the spirit of Shen et al. (2018, section 2.5) and El-Badry
et al. (2023, section 8.2), but with Gaia's real parallax precision and the
exact selection used by the search:

1. Explosions happen at the Milky Way's Type Ia rate (Li et al. 2011:
   0.54 per century), placed like the stars of the disc (exponential in
   radius and height).
2. Before the explosion the binary orbits with the disc (the potential's
   circular speed); the survivor leaves in a random direction at its launch
   speed and coasts in a straight line (over a few Myr the Galaxy's pull
   changes 2,000 km/s by well under 1%).
3. It stays visible for a time tau at absolute magnitude M_G, then fades
   out of reach. Neither number is known from theory -- they are what we
   scan over.
4. Dust: an exponential layer (scale height 120 pc, A_V = 1 mag/kpc in the
   midplane), which reddens and dims stars toward the plane.
5. Gaia measures each star with the parallax error of a real Gaia source of
   the same magnitude (drawn from a random all-sky subsample), and the star
   counts as found if it passes the search's cuts.

Expected number found = rate x tau x (fraction of live survivors that pass).
Dividing the number actually found by this gives the fraction of Type Ia
supernovae that must leave a survivor like the observed ones.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table

from config import (A_G_PER_A_V, DISC_SCALE_HEIGHT_KPC, DISC_SCALE_LENGTH_KPC,
                    DUST_AV_PER_KPC, DUST_SCALE_HEIGHT_KPC, GAIA_DIR, K_KMS, R0_KPC,
                    V_SUN_KMS, Z_SUN_KPC)
from kinematics import E_BPRP_PER_EBV, KPC_PER_KMS_MYR, rotation_velocity

#: E(BP-RP) / A_V for R_V = 3.1 (Casagrande & VandenBerg 2018 coefficients).
E_BPRP_PER_AV = E_BPRP_PER_EBV / 3.1


class ParallaxErrors:
    """Draw parallax errors of real Gaia DR3 sources at a given G."""

    def __init__(self, table: Table | None = None):
        t = table if table is not None else Table.read(GAIA_DIR / "precision.fits")
        ok = np.isfinite(t["parallax_error"]) & np.isfinite(t["phot_g_mean_mag"])
        g = np.asarray(t["phot_g_mean_mag"])[ok]
        e = np.asarray(t["parallax_error"])[ok]
        order = np.argsort(g)
        self.g, self.e = g[order], e[order]

    def draw(self, G, rng, window: float = 0.1) -> np.ndarray:
        G = np.asarray(G, float)
        lo = np.searchsorted(self.g, G - window)
        hi = np.searchsorted(self.g, G + window)
        hi = np.maximum(hi, lo + 1)
        idx = np.minimum(lo + (rng.random(G.size) * (hi - lo)).astype(int), self.g.size - 1)
        return self.e[idx]

    def median(self, G_bins):
        out = []
        for a, b in zip(G_bins[:-1], G_bins[1:]):
            s = (self.g >= a) & (self.g < b)
            out.append(np.median(self.e[s]) if s.sum() > 20 else np.nan)
        return np.array(out)


def sample_births(n, rng):
    """Disc positions (3, n) kpc in the astropy Galactocentric frame, with
    surface density ~ exp(-R/h_R) (so p(R) ~ R exp(-R/h_R)) and density
    ~ exp(-|z|/h_z)."""
    R = rng.gamma(2.0, DISC_SCALE_LENGTH_KPC, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    z = rng.laplace(0.0, DISC_SCALE_HEIGHT_KPC, n)
    return np.array([R * np.cos(phi), R * np.sin(phi), z])


def isotropic(n, rng):
    u = rng.uniform(-1, 1, n)
    phi = rng.uniform(0, 2 * np.pi, n)
    s = np.sqrt(1 - u ** 2)
    return np.array([s * np.cos(phi), s * np.sin(phi), u])


def dust_av(sun, pos, n_steps: int = 48):
    """A_V (mag) from the Sun to positions pos (3, n) through an exponential
    dust layer, by midpoint quadrature along each line of sight."""
    d = np.linalg.norm(pos - sun[:, None], axis=0)
    f = (np.arange(n_steps) + 0.5) / n_steps
    z = sun[2] + np.outer(f, pos[2] - sun[2])                   # (steps, n)
    return DUST_AV_PER_KPC * d * np.exp(-np.abs(z) / DUST_SCALE_HEIGHT_KPC).mean(axis=0)


def simulate(*, v_range, M_G0, color0, tau_myr, n, rng, plx_errors: ParallaxErrors,
             vcirc_fn, selection="union"):
    """Simulate n survivors alive today (ages uniform in [0, tau]) and return
    the fraction that the search would find, plus the found ones' properties."""
    age = rng.uniform(0, tau_myr, n)
    x0 = sample_births(n, rng)
    v = rotation_velocity(x0, vcirc_fn) + rng.uniform(*v_range, n) * isotropic(n, rng)
    x = x0 + v * age * KPC_PER_KMS_MYR
    sun = np.array([-R0_KPC, 0.0, Z_SUN_KPC])
    rel = x - sun[:, None]
    d = np.linalg.norm(rel, axis=0)
    los = rel / d
    vrel = v - np.array(V_SUN_KMS)[:, None]
    rv = (vrel * los).sum(axis=0)
    vtan = np.sqrt(np.maximum((vrel ** 2).sum(axis=0) - rv ** 2, 0))
    pm = vtan / (K_KMS * d)
    av = dust_av(sun, x)
    G = M_G0 + A_G_PER_A_V * av + 5 * np.log10(100.0 * d)
    bprp = color0 + E_BPRP_PER_AV * av
    s_plx = plx_errors.draw(G, rng)
    plx = 1.0 / d + s_plx * rng.standard_normal(n)
    with np.errstate(divide="ignore", invalid="ignore"):
        v_face = K_KMS * pm / plx
        v_low = K_KMS * pm / (plx + s_plx)
    found = ((G < 20) & (pm > 50) & (bprp < 0.5)
             & ((v_face > 600) | (plx < 0)) & ((v_low > 400) | (plx + s_plx < 0)))
    if selection == "elbadry":
        found &= plx / s_plx < 5
    return dict(frac=found.mean(), d=d[found], G=G[found], b=np.degrees(np.arcsin(los[2]))[found],
                age=age[found], n_found=int(found.sum()))


def expected_found(rate_per_yr, tau_myr, frac):
    """Survivors Gaia should find if a fraction 1 of SNe Ia makes one."""
    return rate_per_yr * tau_myr * 1e6 * frac


def poisson_interval(k: int, cl: float = 0.6827):
    """Central (Garwood) interval on a Poisson mean given k events."""
    from scipy.stats import chi2
    a = (1 - cl) / 2
    lo = 0.0 if k == 0 else chi2.ppf(a, 2 * k) / 2
    hi = chi2.ppf(1 - a, 2 * k + 2) / 2
    return lo, hi
