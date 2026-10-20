"""When the photons let go: recombination redshifts from CAMB, rulers from our integrals.

Two moments around recombination set the two rulers this project uses:

* **z_d, the drag epoch** (~1060): when photons stop pushing baryons around.
  The sound waves stall, leaving a preferred galaxy separation r_d ~ 147 Mpc
  -- the BAO ruler.
* **z*, photon decoupling** (~1090): when the universe becomes transparent.
  The CMB's acoustic peaks show the sound horizon at that moment, r_*, at the
  angle theta_* = r_* / D_M(z*) ~ 0.0104 rad (0.6 degrees).

Finding z_d and z* requires solving the atomic physics of hydrogen and helium
recombination, which is CAMB's job (RECFAST, with the BBN helium fraction).
We take only those two redshifts from CAMB, tabulated on a grid of
(omega_b, omega_c) -- the only parameters they depend on -- and compute the
sound horizons and distances ourselves (``background.py``). ``test_physics.py``
checks the result against CAMB's own r_d, r_* and theta_*.

Dark energy is negligible at z > 600 for every model the data allow
(w0 + wa < -0.5 gives rho_de/rho_m < 1e-4 there), so the rulers depend on
omega_b and omega_c alone and can be tabulated too.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.interpolate import RectBivariateSpline

import background as bg
from config import BBN_PREDICTOR, DERIVED_DIR, M_NU_EV, N_EFF, T_CMB

GRID_PATH = DERIVED_DIR / "recombination_grid.npz"

#: Grid of physical densities. Planck allows omega_b = 0.0222 +- 0.00015 and
#: omega_c = 0.119 +- 0.0012; the grid is ~25 sigma wide in each so that even
#: CMB-free fits stay inside it.
GRID_OMBH2 = np.linspace(0.0180, 0.0270, 19)
GRID_OMCH2 = np.linspace(0.085, 0.160, 26)


def camb_params(ombh2, omch2, H0=67.36, w0=-1.0, wa=0.0):
    """CAMB set up exactly as in DESI's runs (PPF dark energy, one 0.06 eV
    neutrino, N_eff = 3.044, PArthENoPE BBN helium, tau and A_s irrelevant here)."""
    import camb

    kw = dict(H0=H0, ombh2=ombh2, omch2=omch2, mnu=M_NU_EV, nnu=N_EFF,
              num_massive_neutrinos=1, TCMB=T_CMB, tau=0.0544, As=2.08e-9, ns=0.9649,
              bbn_predictor=BBN_PREDICTOR)
    if w0 != -1.0 or wa != 0.0:
        kw.update(dark_energy_model="ppf", w=w0, wa=wa)
    return camb.set_params(**kw)


def camb_derived(ombh2, omch2, H0=67.36, w0=-1.0, wa=0.0) -> dict:
    """CAMB's recombination redshifts, rulers and theta_* for one universe."""
    import camb

    pars = camb_params(ombh2, omch2, H0, w0, wa)
    res = camb.get_background(pars)
    d = res.get_derived_params()
    return {"zstar": d["zstar"], "rstar": d["rstar"], "zdrag": d["zdrag"],
            "rdrag": d["rdrag"], "thetastar": d["thetastar"] / 100.0,
            "DMstar": d["DAstar"] * 1000.0, "YHe": pars.YHe, "omnuh2": pars.omnuh2,
            "_results": res}


def build_grid(*, force: bool = False, verbose: bool = True) -> dict:
    """Tabulate z*, z_d (CAMB) and our own r_*, r_d on the (omega_b, omega_c) grid."""
    if GRID_PATH.exists() and not force:
        with np.load(GRID_PATH) as z:
            return {k: z[k] for k in z.files}
    t0 = time.monotonic()
    shape = (GRID_OMBH2.size, GRID_OMCH2.size)
    out = {k: np.empty(shape) for k in ("zstar", "zdrag", "rstar_camb", "rdrag_camb")}
    for i, ob in enumerate(GRID_OMBH2):
        for j, oc in enumerate(GRID_OMCH2):
            d = camb_derived(ob, oc)
            out["zstar"][i, j], out["zdrag"][i, j] = d["zstar"], d["zdrag"]
            out["rstar_camb"][i, j], out["rdrag_camb"][i, j] = d["rstar"], d["rdrag"]
    # Our own sound horizons at CAMB's redshifts (LCDM; dark energy is negligible here).
    ob, oc = np.meshgrid(GRID_OMBH2, GRID_OMCH2, indexing="ij")
    cosmo = bg.Cosmo(h=0.6736, ombh2=ob.ravel(), omch2=oc.ravel(), w0=-1.0, wa=0.0)
    sh = bg.SoundHorizonGrid()
    out["rstar"] = sh.r_s(out["zstar"].ravel(), cosmo).reshape(shape)
    out["rdrag"] = sh.r_s(out["zdrag"].ravel(), cosmo).reshape(shape)
    out["ombh2"], out["omch2"] = GRID_OMBH2, GRID_OMCH2
    np.savez(GRID_PATH, **out)
    if verbose:
        print(f"  recombination grid: {shape[0] * shape[1]} CAMB runs in "
              f"{time.monotonic() - t0:.0f} s; our r_d / CAMB's: "
              f"{np.abs(out['rdrag'] / out['rdrag_camb'] - 1).max():.1e} max deviation")
    return out


class EarlyUniverse:
    """Smooth interpolation of z*, z_d, r_*, r_d over (omega_b, omega_c)."""

    def __init__(self, grid: dict | None = None):
        g = build_grid() if grid is None else grid
        self.ob_range = (float(g["ombh2"][0]), float(g["ombh2"][-1]))
        self.oc_range = (float(g["omch2"][0]), float(g["omch2"][-1]))
        self._splines = {k: RectBivariateSpline(g["ombh2"], g["omch2"], g[k], kx=3, ky=3)
                         for k in ("zstar", "zdrag", "rstar", "rdrag")}

    def inside(self, ombh2, omch2) -> np.ndarray:
        return ((ombh2 >= self.ob_range[0]) & (ombh2 <= self.ob_range[1])
                & (omch2 >= self.oc_range[0]) & (omch2 <= self.oc_range[1]))

    def __call__(self, key: str, ombh2, omch2) -> np.ndarray:
        ob = np.clip(np.atleast_1d(ombh2), *self.ob_range)
        oc = np.clip(np.atleast_1d(omch2), *self.oc_range)
        return self._splines[key].ev(ob, oc)

    def rd(self, cosmo: bg.Cosmo) -> np.ndarray:
        return self("rdrag", cosmo.ombh2, cosmo.omch2)

    def theta_star(self, cosmo: bg.Cosmo, grid: bg.DistanceGrid, dc_nodes=None,
                   integrand=None) -> np.ndarray:
        """theta_* = r_s(z*) / D_M(z*), the CMB's acoustic angle (radians)."""
        zs = self("zstar", cosmo.ombh2, cosmo.omch2)
        dm = grid.comoving_at(zs, cosmo, dc_nodes=dc_nodes, integrand=integrand)
        return self("rstar", cosmo.ombh2, cosmo.omch2) / dm
