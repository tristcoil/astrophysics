"""The expansion history of a real universe -- with curvature -- and the distances
that supernovae, BAO and the CMB measure.

``friedmann.py`` treats a universe as four numbers (Or, Om, Ok, OL). To compare
with data we need the full physical model that DESI and Planck use, so a
universe here is described by *physical* densities, omega = Omega h^2:

* ``h``       -- H0 / (100 km/s/Mpc)
* ``ombh2``   -- ordinary matter (baryons), omega_b
* ``omch2``   -- cold dark matter, omega_c
* ``omk``     -- the curvature parameter Omega_k (positive: open; negative: closed)

and everything else is fixed physics: the CMB photons (2.7255 K), the neutrino
background (N_eff = 3.044, one species with 0.06 eV of mass, whose density is
the exact Fermi-Dirac integral: radiation-like early, matter-like late), and a
cosmological constant taking whatever is left of the critical density:
Omega_Lambda = 1 - Omega_k - Omega_m - Omega_r (CAMB's convention, with the
massive neutrino counted in Omega_m).

Distances (comoving, Mpc). With P(a) = a^4 (H / 100 km/s/Mpc)^2,

    D_C(z) = (c / 100) int_a^1 da' / sqrt(P(a'))      line of sight
    D_M(z) = D_C                         flat space
           = R sinh(D_C / R)             open space (Omega_k > 0)
           = R sin(D_C / R)              closed space (Omega_k < 0),
                                         R = c / (H0 sqrt|Omega_k|), the curvature radius
    D_H(z) = c / H(z),   D_V = (z D_M^2 D_H)^(1/3)

Curvature bends light rays: in closed space they converge again, so a ruler of
fixed size looks *bigger* than in flat space at the same distance; in open
space, smaller. That is how the CMB's acoustic angle measures curvature.

The sound horizon, the BAO and CMB ruler, is the distance a sound wave in the
photon-baryon plasma travelled before the photons let go:

    r_s(z) = (c / 100) int_0^a c_s/c da' / sqrt(P(a')),   c_s = c / sqrt(3 (1 + R_b)),
    R_b = 3 rho_b / (4 rho_gamma) (baryons slow the sound down).

Integrals use Gauss-Legendre quadrature: in u = sqrt(a) for D_C (smooth even
deep in the matter era) and in a for r_s. The redshifts of photon decoupling
(z*) and the drag epoch (z_d) come from CAMB's recombination code, tabulated
on a grid of (omega_b, omega_c); everything else is computed here and checked
against CAMB, curvature included, in ``test_physics.py``.
"""

from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline, RectBivariateSpline

import config as cfg
from friedmann import OMEGA_GAMMA

C_H100 = cfg.C_KM_S / 100.0  # c / (100 km/s/Mpc) in Mpc

# --- Radiation and neutrinos (as in D5, CAMB's conventions) ---------------------

NU_PER_SPECIES = 7.0 / 8.0 * (4.0 / 11.0) ** (4.0 / 3.0)
#: CAMB shares N_eff equally: the massive eigenstate carries N_eff / 3 species.
G_MASSIVE = cfg.N_EFF / 3.0
N_MASSLESS = cfg.N_EFF - G_MASSIVE
#: Photons + the two massless neutrino species (always radiation).
OMEGA_R = OMEGA_GAMMA * (1.0 + N_MASSLESS * NU_PER_SPECIES)
#: Today's density of the 0.06 eV neutrino as CAMB 2.0.4 computes it.
OMEGA_NU = 6.448993754203254e-4 if cfg.M_NU_EV == 0.06 else \
    cfg.M_NU_EV / 94.06410581217612 * G_MASSIVE**0.75
#: DESI/Cobaya's mapping from sampled Omega_m to omega_c.
OMEGA_NU_MAPPING = cfg.M_NU_EV / cfg.NU_MASS_TO_OMEGA

_FD_F0 = 7.0 * np.pi**4 / 120.0
_FD_NR = 1.5 * 1.2020569031595942
_GL_X, _GL_W = np.polynomial.laguerre.laggauss(80)


def _fd_energy_ratio(y) -> np.ndarray:
    """rho / rho_relativistic for a thermal fermion with m/kT = y."""
    y = np.atleast_1d(np.asarray(y, float))
    x = _GL_X[None, :]
    g = x**2 * np.sqrt(x**2 + y[:, None] ** 2) / (1.0 + np.exp(-x))
    return (g @ _GL_W) / _FD_F0


#: m / kT of the massive neutrino today, fixed by its density today.
Y_NU0 = OMEGA_NU / (G_MASSIVE * NU_PER_SPECIES * OMEGA_GAMMA * _FD_NR / _FD_F0)
_NU_REL = G_MASSIVE * NU_PER_SPECIES * OMEGA_GAMMA

_LNA_TAB = np.linspace(np.log(1e-10), 0.0, 2401)
_NU_TAB = CubicSpline(_LNA_TAB, _NU_REL * _fd_energy_ratio(Y_NU0 * np.exp(_LNA_TAB)))


def nu_a4(a) -> np.ndarray:
    """a^4 omega_nu(a) for the massive neutrino: constant while relativistic,
    growing as a once it is slow (around z ~ 100 for 0.06 eV)."""
    a = np.asarray(a, float)
    return _NU_TAB(np.log(np.clip(a, 1e-10, 1.0)))


# --- A batch of universes ---------------------------------------------------------


@dataclass
class Model:
    """Parameters of n universes (1-D arrays of equal length; scalars broadcast)."""

    h: np.ndarray
    ombh2: np.ndarray
    omch2: np.ndarray
    omk: np.ndarray

    def __post_init__(self):
        arrs = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float)) for v in
                                     (self.h, self.ombh2, self.omch2, self.omk)))
        self.h, self.ombh2, self.omch2, self.omk = (np.array(a, float) for a in arrs)

    @property
    def n(self) -> int:
        return self.h.size

    @property
    def omega_cb(self):
        return self.ombh2 + self.omch2

    @property
    def omega_k(self):
        return self.omk * self.h**2

    @property
    def omega_de(self):
        return self.h**2 - self.omega_cb - OMEGA_NU - OMEGA_R - self.omega_k

    @property
    def Om(self):
        """Omega_m including the massive neutrino (CAMB's and DESI's definition)."""
        return (self.omega_cb + OMEGA_NU) / self.h**2

    @property
    def OL(self):
        return self.omega_de / self.h**2

    @property
    def Or(self):
        return OMEGA_R / self.h**2

    @classmethod
    def from_plane(cls, Om, OL, h, ombh2) -> "Model":
        """From (Omega_m, Omega_Lambda) at given h and omega_b, mapping Omega_m to
        omega_c exactly as DESI's chains do; Omega_k follows from the closure."""
        Om, OL, h, ombh2 = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float))
                                                 for v in (Om, OL, h, ombh2)))
        omch2 = Om * h**2 - OMEGA_NU_MAPPING - ombh2
        omega_m = omch2 + ombh2 + OMEGA_NU
        omk = 1.0 - OL - omega_m / h**2 - OMEGA_R / h**2
        return cls(h=h, ombh2=ombh2, omch2=omch2, omk=omk)

    @classmethod
    def from_omk(cls, Om, Ok, h, ombh2) -> "Model":
        """From (Omega_m, Omega_k) at given h and omega_b (DESI's mapping to omega_c)."""
        Om, Ok, h, ombh2 = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float))
                                                 for v in (Om, Ok, h, ombh2)))
        return cls(h=h, ombh2=ombh2, omch2=Om * h**2 - OMEGA_NU_MAPPING - ombh2, omk=Ok)

    @classmethod
    def from_physical(cls, Om, Ok, ombh2, omega_bc) -> "Model":
        """From (Omega_m, Omega_k) and the physical densities the CMB measures:
        h follows from Omega_m h^2 = omega_b + omega_c + omega_nu (DESI's mapping)."""
        Om, Ok, ombh2, omega_bc = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float))
                                                        for v in (Om, Ok, ombh2, omega_bc)))
        h = np.sqrt((omega_bc + OMEGA_NU_MAPPING) / Om)
        return cls.from_omk(Om, Ok, h, ombh2)

    def take(self, idx) -> "Model":
        return Model(self.h[idx], self.ombh2[idx], self.omch2[idx], self.omk[idx])


def p_phys(a, m: Model) -> np.ndarray:
    """a^4 (H/100)^2 for each model (rows) at each a (columns or trailing axes)."""
    a = np.asarray(a, float)
    shape = (m.n,) + (1,) * a.ndim
    r = lambda v: v.reshape(shape)  # noqa: E731
    return (OMEGA_R + r(m.omega_cb) * a + nu_a4(a) + r(m.omega_k) * a**2
            + r(m.omega_de) * a**4)


def hubble(z, m: Model) -> np.ndarray:
    """H(z) in km/s/Mpc, shape (n_models, n_z)."""
    a = 1.0 / (1.0 + np.atleast_1d(np.asarray(z, float)))
    return 100.0 * np.sqrt(p_phys(a, m)) / a**2


_LEG = {k: np.polynomial.legendre.leggauss(k) for k in (32, 48, 64, 96, 128)}


def comoving_distance(z, m: Model, nodes: int = 96) -> np.ndarray:
    """D_C(z) in Mpc (Gauss-Legendre in u = sqrt(a)); shape (n_models, n_z)."""
    z = np.atleast_1d(np.asarray(z, float))
    x, w = _LEG[nodes]
    u_lo = np.sqrt(1.0 / (1.0 + z))
    half = 0.5 * (1.0 - u_lo)
    u = u_lo[:, None] + half[:, None] * (x[None, :] + 1.0)        # (n_z, k)
    f = 2.0 * u[None] / np.sqrt(p_phys(u * u, m))                   # (n_m, n_z, k)
    return C_H100 * half[None, :] * (f @ w)


def comoving_distance_each(z_each, m: Model, nodes: int = 96) -> np.ndarray:
    """D_C to one (model-dependent) redshift per model, e.g. z*: shape (n_models,)."""
    z_each = np.broadcast_to(np.asarray(z_each, float), (m.n,))
    x, w = _LEG[nodes]
    u_lo = np.sqrt(1.0 / (1.0 + z_each))
    half = 0.5 * (1.0 - u_lo)
    u = u_lo[:, None] + half[:, None] * (x[None, :] + 1.0)        # (n_m, k)
    a = u * u
    p = (OMEGA_R + m.omega_cb[:, None] * a + nu_a4(a) + m.omega_k[:, None] * a**2
         + m.omega_de[:, None] * a**4)
    return C_H100 * half * ((2.0 * u / np.sqrt(p)) @ w)


def transverse_distance(dc, m: Model) -> np.ndarray:
    """D_M from D_C: rays spread as sinh in open space and sin in closed space.

    In closed space D_M turns over once light has come more than a quarter of
    the way around the universe; past the antipode no image could form at all
    (D_M <= 0), which is flagged NaN.
    """
    dc = np.asarray(dc, float)
    ok = m.omk.reshape((m.n,) + (1,) * (dc.ndim - 1))
    dh0 = (C_H100 / m.h).reshape(ok.shape)
    sk = np.sqrt(np.abs(ok))
    with np.errstate(invalid="ignore", divide="ignore"):
        x = sk * dc / dh0
        open_ = dh0 * np.sinh(x) / np.where(sk > 0, sk, 1.0)
        closed = dh0 * np.sin(x) / np.where(sk > 0, sk, 1.0)
        # second-order series where |Omega_k| is tiny (avoids 0/0)
        series = dc * (1.0 + ok * (dc / dh0) ** 2 / 6.0)
        dm = np.where(np.abs(ok) < 1e-9, series, np.where(ok > 0, open_, closed))
    return np.where(dm > 0, dm, np.nan)


def sound_horizon_each(z_each, m: Model, nodes: int = 64) -> np.ndarray:
    """r_s(z) in Mpc, one redshift per model: shape (n_models,)."""
    z_each = np.broadcast_to(np.asarray(z_each, float), (m.n,))
    x, w = _LEG[nodes]
    a_hi = 1.0 / (1.0 + z_each)
    a = 0.5 * a_hi[:, None] * (x[None, :] + 1.0)                    # (n_m, k)
    p = (OMEGA_R + m.omega_cb[:, None] * a + nu_a4(a) + m.omega_k[:, None] * a**2
         + m.omega_de[:, None] * a**4)
    Rb = 0.75 * (m.ombh2[:, None] / OMEGA_GAMMA) * a
    f = 1.0 / (np.sqrt(3.0 * (1.0 + Rb)) * np.sqrt(p))
    return C_H100 * 0.5 * a_hi * (f @ w)


def age_gyr(m: Model, nodes: int = 128) -> np.ndarray:
    """Time since the Big Bang, in Gyr: t0 = int_0^1 a da / (100 sqrt(P)), in u = sqrt(a).

    (Valid for universes with a Big Bang, which is every model the data allow.)
    """
    x, w = _LEG[nodes]
    u = 0.5 * (x + 1.0)
    p = p_phys(u * u, m)
    return cfg.HUBBLE_TIME_GYR_X_H0 / 100.0 * 0.5 * ((2.0 * u**3 / np.sqrt(p)) @ w)


def bao_predictions(kinds, z, dm, hz, rd) -> np.ndarray:
    """D_M/r_d, D_H/r_d or D_V/r_d for each entry of a BAO data vector."""
    dh = cfg.C_KM_S / hz
    dv = np.cbrt(np.asarray(z)[None, :] * dm**2 * dh)
    out = np.empty_like(dm)
    for j, kind in enumerate(kinds):
        out[:, j] = {"DM_over_rs": dm[:, j], "DH_over_rs": dh[:, j],
                     "DV_over_rs": dv[:, j]}[kind]
    return out / np.asarray(rd)[:, None]


# --- Supernovae: many redshifts at once ---------------------------------------------


class SupernovaDistances:
    """D_M at every supernova redshift for many models, by a cumulative trapezoid
    in x = ln(1+z) on nodes that include each supernova redshift exactly (step
    <= 0.002, accurate to ~1e-6; the method of projects D1 and D5)."""

    def __init__(self, z, step: float = 0.002):
        z = np.asarray(z, float)
        x_top = float(np.log1p(z.max()))
        fill = np.linspace(0.0, x_top, int(np.ceil(x_top / step)) + 1)
        x_nodes, inv = np.unique(np.concatenate([fill, np.log1p(z)]), return_inverse=True)
        self.x = x_nodes
        self.a = np.exp(-x_nodes)
        self.dx = np.diff(x_nodes)
        self.index = inv[fill.size:]
        self.nu = nu_a4(self.a)

    def comoving(self, m: Model) -> np.ndarray:
        a = self.a[None, :]
        p = (OMEGA_R + m.omega_cb[:, None] * a + self.nu[None, :] + m.omega_k[:, None] * a**2
             + m.omega_de[:, None] * a**4)
        bad = np.any(p <= 0, axis=1)
        # dD_C/dx = (c/100) (1+z) / sqrt(P / a^4) = (c/100) a / sqrt(P)
        f = C_H100 * a / np.sqrt(np.where(p > 0, p, np.nan))
        steps = 0.5 * (f[:, 1:] + f[:, :-1]) * self.dx[None, :]
        dc = np.concatenate([np.zeros((m.n, 1)), np.cumsum(steps, axis=1)], axis=1)
        dc[bad] = np.nan
        return dc[:, self.index]

    def transverse(self, m: Model) -> np.ndarray:
        return transverse_distance(self.comoving(m), m)


# --- Recombination: CAMB's z* and z_d on a grid of physical densities -----------------

GRID_PATH = cfg.DERIVED_DIR / "recombination_grid.npz"
GRID_OMBH2 = np.linspace(0.0180, 0.0270, 19)
GRID_OMCH2 = np.linspace(0.085, 0.160, 26)


def camb_params(ombh2, omch2, H0=67.36, omk=0.0):
    """CAMB set up as in DESI's runs (one 0.06 eV neutrino, N_eff = 3.044,
    PArthENoPE BBN helium); tau, A_s and n_s do not affect the background."""
    import camb

    return camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, omk=omk, mnu=cfg.M_NU_EV,
                           nnu=cfg.N_EFF, num_massive_neutrinos=1, TCMB=cfg.T_CMB,
                           tau=0.0544, As=2.08e-9, ns=0.9649,
                           bbn_predictor=cfg.BBN_PREDICTOR)


def build_recombination_grid(*, force: bool = False, verbose: bool = True) -> dict:
    """z* and z_d from CAMB over (omega_b, omega_c); they depend on nothing else
    (the expansion rate around z ~ 1100 is set by the physical densities alone)."""
    if GRID_PATH.exists() and not force:
        with np.load(GRID_PATH) as z:
            return {k: z[k] for k in z.files}
    import camb

    t0 = time.monotonic()
    shape = (GRID_OMBH2.size, GRID_OMCH2.size)
    zstar, zdrag = np.empty(shape), np.empty(shape)
    for i, ob in enumerate(GRID_OMBH2):
        for j, oc in enumerate(GRID_OMCH2):
            d = camb.get_background(camb_params(ob, oc)).get_derived_params()
            zstar[i, j], zdrag[i, j] = d["zstar"], d["zdrag"]
    out = {"ombh2": GRID_OMBH2, "omch2": GRID_OMCH2, "zstar": zstar, "zdrag": zdrag}
    np.savez(GRID_PATH, **out)
    if verbose:
        print(f"  recombination grid: {zstar.size} CAMB runs in {time.monotonic() - t0:.0f} s")
    return out


class Recombination:
    """Smooth interpolation of CAMB's z* and z_d over (omega_b, omega_c)."""

    def __init__(self, grid: dict | None = None):
        g = build_recombination_grid() if grid is None else grid
        self.ob_range = (float(g["ombh2"][0]), float(g["ombh2"][-1]))
        self.oc_range = (float(g["omch2"][0]), float(g["omch2"][-1]))
        self._spl = {k: RectBivariateSpline(g["ombh2"], g["omch2"], g[k], kx=3, ky=3)
                     for k in ("zstar", "zdrag")}

    def __call__(self, key: str, ombh2, omch2) -> np.ndarray:
        return self._spl[key].ev(np.clip(ombh2, *self.ob_range), np.clip(omch2, *self.oc_range))

    def theta_star(self, m: Model) -> np.ndarray:
        """The CMB's acoustic angle theta_* = r_s(z*) / D_M(z*), in radians."""
        zs = self("zstar", m.ombh2, m.omch2)
        dm = transverse_distance(comoving_distance_each(zs, m), m)
        return sound_horizon_each(zs, m) / dm

    def r_drag(self, m: Model) -> np.ndarray:
        """The BAO ruler: the sound horizon at the drag epoch, in Mpc."""
        return sound_horizon_each(self("zdrag", m.ombh2, m.omch2), m)
