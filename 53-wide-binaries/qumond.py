"""The MOND force between two stars in the Galaxy's gravitational field.

In QUMOND (Milgrom 2010), the true gravitational field g follows from the
Newtonian one g_N by

    div g = div [ nu(|g_N| / a0) g_N ],

with the "simple" interpolating function nu(y) = 1/2 + sqrt(1/4 + 1/y)
(Famaey & Binney 2005), the choice of Banik & Zhao (2018) and Banik et al.
(2024). For a wide binary in the Solar neighbourhood, g_N is the binary's own
Newtonian field plus the Galaxy's Newtonian-equivalent field g_Ne: MOND is
non-linear, so the Galaxy's pull changes the binary's internal dynamics even
though it is uniform across the pair (the "external field effect", EFE).

Following Banik & Zhao (2018), the binary is a test particle orbiting a point
mass M (the total mass), and lengths are measured in the MOND radius
r_M = sqrt(GM/a0), accelerations in a0. Then the only parameter is g_Ne/a0.

The extra source term div[(nu - 1) g_N] acts like a "phantom" mass density.
Away from the star it has a closed form,

    rho_ph = -(1/4 pi G) nu'(y)/y * [3 (g_N . r^)^2 - |g_N|^2] / r^3,

(derived in the README), which falls as 1/r near the star and as a pure
quadrupole (3 cos^2 theta - 1)/r^3 far from it. We expand rho_ph in Legendre
polynomials of cos(theta) (theta measured from the external field) and solve
Poisson's equation multipole by multipole, which gives the field on a grid
without any 3D mesh. By Gauss's law the angle-averaged radial boost is simply
eta(r) = 1 + M_ph(<r)/M, and far out it must reach the analytic EFE value
nu_e (1 + K_e/3) (Banik & Zhao 2018, eq. 16): a built-in check.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.polynomial import legendre as npleg

from config import A0_SI, DERIVED_DIR


def nu_simple(y):
    """Simple interpolating function nu(y), y = g_N / a0."""
    y = np.asarray(y, float)
    return 0.5 + np.sqrt(0.25 + 1.0 / y)


def dnu_simple(y):
    y = np.asarray(y, float)
    return -0.5 / (y * y * np.sqrt(0.25 + 1.0 / y))


def g_from_gN(y):
    """Spherical-symmetry MOND acceleration (in a0) for Newtonian y."""
    return nu_simple(y) * y


def gNe_from_ge(ge: float) -> float:
    """Invert g_e = nu(g_Ne) g_Ne for the simple function (closed form):
    g_Ne = g_e^2 / (1 + g_e) in units of a0."""
    return ge * ge / (1.0 + ge)


def asymptotic_eta(yNe: float) -> float:
    """eta = nu_e (1 + K_e / 3), K_e = dln nu / dln y at the external field."""
    nu_e = float(nu_simple(yNe))
    K_e = float(yNe * dnu_simple(yNe) / nu_e)
    return nu_e * (1.0 + K_e / 3.0)


def galactic_field(vc_kms: float = 232.8, R0_kpc: float = 8.2) -> dict:
    """The Galaxy's pull on the Solar neighbourhood, g_e = v_c^2 / R0
    (Banik et al. 2024 adopt McMillan 2017's 232.8 km/s and 8.2 kpc)."""
    ge = (vc_kms * 1e3) ** 2 / (R0_kpc * 3.0856775814913673e19) / A0_SI
    yNe = gNe_from_ge(ge)
    return dict(ge_a0=ge, gNe_a0=yNe, nu_e=float(nu_simple(yNe)),
                K_e=float(yNe * dnu_simple(yNe) / nu_simple(yNe)), eta_inf=asymptotic_eta(yNe))


def phantom_density(r, mu, yNe: float):
    """rho_ph (units M / r_M^3 with G = 1) at radius r (r_M) and mu = cos(theta)."""
    r = np.asarray(r, float)[..., None] if np.ndim(r) == 1 and np.ndim(mu) == 1 else r
    gr_dot = yNe * mu - 1.0 / r**2                       # g_N . r^
    y2 = yNe**2 - 2.0 * yNe * mu / r**2 + 1.0 / r**4      # |g_N|^2
    y = np.sqrt(np.maximum(y2, 1e-300))
    bracket = (3.0 * gr_dot**2 - y2) / r**3
    return -(dnu_simple(y) / y) * bracket / (4.0 * np.pi)


@dataclass
class QumondField:
    """Tabulated QUMOND field of a unit point mass in an external field.

    ratio_r(r, mu): g_r / (-1/r^2), the radial field over the Newtonian
    point-mass field; t_theta(r, mu): r^2 g_theta (both include the star's own
    Newtonian field; the uniform external field itself is excluded, since the
    binary's barycentre falls freely in it). Grids: ln r and mu.
    """
    yNe: float
    lnr: np.ndarray
    mu: np.ndarray
    ratio_r: np.ndarray
    t_theta: np.ndarray
    u_theta: np.ndarray      # r^2 g_theta / sin(theta): finite on the axis
    eta: np.ndarray          # angle-averaged radial boost on lnr grid
    m_ph: np.ndarray         # enclosed phantom monopole mass (units of M)
    self_g: float            # phantom field at the star itself (a0), along +z
    lmax: int


def solve(yNe: float, *, lmax: int = 192, nmu: int = 512, rmin: float = 1e-5,
          rmax: float = 1e5, nr: int = 6000, lnr_out=None, mu_out=None) -> QumondField:
    """Solve for the QUMOND field multipole by multipole (see module docstring)."""
    lnr = np.linspace(np.log(rmin), np.log(rmax), nr)
    h = lnr[1] - lnr[0]
    r = np.exp(lnr)
    x, w = npleg.leggauss(nmu)
    rho = phantom_density(r[:, None], x[None, :], yNe)            # (nr, nmu)
    ls = np.arange(lmax + 1)
    P = np.array([npleg.legval(x, np.eye(lmax + 1)[l]) for l in ls])  # (L, nmu)
    rho_l = (rho * w[None, :]) @ P.T * ((2 * ls + 1) / 2.0)[None, :]  # (nr, L)

    q = np.exp(-h)
    A = np.zeros_like(rho_l)
    B = np.zeros_like(rho_l)
    f = rho_l * (r**2)[:, None]
    # inner boundary: near the star rho_ph -> 1/(2 pi r), isotropic
    A[0, 0] = r[0] / (4.0 * np.pi)
    qa = q ** (ls + 1)
    qb = q ** ls
    for i in range(1, nr):
        A[i] = qa * A[i - 1] + 0.5 * h * (f[i - 1] * qa + f[i])
    # outer boundary: the far field is a pure quadrupole rho_2 ~ C / r^3
    C2 = rho_l[-1, 2] * r[-1] ** 3
    B[-1, 2] = C2 / (3.0 * r[-1])
    for i in range(nr - 2, -1, -1):
        B[i] = qb * B[i + 1] + 0.5 * h * (f[i] + f[i + 1] * qb)

    pref = 4.0 * np.pi / (2 * ls + 1)
    m_ph = 4.0 * np.pi * A[:, 0] * r          # monopole: 4 pi int rho_0 s^2 ds
    eta = 1.0 + m_ph
    # self field at the star: only l = 1 contributes, g_z(0) = (4 pi/3) int rho_1 ds
    self_g = float(4.0 * np.pi / 3.0 * np.trapezoid(rho_l[:, 1] * r, lnr))

    lnr_o = np.linspace(np.log(1e-3), np.log(1e3), 721) if lnr_out is None else lnr_out
    mu_o = np.cos(np.linspace(0, np.pi, 361)) if mu_out is None else mu_out
    idx = np.searchsorted(lnr, lnr_o)
    idx = np.clip(idx, 1, nr - 1)
    # linear interpolation of A, B in ln r onto the output radii
    t = ((lnr_o - lnr[idx - 1]) / h)[:, None]
    Ao = (1 - t) * A[idx - 1] + t * A[idx]
    Bo = (1 - t) * B[idx - 1] + t * B[idx]
    ro = np.exp(lnr_o)[:, None]
    dphi_dr = -pref[None, :] * (-(ls + 1)[None, :] * Ao + ls[None, :] * Bo) / ro
    phi_l = -pref[None, :] * (Ao + Bo)
    Po = np.array([npleg.legval(mu_o, np.eye(lmax + 1)[l]) for l in ls])          # (L, nmu_o)
    dPo = np.array([npleg.legval(mu_o, npleg.legder(np.eye(lmax + 1)[l])) for l in ls])
    g_r_ph = -(dphi_dr @ Po)                                   # (nr_o, nmu_o)
    sin_t = np.sqrt(np.clip(1 - mu_o**2, 0, 1))[None, :]
    g_t_ph = (phi_l @ dPo) * sin_t / ro                        # -(1/r) dPhi/dtheta
    ratio_r = 1.0 - g_r_ph * ro**2                             # (-1/r^2 + g_r_ph) / (-1/r^2)
    t_theta = g_t_ph * ro**2
    u_theta = (phi_l @ dPo) * ro
    eta_o = np.interp(lnr_o, lnr, eta)
    return QumondField(yNe=yNe, lnr=lnr_o, mu=mu_o, ratio_r=ratio_r, t_theta=t_theta,
                       u_theta=u_theta, eta=eta_o, m_ph=np.interp(lnr_o, lnr, m_ph), self_g=self_g, lmax=lmax)


def cached_field(yNe: float, lmax: int = 192) -> QumondField:
    """solve(), cached on disk (data/derived/qumond_<yNe>_<lmax>.npz)."""
    path = DERIVED_DIR / f"qumond_yNe{yNe:.5f}_l{lmax}.npz"
    if path.exists():
        z = np.load(path)
        return QumondField(yNe=float(z["yNe"]), lnr=z["lnr"], mu=z["mu"], ratio_r=z["ratio_r"],
                           t_theta=z["t_theta"], u_theta=z["u_theta"], eta=z["eta"], m_ph=z["m_ph"],
                           self_g=float(z["self_g"]), lmax=int(z["lmax"]))
    f = solve(yNe, lmax=lmax)
    np.savez_compressed(path, yNe=f.yNe, lnr=f.lnr, mu=f.mu, ratio_r=f.ratio_r,
                        t_theta=f.t_theta, u_theta=f.u_theta, eta=f.eta, m_ph=f.m_ph, self_g=f.self_g, lmax=f.lmax)
    return f


if __name__ == "__main__":
    import time
    gf = galactic_field()
    print({k: round(v, 4) for k, v in gf.items()})
    for lmax in (64, 128, 192):
        t0 = time.time()
        f = solve(gf["gNe_a0"], lmax=lmax)
        i = [np.argmin(abs(np.exp(f.lnr) - x)) for x in (0.1, 0.3, 1.0, 3.0, 10.0, 100.0)]
        print(lmax, f"{time.time() - t0:.1f}s eta:", np.round(f.eta[i], 4), "self_g", round(f.self_g, 4),
              "ratio_r(r=1, mu=1,0,-1):", np.round(f.ratio_r[i[2], [0, 180, 360]], 4),
              "t_theta(r=1, mu=0):", round(float(f.t_theta[i[2], 180]), 4))
