"""Distances, peculiar-velocity estimators, and the Lambda-CDM velocity statistics.

Everything here is physics, not data handling, and is unit-tested in
``test_physics.py``.

Distances are comoving and in h^-1 Mpc unless a name says otherwise: in those
units the Hubble constant is exactly 100 km/s per h^-1 Mpc, so one h^-1 Mpc
corresponds to 100 km/s of expansion at low redshift.

**Peculiar velocities.** A galaxy's redshift mixes the expansion with its own
motion: cz ~ H0 d + v. Given an independent distance d, the "direct" estimate
v = cz - H0 d (Davis & Scrimgeour 2014 add the relativistic corrections) is
unbiased only if distance errors are Gaussian -- but they are Gaussian in the
*distance modulus*, i.e. in log d. A 20% error then scatters v much further
toward negative values (galaxies placed too far) than positive ones. Watkins
& Feldman (2015) instead use v = cz/(1+z) ln(cz/H0 d), which is linear in log d
and so keeps the errors Gaussian; to first order it equals the direct estimate.
"""

from __future__ import annotations

import functools

import numpy as np
from scipy import integrate, special

import config
from config import C_KMS

H100 = 100.0  # km/s per h^-1 Mpc


# --------------------------------------------------------------------------
# Background expansion
# --------------------------------------------------------------------------

def E(z, Om: float):
    """H(z)/H0 for a flat universe of matter and a cosmological constant."""
    z = np.asarray(z, float)
    return np.sqrt(Om * (1 + z) ** 3 + (1 - Om))


@functools.lru_cache(maxsize=8)
def _dc_table(Om: float, zmax: float = 0.5, n: int = 20001):
    z = np.linspace(0.0, zmax, n)
    integrand = 1.0 / E(z, Om)
    dc = integrate.cumulative_trapezoid(integrand, z, initial=0.0) * C_KMS / H100
    return z, dc


def comoving_distance(z, Om: float = config.PLANCK18.Om0):
    """Line-of-sight comoving distance in h^-1 Mpc (z < 0.5; negative z allowed nearby)."""
    zt, dct = _dc_table(Om)
    z = np.asarray(z, float)
    out = np.interp(np.abs(z), zt, dct)
    return np.where(z < 0, -out, out)


def z_from_comoving(r, Om: float = config.PLANCK18.Om0):
    """Inverse of ``comoving_distance``."""
    zt, dct = _dc_table(Om)
    r = np.asarray(r, float)
    out = np.interp(np.abs(r), dct, zt)
    return np.where(r < 0, -out, out)


def growth_rate(Om: float) -> float:
    """f = d ln D / d ln a today, for flat Lambda-CDM, from the growth equation.

    Solves D'' + (3/a + E'/E) D' - 1.5 Om D / (a^5 E^2) = 0 from deep in matter
    domination (D = a). The familiar approximation is f ~ Om^0.55.
    """
    def Ea(a):
        return np.sqrt(Om / a**3 + 1 - Om)

    def dEda(a):
        return -1.5 * Om / a**4 / Ea(a)

    def rhs(a, y):
        D, dD = y
        return [dD, -(3.0 / a + dEda(a) / Ea(a)) * dD + 1.5 * Om / (a**5 * Ea(a) ** 2) * D]

    a0 = 1e-3
    sol = integrate.solve_ivp(rhs, (a0, 1.0), [a0, 1.0], rtol=1e-10, atol=1e-13)
    D, dD = sol.y[:, -1]
    return float(dD / D)


# --------------------------------------------------------------------------
# Cosmicflows-4 conventions
# --------------------------------------------------------------------------

def cf4_curvature_factor(z, q0: float = config.CF4_Q0, j0: float = config.CF4_J0):
    """f(z) = 1 + (1 - q0) z / 2 - (1 - q0 - 3 q0^2 + j0) z^2 / 6.

    Multiplying cz by f gives, to second order in z, the luminosity distance
    times H0 in an expanding universe (Tully et al. 2023, their Eq. 3; Davis &
    Scrimgeour 2014). Without it the Hubble constant would appear to drift with
    distance.
    """
    z = np.asarray(z, float)
    return 1.0 + 0.5 * (1.0 - q0) * z - (1.0 / 6.0) * (1.0 - q0 - 3.0 * q0**2 + j0) * z**2


def vpec_direct(cz, d_mpc, H0: float, *, f=None):
    """Davis & Scrimgeour (2014), as CF4 Eq. 9: (f cz - H0 d) / (1 + H0 d / c)."""
    cz = np.asarray(cz, float)
    f = cf4_curvature_factor(cz / C_KMS) if f is None else f
    hd = H0 * np.asarray(d_mpc, float)
    return (f * cz - hd) / (1.0 + hd / C_KMS)


def vpec_wf(cz, d_mpc, H0: float, *, f=None, log=np.log):
    """Watkins & Feldman (2015): cz_mod / (1 + cz_mod / c) * ln(cz_mod / H0 d).

    ``log`` is the natural logarithm, as Watkins & Feldman's derivation requires
    (-cz ln(1 - v/cz) ~ v). Passing ``np.log10`` reproduces CF4's published
    Table 4 column, which is smaller by a factor ln 10 = 2.303.
    """
    cz = np.asarray(cz, float)
    f = cf4_curvature_factor(cz / C_KMS) if f is None else f
    czm = f * cz
    with np.errstate(divide="ignore", invalid="ignore"):
        return czm / (1.0 + czm / C_KMS) * log(czm / (H0 * np.asarray(d_mpc, float)))


def vpec_cf4_ramp(cz_cmb, v_ls, d_mpc, H0: float, *, log=np.log):
    """CF4's adopted peculiar velocity (their Eq. 11): the direct estimate nearby,
    Watkins & Feldman's beyond, blended linearly up to a Local Sheet velocity of 3,000 km/s."""
    w = np.clip(np.asarray(v_ls, float) / config.CF4_RAMP_KMS, 0.0, 1.0)
    return vpec_direct(cz_cmb, d_mpc, H0) * (1 - w) + vpec_wf(cz_cmb, d_mpc, H0, log=log) * w


# --------------------------------------------------------------------------
# Nusser & Tully (2026): the linear log-distance velocity
# --------------------------------------------------------------------------

def luminosity_distance_h(z_cmb, z_helio, Om: float = config.PLANCK18.Om0):
    """D_L = (1 + z_helio) D_C(z_cmb), in h^-1 Mpc (Nusser & Tully 2026, Eq. A5)."""
    return (1.0 + np.asarray(z_helio, float)) * comoving_distance(z_cmb, Om)


def dm_from_redshift(z_cmb, z_helio, H0: float, Om: float = config.PLANCK18.Om0):
    """Distance modulus a galaxy would have if it had no peculiar velocity."""
    dl_mpc = luminosity_distance_h(z_cmb, z_helio, Om) / (H0 / 100.0)
    return 5.0 * np.log10(dl_mpc) + 25.0


def nt_K(z_cmb, z_helio, Om: float = config.PLANCK18.Om0):
    """The factor K_i that turns a distance-modulus offset into a radial velocity.

    u = -K a dmu with a = ln(10)/5 (Nusser & Tully 2026, Eqs. 6 and A6):
    K = A_U D_L,z with A_U = H(z) / ((1 + z_helio) B), B = 1 - H(z) D_C / (c (1 + z)).
    K tends to cz at low redshift.
    """
    z = np.asarray(z_cmb, float)
    H = H100 * E(z, Om)  # km/s per h^-1 Mpc
    dc = comoving_distance(z, Om)
    B = 1.0 - H * dc / (C_KMS * (1.0 + z))
    A = H / ((1.0 + np.asarray(z_helio, float)) * B)
    return A * luminosity_distance_h(z, z_helio, Om)


A_LN10_5 = np.log(10.0) / 5.0


# --------------------------------------------------------------------------
# Linear power spectrum and velocity correlations
# --------------------------------------------------------------------------

@functools.lru_cache(maxsize=2)
def linear_power(cosmo: config.Cosmology = config.PLANCK18):
    """z = 0 linear matter power spectrum from CAMB, normalised to sigma8.

    Returns k [h/Mpc] and P(k) [(h^-1 Mpc)^3] on a log grid from 1e-4 to 10 h/Mpc.
    """
    import camb

    pars = camb.set_params(H0=cosmo.H0, ombh2=cosmo.Ob0h2, omch2=cosmo.Oc0h2, mnu=cosmo.mnu, omk=0.0,
                           tau=cosmo.tau, As=2.1e-9, ns=cosmo.ns, WantTransfer=True, kmax=20.0)
    pars.set_matter_power(redshifts=[0.0], kmax=20.0)
    results = camb.get_results(pars)
    kh, _, pk = results.get_matter_power_spectrum(minkh=config.PK_KMIN, maxkh=config.PK_KMAX,
                                                  npoints=config.PK_NK)
    pk = pk[0]
    s8 = sigma_r(kh, pk, 8.0)
    pk = pk * (cosmo.sigma8 / s8) ** 2
    return kh, pk


def tophat_window(x):
    """Fourier transform of a spherical top hat: 3 (sin x - x cos x) / x^3."""
    x = np.asarray(x, float)
    out = np.ones_like(x)
    big = x > 1e-3
    xb = x[big]
    out[big] = 3.0 * (np.sin(xb) - xb * np.cos(xb)) / xb**3
    out[~big] = 1.0 - x[~big] ** 2 / 10.0
    return out


def sigma_r(k, pk, R: float) -> float:
    """rms density contrast in top-hat spheres of radius R (h^-1 Mpc)."""
    integrand = pk * tophat_window(k * R) ** 2 * k**3 / (2 * np.pi**2)
    return float(np.sqrt(integrate.trapezoid(integrand, np.log(k))))


def bulk_flow_rms(R, window: str = "tophat", cosmo: config.Cosmology = config.PLANCK18):
    """Expected rms 3D bulk flow (km/s) of an ideal, perfectly sampled sphere of radius R.

    <V^2> = (H0 f)^2 / (2 pi^2) int P(k) W(kR)^2 dk, with H0 = 100 in h^-1 Mpc units.
    ``window`` is "tophat" (a sphere of radius R) or "gaussian" (weight exp(-r^2/2R^2)).
    """
    k, pk = linear_power(cosmo)
    f = growth_rate(cosmo.Om0)
    R = np.atleast_1d(np.asarray(R, float))
    out = []
    for r in R:
        if window == "tophat":
            W = tophat_window(k * r)
        elif window == "gaussian":
            W = np.exp(-0.5 * (k * r) ** 2)
        else:
            raise ValueError(window)
        out.append(np.sqrt((H100 * f) ** 2 / (2 * np.pi**2) * integrate.trapezoid(pk * W**2 * k, np.log(k))))
    return np.array(out)


def _psi_kgrid():
    """A k grid fine enough for the oscillating Bessel kernels out to 700 h^-1 Mpc:
    logarithmic below 0.05 h/Mpc, then uniform with dk = 2e-4 (d(ks) = 0.14 at s = 700)."""
    k1 = np.geomspace(config.PK_KMIN, 0.05, 3000, endpoint=False)
    k2 = np.arange(0.05, config.PK_KMAX, 2e-4)
    return np.concatenate([k1, k2])


@functools.lru_cache(maxsize=2)
def velocity_correlations(cosmo: config.Cosmology = config.PLANCK18):
    """Linear-theory velocity correlation functions Psi_par(s), Psi_perp(s) in (km/s)^2.

    Psi_perp(s) = (H0 f)^2/(2 pi^2) int P(k) j1(ks)/(ks) dk
    Psi_par(s)  = (H0 f)^2/(2 pi^2) int P(k) [j0(ks) - 2 j1(ks)/(ks)] dk
    (Gorski 1989; Nusser & Tully 2026, Appendix B). Cached on disk because the
    oscillating integrals take ~10 s.
    """
    cache = config.DERIVED_DIR / f"psi_{cosmo.name.replace(' ', '_')}.npz"
    s = np.linspace(0.0, config.PSI_SMAX, config.PSI_NS)
    if cache.exists():
        d = np.load(cache)
        if d["s"].shape == s.shape and np.allclose(d["s"], s):
            return d["s"], d["par"], d["perp"]
    kl, pkl = linear_power(cosmo)
    k = _psi_kgrid()
    pk = np.exp(np.interp(np.log(k), np.log(kl), np.log(pkl)))
    f = growth_rate(cosmo.Om0)
    norm = (H100 * f) ** 2 / (2 * np.pi**2)
    wk = np.gradient(k)  # trapezoid-like weights on the non-uniform grid
    par = np.empty_like(s)
    perp = np.empty_like(s)
    for i0 in range(0, len(s), 200):
        ss = s[i0:i0 + 200, None]
        x = k[None, :] * ss
        with np.errstate(invalid="ignore", divide="ignore"):
            j1x = np.where(x > 1e-4, special.spherical_jn(1, x) / x, 1.0 / 3.0 - x**2 / 30.0)
        j0 = special.spherical_jn(0, x)
        perp[i0:i0 + 200] = norm * np.sum(pk * j1x * wk, axis=1)
        par[i0:i0 + 200] = norm * np.sum(pk * (j0 - 2.0 * j1x) * wk, axis=1)
    np.savez(cache, s=s, par=par, perp=perp)
    return s, par, perp


def sigma_v_1d(cosmo: config.Cosmology = config.PLANCK18) -> float:
    """One-dimensional rms peculiar velocity of a point in linear theory."""
    s, par, perp = velocity_correlations(cosmo)
    return float(np.sqrt(par[0]))
