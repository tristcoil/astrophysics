"""The physics: from an Einstein ring to a mass, and from a mass to the speeds of stars.

Everything here is a formula with a known answer, checked in ``test_physics.py``.

Units. Angles on the sky are in arcseconds unless a name says otherwise,
distances in Mpc, masses in solar masses, speeds in km/s.

Three ideas carry the whole project:

1. **The ring weighs the lens.** Light passing a mass M at impact parameter
   R is bent by 4GM/(c^2 R). For a lens and a source almost exactly aligned,
   the source's light reaches us around a ring of angular radius theta_E, the
   Einstein radius, and the projected mass inside that ring is

       M_E = pi * Sigma_cr * (theta_E * D_l)^2,
       Sigma_cr = c^2 D_s / (4 pi G D_l D_ls),

   independent of how the mass is distributed inside the ring.

2. **An isothermal galaxy has one speed.** If the total density falls as
   1/r^2 -- the "singular isothermal sphere", whose rotation curve is flat --
   then theta_E = 4 pi (sigma/c)^2 D_ls / D_s. Lensing therefore predicts a
   velocity dispersion sigma_SIE, which can be compared with the stars' own.
   Only the *ratio* of distances enters, so H0 drops out.

3. **Stars feel all the mass.** The spherical Jeans equation turns a mass
   profile and a light profile into the velocity dispersion a spectrograph
   would measure through the SDSS fibre. With the mass inside the ring fixed by
   lensing, the measured dispersion pins down how steeply the density falls.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy import constants as const
from astropy.cosmology import FlatLambdaCDM
from scipy import special

C_KMS = const.c.to("km/s").value
G_KPC = const.G.to("kpc km2 / (solMass s2)").value  # G in kpc (km/s)^2 / Msun
ARCSEC = np.pi / (180 * 3600)  # radians per arcsecond
MPC_KPC = 1.0e3


# --------------------------------------------------------------------------
# Distances
# --------------------------------------------------------------------------

@lru_cache(maxsize=None)
def cosmology(h0: float, omega_m: float) -> FlatLambdaCDM:
    """Flat Lambda-CDM without radiation, as in the SLACS papers."""
    return FlatLambdaCDM(H0=h0, Om0=omega_m)


def distances(z_l, z_s, h0: float, omega_m: float) -> dict:
    """Angular-diameter distances (Mpc) to the lens, to the source, and between them.

    In a flat universe D_ls = (D_M(z_s) - D_M(z_l)) / (1 + z_s), with D_M the
    transverse comoving distance; astropy implements exactly that."""
    cosmo = cosmology(h0, omega_m)
    z_l = np.asarray(z_l, dtype=float)
    z_s = np.asarray(z_s, dtype=float)
    d_l = cosmo.angular_diameter_distance(z_l).value
    d_s = cosmo.angular_diameter_distance(z_s).value
    d_ls = cosmo.angular_diameter_distance(z_l, z_s).value
    return {"d_l": d_l, "d_s": d_s, "d_ls": d_ls}


def kpc_per_arcsec(d_l_mpc):
    return np.asarray(d_l_mpc) * MPC_KPC * ARCSEC


# --------------------------------------------------------------------------
# From Einstein radius to mass and velocity dispersion
# --------------------------------------------------------------------------

def sigma_crit(d_l, d_s, d_ls):
    """Critical surface density in Msun/kpc^2: c^2 D_s / (4 pi G D_l D_ls)."""
    d_l, d_s, d_ls = (np.asarray(x) * MPC_KPC for x in (d_l, d_s, d_ls))
    return C_KMS**2 * d_s / (4 * np.pi * G_KPC * d_l * d_ls)


def einstein_mass(theta_e_arcsec, d_l, d_s, d_ls):
    """Projected mass (Msun) inside the Einstein radius: pi Sigma_cr R_E^2.

    For the elliptical SIE with intermediate-axis normalisation (Kormann et al.
    1994, used by SLACS) this is the mass inside the critical curve."""
    r_e = np.asarray(theta_e_arcsec) * kpc_per_arcsec(d_l)
    return np.pi * sigma_crit(d_l, d_s, d_ls) * r_e**2


def sigma_sie(theta_e_arcsec, d_s, d_ls):
    """Velocity-dispersion parameter of an isothermal lens: theta_E = 4 pi (sigma/c)^2 D_ls/D_s."""
    theta = np.asarray(theta_e_arcsec) * ARCSEC
    return C_KMS * np.sqrt(theta * np.asarray(d_s) / (4 * np.pi * np.asarray(d_ls)))


def theta_e_sis(sigma_kms, d_s, d_ls):
    """Einstein radius (arcsec) a singular isothermal sphere of dispersion sigma would make."""
    return 4 * np.pi * (np.asarray(sigma_kms) / C_KMS) ** 2 * np.asarray(d_ls) / np.asarray(d_s) / ARCSEC


def jorgensen_sigma(sigma_ap, aperture_radius, r_eff, fraction):
    """Correct a dispersion measured in a circular aperture to one inside R_e*fraction.

    Jorgensen, Franx & Kjaergaard (1995): sigma_ap / sigma_e8 = (R_ap / (R_e/8))^-0.04.
    ``fraction`` = 1/8 gives sigma_e8, 1/2 gives sigma_e2."""
    from config import JORGENSEN_SLOPE
    return np.asarray(sigma_ap) * (np.asarray(aperture_radius)
                                   / (np.asarray(r_eff) * fraction)) ** (-JORGENSEN_SLOPE)


# --------------------------------------------------------------------------
# The singular isothermal ellipsoid (SIE)
# --------------------------------------------------------------------------

def sie_convergence(x, y, b, q):
    """kappa = (b/2) / sqrt(q x^2 + y^2/q), major axis along x.

    In this "intermediate-axis" normalisation (Kormann et al. 1994) the mass
    inside any isodensity ellipse is independent of q, and the critical curve is
    the ellipse q x^2 + y^2/q = b^2 (semi-axes b/sqrt(q) and b*sqrt(q))."""
    return 0.5 * b / np.sqrt(q * x**2 + y**2 / q)


def sie_deflection(x, y, b, q):
    """Deflection (alpha_x, alpha_y) of an SIE centred at the origin, major axis along x.

    Kormann et al. (1994); Keeton & Kochanek (1998) with core s = 0:
        alpha_x = b sqrt(q)/q' * arctan (q' x / psi)
        alpha_y = b sqrt(q)/q' * arctanh(q' y / psi),   psi = sqrt(q^2 x^2 + y^2), q' = sqrt(1 - q^2).
    For q -> 1 it becomes the isothermal sphere, alpha = b * r_hat."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    q = float(min(q, 0.99999))
    qp = np.sqrt(1.0 - q * q)
    psi = np.sqrt(q * q * x * x + y * y)
    psi = np.where(psi == 0, 1e-12, psi)
    scale = b * np.sqrt(q) / qp
    ax = scale * np.arctan(qp * x / psi)
    ay = scale * np.arctanh(np.clip(qp * y / psi, -1 + 1e-15, 1 - 1e-15))
    return ax, ay


def rotate(x, y, angle):
    """Rotate coordinates by -angle (into a frame whose x-axis lies at ``angle``)."""
    c, s = np.cos(angle), np.sin(angle)
    return c * x + s * y, -s * x + c * y


def sie_deflection_sky(x, y, b, q, phi, x0=0.0, y0=0.0):
    """SIE deflection for a lens centred at (x0, y0) with its major axis at angle phi
    (radians, counter-clockwise from +x)."""
    xr, yr = rotate(np.asarray(x) - x0, np.asarray(y) - y0, phi)
    axr, ayr = sie_deflection(xr, yr, b, q)
    c, s = np.cos(phi), np.sin(phi)
    return c * axr - s * ayr, s * axr + c * ayr


def shear_deflection(x, y, gamma, phi_gamma, x0=0.0, y0=0.0):
    """External shear of strength gamma along angle phi_gamma."""
    dx, dy = np.asarray(x) - x0, np.asarray(y) - y0
    g1, g2 = gamma * np.cos(2 * phi_gamma), gamma * np.sin(2 * phi_gamma)
    return g1 * dx + g2 * dy, g2 * dx - g1 * dy


def sie_critical_curve(b, q, phi, x0=0.0, y0=0.0, n=720):
    """The tangential critical curve: for any isothermal lens det(A) = 1 - 2 kappa,
    so it is the kappa = 1/2 contour, the ellipse q x^2 + y^2/q = b^2."""
    t = np.linspace(0, 2 * np.pi, n)
    xr, yr = b / np.sqrt(q) * np.cos(t), b * np.sqrt(q) * np.sin(t)
    c, s = np.cos(phi), np.sin(phi)
    return x0 + c * xr - s * yr, y0 + s * xr + c * yr


def sie_caustics(b, q, phi, x0=0.0, y0=0.0, n=720):
    """Tangential caustic (critical curve mapped to the source plane) and the
    'cut' (the image of the central singularity), both in source-plane coordinates."""
    xc, yc = sie_critical_curve(b, q, phi, x0, y0, n)
    ax, ay = sie_deflection_sky(xc, yc, b, q, phi, x0, y0)
    caustic = (xc - ax, yc - ay)
    t = np.linspace(0, 2 * np.pi, n)
    eps = 1e-6 * b
    xs, ys = x0 + eps * np.cos(t), y0 + eps * np.sin(t)
    ax, ay = sie_deflection_sky(xs, ys, b, q, phi, x0, y0)
    cut = (xs - ax, ys - ay)
    return caustic, cut


# --------------------------------------------------------------------------
# Light profiles
# --------------------------------------------------------------------------

def sersic_b(n):
    """b_n such that half the light of a Sersic profile lies inside R_e."""
    return special.gammaincinv(2 * np.asarray(n), 0.5)


def sersic_light_fraction(r, r_eff, n=4.0):
    """Fraction of a Sersic profile's total light inside projected radius r
    (n = 4 is de Vaucouleurs)."""
    bn = sersic_b(n)
    return special.gammainc(2 * n, bn * (np.asarray(r) / r_eff) ** (1.0 / n))


def sersic_profile(r, r_eff, n):
    """Surface brightness relative to its value at R_e."""
    return np.exp(-sersic_b(n) * ((np.asarray(r) / r_eff) ** (1.0 / n) - 1.0))


def _hernquist_x(s):
    """Hernquist (1990) eq. 33: X(s) = arccosh(1/s)/sqrt(1-s^2) (s<1), arccos(1/s)/sqrt(s^2-1) (s>1)."""
    s = np.asarray(s, dtype=float)
    out = np.ones_like(s)
    lo = s < 1 - 1e-6
    hi = s > 1 + 1e-6
    out[lo] = np.arccosh(1 / s[lo]) / np.sqrt(1 - s[lo] ** 2)
    out[hi] = np.arccos(1 / s[hi]) / np.sqrt(s[hi] ** 2 - 1)
    mid = ~(lo | hi)
    out[mid] = 1.0 - 2.0 * (s[mid] - 1.0) / 3.0  # Taylor expansion about s = 1
    return out


def hernquist_density(r, a):
    """3D density for unit total mass: a / (2 pi r (r + a)^3)."""
    r = np.asarray(r, dtype=float)
    return a / (2 * np.pi * r * (r + a) ** 3)


def hernquist_mass3d(r, a):
    """Fraction of a Hernquist model's mass inside radius r: r^2 / (r + a)^2."""
    r = np.asarray(r, dtype=float)
    return r**2 / (r + a) ** 2


def hernquist_surface(R, a):
    """Projected surface density for unit total mass (Hernquist 1990, eq. 32)."""
    s = np.asarray(R, dtype=float) / a
    s = np.where(np.abs(s - 1) < 1e-4, 1 + 1e-4, s)
    return ((2 + s**2) * _hernquist_x(s) - 3) / (2 * np.pi * a**2 * (1 - s**2) ** 2)


def hernquist_mass2d(R, a):
    """Fraction of a Hernquist model's mass projected inside radius R (Hernquist 1990, eq. 37)."""
    s = np.asarray(R, dtype=float) / a
    s = np.where(np.abs(s - 1) < 1e-4, 1 + 1e-4, s)
    return s**2 * (_hernquist_x(s) - 1) / (1 - s**2)


# --------------------------------------------------------------------------
# Mass profiles
# --------------------------------------------------------------------------

def powerlaw_mass3d(r, gamma, theta_e, m_e=1.0):
    """3D mass inside r of rho ~ r^-gamma, normalised so that the *projected* mass
    inside theta_e is m_e:  M_3D(r) = m_e (2/sqrt(pi)) Gamma(gamma/2)/Gamma((gamma-1)/2) (r/theta_e)^(3-gamma)."""
    k = 2 / np.sqrt(np.pi) * special.gamma(gamma / 2) / special.gamma((gamma - 1) / 2)
    return m_e * k * (np.asarray(r) / theta_e) ** (3 - gamma)


def nfw_mass3d(r, r_s):
    """NFW 3D mass in units of 4 pi rho_s r_s^3: ln(1+x) - x/(1+x)."""
    x = np.asarray(r, dtype=float) / r_s
    return np.log1p(x) - x / (1 + x)


def nfw_mass2d(R, r_s):
    """NFW projected mass in units of 4 pi rho_s r_s^3: ln(x/2) + F(x)
    (Bartelmann 1996; Golse & Kneib 2002)."""
    x = np.asarray(R, dtype=float) / r_s
    f = np.ones_like(x)
    lo, hi = x < 1 - 1e-6, x > 1 + 1e-6
    f[lo] = np.arccosh(1 / x[lo]) / np.sqrt(1 - x[lo] ** 2)
    f[hi] = np.arccos(1 / x[hi]) / np.sqrt(x[hi] ** 2 - 1)
    return np.log(x / 2) + f


# --------------------------------------------------------------------------
# Spherical Jeans equation, seen through a fibre
# --------------------------------------------------------------------------

def fibre_weight(R, aperture_radius, seeing_fwhm):
    """Fraction of the light emitted at projected radius R that lands inside a
    circular fibre after Gaussian seeing: W(R) = 1 - Q_1(R/s, R_ap/s), with Q_1
    Marcum's Q-function, written here as an integral over the fibre."""
    s = seeing_fwhm / (2 * np.sqrt(2 * np.log(2)))
    R = np.atleast_1d(np.asarray(R, dtype=float))
    rr = np.linspace(0, aperture_radius, 400)[None, :]
    arg = R[:, None] * rr / s**2
    # exp(-(r^2 + R^2)/2s^2) I0(rR/s^2) = exp(-(r - R)^2/2s^2) * i0e(rR/s^2)
    integrand = rr / s**2 * np.exp(-((rr - R[:, None]) ** 2) / (2 * s**2)) * special.i0e(arg)
    return np.trapezoid(integrand, rr, axis=1)


class JeansModel:
    """Line-of-sight velocity dispersion inside a seeing-blurred circular fibre,
    for stars following a Hernquist light profile in an arbitrary spherical mass
    distribution, with constant orbital anisotropy beta.

    Radii are in arcsec. Masses enter only through v2(r) = G M(r) / r in
    (km/s)^2, so a mass model is a function r -> G M(r)/r.

        nu sigma_r^2 (r) = r^(-2 beta) int_r^inf r'^(2 beta) nu(r') G M(r') / r'^2 dr'
        I sigma_p^2 (R)  = 2 int_R^inf (1 - beta R^2/r^2) nu sigma_r^2 r dr / sqrt(r^2 - R^2)
        sigma_fibre^2    = int I sigma_p^2 W 2 pi R dR / int I W 2 pi R dR
    """

    def __init__(self, r_eff, aperture_radius, seeing_fwhm, beta=0.0, n_r=400, n_u=200,
                 tracer=None):
        from config import HERNQUIST_RE_OVER_A
        self.a = r_eff / HERNQUIST_RE_OVER_A
        self.beta = beta
        rmax = max(200.0 * self.a, 30 * aperture_radius)
        self.r = np.geomspace(1e-4 * self.a, rmax, n_r)  # 3D radius grid
        # The tracer: Hernquist by default; ``tracer`` = (density(r), surface(R))
        # replaces it (the tests use power laws with known answers).
        density, surface = tracer or (lambda r: hernquist_density(r, self.a),
                                      lambda R: hernquist_surface(R, self.a))
        self.nu = density(self.r)
        # Projected radii covering the fibre plus the seeing halo.
        self.R = np.geomspace(1e-3 * self.a, aperture_radius + 4 * seeing_fwhm, 120)
        self.W = fibre_weight(self.R, aperture_radius, seeing_fwhm)
        self.I = surface(self.R)
        # Projection by u = sqrt(r^2 - R^2): r dr / sqrt(r^2 - R^2) = du.
        u = np.concatenate([[0.0], np.geomspace(1e-4 * self.a, rmax, n_u)])
        self.u = u
        self.rr = np.sqrt(self.R[:, None] ** 2 + u[None, :] ** 2)

    def nu_sigma_r2(self, v2):
        """nu sigma_r^2 on the radial grid, for G M(r)/r = v2(r), by integrating inward."""
        r, nu, b = self.r, self.nu, self.beta
        integrand = r ** (2 * b) * nu * v2(r) / r  # r'^2b nu G M / r'^2 = r'^2b nu (GM/r') / r'
        # cumulative integral from r to the outer edge (tail beyond is negligible)
        seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)
        tail = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])
        return r ** (-2 * b) * tail

    def sigma_fibre(self, v2):
        """Luminosity-weighted dispersion (km/s) inside the fibre."""
        nsr2 = self.nu_sigma_r2(v2)
        log_r = np.log(self.r)
        f = np.interp(np.log(self.rr), log_r, nsr2, right=0.0)
        aniso = 1 - self.beta * (self.R[:, None] / self.rr) ** 2
        i_sp2 = 2 * np.trapezoid(aniso * f, self.u, axis=1)
        w = self.W * self.R  # 2 pi cancels
        num = np.trapezoid(i_sp2 * w, self.R)
        den = np.trapezoid(self.I * w, self.R)
        return np.sqrt(num / den)


def powerlaw_v2(gamma, theta_e, sigma_sis):
    """G M(r)/r in (km/s)^2 for rho ~ r^-gamma with the lensing mass inside theta_e.

    G M_E / R_E = pi sigma_SIS^2 (with sigma_SIS from the Einstein radius), so
    G M(r)/r = pi sigma_SIS^2 (theta_e / r) * M_3D(r)/M_E. For gamma = 2 this is
    2 sigma_SIS^2, the isothermal sphere's v_c^2."""
    def v2(r):
        return np.pi * sigma_sis**2 * (theta_e / r) * powerlaw_mass3d(r, gamma, theta_e)
    return v2


def two_component_v2(f_star, theta_e, sigma_sis, a_hernquist, r_s):
    """G M(r)/r for stars (Hernquist, scale a) plus an NFW halo (scale r_s), with a
    fraction f_star of the projected mass inside theta_e in stars and the rest dark.
    All radii in arcsec."""
    gm_e_over = np.pi * sigma_sis**2 * theta_e  # G M_E in (km/s)^2 arcsec
    m_star_total = f_star / hernquist_mass2d(theta_e, a_hernquist)
    m_dm_unit = (1 - f_star) / nfw_mass2d(theta_e, r_s)

    def v2(r):
        m = m_star_total * hernquist_mass3d(r, a_hernquist) + m_dm_unit * nfw_mass3d(r, r_s)
        return gm_e_over * m / r
    return v2
