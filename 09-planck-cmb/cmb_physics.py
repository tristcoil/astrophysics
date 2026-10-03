"""CMB physics: the formulas behind every number this project produces.

Constants come from ``astropy.constants`` (CODATA 2018) rather than being typed
in, and the cosmology functions take ``astropy.cosmology`` objects so the
expansion history -- photons, massive neutrinos and all -- is computed by a
well-tested library rather than reimplemented here.

Four groups of formulas:

* **Blackbody radiation.** The CMB's spectrum, which FIRAS measured.
* **Angular power spectra.** How a map's fluctuations are summarised by C_l,
  and the conventions (D_l, cosmic variance, binning) needed to read them.
* **Beams and masks.** What a telescope and a cut sky do to a spectrum, and how
  to undo it.
* **The acoustic scale.** Why the peaks sit where they do: a sound horizon
  (set by physics before recombination) seen across the universe (set by its
  geometry and expansion).

A note on the multipole l. A spherical harmonic of degree l has l wavelengths
around a great circle, so structure at multipole l has an angular size of about
``180 deg / l`` -- l = 2 is the scale of the whole sky, l = 200 is about one
degree, l = 2000 about five arcminutes.
"""

from __future__ import annotations

import numpy as np
from astropy import constants as const
from astropy import units as u
from scipy import integrate

H = const.h.si.value
K_B = const.k_B.si.value
C = const.c.si.value
HBAR = const.hbar.si.value

#: 1 MJy/sr in SI units (W m^-2 Hz^-1 sr^-1).
MJY_PER_SR = 1e-20

# --------------------------------------------------------------------------
# Blackbody radiation
# --------------------------------------------------------------------------


def wavenumber_to_ghz(wavenumber_per_cm):
    """Convert spectroscopists' cm^-1 (what FIRAS reports) to GHz.

    A wavenumber is the number of wavelengths per centimetre, so the frequency
    is ``c`` times it: 1 cm^-1 = 29.9792458 GHz.
    """
    return np.asarray(wavenumber_per_cm, float) * C * 100.0 / 1e9


def planck_intensity(nu_ghz, temperature_k):
    """Blackbody specific intensity B_nu(T), in MJy/sr.

    ``B_nu(T) = (2 h nu^3 / c^2) / (exp(h nu / k T) - 1)``

    MJy/sr (megajansky per steradian) is the radio astronomer's unit of surface
    brightness; the CMB peaks near 385 MJy/sr.
    """
    nu = np.asarray(nu_ghz, float) * 1e9
    x = H * nu / (K_B * temperature_k)
    return 2.0 * H * nu**3 / C**2 / np.expm1(x) / MJY_PER_SR


def planck_intensity_dT(nu_ghz, temperature_k):
    """dB_nu/dT in MJy/sr per kelvin.

    Converts a small temperature change into a change in brightness. This is
    why CMB maps are quoted in "K_CMB": a map pixel reading +100 uK means the
    sky there is exactly as bright as a blackbody 100 uK hotter, at every
    frequency.
    """
    nu = np.asarray(nu_ghz, float) * 1e9
    x = H * nu / (K_B * temperature_k)
    ex = np.exp(x)
    return (2.0 * H * nu**3 / C**2) * (x / temperature_k) * ex / (ex - 1.0) ** 2 / MJY_PER_SR


def y_distortion(nu_ghz, temperature_k):
    """Spectral shape of a Compton-y distortion, per unit y, in MJy/sr.

    When CMB photons scatter off hot electrons they are nudged to higher
    energies: the spectrum stays a Planck-like shape but is no longer a
    blackbody. To first order the change is::

        dI = y * T * dB/dT * (x coth(x/2) - 4),    x = h nu / k T

    -- a deficit below ~217 GHz and an excess above. Any energy released into
    the early universe after about a year leaves a distortion like this, which
    is why FIRAS's null result (|y| < 1.5e-5) is such a strong constraint.
    """
    nu = np.asarray(nu_ghz, float) * 1e9
    x = H * nu / (K_B * temperature_k)
    return temperature_k * planck_intensity_dT(nu_ghz, temperature_k) * (
        x / np.tanh(x / 2.0) - 4.0)


def wien_peak_ghz(temperature_k: float) -> float:
    """Frequency at which B_nu peaks: ``h nu = 2.821439 k T`` (Wien's law, per Hz)."""
    return 2.8214393721220787 * K_B * temperature_k / H / 1e9


def photon_number_density(temperature_k: float) -> float:
    """Blackbody photons per cubic centimetre: ``(2 zeta(3) / pi^2) (k T / hbar c)^3``."""
    zeta3 = 1.2020569031595942
    n = 2.0 * zeta3 / np.pi**2 * (K_B * temperature_k / (HBAR * C)) ** 3
    return float(n / 1e6)


# --------------------------------------------------------------------------
# Densities
# --------------------------------------------------------------------------

#: Critical density for h = 1: rho_c = 3 (100 km/s/Mpc)^2 / (8 pi G), in kg/m^3.
RHO_CRIT_H2 = float((3 * (100 * u.km / u.s / u.Mpc) ** 2 / (8 * np.pi * const.G)).si.value)

#: Temperature of the relic neutrinos relative to the photons, (4/11)^(1/3):
#: electron-positron annihilation heated the photons after neutrinos decoupled.
NEUTRINO_TEMPERATURE_RATIO = (4.0 / 11.0) ** (1.0 / 3.0)


def omega_gamma_h2(temperature_k: float) -> float:
    """Physical photon density Omega_gamma h^2 of a blackbody at T.

    Energy density ``a T^4`` with the radiation constant ``a = 4 sigma / c``;
    divided by c^2 to make it a mass density, and by the critical density.
    """
    a_rad = 4.0 * const.sigma_sb.si.value / C
    return a_rad * temperature_k**4 / C**2 / RHO_CRIT_H2


def omega_radiation_h2(temperature_k: float, neff: float = 3.046) -> float:
    """Photons plus relativistic neutrinos: ``omega_gamma (1 + 7/8 (4/11)^(4/3) N_eff)``.

    The 7/8 is the fermion-to-boson energy ratio; N_eff = 3.046 rather than 3
    because neutrino decoupling was not quite complete when electrons and
    positrons annihilated.
    """
    return omega_gamma_h2(temperature_k) * (
        1.0 + 7.0 / 8.0 * NEUTRINO_TEMPERATURE_RATIO**4 * neff)


def z_equality(omega_m_h2: float, temperature_k: float = 2.7255, neff: float = 3.046) -> float:
    """Redshift of matter-radiation equality: ``1 + z_eq = omega_m / omega_r``.

    Before z_eq the universe's expansion was driven by radiation, after it by
    matter. Its value decides how much the dark-matter potential wells decayed
    while the acoustic oscillations were underway -- which is what the third
    peak's height measures. ``omega_m_h2`` should count only matter that was
    non-relativistic at z_eq (so it excludes the neutrinos).
    """
    return omega_m_h2 / omega_radiation_h2(temperature_k, neff) - 1.0


# --------------------------------------------------------------------------
# Angular power spectra
# --------------------------------------------------------------------------


def dl_from_cl(ell, cl):
    """``D_l = l (l + 1) C_l / 2 pi``.

    C_l is the variance per spherical-harmonic mode at multipole l. Plotting
    D_l instead is a convention: for a scale-invariant spectrum of primordial
    fluctuations it is flat on large scales, so departures from flat are the
    physics. D_l is also, roughly, the contribution per log-interval of l to the
    total variance of the map.
    """
    ell = np.asarray(ell, float)
    return ell * (ell + 1.0) * np.asarray(cl, float) / (2.0 * np.pi)


def cl_from_dl(ell, dl):
    """Inverse of ``dl_from_cl``; returns 0 where l = 0."""
    ell = np.asarray(ell, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = 2.0 * np.pi * np.asarray(dl, float) / (ell * (ell + 1.0))
    return np.where(ell > 0, out, 0.0)


def cosmic_variance_fraction(ell, fsky: float = 1.0):
    """Irreducible fractional uncertainty on C_l: ``sqrt(2 / ((2l + 1) f_sky))``.

    At multipole l the sky holds only 2l + 1 independent modes. Measuring the
    variance of a Gaussian from that many samples has this fractional error even
    with a perfect instrument -- we have one sky, and at l = 2 it offers five
    numbers. This, not noise, is why the largest scales are so uncertain.
    """
    ell = np.asarray(ell, float)
    return np.sqrt(2.0 / ((2.0 * ell + 1.0) * fsky))


def knox_error(ell, cl, nl=0.0, *, fsky: float = 1.0, delta_ell: int = 1, cross_nl=None):
    """Gaussian uncertainty of a band power (Knox 1995).

    For an auto-spectrum with noise ``nl``::

        sigma(C_l) = sqrt(2 / ((2l + 1) f_sky delta_l)) * (C_l + N_l)

    For the cross-spectrum of two maps with independent noise ``nl`` and
    ``cross_nl``, the variance is instead
    ``(C_l^2 + (C_l + N1)(C_l + N2)) / ((2l + 1) f_sky delta_l)`` -- the noise
    does not bias the mean, but it still adds scatter.
    """
    ell = np.asarray(ell, float)
    modes = (2.0 * ell + 1.0) * fsky * delta_ell
    cl = np.asarray(cl, float)
    if cross_nl is None:
        return np.sqrt(2.0 / modes) * (cl + nl)
    return np.sqrt((cl**2 + (cl + nl) * (cl + cross_nl)) / modes)


def multipole_to_degrees(ell):
    """Characteristic angular scale of multipole l: ``180 deg / l``."""
    return 180.0 / np.asarray(ell, float)


def planck_bin_edges(lmin: int = 30, lmax: int = 2508, width: int = 30):
    """Bin edges of Planck's binned high-l spectrum: [30, 59], [60, 89], ..., [2490, 2508]."""
    lo = np.arange(lmin, lmax + 1, width)
    hi = np.minimum(lo + width - 1, lmax)
    return lo, hi


def planck_bin(ell, dl, *, lmin: int = 30, lmax: int = 2508, width: int = 30):
    """Bin a D_l spectrum exactly as the Planck Legacy Archive does.

    Planck averages C_l (not D_l) with weights l(l + 1) within each bin, places
    the band power at the weighted-mean multipole l_eff, and converts back::

        C_b = sum(l(l+1) C_l) / sum(l(l+1)),    l_eff = sum(l(l+1) l) / sum(l(l+1))
        D_b = l_eff (l_eff + 1) C_b / 2 pi

    This recipe was not documented alongside the file; it was reverse-engineered
    here and reproduces Planck's own binned best-fit column to 3 parts in 10^9,
    against 5% for the obvious guess of averaging D_l directly.

    Returns ``(l_eff, D_b)``. ``ell`` must be a contiguous run of integers
    covering the requested range.
    """
    ell = np.asarray(ell)
    dl = np.asarray(dl, float)
    lookup = np.full(int(ell.max()) + 1, np.nan)
    lookup[ell.astype(int)] = dl
    lo, hi = planck_bin_edges(lmin, lmax, width)
    l_eff, d_b = np.empty(len(lo)), np.empty(len(lo))
    for i, (a, z) in enumerate(zip(lo, hi)):
        ls = np.arange(a, z + 1, dtype=float)
        w = ls * (ls + 1.0)
        c_l = cl_from_dl(ls, lookup[a:z + 1])
        l_eff[i] = (w * ls).sum() / w.sum()
        d_b[i] = dl_from_cl(l_eff[i], (w * c_l).sum() / w.sum())
    return l_eff, d_b


def map_rms_from_cl(ell, cl, beam=None) -> float:
    """Standard deviation of a map with power spectrum C_l (and optional beam B_l).

    Each multipole contributes ``(2l + 1) C_l B_l^2 / 4 pi`` to the pixel
    variance: 2l + 1 modes, spread over 4 pi steradians. This is the bridge
    between the map and the spectrum -- the spectrum predicts how rough the map
    should look.
    """
    ell = np.asarray(ell, float)
    cl = np.asarray(cl, float)
    b2 = 1.0 if beam is None else np.asarray(beam, float) ** 2
    return float(np.sqrt(np.sum((2.0 * ell + 1.0) * cl * b2) / (4.0 * np.pi)))


# --------------------------------------------------------------------------
# Beams and masks
# --------------------------------------------------------------------------


def gaussian_beam(ell, fwhm_arcmin: float):
    """Transfer function of a Gaussian beam: ``B_l = exp(-l(l+1) sigma^2 / 2)``.

    A telescope blurs the sky with its beam; in harmonic space the blur is a
    multiplication of each a_lm by B_l, so the measured spectrum is ``B_l^2 C_l``.
    ``sigma = FWHM / sqrt(8 ln 2)``.
    """
    sigma = np.radians(fwhm_arcmin / 60.0) / np.sqrt(8.0 * np.log(2.0))
    ell = np.asarray(ell, float)
    return np.exp(-0.5 * ell * (ell + 1.0) * sigma**2)


def _legendre_sum(coeffs_list, mu):
    """Evaluate ``sum_l c_l P_l(mu)`` for several coefficient arrays at once.

    Uses the three-term recurrence ``(l+1) P_{l+1} = (2l+1) mu P_l - l P_{l-1}``,
    which is numerically stable upward for |mu| <= 1, and never stores more than
    two Legendre polynomials at a time.
    """
    lmax = max(len(c) for c in coeffs_list) - 1
    out = [np.zeros_like(mu) for _ in coeffs_list]
    p_prev, p = np.zeros_like(mu), np.ones_like(mu)
    for ell in range(lmax + 1):
        for c, acc in zip(coeffs_list, out):
            if ell < len(c) and c[ell] != 0.0:
                acc += c[ell] * p
        p_prev, p = p, ((2 * ell + 1) * mu * p - ell * p_prev) / (ell + 1)
    return out


def _legendre_project(values_list, mu, weights, lmax: int):
    """``2 pi * integral f(mu) P_l(mu) dmu`` for l = 0..lmax, by Gauss-Legendre quadrature."""
    out = [np.zeros(lmax + 1) for _ in values_list]
    wf = [weights * v for v in values_list]
    p_prev, p = np.zeros_like(mu), np.ones_like(mu)
    for ell in range(lmax + 1):
        for f, acc in zip(wf, out):
            acc[ell] = 2.0 * np.pi * np.dot(f, p)
        p_prev, p = p, ((2 * ell + 1) * mu * p - ell * p_prev) / (ell + 1)
    return out


def deconvolve_mask(pseudo_cls, mask_cl, *, lmax: int, n_nodes: int | None = None):
    """Remove a mask's mode coupling exactly, by working with correlation functions.

    Multiplying the sky by a mask W mixes multipoles: power leaks from bright
    scales into faint ones, filling the spectrum's troughs and shaving its
    peaks. In harmonic space the mixing is a matrix; in angle space it is a
    simple product. For every pair of directions a separation theta apart, the
    masked map's correlation is the sky's correlation times the mask's::

        xi_masked(theta) = xi_sky(theta) * xi_W(theta)

    so we transform the pseudo-spectrum to a correlation function, divide by the
    mask's correlation function, and transform back. This is the idea behind the
    SpICE estimator (Szapudi et al. 2001, ApJ 548, L115; Chon et al. 2004,
    MNRAS 350, 914). It is exact for a statistically isotropic sky as long as
    xi_W never reaches zero -- true for a mask that keeps most of the sky, where
    even antipodal pairs are usually both unmasked.

    ``pseudo_cls`` is a list of pseudo-spectra (auto or cross) of maps that were
    all multiplied by the same mask; ``mask_cl`` is the power spectrum of that
    mask, ideally to about twice the map's lmax. Returns the corrected spectra
    for l = 0..lmax, not yet corrected for the beam.
    """
    from scipy.special import roots_legendre

    if n_nodes is None:
        # n Gauss-Legendre nodes integrate polynomials of degree 2n - 1 exactly.
        # The integrand pairs xi (degree = the pseudo-spectra's lmax) with P_l
        # (degree lmax); dividing by xi_W adds structure on the mask's scales.
        n_nodes = (max(len(c) for c in pseudo_cls) + lmax) // 2 + len(mask_cl) // 2 + 2
    mu, w = roots_legendre(n_nodes)
    mask_cl = np.asarray(mask_cl, float)
    ell_w = np.arange(len(mask_cl))
    xi_w = _legendre_sum([(2 * ell_w + 1) / (4 * np.pi) * mask_cl], mu)[0]
    if np.any(xi_w <= 0):
        raise ValueError("mask correlation function reaches zero; SpICE inversion undefined")
    coeffs = []
    for c in pseudo_cls:
        c = np.asarray(c, float)
        coeffs.append((2 * np.arange(len(c)) + 1) / (4 * np.pi) * c)
    xis = _legendre_sum(coeffs, mu)
    return _legendre_project([x / xi_w for x in xis], mu, w, lmax)


def mask_moments(mask) -> dict:
    """Sky-fraction statistics of a (possibly apodised) mask.

    Multiplying a map by a mask W and taking its spectrum gives a "pseudo"
    spectrum which, for scales much smaller than the mask's features, is the
    true spectrum times ``w2 = <W^2>``, the mean of W^2 over the whole sphere.
    Dividing by w2 is the simplest correction for a cut sky (and the one used
    here); the error-bar formula wants ``f_sky,eff = w2^2 / w4`` instead,
    which accounts for the tapered edge being partly down-weighted.
    """
    w = np.asarray(mask, float)
    w2, w4 = float(np.mean(w**2)), float(np.mean(w**4))
    return {
        "fsky": float(np.mean(w > 0)),
        "w1": float(np.mean(w)),
        "w2": w2,
        "w4": w4,
        "fsky_eff": w2**2 / w4 if w4 > 0 else 0.0,
    }


# --------------------------------------------------------------------------
# The acoustic scale
# --------------------------------------------------------------------------


def z_star_hu_sugiyama(ombh2: float, omh2: float) -> float:
    """Redshift of photon decoupling, Hu & Sugiyama (1996, ApJ 471, 542) eq. E1.

    The redshift at which the optical depth to Thomson scattering reaches one:
    the "surface" of last scattering. A fitting formula to detailed
    recombination calculations, accurate to about 0.1-0.2% for realistic
    densities -- plenty for us, since the answer (~1090) enters the angular scale
    only weakly.
    """
    g1 = 0.0783 * ombh2 ** (-0.238) / (1.0 + 39.5 * ombh2**0.763)
    g2 = 0.560 / (1.0 + 21.1 * ombh2**1.81)
    return 1048.0 * (1.0 + 0.00124 * ombh2 ** (-0.738)) * (1.0 + g1 * omh2**g2)


def baryon_loading(z, ombh2: float, temperature_k: float = 2.7255):
    """``R = 3 rho_b / (4 rho_gamma)``: the baryons' inertia in the photon fluid.

    The fluid that oscillated before recombination was photons (which provide
    the pressure) plus baryons (which add mass but no pressure). R measures the
    baryon share; it grows as ``1/(1+z)`` because radiation dilutes faster.
    """
    return 0.75 * ombh2 / omega_gamma_h2(temperature_k) / (1.0 + np.asarray(z, float))


def sound_speed(z, ombh2: float, temperature_k: float = 2.7255):
    """Sound speed of the photon-baryon fluid, in km/s: ``c / sqrt(3 (1 + R))``.

    Pure radiation would carry sound at c/sqrt(3), 58% of light speed. Baryons
    make the fluid heavier and slow it down.
    """
    return C / 1e3 / np.sqrt(3.0 * (1.0 + baryon_loading(z, ombh2, temperature_k)))


def sound_horizon(cosmo, z: float, ombh2: float) -> float:
    """Comoving distance a sound wave travels before redshift z, in Mpc.

    ``r_s(z) = integral from z to infinity of c_s(z') dz' / H(z')``

    This is the ruler the CMB peaks measure: every acoustic mode that had
    exactly one (or two, three, ...) compression by recombination shows up as
    a peak, and the fundamental one's wavelength is set by r_s. It depends only
    on the physical baryon, matter and radiation densities -- that is, on early-
    universe physics we can calculate -- which is what makes it a standard ruler.

    ``cosmo`` is an ``astropy.cosmology`` object, which supplies H(z) with
    photons and (massive) neutrinos included. The integral is done in
    ``ln(1+z)`` because it spans five orders of magnitude in redshift.
    """
    h0 = cosmo.H0.to(u.km / u.s / u.Mpc).value
    temperature = cosmo.Tcmb0.to(u.K).value

    def integrand(lnz1):
        zz = np.expm1(lnz1)
        return sound_speed(zz, ombh2, temperature) / (h0 * cosmo.efunc(zz)) * (1.0 + zz)

    val, _ = integrate.quad(integrand, np.log1p(z), np.log1p(1e8), limit=200,
                            epsabs=0.0, epsrel=1e-9)
    return float(val)


def angular_scale(cosmo, z_star: float, ombh2: float) -> dict:
    """The angle the sound horizon subtends at last scattering.

    ``theta_* = r_s(z_*) / D_M(z_*)``

    D_M is the comoving angular-diameter distance -- it is where geometry
    enters. In a closed universe, light from last scattering converges on its
    way to us and the ruler looks bigger (peaks move to lower l); in an open one
    it looks smaller. This is the measurement behind "the universe is flat".

    Returns r_s [Mpc], D_M [Mpc], theta_* [rad] and ``l_A = pi / theta_*``, the
    acoustic multipole that sets the peak spacing.
    """
    r_s = sound_horizon(cosmo, z_star, ombh2)
    d_m = float(cosmo.comoving_transverse_distance(z_star).to(u.Mpc).value)
    theta = r_s / d_m
    return {"r_s_mpc": r_s, "d_m_mpc": d_m, "theta_star": theta, "ell_A": np.pi / theta}


def make_cosmology(*, h: float, ombh2: float, omch2: float, omk: float = 0.0,
                   mnu_ev: float = 0.06, neff: float = 3.046, tcmb: float = 2.7255):
    """An ``astropy.cosmology`` object in the Planck 2018 baseline conventions.

    One massive neutrino of ``mnu_ev`` plus two massless ones sharing N_eff.
    Dark energy is a cosmological constant fixed by ``Omega_Lambda = 1 - Omega_m
    - Omega_k - Omega_radiation`` (astropy handles the bookkeeping). Omega_m for
    astropy excludes neutrinos, which it tracks separately as they turn from
    radiation-like to matter-like.
    """
    from astropy.cosmology import LambdaCDM

    om0 = (ombh2 + omch2) / h**2
    base = LambdaCDM(H0=100 * h, Om0=om0, Ode0=0.7, Tcmb0=tcmb, Neff=neff,
                     m_nu=[mnu_ev, 0.0, 0.0] * u.eV, Ob0=ombh2 / h**2)
    # Choose Omega_Lambda so the total gives the requested curvature exactly.
    ode0 = 1.0 - omk - om0 - base.Ogamma0 - base.Onu0
    return LambdaCDM(H0=100 * h, Om0=om0, Ode0=ode0, Tcmb0=tcmb, Neff=neff,
                     m_nu=[mnu_ev, 0.0, 0.0] * u.eV, Ob0=ombh2 / h**2)
