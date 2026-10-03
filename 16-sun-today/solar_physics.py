"""The physics behind every number in this project, in small testable pieces.

Nothing here reads data; ``test_physics.py`` checks each function against a
textbook value or an exact limiting case.
"""

from __future__ import annotations

import numpy as np
from astropy import constants as const

import config as C

K_B = const.k_B.si.value  # J/K
M_H = const.m_p.si.value + const.m_e.si.value  # hydrogen atom, kg
H_PLANCK = const.h.si.value
C_LIGHT = const.c.si.value

# --- The corona as a hydrostatic atmosphere ------------------------------------------------


def mean_particle_mass(he_abundance: float = C.CORONA_HE_ABUNDANCE) -> float:
    """Mean mass per free particle of fully ionised H + He plasma, in units of m_H.

    Per hydrogen nucleus there are 1 proton, y helium nuclei and 1 + 2y electrons
    (y = n_He/n_H), weighing (1 + 4y) m_H in all: mu = (1 + 4y) / (2 + 3y).
    """
    y = he_abundance
    return (1 + 4 * y) / (2 + 3 * y)


def surface_gravity() -> float:
    """g at the photosphere, m/s^2 (274 m/s^2: 28 times Earth's)."""
    return C.GM_SUN / C.R_SUN_M**2


def pressure_scale_height_m(temperature_k, he_abundance: float = C.CORONA_HE_ABUNDANCE):
    """H = kT / (mu m_H g): the height over which pressure drops by a factor e in an
    isothermal atmosphere at the base gravity. About 50 Mm for 1 MK."""
    return K_B * np.asarray(temperature_k) / (mean_particle_mass(he_abundance) * M_H
                                              * surface_gravity())


def temperature_from_scale_height(h_m, he_abundance: float = C.CORONA_HE_ABUNDANCE):
    """Inverse of :func:`pressure_scale_height_m`."""
    return np.asarray(h_m) * mean_particle_mass(he_abundance) * M_H * surface_gravity() / K_B


def hydrostatic_density(r_rsun, temperature_k):
    """n(r)/n(R_sun) for an isothermal atmosphere in the Sun's 1/r^2 gravity:
    exp[-(R/H)(1 - R/r)]. Gravity weakening with height makes it fall slower than a
    plain exponential; far out it levels off (that is why the corona must expand
    into the solar wind)."""
    lam = C.R_SUN_M / pressure_scale_height_m(temperature_k)
    return np.exp(-lam * (1.0 - 1.0 / np.asarray(r_rsun, dtype=float)))


def offlimb_intensity(p_rsun, temperature_k, *, s_max_rsun: float = 3.0, n_s: int = 3001):
    """Optically thin line emission seen at impact parameter p > 1 (in R_sun).

    EUV lines are collisionally excited, so emission per volume goes as n_e^2; the
    image records its integral along the line of sight, r = sqrt(p^2 + s^2). The
    normalisation is arbitrary (only the shape with p is used).
    """
    p = np.atleast_1d(np.asarray(p_rsun, dtype=float))
    s = np.linspace(0.0, s_max_rsun, n_s)
    r = np.sqrt(p[:, None] ** 2 + s[None, :] ** 2)
    integrand = hydrostatic_density(r, temperature_k) ** 2
    return 2.0 * np.trapezoid(integrand, s, axis=1)


def fit_offlimb_temperature(height_mm, intensity, *, t_grid_mk=None):
    """The isothermal hydrostatic temperature that best matches an off-limb profile.

    Fits log(I) = log(A) + log(model(p; T)) with A free, by scanning T and solving
    for A in closed form at each T; the minimum is refined with a parabola. Returns
    (T in K, rms residual in dex, best-fit model on the input heights).
    """
    h = np.asarray(height_mm, dtype=float)
    y = np.log10(np.asarray(intensity, dtype=float))
    good = np.isfinite(y)
    h, y = h[good], y[good]
    p = 1.0 + h * 1e6 / C.R_SUN_M
    if t_grid_mk is None:
        t_grid_mk = np.geomspace(0.2, 6.0, 300)
    chi = np.empty(len(t_grid_mk))
    for i, tmk in enumerate(t_grid_mk):
        model = np.log10(offlimb_intensity(p, tmk * 1e6, n_s=801))
        a = np.mean(y - model)
        chi[i] = np.sum((y - model - a) ** 2)
    i = int(np.argmin(chi))
    t_best = t_grid_mk[i]
    if 0 < i < len(chi) - 1:  # parabolic refinement in log T
        x = np.log(t_grid_mk[i - 1:i + 2])
        c2, c1, _ = np.polyfit(x, chi[i - 1:i + 2], 2)
        if c2 > 0:
            t_best = float(np.exp(-c1 / (2 * c2)))
    model = np.log10(offlimb_intensity(p, t_best * 1e6, n_s=801))
    a = np.mean(y - model)
    rms = float(np.sqrt(np.mean((y - model - a) ** 2)))
    return t_best * 1e6, rms, 10 ** (model + a)


def emission_scale_height_mm(height_mm, intensity) -> float:
    """The plain exponential e-folding height of the intensity (Mm), from a straight
    line fit to log I(h). For n^2 emission it is half the density scale height."""
    h = np.asarray(height_mm, dtype=float)
    y = np.log(np.asarray(intensity, dtype=float))
    good = np.isfinite(y)
    slope = np.polyfit(h[good], y[good], 1)[0]
    return -1.0 / slope


# --- Rotation -------------------------------------------------------------------------------


def rotation_law(lat_deg, a, b, c=0.0):
    """omega(lat) = A + B sin^2(lat) + C sin^4(lat): the standard differential-rotation
    profile (deg/day when A, B, C are)."""
    s2 = np.sin(np.radians(lat_deg)) ** 2
    return a + b * s2 + c * s2**2


def fit_rotation_law(lat_deg, omega, sigma, *, fit_c: bool = True):
    """Weighted least-squares A, B(, C) with their covariance matrix."""
    lat = np.asarray(lat_deg, dtype=float)
    w = 1.0 / np.asarray(sigma, dtype=float) ** 2
    s2 = np.sin(np.radians(lat)) ** 2
    cols = [np.ones_like(s2), s2] + ([s2**2] if fit_c else [])
    x = np.column_stack(cols)
    xtw = x.T * w
    cov = np.linalg.inv(xtw @ x)
    coef = cov @ (xtw @ np.asarray(omega, dtype=float))
    resid = np.asarray(omega) - x @ coef
    chi2 = float(np.sum(w * resid**2))
    dof = len(lat) - len(coef)
    return coef, cov, chi2, dof


def synodic_period_days(sidereal_deg_per_day: float, earth_deg_per_day: float) -> float:
    """How long one turn takes as seen from Earth, which orbits the same way."""
    return 360.0 / (sidereal_deg_per_day - earth_deg_per_day)


def earth_orbital_rate_deg_per_day(t) -> float:
    """Earth's angular speed around the Sun at time t (deg/day), from the change in
    the Sun's apparent ecliptic longitude over +-12 h (0.953 in July, 1.019 in
    January: Kepler's second law)."""
    import astropy.units as u
    from astropy.time import Time
    from sunpy.coordinates import sun

    t = Time(t)
    l1 = sun.apparent_longitude(t + 0.5 * u.day).to_value(u.deg)
    l0 = sun.apparent_longitude(t - 0.5 * u.day).to_value(u.deg)
    return float((l1 - l0 + 180.0) % 360.0 - 180.0)


# --- Radiation from the photosphere ------------------------------------------------------


def planck_lambda(wavelength_m, temperature_k):
    """Planck spectral radiance B_lambda(T), W m^-2 sr^-1 m^-1."""
    lam = np.asarray(wavelength_m, dtype=float)
    t = np.asarray(temperature_k, dtype=float)
    return 2 * H_PLANCK * C_LIGHT**2 / lam**5 / np.expm1(H_PLANCK * C_LIGHT / (lam * K_B * t))


def planck_temperature(wavelength_m, radiance):
    """The temperature whose Planck radiance at this wavelength is ``radiance``."""
    lam = np.asarray(wavelength_m, dtype=float)
    x = 2 * H_PLANCK * C_LIGHT**2 / (lam**5 * np.asarray(radiance, dtype=float))
    return H_PLANCK * C_LIGHT / (lam * K_B * np.log1p(x))


def allen_limb_coefficients(wavelength_a: float) -> tuple[float, float]:
    """(u, v) of Allen's limb-darkening law at a continuum wavelength (4000-15000 A).

    The fifth-order polynomial fits in wavelength to the tables of Allen's
    Astrophysical Quantities that SolarSoft uses (darklimb_u.pro, darklimb_v.pro;
    D. Alexander 1996).
    """
    if not 4000 <= wavelength_a <= 15000:
        raise ValueError("Allen's coefficients are tabulated for 4000-15000 A only")
    au = (-8.9829751, 0.0069093916, -1.8144591e-6, 2.2540875e-10, -1.3389747e-14, 3.0453572e-19)
    av = (9.2891180, -0.0062212632, 1.5788029e-6, -1.9359644e-10, 1.1444469e-14, -2.599494e-19)
    powers = np.array([wavelength_a**k for k in range(6)])
    return float(np.dot(au, powers)), float(np.dot(av, powers))


def allen_limb_darkening(mu, wavelength_a: float):
    """I(mu)/I(1) = 1 - u - v + u mu + v mu^2."""
    u_, v_ = allen_limb_coefficients(wavelength_a)
    mu = np.asarray(mu, dtype=float)
    return 1 - u_ - v_ + u_ * mu + v_ * mu**2


def eddington_limb_darkening(mu):
    """Grey atmosphere in the Eddington approximation: I(mu)/I(1) = (2 + 3 mu)/5."""
    return (2 + 3 * np.asarray(mu, dtype=float)) / 5


def grey_temperature(tau):
    """Eddington grey-atmosphere temperature profile T^4 = (3/4) T_eff^4 (tau + 2/3)."""
    return C.T_EFF_K * (0.75 * (np.asarray(tau, dtype=float) + 2.0 / 3.0)) ** 0.25


def eddington_barbier_temperature(mu, intensity_ratio, wavelength_m, *,
                                  ratio_at_anchor: float, t_anchor: float = C.T_EFF_K):
    """Temperature at optical depth tau = mu, from the limb-darkening curve.

    Eddington-Barbier: the intensity leaving at angle arccos(mu) equals the source
    function at optical depth mu along the vertical, I(mu) = B(T(tau = mu)) (exact
    when B is linear in tau). Relative intensities fix relative Planck radiances; one
    anchor sets the scale: we take T(tau = 2/3) = T_eff, where a grey atmosphere
    radiates its effective temperature. ``ratio_at_anchor`` is I(2/3)/I(1).
    """
    b_anchor = planck_lambda(wavelength_m, t_anchor)
    radiance = b_anchor * np.asarray(intensity_ratio, dtype=float) / ratio_at_anchor
    return planck_temperature(wavelength_m, radiance)


# --- Flares -----------------------------------------------------------------------------------


def thermal_energy_erg(temperature_k, em_cm3, volume_cm3):
    """Thermal energy of an isothermal, fully ionised plasma, 3 n_e k T V.

    Electrons and ions each carry (3/2) k T; with n_e = n_i = sqrt(EM / V) this is
    3 k T sqrt(EM V). In erg (the unit flare papers use).
    """
    k_cgs = K_B * 1e7
    return 3.0 * k_cgs * np.asarray(temperature_k) * np.sqrt(np.asarray(em_cm3) * volume_cm3)


def radiative_cooling_time_s(temperature_k, density_cm3):
    """Time for an optically thin plasma to radiate away its thermal energy,
    3 k T / (n Lambda(T)), with the power-law loss function Lambda = 1.2e-19 T^-1/2
    erg cm^3 s^-1 used for 1-10 MK flare plasma by Cargill, Mariska & Antiochos (1995,
    after Rosner et al. 1978)."""
    k_cgs = K_B * 1e7
    t = np.asarray(temperature_k, dtype=float)
    lam = 1.2e-19 * t**-0.5
    return 3.0 * k_cgs * t / (np.asarray(density_cm3) * lam)


# --- Areas on the sphere ----------------------------------------------------------------------


def msh_per_pixel(pixel_mm: float, mu):
    """Area on the solar surface covered by one image pixel, in millionths of a solar
    hemisphere (the unit sunspot areas are quoted in): foreshortening stretches a
    pixel at angle arccos(mu) over 1/mu times more surface."""
    hemisphere_mm2 = 2 * np.pi * (C.R_SUN_M / 1e6) ** 2
    return pixel_mm**2 / np.asarray(mu, dtype=float) / hemisphere_mm2 * 1e6
