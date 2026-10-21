"""The physics of one star orbiting one point mass, as seen from Earth.

Everything the orbit fit needs is here, written out rather than taken from a
package, so that each step can be read and tested:

* **Kepler's equation** gives where the star is along its ellipse at any time.
* **Three rotations** (the Thiele-Innes constants) turn that ellipse into what
  we see: an offset on the sky in milliarcseconds, plus a line-of-sight
  position and velocity.
* **Kepler's third law** ties the orbit's size to its period and the central
  mass: ``a**3 = G M P**2 / (4 pi**2)``. Astrometry sees ``a`` as an *angle*;
  radial velocities see the orbit's *speed* in km/s. Needing both to agree is
  what fixes the distance R0.

Three effects beyond Newton's point-mass orbit are included, each with a switch
so the fits can compare "Newton" with "Einstein":

* **Roemer delay** -- light from the far side of the orbit takes up to ~8 days
  longer to reach us than light from the near side. Plain light-travel time,
  not relativity, but Do et al. (2019) leave it out of their Newtonian model,
  so it has its own switch.
* **Relativistic redshift** -- time dilation (the transverse Doppler effect,
  ``v**2 / 2c``) plus the gravitational redshift (``GM / rc``), scaled by a
  factor ``f`` that is 0 for Newton and 1 for general relativity. For S2 it
  peaks at about 200 km/s at pericentre.
* **Schwarzschild precession** -- in general relativity the ellipse turns
  forward by ``6 pi GM / (c**2 a (1 - e**2))`` per orbit, about 12 arcminutes
  for S2. Modelled the way pulsar timing does it (Damour & Deruelle 1986): the
  argument of pericentre advances in proportion to the true anomaly, which is
  the exact first-post-Newtonian orbit shape. Scaled by ``f_sp`` (0 = Newton,
  1 = GR).

Coordinates follow Do et al. (2019) and the visual-binary convention: X points
East, Y North, Z *away* from us, along the line of sight to the Galactic Centre.
A positive radial velocity is a receding star. With these conventions the
published S2 elements (i = 134 deg, Omega = 227 deg, omega = 66 deg) describe
the measured orbit: i > 90 deg means the star moves clockwise on the sky.

All functions broadcast: parameters may be arrays of shape ``(W, 1)`` (one row
per MCMC walker) against times of shape ``(N,)``, giving ``(W, N)`` outputs.

Units: time in Julian years (decimal year, J2000-based), GM in units of
10**6 GM_sun, R0 in kpc, angles on the sky in mas, velocities in km/s,
lengths in AU. One AU seen from one kpc subtends exactly one milliarcsecond,
which is why these units are convenient.
"""

from __future__ import annotations

import numpy as np

# --------------------------------------------------------------------------
# Constants (IAU 2012/2015 nominal values; CODATA c)
# --------------------------------------------------------------------------

C_KMS = 299_792.458  # speed of light, km/s
AU_KM = 149_597_870.700  # astronomical unit, km (IAU 2012, exact)
YEAR_S = 365.25 * 86_400.0  # Julian year, s
GM_SUN_KM3S2 = 1.327_124_400_18e11  # heliocentric gravitational constant, km^3/s^2
GM_SUN_AU3YR2 = GM_SUN_KM3S2 * YEAR_S**2 / AU_KM**3  # = 39.4769..., "4 pi^2" in AU, yr
KMS_PER_AUYR = AU_KM / YEAR_S  # 1 AU/yr = 4.74047 km/s
C_AUYR = C_KMS / KMS_PER_AUYR  # speed of light, 63,241 AU/yr
PC_AU = 648_000.0 / np.pi  # one parsec in AU (exact, by definition)
M_SUN_KG = 1.988_409_870_698_051e30  # IAU nominal GM_sun / CODATA 2018 G
G_SI = 6.674_30e-11
PC_M = PC_AU * AU_KM * 1e3
MICROARCSEC_PER_RAD = 180.0 / np.pi * 3600.0 * 1e6

#: Reference epoch of the frame-drift terms (x0 + vx0 * (t - T_REF)).
#: Do et al. (2019) use J2000; we do the same for every data set.
T_REF = 2000.0


def mjd_to_year(mjd):
    """Modified Julian Date -> decimal Julian year (J2000.0 = 2000.0)."""
    return 2000.0 + (np.asarray(mjd, dtype=float) - 51_544.5) / 365.25


def year_to_mjd(year):
    return (np.asarray(year, dtype=float) - 2000.0) * 365.25 + 51_544.5


# --------------------------------------------------------------------------
# Kepler's equation and the orbit in its own plane
# --------------------------------------------------------------------------


def solve_kepler(mean_anomaly, e, *, tol: float = 1e-13, max_iter: int = 30):
    """Solve M = E - e sin E for the eccentric anomaly E (radians).

    Vectorised Newton-Raphson from Danby's starting guess, which converges in a
    handful of steps even at S2's e = 0.886. Returns E in the same branch as M
    (M is wrapped to [-pi, pi) first, then the multiple of 2 pi is restored).
    """
    mean_anomaly = np.asarray(mean_anomaly, dtype=float)
    e = np.asarray(e, dtype=float)
    turns = np.floor((mean_anomaly + np.pi) / (2.0 * np.pi))
    m = mean_anomaly - 2.0 * np.pi * turns  # in [-pi, pi)
    ecc = np.broadcast_to(e, m.shape)
    big_e = m + 0.85 * ecc * np.sign(np.sin(m))
    for _ in range(max_iter):
        f = big_e - ecc * np.sin(big_e) - m
        fp = 1.0 - ecc * np.cos(big_e)
        step = f / fp
        big_e = big_e - step
        if np.all(np.abs(step) < tol):
            break
    return big_e + 2.0 * np.pi * turns, turns


def semi_major_axis_au(gm6, period_yr):
    """Kepler's third law: a = (GM P^2 / 4 pi^2)^(1/3), in AU.

    ``gm6`` is GM in units of 10^6 GM_sun, so for S2 (P = 16.05 yr,
    GM = 4.1) this is about 1,030 AU.
    """
    gm = np.asarray(gm6) * 1e6 * GM_SUN_AU3YR2
    return (gm * np.asarray(period_yr) ** 2 / (4.0 * np.pi**2)) ** (1.0 / 3.0)


def mass_from_orbit(a_au, period_yr):
    """Kepler's third law inverted: the central mass (in M_sun) from a and P."""
    return 4.0 * np.pi**2 * np.asarray(a_au) ** 3 / (GM_SUN_AU3YR2 * np.asarray(period_yr) ** 2)


def precession_per_orbit(gm6, a_au, e):
    """Schwarzschild precession of the pericentre per orbit, radians.

    Delta omega = 6 pi G M / (c^2 a (1 - e^2)).  About 0.0034 rad = 11.6' for S2.
    """
    rg_au = np.asarray(gm6) * 1e6 * GM_SUN_AU3YR2 / C_AUYR**2  # GM/c^2 in AU
    return 6.0 * np.pi * rg_au / (np.asarray(a_au) * (1.0 - np.asarray(e) ** 2))


def orbit_state(t, gm6, period, t0, e, inc_deg, omega_deg, Omega_deg, *, f_sp=0.0):
    """Position (AU) and velocity (km/s) of the star relative to the mass.

    Returns ``(X, Y, Z, VX, VY, VZ, r)`` with X East, Y North, Z away from us.
    ``t0`` is the time of the reference pericentre passage; with precession
    switched on (``f_sp`` > 0) the elements refer to that passage and the
    ellipse turns by ``f_sp * Delta omega`` per orbit either side of it.
    """
    t = np.asarray(t, dtype=float)
    a = semi_major_axis_au(gm6, period)
    n = 2.0 * np.pi / period  # mean motion, rad/yr
    big_e, turns = solve_kepler(n * (t - t0), e)
    cos_e, sin_e = np.cos(big_e), np.sin(big_e)
    sq = np.sqrt(1.0 - e**2)
    one_minus = 1.0 - e * cos_e
    r = a * one_minus
    # True anomaly, unwrapped so that it keeps counting through each orbit.
    e_wrapped = big_e - 2.0 * np.pi * turns
    nu = 2.0 * np.arctan2(np.sqrt(1.0 + e) * np.sin(e_wrapped / 2.0),
                          np.sqrt(1.0 - e) * np.cos(e_wrapped / 2.0))
    nu_unwrapped = nu + 2.0 * np.pi * turns

    k = f_sp * precession_per_orbit(gm6, a, e) / (2.0 * np.pi)  # dimensionless
    inc = np.radians(inc_deg)
    big_om = np.radians(Omega_deg)
    u = np.radians(omega_deg) + nu + k * nu_unwrapped  # argument of latitude

    cu, su = np.cos(u), np.sin(u)
    co, so = np.cos(big_om), np.sin(big_om)
    ci, si = np.cos(inc), np.sin(inc)

    # Unit vectors in the sky frame for "radial" and "along-track" directions.
    north_r = cu * co - su * so * ci
    east_r = cu * so + su * co * ci
    away_r = su * si
    north_t = -su * co - cu * so * ci
    east_t = -su * so + cu * co * ci
    away_t = cu * si

    rdot = n * a * e * sin_e / one_minus  # AU/yr
    r_udot = (1.0 + k) * n * a * sq / one_minus  # AU/yr

    x_east = r * east_r
    y_north = r * north_r
    z_away = r * away_r
    vx = (rdot * east_r + r_udot * east_t) * KMS_PER_AUYR
    vy = (rdot * north_r + r_udot * north_t) * KMS_PER_AUYR
    vz = (rdot * away_r + r_udot * away_t) * KMS_PER_AUYR
    return x_east, y_north, z_away, vx, vy, vz, r


def relativistic_redshift_kms(gm6, r_au, speed_kms):
    """Transverse Doppler + gravitational redshift, in km/s of apparent velocity.

    To first order in 1/c^2 a clock moving at speed v deep in a potential GM/r
    ticks slow by v^2/2c^2 + GM/(r c^2); multiplied by c this shows up as an
    extra, always-positive "velocity". For S2 at pericentre (r = 120 AU,
    v = 7,700 km/s) the two terms are ~100 km/s each.
    """
    gm_km = np.asarray(gm6) * 1e6 * GM_SUN_KM3S2
    return speed_kms**2 / (2.0 * C_KMS) + gm_km / (np.asarray(r_au) * AU_KM * C_KMS)


def roemer_emission_time(t_obs, gm6, period, t0, e, inc_deg, omega_deg, Omega_deg, *, f_sp=0.0):
    """Emission time of light observed at ``t_obs``: t_em = t_obs - Z(t_obs)/c.

    One fixed-point iteration, as in Do et al. (2019); a second changes the
    answer by < 20 minutes. Only the *variation* of the delay matters -- its
    constant part (R0/c, ~26,000 years) is absorbed into t0.
    """
    z = orbit_state(t_obs, gm6, period, t0, e, inc_deg, omega_deg, Omega_deg, f_sp=f_sp)[2]
    return t_obs - z / C_AUYR


def predict(params: dict, t_astro, t_rv, *, roemer: bool = True):
    """Model observables for a parameter dictionary.

    ``params`` holds (arrays or scalars): gm6, r0 (kpc), period, t0, e, inc,
    omega, Omega (deg), f (redshift factor), f_sp (precession factor).
    Returns sky offsets (x_east, y_north) in mas *relative to the mass* at
    ``t_astro`` and the line-of-sight velocity in km/s (without any vz0 offset)
    at ``t_rv``. Frame offsets and drifts are added by the caller.
    """
    el = (params["gm6"], params["period"], params["t0"], params["e"], params["inc"],
          params["omega"], params["Omega"])
    f_sp = params.get("f_sp", 0.0)
    f_rs = params.get("f", 0.0)

    out = {}
    if t_astro is not None:
        t_em = roemer_emission_time(t_astro, *el, f_sp=f_sp) if roemer else t_astro
        x, y, *_ = orbit_state(t_em, *el, f_sp=f_sp)
        out["x"] = x / params["r0"]
        out["y"] = y / params["r0"]
    if t_rv is not None:
        t_em = roemer_emission_time(t_rv, *el, f_sp=f_sp) if roemer else t_rv
        _, _, _, vx, vy, vz, r = orbit_state(t_em, *el, f_sp=f_sp)
        rv = vz
        if np.any(np.asarray(f_rs) != 0):
            rv = rv + f_rs * relativistic_redshift_kms(params["gm6"], r, np.sqrt(vx**2 + vy**2 + vz**2))
        out["rv"] = rv
    return out


# --------------------------------------------------------------------------
# Brute-force check: integrate the post-Newtonian equations of motion
# --------------------------------------------------------------------------


def integrate_1pn(t_eval, gm6, period, t_osc, e, inc_deg, omega_deg, Omega_deg, *,
                  t_start=2000.0, pn=1.0, rtol=1e-12):
    """Numerically integrate S2's orbit with the 1PN test-particle equations.

    This is what Do et al. (2019, eq. S2) do inside their fit:

        d2R/dt2 = -GM R / R^3 + pn * GM / (c^2 R^3) * [(4 GM / R - v^2) R + 4 (R.V) V]

    (Schwarzschild metric, harmonic coordinates), starting from the
    *osculating* Keplerian elements at ``t_start`` (J2000 for them; ``t_osc``
    is the pericentre time of that osculating ellipse). Far too slow for an
    MCMC in Python, it serves as the reference that the fast closed-form model
    (``orbit_state`` with ``f_sp``) is tested against. Returns the state at
    ``t_eval``: X, Y, Z (AU), VX, VY, VZ (km/s).
    """
    from scipy.integrate import solve_ivp

    gm = gm6 * 1e6 * GM_SUN_AU3YR2
    x, y, z, vx, vy, vz, _ = orbit_state(np.array([t_start]), gm6, period, t_osc, e,
                                         inc_deg, omega_deg, Omega_deg)
    y0 = np.array([x[0], y[0], z[0], vx[0] / KMS_PER_AUYR, vy[0] / KMS_PER_AUYR,
                   vz[0] / KMS_PER_AUYR])
    c2 = C_AUYR**2

    def rhs(_, s):
        r_vec, v_vec = s[:3], s[3:]
        r = np.sqrt(r_vec @ r_vec)
        v2 = v_vec @ v_vec
        acc = -gm * r_vec / r**3
        acc += pn * gm / (c2 * r**3) * ((4.0 * gm / r - v2) * r_vec + 4.0 * (r_vec @ v_vec) * v_vec)
        return np.concatenate([v_vec, acc])

    t_eval = np.atleast_1d(np.asarray(t_eval, dtype=float))
    out = np.zeros((6, t_eval.size))
    for sel, t_end in ((t_eval >= t_start, max(t_eval.max(), t_start) + 1e-9),
                       (t_eval < t_start, min(t_eval.min(), t_start) - 1e-9)):
        if not np.any(sel):
            continue
        ts = t_eval[sel]
        order = np.argsort(ts) if t_end > t_start else np.argsort(-ts)
        sol = solve_ivp(rhs, (t_start, t_end), y0, method="DOP853", t_eval=ts[order],
                        rtol=rtol, atol=1e-12 * np.abs(y0).max())
        block = np.empty((6, ts.size))
        block[:, order] = sol.y
        out[:, sel] = block
    out[3:] *= KMS_PER_AUYR
    return out


# --------------------------------------------------------------------------
# Derived quantities
# --------------------------------------------------------------------------


def pericentre_au(gm6, period, e):
    return semi_major_axis_au(gm6, period) * (1.0 - np.asarray(e))


def apocentre_au(gm6, period, e):
    return semi_major_axis_au(gm6, period) * (1.0 + np.asarray(e))


def pericentre_speed_kms(gm6, period, e):
    """Vis-viva at pericentre: v^2 = GM (1 + e) / (a (1 - e))."""
    a = semi_major_axis_au(gm6, period)
    gm = np.asarray(gm6) * 1e6 * GM_SUN_AU3YR2
    return np.sqrt(gm * (1.0 + e) / (a * (1.0 - e))) * KMS_PER_AUYR


def schwarzschild_radius_au(gm6):
    """R_s = 2GM/c^2 -- about 0.08 AU (17 solar radii) for 4 million suns."""
    return 2.0 * np.asarray(gm6) * 1e6 * GM_SUN_AU3YR2 / C_AUYR**2


def shadow_diameter_uas(gm6, r0_kpc):
    """Angular diameter of a non-spinning black hole's shadow, micro-arcsec.

    Light rays passing within an impact parameter sqrt(27) GM/c^2 are captured,
    so the shadow's diameter is 2 sqrt(27) GM/c^2 ~ 10.4 GM/c^2 (Bardeen 1973).
    Spin changes this by at most ~4% for any viewing angle (Johannsen & Psaltis
    2010), which the caller treats as a systematic.
    """
    theta_g = gravitational_angle_uas(gm6, r0_kpc)
    return 2.0 * np.sqrt(27.0) * theta_g


def gravitational_angle_uas(gm6, r0_kpc):
    """theta_g = GM / (c^2 D): the black hole's natural angular unit, micro-arcsec."""
    rg_au = np.asarray(gm6) * 1e6 * GM_SUN_AU3YR2 / C_AUYR**2
    return rg_au / (np.asarray(r0_kpc) * 1e3 * PC_AU) * MICROARCSEC_PER_RAD
