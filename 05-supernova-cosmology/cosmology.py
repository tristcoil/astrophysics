"""The physics: how far away something at redshift z is, in a given universe.

A universe here is a Friedmann-Lemaitre-Robertson-Walker model described by a
handful of numbers:

* ``H0`` -- the Hubble constant, today's expansion rate (km/s/Mpc).
* ``Om`` -- Omega_m, today's matter density as a fraction of the critical
  density (the density that makes space flat).
* ``OL`` -- Omega_Lambda (or Omega_DE), the same for dark energy.
* ``w``  -- dark energy's equation of state, pressure / (density c^2). A
  cosmological constant has exactly w = -1, and its density never dilutes.
* ``Ok`` = 1 - Om - OL -- the curvature term. Positive: open (saddle-shaped)
  space; negative: closed (sphere-like); zero: flat.

Everything follows from the expansion rate at redshift z, relative to today::

    E(z)^2 = H(z)^2 / H0^2 = Om (1+z)^3 + Ok (1+z)^2 + OL (1+z)^(3(1+w))

Matter dilutes as volume, (1+z)^3; curvature as area; a cosmological constant
not at all. Radiation is left out: today it is ~9e-5 of the critical density and
changes distances out to z = 2.3 by well under 0.1% (``test_physics.py`` checks
this against astropy with radiation switched on).

Distances are computed in units of the Hubble distance c/H0 unless stated
otherwise; H0 then only ever appears as an overall scale.

A supernova's *luminosity distance* d_L is defined so that flux = L / (4 pi d_L^2)
holds even in an expanding universe. It is what a standard candle measures::

    distance modulus  mu = m - M = 5 log10(d_L / 10 pc) = 5 log10(d_L / Mpc) + 25
"""

from __future__ import annotations

import numpy as np

from config import C_KM_S

#: Seconds per gigayear (Julian years) and kilometres per megaparsec, for ages.
SEC_PER_GYR = 1e9 * 365.25 * 86400.0
KM_PER_MPC = 3.0856775814913673e19


def hubble_distance_mpc(H0: float) -> float:
    """c/H0 in Mpc: the natural distance scale of a universe expanding at H0."""
    return C_KM_S / H0


def hubble_time_gyr(H0: float) -> float:
    """1/H0 in Gyr: how old the universe would be had it always expanded at H0."""
    return KM_PER_MPC / H0 / SEC_PER_GYR


def e_squared(z, Om, OL, w=-1.0):
    """E(z)^2 = (H(z)/H0)^2. Broadcasts over all arguments."""
    zp1 = 1.0 + np.asarray(z, dtype=float)
    Om, OL, w = np.asarray(Om, float), np.asarray(OL, float), np.asarray(w, float)
    Ok = 1.0 - Om - OL
    return Om * zp1**3 + Ok * zp1**2 + OL * zp1 ** (3.0 * (1.0 + w))


def integration_nodes(z_data, max_step: float = 0.003) -> tuple[np.ndarray, np.ndarray]:
    """Nodes for integrating from z = 0 out to every redshift in ``z_data``.

    The nodes are every data redshift plus a grid uniform in ln(1+z) that caps
    the spacing at ``max_step``, so the integral lands exactly on each supernova
    with no interpolation. Returns the nodes (as redshifts) and, for each data
    point, its node index.
    """
    x_data = np.log1p(np.asarray(z_data, dtype=float))
    top = float(x_data.max())
    fill = np.linspace(0.0, top, int(np.ceil(top / max_step)) + 1)
    x_nodes, inverse = np.unique(np.concatenate([fill, x_data]), return_inverse=True)
    return np.expm1(x_nodes), inverse[fill.size:]


def comoving_distance_nodes(nodes, Om, OL, w=-1.0) -> np.ndarray:
    """Line-of-sight comoving distance D_C = int_0^z dz'/E(z'), in c/H0 units.

    ``Om``, ``OL``, ``w`` may be arrays of equal length (a batch of universes);
    the result then has shape (n_models, n_nodes). Models where E^2 <= 0 somewhere
    in the range (a "bounce": the universe never reached that redshift) are NaN.

    The integral is taken in x = ln(1+z), where it becomes int (1+z)/E dx: that
    integrand is far smoother than 1/E(z) (constant for an empty universe, a
    gentle exponential for matter), so the trapezoid rule on nodes spaced
    <= 0.003 in x is accurate to ~1e-6 -- about 2e-6 mag in distance modulus,
    checked against closed forms, astropy and adaptive quadrature in the tests.
    """
    Om = np.atleast_1d(np.asarray(Om, float))[:, None]
    OL = np.atleast_1d(np.asarray(OL, float))[:, None]
    w = np.atleast_1d(np.asarray(w, float))[:, None]
    e2 = e_squared(nodes[None, :], Om, OL, w)
    bad = np.any(e2 <= 0.0, axis=1)
    integrand = (1.0 + nodes)[None, :] / np.sqrt(np.where(e2 > 0.0, e2, np.nan))
    dx = np.diff(np.log1p(nodes))
    steps = 0.5 * (integrand[:, 1:] + integrand[:, :-1]) * dx
    dc = np.concatenate([np.zeros((Om.shape[0], 1)), np.cumsum(steps, axis=1)], axis=1)
    dc[bad] = np.nan
    return dc


def transverse_comoving(dc, Ok):
    """D_M from D_C: the distance that converts angles to sizes.

    In curved space the light rays spread differently: sinh for open space,
    sin for closed. In closed space D_M turns over once light has come more than a
    quarter of the way around the universe; past the antipode it would be
    negative, which no real observation can produce, so those are NaN.
    """
    dc = np.asarray(dc, float)
    Ok = np.asarray(Ok, float)
    if Ok.ndim == 1 and dc.ndim == 2:
        Ok = Ok[:, None]
    sk = np.sqrt(np.abs(Ok))
    with np.errstate(invalid="ignore", divide="ignore"):
        open_ = np.sinh(sk * dc) / np.where(sk > 0, sk, 1.0)
        closed = np.sin(sk * dc) / np.where(sk > 0, sk, 1.0)
        dm = np.where(Ok > 1e-12, open_, np.where(Ok < -1e-12, closed, dc))
        dm = np.where(dm > 0, dm, np.nan) if dm.ndim else dm
    return dm


def luminosity_distance_shape(z_hd, z_hel, Om, OL, w=-1.0, *, nodes=None, index=None):
    """Luminosity distance in units of c/H0, for a batch of universes.

    Uses the same redshift convention as the published Pantheon+ likelihood:
    the distance is integrated to ``z_hd`` (the CMB-frame redshift corrected for
    the host galaxy's peculiar motion -- the part of the redshift due to
    cosmic expansion), and one factor of (1 + z) uses the *heliocentric* redshift
    ``z_hel``, since that is the redshift that actually stretches the photons we
    receive. Returns shape (n_models, n_sn).
    """
    if nodes is None or index is None:
        nodes, index = integration_nodes(z_hd)
    Om_a = np.atleast_1d(np.asarray(Om, float))
    OL_a = np.atleast_1d(np.asarray(OL, float))
    dc = comoving_distance_nodes(nodes, Om_a, OL_a, w)[:, index]
    dm = transverse_comoving(dc, 1.0 - Om_a - OL_a)
    return (1.0 + np.asarray(z_hel, float))[None, :] * dm


def distance_modulus(d_l_mpc):
    """mu = 5 log10(d_L / 10 pc). A distance of 10 pc gives 0; 1 Mpc gives 25."""
    return 5.0 * np.log10(np.asarray(d_l_mpc, float)) + 25.0


def distance_modulus_curve(z, Om, OL, w=-1.0, H0=70.0, z_hel=None):
    """mu(z) for one universe, on an arbitrary sorted-or-not grid of redshifts."""
    z = np.asarray(z, float)
    z_hel = z if z_hel is None else np.asarray(z_hel, float)
    shape = luminosity_distance_shape(z, z_hel, Om, OL, w)[0]
    return distance_modulus(shape * hubble_distance_mpc(H0))


# --- Closed forms, for the textbook universes and for testing ------------------


def dl_milne(z):
    """Empty universe (Om = OL = 0, Ok = 1): d_L = z (1 + z/2), in c/H0 units.

    Nothing in it to speed up or slow down the expansion: it coasts. This is the
    dividing line between acceleration and deceleration on a Hubble diagram.
    """
    z = np.asarray(z, float)
    return z * (1.0 + 0.5 * z)


def dl_einstein_de_sitter(z):
    """Flat, matter-only universe (Om = 1): d_L = 2 (1+z) (1 - 1/sqrt(1+z))."""
    z = np.asarray(z, float)
    return 2.0 * (1.0 + z) * (1.0 - 1.0 / np.sqrt(1.0 + z))


def dl_de_sitter(z):
    """Flat, Lambda-only universe (OL = 1): d_L = z (1 + z). Pure acceleration."""
    z = np.asarray(z, float)
    return z * (1.0 + z)


# --- Dynamics: deceleration, transition, age -----------------------------------


def deceleration_parameter(z, Om, OL, w=-1.0):
    """q(z) = -a (d^2 a/dt^2) / (da/dt)^2. Negative means accelerating.

    Matter pulls (+Om/2 per unit), a component with w < -1/3 pushes. Curvature
    affects E(z) but exerts no force, so it is absent from the numerator.
    """
    zp1 = 1.0 + np.asarray(z, float)
    num = 0.5 * Om * zp1**3 + 0.5 * (1.0 + 3.0 * w) * OL * zp1 ** (3.0 * (1.0 + w))
    return num / e_squared(z, Om, OL, w)


def q0(Om, OL, w=-1.0):
    """Today's deceleration parameter: Om/2 + (1 + 3w) OL / 2 (= Om/2 - OL for Lambda)."""
    return 0.5 * np.asarray(Om, float) + 0.5 * (1.0 + 3.0 * np.asarray(w, float)) * np.asarray(OL, float)


def transition_redshift(Om, OL):
    """Redshift where a Lambda universe switched from slowing down to speeding up.

    Setting q = 0: Om (1+z)^3 / 2 = OL, so z_t = (2 OL / Om)^(1/3) - 1.
    Valid for any curvature (curvature exerts no force). NaN if never accelerating.
    """
    Om, OL = np.asarray(Om, float), np.asarray(OL, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        zt = np.cbrt(2.0 * OL / Om) - 1.0
    return np.where((OL > 0) & (Om > 0), zt, np.nan)


def _one_over_aE(a, Om, OL, w):
    """1 / (a E(a)): dt = da / (a H). Written in a for use in time integrals."""
    Ok = 1.0 - Om - OL
    e2 = Om * a**-3 + Ok * a**-2 + OL * a ** (-3.0 * (1.0 + w))
    return 1.0 / (a * np.sqrt(e2))


def cosmic_time_of_a(a, Om, OL, w=-1.0, H0=70.0, n: int = 20000) -> np.ndarray:
    """Time since the Big Bang at scale factor a (a = 1 today, a = 1/(1+z)), in Gyr.

    Integrates dt = da / (a H(a)) with the substitution a = u^2, which removes
    the square-root behaviour near a = 0 so a simple trapezoid rule converges.
    ``a`` may extend beyond 1 (the future). Assumes a Big Bang exists.
    """
    a = np.asarray(a, float)
    u = np.linspace(0.0, np.sqrt(a.max()), n)
    uu = u[1:]
    integrand = np.concatenate([[0.0], 2.0 * uu * _one_over_aE(uu**2, Om, OL, w)])
    t = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(u))])
    return np.interp(np.sqrt(a), u, t) * hubble_time_gyr(H0)


def age_gyr(Om, OL, w=-1.0, H0=70.0) -> float:
    """Age of the universe today, in Gyr."""
    return float(cosmic_time_of_a(np.array([1.0]), Om, OL, w, H0)[0])


def lookback_time_gyr(z, Om, OL, w=-1.0, H0=70.0):
    """How long ago light now arriving from redshift z was emitted, in Gyr."""
    a = 1.0 / (1.0 + np.asarray(z, float))
    t = cosmic_time_of_a(np.concatenate([np.atleast_1d(a), [1.0]]), Om, OL, w, H0)
    return t[-1] - t[:-1] if np.ndim(z) else float(t[-1] - t[0])


# --- The geography of the (Om, OL) plane --------------------------------------


def no_big_bang_boundary(Om):
    """Omega_Lambda above which a universe has no Big Bang (it "bounced").

    With enough vacuum energy and positive curvature, running the expansion
    backwards reaches a minimum size and re-expands: there was no singular
    beginning. Carroll, Press & Turner (1992, ARA&A 30, 499):

        OL = 4 Om { f[ (1/3) f^-1( (1 - Om) / Om ) ] }^3,
        f = cosh for Om < 1/2,  f = cos for Om > 1/2.

    Models on or above this line are excluded as soon as we see any redshift
    beyond the bounce.
    """
    Om = np.asarray(Om, float)
    x = (1.0 - Om) / Om
    with np.errstate(invalid="ignore"):
        low = 4.0 * Om * np.cosh(np.arccosh(np.maximum(x, 1.0)) / 3.0) ** 3
        high = 4.0 * Om * np.cos(np.arccos(np.clip(x, -1.0, 1.0)) / 3.0) ** 3
    return np.where(Om < 0.5, low, high)


def bounces(Om, OL) -> np.ndarray:
    """True where a (Om, OL) universe has no Big Bang, by direct analysis.

    Multiply E(a)^2 by a^3: g(a) = Om + Ok a + OL a^3 must stay positive for all
    0 < a <= 1. g(0) = Om >= 0 and g(1) = 1; its only interior minimum is at
    a* = sqrt(-Ok / (3 OL)), where g(a*) = Om + (2/3) Ok a*. Independent of the
    closed-form boundary above, so each checks the other.
    """
    Om, OL = np.broadcast_arrays(np.asarray(Om, float), np.asarray(OL, float))
    Ok = 1.0 - Om - OL
    with np.errstate(invalid="ignore", divide="ignore"):
        a_star = np.sqrt(-Ok / (3.0 * OL))
        g_min = Om + (2.0 / 3.0) * Ok * a_star
    return (Ok < 0) & (OL > 0) & (a_star < 1.0) & (g_min <= 0.0)
