"""The geometry and kinematics of a rotating disc seen from inside it.

Everything here is plain geometry plus the definitions of the velocity frames;
no survey data are read. The conventions (checked in ``test_physics.py``):

* Face-on coordinates: Galactic centre at the origin, the Sun at
  (x, y) = (0, R0), +x toward Galactic longitude 90 deg (the way the Sun is
  moving), +z toward the north Galactic pole. Seen from the north the Galaxy
  turns clockwise.
* Galactocentric azimuth ``beta`` (Reid et al. 2014, 2019): 0 toward the Sun,
  increasing in the direction of rotation, so beta = atan2(x, y).
* The survey's velocities are LSRK: heliocentric velocities plus the
  *standard solar motion* (20 km/s toward RA 18h, Dec +30 deg, B1900), a
  historical convention. The Sun's true motion is different, so we never
  assume, as the classic formulas do, that the LSRK frame moves on a circular
  orbit. Instead the observer moves with (U', V', W') = true solar motion minus
  the standard one, where V' includes the Sun's full orbital speed.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

import config as cfg

G_KPC = 4.30091727e-6  # gravitational constant in kpc (km/s)^2 / Msun


# --------------------------------------------------------------------------
# Velocity frames
# --------------------------------------------------------------------------

@lru_cache(maxsize=1)
def standard_solar_motion() -> np.ndarray:
    """(U, V, W) of the LSRK definition, km/s: 20 km/s toward (18h, +30d) B1900."""
    import astropy.units as u
    from astropy.coordinates import FK4, SkyCoord

    ra, dec = cfg.STANDARD_SOLAR_APEX_B1900
    apex = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=FK4(equinox="B1900")).galactic
    lr, br = apex.l.rad, apex.b.rad
    return cfg.STANDARD_SOLAR_MOTION_KMS * np.array(
        [np.cos(br) * np.cos(lr), np.cos(br) * np.sin(lr), np.sin(br)])


@dataclass(frozen=True)
class Observer:
    """Where the Sun is and how the LSRK frame moves relative to circular orbits.

    ``u, v, w`` are the LSRK frame's velocity in the Galactic rest frame, in the
    heliocentric U (toward the centre), V (direction of rotation) and W (north)
    directions. For the textbook assumption "the LSR moves on a circular orbit"
    these are (0, Theta0, 0).
    """

    r0: float
    u: float
    v: float
    w: float

    @classmethod
    def textbook(cls, r0: float, theta0: float) -> "Observer":
        return cls(r0, 0.0, theta0, 0.0)

    @classmethod
    def reid19(cls, r0: float = cfg.R0_KPC, omega_sun: float = cfg.OMEGA_SUN,
               u_sun: float = cfg.U_SUN_KMS, w_sun: float = cfg.W_SUN_KMS) -> "Observer":
        """The LSRK frame's true motion, from the Sun's measured motion.

        The Sun moves at Omega_sun * R0 in the direction of rotation (the reflex
        of Sgr A*'s apparent motion), plus U_sun and W_sun; the LSRK frame is
        the Sun minus the standard solar motion.
        """
        std = standard_solar_motion()
        return cls(r0, u_sun - std[0], omega_sun * r0 - std[1], w_sun - std[2])

    def frame_velocity_along(self, l_deg, b_deg=0.0):
        """Line-of-sight component of the LSRK frame's own Galactic motion."""
        l, b = np.radians(l_deg), np.radians(b_deg)
        return (self.u * np.cos(l) + self.v * np.sin(l)) * np.cos(b) + self.w * np.sin(b)


# --------------------------------------------------------------------------
# Rotation curves
# --------------------------------------------------------------------------

def theta_flat(r, theta0: float = cfg.THETA0_KMS):
    return np.full_like(np.asarray(r, dtype=float), theta0)


def theta_reid19(r, r0: float = cfg.R0_KPC, a2: float = 0.96, a3: float = 1.62):
    """Reid et al. (2019) fit A5: the Persic, Salucci & Stel (1996) "universal"
    rotation curve in its two-parameter form, transcribed from their Appendix B
    FORTRAN routine ``Univ_RC_from_note``."""
    r = np.asarray(r, dtype=float)
    lam = (a3 / 1.5) ** 5
    rho = r / (a2 * r0)
    log_lam = np.log10(lam)
    term1 = 200.0 * lam ** 0.41
    term2 = np.sqrt(0.80 + 0.49 * log_lam
                    + 0.75 * np.exp(-0.4 * lam) / (0.47 + 2.25 * lam ** 0.4))
    term3 = (0.72 + 0.44 * log_lam) * 1.97 * rho ** 1.22 / (rho ** 2 + 0.61) ** 1.43
    term4 = 1.6 * np.exp(-0.4 * lam) * rho ** 2 / (rho ** 2 + 2.25 * lam ** 0.4)
    return term1 / term2 * np.sqrt(term3 + term4)


#: Clemens (1985, ApJ 295, 422) CO rotation curve for R0 = 8.5 kpc and
#: Theta0 = 220 km/s: a piecewise polynomial in R (kpc), coefficients highest
#: power first, with breaks at 0.09, 0.45 and 1.6 R0.
_CLEMENS = [
    (0.09, [-17731.0, 54904.0, -68287.3, 43980.1, -15809.8, 3069.81, 0.0]),
    (0.45, [-2.110625, 25.073006, -110.73531, 231.87099, -248.1467, 325.0912]),
    (1.60, [0.00129348, -0.08050808, 2.0697271, -28.4080026, 224.562732,
            -1024.068760, 2507.60391, -2342.6564]),
    (np.inf, [234.88]),
]


def theta_clemens85(r):
    r = np.asarray(r, dtype=float)
    out = np.empty_like(r)
    lo = 0.0
    for hi, coeffs in _CLEMENS:
        sel = (r / cfg.R0_IAU_KPC >= lo) & (r / cfg.R0_IAU_KPC < hi)
        out[sel] = np.polyval(coeffs, r[sel])
        lo = hi
    return out


def theta_mgd(r, *, joint2016: bool = True):
    """McClure-Griffiths & Dickey linear fits, Theta/Theta0 = k0 R/R0 + k1, for
    R0 = 8.5 kpc, Theta0 = 220 km/s, valid for 3 < R < 8 kpc: the 2016 joint
    first + fourth quadrant fit (0.171, 0.889) or the 2007 one (0.186, 0.887)."""
    k0, k1 = (0.171, 0.889) if joint2016 else (0.186, 0.887)
    return cfg.THETA0_IAU_KMS * (k0 * np.asarray(r) / cfg.R0_IAU_KPC + k1)


# --------------------------------------------------------------------------
# Geometry of a sightline
# --------------------------------------------------------------------------

def xy_from_ld(l_deg, d, r0: float):
    """Face-on position of a point at longitude l and in-plane distance d."""
    l = np.radians(l_deg)
    return d * np.sin(l), r0 - d * np.cos(l)


def ld_from_xy(x, y, r0: float):
    """Longitude (deg, 0-360) and in-plane distance of a face-on position."""
    dx, dy = x, r0 - y
    return np.degrees(np.arctan2(dx, dy)) % 360.0, np.hypot(dx, dy)


def radius_from_ld(l_deg, d, r0: float):
    return np.sqrt(r0 ** 2 + d ** 2 - 2 * r0 * d * np.cos(np.radians(l_deg)))


def v_lsr(l_deg, b_deg, d, theta_fn, obs: Observer):
    """LSRK velocity of gas in circular rotation at in-plane distance d."""
    r = radius_from_ld(l_deg, d, obs.r0)
    l, b = np.radians(l_deg), np.radians(b_deg)
    v_gas = theta_fn(r) * obs.r0 / r * np.sin(l) * np.cos(b)
    return v_gas - obs.frame_velocity_along(l_deg, b_deg)


def theta_from_terminal(l_deg, v_term, obs: Observer):
    """Circular speed at the tangent point R = R0 |sin l| from a terminal velocity.

    At the tangent point the orbit runs along the line of sight, so the gas
    velocity projects fully: v_term = sign(sin l) Theta(R_t) - (frame motion).
    """
    s = np.sign(np.sin(np.radians(l_deg)))
    return s * (np.asarray(v_term) + obs.frame_velocity_along(l_deg))


def tangent_radius(l_deg, r0: float):
    return r0 * np.abs(np.sin(np.radians(l_deg)))


def outer_distance(l_deg, b_deg, v, theta_fn, obs: Observer, r_max: float = 60.0):
    """Kinematic distance of gas *outside* the solar circle (R > R0).

    Outside the solar circle every radius is met once along a sightline, so a
    velocity maps to a single distance. Returns (R, d) with NaN where the
    velocity is not one an outer-disc orbit can produce on this sightline.
    """
    l_deg, b_deg, v = np.broadcast_arrays(np.asarray(l_deg, float),
                                          np.asarray(b_deg, float),
                                          np.asarray(v, float))
    l, b = np.radians(l_deg), np.radians(b_deg)
    vf = obs.frame_velocity_along(l_deg, b_deg)
    # w(R) = Theta(R) R0 / R, the angular-speed term; must equal (v + vf)/(sin l cos b)
    with np.errstate(divide="ignore", invalid="ignore"):
        target = (v + vf) / (np.sin(l) * np.cos(b))
    r_grid = np.linspace(obs.r0, r_max, 4000)
    w_grid = theta_fn(r_grid) * obs.r0 / r_grid  # decreasing for any sane curve
    if np.any(np.diff(w_grid) >= 0):
        raise ValueError("Theta(R) R0 / R must decrease outward for unique distances")
    r = np.interp(-target, -w_grid, r_grid, left=np.nan, right=np.nan)
    d = obs.r0 * np.cos(l) + np.sqrt(np.clip(r ** 2 - (obs.r0 * np.sin(l)) ** 2, 0, None))
    return r, d


# --------------------------------------------------------------------------
# Spiral arms (Reid et al. 2019, Table 2)
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Arm:
    """A log-periodic spiral with a kink: ln(R/R_kink) = -(beta - beta_kink) tan(psi)."""

    name: str
    beta_kink: float   # deg
    r_kink: float      # kpc
    psi_lt: float      # pitch angle for beta <= beta_kink, deg
    psi_gt: float      # pitch angle for beta > beta_kink, deg
    beta_data: tuple[float, float]  # azimuth range with parallax data, deg
    beta_draw: tuple[float, float]  # range we draw (dashed outside the data)
    #: (beta_below, psi) : below this azimuth, Reid et al. switch to another
    #: pitch angle to follow the arm's CO trace (Sgr-Car beyond its tangency).
    outer_pitch: tuple[float, float] | None = None

    def _kinked(self, beta):
        psi = np.where(beta <= self.beta_kink, self.psi_lt, self.psi_gt)
        return self.r_kink * np.exp(-np.radians(beta - self.beta_kink) * np.tan(np.radians(psi)))

    def radius(self, beta_deg):
        beta = np.asarray(beta_deg, dtype=float)
        r = self._kinked(beta)
        if self.outer_pitch is not None:
            b1, psi1 = self.outer_pitch
            r_b1 = self._kinked(np.asarray(b1, dtype=float))
            r = np.where(beta < b1,
                         r_b1 * np.exp(-np.radians(beta - b1) * np.tan(np.radians(psi1))), r)
        return r

    def xy(self, beta_deg):
        r = self.radius(beta_deg)
        beta = np.radians(beta_deg)
        return r * np.sin(beta), r * np.cos(beta)


#: Reid et al. (2019) Table 2, with R0 = 8.15 kpc. The 3-kpc arms are left out:
#: they are not in circular rotation (the near arm expands toward us at
#: 53 km/s), so a circular-orbit projection would put them in the wrong place in
#: velocity; we draw them from their measured (l, v) tracks instead.
REID19_ARMS = [
    Arm("Norma", 18, 4.46, -1.0, 19.5, (5, 54), (-60, 54)),
    Arm("Scutum-Centaurus", 23, 4.91, 14.1, 12.1, (0, 104), (-55, 104)),
    Arm("Sagittarius-Carina", 24, 6.04, 17.1, 1.0, (-2, 97), (-45, 97),
        outer_pitch=(-33.0, 10.0)),
    Arm("Local", 9, 8.26, 11.4, 11.4, (-8, 34), (-20, 45)),
    Arm("Perseus", 40, 8.87, 10.3, 8.7, (-23, 115), (-40, 130)),
    Arm("Outer", 18, 12.24, 3.0, 9.4, (-16, 71), (-30, 85)),
]

#: The Scutum-Centaurus model continued once around the far side of the
#: Galaxy: the "Outer Scutum-Centaurus" arm that Dame & Thaddeus (2011) found
#: in CO in the first quadrant, 20 kpc away (one maser so far). Azimuth range.
OSC_BETA = (-250.0, -205.0)

#: Dame & Thaddeus (2008, ApJ 683, L143): the expanding 3-kpc arms in the CO
#: longitude-velocity diagram, v = v0 + slope * l (l in degrees, -12 to +8).
THREE_KPC_ARMS = {"Near 3-kpc": (-53.1, 4.16), "Far 3-kpc": (56.0, 4.08)}
THREE_KPC_L_RANGE = (-12.0, 8.0)


def arm_lv(arm: Arm, theta_fn, obs: Observer, n: int = 1500):
    """The arm's track in the longitude-velocity diagram (circular rotation)."""
    beta = np.linspace(*arm.beta_draw, n)
    x, y = arm.xy(beta)
    l, d = ld_from_xy(x, y, obs.r0)
    v = v_lsr(l, 0.0, d, theta_fn, obs)
    in_data = (beta >= arm.beta_data[0]) & (beta <= arm.beta_data[1])
    return l, v, d, in_data


# --------------------------------------------------------------------------
# Mass
# --------------------------------------------------------------------------

def enclosed_mass(r_kpc, v_kms):
    """Mass inside radius r needed to hold a circular orbit of speed v, spherical
    approximation (Msun). For a flattened mass distribution the true mass is
    somewhat smaller, by up to ~30% for a pure exponential disc."""
    return np.asarray(v_kms) ** 2 * np.asarray(r_kpc) / G_KPC
