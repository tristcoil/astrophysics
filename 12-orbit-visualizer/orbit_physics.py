"""Orbital mechanics: from a position and a velocity to everything else.

Horizons gives, for each body at each instant, a position ``r`` and a velocity
``v`` relative to the Sun. Newton's two-body problem says that pair alone fixes
the whole orbit, because two numbers are conserved along it:

* the **specific energy** ``E = v**2/2 - GM/r``. Negative means bound (an
  ellipse); zero is the parabolic borderline; positive means the body escapes on
  a hyperbola and keeps a leftover speed ``v_inf = sqrt(2E)`` at infinity.
* the **specific angular momentum** ``h = r x v``. Its direction is the orbit's
  pole; its size is twice the area swept per unit time -- Kepler's second law
  is literally the statement that ``|h|`` is constant.

From those follow the orbit's shape (eccentricity ``e``), size (semi-major axis
``a = -GM / 2E``), orientation, and -- for escaping bodies -- the directions of
arrival and departure.

Real orbits are not exactly two-body: every planet tugs on every other. So these
functions describe the *osculating* orbit: the conic the body would follow from
this instant if everything but the Sun vanished. Watching how much the
osculating quantities wander along a real trajectory is how we see those tugs.

Units, everywhere unless a name says otherwise: au, days, au/day, radians, and
GM in au**3/day**2. Conversions to km/s live here, once.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

#: The astronomical unit, exact by IAU 2012 definition.
AU_KM = 149_597_870.7
DAY_S = 86_400.0
#: One au/day expressed in km/s (~1731.46).
KMS_PER_AU_DAY = AU_KM / DAY_S
#: Speed of light, exact by definition of the metre.
C_KMS = 299_792.458
#: Speed of light in au/day -- equivalently, one light-day in au (~173.14).
C_AU_DAY = C_KMS / KMS_PER_AU_DAY
LIGHT_DAY_AU = C_AU_DAY

#: The Sun's GM from JPL's DE440 ephemeris (Park et al. 2021), km^3/s^2.
#: Known to 11 significant figures -- far better than G or the Sun's mass
#: separately, which is why orbital mechanics is always done with GM.
GM_SUN_KM3_S2 = 132_712_440_041.279419
GM_SUN = GM_SUN_KM3_S2 * DAY_S**2 / AU_KM**3

#: Obliquity of the ecliptic at J2000 (IAU 1976), the tilt between Horizons'
#: ecliptic frame and the ICRF equator. Horizons documents this exact value.
OBLIQUITY_J2000 = np.radians(84381.448 / 3600.0)


# --------------------------------------------------------------------------
# Vector helpers (all accept a single (3,) vector or an (N, 3) stack)
# --------------------------------------------------------------------------

def dot(a, b):
    return np.sum(np.asarray(a) * np.asarray(b), axis=-1)


def norm(a):
    return np.sqrt(dot(a, a))


def unit(a):
    a = np.asarray(a, dtype=float)
    return a / norm(a)[..., None]


def to_kms(v_au_day):
    """Speed or velocity in au/day -> km/s."""
    return np.asarray(v_au_day) * KMS_PER_AU_DAY


# --------------------------------------------------------------------------
# Conserved quantities and the osculating orbit
# --------------------------------------------------------------------------

def specific_energy(r, v, mu=GM_SUN):
    """Orbital energy per unit mass, au^2/day^2. Negative = bound."""
    return 0.5 * dot(v, v) - mu / norm(r)


def angular_momentum(r, v):
    """Specific angular momentum vector ``r x v``, au^2/day."""
    return np.cross(r, v)


def areal_velocity(r, v):
    """Area swept by the Sun-body line per unit time, au^2/day: ``|r x v| / 2``.

    Kepler's second law says this is constant. It is constant exactly when the
    force on the body points straight at the Sun (a central force has no
    torque), so any drift in it measures the sideways pull of everything else.
    """
    return 0.5 * norm(np.cross(r, v))


def eccentricity_vector(r, v, mu=GM_SUN):
    """Points from the Sun towards perihelion; its length is the eccentricity."""
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    rn = norm(r)[..., None]
    rv = dot(r, v)[..., None]
    vv = dot(v, v)[..., None]
    return ((vv - mu / rn) * r - rv * v) / mu


@dataclass
class Elements:
    """An osculating orbit. Angles in radians, distances in au.

    ``a`` is negative for hyperbolic (unbound) orbits, by the usual convention
    that keeps ``E = -GM / 2a`` valid for every conic.
    """

    a: np.ndarray
    e: np.ndarray
    inc: np.ndarray
    raan: np.ndarray  # longitude of ascending node, Omega
    argp: np.ndarray  # argument of perihelion, omega
    nu: np.ndarray  # true anomaly
    q: np.ndarray  # perihelion distance
    energy: np.ndarray
    mu: float

    @property
    def bound(self):
        return self.energy < 0

    @property
    def period(self):
        """Orbital period in days (NaN for unbound orbits)."""
        a = np.where(self.bound, self.a, np.nan)
        return 2 * np.pi * np.sqrt(a**3 / self.mu)

    @property
    def aphelion(self):
        return np.where(self.bound, self.a * (1 + self.e), np.inf)

    @property
    def v_inf(self):
        """Speed at infinity in au/day (NaN for bound orbits)."""
        return np.sqrt(np.where(self.energy > 0, 2 * self.energy, np.nan))


def elements_from_state(r, v, mu=GM_SUN) -> Elements:
    """Classical orbital elements of the osculating orbit.

    The angles use ``arctan2`` throughout rather than ``arccos`` plus a sign
    fix, which keeps them well defined all the way round the circle.
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    h = np.cross(r, v)
    hn = norm(h)
    hhat = h / hn[..., None]
    zhat = np.array([0.0, 0.0, 1.0])
    node = np.cross(zhat, h)  # points to the ascending node
    evec = eccentricity_vector(r, v, mu)
    e = norm(evec)
    energy = specific_energy(r, v, mu)

    inc = np.arctan2(norm(h[..., :2]), h[..., 2])
    raan = np.mod(np.arctan2(node[..., 1], node[..., 0]), 2 * np.pi)
    argp = np.mod(np.arctan2(dot(np.cross(node, evec), hhat), dot(node, evec)), 2 * np.pi)
    nu = np.mod(np.arctan2(dot(np.cross(evec, r), hhat), dot(evec, r)), 2 * np.pi)
    p = hn**2 / mu  # semi-latus rectum
    q = p / (1 + e)
    a = -mu / (2 * energy)
    return Elements(a=a, e=e, inc=inc, raan=raan, argp=argp, nu=nu, q=q, energy=energy, mu=mu)


def solve_kepler(mean_anomaly, e, tol=1e-14):
    """Eccentric anomaly E from Kepler's equation M = E - e sin E (ellipses).

    Kepler's equation turns "how long since perihelion" into "where on the
    ellipse". It has no closed-form solution; Newton's method converges in a
    handful of steps from the starting guess E = M (or pi for e > 0.8).
    """
    M = np.mod(np.asarray(mean_anomaly, dtype=float), 2 * np.pi)
    E = np.where(e > 0.8, np.pi, M)
    for _ in range(50):
        dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
        E = E - dE
        if np.all(np.abs(dE) < tol):
            break
    return E


def true_anomaly_from_eccentric(E, e):
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))


def _rotation_pqw_to_frame(inc, raan, argp):
    """Matrix taking perifocal (P towards perihelion, Q 90 deg ahead) to the frame."""
    cO, sO = np.cos(raan), np.sin(raan)
    ci, si = np.cos(inc), np.sin(inc)
    cw, sw = np.cos(argp), np.sin(argp)
    return np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si, cw * si, ci],
    ])


def state_from_elements(a, e, inc, raan, argp, nu, mu=GM_SUN):
    """Position and velocity from elements -- the inverse of the function above.

    Works for ellipses (a > 0, e < 1) and hyperbolas (a < 0, e > 1) alike,
    because both are written in terms of the semi-latus rectum
    ``p = a (1 - e**2)``, which is positive for every conic.
    """
    p = a * (1 - e**2)
    rmag = p / (1 + e * np.cos(nu))
    r_pqw = np.array([rmag * np.cos(nu), rmag * np.sin(nu), 0.0])
    v_pqw = np.sqrt(mu / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    R = _rotation_pqw_to_frame(inc, raan, argp)
    return R @ r_pqw, R @ v_pqw


# --------------------------------------------------------------------------
# Speeds
# --------------------------------------------------------------------------

def escape_speed(r, mu=GM_SUN):
    """Speed needed to escape to infinity from distance ``r``: sqrt(2GM/r), au/day."""
    return np.sqrt(2 * mu / np.asarray(r, dtype=float))


def circular_speed(r, mu=GM_SUN):
    return np.sqrt(mu / np.asarray(r, dtype=float))


def hyperbolic_excess_speed(r, v, mu=GM_SUN):
    """Speed left over at infinity, ``sqrt(v**2 - v_esc**2)``, au/day. NaN if bound."""
    energy = specific_energy(r, v, mu)
    return np.sqrt(np.where(energy > 0, 2 * energy, np.nan))


def kepler_period(a, mu=GM_SUN):
    """Kepler's third law in Newton's form: P = 2 pi sqrt(a^3 / GM), days."""
    return 2 * np.pi * np.sqrt(np.asarray(a, dtype=float) ** 3 / mu)


def gm_from_orbit(a, period):
    """Invert Kepler's third law: the central GM implied by an orbit's size and period.

    For a planet this is G(M_sun + m_planet): both bodies orbit their common
    centre of mass, so the relative orbit feels the sum of the two masses.
    """
    return 4 * np.pi**2 * np.asarray(a, dtype=float) ** 3 / np.asarray(period, dtype=float) ** 2


# --------------------------------------------------------------------------
# Why Kepler's third law is not exact: the other planets, smeared into rings
# --------------------------------------------------------------------------

def _ring_coefficients(n_terms=400):
    """c_n = [(2n)! / (4**n (n!)**2)]**2, the Legendre coefficients of a ring."""
    c = np.empty(n_terms)
    c[0] = 1.0
    for n in range(1, n_terms):
        c[n] = c[n - 1] * ((2 * n - 1) / (2 * n)) ** 2
    return c


def ring_gm_factor(radius_ratio):
    """Extra central pull from a ring of mass, in units of the ring's own GM.

    Averaged over its orbit, a planet's gravity is the same as that of a thin
    ring of its mass spread along the orbit. A body circling the Sun in the
    ring's plane at radius ``r`` feels, from a ring of radius ``b``:

    * ``r > b`` (ring inside): an extra *inward* pull. Far away the ring acts
      like its whole mass sitting at the centre (factor 1); closer in the pull
      is stronger than that, because the near side of the ring is closer.
    * ``r < b`` (ring outside): an *outward* pull towards the ring, which
      partly cancels the Sun's -- a negative factor, small unless the ring is
      close.

    Returns ``f`` such that the body's effective central GM gains ``f * GM_ring``.
    ``radius_ratio`` is ``b / r``. Uses the in-plane multipole series of the
    ring potential, ``Phi = -(G m / r_>) * sum c_n (r_< / r_>)**(2n)``,
    differentiated for the radial force; 400 terms converge for any ratio up to
    ~0.95.
    """
    x = float(radius_ratio)
    if not 0 <= x < 0.97 and not x > 1 / 0.97:
        raise ValueError("ring too close to the orbit for the series (ratio near 1)")
    c = _ring_coefficients()
    n = np.arange(len(c))
    if x < 1:  # ring inside the orbit
        return float(np.sum((2 * n + 1) * c * x ** (2 * n)))
    y = 1.0 / x  # ring outside: y = r / b < 1
    return float(-np.sum(2 * n * c * y ** (2 * n + 1)))


def effective_central_gm(index, a, gm_bodies, gm_sun=GM_SUN):
    """Kepler-law GM expected for body ``index`` once the others are included.

    ``a`` are the orbits' semi-major axes and ``gm_bodies`` the bodies' GMs.
    The body's own mass counts in full (the Sun and it orbit each other); each
    other body counts as a ring, via ``ring_gm_factor``. This is the
    orbit-averaged (secular) picture: good when periods are not close to a
    simple ratio, cruder near resonances -- Jupiter-Saturn's 5:2 near miss,
    Uranus-Neptune's near 2:1, and Pluto locked in 3:2 with Neptune.
    """
    total = gm_sun + gm_bodies[index]
    contributions = {}
    for j, (aj, gmj) in enumerate(zip(a, gm_bodies)):
        if j == index or gmj == 0:
            continue
        f = ring_gm_factor(aj / a[index])
        contributions[j] = f * gmj
        total += f * gmj
    return total, contributions


# --------------------------------------------------------------------------
# Unbound orbits: where they come from and where they go
# --------------------------------------------------------------------------

def asymptotes(r, v, mu=GM_SUN):
    """Arrival and departure directions of a hyperbolic orbit.

    Returns ``(radiant, heading, v_inf)``: unit vectors pointing to where the
    body came from (as seen from the Sun, its *radiant*) and where it is going,
    plus the speed at infinity in au/day.

    Along a hyperbola the true anomaly can never exceed
    ``nu_inf = arccos(-1/e)`` -- that is the direction of the asymptote. The
    body approaches along ``-nu_inf`` and leaves along ``+nu_inf``, measured
    from the perihelion direction in the orbital plane. (The asymptotes do not
    pass through the Sun -- they miss it by the impact parameter -- but on the
    sky only their direction matters.)
    """
    r = np.asarray(r, dtype=float)
    v = np.asarray(v, dtype=float)
    evec = eccentricity_vector(r, v, mu)
    e = norm(evec)
    if np.any(e <= 1):
        raise ValueError("asymptotes exist only for hyperbolic orbits (e > 1)")
    P = evec / e[..., None]
    hhat = unit(np.cross(r, v))
    Q = np.cross(hhat, P)
    nu_inf = np.arccos(-1.0 / e)[..., None]
    heading = np.cos(nu_inf) * P + np.sin(nu_inf) * Q
    radiant = np.cos(nu_inf) * P - np.sin(nu_inf) * Q
    return radiant, heading, hyperbolic_excess_speed(r, v, mu)


def flyby_turn_angle(v_inf, r_periapsis, mu):
    """How far a flyby bends a trajectory, in radians.

    For a hyperbola ``e = 1 + r_p v_inf**2 / GM`` and the asymptotes are
    ``2 arcsin(1/e)`` apart: slower or closer passes bend more. Consistent units
    required (e.g. km, km/s, km^3/s^2).
    """
    e = 1.0 + r_periapsis * v_inf**2 / mu
    return 2.0 * np.arcsin(1.0 / e)


def angle_between(a, b):
    """Angle between two vectors, radians (robust for tiny and large angles)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return np.arctan2(norm(np.cross(a, b)), dot(a, b))


# --------------------------------------------------------------------------
# Reference frames
# --------------------------------------------------------------------------

def ecliptic_to_icrs(vec):
    """Rotate J2000 ecliptic coordinates to the ICRS (equatorial) frame.

    A rotation about the shared x-axis (the vernal equinox direction) by the
    obliquity: the ecliptic pole tips over to RA 18h, Dec +66.56 deg.
    """
    vec = np.asarray(vec, dtype=float)
    c, s = np.cos(OBLIQUITY_J2000), np.sin(OBLIQUITY_J2000)
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    return np.stack([x, c * y - s * z, s * y + c * z], axis=-1)


def icrs_to_ecliptic(vec):
    vec = np.asarray(vec, dtype=float)
    c, s = np.cos(OBLIQUITY_J2000), np.sin(OBLIQUITY_J2000)
    x, y, z = vec[..., 0], vec[..., 1], vec[..., 2]
    return np.stack([x, c * y + s * z, -s * y + c * z], axis=-1)


@functools.lru_cache(maxsize=1)
def _icrs_to_galactic_matrix():
    """ICRS -> Galactic rotation, taken from astropy rather than hand-copied.

    Built by transforming the three ICRS axis directions; the columns of the
    matrix are their Galactic Cartesian images.
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    axes = SkyCoord(ra=[0.0, 90.0, 0.0] * u.deg, dec=[0.0, 0.0, 90.0] * u.deg, frame="icrs")
    gal = axes.galactic
    lon, lat = gal.l.radian, gal.b.radian
    cols = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    return cols  # rows: galactic x, y, z; columns: images of ICRS x, y, z


def icrs_to_galactic(vec):
    """Galactic Cartesian: x to the Galactic centre, y along rotation, z to the NGP."""
    return np.asarray(vec, dtype=float) @ _icrs_to_galactic_matrix().T


def spherical_deg(vec):
    """(longitude, latitude) in degrees of a direction, longitude in [0, 360)."""
    vec = np.asarray(vec, dtype=float)
    lon = np.degrees(np.arctan2(vec[..., 1], vec[..., 0])) % 360.0
    lat = np.degrees(np.arctan2(vec[..., 2], np.hypot(vec[..., 0], vec[..., 1])))
    return lon, lat


def radec_from_ecliptic(vec):
    """RA and Dec (degrees, ICRS) of a direction given in ecliptic J2000 coordinates."""
    return spherical_deg(ecliptic_to_icrs(vec))


def galactic_uvw_kms(v_ecliptic_au_day):
    """A velocity in ecliptic J2000 au/day -> Galactic (U, V, W) in km/s.

    U points towards the Galactic centre, V along the direction of Galactic
    rotation, W towards the north Galactic pole -- the convention of the stellar
    kinematics literature.
    """
    return icrs_to_galactic(ecliptic_to_icrs(to_kms(v_ecliptic_au_day)))


# --------------------------------------------------------------------------
# Interpolation and event finding
# --------------------------------------------------------------------------

def hermite(s, f0, f1, d0, d1, h):
    """Cubic Hermite interpolation on one interval.

    ``s`` in [0, 1] is the fractional position, ``f`` the values at the ends,
    ``d`` their time derivatives, ``h`` the interval length. With exact
    derivatives (Horizons gives velocities alongside positions) this is
    fourth-order accurate: halving the step cuts the error sixteenfold.
    """
    s2, s3 = s * s, s * s * s
    return ((2 * s3 - 3 * s2 + 1) * f0 + (s3 - 2 * s2 + s) * h * d0
            + (-2 * s3 + 3 * s2) * f1 + (s3 - s2) * h * d1)


def hermite_slope(s, f0, f1, d0, d1, h):
    s2 = s * s
    return ((6 * s2 - 6 * s) * f0 + (3 * s2 - 4 * s + 1) * h * d0
            + (-6 * s2 + 6 * s) * f1 + (3 * s2 - 2 * s) * h * d1) / h


class Track:
    """A sampled trajectory that can be evaluated at any time in its span."""

    def __init__(self, jd, r, v):
        self.jd = np.asarray(jd, dtype=float)
        self.r = np.asarray(r, dtype=float)
        self.v = np.asarray(v, dtype=float)
        if np.any(np.diff(self.jd) <= 0):
            raise ValueError("track times must increase")

    def _locate(self, t):
        t = np.atleast_1d(np.asarray(t, dtype=float))
        if np.any(t < self.jd[0]) or np.any(t > self.jd[-1]):
            raise ValueError("time outside the track's span")
        i = np.clip(np.searchsorted(self.jd, t, side="right") - 1, 0, len(self.jd) - 2)
        h = self.jd[i + 1] - self.jd[i]
        s = (t - self.jd[i]) / h
        return i, s, h

    def position(self, t):
        i, s, h = self._locate(t)
        return hermite(s[:, None], self.r[i], self.r[i + 1], self.v[i], self.v[i + 1],
                       h[:, None])

    def velocity(self, t):
        # Velocity from linear interpolation of the sampled velocities is
        # plenty for our uses; positions are what need the higher order.
        i, s, h = self._locate(t)
        return self.v[i] * (1 - s[:, None]) + self.v[i + 1] * s[:, None]


def _solve_on_interval(fun, lo=0.0, hi=1.0, tol=1e-13, max_iter=200):
    """Bisection for a sign change of ``fun`` on [lo, hi]."""
    flo = fun(lo)
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        fm = fun(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def crossings(t, f, fdot, level):
    """Times where a sampled function crosses ``level``, with direction (+1 up, -1 down).

    Located by Hermite interpolation using the known derivatives, so the timing
    is limited by the ephemeris, not by the sampling step.
    """
    t = np.asarray(t, dtype=float)
    g = np.asarray(f, dtype=float) - level
    out = []
    for i in np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]:
        h = t[i + 1] - t[i]
        s = _solve_on_interval(lambda s_: hermite(s_, g[i], g[i + 1], fdot[i], fdot[i + 1], h))
        out.append((t[i] + s * h, 1 if g[i + 1] > g[i] else -1))
    return out


def extrema(t, f, fdot, kind="min"):
    """Local minima (or maxima) of a sampled function: [(time, value), ...].

    Found where the derivative changes sign, then refined by solving for the
    zero of the Hermite cubic's slope inside that interval.
    """
    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    fdot = np.asarray(fdot, dtype=float)
    want = (fdot[:-1] < 0) & (fdot[1:] >= 0) if kind == "min" else (fdot[:-1] > 0) & (fdot[1:] <= 0)
    out = []
    for i in np.nonzero(want)[0]:
        h = t[i + 1] - t[i]
        slope = lambda s_: hermite_slope(s_, f[i], f[i + 1], fdot[i], fdot[i + 1], h)
        s = _solve_on_interval(slope)
        out.append((t[i] + s * h, float(hermite(s, f[i], f[i + 1], fdot[i], fdot[i + 1], h))))
    return out


def radial_distance_and_rate(r, v):
    """|r| and d|r|/dt (the radial velocity) along a trajectory."""
    rn = norm(r)
    return rn, dot(r, v) / rn


# --------------------------------------------------------------------------
# Measuring an orbit from its trajectory
# --------------------------------------------------------------------------

def orbital_plane(r, v):
    """Unit normal of the mean orbital plane: the average direction of r x v."""
    return unit(np.mean(unit(np.cross(r, v)), axis=0))


def in_plane_basis(normal):
    """Two orthonormal in-plane axes; the first is the ecliptic x-axis projected in."""
    x = np.array([1.0, 0.0, 0.0])
    ex = unit(x - dot(x, normal) * normal)
    ey = np.cross(normal, ex)
    return ex, ey


def in_plane_angle(r, v, normal):
    """Unwrapped polar angle in the orbital plane, and its time derivative.

    The angle always increases for motion in the sense of ``normal``, which is
    how the orbit's own pole is chosen -- so retrograde orbits need no special
    case.
    """
    ex, ey = in_plane_basis(normal)
    x, y = dot(r, ex), dot(r, ey)
    vx, vy = dot(v, ex), dot(v, ey)
    theta = np.unwrap(np.arctan2(y, x))
    theta_dot = (x * vy - y * vx) / (x * x + y * y)
    return theta, theta_dot


def fit_conic(r, normal):
    """Least-squares conic with the Sun at its focus: 1/r = A + B cos(th) + C sin(th).

    Kepler's first law says an orbit is a conic section with the Sun at one
    focus, ``r = p / (1 + e cos(th - varpi))``. Rewritten as ``1/r`` it is
    *linear* in three unknowns, so every sampled position can be used in one
    ordinary least-squares solve. Then ``p = 1/A``,
    ``e = sqrt(B**2 + C**2)/A``, ``a = p / (1 - e**2)`` and the perihelion
    direction is ``atan2(C, B)``.
    """
    ex, ey = in_plane_basis(normal)
    th = np.arctan2(dot(r, ey), dot(r, ex))
    rn = norm(r)
    design = np.column_stack([np.ones_like(th), np.cos(th), np.sin(th)])
    coef, *_ = np.linalg.lstsq(design, 1.0 / rn, rcond=None)
    A, B, C = coef
    p = 1.0 / A
    e = np.hypot(B, C) / A
    model = design @ coef
    return {"p": p, "e": e, "a": p / (1 - e**2), "varpi": np.arctan2(C, B),
            "rms_frac": float(np.sqrt(np.mean((rn * model - 1.0) ** 2)))}


def measure_orbit(jd, r, v):
    """Period and size of an orbit, measured purely from where the body was.

    * **Sidereal period** -- the time to return to the same direction as seen
      from the Sun: the moments the in-plane angle passes its starting value
      plus whole turns, averaged over every complete orbit in the data.
    * **Semi-major axis** -- a conic fitted, with the Sun at its focus, to all
      positions within those complete orbits.

    Per-orbit values are kept too: their scatter shows how much the other
    planets' pulls make each individual orbit differ from the average.
    """
    jd = np.asarray(jd, dtype=float)
    normal = orbital_plane(r, v)
    theta, theta_dot = in_plane_angle(r, v, normal)
    turns = int(np.floor((theta[-1] - theta[0]) / (2 * np.pi)))
    if turns < 1:
        raise ValueError("fewer than one complete orbit in the data")

    t_cross = [jd[0]]
    for k in range(1, turns + 1):
        level = theta[0] + 2 * np.pi * k
        hits = [c for c in crossings(jd, theta, theta_dot, level) if c[1] > 0]
        t_cross.append(hits[0][0])
    t_cross = np.array(t_cross)
    periods = np.diff(t_cross)

    whole = (jd >= t_cross[0]) & (jd <= t_cross[-1])
    fit = fit_conic(r[whole], normal)
    per_orbit_a = []
    for t0, t1 in zip(t_cross[:-1], t_cross[1:]):
        m = (jd >= t0) & (jd <= t1)
        per_orbit_a.append(fit_conic(r[m], normal)["a"])
    per_orbit_a = np.array(per_orbit_a)

    return {
        "n_orbits": turns,
        "period_days": float((t_cross[-1] - t_cross[0]) / turns),
        "period_scatter_days": float(np.std(periods, ddof=1)) if turns > 1 else None,
        "periods_days": periods,
        "a_au": float(fit["a"]),
        "a_scatter_au": float(np.std(per_orbit_a, ddof=1)) if turns > 1 else None,
        "a_per_orbit_au": per_orbit_a,
        "e": float(fit["e"]),
        "fit_rms_frac": fit["rms_frac"],
        "normal": normal,
        "span_days": float(t_cross[-1] - t_cross[0]),
    }
