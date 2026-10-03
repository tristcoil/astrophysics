"""Coordinate frames, Galactocentric velocities and gravitational potentials.

Everything the three parts of the project share: turning Gaia's observables
(sky position, parallax, proper motion, radial velocity) into positions and
velocities in the Galaxy, and turning those into orbital energies and angular
momenta in a model of the Milky Way's gravity.

Conventions (the same as the papers we compare with):

* Galactocentric Cartesian x, y, z follow astropy: the Galactic Centre at the
  origin, the Sun near x = -R0, y pointing along the direction of rotation at
  the Sun, z towards the North Galactic Pole.
* Cylindrical velocities: v_R positive outwards, v_phi **positive in the
  direction of Galactic rotation** (the disc has v_phi ~ +230 km/s), v_z
  positive northwards.
* Angular momentum L_z is likewise **positive for prograde orbits** (the
  disc), negative for retrograde ones, as in Helmi et al. (2018) and
  Koppelman et al. (2019). Units: kpc km/s.

Every transform is plain numpy (fast for millions of stars) and is checked
against astropy and gala in test_physics.py.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

import numpy as np

from config import (
    AUX_DIR, G_KPC_KMS2_MSUN, GD1_ROTATION, HELMI17, K_KMS, R0_KPC, V_PECULIAR_KMS,
    V_SUN_KMS, Z_SUN_KPC, Helmi17Potential,
)

# --------------------------------------------------------------------------
# Rotation matrices
# --------------------------------------------------------------------------

#: ICRS -> Galactic (astropy's matrix to 16 digits; rows are the Galactic x, y,
#: z axes expressed in ICRS).
ICRS_TO_GAL = np.array([
    [-0.0548756577125916, -0.8734370519556159, -0.4838350736167155],
    [0.4941094371927268, -0.4448297212232952, 0.7469821839866676],
    [-0.8676661375596576, -0.1980763372730005, 0.4559838136873016],
])

#: ICRS -> GD-1 stream frame (Koposov et al. 2010, Appendix).
ICRS_TO_GD1 = np.array(GD1_ROTATION)


def unit_vector(lon_deg, lat_deg):
    """Cartesian unit vectors (3, N) for longitudes/latitudes in degrees."""
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    return np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def tangent_basis(lon_deg, lat_deg):
    """Unit vectors along increasing longitude and latitude, each (3, N)."""
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    e_lon = np.array([-np.sin(lon), np.cos(lon), np.zeros_like(lon)])
    e_lat = np.array([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    return e_lon, e_lat


def rotate_sky(matrix, lon_deg, lat_deg, pm_lon=None, pm_lat=None):
    """Rotate sky positions (and optionally proper motions) into another frame.

    Proper motions are (mu_lon * cos lat, mu_lat) in any consistent unit.
    Returns lon, lat in degrees (lon wrapped to [-180, 180)), and, if proper
    motions were given, the rotated pair plus the 2x2 per-star rotation
    (for transforming covariance matrices).
    """
    v = matrix @ unit_vector(lon_deg, lat_deg)
    lon2 = np.degrees(np.arctan2(v[1], v[0]))
    lat2 = np.degrees(np.arcsin(np.clip(v[2], -1, 1)))
    if pm_lon is None:
        return lon2, lat2
    e_lon, e_lat = tangent_basis(lon_deg, lat_deg)
    f_lon, f_lat = tangent_basis(lon2, lat2)
    # Old tangent-plane basis vectors, expressed in the new frame.
    e_lon_new, e_lat_new = matrix @ e_lon, matrix @ e_lat
    # 2x2 rotation between the two tangent planes (orthonormal).
    m11 = np.sum(f_lon * e_lon_new, axis=0)
    m12 = np.sum(f_lon * e_lat_new, axis=0)
    m21 = np.sum(f_lat * e_lon_new, axis=0)
    m22 = np.sum(f_lat * e_lat_new, axis=0)
    pm_lon2 = m11 * pm_lon + m12 * pm_lat
    pm_lat2 = m21 * pm_lon + m22 * pm_lat
    return lon2, lat2, pm_lon2, pm_lat2, (m11, m12, m21, m22)


def rotate_pm_covariance(m, s_lon, s_lat, corr):
    """Propagate a proper-motion covariance through the 2x2 rotation ``m``.

    Returns the new errors and correlation coefficient."""
    m11, m12, m21, m22 = m
    c11, c22, c12 = s_lon**2, s_lat**2, corr * s_lon * s_lat
    n11 = m11**2 * c11 + 2 * m11 * m12 * c12 + m12**2 * c22
    n22 = m21**2 * c11 + 2 * m21 * m22 * c12 + m22**2 * c22
    n12 = m11 * m21 * c11 + (m11 * m22 + m12 * m21) * c12 + m12 * m22 * c22
    return np.sqrt(n11), np.sqrt(n22), n12 / np.sqrt(n11 * n22)


def gd1_coordinates(ra, dec, pmra=None, pmdec=None):
    """ICRS -> GD-1 stream frame: (phi1, phi2[, pm_phi1_cosphi2, pm_phi2, m])."""
    return rotate_sky(ICRS_TO_GD1, ra, dec, pmra, pmdec)


def gd1_to_icrs(phi1, phi2):
    return tuple(np.mod(x, 360) if i == 0 else x
                 for i, x in enumerate(rotate_sky(ICRS_TO_GD1.T, phi1, phi2)))


def galactic_lb(ra, dec):
    l, b = rotate_sky(ICRS_TO_GAL, ra, dec)
    return np.mod(l, 360), b


# --------------------------------------------------------------------------
# ICRS -> Galactocentric
# --------------------------------------------------------------------------


@functools.lru_cache(maxsize=1)
def _galcen_affine():
    """The Galactocentric transform as x_gc = A @ x_icrs + b, v_gc = A @ v + v_sun.

    Built by pushing unit vectors through astropy's own Galactocentric frame
    once (so the geometry -- Sgr A*'s position, the roll angle, the tilt from
    the Sun's height -- is astropy's exactly), then applied with numpy.
    """
    import astropy.coordinates as coord
    import astropy.units as u

    frame = coord.Galactocentric(
        galcen_distance=R0_KPC * u.kpc, z_sun=Z_SUN_KPC * u.kpc,
        galcen_v_sun=coord.CartesianDifferential(list(V_SUN_KMS) * u.km / u.s))
    pts = np.array([[0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
    c = coord.SkyCoord(coord.CartesianRepresentation(pts * u.kpc), frame="icrs")
    g = c.transform_to(frame).cartesian.xyz.to_value(u.kpc)
    b = g[:, 0]
    A = g[:, 1:] - b[:, None]
    return A, b


@dataclass
class Galactocentric:
    """Positions (kpc) and velocities (km/s) of stars in the Galaxy."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray

    @property
    def R(self):
        return np.hypot(self.x, self.y)

    @property
    def r(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def phi_deg(self):
        """Galactocentric azimuth; the Sun sits at 180 degrees."""
        return np.degrees(np.arctan2(self.y, self.x)) % 360

    @property
    def v_R(self):
        return (self.x * self.vx + self.y * self.vy) / self.R

    @property
    def v_phi(self):
        """Azimuthal velocity, positive along Galactic rotation."""
        return (self.y * self.vx - self.x * self.vy) / self.R

    @property
    def L_z(self):
        """Vertical angular momentum, positive for prograde orbits (kpc km/s)."""
        return self.y * self.vx - self.x * self.vy

    @property
    def L_x(self):
        return self.y * self.vz - self.z * self.vy

    @property
    def L_y(self):
        return self.z * self.vx - self.x * self.vz

    @property
    def L_perp(self):
        """Angular momentum in the plane, sqrt(Lx^2 + Ly^2) (kpc km/s)."""
        return np.hypot(self.L_x, self.L_y)

    @property
    def speed(self):
        return np.sqrt(self.vx**2 + self.vy**2 + self.vz**2)

    def spherical_velocities(self):
        """(v_r, v_theta, v_phi): radial, polar (positive towards the south
        pole, i.e. increasing theta) and azimuthal (prograde positive)."""
        r, R = self.r, self.R
        v_r = (self.x * self.vx + self.y * self.vy + self.z * self.vz) / r
        v_theta = (self.z * (self.x * self.vx + self.y * self.vy) / R - R * self.vz) / r
        return v_r, v_theta, self.v_phi


def galactocentric(ra, dec, distance_kpc, pmra, pmdec, rv) -> Galactocentric:
    """Gaia observables -> Galactocentric positions and velocities."""
    A, b = _galcen_affine()
    u = unit_vector(ra, dec)
    e_ra, e_dec = tangent_basis(ra, dec)
    pos = A @ (u * distance_kpc) + b[:, None]
    vel_icrs = (u * rv + e_ra * (K_KMS * pmra * distance_kpc)
                + e_dec * (K_KMS * pmdec * distance_kpc))
    vel = A @ vel_icrs + np.asarray(V_SUN_KMS)[:, None]
    return Galactocentric(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])


def heliocentric_uvw(ra, dec, distance_kpc, pmra, pmdec, rv):
    """Heliocentric Galactic velocity (U, V, W), km/s: U towards the Galactic
    Centre, V along rotation, W to the North Galactic Pole."""
    u = unit_vector(ra, dec)
    e_ra, e_dec = tangent_basis(ra, dec)
    vel = (u * rv + e_ra * (K_KMS * pmra * distance_kpc) + e_dec * (K_KMS * pmdec * distance_kpc))
    return ICRS_TO_GAL @ vel


def speed_vs_lsr(U, V, W):
    """|v - v_LSR|, the halo criterion of Helmi et al. (2018)."""
    pu, pv, pw = V_PECULIAR_KMS
    return np.sqrt((U + pu) ** 2 + (V + pv) ** 2 + (W + pw) ** 2)


def reflex_correct_pm(ra, dec, distance_kpc, pmra, pmdec):
    """Remove the Sun's own motion from proper motions of stars at a given
    distance: what the proper motion would be if the Sun were at rest in the
    Galactocentric frame. For GD-1 at ~8 kpc the Sun's 250 km/s adds several
    mas/yr, which is why Price-Whelan & Bonaca (2018) correct for it."""
    A, _ = _galcen_affine()
    v_sun_icrs = A.T @ np.asarray(V_SUN_KMS)
    e_ra, e_dec = tangent_basis(ra, dec)
    # Observed motion = (v_star - v_sun) projected; so add the Sun's part back.
    d_pmra = (v_sun_icrs @ e_ra) / (K_KMS * distance_kpc)
    d_pmdec = (v_sun_icrs @ e_dec) / (K_KMS * distance_kpc)
    return pmra + d_pmra, pmdec + d_pmdec


# --------------------------------------------------------------------------
# Potentials
# --------------------------------------------------------------------------


class AnalyticHelmi17:
    """Miyamoto-Nagai disc + Hernquist bulge + NFW halo (Helmi et al. 2017)."""

    name = "Helmi et al. (2017)"

    def __init__(self, p: Helmi17Potential = HELMI17):
        self.p = p
        c = p.halo_concentration
        #: NFW normalisation: M(<r) = M_vir m(r/r_s)/m(c), m(x) = ln(1+x) - x/(1+x).
        self._halo_gm = G_KPC_KMS2_MSUN * p.halo_mass_msun / (np.log1p(c) - c / (1 + c))

    def phi(self, R, z):
        """Potential in (km/s)^2 at cylindrical (R, z) in kpc."""
        p, G = self.p, G_KPC_KMS2_MSUN
        R, z = np.asarray(R, float), np.asarray(z, float)
        disc = -G * p.disc_mass_msun / np.sqrt(R**2 + (p.disc_a_kpc + np.sqrt(z**2 + p.disc_b_kpc**2))**2)
        r = np.sqrt(R**2 + z**2)
        bulge = -G * p.bulge_mass_msun / (r + p.bulge_c_kpc)
        x = r / p.halo_rs_kpc
        with np.errstate(divide="ignore", invalid="ignore"):
            halo = np.where(x > 1e-8, -self._halo_gm * np.log1p(x) / np.where(r > 0, r, 1.0),
                            -self._halo_gm / p.halo_rs_kpc)
        return disc + bulge + halo

    def vertical_profile(self, R, z):
        """Phi(R, z) - Phi(R, 0) on a z grid (exact here)."""
        z = np.asarray(z, float)
        return self.phi(np.full_like(z, R), z) - self.phi(np.array([R]), np.array([0.0]))[0]

    def vcirc(self, R):
        """Circular speed in the midplane, km/s."""
        p, G = self.p, G_KPC_KMS2_MSUN
        R = np.asarray(R, float)
        disc = G * p.disc_mass_msun * R**2 / (R**2 + (p.disc_a_kpc + p.disc_b_kpc)**2) ** 1.5
        bulge = G * p.bulge_mass_msun * R / (R + p.bulge_c_kpc) ** 2
        x = R / p.halo_rs_kpc
        halo = self._halo_gm * (np.log1p(x) - x / (1 + x)) / R
        return np.sqrt(disc + bulge + halo)


class GridPotential:
    """A potential tabulated once on an (R, z) grid and spline-interpolated.

    galpy's McMillan (2017) model takes ~4 ms per evaluation; a million stars
    would take an hour. It is smooth, so a bicubic spline in (log R, |z|) on a
    ~10,000-point grid is accurate to better than 0.01% (test_physics.py
    checks this against direct evaluations) and costs microseconds.
    """

    def __init__(self, name, R_grid, z_grid, phi_grid, R_vc, vc, vertical=None):
        from scipy.interpolate import CubicSpline, RectBivariateSpline

        #: (R, z, Phi(R, z) - Phi(R, 0)) evaluated directly on a fine grid. The
        #: disc's density has a cusp at the midplane (exponential in |z|), so
        #: a spline with 25 pc knots gets Phi right to 1e-7 but blurs its
        #: *curvature* near z = 0 by a few per cent -- and the curvature is
        #: what sets how fast stars oscillate vertically.
        self._vertical = vertical
        self.name = name
        self._R_range = (R_grid[0], R_grid[-1])
        self._z_max = z_grid[-1]
        # Mirror the table to negative z, so the midplane is an interior point
        # of the spline and its curvature there -- which sets how fast stars
        # oscillate through the disc -- comes out right. (Tabulated from z = 0
        # only, the spline's end condition bent the first 25 pc by ~8%.)
        z_full = np.concatenate([-z_grid[:0:-1], z_grid])
        phi_full = np.concatenate([phi_grid[:, :0:-1], phi_grid], axis=1)
        self._spline = RectBivariateSpline(np.log(R_grid), z_full, phi_full, kx=3, ky=3)
        self._vc = CubicSpline(np.log(R_vc), vc)
        self._vc_range = (R_vc[0], R_vc[-1])

    def phi(self, R, z):
        R = np.clip(np.asarray(R, float), *self._R_range)
        z = np.asarray(z, float)
        if np.any(np.abs(z) > self._z_max):
            raise ValueError(f"|z| beyond the tabulated {self._z_max} kpc")
        return self._spline.ev(np.log(R), z)

    def vcirc(self, R):
        R = np.clip(np.asarray(R, float), *self._vc_range)
        return self._vc(np.log(R))

    def vertical_profile(self, R, z):
        z = np.abs(np.asarray(z, float))
        if self._vertical is not None and abs(R - self._vertical[0]) < 1e-9:
            return np.interp(z, self._vertical[1], self._vertical[2])
        return self.phi(np.full_like(z, R), z) - self.phi(np.array([R]), np.array([0.0]))[0]


#: Grid for tabulating slow potentials: log-spaced in R (0.05-120 kpc) and
#: denser near the plane in z, where the disc makes the potential curve most.
_R_GRID = np.geomspace(0.05, 120.0, 120)
_Z_GRID = np.concatenate([np.linspace(0, 1, 41)[:-1], np.linspace(1, 5, 33)[:-1],
                          np.geomspace(5, 60, 20)])
_VC_GRID = np.geomspace(0.1, 120.0, 200)


def mcmillan17() -> GridPotential:
    """McMillan (2017) via galpy, tabulated and cached under data/aux/."""
    path = AUX_DIR / "mcmillan17_grid.npz"
    if path.exists():
        d = np.load(path)
        return GridPotential("McMillan (2017)", d["R"], d["z"], d["phi"], d["R_vc"], d["vc"],
                             (float(d["R_vert"]), d["z_vert"], d["phi_vert"]))
    import warnings

    from galpy.potential import evaluatePotentials, vcirc
    from galpy.potential.mwpotentials import McMillan17

    warnings.filterwarnings("ignore")
    ro, vo = 8.21, 233.1  # McMillan17's native galpy units
    RR, ZZ = np.meshgrid(_R_GRID, _Z_GRID, indexing="ij")
    phi = evaluatePotentials(McMillan17, RR.ravel() / ro, ZZ.ravel() / ro,
                             use_physical=False).reshape(RR.shape) * vo**2
    vc = np.array([vcirc(McMillan17, R / ro, use_physical=False) for R in _VC_GRID]) * vo
    # Fine vertical profile at the Sun's radius: 1 pc steps to 0.2 kpc, then 5 pc.
    z_vert = np.unique(np.concatenate([np.arange(0, 0.2, 0.001), np.arange(0.2, 3.0001, 0.005)]))
    phi_v = evaluatePotentials(McMillan17, np.full_like(z_vert, R0_KPC / ro), z_vert / ro,
                               use_physical=False) * vo**2
    phi_vert = phi_v - phi_v[0]
    np.savez(path, R=_R_GRID, z=_Z_GRID, phi=phi, R_vc=_VC_GRID, vc=vc, R_vert=R0_KPC,
             z_vert=z_vert, phi_vert=phi_vert)
    return GridPotential("McMillan (2017)", _R_GRID, _Z_GRID, phi, _VC_GRID, vc,
                         (R0_KPC, z_vert, phi_vert))


def energy(pot, g: Galactocentric):
    """Specific orbital energy, (km/s)^2."""
    return 0.5 * g.speed**2 + pot.phi(g.R, g.z)


def circularity(pot, E, L_z):
    """L_z / L_z,circ(E): +1 for a prograde circular orbit in the plane, 0 for
    a purely radial one, -1 for a retrograde circular orbit (Koppelman et al.
    2019b). L_z,circ(E) is the angular momentum of the circular orbit that has
    the same energy."""
    R = np.geomspace(0.05, 100.0, 2000)
    vc = pot.vcirc(R)
    E_circ = pot.phi(R, 0.0) + 0.5 * vc**2
    L_circ = R * vc
    order = np.argsort(E_circ)
    Lc = np.interp(E, E_circ[order], L_circ[order])
    return L_z / Lc
