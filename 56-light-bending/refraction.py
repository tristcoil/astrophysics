"""Differential refraction and aberration across a photographic plate.

The air bends light toward the zenith by R ~ 58" tan(z), so stars nearer the
horizon are lifted more than those above them and the field is squeezed along
the vertical (by ~ 58" sec^2 z per radian) and slightly across it. The Earth's
orbital motion tilts every star toward the apex of its motion by ~20" (annual
aberration), again by a slightly different amount at each edge of the field.
Across a 2-degree plate both effects are linear transformations of the star
positions: changes of scale and orientation. Eclipse plates (taken by day in
May) and comparison plates (taken at night months later) suffer different
amounts, and Dyson, Eddington & Davidson computed and removed the difference
"by the usual formulae" before discussing the plate scales.

Here we compute those plate constants numerically, without the 1919
approximations: a grid of points around the field centre is moved (1) by
annual + diurnal aberration with ERFA's rigorous formula and (2) by
refraction along the vertical, using ERFA's refraction constants for the
recorded temperature and pressure; then the linear model
    dx = a x + b y + c,   dy = d x + e y + f
is fitted to the displacements (x, y and the constants in the same units, so
a and e are fractional scale changes).
"""

from __future__ import annotations

import numpy as np

import stars as S


def _grid(half_width_rad: float, n: int = 7):
    u = np.linspace(-half_width_rad, half_width_rad, n)
    return [g.ravel() for g in np.meshgrid(u, u)]


def _tangent_to_sky(xi, eta, ra0, dec0):
    """Inverse gnomonic projection (radians)."""
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / rho, 0.0))
    ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return ra, dec


def _fit_linear(x, y, dx, dy) -> dict:
    A = np.column_stack([x, y, np.ones_like(x)])
    (a, b, c), *_ = np.linalg.lstsq(A, dx, rcond=None)
    (d, e, f), *_ = np.linalg.lstsq(A, dy, rcond=None)
    return {"a": a, "b": b, "c": c, "d": d, "e": e, "f": f}


def plate_effects(ra0_deg: float, dec0_deg: float, utc: str, lon_deg: float, lat_deg: float,
                  height_m: float = 0.0, temp_c: float = 15.0, pressure_hpa: float = 1013.25,
                  rel_humidity: float = 0.5, wavelength_um: float = 0.43, half_width_deg: float = 1.0,
                  frame_equinox: float | None = None) -> dict:
    """Linear plate constants of aberration and refraction for a field centred
    at (ra0, dec0) [ICRS degrees] photographed at ``utc`` from a site.

    Returns fractional constants for 'aberration', 'refraction' and 'total'
    (displacement = apparent - true, in the field's tangent-plane coordinates
    with x east and y north of the mean equator of ``frame_equinox``; ICRS axes
    if None), plus the field's altitude and the parallactic angle.
    """
    import astropy.units as u
    import erfa
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_body_barycentric_posvel
    from astropy.time import Time

    t = Time(utc, scale="utc")
    loc = EarthLocation.from_geodetic(lon_deg * u.deg, lat_deg * u.deg, height_m * u.m)
    ra0, dec0 = np.deg2rad(ra0_deg), np.deg2rad(dec0_deg)
    xi, eta = _grid(np.deg2rad(half_width_deg))
    ra, dec = _tangent_to_sky(xi, eta, ra0, dec0)
    p = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], axis=-1)

    # --- aberration (annual + diurnal), ERFA's relativistic formula ----------
    _, vel = get_body_barycentric_posvel("earth", t)
    v_earth = vel.xyz.to(u.m / u.s).value
    _, gvel = loc.get_gcrs_posvel(t)
    v_obs = v_earth + gvel.xyz.to(u.m / u.s).value
    beta = v_obs / 299792458.0
    bm1 = np.sqrt(1 - beta @ beta)
    # erfa.ab(pnat, v, s, bm1): s (Sun distance, au) only enters the tiny
    # gravitational term; 1 au is fine here.
    p_ab = np.array([erfa.ab(pi, beta, 1.0, bm1) for pi in p])

    def to_xy(vecs):
        r = np.arctan2(vecs[:, 1], vecs[:, 0])
        d = np.arcsin(np.clip(vecs[:, 2] / np.linalg.norm(vecs, axis=1), -1, 1))
        return S.standard_coords(r, d, ra0, dec0)

    x0, y0 = to_xy(p)
    xa, ya = to_xy(p_ab)

    # --- refraction along the vertical ---------------------------------------
    # True topocentric alt/az of the aberrated directions (astropy, no refraction),
    # then lift each point by ERFA's R = A tan z + B tan^3 z.
    sky = SkyCoord(ra=ra * u.rad, dec=dec * u.rad, frame="icrs")
    aa = sky.transform_to(AltAz(obstime=t, location=loc, pressure=0 * u.hPa))
    A, B = erfa.refco(pressure_hpa, temp_c, rel_humidity, wavelength_um)
    z = np.pi / 2 - aa.alt.rad
    dz = A * np.tan(z) + B * np.tan(z) ** 3
    lifted = SkyCoord(az=aa.az, alt=(aa.alt.rad + dz) * u.rad,
                      frame=AltAz(obstime=t, location=loc, pressure=0 * u.hPa)).transform_to("icrs")
    xr_, yr_ = S.standard_coords(lifted.ra.rad, lifted.dec.rad, ra0, dec0)
    xt_, yt_ = S.standard_coords(sky.ra.rad, sky.dec.rad, ra0, dec0)
    # refraction displacement field, measured in the ICRS-axis tangent plane
    xr, yr = xa + (xr_ - xt_), ya + (yr_ - yt_)

    # optionally rotate the axes to the mean equator of another equinox
    rot = 0.0
    if frame_equinox is not None:
        c0 = SkyCoord(ra0 * u.rad, dec0 * u.rad, frame="icrs")
        north = SkyCoord(ra0 * u.rad, (dec0 + 1e-4) * u.rad, frame="icrs")
        e0, e1 = S.to_equinox(c0, frame_equinox), S.to_equinox(north, frame_equinox)
        rot = float(e0.position_angle(e1).rad)   # PA of ICRS north in the equinox frame

    def rotate(x, y):
        return np.cos(rot) * x + np.sin(rot) * y, -np.sin(rot) * x + np.cos(rot) * y

    X0, Y0 = rotate(x0, y0)
    out = {}
    for key, (x1, y1) in {"aberration": (xa, ya), "total": (xr, yr)}.items():
        X1, Y1 = rotate(x1, y1)
        out[key] = _fit_linear(X0, Y0, X1 - X0, Y1 - Y0)
    out["refraction"] = {k: out["total"][k] - out["aberration"][k] for k in out["total"]}
    centre = SkyCoord(ra=ra0 * u.rad, dec=dec0 * u.rad, frame="icrs").transform_to(
        AltAz(obstime=t, location=loc, pressure=0 * u.hPa))
    lst = t.sidereal_time("apparent", longitude=loc.lon).rad
    ha = lst - ra0
    q = np.arctan2(np.sin(ha), np.tan(np.deg2rad(lat_deg)) * np.cos(dec0) - np.sin(dec0) * np.cos(ha))
    out.update({"alt_deg": float(centre.alt.deg), "az_deg": float(centre.az.deg),
                "zd_deg": 90.0 - float(centre.alt.deg), "parallactic_deg": float(np.rad2deg(q)),
                "refco_A_arcsec": float(np.rad2deg(A) * 3600), "refco_B_arcsec": float(np.rad2deg(B) * 3600),
                "v_obs_kms": float(np.linalg.norm(v_obs) / 1e3)})
    return out


def to_revolutions(const: dict, rev_arcsec: float, unit_arcmin: float = 50.0) -> dict:
    """Express fractional plate constants as the 1919 papers did: revolutions
    of the micrometer screw per unit (50') of field coordinate."""
    k = unit_arcmin * 60.0 / rev_arcsec
    return {key: val * k for key, val in const.items() if key in "abde"}


if __name__ == "__main__":
    # sanity check: Sobral at mid-eclipse, provisional weather
    ra0, dec0 = 66.05, 21.85
    r = plate_effects(ra0, dec0, "1919-05-29T12:00:00", -40.35, -3.69, 70.0, temp_c=26.0, pressure_hpa=1005.0,
                      frame_equinox=1919.41)
    for k in ("aberration", "refraction", "total"):
        rv = to_revolutions(r[k], 6.25)
        print(k, {kk: f"{vv:+.4f}" for kk, vv in rv.items()}, "(4-inch rev)")
    print({k: r[k] for k in ("alt_deg", "az_deg", "parallactic_deg", "refco_A_arcsec", "v_obs_kms")})
