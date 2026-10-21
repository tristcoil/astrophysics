"""The 1919 star field in modern terms: Gaia positions, the Sun, GR's deflection.

Table I of Dyson, Eddington & Davidson gives each star's "standard
co-ordinates" x, y (units of 50', x toward the east, y north) and the
displacement expected at Sobral and Principe for 1.75" at the limb. Here we

* propagate the Gaia DR3 positions of the 13 stars to 29 May 1919 (proper
  motion, parallax and radial velocity) and refer them to the mean equator
  and equinox of 1919.4, the frame of Table I;
* find the tangent point, orientation and scale of Table I's coordinates by
  least squares (tests both the identifications and the 1919 star places);
* place the Sun (JPL Horizons, topocentric) in the same coordinates at each
  station's mid-totality, and compute GR's deflection for every star with
  the exact formula delta = (1 + gamma) GM / (c^2 d) cot(chi / 2), where d is
  the Earth-Sun distance and chi the star's angular distance from the Sun.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

import config as C

warnings.filterwarnings("ignore", module="erfa")

UNIT_RAD = np.deg2rad(C.UNIT_TABLE1_ARCMIN / 60.0)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def table1() -> pd.DataFrame:
    """Table I as transcribed (13 stars; x, y in units of 50')."""
    return pd.read_csv(C.TRANSCRIBED_DIR / "dyson1920_table1_stars.csv", comment="#")


def gaia_field() -> pd.DataFrame:
    return pd.read_csv(C.GAIA_DIR / "gaia_dr3_hyades_field.csv")


def gaia_table1() -> pd.DataFrame:
    """Gaia DR3 rows of the 13 Table I stars, indexed by Table I number."""
    g = gaia_field().set_index("Source")
    rows = []
    for no, (src, name, hd) in C.TABLE1_GAIA.items():
        r = g.loc[src].to_dict()
        rows.append({"no": no, "source_id": src, "modern_name": name, "hd": hd, **r})
    return pd.DataFrame(rows).set_index("no")


# ---------------------------------------------------------------------------
# Coordinates
# ---------------------------------------------------------------------------

def _skycoord_gaia(df: pd.DataFrame):
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    plx = np.clip(df["Plx"].to_numpy(float), 0.05, None)
    rv = np.nan_to_num(df["RV"].to_numpy(float))
    return SkyCoord(ra=df["RA_ICRS"].to_numpy(float) * u.deg, dec=df["DE_ICRS"].to_numpy(float) * u.deg,
                    pm_ra_cosdec=np.nan_to_num(df["pmRA"].to_numpy(float)) * u.mas / u.yr,
                    pm_dec=np.nan_to_num(df["pmDE"].to_numpy(float)) * u.mas / u.yr,
                    distance=(1000.0 / plx) * u.pc, radial_velocity=rv * u.km / u.s,
                    obstime=Time(2016.0, format="jyear", scale="tcb"), frame="icrs")


def gaia_at(df: pd.DataFrame, utc: str):
    """Barycentric ICRS directions of Gaia stars at a UTC date (space motion applied)."""
    from astropy.time import Time

    return _skycoord_gaia(df).apply_space_motion(new_obstime=Time(utc, scale="utc"))


def to_equinox(coord, equinox: float = C.EQUINOX_1919):
    """Mean equator and equinox of a Julian epoch (FK5-style, IAU precession)."""
    from astropy.coordinates import FK5
    from astropy.time import Time

    return coord.transform_to(FK5(equinox=Time(equinox, format="jyear")))


def standard_coords(ra, dec, ra0, dec0):
    """Gnomonic (tangent-plane) coordinates, radians; xi east, eta north."""
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return xi, eta


def fit_table1_frame(epoch_utc: str = C.SOBRAL_SECOND_CONTACT_UT) -> dict:
    """Tangent point, rotation and scale that map Gaia (1919, equinox 1919.4)
    onto Table I's x, y. Returns the solution and per-star residuals (arcsec)."""
    from scipy.optimize import least_squares

    t1 = table1().set_index("no")
    g = gaia_table1().loc[t1.index]
    eq = to_equinox(gaia_at(g, epoch_utc))
    ra, de = eq.ra.rad, eq.dec.rad
    xt, yt = t1["x"].to_numpy(float), t1["y"].to_numpy(float)

    def model(p):
        xi, eta = standard_coords(ra, de, p[0], p[1])
        xi, eta = xi / UNIT_RAD, eta / UNIT_RAD
        th, s = p[2], p[3]
        return s * (np.cos(th) * xi - np.sin(th) * eta), s * (np.sin(th) * xi + np.cos(th) * eta)

    def res(p):
        x, y = model(p)
        return np.r_[x - xt, y - yt]

    p0 = [np.deg2rad((4 + 19.5 / 60) * 15), np.deg2rad(21 + 43 / 60), 0.0, 1.0]
    sol = least_squares(res, p0, x_scale=[1e-3, 1e-3, 1e-3, 1e-3])
    x, y = model(sol.x)
    to_as = C.UNIT_TABLE1_ARCMIN * 60.0
    J = sol.jac
    resid = res(sol.x)
    dof = len(resid) - len(sol.x)
    cov = np.linalg.inv(J.T @ J) * (resid @ resid) / dof
    return {
        "ra0_deg": float(np.rad2deg(sol.x[0])), "dec0_deg": float(np.rad2deg(sol.x[1])),
        "rotation_deg": float(np.rad2deg(sol.x[2])), "rotation_err_deg": float(np.rad2deg(np.sqrt(cov[2, 2]))),
        "scale": float(sol.x[3]), "scale_err": float(np.sqrt(cov[3, 3])),
        "params": sol.x, "rms_arcsec": float(np.sqrt(np.mean(resid**2)) * to_as),
        "resid_x_arcsec": (x - xt) * to_as, "resid_y_arcsec": (y - yt) * to_as,
        "x_gaia": x, "y_gaia": y, "no": t1.index.to_numpy(),
    }


def project(frame: dict, coord) -> tuple[np.ndarray, np.ndarray]:
    """Any sky position (SkyCoord, equinox 1919.4) into Table I's x, y units."""
    p = frame["params"]
    xi, eta = standard_coords(coord.ra.rad, coord.dec.rad, p[0], p[1])
    xi, eta = xi / UNIT_RAD, eta / UNIT_RAD
    th, s = p[2], p[3]
    return s * (np.cos(th) * xi - np.sin(th) * eta), s * (np.sin(th) * xi + np.cos(th) * eta)


# ---------------------------------------------------------------------------
# The Sun and the deflection
# ---------------------------------------------------------------------------

def sun_track() -> pd.DataFrame:
    df = pd.read_csv(C.HORIZONS_DIR / "sun_1919_eclipse.csv")
    df["utc"] = pd.to_datetime(df["utc"], format="%Y-%b-%d %H:%M")
    return df


def sun_at(site: str, utc: str) -> dict:
    """Topocentric astrometric Sun (ICRF) interpolated to a UTC time."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    df = sun_track()
    df = df[df["site"] == site]
    t = (df["utc"] - pd.Timestamp(utc)).dt.total_seconds().to_numpy()
    out = {k: float(np.interp(0.0, t, df[k].to_numpy())) for k in
           ("ra_deg", "dec_deg", "el_deg", "az_deg", "ang_diam_arcsec", "delta_au")}
    out["coord"] = SkyCoord(out["ra_deg"] * u.deg, out["dec_deg"] * u.deg, frame="icrs")
    out["radius_arcmin"] = out["ang_diam_arcsec"] / 2 / 60
    return out


def deflection_arcsec(chi_rad, d_au: float, gamma: float = 1.0) -> np.ndarray:
    """Apparent outward displacement of a star at elongation chi from the Sun,
    for an observer at distance d: (1 + gamma) GM / (c^2 d) cot(chi / 2)."""
    return (1 + gamma) * C.GM_SUN / (C.C_LIGHT**2 * d_au * C.AU) / np.tan(np.asarray(chi_rad) / 2) / C.ARCSEC


#: Mid-totality used for the Sun's place at each station. Sobral: 12:00:54 UT,
#: close to the middle of JPL's totality (11:58:10-12:03:28, middle 12:00:49;
#: the report gives only the assumed second contact, 11:58:18 UT, p.300, and no
#: duration; 5 s move the Sun 0.2"); Principe: middle of the exposure series
#: K'-Z, 14:13:09-14:18:03 UT (p.316).
MID_TOTALITY = {"sobral": "1919-05-29T12:00:54", "principe": "1919-05-29T14:15:35"}


def predicted_pattern(frame: dict, site: str, gamma: float = 1.0, limb_value: float | None = None) -> pd.DataFrame:
    """GR's displacement of each Table I star at a station (arcsec, x east, y north),
    plus the Sun's position in Table I units. ``limb_value`` rescales to a given
    deflection at the limb instead (e.g. 1.75 for the 1919 convention)."""
    t1 = table1().set_index("no")
    g = gaia_table1().loc[t1.index]
    utc = MID_TOTALITY[site]
    sun = sun_at(site, utc)
    stars = gaia_at(g, utc)
    chi = stars.separation(sun["coord"]).rad
    defl = deflection_arcsec(chi, sun["delta_au"], gamma)
    if limb_value is not None:
        defl = defl * limb_value / deflection_arcsec(np.deg2rad(sun["radius_arcmin"] / 60), sun["delta_au"], gamma)
    xs, ys = project(frame, to_equinox(stars))
    sx, sy = project(frame, to_equinox(sun["coord"]))
    r = np.hypot(xs - sx, ys - sy)
    return pd.DataFrame({"no": t1.index, "x": xs, "y": ys, "r_units": r, "r_arcmin": r * C.UNIT_TABLE1_ARCMIN,
                         "chi_deg": np.rad2deg(chi), "defl": defl,
                         "dx": defl * (xs - sx) / r, "dy": defl * (ys - sy) / r}).set_index("no").assign(
        sun_x=float(sx[0]) if np.ndim(sx) else float(sx), sun_y=float(sy[0]) if np.ndim(sy) else float(sy),
        sun_radius_arcmin=sun["radius_arcmin"], sun_el_deg=sun["el_deg"])


#: The one printed displacement per station that does not fit the paper's own 1/r pattern
#: (transcribed/dyson1920_table1_stars.csv header; dyson1920_sobral4inch_notes.md A1-A4).
TABLE1_PATTERN_OUTLIER = {"sobral": (13, "dx"), "principe": (8, "dx")}


def sun_from_table1(site: str, drop_outlier: bool = True) -> dict:
    """The Sun's position implied by Table I's predicted displacements: the
    point from which the 13 displacement vectors point radially outward with
    magnitude k / r (k fitted). A check of the 1919 computation. By default the
    one printed value per station that does not fit the pattern is left out."""
    from scipy.optimize import least_squares

    t1 = table1()
    x, y = t1["x"].to_numpy(float), t1["y"].to_numpy(float)
    dx, dy = t1[f"{site}_dx"].to_numpy(float), t1[f"{site}_dy"].to_numpy(float)
    keep = np.ones(2 * len(t1), bool)
    if drop_outlier:
        no, comp = TABLE1_PATTERN_OUTLIER[site]
        i = int(np.flatnonzero(t1["no"].to_numpy() == no)[0])
        keep[i if comp == "dx" else len(t1) + i] = False

    def res(p):
        sx, sy, k = p
        rx, ry = x - sx, y - sy
        r2 = rx**2 + ry**2
        return np.r_[k * rx / r2 - dx, k * ry / r2 - dy][keep]

    sol = least_squares(res, [0.44, -0.26, 0.55])
    resid = res(sol.x)
    return {"sun_x": sol.x[0], "sun_y": sol.x[1], "k_at_unit": sol.x[2],
            "limb_implied": sol.x[2] * C.UNIT_TABLE1_ARCMIN / C.SUN_RADIUS_1919_ARCMIN,
            "rms_arcsec": float(np.sqrt(np.mean(resid**2)))}


def astrometric_xy(frame: dict, utc: str, observer: str = "earth") -> pd.DataFrame:
    """Geocentric *astrometric* positions of the 13 stars at a date (proper
    motion, radial velocity and parallax applied; no aberration and no light
    deflection), in Table I's x, y (units of 50')."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord, get_body_barycentric
    from astropy.time import Time

    t1 = table1().set_index("no")
    g = gaia_table1().loc[t1.index]
    bary = gaia_at(g, utc).cartesian.xyz.to(u.au).value          # 3 x n, barycentric
    earth = get_body_barycentric(observer, Time(utc, scale="utc")).xyz.to(u.au).value
    geo = bary - earth[:, None]
    c = SkyCoord(x=geo[0], y=geo[1], z=geo[2], unit=u.au, representation_type="cartesian", frame="icrs")
    c = SkyCoord(ra=c.spherical.lon, dec=c.spherical.lat, frame="icrs")
    x, y = project(frame, to_equinox(c))
    return pd.DataFrame({"no": t1.index, "x": x, "y": y}).set_index("no")


def epoch_shift(frame: dict, utc_eclipse: str, utc_comparison: str, stars=None) -> pd.DataFrame:
    """How far each star really moved between a comparison and an eclipse
    exposure (proper motion + parallax), in arcsec, after removing the best
    linear plate transformation (which the plate constants absorb anyway).
    Positive = eclipse position minus comparison position."""
    a = astrometric_xy(frame, utc_eclipse)
    b = astrometric_xy(frame, utc_comparison)
    idx = list(a.index if stars is None else stars)
    x, y = b.loc[idx, "x"].to_numpy(), b.loc[idx, "y"].to_numpy()
    dx = (a.loc[idx, "x"] - b.loc[idx, "x"]).to_numpy() * C.UNIT_TABLE1_ARCMIN * 60
    dy = (a.loc[idx, "y"] - b.loc[idx, "y"]).to_numpy() * C.UNIT_TABLE1_ARCMIN * 60
    A = np.column_stack([np.ones_like(x), x, y])
    rx = dx - A @ np.linalg.lstsq(A, dx, rcond=None)[0]
    ry = dy - A @ np.linalg.lstsq(A, dy, rcond=None)[0]
    return pd.DataFrame({"no": idx, "raw_dx": dx, "raw_dy": dy, "dx": rx, "dy": ry}).set_index("no")


def membership() -> pd.DataFrame:
    """Hyades membership flag from Gaia parallax and proper motion (the cluster's
    mean motion ~ (+101, -28) mas/yr at ~21 mas parallax in this field)."""
    g = gaia_table1()
    member = (g["Plx"].between(17, 26)) & (np.hypot(g["pmRA"] - 105, g["pmDE"] + 42) < 25)
    return g.assign(hyades=member)


if __name__ == "__main__":
    f = fit_table1_frame()
    print(f"Table I frame: tangent point RA {f['ra0_deg']:.4f} Dec {f['dec0_deg']:.4f} (equinox 1919.41); "
          f"rotation {f['rotation_deg']:.3f} +- {f['rotation_err_deg']:.3f} deg; scale {f['scale']:.5f} +- {f['scale_err']:.5f}; "
          f"rms {f['rms_arcsec']:.2f}\"")
    for site in ("sobral", "principe"):
        pat = predicted_pattern(f, site, limb_value=1.75)
        imp = sun_from_table1(site)
        print(f"{site}: Sun (Horizons) x={pat['sun_x'].iloc[0]:+.4f} y={pat['sun_y'].iloc[0]:+.4f}; "
              f"implied by Table I x={imp['sun_x']:+.4f} y={imp['sun_y']:+.4f} "
              f"(limb {imp['limb_implied']:.3f}\", rms {imp['rms_arcsec']:.3f}\"); "
              f"Sun radius {pat['sun_radius_arcmin'].iloc[0]:.3f}', elevation {pat['sun_el_deg'].iloc[0]:.1f} deg")
        t1 = table1().set_index("no")
        d = np.hypot(pat["dx"] - t1[f"{site}_dx"], pat["dy"] - t1[f"{site}_dy"])
        print(f"   predicted (1.75\" limb) vs Table I: max |diff| {d.max():.3f}\", rms {np.sqrt(np.mean(d**2)):.3f}\"")
    print(membership()[["modern_name", "Plx", "pmRA", "pmDE", "hyades"]])
