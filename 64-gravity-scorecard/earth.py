"""The ranging station on the rotating, tidally flexing Earth.

* Earth orientation: IERS EOP 20 C04 (daily) interpolated with 4-point Lagrange
  polynomials, plus the diurnal/semidiurnal ocean-tide variations of polar
  motion and UT1 (IERS 2010 ch. 8, ``ORTHO_EOP.F``/``CNMTX.F`` ported below and
  checked against their documented test case). The celestial-to-terrestrial
  matrix is built with ERFA (IAU 2006/2000A, CIO based) including the C04
  celestial pole offsets.
* Station displacements (IERS 2010 ch. 7 as used by Williams & Boggs 2020):
  degree-2 and -3 solid tides with real Love numbers (constant part included),
  the K1 step-2 correction, and the solid pole tide. Ocean loading (< 1 cm at
  inland Apache Point) and atmospheric loading (~2 mm) are not modelled.
"""

from __future__ import annotations

from functools import lru_cache

import erfa
import numpy as np
import pandas as pd

import config as C

ARCSEC = np.pi / 180 / 3600


# ---------------------------------------------------------------------------
# Earth orientation parameters
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def load_eop() -> pd.DataFrame:
    cols = ["yr", "mm", "dd", "hh", "mjd", "x", "y", "ut1_utc", "dX", "dY", "xrt", "yrt", "lod",
            "x_err", "y_err", "ut1_err", "dX_err", "dY_err", "xrt_err", "yrt_err", "lod_err"]
    df = pd.read_csv(C.IERS_DIR / C.IERS_EOP_FILE, comment="#", sep=r"\s+", names=cols)
    # UT1 - TAI is continuous across leap seconds; interpolate that instead of UT1 - UTC
    dat = np.array([erfa.dat(int(y), int(m), int(d), 0.0) for y, m, d in zip(df.yr, df.mm, df.dd)])
    df["ut1_tai"] = df["ut1_utc"] - dat
    return df


def _lagrange4(xq, xs, ys):
    """4-point Lagrange interpolation on a uniform daily grid (IERS recommendation)."""
    i = np.searchsorted(xs, xq) - 2
    i = np.clip(i, 0, len(xs) - 4)
    out = np.zeros_like(xq, dtype=float)
    for j in range(4):
        lj = np.ones_like(xq, dtype=float)
        for k in range(4):
            if k != j:
                lj *= (xq - xs[i + k]) / (xs[i + j] - xs[i + k])
        out += lj * ys[i + j]
    return out


def eop_at(mjd_utc) -> dict:
    """x_p, y_p (arcsec), UT1-TAI (s), dX, dY (arcsec) at UTC MJD, with sub-daily ocean tides."""
    df = load_eop()
    xs = df["mjd"].to_numpy(float)
    q = np.asarray(mjd_utc, float)
    out = {k: _lagrange4(q, xs, df[k].to_numpy(float)) for k in ("x", "y", "ut1_tai", "dX", "dY")}
    dx, dy, dut = ortho_eop(q)
    out["x"] = out["x"] + dx * 1e-6
    out["y"] = out["y"] + dy * 1e-6
    out["ut1_tai"] = out["ut1_tai"] + dut * 1e-6
    return out


# --- IERS 2010 ch. 8: ORTHO_EOP / CNMTX (Ray et al. 1994 ocean-tide EOP) ---

_ORTHOW = np.array([
    -6.77832, -14.86323, 0.47884, -1.45303, 0.16406, 0.42030,
    0.09398, 25.73054, -4.77974, 0.28080, 1.94539, -0.73089,
    14.86283, -6.77846, 1.45234, 0.47888, -0.42056, 0.16469,
    15.30276, -4.30615, 0.07564, 2.28321, -0.45717, -1.62010,
    -1.76335, 1.03364, -0.27553, 0.34569, -0.12343, -0.10146,
    -0.47119, 1.28997, -0.19336, 0.02724, 0.08955, 0.04726]).reshape(3, 12)   # Fortran column order

# (m, HS, PHASE, FREQ) of the 71 lines; n = 2 for all
_CNMTX = np.array([
    (1, -1.94, 9.0899831, 5.18688050), (1, -1.25, 8.8234208, 5.38346657), (1, -6.64, 12.1189598, 5.38439079),
    (1, -1.51, 1.4425700, 5.41398343), (1, -8.02, 4.7381090, 5.41490765), (1, -9.47, 4.4715466, 5.61149372),
    (1, -50.20, 7.7670857, 5.61241794), (1, -1.80, -2.9093042, 5.64201057), (1, -9.54, 0.3862349, 5.64293479),
    (1, 1.52, -3.1758666, 5.83859664), (1, -49.45, 0.1196725, 5.83952086), (1, -262.21, 3.4152116, 5.84044508),
    (1, 1.70, 12.8946194, 5.84433381), (1, 3.43, 5.5137686, 5.87485066), (1, 1.94, 6.4441883, 6.03795537),
    (1, 1.37, -4.2322016, 6.06754801), (1, 7.41, -0.9366625, 6.06847223), (1, 20.62, 8.5427453, 6.07236095),
    (1, 4.14, 11.8382843, 6.07328517), (1, 3.94, 1.1618945, 6.10287781), (1, -7.14, 5.9693878, 6.24878055),
    (1, 1.37, -1.2032249, 6.26505830), (1, -122.03, 2.0923141, 6.26598252), (1, 1.02, -1.7847596, 6.28318449),
    (1, 2.89, 8.0679449, 6.28318613), (1, -7.30, 0.8953321, 6.29946388), (1, 368.78, 4.1908712, 6.30038810),
    (1, 50.01, 7.4864102, 6.30131232), (1, -1.08, 10.7819493, 6.30223654), (1, 2.93, 0.3137975, 6.31759007),
    (1, 5.25, 6.2894282, 6.33479368), (1, 3.95, 7.2198478, 6.49789839), (1, 20.62, -0.1610030, 6.52841524),
    (1, 4.09, 3.1345361, 6.52933946), (1, 3.42, 2.8679737, 6.72592553), (1, 1.69, -4.5128771, 6.75644239),
    (1, 11.29, 4.9665307, 6.76033111), (1, 7.23, 8.2620698, 6.76125533), (1, 1.51, 11.5576089, 6.76217955),
    (1, 2.16, 0.6146566, 6.98835826), (1, 1.38, 3.9101957, 6.98928248), (2, 1.80, 20.6617051, 11.45675174),
    (2, 4.67, 13.2808543, 11.48726860), (2, 16.01, 16.3098310, 11.68477889), (2, 19.32, 8.9289802, 11.71529575),
    (2, 1.30, 5.0519065, 11.73249771), (2, -1.02, 15.8350306, 11.89560406), (2, -4.51, 8.6624178, 11.91188181),
    (2, 120.99, 11.9579569, 11.91280603), (2, 1.13, 8.0808832, 11.93000800), (2, 22.98, 4.5771061, 11.94332289),
    (2, 1.06, 0.7000324, 11.96052486), (2, -1.90, 14.9869335, 12.11031632), (2, -2.18, 11.4831564, 12.12363121),
    (2, -23.58, 4.3105437, 12.13990896), (2, 631.92, 7.6060827, 12.14083318), (2, 1.92, 3.7290090, 12.15803515),
    (2, -4.66, 10.6350594, 12.33834347), (2, -17.86, 3.2542086, 12.36886033), (2, 4.47, 12.7336164, 12.37274905),
    (2, 1.97, 16.0291555, 12.37367327), (2, 17.20, 10.1602590, 12.54916865), (2, 294.00, 6.2831853, 12.56637061),
    (2, -2.46, 2.4061116, 12.58357258), (2, -1.02, 5.0862033, 12.59985198), (2, 79.96, 8.3817423, 12.60077620),
    (2, 23.83, 11.6772814, 12.60170041), (2, 2.59, 14.9728205, 12.60262463), (2, 4.47, 4.0298682, 12.82880334),
    (2, 1.95, 7.3254073, 12.82972756), (2, 1.17, 9.1574019, 13.06071921)])
_SP = np.array([[0.0298, 0.1408, 0.0805, 0.6002, 0.3025, 0.1517],
                [0.0200, 0.0905, 0.0638, 0.3476, 0.1645, 0.0923]])


def _cnmtx(mjd):
    mjd = np.atleast_1d(np.asarray(mjd, float))
    two_pi = 2 * np.pi
    m_arr, hs, phase, freq = _CNMTX.T
    n = 2
    pinm = ((n + m_arr) % 2) * two_pi / 4.0
    anm = {m: {} for m in (1, 2)}
    bnm = {m: {} for m in (1, 2)}
    for k in (-1, 0, 1):
        dt60 = (mjd - k * 2.0) - 37076.5
        alpha = np.mod(phase - pinm, two_pi)[None, :] + np.mod(np.outer(dt60, freq), two_pi)
        for m in (1, 2):
            sel = m_arr == m
            anm[m][k] = (hs[sel] * np.cos(alpha[:, sel])).sum(axis=1)
            bnm[m][k] = -(hs[sel] * np.sin(alpha[:, sel])).sum(axis=1)
    h = []
    for m in (1, 2):
        sp = _SP[m - 1]
        a, b = anm[m], bnm[m]
        ap, am = a[1] + a[-1], a[1] - a[-1]
        bp, bm = b[1] + b[-1], b[1] - b[-1]
        p0 = sp[0] * a[0]
        p1 = sp[1] * a[0] - sp[2] * ap
        p2 = sp[3] * a[0] - sp[4] * ap + sp[5] * bm
        q0 = sp[0] * b[0]
        q1 = sp[1] * b[0] - sp[2] * bp
        q2 = sp[3] * b[0] - sp[4] * bp - sp[5] * am
        h += [p0, q0, p1, q1, p2, q2]
    return np.array(h)          # (12, N)


def ortho_eop(mjd):
    """Ocean-tide sub-daily delta_x, delta_y (micro-arcsec) and delta_UT1 (micro-s)."""
    h = _cnmtx(mjd)
    return tuple(_ORTHOW[k] @ h for k in range(3))


# ---------------------------------------------------------------------------
# Celestial <-> terrestrial rotation
# ---------------------------------------------------------------------------

def c2t_matrix(tt_jd1, tt_jd2, ut1_jd1, ut1_jd2, xp_as, yp_as, dX_as, dY_as):
    """Celestial (GCRS) to terrestrial (ITRS) matrices, shape (N, 3, 3)."""
    X, Y = erfa.xy06(tt_jd1, tt_jd2)
    X = X + dX_as * ARCSEC
    Y = Y + dY_as * ARCSEC
    s = erfa.s06(tt_jd1, tt_jd2, X, Y)
    rc2i = erfa.c2ixys(X, Y, s)
    era = erfa.era00(ut1_jd1, ut1_jd2)
    sp = erfa.sp00(tt_jd1, tt_jd2)
    rpom = erfa.pom00(xp_as * ARCSEC, yp_as * ARCSEC, sp)
    return erfa.c2tcio(rc2i, era, rpom)


# ---------------------------------------------------------------------------
# Station geometry
# ---------------------------------------------------------------------------

def geocentric_to_xyz(radius_m, lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.array([radius_m * np.cos(la) * np.cos(lo), radius_m * np.cos(la) * np.sin(lo), radius_m * np.sin(la)])


def geodetic(xyz):
    """Geodetic latitude, longitude (rad) and height (m) on the GRS80 ellipsoid."""
    lam, phi, h = erfa.gc2gde(C.WGS_A, C.WGS_F, np.asarray(xyz, float))   # ERFA returns (elong, phi, height)
    return phi, lam, h


def enu_basis(xyz):
    """Unit vectors (east, north, up) of the geodetic horizon at a station (ITRS)."""
    phi, lam, _ = geodetic(xyz)
    e = np.array([-np.sin(lam), np.cos(lam), 0.0])
    n = np.array([-np.sin(phi) * np.cos(lam), -np.sin(phi) * np.sin(lam), np.cos(phi)])
    u = np.array([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)])
    return e, n, u


def solid_tide(r_sta, r_moon, r_sun):
    """IERS 2010 eqs. (5) + (6): degree-2 (Moon, Sun) and degree-3 (Moon) solid tides.

    r_sta (N,3) station, r_moon/r_sun (N,3) geocentric positions, all ITRS metres.
    Real nominal Love numbers, constant (permanent) part included.
    """
    L = C.EARTH_LOVE
    rs = np.linalg.norm(r_sta, axis=1, keepdims=True)
    rhat = r_sta / rs
    out = np.zeros_like(r_sta)
    for r_body, gm, deg3 in ((r_moon, C.GM_MOON, True), (r_sun, C.GM_SUN, False)):
        rb = np.linalg.norm(r_body, axis=1, keepdims=True)
        bhat = r_body / rb
        d = np.sum(bhat * rhat, axis=1, keepdims=True)
        f2 = gm / C.GM_EARTH * C.R_EARTH_EQ**4 / rb**3
        out += f2 * (L["h2"] * rhat * (1.5 * d**2 - 0.5) + 3 * L["l2"] * d * (bhat - d * rhat))
        if deg3:
            f3 = gm / C.GM_EARTH * C.R_EARTH_EQ**5 / rb**4
            out += f3 * (L["h3"] * rhat * (2.5 * d**3 - 1.5 * d)
                         + L["l3"] * (7.5 * d**2 - 1.5) * (bhat - d * rhat))
    return out


def k1_correction(r_sta, gmst, omega_node):
    """IERS 2010 step 2, eq. (12a) radial part for K1 and its nodal sideband (metres)."""
    phi, lam, _ = geodetic(r_sta[0])
    _, _, u = enu_basis(r_sta[0])
    dr = 0.0
    for name, theta in (("K1", gmst + np.pi), ("K1_side", gmst + np.pi - omega_node)):
        ip, op = C.K1_STEP2[name]
        dr = dr + (ip * np.sin(theta + lam) + op * np.cos(theta + lam)) * np.sin(2 * phi)
    return np.outer(dr, u)


def pole_tide(r_sta, xp_as, yp_as, year):
    """Solid-Earth pole tide (IERS 2010 eq. 26 with W&B 2020's mean pole and coefficients)."""
    mp = C.POLE_TIDE_MEAN
    m1 = xp_as - (mp["x0"] + mp["xdot"] * (year - 2000.0))
    m2 = -(yp_as - (mp["y0"] + mp["ydot"] * (year - 2000.0)))
    phi, lam, _ = geodetic(r_sta[0])
    theta = np.pi / 2 - phi
    e, n, u = enu_basis(r_sta[0])
    k_r, k_h = C.POLE_TIDE_COEF["radial"], C.POLE_TIDE_COEF["horizontal"]
    s_theta = -k_h * np.cos(2 * theta) * (m1 * np.cos(lam) + m2 * np.sin(lam))   # southward
    s_lam = k_h * np.cos(theta) * (m1 * np.sin(lam) - m2 * np.cos(lam))           # eastward
    s_r = -k_r * np.sin(2 * theta) * (m1 * np.cos(lam) + m2 * np.sin(lam))        # up
    return np.outer(-s_theta, n) + np.outer(s_lam, e) + np.outer(s_r, u)
