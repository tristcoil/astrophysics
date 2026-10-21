"""Gravity Probe B: general relativity's prediction for a gyroscope in its real orbit.

A gyroscope's spin axis S, free of torques, turns as dS/dt = Omega x S with

    Omega_geodetic  = (3/2) (v x a_N) / c^2                 (de Sitter 1916, Fokker)
    Omega_frame     = (G / c^2 r^3) [3 (J . rhat) rhat - J]   (Schiff 1960; Lense-Thirring)
    Omega_Sun       = (3/2) (V x a_Sun) / c^2               (the Earth-Moon system's own
                                                            geodetic precession round the Sun)

where a_N is the Earth's Newtonian gravitational acceleration at the gyroscope
(central term plus the zonal harmonics: this is how the Earth's flattening J2
enters), J the Earth's spin angular momentum, V the gyroscope's velocity
relative to the Sun. GP-B flew drag-free, so no Thomas precession from
non-gravitational forces.

We integrate GP-B's orbit (Stanford's published elements: a = 7027.4 km,
e = 0.0014, i = 90.007 deg, polar, orbital plane containing the guide star
IM Pegasi) over the science phase (28 Aug 2004 - 14 Aug 2005) with the zonal
field J2-J6 of GOCO06s and lunisolar perturbations, integrate Omega along it,
and project the drift of a gyroscope initially pointing at IM Peg onto the
north-south (NS, in the orbital plane) and west-east (WE, across it) directions
that the experiment measured. Signs follow Everitt et al. (2011): ê_WE =
ê_GS x z (local west at the star), ê_NS = ê_WE x ê_GS (local north).
"""

from __future__ import annotations

import datetime as dt

import numpy as np
from numba import njit
from scipy.integrate import solve_ivp

import arcmodel as AM
import config as C
import frames as F
import gravity
import literature as L

AS = np.pi / 180 / 3600
MAS = AS / 1000


def guide_star_unit(mjd: float | None = None) -> np.ndarray:
    """IM Peg (VLBI, Bartel et al. 2015), optionally moved by its proper motion to mjd."""
    g = L.IM_PEG
    ra = np.radians(g["ra_deg"])
    dec = np.radians(g["dec_deg"])
    if mjd is not None:
        dt_yr = (mjd - g["epoch_mjd"]) / 365.25
        dec = dec + g["pmdec_mas_yr"] * MAS * dt_yr
        ra = ra + g["pmra_mas_yr"] * MAS * dt_yr / np.cos(dec)
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def ns_we_basis(gs: np.ndarray):
    """(ê_NS, ê_WE) at the guide star, with Everitt et al.'s definitions."""
    z = np.array([0.0, 0.0, 1.0])
    we = np.cross(gs, z)
    we /= np.linalg.norm(we)
    ns = np.cross(we, gs)
    return ns, we


def elements_to_state(a, e, i, argp, node, M, gm=C.GM_EARTH):
    """Keplerian elements (rad) -> position, velocity."""
    E = M
    for _ in range(30):
        E = E - (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    p = a * (1 - e * e)
    r = p / (1 + e * np.cos(nu))
    rp = np.array([r * np.cos(nu), r * np.sin(nu), 0.0])
    vp = np.sqrt(gm / p) * np.array([-np.sin(nu), e + np.cos(nu), 0.0])
    cO, sO, ci, si, cw, sw = np.cos(node), np.sin(node), np.cos(i), np.sin(i), np.cos(argp), np.sin(argp)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    return R @ rp, R @ vp


@njit(cache=True)
def zonal_accel(r, gm, re, J):
    """Central + zonal (J2..Jn, unnormalised, J_n = -C_n0) acceleration, spin axis = z."""
    x, y, z = r[0], r[1], r[2]
    rn = np.sqrt(x * x + y * y + z * z)
    u = z / rn
    # Legendre P_n(u) and derivatives
    nmax = J.shape[0] - 1
    P = np.zeros(nmax + 2)
    dP = np.zeros(nmax + 2)
    P[0] = 1.0
    P[1] = u
    dP[1] = 1.0
    for n in range(2, nmax + 1):
        P[n] = ((2 * n - 1) * u * P[n - 1] - (n - 1) * P[n - 2]) / n
        dP[n] = n * P[n - 1] + u * dP[n - 1]
    # potential U = gm/r [1 - sum J_n (re/r)^n P_n(u)]; a = grad U
    ar = -gm / rn ** 2
    au = 0.0
    for n in range(2, nmax + 1):
        f = J[n] * (re / rn) ** n
        ar += gm / rn ** 2 * (n + 1) * f * P[n]
        au += -gm / rn ** 2 * f * dP[n]
    # du/dr = (zhat - u rhat)/r
    a = np.empty(3)
    a[0] = ar * x / rn + au * (0.0 - u * x / rn)
    a[1] = ar * y / rn + au * (0.0 - u * y / rn)
    a[2] = ar * z / rn + au * (1.0 - u * z / rn)
    return a


def zonal_coefficients(mjd: float, nmax: int = 6) -> np.ndarray:
    fld = gravity.load(C.GRAVITY_MODEL, nmax, mjd=mjd)
    N = gravity.normalisation(nmax)
    J = np.zeros(nmax + 1)
    for n in range(2, nmax + 1):
        J[n] = -fld.cbar[n, 0] * N[n, 0]
    return J


def predict(n_days: float | None = None, step_s: float = 20.0, nmax: int = 6,
            lunisolar: bool = True) -> dict:
    """Integrate GP-B's orbit and the gyroscope precession over the science phase."""
    o = L.GPB_ORBIT
    start = F.date_to_mjd(L.GPB_SCIENCE[0])
    end = F.date_to_mjd(L.GPB_SCIENCE[1])
    n_days = end - start if n_days is None else n_days
    mjd_mid = start + n_days / 2
    J = zonal_coefficients(mjd_mid, nmax)
    gs0 = guide_star_unit()              # GR prediction: precession relative to inertial space
    ns, we = ns_we_basis(gs0)

    # Orient the orbital plane so that it contains IM Peg (the design condition).
    node = np.radians(o["node_deg"])
    inc = np.radians(o["i_deg"])
    # tweak the node so the guide star lies exactly in the initial plane
    nhat = lambda nd: np.array([np.sin(inc) * np.sin(nd), -np.sin(inc) * np.cos(nd), np.cos(inc)])
    for _ in range(20):
        f = nhat(node) @ gs0
        df = (nhat(node + 1e-7) @ gs0 - f) / 1e-7
        node -= f / df
    # Sun and Moon on an hourly grid (geocentric, ICRF axes; the 0.07 deg between the
    # J2000 pole and the 2005 spin axis is irrelevant here)
    t_grid = np.arange(-2 * 3600.0, (n_days + 1) * 86400.0, 3600.0)
    bodies, vel = AM.body_positions(start + t_grid / 86400.0, with_velocity=True)
    sun, moon, sun_v = bodies[:, 0], bodies[:, 1], vel[:, 0]

    def rhs(t, y, third_bodies=lunisolar):
        r, v = y[:3], y[3:6]
        a = zonal_accel(r, C.GM_EARTH, C.R_EARTH, J)
        if third_bodies:
            x = (t - t_grid[0]) / 3600.0
            k = int(x)
            w = x - k
            for grid, gmb in ((sun, C.GM_SUN), (moon, C.GM_MOON)):
                body = (1 - w) * grid[k] + w * grid[k + 1]
                d = body - r
                a = a + gmb * (d / np.linalg.norm(d) ** 3 - body / np.linalg.norm(body) ** 3)
        return np.concatenate([v, a])

    # The published a = 7027.4 km is a mean value ((apogee + perigee radius)/2 from
    # the 659.1/639.5 km heights). Choose the initial osculating a so that the mean
    # radius over two days matches it (J2 makes osculating a swing by ~ +-9 km).
    a_osc = o["a_m"]
    target = o["a_m"]
    for _ in range(4):
        r0, v0 = elements_to_state(a_osc, o["e"], inc, np.radians(o["argp_deg"]), node, 0.0)
        tt = np.arange(0.0, 2 * 86400.0, 30.0)
        s2 = solve_ivp(lambda t, y: rhs(t, y, False), (0.0, tt[-1]), np.concatenate([r0, v0]),
                       method="DOP853", t_eval=tt, rtol=1e-11, atol=1e-6)
        mean_r = np.mean(np.linalg.norm(s2.y[:3].T, axis=1))
        a_osc += target - mean_r
    r0, v0 = elements_to_state(a_osc, o["e"], inc, np.radians(o["argp_deg"]), node, 0.0)

    t_eval = np.arange(0.0, n_days * 86400.0, step_s)
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), np.concatenate([r0, v0]), method="DOP853",
                    t_eval=t_eval, rtol=1e-11, atol=1e-6)
    R = sol.y[:3].T
    V = sol.y[3:].T
    c2 = C.C_LIGHT ** 2
    aN = np.array([zonal_accel(r, C.GM_EARTH, C.R_EARTH, J) for r in R])
    aN0 = -C.GM_EARTH * R / np.linalg.norm(R, axis=1)[:, None] ** 3
    om_geo = 1.5 * np.cross(V, aN) / c2
    om_geo_point = 1.5 * np.cross(V, aN0) / c2
    rn = np.linalg.norm(R, axis=1)
    rhat = R / rn[:, None]
    Jvec = np.array([0.0, 0.0, C.EARTH_SPIN_PER_MASS * C.GM_EARTH / C.G_NEWTON])
    om_fd = C.G_NEWTON / (c2 * rn[:, None] ** 3) * (3 * (rhat @ Jvec)[:, None] * rhat - Jvec)
    # the Sun: V = velocity relative to the Sun = Earth's heliocentric + orbital
    s_t = np.array([np.interp(t_eval, t_grid, sun[:, k]) for k in range(3)]).T
    sv_t = np.array([np.interp(t_eval, t_grid, sun_v[:, k]) for k in range(3)]).T
    d = s_t - R                       # gyroscope -> Sun
    a_sun = C.GM_SUN * d / np.linalg.norm(d, axis=1)[:, None] ** 3
    V_helio = V - sv_t                # gyroscope velocity relative to the Sun
    om_sun = 1.5 * np.cross(V_helio, a_sun) / c2

    def drift(om):
        """Integrate dS/dt = om x S (first order) and project on NS/WE, mas."""
        dS = np.cumsum(np.cross(om, gs0), axis=0) * step_s
        return dS @ ns / MAS, dS @ we / MAS

    yr = t_eval / C.JULIAN_YEAR_S
    out = {"t_days": t_eval / 86400.0, "node_deg": np.degrees(node)}
    for name, om in (("geodetic", om_geo), ("geodetic_point_mass", om_geo_point),
                     ("frame_dragging", om_fd), ("solar_geodetic", om_sun)):
        dns, dwe = drift(om)
        out[name] = {"ns_mas_yr": float(np.polyfit(yr, dns, 1)[0]),
                     "we_mas_yr": float(np.polyfit(yr, dwe, 1)[0]),
                     "ns_track": dns[::180], "we_track": dwe[::180]}
    # circular-orbit formulas for comparison
    a = o["a_m"]
    out["circular_geodetic_mas_yr"] = float(1.5 * C.GM_EARTH ** 1.5 / (c2 * a ** 2.5)
                                            * C.JULIAN_YEAR_S / MAS)
    out["circular_frame_mas_yr"] = float(C.G_NEWTON * Jvec[2] / (2 * c2 * a ** 3)
                                         * C.JULIAN_YEAR_S / MAS)
    el = F.kepler_elements(R[::360], V[::360])
    out["mean_a_km"] = float(np.mean(el["a"]) / 1e3)
    out["mean_r_km"] = float(np.mean(rn) / 1e3)
    out["initial_osculating_a_km"] = float(a_osc / 1e3)
    out["node_drift_deg"] = float(np.degrees(np.unwrap(el["node"])[-1] - el["node"][0]))
    out["guide_star_offplane_deg"] = np.degrees(np.arcsin(el["h"] @ gs0))[::10]
    return out


def radius_scan(radii_km=(7018.0, 7020.0, 7022.0, 7024.0, 7026.0, 7027.4, 7030.0), n_days=60.0):
    """Geodetic rate (with J2) against the orbit's mean radius: the prediction's sensitivity."""
    out = []
    base = dict(L.GPB_ORBIT)
    for rk in radii_km:
        L.GPB_ORBIT["a_m"] = rk * 1e3
        try:
            p = predict(n_days=n_days, step_s=20.0, lunisolar=False)
        finally:
            L.GPB_ORBIT.update(base)
        out.append({"mean_r_km": p["mean_r_km"], "geodetic_ns": p["geodetic"]["ns_mas_yr"],
                    "point_mass_ns": p["geodetic_point_mass"]["ns_mas_yr"]})
    return out


def run(force: bool = False) -> dict:
    """Full prediction over the science phase (cached)."""
    import json
    path = C.CACHE_DIR / "gpb_prediction.json"
    tracks = C.CACHE_DIR / "gpb_tracks.npz"
    if path.exists() and tracks.exists() and not force:
        return json.loads(path.read_text())
    p = predict(lunisolar=True)
    scan = radius_scan()
    summary = {k: {kk: vv for kk, vv in v.items() if not kk.endswith("track")}
               for k, v in p.items() if isinstance(v, dict)}
    summary.update({k: v for k, v in p.items() if isinstance(v, float)})
    summary["radius_scan"] = scan
    # radius that reproduces GP-B's own -6606.1 (linear in the scan)
    r = np.array([s["mean_r_km"] for s in scan])
    g = np.array([s["geodetic_ns"] for s in scan])
    summary["radius_for_published_km"] = float(np.interp(L.GPB_GR_PREDICTION["ns"], g, r))
    summary["sensitivity_mas_yr_per_km"] = float(np.polyfit(r, g, 1)[0])
    summary["total_ns"] = summary["geodetic"]["ns_mas_yr"] + summary["frame_dragging"]["ns_mas_yr"]
    summary["total_we"] = summary["geodetic"]["we_mas_yr"] + summary["frame_dragging"]["we_mas_yr"]
    path.write_text(json.dumps(summary, indent=1))
    np.savez_compressed(tracks, t_days=p["t_days"][::180],
                        **{f"{k}_{c}": p[k][f"{c}_track"] for k in
                           ("geodetic", "frame_dragging", "solar_geodetic") for c in ("ns", "we")},
                        offplane_arcsec=p["guide_star_offplane_deg"] * 3600)
    return summary
