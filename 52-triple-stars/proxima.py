"""Proxima Centauri's orbit around alpha Centauri A+B.

1. Kervella, Thevenin & Lovis (2017) reproduced from their own inputs: the
   Hipparcos-epoch (J1991.25) positions, parallaxes and proper motions of
   Proxima and of the A+B barycentre, their absolute radial velocities
   (corrected for gravitational redshift and convective blueshift) and the
   three masses. Positions and velocities become 3-D vectors; their
   difference and the total mass give a Keplerian orbit. 100,000 Monte Carlo
   draws of the inputs give the uncertainties, as in the paper.
2. The same with today's data: Gaia DR3 for Proxima (J2016.0) and Akeson et
   al. (2021, ALMA + Hipparcos) for the A+B barycentre (J2019.5), both
   carried to J2016.0; Akeson et al.'s masses (2.7% lighter A+B) and their
   barycentric RV, corrected for alpha Cen A's gravitational redshift the way
   Kervella et al. did.

Light-travel time: Proxima is ~36-45 light-days nearer than alpha Cen, so we
see it later in its orbit; its position is moved back by that interval
(Kervella et al.'s (+0.46", -0.09") correction).
"""

from __future__ import annotations

import numpy as np

import config as C

AU_KM = C.AU / 1e3
PC_AU = 206264.806247
KMS_PER_AUYR = AU_KM / C.YEAR          # 1 au/yr in km/s (4.7405)
MAS = np.pi / 180 / 3600e3
LIGHT_DAY_AU = C.C_LIGHT * C.DAY / C.AU
#: G in km^2 s^-2 au Msun^-1 (v_circ^2 r)
G_KMS_AU = C.GM_SUN / 1e9 / AU_KM


def _hms(s):
    h, m, sec = [float(x) for x in s.split(":")]
    return 15 * (h + m / 60 + sec / 3600)


def _dms(s):
    sign = -1 if s.strip().startswith("-") else 1
    d, m, sec = [abs(float(x)) for x in s.split(":")]
    return sign * (d + m / 60 + sec / 3600)


def _basis(ra, dec):
    """Unit vectors: line of sight (away from us), east, north."""
    ra, dec = np.radians(ra), np.radians(dec)
    u = np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)], -1)
    e = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], -1)
    n = np.stack([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)], -1)
    return u, e, n


def state(ra, dec, plx_mas, pmra, pmdec, rv):
    """Heliocentric position (au) and velocity (km/s) in ICRS axes."""
    u, e, n = _basis(ra, dec)
    d_au = 1000.0 / plx_mas * PC_AU
    vt_e = pmra * 1e-3 / (plx_mas * 1e-3) * KMS_PER_AUYR      # km/s
    vt_n = pmdec * 1e-3 / (plx_mas * 1e-3) * KMS_PER_AUYR
    r = d_au[..., None] * u
    v = rv[..., None] * u + vt_e[..., None] * e + vt_n[..., None] * n
    return r, v


def elements(r, v, mtot, frame):
    """Keplerian elements of the relative orbit. `frame` = (u, e, n) at the
    barycentre: elements are referred to the plane of the sky there, in the
    visual-binary convention (x = north, y = east, z = away from us; i < 90
    deg for counter-clockwise motion on the sky; Omega = position angle of the
    node where the companion recedes; omega from that node)."""
    u, e, n = frame
    X = np.stack([np.sum(r * n, -1), np.sum(r * e, -1), np.sum(r * u, -1)], -1)
    V = np.stack([np.sum(v * n, -1), np.sum(v * e, -1), np.sum(v * u, -1)], -1)
    # (north, east, away) is left-handed; flip y to make it right-handed for the
    # vector algebra, then translate angles back (Omega measured N -> E)
    Xr = X * np.array([1, -1, 1])
    Vr = V * np.array([1, -1, 1])
    mu = G_KMS_AU * mtot                         # km^2 s^-2 au
    rr = np.linalg.norm(Xr, axis=-1)
    vv2 = np.sum(Vr**2, -1)
    a = 1.0 / (2.0 / rr - vv2 / mu)              # au
    h = np.cross(Xr, Vr)
    hn = np.linalg.norm(h, axis=-1)
    evec = np.cross(Vr, h) / mu[..., None] - Xr / rr[..., None]
    ecc = np.linalg.norm(evec, axis=-1)
    inc = np.degrees(np.arccos(np.clip(h[..., 2] / hn, -1, 1)))
    node = np.cross(np.array([0.0, 0.0, 1.0]), h)
    Om_r = np.arctan2(node[..., 1], node[..., 0])
    nn = np.linalg.norm(node, axis=-1)
    cosw = np.sum(node * evec, -1) / (nn * ecc)
    w = np.degrees(np.arccos(np.clip(cosw, -1, 1)))
    w = np.where(evec[..., 2] < 0, 360 - w, w)
    # back to the sky convention: y (east) was flipped, so position angles
    # measured from north toward east are -Omega_r; the motion sense flips too
    Omega = np.degrees(-Om_r) % 360
    inc_sky = 180.0 - inc
    # mean anomaly now and time to periastron
    cosE = (1 - rr / a) / ecc
    sinE = np.sum(Xr * Vr, -1) / (ecc * np.sqrt(mu * a))
    E = np.arctan2(sinE, cosE)
    M = E - ecc * np.sin(E)
    P_s = 2 * np.pi * np.sqrt((a * AU_KM) ** 3 / (mu * AU_KM))  # seconds (mu km^3/s^2)
    P_yr = P_s / C.YEAR
    t_peri_yr = ((-M) % (2 * np.pi)) / (2 * np.pi) * P_yr         # time until next periastron
    return dict(a=a, e=ecc, i=inc_sky, Omega=Omega, omega=w, P_yr=P_yr, d=rr,
                v=np.sqrt(vv2) * 1e3, peri=a * (1 - ecc), apo=a * (1 + ecc),
                t_to_peri_yr=t_peri_yr, bound=a > 0, v_esc=np.sqrt(2 * mu / rr) * 1e3)


def _draw(rng, val, err, n):
    return val + err * rng.standard_normal(n)


def kervella2017(n=100_000, seed=C.RANDOM_SEED) -> dict:
    """Monte Carlo replay of Kervella et al. (2017) from their Tables 1-2."""
    k = C.SYSTEMS["alpha Cen"]["k17"]
    rng = np.random.default_rng(seed)
    AB, Pr = k["AB"], k["Prox"]
    ra_ab, de_ab = _hms(AB["ra"]), _dms(AB["dec"])
    ra_p, de_p = _hms(Pr["ra"]), _dms(Pr["dec"])
    cosd_ab, cosd_p = np.cos(np.radians(de_ab)), np.cos(np.radians(de_p))
    # light-time: Proxima seen ~44.8 days later in its motion: move it back
    dlt_days = (1000 / AB["plx"] - 1000 / Pr["plx"]) * PC_AU / LIGHT_DAY_AU
    dlt_yr = dlt_days / 365.25
    ra_p_c = ra_p - Pr["pmra"] * dlt_yr / 3.6e6 / cosd_p
    de_p_c = de_p - Pr["pmdec"] * dlt_yr / 3.6e6
    out = {}
    for label, nn in (("nominal", 1), ("mc", n)):
        def g(v, e):
            return np.full(nn, v) if nn == 1 else _draw(rng, v, e, nn)
        raA = g(ra_ab, AB["ra_err_mas"] / 3.6e6 / cosd_ab)
        deA = g(de_ab, AB["dec_err_mas"] / 3.6e6)
        raP = g(ra_p_c, Pr["ra_err_mas"] / 3.6e6 / cosd_p)
        deP = g(de_p_c, Pr["dec_err_mas"] / 3.6e6)
        rA, vA = state(raA, deA, g(AB["plx"], AB["plx_err"]), g(AB["pmra"], AB["pmra_err"]),
                       g(AB["pmdec"], AB["pmdec_err"]), g(AB["rv"], AB["rv_err"]))
        rP, vP = state(raP, deP, g(Pr["plx"], Pr["plx_err"]), g(Pr["pmra"], Pr["pmra_err"]),
                       g(Pr["pmdec"], Pr["pmdec_err"]), g(Pr["rv"], Pr["rv_err"]))
        mtot = g(k["m_A"], k["m_A_err"]) + g(k["m_B"], k["m_B_err"]) + g(k["m_C"], k["m_C_err"])
        frame = _basis(raA, deA)
        out[label] = elements(rP - rA, vP - vA, mtot, frame)
    out["light_time_days"] = float(dlt_days)
    out["light_time_shift_arcsec"] = (float(-Pr["pmra"] * dlt_yr / 1e3), float(-Pr["pmdec"] * dlt_yr / 1e3))
    return out


def update_gaia_dr3(n=100_000, seed=C.RANDOM_SEED + 1, epoch=2016.0) -> dict:
    """Gaia DR3 (Proxima) + Akeson et al. 2021 (alpha Cen A+B), at J2016.0."""
    import pandas as pd
    g3 = pd.read_csv(C.RAW / "gaia_dr3_proxima.csv").iloc[0]
    a21 = C.SYSTEMS["alpha Cen"]["a21"]
    k17 = C.SYSTEMS["alpha Cen"]["k17"]
    dv = C.SYSTEMS["alpha Cen"]["dv_grs_A"] / 1e3
    dv_err = C.SYSTEMS["alpha Cen"]["dv_grs_A_err"] / 1e3
    rng = np.random.default_rng(seed)
    # A+B barycentre: carry Akeson's J2019.5 position back to J2016.0
    dt = epoch - a21["epoch"]
    cosd = np.cos(np.radians(a21["dec_deg"]))
    ra_ab = a21["ra_deg"] + a21["pmra"] * dt / 3.6e6 / cosd
    de_ab = a21["dec_deg"] + a21["pmdec"] * dt / 3.6e6
    ra_p, de_p = float(g3["ra"]), float(g3["dec"])
    cosd_p = np.cos(np.radians(de_p))
    dlt_days = (1000 / a21["plx"] - 1000 / g3["parallax"]) * PC_AU / LIGHT_DAY_AU
    dlt_yr = dlt_days / 365.25
    ra_p_c = ra_p - g3["pmra"] * dlt_yr / 3.6e6 / cosd_p
    de_p_c = de_p - g3["pmdec"] * dlt_yr / 3.6e6
    rv_ab = a21["V0"] - dv
    rv_ab_err = np.hypot(a21["V0_err"], dv_err)
    out = {}
    for label, nn in (("nominal", 1), ("mc", n)):
        def g(v, e):
            return np.full(nn, v) if nn == 1 else _draw(rng, v, e, nn)
        # position errors: Akeson's plus the PM error carried over 3.5 yr
        e_ra = np.hypot(a21["ra_err_mas"], a21["pmra_err"] * abs(dt))
        e_de = np.hypot(a21["dec_err_mas"], a21["pmdec_err"] * abs(dt))
        raA = g(ra_ab, e_ra / 3.6e6 / cosd)
        deA = g(de_ab, e_de / 3.6e6)
        raP = g(ra_p_c, g3["ra_error"] / 3.6e6 / cosd_p)
        deP = g(de_p_c, g3["dec_error"] / 3.6e6)
        rA, vA = state(raA, deA, g(a21["plx"], a21["plx_err"]), g(a21["pmra"], a21["pmra_err"]),
                       g(a21["pmdec"], a21["pmdec_err"]), g(rv_ab, rv_ab_err))
        rP, vP = state(raP, deP, g(g3["parallax"], g3["parallax_error"]), g(g3["pmra"], g3["pmra_error"]),
                       g(g3["pmdec"], g3["pmdec_error"]), g(k17["Prox"]["rv"], k17["Prox"]["rv_err"]))
        mtot = (g(a21["m_A"], a21["m_A_err"]) + g(a21["m_B"], a21["m_B_err"])
                + g(k17["m_C"], k17["m_C_err"]))
        out[label] = elements(rP - rA, vP - vA, mtot, _basis(raA, deA))
    out["light_time_days"] = float(dlt_days)
    out["rv_ab_abs"] = float(rv_ab)
    return out


def summarize(res: dict) -> dict:
    mc = res["mc"]
    b = mc["bound"]
    s = {"p_bound": float(np.mean(b))}
    for k in ("a", "e", "i", "Omega", "omega", "P_yr", "d", "v", "v_esc", "peri", "apo", "t_to_peri_yr"):
        x = mc[k][b] if k not in ("d", "v", "v_esc") else mc[k]
        q = np.percentile(x, [16, 50, 84])
        s[k] = dict(median=float(q[1]), lo=float(q[1] - q[0]), hi=float(q[2] - q[1]),
                    nominal=float(np.atleast_1d(res["nominal"][k])[0]))
    return s


if __name__ == "__main__":
    import json
    for name, f in (("Kervella 2017 inputs", kervella2017), ("Gaia DR3 + Akeson 2021", update_gaia_dr3)):
        r = f()
        print(name, json.dumps(summarize(r), indent=0)[:1500])
