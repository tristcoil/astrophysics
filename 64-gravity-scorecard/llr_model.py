"""Predict every APOLLO laser round trip to the Moon from DE440 and general relativity.

This follows the JPL lunar laser range model (Williams & Boggs 2020, IOM
335N-20-01; section numbers below refer to it) step by step:

  T1  the laser fires (UTC -> TAI -> TT -> TDB, with the station's own
      relativistic clock term)                                       sec. 3.7
      the station: ITRS coordinates + plate motion + solid tides + pole
      tide, rotated to GCRS (precession-nutation, UT1, polar motion with
      sub-daily ocean-tide terms), then Lorentz-contracted and rescaled to
      barycentric TDB-compatible coordinates                       sec. 3.3-3.10
  T2  the pulse bounces off a reflector: DE440 Moon + the reflector's
      principal-axis coordinates rotated by DE440's integrated librations,
      displaced by Earth- and Sun-raised tides on the Moon, Lorentz-
      contracted and rescaled                                        sec. 3.11-3.14
      up-leg light time = distance/c + atmosphere (Mendes-Pavlis + FCULa)
      + Shapiro delay of the Sun and the Earth (1 + gamma) + corner cube,
      iterated to < 1 ps                                             sec. 3.15-3.20
  T3  the pulse returns: the same for the down leg                   sec. 3.21
      round trip in TDB converted back to station TAI seconds        sec. 3.22

General relativity enters three times: in the DE440 orbits and librations (the
ephemeris is integrated with the PPN equations of motion, beta = gamma = 1, and
geodetic precession), in the Shapiro delay (~7.6 m one way), and in the time and
space transformations (~1 m effects). Not modelled (each < 1 cm at Apache Point,
W&B 2020): ocean and atmospheric loading, geocentre motion, reflector thermal
expansion, lunar tidal dissipation, the Moon's and Jupiter's Shapiro delays.
"""

from __future__ import annotations

import erfa
import numpy as np
import pandas as pd

import config as C
import earth
from ephem import EARTH, MOON, SUN, default as default_ephemeris

C_KM = C.C_KM
C2_KM = C_KM**2


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _dot(a, b):
    return np.einsum("ij,ij->i", a, b)


def _mv(m, v):
    """Stacked matrix @ vector: (N,3,3) x (N,3)."""
    return np.einsum("nij,nj->ni", m, v)


def _mtv(m, v):
    """Stacked transpose(matrix) @ vector."""
    return np.einsum("nji,nj->ni", m, v)


def tide_on_body(r_site, perturbers, gm_self, radius, love):
    """Degree-2 (+3) elastic tide displacement (eq. 10 of W&B 2020, IERS eqs. 5-6).

    r_site (N,3) site vector from the body centre; perturbers: list of
    (r_body (N,3), GM, include_degree3) in the same frame and units (m).
    """
    rs = np.linalg.norm(r_site, axis=1, keepdims=True)
    rhat = r_site / rs
    out = np.zeros_like(r_site)
    for r_b, gm, deg3 in perturbers:
        rb = np.linalg.norm(r_b, axis=1, keepdims=True)
        bhat = r_b / rb
        d = np.sum(bhat * rhat, axis=1, keepdims=True)
        f2 = gm / gm_self * radius**4 / rb**3
        out += f2 * (love["h2"] * rhat * (1.5 * d**2 - 0.5) + 3 * love["l2"] * d * (bhat - d * rhat))
        if deg3:
            f3 = gm / gm_self * radius**5 / rb**4
            out += f3 * (love["h3"] * rhat * (2.5 * d**3 - 1.5 * d)
                         + love["l3"] * (7.5 * d**2 - 1.5) * (bhat - d * rhat))
    return out


def water_vapour_hpa(rh_pct, t_c, p_hpa):
    """Water-vapour pressure from relative humidity (Giacomo 1982, as in Mendes & Pavlis 2004)."""
    tk = t_c + 273.15
    e_sat = 0.01 * np.exp(1.2378847e-5 * tk**2 - 1.9121316e-2 * tk + 33.93711047 - 6.3431645e3 / tk)
    f_w = 1.00062 + 3.14e-6 * p_hpa + 5.6e-7 * t_c**2
    return rh_pct / 100.0 * f_w * e_sat


def zenith_delay_m(lat_rad, h_m, p_hpa, e_hpa, lam_um):
    """Mendes & Pavlis (2004) zenith delay (IERS 2010 eqs. 9.3-9.7, FCULZD_HPA.F)."""
    sigma = 1.0 / lam_um
    f_s = 1 - 0.00266 * np.cos(2 * lat_rad) - 0.00028e-3 * h_m
    corr = 1.0 + 0.534e-6 * (C.CO2_PPM - 450.0)
    k0, k1, k2, k3 = 238.0185, 19990.975, 57.362, 579.55174
    f_h = 0.01 * corr * (k1 * (k0 + sigma**2) / (k0 - sigma**2) ** 2 + k3 * (k2 + sigma**2) / (k2 - sigma**2) ** 2)
    w0, w1, w2, w3 = 295.235, 2.6422, -0.032380, 0.004028
    f_nh = 0.003101 * (w0 + 3 * w1 * sigma**2 + 5 * w2 * sigma**4 + 7 * w3 * sigma**6)
    zhd = 2.416579e-3 * f_h * p_hpa / f_s
    zwd = 1e-4 * (5.316 * f_nh - 3.759 * f_h) * e_hpa / f_s
    return zhd, zwd


def fcula(lat_rad, h_m, t_c, elev_rad):
    """FCULa mapping function (Mendes et al. 2002; IERS 2010 eq. 9.9, FCUL_A.F)."""
    cphi = np.cos(lat_rad)
    a1 = 0.121008e-02 + 0.17295e-05 * t_c + 0.3191e-04 * cphi - 0.18478e-07 * h_m
    a2 = 0.304965e-02 + 0.2346e-05 * t_c - 0.1035e-03 * cphi - 0.1856e-07 * h_m
    a3 = 0.68777e-01 + 0.1972e-04 * t_c - 0.3458e-02 * cphi + 0.1060e-06 * h_m
    s = np.sin(elev_rad)
    return (1.0 + a1 / (1.0 + a2 / (1.0 + a3))) / (s + a1 / (s + a2 / (s + a3)))


def shapiro_s(gm, r1, r2, rho, gamma=C.GAMMA_PPN):
    """One-way Shapiro delay (s) of a body with GM (km^3/s^2); distances in km (W&B eq. 19)."""
    return (1 + gamma) * gm / C_KM**3 * np.log((r1 + r2 + rho) / (r1 + r2 - rho))


# ---------------------------------------------------------------------------
# Times
# ---------------------------------------------------------------------------

class Epoch:
    """Station-clock instants as UTC calendar day + seconds, with every scale we need.

    Barycentric TDB is kept as (day since J2000, seconds) for the exact ephemeris.
    """

    def __init__(self, year, month, day, utc_sec):
        self.year, self.month, self.day = (np.asarray(a, int) for a in (year, month, day))
        self.utc_sec = np.asarray(utc_sec, float)                    # seconds since 0h UTC of that date
        mjd0, mjd = erfa.cal2jd(self.year, self.month, self.day)
        self.mjd = np.asarray(mjd, float)                            # MJD of 0h of the UTC date
        self.dat = erfa.dat(self.year, self.month, self.day, np.clip(self.utc_sec / C.DAY, 0, 0.999999))
        self.tt_sec = self.utc_sec + self.dat + 32.184
        self.mjd_utc = self.mjd + self.utc_sec / C.DAY

    def shifted(self, dt_s):
        """Same calendar day, dt_s seconds later on the station clock (UTC continuous inside a day)."""
        out = object.__new__(Epoch)
        out.__dict__.update(self.__dict__)
        out.utc_sec = self.utc_sec + dt_s
        out.tt_sec = self.tt_sec + dt_s
        out.mjd_utc = self.mjd + out.utc_sec / C.DAY
        return out

    def tt_jd(self):
        return 2400000.5 + self.mjd, self.tt_sec / C.DAY

    def tdb_minus_tt(self, ut1_frac, elong, u_km, v_km):
        jd1, jd2 = self.tt_jd()
        return erfa.dtdb(jd1, jd2, ut1_frac, elong, u_km, v_km)

    def j2000(self, extra_s=0.0):
        """(day, seconds) since J2000 TT = 2000-01-01 12:00, plus extra seconds."""
        day = (self.mjd - 51545).astype(np.int64)
        return day, 43200.0 + self.tt_sec + extra_s


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

def station_itrs0(params: dict | None = None) -> np.ndarray:
    g = C.APOLLO_GEOCENTRIC
    r = earth.geocentric_to_xyz(g["radius_km"] * 1e3, g["lat_deg"], g["lon_deg"])
    if params and "station_xyz" in params:
        r = np.asarray(params["station_xyz"], float)
    return r


def compute(df: pd.DataFrame, params: dict | None = None, *, gamma: float = C.GAMMA_PPN,
            n_iter: int = 6, ephemeris=None, components: bool = True) -> pd.DataFrame:
    """Computed round trips (s, TAI) and the pieces of the model for each normal point.

    params: optional 'station_xyz' (ITRS m at C.LLR_REFERENCE_EPOCH) and 'station_vel' (m/yr).
    """
    params = params or {}
    eph = ephemeris or default_ephemeris()
    n = len(df)
    t1 = Epoch(df["year"], df["month"], df["day"],
               df["hour"] * 3600.0 + df["minute"] * 60.0 + df["utc_sec_int"] + df["utc_sec_frac"])
    rt_obs = df["rt_s"].to_numpy(float)

    r0 = station_itrs0(params)
    vel = np.asarray(params.get("station_vel", (0.0, 0.0, 0.0)), float)
    phi_s, lam_s, h_s = earth.geodetic(r0)
    e_hat, n_hat, u_hat = earth.enu_basis(r0)
    u_km = np.hypot(r0[0], r0[1]) / 1e3
    v_km = r0[2] / 1e3
    dec_year = df["year"].to_numpy(float) + (t1.mjd - erfa.cal2jd(df["year"], 1, 1)[1]) / 365.25
    refl_pa = np.array([C.REFLECTORS_PA[r] for r in df["reflector"]])              # (N,3) m
    cc_delay = np.array([C.CC_INDEX * C.CC_DEPTH[r] for r in df["reflector"]]) / C.C_LIGHT
    p_hpa = df["pressure_hpa"].to_numpy(float)
    t_c = df["temp_c"].to_numpy(float)
    e_hpa = water_vapour_hpa(df["rh_pct"].to_numpy(float), t_c, p_hpa)
    lam_um = df["wavelength_um"].to_numpy(float)
    zhd, zwd = zenith_delay_m(phi_s, h_s, p_hpa, e_hpa, lam_um)

    def station_state(ep: Epoch):
        """Station GCRS position/velocity (m, m/s), its ITRS->GCRS matrix and TDB - TT (s)."""
        eop = earth.eop_at(ep.mjd_utc)
        ut1_sec = ep.utc_sec + ep.dat + eop["ut1_tai"]
        ut1_jd1, ut1_jd2 = 2400000.5 + ep.mjd, ut1_sec / C.DAY
        tt1, tt2 = ep.tt_jd()
        rc2t = earth.c2t_matrix(tt1, tt2, ut1_jd1, ut1_jd2, eop["x"], eop["y"], eop["dX"], eop["dY"])
        tdb_tt = ep.tdb_minus_tt(np.mod(ut1_sec / C.DAY, 1.0), lam_s, u_km, v_km)
        day, sec = ep.j2000(tdb_tt + params.get("eph_dt", 0.0))
        # geocentric Moon and Sun (km, ICRF ~ GCRS axes) for the solid tides
        pm, _ = eph.posvel(MOON, day, sec, center=EARTH)
        ps, _ = eph.posvel(SUN, day, sec, center=EARTH)
        moon_t = _mv(rc2t, pm.T * 1e3)
        sun_t = _mv(rc2t, ps.T * 1e3)
        years = dec_year + (ep.utc_sec - t1.utc_sec) / C.JULIAN_YEAR
        r_itrs = np.tile(r0, (n, 1)) + np.outer(years - C.LLR_REFERENCE_EPOCH, vel)
        d_tide = earth.solid_tide(r_itrs, moon_t, sun_t)
        gmst = erfa.gmst06(ut1_jd1, ut1_jd2, tt1, tt2)
        t_cent = ((tt1 - 2451545.0) + tt2) / 36525.0
        d_k1 = earth.k1_correction(r_itrs, gmst, erfa.faom03(t_cent))
        d_pole = earth.pole_tide(r_itrs, eop["x"], eop["y"], years)
        r_full = r_itrs + d_tide + d_k1 + d_pole
        r_gcrs = _mtv(rc2t, r_full)
        v_gcrs = _mtv(rc2t, np.cross(np.array([0.0, 0.0, C.OMEGA_EARTH]), r_full))
        up_gcrs = _mtv(rc2t, np.tile(u_hat, (n, 1)))
        out = dict(day=day, sec=sec, rc2t=rc2t, r=r_gcrs, v=v_gcrs, up=up_gcrs, tdb_tt=tdb_tt,
                   tide=d_tide, k1=d_k1, pole=d_pole)
        return out

    def to_bcrs_station(st, day, sec):
        """Barycentric station position (km) at TDB (day, sec): W&B eqs. (11)-(12)."""
        pe, ve = eph.posvel(EARTH, day, sec)
        ps, _ = eph.posvel(SUN, day, sec)
        pe, ve, ps = pe.T, ve.T, ps.T
        u_e = C.GM_SUN / np.linalg.norm(pe - ps, axis=1)                 # km^2/s^2
        rs = st["r"] / 1e3
        rs_b = rs * (1 - C.L_C - u_e / C2_KM)[:, None] - (_dot(ve, rs) / (2 * C2_KM))[:, None] * ve
        return pe + rs_b, pe, ve, ps

    # --- T1 ---
    st1 = station_state(t1)
    x1, pe1, ve1, ps1 = to_bcrs_station(st1, st1["day"], st1["sec"])

    def reflector_bcrs(tau):
        """Reflector position (km, BCRS) at T2 = T1(TDB) + tau and its pieces."""
        day, sec = st1["day"], st1["sec"] + tau
        pm, vm = eph.posvel(MOON, day, sec)
        pe, _ = eph.posvel(EARTH, day, sec)
        ps, _ = eph.posvel(SUN, day, sec)
        pm, vm, pe, ps = pm.T, vm.T, pe.T, ps.T
        m_pa = eph.moon_pa_to_icrf(day, sec)                       # PA -> ICRF
        earth_pa = _mtv(m_pa, (pe - pm) * 1e3)
        sun_pa = _mtv(m_pa, (ps - pm) * 1e3)
        tide_pa = tide_on_body(refl_pa, [(earth_pa, C.GM_EARTH, True), (sun_pa, C.GM_SUN, True)],
                               C.GM_MOON, C.R_MOON, C.MOON_LOVE)
        rr = _mv(m_pa, refl_pa + tide_pa) / 1e3                     # km, ICRF
        u_m = C.GM_SUN / np.linalg.norm(pm - ps, axis=1)
        rr_b = rr * (1 - C.L_M - u_m / C2_KM)[:, None] - (_dot(vm, rr) / (2 * C2_KM))[:, None] * vm
        return pm + rr_b, dict(pm=pm, pe=pe, ps=ps, m_pa=m_pa, tide_pa=tide_pa)

    def leg_delays(x_sta, x_ref, up, pe_sta, ps_sta, pe_ref, ps_ref):
        rho = np.linalg.norm(x_ref - x_sta, axis=1)
        uvec = (x_ref - x_sta) / rho[:, None]
        sin_el = _dot(uvec, up)
        elev = np.arcsin(np.clip(sin_el, -1, 1))
        mapf = fcula(phi_s, h_s, t_c, elev)
        atm = (zhd + zwd) * mapf / C.C_LIGHT
        sh_sun = shapiro_s(C.GM_SUN, np.linalg.norm(x_sta - ps_sta, axis=1),
                           np.linalg.norm(x_ref - ps_ref, axis=1), rho, gamma)
        sh_earth = shapiro_s(C.GM_EARTH, np.linalg.norm(x_sta - pe_sta, axis=1),
                             np.linalg.norm(x_ref - pe_ref, axis=1), rho, gamma)
        return rho, uvec, elev, atm, sh_sun, sh_earth

    # --- up leg ---
    tau12 = rt_obs / 2.0
    for _ in range(n_iter):
        x2, aux2 = reflector_bcrs(tau12)
        rho12, u12, el1, atm1, shs1, she1 = leg_delays(x1, x2, st1["up"], pe1, ps1, aux2["pe"], aux2["ps"])
        tau12_new = rho12 / C_KM + atm1 + shs1 + she1 + cc_delay
        step = np.abs(tau12_new - tau12).max()
        tau12 = tau12_new
        if step < 1e-13:
            break

    # --- down leg ---
    tau23 = tau12.copy()
    for _ in range(n_iter):
        t3 = t1.shifted(tau12 + tau23)                       # station-clock seconds ~ TDB seconds (1e-8)
        st3 = station_state(t3)
        day3, sec3 = st1["day"], st1["sec"] + tau12 + tau23
        x3, pe3, ve3, ps3 = to_bcrs_station(st3, day3, sec3)
        rho23, u23, el3, atm3, shs3, she3 = leg_delays(x3, x2, st3["up"], pe3, ps3, aux2["pe"], aux2["ps"])
        tau23_new = rho23 / C_KM + atm3 + shs3 + she3 + cc_delay
        step = np.abs(tau23_new - tau23).max()
        tau23 = tau23_new
        if step < 1e-13:
            break

    rt_tdb = tau12 + tau23
    # station-clock interval: remove the change of TDB - TT between T1 and T3 (eq. 5)
    t3 = t1.shifted(rt_tdb - (st3["tdb_tt"] - st1["tdb_tt"]))
    st3_final = station_state(t3)
    rt_tai = rt_tdb - (st3_final["tdb_tt"] - st1["tdb_tt"])

    out = pd.DataFrame({
        "rt_obs": rt_obs, "rt_comp": rt_tai,
        "resid_s": rt_obs - rt_tai,
        "resid_m": (rt_obs - rt_tai) * C.C_LIGHT / 2.0,
        "elev_deg": np.degrees(0.5 * (el1 + el3)),
    })
    if components:
        e_g1 = _mtv(st1["rc2t"], np.tile(e_hat, (n, 1)))
        n_g1 = _mtv(st1["rc2t"], np.tile(n_hat, (n, 1)))
        u_g1 = st1["up"]
        e_g3 = _mtv(st3["rc2t"], np.tile(e_hat, (n, 1)))
        n_g3 = _mtv(st3["rc2t"], np.tile(n_hat, (n, 1)))
        u_g3 = st3["up"]
        # one-way range partials (m per m) for station shifts east/north/up
        out["d_east"] = -0.5 * (_dot(u12, e_g1) + _dot(u23 if False else (x2 - x3) / rho23[:, None], e_g3))
        out["d_north"] = -0.5 * (_dot(u12, n_g1) + _dot((x2 - x3) / rho23[:, None], n_g3))
        out["d_up"] = -0.5 * (_dot(u12, u_g1) + _dot((x2 - x3) / rho23[:, None], u_g3))
        # partials for reflector PA coordinates (m per m)
        u_mean = 0.5 * (u12 + (x2 - x3) / rho23[:, None])
        for k, ax in enumerate("xyz"):
            out[f"d_refl_{ax}"] = _dot(u_mean, aux2["m_pa"][:, :, k])
        out["atm_m"] = 0.5 * (atm1 + atm3) * C.C_LIGHT
        out["shapiro_sun_m"] = 0.5 * (shs1 + shs3) * C.C_LIGHT
        out["shapiro_earth_m"] = 0.5 * (she1 + she3) * C.C_LIGHT
        out["tdb_term_m"] = -(st3_final["tdb_tt"] - st1["tdb_tt"]) * C.C_LIGHT / 2.0
        out["tide_up_m"] = _dot(st1["tide"], np.tile(u_hat, (n, 1)))
        out["moon_tide_m"] = _dot(aux2["tide_pa"], refl_pa / np.linalg.norm(refl_pa, axis=1, keepdims=True))
        out["range_km"] = 0.5 * (rho12 + rho23)
        # mean elongation D of the Moon from the Sun (Delaunay argument) at the bounce
        t_cent = (st1["day"] + (st1["sec"] + tau12) / C.DAY) / 36525.0
        out["D_rad"] = np.mod(erfa.fad03(t_cent), 2 * np.pi)
        out["F_rad"] = np.mod(erfa.faf03(t_cent), 2 * np.pi)
        out["l_rad"] = np.mod(erfa.fal03(t_cent), 2 * np.pi)
        out["lp_rad"] = np.mod(erfa.falp03(t_cent), 2 * np.pi)
        out["years"] = dec_year
    return out
