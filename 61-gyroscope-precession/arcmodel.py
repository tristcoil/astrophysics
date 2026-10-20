"""Build the per-arc inputs of the compiled orbit model and run it.

An "arc" is a stretch of time (here a week plus the following days we propagate
into) over which we tabulate everything that does not depend on the satellite:

* Earth orientation on a 1-hour grid: the precession-nutation matrix Q
  (GCRS <- CIRS, with the C04 pole offsets), polar motion W (TIRS <- ITRS) and
  UT1 for the Earth rotation angle;
* the Sun, Moon, Venus, Mars, Jupiter and Saturn (DE440) and the de Sitter
  precession vector on a 30-minute grid;
* the tidal changes of the gravity field (solid, ocean, pole) on a 5-minute grid;
* the static field (GOCO06s at the arc's mid-epoch) to degree 30.

Times inside the model are TT seconds since the arc epoch ``mjd0`` (MJD, TT).
"""

from __future__ import annotations

import functools
from dataclasses import dataclass, field

import erfa
import numpy as np
from jplephem.spk import SPK

import config as C
import forcemodel as FM
import frames as F
import gravity
import tides

MJD0 = 2_400_000.5
ROT_STEP = 3600.0
BODY_STEP = 1800.0
TIDE_STEP = 300.0

BODY_GM = {  # m^3/s^2, DE440 header
    "sun": C.GM_SUN, "moon": C.GM_MOON,
    "venus": 3.24858592000e14, "mars": 4.28283752140e13,
    "jupiter": 1.26712764100e17, "saturn": 3.79405852000e16,
}
BODY_SEGMENTS = {"sun": 10, "venus": 2, "mars": 4, "jupiter": 5, "saturn": 6}


@functools.lru_cache(maxsize=1)
def ephemeris() -> SPK:
    path = C.EPHEM_DIR / "de440.bsp"
    if not path.exists():
        path = C.EPHEM_DIR / "de440s.bsp"
    return SPK.open(str(path))


def body_positions(mjd_tt: np.ndarray, with_velocity: bool = False):
    """Geocentric GCRS (ICRF-aligned) positions (m) of the bodies, shape (N, nb, 3)."""
    k = ephemeris()
    jd_tt = mjd_tt + MJD0
    dtdb = erfa.dtdb(MJD0, mjd_tt, 0.0, 0.0, 0.0, 0.0)
    jd_tdb = jd_tt + dtdb / C.SECONDS_PER_DAY
    p_emb, v_emb = k[0, 3].compute_and_differentiate(jd_tdb)
    p_e, v_e = k[3, 399].compute_and_differentiate(jd_tdb)
    p_m, v_m = k[3, 301].compute_and_differentiate(jd_tdb)
    earth = p_emb + p_e
    earth_v = v_emb + v_e
    out = []
    vel = []
    for name in ("sun", "moon", "venus", "mars", "jupiter", "saturn"):
        if name == "moon":
            p, v = p_m - p_e, v_m - v_e
        else:
            pb, vb = k[0, BODY_SEGMENTS[name]].compute_and_differentiate(jd_tdb)
            p, v = pb - earth, vb - earth_v
        out.append(p.T * 1e3)
        vel.append(v.T * 1e3 / C.SECONDS_PER_DAY)
    pos = np.stack(out, axis=1)
    if with_velocity:
        return pos, np.stack(vel, axis=1)
    return pos


@dataclass
class ArcModel:
    mjd0: float                     # arc epoch, MJD (TT)
    t_end: float                    # s after mjd0 covered by the grids
    M_f: np.ndarray
    M_i: np.ndarray
    Qg: np.ndarray
    Wg: np.ndarray
    tug: np.ndarray
    Cs: np.ndarray
    Ss: np.ndarray
    dCg: np.ndarray
    dSg: np.ndarray
    bodies: np.ndarray
    gmb: np.ndarray
    dsg: np.ndarray
    info: dict = field(default_factory=dict)

    def run(self, y0, prm, t_start=0.0, t_stop=None, h=C.INTEGRATION_STEP_S, stride=10,
            with_var=False, mu_lt=None, mu_ds=None):
        """Integrate from t_start to t_stop (s); returns times and states (every stride steps)."""
        t_stop = self.t_end if t_stop is None else t_stop
        nsteps = int(round(abs(t_stop - t_start) / h))      # h > 0; direction from t_stop
        M_f = self.M_f.copy()
        if mu_lt is not None:
            M_f[8] = mu_lt
        if mu_ds is not None:
            M_f[9] = mu_ds
        y = np.asarray(y0, float)
        if with_var and y.size == 6:
            y = np.concatenate([y, np.eye(6).T.ravel(), np.zeros(6 * FM.NPAR)])
        return FM.propagate(y, float(t_start), float(np.sign(t_stop - t_start) * abs(h)), nsteps,
                            stride, M_f, self.M_i, self.Qg, self.Wg, self.tug, self.Cs, self.Ss,
                            self.dCg, self.dSg, self.bodies, self.gmb, self.dsg,
                            np.asarray(prm, float), with_var)


def unnormalise(cbar: np.ndarray, N: np.ndarray) -> np.ndarray:
    n = cbar.shape[-1]
    return cbar * N[:n, :n]


def build(mjd0: float, span_days: float, sat: str, *, gravity_model: str = C.GRAVITY_MODEL,
          nmax: int = C.GRAVITY_DEGREE, ntide: int = C.OCEAN_TIDE_DEGREE,
          ocean: bool = True, solid: bool = True, pole: bool = True,
          mu_lt: float = 0.0, mu_ds: float = 1.0, t_before: float = 0.0,
          dc20_table: tuple | None = None, srp: bool = True, om1_flip: bool = False,
          field_mjd: float | None = None) -> ArcModel:
    """Tabulate the model for TT times mjd0 - t_before .. mjd0 + span_days.

    dc20_table: optional (mjd, dCbar20) arrays; the correction to the normalised
    C20 is interpolated to the arc's middle and added to the static field (used to
    let the Earth's flattening follow what the satellites themselves measure).
    """
    t_lo = -t_before - 2 * ROT_STEP
    t_hi = span_days * C.SECONDS_PER_DAY + 2 * ROT_STEP

    # --- Earth orientation grid ----------------------------------------------------
    tr = np.arange(t_lo, t_hi + ROT_STEP, ROT_STEP)
    mjd_tt = mjd0 + tr / C.SECONDS_PER_DAY
    mjd_utc = F.tt_mjd_to_utc_mjd(mjd_tt)
    e = F.eop_at(mjd_utc)
    x, y = erfa.xy06(MJD0, mjd_tt)
    x, y = x + e["dx"], y + e["dy"]
    s = erfa.s06(MJD0, mjd_tt, x, y)
    Qg = np.transpose(erfa.c2ixys(x, y, s), (0, 2, 1))               # GCRS <- CIRS
    sp = erfa.sp00(MJD0, mjd_tt)
    Wg = np.transpose(erfa.pom00(e["xp"], e["yp"], sp), (0, 2, 1))   # TIRS <- ITRS
    mjd_tai = mjd_utc + F.tai_minus_utc(mjd_utc) / C.SECONDS_PER_DAY
    tug = (mjd_tai + e["ut1_tai"] / C.SECONDS_PER_DAY) + MJD0 - 2_451_545.0

    # --- bodies and the de Sitter precession vector ------------------------------
    tb = np.arange(t_lo - 2 * BODY_STEP, t_hi + 3 * BODY_STEP, BODY_STEP)
    pos, vel = body_positions(mjd0 + tb / C.SECONDS_PER_DAY, with_velocity=True)
    sun, sun_v = pos[:, 0], vel[:, 0]
    ds = 1.5 * C.GM_SUN / (C.C_LIGHT ** 2 * np.linalg.norm(sun, axis=1) ** 3)[:, None] \
        * np.cross(sun, sun_v)
    gmb = np.array([BODY_GM[b] for b in ("sun", "moon", "venus", "mars", "jupiter", "saturn")])

    # --- static field at the arc's middle ---------------------------------------------
    fld = gravity.load(gravity_model, nmax, mjd=mjd0 + span_days / 2 if field_mjd is None else field_mjd)
    dc20 = 0.0
    if dc20_table is not None:
        dc20 = float(np.interp(mjd0 + span_days / 2, np.asarray(dc20_table[0], float),
                               np.asarray(dc20_table[1], float)))
        fld.cbar[2, 0] += dc20
    Nn = gravity.normalisation(max(nmax, ntide) + 1)
    Cs = unnormalise(fld.cbar, Nn)
    Ss = unnormalise(fld.sbar, Nn)

    # --- tides on a 5-minute grid -------------------------------------------------------
    tt = np.arange(t_lo, t_hi + TIDE_STEP, TIDE_STEP)
    mjd_tt_t = mjd0 + tt / C.SECONDS_PER_DAY
    mjd_utc_t = F.tt_mjd_to_utc_mjd(mjd_tt_t)
    et = F.eop_at(mjd_utc_t)
    mjd_ut1_t = mjd_utc_t + (F.tai_minus_utc(mjd_utc_t) + et["ut1_tai"]) / C.SECONDS_PER_DAY
    beta = tides.doodson_arguments(mjd_tt_t, mjd_ut1_t)
    dC = np.zeros((tt.size, ntide + 1, ntide + 1))
    dS = np.zeros_like(dC)
    if ocean:
        oc, os_ = tides.ocean_tide_coefficients(beta, ntide, om1_flip=om1_flip)
        dC += oc
        dS += os_
    if solid:
        R = F.itrs_to_gcrs_matrices(mjd_utc_t)
        pb = body_positions(mjd_tt_t)
        sun_i = np.einsum("nji,nj->ni", R, pb[:, 0])
        moon_i = np.einsum("nji,nj->ni", R, pb[:, 1])
        sc, ss = tides.solid_tide_coefficients(sun_i, moon_i, beta,
                                               zero_tide=fld.tide_system.startswith("zero"))
        dC[:, :5, :5] += sc
        dS[:, :5, :5] += ss
    if pole:
        pc, ps = tides.pole_tide_coefficients(mjd_utc_t, et["xp"], et["yp"])
        dC[:, 2, 1] += pc
        dS[:, 2, 1] += ps
    dCg = dC * Nn[:ntide + 1, :ntide + 1]
    dSg = dS * Nn[:ntide + 1, :ntide + 1]

    info = C.SATELLITES[sat]
    area = np.pi * (info["diameter_m"] / 2) ** 2
    M_f = np.array([fld.gm, fld.radius, C.C_LIGHT, C.EARTH_SPIN_PER_MASS,
                    area / info["mass_kg"] if srp else 0.0,
                    C.SOLAR_PRESSURE_1AU, C.AU, C.R_SUN, mu_lt, mu_ds,
                    tr[0], ROT_STEP, tt[0], TIDE_STEP, tb[0], BODY_STEP])
    M_i = np.array([nmax, ntide, 6], dtype=np.int64)
    return ArcModel(mjd0=mjd0, t_end=span_days * C.SECONDS_PER_DAY, M_f=M_f, M_i=M_i,
                    Qg=np.ascontiguousarray(Qg), Wg=np.ascontiguousarray(Wg), tug=tug,
                    Cs=np.ascontiguousarray(Cs), Ss=np.ascontiguousarray(Ss),
                    dCg=np.ascontiguousarray(dCg), dSg=np.ascontiguousarray(dSg),
                    bodies=np.ascontiguousarray(pos), gmb=gmb, dsg=np.ascontiguousarray(ds),
                    info={"gravity": gravity_model, "tide_system": fld.tide_system,
                          "sat": sat, "ocean": ocean, "solid": solid, "pole": pole,
                          "dc20": dc20, "mu_lt": mu_lt, "mu_ds": mu_ds, "om1_flip": om1_flip})
