"""Cassini's two-way Doppler during the 2002 solar conjunction, modelled from scratch.

Every second, a Deep Space Network antenna sent Cassini a radio tone of
known frequency f_T (X band, 7.175 GHz, or Ka band, 34.3 GHz); the
spacecraft's transponder multiplied it by a fixed ratio M2 (880/749 at X
band; 14/15 for the Ka-band translator, KaTS) and sent it back, and the
antenna measured the returning frequency f_R. This module predicts f_R from
first principles and forms the residual fractional frequency shift

    y = (f_R,model - f_R,observed) / f_R,model,

with the model's Shapiro term switched off, so that y is what relativity
(and the noise) must explain: y = d(round-trip Shapiro delay)/dt in GR.

**Data.** 1-s closed-loop Doppler observables from the Orbit Data Files (ODF;
``cassini_odf``) and, where the ODFs lack them (the Ka/Ka passes of June 14
and 17), from the Archival Tracking Data Files (``cassini_tdf``). The
received sky frequency is f_R = K f_ref - F (K = 880/221 at X band, 3344/221
at Ka band). The transmitted frequency comes from the ramp tables.

**Model** (all in the Barycentric Celestial Reference System and TDB):

* Ephemerides: JPL Horizons' SSB vectors of the Earth, the Sun and Cassini
  (Cassini: the navigation team's reconstruction; or their SCE1 short-arc
  SPK, ``cassini_spk``), Hermite-interpolated on the exact 2-minute grid
  (the CSV's Julian dates are rounded to 0.86 ms, which matters).
* Stations: DSN 810-005 positions (ITRF93), displaced by the solid-Earth
  tides (IERS 2010, up to 0.3 m twice a day), rotated to the GCRS by astropy
  with IERS Earth orientation, then x_B = x_E + (1 - U/c^2 - L_C) x_G -
  (v_E . x_G) v_E / 2c^2 (Moyer 2000, eq. 4-3).
* Time: ODF tags are UTC; TT = UTC + 64.184 s (TAI - UTC = 32 s all of
  2002); TDB - TT is ERFA's geocentric series plus the topocentric term
  v_E . x_G / c^2. Station clocks tick at dTT/dTDB = 1 - d(TDB - TT)/dt, so
  the received frequency carries the factor (dtau/dt at t1) / (dtau/dt at t3).
* Light time: the signal received at t3 left Cassini at t2 and the station
  at t1, from c (t3 - t2) = |x_S(t3) - x_C(t2)| + c D_down and likewise up,
  where the leg delays D hold the Shapiro delay of the Sun (1PN, plus the
  2PN terms of Klioner & Zschocke 2010, with the Sun at the moment of the
  ray's closest approach), the troposphere (AMC or DSN calibrations, Niell
  mapping) and the Earth's ionosphere (DSN calibrations).
* Doppler: differentiating the light-time equations exactly,

      dt2/dt3 = (1 - n_d.v_S3/c - dD_d/dt3) / (1 - n_d.v_C/c)
      dt1/dt2 = (1 - n_u.v_C/c - dD_u/dt2) / (1 - n_u.v_S1/c)

  with n_u (station to Cassini) and n_d (Cassini to station) unit vectors.
  f_R = M2 f_T(t1) (dt1/dt3) (dtau/dt)(t1) / (dtau/dt)(t3). The delay rates
  come from full light-time solutions 30 s either side.

Because the Shapiro delay enters linearly, the model is evaluated once and
the Doppler formed twice, with the Shapiro rate on (GR) and off; GR's
predicted signal is y_GR = 1 - F(on)/F(off). Fitting y = lambda y_GR +
nuisance gives lambda = (1 + gamma)/2 (``cassini_gamma``).

Not modelled here (handled in the fit or by link combinations): the solar
plasma, errors of the reconstructed trajectory, and residual troposphere.
Neglected (each ~1e-14 or less, mostly diurnal): ocean-tide loading
(~1 cm at Goldstone), sub-daily Earth-orientation variations, antenna
thermal deformation, the ITRF93-to-current-frame difference (cm), the
planets' and the Earth's own Shapiro delays.

The pipeline: ``load_observables`` -> ``residuals`` (1-s model, ~2 min,
cached) -> ``normal_points`` (cleaned 60-s points of every link and of the
X/X + Ka/Ka combination, cached) -> ``cassini_gamma.run`` (fits, variants,
outputs). ``run_cassini()`` runs it all; ``inventory()`` lists the data.
"""

from __future__ import annotations

import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import cassini_media as MEDIA
import cassini_odf as O
import cassini_tdf as TDF
import config as C

warnings.filterwarnings("ignore", module="erfa")

CASSINI_DIR = C.RAW_DIR / "cassini_sce1"
OBS_CACHE = C.CACHE_DIR / "cassini_obs_1s.parquet"
RAMP_CACHE = C.CACHE_DIR / "cassini_ramps.parquet"
RESID_CACHE = C.CACHE_DIR / "cassini_resid_1s.parquet"
INVENTORY_PATH = C.CACHE_DIR / "cassini_inventory.csv"

C_LIGHT = C.C_LIGHT
GM_SUN = C.GM_SUN
#: IAU 2000 L_C = 1.48082686741e-8 (Moyer 2000 eq. 4-3: station coordinates
#: in the GCRS (TT units) to the BCRS (TDB units))
L_C = 1.48082686741e-8
#: TAI - UTC through all of 2002, and TT - TAI
TAI_UTC = 32.0
TT_TAI = 32.184
#: the time origin: 2002-06-04 00:00 TDB = the first Horizons epoch
T0_JD = 2452429.5
T0_1950 = (T0_JD - 2433282.5) * 86400.0        # 1,654,300,800 s
#: spacing of the Horizons vectors (s)
HORIZONS_STEP_S = 120.0

#: DSN antenna positions, ITRF93 epoch 2003.0 (DSN 810-005 module 301, rev. N,
#: Table 2; +-3 cm; the 0.5-yr plate motion to 2002.5 is ~1 cm). ODF station 15
#: is the table's DSS 152 (34-m HEF), 45 is DSS 452, 34 is DSS 341.
DSN_XYZ = {
    14: (-2353621.420, -4641341.472, 3677052.318),
    15: (-2353538.958, -4641649.429, 3676669.984),
    24: (-2354906.711, -4646840.095, 3669242.325),
    25: (-2355022.014, -4646953.204, 3669040.567),
    26: (-2354890.797, -4647166.328, 3668871.755),
    34: (-4461147.093, 2682439.239, -3674393.133),
    43: (-4460894.917, 2682361.507, -3674748.152),
    45: (-4460935.578, 2682765.661, -3674380.982),
    54: (4849434.488, -360723.8999, 4114618.835),
    55: (4849525.256, -360606.0932, 4114495.084),
    63: (4849092.518, -360180.3480, 4115109.251),
    65: (4849339.634, -360427.6637, 4114750.733),
}

#: Transponder turnaround ratios (uplink band, downlink band) -> M2. The Deep
#: Space Transponder turns X around at 880/749 (and makes a Ka downlink from
#: the X uplink at 3344/749, which the closed-loop data never recorded); the
#: Ka-band translator (KaTS) turns 34.3 GHz around at 14/15.
M2 = {(2, 2): 880 / 749, (2, 3): 3344 / 749, (3, 3): 14 / 15}
#: Observable reference: f_R = K x (record reference frequency) - observable.
KREF = {2: 880 / 221, 3: 3344 / 221}
BAND_NAME = {1: "S", 2: "X", 3: "Ka"}
#: sky-level uplink frequency ranges of each band (Hz)
BAND_RANGE = {2: (7.0e9, 7.3e9), 3: (34.0e9, 34.8e9)}


def link_name(ul: int, dl: int) -> str:
    return f"{BAND_NAME[ul]}/{BAND_NAME[dl]}"


# ---------------------------------------------------------------------------
# 1. Data: 1-s two-way Doppler observables and the ramp tables
# ---------------------------------------------------------------------------

def _odf_paths() -> list[str]:
    return sorted(glob.glob(str(CASSINI_DIR / "odf" / "*.odf")))


def load_ramps(refresh: bool = False) -> pd.DataFrame:
    """Sky-level transmitter frequency segments of every station (ODF ramp
    groups of all files, duplicates removed)."""
    if RAMP_CACHE.exists() and not refresh:
        return pd.read_parquet(RAMP_CACHE)
    rp = []
    for p in _odf_paths():
        r = O.read_odf(p)["ramps"]
        if len(r):
            rp.append(r.assign(file=Path(p).name))
    rp = pd.concat(rp, ignore_index=True)
    rp = rp[rp["sky_level"]].drop_duplicates(subset=["station", "t_start1950", "t_end1950", "f_start_hz", "rate_hz_s"])
    rp = rp.sort_values(["station", "t_start1950"]).reset_index(drop=True)
    rp.to_parquet(RAMP_CACHE)
    return rp


def load_observables(refresh: bool = False) -> pd.DataFrame:
    """Every good 1-s two-way Doppler point (all stations and bands).

    Columns: t1950 (UTC s since 1950, mid-count), station, ul_band, dl_band,
    observable (Hz), ref_freq_hz, f_r (received sky frequency, Hz), source.
    Ka/Ka points come from the ODFs; the TDF adds the Ka-channel points the
    ODFs lack (all of them Ka/Ka, identified by their frequency)."""
    if OBS_CACHE.exists() and not refresh:
        return pd.read_parquet(OBS_CACHE)
    rows = []
    for p in _odf_paths():
        d = O.read_odf(p)["doppler"]
        d = d[(d["dtype"] == 12) & (d["bad"] == 0) & (d["tc_s"] == 1.0)]
        rows.append(pd.DataFrame({
            "t1950": d["t1950"].to_numpy(), "station": d["rcv"].to_numpy(), "ul_band": d["ul_band"].to_numpy(),
            "dl_band": d["dl_band"].to_numpy(), "observable": d["observable"].to_numpy(),
            "ref_freq_hz": d["ref_freq_hz"].to_numpy(), "source": "odf"}))
    odf = pd.concat(rows, ignore_index=True)
    odf = odf.drop_duplicates(subset=["t1950", "station", "ul_band", "dl_band"])
    extra = []
    for p in sorted(glob.glob(str(CASSINI_DIR / "tdf" / "*xmmm*.tdf"))):
        o = TDF.observables(TDF.read_tdf(p, stations=(25,)))
        o = o[(o["dl_band"] == 3) & (o["mode"] == 2)]
        if len(o):
            extra.append(pd.DataFrame({
                "t1950": o["t1950"].to_numpy(), "station": 25, "ul_band": 3, "dl_band": 3,
                "observable": o["observable"].to_numpy(), "ref_freq_hz": o["ref_freq_hz"].to_numpy(),
                "source": "tdf"}))
    if extra:
        tdf = pd.concat(extra, ignore_index=True).drop_duplicates(subset=["t1950"])
        have = set(odf.loc[(odf["station"] == 25) & (odf["dl_band"] == 3), "t1950"].round(1))
        tdf = tdf[~tdf["t1950"].round(1).isin(have)]
        odf = pd.concat([odf, tdf], ignore_index=True)
    odf["f_r"] = odf["dl_band"].map(KREF) * odf["ref_freq_hz"] - odf["observable"]
    odf = odf.sort_values(["station", "ul_band", "dl_band", "t1950"]).reset_index(drop=True)
    odf.to_parquet(OBS_CACHE)
    return odf


def transmitted_frequency(rp: pd.DataFrame, station: int, band: int,
                          t1950: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sky-level uplink frequency at station times t (UTC s since 1950), NaN
    outside every ramp segment, and a flag set where the uplink was being
    swept (acquisition ramps of ~150 Hz/s): those points are later dropped."""
    lo, hi = BAND_RANGE[band]
    rr = rp[(rp["station"] == station) & rp["f_start_hz"].between(lo, hi)].sort_values("t_start1950")
    t = np.asarray(t1950, dtype=float)
    if not len(rr):
        return np.full(t.shape, np.nan), np.zeros(t.shape, bool)
    ts, te = rr["t_start1950"].to_numpy(), rr["t_end1950"].to_numpy()
    i = np.clip(np.searchsorted(ts, t, side="right") - 1, 0, len(rr) - 1)
    rate = rr["rate_hz_s"].to_numpy()[i]
    f = rr["f_start_hz"].to_numpy()[i] + rate * (t - ts[i])
    inside = (t >= ts[i]) & (t <= te[i])
    return np.where(inside, f, np.nan), rate != 0.0


# ---------------------------------------------------------------------------
# 2. Ephemerides
# ---------------------------------------------------------------------------

class Ephemeris:
    """Horizons' SSB vectors (m, m/s); t = TDB s since T0 = 2002-06-04 00:00
    TDB. Bodies: earth, cassini, sun. Positions: cubic Hermite through the
    2-minute nodes; velocities: a cubic spline through the nodes' velocities.
    The CSV prints Cassini's coordinates (1e9 km) to 1 cm, and the Hermite
    derivative turns that rounding into velocity errors of ~5e-5 m/s inside
    each 2-minute interval (3e-13 in y; ~4e-14 in a 600-s average), whereas
    the printed velocities are good to 1e-11 m/s."""

    def __init__(self, cassini: str = "horizons"):
        from scipy.interpolate import CubicHermiteSpline, CubicSpline

        rd = dict(float_precision="round_trip")
        eph = pd.read_csv(C.HORIZONS_DIR / "cassini_earth_ssb_2002.csv", **rd).merge(
            pd.read_csv(C.HORIZONS_DIR / "sun_ssb_2002.csv", **rd), on="jd_tdb")
        if abs(eph["jd_tdb"].iloc[0] - T0_JD) > 1e-9:
            raise ValueError("Horizons file does not start at T0")
        # The CSV's Julian dates carry 8 decimals (0.86 ms): at 30 km/s that
        # misplaces the Earth by up to 13 m at a node, and a Hermite spline on
        # those times makes 2-minute Doppler wiggles of ~3e-10. The grid is
        # exactly 120 s, so rebuild the times from the row number.
        ts = (eph["jd_tdb"].to_numpy() - T0_JD) * C.DAY
        exact = np.round(ts / HORIZONS_STEP_S) * HORIZONS_STEP_S
        if np.max(np.abs(exact - ts)) > 1e-3 or np.any(np.diff(exact) != HORIZONS_STEP_S):
            raise ValueError("Horizons epochs are not on the 2-minute grid")
        ts = exact
        self.span = (ts[0], ts[-1])
        self.s, self.v = {}, {}
        for b in ("earth", "cassini", "sun"):
            p = eph[[f"{b}_{k}" for k in "xyz"]].to_numpy() * 1e3
            v = eph[[f"{b}_v{k}" for k in "xyz"]].to_numpy() * 1e3
            self.s[b] = CubicHermiteSpline(ts, p, v, axis=0)
            self.v[b] = CubicSpline(ts, v, axis=0)
        self.cassini_source = cassini
        if cassini == "nav":
            import cassini_spk

            self.nav = cassini_spk.NavTrajectory()

    def pos(self, body: str, t) -> np.ndarray:
        if body == "cassini" and self.cassini_source == "nav":
            return self.nav.pos(t) + self.s["sun"](t)
        return self.s[body](t)

    def vel(self, body: str, t) -> np.ndarray:
        if body == "cassini" and self.cassini_source == "nav":
            return self.nav.vel(t) + self.v["sun"](t)
        return self.v[body](t)


# ---------------------------------------------------------------------------
# 3. Stations: positions, local vertical and clocks on a 60-s grid
# ---------------------------------------------------------------------------

#: IERS Conventions 2010 (ch. 7.1.1, step 1) nominal Love and Shida numbers:
#: degree 2 with their latitude dependence, degree 3
LOVE_H0, LOVE_H2, SHIDA_L0, SHIDA_L2 = 0.6078, -0.0006, 0.0847, 0.0002
LOVE_H3, SHIDA_L3 = 0.292, 0.015
GM_EARTH, GM_MOON, R_EARTH_EQ = 3.986004418e14, 4.9028000661e12, 6378136.6


def solid_tide(x_sta: np.ndarray, x_moon: np.ndarray, x_sun: np.ndarray) -> np.ndarray:
    """Displacement (m) of a station by the solid-Earth tides of the Moon and
    the Sun: degree-2 terms with latitude-dependent Love/Shida numbers and
    the Moon's degree-3 term (IERS Conventions 2010, eqs. 7.5-7.6, step 1;
    the frequency-dependent step-2 corrections, < 1.5 cm at the K1 period,
    are left out). All vectors geocentric, in one frame. The station moves
    by up to ~0.3 m twice a day, ~4e-5 m/s, so a two-way Doppler effect of
    ~2e-13 that changes from day to day with the Moon."""
    r = np.linalg.norm(x_sta, axis=1)
    rh = x_sta / r[:, None]
    sinphi2 = rh[:, 2] ** 2
    p2 = (3 * sinphi2 - 1) / 2
    h2 = LOVE_H0 + LOVE_H2 * p2
    l2 = SHIDA_L0 + SHIDA_L2 * p2
    dr = np.zeros_like(x_sta)
    for xb, gm in ((x_moon, GM_MOON), (x_sun, GM_SUN)):
        R = np.linalg.norm(xb, axis=1)
        Rh = xb / R[:, None]
        c = np.sum(Rh * rh, axis=1)
        f2 = gm * R_EARTH_EQ**4 / (GM_EARTH * R**3)
        dr += f2[:, None] * (h2[:, None] * rh * (1.5 * c**2 - 0.5)[:, None]
                             + 3 * l2[:, None] * c[:, None] * (Rh - c[:, None] * rh))
        f3 = gm * R_EARTH_EQ**5 / (GM_EARTH * R**4)
        dr += f3[:, None] * (LOVE_H3 * rh * (2.5 * c**3 - 1.5 * c)[:, None]
                             + SHIDA_L3 * (7.5 * c**2 - 1.5)[:, None] * (Rh - c[:, None] * rh))
    return dr


class Station:
    """A DSN antenna in the BCRS. Positions are computed by astropy (IERS
    Earth orientation) on a 60-s UTC grid, displaced by the solid-Earth
    tides, and spline-interpolated in TDB."""

    GRID_S = 60.0

    def __init__(self, sta: int, eph: Ephemeris, t_start1950: float, t_stop1950: float, tides: bool = True):
        import astropy.units as u
        import erfa
        from astropy.coordinates import EarthLocation
        from astropy.time import Time
        from scipy.interpolate import CubicSpline

        self.sta = sta
        loc = EarthLocation.from_geocentric(*DSN_XYZ[sta], unit=u.m)
        self.lon_deg, self.lat_deg = loc.lon.deg, loc.lat.deg
        self.height_m = loc.height.to(u.m).value
        grid = np.arange(np.floor(t_start1950 / 60) * 60 - 600, t_stop1950 + 660, self.GRID_S)
        tu = Time(O.utc_from_1950(grid), scale="utc")
        p, v = loc.get_gcrs_posvel(tu)
        xg = p.xyz.to(u.m).value.T
        hi = EarthLocation.from_geodetic(loc.lon, loc.lat, loc.height + 1000 * u.m)
        ph, _ = hi.get_gcrs_posvel(tu)
        up = (ph.xyz.to(u.m).value.T - xg)
        up /= np.linalg.norm(up, axis=1)[:, None]
        tt = tu.tt
        dtdb_geo = erfa.dtdb(tt.jd1, tt.jd2, 0.0, 0.0, 0.0, 0.0)
        t_tdb = grid + TAI_UTC + TT_TAI + dtdb_geo - T0_1950
        if tides:
            from astropy.coordinates import get_body_barycentric

            moon = get_body_barycentric("moon", tu, ephemeris="builtin").xyz.to(u.m).value.T
            x_moon = moon - eph.pos("earth", t_tdb)
            x_sun = eph.pos("sun", t_tdb) - eph.pos("earth", t_tdb)
            self.tide = solid_tide(xg, x_moon, x_sun)
            xg = xg + self.tide
        ve = eph.vel("earth", t_tdb)
        topo = np.sum(ve * xg, axis=1) / C_LIGHT**2
        t_tdb = t_tdb + topo
        self.tdb_grid, self.utc_grid = t_tdb, grid
        self._xg = CubicSpline(t_tdb, xg, axis=0)
        self._up = CubicSpline(t_tdb, up, axis=0)
        self._dtdb = CubicSpline(t_tdb, dtdb_geo + topo)
        self._utc = CubicSpline(t_tdb, grid)
        self.eph = eph

    def utc(self, t_tdb) -> np.ndarray:
        return self._utc(t_tdb)

    def tdb_from_utc(self, t1950) -> np.ndarray:
        return np.interp(t1950, self.utc_grid, self.tdb_grid)

    def bcrs(self, t_tdb) -> tuple[np.ndarray, np.ndarray]:
        """Station BCRS position (m) and velocity (m/s, per TDB second)."""
        xg, vg = self._xg(t_tdb), self._xg(t_tdb, 1)
        xe, ve = self.eph.pos("earth", t_tdb), self.eph.vel("earth", t_tdb)
        r_sun = np.linalg.norm(xe - self.eph.pos("sun", t_tdb), axis=1)
        scale = 1 - (GM_SUN / r_sun) / C_LIGHT**2 - L_C
        vx = np.sum(ve * xg, axis=1)[:, None]
        x = xe + scale[:, None] * xg - vx * ve / (2 * C_LIGHT**2)
        v = ve + scale[:, None] * vg - np.sum(ve * vg, axis=1)[:, None] * ve / (2 * C_LIGHT**2)
        return x, v

    def up(self, t_tdb) -> np.ndarray:
        u_ = self._up(t_tdb)
        return u_ / np.linalg.norm(u_, axis=1)[:, None]

    def clock_rate(self, t_tdb) -> np.ndarray:
        """d(TDB - TT)/dTDB at the station (the station clock keeps TT)."""
        return self._dtdb(t_tdb, 1)


# ---------------------------------------------------------------------------
# 4. Delays along each leg
# ---------------------------------------------------------------------------

def shapiro_delays(x_from: np.ndarray, x_to: np.ndarray, t_from: np.ndarray, t_to: np.ndarray,
                   eph: Ephemeris) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The Sun's gravitational delay (s) on a leg between BCRS points, with
    the Sun where it was when the ray passed closest (Bertotti, Ashby & Iess
    2008). Returns (1PN delay for gamma = 1, 2PN delay in GR, impact
    parameter b in m). 1PN: (1 + gamma) GM/c^3 ln[(r1 + r2 + R)/(r1 + r2 - R)],
    written in the stable form (r1 + r2 + R)^2 / (r1 r2 |n1 + n2|^2); 2PN
    (Klioner & Zschocke 2010, beta = gamma = epsilon = 1): (GM/c^2)^2 R /
    (c r1 r2) [(15/4) arccos(n1.n2)/|n1 x n2| - 4/(1 + n1.n2)], whose second
    ("enhanced") part is ~6 ns per leg at 1.6 solar radii."""
    d = x_to - x_from
    sun0 = eph.pos("sun", 0.5 * (t_from + t_to))
    s = np.clip(-np.sum((x_from - sun0) * d, axis=1) / np.sum(d * d, axis=1), 0.0, 1.0)
    sun = eph.pos("sun", t_from + s * (t_to - t_from))
    r1v, r2v = x_from - sun, x_to - sun
    r1, r2 = np.linalg.norm(r1v, axis=1), np.linalg.norm(r2v, axis=1)
    n1, n2 = r1v / r1[:, None], r2v / r2[:, None]
    R = np.linalg.norm(d, axis=1)
    s12 = np.sum((n1 + n2) ** 2, axis=1)                  # = 2 (1 + n1.n2)
    m = GM_SUN / C_LIGHT**2
    d1 = 2 * m / C_LIGHT * np.log((r1 + r2 + R) ** 2 / (r1 * r2 * s12))
    cross = np.linalg.norm(np.cross(n1, n2), axis=1)
    ang = np.arctan2(cross, np.sum(n1 * n2, axis=1))
    d2 = m**2 * R / (C_LIGHT * r1 * r2) * (3.75 * ang / np.maximum(cross, 1e-300) - 8.0 / s12)
    b = np.linalg.norm(r1v + s[:, None] * d, axis=1)
    return d1, d2, b


def _elevation(up: np.ndarray, x_sta: np.ndarray, v_sta: np.ndarray, x_sc: np.ndarray) -> np.ndarray:
    """Apparent elevation of the spacecraft at the station (aberration of the
    station's velocity included; refraction is in the mapping function)."""
    n = x_sc - x_sta
    n /= np.linalg.norm(n, axis=1)[:, None]
    n = n + v_sta / C_LIGHT
    n /= np.linalg.norm(n, axis=1)[:, None]
    return np.arcsin(np.clip(np.sum(n * up, axis=1), -1.0, 1.0))


class Media:
    """Troposphere and ionosphere delays for the stations in use."""

    def __init__(self, stations: dict, tropo: str = "amc", iono: bool = True, amc_unit: int = 1):
        self.tropo = {s: MEDIA.Troposphere(s, st.lat_deg, st.height_m, tropo, amc_unit=amc_unit)
                      for s, st in stations.items()} if tropo != "none" else {}
        self.iono = MEDIA.Ionosphere() if iono else None

    def delay(self, sta: int, t1950: np.ndarray, el: np.ndarray, freq: np.ndarray) -> np.ndarray:
        d = np.zeros(len(t1950))
        if sta in self.tropo:
            d += self.tropo[sta].delay(t1950, el)
        if self.iono is not None:
            d += self.iono.phase_delay(sta, t1950, freq)
        return d


# ---------------------------------------------------------------------------
# 5. The light-time solution and the predicted received frequency
# ---------------------------------------------------------------------------

def _solve(t3: np.ndarray, st: Station, eph: Ephemeris, media: Media | None,
           f_up: np.ndarray, f_down: np.ndarray, n_iter: int = 4) -> dict:
    """Round-trip light-time solution for reception times t3 (TDB)."""
    S3, V3 = st.bcrs(t3)
    up3 = st.up(t3)
    t2 = t3 - np.linalg.norm(S3 - eph.pos("cassini", t3), axis=1) / C_LIGHT
    u3 = st.utc(t3)
    for _ in range(n_iter):
        C2 = eph.pos("cassini", t2)
        sh1, sh2, b_d = shapiro_delays(C2, S3, t2, t3, eph)
        el3 = _elevation(up3, S3, V3, C2)
        med = media.delay(st.sta, u3, el3, f_down) if media else 0.0
        t2 = t3 - np.linalg.norm(S3 - C2, axis=1) / C_LIGHT - sh1 - sh2 - med
    C2, VC = eph.pos("cassini", t2), eph.vel("cassini", t2)
    D_d = dict(sh=sh1, sh2=sh2, med=med)
    t1 = t2 - np.linalg.norm(C2 - st.bcrs(t2)[0], axis=1) / C_LIGHT
    for _ in range(n_iter):
        S1, V1 = st.bcrs(t1)
        u1 = st.utc(t1)
        su1, su2, b_u = shapiro_delays(S1, C2, t1, t2, eph)
        el1 = _elevation(st.up(t1), S1, V1, C2)
        medu = media.delay(st.sta, u1, el1, f_up) if media else 0.0
        t1 = t2 - np.linalg.norm(C2 - S1, axis=1) / C_LIGHT - su1 - su2 - medu
    S1, V1 = st.bcrs(t1)
    return dict(t1=t1, t2=t2, t3=t3, S1=S1, V1=V1, C2=C2, VC=VC, S3=S3, V3=V3, el1=el1, el3=el3,
                u1=st.utc(t1), b_up=b_u, b_down=b_d,
                D_d=D_d, D_u=dict(sh=su1, sh2=su2, med=medu))


def predict(t1950: np.ndarray, st: Station, eph: Ephemeris, media: Media | None, ul_band: int, dl_band: int,
            rp: pd.DataFrame, h: float = 30.0) -> pd.DataFrame:
    """Predicted received frequency for one station and link at UTC times t."""
    t3 = st.tdb_from_utc(np.asarray(t1950, dtype=float))
    m2 = M2[(ul_band, dl_band)]
    # leg frequencies for the dispersive (ionosphere) delay: nominal values suffice
    f_up = np.full(len(t3), {2: 7.175e9, 3: 34.317e9}[ul_band])
    f_down = m2 * f_up
    L = _solve(t3, st, eph, media, f_up, f_down)
    Lp = _solve(t3 + h, st, eph, media, f_up, f_down)
    Lm = _solve(t3 - h, st, eph, media, f_up, f_down)

    def rate(key_leg, part, tvar):
        dp, dm = Lp[key_leg][part], Lm[key_leg][part]
        if np.isscalar(dp):
            return np.zeros(len(t3))
        return (dp - dm) / (Lp[tvar] - Lm[tvar])

    # delay rates: downlink per unit t3, uplink per unit t2
    rd_sh = rate("D_d", "sh", "t3")
    rd_other = rate("D_d", "sh2", "t3") + rate("D_d", "med", "t3")
    ru_sh = rate("D_u", "sh", "t2")
    ru_other = rate("D_u", "sh2", "t2") + rate("D_u", "med", "t2")
    nu = L["C2"] - L["S1"]
    nu /= np.linalg.norm(nu, axis=1)[:, None]
    nd = L["S3"] - L["C2"]
    nd /= np.linalg.norm(nd, axis=1)[:, None]
    c = C_LIGHT
    a_num = 1 - np.sum(nu * L["VC"], 1) / c
    a_den = 1 - np.sum(nu * L["V1"], 1) / c
    b_num = 1 - np.sum(nd * L["V3"], 1) / c
    b_den = 1 - np.sum(nd * L["VC"], 1) / c
    F_off = (a_num - ru_other) / a_den * (b_num - rd_other) / b_den
    F_on = (a_num - ru_other - ru_sh) / a_den * (b_num - rd_other - rd_sh) / b_den
    clock = (1 - st.clock_rate(L["t1"])) / (1 - st.clock_rate(t3))
    fT, sweep = transmitted_frequency(rp, st.sta, ul_band, L["u1"])
    f_off = m2 * fT * F_off * clock
    out = pd.DataFrame({
        "t_tdb": t3, "t1_tdb": L["t1"], "t1_utc": L["u1"], "f_t": fT, "uplink_sweep": sweep,
        "f_model": f_off, "dop": F_off * clock, "y_gr": 1 - F_on / F_off,
        "el_up_deg": np.rad2deg(L["el1"]), "el_down_deg": np.rad2deg(L["el3"]),
        "b_up_rsun": L["b_up"] / C.R_SUN, "b_down_rsun": L["b_down"] / C.R_SUN,
        "delay_2way_us": (L["D_u"]["sh"] + L["D_d"]["sh"]) * 1e6,
        "y_2pn": rate("D_u", "sh2", "t2") + rate("D_d", "sh2", "t3"),
        "y_media": rate("D_u", "med", "t2") + rate("D_d", "med", "t3"),
        "y_clock": 1 - clock,
    })
    # partials of the fractional residual with respect to Cassini's BCRS
    # position and velocity at t2 (for the trajectory-correction fit)
    VC, V1, V3 = L["VC"], L["V1"], L["V3"]
    Ru = np.linalg.norm(L["C2"] - L["S1"], axis=1)[:, None]
    Rd = np.linalg.norm(L["S3"] - L["C2"], axis=1)[:, None]
    dFdV = (nd - nu) / c
    wu = VC - V1
    wd = VC - V3
    dFdC = (-(wu - np.sum(nu * wu, 1)[:, None] * nu) / Ru - (wd - np.sum(nd * wd, 1)[:, None] * nd) / Rd) / c
    for k, ax in enumerate("xyz"):
        out[f"dy_dv_{ax}"] = -dFdV[:, k]
        out[f"dy_dr_{ax}"] = -dFdC[:, k]
    out["t2_tdb"] = L["t2"]
    return out


# ---------------------------------------------------------------------------
# 6. Residuals for all data
# ---------------------------------------------------------------------------

#: Model configurations. "main" uses Horizons' trajectory over the whole of
#: SCE1 (June 6 - July 5) and is the general-purpose model; "nav" uses the
#: navigation team's SCE1 trajectory (June 9 - July 5), the base of the gamma
#: fit (see cassini_gamma); the others are its robustness variants
#: (troposphere calibration, radiometer unit, ionosphere).
MODEL_CONFIGS = {
    "main": dict(tropo="amc", amc_unit=1, iono=True, cassini="horizons"),
    "nav": dict(tropo="amc", amc_unit=1, iono=True, cassini="nav"),
    "nav_tropo_dsn": dict(tropo="dsn", amc_unit=1, iono=True, cassini="nav"),
    "nav_awvr2": dict(tropo="amc", amc_unit=2, iono=True, cassini="nav"),
    "nav_no_iono": dict(tropo="amc", amc_unit=1, iono=False, cassini="nav"),
    "nav_no_tides": dict(tropo="amc", amc_unit=1, iono=True, cassini="nav", tides=False),
}
#: configurations whose 1-s residuals are cached (~0.4 GB each)
CACHED_RESIDUALS = ("main", "nav")


def resid_cache_path(config: str) -> Path:
    return RESID_CACHE if config == "main" else C.CACHE_DIR / f"cassini_resid_1s_{config}.parquet"


def residuals(config: str = "main", refresh: bool = False, chunk: int = 20000, log=print) -> pd.DataFrame:
    """Model every 1-s observable. Returns the observables with the predicted
    frequency, the residual fractional frequency y_obs (Shapiro off), GR's
    signal y_GR, geometry and the trajectory partials. Configurations in
    CACHED_RESIDUALS are cached at 1 s; the others keep only their normal
    points."""
    cache = resid_cache_path(config) if config in CACHED_RESIDUALS else None
    if cache is not None and cache.exists() and not refresh:
        return pd.read_parquet(cache)
    cfg = MODEL_CONFIGS[config]
    obs = load_observables()
    rp = load_ramps()
    eph = Ephemeris(cfg["cassini"])
    t0, t1 = obs["t1950"].min() - 4 * 3600, obs["t1950"].max() + 600
    stations = {s: Station(s, eph, t0, t1, tides=cfg.get("tides", True))
                for s in sorted(obs["station"].unique()) if s in DSN_XYZ}
    media = Media(stations, tropo=cfg["tropo"], iono=cfg["iono"], amc_unit=cfg["amc_unit"])
    out = []
    for (sta, ul, dl), g in obs.groupby(["station", "ul_band", "dl_band"]):
        if sta not in stations or (ul, dl) not in M2:
            continue
        for i in range(0, len(g), chunk):
            gg = g.iloc[i:i + chunk]
            pr = predict(gg["t1950"].to_numpy(), stations[sta], eph, media, int(ul), int(dl), rp)
            pr.index = gg.index
            out.append(pd.concat([gg, pr], axis=1))
        log(f"  [{config}] modelled station {sta} {link_name(ul, dl)}: {len(g)} points")
    res = pd.concat(out).sort_values(["station", "ul_band", "dl_band", "t1950"]).reset_index(drop=True)
    res["uplink_estimated"] = False
    res = _estimate_missing_ka_uplink(res, log)
    res["y_obs"] = (res["f_model"] - res["f_r"]) / res["f_model"]
    res["utc"] = O.utc_from_1950(res["t1950"].to_numpy())
    if cache is not None:
        res.to_parquet(cache)
    return res


def _estimate_missing_ka_uplink(res: pd.DataFrame, log=print) -> pd.DataFrame:
    """Two Ka/Ka passes (2002 June 14 and 17) exist only in the TDFs, whose
    Ka ramp tables the archive lacks (no k252 ODF/TDF on those days), so the
    Ka uplink frequency is unknown. The DSN reset it every pass and the
    daily values do not follow a smooth curve (a 5th-degree fit leaves 40 Hz,
    1e-9), so each such pass gets the constant f_T that matches its
    simultaneous X/X data in the median; such passes must carry a free bias
    in any fit (flag ``uplink_estimated``)."""
    ka = (res["station"] == 25) & (res["ul_band"] == 3) & res["f_t"].isna()
    if not ka.any():
        return res
    xx = res[(res["station"] == 25) & (res["ul_band"] == 2) & (res["dl_band"] == 2)]
    xx_y = pd.Series(((xx["f_model"] - xx["f_r"]) / xx["f_model"]).to_numpy(), index=xx["t1950"].round(1).to_numpy())
    xx_y = xx_y[~xx_y.index.duplicated()]
    k = res[ka]
    gap = k["t1950"].diff().fillna(1e9) > 1200
    for _, g in k.groupby(gap.cumsum()):
        yx = xx_y.reindex(g["t1950"].round(1).to_numpy()).to_numpy()
        ok = np.isfinite(yx)
        if ok.sum() < 600:
            continue
        ft = g["f_r"].to_numpy()[ok] / (M2[(3, 3)] * g["dop"].to_numpy()[ok] * (1 - yx[ok]))
        f_est = float(np.median(ft))
        idx = g.index
        res.loc[idx, "f_t"] = f_est
        res.loc[idx, "f_model"] = M2[(3, 3)] * f_est * res.loc[idx, "dop"]
        res.loc[idx, "uplink_estimated"] = True
        log(f"  Ka uplink of the pass at {O.utc_from_1950(g['t1950'].iloc[0])} estimated from X/X: {f_est:.1f} Hz")
    return res


# ---------------------------------------------------------------------------
# 7. Passes, cleaning and normal points
# ---------------------------------------------------------------------------

#: minimum elevation (deg) on both legs. The AWVRs retrieve only above 15 deg.
EL_MIN = 15.0
PASS_GAP_S = 1200.0


def clean(res: pd.DataFrame, el_min: float = EL_MIN, clip: float = 5.0, window: int = 61) -> pd.DataFrame:
    """Usable 1-s points with pass labels. Drops points without an uplink
    frequency, those whose uplink was being swept, below ``el_min`` on either
    leg, and outliers beyond ``clip`` robust sigmas from a running median of
    the residual (y_obs - y_GR; the median follows the signal, so the cut is
    signal-blind). Pass = one station's one link without a gap > 20 min;
    track = one station's day of tracking (its links together), the unit of
    the jackknife."""
    r = res[np.isfinite(res["y_obs"]) & ~res["uplink_sweep"] & (res["el_down_deg"] >= el_min)
            & (res["el_up_deg"] >= el_min)].copy()
    r = r.sort_values(["station", "ul_band", "dl_band", "t1950"])
    key = r["station"] * 100 + r["ul_band"] * 10 + r["dl_band"]
    new = (key != key.shift()) | (r["t1950"].diff() > PASS_GAP_S)
    r["pass_id"] = new.cumsum().astype(int)
    r["link"] = [link_name(a, b) for a, b in zip(r["ul_band"], r["dl_band"])]
    rr = r["y_obs"] - r["y_gr"]
    # incoherent or mislabelled stretches (see MAX_PASS_OFFSET) go first
    coherent = rr.groupby(r["pass_id"]).transform("median").abs() < MAX_PASS_OFFSET
    r, rr = r[coherent].copy(), rr[coherent]
    med = rr.groupby(r["pass_id"]).transform(lambda s: s.rolling(window, center=True, min_periods=10).median())
    dev = rr - med
    mad = dev.abs().groupby(r["pass_id"]).transform("median") * 1.4826
    keep = dev.abs() <= clip * mad.clip(lower=1e-15)
    r = r[keep].copy()
    # a track: one station's passes that overlap (or follow within 4 h),
    # labelled by the station and the UTC date of its first point
    r["track"] = ""
    for sta, g in r.groupby("station"):
        spans = g.groupby("pass_id")["t1950"].agg(["min", "max"]).sort_values("min")
        label, cur_end, cur = {}, -1.0, None
        for pid, (t0, t1) in spans.iterrows():
            if cur is None or t0 > cur_end + 4 * 3600:
                cur = f"{sta}_{str(O.utc_from_1950(np.array([t0]))[0].astype('datetime64[h]'))[:13]}"
            cur_end = max(cur_end, t1)
            label[pid] = cur
        r.loc[g.index, "track"] = g["pass_id"].map(label)
    return r


#: a pass whose median residual is this far from zero is not a coherent
#: two-way link of the stated bands (found: a Ka/Ka stretch of June 7 offset
#: by -5.7e-6 with 1e-8 scatter, a non-coherent downlink; one of June 8 offset
#: by -1.57e-4, the X-up/Ka-down link mislabelled Ka/Ka; an X/X stretch of
#: June 17 offset by +3.2e-4). Normal residuals are ~1e-10.
MAX_PASS_OFFSET = 1e-8


def _screen_normal_points(np_: pd.DataFrame, k: float = 5.0, window: int = 31) -> pd.DataFrame:
    """Drop incoherent passes and outlying normal points (beyond k robust
    sigmas from a running median of y_obs - y_GR within the pass)."""
    r = np_["y_obs"] - np_["y_gr"]
    med_pass = r.groupby(np_["pass_id"]).transform("median")
    ok = med_pass.abs() < MAX_PASS_OFFSET
    np_, r = np_[ok].copy(), r[ok]
    np_ = np_.sort_values(["pass_id", "t1950"])
    r = r.loc[np_.index]
    run = r.groupby(np_["pass_id"]).transform(lambda x: x.rolling(window, center=True, min_periods=5).median())
    dev = r - run
    mad = dev.abs().groupby(np_["pass_id"]).transform("median") * 1.4826
    return np_[(dev.abs() <= k * mad.clip(lower=1e-16)) & (dev.abs() < MAX_PASS_OFFSET)]


_NP_MEAN = ["t1950", "t_tdb", "t2_tdb", "y_obs", "y_gr", "y_2pn", "y_media", "f_t", "f_model", "el_up_deg",
            "el_down_deg", "b_up_rsun", "b_down_rsun", "delay_2way_us",
            "dy_dv_x", "dy_dv_y", "dy_dv_z", "dy_dr_x", "dy_dr_y", "dy_dr_z"]


def np_cache_path(config: str = "main", el_min: float = EL_MIN) -> Path:
    tag = "" if (config == "main" and el_min == EL_MIN) else f"_{config}_el{el_min:g}"
    return C.CACHE_DIR / f"cassini_normal_points{tag}.parquet"


def normal_points(config: str = "main", refresh: bool = False, el_min: float = EL_MIN, bin_s: float = 60.0,
                  min_frac: float = 2 / 3, res: pd.DataFrame | None = None, log=print) -> pd.DataFrame:
    """Average the cleaned 1-s points into ``bin_s`` normal points per pass
    (bins on a fixed UTC grid, so different links share bin edges; a bin
    needs ``min_frac`` of its seconds), for every link and for the X/X +
    Ka/Ka plasma-free combination (``cassini_plasma.two_link_1s``). The 1-s
    scatter in the bin gives a white-noise error ``y_err_white``. Cached per
    model configuration and elevation cut; rebuilt from the 1-s residuals
    (and those from the ODFs) when missing."""
    cache = np_cache_path(config, el_min)
    if cache.exists() and not refresh and res is None:
        return pd.read_parquet(cache)
    import cassini_plasma

    res = residuals(config, log=log) if res is None else res
    cl = clean(res, el_min=el_min)
    cl = pd.concat([cl, cassini_plasma.two_link_1s(cl)], ignore_index=True)
    cl["tbin"] = np.floor(cl["t1950"] / bin_s).astype(np.int64)
    g = cl.groupby(["station", "ul_band", "dl_band", "link", "pass_id", "track", "tbin"])
    agg = g[_NP_MEAN].mean()
    agg["n"] = g.size()
    agg["y_sd"] = g["y_obs"].std()
    if "y_plasma_x" in cl:
        agg["y_plasma_x"] = g["y_plasma_x"].mean()
    agg["uplink_estimated"] = g["uplink_estimated"].any()
    agg["source"] = g["source"].first()
    np_ = agg.reset_index()
    np_ = np_[np_["n"] >= min_frac * bin_s].copy()
    np_ = _screen_normal_points(np_)
    np_["y_err_white"] = np_["y_sd"] / np.sqrt(np_["n"])
    np_["utc"] = O.utc_from_1950(np_["t1950"].to_numpy())
    np_ = np_.sort_values(["t1950", "station", "link"]).reset_index(drop=True)
    np_.to_parquet(cache)
    return np_


# ---------------------------------------------------------------------------
# 8. Inventory of the archive's closed-loop Doppler
# ---------------------------------------------------------------------------

def inventory(refresh: bool = False) -> pd.DataFrame:
    """Every closed-loop Doppler data set of SCE1: one row per day, data type
    (1-, 2- or 3-way), receiving and transmitting station and link, with the
    number of 1-s points and their time span. Two-way Ka/Ka points that only
    the TDFs hold are listed with source 'tdf'."""
    if INVENTORY_PATH.exists() and not refresh:
        return pd.read_csv(INVENTORY_PATH)
    rows = []
    for p in _odf_paths():
        d = O.read_odf(p)["doppler"]
        d = d[d["dtype"].isin([11, 12, 13]) & (d["bad"] == 0)]
        for (dt, rcv, xmt, ul, dl), g in d.groupby(["dtype", "rcv", "xmt", "ul_band", "dl_band"]):
            rows.append({"file_day": int(Path(p).name[12:15]), "kind": {11: "1-way", 12: "2-way", 13: "3-way"}[dt],
                         "rx_station": int(rcv), "tx_station": int(xmt) if dt != 11 else -1,
                         "link": (BAND_NAME[dl] + " (1-way)") if dt == 11 else link_name(ul, dl), "source": "odf",
                         "n_1s": len(g), "start_utc": str(g["utc"].min())[:19], "end_utc": str(g["utc"].max())[:19]})
    obs = load_observables()
    t = obs[obs["source"] == "tdf"]
    gap = t["t1950"].diff().fillna(1e9) > PASS_GAP_S
    for _, g in t.groupby(gap.cumsum()):
        u = O.utc_from_1950(g["t1950"].to_numpy())
        rows.append({"file_day": int((u.min().astype("datetime64[D]") - np.datetime64("2002-01-01")).astype(int) + 1),
                     "kind": "2-way", "rx_station": 25, "tx_station": 25, "link": "Ka/Ka", "source": "tdf",
                     "n_1s": len(g), "start_utc": str(u.min())[:19], "end_utc": str(u.max())[:19]})
    inv = pd.DataFrame(rows).sort_values(["file_day", "kind", "rx_station", "link", "start_utc"]).reset_index(drop=True)
    inv["start_utc"] = inv["start_utc"].str.replace("T", " ")
    inv["end_utc"] = inv["end_utc"].str.replace("T", " ")
    inv["hours"] = ((pd.to_datetime(inv["end_utc"]) - pd.to_datetime(inv["start_utc"])).dt.total_seconds() / 3600).round(2)
    inv.to_csv(INVENTORY_PATH, index=False)
    return inv


def dss25_tracks(np_: pd.DataFrame) -> pd.DataFrame:
    """Per DSS 25 track: hours of usable X/X, Ka/Ka and simultaneous (the
    plasma-free combination) 60-s normal points, and the impact parameter."""
    d = np_[np_["station"] == 25]
    g = d.groupby(["track", "link"]).size().unstack(fill_value=0) / 60.0
    out = g.rename(columns=lambda c: f"hours_{c}")
    out["b_min_rsun"] = d.groupby("track")["b_down_rsun"].min()
    out["b_max_rsun"] = d.groupby("track")["b_down_rsun"].max()
    out["first_utc"] = d.groupby("track")["utc"].min().astype(str).str[:16]
    return out.reset_index()


def run_cassini(rebuild: bool = False, variants: bool = True, log=print) -> dict:
    """Rebuild the whole Cassini analysis (for the project's run.py): normal
    points from the cache (or from the ODFs/TDFs if missing, ~3 min), the
    gamma fits with every variant, and the outputs in data/cache
    (cassini_fit.json, cassini_measured_series.parquet). See cassini_gamma.run."""
    import cassini_gamma

    return cassini_gamma.run(rebuild=rebuild, variants=variants, log=log)
