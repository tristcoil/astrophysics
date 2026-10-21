"""Eclipses that switch off and on: SS Lac, HS Hya and V907 Sco.

A third star makes the eclipsing pair's orbit precess about the system's
total angular momentum, so the orbit's tilt to our line of sight changes and
the eclipses come and go. We fit one physical model to a century of Harvard
plates (DASCH) and to TESS:

* The angular-momentum geometry: mutual inclination J of the inner and outer
  orbits, the inclination I of the total angular momentum to the line of
  sight, and the phase h0 of the precession at a reference epoch.
* The third star's mass m3. With the inner masses and the periods it fixes,
  through quadrupole secular theory (Soderhjelm 1975; Zasche et al. 2023 eq.
  2), how fast the orbits precess and how the total angular momentum is shared
  between them; the outer orbit's spectroscopic amplitude ties m3 to the
  outer orbit's tilt.
* The eclipse clock: a modern time of primary minimum and the mean period.
  Old eclipses are predicted from the full 3-D orbit, so the drift of eclipse
  phases caused by the turning orbit is part of the model.
* The plate photometry: a zero point, extra scatter, and an outlier
  population (plates have them), plus per-series zero points fixed from the
  out-of-eclipse median.

TESS enters through per-sector chi-square profiles of the eclipse impact
parameter (which also covers sectors with no eclipses: the profile is flat
beyond grazing) and through the times of minimum.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import config as C
import dasch as D
import eclipse_model as E
import limb_darkening as LD
import tess_lc as T

DEG = np.pi / 180.0
BTJD0 = 2457000.0


# ---------------------------------------------------------------------------
# System specifications (sources in config.py)
# ---------------------------------------------------------------------------

def _msini(K1, K2, e, P):
    """m1 sin^3 i and m2 sin^3 i (Msun) from an SB2 orbit (K in km/s, P in d)."""
    f = 1.0361e-7 * (1 - e * e) ** 1.5 * (K1 + K2) ** 2 * P
    return f * K2, f * K1


@dataclass
class Spec:
    name: str
    P: float                       # mean period guess (d)
    T0: float                      # a primary-minimum epoch guess (JD)
    e: float = 0.0
    omega: float = 90.0            # argument of periastron of the secondary (deg), guess
    m1: float | None = None
    m2: float | None = None
    msini: tuple | None = None     # (m1 sin^3 i, m2 sin^3 i, RV epoch JD)
    R1: float | None = None        # Rsun (None: fitted)
    R2: float | None = None
    teff1: float = 6000.0
    teff2: float = 6000.0
    P_out: float = 100.0
    e_out: float = 0.0
    outer_rv: dict = field(default_factory=dict)   # kind 'K' (km/s) or 'a2proj' (Rsun)
    m3_range: tuple = (0.08, 3.0)
    L3_tess_intrinsic: float = 0.0
    tess_kind: str = "spoc"
    t_ref: float = 2459000.0       # reference epoch of the kinematic parameters (JD)
    fit_radii: bool = False
    fit_apsidal: bool = False


SPECS = {
    "SS Lac": Spec(
        "SS Lac", P=14.4161471, T0=2450381.0742, e=0.136, omega=178.28,
        msini=(*_msini(74.69, 76.77, 0.136, 14.4161471), 2451000.0),
        teff1=9500.0, teff2=9500.0, P_out=678.8, e_out=0.159,
        outer_rv=dict(kind="K", value=5.41, err=0.3, epoch=2451000.0, omega=256.7),
        L3_tess_intrinsic=0.01, tess_kind="cut", t_ref=2422000.0, fit_radii=True, fit_apsidal=True),
    "HS Hya": Spec(
        "HS Hya", P=1.568024, T0=2441374.5934, m1=1.31, m2=1.27, R1=1.2747, R2=1.2161,
        teff1=6500.0, teff2=6400.0, P_out=190.53, e_out=0.246,
        outer_rv=dict(kind="K", value=9.02, err=0.5, epoch=2448047.0, omega=111.2),
        L3_tess_intrinsic=0.005, tess_kind="spoc", t_ref=2458550.0),
    "V907 Sco": Spec(
        "V907 Sco", P=3.7762049, T0=2459362.3272, e=0.0, m1=2.74, m2=2.56, R1=2.68, R2=2.53,
        teff1=10000.0, teff2=10004.0, P_out=142.01, e_out=0.467,
        outer_rv=dict(kind="a2proj", value=34.14, err=2.04, epoch=2459760.0, omega=201.4),
        L3_tess_intrinsic=0.015, tess_kind="cut", t_ref=2459383.9),
}


# ---------------------------------------------------------------------------
# Secular geometry
# ---------------------------------------------------------------------------

def secular(m1, m2, m3, P_in, P_out, e_in, e_out, J):
    """Angular-momentum split and nodal precession rate (quadrupole order).

    Returns eps_in, eps_out (radians: tilt of each orbit to the invariable
    plane), nu (rad/day, the precession frequency; negative = regression),
    gamma = L_in/L_out, a_in, a_out (au)."""
    M_in = m1 + m2
    M = M_in + m3
    a_in = C.semimajor_axis_au(M_in, P_in)
    a_out = C.semimajor_axis_au(M, P_out)
    G = C.G_AU_DAY
    L_in = m1 * m2 / M_in * np.sqrt(G * M_in * a_in * (1 - e_in**2))
    L_out = M_in * m3 / M * np.sqrt(G * M * a_out * (1 - e_out**2))
    gam = L_in / L_out
    eps_in = np.arctan2(np.sin(J), gam + np.cos(J))
    eps_out = J - eps_in
    nu = (-0.75 * (m3 / M) * (2 * np.pi * P_in / P_out**2) * (1 - e_out**2) ** -1.5
          * np.cos(J) * np.sqrt(1 + gam**2 + 2 * gam * np.cos(J)))
    return eps_in, eps_out, nu, gam, a_in, a_out


def outer_inclination(I, eps_out, h):
    """Tilt of the outer orbit to the sky: its normal sits opposite the inner
    one on the precession cone (node h + pi)."""
    nz = np.cos(I) * np.cos(eps_out) + np.sin(I) * np.sin(eps_out) * np.cos(h)
    return np.arccos(np.clip(nz, -1, 1))


def _conj_lambda(h, eps, I, g, e):
    """Mean longitude at which the secondary is in front of the primary
    (primary eclipse) for the given orientation."""
    o = np.array([0.0, np.sin(I), np.cos(I)])          # observer in the invariable frame
    ch, sh = np.cos(h), np.sin(h)
    w1 = np.array([ch * o[0] + sh * o[1], -sh * o[0] + ch * o[1], o[2]])   # R_z(h)^T o
    ce, se = np.cos(eps), np.sin(eps)
    w = np.array([w1[0], ce * w1[1] + se * w1[2], -se * w1[1] + ce * w1[2]])  # R_x(eps)^T
    u = np.arctan2(w[1], w[0])                       # argument of latitude at conjunction
    nu_true = u - g
    if e > 0:
        Ecc = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(nu_true / 2))
        M = Ecc - e * np.sin(Ecc)
    else:
        M = nu_true
    return h + g + M


# ---------------------------------------------------------------------------
# The model
# ---------------------------------------------------------------------------

class Model:
    """Parameter vector -> predictions for one system."""

    def __init__(self, name: str):
        self.name = name
        self.spec = SPECS[name]
        sp = self.spec
        # names of the free parameters, in order
        self.names = ["J", "I", "h0", "m3", "T_ref", "P", "m0", "sig_jit", "f_out", "mu_out", "sig_out"]
        if sp.fit_radii:
            self.names += ["r_sum", "k"]
        if sp.fit_apsidal:
            self.names += ["g0", "gdot"]
        self.ndim = len(self.names)
        self.idx = {n: k for k, n in enumerate(self.names)}
        self.ld_B = [LD.quadratic(sp.teff1, 4.0, "B"), LD.quadratic(sp.teff2, 4.0, "B")]
        self.ld_T = [LD.quadratic(sp.teff1, 4.0, "TESS"), LD.quadratic(sp.teff2, 4.0, "TESS")]
        self.sb_B = E.planck_ratio(sp.teff1, sp.teff2, E.LAMBDA_NM["B"])
        self.sb_T = E.planck_ratio(sp.teff1, sp.teff2, E.LAMBDA_NM["TESS"])
        self.ld_V = [LD.quadratic(sp.teff1, 4.0, "V"), LD.quadratic(sp.teff2, 4.0, "V")]
        self.sb_V = E.planck_ratio(sp.teff1, sp.teff2, E.LAMBDA_NM["V"])

    # -- parameter handling ------------------------------------------------
    def unpack(self, theta):
        return {n: theta[k] for k, n in enumerate(self.names)}

    def masses(self, p, i_rv=None):
        sp = self.spec
        if sp.msini is not None:
            s3 = np.sin(i_rv) ** 3 if i_rv is not None else 1.0
            return sp.msini[0] / s3, sp.msini[1] / s3
        return sp.m1, sp.m2

    def geometry(self, p):
        """Everything derived from the parameters (dict)."""
        sp = self.spec
        J, I, h0, m3 = p["J"] * DEG, p["I"] * DEG, p["h0"] * DEG, p["m3"]
        # inner masses: for SS Lac they follow from m sin^3 i at the RV epoch,
        # which itself depends on the geometry: iterate twice
        m1, m2 = self.masses(p)
        for _ in range(3):
            eps_in, eps_out, nu, gam, a_in, a_out = secular(m1, m2, m3, p["P"], sp.P_out, sp.e, sp.e_out, J)
            if sp.msini is None:
                break
            h_rv = h0 + nu * (sp.msini[2] - sp.t_ref)
            i_rv = E.inclination_to_sky(0.0, 0.0, h_rv, 0.0, eps_in, I)
            m1, m2 = self.masses(p, i_rv)
        if sp.fit_radii:
            r1 = p["r_sum"] / (1 + p["k"])
            r2 = p["r_sum"] - r1
        else:
            a_rsun = a_in / C.R_SUN_AU
            r1, r2 = sp.R1 / a_rsun, sp.R2 / a_rsun
        g0 = p["g0"] * DEG if sp.fit_apsidal else np.radians(sp.omega)
        gdot = p["gdot"] * DEG / 365.25 if sp.fit_apsidal else 0.0
        lam_ref = _conj_lambda(h0 + nu * (p["T_ref"] - sp.t_ref), eps_in, I,
                               g0 + gdot * (p["T_ref"] - sp.t_ref), sp.e)
        n_lam = 2 * np.pi / p["P"]
        lam0 = lam_ref - n_lam * (p["T_ref"] - sp.t_ref)
        return dict(m1=m1, m2=m2, m3=m3, J=J, I=I, h0=h0, eps_in=eps_in, eps_out=eps_out, nu=nu,
                    gamma=gam, a_in=a_in, a_out=a_out, r1=r1, r2=r2, g0=g0, gdot=gdot, lam0=lam0,
                    n_lam=n_lam, P_node_yr=2 * np.pi / abs(nu) / 365.25)

    def positions(self, g, t):
        sp = self.spec
        return E.sky_positions(np.asarray(t, float), sp.t_ref, g["lam0"], g["n_lam"], g["h0"], g["nu"],
                               g["g0"], g["gdot"], sp.e, g["eps_in"], g["I"])

    def flux(self, g, t, band="B", L3=0.0):
        sp = self.spec
        xyz = self.positions(g, t)
        r1, r2 = g["r1"], g["r2"]
        (a1, b1), (a2, b2) = {"B": self.ld_B, "V": self.ld_V}.get(band, self.ld_T)
        sb = {"B": self.sb_B, "V": self.sb_V}.get(band, self.sb_T)
        L1 = r1**2 * (1 - a1 / 3 - b1 / 6)
        L2 = r2**2 * (1 - a2 / 3 - b2 / 6) * sb
        L3abs = L3 * (L1 + L2) / max(1e-12, 1 - L3)
        return E.binary_flux(xyz, r1, r2, L1, L2, L3abs, a1, b1, a2, b2)

    def inclination(self, g, t):
        return E.inclination_to_sky(np.asarray(t, float), self.spec.t_ref, g["h0"], g["nu"], g["eps_in"], g["I"])

    def impact(self, g, t):
        """Projected separation at conjunction in units of r1 + r2 (primary,
        secondary eclipse): > 1 means no eclipse."""
        sp = self.spec
        i = self.inclination(g, t)
        gg = g["g0"] + g["gdot"] * (np.asarray(t, float) - sp.t_ref)
        e = sp.e
        d_pri = (1 - e * e) / (1 + e * np.sin(gg)) if e > 0 else 1.0
        d_sec = (1 - e * e) / (1 - e * np.sin(gg)) if e > 0 else 1.0
        rs = g["r1"] + g["r2"]
        return np.abs(np.cos(i)) * d_pri / rs, np.abs(np.cos(i)) * d_sec / rs

    def primary_minimum_near(self, g, t_guess):
        """Time of the primary minimum (closest sky approach with the
        secondary in front) nearest t_guess."""
        P = 2 * np.pi / g["n_lam"]
        tt = t_guess + np.linspace(-0.6, 0.6, 1201) * P
        xyz = self.positions(g, tt)
        rho = np.hypot(xyz[:, 0], xyz[:, 1])
        rho[xyz[:, 2] < 0] = np.inf
        k = int(np.argmin(rho))
        t1 = tt[max(k - 2, 0)]
        t2 = tt[min(k + 2, len(tt) - 1)]

        def f(x):
            q = self.positions(g, np.array([x]))[0]
            return np.hypot(q[0], q[1]) if q[2] > 0 else 9.0
        r = minimize_scalar(f, bounds=(t1, t2), method="bounded", options={"xatol": 1e-6})
        return float(r.x)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def dasch_data(name: str) -> pd.DataFrame:
    """Cleaned DASCH points with per-series zero points removed (series with at
    least 30 points; the rest keep the global zero point)."""
    d, _ = D.load(name)
    med = np.median(d["magcal_magdep"])
    d["mag"] = d["magcal_magdep"]
    for s, g in d.groupby("series"):
        if len(g) >= 30:
            d.loc[g.index, "mag"] = g["magcal_magdep"] - (np.median(g["magcal_magdep"]) - med)
    d["err"] = d["err"].fillna(np.nanmedian(d["err"])).clip(0.03, 1.0)
    return d


def tess_lightcurves(name: str) -> list[dict]:
    """Normalised TESS light curves per sector with their third light (our
    Gaia-based dilution for cutouts, SPOC's crowding for 2-min data, plus the
    third star's own light)."""
    sp = SPECS[name]
    s = C.SYSTEMS[name]
    out = []
    if sp.tess_kind == "spoc":
        for sec, df, hdr in T.spoc_lightcurves(name):
            f = df["flux"] / np.nanmedian(df["flux"])
            out.append(dict(sector=sec, jd=df["btjd"].to_numpy() + BTJD0, flux=f.to_numpy(),
                            err=(df["err"] / np.nanmedian(df["flux"])).to_numpy(),
                            L3=sp.L3_tess_intrinsic, crowd=hdr.get("CROWDSAP"), exptime=2.0 / 1440))
    else:
        for cut in T.cutouts(name):
            df, meta = T.aperture_photometry(cut, s["ra"], s["dec"], r_ap=2.0)
            dil = T.dilution(name, cut, meta)
            f = df["flux"] / np.nanmedian(df["flux"])
            share = dil["target_share"]
            L3 = 1 - share * (1 - sp.L3_tess_intrinsic)
            out.append(dict(sector=cut.sector, jd=df["btjd"].to_numpy() + BTJD0, flux=f.to_numpy(),
                            err=(df["err"] / np.nanmedian(df["flux"])).to_numpy(), L3=L3,
                            crowd=share, exptime=cut.exptime_days))
    return out


#: third light in V (ground photometry): the third star, and for V907 Sco the
#: 9" companion inside ASAS's 15"-pixel apertures
L3_V = {"SS Lac": 0.01, "HS Hya": 0.01, "V907 Sco": 0.05}


def ground_lightcurves(name: str) -> list[dict]:
    """Hipparcos epoch photometry (Hp, 1989-1993) and ASAS-3 (V, 2000-2009,
    one 'sector' per observing season) in the same form as the TESS sectors."""
    import fetch_data as FD
    out = []
    hip = FD.HIP_IDS.get(name)
    if hip:
        rows = []
        for line in (C.RAW / "hipparcos" / f"hip{hip:06d}_epoch.txt").read_text().splitlines()[1:]:
            parts = [x.strip() for x in line.split("|")]
            try:
                t, m, e, flag = float(parts[0]), float(parts[1]), float(parts[2]), int(parts[3] or 0)
            except (ValueError, IndexError):
                continue
            if flag == 0 and e < 0.05:
                rows.append((t + 2440000.0, m, e))
        a = np.array(rows)
        if len(a) > 20:
            f = 10 ** (-0.4 * (a[:, 1] - np.median(a[:, 1])))
            out.append(dict(sector=-1, label="Hipparcos", jd=a[:, 0], flux=f, err=f * a[:, 2] * 0.921,
                            L3=L3_V[name], crowd=None, exptime=0.0, band="V"))
    path = C.RAW / "asas" / f"{name.lower().replace(' ', '_').replace('-', '_')}_asas3.txt"
    if path.exists():
        blocks, cur = [], None
        for line in path.read_text().splitlines():
            if line.startswith("#dataset"):
                cur = []
                blocks.append(cur)
            elif line.strip() and not line.startswith("#") and cur is not None:
                cur.append(line.split())
        pts = []
        for b in blocks:
            arr = [r for r in b if len(r) >= 13 and r[11] in ("A", "B")]
            if len(arr) < 20:
                continue
            t = np.array([float(r[0]) + 2450000.0 for r in arr])
            mags = np.array([[float(x) for x in r[1:6]] for r in arr])
            errs = np.array([[float(x) for x in r[6:11]] for r in arr])
            good = (mags < 29) & (errs < 0.1)
            # aperture with the smallest robust scatter
            sc = [1.4826 * np.median(np.abs(mags[good[:, k], k] - np.median(mags[good[:, k], k])))
                  if good[:, k].sum() > 10 else np.inf for k in range(5)]
            k = int(np.argmin(sc))
            g = good[:, k]
            m0 = np.median(mags[g, k])
            pts.append(np.column_stack([t[g], mags[g, k] - m0, errs[g, k]]))
        if pts:
            a = np.vstack(pts)
            a = a[np.argsort(a[:, 0])]
            yrs = C.jd_to_year(a[:, 0])
            for season in np.unique(np.floor(yrs - 0.5)):     # seasons centred on each year's visibility
                sel = np.floor(yrs - 0.5) == season
                if sel.sum() < 40:
                    continue
                f = 10 ** (-0.4 * a[sel, 1])
                out.append(dict(sector=-int(season), label=f"ASAS-3 {int(season) + 1}", jd=a[sel, 0], flux=f,
                                err=f * a[sel, 2] * 0.921, L3=L3_V[name], crowd=None,
                                exptime=3.0 / 1440, band="V"))
    return out


# ---------------------------------------------------------------------------
# TESS: chi-square profiles of the impact parameter, sector by sector
# ---------------------------------------------------------------------------

#: how far (fraction of the period) the TESS eclipse may sit from the model's
#: guess: wide for SS Lac, whose last timed eclipses are a century old
DT_RANGE = {"SS Lac": 0.12, "HS Hya": 0.03, "V907 Sco": 0.03}


def _sector_windows(model: Model, g, lc, half_width):
    """Indices of points around predicted primary and secondary minima."""
    jd = lc["jd"]
    P = 2 * np.pi / g["n_lam"]
    t_mid = np.median(jd)
    t_pri = model.primary_minimum_near(g, t_mid)
    xyz = model.positions(g, t_pri + np.linspace(0.2, 0.8, 601) * P)
    rho = np.hypot(xyz[:, 0], xyz[:, 1])
    rho[xyz[:, 2] > 0] = np.inf
    dt_sec = (0.2 + 0.6 * np.argmin(rho) / 600) * P
    wins = []
    n0 = int(np.floor((jd.min() - t_pri) / P)) - 1
    n1 = int(np.ceil((jd.max() - t_pri) / P)) + 1
    for n in range(n0, n1 + 1):
        for kind, tc in (("pri", t_pri + n * P), ("sec", t_pri + n * P + dt_sec)):
            sel = np.abs(jd - tc) < half_width
            if sel.sum() >= 10 and (jd[sel].min() < tc - 0.3 * half_width) and (jd[sel].max() > tc + 0.3 * half_width):
                wins.append((kind, tc, np.where(sel)[0]))
    return wins, t_pri, dt_sec


class SectorFit:
    """One TESS sector against the system's eclipse model seen at a fixed tilt.

    flux = eclipse(t) x [slow trend + orbital-phase harmonics]: the trend is a
    cubic per spacecraft orbit (each sector has two, split at the data gap),
    the harmonics (k = 1..4) carry the ellipsoidal and reflection light. Both
    enter linearly and are solved exactly at every trial. Primary and
    secondary eclipses have their own time offsets (small eccentricities and
    light-travel time move them apart).
    """

    def __init__(self, name, lc, theta_guess, nsub=5, nharm=4, baseline="cubic"):
        self.model = Model(name)
        sp = self.model.spec
        p = self.model.unpack(theta_guess)
        self.g0 = self.model.geometry(p)
        self.P = 2 * np.pi / self.g0["n_lam"]
        self.rs = self.g0["r1"] + self.g0["r2"]
        self.lc = lc
        t = lc["jd"]
        self.t, self.f = t, lc["flux"]
        sig = 1.4826 * np.median(np.abs(np.diff(self.f))) / np.sqrt(2)
        self.sigma_point = float(sig)
        self.e = np.maximum(lc["err"], sig)
        nsub = 3 if lc["exptime"] < 0.005 else nsub
        nsub = 1 if lc["exptime"] == 0 else nsub
        self.band = lc.get("band", "TESS")
        self.sub = (np.arange(nsub) - (nsub - 1) / 2) / max(nsub, 1) * lc["exptime"]
        self.t_pri = self.model.primary_minimum_near(self.g0, np.median(t))
        gg = self.g0["g0"] + self.g0["gdot"] * (self.t_pri - sp.t_ref)
        self.d_pri = (1 - sp.e**2) / (1 + sp.e * np.sin(gg)) if sp.e > 0 else 1.0
        self.dt_range = DT_RANGE[name] * self.P
        # linear baseline basis
        gap = np.argmax(np.diff(t))
        seg = (t > t[gap]).astype(int) if np.diff(t)[gap] > 0.5 else np.zeros(len(t), int)
        if baseline == "constant":
            seg = np.zeros(len(t), int)
        cols = []
        for sgi in np.unique(seg):
            sel = seg == sgi
            tm = 0.5 * (t[sel].min() + t[sel].max())
            hw = max(0.5 * (t[sel].max() - t[sel].min()), 1e-3)
            for q in range(4 if baseline == "cubic" else 1):
                c = np.zeros(len(t))
                c[sel] = ((t[sel] - tm) / hw) ** q
                cols.append(c)
        ph = 2 * np.pi * (t - self.t_pri) / self.P
        for k in range(1, nharm + 1):
            cols += [np.cos(k * ph), np.sin(k * ph)]
        self.basis = np.vstack(cols).T
        self.nbasis = self.basis.shape[1]
        # only points within a few eclipse durations need the occultation code
        self.near = np.ones(len(t), bool)

    def _geom(self, x):
        sp = self.model.spec
        i = np.arccos(min(1.0, x * self.rs / self.d_pri))
        gx = dict(self.g0)
        gx.update(eps_in=0.0, I=i, h0=0.0, nu=0.0)
        gx["lam0"] = (_conj_lambda(0.0, 0.0, i, gx["g0"] + gx["gdot"] * (self.t_pri - sp.t_ref), sp.e)
                      - gx["n_lam"] * (self.t_pri - sp.t_ref))
        return gx

    def eclipse_flux(self, x, dt_pri, dt_sec):
        """Eclipse-only model flux at the data times; secondary eclipses are
        shifted by dt_sec instead of dt_pri. Only points near a conjunction
        go through the occultation code (and the exposure sub-sampling)."""
        gx = self._geom(x)
        xyz0 = self.model.positions(gx, self.t - dt_pri)
        sec_side = xyz0[:, 2] < 0
        dt = np.where(sec_side, dt_sec, dt_pri)
        xyz = self.model.positions(gx, self.t - dt)
        rho = np.hypot(xyz[:, 0], xyz[:, 1])
        pad = 2 * np.pi * (abs(self.sub).max() + 1e-9) / self.P + 0.02
        near = rho < self.rs + pad
        out = np.ones(len(self.t))
        if near.any():
            tn = self.t[near] - dt[near]
            tt = (tn[:, None] + self.sub[None, :]).ravel()
            out[near] = self.model.flux(gx, tt, band=self.band, L3=self.lc["L3"]).reshape(len(tn), -1).mean(axis=1)
        return out

    def chi2(self, x, dt_pri, dt_sec=None, return_model=False):
        m = self.eclipse_flux(x, dt_pri, dt_pri if dt_sec is None else dt_sec)
        w = 1 / self.e
        X = self.basis * m[:, None]
        coef, *_ = np.linalg.lstsq(X * w[:, None], self.f * w, rcond=None)
        fit = X @ coef
        r = (self.f - fit) * w
        if return_model:
            return float(r @ r), fit, m
        return float(r @ r)

    def profile(self):
        from scipy.optimize import minimize
        # 1) where (in time) is the eclipse signal? common offset, a few tilts
        dts = np.linspace(-self.dt_range, self.dt_range, 81)
        best = (np.inf, 0.0, 1.2)
        for x in (0.3, 0.7, 0.9, 0.97, 0.995):
            c = [self.chi2(x, d) for d in dts]
            k = int(np.argmin(c))
            if c[k] < best[0]:
                best = (c[k], dts[k], x)
        chi_none = self.chi2(1.2, 0.0)
        n = len(self.t)
        scale0 = max(1.0, chi_none / (n - self.nbasis))
        detected = bool(chi_none - best[0] > 30 * scale0)
        if detected:
            # a real eclipse is in every part of the sector: require > 3 sigma
            # in each spacecraft orbit separately (guards against red noise
            # found by the wide timing search)
            _, fit_b, _ = self.chi2(best[2], best[1], return_model=True)
            _, fit_n, _ = self.chi2(1.2, 0.0, return_model=True)
            w2 = 1 / self.e**2
            dchi = ((self.f - fit_n) ** 2 - (self.f - fit_b) ** 2) * w2
            gap = np.argmax(np.diff(self.t))
            halves = [slice(0, gap + 1), slice(gap + 1, None)] if np.diff(self.t)[gap] > 0.5 else \
                [slice(0, n // 2), slice(n // 2, None)]
            split = [float(dchi[h].sum()) for h in halves]
            self.split = split
            detected = all(d > 9 * scale0 for d in split)
        else:
            self.split = None
        if detected:
            # 2) refine tilt and both offsets
            fun = lambda q: self.chi2(q[0], q[1], q[2])
            r = minimize(fun, [best[2], best[1], best[1]], method="Nelder-Mead",
                         options=dict(xatol=1e-6, fatol=1e-3, maxiter=600))
            xb, dtp, dts_ = r.x
        else:
            xb, dtp, dts_ = best[2], best[1], best[1]
        # 3) profile in x: timing fixed at the best fit when detected, else the
        # most favourable common offset (conservative limits)
        xs = np.unique(np.concatenate([np.linspace(0.0, 0.9, 19), np.linspace(0.9, 1.12, 45),
                                       np.clip(xb + np.linspace(-0.03, 0.03, 31), 0, 1.12)]))
        CH = []
        for x in xs:
            if detected:
                CH.append(self.chi2(x, dtp, dts_))
            else:
                CH.append(min(self.chi2(x, d) for d in dts[::4]))
        CH = np.array(CH)
        kb = int(np.argmin(CH))
        scale = max(1.0, CH[kb] / (n - self.nbasis - 3))
        return dict(sector=self.lc["sector"], jd_mid=float(np.median(self.t)), x=xs.tolist(),
                    chi2=CH.tolist(), x_best=float(xs[kb]), chi2_min=float(CH[kb]), chi2_noecl=float(chi_none),
                    scale=float(scale), detected=detected, n_points=int(n),
                    t_pri_model=float(self.t_pri), dt_pri=float(dtp), dt_sec=float(dts_),
                    t_pri_obs=float(self.t_pri + dtp) if detected else None, L3=float(self.lc["L3"]),
                    split_dchi2=self.split,
                    sigma_point=self.sigma_point, P=float(self.P))
