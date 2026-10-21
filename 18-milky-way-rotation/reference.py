"""Readers for the comparison data, and the maser rotation speeds.

None of these feed the HI4PI measurement; they are what we grade it against:

* McClure-Griffiths & Dickey (2016): first-quadrant HI terminal velocities from
  a different telescope (the VLA plus Green Bank, 1' resolution).
* Reid et al. (2019): parallaxes and proper motions of masers in young massive
  star-forming regions. From these we compute each star's actual orbital
  speed, a rotation curve measured with distances instead of tangent points.
* Eilers et al. (2019): the circular speed from 23,000 red giant stars (Gaia,
  APOGEE, WISE), 5 to 25 kpc from the centre.
* SPARC (Lelli et al. 2016): rotation curves of 175 other disc galaxies.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

import config as cfg
from config import DATA_DIR

#: Reid et al. (2019) footnote 11: outliers (> 3 sigma in a motion component)
#: removed from their final rotation model, fit A5.
REID19_OUTLIERS = {
    "G015.66-00.49", "G029.86-00.04", "G029.98+00.10", "G030.19-00.16",
    "G030.22-00.18", "G030.41-00.23", "G030.81-00.05", "G032.74-00.07",
    "G078.12+03.63", "G108.59+00.49", "G111.54+00.77",
}

REID19_ARM_NAMES = {
    "GC": "Galactic centre", "Con": "Connecting arm", "3kN": "3-kpc arm (near)",
    "3kF": "3-kpc arm (far)", "Nor": "Norma", "Out": "Outer", "ScF": "Scutum (far)",
    "ScN": "Scutum (near)", "CtN": "Centaurus (near)", "OSC": "Outer Scutum-Centaurus",
    "SgN": "Sagittarius (near)", "SgF": "Sagittarius (far)", "CrN": "Carina (near)",
    "Loc": "Local", "Per": "Perseus", "LoS": "Local spur", "AqS": "Aquarius spur",
    "???": "unassigned",
}


def _require(name: str) -> Path:
    path = DATA_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"{path} not found -- run `python fetch_data.py` first")
    return path


# --------------------------------------------------------------------------
# Published terminal velocities
# --------------------------------------------------------------------------

def mgd16_q1() -> dict:
    """McClure-Griffiths & Dickey (2016) Table 1: l (deg) and the measured
    terminal velocity (km/s). The file's third column, undocumented in the CDS
    ReadMe, is the velocity their joint rotation-curve fit predicts (it equals
    220 (0.889 - 0.829 sin l) to 0.01 km/s); we do not use it."""
    arr = np.loadtxt(_require("mgd16_q1_terminal.dat"))
    return {"l": arr[:, 0], "v": arr[:, 1], "v_model": arr[:, 2]}


# --------------------------------------------------------------------------
# Masers (Reid et al. 2019)
# --------------------------------------------------------------------------

@dataclass
class Masers:
    name: np.ndarray
    arm: np.ndarray
    l: np.ndarray
    b: np.ndarray
    plx: np.ndarray
    e_plx: np.ndarray
    pm_e: np.ndarray
    e_pm_e: np.ndarray
    pm_n: np.ndarray
    e_pm_n: np.ndarray
    vlsr: np.ndarray
    e_vlsr: np.ndarray
    ra: np.ndarray
    dec: np.ndarray


def reid19_masers() -> Masers:
    """Parse the fixed-width CDS table (byte positions from its ReadMe)."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    rows = [line.rstrip("\n") for line in _require("reid19_masers.dat").read_text().splitlines()
            if line.strip()]

    def col(a, b):  # 1-based inclusive byte range
        return [r[a - 1:b] for r in rows]

    def fcol(a, b):
        return np.array([float(x) for x in col(a, b)])

    ra = 15 * (fcol(32, 33) + fcol(35, 36) / 60 + fcol(38, 44) / 3600)
    sign = np.where(np.array(col(47, 47)) == "-", -1.0, 1.0)
    dec = sign * (fcol(48, 49) + fcol(51, 52) / 60 + fcol(54, 59) / 3600)
    gal = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs").galactic
    return Masers(
        name=np.array([s.strip() for s in col(1, 13)]),
        arm=np.array([s.strip() for s in col(114, 116)]),
        l=gal.l.deg, b=gal.b.deg,
        plx=fcol(63, 67), e_plx=fcol(69, 73),
        pm_e=fcol(77, 82), e_pm_e=fcol(84, 87),
        pm_n=fcol(91, 96), e_pm_n=fcol(98, 101),
        vlsr=fcol(105, 108), e_vlsr=fcol(110, 111),
        ra=ra, dec=dec,
    )


def maser_rotation(m: Masers, *, n_mc: int = 2000, seed: int = 1) -> dict:
    """Each maser's Galactocentric radius and orbital speed, with Monte Carlo
    errors from its parallax, proper motion and velocity.

    The full 3D velocity comes from the proper motion (sky-plane) and the
    line-of-sight velocity, converted from LSRK back to heliocentric by
    removing the standard solar motion. The Galactic frame uses Reid et al.'s
    own fit A5 (R0, the Sun's velocity, and its height above the plane,
    5.5 pc).
    """
    import astropy.units as u
    from astropy.coordinates import CartesianDifferential, Galactocentric, SkyCoord

    import galaxy

    std = galaxy.standard_solar_motion()
    lr, br = np.radians(m.l), np.radians(m.b)
    std_los = (std[0] * np.cos(lr) + std[1] * np.sin(lr)) * np.cos(br) + std[2] * np.sin(br)
    frame = Galactocentric(
        galcen_distance=cfg.R0_KPC * u.kpc, z_sun=5.5 * u.pc,
        galcen_v_sun=CartesianDifferential(
            [cfg.U_SUN_KMS, cfg.OMEGA_SUN * cfg.R0_KPC, cfg.W_SUN_KMS] * u.km / u.s))

    rng = np.random.default_rng(seed)
    n = m.l.size
    draws = lambda mu, sig: mu[None, :] + sig[None, :] * rng.standard_normal((n_mc, n))
    plx = draws(m.plx, m.e_plx)
    plx = np.where(plx > 0.01, plx, np.nan)
    pme, pmn = draws(m.pm_e, m.e_pm_e), draws(m.pm_n, m.e_pm_n)
    vhel = draws(m.vlsr, m.e_vlsr.astype(float)) - std_los[None, :]
    # the nominal values first, then the Monte Carlo draws, in one transform
    plx = np.vstack([m.plx[None, :], plx])
    pme = np.vstack([m.pm_e[None, :], pme])
    pmn = np.vstack([m.pm_n[None, :], pmn])
    vhel = np.vstack([(m.vlsr - std_los)[None, :], vhel])
    ra = np.broadcast_to(m.ra, plx.shape)
    dec = np.broadcast_to(m.dec, plx.shape)
    c = SkyCoord(ra=ra.ravel() * u.deg, dec=dec.ravel() * u.deg,
                 distance=(1.0 / plx.ravel()) * u.kpc,
                 pm_ra_cosdec=pme.ravel() * u.mas / u.yr, pm_dec=pmn.ravel() * u.mas / u.yr,
                 radial_velocity=vhel.ravel() * u.km / u.s, frame="icrs")
    gc = c.transform_to(frame)
    x, y = gc.x.to_value(u.kpc), gc.y.to_value(u.kpc)
    vx, vy = gc.v_x.to_value(u.km / u.s), gc.v_y.to_value(u.km / u.s)
    R = np.hypot(x, y).reshape(plx.shape)
    # Rotation is clockwise seen from the north pole, i.e. negative angular
    # momentum about +z in astropy's frame; the orbital speed is -L_z / R.
    vrot = (-(x * vy - y * vx)).reshape(plx.shape) / R
    vrad = ((x * vx + y * vy)).reshape(plx.shape) / R  # outward
    return {
        "R": R[0], "vrot": vrot[0], "vrad": vrad[0],
        "R_err": np.nanstd(R[1:], axis=0), "vrot_err": np.nanstd(vrot[1:], axis=0),
        "frac_plx_err": m.e_plx / m.plx,
    }


# --------------------------------------------------------------------------
# Eilers et al. (2019)
# --------------------------------------------------------------------------

def eilers19() -> dict:
    """R (kpc), v_c, and asymmetric errors (km/s); they assumed R0 = 8.122 kpc."""
    arr = np.loadtxt(_require("eilers19_vc.txt"))
    return {"R": arr[:, 0], "vc": arr[:, 1], "err_lo": arr[:, 2], "err_hi": arr[:, 3]}


# --------------------------------------------------------------------------
# SPARC
# --------------------------------------------------------------------------

def _mrt_rows(path: Path, n_fields: int) -> list[list[str]]:
    """Whitespace-split MRT data rows (SPARC's Table 1 rows sit one byte off
    their declared positions, as project D4 found, so byte slicing fails)."""
    lines = path.read_text().splitlines()
    rules = [i for i, line in enumerate(lines) if line.startswith("-----")]
    rows = [line.split() for line in lines[rules[-1] + 1:] if line.strip()]
    bad = [r for r in rows if len(r) != n_fields]
    if bad:
        raise ValueError(f"{path.name}: {len(bad)} rows lack {n_fields} fields, e.g. {bad[0]}")
    return rows


def sparc() -> dict:
    """{name: {"v_flat", "quality", "L36", "R", "Vobs", "eVobs"}} for 175 galaxies."""
    t1 = _mrt_rows(_require("SPARC_Lelli2016c.mrt"), 19)
    t2 = _mrt_rows(_require("MassModels_Lelli2016c.mrt"), 10)
    out = {r[0]: {"v_flat": float(r[15]), "quality": int(r[17]),
                  "L36": float(r[7]) * 1e9, "R": [], "Vobs": [], "eVobs": []} for r in t1}
    for r in t2:
        g = out[r[0]]
        g["R"].append(float(r[2]))
        g["Vobs"].append(float(r[3]))
        g["eVobs"].append(float(r[4]))
    for g in out.values():
        for k in ("R", "Vobs", "eVobs"):
            g[k] = np.array(g[k])
        if g["v_flat"] == 0:
            g["v_flat"] = np.nan
    return out
