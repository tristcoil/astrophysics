"""Neutron stars are born moving: pulsar speeds, and pulsars traced home.

A core-collapse explosion is never perfectly symmetric. If more ejecta fly
out on one side, momentum conservation sends the neutron star the other way:
a "kick" of hundreds of km/s. We measure it from the ATNF Pulsar Catalogue:

* **Speeds.** Proper motion (mas/yr) times distance gives the speed across
  the sky: v = 4.74 km/s x mu[mas/yr] x d[kpc]. Part of what we see is the
  Galaxy's rotation and the Sun's own motion; following Hobbs et al. (2005)
  we remove both (flat rotation curve, R0 = 8.5 kpc, V0 = 225 km/s, solar
  peculiar motion 15.6 km/s towards l = 48.8, b = 26.3 deg).
* **Traced home.** For pulsars the catalogue associates with a supernova
  remnant, does the proper motion point away from the remnant's centre, and
  does the time to get there match the pulsar's age?
"""

from __future__ import annotations

import re
import warnings

import numpy as np
import pandas as pd

import config
import literature as L
import remnants


def _unit_lb(l_deg, b_deg):
    l, b = np.radians(l_deg), np.radians(b_deg)
    return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])


def galactic_rotation_pm(l_deg: float, b_deg: float, d_kpc: float, *, R0: float = config.HOBBS_R0_KPC,
                         V0: float = config.HOBBS_V0_KMS, sun_pec=config.HOBBS_SUN_PECULIAR) -> tuple[float, float]:
    """Proper motion (mu_l*, mu_b in mas/yr) a star at (l, b, d) would show if
    it simply followed the Galaxy's flat rotation, seen from the moving Sun.

    Heliocentric Galactic axes: x to the Galactic centre, y towards l = 90,
    z to the north Galactic pole. Rotation is clockwise seen from the north,
    so the Sun moves towards +y."""
    pos = d_kpc * _unit_lb(l_deg, b_deg)
    X, Y = pos[0] - R0, pos[1]  # galactocentric, in the plane
    R = np.hypot(X, Y)
    v_star = V0 * np.array([Y / R, -X / R, 0.0])
    v_sun = np.array([0.0, V0, 0.0]) + sun_pec[0] * _unit_lb(sun_pec[1], sun_pec[2])
    dv = v_star - v_sun
    l, b = np.radians(l_deg), np.radians(b_deg)
    e_l = np.array([-np.sin(l), np.cos(l), 0.0])
    e_b = np.array([-np.sin(b) * np.cos(l), -np.sin(b) * np.sin(l), np.cos(b)])
    k = config.KMS_PER_MASYR_KPC * d_kpc
    return float(dv @ e_l / k), float(dv @ e_b / k)


def peculiar_transverse_velocity(ra: float, dec: float, d_kpc: float, pmra: float, pmdec: float) -> dict:
    """Speed across the sky relative to the star's local Galactic rotation."""
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, pm_ra_cosdec=pmra * u.mas / u.yr, pm_dec=pmdec * u.mas / u.yr,
                 frame="icrs")
    g = c.galactic
    ml, mb = galactic_rotation_pm(g.l.deg, g.b.deg, d_kpc)
    pl = g.pm_l_cosb.to_value(u.mas / u.yr) - ml
    pb = g.pm_b.to_value(u.mas / u.yr) - mb
    # Back to equatorial, to quote a position angle.
    back = SkyCoord(l=g.l, b=g.b, pm_l_cosb=pl * u.mas / u.yr, pm_b=pb * u.mas / u.yr, frame="galactic").icrs
    pra, pdec = back.pm_ra_cosdec.to_value(u.mas / u.yr), back.pm_dec.to_value(u.mas / u.yr)
    k = config.KMS_PER_MASYR_KPC * d_kpc
    return {"v_pec_kms": float(np.hypot(pl, pb) * k), "v_obs_kms": float(np.hypot(pmra, pmdec) * k),
            "pmra_pec": float(pra), "pmdec_pec": float(pdec),
            "pa_deg": float(np.degrees(np.arctan2(pra, pdec)) % 360.0)}


# --------------------------------------------------------------------------
# The catalogue
# --------------------------------------------------------------------------

def _raw_assoc() -> dict:
    """PSRJ -> the full ASSOC string (psrqpy keeps only the first entry)."""
    out, name = {}, None
    for line in config.ATNF_DB.read_text(errors="replace").splitlines():
        if line.startswith("PSRJ"):
            name = line.split()[1]
        elif line.startswith("ASSOC") and name:
            out[name] = line.split(None, 1)[1].split()[0] if len(line.split()) > 1 else ""
    return out


def load_pulsars() -> pd.DataFrame:
    """Pulsars with a measured proper motion and a distance, with speeds."""
    import psrqpy

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat = psrqpy.QueryATNF(loadfromdb=str(config.ATNF_DB)).catalogue
    assoc = _raw_assoc()
    num = lambda c: pd.to_numeric(cat[c], errors="coerce")
    t = pd.DataFrame({
        "jname": cat["JNAME"], "ra": num("RAJD"), "dec": num("DECJD"), "gl": num("GL"), "gb": num("GB"),
        "pmra": num("PMRA"), "pmdec": num("PMDEC"), "dist": num("DIST"), "dist_dm": num("DIST_DM"),
        "px": num("PX"), "p0": num("P0"), "p1": num("P1"), "age": num("AGE"),
        "type": cat["TYPE"].astype(str), "binary": cat["BINARY"].astype(str),
    })
    t["assoc"] = t["jname"].map(assoc).fillna("")
    t = t[np.isfinite(t.pmra) & np.isfinite(t.pmdec) & np.isfinite(t.dist) & (t.dist > 0)].copy()
    t["mu"] = np.hypot(t.pmra, t.pmdec)
    t["v_obs"] = config.KMS_PER_MASYR_KPC * t.mu * t.dist
    pec = [peculiar_transverse_velocity(r.ra, r.dec, r.dist, r.pmra, r.pmdec) for r in t.itertuples()]
    t["v_pec"] = [p["v_pec_kms"] for p in pec]
    t["pmra_pec"] = [p["pmra_pec"] for p in pec]
    t["pmdec_pec"] = [p["pmdec_pec"] for p in pec]
    # Hobbs et al.'s classes: recycled = P < 0.1 s and Pdot < 1e-17; pulsars
    # in globular clusters or other galaxies are left out of the statistics.
    t["recycled"] = (t.p0 < 0.1) & (t.p1 < 1e-17)
    t["cluster"] = t.assoc.str.contains(r"\bGC:", regex=True)
    t["extragalactic"] = t.assoc.str.contains("EXGAL")
    t["young"] = (~t.recycled) & (t.age < config.YOUNG_PULSAR_AGE_YR)
    snr = t.assoc.str.extract(r"SNR:([^,\[]+)")[0]
    t["snr"] = snr.fillna("")
    return t.reset_index(drop=True)


def rayleigh_sigma(v: np.ndarray) -> tuple[float, float]:
    """Maximum-likelihood 1-D sigma of a Maxwellian, from 2-D (projected)
    speeds: an isotropic 3-D Maxwellian projects to a Rayleigh distribution
    with the same sigma, whose ML estimate is sqrt(sum v^2 / 2N)."""
    v = np.asarray(v)
    s = np.sqrt(np.sum(v**2) / (2 * len(v)))
    return float(s), float(s / np.sqrt(4 * len(v)))


def speed_statistics(t: pd.DataFrame, seed: int = 2) -> dict:
    rng = np.random.default_rng(seed)
    galactic = t[~t.cluster & ~t.extragalactic]
    out = {}
    for name, sel in (("normal", ~galactic.recycled), ("young", galactic.young),
                      ("recycled", galactic.recycled)):
        s = galactic[sel]
        rec = {"n": int(len(s))}
        for col in ("v_obs", "v_pec"):
            v = s[col].to_numpy()
            boots = [rng.choice(v, len(v)).mean() for _ in range(2000)]
            rec[f"{col}_mean"] = float(v.mean())
            rec[f"{col}_mean_err"] = float(np.std(boots))
            rec[f"{col}_median"] = float(np.median(v))
            rec[f"{col}_sigma"], rec[f"{col}_sigma_err"] = rayleigh_sigma(v)
        out[name] = rec
    return out


#: Catalogue nicknames -> Green's names (ATNF writes either).
SNR_NICKNAMES = {
    "Crab": "G184.6-5.8", "Vela": "G263.9-3.3", "S147": "G180.0-1.7", "CTB80": "G69.0+2.7",
    "3C58": "G130.7+3.1", "CTA1": "G119.5+10.2", "W44": "G34.7-0.4", "PUPPIS_A": "G260.4-3.4",
    "CTB109": "G109.1-1.0", "CTB_37A": "G348.5+0.1", "Kes73": "G27.4+0.0", "Kes75": "G29.7-0.3",
    "W28": "G6.4-0.1", "W41": "G23.3-0.3", "RX_J0852.0-4622": "G266.2-1.2",
}


def trace_home(t: pd.DataFrame, green: pd.DataFrame) -> pd.DataFrame:
    """Pulsars with a remnant association and a proper motion: their offset
    from the remnant's catalogued centre, and whether they move away from it.

    The angle is between the proper motion and the line from the remnant's
    centre to the pulsar: near 0 deg if the pulsar was born at the centre and
    has flown straight out. The travel time is offset / proper motion."""
    rows = []
    for r in t[t.snr != ""].itertuples():
        raw = r.snr.replace("(?)", "")
        doubtful = "(?)" in r.snr
        gname = SNR_NICKNAMES.get(raw)
        if gname is None:
            m = re.match(r"(G\d+\.\d+[+-]\d+\.\d+)", raw)
            gname = m.group(1) if m else None
        if gname is None:
            continue
        g = green[green.name == gname]
        if g.empty:
            continue
        g = g.iloc[0]
        dra = (r.ra - g.ra) * np.cos(np.radians(r.dec)) * 60.0  # arcmin, east
        ddec = (r.dec - g.dec) * 60.0
        offset = float(np.hypot(dra, ddec))
        pa_off = float(np.degrees(np.arctan2(dra, ddec)) % 360.0)
        pa_pm = float(np.degrees(np.arctan2(r.pmra, r.pmdec)) % 360.0)
        angle = abs((pa_pm - pa_off + 180.0) % 360.0 - 180.0)
        travel = offset * 60e3 / r.mu  # yr
        rows.append({"jname": r.jname, "snr": raw, "green": gname, "doubtful": doubtful,
                     "snr_size_arcmin": g.size_arcmin, "offset_arcmin": offset, "offset_over_radius":
                     offset / (0.5 * g.size_arcmin), "pa_offset": pa_off, "pa_pm": pa_pm, "angle": angle,
                     "mu": r.mu, "dist": r.dist, "v_pec": r.v_pec, "travel_yr": travel, "char_age_yr": r.age,
                     "pmra": r.pmra, "pmdec": r.pmdec, "ra": r.ra, "dec": r.dec, "snr_ra": g.ra,
                     "snr_dec": g.dec})
    return pd.DataFrame(rows)


def analyse() -> dict:
    t = load_pulsars()
    green = remnants.load_green()
    homes = trace_home(t, green)
    stats = speed_statistics(t)
    crab = t[t.jname == "J0534+2200"].iloc[0]
    return {
        "n_with_pm": int(len(t)), "stats": stats, "hobbs2005": L.HOBBS2005,
        "homes": homes.to_dict("records"),
        "crab_atnf": {"pmra": float(crab.pmra), "pmdec": float(crab.pmdec), "v_obs": float(crab.v_obs),
                      "v_pec": float(crab.v_pec)},
    }
