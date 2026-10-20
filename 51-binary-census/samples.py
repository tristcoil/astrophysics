"""The samples: Gaia DR3's 100 pc census, its non-single stars and wide pairs.

The 100 pc sample is every Gaia DR3 source with parallax > 10 mas. We do *not*
apply the usual RUWE < 1.4 cut of single-star work (B1 used one): a high RUWE is
itself a sign of an unresolved companion, and cutting on it would delete the
very binaries we are counting.

Each star is labelled with what Gaia says about companions:

* ``nss_type``: the solution type in DR3's two-body orbit table (SB1, SB2,
  Orbital, AstroSpectroSB1, EclipsingBinary...), with its period;
* ``accel_type``: an astrometric acceleration solution (curvature of the
  proper motion: an orbit longer than the 34-month baseline);
* ``trend_type``: a radial-velocity trend (same, in velocity);
* ``wide``: a resolved companion in El-Badry et al.'s (2021) catalogue with a
  chance-alignment ratio R < 0.1.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy.table import Table

import config as C
import stellar

# Solution types grouped by what measured them.
SPECTRO_TYPES = ("SB1", "SB1C", "SB2", "SB2C", "AstroSpectroSB1", "EclipsingSpectro")
ASTRO_TYPES = ("Orbital", "OrbitalAlternative", "OrbitalAlternativeValidated",
               "OrbitalTargetedSearch", "OrbitalTargetedSearchValidated")
ECLIPSING_TYPES = ("EclipsingBinary",)

#: Raghavan et al.'s spectral range, F6V-K3V, as BP-RP edges halfway between
#: neighbouring types in Mamajek's sequence (F5V 0.587 | F6V 0.640 ... K3V
#: 1.21 | K4V 1.34).
SOLAR_BPRP = (0.6135, 1.275)

#: Main-sequence band about Mamajek's dwarf sequence, in magnitudes of M_G:
#: an equal-mass binary sits 0.75 mag above it and a triple of equal stars 1.19
#: mag; metal-poor dwarfs sit up to ~0.5 mag below.
MS_BAND = (-1.25, 0.75)


def _f(col) -> np.ndarray:
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def _i(col, fill=-1) -> np.ndarray:
    return np.ma.filled(np.ma.asarray(col), fill)


def _str(col) -> np.ndarray:
    a = np.ma.filled(np.ma.asarray(col).astype(str), "")
    return np.char.strip(a)


@lru_cache(maxsize=None)
def gaia_table(name: str) -> Table:
    return Table.read(C.GAIA_DIR / f"{name}.fits")


def add_photometry(t: Table) -> Table:
    """Distances, absolute magnitudes and height above the dwarf sequence."""
    plx = _f(t["parallax"])
    t["dist_pc"] = 1000.0 / plx
    t["abs_g"] = _f(t["phot_g_mean_mag"]) + 5 * np.log10(plx / 100.0)
    t["bp_rp"] = _f(t["bp_rp"])
    t["dg_ms"] = t["abs_g"] - stellar.ms_ridge_mg(t["bp_rp"])
    t["mass_colour"] = stellar.mass_from_bprp(t["bp_rp"])
    return t


@lru_cache(maxsize=1)
def load_100pc() -> Table:
    """Every DR3 source with parallax > 10 mas, labelled with its NSS solutions."""
    t = gaia_table("gaia_100pc").copy()
    add_photometry(t)
    poe = _f(t["parallax_over_error"])
    t["good_astrometry"] = poe > 10
    t["ms_band"] = (t["dg_ms"] > MS_BAND[0]) & (t["dg_ms"] < MS_BAND[1])
    attach_nss(t)
    attach_wide(t)
    return t


def attach_nss(t: Table) -> None:
    """Copy each star's DR3 non-single-star solutions onto the sample."""
    sid = np.asarray(t["source_id"], dtype=np.int64)
    orbits = gaia_table("nss_orbits")
    rows = preferred_orbit_rows(orbits, sid)
    has = rows >= 0
    t["nss_type"] = np.where(has, _str(orbits["nss_solution_type"])[np.where(has, rows, 0)], "")
    for col in ("period", "eccentricity", "eccentricity_error", "semi_amplitude_primary",
                "semi_amplitude_primary_error", "semi_amplitude_secondary", "mass_ratio",
                "significance", "goodness_of_fit", "a_thiele_innes", "b_thiele_innes",
                "f_thiele_innes", "g_thiele_innes", "nss_parallax", "center_of_mass_velocity"):
        vals = _f(orbits[col])
        t[f"nss_{col}" if not col.startswith("nss_") else col] = np.where(has, vals[np.where(has, rows, 0)],
                                                                          np.nan)
    for name, label in (("nss_accel", "accel_type"), ("nss_trend", "trend_type"), ("nss_vim", "vim_type")):
        other = gaia_table(name)
        o_sid = np.asarray(other["source_id"], dtype=np.int64)
        pos = {s: i for i, s in enumerate(o_sid)}
        rows = np.array([pos.get(s, -1) for s in sid])
        types = _str(other["nss_solution_type"])
        t[label] = np.where(rows >= 0, types[np.where(rows >= 0, rows, 0)], "")
    t["channel"] = np.select(
        [np.isin(t["nss_type"], SPECTRO_TYPES), np.isin(t["nss_type"], ASTRO_TYPES),
         np.isin(t["nss_type"], ECLIPSING_TYPES), (t["accel_type"] != "") | (t["trend_type"] != "")],
        ["spectroscopic", "astrometric", "eclipsing", "acceleration"], default="")


#: When DR3 gives a source two orbits (5,290 sources have separate Orbital and
#: SB1 solutions), we keep one, in this order: combined astrometric+
#: spectroscopic, double-lined, single-lined, astrometric, eclipsing.
SOLUTION_PREFERENCE = ("AstroSpectroSB1", "EclipsingSpectro", "SB2", "SB2C", "SB1", "SB1C", "Orbital",
                       "OrbitalTargetedSearchValidated", "OrbitalAlternativeValidated", "OrbitalTargetedSearch",
                       "OrbitalAlternative", "EclipsingBinary")


def preferred_orbit_rows(orbits: Table, source_ids) -> np.ndarray:
    """Row of each source's preferred two-body solution (-1 if none)."""
    o_sid = np.asarray(orbits["source_id"], dtype=np.int64)
    rank = {t: i for i, t in enumerate(SOLUTION_PREFERENCE)}
    types = _str(orbits["nss_solution_type"])
    r = np.array([rank.get(x, 99) for x in types])
    order = np.lexsort((r, o_sid))  # by source, then preference
    first = np.concatenate([[True], np.diff(o_sid[order]) != 0])
    best = order[first]
    pos = dict(zip(o_sid[best], best))
    return np.array([pos.get(int(s), -1) for s in np.asarray(source_ids, dtype=np.int64)])


def second_orbit_periods(orbits: Table) -> dict:
    """For sources with two DR3 orbits: {source_id: (types, periods)}."""
    o_sid = np.asarray(orbits["source_id"], dtype=np.int64)
    u, cnt = np.unique(o_sid, return_counts=True)
    dup = set(u[cnt > 1])
    types = _str(orbits["nss_solution_type"])
    per = _f(orbits["period"])
    out: dict = {}
    for i in np.where(np.isin(o_sid, list(dup)))[0]:
        out.setdefault(int(o_sid[i]), []).append((types[i], per[i]))
    return out


@lru_cache(maxsize=1)
def wide_pairs() -> Table:
    """El-Badry et al. (2021) pairs within 200 pc with R < 0.1 (both members
    have DR3 parallaxes; star 1 is always the brighter in G)."""
    e = Table.read(C.VIZIER_DIR / "elbadry21_200pc.fits")
    e = e[_f(e["R"]) < C.R_CHANCE_MAX]
    out = Table()
    out["source_id1"] = np.asarray(e["Source1"], dtype=np.int64)
    out["source_id2"] = np.asarray(e["Source2"], dtype=np.int64)
    for a, b in (("plx1", "Plx1"), ("plx2", "Plx2"), ("g1", "Gmag1"), ("g2", "Gmag2"), ("bprp1", None),
                 ("bprp2", None), ("sep_au", "sepAU"), ("r_chance", "R"), ("ruwe1", "RUWE1"),
                 ("ruwe2", "RUWE2")):
        if b is not None:
            out[a] = _f(e[b])
    out["bprp1"] = _f(e["BPmag1"]) - _f(e["RPmag1"])
    out["bprp2"] = _f(e["BPmag2"]) - _f(e["RPmag2"])
    out["theta_arcsec"] = _f(e["theta"]) * 3600.0
    out["bin_type"] = _str(e["BinType"])
    # absolute magnitudes on the more precise of the two parallaxes
    ep1, ep2 = _f(e["e_Plx1"]), _f(e["e_Plx2"])
    plx = np.where(ep1 <= ep2, out["plx1"], out["plx2"])
    out["plx"] = plx
    out["abs_g1"] = out["g1"] + 5 * np.log10(plx / 100)
    out["abs_g2"] = out["g2"] + 5 * np.log10(plx / 100)
    out["dg"] = out["g2"] - out["g1"]
    return out


def attach_wide(t: Table) -> None:
    """Number of resolved companions of each star, and the closest one."""
    w = wide_pairs()
    sid = np.asarray(t["source_id"], dtype=np.int64)
    index = {s: i for i, s in enumerate(sid)}
    n = np.zeros(len(t), int)
    sep = np.full(len(t), np.nan)
    comp_g = np.full(len(t), np.nan)
    comp_abs_g = np.full(len(t), np.nan)
    role = np.full(len(t), "", dtype="U9")
    for k in range(len(w)):
        for me, other, rl in ((w["source_id1"][k], w["g2"][k], "primary"),
                              (w["source_id2"][k], w["g1"][k], "secondary")):
            i = index.get(int(me))
            if i is None:
                continue
            n[i] += 1
            if not np.isfinite(sep[i]) or w["sep_au"][k] < sep[i]:
                sep[i] = w["sep_au"][k]
                comp_g[i] = other
                comp_abs_g[i] = w["abs_g2"][k] if rl == "primary" else w["abs_g1"][k]
                role[i] = rl
    t["n_wide"] = n
    t["wide_sep_au"] = sep
    t["wide_comp_g"] = comp_g
    t["wide_comp_abs_g"] = comp_abs_g
    t["wide_role"] = role


def solar_type_mask(t: Table) -> np.ndarray:
    """Raghavan-like F6-K3 dwarfs: colour window, main-sequence band, good parallax."""
    c = t["bp_rp"]
    return (c >= SOLAR_BPRP[0]) & (c < SOLAR_BPRP[1]) & t["ms_band"] & t["good_astrometry"]


def fit_grvs_colour(t: Table, deg: int = 4) -> np.ndarray:
    """G_RVS - G_RP as a polynomial in BP-RP, fitted to single-looking dwarfs of
    the 100 pc sample that have a DR3 grvs_mag (Sartoretti et al. 2023 give a
    similar published relation; fitting our own keeps it consistent)."""
    grvs = _f(t["grvs_mag"])
    grp = _f(t["phot_rp_mean_mag"])
    c = t["bp_rp"]
    ok = (np.isfinite(grvs) & np.isfinite(grp) & np.isfinite(c) & t["ms_band"]
          & (_f(t["ruwe"]) < 1.4) & (c > 0) & (c < 4))
    x, y = c[ok], (grvs - grp)[ok]
    coeffs = np.polyfit(x, y, deg)
    for _ in range(3):  # clip the few blends and variables
        r = y - np.polyval(coeffs, x)
        keep = np.abs(r) < 4 * 1.4826 * np.median(np.abs(r - np.median(r)))
        coeffs = np.polyfit(x[keep], y[keep], deg)
    return coeffs
