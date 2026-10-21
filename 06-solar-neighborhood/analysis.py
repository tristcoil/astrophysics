"""The pipeline: from archive tables to a classified 3D map, and what it shows.

1. **Clean Gaia.** Drop spurious astrometric solutions, then compute each star's
   distance, absolute magnitude, Galactic position and temperature.
2. **Fill Gaia's blind spot.** Find the stars Gaia has no usable measurement for
   -- mostly the brightest in the sky -- and take those from Hipparcos.
3. **Classify** every star from its place in the HR diagram.
4. **Measure** the census by type, what the naked eye actually sees, how
   completeness falls with distance, the nearest star cluster (found blind, not
   looked up), and whether we sit above or below the Galactic disc.
5. **Validate** against three independent references: the Gaia team's own
   nearby-star catalogue, the 10 pc census, and SIMBAD spectral types.

Every figure is drawn from the artefacts this module returns, never by
re-running any of it, so what is plotted is exactly what was measured.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone

import numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree

import config as C
import stellar_physics as sp

#: Classes in display order. Main-sequence dwarfs run hot to cool, then the
#: three populations that are not on the main sequence at all.
CLASS_ORDER = ["B", "A", "F", "G", "K", "M", "L/T", "WD", "Giant", "No colour"]
CLASS_LABELS = {
    "B": "B dwarfs", "A": "A dwarfs", "F": "F dwarfs", "G": "G dwarfs", "K": "K dwarfs",
    "M": "M dwarfs", "L/T": "L/T dwarfs", "WD": "White dwarfs", "Giant": "Giants & subgiants",
    "No colour": "No Gaia colour",
}

#: Hipparcos epoch 1991.25 to Gaia DR3 epoch 2016.0.
HIPPARCOS_TO_GAIA_YEARS = 2016.0 - 1991.25


def _f(col) -> np.ndarray:
    """A column as float, with masked entries as NaN."""
    return np.ma.filled(np.ma.asarray(col).astype(float), np.nan)


def _i(col, fill: int = -1) -> np.ndarray:
    """A column as int64, with masked entries as ``fill``."""
    return np.ma.filled(np.ma.asarray(col).astype(np.int64), fill)


def _s(col) -> np.ndarray:
    return np.array([str(v).strip() for v in np.ma.filled(np.ma.asarray(col).astype(str), "")])


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

def load_inputs() -> dict:
    """Read every cached dataset. Raises if fetch_data.py has not been run."""
    missing = [p.name for p in (C.GAIA_PATH, C.HIPPARCOS_PATH, C.GCNS_PATH, C.NAKED_EYE_PATH,
                                C.CENSUS_10PC_PATH, C.DWARF_SEQUENCE_PATH, C.FAMOUS_PATH,
                                C.NAMES_PATH, C.EXOPLANET_HOSTS_PATH) if not p.exists()]
    if missing:
        raise FileNotFoundError(f"missing cached data {missing}; run `python fetch_data.py` first")
    return {
        "gaia": Table.read(C.GAIA_PATH),
        "hipparcos": Table.read(C.HIPPARCOS_PATH),
        "gcns": Table.read(C.GCNS_PATH),
        "naked_eye": Table.read(C.NAKED_EYE_PATH),
        "census": Table.read(C.CENSUS_10PC_PATH),
        "sequence": sp.DwarfSequence.from_text(C.DWARF_SEQUENCE_PATH.read_text()),
        "famous": json.loads(C.FAMOUS_PATH.read_text()),
        "names": Table.read(C.NAMES_PATH),
        "exoplanets": Table.read(C.EXOPLANET_HOSTS_PATH),
    }


# --------------------------------------------------------------------------
# Classification
# --------------------------------------------------------------------------

def classify(teff, m_g, bp_rp, seq: sp.DwarfSequence) -> tuple[np.ndarray, np.ndarray]:
    """Assign every star a class from its place in the HR diagram.

    In order of precedence:

    * **White dwarfs** lie far below the main sequence, below a line through the
      empty gap separating the two sequences.
    * **L/T dwarfs** are anything fainter than M_G = 16.5, where Gaia's colours
      stop being trustworthy; most are brown dwarfs.
    * **Giants and subgiants** lie well above the main sequence: more than
      ``GIANT_OFFSET_MAG`` brighter than a dwarf of the same temperature.
    * Everything else is a **main-sequence dwarf**, and gets the spectral letter
      of its temperature.

    Returns (class, subtype, teff); subtype is the nearest tabulated dwarf
    subtype and is only meaningful for the main-sequence classes. The returned
    temperature equals the input except for colourless late-M dwarfs, which get
    one from their absolute magnitude.
    """
    teff, m_g, bp_rp = (np.asarray(a, dtype=float) for a in (teff, m_g, bp_rp))
    cls = np.full(len(m_g), "No colour", dtype="<U9")
    has_teff = np.isfinite(teff)

    wd = (np.isfinite(bp_rp) & (bp_rp < C.WD_MAX_BP_RP)
          & (m_g > C.WD_LINE_SLOPE * bp_rp + C.WD_LINE_INTERCEPT))
    ultracool = (m_g > C.ULTRACOOL_MG) & ~wd
    # No colour, but too faint to be anything except a late-M dwarf.
    late_m = ~has_teff & ~np.isfinite(bp_rp) & (m_g > C.LATE_M_MG) & ~ultracool
    teff = np.where(late_m, seq.teff_from_mg(m_g), teff)
    has_teff = np.isfinite(teff)
    ms_mg = seq.main_sequence_mg_from_teff(teff)
    giant = (has_teff & ~wd & ~ultracool & np.isfinite(ms_mg)
             & (m_g < ms_mg - C.GIANT_OFFSET_MAG) & (m_g < C.GIANT_MAX_MG))
    dwarf = has_teff & ~wd & ~ultracool & ~giant

    letters = seq.spectral_class(teff)
    cls[dwarf] = letters[dwarf]
    cls[wd] = "WD"
    cls[ultracool] = "L/T"
    cls[giant] = "Giant"
    subtype = np.where(dwarf, seq.spectral_subtype(teff), "")
    return cls, subtype, teff


def display_teff(cls, teff, m_g, seq: sp.DwarfSequence) -> np.ndarray:
    """Temperature used only to colour each star on screen.

    Dwarf-sequence temperatures are the physical estimate. Two classes need a
    fallback: L/T dwarfs, whose Gaia colours are unreliable, take the dwarf
    temperature that matches their absolute magnitude; stars with no colour at
    all get none (and are drawn grey).
    """
    t = np.array(teff, dtype=float)
    ultracool = cls == "L/T"
    t[ultracool] = seq.teff_from_mg(m_g[ultracool])
    lo, hi = np.nanmin(seq.teff[np.isfinite(seq.bp_rp)]), np.nanmax(seq.teff[np.isfinite(seq.bp_rp)])
    # Hot white dwarfs are bluer than the hottest tabulated dwarf colour.
    hot_wd = (cls == "WD") & ~np.isfinite(t)
    t[hot_wd] = hi
    return np.clip(t, 300.0, None)


# --------------------------------------------------------------------------
# 1. Gaia
# --------------------------------------------------------------------------

def build_gaia(gaia: Table, seq: sp.DwarfSequence) -> tuple[Table, dict]:
    """Apply the fidelity cut and derive everything we need per Gaia star."""
    fidelity = _f(gaia["fidelity_v2"])
    good = fidelity > C.MIN_FIDELITY
    l_all, b_all = _f(gaia["l"]), _f(gaia["b"])
    ruwe = _f(gaia["ruwe"])

    def bulge_fraction(mask):
        """Share of sources toward the crowded inner Galaxy (|l| < 30, |b| < 10)."""
        l180 = np.where(l_all > 180, l_all - 360, l_all)
        inner = (np.abs(l180) < 30) & (np.abs(b_all) < 10)
        return float(inner[mask].mean())

    # The patch |l| < 30, |b| < 10 is 2.9% of the sky. A real 100 pc sample is
    # isotropic, so it should hold ~2.9% of the stars -- whatever piles up
    # beyond that is spurious.
    l_edges = np.radians([-30, 30])
    sky_share = float((l_edges[1] - l_edges[0]) * 2 * np.sin(np.radians(10)) / (4 * np.pi))
    diag = {
        "n_candidates": len(gaia),
        "n_rejected_fidelity": int((~good).sum()),
        "n_kept": int(good.sum()),
        "inner_galaxy_sky_share": sky_share,
        "inner_galaxy_share_rejected": bulge_fraction(~good),
        "inner_galaxy_share_kept": bulge_fraction(good),
        "kept_with_ruwe_above_1p4": int((good & (ruwe >= 1.4)).sum()),
        "kept_fraction_ruwe_above_1p4": float((ruwe[good] >= 1.4).mean()),
    }

    g = gaia[good]
    plx = _f(g["parallax"])
    dist = sp.distance_pc(plx)
    x, y, z, l, b = sp.galactic_xyz(_f(g["ra"]), _f(g["dec"]), dist)
    # Cross-check our astropy conversion against Gaia's own catalogue l, b.
    dl = (l - _f(g["l"]) + 180) % 360 - 180
    sep = np.hypot(dl * np.cos(np.radians(b)), b - _f(g["b"])) * 3600
    diag["lb_max_disagreement_arcsec"] = float(np.nanmax(sep))

    gmag = _f(g["phot_g_mean_mag"])
    bp_rp = _f(g["bp_rp"])
    m_g = sp.absolute_magnitude(gmag, plx)

    # Contaminated colours: flagged by the excess factor AND impossibly far below
    # the main sequence (see config). Their temperature comes from brightness.
    cstar, sig = sp.corrected_excess_factor(_f(g["phot_bp_rp_excess_factor"]), bp_rp, gmag)
    in_wd_region = (bp_rp < C.WD_MAX_BP_RP) & (m_g > C.WD_LINE_SLOPE * bp_rp + C.WD_LINE_INTERCEPT)
    contaminated = ((np.abs(cstar) > C.CONTAMINATION_SIGMA * sig)
                    & (m_g > seq.main_sequence_mg(bp_rp) + C.CONTAMINATED_BELOW_MS_MAG) & ~in_wd_region)
    bp_rp_use = np.where(contaminated, np.nan, bp_rp)
    teff = seq.teff_from_bp_rp(bp_rp_use)
    teff[contaminated] = seq.teff_from_mg(m_g[contaminated])
    cls, subtype, teff = classify(teff, m_g, bp_rp_use, seq)
    v_est = gmag - seq.g_minus_v_from_bp_rp(bp_rp_use)
    below_ms = (m_g > seq.main_sequence_mg(bp_rp) + C.CONTAMINATED_BELOW_MS_MAG) & ~in_wd_region
    diag["n_below_main_sequence"] = int(below_ms.sum())
    diag["n_colour_contaminated"] = int(contaminated.sum())

    cat = Table()
    cat["origin"] = np.full(len(g), "Gaia DR3", dtype="<U9")
    cat["source_id"] = _i(g["source_id"])
    cat["hip"] = np.full(len(g), -1, dtype=np.int64)
    cat["ra"], cat["dec"] = _f(g["ra"]), _f(g["dec"])
    cat["parallax"], cat["parallax_error"] = plx, _f(g["parallax_error"])
    cat["dist_pc"] = dist
    cat["x"], cat["y"], cat["z"], cat["l"], cat["b"] = x, y, z, l, b
    cat["pmra"], cat["pmdec"] = _f(g["pmra"]), _f(g["pmdec"])
    cat["rv"] = _f(g["radial_velocity"])
    cat["g_mag"], cat["bp_rp"], cat["m_g"] = gmag, bp_rp, m_g
    cat["v_mag"] = v_est
    cat["v_measured"] = np.zeros(len(g), dtype=bool)
    cat["teff"] = teff
    cat["cls"], cat["subtype"] = cls, subtype
    cat["ruwe"], cat["fidelity"] = _f(g["ruwe"]), _f(g["fidelity_v2"])
    cat["colour_contaminated"] = contaminated
    cat["in_gcns"] = ~np.ma.getmaskarray(np.ma.asarray(g["gcns_prob"]))
    cat["gcns_wd_prob"] = _f(g["wd_prob"])
    return cat, diag


# --------------------------------------------------------------------------
# 2. Hipparcos, where Gaia cannot help
# --------------------------------------------------------------------------

def build_hipparcos_supplement(hip: Table, gaia_raw: Table, gaia_cat: Table,
                               seq: sp.DwarfSequence) -> tuple[Table, dict, dict]:
    """Add the stars Gaia has no usable measurement for, from Hipparcos.

    A Hipparcos star counts as already covered by Gaia if any of these hold:

    * its official Gaia cross-match partner is in our cleaned Gaia sample;
    * its partner has a precise Gaia parallax placing it *beyond* 100 pc
      (Gaia is ~100x more precise, so Gaia wins that argument);
    * the official cross-match left it unpaired, but a Gaia source of similar
      brightness sits within 3 arcsec of where Hipparcos says the star should
      be in 2016.

    Whatever is left is added from Hipparcos, with the same 10% parallax cut.

    Also returns a map from Gaia source_id to Hipparcos V for the covered stars:
    for naked-eye questions a measured V beats one estimated from Gaia colours.
    """
    plx, e_plx = _f(hip["plx"]), _f(hip["e_plx"])
    hip_poe = plx / e_plx
    gid = _i(hip["gaia_source_id"])
    g_plx, g_poe = _f(hip["gaia_parallax"]), _f(hip["gaia_poe"])
    hp = _f(hip["hp_mag"])

    good_ids = np.asarray(gaia_cat["source_id"])
    in_sample = np.isin(gid, good_ids) & (gid > 0)
    gaia_says_beyond = (gid > 0) & (g_poe > C.MIN_PARALLAX_OVER_ERROR) & (g_plx <= C.MIN_PARALLAX_MAS)
    covered = in_sample | gaia_says_beyond

    # Positional rescue for the stars the official cross-match left unpaired or
    # paired with a Gaia source that has no usable astrometry.
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    todo = np.flatnonzero(~covered)
    ra16, dec16 = sp.propagate_position(_f(hip["ra"])[todo], _f(hip["dec"])[todo],
                                        _f(hip["pm_ra"])[todo], _f(hip["pm_de"])[todo],
                                        HIPPARCOS_TO_GAIA_YEARS)
    gaia_coords = SkyCoord(np.asarray(gaia_cat["ra"]) * u.deg, np.asarray(gaia_cat["dec"]) * u.deg)
    idx, sep, _ = SkyCoord(ra16 * u.deg, dec16 * u.deg).match_to_catalog_sky(gaia_coords)
    same_star = (sep.arcsec < C.HIPPARCOS_MATCH_RADIUS_ARCSEC) & \
                (np.abs(np.asarray(gaia_cat["g_mag"])[idx] - hp[todo]) < 1.5)
    rescued = np.zeros(len(hip), dtype=bool)
    rescued[todo[same_star]] = True
    covered |= rescued

    within = (plx > C.MIN_PARALLAX_MAS) & (hip_poe > C.MIN_HIPPARCOS_PARALLAX_OVER_ERROR)
    add = within & ~covered

    # Measured V for Gaia stars that Hipparcos also observed.
    vmag_all = _f(hip["vmag"])
    gaia_v = {}
    for sid, v in zip(gid[in_sample], vmag_all[in_sample]):
        if np.isfinite(v):
            gaia_v[int(sid)] = float(v)
    for k, j in zip(todo[same_star], idx[same_star]):
        if np.isfinite(vmag_all[k]):
            gaia_v[int(gaia_cat["source_id"][j])] = float(vmag_all[k])
    gaia_hip = {int(s): int(h) for s, h in zip(gid[in_sample], _i(hip["hip"])[in_sample])}
    for k, j in zip(todo[same_star], idx[same_star]):
        gaia_hip[int(gaia_cat["source_id"][j])] = int(hip["hip"][k])
    bv_all = _f(hip["b_v"])
    hip_row = {int(h): k for k, h in enumerate(_i(hip["hip"]))}
    gaia_bv = {s: float(bv_all[hip_row[h]]) for s, h in gaia_hip.items() if np.isfinite(bv_all[hip_row[h]])}

    h = hip[add]
    plx_h, e_plx_h = plx[add].copy(), e_plx[add].copy()

    # Physical pairs share one distance. Where two or more added components of
    # the same Hipparcos double (CCDM number) have parallaxes consistent within
    # PAIR_CONSISTENCY_SIGMA, give each the error-weighted mean. The consistency
    # test keeps chance line-of-sight doubles apart. This matters most for
    # Alpha Centauri, whose separate 2007 parallaxes (754.8 and 796.9 mas) would
    # otherwise put B nearer to us than Proxima.
    ccdm = _s(h["ccdm"])
    n_pairs = 0
    for ident in set(ccdm) - {""}:
        k = np.flatnonzero(ccdm == ident)
        if len(k) < 2:
            continue
        w = 1 / e_plx_h[k] ** 2
        mean = float((plx_h[k] * w).sum() / w.sum())
        chi2 = float((((plx_h[k] - mean) / e_plx_h[k]) ** 2).sum())
        if chi2 <= C.PAIR_CONSISTENCY_SIGMA**2 * (len(k) - 1):
            plx_h[k], e_plx_h[k] = mean, 1 / np.sqrt(w.sum())
            n_pairs += 1
    dist = sp.distance_pc(plx_h)
    ra, dec = _f(h["ra"]), _f(h["dec"])
    x, y, z, l, b = sp.galactic_xyz(ra, dec, dist)
    b_v = _f(h["b_v"])
    vmag = _f(h["vmag"])
    teff = seq.teff_from_b_v(b_v)
    g_est = vmag + seq.g_minus_v_from_b_v(b_v)
    # The very hottest stars are bluer than the G-V column reaches; G ~ V there.
    g_est = np.where(np.isfinite(g_est), g_est, vmag)
    bp_rp_est = seq.bp_rp_from_b_v(b_v)
    m_g = sp.absolute_magnitude(g_est, plx_h)
    cls, subtype, teff = classify(teff, m_g, bp_rp_est, seq)

    cat = Table()
    cat["origin"] = np.full(len(h), "Hipparcos", dtype="<U9")
    cat["source_id"] = _i(h["gaia_source_id"])
    cat["hip"] = _i(h["hip"])
    cat["ra"], cat["dec"] = ra, dec
    cat["parallax"], cat["parallax_error"] = plx_h, e_plx_h
    cat["dist_pc"] = dist
    cat["x"], cat["y"], cat["z"], cat["l"], cat["b"] = x, y, z, l, b
    cat["pmra"], cat["pmdec"] = _f(h["pm_ra"]), _f(h["pm_de"])
    cat["rv"] = np.full(len(h), np.nan)
    cat["g_mag"], cat["bp_rp"], cat["m_g"] = g_est, bp_rp_est, m_g
    cat["v_mag"] = vmag
    cat["v_measured"] = np.isfinite(vmag)
    cat["teff"] = teff
    cat["cls"], cat["subtype"] = cls, subtype
    cat["ruwe"] = np.full(len(h), np.nan)
    cat["fidelity"] = np.full(len(h), np.nan)
    cat["colour_contaminated"] = np.zeros(len(h), dtype=bool)
    cat["in_gcns"] = np.zeros(len(h), dtype=bool)
    cat["gcns_wd_prob"] = np.full(len(h), np.nan)

    # How good are Hipparcos-based G and BP-RP estimates? Test them on stars both
    # missions measured well, avoiding the saturated bright end.
    both = in_sample & (vmag_all > 6) & (vmag_all < 10)
    g_true = _f(hip["gaia_g"])[both]
    g_pred = vmag_all[both] + seq.g_minus_v_from_b_v(_f(hip["b_v"])[both])
    resid = g_pred - g_true
    resid = resid[np.isfinite(resid)]

    vmag_add = vmag[np.isfinite(vmag)]
    diag = {
        "n_hipparcos_plx_gt_10": len(hip),
        "n_within_100pc_precise": int(within.sum()),
        "covered_by_crossmatch": int((within & in_sample).sum()),
        "gaia_places_beyond_100pc": int((within & gaia_says_beyond & ~in_sample).sum()),
        "covered_by_position_match": int((within & rescued).sum()),
        "added_from_hipparcos": int(add.sum()),
        "added_brighter_than_v4": int((vmag_add < 4).sum()),
        "added_v4_to_7": int(((vmag_add >= 4) & (vmag_add < 7)).sum()),
        "added_fainter_than_v7": int((vmag_add >= 7).sum()),
        "added_with_gaia_counterpart_no_parallax": int((add & (gid > 0) & ~np.isfinite(g_plx)).sum()),
        "added_physical_pairs_merged_distance": n_pairs,
        "added_in_hipparcos_multiple_systems": int(np.isin(_s(h["multflag"]), ["C", "G", "O", "V", "X"]).sum()),
        "g_estimate_bias_mag": float(np.median(resid)),
        "g_estimate_scatter_mag": float(1.4826 * np.median(np.abs(resid - np.median(resid)))),
        "g_estimate_n_stars": int(len(resid)),
    }
    return cat, diag, {"gaia_v": gaia_v, "gaia_hip": gaia_hip, "gaia_bv": gaia_bv}


# --------------------------------------------------------------------------
# Assembling the catalogue
# --------------------------------------------------------------------------

def _preferred_name(common: str, main_id: str) -> str:
    """Pick the most human name SIMBAD knows for a star.

    Prefers a common name ('Barnard's star') over a catalogue number, and among
    common names the one with fewest digits, then the longest -- which picks
    'Proxima Centauri' over 'Proxima' and 'Luhman 16' over 'WISE J1049-5319AB'.
    """
    names = [n for n in common.split("|") if n]
    if not names:
        return main_id
    return sorted(names, key=lambda n: (sum(ch.isdigit() for ch in n), -len(n)))[0]


def attach_names(cat: Table, names: Table) -> None:
    """Give every star within the named radius its SIMBAD name and spectral type.

    Names are stored exactly as SIMBAD spells them ('* alf CMa'), which is plain
    ASCII and so survives a FITS file; ``stellar_physics.tidy_simbad_name``
    turns them into 'α CMa' at display time.
    """
    # Read through _s(): FITS stores an empty string as a masked value, whose
    # str() is the literal text '--' -- which would otherwise become a name.
    lookup = {}
    for ident, main_id, common, spt, otype in zip(_s(names["ident"]), _s(names["main_id"]),
                                                  _s(names["common_names"]), _s(names["sp_type"]),
                                                  _s(names["otype"])):
        lookup[ident] = (main_id, common, spt, otype)
    n = len(cat)
    name = np.full(n, "", dtype="<U40")
    main = np.full(n, "", dtype="<U40")
    sptype = np.full(n, "", dtype="<U24")
    near = np.flatnonzero(np.asarray(cat["dist_pc"]) <= C.NAMED_RADIUS_PC)
    for i in near:
        keys = []
        if cat["origin"][i] == "Gaia DR3":
            keys.append(f"Gaia DR3 {cat['source_id'][i]}")
        if cat["hip"][i] > 0:
            keys.append(f"HIP {cat['hip'][i]}")
        for key in keys:
            if key in lookup:
                m, common, spt, _ = lookup[key]
                name[i] = _preferred_name(common, m)[:40]
                main[i] = m[:40]
                sptype[i] = spt[:24]
                break
    cat["name"], cat["simbad_id"], cat["simbad_sptype"] = name, main, sptype


def build_catalogue(inputs: dict) -> tuple[Table, dict]:
    seq = inputs["sequence"]
    gaia_cat, gaia_diag = build_gaia(inputs["gaia"], seq)
    hip_cat, hip_diag, links = build_hipparcos_supplement(inputs["hipparcos"], inputs["gaia"], gaia_cat, seq)

    # Measured Hipparcos V (and HIP number) for Gaia stars Hipparcos also saw.
    sids = np.asarray(gaia_cat["source_id"])
    v_meas = np.array([links["gaia_v"].get(int(s), np.nan) for s in sids])
    has_v = np.isfinite(v_meas)
    gaia_cat["v_mag"][has_v] = v_meas[has_v]
    gaia_cat["v_measured"][has_v] = True
    gaia_cat["hip"] = np.array([links["gaia_hip"].get(int(s), -1) for s in sids], dtype=np.int64)

    # The brightest Gaia stars: temperature from Hipparcos B-V (see BRIGHT_COLOUR_G).
    bv = np.array([links["gaia_bv"].get(int(s), np.nan) for s in sids])
    bright = (np.asarray(gaia_cat["g_mag"]) < C.BRIGHT_COLOUR_G) & np.isfinite(bv)
    if bright.any():
        teff_bv = seq.teff_from_b_v(bv[bright])
        cls, subtype, teff_bv = classify(teff_bv, np.asarray(gaia_cat["m_g"])[bright],
                                         seq.bp_rp_from_b_v(bv[bright]), seq)
        gaia_cat["teff"][bright] = teff_bv
        gaia_cat["cls"][bright] = cls
        gaia_cat["subtype"][bright] = subtype
    build_diag_bright = int(bright.sum())

    from astropy.table import vstack

    cat = vstack([gaia_cat, hip_cat])
    cat["display_teff"] = display_teff(np.asarray(cat["cls"]), np.asarray(cat["teff"]),
                                       np.asarray(cat["m_g"]), seq)
    attach_names(cat, inputs["names"])
    gaia_diag["bright_stars_with_hipparcos_colour"] = build_diag_bright
    return cat, {"gaia": gaia_diag, "hipparcos": hip_diag}


# --------------------------------------------------------------------------
# 4. Measurements
# --------------------------------------------------------------------------

def census(cat: Table) -> dict:
    """How many of each kind of star, within several radii, and which are visible."""
    d = np.asarray(cat["dist_pc"])
    cls = np.asarray(cat["cls"])
    v = np.asarray(cat["v_mag"])
    out = {"by_radius": {}, "naked_eye": {}}
    for radius in (10, 20, 50, 100):
        sel = d <= radius
        counts = {c: int((sel & (cls == c)).sum()) for c in CLASS_ORDER}
        out["by_radius"][str(radius)] = {"total": int(sel.sum()), "counts": counts}
        eye = sel & (v < C.NAKED_EYE_V)
        out["naked_eye"][str(radius)] = {
            "total": int(eye.sum()),
            "counts": {c: int((eye & (cls == c)).sum()) for c in CLASS_ORDER},
        }
    stars20 = out["by_radius"]["20"]
    stellar = stars20["total"] - stars20["counts"]["L/T"] - stars20["counts"]["No colour"]
    out["m_dwarf_share_of_stars_20pc"] = stars20["counts"]["M"] / stellar
    return out


def naked_eye_sky(naked_eye: Table) -> dict:
    """What fraction of the stars you can see actually lie inside the map?"""
    v = _f(naked_eye["vmag"])
    plx = _f(naked_eye["plx"])
    visible = v < C.NAKED_EYE_V
    d = sp.distance_pc(np.where(plx > 0, plx, np.nan))
    inside = visible & (plx > C.MIN_PARALLAX_MAS)
    finite = visible & np.isfinite(d)
    return {
        "n_visible": int(visible.sum()),
        "n_inside_map": int(inside.sum()),
        "fraction_inside_map": float(inside.sum() / visible.sum()),
        "median_distance_pc": float(np.nanmedian(d[finite])),
        "fraction_beyond_300pc": float((d[finite] > 300).mean()),
        "note": "Hipparcos parallaxes; beyond ~300 pc their precision is poor, so the "
                "median is robust but individual large distances are not.",
    }


def gaia_limits(cat: Table) -> dict:
    """Gaia's effective bright and faint limits, read off the sample itself."""
    g = np.asarray(cat["g_mag"])[np.asarray(cat["origin"]) == "Gaia DR3"]
    g = g[np.isfinite(g)]
    return {"g_faint": float(np.percentile(g, 99.9)), "g_bright": float(g.min()),
            "g_faint_definition": "99.9th percentile of G in the cleaned Gaia sample",
            "g_bright_definition": "brightest G in the cleaned Gaia sample"}


def density_by_class(cat: Table, cluster_mask: np.ndarray) -> dict:
    """Number density of each class in distance shells.

    If a class were detected completely and stars were spread uniformly, its
    density would be flat with distance. Two things bend it: a *gentle* decline
    for every class because the Galactic disc is thin (stars thin out above and
    below the plane, and a bigger sphere reaches further out of it), and a
    *cliff* where a class grows too faint for Gaia to see.
    """
    edges = np.arange(0.0, C.MAX_DISTANCE_PC + 1, 10.0)
    vol = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    d = np.asarray(cat["dist_pc"])
    cls = np.asarray(cat["cls"])
    keep = ~cluster_mask
    out = {"edges_pc": edges.tolist(), "classes": {}}
    for c in CLASS_ORDER:
        n, _ = np.histogram(d[keep & (cls == c)], bins=edges)
        out["classes"][c] = {"n": n.tolist(), "density_per_pc3": (n / vol).tolist(),
                             "err_per_pc3": (np.sqrt(n) / vol).tolist()}
    return out


def completeness_distances(seq: sp.DwarfSequence, g_faint: float) -> dict:
    """How far Gaia can see each dwarf subtype, from the dwarf sequence alone."""
    out = {}
    for spt, mg in zip(seq.spt, seq.m_g):
        if np.isfinite(mg) and spt[0] in "MLT" and spt.endswith("V") and ("." not in spt):
            out[spt] = float(sp.limiting_distance_pc(mg, g_faint))
    return out


# --------------------------------------------------------------------------
# Finding structure: the nearest star cluster, found blind
# --------------------------------------------------------------------------

def find_overdensities(xyz: np.ndarray, *, n_peaks: int = 2, seed: int = 1) -> tuple[dict, np.ndarray]:
    """Find the densest clumps of stars in 3D, and test whether they are real.

    For every star, count its neighbours within 3 pc. The star with the most
    neighbours marks the densest spot; a mean-shift walk (move to the centroid
    of the neighbours, repeat) settles on the clump's centre. Stars within the
    cluster radius are then masked and the search repeats.

    Is a peak real? A random scatter of 350,000 stars always has *some* crowded
    spot, so the peak is compared with a null: the same stars at the same
    distances, but pointed in random directions. That keeps every distance
    (and so all of Gaia's distance-dependent selection) while destroying any
    real clustering. A real cluster must beat the best clump the null can make.
    """
    r = C.CLUSTER_SEARCH_RADIUS_PC
    tree = cKDTree(xyz)
    counts = tree.query_ball_point(xyz, r, return_length=True) - 1
    dist = np.linalg.norm(xyz, axis=1)

    rng = np.random.default_rng(seed)
    null_max, null_hist = [], None
    top = int(counts.max()) + 1
    for _ in range(3):
        u = rng.normal(size=xyz.shape)
        fake = u / np.linalg.norm(u, axis=1, keepdims=True) * dist[:, None]
        null_counts = cKDTree(fake).query_ball_point(fake, r, return_length=True) - 1
        null_max.append(int(null_counts.max()))
        h = np.bincount(null_counts, minlength=top)[:top]
        null_hist = h if null_hist is None else null_hist + h
    null_hist = null_hist / 3.0

    alive = np.ones(len(xyz), dtype=bool)
    peaks = []
    for _ in range(n_peaks):
        i = int(np.argmax(np.where(alive, counts, -1)))
        centre = xyz[i].copy()
        for _step in range(50):
            members = tree.query_ball_point(centre, r)
            new = xyz[members].mean(axis=0)
            if np.linalg.norm(new - centre) < 1e-3:
                break
            centre = new
        d_c = float(np.linalg.norm(centre))
        n_core = len(tree.query_ball_point(centre, r))
        n_wide = len(tree.query_ball_point(centre, C.CLUSTER_MASK_RADIUS_PC))
        # Field density at this distance from the Sun, away from the clump.
        shell = (np.abs(dist - d_c) < 5.0) & (np.linalg.norm(xyz - centre, axis=1) > 2 * C.CLUSTER_MASK_RADIUS_PC)
        shell_vol = 4 / 3 * np.pi * ((d_c + 5) ** 3 - (d_c - 5) ** 3)
        # Remove the excluded balls' volume from the shell volume, approximately.
        shell_vol -= 4 / 3 * np.pi * (2 * C.CLUSTER_MASK_RADIUS_PC) ** 3 * (10.0 / (4 * C.CLUSTER_MASK_RADIUS_PC))
        field_density = shell.sum() / shell_vol
        exp_core = field_density * 4 / 3 * np.pi * r**3
        exp_wide = field_density * 4 / 3 * np.pi * C.CLUSTER_MASK_RADIUS_PC**3
        lon = float(np.degrees(np.arctan2(centre[1], centre[0])) % 360)
        lat = float(np.degrees(np.arcsin(centre[2] / d_c)))
        peaks.append({
            "centre_xyz_pc": centre.tolist(), "distance_pc": d_c, "l_deg": lon, "b_deg": lat,
            "peak_neighbours": int(counts[i]), "n_within_search_radius": n_core,
            "n_within_mask_radius": n_wide, "field_expected_search_radius": float(exp_core),
            "field_expected_mask_radius": float(exp_wide),
            "excess_within_mask_radius": float(n_wide - exp_wide),
            "beats_null": bool(counts[i] > max(null_max)),
        })
        alive &= np.linalg.norm(xyz - centre, axis=1) > C.CLUSTER_MASK_RADIUS_PC
    hist = {"n_neighbours": list(range(top)), "real": np.bincount(counts, minlength=top).tolist(),
            "null_mean_of_3": null_hist.tolist()}
    return {"search_radius_pc": r, "null_max_neighbours": null_max, "peaks": peaks,
            "histogram": hist}, counts


def identify_peaks(peaks: list, famous: dict) -> None:
    """Label each detected clump with the catalogued cluster it coincides with."""
    for peak in peaks:
        peak["match"] = None
        for label, simbad_name in C.REFERENCE_CLUSTERS.items():
            ref = famous.get(simbad_name)
            if not ref or not ref["plx_mas"]:
                continue
            d_ref = 1000.0 / ref["plx_mas"]
            x, y, z, l, b = sp.galactic_xyz(ref["ra"], ref["dec"], d_ref)
            offset = float(np.linalg.norm(np.array(peak["centre_xyz_pc"]) - np.array([x, y, z]).ravel()))
            if offset < C.CLUSTER_MASK_RADIUS_PC:
                peak["match"] = {
                    "name": label, "simbad_id": ref["main_id"], "literature_distance_pc": d_ref,
                    "literature_plx_mas": ref["plx_mas"], "reference": ref["plx_bibcode"],
                    "offset_pc": offset,
                    "distance_difference_percent": 100 * (peak["distance_pc"] - d_ref) / d_ref,
                }


# --------------------------------------------------------------------------
# Are we above the Galactic disc?
# --------------------------------------------------------------------------

#: Literature values used only to draw a *prediction* to compare against, never
#: fitted: the Sun's height above the mid-plane (Bennett & Bovy 2019, from Gaia
#: DR2 star counts over a far larger volume) and the thin disc's scale height
#: (Bland-Hawthorn & Gerhard 2016 review: 300 +- 50 pc).
LITERATURE_SUN_HEIGHT_PC = 20.8
LITERATURE_SCALE_HEIGHT_PC = (250.0, 300.0, 350.0)


def isothermal_asymmetry(z, sun_height, scale_height):
    """South/north density ratio at height |z| for an isothermal disc.

    The self-gravitating isothermal sheet has n(z) ~ sech^2(z / 2h). If the Sun
    sits a height z_sun above its mid-plane, a point |z| below us is |z| - z_sun
    from the plane and a point |z| above is |z| + z_sun, so the ratio of the two
    densities is [cosh((|z| + z_sun)/2h) / cosh((|z| - z_sun)/2h)]^2.
    """
    z = np.asarray(z, dtype=float)
    return (np.cosh((z + sun_height) / (2 * scale_height)) / np.cosh((z - sun_height) / (2 * scale_height))) ** 2


def _slab_volume(radius, lo, hi):
    """Volume of a sphere (centred on the Sun) between heights lo and hi, 0 <= lo < hi."""
    def f(z):
        return radius**2 * z - z**3 / 3
    return np.pi * (f(min(hi, radius)) - f(min(lo, radius)))


def _ball_slab_overlap(centre_z, radius, lo, hi):
    """Volume of a ball of this radius, centred at height centre_z, between heights lo and hi."""
    u1, u2 = max(lo, centre_z - radius) - centre_z, min(hi, centre_z + radius) - centre_z
    if u2 <= u1:
        return 0.0
    return np.pi * ((radius**2 * u2 - u2**3 / 3) - (radius**2 * u1 - u1**3 / 3))


def north_south(cat: Table, cluster_mask: np.ndarray, clusters: list[dict]) -> dict:
    """Compare star counts below the Sun with those above, height by height.

    The obvious approach -- fit a density profile to star counts versus height --
    fails inside a 100 pc sphere, because counts also fall with *distance* for
    a reason unrelated to the disc: at larger distances more binary companions
    blur into single Gaia sources (in the Galactic plane, where the disc's
    vertical structure plays no role, density still drops ~6% between 15 and
    95 pc). A slab at height -z and a slab at +z, though, have identical volume
    and identical distance distributions inside a sphere, so any selection that
    depends only on distance cancels exactly in their ratio.

    Only dwarfs brighter than M_G = 13 are counted: Gaia sees them to far
    beyond 100 pc, so incompleteness cannot distort the counts either.
    """
    x, y, z = (np.asarray(cat[k]) for k in ("x", "y", "z"))
    m_g = np.asarray(cat["m_g"])
    cls = np.asarray(cat["cls"])
    sel = np.isin(cls, ["A", "F", "G", "K", "M"]) & (m_g < 13.0) & ~cluster_mask
    edges = np.arange(0.0, C.MAX_DISTANCE_PC + 0.1, 10.0)
    centres = 0.5 * (edges[1:] + edges[:-1])
    north = np.array([(sel & (z >= lo) & (z < hi)).sum() for lo, hi in zip(edges[:-1], edges[1:])])
    south = np.array([(sel & (-z >= lo) & (-z < hi)).sum() for lo, hi in zip(edges[:-1], edges[1:])])
    # The masked cluster balls remove volume from particular slabs (the Hyades
    # from the south, Coma Berenices from the northern cap), so compare
    # densities over the unmasked volume, not raw counts. Both balls lie wholly
    # inside the 100 pc sphere, so the overlap is exact.
    vol_n, vol_s = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        full = _slab_volume(C.MAX_DISTANCE_PC, lo, hi)
        cut_n = sum(_ball_slab_overlap(c["centre_xyz_pc"][2], C.CLUSTER_MASK_RADIUS_PC, lo, hi) for c in clusters)
        cut_s = sum(_ball_slab_overlap(-c["centre_xyz_pc"][2], C.CLUSTER_MASK_RADIUS_PC, lo, hi) for c in clusters)
        vol_n.append(full - cut_n)
        vol_s.append(full - cut_s)
    vol_n, vol_s = np.array(vol_n), np.array(vol_s)
    ratio = (south / vol_s) / (north / vol_n)
    err = ratio * np.sqrt(1 / south + 1 / north)

    # One summary number: the pooled excess away from the plane of the Sun
    # itself (|z| > 10 pc), up to 90 pc where the slabs still hold >1000 stars.
    pool = (centres > 10) & (centres < 90)
    w = 1 / err[pool] ** 2
    mean = float((ratio[pool] * w).sum() / w.sum())
    mean_err = float(1 / np.sqrt(w.sum()))

    # In-plane density versus distance: the selection effect that makes a
    # direct profile fit fail, measured where the disc cannot cause it.
    r_cyl = np.hypot(x, y)
    ring_edges = np.arange(10.0, C.MAX_DISTANCE_PC + 0.1, 10.0)
    ring_density = []
    for lo, hi in zip(ring_edges[:-1], ring_edges[1:]):
        n = (sel & (np.abs(z) < 15) & (r_cyl >= lo) & (r_cyl < hi)).sum()
        ring_density.append(n / (np.pi * (hi**2 - lo**2) * 30.0))

    return {
        "bin_centres_pc": centres.tolist(), "n_north": north.tolist(), "n_south": south.tolist(),
        "volume_north_pc3": vol_n.tolist(), "volume_south_pc3": vol_s.tolist(),
        "ratio_south_north": ratio.tolist(), "ratio_err": err.tolist(),
        "pooled_ratio_10_90pc": mean, "pooled_ratio_err": mean_err,
        "pooled_significance_sigma": (mean - 1) / mean_err,
        "prediction": {
            "sun_height_pc": LITERATURE_SUN_HEIGHT_PC,
            "scale_heights_pc": list(LITERATURE_SCALE_HEIGHT_PC),
            "ratio": {str(h): isothermal_asymmetry(centres, LITERATURE_SUN_HEIGHT_PC, h).tolist()
                      for h in LITERATURE_SCALE_HEIGHT_PC},
        },
        "in_plane_ring_edges_pc": ring_edges.tolist(),
        "in_plane_density_per_pc3": ring_density,
        "in_plane_density_drop_percent": 100 * (1 - ring_density[-1] / ring_density[0]),
        "n_stars": int(sel.sum()),
    }


# --------------------------------------------------------------------------
# 5. Validation
# --------------------------------------------------------------------------

def compare_gcns(cat: Table, gcns: Table, gaia_raw: Table) -> dict:
    """Our simple cuts against the Gaia team's carefully curated 100 pc catalogue."""
    ours = cat[np.asarray(cat["origin"]) == "Gaia DR3"]
    ours_ids = np.asarray(ours["source_id"])
    in_gcns = np.asarray(ours["in_gcns"])
    d50 = _f(gcns["dist_50"]) * 1000.0
    within = d50 < C.MAX_DISTANCE_PC
    g_ids = _i(gcns["source_id"])[within]
    found = np.isin(g_ids, ours_ids)
    plx = _f(gcns["parallax"])[within]
    poe = plx / _f(gcns["parallax_error"])[within]
    gmag = _f(gcns["phot_g_mean_mag"])[within]
    miss = ~found
    rejected = np.isin(g_ids, _i(gaia_raw["source_id"])[~(_f(gaia_raw["fidelity_v2"]) > C.MIN_FIDELITY)])
    # White dwarfs: our HR-diagram cut against GCNS's probabilistic classification.
    wdp = np.asarray(ours["gcns_wd_prob"])
    has_p = np.isfinite(wdp)
    ours_wd = np.asarray(ours["cls"]) == "WD"
    theirs_wd = wdp > 0.5
    return {
        "ours_in_gcns_fraction": float(in_gcns.mean()),
        "ours_not_in_gcns": int((~in_gcns).sum()),
        "gcns_within_100pc": int(within.sum()),
        "gcns_within_100pc_found": int(found.sum()),
        "gcns_within_100pc_found_fraction": float(found.mean()),
        "gcns_missed_median_g": float(np.nanmedian(gmag[miss])),
        "gcns_missed_with_poe_below_10": int((miss & (poe <= C.MIN_PARALLAX_OVER_ERROR)).sum()),
        "gcns_missed_with_plx_below_10": int((miss & (plx <= C.MIN_PARALLAX_MAS)).sum()),
        "gcns_missed_rejected_by_fidelity": int((miss & rejected).sum()),
        "wd_both": int((has_p & ours_wd & theirs_wd).sum()),
        "wd_ours_only": int((has_p & ours_wd & ~theirs_wd).sum()),
        "wd_gcns_only": int((has_p & ~ours_wd & theirs_wd).sum()),
    }


MISSED_REASONS = (
    "too faint for Gaia",
    "binary: Gaia parallax spoiled by orbital motion",
    "binary: companion not resolved by Gaia",
    "single: no usable Gaia parallax",
    "other",
)


def _census_kind(objtype: str) -> str:
    t = objtype.rstrip("?")
    return {"LM": "star", "*": "star", "WD": "white dwarf", "BD": "brown dwarf"}.get(t, "other")


def compare_census_10pc(cat: Table, census_10: Table, limits: dict) -> dict:
    """Gaia against the best inventory of the 10 pc neighbourhood from all sources."""
    kind = np.array([_census_kind(str(t)) for t in _s(census_10["ObjType"])])
    keep = kind != "other"  # drop the exoplanets: this is about stars
    t = census_10[keep]
    kind = kind[keep]
    edr3 = _i(t["GaiaEDR3"])
    # HIP is stored as text, e.g. 'HIP 70890'.
    hipn = np.array([int("".join(ch for ch in h if ch.isdigit()) or -1) for h in _s(t["HIP"])])
    plx = _f(t["plx"])
    gmag = _f(t["Gmag"])
    gest = _f(t["Gest"])
    g_any = np.where(np.isfinite(gmag), gmag, gest)

    gaia_ids = set(int(s) for s in np.asarray(cat["source_id"])[np.asarray(cat["origin"]) == "Gaia DR3"])
    hip_ids = set(int(h) for h in np.asarray(cat["hip"])[np.asarray(cat["origin"]) == "Hipparcos"])
    in_gaia = np.array([int(s) in gaia_ids for s in edr3])
    in_hip = np.array([(int(h) in hip_ids) and not g for h, g in zip(hipn, in_gaia)])
    missed = ~in_gaia & ~in_hip

    # Why was each missing object missed? Faintness is checked first; after
    # that, what matters is whether it has companions. The census groups
    # objects into systems, so an object sharing its system number with others
    # is part of a binary or multiple.
    system = _i(t["Sys"])
    members = {sid: int((system == sid).sum()) for sid in set(system)}
    multiple = np.array([members[sid] > 1 for sid in system])
    too_faint = (g_any > limits["g_faint"]) | ((kind == "brown dwarf") & ~np.isfinite(g_any))
    reason = np.full(len(t), "", dtype=f"<U{max(map(len, MISSED_REASONS))}")
    reason[missed & too_faint] = MISSED_REASONS[0]
    reason[missed & (reason == "") & multiple & (edr3 > 0)] = MISSED_REASONS[1]
    reason[missed & (reason == "") & multiple] = MISSED_REASONS[2]
    reason[missed & (reason == "") & (edr3 > 0)] = MISSED_REASONS[3]
    reason[missed & (reason == "")] = MISSED_REASONS[4]

    by_kind = {}
    for k in ("star", "white dwarf", "brown dwarf"):
        sel = kind == k
        by_kind[k] = {
            "total": int(sel.sum()), "gaia": int((sel & in_gaia).sum()),
            "hipparcos": int((sel & in_hip).sum()), "missed": int((sel & missed).sum()),
            "missed_reasons": {r: int((sel & (reason == r)).sum()) for r in MISSED_REASONS},
        }
    ours_10 = int((np.asarray(cat["dist_pc"]) <= 10).sum())
    return {
        "n_objects": int(len(t)), "by_kind": by_kind,
        "our_catalogue_within_10pc": ours_10,
        "points": {
            "dist_pc": (1000.0 / plx).tolist(), "g": g_any.tolist(),
            "m_g": (g_any + 5 * np.log10(plx) - 10).tolist(),
            "kind": kind.tolist(), "status": np.where(in_gaia, "gaia", np.where(in_hip, "hipparcos", "missed")).tolist(),
            "reason": reason.tolist(),
            "name": [o or n for o, n in zip(_s(t["ObjName"]), _s(t["Name"]))],
        },
    }


def _simbad_letter(sptype: str) -> tuple[str, bool]:
    """SIMBAD spectral type -> (our class key, is_evolved)."""
    s = sptype.strip()
    if not s:
        return "", False
    if s[0] == "D":  # white dwarf: DA, DB, DC, DQ, DZ...
        return "WD", False
    for prefix in ("esd", "usd", "sd", "d"):
        if s.startswith(prefix) and len(s) > len(prefix) and s[len(prefix)] in "OBAFGKMLTY":
            s = s[len(prefix):]
            break
    if s[0] not in "OBAFGKMLTY":
        return "", False
    letter = s[0]
    evolved = any(tag in s for tag in ("III", "II", "IV")) and "V" not in s.replace("IV", "")
    if letter in "LTY":
        letter = "L/T"
    return letter, evolved


def compare_simbad_types(cat: Table) -> dict:
    """Our photometric class against published spectral types, within 20 pc."""
    sel = (np.asarray(cat["dist_pc"]) <= C.NAMED_RADIUS_PC) & (np.asarray(cat["simbad_sptype"]) != "")
    ours = np.asarray(cat["cls"])[sel]
    theirs = [_simbad_letter(s) for s in np.asarray(cat["simbad_sptype"])[sel]]
    letters = np.array([t[0] for t in theirs])
    evolved = np.array([t[1] for t in theirs])
    usable = letters != ""
    dwarfs = usable & np.isin(ours, ["B", "A", "F", "G", "K", "M"]) & ~evolved
    agree = dwarfs & (ours == letters)
    order = "OBAFGKM"
    near = dwarfs & np.array([o in order and t in order and abs(order.index(o) - order.index(t)) <= 1
                              for o, t in zip(ours, letters)])
    wd_theirs = usable & (letters == "WD")
    giant_theirs = usable & evolved
    return {
        "n_compared": int(usable.sum()),
        "dwarf_letter_agreement": float(agree.sum() / dwarfs.sum()),
        "dwarf_letter_within_one_class": float(near.sum() / dwarfs.sum()),
        "n_dwarfs_compared": int(dwarfs.sum()),
        "simbad_wd_we_call_wd": f"{int((wd_theirs & (ours == 'WD')).sum())}/{int(wd_theirs.sum())}",
        "simbad_evolved_we_call_giant": f"{int((giant_theirs & (ours == 'Giant')).sum())}/{int(giant_theirs.sum())}",
    }


# --------------------------------------------------------------------------
# Famous stars and planet hosts
# --------------------------------------------------------------------------

def locate_famous(cat: Table, famous: dict) -> dict:
    """Place every named star: from our catalogue if it is there, else from SIMBAD."""
    by_gaia = {int(s): i for i, s in enumerate(np.asarray(cat["source_id"]))
               if cat["origin"][i] == "Gaia DR3"}
    by_hip = {int(h): i for i, h in enumerate(np.asarray(cat["hip"])) if h > 0}

    def entry(star, inside: bool):
        s = famous[star.simbad]
        i = None
        if s["gaia_dr3"] is not None and s["gaia_dr3"] in by_gaia:
            i = by_gaia[s["gaia_dr3"]]
        elif s["hip"] is not None and s["hip"] in by_hip:
            i = by_hip[s["hip"]]
        rec = {"label": star.label, "note": star.note, "simbad_id": s["main_id"],
               "sp_type": s["sp_type"], "vmag": s["vmag"]}
        if i is not None and inside:
            row = cat[i]
            rec.update({"source": str(row["origin"]), "dist_pc": float(row["dist_pc"]),
                        "xyz_pc": [float(row["x"]), float(row["y"]), float(row["z"])],
                        "cls": str(row["cls"]), "teff": float(row["display_teff"]),
                        "m_g": float(row["m_g"]), "catalogue_index": int(i)})
            return rec
        if not s["plx_mas"]:
            rec.update({"source": "no parallax", "dist_pc": None})
            return rec
        d = 1000.0 / s["plx_mas"]
        x, y, z, l, b = sp.galactic_xyz(s["ra"], s["dec"], d)
        rec.update({"source": f"SIMBAD ({s['plx_bibcode']})", "dist_pc": float(d),
                    "xyz_pc": [float(np.ravel(x)[0]), float(np.ravel(y)[0]), float(np.ravel(z)[0])],
                    "cls": None, "teff": None, "catalogue_index": None})
        return rec

    return {"inside": [entry(s, True) for s in C.FAMOUS_STARS],
            "beyond": [entry(s, False) for s in C.BEYOND_STARS]}


def planet_hosts(cat: Table, planets: Table) -> list[dict]:
    """One record per planet-hosting star within the map, placed on our catalogue."""
    by_gaia = {int(s): i for i, s in enumerate(np.asarray(cat["source_id"]))
               if cat["origin"][i] == "Gaia DR3"}
    by_hip = {int(h): i for i, h in enumerate(np.asarray(cat["hip"])) if h > 0}
    hosts: dict[str, dict] = {}
    for row in planets:
        host = str(row["hostname"])
        h = hosts.setdefault(host, {"host": host, "planets": [], "methods": set(),
                                    "gaia": str(row["gaia_dr3_id"]), "hip": str(row["hip_name"]),
                                    "ra": float(row["ra"]), "dec": float(row["dec"]),
                                    "sy_dist": float(row["sy_dist"])})
        h["planets"].append(str(row["pl_name"]))
        h["methods"].add(str(row["discoverymethod"]))
    out = []
    for h in hosts.values():
        i = None
        if h["gaia"].startswith("Gaia DR3 "):
            i = by_gaia.get(int(h["gaia"].split()[-1]))
        if i is None and h["hip"].startswith("HIP "):
            digits = h["hip"].split()[1].rstrip("ABCab")
            if digits.isdigit():
                i = by_hip.get(int(digits))
        if i is not None:
            xyz = [float(cat["x"][i]), float(cat["y"][i]), float(cat["z"][i])]
            d = float(cat["dist_pc"][i])
        else:
            x, y, z, l, b = sp.galactic_xyz(h["ra"], h["dec"], h["sy_dist"])
            xyz = [float(np.ravel(x)[0]), float(np.ravel(y)[0]), float(np.ravel(z)[0])]
            d = h["sy_dist"]
        if d > C.MAX_DISTANCE_PC:
            continue
        out.append({"host": h["host"], "n_planets": len(h["planets"]),
                    "planets": sorted(h["planets"]), "methods": sorted(h["methods"]),
                    "dist_pc": d, "xyz_pc": xyz, "matched_to_catalogue": i is not None,
                    "catalogue_index": i})
    return sorted(out, key=lambda r: r["dist_pc"])


# --------------------------------------------------------------------------

def analyse(inputs: dict | None = None) -> tuple[dict, dict]:
    """Run the whole pipeline. Returns (results for JSON, artefacts for figures)."""
    if inputs is None:
        inputs = load_inputs()
    seq = inputs["sequence"]
    cat, build_diag = build_catalogue(inputs)
    print(f"  catalogue: {len(cat):,} stars "
          f"({build_diag['gaia']['n_kept']:,} Gaia + {build_diag['hipparcos']['added_from_hipparcos']:,} Hipparcos)")

    xyz = np.column_stack([cat["x"], cat["y"], cat["z"]])
    structure, neighbours = find_overdensities(xyz)
    cat["n_neighbours"] = neighbours
    identify_peaks(structure["peaks"], inputs["famous"])
    real_peaks = [p for p in structure["peaks"] if p["beats_null"]]
    cluster_mask = np.zeros(len(cat), dtype=bool)
    for p in real_peaks:
        cluster_mask |= np.linalg.norm(xyz - np.array(p["centre_xyz_pc"]), axis=1) < C.CLUSTER_MASK_RADIUS_PC
    cat["in_cluster"] = cluster_mask

    limits = gaia_limits(cat)
    results = {
        "_generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sample": {
            "max_distance_pc": C.MAX_DISTANCE_PC, "min_parallax_over_error": C.MIN_PARALLAX_OVER_ERROR,
            "min_fidelity": C.MIN_FIDELITY, "n_total": len(cat),
            "n_gaia": int((np.asarray(cat["origin"]) == "Gaia DR3").sum()),
            "n_hipparcos": int((np.asarray(cat["origin"]) == "Hipparcos").sum()),
        },
        "build": build_diag,
        "gaia_limits": limits,
        "census": census(cat),
        "naked_eye_sky": naked_eye_sky(inputs["naked_eye"]),
        "density_by_class": density_by_class(cat, cluster_mask),
        "completeness_distance_pc": completeness_distances(seq, limits["g_faint"]),
        "structure": structure,
        "north_south": north_south(cat, cluster_mask, real_peaks),
        "validation": {
            "gcns": compare_gcns(cat, inputs["gcns"], inputs["gaia"]),
            "census_10pc": compare_census_10pc(cat, inputs["census"], limits),
            "simbad_types": compare_simbad_types(cat),
        },
        "famous": locate_famous(cat, inputs["famous"]),
        "planet_hosts": planet_hosts(cat, inputs["exoplanets"]),
    }
    artefacts = {"catalogue": cat, "sequence": seq}
    return results, artefacts


def load_catalogue() -> Table:
    """Read the saved catalogue, with FITS's masked empty strings turned back into ''."""
    cat = Table.read(C.CATALOGUE_PATH)
    for name in cat.colnames:
        if cat[name].dtype.kind in "SU":
            cat[name] = _s(cat[name])
    return cat


#: Columns kept at full precision in the saved catalogue; everything else is
#: for display and is stored as 32-bit, which halves the file.
FULL_PRECISION = {"ra", "dec", "parallax", "parallax_error", "pmra", "pmdec"}


def save_results(results: dict, cat: Table | None = None) -> None:
    C.RESULTS_PATH.write_text(json.dumps(results, indent=2, default=float))
    if cat is not None:
        out = cat.copy(copy_data=False)
        for name in out.colnames:
            if out[name].dtype == np.float64 and name not in FULL_PRECISION:
                out[name] = out[name].astype(np.float32)
        out.write(C.CATALOGUE_PATH, overwrite=True)
