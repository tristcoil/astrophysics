"""Repeat the search for hypervelocity white dwarfs in Gaia DR3, and ask how
many candidates are real.

A star's speed across the sky is v_tan = 4.74 x proper motion / parallax
(km/s, with mas/yr and mas). Ordinary stars near the Sun move at 20-100 km/s;
halo stars reach ~500 km/s; nothing bound to the Galaxy exceeds its escape
speed, ~550-600 km/s here. So a cut at 600 km/s should leave only runaways --
except that a parallax measured too *small* inflates v_tan, and among 1.5
billion parallaxes some are badly wrong.

Four tests of how many candidates are parallax errors:

1. **The mirror test.** A negative parallax is unphysical: it can only be a
   measurement error. The number of sources with parallax < -k sigma that
   would *also* pass the speed cut measures how common wild errors are.
2. **The parallax-shift test.** How many sigma would each parallax have to
   be off for the star to be an ordinary (< 600 km/s) star -- with Gaia's
   quoted errors, and with the zero point corrected and realistic errors
   (Maiz Apellaniz et al. 2021). The robust survivors are then tested
   against the escape speed for *any* radial velocity.
3. **Follow-up spectra.** El-Badry et al. (2023) took spectra of 22 of the 25
   sources their query returned. We rerun their query and check, candidate by
   candidate, what the parallax alone says about each.
4. **Later data releases.** Shen et al. (2018) followed up seven DR2
   candidates; four were ordinary stars. Gaia DR3 measured all seven again
   with 34 instead of 22 months of data.
"""

from __future__ import annotations

import numpy as np
from astropy.table import Table, join

from config import (FIDELITY_MIN, GAIA_DIR, IPD_HARMONIC_MAX, IPD_MULTI_PEAK_MAX, K_KMS,
                    PLX_SNR_MIN, RUWE_MAX, VIS_PERIODS_MIN, V_TAN_MIN)
from kinematics import distance_grid, log_prior
from literature import CENSUS, ELBADRY_TABLE1, SHEN18_TABLE1


def _with_fidelity(t: Table) -> Table:
    path = GAIA_DIR / "fidelity.fits"
    if not path.exists():
        t["fidelity_v2"] = np.nan
        t["dist_bright_neighbour"] = np.nan
        return t
    f = Table.read(path)["source_id", "fidelity_v2", "dist_nearest_neighbor_at_least_m2_brighter"]
    f.rename_column("dist_nearest_neighbor_at_least_m2_brighter", "dist_bright_neighbour")
    out = join(t, f, keys="source_id", join_type="left")
    for c in ("fidelity_v2", "dist_bright_neighbour"):
        if hasattr(out[c], "filled"):
            out[c] = out[c].filled(np.nan)
    return out


def derive(t: Table) -> Table:
    """Face-value tangential velocity and its 1-sigma 'lower limit' proxy."""
    plx, err, pm = (np.asarray(t[c], float) for c in ("parallax", "parallax_error", "pm"))
    with np.errstate(divide="ignore", invalid="ignore"):
        t["v_tan_face"] = K_KMS * pm / np.abs(plx)
        t["v_tan_lower"] = K_KMS * pm / (np.abs(plx) + err)
        t["abs_M_G"] = np.asarray(t["phot_g_mean_mag"], float) + 5 * np.log10(np.abs(plx) / 100.0)
    t["plx_snr"] = plx / err
    return t


def load_search() -> Table:
    return derive(_with_fidelity(Table.read(GAIA_DIR / "search.fits")))


def quality_masks(t: Table) -> dict[str, np.ndarray]:
    ruwe = np.asarray(t["ruwe"], float) < RUWE_MAX
    fid = np.asarray(t["fidelity_v2"], float)
    fidelity = np.where(np.isfinite(fid), fid > FIDELITY_MIN, False)
    ipd = ((np.asarray(t["ipd_frac_multi_peak"]) <= IPD_MULTI_PEAK_MAX)
           & (np.asarray(t["ipd_gof_harmonic_amplitude"], float) < IPD_HARMONIC_MAX))
    vis = np.asarray(t["visibility_periods_used"]) >= VIS_PERIODS_MIN
    return {"no cut": np.ones(len(t), bool), f"RUWE < {RUWE_MAX}": ruwe,
            "clean images": ipd, f"fidelity > {FIDELITY_MIN}": fidelity,
            "all of these": ruwe & ipd & fidelity & vis}


def mirror_test(t: Table, thresholds=(3, 4, 5, 7, 10, 15, 20)) -> dict:
    """Counts of positive vs negative parallaxes above each significance, for
    each quality cut. The negative count estimates the false positives."""
    snr = np.asarray(t["plx_snr"], float)
    out = {}
    for name, m in quality_masks(t).items():
        rows = []
        for k in thresholds:
            npos = int(np.sum(m & (snr > k)))
            nneg = int(np.sum(m & (snr < -k)))
            rows.append(dict(threshold=k, n_pos=npos, n_neg=nneg,
                             false_fraction=(nneg / npos) if npos else np.nan))
        out[name] = rows
    return out


def prob_fast(plx, plx_err, pm, v_min=V_TAN_MIN, prior="flat_MG") -> float:
    """Posterior probability that v_tan > v_min, given the parallax alone
    (proper motions are precise enough to treat as exact)."""
    d = distance_grid()
    lnp = log_prior(d, prior) - 0.5 * ((plx - 1.0 / d) / plx_err) ** 2 + np.log(d)
    w = np.exp(lnp - lnp.max())
    v = K_KMS * pm * d
    return float(w[v > v_min].sum() / w.sum())


def elbadry_check() -> dict:
    """Rerun of El-Badry et al.'s query against their Table 1."""
    t = derive(_with_fidelity(Table.read(GAIA_DIR / "elbadry_query.fits")))
    ours = {int(i) for i in t["source_id"]}
    theirs = {i for i, _ in ELBADRY_TABLE1}
    verdict = dict(ELBADRY_TABLE1)
    rows = []
    for r in t:
        sid = int(r["source_id"])
        v = verdict.get(sid, "not in their table")
        cls = ("runaway" if ("D6" in v or "LP 40-365" in v) else
               "false positive" if ("DA" in v or "main-sequence" in v) else "unclassified")
        p600 = prob_fast(float(r["parallax"]), float(r["parallax_error"]), float(r["pm"]))
        p1000 = prob_fast(float(r["parallax"]), float(r["parallax_error"]), float(r["pm"]), 1000.0)
        # the parallax a 300 km/s star (a fast halo star) would need, and how
        # many sigma the measurement would have to be off
        plx300 = K_KMS * float(r["pm"]) / 300.0
        rows.append(dict(source_id=sid, verdict=v, cls=cls, G=float(r["phot_g_mean_mag"]),
                         bp_rp=float(r["bp_rp"]), plx=float(r["parallax"]),
                         plx_err=float(r["parallax_error"]), pm=float(r["pm"]),
                         v_tan_face=float(r["v_tan_face"]),
                         v_tan_lower=float(K_KMS * r["pm"] / (r["parallax"] + r["parallax_error"])),
                         p_fast600=p600, p_fast1000=p1000,
                         sigma_to_300=(plx300 - float(r["parallax"])) / float(r["parallax_error"]),
                         fidelity=float(r["fidelity_v2"]),
                         dist_bright_neighbour=float(r["dist_bright_neighbour"])))
    rows.sort(key=lambda r: -r["v_tan_lower"] if r["plx"] + r["plx_err"] > 0 else 1e9)
    return dict(n_returned=len(t), identical_to_table1=(ours == theirs),
                missing=sorted(theirs - ours), extra=sorted(ours - theirs), rows=rows)


def union_selection(search: Table) -> Table:
    """El-Badry et al.'s cuts *without* their parallax_over_error < 5 clause:
    their query plus the precise-parallax sources it set aside."""
    eb = derive(_with_fidelity(Table.read(GAIA_DIR / "elbadry_query.fits")))
    s = search
    precise = ((np.asarray(s["plx_snr"]) >= 5) & (np.asarray(s["phot_g_mean_mag"]) < 20)
               & (np.asarray(s["ruwe"]) < 1.4) & (np.asarray(s["pm"]) > 50)
               & (np.asarray(s["bp_rp"], float) < 0.5)
               & (np.asarray(s["v_tan_face"]) > 600) & (np.asarray(s["v_tan_lower"]) > 400))
    keep = [c for c in s.colnames if c in eb.colnames]
    from astropy.table import vstack
    out = vstack([eb[keep], s[precise][keep]], metadata_conflicts="silent")
    out["subset"] = ["plx/err < 5 (El-Badry query)"] * len(eb) + ["plx/err >= 5"] * int(precise.sum())
    return out


def census_flags(search: Table, union: Table) -> list[dict]:
    """Which known runaways each selection recovers."""
    s_ids = {int(i) for i in search["source_id"][np.asarray(search["plx_snr"]) > PLX_SNR_MIN]}
    u_ids = {int(i) for i in union["source_id"]}
    eb_ids = {i for i, _ in ELBADRY_TABLE1}
    return [dict(name=r.name, cls=r.cls, in_search=r.source_id in s_ids,
                 in_union=r.source_id in u_ids, in_elbadry=r.source_id in eb_ids)
            for r in CENSUS]


def shen_dr2_vs_dr3() -> list[dict]:
    """Shen et al.'s DR2 candidates re-measured by DR3 (best DR3 match per
    DR2 source: smallest angular distance)."""
    dr3 = Table.read(GAIA_DIR / "shen18_dr3.fits")
    dr2 = {int(r["source_id"]): r for r in Table.read(GAIA_DIR / "shen18_dr2.fits")}
    out = []
    for name, sid2, plx2_pub, err2_pub in SHEN18_TABLE1:
        m = dr3[np.asarray(dr3["dr2_source_id"]) == sid2]
        m = m[np.argsort(np.asarray(m["angular_distance"]))][0]
        r2 = dr2[sid2]
        pm2 = float(np.hypot(r2["pmra"], r2["pmdec"]))
        pm3 = float(np.hypot(m["pmra"], m["pmdec"]))
        out.append(dict(name=name, dr2_id=sid2, dr3_id=int(m["source_id"]),
                        plx_dr2=float(r2["parallax"]), err_dr2=float(r2["parallax_error"]),
                        plx_dr3=float(m["parallax"]), err_dr3=float(m["parallax_error"]),
                        pm_dr2=pm2, pm_dr3=pm3,
                        vtan_dr2=K_KMS * pm2 / float(r2["parallax"]),
                        vtan_dr3=K_KMS * pm3 / float(m["parallax"]),
                        vtan3sig_dr2=K_KMS * pm2 / (float(r2["parallax"]) + 3 * float(r2["parallax_error"])),
                        vtan3sig_dr3=K_KMS * pm3 / (float(m["parallax"]) + 3 * float(m["parallax_error"])),
                        shift_sigma=(float(m["parallax"]) - float(r2["parallax"])) / float(r2["parallax_error"]),
                        ruwe_dr3=float(m["ruwe"])))
    return out


# --------------------------------------------------------------------------
# 4. The parallax-shift test, with realistic parallax errors
# --------------------------------------------------------------------------
# A candidate is only interesting if its parallax would have to be wrong by
# many sigma for it to be an ordinary (< 600 km/s) star. Two corrections make
# "sigma" honest: Gaia DR3 parallaxes are offset by a few hundredths of a mas
# (Lindegren et al. 2021 zero point), and the catalogue errors are too small:
# Maiz Apellaniz et al. (2021, globular clusters) found external errors
# sigma_ext = sqrt((k sigma)^2 + sigma_s^2), sigma_s = 10.3 microarcsec, with
# k = 1.1-1.7 for G > 13, 1.7-2.0 for 11 < G < 13, > 2 for G < 9.2.

SIGMA_SYS_MAS = 0.0103
#: J1903-0023: a candidate found in Gaia DR3 that spectroscopy (Bhat et al.
#: 2025) showed to be a bound, eclipsing halo binary with a misleading parallax.
KNOWN_IMPOSTORS = {4265540383431508736: "J1903-0023 (bound eclipsing binary; Bhat et al. 2025)",
                   3841458366321558656: "halo hypervelocity star (Carballo-Bello et al. 2025)"}


def k_factor(G):
    """Representative mid-range values of Maiz Apellaniz et al.'s k(G)."""
    G = np.asarray(G, float)
    return np.select([G < 9.2, G < 11, G < 13], [2.2, 1.4, 1.85], 1.3)


def corrected_parallax(t: Table) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Zero-point-corrected parallax, external error, and the zero point."""
    import warnings
    from zero_point import zpt
    zpt.load_tables()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        zp = zpt.get_zpt(np.asarray(t["phot_g_mean_mag"], float),
                         np.asarray(t["nu_eff_used_in_astrometry"], float),
                         np.asarray(t["pseudocolour"], float), np.asarray(t["ecl_lat"], float),
                         np.asarray(t["astrometric_params_solved"], int), _warnings=False)
    zp = np.where(np.isfinite(zp), zp, -0.017)
    e = np.asarray(t["parallax_error"], float)
    e_ext = np.hypot(k_factor(t["phot_g_mean_mag"]) * e, SIGMA_SYS_MAS)
    return np.asarray(t["parallax"], float) - zp, e_ext, zp


def shift_needed(pm, plx, err, v=V_TAN_MIN):
    """How many sigma the true parallax must exceed the measured one for the
    star to move slower than v across the sky."""
    return (K_KMS * np.asarray(pm) / v - plx) / err


def min_galactocentric_speed(t: Table, plx) -> tuple[np.ndarray, np.ndarray]:
    """Slowest Galactocentric speed the star can have whatever its (unknown)
    radial velocity: the radial velocity can cancel the line-of-sight part of
    the Sun's motion, never the transverse part. Also returns positions."""
    from kinematics import to_galactocentric
    ra, dec = np.asarray(t["ra"], float), np.asarray(t["dec"], float)
    d = 1.0 / np.asarray(plx, float)
    zero = np.zeros(len(t))
    x, v0 = to_galactocentric(ra, dec, d, np.asarray(t["pmra"], float), np.asarray(t["pmdec"], float), zero)
    xs, _ = to_galactocentric(ra, dec, d * 0 + 1e-9, zero, zero, zero)
    rh = (x - xs) / np.linalg.norm(x - xs, axis=0)
    return np.linalg.norm(v0 - (v0 * rh).sum(axis=0) * rh, axis=0), x


def shift_test(search: Table, ks=(0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 15)) -> tuple[dict, Table]:
    """Counts of candidates by the parallax shift they need, with catalogue
    errors and with realistic ones; and the robust survivors, tested against
    the escape speed whatever their radial velocity."""
    import orbits
    pos = search[np.asarray(search["plx_snr"]) > PLX_SNR_MIN]
    clean = quality_masks(pos)["all of these"]
    pm = np.asarray(pos["pm"], float)
    z_cat = shift_needed(pm, np.asarray(pos["parallax"], float), np.asarray(pos["parallax_error"], float))
    plx_c, e_ext, zp = corrected_parallax(pos)
    z_ext = shift_needed(pm, plx_c, e_ext)
    counts = {
        "catalogue errors": {k: int(np.sum(z_cat > k)) for k in ks},
        "catalogue errors, clean astrometry": {k: int(np.sum((z_cat > k) & clean)) for k in ks},
        "zero point + external errors, clean": {k: int(np.sum((z_ext > k) & clean)) for k in ks},
    }
    robust = pos[(z_ext > 3) & clean]
    rz = z_ext[(z_ext > 3) & clean]
    rp = plx_c[(z_ext > 3) & clean]
    vmin, x = min_galactocentric_speed(robust, rp)
    vesc = orbits.escape_speed(x)
    robust["z_ext"] = rz
    robust["plx_corr"] = rp
    robust["v_gc_min"] = vmin
    robust["v_esc"] = vesc
    robust["abs_M_G_corr"] = np.asarray(robust["phot_g_mean_mag"], float) + 5 * np.log10(rp / 100.0)
    names = {r.source_id: r.name for r in CENSUS} | KNOWN_IMPOSTORS
    robust["known_as"] = [names.get(int(i), "") for i in robust["source_id"]]
    robust.sort("z_ext", reverse=True)
    unbound = robust[robust["v_gc_min"] > robust["v_esc"]]
    res = dict(n_candidates=len(pos), n_clean=int(clean.sum()), counts=counts,
               frac_within_1sigma=float(np.mean(z_cat <= 1)),
               zero_point_median=float(np.median(zp)),
               n_robust=len(robust), n_unbound_any_rv=len(unbound),
               unbound=[dict(source_id=int(r["source_id"]), known_as=r["known_as"],
                             G=float(r["phot_g_mean_mag"]), bp_rp=float(r["bp_rp"]),
                             M_G=float(r["abs_M_G_corr"]), plx=float(r["plx_corr"]),
                             z=float(r["z_ext"]), v_gc_min=float(r["v_gc_min"]),
                             v_esc=float(r["v_esc"])) for r in unbound])
    return res, robust
