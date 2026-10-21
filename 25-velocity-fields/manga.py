"""MaNGA: stars against gas, in ten thousand galaxies.

Two questions from the optical side of the velocity-field story:

1. **Do stars rotate slower than gas?** Gas clouds collide and lose their
   random motions, so gas settles onto nearly circular orbits. Stars never
   collide; a stellar disc keeps random velocities of tens of km/s, and that
   "pressure" supports it partly instead of rotation. So stars should lag the
   gas by an amount the Jeans equation predicts from their velocity
   dispersion: the *asymmetric drift* (Binney & Tremaine 2008, eq. 4.228).
   We measure both rotation curves in a few clean MaNGA spirals and check
   that the stars' dispersion accounts for the lag.

2. **How often does gas rotate against its stars?** Gas made inside a galaxy
   rotates with its stars. Gas that is captured -- from a merger, a passing
   companion or the cosmic web -- keeps the spin it arrived with. We fit the
   kinematic position angle of stars and of H-alpha gas in every MaNGA DR17
   galaxy (``kinpa.py``, Krajnovic et al.'s method with Zhou et al. 2022's
   choices) and count the misaligned (> 30 deg) and counter-rotating
   (> 150 deg) ones, by colour.
"""

from __future__ import annotations

import json
import pickle
import time
from multiprocessing import Pool

import numpy as np

import kinpa
import manga_io as M
from config import (COUNTER_ROTATING_DEG, DERIVED_DIR, INTRINSIC_AXIS_RATIO, MANGA_BINS_DIR,
                    MISALIGNED_DEG, PA_FIT_MAX_RE, PA_FIT_MIN_SNR, PA_MAX_ERROR_DEG, SIGMA_Z_OVER_R)

PA_RESULTS = DERIVED_DIR / "manga_pa.json"

#: Minimum number of Voronoi bins with good velocities to fit an axis.
MIN_BINS = 10
#: H-alpha bins must reach this S/N in flux for their velocity to be used.
HA_MIN_SNR = 3.0
#: Zhou et al. (2022): a galaxy is "lineless" if H-alpha S/N < 3 in > 90% of
#: its area within 1.5 Re.
LINELESS_FRACTION = 0.9
#: The "regular rotator" cut, our objective stand-in for visual inspection
#: (see ``census``).
REGULAR_MIN_BINS = 30
REGULAR_MAX_PA_ERR = 20.0
REGULAR_MIN_VAMP = 30.0
REGULAR_AMP_OVER_RESID = 3.0


# --------------------------------------------------------------------------
# Catalogues
# --------------------------------------------------------------------------

def catalogue() -> dict[str, dict]:
    """Per-galaxy properties from the DAP and DRP summary tables, by plate-IFU."""
    dap = M.read_dapall()
    drp = M.read_drpall()
    drp_index = {str(p.decode() if isinstance(p, bytes) else p).strip(): k
                 for k, p in enumerate(drp["plateifu"])}
    names = list(drp.dtype.names)
    out = {}
    ha = 23  # H-alpha channel in the DAPall emission-line arrays (C24 'Ha-6564')
    for row in dap:
        pifu = f"{int(row['PLATE'])}-{int(row['IFUDESIGN'])}"
        k = drp_index.get(pifu)
        rec = {"z": float(row["Z"]), "re_arcsec": float(row["NSA_ELPETRO_TH50_R"]),
               "ba": float(row["NSA_ELPETRO_BA"]), "phi": float(row["NSA_ELPETRO_PHI"]),
               "sersic_n": float(row["NSA_SERSIC_N"]),
               "ew_ha_1re": float(row["EMLINE_GEW_1RE"][ha]), "sfr_tot": float(row["SFR_TOT"]),
               "stellar_sigma_1re": float(row["STELLAR_SIGMA_1RE"]),
               "mngtarg1": int(row["MNGTARG1"]), "ra": float(row["OBJRA"]), "dec": float(row["OBJDEC"])}
        if k is not None:
            d = drp[k]
            absmag = np.asarray(d["nsa_elpetro_absmag"], float)
            rec["nuv_r"] = float(absmag[1] - absmag[4]) if absmag[1] != 0 and absmag[4] != 0 else float("nan")
            rec["abs_r"] = float(absmag[4])
            mass = float(d["nsa_elpetro_mass"])
            rec["log_mstar"] = float(np.log10(mass / 0.7 ** 2)) if mass > 0 else float("nan")  # h=1 -> h=0.7
            rec["ifu_size"] = int(str(d["ifudsgn"]).strip("b' ")[:-2] or 0) if "ifudsgn" in names else 0
        out[pifu] = rec
    return out


def colour_class(nuv_r: float) -> str:
    """Blue cloud / green valley / red sequence by NUV - r (Salim 2014's divisions)."""
    if not np.isfinite(nuv_r):
        return "unknown"
    return "blue" if nuv_r < 4 else ("green" if nuv_r < 5 else "red")


# --------------------------------------------------------------------------
# Position angles, one galaxy at a time
# --------------------------------------------------------------------------

def _component(b, which: str):
    """Bins usable for one component's position angle, and their velocities."""
    base = (b["r_re"] <= PA_FIT_MAX_RE) & (b["snr"] >= PA_FIT_MIN_SNR)
    if which == "star":
        sel = base & (b["star_vel_mask"] == 0) & (b["star_vel_ivar"] > 0)
        return sel, b["star_vel"], b["star_vel_ivar"]
    ha_snr = b["ha_flux"] * np.sqrt(np.clip(b["ha_flux_ivar"], 0, None))
    sel = base & (b["ha_vel_mask"] == 0) & (b["ha_vel_ivar"] > 0) & (ha_snr >= HA_MIN_SNR)
    return sel, b["ha_vel"], b["ha_vel_ivar"]


def pa_one(plateifu: str) -> dict:
    """Kinematic axes and receding sides of stars and gas for one galaxy."""
    try:
        b = M.load_bins(plateifu)
    except (OSError, ValueError) as err:
        return {"plateifu": plateifu, "status": f"unreadable: {err}"}
    inner = b["r_re"] <= PA_FIT_MAX_RE
    area = b["n_spaxels"].astype(float)
    ha_snr = b["ha_flux"] * np.sqrt(np.clip(b["ha_flux_ivar"], 0, None))
    weak = inner & ~((ha_snr >= 3) & (b["ha_flux_mask"] == 0))
    out = {"plateifu": plateifu, "n_bins": int(len(b["bin_id"])),
           "lineless": bool(area[weak].sum() > LINELESS_FRACTION * max(area[inner].sum(), 1))}
    for which in ("star", "gas"):
        sel, vel, ivar = _component(b, which)
        out[f"{which}_n"] = int(sel.sum())
        if sel.sum() < MIN_BINS:
            continue
        x, y = b["x_arcsec"][sel], b["y_arcsec"][sel]
        v = vel[sel] - np.median(vel[sel])
        err = 1 / np.sqrt(ivar[sel])
        try:
            axis, aerr, vsys, _ = kinpa.fit_kinematic_pa(x, y, v, err)
        except Exception as e:  # degenerate geometry (all bins on a line, ...)
            out[f"{which}_error"] = str(e)
            continue
        out[f"{which}_axis"] = axis
        out[f"{which}_pa_err"] = aerr
        out[f"{which}_pa"] = kinpa.receding_pa(x, y, v - vsys, axis)
        # how much of the field the rotation explains: amplitude of the
        # symmetric (rotating) part relative to the measurement errors
        from scipy.interpolate import LinearNDInterpolator
        from scipy.spatial import Delaunay
        interp = LinearNDInterpolator(Delaunay(np.column_stack([x, y])), v - vsys)
        vsym = kinpa.symmetrized(interp, x, y, v - vsys, axis)
        out[f"{which}_vamp"] = float(np.nanpercentile(np.abs(vsym), 90))
        out[f"{which}_resid_rms"] = float(np.sqrt(np.nanmean((v - vsys - vsym) ** 2)))
        out[f"{which}_median_err"] = float(np.median(err))
    if "star_pa" in out and "gas_pa" in out:
        out["dpa"] = kinpa.angle_between(out["star_pa"], out["gas_pa"])
    out["status"] = "ok"
    return out


def fit_all(*, rebuild=False, workers=6, quiet=False) -> list[dict]:
    """Position angles for every streamed galaxy (cached; parallel)."""
    if PA_RESULTS.exists() and not rebuild:
        return json.loads(PA_RESULTS.read_text())
    names = sorted(p.stem for p in MANGA_BINS_DIR.glob("*-*.npz") if ".part" not in p.name)
    t0 = time.time()
    with Pool(workers) as pool:
        results = pool.map(pa_one, names, chunksize=20)
    PA_RESULTS.write_text(json.dumps(results))
    if not quiet:
        print(f"    kinematic position angles for {len(results):,} galaxies in {time.time() - t0:.0f} s")
    return results


# --------------------------------------------------------------------------
# The census
# --------------------------------------------------------------------------

def census(results: list[dict], cat: dict[str, dict]) -> dict:
    """Counts and fractions of misaligned and counter-rotating gas.

    The funnel follows Zhou et al. (2022): galaxies with DAP output -> with
    emission lines (not "lineless") -> both axes measured with error <= 60
    deg ("measured"; their automatic cuts) -> misaligned / counter-rotating.
    Zhou et al. then inspected their candidates by eye and removed mergers,
    irregulars and multiple galaxies (37% of them). Our objective stand-in is
    the "regular" sample: both velocity fields dominated by ordered rotation
    (at least 30 bins, axis error <= 20 deg, rotation amplitude >= 30 km/s
    and at least three times the scatter left after removing it). A noisy or
    disordered field gives a near-random axis, which lands in the misaligned
    range most of the time; this cut removes such galaxies whichever way
    their axes point.
    """
    rows = []
    for r in results:
        c = cat.get(r["plateifu"], {})
        rows.append({**r, **{k: c.get(k) for k in ("nuv_r", "log_mstar", "ew_ha_1re", "sersic_n", "ba",
                                                     "z", "mngtarg1", "re_arcsec")}})
    n_all = len(rows)
    emission = [r for r in rows if r.get("status") == "ok" and not r.get("lineless")]
    measured = [r for r in emission if "dpa" in r and r["star_pa_err"] <= PA_MAX_ERROR_DEG
                and r["gas_pa_err"] <= PA_MAX_ERROR_DEG]

    def regular(r, comp):
        return (r[f"{comp}_n"] >= REGULAR_MIN_BINS and r[f"{comp}_pa_err"] <= REGULAR_MAX_PA_ERR
                and r[f"{comp}_vamp"] >= REGULAR_MIN_VAMP
                and r[f"{comp}_vamp"] >= REGULAR_AMP_OVER_RESID * r[f"{comp}_resid_rms"])

    clean = [r for r in measured if regular(r, "star") and regular(r, "gas")]
    out = {"n_streamed": n_all, "n_emission": len(emission), "n_measured": len(measured),
           "n_regular": len(clean), "by_sample": {}}
    for label, sample in (("measured", measured), ("regular", clean)):
        dpa = np.array([r["dpa"] for r in sample])
        stats = {"n": len(sample), "n_misaligned": int(np.sum(dpa > MISALIGNED_DEG)),
                 "n_counter": int(np.sum(dpa > COUNTER_ROTATING_DEG)),
                 "n_perpendicular": int(np.sum((dpa > 60) & (dpa < 120))),
                 "frac_misaligned": float(np.mean(dpa > MISALIGNED_DEG)),
                 "frac_counter": float(np.mean(dpa > COUNTER_ROTATING_DEG)),
                 "hist": np.histogram(dpa, bins=np.arange(0, 181, 10))[0].tolist(),
                 "by_colour": {}}
        for cls in ("blue", "green", "red"):
            sub = np.array([r["dpa"] for r in sample if colour_class(r.get("nuv_r") or np.nan) == cls])
            if len(sub):
                stats["by_colour"][cls] = {"n": int(len(sub)), "misaligned": int(np.sum(sub > MISALIGNED_DEG)),
                                           "counter": int(np.sum(sub > COUNTER_ROTATING_DEG)),
                                           "frac_misaligned": float(np.mean(sub > MISALIGNED_DEG)),
                                           "frac_counter": float(np.mean(sub > COUNTER_ROTATING_DEG))}
        # fraction against stellar mass
        lm = np.array([r.get("log_mstar") or np.nan for r in sample])
        edges = np.arange(8.5, 12.01, 0.5)
        stats["mass_bins"] = edges.tolist()
        stats["frac_misaligned_by_mass"] = [
            float(np.mean(dpa[(lm >= lo) & (lm < hi)] > MISALIGNED_DEG)) if np.sum((lm >= lo) & (lm < hi)) >= 20
            else None for lo, hi in zip(edges[:-1], edges[1:])]
        stats["n_by_mass"] = [int(np.sum((lm >= lo) & (lm < hi))) for lo, hi in zip(edges[:-1], edges[1:])]
        out["by_sample"][label] = stats
    def order(r):  # cleanest first: rotation amplitude over residual scatter, the worse of the two fields
        return -min(r["star_vamp"] / max(r["star_resid_rms"], 1e-3), r["gas_vamp"] / max(r["gas_resid_rms"], 1e-3))

    out["counter_rotators"] = [
        {k: r.get(k) for k in ("plateifu", "dpa", "star_pa", "gas_pa", "star_pa_err", "gas_pa_err", "nuv_r",
                               "log_mstar", "star_vamp", "gas_vamp", "star_resid_rms", "gas_resid_rms", "star_n",
                               "gas_n", "ew_ha_1re", "sersic_n")}
        for r in sorted((r for r in clean if r["dpa"] > COUNTER_ROTATING_DEG), key=order)]
    return out


# --------------------------------------------------------------------------
# Asymmetric drift in a few clean spirals
# --------------------------------------------------------------------------

def inclination_from_ba(ba: float, q0: float = INTRINSIC_AXIS_RATIO) -> float:
    """Hubble (1926): cos^2 i = ((b/a)^2 - q0^2) / (1 - q0^2) for an oblate disc of thickness q0."""
    c2 = (ba ** 2 - q0 ** 2) / (1 - q0 ** 2)
    return float(np.degrees(np.arccos(np.sqrt(np.clip(c2, 0, 1)))))


def ring_curves(b, pa_deg, incl_deg, r_edges_re, re_arcsec):
    """Rotation of stars and gas, and stellar dispersion, in rings (fixed geometry)."""
    pa, inc = np.radians(pa_deg), np.radians(incl_deg)
    x, y = b["x_arcsec"], b["y_arcsec"]
    # kinpa's convention: the receding major axis points along (-sin PA, cos PA)
    along = -x * np.sin(pa) + y * np.cos(pa)
    across = (-x * np.cos(pa) - y * np.sin(pa)) / np.cos(inc)
    r = np.hypot(along, across) / re_arcsec
    cos_t = np.where(r > 0, along / np.where(np.hypot(along, across) > 0, np.hypot(along, across), 1), 0)
    ha_snr = b["ha_flux"] * np.sqrt(np.clip(b["ha_flux_ivar"], 0, None))
    good_g = (b["ha_vel_mask"] == 0) & (b["ha_vel_ivar"] > 0) & (ha_snr >= 5)
    good_s = (b["star_vel_mask"] == 0) & (b["star_vel_ivar"] > 0) & (b["snr"] >= 10)
    good_d = good_s & (b["star_sigma_mask"] == 0) & (b["star_sigma_ivar"] > 0)
    vs_g = np.median(b["ha_vel"][good_g & (r < 0.5)]) if np.any(good_g & (r < 0.5)) else np.median(b["ha_vel"][good_g])
    vs_s = np.median(b["star_vel"][good_s & (r < 0.5)]) if np.any(good_s & (r < 0.5)) else np.median(b["star_vel"][good_s])
    sig2 = b["star_sigma"] ** 2 - b["star_sigmacorr0"] ** 2  # DAP: quadrature correction for resolution
    rows = []
    for lo, hi in zip(r_edges_re[:-1], r_edges_re[1:]):
        ring = (r >= lo) & (r < hi) & (np.abs(cos_t) > 0.3)  # keep away from the minor axis
        row = {"r_re": 0.5 * (lo + hi)}
        for name, good, vel, vsys in (("gas", good_g, b["ha_vel"], vs_g), ("star", good_s, b["star_vel"], vs_s)):
            s = ring & good
            if s.sum() >= 6:
                proj = np.sin(inc) * cos_t[s]
                row[f"v_{name}"] = float(np.sum((vel[s] - vsys) * proj) / np.sum(proj ** 2))
                row[f"n_{name}"] = int(s.sum())
        s = (r >= lo) & (r < hi) & good_d & (sig2 > 0)
        if s.sum() >= 6:
            row["sigma_los"] = float(np.sqrt(np.median(sig2[s])))
        rows.append(row)
    return rows


def asymmetric_drift(rows, incl_deg, re_arcsec, *, sigma_z_over_r=SIGMA_Z_OVER_R):
    """Correct stellar rotation for asymmetric drift (Binney & Tremaine eq. 4.228).

    Assumes an exponential disc (scale length h = Re / 1.678) whose radial
    dispersion falls as exp(-R / 2h) (constant scale height), and
    sigma_phi^2 = sigma_R^2 / 2 (the epicyclic ratio for a flat rotation
    curve). Then V_c^2 = V_star^2 + sigma_R^2 (2 R / h - 1/2). The radial
    dispersion follows from the line-of-sight one averaged around the ring:
    <sigma_los^2> = sigma_z^2 cos^2 i + (sigma_R^2 + sigma_phi^2) sin^2 i / 2.
    """
    inc = np.radians(incl_deg)
    k = sigma_z_over_r ** 2 * np.cos(inc) ** 2 + 0.75 * np.sin(inc) ** 2
    for row in rows:
        r_over_h = row["r_re"] * 1.678
        bracket = 2 * r_over_h - 0.5
        if "v_star" in row and "sigma_los" in row:
            sigma_r = row["sigma_los"] / np.sqrt(k)
            row["sigma_r"] = float(sigma_r)
            row["v_star_corrected"] = float(np.sqrt(max(row["v_star"] ** 2 + sigma_r ** 2 * bracket, 0)))
        if "v_star" in row and "v_gas" in row and bracket > 0:
            # run the Jeans equation backwards: the radial dispersion the observed lag implies
            d2 = row["v_gas"] ** 2 - row["v_star"] ** 2
            row["sigma_r_implied"] = float(np.sqrt(d2 / bracket)) if d2 > 0 else 0.0
    return rows


#: MNGTARG1 bit of MaNGA's secondary sample, observed out to 2.5 effective radii.
SECONDARY_BIT = 1 << 11


def pick_drift_galaxies(results, cat, n=3) -> list[str]:
    """Clean, inclined, aligned, gas-rich discs, observed far out.

    The asymmetric-drift correction assumes an exponential disc of stars, so
    the showcase avoids big bulges: disc-dominated light (Sersic n < 2),
    moderate mass (log M* 9.8-10.7), low stellar dispersion within 1 Re
    (< 80 km/s), inclination ~55-70 deg, stars and gas aligned within 10 deg,
    both axes good to 5 deg, from the secondary sample (coverage to 2.5 Re)
    with bins out to at least 1.8 Re.
    """
    cands = []
    for r in results:
        c = cat.get(r["plateifu"], {})
        if r.get("status") != "ok" or "dpa" not in r or r["dpa"] > 10:
            continue
        if not (0.35 <= (c.get("ba") or 0) <= 0.6 and (c.get("sersic_n") or 9) < 2.0):
            continue
        if not (9.8 <= (c.get("log_mstar") or 0) <= 10.7) or not (0 < (c.get("stellar_sigma_1re") or 999) < 80):
            continue
        if not (c.get("mngtarg1", 0) & SECONDARY_BIT):
            continue
        if r["star_n"] < 120 or r["gas_n"] < 120 or r["star_pa_err"] > 5 or r["gas_pa_err"] > 5:
            continue
        cands.append((min(r["star_n"], r["gas_n"]), r["plateifu"]))
    picks = []
    for _, pifu in sorted(cands, reverse=True):
        b = M.load_bins(pifu)
        good = (b["star_vel_mask"] == 0) & (b["snr"] >= 10)
        if np.sum(good & (b["r_re"] >= 1.5)) >= 20 and b["r_re"][good].max() >= 1.8:
            picks.append(pifu)
        if len(picks) == n:
            break
    return picks


def drift_showcase(results, cat, *, n=3) -> list[dict]:
    out = []
    for pifu in pick_drift_galaxies(results, cat, n):
        r = next(x for x in results if x["plateifu"] == pifu)
        c = cat[pifu]
        b = M.load_bins(pifu)
        incl = inclination_from_ba(c["ba"])
        edges = np.arange(0.0, 2.01, 0.25)
        rows = asymmetric_drift(ring_curves(b, r["gas_pa"], incl, edges, c["re_arcsec"]), incl, c["re_arcsec"])
        out.append({"plateifu": pifu, "incl": incl, "re_arcsec": c["re_arcsec"], "log_mstar": c.get("log_mstar"),
                    "nuv_r": c.get("nuv_r"), "gas_pa": r["gas_pa"], "star_pa": r["star_pa"], "rows": rows})
    return out


def run(*, rebuild=False, quiet=False) -> dict:
    path = DERIVED_DIR / "manga_summary.pkl"
    if path.exists() and not rebuild:
        with open(path, "rb") as fh:
            return pickle.load(fh)
    cat = catalogue()
    results = fit_all(rebuild=rebuild, quiet=quiet)
    summary = {"census": census(results, cat), "drift": drift_showcase(results, cat)}
    with open(path, "wb") as fh:
        pickle.dump(summary, fh)
    return summary
