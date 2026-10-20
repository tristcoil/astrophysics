"""Every measurement in the project, from cached data to data/results.json.

Steps (each cached under data/cache/ so a re-run is quick):

1. The disc pipeline (discs.py) for the 18 single DSHARP discs and HL Tau.
2. Comparisons with DSHARP II: orientations, ring/gap radii, widths and depths, R_dust;
   with DSHARP I: flux, peak, noise; with ALMA 2015: HL Tau's orientation.
3. Planet masses from our gap widths with DSHARP VII's relation, and their check
   against DSHARP VII's own table; the planets among the known exoplanets.
4. The snow-line test: midplane temperature at every gap we find.
5. Spiral pitch angles (tracing and the log-spiral Fourier transform) vs DSHARP III.
6. HD 163296's CO kinematics: stellar mass and emission surface vs Teague et al. 2019,
   the rotation curve, and the velocity residuals at Pinte et al.'s two kinks.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
import zlib

import numpy as np
from scipy import ndimage

import alma_images as A
import discs as D
import geometry as G
import kinematics as K
import literature as L
import planets as P
import profiles as PR
import spirals as S
from config import (ALMA2015_HLTAU_GEOMETRY, CACHE_DIR, CO_CUBE_PATH, EXOPLANET_PATH, HUANG_SPIRAL_RANGES,
                    HUANG_SPIRALS_LOG, KURTOVIC_GEOMETRY, M_EARTH_CGS, PINTE2020_KINKS, RESULTS_PATH, SINGLE_DISCS,
                    SNOWLINES_K, TARGETS, TEAGUE2019_HD163296, ZHANG_ALPHAS, ensure_dirs)
from distances import gaia_distances

DISC_KEYS = [t.key for t in SINGLE_DISCS] + ["HLTau"]


def seed(name: str) -> int:
    """Deterministic per-task seed (Python's str hash is randomised per process)."""
    return zlib.crc32(name.encode())


# ---------------------------------------------------------------------------
# 1. Discs
# ---------------------------------------------------------------------------
def _work(key: str):
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    from threadpoolctl import threadpool_limits

    with threadpool_limits(1):
        t0 = time.time()
        D.analyse(key, np.random.default_rng(seed(key)), verbose=True)
        return key, time.time() - t0


def run_discs(refit: bool = False, processes: int = 4) -> dict:
    todo = [k for k in DISC_KEYS if refit or D.load(k) is None]
    if todo:
        print(f"disc pipeline for {len(todo)} discs ({processes} processes)")
        with mp.get_context("fork").Pool(processes) as pool:
            for key, dt in pool.imap_unordered(_work, todo):
                pass
    return {k: D.load(k) for k in DISC_KEYS}


# ---------------------------------------------------------------------------
# 2. Comparisons with DSHARP I/II
# ---------------------------------------------------------------------------
def compare_geometry(res: dict) -> dict:
    geo = L.dsharp2_geometry()
    rows = {}
    for key, o in res.items():
        g = o["geometry"]
        if key == "HLTau":
            pub = dict(inc=ALMA2015_HLTAU_GEOMETRY["inc_deg"], pa=ALMA2015_HLTAU_GEOMETRY["pa_deg"], method="E",
                       source="ALMA Partnership 2015")
            pub_x = pub_y = (None, None)
            rd_pub = None
        else:
            h = geo[key]
            pub = dict(inc=(h.inc.value, h.inc.err), pa=(h.pa.value, h.pa.err), method=h.method, source="DSHARP II")
            pub_x, pub_y = (h.dx_mas.value, h.dx_mas.err), (h.dy_mas.value, h.dy_mas.err)
            rd_pub = h.r_dust_au.value
        d_inc = g["inc"] - pub["inc"][0]
        d_pa = float(G.angle_diff_180(g["pa"], pub["pa"][0]))
        s_inc = np.hypot(g["inc_err"] or 0, pub["inc"][1] or 0)
        s_pa = np.hypot(g["pa_err"] or 0, pub["pa"][1] or 0)
        rows[key] = dict(
            name=o["name"], ours=dict(method=g["method"], inc=g["inc"], inc_err=g["inc_err"], pa=g["pa"] % 180.0,
                                      pa_err=g["pa_err"], x0_mas=g["x0"] * 1e3, y0_mas=g["y0"] * 1e3,
                                      n_ellipses=g.get("n_features", 0)),
            gauss=dict(inc=o["geometry_gauss"]["inc"], inc_err=o["geometry_gauss"]["inc_err"],
                       pa=o["geometry_gauss"]["pa"], pa_err=o["geometry_gauss"]["pa_err"]),
            published=dict(method=pub["method"], inc=pub["inc"][0], inc_err=pub["inc"][1], pa=pub["pa"][0],
                           pa_err=pub["pa"][1], x0_mas=pub_x[0], y0_mas=pub_y[0], source=pub["source"]),
            d_inc=d_inc, d_inc_sigma=d_inc / s_inc if s_inc > 0 else None,
            d_pa=d_pa, d_pa_sigma=d_pa / s_pa if s_pa > 0 else None,
            r_dust_au_dsharp=o["r_dust_au_dsharp"], r_dust_pub=rd_pub, r_dust_au=o["r_dust_au"],
            d_x_mas=(g["x0"] * 1e3 - pub_x[0]) if pub_x[0] is not None else None,
            d_y_mas=(g["y0"] * 1e3 - pub_y[0]) if pub_y[0] is not None else None,
        )
    return rows


def compare_images(res: dict) -> dict:
    tab = L.dsharp1_images()
    out = {}
    for key, o in res.items():
        if key not in tab:
            continue
        t = tab[key]
        out[key] = dict(flux_mjy=o["flux"] * 1e3, flux_pub=t.flux_mjy, peak_mjy=o["peak"] * 1e3, peak_pub=t.peak_mjy,
                        rms_ujy=o["noise"] * 1e6, rms_pub=t.rms_ujy,
                        beam_mas=(o["beam"][0] * 1e3, o["beam"][1] * 1e3), beam_pub=t.beam_mas)
    return out


def compare_features(res: dict) -> dict:
    """Ring/gap radii vs DSHARP II, by method of the published measurement."""
    rows = []
    for key, o in res.items():
        beam = np.sqrt(o["beam"][0] * o["beam"][1])
        for m in o["matches"]:
            row = dict(disc=key, label=m["label"], kind=m["kind"], method=m["method"], approx=m["approx"],
                       r_pub_mas=m["r_pub"] * 1e3, r_pub_err_mas=(m["r_pub_err"] or np.nan) * 1e3,
                       found="ours" in m, beam_mas=beam * 1e3)
            if "ours" in m:
                row.update(ours=m["ours"], r_ours_mas=m["r_ours"] * 1e3, r_ours_err_mas=m["r_ours_err"] * 1e3,
                           sigma=m["sigma"], our_method=m["our_method"],
                           dr_mas=(m["r_ours"] - m["r_pub"]) * 1e3)
            rows.append(row)
    found = [r for r in rows if r["found"]]
    by_method = {}
    for meth in ("E", "R", "V", "-"):
        sel = [r for r in rows if r["method"] == meth]
        if not sel:
            continue
        f = [r for r in sel if r["found"]]
        dr = np.array([r["dr_mas"] for r in f])
        by_method[meth] = dict(n=len(sel), found=len(f),
                               median_abs_dr_mas=float(np.median(np.abs(dr))) if f else None,
                               median_abs_dr_beam=float(np.median([abs(r["dr_mas"]) / r["beam_mas"] for r in f]))
                               if f else None)
    extras = []
    for key, o in res.items():
        matched = {m.get("ours") for m in o["matches"]}
        for f in o["features"]:
            if f["label"] not in matched:
                extras.append(dict(disc=key, label=f["label"], sigma=f["sigma"], secure=f["secure"]))
    return dict(rows=rows, by_method=by_method, n_published=len(rows), n_found=len(found), extras=extras)


def compare_widths(res: dict) -> list[dict]:
    """Our gap/ring widths and depths (DSHARP II definitions) vs DSHARP II's."""
    out = []
    feats = {f.name: f for f in L.dsharp2_features()}
    for key, o in res.items():
        d_ds = o["d_dsharp"]
        ours_to_pub = {m["ours"]: m["label"] for m in o["matches"] if "ours" in m}
        for p in o["pairs"]:
            gap_pub = ours_to_pub.get(p["gap"])
            ring_pub = ours_to_pub.get(p["ring"])
            e = p.get("errors", {})
            row = dict(disc=key, gap=p["gap"], ring=p["ring"], gap_pub=gap_pub, ring_pub=ring_pub, case=p["case"],
                       gap_width_au=p["gap_width"] * d_ds if p["gap_width"] is not None else None,
                       gap_width_err=(e.get("gap_width") or np.nan) * d_ds,
                       ring_width_au=p["ring_width"] * d_ds if p["ring_width"] is not None else None,
                       ring_width_err=(e.get("ring_width") or np.nan) * d_ds,
                       depth=p["depth"], depth_err=e.get("depth"), delta=p["delta"], delta_err=e.get("delta"))
            if gap_pub and f"{key}:{gap_pub}" in feats:
                fp = feats[f"{key}:{gap_pub}"]
                row.update(gap_width_pub=fp.width_au.value, gap_width_pub_err=fp.width_au.err,
                           gap_width_pub_upper=fp.width_au.upper_limit, depth_pub=fp.depth.value,
                           depth_pub_err=fp.depth.err)
            if ring_pub and f"{key}:{ring_pub}" in feats:
                fr = feats[f"{key}:{ring_pub}"]
                row.update(ring_width_pub=fr.width_au.value, ring_width_pub_err=fr.width_au.err)
            out.append(row)
    return out


# ---------------------------------------------------------------------------
# 3. Planets
# ---------------------------------------------------------------------------
def _scaled_star(key: str, d_new: float) -> tuple[float, float]:
    """(L*, M*) from DSHARP I, with L* rescaled from the DR2 to our distance (L ~ d^2)."""
    st = L.dsharp1_stars()[key]
    return st.lstar * (d_new / st.distance.value) ** 2, st.mstar


def planet_masses(res: dict) -> dict:
    """Our widths -> planet masses, for DSHARP VII's 12 fitted gaps.

    DSHARP VII's AS 209 "99 au" gap is the pair D90 + D105 taken as one gap: they read the
    thin ring B97 between them as material left in the planet's horseshoe orbit (their
    sec. 3.3.2), so the gap runs from inside D90 to the outer edge of D105, its depth and
    I_edge set by D105 and B120. We measure it the same way by dropping B97 and D90 from
    the feature list. Stokes numbers (how strongly the grains feel the gas) come from
    DSHARP VII's own table: they read them off simulations (their Fig. 13) that are not
    published as numbers."""
    stars = L.dsharp1_stars()
    rows = []
    for zp in [p for p in L.dsharp7_planets() if p.from_fit]:
        o = res[zp.key]
        d_ds, d_new = o["d_dsharp"], o["d_gaia"]
        prof = D.load_profile(zp.key)
        feats = [(f["label"], f["kind"], f["r"]) for f in o["features"]]
        merged = zp.key == "AS209" and abs(zp.r_gap_au - 99) < 2
        if merged:
            feats = [f for f in feats if not (f[1] in "BD" and 85 <= f[2] * d_ds <= 100)]
        gaps = [f for f in feats if f[1] == "D"]
        target_r = (105.5 if merged else zp.r_gap_au) / d_ds
        g = min(gaps, key=lambda f: abs(f[2] - target_r)) if gaps else None
        if g is None or abs(g[2] - target_r) * d_ds > 6.0:
            rows.append(dict(disc=zp.key, r_gap_pub=zp.r_gap_au, found=False))
            continue
        pairs, errs = PR.pair_measures_with_errors(prof, feats, 0.1 / d_ds, np.sqrt(o["beam"][0] * o["beam"][1]),
                                                   np.random.default_rng(seed(zp.key + str(zp.r_gap_au))))
        pr = next((p for p in pairs if p.gap == g[0]), None)
        if pr is None or not np.isfinite(pr.delta):
            rows.append(dict(disc=zp.key, r_gap_pub=zp.r_gap_au, found=False))
            continue
        delta, delta_err = float(pr.delta), float(errs[(pr.gap, pr.ring)]["delta"])
        r_gap_ds = zp.r_gap_au if merged else g[2] * d_ds  # DSHARP VII puts AS 209's planet at 99 au
        # (a) DSHARP VII's own inputs (DSHARP I star, DR2 distance) with our width
        st = stars[zp.key]
        hr_ds = float(P.aspect_ratio(r_gap_ds, st.lstar, st.mstar))
        # (b) the same at our Gaia DR3 distance
        lstar, mstar = _scaled_star(zp.key, d_new)
        r_gap_new = r_gap_ds * d_new / d_ds
        hr_new = float(P.aspect_ratio(r_gap_new, lstar, mstar))
        masses, masses_new = {}, {}
        for model, st_units in zip(("DSD1", "1mm", "DSD2"), zp.st_used):
            if st_units is None:
                continue
            for alpha in ZHANG_ALPHAS:
                if (model, alpha) not in zp.masses:
                    continue  # DSHARP VII found no solution for this combination
                m, lo, hi = P.planet_mass_mjup(delta, hr_ds, alpha, st_units, model, st.mstar)
                masses[f"{model}_{alpha:g}"] = dict(m=float(m), lo=float(lo), hi=float(hi), pub=zp.masses[(model, alpha)])
                m2, lo2, hi2 = P.planet_mass_mjup(delta, hr_new, alpha, st_units, model, mstar)
                masses_new[f"{model}_{alpha:g}"] = dict(m=float(m2), lo=float(lo2), hi=float(hi2))
        # the width's own error, propagated to the preferred mass (q ~ Delta^(1/B))
        st_units = zp.st_used[0]
        _, a_fit, b_fit, *_ = P.width_fit(st_units, "DSD1")
        frac = (1.0 / b_fit) * delta_err / delta
        pref = masses.get("DSD1_0.001")
        rows.append(dict(disc=zp.key, name=o["name"], r_gap_pub=zp.r_gap_au, our_gap=g[0], merged=merged,
                         found=True, r_gap_au_dsharp=r_gap_ds, r_gap_au=r_gap_new, delta=delta, delta_err=delta_err,
                         delta_pub=zp.width, h_over_r=hr_ds, h_over_r_pub=zp.h_over_r, h_over_r_dr3=hr_new,
                         mstar=mstar, st_used=zp.st_used, masses=masses, masses_dr3=masses_new,
                         preferred=dict(m=pref["m"], width_frac_err=frac, pub=pref["pub"]) if pref else None))
    found = [r for r in rows if r["found"] and r["preferred"]]
    ratio = np.array([r["preferred"]["m"] / r["preferred"]["pub"] for r in found])
    return dict(rows=rows, n=len(rows), n_found=len(found),
                median_ratio=float(np.median(ratio)) if found else None,
                log_rms=float(np.sqrt(np.mean(np.log10(ratio) ** 2))) if found else None)


def known_exoplanets() -> dict:
    """Planets with a measured mass (or m sin i) and semi-major axis, as in A2 (masses filled
    in from the archive's mass-radius relation are left out)."""
    import pandas as pd

    df = pd.read_csv(EXOPLANET_PATH)
    calc = df["pl_bmasse_reflink"].fillna("").str.contains("CALCULATED_VALUE", regex=False)
    ok = (df["pl_bmassj"] > 0) & (df["pl_orbsmax"] > 0) & (df["pl_bmassprov"] != "M-R relationship") & ~calc
    ok &= df["pl_bmasselim"].fillna(0).eq(0) & df["pl_orbsmaxlim"].fillna(0).eq(0)
    sub = df[ok]
    groups = {}
    for meth, g in sub.groupby("discoverymethod"):
        groups[meth] = dict(a=g["pl_orbsmax"].tolist(), m=g["pl_bmassj"].tolist(), names=g["pl_name"].tolist(),
                            age_myr=(g["st_age"] * 1e3).tolist())
    young = sub[(sub["st_age"] < 0.03)]  # younger than 30 Myr
    return dict(n=int(ok.sum()), n_total=int(len(df)), by_method=groups,
                young=dict(names=young["pl_name"].tolist(), a=young["pl_orbsmax"].tolist(),
                           m=young["pl_bmassj"].tolist(), method=young["discoverymethod"].tolist()))


# ---------------------------------------------------------------------------
# 4. Snow lines
# ---------------------------------------------------------------------------
def snowline_test(res: dict) -> dict:
    """Midplane temperature (DSHARP II eq. 1, phi = 0.02) at every secure gap and ring we
    find. If gaps were carved by ice lines, they would cluster at the freeze-out
    temperatures of CO, N2 or CO2, whatever the star's luminosity."""
    rows = []
    for key, o in res.items():
        if key == "HLTau":
            continue
        lstar, _ = _scaled_star(key, o["d_gaia"])
        for f in o["features"]:
            if not f["secure"]:
                continue
            t = float(P.midplane_temperature(f["r_au"], lstar))
            rows.append(dict(disc=key, label=f["label"], kind=f["kind"], r_au=f["r_au"], t_mid=t))
    gaps = np.array([r["t_mid"] for r in rows if r["kind"] == "D"])
    frac = {}
    for name, (lo, hi) in SNOWLINES_K.items():
        frac[name] = float(np.mean((gaps >= lo) & (gaps <= hi))) if gaps.size else None
    # how often would gaps fall into those ranges if their temperatures were spread evenly in log T
    # over the observed range? (a crude null hypothesis)
    lo_t, hi_t = np.log(gaps.min()), np.log(gaps.max())
    null = {name: float(max(0.0, min(np.log(hi), hi_t) - max(np.log(lo), lo_t)) / (hi_t - lo_t))
            for name, (lo, hi) in SNOWLINES_K.items()}
    return dict(rows=rows, n_gaps=int(gaps.size), fraction_in_range=frac, uniform_log_t_expectation=null,
                t_range=(float(gaps.min()), float(gaps.max())))


# ---------------------------------------------------------------------------
# 5. Spirals
# ---------------------------------------------------------------------------
def _m2_phase(pol: S.PolarImage, r_lo: float, r_hi: float, p: float, m: int = 2) -> float:
    u = np.linspace(np.log(r_lo), np.log(r_hi), 400)
    th = np.radians(pol.theta)
    res = np.nan_to_num(pol.residual)
    ru = np.array([np.interp(np.exp(u), pol.r, res[:, j]) for j in range(th.size)]).T
    w = np.exp(2 * u)[:, None] * np.clip(ru, 0, None)
    return float(-np.angle(np.sum(w * np.exp(-1j * (m * th[None, :] + p * u[:, None])))))


def spiral_analysis(res: dict) -> dict:
    out = {}
    for key in ("Elias27", "IMLup", "WaOph6"):
        o = res[key]
        g = o["geometry"]
        d = o["d_dsharp"]
        geom = G.DiscGeometry(g["x0"], g["y0"], g["inc"], g["pa"])
        img = A.load_continuum(key).cutout(o["half_width"], *o["centre_guess"])
        ranges = HUANG_SPIRAL_RANGES[key]
        lo, hi = ranges[0][0] / d, ranges[-1][1] / d
        pol = S.polar_image(img, geom, hi * 1.2, img.beam_mean / 4, 1.0)
        p_grid, amp = S.log_spiral_spectrum(pol, lo, hi, m=2)
        p_max, b, mu = S.pitch_from_spectrum(p_grid, amp)
        # half-maximum width of the Fourier peak -> a rough uncertainty on the pitch
        k = int(np.argmin(np.abs(p_grid - p_max)))
        half = amp[k] / 2
        left = k
        while left > 0 and amp[left] > half:
            left -= 1
        right = k
        while right < amp.size - 1 and amp[right] > half:
            right += 1
        mu_lo = np.degrees(np.arctan(2 / max(abs(p_grid[right]), abs(p_grid[left]))))
        mu_hi = np.degrees(np.arctan(2 / min(abs(p_grid[right]), abs(p_grid[left]))))
        phase = _m2_phase(pol, lo, hi, p_max)
        u_mid = 0.5 * (np.log(lo) + np.log(hi))
        # trace: Elias 27's outer arms are measured on a smoothed copy (~80 mas, as DSHARP III's
        # tapered "mid-res" image)
        pts = []
        for (r1, r2, dth) in ranges:
            im_use = img
            if key == "Elias27" and r1 >= 110:
                extra = np.sqrt(max(0.080**2 - img.beam_mean**2, 0.0)) / 2.3548 / img.pix
                # convolving Jy/beam with a unit-sum kernel: rescale by the ratio of beam areas
                im_use = A.Image(img.key, ndimage.gaussian_filter(img.data, extra) * (0.080**2 / img.beam_mean**2),
                                 img.x, img.y, img.pix, 0.080, 0.080, 0.0, img.freq_hz, img.ra0_deg, img.dec0_deg,
                                 img.header)
            pol_use = S.polar_image(im_use, geom, hi * 1.2, img.beam_mean / 4, 1.0)
            pts += S.ray_maxima(pol_use, r1 / d, r2 / d, dth, 3 * o["noise"])
        arms = {}
        for arm in (0, 1):
            th0 = np.degrees((phase + 2 * np.pi * arm - p_max * u_mid) / 2)
            th0 = (th0 + 180) % 360 - 180
            ap = S.link_arm(pts, b, th0, np.exp(u_mid), 0.3 * abs(b) * np.pi)
            sig = S.radial_sigma(ap, geom, img.beam_mean)
            segs = [("", np.ones(ap.r.size, bool))]
            if key == "IMLup":  # DSHARP III fit the inner (a) and outer (b) parts separately, split at 75 au
                segs = [("a", ap.r * d < 75), ("b", ap.r * d >= 75)]
            for suffix, sel in segs:
                if sel.sum() < 5:
                    continue
                sub = S.ArmPoints(ap.theta[sel], ap.r[sel], ap.r_err[sel])
                f = S.fit_log_spiral(sub, sig[sel])
                arms[f"arm{arm + 1}{suffix}"] = dict(r0_au=f.r0 * d, b=f.b, b_err=f.b_err, pitch=f.pitch,
                                                     pitch_err=f.pitch_err, n=f.n, theta_range=f.theta_range,
                                                     chi2_red=f.extra["chi2_red"],
                                                     theta=sub.theta.tolist(), r_au=(sub.r * d).tolist())
        # pair our arms with DSHARP III's by the sign-consistent pitch and R0 at matching windings
        out[key] = dict(fourier=dict(p_max=p_max, b=b, pitch=mu, pitch_range=(mu_lo, mu_hi),
                                     p_grid=p_grid[::4].tolist(), amp=amp[::4].tolist()),
                        arms=arms, published={k[1]: dict(r0_au=v[0], b=v[2], b_err=v[3], pitch=v[4], pitch_err=v[5])
                                              for k, v in HUANG_SPIRALS_LOG.items() if k[0] == key},
                        geometry=dict(inc=g["inc"], pa=g["pa"], method=g["method"]))
    return out


# ---------------------------------------------------------------------------
# 6. HD 163296 kinematics
# ---------------------------------------------------------------------------
KIN_CACHE = CACHE_DIR / "kinematics"


def kinematics(res: dict, refit: bool = False) -> dict:
    KIN_CACHE.mkdir(parents=True, exist_ok=True)
    cache = KIN_CACHE / "hd163296.json"
    if cache.exists() and not refit:
        return json.loads(cache.read_text())
    if not CO_CUBE_PATH.exists():
        return {"skipped": "CO cube not downloaded (fetch_data.py --skip-co)"}
    rng = np.random.default_rng(seed("HD163296 CO"))
    cube = A.read_cube(CO_CUBE_PATH, "HD163296", 5.0)
    rms = K.channel_rms(cube)
    g_cont = res["HD163296"]["geometry"]
    # line centres within the disc's velocity range: the line-of-sight speeds stay below ~10 km/s
    # of the systemic velocity outside 0.3", so search +/- 9 km/s around the cube's flux-weighted mean
    spec = np.nansum(cube.data, axis=(1, 2))
    v_mean = float(np.sum(spec * cube.velocity * (spec > 0)) / np.sum(spec * (spec > 0)))
    lm = K.line_centres(cube, rms, v_mean - 9.0, v_mean + 9.0)
    np.savez_compressed(KIN_CACHE / "linecentre.npz", v0=lm.v0, dv0=lm.dv0, peak=lm.peak, x=lm.x, y=lm.y, rms=rms)
    vsys0 = float(np.nanmedian(lm.v0))
    d_fit = TEAGUE2019_HD163296["distance_pc"]  # fit at Teague et al.'s distance; M scales with d
    # which half of the major axis is red, and which half of the minor axis is near?
    # Four ways the disc can sit on the sky; each is fitted (surface without taper, robust loss)
    # and judged by the median absolute velocity residual
    xx_s, yy_s = np.meshgrid(lm.x, lm.y)
    trials = []
    for pa in (g_cont["pa"] % 180.0, g_cont["pa"] % 180.0 + 180.0):
        for near in (1, -1):
            guess = K.DiscModel(g_cont["x0"], g_cont["y0"], g_cont["inc"], pa, 2.0, vsys0, 0.25, 1.0, 0.0, 4.0,
                                distance=d_fit, near=near)
            f = K.fit_rotation(lm, guess, (0.4, 3.0), 0.1, rng, mcmc=False, fix_taper=True)
            mv_t, r_t = f.model.velocity(xx_s, yy_s)
            sel_t = (r_t > 0.4) & (r_t < 3.0) & np.isfinite(lm.v0)
            mad = float(np.median(np.abs(lm.v0 - mv_t)[sel_t]))
            trials.append((mad, pa, near, f))
    trials.sort(key=lambda t: t[0])
    best = trials[0][3]
    orientation = [dict(pa_start=t[1], near=t[2], mad_ms=t[0] * 1e3, pa=t[3].model.pa % 360, mstar=t[3].model.mstar,
                        near_side_pa=t[3].model.near_side_pa) for t in trials]
    start = K.DiscModel(**{**best.model.__dict__, "z1": 0.01})
    inc_prior = (g_cont["inc"], max(g_cont["inc_err"], 0.1))
    main = K.fit_rotation(lm, start, (0.4, 3.0), 0.1, rng, inc_prior=inc_prior, mcmc=True)
    m = main.model
    variants = {}
    for name, kw in (("no inclination prior", dict(r=(0.4, 3.0), prior=None)),
                     ("r < 3.5\"", dict(r=(0.4, 3.5), prior=inc_prior)),
                     ("r < 2.5\"", dict(r=(0.4, 2.5), prior=inc_prior)),
                     ("r > 0.6\"", dict(r=(0.6, 3.0), prior=inc_prior)),
                     ("no taper (z1 = 0)", dict(r=(0.4, 3.0), prior=inc_prior, fix=True))):
        f = K.fit_rotation(lm, m, kw["r"], 0.1, rng, inc_prior=kw["prior"], mcmc=False, fix_taper=kw.get("fix", False))
        variants[name] = dict(mstar=f.model.mstar, inc=f.model.inc, z0=f.model.z0, phi=f.model.phi, chi2_red=f.chi2_red)
    sys_m = float(np.std([v["mstar"] for v in variants.values()] + [m.mstar]))
    # rotation curve on the fitted surface
    rc = K.rotation_curve(lm, m, np.arange(0.3, 3.6, 0.1), 0.1, rng)
    # Pressure minima: v_phi^2 = v_K^2 + (r/rho) dP/dr, so the gas turns from sub- to
    # super-Keplerian (dv/v crossing zero upwards) where the pressure is lowest -- a gas gap.
    # Restricted to 0.5"-2.9", away from the beam-smeared centre and the tapering outer edge.
    dvv = rc["v_phi"] / rc["v_kep"] - 1.0
    crossings = []
    for j in range(dvv.size - 1):
        if 0.5 <= rc["r"][j] <= 2.9 and dvv[j] < 0 <= dvv[j + 1]:
            r0 = rc["r"][j] + (0 - dvv[j]) * (rc["r"][j + 1] - rc["r"][j]) / (dvv[j + 1] - dvv[j])
            crossings.append(float(r0))
    # residual map and the two kinks of Pinte et al. 2020
    xx, yy = np.meshgrid(lm.x, lm.y)
    mv, r_map = m.velocity(xx, yy)
    resid = lm.v0 - mv
    np.savez_compressed(KIN_CACHE / "model.npz", vmodel=mv, resid=resid, r=r_map)
    kinks = {}
    d_gaia = gaia_distances()["HD163296"]["distance"]
    for name in ("HD163296 #1", "HD163296 #2"):
        sep, pa_k, snr, r_gap = PINTE2020_KINKS[name]
        xk, yk = sep * np.sin(np.radians(pa_k)), sep * np.cos(np.radians(pa_k))
        r_mid, phi_mid, _ = m.disc_coords(xk, yk, surface=False)
        v_planet = float(m.vsys + m.v_kep(r_mid) * np.cos(phi_mid) * np.sin(np.radians(m.inc)))
        z_s = float(m.height(r_mid))
        xs, ys = m.sky(r_mid, phi_mid, z_s)  # the same orbit position on the front emission surface
        # beam-averaged residual there vs. along the same annulus of the surface
        r_s, phi_s, _ = m.disc_coords(xx, yy)
        ann = (np.abs(r_s - r_mid) < 0.05) & np.isfinite(resid)
        beam = 0.1
        near = np.hypot(xx - xs, yy - ys) < beam / 2
        # sample the annulus in beam-sized azimuth bins
        nb = max(8, int(2 * np.pi * r_mid / beam))
        edges = np.linspace(-np.pi, np.pi, nb + 1)
        vals = [np.nanmean(resid[ann & (phi_s >= a) & (phi_s < b)]) for a, b in zip(edges[:-1], edges[1:])
                if np.sum(ann & (phi_s >= a) & (phi_s < b)) > 3]
        vals = np.array(vals)
        here = float(np.nanmean(resid[near])) if np.any(near & np.isfinite(resid)) else None
        kinks[name] = dict(sep=sep, pa=pa_k, r_mid_arcsec=float(r_mid), r_mid_au=float(r_mid * d_gaia),
                           phi_mid_deg=float(np.degrees(phi_mid)), v_planet=v_planet, surface_xy=(float(xs), float(ys)),
                           resid_here=here, annulus_mean=float(np.nanmean(vals)), annulus_std=float(np.nanstd(vals)),
                           z_score=(here - float(np.nanmean(vals))) / float(np.nanstd(vals)) if here is not None else None,
                           v_kep=float(m.v_kep(r_mid)), snr_pinte=snr, gap_au=r_gap)
    out = dict(
        rms_mjy=rms * 1e3, n_valid=int(np.isfinite(lm.v0).sum()), v_mean=v_mean,
        orientation_trials=orientation, near_side_pa=m.near_side_pa, rotation_on_sky=m.rotation_on_sky,
        fit=dict(x0_mas=m.x0 * 1e3, y0_mas=m.y0 * 1e3, inc=m.inc, pa=m.pa % 360, mstar=m.mstar, vsys=m.vsys,
                 z0_mas=m.z0 * 1e3, phi=m.phi, z1_mas=m.z1 * 1e3, psi=m.psi, distance=d_fit,
                 errors={k: (v * 1e3 if k in ("x0", "y0", "z0", "z1") else v) for k, v in main.errors.items()},
                 jitter_ms=main.jitter * 1e3, chi2_red=main.chi2_red, n_points=main.n_points,
                 n_clipped=main.extra["n_clipped"], inc_prior=inc_prior, r_range=main.r_range),
        mstar_systematic=sys_m, variants=variants,
        mstar_at_gaia=dict(distance=d_gaia, mstar=m.mstar * d_gaia / d_fit),
        published=TEAGUE2019_HD163296,
        rotation_curve={k: v.tolist() for k, v in rc.items()},
        pressure_minima_au=[c * gaia_distances()["HD163296"]["distance"] for c in crossings],
        pressure_minima_au_101pc=[c * 101.0 for c in crossings],
        kinks=kinks,
        resid_mad_ms=float(np.nanmedian(np.abs(resid[(r_map > 0.4) & (r_map < 3.0)])) * 1e3),
    )
    cache.write_text(json.dumps(D._jsonable(out), indent=1))
    return out


# ---------------------------------------------------------------------------
# Dust masses (optically thin, order-of-magnitude)
# ---------------------------------------------------------------------------
def dust_masses(res: dict) -> dict:
    """M_dust = F_nu d^2 / (kappa_nu B_nu(T)), the standard optically thin estimate with
    T = 20 K and kappa = 2.3 cm^2/g at 230 GHz scaled as nu^0.4 (Andrews et al. 2013's
    convention, used in most surveys), and with the DSHARP opacity (0.43 cm^2/g). Optically
    thick inner regions make both lower limits."""
    from config import PC_CGS

    out = {}
    for t in TARGETS:
        key = t.key
        o = res.get(key)
        if o is None:
            flux = L.dsharp1_images()[key].flux_mjy * 1e-3
            nu = L.dsharp1_images()[key].freq_ghz * 1e9
        else:
            flux, nu = o["flux"], o["freq_hz"]
        d = gaia_distances()[key]["distance"] * PC_CGS
        b20 = float(A.planck(nu, 20.0))
        kappa_std = 2.3 * (nu / 230e9) ** 0.4
        m_std = flux * 1e-23 * d**2 / (kappa_std * b20) / M_EARTH_CGS
        m_dsharp = flux * 1e-23 * d**2 / (0.43 * b20) / M_EARTH_CGS
        out[key] = dict(flux_mjy=flux * 1e3, m_dust_earth=m_std, m_dust_earth_dsharp_opacity=m_dsharp)
    return out


def run(refit: bool = False) -> dict:
    ensure_dirs()
    t0 = time.time()
    res = run_discs(refit)
    print(f"[{time.time() - t0:5.0f}s] discs done")
    out = dict(
        distances=gaia_distances(),
        geometry=compare_geometry(res),
        images=compare_images(res),
        features=compare_features(res),
        widths=compare_widths(res),
        planets=planet_masses(res),
        exoplanets=known_exoplanets(),
        snowlines=snowline_test(res),
        dust_masses=dust_masses(res),
        multiples=KURTOVIC_GEOMETRY,
    )
    print(f"[{time.time() - t0:5.0f}s] comparisons, planets, snow lines done")
    out["spirals"] = spiral_analysis(res)
    print(f"[{time.time() - t0:5.0f}s] spirals done")
    out["kinematics"] = kinematics(res, refit)
    print(f"[{time.time() - t0:5.0f}s] kinematics done")
    out["discs"] = {k: {kk: v for kk, v in o.items() if kk not in ("matches",)} for k, o in res.items()}
    RESULTS_PATH.write_text(json.dumps(D._jsonable(out), indent=1))
    return out


if __name__ == "__main__":
    run()
