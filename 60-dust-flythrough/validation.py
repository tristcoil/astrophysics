"""Checks of our map processing against the published paper and its posterior samples.

1. ``fig4_check``: read the published Figure 4 of Edenhofer et al. (2024)
   (the face-on view, from the arXiv source) back into numbers. Its colour
   bar is linear: the tick marks give rows -> A_V, the bar's own pixels give
   grey -> A_V, and the panel's frame gives pixels -> (x, y). The result is
   compared pixel by pixel with our projection.
2. ``holder_test``: the dust density between voxel centres is defined by
   interpolating ln(density) (the authors' convention). Interpolating the
   logarithm of the posterior *mean* is not the same as averaging the
   interpolated *samples*: by Hölder's inequality the latter is never larger,
   and the gap grows with the scatter between samples. We measure the gap on
   four shell pairs of the 12 released samples (fetched by byte range).
"""

from __future__ import annotations

import json

import healpy as hp
import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates

import config as C
import grids
from dustmap import DustSphere

FIG4 = C.LIT_DIR / "2308.01295" / "res" / "galactic.jpg"
#: Layout of the face-on panel and its colour bar in that JPEG (pixels),
#: found from the dark frame lines (see ``_frame``): panel frame columns
#: 341-344 and 2325-2328, rows 38-41 and 2022-2025; colour-bar frame columns
#: 2449-2453 and 2525-2528.
PANEL = {"c0": 344.5, "c1": 2324.5, "r0": 41.5, "r1": 2021.5, "half": 1250.0}
BAR = {"left": 2451, "right": 2528, "inner": (2462, 2515), "ticks": (1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2)}


def _bar_calibration(img: np.ndarray) -> dict:
    """Rows -> A_V from the tick marks; grey -> A_V from the bar itself."""
    strip = img[:, BAR["right"] + 3: BAR["right"] + 14]
    dark = np.flatnonzero((strip < 100).mean(axis=1) > 0.6)
    groups = np.split(dark, np.flatnonzero(np.diff(dark) > 2) + 1)
    rows = np.array([g.mean() for g in groups if g.size >= 2])[: len(BAR["ticks"])]
    slope, icpt = np.polyfit(rows, BAR["ticks"], 1)
    resid = np.polyval((slope, icpt), rows) - np.array(BAR["ticks"])
    frame = np.flatnonzero(img[:, BAR["left"]] < 100)
    frame = frame[(frame > 200) & (frame < 2100)]
    top, bottom = frame.min(), frame.max()
    r = np.arange(top + 6, bottom - 4)
    grey = np.median(img[r][:, BAR["inner"][0]:BAR["inner"][1]], axis=1)
    value = slope * r + icpt
    order = np.argsort(grey)
    return {"slope": slope, "intercept": icpt, "tick_rms": float(np.std(resid)),
            "vmin": float(slope * bottom + icpt), "vmax": float(slope * top + icpt),
            "grey": grey[order], "value": value[order]}


def fig4_check() -> dict:
    img = np.asarray(Image.open(FIG4).convert("L")).astype(float)
    cal = _bar_calibration(img)
    p = PANEL
    rows = np.arange(int(np.ceil(p["r0"])), int(p["r1"]) + 1)
    cols = np.arange(int(np.ceil(p["c0"])), int(p["c1"]) + 1)
    panel = img[np.ix_(rows, cols)]
    grey_u, inv = np.unique(cal["grey"], return_inverse=True)
    val_u = np.bincount(inv, weights=cal["value"]) / np.bincount(inv)
    theirs = np.interp(panel, grey_u, val_u)
    px = (cols + 0.5 - p["c0"]) / (p["c1"] - p["c0"]) * 2 * p["half"] - p["half"]
    py = p["half"] - (rows + 0.5 - p["r0"]) / (p["r1"] - p["r0"]) * 2 * p["half"]
    PX, PY = np.meshgrid(px, py)
    g = grids.build_main()
    xs, h = g["xs"], float(g["h"])
    ours = map_coordinates(C.ZGR23_TO_AV * g["proj_xy"], [(PX - xs[0]) / h, (PY - xs[0]) / h], order=1)
    R = np.hypot(PX, PY)
    ok = (R < 1200) & (theirs < 0.92 * cal["vmax"]) & (ours < 0.92 * cal["vmax"])
    sig = ok & (theirs > 0.2)
    fit = np.polyfit(theirs[ok], ours[ok], 1)
    by_r = []
    for lo in range(0, 1200, 200):
        s = sig & (R >= lo) & (R < lo + 200)
        by_r.append({"r_lo": lo, "r_hi": lo + 200, "ratio": float(np.median(ours[s] / theirs[s])), "n": int(s.sum())})
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    q999 = float(np.quantile(C.ZGR23_TO_AV * g["proj_xy"][np.hypot(X, Y) < C.R_OUTER_MAIN], 0.999))
    out = {"bar_vmin": cal["vmin"], "bar_vmax": cal["vmax"], "tick_rms": cal["tick_rms"],
           "our_q999": q999, "pixels": int(ok.sum()), "pearson_r": float(np.corrcoef(theirs[ok], ours[ok])[0, 1]),
           "slope": float(fit[0]), "offset": float(fit[1]),
           "median_ratio": float(np.median(ours[sig] / theirs[sig])), "ratio_by_distance": by_r}
    np.savez_compressed(C.CACHE_DIR / "fig4_pixels.npz", theirs=theirs.astype(np.float32),
                        ours=ours.astype(np.float32), px=px, py=py)
    return out


# ---------------------------------------------------------------------------
# Posterior samples: mean-of-interpolations vs interpolation-of-mean
# ---------------------------------------------------------------------------

def _sample_shell(s: int, shell: int) -> np.ndarray:
    path = C.EDENHOFER_DIR / "sample_shells" / f"s{s:02d}_r{shell:03d}.f32"
    return np.fromfile(path, dtype=">f4").astype(np.float32)


def holder_test(n_points: int = 200_000, seed: int = 20261019) -> dict:
    sphere = DustSphere("main")
    rng = np.random.default_rng(seed)
    rows = []
    for i in C.SAMPLE_SHELL_PAIRS:
        lo = np.stack([_sample_shell(s, i) for s in range(12)])
        hi = np.stack([_sample_shell(s, i + 1) for s in range(12)])
        mean_lo, mean_hi = sphere.mean[i], sphere.mean[i + 1]
        # (1) the released mean is the average of the samples at the voxel centres
        dev = max(np.max(np.abs(lo.mean(0) / mean_lo - 1)), np.max(np.abs(hi.mean(0) / mean_hi - 1)))
        # (2) random points between the two spheres
        l = rng.uniform(0, 360, n_points)
        b = np.degrees(np.arcsin(rng.uniform(-1, 1, n_points)))
        t = rng.uniform(0, 1, n_points)
        pix, w = hp.get_interp_weights(sphere.nside, l, b, nest=True, lonlat=True)
        def geo(a_lo, a_hi):
            la = (np.log(a_lo[..., pix]) * w).sum(-2)
            lb = (np.log(a_hi[..., pix]) * w).sum(-2)
            return np.exp((1 - t) * la + t * lb)
        from_mean = geo(mean_lo, mean_hi)
        from_samples = geo(lo, hi).mean(0)
        # the same with the samples' log-scatter, for the log-normal estimate
        ratio = from_mean / from_samples
        rel_scatter = np.median(lo.std(0) / lo.mean(0))
        rows.append({"shell": i, "r_pc": float(sphere.radii[i]), "max_mean_vs_samples_at_centres": float(dev),
                     "ratio_median": float(np.median(ratio)), "ratio_p16": float(np.quantile(ratio, 0.16)),
                     "ratio_p84": float(np.quantile(ratio, 0.84)),
                     "ratio_of_sums": float(from_mean.sum() / from_samples.sum()),
                     "voxel_rel_scatter_median": float(rel_scatter)})
    return {"pairs": rows}


# ---------------------------------------------------------------------------
# The dereddening tool against independent extinctions
# ---------------------------------------------------------------------------

def extinction_checks(seed: int = 7) -> dict:
    """Our integrated extinction vs dustmaps' integrated mode, cluster photometry and SFD."""
    import healpy as hp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from dustmaps.config import config as dmc
    import catalogs
    from deredden import Extinction, coefficient

    dmc["data_dir"] = str(C.DUSTMAPS_DIR)
    ext = Extinction()
    s = ext.sphere
    out = {}

    # (a) dustmaps' integrated mode: it puts the sum to each voxel's outer edge at the voxel centre
    rng = np.random.default_rng(seed)
    n = 100_000
    l = rng.uniform(0, 360, n)
    b = np.degrees(np.arcsin(rng.uniform(-1, 1, n)))
    d = np.exp(rng.uniform(np.log(80), np.log(1240), n))
    ours = ext.zgr23(l, b, d)
    link = C.DUSTMAPS_DIR / "edenhofer_2023"
    link.mkdir(exist_ok=True)
    target = link / C.EDENHOFER_FILES["main"][0]
    if not target.exists():
        target.symlink_to(s.path)
    from dustmaps.edenhofer2023 import Edenhofer2023Query
    q = Edenhofer2023Query(integrated=True)
    theirs = q.query(SkyCoord(l=l * u.deg, b=b * u.deg, distance=d * u.pc, frame="galactic"))
    del q
    k = np.clip(np.searchsorted(s.bounds, d, side="right") - 1, 0, s.nshell - 1)
    half = 0.5 * np.diff(s.bounds)[k]
    ours_shift = ext.zgr23(l, b, d + half)
    ok = ours > 1e-3
    out["dustmaps_integrated"] = {
        "median_ratio_dustmaps_over_ours": float(np.median(theirs[ok] / ours[ok])),
        "p99_ratio": float(np.quantile(theirs[ok] / ours[ok], 0.99)),
        "median_ratio_vs_ours_half_voxel_further": float(np.median(theirs[ok] / ours_shift[ok])),
        "p99_ratio_shifted": float(np.quantile(theirs[ok] / ours_shift[ok], 0.99)),
        "max_abs_diff_av_mag": float(C.ZGR23_TO_AV * np.max(np.abs(theirs - ours))),
        "max_abs_diff_shifted_av_mag": float(C.ZGR23_TO_AV * np.max(np.abs(theirs - ours_shift)))}

    # (b) Hunt & Reffert (2023) clusters: photometric A_V from their isochrone network
    hr = catalogs.hunt_reffert()
    sel = hr["reliable"] & (hr["dist50"] > 100) & (hr["dist50"] < 1200) & (hr["Type"] == "o")
    c = hr[sel].copy()
    e = ext.zgr23(c["GLON"].to_numpy(), c["GLAT"].to_numpy(), c["dist50"].to_numpy())
    c["A550_map"] = coefficient(550) * e
    c["AV_map_540"] = C.ZGR23_TO_AV * e
    dav = c["A550_map"] - c["AV50"]
    big = c["AV50"] > 0.3
    out["hunt_reffert"] = {
        "n_clusters": int(len(c)), "median_diff_mag": float(np.median(dav)),
        "robust_sd_mag": float(1.4826 * np.median(np.abs(dav - np.median(dav)))),
        "median_ratio_AV_gt_0p3": float(np.median(c["A550_map"][big] / c["AV50"][big])),
        "spearman": float(c[["A550_map", "AV50"]].corr(method="spearman").iloc[0, 1]),
        "by_distance": [{"lo": lo, "hi": hi,
                         "median_diff": float(np.median(dav[(c.dist50 >= lo) & (c.dist50 < hi)])),
                         "n": int(((c.dist50 >= lo) & (c.dist50 < hi)).sum())}
                        for lo, hi in ((100, 300), (300, 600), (600, 900), (900, 1200))]}
    c[["Name", "GLON", "GLAT", "dist50", "AV16", "AV50", "AV84", "diffAV50", "A550_map", "AV_map_540",
       "age_myr"]].to_parquet(C.CACHE_DIR / "extinction_clusters.parquet")

    # (c) SFD (1998) at high latitude, where the whole dust column lies within the map
    from dustmaps.sfd import SFDQuery
    sfd = SFDQuery()
    nside = 64
    lp, bp = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), lonlat=True)
    hi_b = np.abs(bp) > 30
    col = ext.zgr23(lp[hi_b], bp[hi_b], np.full(hi_b.sum(), 1248.0))
    ebv_sfd = sfd(SkyCoord(l=lp[hi_b] * u.deg, b=bp[hi_b] * u.deg, frame="galactic"))
    ebv_map = (coefficient(440) - coefficient(550)) * col
    sf11 = 0.86 * ebv_sfd
    good = ebv_sfd < 0.3
    slope = np.polyfit(sf11[good], ebv_map[good], 1)
    out["sfd_high_latitude"] = {
        "n_pixels": int(good.sum()), "median_ratio_map_over_sf11": float(np.median(ebv_map[good] / sf11[good])),
        "fit_map_vs_sf11": [float(v) for v in slope],
        "pearson": float(np.corrcoef(sf11[good], ebv_map[good])[0, 1])}
    np.savez_compressed(C.CACHE_DIR / "extinction_sfd.npz", l=lp[hi_b], b=bp[hi_b], ebv_map=ebv_map, sf11=sf11)
    return out


def run() -> dict:
    out = {"fig4": fig4_check(), "holder": holder_test(), "extinction": extinction_checks()}
    (C.CACHE_DIR / "validation.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    res = run()
    print(json.dumps(res, indent=1))
