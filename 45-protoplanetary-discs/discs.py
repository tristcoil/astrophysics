"""The per-disc pipeline: orientation, profile, gaps and rings, widths, size -- all our own.

For each disc, from its DSHARP image alone (published values are used only afterwards,
to grade the results):

1. Find the disc (peak of a smoothed image) and its size (curve of growth in circles);
   measure the noise far from it.
2. A first orientation from a 2-D Gaussian fit (DSHARP II's "G" method).
3. A first radial profile; its significant minima and maxima are the candidate gaps and rings.
4. For each high-contrast, well-resolved feature: trace its ridge or trough along spokes
   and fit an ellipse (DSHARP II's "E" method). The weighted mean of these ellipses is the
   disc's orientation; the tracing and fitting are repeated once in that orientation.
5. The final profile (averaged over the same azimuths as DSHARP II: wedges along the
   major axis for highly inclined discs, minus asymmetric regions), its gaps and rings,
   their widths and depths, and the disc's dust radius R_dust (95% of the flux).

Radial bins are 1 au wide at DSHARP's distance (as in DSHARP II), so that the "intensity
of the bin containing r" in the width and depth definitions means the same thing.
"""

from __future__ import annotations

import json
from dataclasses import asdict

import numpy as np
from scipy import ndimage

import alma_images as A
import geometry as G
import literature as L
import profiles as P
from config import CACHE_DIR, HLTAU, TARGET_BY_KEY, Target
from distances import gaia_distances

DISC_CACHE = CACHE_DIR / "discs"

# Feature selection for ellipse fits
E_MIN_SIGMA = 10.0  # the extremum must stand out by this many profile errors
E_MIN_R_BEAMS = 2.5  # and lie at least this many beams from the centre
E_MIN_POINTS = 12
E_MIN_FRACTION = 0.5  # of the spokes tried must give a valid extremum (ALMA 2015 used 55%)
E_MIN_MINOR_SEP_BEAMS = 0.6  # projected distance to the neighbouring features along the minor axis
SECURE_SIGMA = 5.0  # features at 3-5 sigma are reported as marginal
MIN_SIGMA = 3.0
FEATURE_RMAX_RDUST = 1.15  # features are searched for out to this multiple of R_dust
PA_KNOWN_DEG = 10.0  # Gaussian PA error below which the first-pass profile may use a wedge


def dsharp_distance(key: str) -> float:
    """The distance DSHARP II used (DSHARP I Table 1; 147 pc for HL Tau)."""
    if key == "HLTau":
        return 147.0
    return L.dsharp1_stars()[key].distance.value


def target(key: str) -> Target:
    return HLTAU if key == "HLTau" else TARGET_BY_KEY[key]


def extent_mask(img: A.Image, xc: float, yc: float, nsigma: float = 4.0, smooth_beams: float = 1.5):
    """Pixels belonging to the disc: where the image, smoothed to ``smooth_beams`` beams, is
    above ``nsigma`` times its own noise, keeping the region around (xc, yc) and every other
    region of at least two beams whose centroid lies near the centre -- a ring cut off by an
    empty gap is an annulus centred on the star, while a background source or a noise
    spike is not."""
    sm = ndimage.gaussian_filter(img.data, smooth_beams * img.beam_mean / img.pix / 2.3548)
    xx, yy = img.grids()
    far = np.hypot(xx - xc, yy - yc) > 0.6 * (img.x.max() - img.x.min()) / 2
    noise = A.robust_std(sm[far]) if far.sum() > 1000 else A.robust_std(sm)
    lab, n = ndimage.label(sm > nsigma * noise)
    j = np.argmin(np.abs(img.y - yc))
    i = np.argmin(np.abs(img.x - xc))
    keep = np.zeros(n + 1, dtype=bool)
    if lab[j, i] > 0:
        keep[lab[j, i]] = True
    area = np.bincount(lab.ravel(), minlength=n + 1)
    cx = np.bincount(lab.ravel(), weights=xx.ravel(), minlength=n + 1) / np.maximum(area, 1)
    cy = np.bincount(lab.ravel(), weights=yy.ravel(), minlength=n + 1) / np.maximum(area, 1)
    big = area >= 2 * img.beam_pix
    near = np.hypot(cx - xc, cy - yc) < max(0.15, 2 * img.beam_mean)
    keep |= big & near
    keep[0] = False
    mask = keep[lab]
    if not mask.any():
        return np.hypot(xx - xc, yy - yc) < 3 * img.beam_mean
    return mask


def locate(img: A.Image) -> tuple[float, float, float, float]:
    """Disc centre (peak of the image smoothed over ~3 beams), the radius of its 4-sigma
    extent and the flux inside it (circles; refined later with the disc's orientation)."""
    sm = ndimage.gaussian_filter(img.data, 3.0 * img.beam_mean / img.pix / 2.3548)
    j, i = np.unravel_index(np.nanargmax(sm), sm.shape)
    xc, yc = float(img.x[i]), float(img.y[j])
    xx, yy = img.grids()
    rr = np.hypot(xx - xc, yy - yc)
    r_ext = float(rr[extent_mask(img, xc, yc)].max())
    total = float(img.data[rr <= 1.1 * r_ext].sum() / img.beam_pix)
    return xc, yc, r_ext, total


def curve_of_growth(img: A.Image, geom: G.DiscGeometry, r_max: float, step: float):
    """Flux inside deprojected radius r (ellipses with the disc's orientation)."""
    xx, yy = img.grids()
    r, _ = geom.polar(xx, yy)
    edges = np.arange(0.0, r_max + step, step)
    f = np.cumsum(np.histogram(r, bins=edges, weights=img.data)[0]) / img.beam_pix
    return edges[1:], f


def r_dust(img: A.Image, geom: G.DiscGeometry, r_max: float, step: float, frac: float = 0.95):
    """DSHARP II's R_dust: the radius enclosing 95% of the disc's total flux, in elliptical
    apertures matching its orientation. DSHARP II takes the total as "the asymptotic value of
    the flux curve"; here it is the median of the curve of growth between 1.3 and 1.8 times
    the disc's 4-sigma extent (deprojected). Closer in, faint outer emission below the 4-sigma
    contour is still adding flux (Elias 27's curve rises by 7% between 1.0 and 1.4 times the
    extent); the maximum of the curve would instead chase noise far outside the disc."""
    rr, f = curve_of_growth(img, geom, r_max, step)
    xx, yy = img.grids()
    r, _ = geom.polar(xx, yy)
    r_ext = float(r[extent_mask(img, geom.x0, geom.y0)].max())
    sel = (rr >= 1.3 * r_ext) & (rr <= min(1.8 * r_ext, r_max))
    if sel.sum() < 3:
        sel = (rr >= r_ext) & (rr <= r_max)
    total = float(np.median(f[sel])) if sel.sum() else float(f[rr <= r_max].max())
    return float(rr[np.argmax(f >= frac * total)]), total


def _windows(ext: list[P.Extremum], k: int, beam: float) -> tuple[float, float]:
    """Spoke range for tracing feature k: up to 80% of the way to its neighbours."""
    e = ext[k]
    lo = e.r - 0.8 * (e.r - ext[k - 1].r) if k > 0 else max(e.r - 1.5 * beam, 0.5 * e.r)
    hi = e.r + 0.8 * (ext[k + 1].r - e.r) if k + 1 < len(ext) else e.r + 1.5 * beam
    return lo, hi


def ellipse_fits(img: A.Image, geom: G.DiscGeometry, ext: list[P.Extremum], noise: float, t: Target,
                 rng: np.random.Generator) -> list[tuple[P.Extremum, G.EllipseFit]]:
    out = []
    d_ds = dsharp_distance(t.key)
    for k, e in enumerate(ext):
        if e.prominence < E_MIN_SIGMA or e.r < E_MIN_R_BEAMS * img.bmaj:
            continue
        if t.spiral_range_au is not None:
            # spiral arms bias a ring's ridge; inside the arms only gaps with the arm crossings
            # masked are traced (as DSHARP II did for Elias 27's D69)
            lo_au, hi_au = t.spiral_range_au
            if lo_au <= e.r * d_ds <= hi_au and not (e.kind == "D" and t.ridge_exclude):
                continue
        lo, hi = _windows(ext, k, img.bmaj)
        if hi - lo < 0.75 * img.bmaj:
            continue
        # along the minor axis everything is foreshortened by cos i: a feature whose
        # neighbours come closer than ~half a beam there blends with them, and its traced
        # ellipse comes out too round (HD 142666's rings at i = 62 deg)
        sep = min(e.r - ext[k - 1].r if k > 0 else np.inf, ext[k + 1].r - e.r if k + 1 < len(ext) else np.inf)
        if np.isfinite(sep) and sep * np.cos(np.radians(geom.inc)) < E_MIN_MINOR_SEP_BEAMS * img.beam_mean:
            continue
        n_spokes = max(8, int(round(2 * np.pi * e.r / img.bmaj)))
        pts = G.extract_ridge(img, geom, e.kind, lo, hi, 360.0 / n_spokes, noise, t.ridge_exclude)
        if pts.r.size < E_MIN_POINTS or pts.r.size < E_MIN_FRACTION * pts.n_tried:
            continue
        fit = G.fit_ellipse(pts, geom, e.r, (0.0, 180.0), rng)
        if not (np.isfinite(fit.errors["inc"]) and fit.errors["inc"] > 0):
            continue
        out.append((e, fit))
    return out


def _extrema(prof: P.Profile, img: A.Image, noise: float, r_feat: float) -> list[P.Extremum]:
    ext = P.find_extrema(prof, img.beam_mean, 0.75 * img.beam_mean, MIN_SIGMA, noise=noise)
    return [e for e in ext if e.r <= r_feat]


def label(kind: str, r_arcsec: float, d_pc: float) -> str:
    return f"{kind}{int(round(r_arcsec * d_pc))}"


def match_features(ours: list[dict], key: str, beam: float, d_ds: float) -> list[dict]:
    """Pair each DSHARP II feature of this disc with our nearest feature of the same kind."""
    out = []
    used = set()
    for f in [f for f in L.dsharp2_features() if f.key == key and f.r_mas.value is not None]:
        r_pub = f.r_mas.value / 1e3
        tol = (1.5 if f.r_mas.approx else 0.75) * beam
        cands = [(abs(o["r"] - r_pub), i) for i, o in enumerate(ours) if o["kind"] == f.kind and i not in used]
        cands = [c for c in cands if c[0] <= max(tol, 2.0 / d_ds)]
        entry = dict(label=f.label, kind=f.kind, method=f.method, r_pub=r_pub,
                     r_pub_err=(f.r_mas.err / 1e3) if f.r_mas.err else (1.0 / d_ds if not f.r_mas.approx else np.nan),
                     approx=f.r_mas.approx, width_pub=f.width_au.value, width_pub_err=f.width_au.err,
                     width_upper=f.width_au.upper_limit, depth_pub=f.depth.value, depth_pub_err=f.depth.err)
        if cands:
            _, i = min(cands)
            used.add(i)
            entry.update(ours=ours[i]["label"], r_ours=ours[i]["r"], r_ours_err=ours[i]["r_err"],
                         sigma=ours[i]["sigma"], our_method=ours[i]["method"])
        out.append(entry)
    return out


def analyse(key: str, rng: np.random.Generator, verbose: bool = True) -> dict:
    t = target(key)
    d_ds = dsharp_distance(key)
    d_new = gaia_distances()[key]["distance"]
    full = A.load_continuum(key)
    xc, yc, r_ext_c, _ = locate(full)
    noise = A.measure_rms(full, 1.5 * r_ext_c + 0.5, xc, yc)
    peak = float(full.data.max())
    half = min(1.9 * r_ext_c + 0.5, 0.5 * (full.x.max() - full.x.min()) - 0.1)
    img = full.cutout(half, xc, yc)
    del full
    bw = 1.0 / d_ds
    r_max = half - 0.05

    # step 2: Gaussian orientation
    gf = G.fit_gaussian(img, xc, yc, r_ext_c)
    geom_g = G.DiscGeometry(gf.x0, gf.y0, gf.inc, gf.pa, "G", gf.errors["x0"], gf.errors["y0"], gf.errors["inc"],
                            gf.errors["pa"])
    # steps 3-4: profile -> features -> ellipses -> orientation, twice
    geom = geom_g
    rd0, _ = r_dust(img, geom_g, r_max, bw)
    r_feat = min(FEATURE_RMAX_RDUST * rd0, r_max)  # beyond the dust disc's edge the profile is noise
    # a wedge along the major axis needs a known major axis: for nearly face-on discs whose
    # Gaussian PA is undetermined (Sz 114) the first pass averages over all azimuths
    wedge0 = t.wedge_deg if (np.isfinite(gf.errors["pa"]) and gf.errors["pa"] < PA_KNOWN_DEG) else None
    prof = P.radial_profile(img, geom, bw, r_max, wedge0, t.exclude_theta)
    ext = _extrema(prof, img, noise, r_feat)
    efits: list = []
    geom_e = None
    for _ in range(2):
        efits = ellipse_fits(img, geom, ext, noise, t, rng)
        if not efits:
            break
        geom_e = G.weighted_geometry([f for _, f in efits], (0.0, 180.0))
        geom = geom_e
        prof = P.radial_profile(img, geom, bw, r_max, t.wedge_deg, t.exclude_theta)
        rd0, _ = r_dust(img, geom, r_max, bw)
        r_feat = min(FEATURE_RMAX_RDUST * rd0, r_max)
        ext = _extrema(prof, img, noise, r_feat)
    # step 5: final features
    feats = []
    for e in ext:
        fit = next((f for ee, f in efits if ee.kind == e.kind and abs(f.a - e.r) < 0.5 * img.bmaj), None)
        r, r_err, method = (fit.a, fit.errors["a"], "E") if fit is not None else (e.r, bw / np.sqrt(12.0), "R")
        feats.append(dict(kind=e.kind, r=float(r), r_err=float(r_err), method=method, sigma=e.prominence,
                          contrast=e.contrast, intensity=e.intensity, secure=bool(e.prominence >= SECURE_SIGMA),
                          label=label(e.kind, r, d_ds), r_au=float(r * d_new), r_au_dsharp=float(r * d_ds)))
    feature_list = [(f["label"], f["kind"], f["r"]) for f in feats]
    pairs, perr = P.pair_measures_with_errors(prof, feature_list, 0.1 / d_ds, img.beam_mean, rng)
    rd, flux = r_dust(img, geom, r_max, bw)
    out = dict(
        key=key, name=t.name, d_dsharp=d_ds, d_gaia=d_new, centre_guess=(xc, yc), noise=noise, peak=peak,
        flux=flux, r_dust=rd, r_dust_au=rd * d_new, r_dust_au_dsharp=rd * d_ds,
        beam=(img.bmaj, img.bmin, img.bpa), freq_hz=img.freq_hz, half_width=half,
        geometry_gauss=geom_g.as_dict(), gauss_fwhm=(gf.fwhm_maj, gf.fwhm_min),
        geometry_ellipse=geom_e.as_dict() if geom_e is not None else None,
        geometry=geom.as_dict(),
        ellipse_fits=[dict(label=label(e.kind, f.a, d_ds), kind=e.kind, a=f.a, inc=f.inc, pa=f.pa, x0=f.x0, y0=f.y0,
                           sigma=f.sigma, errors=f.errors, n_points=f.n_points, n_tried=f.n_tried) for e, f in efits],
        features=feats,
        pairs=[dict(asdict(m), errors=perr.get((m.gap, m.ring), {})) for m in pairs],
        matches=match_features(feats, key, img.bmaj, d_ds) if key in {f.key for f in L.dsharp2_features()} else [],
        wedge_deg=t.wedge_deg, exclude_theta=t.exclude_theta,
    )
    DISC_CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(DISC_CACHE / f"{key}_profile.npz", r=prof.r, intensity=prof.intensity, error=prof.error,
                        scatter=prof.scatter, nbeam=prof.nbeam, bin_width=prof.bin_width)
    (DISC_CACHE / f"{key}.json").write_text(json.dumps(_jsonable(out), indent=1))
    if verbose:
        g = geom
        print(f"  {key:9s} {g.method} inc {g.inc:6.2f}+/-{g.inc_err:4.2f} PA {g.pa:6.2f}+/-{g.pa_err:4.2f} "
              f"centre ({g.x0 * 1e3:+6.1f},{g.y0 * 1e3:+6.1f}) mas; {len(feats)} features "
              f"({sum(f['secure'] for f in feats)} secure), {len(efits)} ellipses; R_dust {rd * d_ds:5.1f} au")
    return out


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def load(key: str) -> dict | None:
    p = DISC_CACHE / f"{key}.json"
    return json.loads(p.read_text()) if p.exists() else None


def load_profile(key: str) -> P.Profile:
    z = np.load(DISC_CACHE / f"{key}_profile.npz")
    return P.Profile(z["r"], z["intensity"], z["error"], z["scatter"], np.zeros_like(z["r"]), z["nbeam"],
                     float(z["bin_width"]))
