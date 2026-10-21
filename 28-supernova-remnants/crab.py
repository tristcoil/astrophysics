"""The Crab: sixty years of filament motions, rewound to the explosion.

Inputs are the registered images from ``plates.py``. The main measurement
compares the POSS-I red plate (1951 November 7) with the Pan-STARRS r-band
stack (2009-2014, effective date 2011.9): both bands are dominated by the
filaments' hydrogen-alpha, [N II] and [S II] lines. Five other band-matched
pairs check it. Pairs that mix bands (e.g. a red plate against the
near-infrared i band, which sees mostly the smooth synchrotron nebula, not the
filaments) fail to give coherent vectors and are not used.
"""

from __future__ import annotations

import hashlib
import json

import numpy as np
from scipy import ndimage

import astrometry as A
import config
import kicks
import literature as L
import motions as M
import plates

#: Band-matched pairs, (first image, second image). Red: the filaments' H-alpha,
#: [N II] and [S II]; blue: [O III], H-beta and, on 103aO, [O II] -- plus a
#: larger share of the smooth synchrotron light.
PAIRS = {
    "red 1951-2012": ("crab_poss1_red", "crab_ps1_r"),
    "red 1951-1991": ("crab_poss1_red", "crab_poss2_red"),
    "red 1991-2012": ("crab_poss2_red", "crab_ps1_r"),
    "blue 1951-2012": ("crab_poss1_blue", "crab_ps1_g"),
    "blue 1951-1992": ("crab_poss1_blue", "crab_poss2_blue"),
    "blue 1992-2012": ("crab_poss2_blue", "crab_ps1_g"),
}
MAIN = "red 1951-2012"

#: Processing choices varied one at a time on the main pair, to see how much
#: the date depends on them (the spread is quoted as a systematic error).
VARIANTS = {
    "patch 14\"": {"patch": 14.0},
    "patch 30\"": {"patch": 30.0, "step": 8.0},
    "high-pass 8\"": {"highpass": 8.0},
    "high-pass 20\"": {"highpass": 20.0},
    "blur to 6\"": {"fwhm": 6.0},
    "CCD stretch x0.2": {"asinh_scale": 10.0},
    "CCD stretch x6": {"asinh_scale": 300.0},
    "stricter matches": {"min_peak": 0.65},
    "larger region": {"region_pad": 6.0},
}

#: Number of position-angle wedges, as in Martin et al. (2025).
N_WEDGES = 20
#: Asinh softening for the Pan-STARRS stacks (counts), before blurring.
PS1_ASINH = 50.0


def nebula_region(regs: dict, grid: A.Grid, gaia: A.GaiaStars) -> np.ndarray:
    """Where the nebula is: the Pan-STARRS r image with stars filled in,
    blurred to 5", above 10 sigma of the sky; the largest connected piece,
    holes filled."""
    cache = config.DERIVED_DIR / "crab_region.npy"
    if cache.exists():
        return np.load(cache)
    r = regs["crab_ps1_r"]
    radii = M.star_radii(r, grid, gaia)
    stars = M.star_mask(grid, gaia, r.epoch, r.fwhm_arcsec, radius_model=radii)
    filled = M.inpaint(r.image, stars | ~np.isfinite(r.image), 2.0 / grid.step)
    filled = ndimage.gaussian_filter(filled, config.COMMON_FWHM_ARCSEC / 2.3548 / grid.step)
    edge = filled[: int(60 / grid.step), :]  # a 1' strip along the south edge: sky
    sky, sigma = np.median(edge), A.robust_sigma(edge)
    region = ndimage.binary_opening(filled > sky + 10 * sigma, iterations=3)
    lab, n = ndimage.label(region)
    sizes = ndimage.sum(region, lab, range(1, n + 1))
    region = ndimage.binary_fill_holes(lab == (np.argmax(sizes) + 1))
    np.save(cache, region)
    return region


class Crab:
    """Registered images, the nebula region, and cached prepared images."""

    def __init__(self, refresh: bool = False):
        self.grid = plates.grid_for(config.CRAB)
        self.gaia = A.read_gaia("crab")
        self.regs = plates.register_all("crab", refresh=refresh, grid=self.grid)
        self.region = nebula_region(self.regs, self.grid, self.gaia)
        self._prep: dict = {}

    def prepared(self, key: str, *, fwhm=config.COMMON_FWHM_ARCSEC, highpass=config.HIGHPASS_ARCSEC,
                 asinh_scale=PS1_ASINH) -> M.Prepared:
        tag = (key, fwhm, highpass, asinh_scale if "ps1" in key else None)
        if tag not in self._prep:
            self._prep[tag] = M.prepare(self.regs[key], self.grid, self.gaia, self.region, fwhm=fwhm,
                                        highpass=highpass, asinh_scale=asinh_scale if "ps1" in key else None)
        return self._prep[tag]

    def vectors(self, a: str, b: str, *, patch=config.PATCH_ARCSEC, step=config.PATCH_STEP_ARCSEC,
                fwhm=config.COMMON_FWHM_ARCSEC, highpass=config.HIGHPASS_ARCSEC, asinh_scale=PS1_ASINH,
                region_pad=-4.0, min_peak=0.5, refresh=False) -> M.Vectors:
        """Tracked, quality-cut vectors from image a to image b (cached).

        Patch centres lie inside the nebula region shrunk by 4" (or grown,
        for a positive ``region_pad``)."""
        tag = f"{a}__{b}__p{patch}_s{step}_f{fwhm}_h{highpass}_a{asinh_scale}_r{region_pad}_m{min_peak}"
        cache = config.DERIVED_DIR / f"vec_{hashlib.sha1(tag.encode()).hexdigest()[:12]}.npz"
        if cache.exists() and not refresh:
            z = np.load(cache, allow_pickle=True)
            if str(z["tag"]) == tag:
                return M.Vectors((a, b), *(z[f] for f in ("xi", "eta", "dxi", "deta", "t1", "t2", "peak",
                                                           "second", "back_err")))
        pa = self.prepared(a, fwhm=fwhm, highpass=highpass, asinh_scale=asinh_scale)
        pb = self.prepared(b, fwhm=fwhm, highpass=highpass, asinh_scale=asinh_scale)
        it = int(round(abs(region_pad) / self.grid.step))
        pts = (ndimage.binary_erosion(self.region, iterations=it) if region_pad < 0
               else ndimage.binary_dilation(self.region, iterations=it))
        v = M.track(pa, pb, self.grid, pts, patch=patch, step=step)
        v = v.subset(M.quality_cut(v, min_peak=min_peak))
        np.savez_compressed(cache, tag=tag, **v.to_dict())
        return v


def _pair_summary(v: M.Vectors) -> tuple[dict, M.ExpansionFit]:
    fit = M.fit_expansion(v)
    used = v.subset(fit.inliers)
    t_ca, _, _ = M.closest_approach(used)
    mu = np.hypot(used.mu_xi, used.mu_eta) * 1000.0
    return {
        "n_vectors": int(len(v.xi)), "n_used": fit.n,
        "epoch_1": float(np.median(v.t1)), "epoch_2": float(np.median(v.t2)),
        "baseline_yr": float(np.median(v.dt)), "t_mid": float(np.median(used.tmid)),
        "T_fit": fit.T, "T_fit_err": float(fit.errors[2]),
        "centre_xi": fit.xi_c, "centre_eta": fit.eta_c,
        "centre_err": [float(fit.errors[0]), float(fit.errors[1])],
        "T_closest_approach": t_ca, "rms_mu_mas_yr": fit.rms_mu * 1000.0,
        "mu_median_mas_yr": float(np.median(mu)), "mu_p95_mas_yr": float(np.percentile(mu, 95)),
    }, fit


def wedge_dates(v: M.Vectors, xi_c: float, eta_c: float, n_wedges: int = N_WEDGES) -> dict:
    """Convergence date in position-angle wedges about a fixed centre.

    Within one wedge the centre is poorly constrained (all vectors point the
    same way), so the global centre is kept and only the slope of radial
    proper motion against radius is fitted, through the origin, with clipping."""
    T, r, pa = M.per_vector_dates(v, xi_c, eta_c)
    dx, dy = v.xi - xi_c, v.eta - eta_c
    mu_r = (v.mu_xi * dx + v.mu_eta * dy) / r
    out = {"pa": [], "T": [], "T_err": [], "n": []}
    width = 360.0 / n_wedges
    for k in range(n_wedges):
        sel = (pa >= k * width) & (pa < (k + 1) * width) & (r > 30)
        if sel.sum() < 15:
            continue
        x, y, tm = r[sel], mu_r[sel], v.tmid[sel]
        keep = np.ones(len(x), bool)
        for _ in range(5):
            slope = np.sum(x[keep] * y[keep]) / np.sum(x[keep] ** 2)
            res = y - slope * x
            s = 1.4826 * np.median(np.abs(res[keep]))
            keep = np.abs(res) < 4 * s
        slope = np.sum(x[keep] * y[keep]) / np.sum(x[keep] ** 2)
        res = y[keep] - slope * x[keep]
        slope_err = np.sqrt(np.sum(res**2) / max(1, keep.sum() - 1) / np.sum(x[keep] ** 2))
        age = 1.0 / slope
        out["pa"].append((k + 0.5) * width)
        out["T"].append(float(np.median(tm) - age))
        out["T_err"].append(float(age**2 * slope_err))
        out["n"].append(int(keep.sum()))
    return out


def radial_dates(v: M.Vectors, xi_c: float, eta_c: float, edges=(0, 60, 100, 140, 180, 260)) -> dict:
    """Median per-vector convergence date in rings about the centre."""
    T, r, _ = M.per_vector_dates(v, xi_c, eta_c)
    out = {"r_in": [], "r_out": [], "T_median": [], "n": []}
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (r >= lo) & (r < hi) & np.isfinite(T)
        if sel.sum() < 10:
            continue
        out["r_in"].append(lo), out["r_out"].append(hi)
        out["T_median"].append(float(np.median(T[sel]))), out["n"].append(int(sel.sum()))
    return out


def pulsar(fit: M.ExpansionFit, gaia: A.GaiaStars, t_mid: float) -> dict:
    """The Crab pulsar's Gaia DR3 motion, wound back to the explosion.

    Our coordinate origin is the pulsar's own position at Gaia's epoch
    2016.0, so its position at any date is just -mu (2016 - date)."""
    i = int(np.flatnonzero(gaia.source_id == config.CRAB_PULSAR_GAIA)[0])
    pmra, pmdec, plx = float(gaia.pmra[i]), float(gaia.pmdec[i]), float(gaia.parallax[i])
    at = lambda year: (pmra * (year - config.GAIA_EPOCH) / 1000.0, pmdec * (year - config.GAIA_EPOCH) / 1000.0)
    p1054, pT = at(1054.5), at(fit.T)
    off = lambda p: float(np.hypot(p[0] - fit.xi_c, p[1] - fit.eta_c))
    mu = float(np.hypot(pmra, pmdec))
    d = config.CRAB.distance_kpc
    pec = kicks.peculiar_transverse_velocity(config.CRAB.ra, config.CRAB.dec, d, pmra, pmdec)
    # Where is the pulsar *today* relative to the centre, and how fast must
    # it have moved to get there since 1054?
    now = at(t_mid)
    sep_now = off(now)
    return {
        "gaia_source_id": config.CRAB_PULSAR_GAIA, "pmra_mas_yr": pmra, "pmdec_mas_yr": pmdec,
        "parallax_mas": plx, "parallax_err_mas": 0.079,
        "mu_mas_yr": mu, "pa_deg": float(np.degrees(np.arctan2(pmra, pmdec)) % 360.0),
        "v_transverse_kms": mu * config.KMS_PER_MASYR_KPC * d,
        "v_peculiar_kms": pec["v_pec_kms"], "v_peculiar_pa_deg": pec["pa_deg"],
        "position_1054": list(p1054), "position_at_T": list(pT),
        "offset_1054_from_centre_arcsec": off(p1054), "offset_T_from_centre_arcsec": off(pT),
        "separation_from_centre_now_arcsec": sep_now,
        "speed_needed_since_1054_mas_yr": sep_now / (t_mid - 1054.5) * 1000.0,
    }


def acceleration(T: float, t_mid: float, t_true: float = L.CRAB_SN_YEAR) -> dict:
    """What the gap between the convergence date and 1054 means.

    If the filaments had coasted since 1054, today's speed v would equal their
    average speed r / age and the rewind would land on 1054. Landing later
    means v > r / age: they have sped up. The ratio of ages
    q = (t_mid - T) / (t_mid - 1054) is the average speed over the current
    speed. For a constant push since the explosion, r = v0 t + a t^2 / 2 and
    the fraction of today's speed that was gained is 2 (1 - q). A nebula
    driven by a pulsar of constant power grows as R ~ t^(6/5) (Chevalier &
    Fransson 1992), giving q = 5/6."""
    q = (t_mid - T) / (t_mid - t_true)
    return {
        "apparent_age_yr": t_mid - T, "true_age_yr": t_mid - t_true, "q": q,
        "speed_now_over_mean": 1.0 / q, "fraction_gained_constant_push": 2.0 * (1.0 - q),
        "T_if_R_t65": t_mid - (5.0 / 6.0) * (t_mid - t_true),
    }


def analyse(refresh: bool = False) -> dict:
    cache = config.DERIVED_DIR / "crab_results.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())
    crab = Crab()
    res = {"pairs": {}, "plates": {}}
    for key, reg in crab.regs.items():
        res["plates"][key] = {"label": reg.plate.label, "epoch": reg.epoch, "epoch_span": reg.epoch_span,
                              "fwhm_arcsec": reg.fwhm_arcsec, "solutions": reg.solutions}
    fits = {}
    for name, (a, b) in PAIRS.items():
        v = crab.vectors(a, b, refresh=refresh)
        res["pairs"][name], fits[name] = _pair_summary(v)
        print(f"  {name}: T = {res['pairs'][name]['T_fit']:.1f} +- {res['pairs'][name]['T_fit_err']:.1f} "
              f"({res['pairs'][name]['n_used']} vectors)")

    v = crab.vectors(*PAIRS[MAIN])
    fit = fits[MAIN]
    used = v.subset(fit.inliers)
    boot = M.block_bootstrap(used, lambda s: (lambda f: (f.T, f.xi_c, f.eta_c))(M.fit_expansion(s)))
    variants = {}
    for name, kw in VARIANTS.items():
        vv = crab.vectors(*PAIRS[MAIN], **kw)
        f = M.fit_expansion(vv)
        variants[name] = {"T": f.T, "centre": [f.xi_c, f.eta_c], "n": f.n}
    var_T = np.array([x["T"] for x in variants.values()])
    t_mid = res["pairs"][MAIN]["t_mid"]
    res["main"] = {
        "pair": MAIN, "T": fit.T, "T_err_formal": float(fit.errors[2]),
        "T_err_bootstrap": float(np.std(boot[:, 0])),
        "T_err_processing": float(np.std(np.r_[var_T, fit.T])),
        "centre": [fit.xi_c, fit.eta_c],
        "centre_err_bootstrap": [float(np.std(boot[:, 1])), float(np.std(boot[:, 2]))],
        "centre_err_processing": [float(np.std([x["centre"][0] for x in variants.values()])),
                                  float(np.std([x["centre"][1] for x in variants.values()]))],
        "variants": variants,
        "expansion_rate_pct_yr": 100.0 / (t_mid - fit.T),
        "velocity_median_kms": res["pairs"][MAIN]["mu_median_mas_yr"] * config.KMS_PER_MASYR_KPC * config.CRAB.distance_kpc,
        "velocity_p95_kms": res["pairs"][MAIN]["mu_p95_mas_yr"] * config.KMS_PER_MASYR_KPC * config.CRAB.distance_kpc,
    }
    m = res["main"]
    m["T_err_total"] = float(np.hypot(m["T_err_bootstrap"], m["T_err_processing"]))
    red = [k for k in PAIRS if k.startswith("red")]
    m["T_red_pairs"] = {k: res["pairs"][k]["T_fit"] for k in red}
    res["wedges"] = wedge_dates(used, fit.xi_c, fit.eta_c)
    res["radial"] = radial_dates(used, fit.xi_c, fit.eta_c)
    res["pulsar"] = pulsar(fit, crab.gaia, t_mid)
    res["acceleration"] = acceleration(fit.T, t_mid)
    # Martin et al.'s centre is relative to Duncan's star; put it in our frame.
    i = int(np.flatnonzero(crab.gaia.source_id == config.DUNCAN_STAR_GAIA)[0])
    dxi, deta = A.to_tangent(crab.gaia.ra[i], crab.gaia.dec[i], config.CRAB.ra, config.CRAB.dec)
    res["martin2025_centre"] = [float(dxi + L.MARTIN2025["centre_from_duncan"][0]),
                                float(deta + L.MARTIN2025["centre_from_duncan"][1])]
    res["duncan_star"] = [float(dxi), float(deta)]
    cache.write_text(json.dumps(res, indent=1, default=float))
    return res
