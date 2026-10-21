"""Cas A in X-rays: how fast the blast wave and the debris expand, 2000-2019.

Chandra has imaged Cas A almost every year since 2000. Fifteen epochs of the
whole remnant are public (``config.CHANDRA_EPOCHS``). The method follows
Vink, Patnaude & Castro (2022), who did this with the same observations:

* The deep 2004 image (eight observations, 0.97 Ms, the "Very Large
  Project") is the **model**. To compare it with another epoch, stretch it
  about the expansion centre by a factor s, shift it by the pointing error
  (dx, dy) of that observation, scale its brightness, and ask how well it
  predicts the counts actually recorded, pixel by pixel, with Poisson
  statistics (the Cash statistic).
* First the whole remnant fixes each epoch's pointing correction and a
  global stretch; then every 20-degree sector of the forward-shock rim and
  of the bright ejecta ring gets its own stretch with the pointing held.
* Stretch against time is a straight line whose slope is the **expansion
  rate** a = (dR/dt)/R. Its inverse is the age the remnant would have if it
  had always expanded this fast; the true age over that is the **expansion
  parameter** m = a t (m = 1 for free expansion, 0.4 for a Sedov-Taylor blast
  wave).

Differences from Vink et al.: they extracted 4.2-6 keV images (synchrotron
light that traces the shock itself) from the event lists; we use the
archive's ready-made all-energy "centre images" (0.7 MB per observation
instead of ~500 MB of events), as DeLaney & Rudnick (2003) did for the first
such measurement. Our forward-shock sectors therefore also see some ejecta.
"""

from __future__ import annotations

import glob
import json
import warnings
from dataclasses import dataclass

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from matplotlib.path import Path as MplPath
from scipy import ndimage, optimize

import astrometry as A
import config

warnings.filterwarnings("ignore", category=FITSFixedWarning)

#: Thorstensen et al. (2001) expansion centre from the fastest optical knots:
#: the centre every stretch is made about.
CENTRE_RA, CENTRE_DEC = config.CASA.ra, config.CASA.dec
#: Half-size of the box cut around Cas A (pixels of 0.492"): 3.9'.
HALF_BOX = 480
#: Thorstensen et al.'s date, used to turn rates into expansion parameters.
EXPLOSION_YEAR = 1671.3

#: Regions, radii in arcmin about the centre, as Vink et al. (2022): the
#: forward-shock annulus of their Fig. 1 and their reverse-shock annulus.
#: In our all-energy images the inner one is dominated by the bright ring of
#: shocked ejecta just outside the reverse shock.
FORWARD_SHOCK = (2.4, 3.0)
EJECTA_RING = (1.67, 2.18)
WHOLE = (0.0, 3.2)
SECTOR_DEG = 20.0


def chandra_year(t_sec: float) -> float:
    """Chandra time (s since 1998-01-01 TT, MJD 50814) as a Julian year."""
    return 2000.0 + (50814.0 + t_sec / 86400.0 - 51544.5) / 365.25


@dataclass
class Epoch:
    label: str
    obsids: tuple
    year: float  # exposure-weighted mid-time
    exposure: float  # s
    counts: np.ndarray  # summed counts on the common pixel grid (box around Cas A)
    footprint: np.ndarray  # True where every observation of this epoch has full exposure


@dataclass
class Frame:
    """The common pixel grid: the archive's 'sky' pixels of the reference
    observation, cut to a box around Cas A. North up, east left, 0.492"."""

    wcs: WCS
    x0: int  # box corner in the reference image
    y0: int
    shape: tuple

    def radec_to_box(self, ra, dec):
        x, y = self.wcs.world_to_pixel_values(ra, dec)
        return np.asarray(x) - self.x0, np.asarray(y) - self.y0


def _files(obsid: int) -> tuple[str, str]:
    img = glob.glob(str(config.CHANDRA_DIR / f"{obsid:05d}_*cntr_img2.fits.gz"))
    fov = glob.glob(str(config.CHANDRA_DIR / f"{obsid:05d}_*fov1.fits.gz"))
    if len(img) != 1 or len(fov) != 1:
        raise FileNotFoundError(f"ObsID {obsid}: run fetch_data.py --only chandra")
    return img[0], fov[0]


def _footprint(fov_path: str, img_wcs: WCS, shape) -> np.ndarray:
    """The S3 chip's outline (from the field-of-view file), rasterised, then
    eroded by 16" so the dither-blurred edge with partial exposure is dropped."""
    with fits.open(fov_path) as h:
        tab = h[1]
        row = tab.data[0]
        n = int(np.sum(np.isfinite(row["X"]) & (row["X"] != 0)))
        cx, cy = row["X"][:n], row["Y"][:n]
        cols = tab.columns
        w = WCS(naxis=2)
        w.wcs.ctype = [cols["X"].coord_type, cols["Y"].coord_type]
        w.wcs.crpix = [cols["X"].coord_ref_point, cols["Y"].coord_ref_point]
        w.wcs.crval = [cols["X"].coord_ref_value, cols["Y"].coord_ref_value]
        w.wcs.cdelt = [cols["X"].coord_inc, cols["Y"].coord_inc]
    ra, dec = w.pixel_to_world_values(cx - 1, cy - 1)  # FITS 1-based -> 0-based
    px, py = img_wcs.world_to_pixel_values(ra, dec)
    yy, xx = np.mgrid[: shape[0], : shape[1]]
    inside = MplPath(np.c_[px, py]).contains_points(np.c_[xx.ravel(), yy.ravel()]).reshape(shape)
    return ndimage.binary_erosion(inside, iterations=int(16 / config.CHANDRA_PIXEL_ARCSEC))


def load(refresh: bool = False) -> tuple[Frame, list[Epoch]]:
    """Every epoch summed onto the reference frame (cached in data/derived)."""
    cache = config.DERIVED_DIR / "chandra_epochs.npz"
    ref_img, _ = _files(config.CHANDRA_EPOCHS[config.CHANDRA_REFERENCE_EPOCH][0])
    ref_wcs = WCS(fits.getheader(ref_img))
    xc, yc = ref_wcs.world_to_pixel_values(CENTRE_RA, CENTRE_DEC)
    x0, y0 = int(round(float(xc))) - HALF_BOX, int(round(float(yc))) - HALF_BOX
    frame = Frame(ref_wcs, x0, y0, (2 * HALF_BOX, 2 * HALF_BOX))
    if cache.exists() and not refresh:
        z = np.load(cache, allow_pickle=True)
        meta = json.loads(str(z["meta"]))
        # Counts are cached as float32 to save space but must be float64 here:
        # the Cash statistic sums ~1e7 over half a million pixels, and float32
        # rounding would swamp the Delta C ~ 1 the fits have to resolve.
        epochs = [Epoch(m["label"], tuple(m["obsids"]), m["year"], m["exposure"],
                        z[f"c{i}"].astype(np.float64), z[f"f{i}"]) for i, m in enumerate(meta)]
        return frame, epochs

    epochs = []
    for label, obsids in config.CHANDRA_EPOCHS.items():
        total = np.zeros(frame.shape)
        foot = np.ones(frame.shape, bool)
        t_weighted, exp_total = 0.0, 0.0
        for obsid in obsids:
            img_path, fov_path = _files(obsid)
            with fits.open(img_path) as h:
                data = h[0].data.astype(float)
                hdr = h[0].header
            w = WCS(hdr)
            fp = _footprint(fov_path, w, data.shape)
            # Where does our box sit in this image? Same nominal pointing ->
            # identical pixel grid; otherwise (ObsID 14482) shift by the
            # offset between the grids (order-1 interpolation, sums conserved).
            x_here, y_here = w.world_to_pixel_values(*ref_wcs.pixel_to_world_values(x0, y0))
            dx, dy = float(x_here) - x0, float(y_here) - y0
            if abs(dx) < 1e-3 and abs(dy) < 1e-3:
                cut = data[y0 : y0 + frame.shape[0], x0 : x0 + frame.shape[1]]
                fcut = fp[y0 : y0 + frame.shape[0], x0 : x0 + frame.shape[1]]
            else:
                yy, xx = np.mgrid[: frame.shape[0], : frame.shape[1]]
                cut = ndimage.map_coordinates(data, [yy + y0 + dy, xx + x0 + dx], order=1, cval=0.0)
                fcut = ndimage.map_coordinates(fp.astype(float), [yy + y0 + dy, xx + x0 + dx], order=0) > 0.5
            if cut.shape != frame.shape:
                raise RuntimeError(f"ObsID {obsid}: box falls off the image")
            total += cut
            foot &= fcut
            mid = chandra_year(0.5 * (hdr["TSTART"] + hdr["TSTOP"]))
            t_weighted += mid * hdr["EXPOSURE"]
            exp_total += hdr["EXPOSURE"]
        epochs.append(Epoch(label, tuple(obsids), t_weighted / exp_total, exp_total, total, foot))

    meta = [{"label": e.label, "obsids": list(e.obsids), "year": e.year, "exposure": e.exposure} for e in epochs]
    np.savez_compressed(cache, meta=json.dumps(meta),
                        **{f"c{i}": e.counts.astype(np.float32) for i, e in enumerate(epochs)},
                        **{f"f{i}": e.footprint for i, e in enumerate(epochs)})
    return frame, epochs


# --------------------------------------------------------------------------
# The stretch model
# --------------------------------------------------------------------------

class StretchModel:
    """The reference image, stretched about the centre and shifted.

    model(x) = A * R(c + (x - c - d) / s)

    R is the reference image smoothed by a 0.5-pixel Gaussian (it has ~20x
    the counts of a 50 ks epoch, so its own noise matters little). A is
    solved for analytically: for Poisson data the best scale is
    sum(counts) / sum(R'), over the pixels being fitted."""

    def __init__(self, reference: np.ndarray, centre_xy: tuple[float, float]):
        self.R = ndimage.gaussian_filter(reference.astype(np.float64), 0.5)
        self.cx, self.cy = centre_xy
        yy, xx = np.mgrid[: reference.shape[0], : reference.shape[1]]
        self.xx, self.yy = xx.astype(float), yy.astype(float)

    def predict(self, s: float, dx: float, dy: float, sel: np.ndarray) -> np.ndarray:
        x = self.cx + (self.xx[sel] - self.cx - dx) / s
        y = self.cy + (self.yy[sel] - self.cy - dy) / s
        return ndimage.map_coordinates(self.R, [y, x], order=3, mode="nearest")

    def cash(self, counts: np.ndarray, sel: np.ndarray, s: float, dx: float, dy: float) -> float:
        """Cash (1979) statistic, C = 2 sum(m - n ln m), with A at its optimum."""
        r = np.clip(self.predict(s, dx, dy, sel), 1e-6, None)
        n = counts[sel].astype(np.float64)
        a = n.sum() / r.sum()
        m = a * r
        return float(2.0 * np.sum(m - n * np.log(m)))


def region_mask(frame: Frame, inner_arcmin: float, outer_arcmin: float, pa_range=None) -> np.ndarray:
    """Annulus (and optionally a position-angle sector, deg east of north) about the centre."""
    cx, cy = frame.radec_to_box(CENTRE_RA, CENTRE_DEC)
    yy, xx = np.mgrid[: frame.shape[0], : frame.shape[1]]
    dx, dy = xx - float(cx), yy - float(cy)
    r = np.hypot(dx, dy) * config.CHANDRA_PIXEL_ARCSEC / 60.0
    sel = (r >= inner_arcmin) & (r < outer_arcmin)
    if pa_range is not None:
        pa = np.degrees(np.arctan2(-dx, dy)) % 360.0  # east is -x in the sky frame
        lo, hi = pa_range
        sel &= (pa >= lo) & (pa < hi)
    return sel


def fit_global(model: StretchModel, ep: Epoch, sel: np.ndarray) -> dict:
    """Stretch and pointing shift of one epoch against the reference, whole remnant."""
    use = sel & ep.footprint

    def f(p):
        return model.cash(ep.counts, use, 1.0 + p[0] * 1e-3, p[1], p[2])

    best = optimize.minimize(f, [0.0, 0.0, 0.0], method="Nelder-Mead",
                             options={"xatol": 1e-4, "fatol": 0.05, "maxiter": 2000})
    p = best.x
    # Curvature errors: Delta C = 1 at one sigma, from a numerical Hessian.
    h = np.array([0.02, 0.01, 0.01])
    H = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ei, ej = np.eye(3)[i] * h[i], np.eye(3)[j] * h[j]
            H[i, j] = (f(p + ei + ej) - f(p + ei - ej) - f(p - ei + ej) + f(p - ei - ej)) / (4 * h[i] * h[j])
    cov = np.linalg.inv(H / 2.0)
    return {"s": 1.0 + p[0] * 1e-3, "s_err": np.sqrt(cov[0, 0]) * 1e-3, "dx": p[1], "dy": p[2],
            "dx_err": np.sqrt(cov[1, 1]), "dy_err": np.sqrt(cov[2, 2]), "n_pix": int(use.sum()),
            "counts": float(ep.counts[use].sum())}


def fit_stretch(model: StretchModel, ep: Epoch, sel: np.ndarray, dx: float, dy: float) -> tuple[float, float]:
    """Stretch of one region of one epoch with the pointing held fixed:
    a one-dimensional minimisation, with the error where C rises by 1."""
    use = sel & ep.footprint
    if use.sum() < 200:
        return np.nan, np.nan
    f = lambda e: model.cash(ep.counts, use, 1.0 + e * 1e-3, dx, dy)
    res = optimize.minimize_scalar(f, bounds=(-40.0, 60.0), method="bounded", options={"xatol": 1e-3})
    e0 = res.x
    h = 0.05
    curv = (f(e0 + h) - 2 * f(e0) + f(e0 - h)) / h**2
    err = np.sqrt(2.0 / curv) if curv > 0 else np.nan
    return 1.0 + e0 * 1e-3, err * 1e-3


def rate_fit(t: np.ndarray, s: np.ndarray, s_err: np.ndarray, t_ref: float, *, quadratic: bool = False) -> dict:
    """Straight line (or parabola) through stretch against time, forced
    through s = 1 at the reference epoch. Errors are inflated by the scatter
    about the line when it exceeds the statistical errors (it does: pointing
    and calibration differences between epochs are not in the Cash errors)."""
    ok = np.isfinite(s) & np.isfinite(s_err) & (s_err > 0)
    dt = t[ok] - t_ref
    y = s[ok] - 1.0
    w = 1.0 / s_err[ok] ** 2
    X = np.c_[dt, 0.5 * dt**2] if quadratic else dt[:, None]
    W = np.diag(w)
    cov = np.linalg.inv(X.T @ W @ X)
    beta = cov @ X.T @ W @ y
    resid = y - X @ beta
    chi2 = float(np.sum(w * resid**2))
    dof = max(1, len(y) - X.shape[1])
    scale = max(1.0, chi2 / dof)
    err = np.sqrt(np.diag(cov) * scale)
    out = {"a": float(beta[0]), "a_err": float(err[0]), "chi2_dof": chi2 / dof, "n": int(ok.sum())}
    if quadratic:
        out.update({"b": float(beta[1]), "b_err": float(err[1])})
    return out


def run(refresh: bool = False) -> dict:
    """All fits; returns a JSON-able dict (cached)."""
    cache = config.DERIVED_DIR / "chandra_fits.json"
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())
    frame, epochs = load()
    ref = next(e for e in epochs if e.label == config.CHANDRA_REFERENCE_EPOCH)
    cx, cy = frame.radec_to_box(CENTRE_RA, CENTRE_DEC)
    model = StretchModel(ref.counts, (float(cx), float(cy)))
    whole = region_mask(frame, *WHOLE)

    out = {"reference": ref.label, "t_ref": ref.year, "epochs": []}
    for ep in epochs:
        if ep is ref:
            g = {"s": 1.0, "s_err": 0.0, "dx": 0.0, "dy": 0.0, "dx_err": 0.0, "dy_err": 0.0}
        else:
            g = fit_global(model, ep, whole)
        rec = {"label": ep.label, "obsids": list(ep.obsids), "year": ep.year, "exposure_ks": ep.exposure / 1e3,
               "global": g, "sectors": {}}
        if ep is not ref:
            for name, (r0, r1) in (("forward_shock", FORWARD_SHOCK), ("ejecta_ring", EJECTA_RING)):
                full = fit_stretch(model, ep, region_mask(frame, r0, r1), g["dx"], g["dy"])
                rec["sectors"][f"{name}:all"] = full
                for lo in np.arange(0.0, 360.0, SECTOR_DEG):
                    sel = region_mask(frame, r0, r1, (lo, lo + SECTOR_DEG))
                    rec["sectors"][f"{name}:{lo + SECTOR_DEG / 2:.0f}"] = fit_stretch(model, ep, sel, g["dx"], g["dy"])
        out["epochs"].append(rec)
        print(f"  {ep.label}: s = {g['s']:.5f} +- {g['s_err']:.5f}, pointing ({g['dx']:+.2f}, {g['dy']:+.2f}) px",
              flush=True)
    cache.write_text(json.dumps(out, indent=1, default=float))
    return out


# --------------------------------------------------------------------------
# The neutron star
# --------------------------------------------------------------------------

#: The central compact object, CXOU J232327.9+584842 (Fesen et al. 2006):
#: a starting guess only; every epoch is centroided.
CCO_RA, CCO_DEC = 350.866429, 58.811808


def cco_centroid(counts: np.ndarray, x0: float, y0: float, *, sigma: float = 1.2, iterations: int = 6):
    """Gaussian-weighted centroid of the point source near (x0, y0), after
    subtracting the remnant's local surface brightness (median of a 4-7 px
    annulus). Returns x, y (box pixels) and the net counts used."""
    x, y = x0, y0
    r_out = 7
    for _ in range(iterations):
        ix, iy = int(round(x)), int(round(y))
        st = counts[iy - r_out : iy + r_out + 1, ix - r_out : ix + r_out + 1].astype(float)
        yy, xx = np.mgrid[-r_out : r_out + 1, -r_out : r_out + 1]
        rr = np.hypot(xx, yy)
        st = st - np.median(st[(rr >= 4) & (rr <= 7)])
        w = np.exp(-0.5 * ((xx - (x - ix)) ** 2 + (yy - (y - iy)) ** 2) / sigma**2) * np.clip(st, 0, None)
        x = ix + (w * xx).sum() / w.sum()
        y = iy + (w * yy).sum() / w.sum()
    net = float(np.clip(st, 0, None)[rr <= 3].sum())
    return x, y, net


def neutron_star(fits_result: dict | None = None) -> dict:
    """Positions of the neutron star at every epoch and its proper motion.

    Two frames: (1) the archive's own astrometry, as Holland-Ashford et al.'s
    method 1; (2) corrected by each epoch's pointing shift from the stretch
    fit, i.e. relative to the expanding remnant. The second removes the
    ~0.3-0.6" pointing scatter of Chandra, but any lopsided expansion that
    the shift absorbed comes along with it."""
    fits_result = fits_result or run()
    frame, epochs = load()
    gx, gy = frame.radec_to_box(CCO_RA, CCO_DEC)
    by_label = {e["label"]: e for e in fits_result["epochs"]}
    rows = []
    for ep in epochs:
        x, y, net = cco_centroid(ep.counts, float(gx), float(gy))
        g = by_label[ep.label]["global"]
        rows.append({"label": ep.label, "year": ep.year, "x": x, "y": y, "net_counts": net,
                     "x_corr": x - g["dx"], "y_corr": y - g["dy"]})
    t = np.array([r["year"] for r in rows])
    pix = config.CHANDRA_PIXEL_ARCSEC * 1000.0  # mas per pixel
    out = {"epochs": rows}
    for tag, kx, ky in (("archive", "x", "y"), ("corrected", "x_corr", "y_corr")):
        X = np.array([r[kx] for r in rows])
        Y = np.array([r[ky] for r in rows])
        A_ = np.c_[np.ones_like(t), t - t.mean()]
        cx, *_ = np.linalg.lstsq(A_, X, rcond=None)
        cy, *_ = np.linalg.lstsq(A_, Y, rcond=None)
        rx, ry = X - A_ @ cx, Y - A_ @ cy
        dof = len(t) - 2
        sx = np.sqrt(np.sum(rx**2) / dof) * np.sqrt(np.linalg.inv(A_.T @ A_)[1, 1])
        sy = np.sqrt(np.sum(ry**2) / dof) * np.sqrt(np.linalg.inv(A_.T @ A_)[1, 1])
        mu_e, mu_n = -cx[1] * pix, cy[1] * pix  # east is -x on the sky grid
        se, sn = sx * pix, sy * pix
        mu = np.hypot(mu_e, mu_n)
        mu_err = np.hypot(mu_e * se, mu_n * sn) / mu
        d = config.CASA.distance_kpc
        out[tag] = {"mu_east_mas_yr": float(mu_e), "mu_north_mas_yr": float(mu_n),
                    "mu_east_err": float(se), "mu_north_err": float(sn),
                    "mu_mas_yr": float(mu), "mu_err": float(mu_err),
                    "pa_deg": float(np.degrees(np.arctan2(mu_e, mu_n)) % 360.0),
                    "pa_err_deg": float(np.degrees(mu_err / mu)),
                    "v_kms": float(mu * config.KMS_PER_MASYR_KPC * d),
                    "v_err_kms": float(mu_err * config.KMS_PER_MASYR_KPC * d),
                    "rms_residual_mas": float(np.sqrt(np.mean(rx**2 + ry**2)) * pix)}
    # Offset from the optical expansion centre, at the reference epoch.
    ref = next(r for r in rows if r["label"] == config.CHANDRA_REFERENCE_EPOCH)
    ra, dec = frame.wcs.pixel_to_world_values(ref["x"] + frame.x0, ref["y"] + frame.y0)
    xi, eta = A.to_tangent(float(ra), float(dec), CENTRE_RA, CENTRE_DEC)
    age = ref["year"] - EXPLOSION_YEAR
    off = float(np.hypot(xi, eta))
    out["offset"] = {"east_arcsec": float(xi), "north_arcsec": float(eta), "arcsec": off,
                     "pa_deg": float(np.degrees(np.arctan2(xi, eta)) % 360.0), "age_yr": age,
                     "implied_mu_mas_yr": off / age * 1000.0,
                     "implied_v_kms": off / age * 1000.0 * config.KMS_PER_MASYR_KPC * config.CASA.distance_kpc,
                     "ra": float(ra), "dec": float(dec)}
    return out


# --------------------------------------------------------------------------
# Summaries: rates, expansion parameters, speeds
# --------------------------------------------------------------------------

def shock_radius(reference: np.ndarray, frame: Frame, pa_lo: float, pa_hi: float) -> float:
    """Forward-shock radius (arcmin) in one sector of the deep reference image:
    where the radial profile, beyond its outermost bright rim, falls half-way
    from the rim to the background outside."""
    cx, cy = frame.radec_to_box(CENTRE_RA, CENTRE_DEC)
    yy, xx = np.mgrid[: frame.shape[0], : frame.shape[1]]
    dx, dy = xx - float(cx), yy - float(cy)
    r = np.hypot(dx, dy) * config.CHANDRA_PIXEL_ARCSEC / 60.0
    pa = np.degrees(np.arctan2(-dx, dy)) % 360.0
    sel = (pa >= pa_lo) & (pa < pa_hi)
    edges = np.arange(2.0, 3.6, 0.02)
    prof = np.array([np.mean(reference[sel & (r >= a) & (r < b)]) for a, b in zip(edges[:-1], edges[1:])])
    mid = 0.5 * (edges[:-1] + edges[1:])
    prof = ndimage.uniform_filter1d(prof, 3)
    bg = np.median(prof[mid > 3.3])
    rim = np.argmax(np.where(mid < 3.2, prof, -np.inf))
    half = bg + 0.5 * (prof[rim] - bg)
    beyond = np.flatnonzero((np.arange(len(prof)) > rim) & (prof < half))
    return float(mid[beyond[0]]) if len(beyond) else np.nan


def summarise(fits_result: dict | None = None) -> dict:
    """Expansion rate, timescale, expansion parameter and speed per region."""
    fits_result = fits_result or run()
    frame, epochs = load()
    ref = next(e for e in epochs if e.label == config.CHANDRA_REFERENCE_EPOCH)
    t_ref = fits_result["t_ref"]
    age = t_ref - EXPLOSION_YEAR
    eps = [e for e in fits_result["epochs"] if e["label"] != fits_result["reference"]]
    t = np.array([e["year"] for e in eps])
    kpc_km = 3.0857e16
    arcmin_rad = np.pi / 180.0 / 60.0
    out = {"t_ref": t_ref, "age_at_ref": age, "regions": {}}

    def one(s, s_err, radius_arcmin):
        lin = rate_fit(t, s, s_err, t_ref)
        quad = rate_fit(t, s, s_err, t_ref, quadratic=True)
        a = lin["a"]  # per year
        rec = {"rate_pct_yr": 100 * a, "rate_err": 100 * lin["a_err"], "chi2_dof": lin["chi2_dof"],
               "tau_exp_yr": 1.0 / a if a > 0 else np.nan, "m": a * age,
               "m_err": lin["a_err"] * age, "radius_arcmin": radius_arcmin,
               "v_kms": a * radius_arcmin * arcmin_rad * config.CASA.distance_kpc * kpc_km / 3.15576e7,
               "b_per_yr2": quad.get("b"), "b_err": quad.get("b_err"),
               "t": t.tolist(), "s": list(map(float, s)), "s_err": list(map(float, s_err))}
        rec["v_err_kms"] = rec["v_kms"] * lin["a_err"] / abs(a) if a else np.nan
        return rec

    g = np.array([e["global"]["s"] for e in eps])
    ge = np.array([e["global"]["s_err"] for e in eps])
    out["regions"]["whole"] = one(g, ge, np.nan)
    for name in ("forward_shock", "ejecta_ring"):
        keys = [k.split(":")[1] for k in eps[0]["sectors"] if k.startswith(name + ":")]
        out["regions"][name] = {}
        for k in keys:
            s = np.array([e["sectors"][f"{name}:{k}"][0] for e in eps])
            se = np.array([e["sectors"][f"{name}:{k}"][1] for e in eps])
            if k == "all":
                radius = 2.8 if name == "forward_shock" else 1.9
            else:
                pa = float(k)
                radius = (shock_radius(ref.counts, frame, pa - SECTOR_DEG / 2, pa + SECTOR_DEG / 2)
                          if name == "forward_shock" else 1.9)
            out["regions"][name][k] = one(s, se, radius)
    return out
