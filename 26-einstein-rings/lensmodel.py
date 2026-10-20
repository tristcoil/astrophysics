"""Our own lens models of the Hubble images.

For each grade-A lens we fit, to the F814W image, three things at once:

* **the lens galaxy's light**: a multi-Gaussian expansion (14 concentric
  elliptical Gaussians sharing a centre, axis ratio and orientation, with free
  amplitudes), which can follow any smooth elliptical light profile;
* **the lens's mass**: a singular isothermal ellipsoid (SIE) centred on the
  light, with free Einstein radius b, axis ratio and orientation -- the model
  SLACS used, so that our Einstein radii can be compared with theirs;
* **the background galaxy**: an elliptical Sersic profile in the *source*
  plane. Each image pixel is traced back through the lens (beta = theta -
  alpha(theta)) to where its light left the source.

The model image is blurred by Hubble's point-spread function and compared with
the data pixel by pixel. Amplitudes enter linearly and are solved exactly at
every step; the 9 shape parameters are found by least squares from several
starting points.

Nothing here reads the published lens models. The starting Einstein radius
comes from the image itself: the radius at which the arcs, traced back to the
source plane, converge most tightly onto one spot.
"""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass, field

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage, optimize
from scipy.signal import fftconvolve

import lensing as L
from config import (HST_DIR, HST_PIXEL_ARCSEC, MODEL_BOX_PIXELS, MODEL_DIR, PSF_FWHM_ARCSEC,
                    PSF_MOFFAT_BETA, SUPERSAMPLE)

FILTER = "acs_wfc_f814w"
MGE_SIGMAS = np.geomspace(0.012, 3.0, 14)  # arcsec
FIT_RADIUS = 2.9  # arcsec: the circle we fit inside the 6.4" box
CORE_MASK = 0.25  # arcsec: the lens galaxy's centre, where light models are poorest
#                    (SLACS rings are all >= 0.69"; nothing lensed is lost)
B_MIN = 0.4  # arcsec: smallest Einstein radius allowed, safely outside the core mask
MGE_SPLIT = 0.35  # arcsec: Gaussians narrower than this share one shape (the bulge), wider ones another
OUTLIER_SNR = 3.0  # an outlier sticks out by 3 sigma (3x3 sum) ...
OUTLIER_REL = 0.03  # ... and by at least 3% of the lens light there
OUTLIER_MAX_FRAC = 0.35  # never mask more than this share of the fitted area as outliers
ARC_ZONE = (0.35, 2.3)  # arcsec: where arcs are looked for (SLACS rings are 0.7-1.8")
MIN_ARC_PIXELS = 8  # smallest connected patch counted as lensed light
MASS_E_MAX = 0.48  # mass axis ratio q >= 0.35 (published SLACS values: q >= 0.51)
LIGHT_E_MAX = 0.65  # light axis ratio q >= 0.21 (an edge-on disc is about 0.2)
SOURCE_R_MAX = 0.5  # arcsec: lensed galaxies at z ~ 0.6 are a few kpc across


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------

@dataclass
class Cutout:
    name: str
    sci: np.ndarray  # electrons/s
    var: np.ndarray  # noise variance per pixel (background + Poisson)
    x: np.ndarray  # arcsec, increasing to the right (west)
    y: np.ndarray  # arcsec, increasing upward (north)
    north_angle: float  # direction of north in image coordinates (radians from +x)
    parity: int  # +1 if east is counter-clockwise from north
    sigma_bkg: float
    exptime: float
    centre_pix: tuple[int, int]


def load_cutout(name: str, ra: float, dec: float, box: int = MODEL_BOX_PIXELS,
                filt: str = FILTER) -> Cutout:
    path = HST_DIR / name / f"{filt}.fits"
    with fits.open(path) as hdul, warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sci = hdul["SCI"].data.astype(float)
        wht = hdul["WHT"].data.astype(float)
        wcs = WCS(hdul["SCI"].header)
        exptime = float(hdul[0].header.get("EXPTIME", np.nan))
    bad = ~np.isfinite(sci) | ~np.isfinite(wht) | (wht <= 0)
    sci[bad] = 0.0
    wht[bad] = 0.0

    # Centre on the lens galaxy: the light peak near the catalogue position.
    tx, ty = wcs.world_to_pixel_values(ra, dec)
    tx, ty = int(round(float(tx))), int(round(float(ty)))
    win = ndimage.gaussian_filter(sci, 1.0)[ty - 6:ty + 7, tx - 6:tx + 7]
    dy, dx = np.unravel_index(np.argmax(win), win.shape)
    cx, cy = tx - 6 + dx, ty - 6 + dy
    h = box // 2
    if cy - h < 0 or cx - h < 0 or cy + h > sci.shape[0] or cx + h > sci.shape[1]:
        raise ValueError(f"{name}: lens too close to the cutout edge")
    s = sci[cy - h:cy + h, cx - h:cx + h]
    w = wht[cy - h:cy + h, cx - h:cx + h]

    # Directions of north and east on the image, from the WCS.
    p0 = np.array(wcs.world_to_pixel_values(ra, dec))
    pn = np.array(wcs.world_to_pixel_values(ra, dec + 1 / 3600))
    pe = np.array(wcs.world_to_pixel_values(ra + 1 / 3600 / np.cos(np.radians(dec)), dec))
    north = np.arctan2(*(pn - p0)[::-1])
    east = np.arctan2(*(pe - p0)[::-1])
    parity = 1 if np.sin(east - north) > 0 else -1

    # Background noise from the full cutout's pixels well away from the galaxy,
    # after removing the slowly varying light of the galaxy's outskirts.
    yy, xx = np.mgrid[:sci.shape[0], :sci.shape[1]]
    far = (np.hypot(xx - cx, yy - cy) > 60) & (wht > 0)
    smooth = ndimage.median_filter(sci, size=15)
    resid = (sci - smooth)[far]
    sigma_bkg = 1.4826 * np.median(np.abs(resid - np.median(resid)))

    idx = np.arange(box) - h  # pixel centres; the galaxy peak is at index h
    x = idx[None, :] * HST_PIXEL_ARCSEC * np.ones((box, 1))
    y = idx[:, None] * HST_PIXEL_ARCSEC * np.ones((1, box))
    # Poisson noise of the source photons: WHT is the effective exposure time.
    var = sigma_bkg**2 + np.clip(s, 0, None) / np.maximum(w, 1.0)
    var[w <= 0] = np.inf
    return Cutout(name, s, var, x, y, float(north), parity, float(sigma_bkg), exptime, (cx, cy))


def pa_to_angle(pa_deg: float, cut: Cutout) -> float:
    """Sky position angle (degrees east of north) -> image angle (radians from +x)."""
    return cut.north_angle + cut.parity * np.radians(pa_deg)


def angle_to_pa(angle: float, cut: Cutout) -> float:
    """Image angle -> sky position angle in [0, 180) degrees (an axis, not a direction)."""
    return float(np.degrees((angle - cut.north_angle) * cut.parity) % 180.0)


# --------------------------------------------------------------------------
# Model pieces
# --------------------------------------------------------------------------

def moffat_psf(fwhm=PSF_FWHM_ARCSEC, beta=PSF_MOFFAT_BETA, pix=HST_PIXEL_ARCSEC, half=12):
    alpha = fwhm / (2 * np.sqrt(2 ** (1 / beta) - 1))
    g = np.arange(-half, half + 1) * pix
    r2 = g[None, :] ** 2 + g[:, None] ** 2
    k = (1 + r2 / alpha**2) ** (-beta)
    return k / k.sum()


def _supersample(cut: Cutout, n=SUPERSAMPLE):
    """Sub-pixel coordinates: each pixel split into n x n rays."""
    off = (np.arange(n) + 0.5) / n - 0.5
    ox, oy = np.meshgrid(off, off)
    xs = cut.x[..., None] + ox.ravel() * HST_PIXEL_ARCSEC
    ys = cut.y[..., None] + oy.ravel() * HST_PIXEL_ARCSEC
    return xs, ys  # shape (box, box, n*n)


def e_to_qphi(e1, e2, e_max=0.95):
    e = min(np.hypot(e1, e2), e_max)
    return (1 - e) / (1 + e), 0.5 * np.arctan2(e2, e1)


def qphi_to_e(q, phi):
    e = (1 - q) / (1 + q)
    return e * np.cos(2 * phi), e * np.sin(2 * phi)


def sersic_source(bx, by, xs, ys, r_s, n_s, q_s, phi_s):
    xr, yr = L.rotate(bx - xs, by - ys, phi_s)
    r = np.sqrt(q_s * xr**2 + yr**2 / q_s)
    return np.exp(-L.sersic_b(n_s) * ((r / r_s) ** (1 / n_s) - 1))


def cosmic_rays(cut: Cutout, nsigma: float = 6.0):
    """Pixels much sharper than Hubble's point-spread function: cosmic-ray hits,
    which single 420-s snapshot exposures cannot reject. A pixel is flagged if
    it stands > 6 sigma above the median of its 3x3 neighbourhood and that
    excess is most of its height above the wider (7x7) surroundings -- true of
    a one-pixel hit, not of a star or an arc, whose light is spread over the PSF."""
    s = cut.sci
    sig = np.sqrt(np.where(np.isfinite(cut.var), cut.var, np.inf))
    d3 = s - ndimage.median_filter(s, size=3)
    d7 = s - ndimage.median_filter(s, size=7)
    hit = (d3 / sig > nsigma) & (d3 > 0.6 * np.maximum(d7, 1e-30))
    return ndimage.binary_dilation(hit, iterations=1)


def manual_mask(cut: Cutout):
    """Neighbouring galaxies masked by hand (config.MANUAL_MASKS): circles given
    as (arcsec east, arcsec north, radius) from the lens centre."""
    from config import MANUAL_MASKS

    m = np.zeros(cut.sci.shape, dtype=bool)
    for de, dn, rad in MANUAL_MASKS.get(cut.name, []):
        east = cut.north_angle + cut.parity * np.pi / 2
        x = dn * np.cos(cut.north_angle) + de * np.cos(east)
        y = dn * np.sin(cut.north_angle) + de * np.sin(east)
        m |= np.hypot(cut.x - x, cut.y - y) < rad
    return m


class LensFitter:
    """Holds one cutout and fits the light + SIE + source model to it."""

    def __init__(self, cut: Cutout, psf=None, hand_masks: bool = False):
        self.cut = cut
        self.psf = moffat_psf() if psf is None else psf
        self.xs, self.ys = _supersample(cut)
        r = np.hypot(cut.x, cut.y)
        self.r = r
        self.sharp = cosmic_rays(cut)
        self.manual = manual_mask(cut) if hand_masks else np.zeros(cut.sci.shape, dtype=bool)
        self.base_mask = ((r < FIT_RADIUS) & (r > CORE_MASK) & np.isfinite(cut.var)
                          & ~self.sharp & ~self.manual)
        self.mask = self.base_mask.copy()
        self.light_shape = None  # (x0, y0, q, phi)
        self.light_basis = None

    # ---- lens light --------------------------------------------------------
    def _mge_basis(self, x0, y0, shapes):
        """Concentric elliptical Gaussians. ``shapes`` = ((q_in, phi_in), (q_out, phi_out)):
        Gaussians narrower than MGE_SPLIT take the inner shape, wider ones the outer,
        so a round bulge inside a flat disc (or isophotes that twist) can be followed."""
        cols = []
        rr2 = []
        for q, phi in shapes:
            xr, yr = L.rotate(self.xs - x0, self.ys - y0, phi)
            rr2.append(q * xr**2 + yr**2 / q)
        for s in MGE_SIGMAS:
            r2 = rr2[0] if s < MGE_SPLIT else rr2[1]
            img = np.exp(-0.5 * r2 / s**2).mean(axis=-1)
            img = fftconvolve(img, self.psf, mode="same")
            cols.append(img / max(img.sum(), 1e-30))
        return np.array(cols)  # (K, box, box)

    def _solve_linear(self, columns, mask):
        """Weighted linear least squares for the amplitudes; returns (model, amps, chi2)."""
        w = 1.0 / self.cut.var[mask]
        a = columns[:, mask].T * np.sqrt(w)[:, None]
        b = self.cut.sci[mask] * np.sqrt(w)
        amps, *_ = np.linalg.lstsq(a, b, rcond=None)
        model = np.tensordot(amps, columns, axes=1)
        chi2 = float(np.sum((self.cut.sci[mask] - model[mask]) ** 2 * w))
        return model, amps, chi2

    def fit_light(self, n_iter=3):
        """Fit the lens galaxy's light alone, masking anything that sticks out
        above it (arcs, companions) and refitting."""
        ones = np.ones_like(self.cut.sci)[None]

        def shapes_of(p):
            return (e_to_qphi(p[2], p[3], LIGHT_E_MAX), e_to_qphi(p[4], p[5], LIGHT_E_MAX))

        def cost(p):
            cols = np.concatenate([self._mge_basis(p[0], p[1], shapes_of(p)), ones])
            return self._solve_linear(cols, self.mask)[2]

        p = np.array([0.0, 0.0, 0.05, 0.0, 0.05, 0.0])
        for _ in range(n_iter):
            res = optimize.minimize(cost, p, method="Nelder-Mead",
                                    options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 1200})
            p = res.x
            cols = np.concatenate([self._mge_basis(p[0], p[1], shapes_of(p)), ones])
            model, amps, _ = self._solve_linear(cols, self.mask)
            outlier = self._outliers(model)
            self.mask = self.base_mask & ~outlier
        x0, y0 = p[0], p[1]
        (q_in, phi_in), (q, phi) = shapes_of(p)
        self.light_inner = (q_in, phi_in)
        self.light_shape = (x0, y0, q, phi)
        self.light_basis = cols
        self.light_model, self.light_amps, _ = self._solve_linear(cols, self.mask)
        self.light_residual = self.cut.sci - self.light_model
        self.outlier_mask = outlier
        return self.light_shape

    def _outliers(self, model):
        """Pixels that stick out above the lens-light model: by OUTLIER_SNR in a
        3x3 sum *and* by OUTLIER_REL of the model, so imperfections near a bright
        centre are not mistaken for arcs. If that would mask too much of the
        galaxy, the S/N bar is raised until it does not."""
        resid = self.cut.sci - model
        snr = ndimage.uniform_filter(np.nan_to_num(resid / np.sqrt(self.cut.var)), 3) * 3
        rel = ndimage.uniform_filter(resid, 3) > OUTLIER_REL * np.maximum(model, 0)
        area = self.base_mask.sum()
        thr = OUTLIER_SNR
        while True:
            out = ndimage.binary_dilation((snr > thr) & rel, iterations=2)
            if (out & self.base_mask).sum() <= OUTLIER_MAX_FRAC * area or thr > 50:
                return out
            thr *= 1.5

    # ---- first guess from the arcs ----------------------------------------
    def arc_pixels(self):
        """Pixels significantly brighter than the lens light: the arcs (and any
        companions, which the radius limits and later masking deal with)."""
        snr = ndimage.uniform_filter(np.nan_to_num(self.light_residual / np.sqrt(self.cut.var)), 3) * 3
        rel = (ndimage.uniform_filter(self.light_residual, 3)
               > OUTLIER_REL * np.maximum(self.light_model, 0))
        sel = ((snr > 3.0) & rel & (self.r > ARC_ZONE[0]) & (self.r < ARC_ZONE[1])
               & ~self.sharp & ~self.manual)
        labels, n = ndimage.label(sel)
        if n:
            sizes = ndimage.sum(sel, labels, index=np.arange(1, n + 1))
            keep = np.isin(labels, 1 + np.flatnonzero(sizes >= MIN_ARC_PIXELS))
            sel &= keep
        return sel

    def initial_guesses(self, n_max: int = 3):
        """Einstein radii at which the arc pixels, traced back through an SIE with
        the light's shape, land closest together in the source plane: the
        deepest local minima of that spread, best first."""
        x0, y0, q_l, phi_l = self.light_shape
        sel = self.arc_pixels()
        if sel.sum() < 10:
            raise RuntimeError("no arcs detected above the lens light")
        f = np.clip(self.light_residual[sel], 0, None)
        px, py = self.cut.x[sel], self.cut.y[sel]
        q_m = float(np.clip(q_l, 0.4, 0.95))
        grid = np.arange(B_MIN, 2.31, 0.01)
        scores, centres = [], []
        for b in grid:
            ax, ay = L.sie_deflection_sky(px, py, b, q_m, phi_l, x0, y0)
            bx, by = px - ax, py - ay
            mx, my = np.average(bx, weights=f), np.average(by, weights=f)
            spread = np.average((bx - mx) ** 2 + (by - my) ** 2, weights=f)
            # spread is compared with the arcs' own extent, so a lens that
            # simply shrinks everything to a point is not preferred
            scores.append(spread / np.average((px - x0) ** 2 + (py - y0) ** 2, weights=f))
            centres.append((mx, my))
        scores = np.array(scores)
        is_min = np.r_[True, scores[1:-1] <= np.minimum(scores[:-2], scores[2:]), True]
        order = [i for i in np.argsort(scores) if is_min[i]]
        picks = []
        for i in order:
            if all(abs(grid[i] - grid[j]) > 0.15 for j in picks):
                picks.append(i)
            if len(picks) == n_max:
                break
        self.guess_curve = (grid, scores)
        return [{"b": float(grid[i]), "q_m": q_m, "phi_m": float(phi_l),
                 "xs": float(centres[i][0]), "ys": float(centres[i][1]),
                 "score": float(scores[i])} for i in picks]

    # ---- the full model -----------------------------------------------------
    def lensed_source(self, b, q_m, phi_m, xs, ys, r_s, n_s, q_s, phi_s, shear=None):
        x0, y0, _, _ = self.light_shape
        ax, ay = L.sie_deflection_sky(self.xs, self.ys, b, q_m, phi_m, x0, y0)
        if shear is not None:
            gx, gy = L.shear_deflection(self.xs, self.ys, shear[0], shear[1], x0, y0)
            ax, ay = ax + gx, ay + gy
        img = sersic_source(self.xs - ax, self.ys - ay, xs, ys, r_s, n_s, q_s, phi_s).mean(axis=-1)
        return fftconvolve(img, self.psf, mode="same")

    def unpack(self, p):
        b, e1m, e2m, xs, ys, log_rs, n_s, e1s, e2s = p[:9]
        q_m, phi_m = e_to_qphi(e1m, e2m, MASS_E_MAX)
        q_s, phi_s = e_to_qphi(e1s, e2s)
        out = dict(b=b, q_m=q_m, phi_m=phi_m, xs=xs, ys=ys, r_s=10**log_rs, n_s=n_s,
                   q_s=q_s, phi_s=phi_s)
        if len(p) > 9:
            g1, g2 = p[9:11]
            out["shear"] = (float(np.hypot(g1, g2)), 0.5 * float(np.arctan2(g2, g1)))
        return out

    def model_image(self, p):
        par = self.unpack(p)
        src = self.lensed_source(**par)
        cols = np.concatenate([self.light_basis[:-1], src[None], self.light_basis[-1:]])
        model, amps, chi2 = self._solve_linear(cols, self.fit_mask)
        return model, amps, chi2, src

    def residuals(self, p):
        model, amps, _, _ = self.model_image(p)
        r = (self.cut.sci - model)[self.fit_mask] / np.sqrt(self.cut.var[self.fit_mask])
        if amps[-2] < 0:  # a negative source is not a source
            r = r * (1 + abs(amps[-2]))
        return r

    def fit(self, with_shear: bool = False, starts_per_guess: int = 3, seed: int = 0):
        """Fit from each candidate Einstein radius and keep the best chi^2.

        All candidates are judged on the same pixels: everything inside the
        arcs' reach, minus cosmic rays, hand-masked neighbours, and anything
        that sticks out above the lens light beyond the arc zone."""
        if self.light_shape is None:
            self.fit_light()
        guesses = self.initial_guesses()
        reach = min(FIT_RADIUS, 1.7 * max(g["b"] for g in guesses) + 0.5)
        far_outliers = self.outlier_mask & (self.r > ARC_ZONE[1])
        self.fit_mask = self.base_mask & (self.r < reach) & ~far_outliers
        rng = np.random.default_rng(seed)
        em = MASS_E_MAX / np.sqrt(2)  # box bounds inside the |e| cap
        lo = [B_MIN, -em, -em, -1.5, -1.5, np.log10(0.01), 0.3, -0.8, -0.8]
        hi = [3.0, em, em, 1.5, 1.5, np.log10(SOURCE_R_MAX), 6.0, 0.8, 0.8]
        if with_shear:
            lo += [-0.3, -0.3]
            hi += [0.3, 0.3]
        best, best_guess, tried = None, None, []
        for g in guesses:
            e1m, e2m = qphi_to_e(g["q_m"], g["phi_m"])
            best_here = None
            for k in range(starts_per_guess):
                jitter = 0 if k == 0 else 1
                p0 = [g["b"] * (1 + 0.04 * jitter * rng.normal()), e1m, e2m,
                      g["xs"] + 0.06 * jitter * rng.normal(), g["ys"] + 0.06 * jitter * rng.normal(),
                      np.log10([0.08, 0.2, 0.04][k % 3]), 1.0, 0.0, 0.0]
                if with_shear:
                    p0 += [0.0, 0.0]
                p0 = np.clip(p0, np.array(lo) + 1e-6, np.array(hi) - 1e-6)
                try:
                    res = optimize.least_squares(self.residuals, p0, bounds=(lo, hi),
                                                 x_scale="jac", diff_step=1e-3, max_nfev=400)
                except (ValueError, np.linalg.LinAlgError):
                    continue
                if best_here is None or res.cost < best_here.cost:
                    best_here = res
            if best_here is None:
                continue
            tried.append({"b_start": g["b"], "b_end": float(best_here.x[0]),
                          "chi2": float(2 * best_here.cost)})
            if best is None or best_here.cost < best.cost:
                best, best_guess = best_here, g
        if best is None:
            raise RuntimeError("every start failed")
        self.best = best
        self.tried = tried
        return self.summary(best, best_guess)

    def summary(self, res, guess):
        p = res.x
        par = self.unpack(p)
        model, amps, chi2, src = self.model_image(p)
        n = int(self.fit_mask.sum())
        dof = n - len(p) - len(amps)
        # Formal errors from the Jacobian, inflated by the reduced chi^2.
        try:
            cov = np.linalg.inv(res.jac.T @ res.jac) * max(chi2 / dof, 1.0)
            err_b = float(np.sqrt(cov[0, 0]))
        except np.linalg.LinAlgError:
            err_b = float("nan")
        src_amp = float(amps[-2])
        # Signal-to-noise of the lensed source in the data.
        sn_source = float(src_amp * np.sqrt(np.sum(src[self.fit_mask] ** 2 / self.cut.var[self.fit_mask])))
        out = {
            "name": self.cut.name,
            "b": float(par["b"]), "b_err_formal": err_b,
            "q": float(par["q_m"]), "pa": angle_to_pa(par["phi_m"], self.cut),
            "source_x": float(par["xs"] - self.light_shape[0]),
            "source_y": float(par["ys"] - self.light_shape[1]),
            "source_r": float(par["r_s"]), "source_n": float(par["n_s"]),
            "source_q": float(par["q_s"]),
            "source_sn": sn_source,
            "chi2": chi2, "n_pix": n, "chi2_red": chi2 / dof,
            "light_q": float(self.light_shape[2]),
            "light_q_inner": float(self.light_inner[0]),
            "light_pa_inner": angle_to_pa(self.light_inner[1], self.cut),
            "light_pa": angle_to_pa(self.light_shape[3], self.cut),
            "light_centre": [float(self.light_shape[0]), float(self.light_shape[1])],
            "b_initial": guess["b"],
            "candidates": getattr(self, "tried", []),
            "n_cosmic_ray_pixels": int(self.sharp.sum()),
            "manual_mask_pixels": int(self.manual.sum()),
            "fit_radius": float(np.max(self.r[self.fit_mask])),
            "sigma_bkg": self.cut.sigma_bkg, "exptime": self.cut.exptime,
            "params": [float(v) for v in p],
        }
        if "shear" in par:
            out["shear"] = par["shear"][0]
            out["shear_pa"] = angle_to_pa(par["shear"][1], self.cut)
        self.model, self.source_image = model, src * src_amp
        self.light_only = np.tensordot(np.delete(amps, -2), np.delete(
            np.concatenate([self.light_basis[:-1], src[None], self.light_basis[-1:]]), -2, axis=0), axes=1)
        return out

    def save_images(self, path):
        """Store the images the figures need: data, lens light, model, residual."""
        np.savez_compressed(
            path, data=self.cut.sci.astype(np.float32), model=self.model.astype(np.float32),
            light=self.light_only.astype(np.float32), var=self.cut.var.astype(np.float32),
            fit_mask=self.fit_mask, source=self.source_image.astype(np.float32),
            light_shape=np.array(self.light_shape), params=self.best.x,
            light_full=self.light_model.astype(np.float32), sharp=self.sharp, manual=self.manual,
            north_angle=self.cut.north_angle, parity=self.cut.parity,
            centre_pix=np.array(self.cut.centre_pix))


#: Three sets of fits: "" the automatic ones (the headline), "shear" with an
#: external shear added, "hand_masked" refits of the lenses in config.MANUAL_MASKS.
VARIANTS = ("", "shear", "hand_masked")


def model_path(name: str, variant: str = "", suffix: str = ".json"):
    if variant == "hand_masked":
        (MODEL_DIR / "hand_masked").mkdir(exist_ok=True)
        return MODEL_DIR / "hand_masked" / f"{name}{suffix}"
    return MODEL_DIR / f"{name}{'_' + variant if variant else ''}{suffix}"


def fit_lens(name: str, ra: float, dec: float, variant: str = "", save: bool = True) -> dict:
    cut = load_cutout(name, ra, dec)
    fitter = LensFitter(cut, hand_masks=(variant == "hand_masked"))
    fitter.fit_light()
    result = fitter.fit(with_shear=(variant == "shear"))
    if save:
        fitter.save_images(model_path(name, variant, ".npz"))
        model_path(name, variant).write_text(json.dumps(result, indent=1))
    return result


# --------------------------------------------------------------------------
# Every lens
# --------------------------------------------------------------------------

def _single_thread():
    """One BLAS thread per worker: the pool already uses the cores."""
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def _fit_one(args):
    name, ra, dec, variant = args
    try:
        return fit_lens(name, ra, dec, variant=variant)
    except Exception as err:  # recorded, reported, and the batch carries on
        out = {"name": name, "failed": f"{type(err).__name__}: {err}"}
        model_path(name, variant).write_text(json.dumps(out, indent=1))
        return out


def cached_result(name: str, variant: str = "") -> dict | None:
    path = model_path(name, variant)
    return json.loads(path.read_text()) if path.exists() else None


def fit_all(lenses, variant: str = "", processes: int = 6, refit: bool = False,
            quiet: bool = False) -> list[dict]:
    """Fit every lens that has an F814W image, in parallel; results are cached
    as JSON (and images as .npz) under data/lens_models/. The "hand_masked"
    variant refits only the lenses that config.MANUAL_MASKS names."""
    from multiprocessing import get_context

    from config import MANUAL_MASKS

    todo = []
    for row in lenses.itertuples():
        if not (HST_DIR / row.name / f"{FILTER}.fits").exists():
            continue
        if variant == "hand_masked" and row.name not in MANUAL_MASKS:
            continue
        if refit or cached_result(row.name, variant) is None:
            todo.append((row.name, float(row.ra), float(row.dec), variant))
    if todo:
        if not quiet:
            print(f"  fitting {len(todo)} lenses ({variant or 'automatic'}) "
                  f"on {processes} processes ...", flush=True)
        with get_context("fork").Pool(processes, initializer=_single_thread) as pool:
            for res in pool.imap_unordered(_fit_one, todo):
                if not quiet:
                    msg = res.get("failed") or f"b = {res['b']:.3f}\""
                    print(f"    {res['name']}: {msg}", flush=True)
    out = []
    for row in lenses.itertuples():
        res = cached_result(row.name, variant)
        if res is not None:
            out.append(res)
    return out


if __name__ == "__main__":
    import argparse

    import catalog

    parser = argparse.ArgumentParser(description="Fit lens models to every grade-A lens")
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--variant", choices=VARIANTS, default="",
                        help="'shear' adds external shear; 'hand_masked' refits the hand-masked lenses")
    parser.add_argument("--refit", action="store_true", help="ignore cached fits")
    parser.add_argument("names", nargs="*", help="only these lenses")
    args = parser.parse_args()
    lenses = catalog.read_grade_a()
    if args.names:
        lenses = lenses[lenses["name"].isin(args.names)]
    fit_all(lenses, variant=args.variant, processes=args.processes, refit=args.refit)
