"""A small regularised-maximum-likelihood (RML) imager that uses only closure quantities.

With ~250 Fourier samples and 4096 pixels, infinitely many images fit the data.
RML picks the one that balances fitting the data against simple preferences
("regularisers"): stay close to a featureless prior blob (maximum entropy), be
smooth (total squared variation) and have sharp-but-few edges (total variation).
This is the approach of the EHT's eht-imaging and SMILI pipelines (M87 Paper IV,
section 6). Because we fit only closure phases and log closure amplitudes, the
unknown station gains never enter, and no self-calibration is needed.

The image is I = F * softmax(z): always positive, with a fixed total flux F
(closure quantities cannot measure flux). Gradients are analytic, so L-BFGS can
optimise all pixels at once.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import minimize

from config import LMT_GAIN_SIGMA, RAD_PER_UAS, SEFD_ERROR_BUDGET
from uvdata import Prepared, debiased_amplitude


def gain_sigma(station: str) -> float:
    """Prior width of a station's log amplitude gain (see config.SEFD_ERROR_BUDGET)."""
    if station == "LM":
        return LMT_GAIN_SIGMA
    return 0.5 * np.log1p(SEFD_ERROR_BUDGET.get(station, 0.15))


def pixel_grid(npix: int, fov_uas: float):
    """Pixel centres, ehtim's layout: column 0 = east (+x), row 0 = north (+y)."""
    psize = fov_uas / npix
    xs = (np.arange(0, -npix, -1) * psize + 0.5 * fov_uas - 0.5 * psize)
    X, Y = np.meshgrid(xs, xs)
    return X, Y, psize


def gaussian_image(npix: int, fov_uas: float, fwhm_uas: float, flux: float = 1.0):
    X, Y, _ = pixel_grid(npix, fov_uas)
    s = fwhm_uas / np.sqrt(8 * np.log(2))
    img = np.exp(-(X**2 + Y**2) / (2 * s**2))
    return flux * img / img.sum()


@dataclass
class _Block:
    A: np.ndarray  # DFT matrix (n_vis x n_pix)
    D: np.ndarray  # closure design matrix (n_closure x n_vis)
    obs: np.ndarray
    W: np.ndarray  # whitener
    kind: str
    n: int


@dataclass
class _AmpBlock:
    A: np.ndarray  # DFT matrix of the baselines used
    obs: np.ndarray  # debiased amplitudes (Jy)
    sigma: np.ndarray
    g1: np.ndarray  # index of each end's gain parameter
    g2: np.ndarray
    n: int


class ClosureImager:
    """RML imaging from closure phases and log closure amplitudes of one or more data sets.

    Data sets (e.g. the two bands of one night) are combined, but closure
    quantities are never formed across them.
    """

    def __init__(self, preps: list[Prepared], *, npix: int = 64, fov_uas: float = 128.0, flux: float = 0.6,
                 prior_fwhm_uas: float = 40.0, amplitudes: bool = False):
        self.npix, self.fov, self.flux = npix, fov_uas, flux
        self.X, self.Y, self.psize = pixel_grid(npix, fov_uas)
        x = self.X.ravel() * RAD_PER_UAS
        y = self.Y.ravel() * RAD_PER_UAS
        self.prior = gaussian_image(npix, fov_uas, prior_fwhm_uas, 1.0).ravel()
        self.blocks: list[_Block] = []
        for prep in preps:
            A = np.exp(2j * np.pi * (np.outer(prep.avg.u, x) + np.outer(prep.avg.v, y)))
            for cs in (prep.cphase, prep.lcamp):
                if len(cs):
                    # keep only the visibilities this closure set uses
                    used = np.flatnonzero(np.any(cs.design != 0, axis=0))
                    self.blocks.append(_Block(A[used], cs.design[:, used], cs.values, cs.whitener(), cs.kind, len(cs)))
        self.n_cp = sum(b.n for b in self.blocks if b.kind == "cphase")
        self.n_la = sum(b.n for b in self.blocks if b.kind == "lcamp")
        # Optional visibility amplitudes, each station's gain free per scan within its prior.
        # Intra-site baselines are left out: they see the arcsecond-scale jet, not the image.
        self.amp_blocks: list[_AmpBlock] = []
        self.gain_keys: list[tuple] = []
        if amplitudes:
            for k, prep in enumerate(preps):
                d = prep.avg
                amp = debiased_amplitude(np.abs(d.vis), d.sigma)
                use = (~d.intra_site()) & (amp > 0)
                index = {}
                g1, g2 = [], []
                for i in np.flatnonzero(use):
                    for st, lst in ((d.t1[i], g1), (d.t2[i], g2)):
                        key = (k, int(d.scan[i]), str(st))
                        if key not in index:
                            index[key] = len(self.gain_keys)
                            self.gain_keys.append(key)
                        lst.append(index[key])
                A = np.exp(2j * np.pi * (np.outer(d.u[use], x) + np.outer(d.v[use], y)))
                self.amp_blocks.append(_AmpBlock(A, amp[use], d.sigma[use], np.array(g1), np.array(g2), int(use.sum())))
        self.n_amp = sum(b.n for b in self.amp_blocks)
        self.n_gain = len(self.gain_keys)
        self.gain_prior_sigma = np.array([gain_sigma(k[2]) for k in self.gain_keys])
        self.gains = np.zeros(self.n_gain)
        # reference values that make the regularisers dimensionless
        p0 = self.prior.reshape(npix, npix)
        self._tsv0 = self._tsv(p0)[0]
        self._tv0 = self._tv(p0)[0]

    # ---- regularisers (on the normalised image p, sum p = 1) ----
    def _mem(self, p):
        q = self.prior
        val = np.sum(p * np.log(p / q))
        return val, np.log(p / q) + 1.0

    @staticmethod
    def _tsv(img2d):
        dx = np.diff(img2d, axis=1)
        dy = np.diff(img2d, axis=0)
        val = np.sum(dx**2) + np.sum(dy**2)
        g = np.zeros_like(img2d)
        g[:, :-1] -= 2 * dx
        g[:, 1:] += 2 * dx
        g[:-1, :] -= 2 * dy
        g[1:, :] += 2 * dy
        return val, g

    @staticmethod
    def _tv(img2d, eps=1e-10):
        dx = np.zeros_like(img2d)
        dy = np.zeros_like(img2d)
        dx[:, :-1] = np.diff(img2d, axis=1)
        dy[:-1, :] = np.diff(img2d, axis=0)
        mag = np.sqrt(dx**2 + dy**2 + eps)
        val = np.sum(mag)
        g = np.zeros_like(img2d)
        gx, gy = dx / mag, dy / mag
        g -= gx + gy
        g[:, 1:] += gx[:, :-1]
        g[1:, :] += gy[:-1, :]
        return val, g

    def _centroid(self, p):
        cx, cy = np.sum(p * self.X.ravel()), np.sum(p * self.Y.ravel())
        scale = (self.fov / 8.0) ** 2
        val = (cx**2 + cy**2) / scale
        g = 2 * (cx * self.X.ravel() + cy * self.Y.ravel()) / scale
        return val, g

    # ---- data terms ----
    def data_terms(self, img_flat):
        """Reduced chi^2 of closure phases and log closure amplitudes, with gradients w.r.t. pixels."""
        chi = {"cphase": 0.0, "lcamp": 0.0}
        grad = {"cphase": np.zeros_like(img_flat), "lcamp": np.zeros_like(img_flat)}
        for b in self.blocks:
            V = b.A @ img_flat
            if b.kind == "cphase":
                r = np.angle(np.exp(1j * (b.obs - b.D @ np.angle(V))))
            else:
                r = b.obs - b.D @ np.log(np.abs(V))
            w = b.W @ r
            chi[b.kind] += w @ w
            g_model = -2.0 * (b.D.T @ (b.W.T @ w))  # d chi2 / d (phase or ln amp) per visibility
            if b.kind == "cphase":
                grad[b.kind] += np.imag(b.A.T @ (g_model / V))
            else:
                grad[b.kind] += np.real(b.A.T @ (g_model / V))
        n = {"cphase": max(self.n_cp, 1), "lcamp": max(self.n_la, 1)}
        return ({k: chi[k] / n[k] for k in chi}, {k: grad[k] / n[k] for k in grad})

    def amp_term(self, img_flat, gains):
        """(chi^2 of amplitudes + gain prior) / N_amp, and its gradients."""
        val = 0.0
        g_img = np.zeros_like(img_flat)
        g_gain = np.zeros(self.n_gain)
        for b in self.amp_blocks:
            V = b.A @ img_flat
            absV = np.abs(V)
            scale = np.exp(gains[b.g1] + gains[b.g2])
            model = scale * absV
            r = (b.obs - model) / b.sigma
            val += r @ r
            dmodel = -2.0 * r / b.sigma  # d chi2 / d model amplitude
            g_img += np.real(b.A.T @ (dmodel * scale * np.conj(V) / absV))
            np.add.at(g_gain, b.g1, dmodel * model)
            np.add.at(g_gain, b.g2, dmodel * model)
        prior = gains / self.gain_prior_sigma
        val += prior @ prior
        g_gain += 2.0 * gains / self.gain_prior_sigma**2
        n = max(self.n_amp, 1)
        return val / n, g_img / n, g_gain / n

    def objective(self, x, weights):
        npx = self.npix * self.npix
        z = x[:npx] - x[:npx].max()
        gains = x[npx:] if self.n_gain else np.zeros(0)
        e = np.exp(z)
        p = e / e.sum()
        img = self.flux * p
        red, rgrad = self.data_terms(img)
        val = weights["cphase"] * red["cphase"] + weights["lcamp"] * red["lcamp"]
        g_img = weights["cphase"] * rgrad["cphase"] + weights["lcamp"] * rgrad["lcamp"]
        g_gain = np.zeros(self.n_gain)
        if self.n_amp and weights.get("amp", 0):
            v, gi, gg = self.amp_term(img, gains)
            val += weights["amp"] * v
            g_img += weights["amp"] * gi
            g_gain += weights["amp"] * gg
        g_p = self.flux * g_img
        p2 = p.reshape(self.npix, self.npix)
        if weights.get("mem", 0):
            v, g = self._mem(p)
            val += weights["mem"] * v
            g_p += weights["mem"] * g
        if weights.get("tsv", 0):
            v, g = self._tsv(p2)
            val += weights["tsv"] * v / self._tsv0
            g_p += weights["tsv"] * g.ravel() / self._tsv0
        if weights.get("tv", 0):
            v, g = self._tv(p2)
            val += weights["tv"] * v / self._tv0
            g_p += weights["tv"] * g.ravel() / self._tv0
        if weights.get("centroid", 0):
            v, g = self._centroid(p)
            val += weights["centroid"] * v
            g_p += weights["centroid"] * g
        g_z = p * (g_p - np.sum(p * g_p))  # chain rule through the softmax
        return val, np.concatenate([g_z, g_gain])

    def solve(self, init_img, weights, maxiter=600):
        z0 = np.log(np.clip(init_img.ravel() / init_img.sum(), 1e-12, None))
        x0 = np.concatenate([z0, self.gains])
        res = minimize(self.objective, x0, args=(weights,), jac=True, method="L-BFGS-B",
                       options={"maxiter": maxiter, "maxcor": 20})
        npx = self.npix * self.npix
        z = res.x[:npx] - res.x[:npx].max()
        self.gains = res.x[npx:].copy()
        p = np.exp(z) / np.exp(z).sum()
        return (self.flux * p).reshape(self.npix, self.npix)

    def reduced_chi2(self, img) -> dict:
        red, _ = self.data_terms(img.ravel())
        if self.n_amp:
            v, _, _ = self.amp_term(img.ravel(), self.gains)
            red["amp"] = v
        return red

    def image(self, *, rounds: int = 4, weights=None, blur_uas: float = 20.0, init=None, data_growth: float = 1.0):
        """Multi-round imaging in the style of the EHT pipelines.

        Start from the Gaussian prior; after each round blur the result to the
        array's resolution and image again (optionally with a heavier data weight),
        which helps the optimiser out of poor local minima.
        """
        w = dict(cphase=1.0, lcamp=1.0, amp=1.0, mem=1.0, tsv=1.0, tv=0.3, centroid=0.1)
        if weights:
            w.update(weights)
        img = (init if init is not None else self.prior.reshape(self.npix, self.npix) * self.flux)
        history = []
        for k in range(rounds):
            if k > 0:
                img = gaussian_filter(img, blur_uas / np.sqrt(8 * np.log(2)) / self.psize)
                w["cphase"] *= data_growth
                w["lcamp"] *= data_growth
            img = self.solve(img, w)
            history.append(self.reduced_chi2(img))
        self.history = history
        return img
