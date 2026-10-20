"""Our own time-delay estimator: one intrinsic quasar curve, shifted, plus microlensing.

Every image of a lensed quasar shows the *same* intrinsic flickering, shifted
by that image's light-travel delay. On top of that, stars in the lens galaxy
act as microlenses and slowly brighten or dim each image independently
("extrinsic" variability, typically 0.1-0.5 mag over years). The model is

    m_i(t) = s(t - tau_i) + l_i(t) + o_{i, telescope}

* ``s``: the intrinsic light curve, a cubic B-spline with uniform knots every
  ``knot_int`` days (~2-4 weeks) and a curvature penalty (a P-spline; Eilers &
  Marx 1996) so that it stays smooth where no image samples it (season gaps).
* ``l_i``: microlensing of image i, a much stiffer B-spline (knots every
  ``knot_ml`` ~ 150-400 days). The reference image has none: its microlensing
  is absorbed into ``s``, so every other ``l_i`` is microlensing *relative to
  the reference image*.
* ``o_{i,tel}``: one magnitude offset per image and extra telescope
  (instruments differ in zero point, colour response and deblending).

For fixed delays tau the model is linear in all spline coefficients and
offsets, so the best fit is one penalized weighted least-squares solve.
The delays are found by maximizing the Bayesian evidence p(data | tau), with
the spline coefficients integrated out under their Gaussian smoothness prior:

    -2 ln p = chi^2_penalized(tau) + ln det(A'WA + P) + const.

Plain chi^2 has a trap: a delay that pushes one image's points into the
other images' seasonal gaps lets the intrinsic curve fit those points alone,
and chi^2 drops for no physical reason. The log-determinant is the Occam
factor that charges for exactly that freedom (it comes free from the
Cholesky factors of the solve). This is the same idea as PyCS's
free-knot spline technique (Tewes et al. 2013a), but with fixed knots plus a
smoothness penalty instead of optimized knot positions, written from scratch.

Uncertainties come from simulated light curves, as in COSMOGRAIL: the fitted
model is re-sampled at the real epochs with known "true" delays, the real
residuals (which carry the fast microlensing and photometric noise) are added
back after shifting them in time, and the estimator is run blind on each mock.
The spread of (measured - true) gives the random error; its mean, the bias.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
from scipy.interpolate import BSpline
from scipy.linalg import cho_factor, cho_solve, cholesky_banded, cho_solve_banded, LinAlgError
from scipy.optimize import minimize

from lightcurves import LightCurves

K = 3  # cubic splines


def uniform_knots(lo: float, hi: float, step: float) -> np.ndarray:
    """P-spline knot vector: equal spacing, K extra knots beyond each end."""
    n = max(int(np.ceil((hi - lo) / step)), 1)
    h = (hi - lo) / n
    return lo + h * np.arange(-K, n + K + 1)


def basis(x: np.ndarray, knots: np.ndarray) -> sp.csr_matrix:
    """Sparse cubic B-spline design matrix (rows: x, columns: coefficients)."""
    return BSpline.design_matrix(np.asarray(x, float), knots, K).tocsr()


def second_difference(n: int) -> np.ndarray:
    d = np.zeros((n - 2, n))
    i = np.arange(n - 2)
    d[i, i], d[i, i + 1], d[i, i + 2] = 1.0, -2.0, 1.0
    return d


@dataclass
class SplineSettings:
    """Estimator settings (the analogue of PyCS's knot steps)."""

    knot_int: float = 20.0  # days between intrinsic-curve knots
    knot_ml: float = 150.0  # days between microlensing knots
    curv_int: float = 1.0  # mag: prior scale of the intrinsic curve's 2nd difference per knot
    curv_ml: float = 0.3  # mag: same for the microlensing curves
    max_shift: float = 250.0  # days: the intrinsic curve covers t - tau for |tau| up to this

    def label(self) -> str:
        return f"int {self.knot_int:g} d / ML {self.knot_ml:g} d"


@dataclass
class Fit:
    taus: dict[str, float]  # image -> delay relative to the reference image (days)
    chi2: float  # penalized chi^2 at the optimum
    chi2_data: float
    ndata: int
    coef: np.ndarray = field(repr=False, default=None)
    model: dict[str, np.ndarray] = field(repr=False, default=None)  # image -> model at the data epochs
    resid: dict[str, np.ndarray] = field(repr=False, default=None)


class JointSpline:
    """The linear model above for one set of light curves.

    ``ref``: the reference image (delay 0); ``others``: the images whose
    delays are free. Delays are t_image - t_ref (positive = arrives later).
    """

    def __init__(self, lc: LightCurves, ref: str, settings: SplineSettings | None = None,
                 err_scale: dict | None = None, objective: str = "evidence"):
        self.objective = objective  # "evidence" (default) or "chi2"
        self.lc = lc
        self.ref = ref
        self.images = list(lc.images)
        self.others = [im for im in self.images if im != ref]
        self.s = settings or SplineSettings()
        t = lc.t
        self.t = t
        s = self.s
        # Intrinsic curve: must cover t - tau for every trial tau.
        self.kn_int = uniform_knots(t.min() - s.max_shift - 5, t.max() + s.max_shift + 5, s.knot_int)
        self.n_int = len(self.kn_int) - K - 1
        self.kn_ml = uniform_knots(t.min() - 1.0, t.max() + 1.0, s.knot_ml)
        self.n_ml = len(self.kn_ml) - K - 1
        self.ml_basis = basis(t, self.kn_ml)
        tels = lc.telescopes
        self.tels = tels
        self.n_off = len(tels) - 1  # per image; the first telescope has no offset
        self.tel_idx = np.array([tels.index(x) for x in lc.tel])
        # Column layout: [intrinsic | ML of each non-ref image | offsets of each image]
        self.col_ml = {im: self.n_int + j * self.n_ml for j, im in enumerate(self.others)}
        base_off = self.n_int + len(self.others) * self.n_ml
        self.col_off = {im: base_off + j * self.n_off for j, im in enumerate(self.images)}
        self.ncol = base_off + len(self.images) * self.n_off
        # Penalty (Gaussian random-walk prior on curvature).
        P = np.zeros((self.ncol, self.ncol))
        D = second_difference(self.n_int)
        P[: self.n_int, : self.n_int] = D.T @ D / s.curv_int**2
        Dm = second_difference(self.n_ml)
        for im in self.others:
            c = self.col_ml[im]
            P[c: c + self.n_ml, c: c + self.n_ml] = Dm.T @ Dm / s.curv_ml**2
        # Weak ridge on offsets (keeps the system regular if a telescope has
        # no overlap with any other; 1 mag scale, irrelevant otherwise).
        P[base_off:, base_off:] += np.eye(self.ncol - base_off)
        self.P = P
        # Data, weights (errors optionally rescaled per image and telescope).
        self.y = {im: lc.mag[im] for im in self.images}
        self.err_scale = err_scale or {}
        self.w = {}
        for im in self.images:
            sc = np.array([self.err_scale.get((im, tel), 1.0) for tel in lc.tel])
            self.w[im] = 1.0 / (lc.err[im] * sc) ** 2
        # Fixed (tau-independent) parts of the design matrix.
        self._fixed = {}
        for im in self.images:
            rows = []
            if im != ref:
                rows.append((self.col_ml[im], self.ml_basis))
            if self.n_off:
                o = np.zeros((self.lc.n, self.n_off))
                m = self.tel_idx > 0
                o[np.nonzero(m)[0], self.tel_idx[m] - 1] = 1.0
                rows.append((self.col_off[im], sp.csr_matrix(o)))
            self._fixed[im] = rows

    # -- design ------------------------------------------------------------
    def design(self, im: str, tau: float, t=None) -> sp.csr_matrix:
        """Design matrix rows for image ``im`` at its delay ``tau``."""
        t = self.t if t is None else t
        blocks = [(0, basis(t - tau, self.kn_int))]
        if t is self.t:
            blocks += self._fixed[im]
        else:
            if im != self.ref:
                blocks.append((self.col_ml[im], basis(t, self.kn_ml)))
        n = len(t)
        mats = []
        last = 0
        for col, M in sorted(blocks, key=lambda b: b[0]):
            if col > last:
                mats.append(sp.csr_matrix((n, col - last)))
            mats.append(M)
            last = col + M.shape[1]
        if last < self.ncol:
            mats.append(sp.csr_matrix((n, self.ncol - last)))
        return sp.hstack(mats, format="csr")

    # -- normal equations, assembled incrementally ---------------------------
    #
    # Columns split into the intrinsic curve (banded block H: cubic B-splines
    # overlap only their 3 neighbours, and the curvature penalty is
    # pentadiagonal) and the rest (microlensing + offsets, R). The system
    # [[H, X], [X', R]] c = b is solved by a Schur complement: a banded
    # Cholesky of H and a small dense solve for R - X' H^-1 X. Only the
    # blocks of the image whose delay changed are recomputed.

    def _prepare(self):
        ni = self.n_int
        self._F, self._G = {}, self.P[ni:, ni:].copy()
        for im in self.images:
            F = self.design(im, 0.0)[:, ni:].tocsr()
            self._F[im] = F
            self._G += (F.T @ F.multiply(self.w[im][:, None]).tocsr()).toarray()
        Pb = np.zeros((K + 1, ni))
        for k in range(K + 1):
            Pb[K - k, k:] = np.diagonal(self.P[:ni, :ni], k)
        self._Pband = Pb
        self._cache = {im: {} for im in self.images}

    def _intrinsic_blocks(self, im: str, tau: float):
        """(B, H in banded form, X = B'WF) for image ``im`` at delay ``tau`` (cached)."""
        c = self._cache[im]
        key = round(float(tau), 9)
        hit = c.get(key)
        if hit is not None:
            return hit
        B = basis(self.t - tau, self.kn_int)
        Bw = B.multiply(self.w[im][:, None]).tocsr()
        Hs = (B.T @ Bw).tocsr()
        Hb = np.zeros((K + 1, self.n_int))
        for k in range(K + 1):
            Hb[K - k, k:] = Hs.diagonal(k)
        X = (Bw.T @ self._F[im]).toarray()
        if len(c) > 48:  # small FIFO: coordinate descent only revisits recent delays
            c.pop(next(iter(c)))
        c[key] = (B, Hb, X)
        return c[key]

    def solve(self, taus: dict[str, float], y: dict | None = None, full: bool = False):
        """Best-fit coefficients and penalized chi^2 for the given delays.

        Penalized chi^2 at the optimum = y'Wy - b'c, with b = A'Wy and c the
        solution of (A'WA + P) c = b, so no residuals are needed for a scan.
        """
        if not hasattr(self, "_F"):
            self._prepare()
        y = self.y if y is None else y
        ni = self.n_int
        Hb = self._Pband.copy()
        X = np.zeros((ni, self.ncol - ni))
        b = np.zeros(self.ncol)
        yWy = 0.0
        blocks = {}
        for im in self.images:
            tau = 0.0 if im == self.ref else taus[im]
            B, Hbi, Xi = self._intrinsic_blocks(im, tau)
            blocks[im] = B
            Hb += Hbi
            X += Xi
            wy = self.w[im] * y[im]
            b[:ni] += B.T @ wy
            b[ni:] += self._F[im].T @ wy
            yWy += float(y[im] @ wy)
        try:
            cb = cholesky_banded(Hb, lower=False, check_finite=False)
            HiX = cho_solve_banded((cb, False), X, check_finite=False)
            Hib = cho_solve_banded((cb, False), b[:ni], check_finite=False)
            S = self._G - X.T @ HiX
            cf = cho_factor(S, lower=True, check_finite=False)
            c_rest = cho_solve(cf, b[ni:] - X.T @ Hib, check_finite=False)
        except (LinAlgError, ValueError):
            return (np.inf, np.inf, None, None, None) if full else (np.inf, None, None)
        coef = np.concatenate([Hib - HiX @ c_rest, c_rest])
        total = yWy - float(b @ coef)
        # log det(A'WA + P) = log det H + log det(Schur complement)
        self._last_logdet = 2.0 * (np.sum(np.log(cb[K])) + np.sum(np.log(np.diag(cf[0]))))
        if self.objective == "evidence":
            total = total + self._last_logdet
        if not full:
            return total, coef, None
        model, resid, chi2_data = {}, {}, 0.0
        for im in self.images:
            mdl = blocks[im] @ coef[:ni] + self._F[im] @ coef[ni:]
            r = y[im] - mdl
            chi2_data += float(np.sum(self.w[im] * r * r))
            model[im], resid[im] = mdl, r
        return total, chi2_data, coef, model, resid

    def chi2(self, taus: dict[str, float], y=None) -> float:
        return self.solve(taus, y)[0]

    def fit(self, taus: dict[str, float], y=None) -> Fit:
        tot, c2d, coef, model, resid = self.solve(taus, y, full=True)
        n = sum(len(v) for v in self.y.values())
        return Fit(dict(taus), tot, c2d, n, coef, model, resid)

    def intrinsic(self, coef: np.ndarray, t: np.ndarray) -> np.ndarray:
        return basis(t, self.kn_int) @ coef[: self.n_int]

    def microlensing(self, coef: np.ndarray, im: str, t: np.ndarray) -> np.ndarray:
        if im == self.ref:
            return np.zeros_like(np.asarray(t, float))
        c = self.col_ml[im]
        return basis(t, self.kn_ml) @ coef[c: c + self.n_ml]

    def offsets(self, coef: np.ndarray, im: str) -> np.ndarray:
        """Per-epoch telescope offset of image ``im``."""
        if not self.n_off:
            return np.zeros(self.lc.n)
        o = np.concatenate([[0.0], coef[self.col_off[im]: self.col_off[im] + self.n_off]])
        return o[self.tel_idx]

    # -- optimization ------------------------------------------------------
    def optimize(self, start: dict[str, float], width: float = 10.0, step: float = 0.5,
                 rounds: int = 3, y=None, polish: bool = True) -> dict[str, float]:
        """Coordinate descent over 1-D delay scans, then a Nelder-Mead polish.

        Each round scans every free delay over +-``width`` days around the
        current value (others held fixed) and moves to the best grid point;
        the scan half-width shrinks each round. This finds the right basin
        far more reliably than a local optimizer started cold.
        """
        taus = dict(start)
        w = width
        for _ in range(rounds):
            for im in self.others:
                grid = taus[im] + np.arange(-w, w + step / 2, step)
                vals = []
                for g in grid:
                    trial = dict(taus)
                    trial[im] = g
                    vals.append(self.chi2(trial, y))
                vals = np.array(vals)
                j = int(np.argmin(vals))
                taus[im] = float(grid[j])
                # parabolic refinement between neighbours
                if 0 < j < len(grid) - 1 and np.all(np.isfinite(vals[j - 1: j + 2])):
                    a, b, c = vals[j - 1: j + 2]
                    den = a - 2 * b + c
                    if den > 0:
                        taus[im] = float(grid[j] + 0.5 * step * (a - c) / den)
            w = max(w / 2.0, 2 * step)
        if polish:
            x0 = np.array([taus[im] for im in self.others])
            res = minimize(lambda x: self.chi2(dict(zip(self.others, x)), y), x0, method="Nelder-Mead",
                           options={"xatol": 0.02, "fatol": 0.01, "initial_simplex": np.vstack([x0] + [x0 + np.eye(len(x0))[i] * 0.7 for i in range(len(x0))])})
            if res.fun <= self.chi2(taus, y):
                taus = dict(zip(self.others, map(float, res.x)))
        return taus

    def scan(self, im: str, grid: np.ndarray, taus: dict[str, float], profile: bool = False,
             width: float = 6.0) -> np.ndarray:
        """chi^2 as a function of one delay, others fixed (or re-optimized if ``profile``)."""
        out = []
        cur = dict(taus)
        for g in grid:
            trial = dict(cur)
            trial[im] = float(g)
            if profile and len(self.others) > 1:
                rest = [o for o in self.others if o != im]
                sub = JointSplineProfile(self, im, float(g))
                best = sub.optimize({o: cur[o] for o in rest}, width=width)
                trial.update(best)
                cur.update(best)
            out.append(self.chi2(trial))
        return np.array(out)


class JointSplineProfile:
    """Helper: optimize the other delays with one delay held fixed."""

    def __init__(self, model: JointSpline, fixed_im: str, fixed_tau: float):
        self.m, self.fi, self.ft = model, fixed_im, fixed_tau

    def optimize(self, start: dict[str, float], width: float = 6.0) -> dict[str, float]:
        taus = dict(start)
        for _ in range(2):
            for im in taus:
                grid = taus[im] + np.arange(-width, width + 0.25, 0.5)
                vals = [self.m.chi2({**taus, im: g, self.fi: self.ft}) for g in grid]
                taus[im] = float(grid[int(np.argmin(vals))])
            width = max(width / 2, 1.0)
        return taus


# ---------------------------------------------------------------------------
# Error-bar rescaling
# ---------------------------------------------------------------------------


def rescale_errors(lc: LightCurves, ref: str, taus: dict, settings: SplineSettings, iterations: int = 2,
                   min_points: int = 15) -> dict:
    """Per (image, telescope) factors that make the reduced chi^2 of each group ~1.

    Deconvolution-photometry error bars are internal estimates; the scatter
    about a smooth model also contains fast microlensing and systematics. The
    factors are never allowed below 1 (we only inflate errors).
    """
    scale = {}
    for _ in range(iterations):
        m = JointSpline(lc, ref, settings, err_scale=scale)
        fit = m.fit(taus)
        new = {}
        for im in lc.images:
            for tel in lc.telescopes:
                sel = lc.tel == tel
                if sel.sum() < min_points:
                    new[(im, tel)] = scale.get((im, tel), 1.0)
                    continue
                r = fit.resid[im][sel] / (lc.err[im][sel] * scale.get((im, tel), 1.0))
                new[(im, tel)] = scale.get((im, tel), 1.0) * max(1.0, float(np.sqrt(np.mean(r * r))))
        scale = new
    return scale


# ---------------------------------------------------------------------------
# Simulated light curves for the error budget
# ---------------------------------------------------------------------------


def make_mock(model: JointSpline, fit: Fit, true_taus: dict, rng: np.random.Generator) -> dict:
    """Mock magnitudes at the real epochs: model with ``true_taus`` + shifted real residuals.

    The residuals of each image are rolled by a random number of epochs (the
    same roll for all images would preserve their mutual correlation; we roll
    independently because microlensing is independent between images).
    """
    coef = fit.coef
    y = {}
    for im in model.images:
        tau = 0.0 if im == model.ref else true_taus[im]
        mdl = model.design(im, tau) @ coef
        r = fit.resid[im]
        y[im] = mdl + np.roll(r, int(rng.integers(1, len(r))))
    return y


def mock_errors(model: JointSpline, fit: Fit, n: int = 100, spread: float = 5.0, seed: int = 1,
                width: float = 8.0) -> dict:
    """Run the estimator on ``n`` mocks. Returns per-image arrays of (measured - true)."""
    rng = np.random.default_rng(seed)
    out = {im: [] for im in model.others}
    for _ in range(n):
        true = {im: fit.taus[im] + rng.uniform(-spread, spread) for im in model.others}
        y = make_mock(model, fit, true, rng)
        start = {im: fit.taus[im] + rng.uniform(-spread, spread) for im in model.others}
        est = model.optimize(start, width=width, y=y)
        for im in model.others:
            out[im].append(est[im] - true[im])
    return {im: np.array(v) for im, v in out.items()}


def pair_delays(taus: dict[str, float], ref: str) -> dict:
    """Delays for every image pair from delays relative to ``ref``: (X, Y) -> t_Y - t_X."""
    full = {ref: 0.0, **taus}
    ims = list(full)
    return {(a, b): full[b] - full[a] for i, a in enumerate(ims) for b in ims[i + 1:]}


# ---------------------------------------------------------------------------
# One complete measurement
# ---------------------------------------------------------------------------


@dataclass
class Measurement:
    lens: str
    ref: str
    settings: SplineSettings
    taus: dict  # image -> delay vs ref (days)
    err: dict  # image -> total 1-sigma error (random + bias in quadrature)
    rand: dict  # image -> random error (spread of mock estimates)
    bias: dict  # image -> mean (measured - true) on mocks
    cov: np.ndarray  # covariance of the free delays from the mocks (order: others)
    others: list
    chi2: float
    ndata: int
    err_scale: dict
    mock_dev: dict = field(repr=False, default=None)


def measure(lc: LightCurves, ref: str, start: dict, settings: SplineSettings, n_mock: int = 100,
            width: float = 10.0, spread: float = 5.0, seed: int = 1) -> Measurement:
    """Rescale errors, find the delays, then run the mock error analysis.

    Errors are rescaled at the starting delays first (so the wide search
    already weighs the telescopes sensibly), then again at the optimum.
    """
    scale = rescale_errors(lc, ref, start, settings, iterations=1)
    taus = JointSpline(lc, ref, settings, err_scale=scale).optimize(start, width=width)
    scale = rescale_errors(lc, ref, taus, settings)
    model = JointSpline(lc, ref, settings, err_scale=scale)
    taus = model.optimize(taus, width=min(width, 6.0))
    fit = model.fit(taus)
    dev = mock_errors(model, fit, n=n_mock, spread=spread, seed=seed, width=spread + 3.0) if n_mock else None
    others = model.others
    if dev is not None:
        M = np.array([dev[im] for im in others])
        cov = np.atleast_2d(np.cov(M))
        rand = {im: float(np.std(dev[im], ddof=1)) for im in others}
        bias = {im: float(np.mean(dev[im])) for im in others}
        err = {im: float(np.hypot(rand[im], bias[im])) for im in others}
    else:
        cov = np.full((len(others), len(others)), np.nan)
        rand = bias = err = {im: np.nan for im in others}
    return Measurement(lc.key, ref, settings, taus, err, rand, bias, cov, others, fit.chi2_data, fit.ndata, scale, dev)
