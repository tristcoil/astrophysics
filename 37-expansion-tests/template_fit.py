"""Our own time-dilation measurement: one template for all supernovae, b free.

Every measured flux is modelled as

    flux_ij = A_i,band * S(tau_ij, lambda_ij),   tau = t / D(z_i),   lambda = lambda_band / (1 + z_i)

* t is the observer-frame time since the supernova's peak;
* D(z) is the time-dilation factor under test: D = (1 + z)^b for a power law
  (b = 1 expansion, b = 0 tired light), or one free factor per redshift bin;
* lambda is the rest-frame wavelength the band samples. Wavelength stretching
  by (1 + z) is the redshift itself, which tired light also accepts, so this
  carries no time-dilation assumption;
* S is a smooth surface (cubic B-splines in tau and in log lambda), fitted to
  all supernovae at once. It knows nothing about D: it is whatever single
  light-curve family best describes the data once times are divided by D;
* A_i,band is a free amplitude for each supernova and band (absorbing
  distance, dust, colour and calibration).

The model is linear in S for fixed amplitudes and in the amplitudes for fixed
S. For a trial D both are found by Gauss-Newton with the amplitudes profiled
out (variable projection); D's parameters can be fitted in the same way. The
amplitudes make S's overall brightness at each wavelength unobservable, so S
is normalised to S(0, lambda) = 1: one at the supernova's peak date. The
chi-square is compared between trial D's on exactly the same data points with
exactly the same errors: the observer-frame window is fixed, so no choice of D
changes which measurements are used or how they are weighted.

Supernovae also differ in intrinsic light-curve width ("stretch", about 10%),
which the single template cannot follow. Each point's error therefore gets a
floor of ``eps`` times its own light curve's observed peak flux, so that a
well measured nearby supernova does not outvote everything else with
residuals that are really stretch. The floor is computed from the data alone,
the same for every trial D. Uncertainties on b come from bootstrapping
supernovae, which includes the stretch scatter honestly.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import interpolate, optimize, sparse

from config import DES_FILTERS

#: Observer-frame window used by the fit (days from peak). Fixed for every
#: trial D. For D >= 1 every point maps inside the template's tau range.
FIT_WINDOW = (-25.0, 85.0)
TAU_RANGE = (-26.0, 86.0)
#: Rest-frame wavelength range of the template (nm). Below ~300 nm a Type Ia
#: supernova is faint and its light curve diverse; such points are not used.
LAMBDA_RANGE = (300.0, 860.0)
TAU_KNOT_STEP = 3.0
LOGLAM_KNOT_STEP = 0.045
#: Curvature penalty on the spline coefficients (keeps sparsely sampled
#: corners of the surface well behaved; negligible where data are dense).
SMOOTHING = 1e-2
#: Error floor as a fraction of each light curve's observed peak flux.
EPS_DEFAULT = 0.05
MIN_POINTS_PER_CURVE = 3


def _knots(inner, k=3):
    return np.concatenate([[inner[0]] * k, inner, [inner[-1]] * k])


def _basis(x, inner, k=3):
    """Sparse B-spline design matrix (n_points x n_basis) on clamped knots."""
    return interpolate.BSpline.design_matrix(np.clip(x, inner[0], inner[-1]), _knots(inner, k), k)


def _n_basis(inner, k=3):
    return len(inner) + k - 1


def _second_difference(n):
    d = np.zeros((n - 2, n))
    for i in range(n - 2):
        d[i, i:i + 3] = (1.0, -2.0, 1.0)
    return d


def curve_peak_flux(points: pd.DataFrame) -> pd.Series:
    """Each light curve's peak flux, estimated from its own data only.

    The largest value of a three-point running median of the fluxes in time
    order (robust to a single outlier); at least the curve's flux error.
    """
    def peak(g):
        f = g.sort_values("t")["flux"].to_numpy()
        if f.size >= 3:
            f = np.median(np.lib.stride_tricks.sliding_window_view(f, 3), axis=1)
        return max(f.max(), g["err"].median())
    return points.groupby("curve", sort=False)[["t", "flux", "err"]].apply(peak)


class TemplateFitter:
    """Holds the data once; ``fit(dilation)`` returns the chi-square for any D."""

    def __init__(self, points: pd.DataFrame, *, eps=EPS_DEFAULT, smoothing=SMOOTHING,
                 tau_step=TAU_KNOT_STEP, loglam_step=LOGLAM_KNOT_STEP,
                 lambda_range=LAMBDA_RANGE, fit_window=FIT_WINDOW):
        lam = np.array([DES_FILTERS[b][0] for b in points["band"]]) / (1.0 + points["z"].to_numpy())
        keep = ((points["t"] >= fit_window[0]) & (points["t"] <= fit_window[1])
                & (lam >= lambda_range[0]) & (lam <= lambda_range[1])).to_numpy()
        p = points[keep].copy()
        p["loglam"] = np.log10(lam[keep])
        # Light curves (supernova + band) with too few points cannot pin an amplitude.
        p["curve"] = p["cid"].astype(str) + p["band"]
        counts = p.groupby("curve")["t"].transform("size")
        p = p[counts >= MIN_POINTS_PER_CURVE].reset_index(drop=True)
        peak = curve_peak_flux(p)
        p["peak"] = p["curve"].map(peak)
        self.points = p
        self.curve_codes, self.curve_names = pd.factorize(p["curve"])
        self.n_curves = len(self.curve_names)
        self.cids = p["cid"].to_numpy()
        self.z = p["z"].to_numpy()
        self.t = p["t"].to_numpy()
        self.f = p["flux"].to_numpy()
        self.e = p["err"].to_numpy()
        self.loglam = p["loglam"].to_numpy()
        self.eps = eps
        self.sigma2 = self.e ** 2 + (eps * p["peak"].to_numpy()) ** 2

        self.tau_knots = np.arange(TAU_RANGE[0], TAU_RANGE[1] + 1e-9, tau_step)
        lo, hi = np.log10(lambda_range[0]), np.log10(lambda_range[1])
        self.lam_knots = np.linspace(lo, hi, int(round((hi - lo) / loglam_step)) + 1)
        self.n_tau, self.n_lam = _n_basis(self.tau_knots), _n_basis(self.lam_knots)
        self.L = _basis(self.loglam, self.lam_knots).toarray()
        self.n_coef = self.n_tau * self.n_lam
        dt, dl = _second_difference(self.n_tau), _second_difference(self.n_lam)
        self.penalty = smoothing * (np.kron(dt.T @ dt, np.eye(self.n_lam))
                                    + np.kron(np.eye(self.n_tau), dl.T @ dl))
        # Normalisation. Each light curve has its own amplitude, so the
        # template's brightness at each wavelength is unobservable: S -> g(lambda) S
        # with A -> A / g changes (almost) nothing, a flat valley that stalls any
        # optimiser. Remove it exactly by requiring S(tau = 0, lambda) = 1 at
        # every wavelength: the template is normalised to its value at the
        # supernova's (SALT3) peak date. tau = 0 is t = 0 for every D, so the
        # normalisation does not depend on the dilation under test. Because the
        # lambda B-splines sum to one, it is enough that, in every lambda column,
        # sum_a T_a(0) c[a, b] = 1. One coefficient per column (a*) is eliminated:
        #     c = c0 + M u,  u = the free coefficients.
        t0 = _basis(np.array([0.0]), self.tau_knots).toarray()[0]
        nz = np.nonzero(t0)[0]
        astar = int(nz[np.argmax(t0[nz])])
        others = [a for a in nz if a != astar]
        idx = lambda a, b: a * self.n_lam + b
        eliminated = [idx(astar, b) for b in range(self.n_lam)]
        self.free = np.setdiff1d(np.arange(self.n_coef), eliminated)
        col_of = {c: j for j, c in enumerate(self.free)}
        rows, cols, vals = [], [], []
        for c in self.free:
            rows.append(c); cols.append(col_of[c]); vals.append(1.0)
        self.c0 = np.zeros(self.n_coef)
        for b in range(self.n_lam):
            self.c0[idx(astar, b)] = 1.0 / t0[astar]
            for a in others:
                rows.append(idx(astar, b)); cols.append(col_of[idx(a, b)])
                vals.append(-t0[a] / t0[astar])
        self.M = sparse.csr_matrix((vals, (rows, cols)), shape=(self.n_coef, self.free.size))
        self.Pu = np.asarray(self.M.T @ self.penalty @ self.M)
        self.Pc0 = self.M.T @ (self.penalty @ self.c0)
        self.coef = None   # template coefficients of the last fit
        self.u = None      # its free part (warm start between calls)

    # -- the model -----------------------------------------------------------

    def design(self, tau):
        """Sparse design matrix of the tensor-product basis at (tau, log lambda)."""
        T = _basis(tau, self.tau_knots).tocoo()
        # Each row has 4 tau and 4 lambda non-zeros: build the 16 products directly.
        n = len(tau)
        trow = np.argsort(T.row, kind="stable")
        tcol = T.col[trow].reshape(n, -1)
        tval = T.data[trow].reshape(n, -1)
        lcol = np.argsort(-self.L, axis=1)[:, :4]
        lval = np.take_along_axis(self.L, lcol, axis=1)
        cols = (tcol[:, :, None] * self.n_lam + lcol[:, None, :]).reshape(n, -1)
        vals = (tval[:, :, None] * lval[:, None, :]).reshape(n, -1)
        rows = np.repeat(np.arange(n), cols.shape[1])
        return sparse.csr_matrix((vals.ravel(), (rows, cols.ravel())), shape=(n, self.n_coef))

    def _amplitudes(self, s, w):
        num = np.bincount(self.curve_codes, w * self.f * s, minlength=self.n_curves)
        den = np.bincount(self.curve_codes, w * s * s, minlength=self.n_curves)
        return np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)

    def _solve_template(self, XM, Xc0, a, w):
        """Free template coefficients u given the amplitudes: linear least squares."""
        sw = np.sqrt(w)
        B = XM.multiply((a * sw)[:, None]).tocsr()
        y = sw * (self.f - a * Xc0)
        lhs = (B.T @ B).toarray() + self.Pu
        rhs = B.T @ y - self.Pc0
        return np.linalg.solve(lhs, rhs)

    def _state(self, u, theta, phi, w):
        """Model quantities for template u and dilation parameters theta.

        D = exp(phi @ theta) per point, so for phi = ln(1 + z) the single
        parameter is b, and for phi = redshift-bin indicators each parameter is
        the log of that bin's dilation factor.
        """
        tau = self.t * np.exp(-(phi @ theta)) if phi is not None else self._tau_fixed
        X = self.design(tau)
        XM = (X @ self.M).tocsr()
        Xc0 = X @ self.c0
        s = Xc0 + XM @ u
        amp = self._amplitudes(s, w)
        model = amp[self.curve_codes] * s
        chi2 = float(np.sum(w * (self.f - model) ** 2))
        c = self.c0 + self.M @ u
        return {"tau": tau, "X": X, "XM": XM, "Xc0": Xc0, "s": s, "amp": amp,
                "model": model, "chi2": chi2, "obj": chi2 + float(c @ self.penalty @ c), "c": c}

    def _gauss_newton_step(self, st, u, w, phi, free_theta, damping):
        """One Gauss-Newton step with the amplitudes eliminated.

        The amplitudes enter linearly, so they are profiled out (variable
        projection): the Schur complement of the amplitude block gives a
        system in the template coefficients (and dilation parameters) alone.
        """
        sw = np.sqrt(w)
        a = st["amp"][self.curve_codes]
        r = sw * (self.f - a * st["s"])
        B = st["XM"].multiply((sw * a)[:, None]).tocsr()       # -d r / d u
        n = len(self.f)
        G = sparse.csr_matrix((sw * st["s"], (np.arange(n), self.curve_codes)),
                              shape=(n, self.n_curves))         # -d r / d amp
        H = (B.T @ B).toarray() + self.Pu
        g = B.T @ r - self.M.T @ (self.penalty @ st["c"])      # minus half the gradient
        HcA = (B.T @ G).toarray()
        if free_theta:
            h = 1e-3
            dsdtau = ((self.design(st["tau"] + h) - self.design(st["tau"] - h)) @ st["c"]) / (2 * h)
            J = (sw * a * dsdtau * (-st["tau"]))[:, None] * phi   # -d r / d theta
            k = phi.shape[1]
            BJ = B.T @ J
            H = np.block([[H, BJ], [BJ.T, J.T @ J]])
            g = np.concatenate([g, J.T @ r])
            HcA = np.vstack([HcA, np.asarray(G.T @ J).T])
        D = np.bincount(self.curve_codes, w * st["s"] ** 2, minlength=self.n_curves)
        Dinv = np.where(D > 0, 1.0 / np.where(D > 0, D, 1.0), 0.0)
        gA = G.T @ r
        schur = H - (HcA * Dinv[None, :]) @ HcA.T
        rhs = g - HcA @ (Dinv * gA)
        step = np.linalg.solve(schur + damping * np.diag(np.diag(schur)), rhs)
        self._last_schur = schur
        return step

    def fit(self, dilation=None, *, phi=None, theta0=None, free_theta=False, weights=None,
            max_iter=100, tol=1e-12, u0=None):
        """Best template and amplitudes, and optionally the dilation parameters.

        Either give ``dilation`` (one fixed D per point), or ``phi`` (n_points x
        k) with ``theta0``: D = exp(phi @ theta). With ``free_theta`` the
        parameters theta are fitted together with the template.
        ``weights`` multiplies each point's chi-square term (bootstrap
        multiplicities; 0 removes a point). Returns chi2, theta and, when
        theta is free, its covariance from the curvature of the chi-square.
        """
        if dilation is not None:
            self._tau_fixed = self.t / dilation
            phi, theta = None, np.zeros(0)
            free_theta = False
        else:
            theta = np.asarray(theta0, float).copy()
        mult = np.ones_like(self.f) if weights is None else np.asarray(weights, float)
        w = mult / self.sigma2
        if u0 is not None:
            u = u0
        elif self.u is not None:
            u = self.u
        else:
            st = self._state(np.zeros(self.free.size), theta, phi, w)
            amp0 = pd.Series(self.f).groupby(self.curve_codes).max().to_numpy()
            u = self._solve_template(st["XM"], st["Xc0"], np.where(amp0 > 0, amp0, 1.0)[self.curve_codes], w)
        for _ in range(3):
            st = self._state(u, theta, phi, w)
            u = self._solve_template(st["XM"], st["Xc0"], st["amp"][self.curve_codes], w)
        st = self._state(u, theta, phi, w)
        damping = 1e-4
        nu = self.free.size
        for it in range(max_iter):
            step = self._gauss_newton_step(st, u, w, phi, free_theta, damping)
            u_t = u + step[:nu]
            th_t = theta + step[nu:] if free_theta else theta
            st_t = self._state(u_t, th_t, phi, w)
            if st_t["obj"] <= st["obj"]:
                converged = st["obj"] - st_t["obj"] < tol * st["obj"]
                u, theta, st = u_t, th_t, st_t
                damping = max(damping / 10.0, 1e-10)
                if converged:
                    break
            else:
                damping *= 10.0
                if damping > 1e6:
                    break
        self.u = u
        self.coef = st["c"]
        used = mult > 0
        out = {"chi2": st["chi2"], "objective": st["obj"], "n_points": int(used.sum()),
               "iterations": it + 1, "n_curves": int(np.unique(self.curve_codes[used]).size),
               "theta": theta, "coef": self.coef, "amp": st["amp"], "model": st["model"],
               "tau": st["tau"]}
        if free_theta:
            self._gauss_newton_step(st, u, w, phi, True, 0.0)
            cov = np.linalg.inv(self._last_schur)
            out["theta_cov"] = cov[nu:, nu:]
        return out

    def reset(self):
        self.coef = None
        self.u = None

    # -- the tests ------------------------------------------------------------

    def phi_power(self):
        """Basis for D = (1 + z)^b: ln D = b ln(1 + z)."""
        return np.log1p(self.z)[:, None]

    def phi_bins(self, edges):
        """Basis for one dilation factor per redshift bin (the first bin fixed at D = 1).

        The template's own time scale absorbs an overall factor, so only the
        ratios D_k / D_1 are measurable.
        """
        k = np.digitize(self.z, edges) - 1
        return np.stack([(k == j).astype(float) for j in range(1, len(edges) - 1)], axis=1)

    def chi2_power(self, b, **kw):
        """Minimum chi-square with D = (1 + z)^b held fixed."""
        return self.fit(phi=self.phi_power(), theta0=[b], **kw)["chi2"]

    def fit_power(self, b0=1.0, **kw):
        """b fitted together with the template; error from the chi-square curvature."""
        r = self.fit(phi=self.phi_power(), theta0=[b0], free_theta=True, **kw)
        r["b"] = float(r["theta"][0])
        r["b_err_curvature"] = float(np.sqrt(r["theta_cov"][0, 0]))
        return r

    def template_grid(self, tau, lam_nm):
        """Evaluate the current template S on a grid of (tau, rest wavelength)."""
        T = _basis(np.asarray(tau, float), self.tau_knots).toarray()
        L = _basis(np.log10(np.asarray(lam_nm, float)), self.lam_knots).toarray()
        return T @ self.coef.reshape(self.n_tau, self.n_lam) @ L.T
