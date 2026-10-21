"""Fitting Keplerian orbits to radial velocities from several instruments.

The model for a velocity measured at time t by instrument k is

    v(t) = sum over planets of Keplerian(t)          (nonlinear parameters)
         + gamma_k                                    (one zero point per instrument)
         + optional slow trend and activity terms     (linear parameters)

with independent Gaussian noise of variance ``err^2 + s_k^2``. The quoted
error ``err`` is photon noise; ``s_k`` ("jitter") is the extra scatter every
star and instrument adds on top -- granulation, oscillations, spots,
calibration -- which must be *fitted*, or the error bars lie.

Parameters that enter linearly (zero points, trends, activity coefficients)
are integrated out analytically: with flat priors the integral over a linear
Gaussian model is exact, so the sampler only explores the orbits and jitters.
Their values are recovered afterwards from their conditional Gaussian.

Orbits are parameterised by (P, K, sqrt(e) cos w, sqrt(e) sin w, T_conj): the
middle pair gives uniform priors on e and w without piling up at e = 0, and
the time of conjunction keeps the phase meaningful for circular orbits.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field, replace

import numpy as np
from scipy import optimize

import rv_physics as phys
from config import CACHE_DIR, E_MAX, JITTER_MAX

K_MAX = 3000.0  # m/s, upper bound of the uniform prior on semi-amplitudes


# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------

@dataclass
class RVData:
    """Velocities of one star from one or more instruments, time-ordered."""

    star: str
    t: np.ndarray  # BJD (TDB)
    y: np.ndarray  # m/s
    err: np.ndarray  # m/s, quoted (photon-noise) uncertainty
    inst: np.ndarray  # integer instrument index into inst_names
    inst_names: list[str]
    activity: dict[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self):
        order = np.argsort(self.t, kind="stable")
        self.t = np.asarray(self.t, float)[order]
        self.y = np.asarray(self.y, float)[order]
        self.err = np.asarray(self.err, float)[order]
        self.inst = np.asarray(self.inst, int)[order]
        self.activity = {k: np.asarray(v, float)[order] for k, v in self.activity.items()}

    @property
    def n(self) -> int:
        return self.t.size

    @property
    def n_inst(self) -> int:
        return len(self.inst_names)

    @property
    def baseline(self) -> float:
        return float(np.ptp(self.t))

    def counts(self) -> dict[str, int]:
        return {name: int(np.sum(self.inst == k)) for k, name in enumerate(self.inst_names)}

    def subset(self, mask) -> "RVData":
        """Keep a subset of points; instruments left empty are dropped."""
        mask = np.asarray(mask, bool)
        used = [k for k in range(self.n_inst) if np.any(mask & (self.inst == k))]
        remap = {old: new for new, old in enumerate(used)}
        return RVData(
            self.star, self.t[mask], self.y[mask], self.err[mask],
            np.array([remap[k] for k in self.inst[mask]], int),
            [self.inst_names[k] for k in used],
            {k: v[mask] for k, v in self.activity.items()},
        )

    def only(self, *names: str) -> "RVData":
        keep = [self.inst_names.index(n) for n in names if n in self.inst_names]
        return self.subset(np.isin(self.inst, keep))

    def with_values(self, y) -> "RVData":
        return replace(self, y=np.asarray(y, float), activity=dict(self.activity))

    def fingerprint(self) -> str:
        h = hashlib.sha256()
        for arr in (self.t, self.y, self.err, self.inst):
            h.update(np.ascontiguousarray(arr).tobytes())
        h.update(json.dumps(self.inst_names).encode())
        return h.hexdigest()[:16]


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

@dataclass
class PlanetSpec:
    """Starting guess and prior range for one planet."""

    P: float
    K: float
    tc: float
    h: float = 0.0
    k: float = 0.0
    circular: bool = False
    P_bounds: tuple[float, float] | None = None
    label: str = ""

    def bounds(self) -> tuple[float, float]:
        return self.P_bounds if self.P_bounds else (self.P / 1.5, self.P * 1.5)


#: Points closer in time than this (days), from the same instrument, belong to
#: the same night.
NIGHT_GAP_D = 0.4
#: An instrument gets a night-to-night noise term when more than this fraction
#: of its points share a night with another point.
NIGHT_TERM_MIN_FRACTION = 0.5


def night_groups(t, inst, which):
    """Group the points of the instruments ``which`` into nights.

    Returns (order, starts, group_inst): point indices sorted so each night is
    contiguous, the start of each night within ``order``, and each night's
    instrument.
    """
    pts = np.flatnonzero(np.isin(inst, which))
    if pts.size == 0:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0, int)
    order = pts[np.lexsort((t[pts], inst[pts]))]
    new = np.r_[True, (np.diff(inst[order]) != 0) | (np.diff(t[order]) > NIGHT_GAP_D)]
    starts = np.flatnonzero(new)
    return order, starts, inst[order][starts]


def instruments_needing_night_term(data: RVData) -> list[str]:
    out = []
    for k, name in enumerate(data.inst_names):
        t = np.sort(data.t[data.inst == k])
        if t.size < 3:
            continue
        g = np.r_[0, np.cumsum(np.diff(t) > NIGHT_GAP_D)]
        sizes = np.bincount(g)
        if np.sum(sizes[sizes > 1]) / t.size > NIGHT_TERM_MIN_FRACTION:
            out.append(name)
    return out


class KeplerModel:
    """Keplerian orbits + noise terms (sampled) + linear terms (marginalised).

    Noise: each instrument has a white "jitter" added in quadrature to the
    quoted errors. Instruments that took many exposures per night also get a
    *night-to-night* term: points from the same night share a random offset
    of variance ``sigma_night^2``. Stellar granulation and supergranulation,
    and nightly calibration drifts, move all of a night's velocities together;
    without this term, forty exposures on one night would count as forty
    independent measurements and a few nights would dominate the fit. The
    covariance is then block-diagonal, and its inverse and determinant have
    closed forms (Sherman-Morrison), so the likelihood stays exact and fast.

    Linear terms: per-instrument offsets, an optional polynomial trend in
    time, and optional linear decorrelation against activity indicators.
    """

    def __init__(self, data: RVData, planets: list[PlanetSpec], *, trend: int = 0,
                 activity: tuple[str, ...] = (), t_ref: float | None = None,
                 fixed_jitter: dict[str, float] | None = None,
                 night_term: list[str] | str | None = "auto"):
        self.data = data
        self.planets = [replace(p) for p in planets]
        self.trend = trend
        self.activity = tuple(activity)
        self.t_ref = float(np.round(np.median(data.t))) if t_ref is None else float(t_ref)
        self.fixed_jitter = dict(fixed_jitter or {})
        if night_term == "auto":
            night_term = instruments_needing_night_term(data)
        self.night_insts = list(night_term or [])

        names, x0, lo, hi = [], [], [], []
        for j, p in enumerate(self.planets, start=1):
            pb = p.bounds()
            entries = [(f"P{j}", p.P, *pb), (f"K{j}", p.K, 0.0, K_MAX)]
            if not p.circular:
                entries += [(f"h{j}", p.h, -1.0, 1.0), (f"k{j}", p.k, -1.0, 1.0)]
            entries += [(f"tc{j}", p.tc, p.tc - 0.5 * p.P, p.tc + 0.5 * p.P)]
            for n_, v, a, b in entries:
                names.append(n_); x0.append(v); lo.append(a); hi.append(b)
        self.jitter_names = []
        for name in data.inst_names:
            if name in self.fixed_jitter:
                continue
            names.append(f"jit[{name}]")
            x0.append(3.0); lo.append(0.0); hi.append(JITTER_MAX)
            self.jitter_names.append(name)
        for name in self.night_insts:
            names.append(f"night[{name}]")
            x0.append(1.0); lo.append(0.0); hi.append(JITTER_MAX)
        idx = [data.inst_names.index(n) for n in self.night_insts]
        self._g_order, self._g_starts, self._g_inst = night_groups(data.t, data.inst, idx)
        self.names = names
        self.x0 = np.array(x0, float)
        self.lo = np.array(lo, float)
        self.hi = np.array(hi, float)
        self.ndim = len(names)
        self._index = {n_: i for i, n_ in enumerate(names)}
        self.X, self.linear_names = self._design()

    # -- linear part -----------------------------------------------------------
    def _design(self):
        d = self.data
        cols, names = [], []
        for k, name in enumerate(d.inst_names):
            cols.append((d.inst == k).astype(float)); names.append(f"gamma[{name}]")
        tau = (d.t - self.t_ref) / 365.25
        for p in range(1, self.trend + 1):
            cols.append(tau ** p); names.append(f"trend_{p} (m/s/yr^{p})")
        for a in self.activity:
            v = d.activity[a].copy()
            if not np.all(np.isfinite(v)):
                raise ValueError(f"activity index {a} has missing values")
            # Centre within each instrument so the coefficient measures the
            # velocity response to *changes* in activity, not zero points.
            for k in range(d.n_inst):
                m = d.inst == k
                v[m] -= v[m].mean()
            cols.append(v); names.append(f"c[{a}]")
        return np.column_stack(cols), names

    # -- parameters ------------------------------------------------------------
    def _col(self, theta, name, default=None):
        i = self._index.get(name)
        if i is None:
            return np.full(theta.shape[0], default, float)
        return theta[:, i]

    def planet_arrays(self, theta):
        """(W, n_planets) arrays of P, K, e, omega, tc for W parameter vectors."""
        theta = np.atleast_2d(theta)
        P, K, e, w, tc = [], [], [], [], []
        for j, _ in enumerate(self.planets, start=1):
            h = self._col(theta, f"h{j}", 0.0)
            k = self._col(theta, f"k{j}", 0.0)
            ej, wj = phys.he_to_eomega(h, k)
            # The prior rejects e >= E_MAX; clipping here only keeps an
            # optimiser that strays outside it from evaluating e >= 1.
            ej = np.minimum(ej, E_MAX)
            P.append(theta[:, self._index[f"P{j}"]]); K.append(theta[:, self._index[f"K{j}"]])
            e.append(ej); w.append(wj); tc.append(theta[:, self._index[f"tc{j}"]])
        stack = lambda xs: np.stack(xs, axis=1) if xs else np.zeros((theta.shape[0], 0))
        return stack(P), stack(K), stack(e), stack(w), stack(tc)

    def jitter(self, theta):
        """(W, n_inst) jitter for each instrument."""
        theta = np.atleast_2d(theta)
        out = np.empty((theta.shape[0], self.data.n_inst))
        for k, name in enumerate(self.data.inst_names):
            if name in self.fixed_jitter:
                out[:, k] = self.fixed_jitter[name]
            else:
                out[:, k] = theta[:, self._index[f"jit[{name}]"]]
        return out

    def keplerian(self, theta, t=None, planets=None):
        """Sum of Keplerian velocities, shape (W, len(t))."""
        t = self.data.t if t is None else np.asarray(t, float)
        P, K, e, w, tc = self.planet_arrays(theta)
        out = np.zeros((P.shape[0], t.size))
        for j in range(P.shape[1]):
            if planets is not None and j not in planets:
                continue
            out += phys.rv_keplerian(t[None, :], P[:, j:j + 1], K[:, j:j + 1],
                                     e[:, j:j + 1], w[:, j:j + 1], tc[:, j:j + 1])
        return out

    # -- probability -------------------------------------------------------------
    def log_prior(self, theta):
        theta = np.atleast_2d(theta)
        ok = np.all((theta >= self.lo) & (theta <= self.hi), axis=1)
        for j, p in enumerate(self.planets, start=1):
            if not p.circular:
                ok &= (theta[:, self._index[f"h{j}"]] ** 2 + theta[:, self._index[f"k{j}"]] ** 2) < E_MAX
        return np.where(ok, 0.0, -np.inf)

    def night_sigma(self, theta):
        """(W, n_night_instruments) night-to-night noise."""
        theta = np.atleast_2d(theta)
        if not self.night_insts:
            return np.zeros((theta.shape[0], 0))
        return np.stack([theta[:, self._index[f"night[{n}]"]] for n in self.night_insts], axis=1)

    def _night_variance_per_group(self, theta):
        """(W, n_groups) shared-offset variance of each night."""
        sig = self.night_sigma(theta)
        col = {self.data.inst_names.index(n): j for j, n in enumerate(self.night_insts)}
        return sig[:, [col[k] for k in self._g_inst]] ** 2

    def _linear_solve(self, theta):
        """Generalised least squares for the linear parameters, for each of W
        parameter vectors, under the block-diagonal noise covariance.

        With white variances D and, for each night g, a shared-offset variance
        c_g, the inverse covariance is D^-1 minus, per night, the rank-one term
        f_g (D^-1 1)(D^-1 1)^T with f_g = c_g / (1 + c_g sum_g D^-1), and
        log det C = sum log D + sum_g log(1 + c_g sum_g D^-1).

        Returns residual r (before linear terms), log det C, beta-hat, the
        normal matrix A = X^T C^-1 X, and chi^2 at beta-hat.
        """
        theta = np.atleast_2d(theta)
        d = self.data
        r = d.y[None, :] - self.keplerian(theta)
        var = d.err[None, :] ** 2 + self.jitter(theta)[:, d.inst] ** 2
        u = 1.0 / var
        A = np.einsum("wn,nl,nm->wlm", u, self.X, self.X, optimize=True)
        b = np.einsum("wn,nl->wl", u * r, self.X, optimize=True)
        rr = np.sum(u * r * r, axis=1)
        logdet = np.sum(np.log(var), axis=1)
        if self._g_starts.size:
            o, s = self._g_order, self._g_starts
            c = self._night_variance_per_group(theta)
            uo = u[:, o]
            S = np.add.reduceat(uo, s, axis=1)
            Ur = np.add.reduceat(uo * r[:, o], s, axis=1)
            UX = np.add.reduceat(uo[:, :, None] * self.X[o][None, :, :], s, axis=1)
            f = c / (1.0 + c * S)
            A = A - np.einsum("wg,wgl,wgm->wlm", f, UX, UX, optimize=True)
            b = b - np.einsum("wg,wgl,wg->wl", f, UX, Ur, optimize=True)
            rr = rr - np.sum(f * Ur * Ur, axis=1)
            logdet = logdet + np.sum(np.log1p(c * S), axis=1)
        beta = np.linalg.solve(A, b[..., None])[..., 0]
        chi2 = rr - np.sum(b * beta, axis=1)
        return r, logdet, beta, A, chi2

    def log_like(self, theta, *, marginal: bool = True):
        """Gaussian log-likelihood with linear parameters either integrated out
        (``marginal``, for sampling) or set to their best values (profile)."""
        _, logdet, _, A, chi2 = self._linear_solve(theta)
        ll = -0.5 * (chi2 + logdet + self.data.n * np.log(2 * np.pi))
        if marginal:
            sign, logdet_A = np.linalg.slogdet(A)
            ll = ll + 0.5 * self.X.shape[1] * np.log(2 * np.pi) - 0.5 * logdet_A
        return ll

    def log_prob(self, theta, marginal: bool = True):
        theta = np.atleast_2d(theta)
        lp = self.log_prior(theta)
        out = np.full(theta.shape[0], -np.inf)
        good = np.isfinite(lp)
        if np.any(good):
            out[good] = lp[good] + self.log_like(theta[good], marginal=marginal)
        return out

    # -- convenience for a single parameter vector ---------------------------------
    def linear_params(self, theta):
        _, _, beta, A, _ = self._linear_solve(theta)
        cov = np.linalg.inv(A[0])
        return dict(zip(self.linear_names, beta[0])), dict(zip(self.linear_names, np.sqrt(np.diag(cov))))

    def full_model(self, theta, t=None, inst=None, *, include_activity: bool = True):
        """Keplerians + offsets + trend (+ activity terms) for one theta. For
        times other than the data's, activity terms are omitted."""
        theta = np.atleast_2d(theta)[:1]
        _, _, beta, _, _ = self._linear_solve(theta)
        beta = beta[0]
        if t is None:
            X = self.X if include_activity else self.X[:, :len(self.linear_names) - len(self.activity)]
            b = beta if include_activity else beta[:X.shape[1]]
            return self.keplerian(theta)[0] + X @ b
        t = np.asarray(t, float)
        inst = np.zeros(t.size, int) if inst is None else np.asarray(inst, int)
        lin = beta[inst]
        tau = (t - self.t_ref) / 365.25
        for p in range(1, self.trend + 1):
            lin = lin + beta[self.data.n_inst + p - 1] * tau ** p
        return self.keplerian(theta, t)[0] + lin

    def offsets_and_trend(self, theta):
        """The linear part of the model evaluated at the data (no Keplerians)."""
        theta = np.atleast_2d(theta)[:1]
        _, _, beta, _, _ = self._linear_solve(theta)
        return self.X @ beta[0]

    def residuals(self, theta):
        return self.data.y - self.full_model(theta)

    def whitened_residuals(self, theta):
        """Residuals z with z.z = r^T C^-1 r exactly (linear terms at their
        best values), for least-squares optimisation of the orbits.

        Each night block C = D + c 1 1^T is whitened by W = (I - g v v^T) D^-1/2
        with v = D^-1/2 1 and g = (1 - sqrt(1 - f |v|^2)) / |v|^2, since then
        W^T W = C^-1.
        """
        theta = np.atleast_2d(theta)[:1]
        r, _, beta, _, _ = self._linear_solve(theta)
        res = r[0] - self.X @ beta[0]
        d = self.data
        sig = np.sqrt(d.err ** 2 + self.jitter(theta)[0, d.inst] ** 2)
        z = res / sig
        if self._g_starts.size:
            o, s = self._g_order, self._g_starts
            c = self._night_variance_per_group(theta)[0]
            u = 1.0 / sig[o] ** 2
            S = np.add.reduceat(u, s)
            Ur = np.add.reduceat(u * res[o], s)
            f = c / (1.0 + c * S)
            g = np.where(S > 0, (1.0 - np.sqrt(np.clip(1.0 - f * S, 0.0, None))) / S, 0.0)
            sizes = np.diff(np.r_[s, o.size])
            z[o] -= np.repeat(g * Ur, sizes) / sig[o]
        return z

    def effective_errors(self, theta):
        """Per-point errors for periodogram weights. Quoted error and jitter in
        quadrature, plus, for a point in a night of n points with night noise
        s, a share n s^2: the night's n points then carry together the weight
        of one measurement of variance (white/n + s^2), as they should."""
        d = self.data
        var = d.err ** 2 + self.jitter(theta)[0, d.inst] ** 2
        if self._g_starts.size:
            c = self._night_variance_per_group(theta)[0]
            sizes = np.diff(np.r_[self._g_starts, self._g_order.size])
            extra = np.repeat(c * sizes, sizes)
            var[self._g_order] += extra
        return np.sqrt(var)

    def n_params(self) -> int:
        return self.ndim + self.X.shape[1]

    def n_effective(self) -> int:
        """Number of independent measurements: each night of an instrument
        with a night-to-night noise term counts once."""
        return int(self.data.n - self._g_order.size + self._g_starts.size)


# ---------------------------------------------------------------------------
# Maximum a-posteriori fit
# ---------------------------------------------------------------------------

@dataclass
class FitResult:
    model: KeplerModel
    theta: np.ndarray
    log_like: float  # profile (maximised) log-likelihood
    success: bool

    @property
    def bic(self) -> float:
        return self.model.n_params() * np.log(self.model.data.n) - 2 * self.log_like

    def params(self) -> dict[str, float]:
        return dict(zip(self.model.names, self.theta))

    def planets(self) -> list[dict]:
        P, K, e, w, tc = (a[0] for a in self.model.planet_arrays(self.theta))
        return [dict(P=float(P[j]), K=float(K[j]), e=float(e[j]), omega=float(w[j]), tc=float(tc[j]))
                for j in range(P.size)]

    def rms(self) -> float:
        return float(np.std(self.model.residuals(self.theta)))


def _initial_jitter(model: KeplerModel, theta):
    """Per-instrument excess scatter of the residuals over the quoted errors."""
    d = model.data
    r = d.y - model.keplerian(theta)[0]
    out = theta.copy()
    for k, name in enumerate(d.inst_names):
        key = f"jit[{name}]"
        if key not in model._index:
            continue
        m = d.inst == k
        rk = r[m] - np.average(r[m], weights=1 / d.err[m] ** 2)
        excess = np.var(rk) - np.mean(d.err[m] ** 2) if m.sum() > 2 else 1.0
        out[model._index[key]] = float(np.clip(np.sqrt(max(excess, 0.25)), 0.5, JITTER_MAX * 0.9))
    return out


#: Bump when fit_map's algorithm changes, so cached optima are recomputed.
FIT_CACHE_VERSION = 2


def _map_cache_path(model: KeplerModel, theta0, n_rounds: int, tol: float):
    spec = {
        "data": model.data.fingerprint(), "X": hashlib.sha256(np.ascontiguousarray(model.X).tobytes()).hexdigest(),
        "names": model.names, "lo": np.round(model.lo, 10).tolist(), "hi": np.round(model.hi, 10).tolist(),
        "x0": np.round(model.x0, 10).tolist(),
        "theta0": None if theta0 is None else np.round(np.asarray(theta0, float), 10).tolist(),
        "trend": model.trend, "activity": model.activity, "fixed": model.fixed_jitter,
        "night": model.night_insts, "t_ref": model.t_ref, "rounds": n_rounds, "tol": tol,
        "version": FIT_CACHE_VERSION,
    }
    key = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:20]
    d = CACHE_DIR / "map"
    d.mkdir(exist_ok=True)
    return d / f"{key}.npy"


def fit_map(model: KeplerModel, theta0=None, *, n_rounds: int = 12, tol: float = 1e-7,
            cache: bool = True) -> FitResult:
    """Maximum a-posteriori fit, cached on disk (the optimisation is
    deterministic, so the same model, data and start give the same answer)."""
    path = _map_cache_path(model, theta0, n_rounds, tol) if cache else None
    if path is not None:
        _note_cache_use(path)
    if path is not None and path.exists():
        theta = np.load(path)
        ll = float(model.log_like(theta, marginal=False)[0])
        return FitResult(model, theta, ll, bool(np.isfinite(ll)))
    res = _fit_map(model, theta0, n_rounds=n_rounds, tol=tol)
    if path is not None and res.success:
        np.save(path, res.theta)
    return res


def _fit_map(model: KeplerModel, theta0=None, *, n_rounds: int = 12, tol: float = 1e-7) -> FitResult:
    """Maximise the posterior (profile likelihood x priors).

    Alternates two well-conditioned steps until the log-likelihood stops
    improving: (1) a Levenberg-Marquardt least-squares fit of the orbits with
    the noise held fixed, on exactly whitened residuals, so it minimises the
    true chi^2 of the correlated-noise model; (2) a bounded quasi-Newton fit
    of the noise parameters (all in m/s, so well scaled) with the orbits held
    fixed. A joint quasi-Newton step would have to handle a period known to
    one part in a million next to jitters known to 10%, and stalls.
    """
    theta = np.array(model.x0 if theta0 is None else theta0, float)
    theta = np.clip(theta, model.lo + 1e-9, model.hi - 1e-9)
    if theta0 is None:
        theta = _initial_jitter(model, theta)
    noise_idx = [i for i, n in enumerate(model.names) if n.startswith(("jit[", "night["))]
    orb_idx = [i for i in range(model.ndim) if i not in noise_idx]
    neg = lambda th: -model.log_prob(th, marginal=False)[0]
    best = neg(theta)
    success = np.isfinite(best)
    for _ in range(n_rounds):
        if orb_idx:
            def resid(x):
                th = theta.copy(); th[orb_idx] = x
                return model.whitened_residuals(th)
            lo = model.lo[orb_idx]; hi = model.hi[orb_idx]
            x0 = np.clip(theta[orb_idx], lo + 1e-12, hi - 1e-12)
            # Tolerances of 1e-10 are far below any quoted uncertainty; tighter
            # ones make Levenberg-Marquardt crawl for thousands of iterations.
            sol = optimize.least_squares(resid, x0, bounds=(lo, hi), x_scale="jac",
                                         max_nfev=600, ftol=1e-10, xtol=1e-10, gtol=1e-10)
            trial = theta.copy(); trial[orb_idx] = sol.x
            if neg(trial) <= best:
                theta = trial; best = neg(theta)
        if noise_idx:
            # (a) Per instrument, a coarse grid over (jitter, night noise)
            # with everything else fixed. When most nights hold one or two
            # points the two terms trade off along a ridge with more than one
            # local maximum, and a local optimiser keeps whichever it met first.
            theta = _grid_noise(model, theta, neg)
            # (b) Polish jointly in sqrt(sigma): sigma >= 0 then holds
            # automatically and the gradient stays finite at sigma = 0, where
            # a bounded optimiser in sigma itself stalls.
            hi_n = model.hi[noise_idx]

            def negn(x):
                th = theta.copy(); th[noise_idx] = np.minimum(x * x, hi_n)
                return neg(th)
            for start in (np.sqrt(np.maximum(theta[noise_idx], 0.0)),
                          np.full(len(noise_idx), 1.0)):
                sol = optimize.minimize(negn, start, method="L-BFGS-B",
                                        options={"maxiter": 3000, "ftol": 1e-13, "gtol": 1e-9})
                if np.isfinite(sol.fun) and sol.fun < neg(theta):
                    theta[noise_idx] = np.minimum(sol.x ** 2, hi_n)
        new = neg(theta)
        converged = best - new < tol
        best = min(best, new)
        if converged:
            break
    ll = float(model.log_like(theta, marginal=False)[0])
    return FitResult(model, theta, ll, bool(success and np.isfinite(ll)))


NOISE_GRID = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 35.0])


def _grid_noise(model: KeplerModel, theta, neg):
    """Coarse grid search of each instrument's noise terms (see fit_map)."""
    theta = theta.copy()
    for name in model.data.inst_names:
        idx = [model._index[n] for n in (f"jit[{name}]", f"night[{name}]") if n in model._index]
        if not idx:
            continue
        grid = [NOISE_GRID] * len(idx)
        mesh = np.array(np.meshgrid(*grid, indexing="ij")).reshape(len(idx), -1).T
        trials = np.repeat(theta[None, :], len(mesh) + 1, axis=0)
        trials[1:, idx] = mesh  # row 0 keeps the current values
        vals = -model.log_prob(trials, marginal=False)
        theta = trials[int(np.argmin(vals))]
    return theta


def fit_best_of(model_factory, starts, **kwargs) -> FitResult:
    """Fit from several starting points (e.g. eccentricity vectors) and keep
    the best. ``model_factory(start)`` returns (model, theta0)."""
    best = None
    for start in starts:
        model, theta0 = model_factory(start)
        res = fit_map(model, theta0, **kwargs)
        if best is None or res.log_like > best.log_like:
            best = res
    return best


def laplace_errors(model: KeplerModel, theta, *, rel_step: float = 1e-4):
    """Approximate 1-sigma errors from the curvature of the log posterior at
    its maximum (finite-difference Hessian). Quick, and adequate for strong
    signals; the MCMC is the real answer."""
    theta = np.asarray(theta, float)
    n = theta.size
    step = np.maximum(np.abs(theta) * rel_step, 1e-6)
    for i, name in enumerate(model.names):
        if name.startswith("P"):
            step[i] = max(theta[i] * 1e-7, 1e-8)
        if name.startswith("tc"):
            step[i] = 1e-4
    f = lambda th: model.log_prob(th, marginal=False)[0]
    f0 = f(theta)
    H = np.empty((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ei[i] = step[i]
            ej = np.zeros(n); ej[j] = step[j]
            fpp = f(theta + ei + ej); fpm = f(theta + ei - ej)
            fmp = f(theta - ei + ej); fmm = f(theta - ei - ej)
            H[i, j] = H[j, i] = (fpp - fpm - fmp + fmm) / (4 * step[i] * step[j])
    try:
        cov = np.linalg.inv(-H)
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        cov, err = None, np.full(n, np.nan)
    if not np.isfinite(f0):
        err[:] = np.nan
    return err, cov


# ---------------------------------------------------------------------------
# Periodogram-driven planet search
# ---------------------------------------------------------------------------

def planet_from_sinusoid(P, amp, phase, t_ref) -> PlanetSpec:
    """Initial circular orbit from a fitted sinusoid ``A cos(2 pi t/P - phi)``.

    A circular orbit with our conventions is ``-K sin(2 pi (t - tc)/P)``, so
    ``tc = P (phi + pi/2) / (2 pi)``, moved to the orbit nearest ``t_ref``.
    """
    tc = P * (phase + np.pi / 2) / (2 * np.pi)
    tc += P * np.round((t_ref - tc) / P)
    return PlanetSpec(P=float(P), K=float(max(amp, 0.1)), tc=float(tc), circular=True)


@dataclass
class SearchStep:
    n_planets_before: int
    freqs: np.ndarray
    power: np.ndarray
    peak_period: float
    peak_power: float
    fap: float
    top_periods: list[float]
    top_powers: list[float]
    rms_before: float


def residual_periodogram(fit: FitResult, freqs, *, multi_offset: bool = True):
    """GLS of the residuals of a fit, weighted by the fitted noise (quoted error
    and jitter), with instrument zero points refitted at each frequency."""
    m = fit.model
    r = m.residuals(fit.theta)
    err = m.effective_errors(fit.theta)
    inst = m.data.inst if multi_offset else None
    return gls(m.data.t, r, err, freqs, inst, return_fit=True)


def gls(*args, **kwargs):  # re-exported for convenience
    return phys.gls(*args, **kwargs)


def search_planets(data: RVData, freqs, *, n_max: int, fap_threshold: float,
                   trend: int = 0, activity: tuple[str, ...] = (), fixed_jitter=None,
                   exclude=None, P_bounds_factor: float = 1.02, quiet: bool = False):
    """Find planets one at a time: periodogram of the residuals -> highest peak
    -> add a circular orbit there -> refit *everything* jointly -> repeat,
    until the highest peak is no longer significant.

    ``exclude`` is an optional list of (P_lo, P_hi) ranges whose peaks are
    skipped (used only for the explicit alias experiment).
    """
    planets: list[PlanetSpec] = []
    model = KeplerModel(data, planets, trend=trend, activity=activity, fixed_jitter=fixed_jitter)
    fit = fit_map(model)
    steps: list[SearchStep] = []
    for i in range(n_max + 1):
        power, amp, phase = residual_periodogram(fit, freqs)
        use = np.ones(freqs.size, bool)
        for lo, hi in (exclude or []):
            use &= ~((1 / freqs >= lo) & (1 / freqs <= hi))
        pw = np.where(use, power, 0.0)
        top = phys.top_peaks(freqs, pw, n=6)
        best = top[0]
        err_eff = fit.model.effective_errors(fit.theta)
        fap = float(phys.fap_baluev(pw[best], data.t, err_eff, freqs.max(), n_offsets=data.n_inst,
                                    n_eff=fit.model.n_effective()))
        steps.append(SearchStep(len(planets), freqs, power, float(1 / freqs[best]), float(pw[best]),
                                fap, [float(1 / freqs[j]) for j in top], [float(pw[j]) for j in top],
                                fit.rms()))
        if not quiet:
            print(f"    step {i}: {len(planets)} planets, rms {fit.rms():.2f} m/s -> peak "
                  f"{1 / freqs[best]:.4f} d (power {pw[best]:.3f}, FAP {fap:.1e})")
        if fap > fap_threshold or i == n_max:
            break
        P = 1 / freqs[best]
        new = planet_from_sinusoid(P, amp[best], phase[best], model.t_ref)
        new.P_bounds = (P / P_bounds_factor, P * P_bounds_factor) if P < data.baseline / 3 else (P / 2, P * 2)
        planets = [replace(p, P=q["P"], K=q["K"], tc=q["tc"]) for p, q in zip(planets, fit.planets())] + [new]
        model = KeplerModel(data, planets, trend=trend, activity=activity, fixed_jitter=fixed_jitter)
        fit = fit_map(model)
    return steps, fit


def signal_by_instrument(fit: FitResult, period: float) -> dict:
    """Fit a sinusoid at ``period`` to each instrument's residuals separately.

    A planet moves the star the same way whoever is watching, so a real signal
    has the same amplitude and phase in every instrument. The (a, b)
    coefficients of ``a cos + b sin`` are compared with a chi-square test
    against their weighted mean.
    """
    from scipy import stats

    m = fit.model
    d = m.data
    r = m.residuals(fit.theta)
    err = m.effective_errors(fit.theta)
    rows, vecs, covs = [], [], []
    for k, name in enumerate(d.inst_names):
        sel = d.inst == k
        if sel.sum() < 8:
            continue
        ph = 2 * np.pi * d.t[sel] / period
        X = np.column_stack([np.cos(ph), np.sin(ph), np.ones(sel.sum())])
        W = 1.0 / err[sel] ** 2
        A = X.T @ (W[:, None] * X)
        cov = np.linalg.inv(A)
        beta = cov @ (X.T @ (W * r[sel]))
        amp = float(np.hypot(beta[0], beta[1]))
        # amplitude error by linear propagation
        g = np.array([beta[0], beta[1]]) / max(amp, 1e-12)
        amp_err = float(np.sqrt(g @ cov[:2, :2] @ g))
        rows.append({"instrument": name, "n": int(sel.sum()), "amplitude": amp, "amplitude_err": amp_err})
        vecs.append(beta[:2]); covs.append(cov[:2, :2])
    if len(vecs) < 2:
        return {"rows": rows, "chi2": None, "dof": 0, "p_consistent": None}
    inv = [np.linalg.inv(c) for c in covs]
    mean = np.linalg.solve(sum(inv), sum(i @ v for i, v in zip(inv, vecs)))
    chi2 = float(sum((v - mean) @ i @ (v - mean) for v, i in zip(vecs, inv)))
    dof = 2 * (len(vecs) - 1)
    return {"rows": rows, "chi2": chi2, "dof": dof, "p_consistent": float(stats.chi2.sf(chi2, dof)),
            "mean_amplitude": float(np.hypot(*mean))}


# ---------------------------------------------------------------------------
# MCMC
# ---------------------------------------------------------------------------

@dataclass
class Posterior:
    model: KeplerModel
    samples: np.ndarray  # (n_samples, ndim)
    log_prob: np.ndarray
    acceptance: float
    autocorr: np.ndarray | None
    n_steps: int

    def column(self, name):
        return self.samples[:, self.model.names.index(name)]

    def planet_samples(self):
        return self.model.planet_arrays(self.samples)

    def linear_samples(self, rng, chunk: int = 4000):
        """Draw the marginalised linear parameters from their conditional
        Gaussian given each sample."""
        out = []
        for i0 in range(0, len(self.samples), chunk):
            th = self.samples[i0:i0 + chunk]
            _, _, beta, A, _ = self.model._linear_solve(th)
            L = np.linalg.cholesky(np.linalg.inv(A))
            z = rng.standard_normal(beta.shape)
            out.append(beta + np.einsum("wij,wj->wi", L, z))
        return np.concatenate(out)


def _initial_scales(model: KeplerModel, theta) -> np.ndarray:
    """Spread of the walkers' starting ball, parameter by parameter: a small
    fraction of each quantity's natural scale (1% of a periodogram peak's
    width in period, a thousandth of an orbit in phase). A ball this small
    stays inside the best fit's own peak -- for 55 Cnc e the neighbouring
    period sidelobes are only 6e-5 d apart -- and the ensemble expands to
    the posterior's true width within the burn-in."""
    T = model.data.baseline
    out = np.empty(theta.size)
    for i, n in enumerate(model.names):
        v = theta[i]
        if n.startswith("P"):
            out[i] = 0.01 * v * v / T
        elif n.startswith("K"):
            out[i] = max(0.01 * v, 0.01)
        elif n.startswith(("h", "k")) and n[1:].isdigit():
            out[i] = 0.005
        elif n.startswith("tc"):
            j = n[2:]
            out[i] = 1e-3 * theta[model.names.index(f"P{j}")]
        else:  # jitter, night noise
            out[i] = max(0.01 * v, 0.01)
    return out


def _note_cache_use(path) -> None:
    """If RV_CACHE_LOG is set, record every cache file a run touches, so that
    files left over from earlier versions can be pruned."""
    import os
    log = os.environ.get("RV_CACHE_LOG")
    if log:
        with open(log, "a") as fh:
            fh.write(f"{path}\n")


def _cache_key(model: KeplerModel, tag: str, n_steps: int, n_walkers: int, seed: int) -> str:
    spec = {
        "data": model.data.fingerprint(), "names": model.names, "lo": model.lo.tolist(),
        "hi": model.hi.tolist(), "trend": model.trend, "activity": model.activity,
        "fixed": model.fixed_jitter, "steps": n_steps, "walkers": n_walkers, "seed": seed,
        "t_ref": model.t_ref, "night": model.night_insts, "version": 6,
    }
    return tag + "_" + hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:12]


def run_mcmc(fit: FitResult, *, tag: str, n_steps: int, burn: int, thin: int,
             walkers_per_dim: int, seed: int, quiet: bool = False) -> Posterior:
    """Sample the posterior with emcee (differential-evolution moves, which
    cope better than the default stretch move with correlated, many-
    dimensional posteriors). Cached on disk by a hash of data + model."""
    import emcee

    model = fit.model
    ndim = model.ndim
    n_walkers = max(2 * ndim + 2, walkers_per_dim * ndim)
    key = _cache_key(model, tag, n_steps, n_walkers, seed)
    path = CACHE_DIR / f"{key}.npz"
    _note_cache_use(path)
    if path.exists():
        z = np.load(path)
        return Posterior(model, z["samples"], z["log_prob"], float(z["acceptance"]),
                         z["autocorr"] if z["autocorr"].size else None, int(z["n_steps"]))

    rng = np.random.default_rng(seed)
    scale = _initial_scales(model, fit.theta)
    # Start the walkers in a small ball around the best fit. Draws are
    # clipped into the prior box (a jitter of exactly zero sits on its
    # bound), and the ball shrinks whenever draws keep landing outside the
    # priors (e.g. beyond e = E_MAX) -- with 30+ parameters a fixed ball can
    # almost never produce a valid walker.
    lo = model.lo + 1e-9 * (model.hi - model.lo)
    hi = model.hi - 1e-9 * (model.hi - model.lo)
    p0, tries, s = [], 0, scale
    while len(p0) < n_walkers:
        trial = np.clip(fit.theta + s * rng.standard_normal(ndim), lo, hi)
        if np.isfinite(model.log_prob(trial)[0]):
            p0.append(trial)
        tries += 1
        if tries % 100 == 0:
            s = s * 0.5
        if tries > 20000:
            raise RuntimeError(f"could not initialise MCMC walkers for {tag}")
    p0 = np.array(p0)
    moves = [(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)]
    sampler = emcee.EnsembleSampler(n_walkers, ndim, lambda th: model.log_prob(th),
                                    vectorize=True, moves=moves)
    sampler.random_state = np.random.RandomState(int(rng.integers(2**31))).get_state()
    if not quiet:
        print(f"    MCMC {tag}: {n_walkers} walkers x {n_steps} steps, {ndim} parameters")
    sampler.run_mcmc(p0, n_steps, progress=False)
    chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
    lp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    try:
        tau = sampler.get_autocorr_time(discard=burn, quiet=True)
    except Exception:
        tau = np.array([])
    acc = float(np.mean(sampler.acceptance_fraction))
    np.savez_compressed(path, samples=chain, log_prob=lp, acceptance=acc,
                        autocorr=tau, n_steps=n_steps)
    return Posterior(model, chain, lp, acc, tau if tau.size else None, n_steps)


def summarize(x, *, upper_limit: float | None = None) -> dict:
    """Median and 68% interval of a sample; optionally a one-sided upper limit."""
    x = np.asarray(x, float)
    lo, med, hi = np.percentile(x, [15.865, 50, 84.135])
    out = {"median": float(med), "minus": float(med - lo), "plus": float(hi - med)}
    if upper_limit is not None:
        out[f"upper_{int(round(upper_limit * 100))}"] = float(np.percentile(x, 100 * upper_limit))
    return out
