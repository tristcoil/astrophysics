"""The alpha-knee: where Type Ia supernovae start to show.

Massive stars explode within ~3-40 Myr of forming and return alpha elements
(Mg, Si, Ca, ...) and iron in a fixed ratio, so early stars sit on a
*plateau* of [Mg/Fe]. Type Ia supernovae add iron after a delay of ~0.1 Gyr to
many Gyr; once they do, [Mg/Fe] falls as [Fe/H] keeps rising. The [Fe/H] at
the turn -- the knee -- says how much iron a galaxy had made before its first
Type Ia went off, i.e. how fast it formed stars.

Two populations are measured the same way:

* the Milky Way's high-alpha (thick) disc, as the ridge line of its sequence:
  in each [Fe/H] bin a one- or two-Gaussian fit to [Mg/Fe] (chosen by BIC),
  the upper component's mean;
* Gaia-Sausage-Enceladus, star by star.

Both are fitted with a plateau that breaks into a straight decline at
[Fe/H]_knee, with intrinsic scatter, by maximum likelihood on a grid of knee
positions; uncertainties come from bootstrap resampling of the stars.
"""

from __future__ import annotations

import numpy as np

KNEE_GRID = np.round(np.arange(-2.4, 0.01, 0.01), 3)


# ---------------------------------------------------------------------------
# Plateau + decline
# ---------------------------------------------------------------------------


def _design(x, xk):
    return np.column_stack([np.ones_like(x), np.clip(x - xk, 0.0, None)])


def broken_line_fit(x, y, e=None, *, grid=KNEE_GRID, sigma_int=None, n_iter=3):
    """Maximum-likelihood plateau + decline. Returns knee, plateau, slope,
    intrinsic scatter and the log-likelihood profile over the knee grid.

    The likelihood is Gaussian with variance e^2 + s^2; s (the intrinsic
    scatter) is re-estimated a few times and then held fixed while the knee is
    scanned, so the profile compares knees at a common noise level.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    e = np.zeros_like(x) if e is None else np.asarray(e, float)
    grid = grid[(grid > np.percentile(x, 2)) & (grid < np.percentile(x, 98))]
    s = np.std(y) if sigma_int is None else sigma_int
    for _ in range(n_iter if sigma_int is None else 1):
        w = 1.0 / (e ** 2 + s ** 2)
        chi = np.empty(grid.size)
        coefs = []
        for i, xk in enumerate(grid):
            A = _design(x, xk)
            coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
            r = y - A @ coef
            chi[i] = np.sum(w * r ** 2)
            coefs.append(coef)
        best = int(np.argmin(chi))
        r = y - _design(x, grid[best]) @ coefs[best]
        if sigma_int is None:
            s = np.sqrt(max(np.mean(r ** 2 - e ** 2), 1e-6))
    lnl = -0.5 * chi
    # A plain straight line (no plateau), for the significance of the knee.
    A = np.column_stack([np.ones_like(x), x])
    w = 1.0 / (e ** 2 + s ** 2)
    coef, *_ = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)
    chi_line = np.sum(w * (y - A @ coef) ** 2)
    return dict(knee=float(grid[best]), plateau=float(coefs[best][0]),
                slope=float(coefs[best][1]), sigma_int=float(s), grid=grid, lnl=lnl,
                delta_chi2_vs_line=float(chi_line - chi[best]), n=int(x.size))


def bootstrap_knee(x, y, e=None, *, n_boot=200, seed=0):
    """Knee and plateau with bootstrap errors (stars resampled)."""
    rng = np.random.default_rng(seed)
    best = broken_line_fit(x, y, e)
    knees, plats, slopes = [], [], []
    for _ in range(n_boot):
        i = rng.integers(0, x.size, x.size)
        f = broken_line_fit(x[i], y[i], None if e is None else e[i], sigma_int=best["sigma_int"])
        knees.append(f["knee"]), plats.append(f["plateau"]), slopes.append(f["slope"])
    knees = np.array(knees)
    best.update(
        knee_lo=float(np.percentile(knees, 16)), knee_hi=float(np.percentile(knees, 84)),
        knee_p2=float(np.percentile(knees, 2.5)), knee_p97=float(np.percentile(knees, 97.5)),
        plateau_err=float(np.std(plats)), slope_err=float(np.std(slopes)), boot_knees=knees)
    return best


# ---------------------------------------------------------------------------
# The high-alpha sequence's ridge line
# ---------------------------------------------------------------------------


def _gmm1d(y, k, n_iter=200, seed=0):
    """Plain EM for a 1-D Gaussian mixture; returns weights, means, sigmas, lnL."""
    rng = np.random.default_rng(seed)
    mu = np.sort(rng.choice(y, k, replace=False)) if k > 1 else np.array([y.mean()])
    if k == 2:
        mu = np.percentile(y, [25, 85])
    sig = np.full(k, y.std() / k + 1e-3)
    w = np.full(k, 1.0 / k)
    for _ in range(n_iter):
        p = w * np.exp(-0.5 * ((y[:, None] - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        tot = p.sum(1, keepdims=True) + 1e-300
        r = p / tot
        nk = r.sum(0) + 1e-12
        w, mu = nk / y.size, (r * y[:, None]).sum(0) / nk
        sig = np.sqrt((r * (y[:, None] - mu) ** 2).sum(0) / nk) + 1e-4
    lnl = float(np.sum(np.log(tot)))
    order = np.argsort(mu)
    return w[order], mu[order], sig[order], lnl


def high_alpha_ridge(feh, xfe, *, edges=np.arange(-1.6, 0.41, 0.1), min_n=40,
                     min_sep=0.08, n_boot=30, seed=0):
    """Ridge line of the high-alpha sequence, bin by bin in [Fe/H].

    A second component is accepted only if BIC prefers it, it carries >= 5% of
    the stars, and the two means differ by > min_sep (otherwise one sequence
    is simply being split in two). The ridge is the upper mean.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (feh >= lo) & (feh < hi) & np.isfinite(xfe)
        y = xfe[s]
        if y.size < min_n:
            continue

        def ridge(yy):
            w1, m1, s1, l1 = _gmm1d(yy, 1)
            w2, m2, s2, l2 = _gmm1d(yy, 2)
            bic1 = -2 * l1 + 2 * np.log(yy.size)
            bic2 = -2 * l2 + 5 * np.log(yy.size)
            two = bic2 < bic1 - 10 and w2.min() > 0.05 and (m2[1] - m2[0]) > min_sep
            return (m2[1], s2[1], w2[1], m2[0], True) if two else (m1[0], s1[0], 1.0, np.nan, False)

        m, sg, wt, mlow, two = ridge(y)
        boots = [ridge(rng.choice(y, y.size))[0] for _ in range(n_boot)]
        rows.append(dict(feh=float(0.5 * (lo + hi)), feh_median=float(np.median(feh[s])),
                         ridge=float(m), ridge_err=float(np.std(boots)), sigma=float(sg),
                         high_weight=float(wt), low_mean=float(mlow), two_components=bool(two),
                         n=int(y.size)))
    return rows


# ---------------------------------------------------------------------------
# Both disc sequences at once
# ---------------------------------------------------------------------------


def _two_seq_unpack(p, xk):
    plat, slope, ls_h, c0, c1, c2, ls_l, x0, lw = p
    return dict(knee=xk, plateau=plat, slope=slope, sigma_high=np.exp(ls_h), low_c=(c0, c1, c2),
                sigma_low=np.exp(ls_l), w_x0=x0, w_width=np.exp(lw))


def _two_seq_nll(p, xk, x, y, e2):
    """Negative log-likelihood and its analytic gradient."""
    plat, slope, ls_h, c0, c1, c2, ls_l, x0, lw = p
    s_h2, s_l2, W = np.exp(2 * ls_h), np.exp(2 * ls_l), np.exp(lw)
    c = np.clip(x - xk, 0.0, None)
    mu_h = plat + slope * c
    mu_l = c0 + c1 * x + c2 * x ** 2
    v_h, v_l = s_h2 + e2, s_l2 + e2
    u = np.clip((x - x0) / W, -50, 50)
    w = 1.0 / (1.0 + np.exp(u))                      # high-alpha share
    dh, dl = y - mu_h, y - mu_l
    ph = np.exp(-0.5 * dh ** 2 / v_h) / np.sqrt(v_h)
    pl = np.exp(-0.5 * dl ** 2 / v_l) / np.sqrt(v_l)
    L = w * ph + (1 - w) * pl + 1e-300
    rh, rl = w * ph / L, (1 - w) * pl / L
    gh, gl = rh * dh / v_h, rl * dl / v_l            # d lnL / d mu
    vh, vl = rh * (dh ** 2 / v_h ** 2 - 1 / v_h) / 2, rl * (dl ** 2 / v_l ** 2 - 1 / v_l) / 2
    gw = (ph - pl) / L * (-w * (1 - w))              # d lnL / d u
    grad = np.array([gh.sum(), (gh * c).sum(), (vh * 2 * s_h2).sum(), gl.sum(), (gl * x).sum(),
                     (gl * x ** 2).sum(), (vl * 2 * s_l2).sum(), (gw * (-1 / W)).sum(),
                     (gw * (-u)).sum()])
    return -np.sum(np.log(L)), -grad


#: Bounds keep the mixture physical: both sequences 0.01-0.15 dex wide, the
#: hand-over from high- to low-alpha somewhere in -1 < [Fe/H] < 0.3 over a
#: width of 0.03-0.6 dex.
_BOUNDS = [(0.0, 0.6), (-2.0, 0.5), (np.log(0.01), np.log(0.15)), (-0.3, 0.3), (-1.0, 1.0),
           (-1.0, 1.0), (np.log(0.01), np.log(0.15)), (-1.0, 0.3), (np.log(0.03), np.log(0.6))]


def _start(x, y):
    """Starting point from the data: the plateau from metal-poor stars, the
    low-alpha sequence from a quadratic through metal-rich medians."""
    plat = float(np.median(y[(x > -1.3) & (x < -0.9)]))
    rich = x > -0.1
    c2, c1, c0 = np.polyfit(x[rich], y[rich], 2)
    return np.array([plat, -0.4, np.log(0.04), c0, c1, c2, np.log(0.04), -0.4, np.log(0.12)])


def two_sequence_fit(feh, xfe, err, *, coarse=np.round(np.arange(-1.3, -0.19, 0.05), 3),
                     fine_step=0.005):
    """The high- and low-alpha disc sequences fitted together star by star.

    p([X/Fe] | [Fe/H]) = w N(plateau + decline beyond the knee, s_h)
                       + (1 - w) N(quadratic low-alpha sequence, s_l),
    with the high-alpha share w a logistic function of [Fe/H] and each star's
    measurement error added in quadrature. The knee is profiled on a coarse
    grid, then a fine one around the best value, each point started both from
    the data-driven guess and from its neighbour's solution (the better one
    kept: the likelihood has local optima). The 1-sigma range is where the
    profile log-likelihood drops by 0.5.
    """
    from scipy.optimize import minimize

    x, y = np.asarray(feh, float), np.asarray(xfe, float)
    e2 = np.asarray(err, float) ** 2
    p0 = _start(x, y)

    def fit_at(xk, starts):
        best = None
        for p in starts:
            r = minimize(_two_seq_nll, p, args=(xk, x, y, e2), jac=True, method="L-BFGS-B",
                         bounds=_BOUNDS)
            if best is None or r.fun < best.fun:
                best = r
        return best

    def profile(grid):
        prof, pars, prev = [], [], None
        for xk in grid:
            r = fit_at(xk, [p0] if prev is None else [p0, prev])
            prev = r.x
            prof.append(-r.fun)
            pars.append(r.x.copy())
        return np.array(prof), pars

    prof_c, pars_c = profile(coarse)
    i = int(np.argmax(prof_c))
    fine = np.round(np.arange(coarse[i] - 0.06, coarse[i] + 0.0601, fine_step), 4)
    prof, pars = profile(fine)
    j = int(np.argmax(prof))
    ok = prof >= prof[j] - 0.5
    out = _two_seq_unpack(pars[j], float(fine[j]))
    out.update(knee_lo=float(fine[ok].min()), knee_hi=float(fine[ok].max()),
               grid_coarse=coarse, grid_fine=fine, profile_coarse=prof_c - prof[j],
               profile_fine=prof - prof[j], n=int(x.size))
    out["low_c"] = [float(c) for c in out["low_c"]]
    return {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in out.items()}


def high_alpha_probability(fit, feh, xfe, err):
    """Posterior probability that each star belongs to the high-alpha sequence."""
    x, y, e2 = np.asarray(feh, float), np.asarray(xfe, float), np.asarray(err, float) ** 2
    mu_h = fit["plateau"] + fit["slope"] * np.clip(x - fit["knee"], 0, None)
    c0, c1, c2 = fit["low_c"]
    mu_l = c0 + c1 * x + c2 * x ** 2
    v_h, v_l = fit["sigma_high"] ** 2 + e2, fit["sigma_low"] ** 2 + e2
    w = 1.0 / (1.0 + np.exp(np.clip((x - fit["w_x0"]) / fit["w_width"], -50, 50)))
    ph = w * np.exp(-0.5 * (y - mu_h) ** 2 / v_h) / np.sqrt(v_h)
    pl = (1 - w) * np.exp(-0.5 * (y - mu_l) ** 2 / v_l) / np.sqrt(v_l)
    return ph / (ph + pl + 1e-300)


def alpha_division(feh, xfe, *, edges=np.arange(-0.8, 0.21, 0.1), min_sep=0.08):
    """Straight line through the valley between the high- and low-alpha
    sequences: in each [Fe/H] bin where two components are resolved, the
    [Mg/Fe] at which they are equally likely. Returns (intercept, slope)."""
    xs, ys = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (feh >= lo) & (feh < hi) & np.isfinite(xfe)
        if s.sum() < 200:
            continue
        w, m, sg, _ = _gmm1d(xfe[s], 2)
        if (m[1] - m[0]) < min_sep or w.min() < 0.05:
            continue
        yy = np.linspace(m[0], m[1], 400)
        p0 = w[0] * np.exp(-0.5 * ((yy - m[0]) / sg[0]) ** 2) / sg[0]
        p1 = w[1] * np.exp(-0.5 * ((yy - m[1]) / sg[1]) ** 2) / sg[1]
        xs.append(float(np.median(feh[s])))
        ys.append(float(yy[np.argmin(np.abs(np.log(p0 / p1)))]))
    slope, intercept = np.polyfit(xs, ys, 1)
    return float(intercept), float(slope), np.array(xs), np.array(ys)


def knee_from_ridge(rows, *, n_boot=500, seed=0, max_jump=0.08):
    """Plateau + decline fitted to the ridge points (weighted by their errors);
    the knee's uncertainty from resampling the points within their errors.
    The ridge is followed only until it jumps down by more than ``max_jump``
    from one bin to the next: there the fit has lost the high-alpha sequence
    under the far more numerous low-alpha stars."""
    keep = [rows[0]]
    for r in rows[1:]:
        if keep[-1]["ridge"] - r["ridge"] > max_jump:
            break
        keep.append(r)
    rows = keep
    x = np.array([r["feh_median"] for r in rows])
    y = np.array([r["ridge"] for r in rows])
    e = np.array([max(r["ridge_err"], 0.003) for r in rows])
    best = broken_line_fit(x, y, e, grid=np.round(np.arange(-1.6, 0.31, 0.005), 3))
    rng = np.random.default_rng(seed)
    knees = [broken_line_fit(x, y + rng.normal(0, e), e, sigma_int=best["sigma_int"],
                             grid=np.round(np.arange(-1.6, 0.31, 0.005), 3))["knee"]
             for _ in range(n_boot)]
    best.update(knee_lo=float(np.percentile(knees, 16)), knee_hi=float(np.percentile(knees, 84)))
    return best
