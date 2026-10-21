"""How many Higgs bosons, at what mass, and how sure can we be?

An unbinned maximum-likelihood fit to the four-lepton masses of the selected
events, channel by channel (4mu, 4e, 2e2mu), in 105 < m4l < 160 GeV:

    density = mu * S_c * f_sig(m; m_H) + theta_ZZ * B_c * f_ZZ(m) + theta_ZX * X_c * f_ZX(m)

* **Signal** shape: a double-sided Crystal Ball fitted to the simulated Higgs
  boson (generated at 125 GeV) in each channel, slid along the mass axis to a
  trial mass m_H and widened by the data/simulation resolution ratio measured
  on the Z peaks (``calibration.py``). mu is the signal strength: 1 is the
  Standard Model rate as normalised by the outreach analysis.
* **ZZ** background: the simulated ZZ events (smoothed), with the outreach
  normalisation and a 10% uncertainty (nuisance theta_ZZ).
* **Z+X** (misidentified leptons; neglected by the outreach analysis): shape
  and size from the same-sign control region, assuming as many opposite-sign
  as same-sign fakes, with a 50% uncertainty (nuisance theta_ZX).

Test statistic for a discovery: q0 = -2 ln [L(mu=0) / L(mu_hat)], with
mu_hat >= 0 (Cowan, Cranmer, Gross & Vitells 2011, Eur. Phys. J. C 71, 1554).
The local significance is sqrt(q0); its distribution is checked with
pseudo-experiments, and the "look-elsewhere effect" of scanning m_H over
110-160 GeV is estimated with the upcrossing method of Gross & Vitells (2010,
Eur. Phys. J. C 70, 525).
"""

from __future__ import annotations

import numpy as np
from iminuit import Minuit
from scipy import stats
from scipy.optimize import brentq

import fourlepton as fl

FIT_RANGE = (105.0, 160.0)
SCAN = np.arange(110.0, 160.0 + 1e-9, 0.25)
MC_MH = 125.0
GRID = np.linspace(FIT_RANGE[0], FIT_RANGE[1], 2201)
CHANNELS = fl.CHANNELS


# ---------------------------------------------------------------------------
# Shapes
# ---------------------------------------------------------------------------


def dscb(x, mu, s, al, nl, ar, nr):
    """Double-sided Crystal Ball (unnormalised)."""
    t = (np.asarray(x, dtype=float) - mu) / s
    out = np.exp(-0.5 * t * t)
    lo, hi = t < -al, t > ar
    Al = (nl / al) ** nl * np.exp(-0.5 * al * al)
    Ar = (nr / ar) ** nr * np.exp(-0.5 * ar * ar)
    out[lo] = Al * (nl / al - al - t[lo]) ** (-nl)
    out[hi] = Ar * (nr / ar - ar + t[hi]) ** (-nr)
    return out


def fit_signal_shape(m: np.ndarray) -> dict:
    """Unbinned fit of a DSCB to simulated Higgs masses within a wide window."""
    lo, hi = 100.0, 150.0
    m = m[(m > lo) & (m < hi)]
    xs = np.linspace(lo, hi, 4001)

    def nll(mu, s, al, nl, ar, nr):
        norm = np.trapezoid(dscb(xs, mu, s, al, nl, ar, nr), xs)
        return -np.sum(np.log(dscb(m, mu, s, al, nl, ar, nr) / norm + 1e-300))

    mn = Minuit(nll, mu=124.8, s=1.5, al=1.0, nl=3.0, ar=1.5, nr=5.0)
    mn.errordef = Minuit.LIKELIHOOD
    mn.limits = [(120, 130), (0.3, 5), (0.2, 5), (1.01, 50), (0.2, 5), (1.01, 50)]
    mn.migrad()
    return dict(zip(["mu", "s", "al", "nl", "ar", "nr"], np.array(mn.values).tolist())) | {
        "n_events": int(len(m))}


def kde(samples: np.ndarray, weights: np.ndarray, bandwidth: float, x: np.ndarray):
    """Weighted Gaussian kernel density, evaluated on x (not normalised)."""
    d = (x[:, None] - samples[None, :]) / bandwidth
    return (np.exp(-0.5 * d * d) * weights[None, :]).sum(axis=1)


class Model:
    """The per-channel densities on a fine mass grid, interpolated at the data."""

    def __init__(self, mc: dict, ss_data: dict, ss_mc: dict, *, resolution_ratio=None,
                 zz_sigma: float = 0.10, zx_sigma: float = 0.50, os_ss_ratio: float = 1.0):
        lo, hi = FIT_RANGE
        self.resolution_ratio = resolution_ratio or {c: 1.0 for c in CHANNELS}
        self.sig_par, self.S, self.B, self.X = {}, {}, {}, {}
        self.f_zz, self.f_zx = {}, {}
        H, ZZ = mc["H"], mc["ZZ"]
        for c in CHANNELS:
            sel = H["channel"] == c
            self.sig_par[c] = fit_signal_shape(H["m4l"][sel])
            win = sel & (H["m4l"] > lo) & (H["m4l"] < hi)
            self.S[c] = float(H["weight"][win].sum())
            zsel = ZZ["channel"] == c
            mz = ZZ["m4l"][zsel]
            wz = ZZ["weight"][zsel]
            near = (mz > lo - 15) & (mz < hi + 15)
            dens = kde(mz[near], wz[near], 3.0, GRID)
            self.f_zz[c] = dens / np.trapezoid(dens, GRID)
            self.B[c] = float(wz[(mz > lo) & (mz < hi)].sum())
        # Z+X: one shape for all channels (a Landau-like Moyal distribution
        # fitted to all same-sign data events), per-channel yields from counts
        m_ss = ss_data["m4l"]
        loc, scale = stats.moyal.fit(m_ss[(m_ss > 70) & (m_ss < 600)])
        self.zx_shape = {"loc": float(loc), "scale": float(scale)}
        frac_in = stats.moyal.cdf(hi, loc, scale) - stats.moyal.cdf(lo, loc, scale)
        dens = stats.moyal.pdf(GRID, loc, scale)
        self.n_ss, self.n_ss_mc, self.os_ss_ratio = {}, {}, os_ss_ratio
        for c in CHANNELS:
            n_ss = float(np.sum(ss_data["channel"] == c))
            n_ss_mc = float(sum(ss_mc[k]["weight"][ss_mc[k]["channel"] == c].sum()
                                for k in ("ZZ", "H")))
            self.n_ss[c], self.n_ss_mc[c] = n_ss, n_ss_mc
            self.X[c] = max(n_ss - n_ss_mc, 0.0) * os_ss_ratio * frac_in
            self.f_zx[c] = dens / np.trapezoid(dens, GRID)
        self.zz_sigma, self.zx_sigma = zz_sigma, zx_sigma

    def f_sig(self, c: str, mh: float, x=GRID):
        p = dict(self.sig_par[c])
        p.pop("n_events", None)
        p["mu"] += mh - MC_MH
        p["s"] *= self.resolution_ratio[c]
        d = dscb(GRID, **p)
        d = d / np.trapezoid(d, GRID)
        return d if x is GRID else np.interp(x, GRID, d)

    def expected(self) -> dict:
        return {c: {"S": self.S[c], "ZZ": self.B[c], "ZX": self.X[c]} for c in CHANNELS}


# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------


class Likelihood:
    def __init__(self, model: Model, data: dict, *, constrain: bool = True):
        lo, hi = FIT_RANGE
        self.model = model
        self.constrain = constrain
        self.m = {c: data["m4l"][(data["channel"] == c) & (data["m4l"] > lo)
                                 & (data["m4l"] < hi)] for c in CHANNELS}
        self.bz = {c: np.interp(self.m[c], GRID, model.f_zz[c]) for c in CHANNELS}
        self.bx = {c: np.interp(self.m[c], GRID, model.f_zx[c]) for c in CHANNELS}
        self._sig_cache = {}

    def sig(self, c, mh):
        key = (c, round(mh, 6))
        if key not in self._sig_cache:
            self._sig_cache[key] = self.model.f_sig(c, mh, self.m[c])
        return self._sig_cache[key]

    def nll(self, mu, mh, tz, tx):
        """-2 ln L (extended), with log-normal constraints on the nuisances."""
        M = self.model
        tot = 0.0
        for c in CHANNELS:
            nu = mu * M.S[c] + tz * M.B[c] + tx * M.X[c]
            dens = mu * M.S[c] * self.sig(c, mh) + tz * M.B[c] * self.bz[c] + tx * M.X[c] * self.bx[c]
            tot += 2 * (nu - np.sum(np.log(np.maximum(dens, 1e-300))))
        if self.constrain:
            tot += (np.log(tz) / M.zz_sigma) ** 2 + (np.log(tx) / M.zx_sigma) ** 2
        return tot

    def fit(self, *, mu=None, mh=None):
        """Minimise over the free parameters; fix mu and/or mh if given."""
        m = Minuit(lambda a, b, c, d: self.nll(a, b, c, d), a=1.0 if mu is None else mu,
                   b=125.0 if mh is None else mh, c=1.0, d=1.0)
        m.errordef = 1.0
        m.limits["a"] = (0.0, 20.0)
        m.limits["b"] = (SCAN[0], SCAN[-1])
        m.limits["c"] = (0.3, 3.0)
        m.limits["d"] = (0.05, 10.0)
        if mu is not None:
            m.fixed["a"] = True
        if mh is not None:
            m.fixed["b"] = True
        if mu is not None and mu == 0.0:
            m.fixed["b"] = True       # m_H is undefined without signal
        m.migrad()
        return {"nll": m.fval, "mu": m.values["a"], "mh": m.values["b"], "tz": m.values["c"],
                "tx": m.values["d"], "minuit": m}

    def scan_q0(self, masses=SCAN) -> np.ndarray:
        b = self.fit(mu=0.0, mh=125.0)["nll"]
        out = []
        for mh in masses:
            s = self.fit(mh=mh)
            out.append(max(b - s["nll"], 0.0) if s["mu"] > 0 else 0.0)
        return np.array(out)

    def profile_mass(self, masses) -> np.ndarray:
        return np.array([self.fit(mh=mh)["nll"] for mh in masses])

    def profile_mu(self, mus) -> np.ndarray:
        return np.array([self.fit(mu=mu)["nll"] for mu in mus])


def interval(x, nll2, level=1.0):
    """Points where -2 Delta ln L crosses `level` around the minimum (linear interp.)."""
    d = nll2 - nll2.min()
    i = int(np.argmin(d))
    lo = hi = np.nan
    for j in range(i, 0, -1):
        if d[j - 1] >= level > d[j]:
            lo = np.interp(level, [d[j], d[j - 1]], [x[j], x[j - 1]])
            break
    for j in range(i, len(d) - 1):
        if d[j + 1] >= level > d[j]:
            hi = np.interp(level, [d[j], d[j + 1]], [x[j], x[j + 1]])
            break
    return x[i], lo, hi


# ---------------------------------------------------------------------------
# Pseudo-experiments: background-only toys with fixed nuisances
# ---------------------------------------------------------------------------


def toy_q0_scans(model: Model, n_toys: int, seed: int, masses=SCAN, mu_sig: float = 0.0,
                 mh_sig: float = 125.0) -> np.ndarray:
    """q0(m_H) for each toy, profiling mu >= 0 with the nuisances fixed at their
    nominal values (adequate for the trial-factor estimate).

    Each toy draws Poisson event counts per channel and masses from the model;
    q0 at every trial mass comes from a few vectorised Newton steps in mu.
    """
    rng = np.random.default_rng(seed)
    sig_dens = {c: np.stack([model.f_sig(c, mh) for mh in masses]) for c in CHANNELS}
    out = np.zeros((n_toys, len(masses)))
    for t in range(n_toys):
        grad_terms = []
        for c in CHANNELS:
            B = model.B[c] + model.X[c]
            fb = (model.B[c] * model.f_zz[c] + model.X[c] * model.f_zx[c]) / max(B, 1e-12)
            n_b = rng.poisson(B)
            mb = _sample(rng, GRID, fb, n_b)
            ms = np.zeros(0)
            if mu_sig > 0:
                n_s = rng.poisson(mu_sig * model.S[c])
                ms = _sample(rng, GRID, model.f_sig(c, mh_sig), n_s)
            m = np.concatenate([mb, ms])
            r = np.stack([np.interp(m, GRID, sd) for sd in sig_dens[c]]) * model.S[c]
            b = np.interp(m, GRID, fb) * B
            grad_terms.append((model.S[c], r, b))
        out[t] = _q0_newton(grad_terms)
    return out


def _sample(rng, grid, dens, n):
    cdf = np.cumsum(dens)
    cdf = cdf / cdf[-1]
    return np.interp(rng.random(n), cdf, grid)


def _q0_newton(terms, iters: int = 30) -> np.ndarray:
    """max over mu >= 0 of 2[ln L(mu) - ln L(0)] for each trial mass."""
    nm = terms[0][1].shape[0]
    mu = np.zeros(nm)
    for _ in range(iters):
        g = np.zeros(nm)
        h = np.zeros(nm)
        for S, r, b in terms:
            den = mu[:, None] * r + b[None, :]
            g += -S + np.sum(r / den, axis=1)
            h += -np.sum((r / den) ** 2, axis=1)
        step = -g / np.minimum(h, -1e-12)
        mu = np.maximum(mu + step, 0.0)
    q = np.zeros(nm)
    for S, r, b in terms:
        q += 2 * (-mu * S + np.sum(np.log1p(mu[:, None] * r / b[None, :]), axis=1))
    return np.maximum(q, 0.0)


def gross_vitells(q_obs: float, toy_scans: np.ndarray, u0: float = 1.0) -> dict:
    """Global p-value from the mean number of upcrossings of level u0."""
    up = np.sum((toy_scans[:, 1:] > u0) & (toy_scans[:, :-1] <= u0), axis=1)
    n_up = float(np.mean(up))
    p_local = 0.5 * stats.chi2.sf(q_obs, 1)
    p_global = p_local + n_up * np.exp(-(q_obs - u0) / 2)
    return {"u0": u0, "mean_upcrossings": n_up, "mean_upcrossings_err": float(np.std(up) /
                                                                            np.sqrt(len(up))),
            "p_local": float(p_local), "p_global": float(min(p_global, 1.0)),
            "z_local": float(stats.norm.isf(p_local)),
            "z_global": float(stats.norm.isf(min(p_global, 0.5)))}


def z_from_p(p: float) -> float:
    return float(stats.norm.isf(p))


def counting_significance(n_obs: int, b: float) -> dict:
    """Poisson p-value of n_obs given background b (one-sided)."""
    p = float(stats.poisson.sf(n_obs - 1, b))
    return {"n_obs": n_obs, "b": b, "p": p, "z": z_from_p(p)}


def mass_at(level_fn, lo, hi):
    return brentq(level_fn, lo, hi)
