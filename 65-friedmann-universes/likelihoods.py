"""Scoring curved universes against supernovae, BAO and the CMB.

Every likelihood returns -2 ln L (a chi^2, up to a constant) for a *batch* of
universes, so a whole grid of the (Omega_m, Omega_Lambda) plane is scored at once.

* **Supernovae (Pantheon+)**: standardised magnitudes with their full 1590 x 1590
  covariance; the unknown offset (supernova luminosity combined with H0) is
  removed exactly, as in D1.
* **BAO (DESI DR2)**: 13 distance ratios D_X / r_d. Alone, BAO measure only
  ratios of distances, so the ruler's size (h r_d) is a free parameter,
  marginalised with DESI's flat prior. With the CMB, r_d is computed.
* **CMB (compressed Planck PR4)**: a Gaussian in (theta_*, omega_b,
  omega_b + omega_c). The CMB measures the physical densities directly; a
  universe's curvature and dark energy then only change theta_* through the
  distance to the last-scattering surface.

**Marginalising the physical densities.** In the plane, Omega_m and Omega_Lambda
are fixed, and the nuisance parameters x = (omega_b, omega_b + omega_c) are not
(h follows from Omega_m h^2 = omega_m). The CMB pins x through the prior, and
theta_* -- measured to 0.025% -- pins one combination of x far more tightly
still, so a grid over x would be hopeless. Instead, at every point of the
plane, Gauss-Newton iterations find the best x, and the Gaussian integral over
x is done analytically around it (the Laplace approximation):

    -2 ln L = chi^2(x_best) + ln det F,     F = J^T C^-1 J,

where J is the Jacobian of the predictions with respect to x. The predictions
are smooth and nearly linear over the tiny range of x the data allow, so this is
accurate to well below the precision needed (``test_physics.py`` checks it
against brute-force integration over x).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import cholesky, solve_triangular

import background as bg
import config as cfg


def whitening(cov: np.ndarray) -> np.ndarray:
    """W with W^T W = C^-1 (so chi^2 = |W r|^2): W = L^-1 for C = L L^T."""
    lower = cholesky(cov, lower=True)
    return solve_triangular(lower, np.eye(cov.shape[0]), lower=True)


def _chunks(n: int, size: int):
    for start in range(0, n, size):
        yield slice(start, min(n, start + size))


# ---------------------------------------------------------------------------
# Supernovae
# ---------------------------------------------------------------------------


class SupernovaLikelihood:
    """Pantheon+ with the offset profiled (identical to marginalising it with a
    flat prior: the model is linear in it, with the same width everywhere)."""

    def __init__(self, data: dict):
        self.z_hd, self.z_hel, self.m = data["z_hd"], data["z_hel"], data["m"]
        self.W = whitening(data["cov"])
        self.u = self.W @ np.ones(self.m.size)
        self.uu = float(self.u @ self.u)
        self.wm = self.W @ self.m
        self.dist = bg.SupernovaDistances(self.z_hd)
        self.n = self.m.size

    def chi2(self, model: bg.Model, chunk: int = 1500) -> tuple[np.ndarray, np.ndarray]:
        """Profiled chi^2 and best offset M + 25 + 5 log10(c/H0 stuff) per model."""
        chi2 = np.empty(model.n)
        off = np.empty(model.n)
        for sl in _chunks(model.n, chunk):
            m = model.take(sl)
            dm = self.dist.transverse(m)
            with np.errstate(invalid="ignore", divide="ignore"):
                mu = 5.0 * np.log10((1.0 + self.z_hel)[None, :] * dm) + 25.0
            bad = ~np.all(np.isfinite(mu), axis=1)
            mu[bad] = 0.0
            r = self.wm[:, None] - self.W @ mu.T
            ur = self.u @ r
            c = np.einsum("ij,ij->j", r, r) - ur**2 / self.uu
            c[bad] = np.inf
            chi2[sl], off[sl] = c, ur / self.uu
        return chi2, off


# ---------------------------------------------------------------------------
# BAO
# ---------------------------------------------------------------------------


class BAOLikelihood:
    """DESI DR2 BAO: a Gaussian in 13 distance ratios."""

    def __init__(self, data: dict):
        self.z, self.value, self.kind, self.cov = data["z"], data["value"], data["kind"], data["cov"]
        self.W = whitening(self.cov)
        self.uz, self.zi = np.unique(self.z, return_inverse=True)
        self.n = self.value.size

    def distances(self, model: bg.Model) -> np.ndarray:
        """The data vector's distances in Mpc with r_d = 1, shape (n_models, 13)."""
        dc = bg.comoving_distance(self.uz, model, nodes=48)
        dm = bg.transverse_distance(dc, model)[:, self.zi]
        hz = bg.hubble(self.uz, model)[:, self.zi]
        return bg.bao_predictions(self.kind, self.z, dm, hz, np.ones(model.n))

    def chi2_calibrated(self, model: bg.Model, rd) -> np.ndarray:
        pred = self.distances(model) / np.asarray(rd)[:, None]
        r = self.W @ (self.value[None, :] - pred).T
        out = np.einsum("ij,ij->j", r, r)
        return np.where(np.isfinite(out), out, np.inf)

    def m2lnl_uncalibrated(self, model: bg.Model, chunk: int = 20000) -> dict:
        """-2 ln L with the ruler marginalised: D_X/r_d = s D_X, s = 1/r_d.

        The chi^2 is a parabola in s; with DESI's flat prior on h r_d (h fixed,
        so flat in r_d = 1/s, a Jacobian 1/s^2) the integral over s is

            -2 ln L = chi^2_min + ln A + 4 ln s_best + const,  A = m^T C^-1 m.

        (The 1/s^2 factor varies by 2% across the 1% width of s, so its
        curvature correction, ~3 sigma_s^2 / s^2 ~ 3e-4, is negligible.)
        Returns -2 ln L, chi^2 at the best ruler, and the best h r_d.
        """
        out = {k: np.empty(model.n) for k in ("m2lnl", "chi2", "hrd")}
        d = self.W @ self.value
        for sl in _chunks(model.n, chunk):
            m = model.take(sl)
            u = self.W @ self.distances(m).T
            A = np.einsum("ij,ij->j", u, u)
            B = d @ u
            s = B / A
            chi2 = d @ d - B**2 / A
            bad = ~np.isfinite(chi2) | (s <= 0)
            rd = 1.0 / s
            # DESI's prior box on h r_d (far from where the likelihood lives)
            hrd = m.h * rd
            bad |= (hrd < cfg.PRIOR_HRD[0]) | (hrd > cfg.PRIOR_HRD[1])
            m2 = chi2 + np.log(A) + 4.0 * np.log(np.abs(s))
            out["m2lnl"][sl] = np.where(bad, np.inf, m2)
            out["chi2"][sl] = np.where(bad, np.inf, chi2)
            out["hrd"][sl] = hrd
        return out


# ---------------------------------------------------------------------------
# The compressed CMB prior
# ---------------------------------------------------------------------------


@dataclass
class CMBPrior:
    """Gaussian in (theta_*, omega_b, omega_b + omega_c)."""

    name: str
    mean: np.ndarray
    cov: np.ndarray
    note: str = ""
    W: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.mean = np.asarray(self.mean, float)
        self.cov = np.asarray(self.cov, float)
        self.W = whitening(self.cov)

    @property
    def sigma(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))

    def as_dict(self) -> dict:
        s = self.sigma
        return {"name": self.name, "mean": self.mean.tolist(), "sigma": s.tolist(),
                "corr": (self.cov / np.outer(s, s)).tolist(), "note": self.note}


def cmb_prior_from_chains(chains: dict) -> CMBPrior:
    """The compressed prior rebuilt from Lemos & Lewis's chains: weighted mean and
    covariance of (theta_*, omega_b, omega_b + omega_c). The third mean is that of
    omega_b + omega_c itself (D5 found that DESI's version used the mean of
    omega_m h^2, which includes the neutrino's 0.00064)."""
    w = chains["weight"]
    X = np.vstack([chains["thetastar"] / 100.0, chains["ombh2"],
                   chains["ombh2"] + chains["omch2"]])
    return CMBPrior("Planck PR4, compressed", np.average(X, axis=1, weights=w),
                    np.cov(X, aweights=w),
                    note=f"{w.size} chain rows after {cfg.PLANCK_EARLY_BURN:.0%} burn-in")


def cmb_prior_from_d3(d3_fit: dict) -> CMBPrior:
    """Our own CMB prior, from project D3's LCDM fit to the Planck 2018 TT spectrum.

    D3 fitted (omega_b, omega_c, 100 theta_MC, ln A_s, n_s) with tau fixed and
    stored errors and correlations. We keep the (omega_b, omega_c, theta_MC)
    block (marginalising amplitude and tilt) and convert CosmoMC's approximate
    theta_MC into the exact theta_* with CAMB's Jacobian at D3's best fit. (The
    same construction as D5's; recomputed here.)
    """
    import camb

    p = d3_fit["parameters"]
    names = d3_fit["correlation"]["names"]
    err = np.array([p[k]["error"] for k in names])
    cov5 = np.array(d3_fit["correlation"]["matrix"]) * np.outer(err, err)
    idx = [names.index(k) for k in ("ombh2", "omch2", "theta100")]
    cov3 = cov5[np.ix_(idx, idx)]
    x0 = np.array([p["ombh2"]["value"], p["omch2"]["value"], p["theta100"]["value"]])

    def theta_star(x):
        pars = camb.set_params(ombh2=x[0], omch2=x[1], cosmomc_theta=x[2] / 100.0, mnu=0.06,
                               nnu=3.046, tau=d3_fit.get("tau_fixed", 0.0522), As=2.1e-9,
                               ns=0.9649)  # D3's settings (CAMB defaults)
        return camb.get_background(pars).get_derived_params()["thetastar"] / 100.0

    steps = np.sqrt(np.diag(cov3)) * 0.5
    grad = np.array([(theta_star(x0 + steps[i] * np.eye(3)[i])
                      - theta_star(x0 - steps[i] * np.eye(3)[i])) / (2 * steps[i])
                     for i in range(3)])
    J = np.array([grad, [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    mean = np.array([theta_star(x0), x0[0], x0[0] + x0[1]])
    return CMBPrior("D3: our Planck 2018 TT fit", mean, J @ cov3 @ J.T,
                    note="TT only, tau fixed, lensing as in LCDM")


# ---------------------------------------------------------------------------
# The CMB (and calibrated BAO) at fixed (Omega_m, Omega_Lambda)
# ---------------------------------------------------------------------------


class PhysicalDensityMarginal:
    """-2 ln L(Omega_m, Omega_k) with (omega_b, omega_b + omega_c) marginalised.

    ``cmb`` is required (it is what calibrates the physical densities); ``bao``
    adds DESI's distances with the ruler r_d computed from x. Supernovae, which
    do not depend on x, are added separately. The plane is parametrised by
    (Omega_m, Omega_k) because Omega_Lambda = 1 - Omega_m - Omega_k - Omega_r
    depends (very slightly) on h through Omega_r; the map between the two is a
    shear with unit Jacobian, so posteriors carry over unchanged.
    """

    def __init__(self, cmb: CMBPrior, rec: bg.Recombination, bao: BAOLikelihood | None = None,
                 iterations: int = 5):
        self.cmb, self.rec, self.bao = cmb, rec, bao
        self.iterations = iterations
        self.step = 0.02 * cmb.sigma[1:]  # finite-difference steps in x

    def predictions(self, Om, Ok, x) -> tuple[np.ndarray, np.ndarray]:
        """Stacked whitened data-minus-model residuals, and the model's h."""
        m = bg.Model.from_physical(Om, Ok, x[:, 0], x[:, 1])
        th = self.rec.theta_star(m)
        pred = np.column_stack([th, x[:, 0], x[:, 1]])
        r = (self.cmb.W @ (self.cmb.mean[None, :] - pred).T).T
        if self.bao is not None:
            rd = self.rec.r_drag(m)
            bp = self.bao.distances(m) / rd[:, None]
            r = np.hstack([r, (self.bao.W @ (self.bao.value[None, :] - bp).T).T])
        return r, m.h

    def __call__(self, Om, Ok, chunk: int = 4000) -> dict:
        Om = np.atleast_1d(np.asarray(Om, float)).ravel()
        Ok = np.atleast_1d(np.asarray(Ok, float)).ravel()
        out = {k: np.empty(Om.size) for k in ("m2lnl", "chi2", "ombh2", "omega_bc", "h", "OL",
                                              "sig_ombh2", "sig_omega_bc", "corr")}
        for sl in _chunks(Om.size, chunk):
            res = self._solve(Om[sl], Ok[sl])
            for k in out:
                out[k][sl] = res[k]
        return out

    def _solve(self, Om, Ok) -> dict:
        n = Om.size
        x = np.tile(self.cmb.mean[1:], (n, 1))
        for it in range(self.iterations + 1):
            r0, h = self.predictions(Om, Ok, x)
            J = np.empty(r0.shape + (2,))
            for j in range(2):
                dx = np.zeros(2)
                dx[j] = self.step[j]
                rj, _ = self.predictions(Om, Ok, x + dx)
                J[..., j] = -(rj - r0) / self.step[j]  # d(model, whitened)/dx
            F = np.einsum("nki,nkj->nij", J, J)
            g = np.einsum("nki,nk->ni", J, r0)
            if it == self.iterations:
                break
            delta = np.linalg.solve(F, g[..., None])[..., 0]
            # damp absurd steps (far from the data the problem is not quadratic)
            scale = np.maximum(1.0, np.max(np.abs(delta) / (20 * self.cmb.sigma[1:]), axis=1))
            x = x + delta / scale[:, None]
        chi2 = np.einsum("nk,nk->n", r0, r0)
        sign, logdet = np.linalg.slogdet(F)
        cov = np.linalg.inv(F)
        bad = ~np.isfinite(chi2) | (sign <= 0) | (h < cfg.PRIOR_H0[0] / 100) | \
            (h > cfg.PRIOR_H0[1] / 100)
        m2 = np.where(bad, np.inf, chi2 + logdet)
        OL = bg.Model.from_physical(Om, Ok, x[:, 0], x[:, 1]).OL
        return {"m2lnl": m2, "chi2": np.where(bad, np.inf, chi2), "ombh2": x[:, 0],
                "omega_bc": x[:, 1], "h": h, "OL": OL, "sig_ombh2": np.sqrt(cov[:, 0, 0]),
                "sig_omega_bc": np.sqrt(cov[:, 1, 1]),
                "corr": cov[:, 0, 1] / np.sqrt(cov[:, 0, 0] * cov[:, 1, 1])}


# ---------------------------------------------------------------------------
# From -2 ln L on a grid to probabilities
# ---------------------------------------------------------------------------

#: Probability inside 1, 2 and 3 sigma of a Gaussian.
SIGMA_MASS = (0.682689, 0.954500, 0.997300)


def grid_posterior(m2lnl: np.ndarray, mask=None) -> np.ndarray:
    """Normalised posterior (sums to 1) on a uniform grid, flat prior inside ``mask``."""
    finite = np.isfinite(m2lnl)
    p = np.zeros_like(m2lnl)
    p[finite] = np.exp(-0.5 * (m2lnl[finite] - m2lnl[finite].min()))
    if mask is not None:
        p = np.where(mask, p, 0.0)
    return p / p.sum()


def credible_levels(p: np.ndarray, masses=SIGMA_MASS) -> list[float]:
    """Density thresholds of the highest-posterior-density regions holding each mass."""
    flat = np.sort(p.ravel())[::-1]
    cum = np.cumsum(flat)
    cum /= cum[-1]
    return [float(flat[min(np.searchsorted(cum, m), flat.size - 1)]) for m in masses]


def weighted_summary(x: np.ndarray, w: np.ndarray) -> dict:
    """Mean, std and the central 68% interval of a weighted sample (or grid)."""
    x, w = np.asarray(x, float).ravel(), np.asarray(w, float).ravel()
    keep = (w > 0) & np.isfinite(x)
    x, w = x[keep], w[keep] / w[keep].sum()
    mean = float(np.sum(w * x))
    std = float(np.sqrt(np.sum(w * (x - mean) ** 2)))
    order = np.argsort(x)
    cdf = np.cumsum(w[order])
    lo, med, hi = np.interp([0.158655, 0.5, 0.841345], cdf, x[order])
    lo95, hi95 = np.interp([0.02275, 0.97725], cdf, x[order])
    return {"mean": mean, "std": std, "median": float(med), "lo68": float(lo),
            "hi68": float(hi), "lo95": float(lo95), "hi95": float(hi95)}
