"""From time-delay distances to H0.

The time-delay distance of a lens,

    D_dt = (1 + z_d) D_d D_s / D_ds,

is a ratio of angular-diameter distances (lens, source, lens-to-source), so
it scales as 1/H0 for any fixed shape of the expansion history. Each lens
team publishes a posterior for D_dt (and, where stellar kinematics were used,
the distance to the lens D_d). Turning those into H0 needs a cosmological
model for the distance-redshift relation: here flat LCDM with
Omega_m uniform in [0.05, 0.5], exactly as Wong et al. (2020).

The per-lens likelihoods follow H0LiCOW's own code (lensutils.py in their
release), reimplemented here:

* B1608+656: analytic skewed log-normal fits for D_dt (Suyu et al. 2010) and
  D_d (Jee et al. 2019), treated as independent.
* RX J1131, PG 1115: Gaussian kernel density (bandwidth 20 Mpc) on an 80 x 80
  histogram of the (D_d, D_dt) samples.
* HE 0435, WFI 2033: the same in 1-D (400 bins) for D_dt alone; WFI 2033's
  samples carry weights and are cut to 0 < D_dt < 8000 Mpc.
* SDSS J1206: scipy's Gaussian KDE (Scott's rule) on its (D_d, D_dt) samples,
  as lenstronomy's KDELikelihood does.

Instead of MCMC, the posterior is evaluated on a dense (H0, Omega_m) grid: in
flat LCDM every distance is (c/H0) times a function of Omega_m, so each grid
column costs one set of integrals.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.stats import gaussian_kde

import config as C

# ---------------------------------------------------------------------------
# Distances in flat LCDM
# ---------------------------------------------------------------------------


def comoving_unit(z, om: float, n: int = 2001) -> np.ndarray:
    """Line-of-sight comoving distance in units of c/H0 (flat LCDM, no radiation).

    No radiation term, matching astropy's FlatLambdaCDM default (Tcmb0 = 0)
    used by H0LiCOW's code.
    """
    z = np.atleast_1d(np.asarray(z, float))
    zz = np.linspace(0.0, z.max(), n)
    integrand = 1.0 / np.sqrt(om * (1 + zz) ** 3 + 1.0 - om)
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(zz))])
    return np.interp(z, zz, cum)


def distances(z_d: float, z_s: float, om, h0) -> tuple[np.ndarray, np.ndarray]:
    """(D_d, D_dt) in Mpc for arrays of Omega_m and H0 (broadcast together)."""
    om = np.atleast_1d(om)
    dd_u, ddt_u = [], []
    for o in om.ravel():
        x_d, x_s = comoving_unit([z_d, z_s], float(o))
        dd = x_d / (1 + z_d)
        ds = x_s / (1 + z_s)
        dds = (x_s - x_d) / (1 + z_s)
        dd_u.append(dd)
        ddt_u.append((1 + z_d) * dd * ds / dds)
    dd_u = np.array(dd_u).reshape(om.shape)
    ddt_u = np.array(ddt_u).reshape(om.shape)
    dh = C.C_KM_S / np.asarray(h0, float)  # Hubble distance, Mpc
    return dd_u * dh, ddt_u * dh


# ---------------------------------------------------------------------------
# Per-lens likelihoods
# ---------------------------------------------------------------------------


def _skewlognormal_logpdf(x, lam, mu, sigma):
    x = np.asarray(x, float)
    out = np.full(x.shape, -np.inf)
    ok = x > lam
    u = np.log(x[ok] - lam)
    out[ok] = -((u - mu) ** 2) / (2 * sigma**2) - np.log(np.sqrt(2 * np.pi) * sigma * (x[ok] - lam))
    return out


def _wq(x, w, qs):
    """Weighted quantiles."""
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    return np.interp(qs, cw / cw[-1], x[o])


def _gauss_matrix(x_eval, centers, h):
    return np.exp(-0.5 * ((x_eval[:, None] - centers[None, :]) / h) ** 2)


@dataclass
class LensLikelihood:
    """log L(D_d, D_dt) for one lens, tabulated on a fine mesh and interpolated."""

    key: str
    z_d: float
    z_s: float
    kind: str  # "1d" (D_dt only), "2d" (D_d, D_dt) or "skewlognormal"
    ddt: np.ndarray | None = None
    dd: np.ndarray | None = None
    weights: np.ndarray | None = None
    params: dict | None = None
    #: None: H0LiCOW's binning (400 bins in 1-D, 80 x 80 in 2-D, over the full sample
    #: range). A number: bins of about that many Mpc over the 0.05-99.95% range, used
    #: for the rescaled posteriors, whose wider spread would otherwise make the
    #: fixed bins coarser than the 20 Mpc kernel (and the density bumpy).
    bin_mpc: float | None = None

    def build(self):
        """Tabulate log-density on a mesh (KDE evaluated exactly at mesh nodes)."""
        if self.kind == "skewlognormal":
            return self
        h = 20.0  # Mpc, H0LiCOW's bandwidth
        w = np.ones_like(self.ddt) if self.weights is None else self.weights
        if self.kind == "1d":
            if self.bin_mpc:
                lo, hi = _wq(self.ddt, w, [0.0005, 0.9995])
                hist, edges = np.histogram(self.ddt, bins=int(np.clip((hi - lo) / self.bin_mpc, 50, 800)),
                                           range=(lo, hi), weights=w)
            else:
                hist, edges = np.histogram(self.ddt, bins=400, weights=w)
            cen = 0.5 * (edges[1:] + edges[:-1])
            keep = hist > 0
            x = np.linspace(edges[0] - 6 * h, edges[-1] + 6 * h, 3000)
            dens = _gauss_matrix(x, cen[keep], h) @ hist[keep]
            dens /= dens.sum() * (x[1] - x[0])
            self._x1 = x
            self._logd1 = np.log(np.clip(dens, 1e-300, None))
        elif self.kind == "2d":
            if self.bin_mpc:
                r_dd = _wq(self.dd, w, [0.0005, 0.9995])
                r_ddt = _wq(self.ddt, w, [0.0005, 0.9995])
                nb = [int(np.clip((r[1] - r[0]) / self.bin_mpc, 40, 400)) for r in (r_dd, r_ddt)]
                hist, e_dd, e_ddt = np.histogram2d(self.dd, self.ddt, bins=nb, range=[r_dd, r_ddt], weights=w)
            else:
                hist, e_dd, e_ddt = np.histogram2d(self.dd, self.ddt, bins=80, weights=w)
            c_dd, c_ddt = 0.5 * (e_dd[1:] + e_dd[:-1]), 0.5 * (e_ddt[1:] + e_ddt[:-1])
            x_dd = np.linspace(e_dd[0] - 6 * h, e_dd[-1] + 6 * h, 400)
            x_ddt = np.linspace(e_ddt[0] - 6 * h, e_ddt[-1] + 6 * h, 400)
            # Separable Gaussian kernel on a regular grid of bin centres:
            # density = Gx . hist . Gy^T (exactly sklearn's KDE on the binned samples).
            dens = _gauss_matrix(x_dd, c_dd, h) @ hist @ _gauss_matrix(x_ddt, c_ddt, h).T
            self._mesh2 = (x_dd, x_ddt, np.log(np.clip(dens, 1e-300, None)))
        elif self.kind == "kde_scipy":
            kde = gaussian_kde(np.vstack([self.dd, self.ddt]))
            pad_dd, pad_ddt = 4 * np.std(self.dd), 4 * np.std(self.ddt)
            x_dd = np.linspace(self.dd.min() - pad_dd, self.dd.max() + pad_dd, 160)
            x_ddt = np.linspace(self.ddt.min() - pad_ddt, self.ddt.max() + pad_ddt, 160)
            G = np.meshgrid(x_dd, x_ddt, indexing="ij")
            dens = kde(np.vstack([G[0].ravel(), G[1].ravel()])).reshape(G[0].shape)
            self._mesh2 = (x_dd, x_ddt, np.log(np.clip(dens, 1e-300, None)))
        return self

    def logl(self, dd, ddt) -> np.ndarray:
        dd, ddt = np.asarray(dd, float), np.asarray(ddt, float)
        if self.kind == "skewlognormal":
            p = self.params
            out = _skewlognormal_logpdf(ddt, p["lam"], p["mu"], p["sigma"])
            if "lam_d" in p:
                out = out + _skewlognormal_logpdf(dd, p["lam_d"], p["mu_d"], p["sigma_d"])
            return out
        if self.kind == "1d":
            out = np.interp(ddt, self._x1, self._logd1, left=-np.inf, right=-np.inf)
            return out
        x_dd, x_ddt, logd = self._mesh2
        f = RegularGridInterpolator((x_dd, x_ddt), logd, bounds_error=False, fill_value=-np.inf)
        return f(np.stack([dd, ddt], axis=-1))


def load_h0licow(key: str) -> LensLikelihood:
    """The published H0LiCOW likelihood of one lens (Wong et al. 2020)."""
    L = C.LENSES[key]
    d = C.DATA / "h0licow"
    zd, zs = L["z_d"], L["z_s"]
    if key == "B1608":
        row = np.loadtxt(d / "B1608_Dd_Ddt_params.dat")
        # columns: zd, zs, lambda_d, mu_d (stored as exp(mu)), sigma_d, lambda_dt, mu_dt, sigma_dt
        # (the header labels the last three "lambda_sft, mu_sft, sigma_sft" but
        #  H0LiCOW's notebook gives D_dt: mu=7.0531390, sigma=0.2282395, lam=4000
        #  and D_d: mu=6.79671, sigma=0.1836, lam=334.2; the file's numbers
        #  4000.0 7.053 0.2282 | 334.2 894.9 0.1836 match with mu_d = ln 894.9)
        params = dict(lam=4000.0, mu=7.0531390, sigma=0.2282395, lam_d=334.2, mu_d=6.79671, sigma_d=0.1836)
        assert abs(row[2] - 4000.0) < 1e-6 and abs(np.log(row[6]) - 6.79671) < 1e-4
        return LensLikelihood(key, zd, zs, "skewlognormal", params=params)
    if key == "RXJ1131":
        a = np.loadtxt(d / "RXJ1131_AO+HST_Dd_Ddt.dat", skiprows=1)
        return LensLikelihood(key, zd, zs, "2d", dd=a[:, 0], ddt=a[:, 1]).build()
    if key == "PG1115":
        a = np.loadtxt(d / "PG1115_AO+HST_Dd_Ddt.dat", skiprows=1)
        return LensLikelihood(key, zd, zs, "2d", dd=a[:, 0], ddt=a[:, 1]).build()
    if key == "HE0435":
        a = np.loadtxt(d / "HE0435_Ddt_AO+HST.dat")
        return LensLikelihood(key, zd, zs, "1d", ddt=a).build()
    if key == "WFI2033":
        a = np.loadtxt(d / "wfi2033_dt_bic.dat", delimiter=",", skiprows=1)
        keep = (a[:, 0] > 0) & (a[:, 0] < 8000)
        return LensLikelihood(key, zd, zs, "1d", ddt=a[keep, 0], weights=a[keep, 1]).build()
    if key == "J1206":
        a = np.loadtxt(d / "J1206_final.csv", delimiter=",", skiprows=1)
        return LensLikelihood(key, zd, zs, "kde_scipy", ddt=a[:, 0], dd=a[:, 1]).build()
    raise KeyError(key)


# ---------------------------------------------------------------------------
# Grid posterior
# ---------------------------------------------------------------------------

H0_GRID = np.arange(30.0, 150.0001, 0.05)  # the prior is uniform in [0, 150]; below 30 the likelihoods vanish
OM_GRID = np.linspace(0.05, 0.5, 91)


@dataclass
class H0Posterior:
    h0: np.ndarray
    om: np.ndarray
    logpost: np.ndarray  # (n_om, n_h0)

    @property
    def p_h0(self) -> np.ndarray:
        lp = self.logpost - np.nanmax(self.logpost)
        p = np.exp(lp).sum(axis=0)
        return p / (p.sum() * (self.h0[1] - self.h0[0]))

    def summary(self) -> tuple[float, float, float]:
        """(median, minus, plus) from the 16th/50th/84th percentiles."""
        p = self.p_h0
        cdf = np.cumsum(p)
        cdf /= cdf[-1]
        q = np.interp([0.15865525, 0.5, 0.84134475], cdf, self.h0)
        return float(q[1]), float(q[1] - q[0]), float(q[2] - q[1])


def grid_logl(lik: LensLikelihood, h0=H0_GRID, om=OM_GRID) -> np.ndarray:
    """log-likelihood of one lens on the (Omega_m, H0) grid."""
    dd, ddt = distances(lik.z_d, lik.z_s, om[:, None], h0[None, :])
    return lik.logl(dd, ddt)


def posterior(logls: list[np.ndarray], h0=H0_GRID, om=OM_GRID) -> H0Posterior:
    """Flat priors (H0 in [0,150], Omega_m in [0.05,0.5]) times the lenses' likelihoods."""
    lp = np.sum(logls, axis=0)
    return H0Posterior(h0, om, lp)


# ---------------------------------------------------------------------------
# Swapping in our own delays
# ---------------------------------------------------------------------------


def delay_scale(our: dict, our_cov: np.ndarray, pub: dict, pub_err: dict) -> dict:
    """How much our delays rescale a published lens model's distances.

    For a fixed lens model every delay is D_dt * (Fermat potential
    difference) / c, so changing the measured delays by a common factor r
    changes the inferred D_dt by exactly r. With several delays, r is the
    generalized least-squares scale between our delay vector d (covariance
    C) and the delays the model was fitted to, p (which stand in for the
    model's Fermat-potential *ratios*): r = p'C^-1 d / p'C^-1 p, with
    fractional error sigma_r / r.

    ``our``/``pub``: pair -> delay (same pairs, same order); ``pub_err``:
    pair -> symmetric 1-sigma error of the published delays.
    Returns r, its fractional error, and the published fractional error
    computed the same way (for comparison with Wong et al. 2020 Table 3).
    """
    pairs = list(pub)
    p = np.array([pub[k] for k in pairs])
    d = np.array([our[k] for k in pairs])
    Ci = np.linalg.inv(our_cov)
    r = float(p @ Ci @ d / (p @ Ci @ p))
    sig_r = float(1.0 / np.sqrt(p @ Ci @ p))
    Cp = np.diag([pub_err[k] ** 2 for k in pairs])
    frac_pub = float(1.0 / np.sqrt(p @ np.linalg.inv(Cp) @ p))
    return dict(r=r, frac_ours=sig_r / r, frac_pub=frac_pub, pairs=pairs)


def swap_delays(lik: LensLikelihood, r: float, frac_ours: float, frac_pub: float) -> LensLikelihood:
    """A copy of a sample-based likelihood with our delays in place of the published ones.

    In log space the published D_dt width is the lens model's width plus
    the delays' width in quadrature. We shift every sample by ln r and
    rescale its offset from the mean so that the delay part of the variance
    becomes ours: s^2 = (V - frac_pub^2 + frac_ours^2) / V. D_d moves with
    D_dt sample by sample, since the kinematics constrain D_d/D_dt, not the
    delays.
    """
    if lik.kind == "skewlognormal":
        raise ValueError("B1608 is not re-measured (radio light curves not public)")
    x = np.log(lik.ddt)
    w = np.ones_like(x) if lik.weights is None else lik.weights
    mu = np.average(x, weights=w)
    V = np.average((x - mu) ** 2, weights=w)
    s = np.sqrt(max(V - frac_pub**2 + frac_ours**2, 1e-6) / V)
    x_new = mu + np.log(r) + (x - mu) * s
    ddt_new = np.exp(x_new)
    dd_new = None if lik.dd is None else lik.dd * ddt_new / lik.ddt
    return LensLikelihood(lik.key, lik.z_d, lik.z_s, lik.kind, ddt=ddt_new, dd=dd_new,
                          weights=lik.weights, bin_mpc=10.0).build()
