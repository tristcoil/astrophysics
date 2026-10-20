"""Forward model: how many clusters should eRASS1 have found, for a given universe?

For each cosmology (Omega_m, sigma_8) and count rate-mass relation, we predict
the number of clusters in bins of measured count rate and photometric redshift:

    N(CR_hat bin, z_hat bin) = Area * integral dz  dV/dz/dOmega
        * integral dlnM  dn/dlnM(M, z) * P(lambda_hat > 3 | M, z)
        * integral dlnCR  P(CR | M, z)                  <- scaling relation, scatter
        * < P(detected | CR, z, sky) * P(CR_hat in bin | CR, sky)
            * [z_hat below the sky's optical depth limit] >_sky
        * P(z_hat in bin | z)

The sky average <...> runs over the footprint's mix of exposure times,
Galactic absorption and optical depth (footprint.survey_classes) and does
not depend on cosmology, so it is computed once ("kernel") and the
likelihood only redoes the cheap part. Measured count rates scatter around the
true ones by sigma_ln = sqrt(a / photons), calibrated on the catalogue's own
error bars (a = 2.5, i.e. 1.6 x Poisson).

The data side: each cluster counts as (1 - PCONT), its probability of being a
real cluster rather than an AGN or a noise fluctuation (Bulbul et al. 2024).
The likelihood is Poisson in the 20 x 14 bins.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import ndtr

import config as C
import cosmology as CO
import footprint as F
import selection as S

@dataclass
class Scaling:
    A_X: float = C.SCALING_PUBLISHED["A_X"][0]
    B_X: float = C.SCALING_PUBLISHED["B_X"][0]
    F_X: float = C.SCALING_PUBLISHED["F_X"][0]
    G_X: float = C.SCALING_PUBLISHED["G_X"][0]
    sigma_X: float = C.SCALING_PUBLISHED["sigma_X"][0]
    mass_scale: float = 1.0  # M_true / M_assumed: a global mass bias (tests)

    def mean_lncr(self, cosmo: CO.Cosmology, lnM, z):
        """<ln CR | M, z> on a (z, M) grid (Ghirardini et al. Eqs. 22-24)."""
        z = np.atleast_1d(z)
        lz = np.log((1 + z) / (1 + C.Z_PIVOT))
        dl = cosmo.luminosity_distance(np.append(z, C.Z_PIVOT))
        ez = cosmo.E(np.append(z, C.Z_PIVOT))
        ev = C.D_X * np.log(dl[:-1] / dl[-1]) + C.E_X * np.log(ez[:-1] / ez[-1]) + self.G_X * lz
        slope = self.B_X + self.F_X * lz
        lnm = np.asarray(lnM)[None, :] - np.log(C.M_PIVOT * self.mass_scale)
        return np.log(C.CR_PIVOT * self.A_X) + slope[:, None] * lnm + ev[:, None]


def bins():
    lo, hi, n = C.CR_HAT_EDGES_LOG10
    cr_edges = np.log(10.0 ** np.linspace(lo, hi, n + 1))
    z0, z1, nz = C.Z_HAT_EDGES
    return cr_edges, np.linspace(z0, z1, nz + 1)


def grids():
    lnm = np.log(np.geomspace(*C.LNM_RANGE))
    z = np.linspace(*C.Z_RANGE)
    lncr = np.log(np.geomspace(*C.LNCR_RANGE))
    return lnm, z, lncr


def richness_completeness(lnm, z):
    """P(lambda_hat > 3 | M, z) from the published richness-mass relation.

    ln lambda ~ N(ln A + b(z) ln(M/M_p) + C ln((1+z)/(1+z_p)), sigma_lambda),
    and the measured richness scatters by Poisson noise, sqrt(lambda).
    """
    r = C.RICHNESS_PUBLISHED
    lz = np.log((1 + z) / (1 + C.Z_PIVOT))[:, None]
    mu = r["lnA_l"] + (r["B_l"] + r["D_l"] * lz) * (lnm[None, :] - np.log(C.M_PIVOT)) + r["C_l"] * lz
    x = np.linspace(-4, 4, 41)
    w = np.exp(-0.5 * x**2)
    w /= w.sum()
    lam = np.exp(mu[..., None] + r["sigma_l"] * x)
    return (ndtr((lam - C.LAMBDA_MIN) / np.sqrt(lam)) * w).sum(-1)


def measurement_sigma(photons, a=2.52):
    return np.sqrt(a / np.maximum(photons, 0.05))


def build_kernel(selection: S.SelectionModel | None = None, classes: dict | None = None,
                 tag: str = "default", meas_a: float = 2.52, use_nh: bool = True,
                 use_zvlim: bool = True, refresh: bool = False) -> dict:
    """Sky-averaged selection x measurement kernel K[cr_hat bin, z, CR] (cached)."""
    path = C.CACHE_DIR / f"selection_kernel_{tag}.npz"
    if path.exists() and not refresh:
        return dict(np.load(path))
    selection = S.SelectionModel() if selection is None else selection
    classes = F.survey_classes(F.build()) if classes is None else classes
    lnm, z, lncr = grids()
    cr_edges, _ = bins()
    cr = np.exp(lncr)
    inner = cr_edges[1:-1]
    K = np.zeros((len(cr_edges) - 1, z.size, lncr.size))
    t, nh, zv, w = classes["t"], classes["nh"], classes["zvlim"], classes["weight"]
    for c in range(w.size):
        p_det = selection(cr[None, :], z[:, None], t[c], nh[c] if use_nh else None)  # (z, CR)
        if use_zvlim:
            p_det = p_det * (C.PHOTOZ_BIAS * z < zv[c])[:, None]
        photons = cr * S.absorbed_fraction(nh[c]) * t[c]
        sig = measurement_sigma(photons, meas_a)
        cdf = ndtr((inner[:, None] - lncr[None, :]) / sig[None, :])  # P(ln CR_hat < edge)
        cdf = np.vstack([np.zeros(lncr.size), cdf, np.ones(lncr.size)])
        pbin = np.diff(cdf, axis=0)  # open-ended first and last bins
        K += w[c] * pbin[:, None, :] * p_det[None, :, :]
    out = {"K": K, "z": z, "lncr": lncr, "area_deg2_mask": classes["area_deg2"]}
    np.savez_compressed(path, **out)
    return out


def photoz_matrix(z, z_edges):
    """P(z_hat in bin j | z): Gaussian, mean b_z z, width sigma_z (1 + z)."""
    mu = C.PHOTOZ_BIAS * z
    sd = C.PHOTOZ_SIGMA * (1 + z)
    cdf = ndtr((z_edges[:, None] - mu[None, :]) / sd[None, :])
    return np.diff(cdf, axis=0)  # (n_zhat, n_z)


class CountsModel:
    def __init__(self, kernel: dict | None = None, area_deg2: float = C.AREA_COSMO_DEG2):
        self.kernel = build_kernel() if kernel is None else kernel
        self.lnm, self.z, self.lncr = grids()
        self.cr_edges, self.z_edges = bins()
        self.area_sr = area_deg2 * (np.pi / 180) ** 2
        self.dz = np.gradient(self.z)
        self.dlnm = np.gradient(self.lnm)
        self.dlncr = np.gradient(self.lncr)
        self.Pz = photoz_matrix(self.z, self.z_edges)
        self.f_lambda = richness_completeness(self.lnm, self.z)
        # integrand weights restricted to the observed redshift window
        self.Kz = np.einsum("izk,jz->ijzk", self.kernel["K"], self.Pz)  # (cr_hat, z_hat, z, CR)

    def dn_dlncr(self, cosmo: CO.Cosmology, sr: Scaling):
        """Clusters per steradian per unit z per unit ln CR (true CR), (z, CR)."""
        M = np.exp(self.lnm)
        hmf = cosmo.dndlnM(M, self.z) * self.f_lambda  # (z, M)
        mu = sr.mean_lncr(cosmo, self.lnm, self.z)  # (z, M)
        g = np.exp(-0.5 * ((self.lncr[None, None, :] - mu[:, :, None]) / sr.sigma_X) ** 2)
        g /= np.sqrt(2 * np.pi) * sr.sigma_X
        dn = np.einsum("zm,zmk->zk", hmf * self.dlnm[None, :], g)
        return dn * cosmo.dV_dz_dOmega(self.z)[:, None]

    def expected(self, cosmo: CO.Cosmology, sr: Scaling | None = None):
        sr = Scaling() if sr is None else sr
        dn = self.dn_dlncr(cosmo, sr) * self.dz[:, None] * self.dlncr[None, :]
        return self.area_sr * np.einsum("ijzk,zk->ij", self.Kz, dn)

    def detected_true(self, cosmo, sr=None):
        """Expected detections per true (z, CR) cell, before binning (for plots)."""
        sr = Scaling() if sr is None else sr
        dn = self.dn_dlncr(cosmo, sr)
        return self.area_sr * dn * self.kernel["K"].sum(0) * (self.Pz.sum(0)[:, None])


def data_counts(df, cr_edges=None, z_edges=None, weighted: bool = True):
    """The observed (weighted) counts in the same bins."""
    if cr_edges is None:
        cr_edges, z_edges = bins()
    inner = cr_edges.copy()
    inner[0], inner[-1] = -np.inf, np.inf
    w = df.weight.values if weighted else None
    n, _, _ = np.histogram2d(np.log(df.cr.values), df.z.values, bins=[inner, z_edges], weights=w)
    return n


def loglike(model_counts, data):
    m = np.maximum(model_counts, 1e-30)
    return float(np.sum(data * np.log(m) - m))
