"""Scoring universes against the data: BAO, supernovae, and the CMB prior.

Each likelihood returns chi^2 = -2 ln L (up to a constant) for a batch of
universes. They are written to match, term by term, the likelihoods DESI used
for its DR2 cosmology (the Cobaya implementations), so that our chi^2 values can
be compared with the per-likelihood chi^2 DESI published at its best fits.

* **BAO**: a Gaussian in the 13 DESI DR2 distance ratios (D_V/r_d, D_M/r_d,
  D_H/r_d) with their covariance.
* **Supernovae**: a Gaussian in distance moduli with the sample's full
  covariance. Every sample has an unknown overall offset (the supernova
  absolute magnitude combined with H0), removed exactly by profiling, which for
  a Gaussian with a flat prior equals marginalising it.
* **CMB**: a Gaussian prior on (theta_*, omega_b, omega_b + omega_c), the
  "compressed" early-universe information from Planck.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import cholesky, solve_triangular

import background as bg
from config import CMB_PRIOR_DESI, PLANCK_EARLY_BURN


def whitening_from_cov(cov: np.ndarray) -> np.ndarray:
    """A with A^T A = C^-1 (so chi^2 = |A r|^2): A = L^-1 for C = L L^T."""
    lower = cholesky(cov, lower=True)
    return solve_triangular(lower, np.eye(cov.shape[0]), lower=True)


def whitening_from_invcov(invcov: np.ndarray) -> np.ndarray:
    """A with A^T A = W for a given inverse covariance W = U U^T: A = U^T."""
    return cholesky(0.5 * (invcov + invcov.T), lower=True).T


# ---------------------------------------------------------------------------
# BAO
# ---------------------------------------------------------------------------


@dataclass
class BAOLikelihood:
    """Gaussian likelihood for a BAO data vector."""

    z: np.ndarray
    value: np.ndarray
    kind: np.ndarray
    cov: np.ndarray
    name: str = "DESI DR2"
    whiten: np.ndarray = field(init=False, repr=False)

    def __post_init__(self):
        self.whiten = whitening_from_cov(self.cov)
        self.unique_z, self.z_index = np.unique(self.z, return_inverse=True)

    @property
    def n(self) -> int:
        return self.value.size

    def model_unscaled(self, dm_unique: np.ndarray, hz_unique: np.ndarray) -> np.ndarray:
        """The data vector's distances in Mpc (r_d = 1), (n_models, n)."""
        dm = dm_unique[:, self.z_index]
        hz = hz_unique[:, self.z_index]
        return bg.bao_predictions(self.kind, self.z, dm, hz, np.ones(dm.shape[0]))

    def chi2(self, model: np.ndarray) -> np.ndarray:
        r = self.whiten @ (self.value[None, :] - model).T
        return np.einsum("ij,ij->j", r, r)

    def chi2_profile_scale(self, unscaled: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """chi^2 minimised over 1/r_d (the model is linear in it). Returns chi^2, best r_d."""
        u = self.whiten @ unscaled.T  # (n, k)
        d = self.whiten @ self.value
        uu = np.einsum("ij,ij->j", u, u)
        ud = d @ u
        s = ud / uu
        return d @ d - ud**2 / uu, 1.0 / s


# ---------------------------------------------------------------------------
# Supernovae
# ---------------------------------------------------------------------------


@dataclass
class SNLikelihood:
    """Distance moduli with full covariance and an analytically removed offset."""

    name: str
    z_hd: np.ndarray
    z_hel: np.ndarray
    mu: np.ndarray
    whiten: np.ndarray = field(repr=False)

    @classmethod
    def from_sample(cls, sample: dict) -> "SNLikelihood":
        A = (whitening_from_cov(sample["cov"]) if "cov" in sample
             else whitening_from_invcov(sample["invcov"]))
        return cls(name=sample["name"], z_hd=sample["z_hd"], z_hel=sample["z_hel"],
                   mu=sample["mu"], whiten=A)

    def __post_init__(self):
        self.u = self.whiten @ np.ones(self.mu.size)
        self.uu = float(self.u @ self.u)
        self.wmu = self.whiten @ self.mu

    @property
    def n(self) -> int:
        return self.mu.size

    def model_mu(self, dm: np.ndarray) -> np.ndarray:
        """5 log10(d_L / Mpc) + 25 with d_L = (1 + z_hel) D_M(z_hd), as in Cobaya."""
        return 5.0 * np.log10((1.0 + self.z_hel)[None, :] * dm) + 25.0

    def chi2(self, dm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Profiled chi^2 and the best offset, for D_M (Mpc) at the SN redshifts."""
        m = self.model_mu(dm)
        r = self.wmu[:, None] - self.whiten @ m.T  # (n, k)
        ur = self.u @ r
        return np.einsum("ij,ij->j", r, r) - ur**2 / self.uu, ur / self.uu

    def residuals(self, dm_one: np.ndarray) -> np.ndarray:
        """Hubble residuals (data - model - best offset) for one model, for plots."""
        m = self.model_mu(dm_one[None, :])[0]
        _, off = self.chi2(dm_one[None, :])
        return self.mu - m - off[0]


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

    def __post_init__(self):
        self.mean = np.asarray(self.mean, float)
        self.cov = np.asarray(self.cov, float)
        self.icov = np.linalg.inv(self.cov)

    def chi2(self, theta_star, ombh2, omega_bc) -> np.ndarray:
        r = np.vstack([theta_star, ombh2, omega_bc]).T - self.mean[None, :]
        return np.einsum("ki,ij,kj->k", r, self.icov, r)

    @property
    def sigma(self) -> np.ndarray:
        return np.sqrt(np.diag(self.cov))


def cmb_prior_desi() -> CMBPrior:
    """Exactly what DESI's chains used (its released likelihood configuration)."""
    return CMBPrior("PR4 compressed (DESI as run)", CMB_PRIOR_DESI["mean"], CMB_PRIOR_DESI["cov"],
                    note="third mean is <omega_m h^2> (incl. the 0.06 eV neutrino), "
                         "compared with omega_b + omega_c")


def cmb_prior_from_chains(chains: dict, *, use_omegam_mean: bool = False) -> CMBPrior:
    """Our recomputation from Lemos & Lewis's chains (weighted mean and covariance).

    ``use_omegam_mean=False`` (default) gives the self-consistent version: the
    mean of omega_b + omega_c itself. ``True`` reproduces DESI's choice.
    """
    w = chains["weight"]
    th = chains["thetastar"] / 100.0
    ob, oc = chains["ombh2"], chains["omch2"]
    X = np.vstack([th, ob, ob + oc])
    mean = np.average(X, axis=1, weights=w)
    if use_omegam_mean:
        mean[2] = np.average(chains["omegamh2"], weights=w)
    cov = np.cov(X, aweights=w)
    name = "PR4 compressed (" + ("DESI's mean" if use_omegam_mean else "consistent") + ")"
    return CMBPrior(name, mean, cov, note=f"{len(w)} chain rows after {PLANCK_EARLY_BURN:.0%} burn-in")


def cmb_prior_from_d3(d3_fit: dict, early) -> CMBPrior:
    """Our own CMB prior, from project D3's LCDM fit to the Planck 2018 TT spectrum.

    D3 fitted (omega_b, omega_c, 100 theta_MC, ln A, n_s) with tau fixed and
    published the covariance as errors + correlation matrix. We keep the
    (omega_b, omega_c, theta_MC) block -- marginalising over the amplitude and
    tilt -- and convert theta_MC (CosmoMC's approximation) to theta_* with the
    Jacobian of CAMB's exact theta_* at D3's best fit.
    """
    import camb

    p = d3_fit["parameters"]
    names = d3_fit["correlation"]["names"]
    err = np.array([p[k]["error"] for k in names])
    corr = np.array(d3_fit["correlation"]["matrix"])
    cov5 = corr * np.outer(err, err)
    idx = [names.index(k) for k in ("ombh2", "omch2", "theta100")]
    cov3 = cov5[np.ix_(idx, idx)]
    x0 = np.array([p["ombh2"]["value"], p["omch2"]["value"], p["theta100"]["value"]])

    def theta_star(x):
        pars = camb.set_params(ombh2=x[0], omch2=x[1], cosmomc_theta=x[2] / 100.0,
                               mnu=0.06, nnu=3.046, tau=d3_fit.get("tau_fixed", 0.0522),
                               As=2.1e-9, ns=0.9649)  # D3's settings (CAMB defaults)
        return camb.get_background(pars).get_derived_params()["thetastar"] / 100.0

    th0 = theta_star(x0)
    steps = np.sqrt(np.diag(cov3)) * 0.5
    grad = np.array([(theta_star(x0 + steps[i] * np.eye(3)[i])
                      - theta_star(x0 - steps[i] * np.eye(3)[i])) / (2 * steps[i])
                     for i in range(3)])
    # (theta_*, omega_b, omega_b + omega_c) = J (omega_b, omega_c, theta_MC)
    J = np.array([grad, [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]])
    mean = np.array([th0, x0[0], x0[0] + x0[1]])
    return CMBPrior("D3: our Planck 2018 TT fit", mean, J @ cov3 @ J.T,
                    note="TT only, tau fixed, lensing as in LCDM")
