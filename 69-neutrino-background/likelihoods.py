"""Scoring universes against data.

Every likelihood returns a vector of whitened residuals r for one universe
(χ² = |r|²), so a data combination is just the concatenation and fitting is
least squares (see fitting.py).

* :class:`BinnedTT` — Planck's 2018 binned temperature spectrum with
  independent errors per bin: exactly project D3's likelihood (83 bands,
  ℓ = 30–2508), the spec's "D3's TT spectrum".
* :class:`PlikLite` — Planck's own foreground-marginalised high-ℓ likelihood
  ("plik_lite", Planck 2018 V): TT, TE, EE bandpowers with their full
  covariance and the overall calibration y_P = A_planck.
* :class:`GaussianPrior` — one parameter, e.g. the optical depth τ from
  Planck's low-ℓ polarisation, standing in for the low-ℓ likelihoods.
* :class:`PlanckLensing2018`, :class:`ACTLensing` — CMB lensing
  reconstruction band powers with the linear corrections for the
  CMB-dependence of the reconstruction.
* :class:`BAO` — DESI DR2's 13 distance ratios with their covariance;
  :class:`BAO2018` — the 6dF + MGS + BOSS DR12 set Planck 2018 used.
* :class:`Supernovae` — a Type Ia sample with full covariance; the unknown
  absolute magnitude is profiled out analytically.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import cholesky, solve_triangular

import config as cfg


def whitening(cov: np.ndarray) -> np.ndarray:
    """A with Aᵀ A = C⁻¹ (so χ² = |A r|²): A = L⁻¹ for C = L Lᵀ."""
    lower = cholesky(cov, lower=True)
    return solve_triangular(lower, np.eye(cov.shape[0]), lower=True)


def whitening_from_inverse(icov: np.ndarray) -> np.ndarray:
    """A with Aᵀ A = W for a given inverse covariance W = U Uᵀ: A = Uᵀ."""
    return cholesky(0.5 * (icov + icov.T), lower=True).T


def cl_from_dl(ell, dl):
    ell = np.asarray(ell, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(ell > 0, 2 * np.pi * np.asarray(dl, float) / (ell * (ell + 1)), 0.0)


# ---------------------------------------------------------------------------
# Planck temperature, D3's way
# ---------------------------------------------------------------------------


class BinnedTT:
    """Planck 2018 binned TT (COM_PowerSpect_CMB-TT-binned_R3.01), diagonal errors.

    Planck bins C_ℓ (not D_ℓ) with weights ℓ(ℓ+1) in bins of 30 from ℓ = 30 and
    places each band at the weighted-mean multipole; D3 reverse-engineered this
    and reproduced Planck's binned best-fit column to 3 parts in 10⁹. The model
    is binned identically.
    """

    name = "Planck TT (D3)"
    nuisance = ()

    def __init__(self, table: np.ndarray, lmin: int = 30, lmax: int = 2508, width: int = 30):
        self.l_eff, self.dl, self.err = table[:, 0], table[:, 1], table[:, 2]
        self.bestfit_binned = table[:, 4] if table.shape[1] > 4 else None
        lo = np.arange(lmin, lmax + 1, width)
        hi = np.minimum(lo + width - 1, lmax)
        self.lo, self.hi = lo, hi
        if len(lo) != len(self.dl):
            raise ValueError(f"{len(self.dl)} bands, expected {len(lo)}")
        # Binning matrix acting on D_ℓ for ℓ = 0..lmax: D_b = l_eff(l_eff+1)/2π × Σ w C_ℓ / Σ w.
        self.W = np.zeros((len(lo), lmax + 1))
        for i, (a, z) in enumerate(zip(lo, hi)):
            ls = np.arange(a, z + 1, dtype=float)
            w = ls * (ls + 1.0)
            leff = (w * ls).sum() / w.sum()
            # C_ℓ = 2π D_ℓ / (ℓ(ℓ+1)) and w = ℓ(ℓ+1): Σ w C_ℓ = 2π Σ D_ℓ
            self.W[i, a:z + 1] = 2 * np.pi / w.sum() * leff * (leff + 1) / (2 * np.pi)
        self.lmax = lmax

    def model(self, theory: dict) -> np.ndarray:
        return self.W @ theory["tt"][: self.lmax + 1]

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        return (self.model(theory) - self.dl) / self.err


# ---------------------------------------------------------------------------
# Planck plik_lite (TT, TE, EE with covariance)
# ---------------------------------------------------------------------------


class PlikLite:
    """Planck 2018 plik_lite: CMB-only bandpowers, foregrounds marginalised.

    The data vector is binned C_ℓ (TT, then TE, then EE), each band a weighted
    sum of C_ℓ over [blmin, blmax] with weights ``bweight``. The model C_ℓ is
    divided by A_planck² (Planck's calibration), whose Gaussian prior
    1.0000 ± 0.0025 is added as one more residual. ``use`` picks spectra and
    ``lmax_tt``/``lmax_pol`` cut multipoles (P-ACT keeps TT < 1000, TE/EE < 600).
    """

    nuisance = ("A_planck",)

    def __init__(self, data: dict, *, use=("TT", "TE", "EE"), lmin: int = 30,
                 lmax_tt: int = 2508, lmax_pol: int = 1996, name: str | None = None):
        self.name = name or "Planck plik_lite " + "".join(use)
        nbins = data["nbins"]  # {"TT": 215, "TE": 199, "EE": 199}
        blmin, blmax, bw = data["blmin"], data["blmax"], data["bweight"]
        lmin_data = data["plmin"]
        offset = {"TT": 0, "TE": nbins["TT"], "EE": nbins["TT"] + nbins["TE"]}
        keep, rows = [], []
        self.spec_index = []
        for s in ("TT", "TE", "EE"):
            if s not in use:
                continue
            lmax_s = lmax_tt if s == "TT" else lmax_pol
            for b in range(nbins[s]):
                lo = blmin[b] + lmin_data
                hi = blmax[b] + lmin_data
                if lo < lmin or hi > lmax_s:
                    continue
                keep.append(offset[s] + b)
                row = np.zeros(hi + 1)
                row[lo:hi + 1] = bw[blmin[b]:blmax[b] + 1]
                rows.append((s, row))
                self.spec_index.append(s)
        self.keep = np.array(keep)
        self.lmax = max(len(r) for _, r in rows) - 1
        self.Wt = {s: np.zeros((len(rows), self.lmax + 1)) for s in ("TT", "TE", "EE")}
        for i, (s, row) in enumerate(rows):
            self.Wt[s][i, : len(row)] = row
        self.data = data["cl"][self.keep]
        cov = data["cov"][np.ix_(self.keep, self.keep)]
        self.A = whitening(cov)
        self.ell = np.arange(self.lmax + 1)
        self.n = len(self.keep)

    def model(self, theory: dict) -> np.ndarray:
        """Binned C_ℓ (µK²) before calibration."""
        out = np.zeros(self.n)
        for s, key in (("TT", "tt"), ("TE", "te"), ("EE", "ee")):
            cl = cl_from_dl(self.ell, theory[key][: self.lmax + 1])
            out += self.Wt[s] @ cl
        return out

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        a = p.get("A_planck", 1.0)
        r = self.A @ (self.model(theory) / a**2 - self.data)
        return np.concatenate([r, [(a - 1.0) / cfg.A_PLANCK_PRIOR[1]]])


# ---------------------------------------------------------------------------
# CMB lensing
# ---------------------------------------------------------------------------


class PlanckLensing2018:
    """Planck 2018 lensing reconstruction ("consext8", L = 8-400), as Planck 2018 VI.

    Each band power is a window over [L(L+1)]² C_L^φφ/2π plus a linear
    correction for how the reconstruction depends on the CMB spectra (the
    calibrated D_ℓ^TT, EE, TE and the lensing spectrum itself), minus that
    correction on Planck's fiducial model.
    """

    name = "Planck 2018 lensing"
    nuisance = ("A_planck",)

    def __init__(self, data: dict):
        self.data = data["data"]
        self.A = whitening(data["cov"])
        self.W, self.M, self.fid = data["W"], data["M"], data["fid_correction"]
        self.lmax = data["lmax_theory"]

    def model(self, theory: dict, a_planck: float = 1.0) -> np.ndarray:
        n = self.lmax + 1
        pp = theory["pp"][:n]
        cal = 1.0 / a_planck**2
        out = self.W @ pp + self.M["pp"] @ pp - self.fid
        for x in ("tt", "ee", "te"):
            out = out + cal * (self.M[x] @ theory[x][:n])
        return out

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        return self.A @ (self.model(theory, p.get("A_planck", 1.0)) - self.data)


class ACTLensing:
    """ACT DR6 + Planck PR4 lensing (act_dr6_lenslike v1.2, ``actplanck_baseline``).

    Band powers of C_L^κκ = [L(L+1)]² C_L^φφ/4: ten ACT bins (L = 40-763) and
    nine Planck PR4 bins (L = 8-400). The theory is corrected linearly for the
    CMB-dependence of the reconstruction's normalisation and N1 bias; those
    corrections were projected onto the bins once (fetch_data.reduce_act_lensing),
    so each band is K·C^κκ + Σ_X M_X·C_X + const with lensed C_X in µK².
    The inverse covariance carries the Hartlap factor (400 − 19 − 2)/(400 − 1).
    """

    name = "ACT DR6 + Planck PR4 lensing"
    nuisance = ()

    def __init__(self, data: dict, n_sims: int = 400):
        self.d = data
        self.data = data["data"]
        nb = self.data.size
        self.hartlap = (n_sims - nb - 2.0) / (n_sims - 1.0)
        self.A = np.sqrt(self.hartlap) * whitening(data["cov"])
        self.n = data["act_K"].shape[1]
        ell = np.arange(self.n, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            self.dl_to_cl = np.where(ell > 0, 2 * np.pi / (ell * (ell + 1)), 0.0)

    def model(self, theory: dict) -> np.ndarray:
        n = self.n
        ckk = theory["pp"][:n] * 2 * np.pi / 4.0
        cls = {x: theory[x][:n] * self.dl_to_cl for x in ("tt", "ee", "bb", "te")}
        out = []
        for t in ("act", "planck"):
            b = self.d[f"{t}_K"] @ ckk + self.d[f"{t}_const"]
            for x in ("tt", "ee", "bb", "te"):
                b = b + self.d[f"{t}_M_{x}"] @ cls[x]
            out.append(b)
        return np.concatenate(out)

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        return self.A @ (self.model(theory) - self.data)


# ---------------------------------------------------------------------------
# Planck's low-ℓ likelihoods (ℓ = 2-29)
# ---------------------------------------------------------------------------


class PlanckLowlTT:
    """Planck 2018 low-ℓ temperature (Commander), ℓ = 2-29.

    A Gaussianised Blackwell-Rao likelihood (Planck 2018 V): each D_ℓ maps to a
    variable x_ℓ(D_ℓ) through a tabulated monotonic transform, and the x_ℓ are
    Gaussian with mean μ and covariance C:

        −2 ln L = (x − μ)ᵀ C⁻¹ (x − μ) − 2 Σ ln(dx_ℓ/dD_ℓ) + const,

    normalised (as Planck's code does) to zero at the reference spectrum. As
    residuals: the whitened x − μ, plus one term √(K − 2 Σ ln dx/dD) for the
    Jacobian (K a constant that keeps it positive and is removed again in
    ``chi2_constant``). Ported from Cobaya's translation of Planck's clik code.
    """

    name = "Planck low-ℓ TT"
    nuisance = ("A_planck",)
    lmin, lmax = 2, 29

    def __init__(self, data: dict):
        sl = slice(self.lmin - 2, self.lmax - 1)
        self.A = whitening(data["cov"][sl, sl])
        self.mu = data["mu"][sl]
        self.cl_grid = data["cl_grid"][:, sl]
        self.x_grid = data["x_grid"][:, sl]
        n = self.lmax - self.lmin + 1
        self.bounds = np.zeros((n, 2))
        for i in range(n):  # where the transform saturates at x = ±5 (Planck's prior edges)
            j = 0
            while abs(self.x_grid[j, i] + 5) < 1e-4:
                j += 1
            self.bounds[i, 0] = self.cl_grid[j + 2, i]
            j = self.x_grid.shape[0] - 1
            while abs(self.x_grid[j, i] - 5) < 1e-4:
                j -= 1
            self.bounds[i, 1] = self.cl_grid[j - 2, i]
        self._splines = None
        ref = np.zeros(self.lmax + 1)
        ref[self.lmin:] = data["mu_sigma"][sl]
        x, lj = self._transform(ref[self.lmin:])
        r = self.A @ (x - self.mu)
        self.offset = float(lj - 0.5 * r @ r)  # ln L at the reference spectrum
        self.K = 400.0
        self.chi2_constant = -self.K

    def _build(self):
        from scipy.interpolate import InterpolatedUnivariateSpline

        sp = [InterpolatedUnivariateSpline(self.cl_grid[:, i], self.x_grid[:, i])
              for i in range(self.cl_grid.shape[1])]
        self._splines = [(s, s.derivative()) for s in sp]

    def __getstate__(self):
        state = dict(self.__dict__)
        state["_splines"] = None  # FITPACK splines are rebuilt in each process
        return state

    def _transform(self, d):
        if self._splines is None:
            self._build()
        d = np.clip(d, self.bounds[:, 0], self.bounds[:, 1])
        x = np.array([s(v) for (s, _), v in zip(self._splines, d)], dtype=float).ravel()
        dx = np.array([ds(v) for (_, ds), v in zip(self._splines, d)], dtype=float).ravel()
        return x, float(np.sum(np.log(np.maximum(dx, 1e-300))))

    def chi2(self, theory: dict, a_planck: float = 1.0) -> float:
        """−2 ln L in Planck's normalisation (zero at the reference spectrum)."""
        x, lj = self._transform(theory["tt"][self.lmin:self.lmax + 1] / a_planck**2)
        r = self.A @ (x - self.mu)
        return float(-2 * (lj - 0.5 * r @ r - self.offset))

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        x, lj = self._transform(theory["tt"][self.lmin:self.lmax + 1] / p.get("A_planck", 1.0)**2)
        r = self.A @ (x - self.mu)
        jac = self.K + 2 * self.offset - 2 * lj
        return np.concatenate([r, [np.sqrt(max(jac, 0.0))]])


class PlanckLowlEE:
    """Planck 2018 low-ℓ E-mode polarisation (SimAll), ℓ = 2-29.

    Planck tabulated ln L_ℓ(D_ℓ^EE) for each multipole on a grid of 1e-4 µK²;
    the total is the sum over ℓ. This is what pins the optical depth τ. Planck's
    code reads the table at the grid point below D_ℓ (a step function); we
    interpolate linearly through the middle of each step so the fitter sees a
    smooth likelihood. As residuals:
    sign(D_ℓ − D_ℓ,peak) √(2 (ln L_ℓ,max − ln L_ℓ(D_ℓ))), whose squares sum to
    −2 ln L up to the constant ``chi2_constant``.
    """

    name = "Planck low-ℓ EE"
    nuisance = ("A_planck",)
    lmin, lmax, step = 2, 29, 1e-4

    def __init__(self, table: np.ndarray):
        self.table = table
        self.grid = np.arange(table.shape[0]) * self.step
        self.peak = self.grid[np.argmax(table, axis=0)]
        self.lmax_val = table.max(axis=0)
        self.chi2_constant = float(-2 * self.lmax_val.sum())

    def lnl(self, d: np.ndarray, *, truncate: bool = False) -> np.ndarray:
        if truncate:  # exactly Planck's (and Cobaya's) lookup
            idx = np.clip((d / self.step).astype(int), 0, self.table.shape[0] - 1)
            return self.table[idx, np.arange(d.size)]
        # the step function is centred half a grid step above each table row
        return np.array([np.interp(v - 0.5 * self.step, self.grid, self.table[:, i])
                         for i, v in enumerate(d)])

    def chi2(self, theory: dict, a_planck: float = 1.0, *, truncate: bool = True) -> float:
        d = theory["ee"][self.lmin:self.lmax + 1] / a_planck**2
        return float(-2 * self.lnl(d, truncate=truncate).sum())

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        d = theory["ee"][self.lmin:self.lmax + 1] / p.get("A_planck", 1.0)**2
        drop = np.maximum(2 * (self.lmax_val - self.lnl(d)), 0.0)
        return np.sign(d - self.peak) * np.sqrt(drop)


# ---------------------------------------------------------------------------
# Priors
# ---------------------------------------------------------------------------


class GaussianPrior:
    nuisance = ()

    def __init__(self, param: str, mean: float, sigma: float, name: str | None = None):
        self.param, self.mean, self.sigma = param, mean, sigma
        self.name = name or f"prior on {param}"

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        return np.array([(p[self.param] - self.mean) / self.sigma])


# ---------------------------------------------------------------------------
# BAO
# ---------------------------------------------------------------------------


class BAO:
    """DESI DR2 BAO: D_V/r_d, D_M/r_d and D_H/r_d at seven redshifts, Gaussian."""

    nuisance = ()

    def __init__(self, data: dict, name: str = "DESI DR2 BAO"):
        self.name = name
        self.z, self.value, self.kind = data["z"], data["value"], data["kind"]
        self.A = whitening(data["cov"])
        zu = np.asarray(cfg.DESI_Z_UNIQUE)
        self.index = np.array([int(np.argmin(np.abs(zu - z))) for z in self.z])
        if np.max(np.abs(zu[self.index] - self.z)) > 1e-6:
            raise ValueError("BAO redshifts differ from config.DESI_Z_UNIQUE")

    def model(self, theory: dict) -> np.ndarray:
        dm = theory["bao_dm"][self.index]
        dh = cfg.SPEED_OF_LIGHT_KMS / theory["bao_hz"][self.index]
        rd = float(theory["derived"][list(_derived_keys()).index("rdrag")])
        dv = (self.z * dm**2 * dh) ** (1.0 / 3.0)
        out = np.where(self.kind == "DM_over_rs", dm, np.where(self.kind == "DH_over_rs", dh, dv))
        return out / rd

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        return self.A @ (self.model(theory) - self.value)


class BAO2018:
    """The BAO data Planck 2018 VI combined with the CMB: 6dFGS, SDSS DR7 MGS and
    BOSS DR12 consensus (BAO only), implemented as in Cobaya.

    * 6dF: r_s/D_V at z = 0.106, with r_s = 1.027369826 r_d (its own convention);
    * MGS: a tabulated χ² of α = D_V(0.15)/(4.29720761315 r_d), not a Gaussian;
    * DR12: D_M × (147.78/r_d) and H × (r_d/147.78) at z = 0.38, 0.51, 0.61.
    """

    name = "BAO 2018 (6dF + MGS + BOSS DR12)"
    nuisance = ()

    def __init__(self, data: dict):
        self.d = data
        self._spline = None
        m = data["mgs"]
        fine = np.linspace(*m["alpha_bounds"], 20001)
        vals = self.mgs_spline(fine)
        self.mgs_alpha0 = float(fine[np.argmin(vals)])
        self.mgs_chi2_min = float(vals.min())
        self.chi2_constant = self.mgs_chi2_min
        self.A12 = whitening(data["dr12"]["cov"])

    @property
    def mgs_spline(self):
        """Cubic spline through the MGS χ²(α) table (built lazily: SciPy's
        splines do not pickle, and profile workers receive this object pickled)."""
        if self._spline is None:
            from scipy.interpolate import CubicSpline

            m = self.d["mgs"]
            self._spline = CubicSpline(np.linspace(*m["alpha_bounds"], m["chi2"].size), m["chi2"])
        return self._spline

    def __getstate__(self):
        state = dict(self.__dict__)
        state["_spline"] = None
        return state

    @staticmethod
    def _dm_h(theory: dict, z):
        dm = np.interp(z, theory["sn_z"], theory["sn_dm"])
        hz = np.interp(z, theory["sn_z"], theory["sn_hz"])
        return dm, hz

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        rd = float(theory["derived"][list(_derived_keys()).index("rdrag")])
        c = cfg.SPEED_OF_LIGHT_KMS
        out = []
        s6 = self.d["sixdf"]
        dm, hz = self._dm_h(theory, s6["z"])
        dv = (s6["z"] * dm**2 * c / hz) ** (1 / 3)
        out.append((rd * s6["rs_rescale"] / dv - s6["value"]) / s6["err"])
        mg = self.d["mgs"]
        dm, hz = self._dm_h(theory, mg["z"])
        alpha = (mg["z"] * dm**2 * c / hz) ** (1 / 3) / (rd * mg["rs_rescale"])
        lo, hi = mg["alpha_bounds"]
        chi2 = float(self.mgs_spline(min(max(alpha, lo), hi))) - self.mgs_chi2_min
        out.append(np.sign(alpha - self.mgs_alpha0) * np.sqrt(max(chi2, 0.0)))
        d12 = self.d["dr12"]
        dm, hz = self._dm_h(theory, d12["z"])
        rs = d12["rs_fid"]
        model = np.where(d12["kind"] == "DM_over_rs", dm * rs / rd, hz * rd / rs)
        out.extend(self.A12 @ (model - d12["value"]))
        return np.array(out, dtype=float)


def _derived_keys():
    import cmb_model
    return cmb_model.DERIVED_KEYS


# ---------------------------------------------------------------------------
# Supernovae
# ---------------------------------------------------------------------------


class Supernovae:
    """Distance moduli with full covariance and an analytically profiled offset.

    μ_model = 5 log10((1 + z_hel) D_M(z_HD) / Mpc) + 25, as in DESI's (Cobaya's)
    likelihoods. The unknown offset (absolute magnitude combined with H0) is
    removed by projecting the whitened residual orthogonal to the whitened
    vector of ones, which equals marginalising it with a flat prior.
    """

    nuisance = ()

    def __init__(self, sample: dict):
        self.name = sample["label"]
        self.z_hd, self.z_hel, self.mu = sample["z_hd"], sample["z_hel"], sample["mu"]
        self.A = (whitening(sample["cov"]) if "cov" in sample
                  else whitening_from_inverse(sample["invcov"]))
        u = self.A @ np.ones(self.mu.size)
        self.u_hat = u / np.linalg.norm(u)
        self.wmu = self.A @ self.mu

    def residual(self, theory: dict, p: dict) -> np.ndarray:
        dm = np.interp(self.z_hd, theory["sn_z"], theory["sn_dm"])
        model = 5.0 * np.log10((1.0 + self.z_hel) * dm) + 25.0
        r = self.wmu - self.A @ model
        return r - self.u_hat * (self.u_hat @ r)
