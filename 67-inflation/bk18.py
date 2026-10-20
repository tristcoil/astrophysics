"""Our own BICEP/Keck BK18 likelihood, rebuilt from the released CosmoMC data.

BK18 measured the B-mode power in 66 auto- and cross-spectra between eleven
maps of the same patch of sky: BICEP3 at 95 GHz, Keck at 95, 150 and 220 GHz,
WMAP at 23 and 33 GHz, Planck at 30, 44, 143, 217 and 353 GHz. Each map sees
the same CMB B-modes (lensing + whatever gravitational waves there are) but
different amounts of polarized dust (bright at high frequency) and synchrotron
(bright at low frequency), so together they separate the three.

The model for the spectrum between maps i and j in bandpower bin b is

    C_b^ij = sum_l W_bl^ij [ D_l^lens + r D_l^tens
                             + A_d f_i f_j (l/80)^alpha_d
                             + A_s g_i g_j (l/80)^alpha_s
                             + eps sqrt(A_d A_s) (f_i g_j + g_i f_j) (l/80)^((alpha_d+alpha_s)/2) ]
             + N_b^ij,

with D_l = l(l+1)C_l/2pi in uK^2, bandpower window functions W (which include
the filtering of the maps), dust f_i = modified blackbody (beta_d, T_d = 19.6 K)
relative to 353 GHz, synchrotron g_i = power law (beta_s) relative to 23 GHz, both
integrated over each map's bandpass and converted to CMB temperature units, and
the noise bias N. The likelihood is the Hamimeche & Lewis (2008) approximation:
for each bin the 11x11 matrices of model (C) and data (C_hat) are compared via

    X_b = C_fid^(1/2) U g(D) U^T C_fid^(1/2),  C^(-1/2) C_hat C^(-1/2) = U D U^T,
    g(x) = sign(x - 1) sqrt(2 (x - ln x - 1)),

and -2 ln L = X^T M^-1 X with the released bandpower covariance M.

This follows the dataset definition (BK18lf_dust.dataset) and the CosmoMC /
Cobaya implementations; it is validated in ``validate_against_cobaya``.
"""

from __future__ import annotations

import io
import math
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import config as cfg

GHZ_KELVIN = 0.0479924466  # h / k_B in K/GHz
T_CMB = 2.7255
L_PIVOT = 80.0
DATA_PREFIX = "BK18_cosmomc/data/BK18lf_dust/"

#: Parameter order used by our sampler.
PARAMS = ("r", "BBdust", "BBsync", "BBalphadust", "BBalphasync", "BBbetadust",
          "BBbetasync", "BBdustsynccorr")
#: BK18 baseline priors (batch3/BK18lf.ini + BK18_04_BK18lf_freebdust.ini):
#: flat ranges; beta_sync also has a Gaussian prior -3.1 +- 0.3 (WMAP/Planck).
PRIOR_RANGES = {"r": (0.0, 0.5), "BBdust": (0.0, 15.0), "BBsync": (0.0, 50.0),
                "BBalphadust": (-1.0, 0.0), "BBalphasync": (-1.0, 0.0),
                "BBbetadust": (0.80, 2.40), "BBbetasync": (-4.5, -2.0),
                "BBdustsynccorr": (-1.0, 1.0)}
PRIOR_GAUSS = {"BBbetasync": (-3.1, 0.3)}
#: The cosmology BK18 fixed for its own analysis (batch3/BK18lfonly.ini):
#: Planck 2018 Plik best fit, tensors with n_t = 0 at k = 0.05 Mpc^-1.
BK18_COSMOLOGY = {"ombh2": 0.022383, "omch2": 0.12011, "tau": 0.0543, "ns": 0.96605,
                  "theta100": 1.040909, "logA": 3.0448}


# ---------------------------------------------------------------------------
# Reading the release
# ---------------------------------------------------------------------------


def _members(tgz: Path) -> dict[str, bytes]:
    out = {}
    with tarfile.open(tgz, "r:gz") as tar:
        for m in tar.getmembers():
            if m.isfile() and m.name.startswith(DATA_PREFIX):
                out[m.name[len(DATA_PREFIX):]] = tar.extractfile(m).read()
    return out


def _parse_dataset(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _table(blob: bytes) -> np.ndarray:
    return np.loadtxt(io.BytesIO(blob), comments="#")


@dataclass
class Bandpass:
    nu: np.ndarray
    resp: np.ndarray
    dnu: np.ndarray
    th_dust: float
    th_sync: float


def _bandpass(blob: bytes, f_dust: float, f_sync: float) -> Bandpass:
    a = _table(blob)
    nu, resp = a[:, 0], a[:, 1]
    dnu = np.hstack(((nu[1] - nu[0]), (nu[2:] - nu[:-2]) / 2, (nu[-1] - nu[-2])))

    def th(n):  # dB/dT at T_CMB, up to constants
        x = GHZ_KELVIN * n / T_CMB
        return n ** 4 * np.exp(x) / np.expm1(x) ** 2

    th_int = np.sum(dnu * resp * th(nu))
    return Bandpass(nu, resp, dnu, th_int / th(f_dust), th_int / th(f_sync))


def _parse_release(tgz: Path) -> dict:
    """Everything the likelihood needs, read out of BK18_cosmomc.tgz."""
    files = _members(tgz)
    ds = _parse_dataset(files["BK18lf_dust.dataset"].decode())
    maps = ds["map_names"].split()
    nbins = int(ds["nbins"])
    lmin, lmax = int(ds["cl_lmin"]), int(ds["cl_lmax"])
    idx = {m: i for i, m in enumerate(maps)}

    def pairs(order: str):
        return [tuple(idx[x] for x in spec.split("x")) for spec in order.split()]

    pr = pairs(ds["cl_hat_order"])
    for key in ("cl_fiducial_order", "cl_noise_order", "covmat_cl", "bin_window_in_order"):
        if pairs(ds[key]) != pr:
            raise ValueError(f"{key} differs from cl_hat_order; reorder needed")
    if ds.get("cl_hat_includes_noise", "F") != "F" or \
            ds.get("cl_fiducial_includes_noise", "F") != "F":
        raise ValueError("expected noise-free cl_hat and fiducial files")
    ncl = len(pr)

    def binned(name):
        a = _table(files[name])
        if a.shape != (nbins, ncl + 1):
            raise ValueError(f"{name}: shape {a.shape}")
        return a[:, 1:]

    cov = _table(files[ds["covmat_fiducial"]])
    if cov.shape != (nbins * ncl,) * 2:
        raise ValueError("covariance shape")
    # window functions W[b, spectrum, l] for l in [lmin, lmax]; any weight outside
    # that range is dropped, as CosmoMC/Cobaya do
    windows = np.zeros((nbins, ncl, lmax - lmin + 1))
    dropped = 0.0
    for b in range(nbins):
        w = _table(files[ds["bin_window_files"].replace("%u", str(b + 1))])
        ell = w[:, 0].astype(int)
        use = (ell >= lmin) & (ell <= lmax)
        windows[b][:, ell[use] - lmin] = w[use, 1:].T
        dropped = max(dropped, float(np.abs(w[~use, 1:]).sum() / np.abs(w[:, 1:]).sum()))
    f_dust, f_sync = float(ds["fpivot_dust"]), float(ds["fpivot_sync"])
    bps = [_bandpass(files[ds[f"bandpass[{m}]"]], f_dust, f_sync) for m in maps]
    return {
        "maps": np.array(maps), "nbins": nbins, "lmin": lmin, "lmax": lmax,
        "pi": np.array([a for a, _ in pr]), "pj": np.array([b for _, b in pr]),
        "cl_hat": binned(ds["cl_hat_file"]), "noise": binned(ds["cl_noise_file"]),
        "fiducial": binned(ds["cl_fiducial_file"]), "cov": cov, "windows": windows,
        "window_weight_dropped": dropped, "f_dust": f_dust, "f_sync": f_sync,
        "bp_nu": np.concatenate([bp.nu for bp in bps]),
        "bp_w": np.concatenate([bp.dnu * bp.resp for bp in bps]),
        "bp_start": np.cumsum([0] + [bp.nu.size for bp in bps])[:-1],
        "th_dust": np.array([bp.th_dust for bp in bps]),
        "th_sync": np.array([bp.th_sync for bp in bps]),
    }


def load_release(tgz: Path = cfg.BK18_DIR / "BK18_cosmomc.tgz") -> dict:
    """The parsed release, cached in data/derived (re-parsed if the tarball changes)."""
    cache = cfg.DERIVED_DIR / "bk18_parsed.npz"
    stamp = np.array([tgz.stat().st_size, int(tgz.stat().st_mtime)])
    if cache.exists():
        with np.load(cache) as z:
            if np.array_equal(z["stamp"], stamp):
                return {k: z[k] for k in z.files if k != "stamp"}
    d = _parse_release(tgz)
    np.savez_compressed(cache, stamp=stamp, **d)
    return d


class BK18:
    """The BK18 B-mode likelihood (11 maps, 9 bins, 66 spectra)."""

    def __init__(self, tgz: Path = cfg.BK18_DIR / "BK18_cosmomc.tgz"):
        d = load_release(tgz)
        self.maps = [str(m) for m in d["maps"]]
        self.nmap = len(self.maps)
        self.nbins = int(d["nbins"])
        self.lmin, self.lmax = int(d["lmin"]), int(d["lmax"])
        self.ells = np.arange(self.lmin, self.lmax + 1)
        self.f_dust, self.f_sync = float(d["f_dust"]), float(d["f_sync"])
        self.pi, self.pj = d["pi"], d["pj"]
        self.ncl = self.pi.size
        self.cl_hat, self.noise, self.fiducial = d["cl_hat"], d["noise"], d["fiducial"]
        self.cov = d["cov"]
        self.inv_cov = np.linalg.inv(self.cov)
        self.windows = d["windows"]
        self.window_weight_dropped = float(d["window_weight_dropped"])
        self._bp_nu, self._bp_w, self._bp_start = d["bp_nu"], d["bp_w"], d["bp_start"]
        self._th_dust, self._th_sync = d["th_dust"], d["th_sync"]
        # matrices of the data and of the fiducial model, both with noise added
        self.chat = self._matrices(self.cl_hat + self.noise)
        cfid = self._matrices(self.fiducial + self.noise)
        self.cfid_half = np.array([_sqrtm(c) for c in cfid])
        self.lrat = self.ells / L_PIVOT
        self._wflat = self.windows.reshape(self.nbins * self.ncl, -1)
        self.cmb_binned = None

    # -- theory ---------------------------------------------------------------

    def _matrices(self, vec: np.ndarray) -> np.ndarray:
        """(nbins, ncl) spectra -> (nbins, nmap, nmap) symmetric matrices."""
        m = np.zeros((self.nbins, self.nmap, self.nmap))
        m[:, self.pi, self.pj] = vec
        m[:, self.pj, self.pi] = vec
        return m

    def bin(self, dl: np.ndarray) -> np.ndarray:
        """Bin a spectrum D_l (l = lmin..lmax) through every window: (nbins, ncl)."""
        return (self._wflat @ dl).reshape(self.nbins, self.ncl)

    def set_cmb(self, dl_lens: np.ndarray, dl_tens_r1: np.ndarray) -> None:
        """Lensed scalar BB and tensor BB for r = 1, as D_l on l = lmin..lmax (uK^2)."""
        self.dl_lens = dl_lens[self.lmin:self.lmax + 1] if dl_lens.size > self.ells.size \
            else dl_lens
        self.dl_tens = dl_tens_r1[self.lmin:self.lmax + 1] if dl_tens_r1.size > self.ells.size \
            else dl_tens_r1
        self.cmb_binned = (self.bin(self.dl_lens), self.bin(self.dl_tens))

    def dust_scaling(self, beta: float, t_dust: float = 19.6) -> np.ndarray:
        """Modified-blackbody dust in each map relative to 353 GHz (CMB temperature units)."""
        nu0 = self.f_dust
        gb0 = nu0 ** (3 + beta) / np.expm1(GHZ_KELVIN * nu0 / t_dust)
        nu = self._bp_nu
        gb = np.add.reduceat(self._bp_w * nu ** (3 + beta) / np.expm1(GHZ_KELVIN * nu / t_dust),
                             self._bp_start)
        return gb / gb0 / self._th_dust

    def sync_scaling(self, beta: float) -> np.ndarray:
        """Power-law synchrotron in each map relative to 23 GHz (CMB temperature units)."""
        pl0 = self.f_sync ** (2 + beta)
        pl = np.add.reduceat(self._bp_w * self._bp_nu ** (2 + beta), self._bp_start)
        return pl / pl0 / self._th_sync

    def model(self, p: dict, cmb: tuple | None = None) -> np.ndarray:
        """Model bandpowers (nbins, ncl) without noise."""
        lens, tens = cmb if cmb is not None else self.cmb_binned
        fd = self.dust_scaling(p["BBbetadust"], p.get("BBTdust", 19.6))
        fs = self.sync_scaling(p["BBbetasync"])
        ad, asy = p["BBdust"], p["BBsync"]
        dust = self.bin(ad * self.lrat ** p["BBalphadust"])
        sync = self.bin(asy * self.lrat ** p["BBalphasync"])
        out = lens + p["r"] * tens
        out = out + dust * (fd[self.pi] * fd[self.pj]) + sync * (fs[self.pi] * fs[self.pj])
        eps = p["BBdustsynccorr"]
        if eps != 0 and ad > 0 and asy > 0:
            cross = self.bin(eps * math.sqrt(ad * asy)
                             * self.lrat ** (0.5 * (p["BBalphadust"] + p["BBalphasync"])))
            out = out + cross * (fd[self.pi] * fs[self.pj] + fs[self.pi] * fd[self.pj])
        return out

    # -- likelihood -------------------------------------------------------------

    def chi2_from_model(self, model: np.ndarray) -> float:
        """-2 ln L (Hamimeche-Lewis) for model bandpowers (nbins, ncl) without noise."""
        x = _hl_batch(self._matrices(model + self.noise), self.chat, self.cfid_half)
        if x is None:
            return np.inf
        v = x[:, self.pi, self.pj].ravel()
        return float(v @ self.inv_cov @ v)

    def chi2(self, p: dict, cmb: tuple | None = None) -> float:
        return self.chi2_from_model(self.model(p, cmb))

    def log_prior(self, p: dict) -> float:
        for k, (lo, hi) in PRIOR_RANGES.items():
            if not lo <= p[k] <= hi:
                return -np.inf
        lp = 0.0
        for k, (mu, sig) in PRIOR_GAUSS.items():
            lp -= 0.5 * ((p[k] - mu) / sig) ** 2
        return lp

    def log_post(self, theta: np.ndarray) -> float:
        p = dict(zip(PARAMS, theta))
        lp = self.log_prior(p)
        if not np.isfinite(lp):
            return -np.inf
        return lp - 0.5 * self.chi2(p)


def _sqrtm(c: np.ndarray) -> np.ndarray:
    d, u = np.linalg.eigh(c)
    return (u * np.sqrt(d)) @ u.T


def _hl_batch(c: np.ndarray, chat: np.ndarray, cfhalf: np.ndarray):
    """Hamimeche-Lewis transform of every bin's matrices at once (as in CMBlikes):
    X = C_fid^(1/2) g(C^(-1/2) C_hat C^(-1/2)) C_fid^(1/2), g acting on eigenvalues."""
    d, u = np.linalg.eigh(c)
    if np.any(d <= 0):
        return None
    ut = np.swapaxes(u, 1, 2)
    isq = (u / np.sqrt(d)[:, None, :]) @ ut           # C^(-1/2)
    m = isq @ chat @ isq
    e, v = np.linalg.eigh(m)
    if np.any(e <= 0):
        return None
    g = np.sign(e - 1) * np.sqrt(2 * np.maximum(0.0, e - np.log(e) - 1))
    return cfhalf @ ((v * g[:, None, :]) @ np.swapaxes(v, 1, 2)) @ cfhalf


# ---------------------------------------------------------------------------
# CMB B-mode spectra from CAMB
# ---------------------------------------------------------------------------


def camb_bb(cosmo: dict = BK18_COSMOLOGY, *, lmax: int = 1500, r: float = 1.0,
            nt: float | None = 0.0, lens_accuracy: int = 1) -> dict:
    """Lensed scalar BB and tensor BB (for ``r``) as D_l in uK^2, l = 0..lmax."""
    import camb

    pars = camb.CAMBparams()
    pars.set_cosmology(ombh2=cosmo["ombh2"], omch2=cosmo["omch2"], tau=cosmo["tau"],
                       cosmomc_theta=cosmo["theta100"] / 100.0, mnu=0.06, num_massive_neutrinos=1)
    pars.InitPower.set_params(As=math.exp(cosmo["logA"]) * 1e-10, ns=cosmo["ns"], r=r, nt=nt,
                              pivot_scalar=0.05, pivot_tensor=0.05)
    pars.WantTensors = True
    pars.set_for_lmax(lmax + 500, lens_potential_accuracy=lens_accuracy)
    pars.NonLinear = camb.model.NonLinear_lens
    res = camb.get_results(pars)
    cls = res.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)
    return {"lensed_scalar_bb": cls["lensed_scalar"][:lmax + 1, 2],
            "tensor_bb": cls["tensor"][:lmax + 1, 2],
            "tensor_ee": cls["tensor"][:lmax + 1, 1],
            "unlensed_scalar_ee": cls["unlensed_scalar"][:lmax + 1, 1],
            "lensed_scalar_ee": cls["lensed_scalar"][:lmax + 1, 1]}


# ---------------------------------------------------------------------------
# Sampling r and the seven foreground parameters
# ---------------------------------------------------------------------------

_LIKE: BK18 | None = None  # set before forking the worker pool


def _log_post(theta: np.ndarray) -> float:
    return _LIKE.log_post(theta)


def sample(like: BK18, *, nwalkers: int = 48, nsteps: int = 12000, burn: int = 2000,
           thin: int = 10, seed: int = 20261019, processes: int = 6) -> dict:
    """emcee over (r, A_d, A_s, alpha_d, alpha_s, beta_d, beta_s, eps) with BK18's priors."""
    import multiprocessing as mp

    import emcee

    global _LIKE
    _LIKE = like
    rng = np.random.default_rng(seed)
    start = np.array([0.015, 4.4, 1.0, -0.5, -0.5, 1.5, -3.1, 0.0])
    scale = np.array([0.006, 0.3, 0.4, 0.08, 0.08, 0.06, 0.15, 0.15])
    p0 = start + scale * rng.standard_normal((nwalkers, len(PARAMS)))
    for j, k in enumerate(PARAMS):  # keep every walker inside the prior box
        lo, hi = PRIOR_RANGES[k]
        p0[:, j] = np.clip(p0[:, j], lo + 0.02 * (hi - lo), hi - 0.02 * (hi - lo))
    with mp.get_context("fork").Pool(processes) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, len(PARAMS), _log_post, pool=pool)
        sampler.random_state = np.random.RandomState(seed).get_state()  # emcee seeding trap
        sampler.run_mcmc(p0, nsteps, progress=False)
    chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
    logp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    tau = sampler.get_autocorr_time(discard=burn, tol=0)
    return {"chain": chain, "log_post": logp, "tau": tau,
            "acceptance": float(np.mean(sampler.acceptance_fraction)),
            "nwalkers": nwalkers, "nsteps": nsteps, "burn": burn, "thin": thin}


def r_likelihood_curve(r: np.ndarray, *, grid: np.ndarray, bandwidth: float) -> np.ndarray:
    """Smoothed marginal posterior of r (flat prior, so proportional to the likelihood),
    with a reflection boundary at r = 0; normalised to a peak of 1 on ``grid``."""
    d = grid[1] - grid[0]
    edges = np.concatenate([grid - d / 2, [grid[-1] + d / 2]])
    h, _ = np.histogram(np.concatenate([r, -r]), np.concatenate([-edges[::-1], edges[1:]]))
    from scipy.ndimage import gaussian_filter1d
    sm = gaussian_filter1d(h.astype(float), bandwidth / d, mode="constant")
    sm = sm[grid.size:]
    return sm / sm.max()


def published_r_likelihood(path: Path = cfg.BK18_DIR / "BK18_r_likelihood_20210607.txt"):
    a = np.loadtxt(path, comments="#")
    return a[:, 0], a[:, 1]


# ---------------------------------------------------------------------------
# The marginal r posterior without a smoothing kernel (Rao-Blackwell)
# ---------------------------------------------------------------------------


def _conditional_worker(args):
    theta, r_grid = args
    p = dict(zip(PARAMS, theta))
    lens, tens = _LIKE.cmb_binned
    base = _LIKE.model({**p, "r": 0.0})
    chi2 = np.array([_LIKE.chi2_from_model(base + r * tens) for r in r_grid])
    return chi2


def rao_blackwell_r(like: BK18, chain: np.ndarray, r_grid: np.ndarray, *, n_use: int = 3000,
                    processes: int = 6, seed: int = 7) -> dict:
    """p(r | data) = mean over posterior foreground samples of p(r | foregrounds, data).

    For each sampled set of foreground parameters the likelihood is evaluated
    on the whole r grid and normalised (flat prior on r); averaging these exact
    conditionals estimates the marginal posterior with no smoothing kernel, so
    there is no bias at the r = 0 boundary.
    """
    import multiprocessing as mp

    global _LIKE
    _LIKE = like
    rng = np.random.default_rng(seed)
    idx = rng.choice(chain.shape[0], size=min(n_use, chain.shape[0]), replace=False)
    with mp.get_context("fork").Pool(processes) as pool:
        chi2 = np.array(pool.map(_conditional_worker, [(chain[i], r_grid) for i in idx],
                                 chunksize=16))
    logl = -0.5 * chi2
    logl -= logl.max(axis=1, keepdims=True)
    like_r = np.exp(logl)
    norm = np.trapezoid(like_r, r_grid, axis=1)
    cond = like_r / norm[:, None]
    post = cond.mean(axis=0)
    # bootstrap over the foreground samples for the Monte Carlo error
    boot = np.array([cond[rng.integers(0, len(idx), len(idx))].mean(axis=0) for _ in range(200)])
    return {"r_grid": r_grid, "posterior": post, "boot_std": boot.std(axis=0),
            "n_used": len(idx), "boot_r95": np.array([_limit(r_grid, b) for b in boot])}


def _limit(grid: np.ndarray, p: np.ndarray, level: float = 0.95) -> float:
    c = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(grid))])
    return float(np.interp(level, c / c[-1], grid))
