"""Gaia DR3: the Leavitt law's zero point from Milky Way Cepheid parallaxes.

The Magellanic Clouds give the *shape* of the Leavitt law (slope, scatter)
with every star at nearly the same distance, but not its *zero point*: how
luminous a 10-day Cepheid really is. Parallaxes of nearby Milky Way Cepheids
give that, and comparing the two gives the absolute distance to the LMC.

To compare like with like, both sides use Gaia's own photometry, in the
reddening-free Gaia Wesenheit W_G = G - 1.90 (BP - RP), from the light-curve
(intensity) averages of Gaia's Cepheid pipeline (``vari_cepheid``).

The fit is done in **parallax space**: each star's Wesenheit magnitude and
period predict a parallax, compared with the measured one,

    plx_pred = 10^(-0.2 (W - M_W - 10)) mas,   M_W = alpha + beta (log P - 1) + gamma [Fe/H],
    plx_obs - zp_L21 - delta  ~  N(plx_pred, (k sigma_plx)^2 + (0.2 ln10 plx_pred sigma_int)^2),

which keeps low-S/N and even negative parallaxes valid and so avoids the
Lutz-Kelker bias that cutting on parallax S/N or inverting parallaxes brings.
zp_L21 is the Lindegren et al. (2021) zero point; delta is any residual
offset, fitted (Riess et al. 2021 found one for bright Cepheids); k inflates
the catalogue parallax errors; sigma_int is the law's intrinsic scatter.
"""

from __future__ import annotations

import warnings

import numpy as np
from astropy.table import Table

import config as C
import ogle

LN10_5 = 0.2 * np.log(10)


# ---------------------------------------------------------------------------
# Photometry helpers
# ---------------------------------------------------------------------------

def excess_cstar(bp_rp, c):
    """Riello et al. (2021) corrected BP/RP flux excess C* (their eq. 6, Table 2).

    For an isolated source C* scatters around 0; light from neighbours inside
    the BP/RP windows pushes it positive.
    """
    x = np.asarray(bp_rp, float)
    f = np.where(x < 0.5, 1.154360 + 0.033772 * x + 0.032277 * x ** 2,
                 np.where(x < 4.0, 1.162004 + 0.011464 * x + 0.049255 * x ** 2 - 0.005879 * x ** 3,
                          1.057572 + 0.140537 * x))
    return np.asarray(c, float) - f


def sigma_cstar(g):
    """1-sigma scatter of C* for well-behaved isolated sources (Riello et al. 2021, eq. 18)."""
    return 0.0059898 + 8.817481e-12 * np.asarray(g, float) ** 7.618399


def wesenheit_g(g, bp, rp):
    return np.asarray(g, float) - C.R_WG * (np.asarray(bp, float) - np.asarray(rp, float))


def _col(tab, name, fill=np.nan):
    """A table column as a plain array: masked entries filled, byte strings decoded."""
    col = tab[name]
    kind = col.dtype.kind
    if kind == "S":
        col = np.char.decode(np.ma.filled(col, b"") if hasattr(col, "mask") else np.asarray(col), "ascii")
        return np.char.strip(col)
    if kind in "iu" and isinstance(fill, float) and np.isnan(fill):
        col = col.astype(float)
    if hasattr(col, "mask"):
        return np.asarray(np.ma.filled(col, fill))
    return np.asarray(col)


# ---------------------------------------------------------------------------
# Magellanic Cepheids in Gaia
# ---------------------------------------------------------------------------

def magellanic() -> dict:
    """Gaia DR3 counterparts of the OGLE Magellanic Cepheids, aligned with ``ogle.catalog`` rows.

    The counterpart is the nearest Gaia source (median separation 0.09").
    ``n_nb`` counts every Gaia source within 2", a crowding diagnostic.
    """
    x = Table.read(C.GAIA_DIR / "xmatch_ogle_mc_gaiadr3.csv", format="ascii.csv")
    x.sort(["id", "angDist"])
    ids = np.asarray(x["id"]).astype(str)
    xcol = {key: _col(x, col).astype(float) for key, col in
            (("sep", "angDist"), ("G", "Gmag"), ("BP", "BPmag"), ("RP", "RPmag"), ("C", "E(BP/RP)"),
             ("ruwe", "RUWE"), ("plx", "Plx"), ("e_plx", "e_Plx"))}
    xsrc = _col(x, "Source", 0).astype(np.int64)
    # first (nearest) row of each OGLE id, and the number of Gaia sources within the radius
    uid, first, count = np.unique(ids, return_index=True, return_counts=True)
    pos = {s: (f, c) for s, f, c in zip(uid, first, count)}
    v = Table.read(C.GAIA_DIR / "vari_cepheid_x_gaia_source.fits")
    vrow = {int(s): k for k, s in enumerate(_col(v, "source_id", 0))}
    vcol = {key: _col(v, col) for key, col in
            (("G_int", "int_average_g"), ("BP_int", "int_average_bp"), ("RP_int", "int_average_rp"),
             ("gaia_type", "type_best_classification"), ("gaia_mode", "mode_best_classification"))}
    out = {}
    for gal in C.GALAXIES:
        cat = ogle.catalog(gal)
        f_c = np.array([pos.get(s, (-1, 0)) for s in cat["id"]])
        r, n_nb = f_c[:, 0], f_c[:, 1]
        has = r >= 0
        res = {key: np.where(has, col[np.maximum(r, 0)], np.nan) for key, col in xcol.items()}
        res["n_nb"] = n_nb
        res["source_id"] = np.where(has, xsrc[np.maximum(r, 0)], 0)
        # magnitude difference to the second-nearest source, if any
        r2 = np.where(n_nb > 1, r + 1, 0)
        res["nb_dG"] = np.where(n_nb > 1, xcol["G"][r2] - res["G"], np.nan)
        vk = np.array([vrow.get(int(s), -1) for s in res["source_id"]])
        inv = vk >= 0
        for key in ("G_int", "BP_int", "RP_int"):
            res[key] = np.where(inv, vcol[key][np.maximum(vk, 0)].astype(float), np.nan)
        for key in ("gaia_type", "gaia_mode"):
            res[key] = np.where(inv, vcol[key][np.maximum(vk, 0)], "")
        res["W_G"] = wesenheit_g(res["G_int"], res["BP_int"], res["RP_int"])
        res["cstar"] = excess_cstar(res["BP"] - res["RP"], res["C"])
        res["cstar_sig"] = res["cstar"] / sigma_cstar(res["G"])
        out[gal] = res
    return out


# ---------------------------------------------------------------------------
# Milky Way calibrators
# ---------------------------------------------------------------------------

def lindegren_zero_point(g, nu_eff, pseudocolour, ecl_lat, solved):
    """Lindegren et al. (2021) parallax zero point in mas (NaN outside its calibrated range)."""
    from zero_point import zpt
    zpt.load_tables()
    solved = np.asarray(solved, int)
    ok = np.isin(solved, (31, 95))  # 2-parameter solutions have no parallax
    out = np.full(len(solved), np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out[ok] = zpt.get_zpt(np.asarray(g, float)[ok], np.asarray(nu_eff, float)[ok],
                              np.asarray(pseudocolour, float)[ok], np.asarray(ecl_lat, float)[ok],
                              solved[ok], _warnings=False)
    return out


def milky_way() -> dict:
    """OGLE's Galactic Cepheid list joined to Gaia DR3 astrometry and Cepheid photometry."""
    lst = ogle.galactic_list()
    t = Table.read(C.GAIA_DIR / "galactic_cepheids_gaia.fits")
    row = {int(s): k for k, s in enumerate(_col(t, "source_id", 0))}
    k_of = np.array([row.get(int(s), -1) for s in lst["source_id"]])
    have = k_of >= 0
    out = {key: lst[key][have] for key in ("name", "source", "mode", "P_F", "P_1O", "P_2O", "l", "b",
                                           "ra", "dec", "source_id")}
    k = k_of[have]

    def g(name, fill=np.nan):
        return _col(t, name, fill)[k]

    for name in ("parallax", "parallax_error", "ruwe", "phot_g_mean_mag", "phot_bp_mean_mag",
                 "phot_rp_mean_mag", "phot_bp_rp_excess_factor", "nu_eff_used_in_astrometry", "pseudocolour",
                 "ecl_lat", "int_average_g", "int_average_bp", "int_average_rp", "pf", "p1_o",
                 "ipd_frac_multi_peak", "radial_velocity", "metallicity", "metallicity_error",
                 "peak_to_peak_g", "pmra", "pmdec"):
        out[name] = g(name).astype(float)
    out["solved"] = g("astrometric_params_solved", 0).astype(int)
    out["gaia_type"] = g("type_best_classification", "").astype(str)
    out["gaia_mode"] = g("mode_best_classification", "").astype(str)
    out["zp_l21"] = lindegren_zero_point(out["phot_g_mean_mag"], out["nu_eff_used_in_astrometry"],
                                         out["pseudocolour"], out["ecl_lat"], out["solved"])
    out["W_G"] = wesenheit_g(out["int_average_g"], out["int_average_bp"], out["int_average_rp"])
    out["bp_rp"] = out["int_average_bp"] - out["int_average_rp"]
    out["cstar"] = excess_cstar(out["phot_bp_mean_mag"] - out["phot_rp_mean_mag"], out["phot_bp_rp_excess_factor"])
    out["cstar_sig"] = out["cstar"] / sigma_cstar(out["phot_g_mean_mag"])
    # the period of the dominant mode, and that mode, from OGLE's list
    out["P"] = np.where(np.char.startswith(out["mode"].astype(str), "F"), out["P_F"], out["P_1O"])
    out["logP"] = np.log10(out["P"])
    return out


def calibrator_mask(mw: dict, mode: str, *, g_min=C.G_MIN_L21, g_max=C.G_MAX, ruwe_max=C.RUWE_MAX,
                    require_gaia_dcep=False) -> np.ndarray:
    """Single-mode Cepheids of ``mode`` ('F' or '1O') with clean astrometry and photometry."""
    m = (mw["mode"] == mode) & np.isfinite(mw["W_G"]) & np.isfinite(mw["parallax"])
    m &= np.isfinite(mw["zp_l21"]) & (mw["ruwe"] < ruwe_max)
    m &= (mw["phot_g_mean_mag"] > g_min) & (mw["phot_g_mean_mag"] < g_max)
    m &= np.isin(mw["solved"], (31, 95))
    if require_gaia_dcep:
        m &= mw["gaia_type"] == "DCEP"
    return m


# ---------------------------------------------------------------------------
# The parallax-space fit
# ---------------------------------------------------------------------------

#: Parameters of the model, their defaults when held fixed, and starting values.
#:   alpha  absolute W at log P = 1 and [Fe/H] = 0 (mag)
#:   beta   slope (mag/dex), normally fixed to the LMC's
#:   R      Wesenheit coefficient W = G - R (BP - RP)
#:   gamma  metallicity term (mag/dex)
#:   delta  residual parallax offset after Lindegren et al. (mas; subtracted)
#:   lnsig  ln of the intrinsic scatter (mag)
#:   k      inflation of the catalogue parallax errors
DEFAULTS = {"alpha": -5.9, "beta": -3.3, "R": C.R_WG, "gamma": 0.0, "delta": 0.0, "lnsig": np.log(0.1), "k": 1.0}


class ParallaxFit:
    """Likelihood of Milky Way Cepheid parallaxes given a Leavitt law in W_G.

    ``free`` names the fitted parameters; the rest take ``fixed`` or DEFAULTS.
    """

    def __init__(self, g, bp, rp, logp, plx, e_plx, zp, *, feh=None, free=("alpha", "delta", "lnsig"),
                 fixed=None, x0=C.LOGP_PIVOT):
        self.g, self.bp, self.rp = (np.asarray(a, float) for a in (g, bp, rp))
        self.logp, self.plx, self.e_plx, self.zp = (np.asarray(a, float) for a in (logp, plx, e_plx, zp))
        self.feh = np.zeros_like(self.g) if feh is None else np.asarray(feh, float)
        self.free = tuple(free)
        self.fixed = {**DEFAULTS, **(fixed or {})}
        self.x0 = x0
        self.mask = np.ones(len(self.g), bool)

    def params(self, theta):
        p = dict(self.fixed)
        p.update(zip(self.free, theta))
        return p

    def predict(self, p):
        w = self.g - p["R"] * (self.bp - self.rp)
        m = p["alpha"] + p["beta"] * (self.logp - self.x0) + p["gamma"] * self.feh
        return 10 ** (-0.2 * (w - m - 10.0))

    def terms(self, theta):
        p = self.params(theta)
        pred = self.predict(p)
        var = (p["k"] * self.e_plx) ** 2 + (LN10_5 * pred * np.exp(p["lnsig"])) ** 2
        r = self.plx - self.zp - p["delta"] - pred
        return r, var, pred

    def nll(self, theta):
        p = self.params(theta)
        if "k" in self.free and p["k"] <= 0.3:
            return np.inf
        if p["lnsig"] < np.log(1e-3):
            return np.inf
        r, var, _ = self.terms(theta)
        m = self.mask
        return 0.5 * np.sum(r[m] ** 2 / var[m] + np.log(var[m]))

    def start(self):
        return np.array([self.fixed[name] for name in self.free])

    def fit(self, clip=4.0, iters=10):
        """MAP fit with iterative clipping of > ``clip``-sigma parallax outliers."""
        from scipy.optimize import minimize
        theta = self.start()
        for _ in range(iters):
            res = minimize(self.nll, theta, method="Nelder-Mead",
                           options={"xatol": 1e-7, "fatol": 1e-7, "maxiter": 40000, "maxfev": 40000})
            res = minimize(self.nll, res.x, method="Powell", options={"xtol": 1e-8, "ftol": 1e-10})
            theta = res.x
            r, var, _ = self.terms(theta)
            new = np.abs(r / np.sqrt(var)) < clip
            if np.array_equal(new, self.mask):
                break
            self.mask = new
        self.theta = theta
        return theta

    def sample(self, nwalkers=32, nsteps=3000, burn=1000, seed=C.SEED):
        """emcee posterior around the MAP (outlier mask frozen)."""
        import emcee
        ndim = len(self.free)
        rng = np.random.default_rng(seed)
        p0 = self.theta + 1e-4 * rng.standard_normal((nwalkers, ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lambda t: -self.nll(t))
        sampler.random_state = np.random.RandomState(seed).get_state()  # emcee copies the global RNG at init
        sampler.run_mcmc(p0, nsteps, progress=False)
        self.chain = sampler.get_chain(discard=burn, flat=True)
        try:
            self.tau = sampler.get_autocorr_time(discard=burn, quiet=True)
        except Exception:
            self.tau = np.full(ndim, np.nan)
        return self.chain

    def summary(self) -> dict:
        out = {"n": int(self.mask.sum()), "n_clipped": int((~self.mask).sum()), "free": list(self.free)}
        for j, name in enumerate(self.free):
            q = np.percentile(self.chain[:, j], [16, 50, 84])
            out[name] = {"map": float(self.theta[j]), "median": float(q[1]), "lo": float(q[1] - q[0]),
                         "hi": float(q[2] - q[1]), "std": float(np.std(self.chain[:, j]))}
        r, var, pred = self.terms(self.theta)
        m = self.mask
        out["chi2"] = float(np.sum(r[m] ** 2 / var[m]))
        out["sigma_int"] = float(np.exp(self.params(self.theta)["lnsig"]))
        return out


# ---------------------------------------------------------------------------
# The joint model: the LMC's distance as a parameter of the Milky Way fit
# ---------------------------------------------------------------------------

#: Mean reddening of LMC Cepheids, E(BP-RP) ~ E(V-I) ~ 0.10 (Skowron et al. 2021
#: maps). Only used to set the zero of each Milky Way star's reddening estimate
#: in the dust-law scatter term, so its precision hardly matters.
E_BPRP_LMC = 0.10


class LMCShape:
    """The LMC's Leavitt law in Gaia bands for any Wesenheit coefficient R.

    G and (BP - RP) are each fitted against log P on one fixed, clipped sample
    (the clipping is done once in W with R = 1.90), so the Wesenheit relation
    for any R follows exactly: a_W(R) = a_G - R a_c, b_W(R) = b_G - R b_c.
    The magnitudes are first corrected to the LMC centre with our disk fit.
    """

    def __init__(self, logp, g, bp, rp, geo_corr=None, keep=None):
        import plfit
        logp, g, bp, rp = (np.asarray(a, float) for a in (logp, g, bp, rp))
        corr = np.zeros_like(g) if geo_corr is None else np.asarray(geo_corr, float)
        g = g - corr  # distance correction applies equally to every band
        w = g - C.R_WG * (bp - rp)
        fw = plfit.fit_line(logp, w, keep=keep)
        m = fw["mask"]
        fg = plfit.fit_line(logp, g, keep=m, nsigma=1e9)
        fc = plfit.fit_line(logp, bp - rp, keep=m, nsigma=1e9)
        self.n = int(m.sum())
        self.a_g, self.b_g, self.a_c, self.b_c = fg["a"], fg["b"], fc["a"], fc["b"]
        self.sigma_w = fw["sigma"]
        self.e_a = fw["ea"]  # zero-point error of the W relation (nearly R-independent)
        self.fit_w190 = fw

    def a(self, R):
        return self.a_g - R * self.a_c

    def b(self, R):
        return self.b_g - R * self.b_c

    def intrinsic_colour(self, logp):
        """(BP - RP)_0 of an LMC Cepheid at this period (LMC colour minus its mean reddening)."""
        return self.a_c + self.b_c * (np.asarray(logp) - C.LOGP_PIVOT) - E_BPRP_LMC


#: Joint-model parameters, defaults (when fixed) and starting values.
JOINT_DEFAULTS = {"mu": 18.5, "gamma": -0.4, "delta": 0.015, "R": C.R_WG, "lnsig": np.log(0.08),
                  "sigR": 0.1, "k": 1.3, "d1O": 0.0}


class JointFit:
    """Milky Way parallaxes given the LMC's Leavitt law shifted by its distance.

    For each Milky Way Cepheid of mode m (F or 1O):

        M_W(R) = a_m(R) + b_m(R) (log P - 1) - mu_LMC + gamma ([Fe/H] - [Fe/H]_LMC) [+ d1O for 1O]
        W(R)   = G - R (BP - RP)
        plx_pred = 10^(-0.2 (W - M_W - 10))  mas
        plx_obs - zp_L21 - delta ~ N(plx_pred, (k s_plx)^2 + (0.2 ln10 plx_pred)^2 (sig^2 + (sigR E)^2))

    E is the star's reddening, its BP - RP minus the LMC's intrinsic colour at
    the same period; sigR is the sightline-to-sightline scatter of the dust law
    (in units of R), so heavily reddened stars are trusted less. Priors are
    uniform: 0 <= sigR < 1, 0.3 < k < 5, 1 < R < 2.8.
    """

    def __init__(self, mw: dict, sel: np.ndarray, feh: np.ndarray, shapes: dict, feh_lmc: float,
                 free=("mu", "gamma", "delta", "R", "lnsig", "sigR", "k"), fixed=None):
        self.sel = np.asarray(sel, bool)
        s = self.sel
        self.g, self.bp, self.rp = mw["int_average_g"][s], mw["int_average_bp"][s], mw["int_average_rp"][s]
        self.logp, self.feh = mw["logP"][s], np.asarray(feh, float)[s]
        self.plx, self.e_plx, self.zp = mw["parallax"][s], mw["parallax_error"][s], mw["zp_l21"][s]
        self.mode = mw["mode"][s]
        self.is1o = self.mode == "1O"
        self.shapes = shapes
        self.feh_lmc = feh_lmc
        self.free = tuple(free)
        self.fixed = {**JOINT_DEFAULTS, **(fixed or {})}
        self.mask = np.ones(int(s.sum()), bool)
        ext = {}
        for m in ("F", "1O"):
            mm = self.mode == m
            if mm.any():
                ext[m] = self.bp[mm] - self.rp[mm] - shapes[m].intrinsic_colour(self.logp[mm])
        self.E = np.zeros_like(self.g)
        for m, v in ext.items():
            self.E[self.mode == m] = np.clip(v, 0, None)

    def params(self, theta):
        p = dict(self.fixed)
        p.update(zip(self.free, theta))
        return p

    def absolute_w(self, p):
        R = p["R"]
        a = np.where(self.is1o, self.shapes["1O"].a(R) if "1O" in self.shapes else np.nan, self.shapes["F"].a(R))
        b = np.where(self.is1o, self.shapes["1O"].b(R) if "1O" in self.shapes else np.nan, self.shapes["F"].b(R))
        return (a + b * (self.logp - C.LOGP_PIVOT) - p["mu"] + p["gamma"] * (self.feh - self.feh_lmc)
                + np.where(self.is1o, p["d1O"], 0.0))

    #: Correct the photometric parallax for the scatter of the law. A star's
    #: magnitude scatters by eps about the law, so its photometric parallax is
    #: plx_true 10^(-0.2 eps), whose mean over eps is plx_true exp(+(0.2 ln10)^2 s^2 / 2):
    #: too large on average. Uncorrected, the fit pushes mu_LMC up by ~0.23 s^2 per
    #: star (+0.030 mag overall in our injection test, run.py section 7).
    lognormal_correction = True

    def terms(self, theta):
        p = self.params(theta)
        w = self.g - p["R"] * (self.bp - self.rp)
        pred = 10 ** (-0.2 * (w - self.absolute_w(p) - 10.0))
        smag2 = np.exp(2 * p["lnsig"]) + (p["sigR"] * self.E) ** 2
        if self.lognormal_correction:
            pred = pred * np.exp(-0.5 * LN10_5 ** 2 * smag2)
        var = (p["k"] * self.e_plx) ** 2 + (LN10_5 * pred) ** 2 * smag2
        r = self.plx - self.zp - p["delta"] - pred
        return r, var, pred, p

    def lnprob(self, theta):
        p = self.params(theta)
        if not (0.3 < p["k"] < 5 and -9 < p["lnsig"] < 0 and 0 <= p["sigR"] < 1 and 1.0 < p["R"] < 2.8):
            return -np.inf
        r, var, _, _ = self.terms(theta)
        m = self.mask
        return -0.5 * np.sum(r[m] ** 2 / var[m] + np.log(var[m]))

    def start(self):
        return np.array([self.fixed[n] for n in self.free], float)

    def fit(self, clip=4.0, iters=10):
        """MAP fit with iterative clipping of > ``clip``-sigma parallax outliers."""
        from scipy.optimize import minimize
        theta = self.start()
        for _ in range(iters):
            f = lambda t: -self.lnprob(t)
            res = minimize(f, theta, method="Nelder-Mead",
                           options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 60000, "maxfev": 60000})
            res = minimize(f, res.x, method="Powell", options={"xtol": 1e-9, "ftol": 1e-11})
            theta = res.x
            r, var, _, _ = self.terms(theta)
            new = np.abs(r / np.sqrt(var)) < clip
            if np.array_equal(new, self.mask):
                break
            self.mask = new
        self.theta = theta
        return theta

    def sample(self, nwalkers=40, nsteps=9000, burn=2500, seed=C.SEED):
        import emcee
        ndim = len(self.free)
        rng = np.random.default_rng(seed)
        scale = np.array([{"mu": 0.01, "gamma": 0.03, "delta": 0.002, "R": 0.01, "lnsig": 0.05, "sigR": 0.02,
                           "k": 0.03, "d1O": 0.01}[n] for n in self.free])
        p0 = self.theta + scale * rng.standard_normal((nwalkers, ndim))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, self.lnprob)
        sampler.random_state = np.random.RandomState(seed).get_state()
        sampler.run_mcmc(p0, nsteps, progress=False)
        self.chain = sampler.get_chain(discard=burn, flat=True)
        self.tau = sampler.get_autocorr_time(discard=burn, quiet=True)
        return self.chain

    def summary(self) -> dict:
        out = {"n": int(self.mask.sum()), "n_clipped": int((~self.mask).sum()),
               "n_F": int(np.sum(self.mask & ~self.is1o)), "n_1O": int(np.sum(self.mask & self.is1o)),
               "map": dict(zip(self.free, map(float, self.theta)))}
        if hasattr(self, "chain"):
            for j, name in enumerate(self.free):
                q = np.percentile(self.chain[:, j], [16, 50, 84])
                out[name] = {"median": float(q[1]), "lo": float(q[1] - q[0]), "hi": float(q[2] - q[1]),
                             "std": float(np.std(self.chain[:, j]))}
            out["tau_max"] = float(np.nanmax(self.tau))
        return out
