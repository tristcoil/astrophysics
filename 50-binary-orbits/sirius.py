"""Sirius A and B: our own orbit from 154 years of measurements, and the masses.

Data: Bond et al. 2017 (ApJ 840, 70), the critical compilation of every
published measure of Sirius B relative to Sirius A (2,367 rows, 1862-2016:
visual micrometers, photographic plates, CCDs, and their own HST astrometry),
already corrected to the J2000 equator and for the system's space motion.

Three fits, from least to most independent of the paper:

1. **Blind start** (``thiele_innes_search``): no prior knowledge. For each
   trial (P, e, T) the four Thiele-Innes constants follow by linear least
   squares; a grid plus a polish finds the global minimum.
2. **Replication** (``fit_bond_weights``): Bond's own adopted uncertainties
   and rejections, plain weighted least squares in separation and position
   angle. Should reproduce their Table 4.
3. **Our fit** (``fit_mcmc``): our own noise model, fitted together with the
   orbit by MCMC: per-method scatter (micrometer errors that shrink with
   aperture and improve over time; isotropic photographic and CCD errors; a
   free scale on the HST errors) and an outlier mixture in place of sigma
   clipping. Only Bond's documented rejections (a false image, superseded HST
   reductions) are dropped by hand.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
from scipy import optimize

import config as C
import literature as L
import orbits

DATA_FILE = C.LIT_DIR / "bond2017" / "sirius.dat"
PARAM_NAMES = ("P", "T", "e", "a", "i", "Omega", "omega")
METHODS = {"M": "visual micrometer", "P": "photographic", "C": "CCD / mid-infrared", "H": "HST"}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------


def _f(s: str) -> float:
    s = s.strip()
    return float(s) if s else np.nan


def load_measures(path=DATA_FILE) -> pd.DataFrame:
    """Read Bond et al. 2017's sirius.dat (byte layout from its ReadMe)."""
    rows = []
    for line in path.read_text(encoding="latin-1").splitlines():
        line = line.ljust(157)
        rows.append({
            "date": line[0:11].strip(), "by": _f(line[12:21]),
            "pa_meas": _f(line[23:30]), "sep_meas": _f(line[32:40]),
            "pa": _f(line[43:50]), "e_pa": _f(line[52:57]),
            "sep": _f(line[58:66]), "e_sep": _f(line[68:74]),
            "obs": line[76:98].strip(), "method": line[100], "tel": _f(line[102:105]),
            "notes": line[106:114].strip(), "bibcode": line[115:134].strip(), "refs": line[136:157].strip(),
        })
    df = pd.DataFrame(rows)
    df["year"] = orbits.besselian_to_year(df["by"])          # Julian epoch
    df["has_pa"] = df["pa"] > -90.0
    df["has_sep"] = df["sep"] > -90.0
    df.loc[~df["has_pa"], "pa"] = np.nan
    df.loc[~df["has_sep"], "sep"] = np.nan
    # Bond et al.'s own use: uncertainties of 0 mark rejected measures
    df["bond_used"] = ((df["e_pa"] > 0) | (df["e_sep"] > 0)) & (df["has_pa"] | df["has_sep"])
    note_codes = df["notes"].str.findall(r"[RCMNU]\d+")
    df["bond_r1"] = note_codes.apply(lambda c: "R1" in c)
    excl = note_codes.apply(lambda c: any(x in c for x in C.SIRIUS_EXCLUDE_NOTES))
    unresolved = df["method"].eq("M") & ~df["has_pa"] & ~df["has_sep"]
    df["ours_used"] = (df["has_pa"] | df["has_sep"]) & ~excl & ~unresolved
    df["group"] = df["method"]
    return df


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


def model_rho_theta(theta7, t):
    P, T, e, a, i, O, w = theta7
    return orbits.rho_theta(t, P, T, e, a, i, O, w)


def residuals_pa_sep(theta7, df):
    """(PA residual in deg, separation residual in arcsec, model separation) per row."""
    rho, th = model_rho_theta(theta7, df["year"].to_numpy())
    dpa = orbits.wrap_angle_deg(df["pa"].to_numpy() - th)
    dsep = df["sep"].to_numpy() - rho
    return dpa, dsep, rho


# ---------------------------------------------------------------------------
# 1. Blind Thiele-Innes search
# ---------------------------------------------------------------------------


def _ti_chi2(P, e, T, t, x, y, w):
    """Weighted linear solve for (A, F) and (B, G); returns chi2 and the constants."""
    X, Y = orbits.unit_orbit(t, P, T, e)
    M = np.stack([X, Y], axis=1) * w[:, None]
    lhs = M.T @ M
    try:
        AF = np.linalg.solve(lhs, M.T @ (x * w))
        BG = np.linalg.solve(lhs, M.T @ (y * w))
    except np.linalg.LinAlgError:
        return np.inf, None
    rx = (x - (AF[0] * X + AF[1] * Y)) * w
    ry = (y - (BG[0] * X + BG[1] * Y)) * w
    return float(rx @ rx + ry @ ry), (AF[0], BG[0], AF[1], BG[1])


def thiele_innes_search(df, P_grid=np.linspace(20.0, 120.0, 201), e_grid=np.linspace(0.0, 0.95, 39),
                        n_T=48):
    """Global search over (P, e, T) with the Thiele-Innes constants solved linearly.

    Uses the measures with both PA and separation from photographic plates,
    CCDs and HST (micrometer data are too noisy for the coarse grid), with
    isotropic weights: 1/0.1" for plates and CCDs, 1/0.005" for HST.
    """
    m = df["ours_used"] & df["has_pa"] & df["has_sep"] & df["method"].isin(["P", "C", "H"])
    d = df[m]
    t = d["year"].to_numpy()
    th = np.radians(d["pa"].to_numpy())
    x, y = d["sep"].to_numpy() * np.cos(th), d["sep"].to_numpy() * np.sin(th)
    w = np.where(d["method"].to_numpy() == "H", 1 / 0.005, 1 / 0.1)
    best = (np.inf, None)
    t0 = t.min()
    for P in P_grid:
        for e in e_grid:
            for T in t0 + P * np.arange(n_T) / n_T:
                chi2, _ = _ti_chi2(P, e, T, t, x, y, w)
                if chi2 < best[0]:
                    best = (chi2, (P, e, T))
    # polish (P, e, T) with Nelder-Mead
    P0, e0, T0 = best[1]
    res = optimize.minimize(lambda p: _ti_chi2(p[0], np.clip(p[1], 0, 0.99), p[2], t, x, y, w)[0],
                            [P0, e0, T0], method="Nelder-Mead",
                            options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 20000})
    P, e, T = res.x
    e = float(np.clip(e, 0, 0.99))
    chi2, (A, B, F, G) = _ti_chi2(P, e, T, t, x, y, w)
    a, i, O, om = orbits.campbell_from_thiele_innes(A, B, F, G)
    # periastron time folded near the data
    return {"theta": np.array([P, T, e, float(a), float(i), float(O), float(om)]), "chi2": chi2,
            "n_points": int(m.sum()), "grid_best": best[1]}


# ---------------------------------------------------------------------------
# 2. Replication with Bond et al.'s adopted weights
# ---------------------------------------------------------------------------


def _stack_residuals(theta7, d, sig_pa, sig_sep):
    dpa, dsep, _ = residuals_pa_sep(theta7, d)
    r = []
    mpa = np.isfinite(dpa) & (sig_pa > 0)
    msep = np.isfinite(dsep) & (sig_sep > 0)
    r.append(dpa[mpa] / sig_pa[mpa])
    r.append(dsep[msep] / sig_sep[msep])
    return np.concatenate(r)


def fit_bond_weights(df, start, *, mask=None):
    """Least squares in (PA, sep) with Bond et al.'s adopted uncertainties."""
    d = df[df["bond_used"] if mask is None else mask]
    sig_pa = d["e_pa"].to_numpy()
    sig_sep = d["e_sep"].to_numpy()
    res = optimize.least_squares(lambda p: _stack_residuals(p, d, sig_pa, sig_sep), start,
                                 x_scale="jac", method="lm", xtol=1e-14, ftol=1e-14, max_nfev=20000)
    J = res.jac
    cov = np.linalg.inv(J.T @ J)
    n, k = len(res.fun), len(start)
    chi2 = float(res.fun @ res.fun)
    return {"theta": res.x, "err": np.sqrt(np.diag(cov)), "cov": cov, "chi2": chi2, "n": n,
            "chi2_nu": chi2 / (n - k)}


# ---------------------------------------------------------------------------
# 3. Our noise model + MCMC
# ---------------------------------------------------------------------------

#: Noise hyperparameters (all but the betas and the time trend are logs):
#: ln s_H      scale applied to the published HST errors
#: ln sig_P    photographic plates, isotropic (arcsec)
#: ln sig_C    CCD and mid-infrared, isotropic (arcsec)
#: ln sig_Mr   micrometer separation error at D = 0.5 m, epoch 1900 (arcsec)
#: ln sig_Mt   micrometer position-angle error at D = 0.5 m, epoch 1900 (deg)
#: beta        aperture exponent, sigma ~ (D / 0.5 m)^-beta
#: gamma       change of ln sigma per century (micrometers)
#: logit f_M   outlier fraction, micrometer
#: logit f_P   outlier fraction, plates + CCD
#: bias_Msep   micrometer separation bias at rho = 4", growing as (4"/rho)^4: visual
#:             observers overestimate small separations in Sirius A's glare
#: bias_Mpa    constant offset of micrometer position angles (deg)
NOISE_NAMES = ("ln_sH", "ln_sigP", "ln_sigC", "ln_sigMr", "ln_sigMt", "beta", "gamma", "logit_fM", "logit_fP",
               "bias_Msep", "bias_Mpa")
#: Reference separation and steepness of the micrometer separation bias.
MSEP_BIAS_RHO, MSEP_BIAS_POWER = 4.0, 4.0
ALL_NAMES = ("P", "T", "sqrte_cosw", "sqrte_sinw", "a", "cosi", "Omega") + NOISE_NAMES


class SiriusLikelihood:
    """Likelihood of the orbit + noise model for a set of measures (vectorised)."""

    def __init__(self, df, use=None):
        d = df[df["ours_used"] if use is None else use].reset_index(drop=True)
        self.d = d
        self.t = d["year"].to_numpy()
        self.pa = d["pa"].to_numpy()
        self.sep = d["sep"].to_numpy()
        self.has_pa = d["has_pa"].to_numpy()
        self.has_sep = d["has_sep"].to_numpy()
        g = d["method"].to_numpy()
        self.isH, self.isP, self.isC, self.isM = (g == "H"), (g == "P"), (g == "C"), (g == "M")
        self.eH_pa = d["e_pa"].to_numpy()
        self.eH_sep = d["e_sep"].to_numpy()
        D = d["tel"].to_numpy()
        self.lnD = np.log(np.where(D > 0, D, C.MICROMETER_D_REF) / C.MICROMETER_D_REF)
        self.cent = (self.t - 1900.0) / 100.0
        self.present = tuple(m.any() for m in (self.isH, self.isP, self.isC, self.isM))

    @staticmethod
    def unpack(p):
        P, T, sc, ss, a, cosi, O = p[:7]
        e = sc * sc + ss * ss
        w = np.degrees(np.arctan2(ss, sc)) % 360.0
        i = np.degrees(np.arccos(cosi))
        return np.array([P, T, e, a, i, O, w]), p[7:]

    @staticmethod
    def pack(theta7, noise):
        P, T, e, a, i, O, w = theta7
        se = np.sqrt(e)
        return np.concatenate([[P, T, se * np.cos(np.radians(w)), se * np.sin(np.radians(w)), a,
                                np.cos(np.radians(i)), O], noise])

    def sigmas(self, noise, rho_model):
        ln_sH, ln_sP, ln_sC, ln_sMr, ln_sMt, beta, gamma = noise[:7]
        sig_pa = np.empty_like(self.t)
        sig_sep = np.empty_like(self.t)
        # HST: published errors times a common scale
        sH = np.exp(ln_sH)
        sig_pa[self.isH] = self.eH_pa[self.isH] * sH
        sig_sep[self.isH] = self.eH_sep[self.isH] * sH
        # plates and CCDs: isotropic, converted to an angle at the model separation
        for m, ln_s in ((self.isP, ln_sP), (self.isC, ln_sC)):
            s = np.exp(ln_s)
            sig_sep[m] = s
            sig_pa[m] = np.degrees(s / rho_model[m])
        # micrometers: aperture and epoch dependence
        m = self.isM
        scale = np.exp(-beta * self.lnD[m] + gamma * self.cent[m])
        sig_sep[m] = np.exp(ln_sMr) * scale
        sig_pa[m] = np.exp(ln_sMt) * scale
        return sig_pa, sig_sep

    def outlier_frac(self, noise):
        fM = 1.0 / (1.0 + np.exp(-noise[7]))
        fP = 1.0 / (1.0 + np.exp(-noise[8]))
        f = np.zeros_like(self.t)
        f[self.isM] = fM
        f[self.isP | self.isC] = fP
        return f

    @staticmethod
    def _mix_logpdf(r, s, f):
        """log[(1-f) N(r; 0, s) + f N(r; 0, W s)] (stable)."""
        W = C.OUTLIER_WIDTH
        z = r / s
        lg = -0.5 * z * z - np.log(s) - 0.5 * np.log(2 * np.pi)
        lo = -0.5 * (z / W) ** 2 - np.log(W * s) - 0.5 * np.log(2 * np.pi)
        with np.errstate(divide="ignore"):
            a = np.log1p(-f) + lg
            b = np.where(f > 0, np.log(f) + lo, -np.inf)
        return np.logaddexp(a, b)

    def bias(self, noise, rho_model):
        """Systematic offsets (PA in deg, separation in arcsec) expected for each measure."""
        b_sep = np.zeros_like(self.t)
        b_pa = np.zeros_like(self.t)
        m = self.isM
        if len(noise) > 9:
            b_sep[m] = noise[9] * (MSEP_BIAS_RHO / rho_model[m]) ** MSEP_BIAS_POWER
            b_pa[m] = noise[10]
        return b_pa, b_sep

    def components(self, p):
        theta7, noise = self.unpack(p)
        rho, th = orbits.rho_theta(self.t, *theta7)
        b_pa, b_sep = self.bias(noise, rho)
        dpa = orbits.wrap_angle_deg(self.pa - th - b_pa)
        dsep = self.sep - rho - b_sep
        sig_pa, sig_sep = self.sigmas(noise, rho)
        f = self.outlier_frac(noise)
        return dpa, dsep, sig_pa, sig_sep, f, rho

    def loglike(self, p):
        dpa, dsep, sig_pa, sig_sep, f, _ = self.components(p)
        lp = self._mix_logpdf(dpa[self.has_pa], sig_pa[self.has_pa], f[self.has_pa]).sum()
        ls = self._mix_logpdf(dsep[self.has_sep], sig_sep[self.has_sep], f[self.has_sep]).sum()
        return lp + ls

    def logprior(self, p):
        P, T, sc, ss, a, cosi, O = p[:7]
        e = sc * sc + ss * ss
        noise = p[7:]
        if not (1.0 < P < 500.0 and 0.0 <= e < 0.99 and 0.0 < a < 100.0 and -1.0 < cosi < 1.0
                and 0.0 <= O < 180.0):
            return -np.inf
        if not (-3.0 < (T - 2000.0) / P < 3.0):
            return -np.inf
        if np.any(np.abs(noise[:5]) > 12.0) or not (-3.0 < noise[5] < 5.0) or abs(noise[6]) > 5.0:
            return -np.inf
        if np.any(np.abs(noise[7:9]) > 12.0) or np.any(np.abs(noise[9:]) > 2.0):
            return -np.inf
        # flat in (sqrt(e) cos w, sqrt(e) sin w) is flat in e and w; flat in cos i is isotropic
        # orientation; the outlier fractions get a weak Beta(1, 9)-like pull toward small values
        fpen = 0.0
        for lf in noise[7:9]:
            f = 1.0 / (1.0 + np.exp(-lf))
            fpen += 8.0 * np.log1p(-f) + np.log(f) + np.log1p(-f)   # Beta(1,9) density x logit Jacobian
        return fpen

    def __call__(self, p):
        lp = self.logprior(p)
        if not np.isfinite(lp):
            return -np.inf
        ll = self.loglike(p)
        return lp + ll if np.isfinite(ll) else -np.inf

    def outlier_probability(self, p):
        """Posterior probability that each PA and each separation is an outlier."""
        dpa, dsep, sig_pa, sig_sep, f, _ = self.components(p)
        out = {}
        for name, r, s, m in (("pa", dpa, sig_pa, self.has_pa), ("sep", dsep, sig_sep, self.has_sep)):
            W = C.OUTLIER_WIDTH
            z = r / s
            lg = np.log1p(-f) - 0.5 * z * z - np.log(s)
            with np.errstate(divide="ignore"):
                lo = np.where(f > 0, np.log(f), -np.inf) - 0.5 * (z / W) ** 2 - np.log(W * s)
            with np.errstate(invalid="ignore"):
                q = np.exp(lo - np.logaddexp(lg, lo))
            out[name] = np.where(m, q, np.nan)
        return out


def person(obs_code: str) -> str:
    """Observer from Bond et al.'s code ('Burnham_1883a' -> 'Burnham')."""
    import re
    return re.sub(r"_\d{4}[a-z]?$", "", obs_code)


#: Hyperparameters of the hierarchical model (after the 7 orbit parameters):
HIER_NAMES = ("ln_sH", "ln_sigP", "ln_sigC", "ln_sigMr", "ln_sigMt", "beta", "gamma", "bias_Msep",
              "mean_Mpa", "ln_tauMpa", "ln_tauMsep", "ln_tauPpa", "ln_tauPsep")
HIER_SCALE = np.array([0.004, 0.006, 0.0005, 0.0005, 0.0025, 0.0005, 0.07,
                       0.05, 0.05, 0.2, 0.03, 0.03, 0.05, 0.05, 0.02, 0.05, 0.1, 0.1, 0.1, 0.1])


class HierLikelihood(SiriusLikelihood):
    """Gaussian likelihood with each observer's personal equation integrated out.

    Every visual observer (and every photographic programme) has his own
    systematic offsets in position angle and separation, drawn from a normal
    population: PA offsets ~ N(mean_Mpa, tau_Mpa^2) for micrometers and
    N(0, tau_Ppa^2) for plates; separation offsets ~ N(0, tau^2) likewise. The
    offsets are marginalised analytically (Sherman-Morrison for each observer's
    block), so 1,900 measures by 119 people count for what they are worth.
    Outliers must be removed beforehand (``screen``); the per-point noise model
    and the glare bias of small micrometer separations are as in the parent.
    """

    def __init__(self, df, use_pa, use_sep):
        super().__init__(df, use=use_pa | use_sep)
        d = self.d
        # masks aligned with self.d (the rows kept by the parent)
        idx = df.index[use_pa | use_sep]
        self.m_pa = use_pa.loc[idx].to_numpy() & self.has_pa
        self.m_sep = use_sep.loc[idx].to_numpy() & self.has_sep
        persons = d["obs"].map(person).to_numpy()
        gid = np.full(len(d), -1)
        names = []
        for pname, meth in sorted(set(zip(persons, d["method"]))):
            sel = (persons == pname) & (d["method"].to_numpy() == meth)
            if meth in ("M", "P"):
                gid[sel] = len(names)
                names.append(f"{meth}:{pname}")
        self.gid = gid
        self.group_names = names
        self.group_is_M = np.array([n.startswith("M:") for n in names])
        self.ngroups = len(names)

    def unpack(self, p):
        theta7, _ = SiriusLikelihood.unpack(p)
        return theta7, p[7:]

    def pack(self, theta7, hyper):
        return SiriusLikelihood.pack(theta7, hyper)

    def sigmas(self, noise, rho_model):
        return SiriusLikelihood.sigmas(self, noise[:7], rho_model)

    def bias(self, noise, rho_model):
        b_sep = np.zeros_like(self.t)
        m = self.isM
        b_sep[m] = noise[7] * (MSEP_BIAS_RHO / rho_model[m]) ** MSEP_BIAS_POWER
        b_pa = np.where(m, noise[8], 0.0)
        return b_pa, b_sep

    def _block_loglike(self, r, sig, mask, tau_per_group):
        """Sum of log N over all points, with per-group random offsets marginalised."""
        w = np.where(mask, 1.0 / sig**2, 0.0)
        r = np.where(mask, r, 0.0)
        ll = -0.5 * np.sum(w * r * r) - np.sum(np.where(mask, np.log(sig), 0.0)) \
            - 0.5 * np.log(2 * np.pi) * mask.sum()
        g = self.gid
        ing = (g >= 0) & mask
        S1 = np.bincount(g[ing], weights=w[ing], minlength=self.ngroups)
        Sr = np.bincount(g[ing], weights=(w * r)[ing], minlength=self.ngroups)
        t2 = tau_per_group**2
        denom = 1.0 + t2 * S1
        ll += 0.5 * np.sum(t2 * Sr * Sr / denom) - 0.5 * np.sum(np.log(denom))
        return ll

    def loglike(self, p):
        theta7, h = self.unpack(p)
        rho, th = orbits.rho_theta(self.t, *theta7)
        b_pa, b_sep = self.bias(h, rho)
        dpa = orbits.wrap_angle_deg(self.pa - th - b_pa)
        dsep = self.sep - rho - b_sep
        sig_pa, sig_sep = self.sigmas(h, rho)
        tau_pa = np.where(self.group_is_M, np.exp(h[9]), np.exp(h[11]))
        tau_sep = np.where(self.group_is_M, np.exp(h[10]), np.exp(h[12]))
        return (self._block_loglike(np.nan_to_num(dpa), sig_pa, self.m_pa, tau_pa)
                + self._block_loglike(np.nan_to_num(dsep), sig_sep, self.m_sep, tau_sep))

    def logprior(self, p):
        P, T, sc, ss, a, cosi, O = p[:7]
        e = sc * sc + ss * ss
        h = p[7:]
        if not (1.0 < P < 500.0 and 0.0 <= e < 0.99 and 0.0 < a < 100.0 and -1.0 < cosi < 1.0
                and 0.0 <= O < 180.0 and -3.0 < (T - 2000.0) / P < 3.0):
            return -np.inf
        if np.any(np.abs(h[:5]) > 12.0) or not (-3.0 < h[5] < 5.0) or abs(h[6]) > 5.0:
            return -np.inf
        if abs(h[7]) > 2.0 or abs(h[8]) > 5.0 or np.any(h[9:] < -9.0) or np.any(h[9:] > 2.5):
            return -np.inf
        return 0.0

    def observer_offsets(self, p):
        """Posterior mean personal equation of each group, given the parameters p."""
        theta7, h = self.unpack(p)
        rho, th = orbits.rho_theta(self.t, *theta7)
        b_pa, b_sep = self.bias(h, rho)
        dpa = np.nan_to_num(orbits.wrap_angle_deg(self.pa - th - b_pa))
        dsep = np.nan_to_num(self.sep - rho - b_sep)
        sig_pa, sig_sep = self.sigmas(h, rho)
        out = {}
        for lab, r, sig, mask, tau in (
                ("pa", dpa, sig_pa, self.m_pa, np.where(self.group_is_M, np.exp(h[9]), np.exp(h[11]))),
                ("sep", dsep, sig_sep, self.m_sep, np.where(self.group_is_M, np.exp(h[10]), np.exp(h[12])))):
            w = np.where(mask, 1 / sig**2, 0.0)
            ing = (self.gid >= 0) & mask
            S1 = np.bincount(self.gid[ing], weights=w[ing], minlength=self.ngroups)
            Sr = np.bincount(self.gid[ing], weights=(w * r)[ing], minlength=self.ngroups)
            n = np.bincount(self.gid[ing], minlength=self.ngroups)
            mean = tau**2 * Sr / (1 + tau**2 * S1)
            sd = np.sqrt(tau**2 / (1 + tau**2 * S1))
            out[lab] = pd.DataFrame({"group": self.group_names, "n": n, "offset": mean, "err": sd})
        return out


def screen(df, p_mix, like_mix: SiriusLikelihood, threshold=0.5):
    """Outlier flags from the mixture model: keep a PA (separation) if P(outlier) < threshold."""
    q = like_mix.outlier_probability(p_mix)
    use_pa = pd.Series(False, index=df.index)
    use_sep = pd.Series(False, index=df.index)
    idx = df.index[df["ours_used"]]
    use_pa.loc[idx] = like_mix.has_pa & (np.nan_to_num(q["pa"], nan=1.0) < threshold)
    use_sep.loc[idx] = like_mix.has_sep & (np.nan_to_num(q["sep"], nan=1.0) < threshold)
    return use_pa, use_sep, q


def default_hyper(p_mix_noise=None):
    """Starting hyperparameters for the hierarchical model."""
    h = np.array([np.log(1.35), np.log(0.06), np.log(0.15), np.log(0.30), np.log(1.4), 0.3, -0.5,
                  0.18, 0.3, np.log(0.5), np.log(0.05), np.log(0.1), np.log(0.01)])
    if p_mix_noise is not None:
        h[:7] = p_mix_noise[:7]
        h[7] = p_mix_noise[9]
    return h


def default_noise():
    """Starting noise hyperparameters (roughly Bond et al.'s adopted values)."""
    return np.array([np.log(1.3), np.log(0.08), np.log(0.08), np.log(0.33), np.log(1.7), 0.3, 0.0,
                     np.log(0.05 / 0.95), np.log(0.02 / 0.98), 0.0, 0.0])


#: Typical posterior widths of the packed parameters, used to put the optimiser
#: on a well-conditioned scale (orbit: Bond et al.'s errors; noise: ~5%).
PACKED_SCALE = np.array([0.004, 0.006, 0.0005, 0.0005, 0.0025, 0.0005, 0.07,
                         0.05, 0.05, 0.2, 0.03, 0.03, 0.05, 0.05, 0.2, 0.3, 0.02, 0.05])


def _minimize_scaled(fun, p0, scale, idx, method):
    """Minimise fun over the parameters idx of p0, in units of scale."""
    def f(x):
        p = p0.copy()
        p[idx] = p0[idx] + scale[idx] * x
        v = fun(p)
        return v if np.isfinite(v) else 1e300
    opts = ({"maxiter": 40000, "maxfev": 40000, "xatol": 1e-6, "fatol": 1e-8} if method == "Nelder-Mead"
            else {"maxiter": 40000, "xtol": 1e-6, "ftol": 1e-12})
    res = optimize.minimize(f, np.zeros(len(idx)), method=method, options=opts)
    p = p0.copy()
    p[idx] = p0[idx] + scale[idx] * res.x
    return p


def fit_map(like: SiriusLikelihood, theta7_start, noise_start=None, scale=PACKED_SCALE):
    """Maximum a posteriori point: alternate noise-only and orbit-only steps, then polish all.

    A single 16-dimensional simplex wanders off (it can "explain" the HST data by
    inflating their errors instead of moving the orbit), so the search runs on
    scaled parameters and in blocks.
    """
    p = like.pack(theta7_start, default_noise() if noise_start is None else noise_start)
    nlp = lambda q: -like(q)
    orbit_idx, noise_idx = np.arange(7), np.arange(7, len(p))
    for _ in range(3):
        p = _minimize_scaled(nlp, p, scale, noise_idx, "Nelder-Mead")
        p = _minimize_scaled(nlp, p, scale, orbit_idx, "Nelder-Mead")
    for method in ("Powell", "Nelder-Mead"):
        p = _minimize_scaled(nlp, p, scale, np.arange(len(p)), method)
    return p, like(p)


def fit_mcmc(like: SiriusLikelihood, p_map, *, nwalkers=C.MCMC_WALKERS, nsteps=C.MCMC_STEPS,
             burn=C.MCMC_BURN, thin=C.MCMC_THIN, seed=C.SEED, cache=None, label="", scale=None):
    """emcee posterior around the MAP point. Returns flat samples (unpacked to elements)."""
    import emcee

    if cache is not None and cache.exists():
        z = np.load(cache)
        return {"chain": z["chain"], "logp": z["logp"], "acc": float(z["acc"]), "tau": z["tau"]}
    rng = np.random.default_rng(seed)
    ndim = len(p_map)
    scale = (HIER_SCALE if len(p_map) == len(HIER_SCALE) else PACKED_SCALE) * 0.01 if scale is None else scale
    p0 = p_map + scale * rng.standard_normal((nwalkers, ndim))
    for k in range(nwalkers):
        while not np.isfinite(like(p0[k])):
            p0[k] = p_map + scale * rng.standard_normal(ndim)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, like)
    sampler.random_state = np.random.RandomState(seed).get_state()   # see the emcee seeding trap
    t0 = time.time()
    sampler.run_mcmc(p0, nsteps, progress=False)
    print(f"    MCMC {label}: {nwalkers}x{nsteps} in {time.time() - t0:.0f} s, "
          f"acceptance {np.mean(sampler.acceptance_fraction):.2f}")
    try:
        tau = sampler.get_autocorr_time(discard=burn, quiet=True)
    except Exception:  # noqa: BLE001
        tau = np.full(ndim, np.nan)
    chain = sampler.get_chain(discard=burn, thin=thin, flat=True)
    logp = sampler.get_log_prob(discard=burn, thin=thin, flat=True)
    out = {"chain": chain, "logp": logp, "acc": float(np.mean(sampler.acceptance_fraction)), "tau": tau}
    if cache is not None:
        np.savez_compressed(cache, **out)
    return out


def elements_from_chain(like: SiriusLikelihood, chain):
    """Flat chain of packed vectors -> (n, 7) array of Campbell elements."""
    return np.array([like.unpack(p)[0] for p in chain])


def summarize(samples, names=PARAM_NAMES):
    q = np.percentile(samples, [15.865, 50, 84.135], axis=0)
    return {n: {"median": float(q[1, k]), "lo": float(q[1, k] - q[0, k]), "hi": float(q[2, k] - q[1, k]),
                "mean": float(np.mean(samples[:, k])), "std": float(np.std(samples[:, k]))}
            for k, n in enumerate(names)}


# ---------------------------------------------------------------------------
# Masses
# ---------------------------------------------------------------------------


def masses_from_samples(P, a, rng, *, plx=L.SIRIUS_PLX_ADOPTED, aA=L.BOND17_AA):
    """Monte Carlo masses: orbit samples x parallax x Sirius A's orbit size.

    plx, aA: dicts with "value" and "err" in arcsec.
    """
    n = len(P)
    pi = rng.normal(plx["value"], plx["err"], n)
    a_A = rng.normal(aA["value"], aA["err"], n)
    Mtot = orbits.total_mass(a, pi, P)
    fB = a_A / a
    return {"Mtot": Mtot, "MA": Mtot * (1 - fB), "MB": Mtot * fB, "fB": fB, "plx": pi}


def error_budget(theta, cov_Pa, plx, aA):
    """Linear error budget for M_A and M_B (as Bond et al.'s Table 7)."""
    P, a = theta
    Mtot = orbits.total_mass(a, plx["value"], P)
    fB = aA["value"] / a
    MA, MB = Mtot * (1 - fB), Mtot * fB
    sP, sa = np.sqrt(np.diag(cov_Pa))
    rows = {
        "parallax": (3 * plx["err"] / plx["value"] * MA, 3 * plx["err"] / plx["value"] * MB),
        "a (relative orbit)": (abs(MA * 3 * sa / a + Mtot * aA["value"] / a**2 * sa),
                               abs(MB * 3 * sa / a - Mtot * aA["value"] / a**2 * sa)),
        "a_A (Sirius A's orbit)": (Mtot * aA["err"] / a, Mtot * aA["err"] / a),
        "period": (2 * sP / P * MA, 2 * sP / P * MB),
    }
    total = tuple(np.sqrt(sum(r[k] ** 2 for r in rows.values())) for k in (0, 1))
    return {"MA": MA, "MB": MB, "rows": rows, "total": total}


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def fold_T(theta7, ref=2000.0):
    """Move the periastron epoch to the passage closest to ``ref``."""
    th = np.array(theta7, float)
    th[1] = th[1] + th[0] * np.round((ref - th[1]) / th[0])
    return th


def laplace(like, p, idx, scale, n_draw=4000, seed=C.SEED):
    """Gaussian (Laplace) approximation around a MAP point over the parameters idx.

    Returns draws of the 7 Campbell elements. The Hessian of -log p is taken by
    central differences in units of ``scale``.
    """
    k = len(idx)
    h = 0.5
    f0 = -like(p)
    H = np.zeros((k, k))
    def f(dx):
        q = p.copy()
        q[idx] = p[idx] + scale[idx] * dx
        return -like(q)
    E = np.eye(k) * h
    for a in range(k):
        for b in range(a, k):
            if a == b:
                H[a, a] = (f(E[a]) - 2 * f0 + f(-E[a])) / h**2
            else:
                H[a, b] = H[b, a] = (f(E[a] + E[b]) - f(E[a] - E[b]) - f(-E[a] + E[b]) + f(-E[a] - E[b])) / (4 * h * h)
    cov = np.linalg.inv(H)
    rng = np.random.default_rng(seed)
    dx = rng.multivariate_normal(np.zeros(k), cov, n_draw)
    draws = np.repeat(p[None, :], n_draw, axis=0)
    draws[:, idx] = p[idx] + dx * scale[idx]
    return elements_from_chain(like, draws)


def _cached_map(path, fn):
    if path.exists():
        return np.load(path)
    p = fn()
    np.save(path, p)
    return p


def run(verbose=True):
    """Every Sirius fit. Expensive steps are cached in data/cache/."""
    say = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
    df = load_measures()
    out = {"n_rows": len(df), "n_by_method": df["method"].value_counts().to_dict(),
           "n_ours_used": int(df["ours_used"].sum()), "n_bond_used": int(df["bond_used"].sum()),
           "years": (float(df["year"].min()), float(df["year"].max())), "df": df}

    say("  1. Thiele-Innes blind search ...")
    cache = C.CACHE_DIR / "sirius_ti_search.json"
    if cache.exists():
        ti = json.loads(cache.read_text())
    else:
        ti = thiele_innes_search(df)
        ti = {**ti, "theta": ti["theta"].tolist()}
        cache.write_text(json.dumps(ti, default=float))
    ti["theta"] = fold_T(ti["theta"])
    out["ti"] = ti

    say("  2. replication with Bond et al.'s adopted weights ...")
    rep = fit_bond_weights(df, ti["theta"])
    rep["theta"] = fold_T(rep["theta"])
    out["replication"] = rep

    say("  3. robust screening: outlier-mixture model ...")
    lm = SiriusLikelihood(df)
    p_mix = _cached_map(C.CACHE_DIR / "sirius_mixture_map.npy", lambda: fit_map(lm, rep["theta"])[0])
    use_pa, use_sep, q = screen(df, p_mix, lm)
    out["mixture"] = {"p": p_mix, "theta": lm.unpack(p_mix)[0], "noise": dict(zip(NOISE_NAMES, lm.unpack(p_mix)[1])),
                      "q": q, "like": lm}
    flagged = df["ours_used"] & ~((use_pa | ~df["has_pa"]) & (use_sep | ~df["has_sep"]))
    out["screen"] = {"use_pa": use_pa, "use_sep": use_sep, "flagged": flagged,
                     "n_flagged": int(flagged.sum()), "n_flagged_bond_r1": int((flagged & df["bond_r1"]).sum()),
                     "n_bond_r1": int(df["bond_r1"].sum())}

    say("  4. hierarchical model (each observer's personal equation marginalised) ...")
    hl = HierLikelihood(df, use_pa, use_sep)
    h0 = default_hyper(lm.unpack(p_mix)[1])
    p_h = _cached_map(C.CACHE_DIR / "sirius_hier_map.npy",
                      lambda: fit_map(hl, rep["theta"], h0, scale=HIER_SCALE)[0])
    mc = fit_mcmc(hl, p_h, cache=C.CACHE_DIR / "sirius_hier_mcmc.npz", label="hierarchical, all data")
    el = elements_from_chain(hl, mc["chain"])
    out["hier"] = {"p_map": p_h, "theta_map": hl.unpack(p_h)[0], "elements": el, "chain": mc["chain"],
                   "acc": mc["acc"], "tau": mc["tau"], "summary": summarize(el),
                   "hyper_summary": summarize(mc["chain"][:, 7:], HIER_NAMES), "like": hl,
                   "offsets": hl.observer_offsets(np.median(mc["chain"], axis=0))}

    say("  5. alternative data choices (systematics) ...")
    alts = {}
    nonM = ~df["method"].eq("M")
    for lab, mask, hidx in (("plates + CCD + HST", nonM, np.r_[7, 8, 9, 18, 19]),
                            ("HST only", df["method"].eq("H"), np.r_[7])):
        hla = HierLikelihood(df, use_pa & mask, use_sep & mask)
        idx = np.r_[0:7, hidx]
        def _fit(hla=hla, idx=idx):
            q0 = hla.pack(rep["theta"], h0)
            for _ in range(4):
                q0 = _minimize_scaled(lambda z: -hla(z), q0, HIER_SCALE, idx, "Nelder-Mead")
            return _minimize_scaled(lambda z: -hla(z), q0, HIER_SCALE, idx, "Powell")
        tag = lab.replace(" ", "").replace("+", "_")
        pa_ = _cached_map(C.CACHE_DIR / f"sirius_alt_{tag}.npy", _fit)
        draws = laplace(hla, pa_, np.arange(7), HIER_SCALE)
        alts[lab] = {"theta": hla.unpack(pa_)[0], "summary": summarize(draws), "draws": draws,
                     "n_pa": int(hla.m_pa.sum()), "n_sep": int(hla.m_sep.sum())}
    out["alternatives"] = alts
    return out
