"""Fit the precessing-orbit model to a century of plates and to TESS.

    python switch_fit.py profiles "V907 Sco"   # TESS chi-square profiles (cached)
    python switch_fit.py fit "V907 Sco"        # nested sampling (resumable)

Likelihood
----------
* DASCH: every cleaned plate measurement. Model magnitude m0 - 2.5 log F(t),
  with F from the eclipse model (B band); a two-component mixture: Gaussian
  with the plate error plus a fitted extra scatter, and a broad outlier
  population (DASCH has many: blends, plate defects; ~10-20% for some stars).
* TESS: each sector's chi-square profile of the impact parameter x (scaled so
  its best fit has chi2/dof = 1), evaluated at the model's x for that epoch.
* TESS times of primary minimum (sectors with detected eclipses), with a
  5-minute floor for the light-travel-time and dynamical timing wobbles the
  model leaves out.
* The outer orbit's spectroscopic amplitude, which links the third star's mass
  to the outer orbit's tilt at the epoch of the RVs.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.special import logsumexp

import config as C
import switching as S

CACHE = C.CACHE / "switching"
CACHE.mkdir(parents=True, exist_ok=True)
DEG = np.pi / 180

#: third light in the plates' blue band (third star; for V907 Sco also the
#: 9" companion, which blends with it on the plates)
L3_B = {"SS Lac": 0.01, "HS Hya": 0.01, "V907 Sco": 0.04}
#: timing floor for TESS minima (days): unmodelled light-time + dynamical ETVs
TIMING_FLOOR = 5.0 / 1440
#: systematic floor on the TESS tilts for the secular fit (deg): its outer-
#: period terms are missing (largest for compact, eccentric V907 Sco)
SIGMA_I_SYS_SECULAR = {"SS Lac": 0.2, "HS Hya": 0.2, "V907 Sco": 0.5}


def _slug(name):
    return name.lower().replace(" ", "_").replace("-", "_")


# ---------------------------------------------------------------------------
# Starting points: the published ephemerides, re-anchored on TESS
# ---------------------------------------------------------------------------

def guess(name: str) -> dict:
    sp = S.SPECS[name]
    if name == "HS Hya":
        n = np.round((2458555.6 - sp.T0) / sp.P)
        T_ref = sp.T0 + (n + 0.12) * sp.P + 0.0038      # TESS S9 minimum (this work)
        th = dict(J=55, I=80, h0=0, m3=0.56)
    elif name == "V907 Sco":
        T_ref = 2459362.3272 + 0.0053                    # Zasche et al. epoch, S39 offset
        th = dict(J=26, I=71, h0=0, m3=1.06)
    else:   # SS Lac: the plates put the eclipses 0.21 and 0.79 in phase after periastron;
        # its clock is anchored in 1921, when it eclipsed
        n = np.round((sp.t_ref - sp.T0) / sp.P)
        T_ref = sp.T0 + (n + 0.21) * sp.P
        th = dict(J=27, I=84, h0=0, m3=0.8, r_sum=0.144, k=1.0, g0=sp.omega, gdot=0.4)
    th.update(T_ref=T_ref, P=sp.P, m0=9.0, sig_jit=0.1, f_out=0.08, mu_out=0.0, sig_out=0.8)
    return th


def theta_array(model: S.Model, d: dict) -> np.ndarray:
    return np.array([d[n] for n in model.names], float)


# ---------------------------------------------------------------------------
# TESS profiles (cached)
# ---------------------------------------------------------------------------

def profiles(name: str, force=False) -> list[dict]:
    path = CACHE / f"profiles_{_slug(name)}.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    m = S.Model(name)
    th = theta_array(m, guess(name))
    out = []
    for lc in S.tess_lightcurves(name):
        t0 = time.time()
        sf = S.SectorFit(name, lc, th)
        pr = sf.profile()
        if pr is None:
            continue
        pr["seconds"] = time.time() - t0
        out.append(pr)
        print(f"{name} S{pr['sector']}: x_best {pr['x_best']:.4f} detected {pr['detected']} "
              f"scale {pr['scale']:.2f} ({pr['seconds']:.0f}s)", flush=True)
        path.write_text(json.dumps(out))
    path.write_text(json.dumps(out))
    return out


# ---------------------------------------------------------------------------
# Likelihood and priors
# ---------------------------------------------------------------------------

#: systematic floor on the Hipparcos/ASAS tilts (passband and third-light
#: assumptions for ground photometry)
SIGMA_I_SYS_GROUND = 0.3


class GlobalFit:
    def __init__(self, name: str, use_ground: bool = False, J_range=(0.0, 180.0)):
        self.name = name
        self.use_ground = use_ground
        self.model = S.Model(name)
        self.sp = self.model.spec
        d = S.dasch_data(name)
        self.t = d["date_jd"].to_numpy()
        self.mag = d["mag"].to_numpy()
        self.err = d["err"].to_numpy()
        self.prof = profiles(name)
        self.set_systematics(SIGMA_I_SYS_SECULAR[name])
        if use_ground:
            gp = json.loads((CACHE / f"ground_profiles_{_slug(name)}.json").read_text())
            keep = self.prof
            self.prof = gp
            self.set_systematics(max(SIGMA_I_SYS_GROUND, SIGMA_I_SYS_SECULAR[name]))
            for pr in gp:
                pr["detected"] = False          # ground timings are not used
            self.prof = keep + gp
        g = guess(name)
        P = self.sp.P
        # prior ranges
        self.bounds = {
            "J": tuple(J_range), "I": (0.0, 90.0), "h0": (0.0, 360.0), "m3": self.sp.m3_range,
            "T_ref": (g["T_ref"] - 0.06 * P, g["T_ref"] + 0.06 * P),
            "P": (P * (1 - 5e-5), P * (1 + 5e-5)),
            "m0": (float(np.median(self.mag)) - 0.5, float(np.median(self.mag)) + 0.5),
            "sig_jit": (0.0, 0.4), "f_out": (0.0, 0.4), "mu_out": (-1.5, 1.5), "sig_out": (0.2, 3.0),
            "r_sum": (0.08, 0.25), "k": (0.7, 1.4), "g0": (0.0, 360.0), "gdot": (-2.0, 2.0),
        }
        if name == "SS Lac":
            # the plates' eclipses fix the 1921 phase only loosely: allow +-0.1 P
            self.bounds["T_ref"] = (g["T_ref"] - 0.1 * P, g["T_ref"] + 0.1 * P)

    def set_systematics(self, sigma_i_deg: float):
        """Smooth each TESS profile with a Gaussian penalty on a systematic
        offset (infimal convolution): c_eff(x) = min_x' [c(x') + (x-x')^2/s^2],
        with s the x-equivalent of sigma_i at that sector's tilt. This stands
        for what a model leaves out (the secular model has no outer-period
        terms; those move V907 Sco's tilt by tenths of a degree)."""
        rs0 = None
        for pr in self.prof:
            x = np.array(pr["x"])
            c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
            xf = np.linspace(0.0, x.max(), 1500)
            cf = np.interp(xf, x, c)
            if sigma_i_deg > 0:
                if rs0 is None:
                    g = self.model.geometry(self.model.unpack(theta_array(self.model, guess(self.name))))
                    rs0 = g["r1"] + g["r2"]
                i_best = np.arccos(min(1.0, pr["x_best"] * rs0))
                sx = max(1e-4, np.sin(i_best) * np.radians(sigma_i_deg) / rs0)
                ce = np.min(cf[None, :] + (xf[:, None] - xf[None, :]) ** 2 / sx**2, axis=1)
            else:
                ce = cf
            pr["_x"], pr["_c"] = xf, ce - ce.min()

    # -- priors ----------------------------------------------------------------
    def prior_transform(self, u):
        x = np.empty_like(u)
        for k, n in enumerate(self.model.names):
            lo, hi = self.bounds[n]
            if n in ("J",):                       # isotropic: cos uniform
                c_lo, c_hi = np.cos(np.radians(lo)), np.cos(np.radians(hi))
                x[k] = np.degrees(np.arccos(c_lo - u[k] * (c_lo - c_hi)))
            elif n == "I":                        # isotropic on [0, 90]
                x[k] = np.degrees(np.arccos(1 - u[k]))
            else:
                x[k] = lo + (hi - lo) * u[k]
        return x

    # -- model pieces ------------------------------------------------------------
    def conjunction_time(self, g, t_guess):
        """Time of primary conjunction nearest t_guess (two fixed-point steps)."""
        sp = self.sp
        t = t_guess
        for _ in range(3):
            h = g["h0"] + g["nu"] * (t - sp.t_ref)
            gg = g["g0"] + g["gdot"] * (t - sp.t_ref)
            lc = S._conj_lambda(h, g["eps_in"], g["I"], gg, sp.e)
            lam = g["lam0"] + g["n_lam"] * (t - sp.t_ref)
            d = (lc - lam + np.pi) % (2 * np.pi) - np.pi
            t = t + d / g["n_lam"]
        return t

    def outer_rv_term(self, g, p):
        rv = self.sp.outer_rv
        if not rv:
            return 0.0
        h = g["h0"] + g["nu"] * (rv["epoch"] - self.sp.t_ref)
        i_out = S.outer_inclination(g["I"], g["eps_out"], h)
        M = g["m1"] + g["m2"] + g["m3"]
        frac = g["m3"] / M
        if rv["kind"] == "a2proj":
            model = g["a_out"] / C.R_SUN_AU * frac * np.sin(i_out)
        else:
            v_orb = 2 * np.pi * g["a_out"] * C.AU / 1e3 / (self.sp.P_out * C.DAY) / np.sqrt(1 - self.sp.e_out**2)
            model = v_orb * frac * np.sin(i_out)
        return -0.5 * ((model - rv["value"]) / rv["err"]) ** 2

    def lnlike(self, theta, parts=False):
        p = self.model.unpack(theta)
        try:
            g = self.model.geometry(p)
        except (FloatingPointError, ValueError):
            return -1e300
        if not np.isfinite(g["nu"]) or g["r1"] <= 0 or g["r2"] <= 0:
            return -1e300
        # plates
        F = self.model.flux(g, self.t, band="B", L3=L3_B[self.name])
        mod = p["m0"] - 2.5 * np.log10(F)
        s2 = self.err**2 + p["sig_jit"] ** 2
        r = self.mag - mod
        lp_in = -0.5 * r * r / s2 - 0.5 * np.log(2 * np.pi * s2)
        so2 = p["sig_out"] ** 2 + self.err**2
        lp_out = -0.5 * (r - p["mu_out"]) ** 2 / so2 - 0.5 * np.log(2 * np.pi * so2)
        f = p["f_out"]
        ll_d = float(np.sum(np.logaddexp(np.log1p(-f) + lp_in, np.log(f + 1e-300) + lp_out)))
        # TESS
        ll_t = 0.0
        ll_time = 0.0
        for pr in self.prof:
            x = self.model.impact(g, pr["jd_mid"])[0]
            ll_t += -0.5 * float(np.interp(x, pr["_x"], pr["_c"]))
            if pr["detected"] and pr["t_pri_obs"] is not None:
                tm = self.conjunction_time(g, pr["t_pri_obs"])
                ll_time += -0.5 * ((tm - pr["t_pri_obs"]) / TIMING_FLOOR) ** 2
        ll_rv = self.outer_rv_term(g, p)
        tot = ll_d + ll_t + ll_time + ll_rv
        if parts:
            return dict(dasch=ll_d, tess=ll_t, timing=ll_time, rv=ll_rv, total=tot)
        return tot if np.isfinite(tot) else -1e300


# ---------------------------------------------------------------------------
# Nested sampling (resumable)
# ---------------------------------------------------------------------------

def _tag_options(tag: str) -> dict:
    """Run variants by tag: 'all' adds Hipparcos/ASAS-3; 'retro'/'pro' limit
    the mutual inclination to retrograde/prograde, 'quietpro'/'quietretro' to
    outside the Kozai-Lidov window (J < 39.2 or > 140.8 deg); 'rerun' uses
    another seed and twice the live points (a check that the sampler found
    every mode)."""
    parts = set(tag.split("_")) if tag else set()
    J = (90.0, 180.0) if "retro" in parts else (0.0, 90.0) if "pro" in parts else (0.0, 180.0)
    # outside the quadrupole Kozai-Lidov window (cos^2 J > 3/5): geometries in
    # which the inner eccentricity cannot be pumped (the secular model holds)
    if "quietpro" in parts:
        J = (0.0, 39.23)
    if "quietretro" in parts:
        J = (140.77, 180.0)
    return dict(use_ground="all" in parts, J_range=J, rerun="rerun" in parts)


def run(name: str, nlive: int = 500, dlogz: float = 0.2, seed: int = C.RANDOM_SEED, tag: str = "",
        nproc: int = 1):
    """tag '' : plates + TESS; tag 'all': plus Hipparcos and ASAS-3 (see _tag_options).
    nproc > 1 evaluates likelihoods in a process pool (also when resuming a
    checkpoint written by a serial run)."""
    import dynesty
    opt = _tag_options(tag)
    if opt["rerun"]:
        nlive, seed = 2 * nlive, seed + 1000
    fit = GlobalFit(name, use_ground=opt["use_ground"], J_range=opt["J_range"])
    sfx = f"_{tag}" if tag else ""
    ckpt = CACHE / f"dynesty_{_slug(name)}{sfx}.save"
    out = CACHE / f"posterior_{_slug(name)}{sfx}.npz"
    if out.exists():
        return dict(np.load(out, allow_pickle=True))
    rng = np.random.default_rng(seed)
    pool = None
    if nproc > 1:
        import multiprocessing as mp
        pool = mp.get_context("fork").Pool(nproc)
    try:
        if ckpt.exists():
            sampler = dynesty.NestedSampler.restore(str(ckpt), pool=pool)
            sampler.run_nested(resume=True, dlogz=dlogz, checkpoint_file=str(ckpt), print_progress=False)
        else:
            sampler = dynesty.NestedSampler(fit.lnlike, fit.prior_transform, fit.model.ndim, nlive=nlive,
                                            sample="rslice", bound="multi", rstate=rng, pool=pool,
                                            queue_size=nproc if pool else None)
            sampler.run_nested(dlogz=dlogz, checkpoint_file=str(ckpt), print_progress=False)
    finally:
        if pool is not None:
            pool.close()
    res = sampler.results
    w = np.exp(res.logwt - res.logz[-1])
    from dynesty.utils import resample_equal
    post = resample_equal(res.samples, w / w.sum(), rstate=np.random.default_rng(seed))
    np.savez(out, samples=res.samples, logwt=res.logwt, logl=res.logl, logz=res.logz[-1],
             logzerr=res.logzerr[-1], posterior=post, names=np.array(fit.model.names),
             ncall=int(np.sum(res.ncall)))
    return dict(np.load(out, allow_pickle=True))


if __name__ == "__main__":
    what, name = sys.argv[1], sys.argv[2]
    if what == "profiles":
        profiles(name, force=len(sys.argv) > 3)
    elif what == "fit":
        t0 = time.time()
        args = [a for a in sys.argv[3:] if not a.startswith("--nproc")]
        nproc = int(next((a.split("=")[1] for a in sys.argv[3:] if a.startswith("--nproc=")), 1))
        r = run(name, tag=args[0] if args else "", nproc=nproc)
        print(name, "logZ", float(r["logz"]), "ncall", int(r["ncall"]), f"{(time.time() - t0) / 60:.1f} min")


# ---------------------------------------------------------------------------
# Era-by-era tilts from the plates (for the figures; Zasche et al.'s Table 2
# style): the global fit's geometry and eclipse clock, only the tilt free
# ---------------------------------------------------------------------------

ERA_EDGES = {"SS Lac": [1890, 1905, 1915, 1925, 1935, 1945, 1955, 1990],
             "HS Hya": [1890, 1910, 1925, 1935, 1945, 1955, 1990],
             "V907 Sco": [1889, 1900, 1910, 1920, 1930, 1945, 1960, 1990]}


def era_inclinations(name: str, theta, edges=None, n_grid=181):
    """Profile likelihood of the inner orbit's tilt in each era of plates.
    Returns, per era: year range, number of plates, best tilt (deg, folded to
    0-90: eclipses cannot tell i from 180 - i) and the 1-sigma range, or an
    upper limit when no eclipse is seen."""
    fit = GlobalFit(name)
    m = fit.model
    p = m.unpack(theta)
    g = m.geometry(p)
    edges = edges or ERA_EDGES[name]
    yrs = C.jd_to_year(fit.t)
    out = []
    i_grid = np.linspace(60.0, 90.0, n_grid)
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (yrs >= lo) & (yrs < hi)
        if sel.sum() < 20:
            continue
        t = fit.t[sel]
        jd_mid = float(np.median(t))
        # orientation frozen at the era's middle, tilt replaced; the clock
        # (which conjunction when) stays the global model's
        lnl = []
        for ideg in i_grid:
            gx = dict(g)
            h_mid = g["h0"] + g["nu"] * (jd_mid - fit.sp.t_ref)
            # rotate the precession phase so that the tilt at mid-era is ideg
            gx["h0"], gx["eps_in"], gx["I"], gx["nu"] = 0.0, 0.0, np.radians(ideg), 0.0
            gg = g["g0"] + g["gdot"] * (jd_mid - fit.sp.t_ref)
            tc = fit.conjunction_time(g, jd_mid)
            gx["g0"], gx["gdot"] = gg, 0.0
            gx["lam0"] = _conj_lambda_fixed(gx, tc, fit.sp)
            F = m.flux(gx, t, band="B", L3=L3_B[name])
            mod = p["m0"] - 2.5 * np.log10(F)
            s2 = fit.err[sel] ** 2 + p["sig_jit"] ** 2
            r = fit.mag[sel] - mod
            lp_in = -0.5 * r * r / s2 - 0.5 * np.log(2 * np.pi * s2)
            so2 = p["sig_out"] ** 2 + fit.err[sel] ** 2
            lp_out = -0.5 * (r - p["mu_out"]) ** 2 / so2 - 0.5 * np.log(2 * np.pi * so2)
            f = p["f_out"]
            lnl.append(float(np.sum(np.logaddexp(np.log1p(-f) + lp_in, np.log(f + 1e-300) + lp_out))))
        lnl = np.array(lnl)
        k = int(np.argmax(lnl))
        ok = lnl > lnl.max() - 0.5
        ok95 = lnl > lnl.max() - 2.0
        no_ecl = lnl[0]            # i = 60 deg: no eclipse for these systems
        out.append(dict(lo=lo, hi=hi, n=int(sel.sum()), year=float(C.jd_to_year(jd_mid)),
                        i_best=float(i_grid[k]), i_lo=float(i_grid[ok].min()), i_hi=float(i_grid[ok].max()),
                        i_hi95=float(i_grid[ok95].max()), dlnl_vs_none=float(lnl.max() - no_ecl)))
    return out


def _conj_lambda_fixed(gx, tc, sp):
    """lam0 such that the frozen orbit has its primary conjunction at tc."""
    lam_c = S._conj_lambda(gx["h0"], gx["eps_in"], gx["I"], gx["g0"], sp.e)
    return lam_c - gx["n_lam"] * (tc - sp.t_ref)


# ---------------------------------------------------------------------------
# Summaries and predictions from a finished posterior
# ---------------------------------------------------------------------------

def _seasons(years, x_min):
    """Start/end years of the intervals with eclipses (impact parameter < 1)."""
    ecl = x_min < 1.0
    edges = np.flatnonzero(np.diff(ecl.astype(int)))
    out = []
    start = years[0] if ecl[0] else None
    for k in edges:
        if ecl[k + 1]:
            start = years[k + 1]
        else:
            out.append((start, years[k]))
            start = None
    if start is not None:
        out.append((start, None))
    return out


#: the posterior quoted for each system (see README section 5): "all" adds
#: Hipparcos/ASAS-3; "quietpro" is the fit outside the Kozai-Lidov window (the
#: unrestricted SS Lac fit lands in it, where the secular model fails)
FINAL_TAG = {"SS Lac": "quietpro", "HS Hya": "all", "V907 Sco": ""}


def posterior_path(name: str, tag: str = "") -> Path:
    return CACHE / f"posterior_{_slug(name)}{'_' + tag if tag else ''}.npz"


def alternative_solution(name: str):
    """The competing solution with the other sense of orbital motion, and
    the log odds of the final one against it. Prograde-only and
    retrograde-only runs cover equal prior masses; failing a prograde-only
    run, a main run over the full isotropic prior that found only prograde
    samples covers half the prior mass (Z_pro = 2 Z_main). The Kozai-quiet
    runs cover equal masses."""
    tag = FINAL_TAG[name]
    if tag.startswith("quiet"):
        alt = "quietretro" if tag == "quietpro" else "quietpro"
        half = 0.0
    else:
        alt = f"{tag}_retro" if tag else "retro"
        half = np.log(2.0)
    pm, pa = posterior_path(name, tag), posterior_path(name, alt)
    if not (pm.exists() and pa.exists()):
        return None
    dm, da = np.load(pm, allow_pickle=True), np.load(pa, allow_pickle=True)
    if half and np.mean(dm["posterior"][:, list(dm["names"]).index("J")] < 90) < 0.999:
        half = 0.0          # the main run sampled both senses itself
    lnZ_final, err = float(dm["logz"]) + half, float(dm["logzerr"])
    # a prograde-only run that found a *different* prograde maximum (runs get
    # trapped) adds its own evidence: distinct modes' evidences add
    pp = posterior_path(name, f"{tag}_pro" if tag else "pro")
    modes = [dict(tag=tag, lnZ_half=lnZ_final, I=float(np.median(dm["posterior"][:, 1])))]
    if half and pp.exists():
        dp = np.load(pp, allow_pickle=True)
        I_p = float(np.median(dp["posterior"][:, 1]))
        if abs(I_p - modes[0]["I"]) > 10.0:
            modes.append(dict(tag=f"{tag}_pro", lnZ_half=float(dp["logz"]), I=I_p))
            lnZ_final = float(np.logaddexp(lnZ_final, float(dp["logz"])))
    lnB = lnZ_final - float(da["logz"])
    return dict(tag=alt, lnB=lnB, lnB_err=float(np.hypot(err, float(da["logzerr"]))), final_modes=modes)


def summarize(name: str, n_draws=2000, seed=C.RANDOM_SEED, tag: str = ""):
    """Posterior summaries: parameters, derived quantities (masses, precession
    period, angles) and the eclipse calendar (past seasons, next return)."""
    d = np.load(posterior_path(name, tag), allow_pickle=True)
    post = d["posterior"]
    names = list(d["names"])
    rng = np.random.default_rng(seed)
    draws = post[rng.choice(len(post), min(n_draws, len(post)), replace=False)]
    m = S.Model(name)
    years = np.arange(1700.0, 2701.0, 0.5)
    jd = C.year_to_jd(years)
    der = {k: [] for k in ("P_node_yr", "eps_in", "eps_out", "gamma", "m1", "m2", "m3", "J", "I",
                           "i_1900", "i_2026", "last_end", "next_start", "season_len")}
    for th in draws:
        g = m.geometry(m.unpack(th))
        xp, xs = m.impact(g, jd)
        xmin = np.minimum(xp, xs)
        seas = _seasons(years, xmin)
        now = 2026.8
        past = [s for s in seas if s[1] is not None and s[1] < now]
        future = [s for s in seas if s[0] is not None and s[0] > now]
        cur = [s for s in seas if (s[0] is None or s[0] <= now) and (s[1] is None or s[1] >= now)]
        der["last_end"].append(past[-1][1] if past else np.nan)
        der["next_start"].append(future[0][0] if future else np.nan)
        lens = [s[1] - s[0] for s in seas if s[0] is not None and s[1] is not None]
        der["season_len"].append(np.median(lens) if lens else np.nan)
        for k in ("P_node_yr", "eps_in", "eps_out", "gamma", "m1", "m2", "m3", "J", "I"):
            v = g[k]
            der[k].append(np.degrees(v) if k in ("eps_in", "eps_out", "J", "I") else v)
        der["i_1900"].append(float(np.degrees(m.inclination(g, C.year_to_jd(1900.0)))))
        der["i_2026"].append(float(np.degrees(m.inclination(g, C.year_to_jd(2026.8)))))
        if cur:
            der.setdefault("current_end", []).append(cur[0][1] if cur[0][1] is not None else np.nan)
    def q(x):
        x = np.asarray(x, float)
        x = x[np.isfinite(x)]
        if len(x) == 0:
            return None
        a = np.percentile(x, [2.5, 16, 50, 84, 97.5])
        return dict(median=float(a[2]), lo=float(a[2] - a[1]), hi=float(a[3] - a[2]), p2_5=float(a[0]),
                    p97_5=float(a[4]), n=int(len(x)))
    out = {"logZ": float(d["logz"]), "logZ_err": float(d["logzerr"]), "ncall": int(d["ncall"]),
           "n_posterior": int(len(post))}
    out["params"] = {n: q(post[:, k]) for k, n in enumerate(names)}
    out["derived"] = {k: q(v) for k, v in der.items()}
    out["frac_prograde"] = float(np.mean(post[:, names.index("J")] < 90))
    out["frac_no_next_before_2700"] = float(np.mean(~np.isfinite(der["next_start"])))
    return out


def ground_profile(name: str, lc: dict, theta) -> dict | None:
    """Profile for Hipparcos/ASAS data: constant + two-harmonic baseline, then
    clip points more than 5 robust sigma from the best model (survey outliers;
    true eclipse points follow the model and stay) and refit."""
    sf = S.SectorFit(name, lc, theta, nharm=2, baseline="constant")
    pr = sf.profile()
    if pr is None:
        return None
    _, fit, _ = sf.chi2(pr["x_best"], pr["dt_pri"], pr["dt_sec"], return_model=True)
    r = lc["flux"] - fit
    sig = 1.4826 * np.median(np.abs(r - np.median(r)))
    keep = np.abs(r) < 5 * sig
    lc2 = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in lc.items()}
    sf2 = S.SectorFit(name, lc2, theta, nharm=2, baseline="constant")
    pr2 = sf2.profile()
    pr2["label"] = lc["label"]
    pr2["n_clipped"] = int((~keep).sum())
    pr2["source"] = "ground"
    return pr2


def ground_profiles(name: str, theta, force=False) -> list[dict]:
    path = CACHE / f"ground_profiles_{_slug(name)}.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    out = []
    for lc in S.ground_lightcurves(name):
        pr = ground_profile(name, lc, theta)
        if pr is not None:
            out.append(pr)
            print(f"{name} {lc['label']}: x_best {pr['x_best']:.3f} detected {pr['detected']} "
                  f"scale {pr['scale']:.2f} clipped {pr['n_clipped']}", flush=True)
    path.write_text(json.dumps(out))
    return out


def ground_check(name: str, tag: str = "", n: int = 400, seed=C.RANDOM_SEED) -> list[dict]:
    """Posterior predictive check against the Hipparcos/ASAS-3 profiles: the
    tilt each ground data set measures (best x and its delta-chi2 < 1 range,
    converted with each draw's radii) against the tilt the fit predicts."""
    path = CACHE / f"ground_profiles_{_slug(name)}.json"
    if not path.exists():
        return []
    gp = json.loads(path.read_text())
    post = np.load(posterior_path(name, tag), allow_pickle=True)["posterior"]
    rng = np.random.default_rng(seed)
    draws = post[rng.choice(len(post), min(n, len(post)), replace=False)]
    m = S.Model(name)
    geo = [m.geometry(m.unpack(th)) for th in draws]
    out = []
    for pr in gp:
        x = np.array(pr["x"])
        c = (np.array(pr["chi2"]) - pr["chi2_min"]) / pr["scale"]
        xp = np.array([m.impact(g, pr["jd_mid"])[0] for g in geo])
        rs = np.array([g["r1"] + g["r2"] for g in geo])
        dchi = np.interp(xp, x, c)

        def tilt(xx):
            return np.degrees(np.arccos(np.clip(np.asarray(xx) * rs, -1, 1)))
        ok = x[c < 1]
        i_obs = tilt(pr["x_best"])
        q = lambda a: [float(v) for v in np.percentile(a, [16, 50, 84])]
        out.append(dict(label=pr["label"], year=float(C.jd_to_year(pr["jd_mid"])), detected=bool(pr["detected"]),
                        x_obs=float(pr["x_best"]), x_obs_range=[float(ok.min()), float(ok.max())],
                        tilt_obs=q(i_obs), tilt_obs_range=[float(np.median(tilt(ok.max()))), float(np.median(tilt(ok.min())))],
                        tilt_pred=q(tilt(xp)), x_pred=q(xp), dchi2_median=float(np.median(dchi)),
                        dchi2_p16_p84=[float(np.percentile(dchi, 16)), float(np.percentile(dchi, 84))]))
    return out
