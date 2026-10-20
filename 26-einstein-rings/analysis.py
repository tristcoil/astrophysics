"""From rings to masses to dark matter: every number in the README comes from here.

The steps, in the order the README tells them:

1. ``lens_table``      Einstein radius + two redshifts -> mass inside the ring,
                       and the isothermal velocity dispersion sigma_SIE.
2. ``isothermal_test`` Is sigma_SIE the stars' own velocity dispersion?
3. ``density_slopes``  Lensing mass + stellar dispersion + Jeans equation ->
                       how steeply each galaxy's total density falls.
4. ``dark_fractions``  Stellar masses (two IMFs) against the lensing mass.
5. ``decompose``       Stars + an NFW dark halo fitted to lensing + dynamics:
                       how heavy must the stars be?
6. ``discovery``       The SDSS spectra: two redshifts in one fibre, and
                       which candidates turned out to be lenses.
7. ``our_einstein_radii`` Our own lens models against the published ones.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import catalog
import lensing as L
from config import (EINSTEIN_RADIUS_ERROR, FIBRE_RADIUS_ARCSEC, GAMMA_RANGE,
                    HERNQUIST_RE_OVER_A, NFW_BREAK_RADIUS_KPC, PLANCK_H0, PLANCK_OMEGA_M,
                    SEEING_FWHM_ARCSEC, SIGMA_ERROR_FLOOR, SLACS_H0, SLACS_OMEGA_M)

EARLY_TYPES = {"E", "S0", "E/S0"}  # SLACS IX morphologies counted as early-type


# --------------------------------------------------------------------------
# 1. The lens table
# --------------------------------------------------------------------------

def lens_table(h0: float = SLACS_H0, omega_m: float = SLACS_OMEGA_M) -> pd.DataFrame:
    """One row per grade-A lens with everything derived from ring + redshifts."""
    g = catalog.read_grade_a().copy()
    d = L.distances(g["z_lens"], g["z_source"], h0, omega_m)
    for k, v in d.items():
        g[k] = v
    g["kpc_per_arcsec"] = L.kpc_per_arcsec(g["d_l"])

    # Einstein radius in arcsec: SLACS V's SIE model where there is one; for
    # the lenses modelled only in SLACS IX, their kpc value converted back with
    # their own cosmology (H0 = 70, Omega_m = 0.3).
    d70 = L.distances(g["z_lens"], g["z_source"], SLACS_H0, SLACS_OMEGA_M)
    theta_ix = g["r_ein_kpc"] / L.kpc_per_arcsec(d70["d_l"])
    g["theta_e"] = g["b_sie"].where(g["b_sie"].notna(), theta_ix)
    g["theta_e_from"] = np.where(g["b_sie"].notna(), "SLACS V",
                                 np.where(theta_ix.notna(), "SLACS IX", ""))

    g["r_ein_kpc_ours"] = g["theta_e"] * g["kpc_per_arcsec"]
    g["m_ein"] = L.einstein_mass(g["theta_e"], g["d_l"], g["d_s"], g["d_ls"])
    g["sigma_cr"] = L.sigma_crit(g["d_l"], g["d_s"], g["d_ls"])
    g["sigma_sie"] = L.sigma_sie(g["theta_e"], g["d_s"], g["d_ls"])

    # Stellar velocity dispersion inside the 3" fibre (SLACS V where it has
    # one, else SLACS IX), and the de Vaucouleurs effective radius in F814W.
    cand = catalog.read_candidates().set_index("name")
    sig_v = cand["sigma"].reindex(g["name"]).to_numpy(dtype=float)
    esig_v = cand["e_sigma"].reindex(g["name"]).to_numpy(dtype=float)
    g["sigma_star"] = np.where(np.isfinite(sig_v), sig_v, g["sigma"])
    g["e_sigma_star"] = np.where(np.isfinite(sig_v), esig_v, g["e_sigma"])
    g["r_eff"] = g["re_arcsec"].where(g["re_arcsec"].notna(), g["re_i"])
    g["sigma_e8"] = L.jorgensen_sigma(g["sigma_star"], FIBRE_RADIUS_ARCSEC, g["r_eff"], 1 / 8)
    g["sigma_e2"] = L.jorgensen_sigma(g["sigma_star"], FIBRE_RADIUS_ARCSEC, g["r_eff"], 1 / 2)

    # Stellar mass inside the ring, for each IMF: total stellar mass (SLACS IX,
    # stellar-population fits to HST colours, H0 = 70) rescaled to our H0 and
    # times the share of de Vaucouleurs light inside the Einstein radius.
    lum_scale = (SLACS_H0 / h0) ** 2  # luminosity distance squared
    g["light_frac_ein"] = L.sersic_light_fraction(g["theta_e"], g["r_eff"], 4.0)
    for imf in ("chab", "salp"):
        mstar = 10 ** g[f"log_mstar_{imf}"] * lum_scale
        g[f"mstar_{imf}"] = mstar
        g[f"mstar_ein_{imf}"] = mstar * g["light_frac_ein"]
        g[f"fstar_ein_{imf}_ours"] = g[f"mstar_ein_{imf}"] / g["m_ein"]
    g["early_type"] = g["morph_type"].isin(EARLY_TYPES)
    return g


def reproduce_slacs_ix(t: pd.DataFrame) -> dict:
    """Our ring sizes, ring masses and stellar fractions against SLACS IX's table."""
    ok = t["r_ein_kpc"].notna() & t["theta_e"].notna()
    dr = t.loc[ok, "r_ein_kpc_ours"] - t.loc[ok, "r_ein_kpc"]
    dm = np.log10(t.loc[ok, "m_ein"]) - t.loc[ok, "log_m_ein"]
    fs = t["fstar_chab"].notna() & ok
    dfc = t.loc[fs, "fstar_ein_chab_ours"] - t.loc[fs, "fstar_chab"]
    dfs = t.loc[fs, "fstar_ein_salp_ours"] - t.loc[fs, "fstar_salp"]
    return {
        "n": int(ok.sum()),
        "r_ein_max_abs_diff_kpc": float(dr.abs().max()),
        "log_m_ein_median_diff": float(dm.median()), "log_m_ein_max_abs_diff": float(dm.abs().max()),
        "n_fstar": int(fs.sum()),
        "fstar_chab_median_diff": float(dfc.median()), "fstar_chab_rms_diff": float(np.sqrt(np.mean(dfc**2))),
        "fstar_salp_median_diff": float(dfs.median()), "fstar_salp_rms_diff": float(np.sqrt(np.mean(dfs**2))),
        "fstar_chab_mean": float(t.loc[fs, "fstar_ein_chab_ours"].mean()),
        "fstar_chab_std": float(t.loc[fs, "fstar_ein_chab_ours"].std()),
        "fstar_salp_mean": float(t.loc[fs, "fstar_ein_salp_ours"].mean()),
        "fstar_salp_std": float(t.loc[fs, "fstar_ein_salp_ours"].std()),
        "m_ein_range": [float(t.loc[ok, "m_ein"].min()), float(t.loc[ok, "m_ein"].max())],
        "r_ein_kpc_range": [float(t.loc[ok, "r_ein_kpc_ours"].min()), float(t.loc[ok, "r_ein_kpc_ours"].max())],
        "theta_e_range": [float(t.loc[ok, "theta_e"].min()), float(t.loc[ok, "theta_e"].max())],
    }


def cosmology_dependence(t70: pd.DataFrame) -> dict:
    """How much the choice of H0 and Omega_m moves each quantity."""
    tp = lens_table(PLANCK_H0, PLANCK_OMEGA_M)
    ok = t70["theta_e"].notna()
    return {
        "sigma_sie_ratio_median": float(np.median(tp.loc[ok, "sigma_sie"] / t70.loc[ok, "sigma_sie"])),
        "m_ein_ratio_median": float(np.median(tp.loc[ok, "m_ein"] / t70.loc[ok, "m_ein"])),
        "fstar_ratio_median": float(np.nanmedian(tp.loc[ok, "fstar_ein_chab_ours"]
                                                 / t70.loc[ok, "fstar_ein_chab_ours"])),
    }


# --------------------------------------------------------------------------
# 2. Is the lensing dispersion the stars' dispersion?
# --------------------------------------------------------------------------

def _proportional_fit(y, x, ey):
    """Least-squares f in y = f x, with errors on y; returns f, its error, rms of y/x - f,
    and the reduced chi^2."""
    w = 1 / ey**2
    f = np.sum(w * x * y) / np.sum(w * x * x)
    ef = 1 / np.sqrt(np.sum(w * x * x))
    chi2 = np.sum(w * (y - f * x) ** 2) / (len(y) - 1)
    return float(f), float(ef), float(np.std(y / x - f)), float(chi2)


def isothermal_test(t: pd.DataFrame) -> dict:
    """Bolton et al. (2008b) section 4: sigma_stars against sigma_SIE for the
    early-type grade-A lenses with SIE models and well-measured dispersions."""
    sel = (t["has_sie_model"] & t["good_sigma"].eq(True)
           & (t["morphology"] == "E") & t["sigma_star"].notna())
    s = t[sel]
    ey = np.maximum(s["e_sigma_star"], SIGMA_ERROR_FLOOR * s["sigma_star"])
    f_sdss = _proportional_fit(s["sigma_star"].to_numpy(), s["sigma_sie"].to_numpy(), ey.to_numpy())
    ey8 = ey * s["sigma_e8"] / s["sigma_star"]
    f_e8 = _proportional_fit(s["sigma_e8"].to_numpy(), s["sigma_sie"].to_numpy(), ey8.to_numpy())

    # Intrinsic scatter: the fraction of sigma_SIE that, added in quadrature to
    # the measurement errors, brings chi^2/dof to 1.
    def chi2_with(frac):
        e = np.sqrt(ey8**2 + (frac * s["sigma_sie"]) ** 2)
        w = 1 / e**2
        f = np.sum(w * s["sigma_e8"] * s["sigma_sie"]) / np.sum(w * s["sigma_sie"] ** 2)
        return np.sum(w * (s["sigma_e8"] - f * s["sigma_sie"]) ** 2) / (len(s) - 1)
    fracs = np.linspace(0, 0.3, 3001)
    c = np.array([chi2_with(fr) for fr in fracs])
    intrinsic = float(fracs[np.argmin(np.abs(c - 1))])

    # Does the ratio depend on Einstein radius / effective radius (it should
    # not, if the density is isothermal at every radius)?
    ratio = s["sigma_e8"] / s["sigma_sie"]
    rho = pd.Series(ratio.to_numpy()).corr(pd.Series((s["theta_e"] / s["r_eff"]).to_numpy()),
                                           method="spearman")
    return {
        "n": int(sel.sum()), "names": list(s["name"]),
        "f_sdss": f_sdss, "f_e8": f_e8, "intrinsic_scatter_frac": intrinsic,
        "spearman_ratio_vs_rein_over_reff": float(rho),
        "sigma_sie_range": [float(s["sigma_sie"].min()), float(s["sigma_sie"].max())],
    }


# --------------------------------------------------------------------------
# 3. Total density slopes from lensing + stellar dynamics
# --------------------------------------------------------------------------

GAMMA_GRID = np.round(np.arange(GAMMA_RANGE[0], GAMMA_RANGE[1] + 1e-9, 0.01), 3)


def _posterior_summary(grid, post):
    cdf = np.cumsum(post)
    cdf /= cdf[-1]
    q = np.interp([0.16, 0.5, 0.84], cdf, grid)
    return float(q[1]), float(q[1] - q[0]), float(q[2] - q[1])


def density_slope(row, beta: float = 0.0, seeing: float = SEEING_FWHM_ARCSEC,
                  grid=GAMMA_GRID) -> dict:
    """Koopmans et al. (2006, 2009): a power-law total density rho ~ r^-gamma,
    normalised to the lensing mass inside the Einstein radius, with the stars a
    Hernquist tracer of the light. The Jeans equation predicts the dispersion
    the SDSS fibre would see for each gamma; the measured one picks gamma."""
    sis = float(L.sigma_sie(row.theta_e, row.d_s, row.d_ls))
    jm = L.JeansModel(row.r_eff, FIBRE_RADIUS_ARCSEC, seeing, beta=beta)
    s_mod = np.array([jm.sigma_fibre(L.powerlaw_v2(g, row.theta_e, sis)) for g in grid])
    # The 2% Einstein-radius error moves sigma by 1%; add it to the measurement error.
    err = np.sqrt(row.e_sigma_star**2 + (0.5 * EINSTEIN_RADIUS_ERROR * s_mod) ** 2)
    like = np.exp(-0.5 * ((row.sigma_star - s_mod) / err) ** 2)
    if like.sum() == 0 or not np.isfinite(like.sum()):
        like = np.zeros_like(grid)
        like[np.argmin(np.abs(row.sigma_star - s_mod))] = 1.0
    med, lo, hi = _posterior_summary(grid, like)
    return {"name": row.name, "gamma": med, "gamma_lo": lo, "gamma_hi": hi,
            "sigma_model_iso": float(np.interp(2.0, grid, s_mod)),
            "like": like, "sigma_model": s_mod}


def hierarchical_gaussian(likes: np.ndarray, grid: np.ndarray, mus: np.ndarray,
                          scatters: np.ndarray) -> dict:
    """Population mean and intrinsic scatter of a quantity measured per object
    with likelihoods ``likes`` (n_obj x n_grid) on ``grid``, for a Gaussian
    population: L(mu, s) = prod_i sum_x like_i(x) N(x | mu, s)."""
    from scipy.special import ndtr

    likes = likes / likes.sum(axis=1, keepdims=True)
    # Population probability of each grid cell, from the Gaussian's CDF at the
    # cell edges: exact however narrow the Gaussian is compared with the grid.
    edges = np.concatenate([[grid[0] - 0.5 * (grid[1] - grid[0])],
                            0.5 * (grid[1:] + grid[:-1]),
                            [grid[-1] + 0.5 * (grid[-1] - grid[-2])]])
    s_eff = np.maximum(scatters, 1e-6)
    z = (edges[None, None, :] - mus[:, None, None]) / s_eff[None, :, None]
    pop = np.diff(ndtr(z), axis=2)  # (n_mu, n_s, n_grid)
    logl = np.zeros((len(mus), len(scatters)))
    for li in likes:
        logl += np.log(np.maximum(np.tensordot(pop, li, axes=([2], [0])), 1e-300))
    i, j = np.unravel_index(np.argmax(logl), logl.shape)
    prob = np.exp(logl - logl.max())
    pm = prob.sum(axis=1)
    ps = prob.sum(axis=0)
    mu_med, mu_lo, mu_hi = _posterior_summary(mus, pm)
    s_med, s_lo, s_hi = _posterior_summary(scatters, ps)
    return {"mu_ml": float(mus[i]), "scatter_ml": float(scatters[j]),
            "mu": mu_med, "mu_err": 0.5 * (mu_lo + mu_hi), "mu_lo": mu_lo, "mu_hi": mu_hi,
            "scatter": s_med, "scatter_lo": s_lo, "scatter_hi": s_hi,
            "mu_grid": mus, "mu_post": pm / pm.max(), "s_grid": scatters, "s_post": ps / ps.max()}


def slope_sample(t: pd.DataFrame) -> pd.DataFrame:
    """Early-type lenses with an Einstein radius, a dispersion and an effective radius."""
    sel = (t["early_type"] & t["theta_e"].notna() & t["sigma_star"].notna()
           & t["e_sigma_star"].notna() & t["r_eff"].notna())
    # SLACS V flags dispersions it does not trust; honour the flag where it exists.
    sel &= ~t["good_sigma"].eq(False)
    return t[sel]


def density_slopes(t: pd.DataFrame, beta: float = 0.0, seeing: float = SEEING_FWHM_ARCSEC) -> dict:
    s = slope_sample(t)
    per = [density_slope(r, beta=beta, seeing=seeing) for r in s.itertuples()]
    likes = np.array([p["like"] for p in per])
    pop = hierarchical_gaussian(likes, GAMMA_GRID, np.arange(1.85, 2.35, 0.002),
                                np.arange(0.0, 0.40, 0.004))
    table = pd.DataFrame([{k: v for k, v in p.items() if k not in ("like", "sigma_model")} for p in per])
    return {"table": table, "population": pop, "likes": likes, "n": len(per), "beta": beta,
            "seeing": seeing}


# --------------------------------------------------------------------------
# 4. How much of the mass inside the ring is dark?
# --------------------------------------------------------------------------

def dark_fractions(t: pd.DataFrame, slopes: dict | None = None) -> dict:
    """Stellar mass inside the Einstein radius against the lensing mass, for
    Chabrier and Salpeter IMFs; and the heaviest IMF each lens allows."""
    s = t[t["theta_e"].notna() & t["mstar_chab"].notna()].copy()
    out = {"n": len(s)}
    for imf in ("chab", "salp"):
        f = s[f"fstar_ein_{imf}_ours"]
        # error: 0.1-dex-level stellar-mass errors (SLACS IX) plus 2% on M_E
        ef = f * np.log(10) * np.hypot(s[f"e_log_mstar_{imf}"], EINSTEIN_RADIUS_ERROR / np.log(10))
        out[imf] = {
            "fstar_mean": float(f.mean()), "fstar_median": float(f.median()),
            "fstar_std": float(f.std()),
            "fdm_mean": float(1 - f.mean()),
            "n_fstar_above_1": int((f > 1).sum()),
            "n_fstar_above_1_by_2sigma": int(((f - 1) / ef > 2).sum()),
            "names_above_1": list(s.loc[f > 1, "name"]),
        }
    # The heaviest stellar population each lens allows, relative to Chabrier:
    # all of the lensing mass in stars, none dark.
    s["alpha_max_chab"] = 1 / s["fstar_ein_chab_ours"]
    salp_over_chab = 10 ** np.median(s["log_mstar_salp"] - s["log_mstar_chab"])
    out["salp_over_chab"] = float(salp_over_chab)
    out["alpha_max_median"] = float(s["alpha_max_chab"].median())
    out["alpha_max_min"] = float(s["alpha_max_chab"].min())
    out["frac_lenses_allowing_salpeter"] = float((s["alpha_max_chab"] > salp_over_chab).mean())
    # Does the dark share grow with mass (Auger et al. 2010)?
    from scipy.stats import spearmanr

    fdm = 1 - s["fstar_ein_chab_ours"]
    rho, p = spearmanr(fdm, s["sigma_star"], nan_policy="omit")
    out["spearman_fdm_vs_sigma"], out["p_fdm_vs_sigma"] = float(rho), float(p)
    rho, p = spearmanr(fdm, s["mstar_chab"], nan_policy="omit")
    out["spearman_fdm_vs_mstar"], out["p_fdm_vs_mstar"] = float(rho), float(p)
    out["table"] = s[["name", "theta_e", "r_eff", "m_ein", "mstar_chab", "mstar_salp",
                      "mstar_ein_chab", "mstar_ein_salp",
                      "fstar_ein_chab_ours", "fstar_ein_salp_ours", "fstar_chab", "fstar_salp",
                      "alpha_max_chab", "sigma_star", "light_frac_ein"]]

    # Dark fraction inside r_e/2 (Auger et al. 2010's aperture), extrapolating
    # the lensing mass with each lens's own density slope.
    if slopes is not None:
        sl = slopes["table"].set_index("name")
        s2 = s[s["name"].isin(sl.index)].copy()
        gam = sl.loc[s2["name"], "gamma"].to_numpy()
        r_half = 0.5 * s2["r_eff"]
        m_re2 = s2["m_ein"] * (r_half / s2["theta_e"]) ** (3 - gam)
        frac_re2 = L.sersic_light_fraction(r_half, s2["r_eff"], 4.0)  # 32%
        for imf in ("chab", "salp"):
            s2[f"fdm_re2_{imf}"] = 1 - s2[f"mstar_{imf}"] * frac_re2 / m_re2
        s2["log_m_re2"] = np.log10(m_re2)
        out["re2_table"] = s2[["name", "log_m_re2", "fdm_re2_chab", "fdm_re2_salp", "sigma_star",
                               "mstar_chab"]]
        out["fdm_re2_chab_mean"] = float(s2["fdm_re2_chab"].mean())
        out["fdm_re2_salp_mean"] = float(s2["fdm_re2_salp"].mean())
    return out


# --------------------------------------------------------------------------
# 5. Stars + dark halo: how heavy are the stars?
# --------------------------------------------------------------------------

FSTAR_GRID = np.linspace(0.0, 1.0, 201)
LOG_ALPHA_GRID = np.round(np.arange(-1.2, 0.81, 0.01), 3)


def decompose(row, beta: float = 0.0, r_break_kpc: float = NFW_BREAK_RADIUS_KPC) -> dict:
    """Treu et al. (2010): stars as a Hernquist profile plus an NFW dark halo with
    a 30-kpc break radius. Lensing fixes their sum inside the Einstein radius;
    the stars' dispersion fixes the split. Returns the posterior of the stellar
    fraction f* inside the ring and of the IMF mismatch parameter
    alpha = (stellar mass from lensing + dynamics) / (stellar mass for a Salpeter IMF)."""
    sis = float(L.sigma_sie(row.theta_e, row.d_s, row.d_ls))
    a = row.r_eff / HERNQUIST_RE_OVER_A
    r_s = r_break_kpc / row.kpc_per_arcsec
    jm = L.JeansModel(row.r_eff, FIBRE_RADIUS_ARCSEC, SEEING_FWHM_ARCSEC, beta=beta)
    s_mod = np.array([jm.sigma_fibre(L.two_component_v2(f, row.theta_e, sis, a, r_s))
                      for f in FSTAR_GRID])
    err = np.sqrt(row.e_sigma_star**2 + (0.5 * EINSTEIN_RADIUS_ERROR * s_mod) ** 2)
    like_f = np.exp(-0.5 * ((row.sigma_star - s_mod) / err) ** 2)
    f_med, f_lo, f_hi = _posterior_summary(FSTAR_GRID, like_f)

    # log alpha = log(f* M_E) - log(M*_Salp inside R_E). The stellar-population
    # mass carries its own ~0.1 dex error, folded in by convolution.
    log_mld = np.log10(np.maximum(FSTAR_GRID, 1e-6) * row.m_ein)
    log_msps = np.log10(row.mstar_ein_salp)
    e_sps = row.e_log_mstar_salp if np.isfinite(row.e_log_mstar_salp) else 0.1
    la = log_mld - log_msps
    like_a = np.zeros_like(LOG_ALPHA_GRID)
    for la_i, w in zip(la[1:], like_f[1:]):  # f* = 0 means no stars: log alpha = -inf
        like_a += w * np.exp(-0.5 * ((LOG_ALPHA_GRID - la_i) / e_sps) ** 2)
    a_med, a_lo, a_hi = _posterior_summary(LOG_ALPHA_GRID, like_a)
    return {"name": row.name, "fstar": f_med, "fstar_lo": f_lo, "fstar_hi": f_hi,
            "log_alpha": a_med, "log_alpha_lo": a_lo, "log_alpha_hi": a_hi,
            "like_alpha": like_a, "like_f": like_f, "sigma_model": s_mod,
            "fstar_best": float(FSTAR_GRID[np.argmax(like_f)])}


def decompositions(t: pd.DataFrame, beta: float = 0.0,
                   r_break_kpc: float = NFW_BREAK_RADIUS_KPC) -> dict:
    s = slope_sample(t)
    s = s[s["mstar_salp"].notna()]
    per = [decompose(r, beta=beta, r_break_kpc=r_break_kpc) for r in s.itertuples()]
    likes = np.array([p["like_alpha"] for p in per])
    pop = hierarchical_gaussian(likes, LOG_ALPHA_GRID, np.arange(-0.3, 0.3, 0.002),
                                np.arange(0.0, 0.30, 0.004))
    table = pd.DataFrame([{k: v for k, v in p.items()
                           if k not in ("like_alpha", "like_f", "sigma_model")} for p in per])
    table = table.merge(s[["name", "sigma_star", "sigma_sie", "mstar_salp"]], on="name")
    # Trend with velocity dispersion: least squares on the medians, weighted.
    e = 0.5 * (table["log_alpha_lo"] + table["log_alpha_hi"])
    x = np.log10(table["sigma_star"])
    w = 1 / e**2
    A = np.vstack([x - np.log10(250), np.ones_like(x)]).T
    cov = np.linalg.inv(A.T @ (A * w.to_numpy()[:, None]))
    coef = cov @ (A.T @ (w * table["log_alpha"]).to_numpy())
    salp_over_chab = float(np.median(s["log_mstar_salp"] - s["log_mstar_chab"]))
    return {"table": table, "population": pop, "n": len(per), "beta": beta,
            "r_break_kpc": r_break_kpc,
            "trend_slope": float(coef[0]), "trend_slope_err": float(np.sqrt(cov[0, 0])),
            "trend_at_250": float(coef[1]),
            "log_salp_over_chab": salp_over_chab}


def conspiracy_profiles(row, fstar: float, r_break_kpc: float = NFW_BREAK_RADIUS_KPC) -> dict:
    """3D densities of the stars, the dark halo and their sum for one lens, in
    units where the isothermal sphere is flat (rho r^2)."""
    a = row.r_eff / HERNQUIST_RE_OVER_A
    r_s = r_break_kpc / row.kpc_per_arcsec
    r = np.geomspace(0.03, 10, 300) * row.r_eff  # arcsec
    m_star = fstar / L.hernquist_mass2d(row.theta_e, a)
    m_dm = (1 - fstar) / L.nfw_mass2d(row.theta_e, r_s)
    rho_star = m_star * L.hernquist_density(r, a)
    x = r / r_s
    rho_dm = m_dm / (4 * np.pi * r_s**3) / (x * (1 + x) ** 2)
    rho = rho_star + rho_dm
    slope = -np.gradient(np.log(rho), np.log(r))
    sel = (r > 0.2 * row.r_eff) & (r < 1.3 * row.r_eff)
    return {"r_over_re": r / row.r_eff, "rho_star": rho_star, "rho_dm": rho_dm, "rho": rho,
            "slope": slope, "mean_slope": float(np.mean(slope[sel])),
            "theta_e_over_re": float(row.theta_e / row.r_eff)}


# --------------------------------------------------------------------------
# 6. Discovery: the SDSS spectra
# --------------------------------------------------------------------------

def discovery() -> dict:
    import spectra
    cand = catalog.read_candidates()
    rows = []
    for r in cand.itertuples():
        sp = spectra.read_spectrum(r.name)
        res = spectra.search_background(sp, r.z_lens)
        rows.append({"name": r.name, "grade": r.grade, "z_lens": r.z_lens, "z_source": r.z_source,
                     "z_found": res.z_best, "score": res.score_best,
                     "sn_oii": res.lines_at_best["[O II]"], "sn_median": sp.sn_median,
                     "vdisp_dr17": sp.vdisp, "e_vdisp_dr17": sp.vdisp_err,
                     "sigma": r.sigma, "e_sigma": r.e_sigma})
    df = pd.DataFrame(rows)
    df["recovered"] = (df["z_found"] - df["z_source"]).abs() < 0.002
    # Predicted Einstein radius from the spectrum alone (isothermal sphere with
    # the measured dispersion): the information SLACS had before Hubble looked.
    d = L.distances(df["z_lens"], df["z_source"], SLACS_H0, SLACS_OMEGA_M)
    df["theta_pred"] = L.theta_e_sis(df["sigma"].astype(float), d["d_s"], d["d_ls"])
    bins = np.array([0.0, 0.8, 1.0, 1.2, 1.4, 3.0])
    has = df["theta_pred"].notna()
    rate = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = has & (df["theta_pred"] >= lo) & (df["theta_pred"] < hi)
        rate.append({"lo": lo, "hi": hi, "n": int(m.sum()),
                     "frac_a": float((df.loc[m, "grade"] == "A").mean()) if m.any() else np.nan,
                     "frac_ab": float(df.loc[m, "grade"].isin(["A", "B"]).mean()) if m.any() else np.nan})
    ok = df["sigma"].notna()
    return {
        "table": df, "n": len(df),
        "n_recovered": int(df["recovered"].sum()),
        "recovered_by_grade": df.groupby("grade")["recovered"].agg(["sum", "count"]).to_dict("index"),
        "missed": df.loc[~df["recovered"], ["name", "grade", "z_lens", "z_source", "z_found",
                                             "score"]].to_dict("records"),
        "confirmation_vs_theta": rate,
        "vdisp_ratio_median": float(np.median(df.loc[ok, "vdisp_dr17"] / df.loc[ok, "sigma"])),
        "vdisp_ratio_p16_p84": [float(v) for v in np.percentile(
            df.loc[ok, "vdisp_dr17"] / df.loc[ok, "sigma"], [16, 84])],
    }


def isothermal_test_with_ours(t: pd.DataFrame, ours: dict) -> dict:
    """The isothermal test again, with our own Einstein radii (unflagged fits)
    in place of the published ones."""
    df = ours["table"]
    mine = df[(df["failed"] == "") & (df["flag"] == "")].set_index("name")["b_ours"]
    t2 = t.copy()
    t2 = t2[t2["name"].isin(mine.index)]
    t2["theta_e"] = mine.reindex(t2["name"]).to_numpy()
    t2["sigma_sie"] = L.sigma_sie(t2["theta_e"], t2["d_s"], t2["d_ls"])
    return isothermal_test(t2)


# --------------------------------------------------------------------------
# 7. Our own Einstein radii
# --------------------------------------------------------------------------

def _fit_flag(res: dict) -> str:
    """A fit's own warning signs, judged without the published model: a
    parameter pinned at its bound means the data pushed it out of range."""
    import lensmodel
    if res.get("failed"):
        return "failed"
    if res["b"] <= lensmodel.B_MIN + 0.01 or res["b"] >= 2.95:
        return "Einstein radius at its bound"
    if res["q"] <= (1 - lensmodel.MASS_E_MAX) / (1 + lensmodel.MASS_E_MAX) + 0.01:
        return "mass flattening at its bound"
    if res["source_r"] >= lensmodel.SOURCE_R_MAX - 0.005:
        return "source size at its bound"
    return ""


def our_einstein_radii(t: pd.DataFrame, variant: str = "") -> dict:
    import lensmodel
    rows = []
    for r in t.itertuples():
        if variant == "combined":  # automatic fits, with the hand-masked refits where they exist
            res = lensmodel.cached_result(r.name, "hand_masked") or lensmodel.cached_result(r.name)
        else:
            res = lensmodel.cached_result(r.name, variant)
        if res is None:
            continue
        rows.append({"name": r.name, "failed": res.get("failed", ""),
                     "b_ours": res.get("b", np.nan), "q_ours": res.get("q", np.nan),
                     "pa_ours": res.get("pa", np.nan), "source_sn": res.get("source_sn", np.nan),
                     "chi2_red": res.get("chi2_red", np.nan), "b_initial": res.get("b_initial", np.nan),
                     "shear": res.get("shear", np.nan),
                     "b_pub": r.b_sie if np.isfinite(r.b_sie) else r.theta_e,
                     "b_pub_from": r.theta_e_from, "q_pub": r.q_sie, "pa_pub": r.pa_sie,
                     "ring": bool(r.ring_subset) if isinstance(r.ring_subset, (bool, np.bool_)) else False,
                     "light_q": res.get("light_q", np.nan), "light_pa": res.get("light_pa", np.nan),
                     "morph": r.morph_type, "flag": _fit_flag(res),
                     "hand_masked": res.get("manual_mask_pixels", 0) > 0})
    df = pd.DataFrame(rows)
    if df.empty:
        return {"table": df, "n": 0}
    ok = (df["failed"] == "") & df["b_pub"].notna()
    d = (df["b_ours"] / df["b_pub"] - 1)[ok]
    out = {"table": df, "n_fitted": int((df["failed"] == "").sum()), "n_failed": int((df["failed"] != "").sum()),
           "n_compared": int(ok.sum()), "n_flagged": int((df["flag"] != "").sum()),
           "flagged": df.loc[df["flag"] != "", ["name", "flag"]].to_dict("records"),
           "hand_masked": list(df.loc[df["hand_masked"], "name"])}
    good = ok & (df["flag"] == "")
    dg = (df["b_ours"] / df["b_pub"] - 1)[good]
    if good.any():
        out.update({"n_unflagged": int(good.sum()), "unflagged_median_frac_diff": float(np.median(dg)),
                    "unflagged_robust_scatter": float(1.4826 * np.median(np.abs(dg - np.median(dg)))),
                    "unflagged_n_within_5pct": int((dg.abs() < 0.05).sum()),
                    "unflagged_n_off_by_10pct": int((dg.abs() > 0.10).sum())})
    if ok.any():
        mad = 1.4826 * np.median(np.abs(d - np.median(d)))
        out.update({
            "median_frac_diff": float(np.median(d)), "robust_scatter": float(mad),
            "rms_frac_diff": float(np.sqrt(np.mean(d**2))),
            "n_within_2pct": int((d.abs() < 0.02).sum()), "n_within_5pct": int((d.abs() < 0.05).sum()),
            "n_off_by_10pct": int((d.abs() > 0.10).sum()),
            "names_off_by_10pct": list(df.loc[ok, "name"][d.abs() > 0.10]),
        })
        both = ok & df["q_pub"].notna()
        out["q_rms_diff"] = float(np.sqrt(np.mean((df.loc[both, "q_ours"] - df.loc[both, "q_pub"]) ** 2)))
        ring = ok & df["ring"]
        if ring.any():
            dr = (df["b_ours"] / df["b_pub"] - 1)[ring]
            out["ring_median_frac_diff"] = float(np.median(dr))
            out["ring_robust_scatter"] = float(1.4826 * np.median(np.abs(dr - np.median(dr))))
            out["n_ring"] = int(ring.sum())
    return out
