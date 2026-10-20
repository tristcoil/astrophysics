"""Cluster masses: from photons to solar masses, and Planck's hydrostatic bias.

A cluster's mass is not measured directly; it is inferred from its count rate
through the scaling relation. With an intrinsic scatter of 0.98 in ln CR (a
factor 2.7), a single count rate is a poor mass estimate on its own, and the
steep mass function matters: at a given count rate there are many more light
clusters scattered up than heavy ones scattered down (Eddington bias). The
posterior for each cluster is therefore

    P(M | CR_hat, z) ∝ dn/dlnM(M, z) * integral dlnCR  P(CR | M, z)
                       * P(detected | CR, z, sky) * P(CR_hat | CR, sky)

with that cluster's own exposure and absorption. The catalogue's M500 PDFs
(Bulbul et al. 2024) were made the same way with the same relation, so
reproducing them checks our pipeline cluster by cluster.

Planck's PSZ2 masses (MSZ) come from an SZ-mass relation calibrated on X-ray
masses that assume hydrostatic equilibrium. If the hot gas is partly held up by
turbulence and bulk motions, those masses come out low by a factor (1 - b).
eRASS1's masses are calibrated with weak lensing, which measures total mass
without that assumption, so MSZ / M_eRASS1 for the same clusters measures (1 - b).
"""

from __future__ import annotations

import numpy as np

import catalog as K
import config as C
import cosmology as CO
import model as MD
import selection as S


def cluster_mass_posteriors(df, cosmo: CO.Cosmology, sr: MD.Scaling | None = None,
                            selection: S.SelectionModel | None = None, meas_a: float = 2.52):
    """Posterior mass grid for every cluster: returns (lnM grid, P (n, nM))."""
    sr = MD.Scaling() if sr is None else sr
    selection = S.SelectionModel() if selection is None else selection
    lnm = np.log(np.geomspace(3e12, 1e16, 160))
    lncr = np.log(np.geomspace(1e-4, 300, 241))
    dlncr = np.gradient(lncr)
    zc = df.z.values
    zgrid = np.linspace(0.05, 0.85, 81)
    hmf_grid = cosmo.dndlnM(np.exp(lnm), zgrid)  # (z, M)
    mu_grid = sr.mean_lncr(cosmo, lnm, zgrid)
    out = np.empty((len(df), lnm.size))
    for i, (z, crhat, t, nh) in enumerate(zip(zc, df.cr.values, df.exp.values, df.nh.values)):
        j = np.clip(np.searchsorted(zgrid, z), 1, zgrid.size - 1)
        w = (z - zgrid[j - 1]) / (zgrid[j] - zgrid[j - 1])
        hmf = (1 - w) * hmf_grid[j - 1] + w * hmf_grid[j]
        mu = (1 - w) * mu_grid[j - 1] + w * mu_grid[j]
        pdet = selection(np.exp(lncr), z, t, nh)
        sig = MD.measurement_sigma(np.exp(lncr) * S.absorbed_fraction(nh) * t, meas_a)
        meas = np.exp(-0.5 * ((np.log(crhat) - lncr) / sig) ** 2) / sig
        kern = pdet * meas * dlncr  # over true CR
        g = np.exp(-0.5 * ((lncr[None, :] - mu[:, None]) / sr.sigma_X) ** 2)
        like = g @ kern
        p = hmf * like
        out[i] = p / np.trapezoid(p, lnm)
    return lnm, out


def summarize_posteriors(lnm, P):
    cdf = np.cumsum(P * np.gradient(lnm)[None, :], axis=1)
    cdf /= cdf[:, -1:]
    q = lambda f: np.exp(np.array([np.interp(f, c, lnm) for c in cdf]))  # noqa: E731
    return q(0.5), q(0.16), q(0.84)


def compare_with_catalogue(df, m_med) -> dict:
    ok = df.in_primary.values & np.isfinite(df.m500.values) & (df.m500.values > 0)
    r = m_med[ok] / df.m500.values[ok]
    return {"n": int(ok.sum()), "median_ratio": float(np.median(r)),
            "p16": float(np.percentile(r, 16)), "p84": float(np.percentile(r, 84)),
            "rms_dex": float(np.std(np.log10(r)))}


def hydrostatic_bias(primary, psz2, df=None, cosmo=None, our_masses=None,
                     snr_cuts=(4.5, 6.0, 7.0, 8.0, 10.0), n_boot=4000, seed=7) -> dict:
    """M_SZ(Planck) / M500(eRASS1) for matched clusters, four ways.

    1. the catalogue's M500 column (Bulbul et al.'s point estimates, after
       their correction that maps point estimates onto the expected mass
       distribution), median ratio;
    2. the median of the catalogue's own mass PDFs;
    3. our X-ray-only posteriors (wrong prior for SZ-selected clusters:
       Planck saw them *because* they are massive);
    4. "forward": predict each cluster's count rate from its Planck mass
       through the published relation and compare with the observed one,
       ln(1-b) = -<ln CR_hat - mu(M_SZ)> / B_X. X-ray selection is irrelevant
       at these masses, so no Eddington correction is involved on that side.
    The spread between them is the honest result: a mass proxy with a factor
    2.7 scatter cannot weigh individual clusters.
    """
    mt = K.match_psz2(primary, psz2)
    ok = (mt.sz_msz > 0) & np.isfinite(mt.x_m500) & (mt.x_z > 0.02)
    mt = mt[ok].reset_index(drop=True)
    Ma, pdf, det = K.primary_mass_pdfs()
    idx = {d: i for i, d in enumerate(det)}
    rows_pdf = np.array([idx[d] for d in mt.x_detuid])
    cdf = np.cumsum(pdf[rows_pdf], axis=1)
    cdf /= cdf[:, -1:]
    mt["x_pdf_median"] = [np.interp(0.5, c, m) for c, m in zip(cdf, Ma[rows_pdf])]
    est = {"catalogue M500 column": (mt.sz_msz / mt.x_m500).values,
           "catalogue PDF medians": (mt.sz_msz / mt.x_pdf_median).values}
    if df is not None and our_masses is not None:
        ours = dict(zip(df.name, our_masses))
        crs = dict(zip(df.name, df.cr))
        zs = dict(zip(df.name, df.z))
        r_ours = np.array([s / ours[n] if n in ours else np.nan for n, s in zip(mt.x_name, mt.sz_msz)])
        est["our X-ray-only posteriors"] = r_ours
        sr = MD.Scaling()
        fwd = np.full(len(mt), np.nan)
        for k, (n, msz) in enumerate(zip(mt.x_name, mt.sz_msz)):
            if n in crs:
                z = zs[n]
                mu = sr.mean_lncr(cosmo, np.array([np.log(msz)]), np.array([z]))[0, 0]
                slope = sr.B_X + sr.F_X * np.log((1 + z) / (1 + C.Z_PIVOT))
                fwd[k] = np.exp(-(np.log(crs[n]) - mu) / slope)
        est["forward: count rate predicted from Planck mass"] = fwd
    rng = np.random.default_rng(seed)
    table = []
    for name, ratio in est.items():
        for cut in snr_cuts:
            sel = (mt.sz_snr.values >= cut) & np.isfinite(ratio)
            r = np.log(ratio[sel])
            if name.startswith("forward"):
                centre = np.exp(np.mean(r))  # the estimator is a mean in ln CR
                boot = np.exp(np.mean(rng.choice(r, (n_boot, r.size)), axis=1))
            else:
                centre = np.exp(np.median(r))
                boot = np.exp(np.median(rng.choice(r, (n_boot, r.size)), axis=1))
            table.append({"estimator": name, "snr_min": cut, "n": int(sel.sum()), "one_minus_b": float(centre),
                          "err": float(np.std(boot)), "scatter_dex": float(np.std(r) / np.log(10))})
    lm = np.log10(mt.x_m500.values)
    slope = np.polyfit(lm - 14.7, np.log10(est["catalogue M500 column"]), 1)[0]
    return {"table": table, "n_matched": int(len(mt)), "slope_dlog_ratio_dlogM": float(slope),
            "matches": mt[["x_name", "x_detuid", "x_z", "x_m500", "x_m500_lo", "x_m500_hi", "x_pdf_median",
                           "sz_name", "sz_snr", "sz_z", "sz_msz", "sz_cosmo"]]}
