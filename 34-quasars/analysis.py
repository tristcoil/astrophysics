"""Every analysis step, in order, collecting the numbers the figures and README use.

    results = analysis.run()     # ~3-5 min once the data are cached

Sections:
 1. census       - who is in DR16Q, and which programme chose them
 2. composite    - our stacked spectrum vs Vanden Berk et al. 2001
 3. photometry   - colours vs redshift, the stellar locus, K-corrections
 4. qlf          - the luminosity function (Richards et al. 2006 reproduced; DR7 uniform)
 5. era          - the number density of luminous quasars through cosmic time
 6. masses       - black hole masses, Eddington ratios, three emission lines compared
 7. growth       - the z > 5.7 quasars and the seeds they need
"""

from __future__ import annotations

import json
import pickle
import time

import numpy as np
import pandas as pd

import catalogue
import composite
import fetch_data
import highz
import masses
import photometry
import physics
import qlf
from config import (COMPOSITE_GRID, CIV_MGII_Z, DERIVED, DR7_TABLE, DR7_UNIFORM_AREA_DEG2,
                    LUMINOUS_M_CUT, QLF_Z_EDGES, R06_AREA_DEG2, RESULTS_JSON, RESULTS_PICKLE,
                    RICHARDS06_TABLES, SHENKELLY12_TABLES, SPEC_LIST, Z_SPLIT_DR7, Z_SPLIT_R06)
from fetch_data import load_vizier_table

M_EDGES = np.round(np.arange(-31.2, -22.5 + 1e-9, 0.3), 2)
SK12_Z_EDGES = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.4, 2.9, 3.5, 4.0, 4.5, 5.0])
#: finer redshift bins for the density history
ERA_Z_EDGES = np.array([0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9,
                        3.2, 3.5, 3.8, 4.1, 4.5, 5.0])


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def robust_std(x) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def kcorr_r06() -> qlf.KCorrection:
    t4 = load_vizier_table(RICHARDS06_TABLES["kcorr"])
    return qlf.KCorrection(t4["z"], t4["KCorr"], "Richards et al. 2006")


# ---------------------------------------------------------------------------
# 1. Census
# ---------------------------------------------------------------------------

def census(cat: pd.DataFrame) -> dict:
    z = cat["zbest"].to_numpy()
    ok = z > 0
    edges = np.arange(0.0, 7.2 + 1e-9, 0.05)
    hist = {p: np.histogram(z[ok & (cat["programme"] == p).to_numpy()], edges)[0]
            for p in catalogue.PROGRAMMES}
    by_prog = (cat[ok].groupby("programme")["zbest"]
               .agg(n="size", median="median", q16=lambda s: s.quantile(0.16),
                    q84=lambda s: s.quantile(0.84)))
    return dict(
        n=len(cat), n_valid_z=int(ok.sum()), z_edges=edges, z_hist=hist,
        by_programme=by_prog.reset_index().to_dict("records"),
        n_z_gt=dict((t, int((z > t).sum())) for t in (2.0, 3.0, 4.0, 5.0, 6.0)),
        z_max=float(z.max()), z_median=float(np.median(z[ok])),
        redshift_source=cat["r_z"].value_counts().to_dict(),
        n_masses=int(np.isfinite(cat["m_fid"]).sum()),
        ra=cat["RA"].to_numpy(np.float32), dec=cat["DEC"].to_numpy(np.float32),
        z=z.astype(np.float32),
    )


# ---------------------------------------------------------------------------
# 2. Composite spectrum
# ---------------------------------------------------------------------------

def build_composite() -> dict:
    sample = pd.read_parquet(SPEC_LIST)
    paths = [fetch_data.spec_path(int(p), int(m), int(f))
             for p, m, f in zip(sample["PLATE"], sample["MJD"], sample["FIBERID"])]
    present = np.array([p.exists() for p in paths])
    sample, paths = sample[present].reset_index(drop=True), [p for p, ok in zip(paths, present) if ok]
    log(f"composite: stacking {len(sample):,} spectra")
    cache = DERIVED / "rest_frame_stack.npy"
    if cache.exists() and np.load(cache, mmap_mode="r").shape == (len(sample), len(COMPOSITE_GRID)):
        flux, bad = np.load(cache), []
    else:
        flux, _, bad = composite.load_stack(sample, paths)
        np.save(cache, flux)
    wave = COMPOSITE_GRID
    med = composite.median_composite(flux, wave, z=sample["Z_SYS"].to_numpy())
    gm = composite.geometric_mean_composite(flux, wave)
    a_gm, pl_gm = composite.power_law_fit(wave, gm["flux"], composite.VB01_UV_WINDOWS)
    a_gm_red, _ = composite.power_law_fit(wave, gm["flux"], composite.VB01_RED_WINDOWS)
    a_med, pl_med = composite.power_law_fit(wave, med["flux"], composite.VB01_UV_WINDOWS)
    a_med_red, _ = composite.power_law_fit(wave, med["flux"], composite.VB01_RED_WINDOWS)

    vb, _ = fetch_data.read_vandenberk()
    # put both median composites on one scale: unit mean over 3020-3100 A
    win = (wave >= 3020) & (wave <= 3100)
    ours = med["flux"] / np.nanmean(med["flux"][win])
    vbw = (vb["wave"] >= 3020) & (vb["wave"] <= 3100)
    vbf = vb["flux"].to_numpy() / vb["flux"][vbw].mean()
    vb_on_grid = np.interp(wave, vb["wave"], vbf, left=np.nan, right=np.nan)
    ratio = ours / vb_on_grid
    cmp = {}
    for lo, hi in ((1250, 2000), (2000, 3500), (3500, 5000), (5000, 7000), (1250, 7000)):
        m = (wave >= lo) & (wave <= hi) & np.isfinite(ratio)
        cmp[f"{lo}-{hi}"] = dict(median_ratio=float(np.median(ratio[m])),
                                 mad_pct=float(100 * np.median(np.abs(ratio[m] / np.median(ratio[m]) - 1))))
    # S/N of our median composite per 1 A pixel (flux / interquantile error)
    snr = ours / (med["err"] / np.nanmean(med["flux"][win]))
    n_used = int(np.isfinite(med["scale"]).sum())

    # The composite depends on whom you stack: the same machinery in bins of luminosity
    # (DR16Q's M_i(z=2)). Faint quasars show their host galaxies' starlight and
    # stronger lines (the Baldwin effect); the most luminous match Vanden Berk et al.
    mi = sample["M_I"].to_numpy()
    lum_bins = [(-30.0, -26.0), (-26.0, -25.0), (-25.0, -24.0), (-24.0, -23.0), (-23.0, -18.0)]
    lum = []
    for lo, hi in lum_bins:
        sel = (mi >= lo) & (mi < hi)
        mb = composite.median_composite(flux[sel], wave, z=sample["Z_SYS"].to_numpy()[sel], n_iter=4)
        f = mb["flux"]
        norm_win = (wave >= 1440) & (wave <= 1460)
        if not np.isfinite(f[norm_win]).any():
            norm_win = win
        fn = f / np.nanmean(f[norm_win])
        fn_ok = np.where(mb["n"] >= 20, fn, np.nan)
        civ_ew = line_ew(wave, fn_ok, (1500.0, 1600.0), ((1440.0, 1460.0), (1690.0, 1710.0)))
        o3_ew = line_ew(wave, fn_ok, (4995.0, 5020.0), ((4740.0, 4760.0), (5080.0, 5100.0)))
        lum.append(dict(m_lo=lo, m_hi=hi, n=int(sel.sum()), flux=fn, n_pix=mb["n"],
                        z_median=float(np.median(sample["Z_SYS"].to_numpy()[sel])),
                        civ_ew=civ_ew, oiii_ew=o3_ew))
    sel = mi < -25.0
    # A template for luminous quasars like those in the SDSS uniform sample: the
    # M_i < -25 median composite where >= 30 of them cover a pixel, the full
    # composite (scaled in 4000-4600 A) beyond. Used for colours and K-corrections.
    med_l = composite.median_composite(flux[sel], wave, z=sample["Z_SYS"].to_numpy()[sel], n_iter=4)
    lum_f = med_l["flux"] / np.nanmean(med_l["flux"][win])
    ov = (wave >= 4000) & (wave <= 4600) & (med_l["n"] >= 30) & np.isfinite(ours)
    tscale = float(np.median(lum_f[ov] / ours[ov]))
    template = np.where(med_l["n"] >= 30, lum_f, ours * tscale)
    a_tmpl, pl_tmpl = composite.power_law_fit(wave, template, composite.VB01_UV_WINDOWS)
    gm_lum = composite.geometric_mean_composite(flux[sel], wave)
    a_gm_lum, _ = composite.power_law_fit(wave, gm_lum["flux"], composite.VB01_UV_WINDOWS)
    try:
        a_med_lum, _ = composite.power_law_fit(wave, lum_f, composite.VB01_UV_WINDOWS)
    except (ValueError, IndexError):
        a_med_lum = np.nan
    return dict(
        wave=wave, median=ours, median_err=med["err"] / np.nanmean(med["flux"][win]),
        n_per_pixel=med["n"], gm=gm["flux"] / np.nanmean(gm["flux"][win]), gm_n=gm["n"],
        gm_n_spectra=gm["n_spectra"], gm_frac_nonpositive=gm["frac_nonpositive"],
        alpha_gm=a_gm, alpha_gm_red=a_gm_red, alpha_median=a_med, alpha_median_red=a_med_red,
        powerlaw_median=pl_med(wave) / np.nanmean(med["flux"][win]),
        powerlaw_gm=pl_gm(wave) / np.nanmean(gm["flux"][win]),
        scale_history=med["history"], n_spectra=len(sample), n_used=n_used, n_unreadable=len(bad),
        z=sample["Z_SYS"].to_numpy(), mi=sample["M_I"].to_numpy(),
        vb_wave=vb["wave"].to_numpy(), vb_flux=vbf, vb_on_grid=vb_on_grid, ratio_to_vb=ratio,
        vb_comparison=cmp, snr_2800=float(np.nanmedian(snr[(wave > 2790) & (wave < 2810)])),
        wave_min=float(wave[med["n"] >= 10].min()), wave_max=float(wave[med["n"] >= 10].max()),
        luminosity_bins=lum, alpha_gm_luminous=a_gm_lum, n_gm_luminous=gm_lum["n_spectra"],
        template=template, template_n=med_l["n"], alpha_template=a_tmpl,
        powerlaw_template=pl_tmpl(wave), template_scale=tscale,
        alpha_median_luminous=a_med_lum, n_luminous=int(sel.sum()),
        gm_coverage=dict(n_1350=int(np.isfinite(gm["flux"]) [(wave >= 1350) & (wave <= 1365)].sum()),
                         n_spectra_at_1350=int(np.median(gm["n"][(wave >= 1350) & (wave <= 1365)])),
                         n_spectra_at_4200=int(np.median(gm["n"][(wave >= 4200) & (wave <= 4230)]))),
    )


def line_ew(wave, flux, line, side_windows) -> float:
    """Rest-frame equivalent width (A) of an emission line above a linear continuum
    drawn through two side windows."""
    pts = []
    for lo, hi in side_windows:
        m = (wave >= lo) & (wave <= hi) & np.isfinite(flux)
        if m.sum() < 5:
            return np.nan
        pts.append((np.mean(wave[m]), np.median(flux[m])))
    (x1, y1), (x2, y2) = pts
    m = (wave >= line[0]) & (wave <= line[1]) & np.isfinite(flux)
    if m.sum() < 0.8 * (line[1] - line[0]):
        return np.nan
    cont = y1 + (y2 - y1) * (wave[m] - x1) / (x2 - x1)
    return float(np.sum((flux[m] - cont) / cont) * 1.0)


# ---------------------------------------------------------------------------
# 3. Photometry: colours, the stellar locus, K-corrections
# ---------------------------------------------------------------------------

#: line-free windows used to trace the total (power law + host galaxy) continuum
#: redward of 4200 A, where Richards et al. remove host light from the K-correction.
RED_CONTINUUM_WINDOWS = ((4200, 4230), (5470, 5520), (6005, 6035), (7160, 7180), (7900, 8000),
                         (8800, 8900))


def quasar_only_spectrum(comp: dict):
    """Power-law continuum + emission lines, host-galaxy light removed (as Richards et al. did).

    Blueward of 4230 A the continuum is the alpha_nu power law through Vanden Berk's
    windows. Redward, the observed continuum rises above it (host starlight);
    we trace it through line-free windows and keep only the lines on top.
    """
    w, f = comp["wave"], comp["template"]
    pl = comp["powerlaw_template"]
    cont = pl.copy()
    pts_w, pts_r = [], []
    for lo, hi in RED_CONTINUUM_WINDOWS:
        m = (w >= lo) & (w <= hi) & np.isfinite(f)
        if m.sum() > 5:
            pts_w.append(np.mean(w[m]))
            pts_r.append(np.mean(f[m]) / np.mean(pl[m]))
    red = w > 4215
    excess = np.interp(w[red], pts_w, pts_r)   # total continuum / power law
    cont[red] = pl[red] * excess
    lines = f - cont
    qso = pl + np.where(np.isfinite(lines), lines, 0.0)
    return cont, qso, dict(window_waves=pts_w, excess_ratio=pts_r)


def photometry_analysis(comp: dict, cat: pd.DataFrame, stars: pd.DataFrame, k06) -> dict:
    w = comp["wave"]
    cont, qso, host = quasar_only_spectrum(comp)
    # an intrinsic spectrum to shift to any redshift: forest-free below Ly alpha
    good = np.isfinite(comp["template"])
    base = np.where(good, comp["template"], comp["powerlaw_template"])
    intrinsic = photometry.intrinsic_blue_extension(w, base, comp["alpha_template"], 1450.0)
    zgrid = np.round(np.arange(0.10, 5.50 + 1e-9, 0.02), 2)
    track = photometry.colour_track(w, intrinsic, zgrid, igm=True)
    track_noigm = photometry.colour_track(w, intrinsic, zgrid, igm=False)

    # observed colours of DR16Q quasars bright enough for clean photometry
    sel = (cat["zbest"] > 0.1) & (cat["i0"] < 19.5) & (cat["i0"] > 15)
    obs = cat.loc[sel, ["zbest", "u-g", "g-r", "r-i", "i-z"]]
    zb = np.arange(0.1, 5.0 + 1e-9, 0.05)
    idx = np.digitize(obs["zbest"], zb) - 1
    obs_med = {}
    for c in ("u-g", "g-r", "r-i", "i-z"):
        vals = obs[c].to_numpy()
        med, lo, hi, n = [], [], [], []
        for k in range(len(zb) - 1):
            v = vals[idx == k]
            v = v[np.isfinite(v)]
            n.append(len(v))
            if len(v) >= 20:
                q = np.percentile(v, [16, 50, 84])
                lo.append(q[0]); med.append(q[1]); hi.append(q[2])
            else:
                lo.append(np.nan); med.append(np.nan); hi.append(np.nan)
        obs_med[c] = dict(z=0.5 * (zb[1:] + zb[:-1]), median=np.array(med), p16=np.array(lo),
                          p84=np.array(hi), n=np.array(n))
    agreement = {}
    for c in ("u-g", "g-r", "r-i", "i-z"):
        syn = np.interp(obs_med[c]["z"], zgrid, track[c])
        d = syn - obs_med[c]["median"]
        for lo, hi in ((0.3, 2.2), (2.2, 3.5), (3.5, 4.8)):
            m = (obs_med[c]["z"] >= lo) & (obs_med[c]["z"] < hi) & np.isfinite(d)
            agreement[f"{c} {lo}-{hi}"] = float(np.median(np.abs(d[m]))) if m.any() else np.nan

    # distance of the quasar colour track from the stellar locus (4D colour space)
    sc = stars[["u", "g", "r", "i", "z"]].to_numpy()
    scol = np.column_stack([sc[:, 0] - sc[:, 1], sc[:, 1] - sc[:, 2], sc[:, 2] - sc[:, 3],
                            sc[:, 3] - sc[:, 4]])
    qcol = np.column_stack([track["u-g"], track["g-r"], track["r-i"], track["i-z"]])
    from scipy.spatial import cKDTree
    def locus_distance(cols):
        """Distance from the quasar track to the 20th-nearest star in the given colours
        (robust to stray stars: "how deep inside the stellar locus")."""
        tree = cKDTree(scol[:, cols])
        q = qcol[:, cols]
        finite = np.all(np.isfinite(q), axis=1)
        out = np.full(len(zgrid), np.nan)
        dd, _ = tree.query(q[finite], k=20)
        out[finite] = dd[:, -1]
        return out

    dist = locus_distance([0, 1, 2, 3])
    # The SDSS selection used ugri colours for its low-redshift branch and griz
    # colours (where u is too faint to measure) for the high-redshift branch.
    dist_ugri = locus_distance([0, 1, 2])
    dist_griz = locus_distance([1, 2, 3])

    sf = qlf.SelectionFunction(load_vizier_table(RICHARDS06_TABLES["selfunc"]))
    sf_marg = sf.marginal(zgrid, 19.1)
    sf_marg_faint = sf.marginal(zgrid, 20.2)

    # K-corrections: ours (alpha = -0.5 continuum + our emission-line term) vs Richards et al.
    kz = np.round(np.arange(0.0, 5.49 + 1e-9, 0.01), 2)
    k_em = photometry.kcorr_emission(w, qso, comp["powerlaw_template"], kz, "i")
    k_ours = photometry.kcorr_continuum(kz) + k_em
    k_r06 = k06(kz)
    dk = k_ours - k_r06
    kcmp = {}
    for lo, hi in ((0.0, 0.5), (0.5, 2.2), (2.2, 4.0), (4.0, 5.1)):
        m = (kz >= lo) & (kz < hi) & np.isfinite(dk)
        kcmp[f"{lo}-{hi}"] = dict(median=float(np.median(dk[m])), max_abs=float(np.max(np.abs(dk[m]))))

    # colours needed to turn g/r survey limits into i limits
    gi = track["g"] - track["i"]
    ri = track["r"] - track["i"]
    return dict(zgrid=zgrid, track=track, track_noigm=track_noigm, observed=obs_med,
                agreement=agreement, locus_distance=dist, locus_distance_ugri=dist_ugri,
                locus_distance_griz=dist_griz, selfunc_marginal=sf_marg,
                selfunc_marginal_faint=sf_marg_faint, kz=kz, k_ours=k_ours, k_r06=k_r06,
                k_em=k_em, k_comparison=kcmp, g_minus_i=gi, r_minus_i=ri,
                intrinsic=intrinsic, continuum=cont, quasar_only=qso, host=host,
                star_colours=scol.astype(np.float32))


# ---------------------------------------------------------------------------
# 4-5. Luminosity function and the quasar era
# ---------------------------------------------------------------------------

def _compare(ours: pd.DataFrame, theirs: pd.DataFrame, z_edges, zcol, mcol, logcol, errcol=None):
    zc = 0.5 * (z_edges[1:] + z_edges[:-1])
    th = theirs.copy()
    th["zmid"] = [zc[np.argmin(np.abs(zc - v))] for v in th[zcol]]
    m = th.merge(ours, left_on=[th["zmid"].round(3), th[mcol].round(2)],
                 right_on=[ours["z_mid"].round(3), ours["m_mid"].round(2)], how="left")
    d = np.log10(m["phi"]) - m[logcol]
    good = np.isfinite(d)
    out = dict(n_bins=int(len(th)), n_matched=int(good.sum()),
               median_dlog=float(np.median(d[good])), robust_std_dlog=robust_std(d[good]),
               frac_within_0p05=float(np.mean(np.abs(d[good]) < 0.05)),
               frac_within_0p1=float(np.mean(np.abs(d[good]) < 0.1)))
    return out, m.assign(dlog=d)


def qlf_analysis(k06) -> dict:
    t5 = load_vizier_table(RICHARDS06_TABLES["sample"])
    t6 = load_vizier_table(RICHARDS06_TABLES["qlf"])
    sf = qlf.SelectionFunction(load_vizier_table(RICHARDS06_TABLES["selfunc"]))

    # (a) Richards et al. 2006 rebuilt from their own sample: their weights (Cor already
    # includes the 95% image-quality factor), their floor, their stated flux limits.
    keep = ~((t5["z"] < Z_SPLIT_R06) & (t5["imag"] > 19.1))
    s5 = t5[keep]
    sv06 = qlf.Survey(area_deg2=R06_AREA_DEG2, kcorr=k06, z_split=Z_SPLIT_R06, cosmology="ref",
                      name="SDSS DR3 (Richards et al. 2006)")
    w5 = 1.0 / np.maximum(s5["Cor"].to_numpy(), 1.0 / 3.0)
    q06 = qlf.binned_qlf(s5["z"], s5["iMAG"], w5, sv06, QLF_Z_EDGES, M_EDGES)
    rep_stats, rep_table = _compare(q06, t6, QLF_Z_EDGES, "z", "iMAG", "logPhi")
    rep_stats["n_excluded_faint_lowz"] = int((~keep).sum())

    # (b) the DR7 uniform sample (Shen et al. 2011 flag 1; Shen & Kelly 2012's 57,959)
    d7 = pd.read_parquet(DR7_TABLE)
    u = d7[(d7["Uni"] == 1) & (d7["z"] >= 0.3) & (d7["z"] <= 5.0)].copy()
    kind = np.where(u["FIRST"] >= 1, "radio", np.where(u["dr7q_Mph"] == 1, "extend", "point"))
    p = sf(u["imag0"].to_numpy(), u["z"].to_numpy(), kind)
    u["selfunc"] = p
    u["weight"] = qlf.weights_from_selection(p)
    sv7 = qlf.Survey(area_deg2=DR7_UNIFORM_AREA_DEG2, kcorr=k06, z_split=Z_SPLIT_DR7,
                     cosmology="ref", name="SDSS DR7 uniform")
    q7 = qlf.binned_qlf(u["z"], u["iMAG"], u["weight"], sv7, QLF_Z_EDGES, M_EDGES)
    q7_sk = qlf.binned_qlf(u["z"], u["iMAG"], u["weight"], sv7, SK12_Z_EDGES, M_EDGES)
    sk3 = load_vizier_table(SHENKELLY12_TABLES["qlf"])
    sk_stats, sk_table = _compare(q7_sk, sk3, SK12_Z_EDGES, "zbin", "magz2", "logPhi")
    q7_sk_noimg = qlf.binned_qlf(u["z"], u["iMAG"], qlf.weights_from_selection(p, image_completeness=1.0),
                                 sv7, SK12_Z_EDGES, M_EDGES)
    sk_stats_noimg, sk_table_noimg = _compare(q7_sk_noimg, sk3, SK12_Z_EDGES, "zbin", "magz2", "logPhi")
    dr3_dr7_stats, dr3_dr7_table = _compare(q7, t6, QLF_Z_EDGES, "z", "iMAG", "logPhi")

    # (c) the quasar era in D3's cosmology: luminous quasars per Gpc^3 vs time
    u["Mi_d3"] = u["iMAG"] - (physics.distance_modulus(u["z"], "d3") - physics.distance_modulus(u["z"], "ref"))
    sv7d3 = qlf.Survey(area_deg2=DR7_UNIFORM_AREA_DEG2, kcorr=k06, z_split=Z_SPLIT_DR7,
                       cosmology="d3", name="SDSS DR7 uniform (D3 cosmology)")
    qera = qlf.binned_qlf(u["z"], u["Mi_d3"], u["weight"], sv7d3, ERA_Z_EDGES, M_EDGES)
    era = {}
    for cut in (-28.5, -27.6, -26.7):
        dens = qlf.cumulative_density(qera, cut)
        cb = complete_bins(sv7d3, cut)
        dens["complete"] = [any(abs(a - zl) < 1e-9 for a, _ in cb) for zl in dens["z_lo"]]
        era[cut] = dens
    peak = era_peak(u, sv7d3, LUMINOUS_M_CUT)
    # how much the peak moves with the choices we made
    variants = {}
    for name, floor in (("floor 0.2", 0.2), ("floor 0.5", 0.5), ("no floor", 1e-3)):
        uu = u.assign(weight=qlf.weights_from_selection(u["selfunc"], floor=floor))
        variants[name] = era_peak(uu, sv7d3, LUMINOUS_M_CUT, n_boot=0)["z_peak"]
    variants["no selection-function correction"] = era_peak(
        u.assign(weight=1.0 / 0.95), sv7d3, LUMINOUS_M_CUT, n_boot=0)["z_peak"]
    variants["parabola in z"] = era_peak(u, sv7d3, LUMINOUS_M_CUT, n_boot=0, in_log=False)["z_peak"]
    variants["z < 4.1 only"] = era_peak(u, sv7d3, LUMINOUS_M_CUT, n_boot=0, z_max=4.1)["z_peak"]
    variants["z > 1.9 only"] = era_peak(u, sv7d3, LUMINOUS_M_CUT, n_boot=0, z_min=1.9)["z_peak"]
    for cut in (-28.2, -27.0):
        variants[f"M_i < {cut}"] = era_peak(u, sv7d3, cut, n_boot=0)["z_peak"]
    peak["variants"] = variants
    # the same density history from Richards et al.'s published bins
    r06_dens = qlf.cumulative_density(
        t6.assign(z_lo=lambda d: [QLF_Z_EDGES[np.argmin(np.abs(0.5 * (QLF_Z_EDGES[1:] + QLF_Z_EDGES[:-1]) - v))] for v in d["z"]],
                  z_hi=lambda d: [QLF_Z_EDGES[np.argmin(np.abs(0.5 * (QLF_Z_EDGES[1:] + QLF_Z_EDGES[:-1]) - v)) + 1] for v in d["z"]],
                  m_lo=lambda d: d["iMAG"] - 0.15, m_hi=lambda d: d["iMAG"] + 0.15,
                  phi=lambda d: 10 ** d["logPhi"], phi_err=lambda d: d["errPhi"] * 1e-9,
                  z_mean=lambda d: d["<z>"], n=lambda d: d["NQ"], fill=lambda d: d["Fill"].astype(float)),
        LUMINOUS_M_CUT)
    r06_dens_d3 = r06_to_d3(t6, LUMINOUS_M_CUT)
    return dict(r06_reproduction=rep_stats, r06_table=rep_table, q06=q06, r06_luminous_density_d3=r06_dens_d3,
                dr7=q7, dr7_sk12_bins=q7_sk, sk12_comparison=sk_stats, sk12_table=sk_table,
                sk12_comparison_without_image_factor=sk_stats_noimg, sk12_table_noimg=sk_table_noimg,
                dr3_vs_dr7=dr3_dr7_stats, dr3_dr7_table=dr3_dr7_table,
                n_dr7_uniform=len(u), dr7_kind_counts=pd.Series(kind).value_counts().to_dict(),
                dr7_floor_frac=float(np.mean(p < 1.0 / 3.0)),
                dr7_weights_sum=float(u["weight"].sum()),
                era=era, era_peak=peak, r06_luminous_density=r06_dens,
                dr7_sample=u[["z", "iMAG", "Mi_d3", "imag0", "selfunc", "weight"]].reset_index(drop=True))


def r06_to_d3(t6: pd.DataFrame, cut: float) -> pd.DataFrame:
    """Richards et al.'s density of quasars brighter than ``cut``, moved into D3's cosmology.

    Their bins are in (h, Om) = (0.7, 0.3). In D3's cosmology (H0 = 66.66) every
    quasar sits a little farther away and is Delta = DM_D3 - DM_ref ~ 0.1 mag
    more luminous, so M_D3 < cut means M_ref < cut + Delta: the bins brighter
    than the cut plus a Delta/0.3 share of the next one. The density then
    scales by the ratio of comoving volume elements, dV_ref/dV_D3.
    """
    zc = 0.5 * (QLF_Z_EDGES[1:] + QLF_Z_EDGES[:-1])
    rows = []
    for zm, g in t6.groupby("z"):
        k = int(np.argmin(np.abs(zc - zm)))
        zmean = float(np.average(g["<z>"], weights=g["NQ"]))
        delta = float(physics.distance_modulus(zmean, "d3") - physics.distance_modulus(zmean, "ref"))
        phi = 10 ** g["logPhi"].to_numpy()
        err = g["errPhi"].to_numpy() * 1e-9
        m = g["iMAG"].to_numpy()
        full = m + 0.15 <= cut + 1e-6
        part = np.isclose(m - 0.15, cut)
        dens = np.sum(phi[full] * 0.3) + np.sum(phi[part] * delta)
        var = np.sum((err[full] * 0.3) ** 2) + np.sum((err[part] * delta) ** 2)
        vol = float(physics.dvdz_per_deg2(zmean, "ref") / physics.dvdz_per_deg2(zmean, "d3"))
        # a real measurement needs at least one whole bin brighter than the cut
        complete = bool(full.any() and part.any() and np.all(g.loc[full | part, "Fill"] == 1))
        rows.append(dict(z_lo=QLF_Z_EDGES[k], z_hi=QLF_Z_EDGES[k + 1], z_mean=zmean, delta_mag=delta,
                         volume_ratio=vol, density=dens * vol, err=np.sqrt(var) * vol, complete=complete))
    return pd.DataFrame(rows)


def era_peak(u: pd.DataFrame, survey: qlf.Survey, cut: float, n_boot: int = 300,
             seed: int = 7, in_log: bool = True, z_min: float = 0.0, z_max: float = 9.0) -> dict:
    """Redshift of peak luminous-quasar density: a parabola in log n vs log(1+z),
    fitted to the complete bins, with a bootstrap over quasars."""
    medges = M_EDGES[M_EDGES <= cut + 1e-9]
    # the luminosity bins must reach the cut exactly
    if abs(medges[-1] - cut) > 1e-6:
        medges = np.append(medges, cut)
    keep_bins = complete_bins(survey, cut)
    keep_bins = [(a, b) for a, b in keep_bins if a >= z_min - 1e-9 and b <= z_max + 1e-9]

    def peak_of(df):
        q = qlf.binned_qlf(df["z"], df["Mi_d3"], df["weight"], survey, ERA_Z_EDGES, medges)
        d = qlf.cumulative_density(q, cut)
        d = d[(d["density"] > 0) & [any(abs(a - zl) < 1e-9 for a, _ in keep_bins) for zl in d["z_lo"]]]
        x = np.log10(1 + d["z_mean"].to_numpy()) if in_log else d["z_mean"].to_numpy()
        y = np.log10(d["density"].to_numpy())
        sig = d["err"].to_numpy() / d["density"].to_numpy() / np.log(10)
        c = np.polyfit(x, y, 2, w=1 / sig)
        if c[0] >= 0:
            return np.nan, np.nan, d
        xp = -c[1] / (2 * c[0])
        zp = 10 ** xp - 1 if in_log else xp
        return float(zp), float(10 ** np.polyval(c, xp)), d

    zpk, npk, d = peak_of(u)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        bs = u.iloc[rng.integers(0, len(u), len(u))]
        boots.append(peak_of(bs)[0])
    boots = np.array(boots)
    boots = boots[np.isfinite(boots)]
    if n_boot == 0:
        return dict(z_peak=zpk, density_peak_per_gpc3=npk * 1e9, bins=d)
    return dict(z_peak=zpk, z_peak_p16=float(np.percentile(boots, 16)),
                z_peak_p84=float(np.percentile(boots, 84)), density_peak_per_gpc3=npk * 1e9,
                age_at_peak_gyr=float(physics.age_myr(zpk) / 1e3),
                age_p16=float(physics.age_myr(np.percentile(boots, 84)) / 1e3),
                age_p84=float(physics.age_myr(np.percentile(boots, 16)) / 1e3),
                z_range=(float(keep_bins[0][0]), float(keep_bins[-1][1])), bins=d, n_boot=len(boots))


def complete_bins(survey: qlf.Survey, cut: float, m_bright: float = -29.6):
    """Redshift bins of ERA_Z_EDGES in which every quasar with m_bright < M < cut would
    have been inside the flux limits (the bright limit i = 15 never bites above
    M ~ -29.6 at the redshifts used: brighter quasars hold < 2% of the density)."""
    out = []
    for a, b in zip(ERA_Z_EDGES[:-1], ERA_Z_EDGES[1:]):
        _, fill = survey.accessible_volume(a, b, m_bright, cut, nz=200)
        if fill > 0.999:
            out.append((float(a), float(b)))
    return out


# ---------------------------------------------------------------------------
# 6. Masses and Eddington ratios
# ---------------------------------------------------------------------------

def limit_lbol(survey_mag: float, z, colour=None, k06=None):
    """L_bol of a quasar exactly at an apparent-magnitude limit, via M_i(z=2) -> L_2500.

    ``colour`` is (X - i) at each z when the limit is in another band X. L_bol =
    5.15 x lambda L_3000 with lambda L_3000 = 0.913 lambda L_2500 for alpha_nu = -0.5.
    """
    i_lim = survey_mag - (0.0 if colour is None else colour)
    mi = physics.absolute_mag(i_lim, z, k06, "d3")
    l2500 = physics.luminosity_2500(mi)  # erg/s/Hz
    lam_l = l2500 + np.log10(2.99792458e18 / 2500.0)
    return lam_l + np.log10(0.913) + np.log10(5.15)


def mass_analysis(cat: pd.DataFrame, phot: dict, k06) -> dict:
    ok = np.isfinite(cat["m_fid"]) & np.isfinite(cat["log_lbol"]) & (cat["zbest"] > 0)
    s = cat.loc[ok, ["zbest", "m_fid", "log_lbol", "log_edd", "line_fid", "programme", "m_hb",
                     "m_mgii", "m_mgii_vo09", "m_civ", "m_civ_c17", "civ_blueshift", "IF_BOSS_SDSS",
                     "fwhm_civ", "fwhm_mgii", "log_l1350", "LOGL3000"]]
    zb = np.arange(0.0, 5.0 + 1e-9, 0.25)
    idx = np.digitize(s["zbest"], zb) - 1
    per_z = []
    for k in range(len(zb) - 1):
        g = s[idx == k]
        if len(g) < 50:
            continue
        per_z.append(dict(z_lo=zb[k], z_hi=zb[k + 1], n=len(g),
                          m_p5=float(g["m_fid"].quantile(0.05)), m_p50=float(g["m_fid"].median()),
                          m_p95=float(g["m_fid"].quantile(0.95)), m_max=float(g["m_fid"].max()),
                          edd_p50=float(g["log_edd"].median()), edd_std=robust_std(g["log_edd"]),
                          frac_super_edd=float(np.mean(g["log_edd"] > 0))))
    per_z = pd.DataFrame(per_z)

    # Kollmeier et al.'s regime: 0.3 < z < 4, 1e45 < L_bol < 1e47
    k = s[(s["zbest"] > 0.3) & (s["zbest"] < 4) & (s["log_lbol"] > 45) & (s["log_lbol"] < 47)]
    e = k["log_edd"].to_numpy()
    kollmeier = dict(n=len(k), median=float(np.median(e)), mean=float(np.mean(e)),
                     std=float(np.std(e)), robust_std=robust_std(e),
                     frac_above_1=float(np.mean(e > 0)), frac_above_0p3=float(np.mean(e > np.log10(0.3))),
                     frac_below_0p01=float(np.mean(e < -2)))
    # how much of the width the mass uncertainty alone explains
    kollmeier["width_after_0p3dex_mass_noise"] = float(np.sqrt(max(kollmeier["robust_std"] ** 2 - 0.3 ** 2, 0)))
    from scipy.stats import norm
    # a lognormal with the observed centre and width would put this much above Eddington...
    kollmeier["gauss_frac_above_1"] = float(norm.sf(-kollmeier["median"] / kollmeier["robust_std"]))
    # ...and with the 0.3 dex of mass error taken out of the width, this much
    w0 = kollmeier["width_after_0p3dex_mass_noise"]
    kollmeier["intrinsic_frac_above_1"] = float(norm.sf(-kollmeier["median"] / w0)) if w0 > 0 else 0.0

    # the three rulers: masses from two lines in the same spectrum
    def pair(a, b, zlo, zhi, boss_only=True):
        m = np.isfinite(s[a]) & np.isfinite(s[b]) & (s["zbest"] > zlo) & (s["zbest"] < zhi)
        if boss_only:
            m &= s["IF_BOSS_SDSS"] == "BOSS"
        d = (s.loc[m, a] - s.loc[m, b]).to_numpy()
        return dict(n=int(m.sum()), median=float(np.median(d)), robust_std=robust_std(d),
                    a=s.loc[m, a].to_numpy(np.float32), b=s.loc[m, b].to_numpy(np.float32))

    hb_mg = pair("m_mgii", "m_hb", 0.4, 0.9)
    mg_c4 = pair("m_civ", "m_mgii", CIV_MGII_Z[0], CIV_MGII_Z[1])
    # bias of C IV vs Mg II as a function of the C IV blueshift, before/after Coatman et al.
    m = (np.isfinite(s["civ_blueshift"]) & np.isfinite(s["m_mgii"]) & np.isfinite(s["m_civ"])
         & (s["zbest"] > CIV_MGII_Z[0]) & (s["zbest"] < CIV_MGII_Z[1]))
    bs = s.loc[m, "civ_blueshift"].to_numpy()
    d_raw = (s.loc[m, "m_civ"] - s.loc[m, "m_mgii"]).to_numpy()
    d_c17 = (s.loc[m, "m_civ_c17"] - s.loc[m, "m_mgii"]).to_numpy()
    edges = np.arange(-1000, 4001, 250)
    bi = np.digitize(bs, edges) - 1
    rows = []
    for j in range(len(edges) - 1):
        sel = bi == j
        if sel.sum() < 30:
            continue
        rows.append(dict(bs_lo=edges[j], bs_hi=edges[j + 1], bs_mid=0.5 * (edges[j] + edges[j + 1]),
                         n=int(sel.sum()), raw=float(np.median(d_raw[sel])), c17=float(np.median(d_c17[sel])),
                         raw_q16=float(np.percentile(d_raw[sel], 16)), raw_q84=float(np.percentile(d_raw[sel], 84)),
                         c17_q16=float(np.percentile(d_c17[sel], 16)), c17_q84=float(np.percentile(d_c17[sel], 84))))
    fin = np.isfinite(d_c17)
    # our own recalibration of C IV onto Mg II, two-fold cross-validated: fit on one
    # random half, judge on the other, then swap
    sub = s.loc[m]
    rng = np.random.default_rng(11)
    half = rng.random(len(sub)) < 0.5
    d_ours = np.full(len(sub), np.nan)
    cals = []
    for train in (half, ~half):
        cal = masses.fit_civ_recalibration(sub["civ_blueshift"].to_numpy()[train],
                                           sub["fwhm_civ"].to_numpy()[train], sub["fwhm_mgii"].to_numpy()[train],
                                           sub["log_l1350"].to_numpy()[train], sub["m_mgii"].to_numpy()[train])
        cals.append(cal)
        test = ~train
        d_ours[test] = masses.recalibrated_civ_mass(sub["civ_blueshift"].to_numpy()[test],
                                                    sub["fwhm_civ"].to_numpy()[test],
                                                    sub["log_l1350"].to_numpy()[test], cal) - sub["m_mgii"].to_numpy()[test]
    cal_all = masses.fit_civ_recalibration(sub["civ_blueshift"].to_numpy(), sub["fwhm_civ"].to_numpy(),
                                           sub["fwhm_mgii"].to_numpy(), sub["log_l1350"].to_numpy(),
                                           sub["m_mgii"].to_numpy())
    fo = np.isfinite(d_ours)
    for row in rows:
        sel = (bs >= row["bs_lo"]) & (bs < row["bs_hi"]) & fo
        row["ours"] = float(np.median(d_ours[sel])) if sel.sum() >= 30 else np.nan
    # the width ratio itself, against Coatman's C IV / H alpha line
    ratio = (sub["fwhm_civ"] / sub["fwhm_mgii"]).to_numpy()
    ratio_bins = []
    for lo, hi in ((-1000, 0), (0, 500), (500, 1000), (1000, 1500), (1500, 2000), (2000, 3000), (3000, 4000)):
        sel = (bs >= lo) & (bs < hi)
        if sel.sum() >= 30:
            ratio_bins.append(dict(lo=lo, hi=hi, mid=float(np.median(bs[sel])), n=int(sel.sum()),
                                   ratio=float(np.median(ratio[sel])),
                                   coatman=float(0.41 * np.median(bs[sel]) / 1e3 + 0.62)))
    # the blueshift trend at fixed luminosity: is it the outflow, or just luminosity?
    fixed_l = []
    for lo, hi in ((44.5, 45.0), (45.0, 45.5), (45.5, 46.0), (46.0, 46.5)):
        sel = (sub["log_l1350"].to_numpy() >= lo) & (sub["log_l1350"].to_numpy() < hi)
        if sel.sum() > 200:
            fixed_l.append(dict(l_lo=lo, l_hi=hi, n=int(sel.sum()),
                                slope=float(np.polyfit(bs[sel] / 1e3, d_raw[sel], 1)[0]),
                                bs_median=float(np.median(bs[sel]))))
    blueshift = dict(n=int(m.sum()), bins=pd.DataFrame(rows), bs=bs.astype(np.float32),
                     d_raw=d_raw.astype(np.float32), d_c17=d_c17.astype(np.float32),
                     d_ours=d_ours.astype(np.float32),
                     corr_raw=float(np.corrcoef(bs, d_raw)[0, 1]),
                     corr_c17=float(np.corrcoef(bs[fin], d_c17[fin])[0, 1]),
                     corr_ours=float(np.corrcoef(bs[fo], d_ours[fo])[0, 1]),
                     slope_raw=float(np.polyfit(bs / 1000, d_raw, 1)[0]),
                     slope_c17=float(np.polyfit(bs[fin] / 1000, d_c17[fin], 1)[0]),
                     slope_ours=float(np.polyfit(bs[fo] / 1000, d_ours[fo], 1)[0]),
                     median_raw=float(np.median(d_raw)), median_c17=float(np.median(d_c17[fin])),
                     median_ours=float(np.median(d_ours[fo])),
                     std_raw=robust_std(d_raw), std_c17=robust_std(d_c17[fin]), std_ours=robust_std(d_ours[fo]),
                     n_c17=int(fin.sum()), calibration=cal_all, calibration_halves=cals,
                     ratio_bins=pd.DataFrame(ratio_bins), fixed_luminosity=pd.DataFrame(fixed_l),
                     bs_median=float(np.median(bs)), bs_p90=float(np.percentile(bs, 90)))

    # recipe offsets vs luminosity: the recipes were matched on luminous quasars
    mm = (np.isfinite(s["m_hb"]) & np.isfinite(s["m_mgii"]) & (s["zbest"] > 0.4) & (s["zbest"] < 0.9)
          & (s["IF_BOSS_SDSS"] == "BOSS"))
    l3 = s.loc[mm, "LOGL3000"].to_numpy()
    rec = []
    for lo, hi in ((43.5, 44.0), (44.0, 44.5), (44.5, 45.0), (45.0, 45.5), (45.5, 46.5)):
        sel = (l3 >= lo) & (l3 < hi)
        if sel.sum() >= 50:
            rec.append(dict(l_lo=lo, l_hi=hi, l_mid=float(np.median(l3[sel])), n=int(sel.sum()),
                            s11_minus_hb=float(np.median((s.loc[mm, "m_mgii"] - s.loc[mm, "m_hb"]).to_numpy()[sel])),
                            vo09_minus_hb=float(np.median((s.loc[mm, "m_mgii_vo09"] - s.loc[mm, "m_hb"]).to_numpy()[sel]))))
    l1 = s.loc[m, "log_l1350"].to_numpy()
    rec_c4 = []
    for lo, hi in ((44.5, 45.0), (45.0, 45.5), (45.5, 46.0), (46.0, 46.5), (46.5, 47.5)):
        sel = (l1 >= lo) & (l1 < hi)
        if sel.sum() >= 50:
            rec_c4.append(dict(l_lo=lo, l_hi=hi, l_mid=float(np.median(l1[sel])), n=int(sel.sum()),
                               civ_minus_s11=float(np.median(d_raw[sel])),
                               bs_median=float(np.median(bs[sel]))))
    recipes = dict(mg_hb=pd.DataFrame(rec), c4_mg=pd.DataFrame(rec_c4),
                   vo09_minus_hb_all=float(np.median((s.loc[mm, "m_mgii_vo09"] - s.loc[mm, "m_hb"]).to_numpy())))

    # the step at z = 2, where the fiducial mass switches from Mg II to C IV
    zs = np.round(np.arange(1.6, 2.45, 0.1), 2)
    step = []
    allm = s[(s["IF_BOSS_SDSS"] == "BOSS")]
    rc = masses.recalibrated_civ_mass(allm["civ_blueshift"], allm["fwhm_civ"], allm["log_l1350"], cal_all)
    for lo in zs[:-1]:
        sel = ((allm["zbest"] >= lo) & (allm["zbest"] < lo + 0.1)).to_numpy()
        both = sel & np.isfinite(allm["m_mgii"]).to_numpy() & np.isfinite(allm["m_civ"]).to_numpy() & np.isfinite(rc)
        step.append(dict(z_mid=float(lo + 0.05), n=int(both.sum()),
                         fiducial=float(np.median(allm["m_fid"].to_numpy()[both])),
                         mgii=float(np.median(allm["m_mgii"].to_numpy()[both])),
                         civ=float(np.median(allm["m_civ"].to_numpy()[both])),
                         civ_recal=float(np.median(rc[both]))))
    step = pd.DataFrame(step)

    # detection limits drawn in the mass-redshift plane
    zl = np.linspace(0.1, 5.0, 200)
    lb_sdss = np.where(zl < 2.9, limit_lbol(19.1, zl, k06=k06), limit_lbol(20.2, zl, k06=k06))
    g_i = np.interp(zl, phot["zgrid"], phot["g_minus_i"])
    r_i = np.interp(zl, phot["zgrid"], phot["r_minus_i"])
    # BOSS/eBOSS: g < 22 or r < 21.85 - the fainter of the two limits applies
    lb_eboss = np.minimum(limit_lbol(22.0, zl, colour=g_i, k06=k06),
                          limit_lbol(21.85, zl, colour=r_i, k06=k06))
    led = np.log10(physics.eddington_per_msun())
    limits = dict(z=zl, lbol_sdss=lb_sdss, lbol_eboss=lb_eboss,
                  m_min_sdss=lb_sdss - led, m_min_eboss=lb_eboss - led)

    return dict(n=len(s), per_z=per_z, kollmeier=kollmeier, hb_mg=hb_mg, mg_c4=mg_c4,
                blueshift=blueshift, limits=limits, recipes=recipes, z2_step=step,
                z=s["zbest"].to_numpy(np.float32), m=s["m_fid"].to_numpy(np.float32),
                lbol=s["log_lbol"].to_numpy(np.float32), edd=s["log_edd"].to_numpy(np.float32),
                line=s["line_fid"].to_numpy(), heaviest=cat.loc[ok].nlargest(5, "m_fid")[
                    ["SDSS_NAME", "zbest", "m_fid", "log_lbol", "log_edd", "line_fid"]].to_dict("records"),
                frac_super_edd_all=float(np.mean(s["log_edd"] > 0)),
                median_edd_all=float(s["log_edd"].median()))


# ---------------------------------------------------------------------------
# 7. Growth
# ---------------------------------------------------------------------------

def growth_analysis() -> dict:
    hz = highz.load()
    seeds = highz.seed_table(hz)
    j0313 = highz.j0313_checks()
    ours_j0313 = hz[hz["name"] == "J0313-1806"].iloc[0]
    sens = highz.efficiency_sensitivity(float(ours_j0313["log_mbh"]), float(ours_j0313["z"]))
    tracks = {r["name"]: highz.growth_track(r["log_mbh"], r["z"]) for _, r in hz.iterrows()}
    frac_heavy = {zs: float(np.mean(seeds[f"log_seed_{zs}"] > 2.0)) for zs in (30, 20, 15)}
    frac_heavy_obs = {zs: float(np.mean(seeds[f"log_seed_obs_{zs}"] > 2.0)) for zs in (30, 20, 15)}
    return dict(hz=hz, seeds=seeds, j0313=j0313, j0313_s11=float(ours_j0313["log_mbh"]),
                j0313_vo09=float(ours_j0313["log_mbh_vo09"]), efficiency=sens, tracks=tracks,
                jwst=highz.jwst_seeds(), frac_seed_above_100=frac_heavy,
                frac_seed_above_100_at_observed_rate=frac_heavy_obs,
                efold_myr=physics.efold_time_myr(), t_100_to_1e9_myr=float(physics.time_to_grow_myr(100, 1e9)),
                edd_median=float(np.median(hz["log_edd"])), mass_median=float(np.median(hz["log_mbh"])),
                vo09_minus_s11_median=float(np.median(hz["log_mbh_vo09"] - hz["log_mbh"])))


# ---------------------------------------------------------------------------

def run() -> dict:
    t0 = time.time()
    k06 = kcorr_r06()
    log("loading catalogue (DR16Q + Wu & Shen 2022)")
    cat = catalogue.load(kcorr=k06)
    res = dict(created=time.strftime("%Y-%m-%d %H:%M"))
    log("1. census")
    res["census"] = census(cat)
    log("2. composite spectrum")
    res["composite"] = build_composite()
    log("3. colours, stellar locus, K-corrections")
    res["photometry"] = photometry_analysis(res["composite"], cat, fetch_data.fetch_stars(), k06)
    log("4-5. luminosity function and the quasar era")
    res["qlf"] = qlf_analysis(k06)
    log("6. black hole masses and Eddington ratios")
    res["masses"] = mass_analysis(cat, res["photometry"], k06)
    log("7. the growth problem")
    res["growth"] = growth_analysis()
    # luminosity-redshift data for the figures
    lz = cat[(cat["zbest"] > 0) & np.isfinite(cat["Mi_z2"])]
    res["lum_z"] = dict(z=lz["zbest"].to_numpy(np.float32), mi=lz["Mi_z2"].to_numpy(np.float32),
                        programme=lz["programme"].to_numpy())
    DERIVED.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_PICKLE, "wb") as fh:
        pickle.dump(res, fh)
    RESULTS_JSON.write_text(json.dumps(headline_numbers(res), indent=2, default=float))
    log(f"analysis done in {time.time() - t0:.0f} s")
    return res


def _scalars(d: dict) -> dict:
    """The JSON-friendly part of a results dict (numbers, strings, small nested dicts)."""
    out = {}
    for k, v in d.items():
        if isinstance(v, (int, float, str, np.integer, np.floating)):
            out[k] = v
        elif isinstance(v, dict):
            sub = _scalars(v)
            if sub:
                out[k] = sub
    return out


def headline_numbers(res: dict) -> dict:
    c, q, m, g, p = res["composite"], res["qlf"], res["masses"], res["growth"], res["photometry"]
    return dict(
        n_quasars=res["census"]["n"], z_max=res["census"]["z_max"],
        composite=dict(n=c["n_spectra"], alpha_gm=c["alpha_gm"], alpha_gm_red=c["alpha_gm_red"],
                       alpha_median=c["alpha_median"], alpha_median_red=c["alpha_median_red"],
                       vb=c["vb_comparison"], snr_2800=c["snr_2800"]),
        kcorr=p["k_comparison"], colours=p["agreement"],
        r06_reproduction=q["r06_reproduction"], sk12=q["sk12_comparison"], dr3_vs_dr7=q["dr3_vs_dr7"],
        era_peak=_scalars(q["era_peak"]),
        kollmeier=m["kollmeier"],
        hb_mg={k: v for k, v in m["hb_mg"].items() if k not in ("a", "b")},
        mg_c4={k: v for k, v in m["mg_c4"].items() if k not in ("a", "b")},
        blueshift=_scalars(m["blueshift"]), z2_step=m["z2_step"].round(3).to_dict("records"),
        recipes=dict(mg_hb=m["recipes"]["mg_hb"].round(3).to_dict("records"),
                     c4_mg=m["recipes"]["c4_mg"].round(3).to_dict("records")),
        sk12_without_image_factor=q["sk12_comparison_without_image_factor"],
        composite_luminous=dict(alpha_gm=c["alpha_gm_luminous"], n_gm=c["n_gm_luminous"],
                                alpha_median=c["alpha_median_luminous"], alpha_template=c["alpha_template"],
                                civ_ew=[b["civ_ew"] for b in c["luminosity_bins"]],
                                oiii_ew=[b["oiii_ew"] for b in c["luminosity_bins"]]),
        growth=dict(j0313=g["j0313"], j0313_s11=g["j0313_s11"], j0313_vo09=g["j0313_vo09"],
                    frac_seed_above_100=g["frac_seed_above_100"],
                    efold_myr=g["efold_myr"], t_100_to_1e9_myr=g["t_100_to_1e9_myr"]),
    )


if __name__ == "__main__":
    run()
