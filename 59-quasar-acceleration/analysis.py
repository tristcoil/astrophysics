"""The measurement: fits of the quasar proper-motion field, their error bars,
and the tests that the result is the Galaxy and not Gaia.

Every function returns a JSON-able dict; run.py caches each one in
data/cache/ so figures and README numbers can be regenerated in seconds.
"""

from __future__ import annotations

import math
import multiprocessing as mp
import time

import numpy as np
import pandas as pd

import config as C
import fit
import vsh

A = np.array(C.A_G)
MAS_TO_UAS = 1000.0


# --------------------------------------------------------------------------
# Samples
# --------------------------------------------------------------------------
def make_sample(df: pd.DataFrame, mask=None, error_scale=None) -> fit.Sample:
    d = df if mask is None else df[np.asarray(mask)]
    scale = 1.0 if error_scale is None else np.asarray(error_scale)[np.asarray(mask)] \
        if mask is not None else np.asarray(error_scale)
    return fit.Sample(
        ra=np.radians(d["ra"].to_numpy()), dec=np.radians(d["dec"].to_numpy()),
        pmra=d["pmra"].to_numpy() * MAS_TO_UAS, pmdec=d["pmdec"].to_numpy() * MAS_TO_UAS,
        pmra_error=d["pmra_error"].to_numpy() * MAS_TO_UAS * scale,
        pmdec_error=d["pmdec_error"].to_numpy() * MAS_TO_UAS * scale,
        pmra_pmdec_corr=d["pmra_pmdec_corr"].to_numpy(),
    )


# --------------------------------------------------------------------------
# Vectors: amplitude, direction, uncertainties
# --------------------------------------------------------------------------
def direction(vec: np.ndarray) -> tuple[float, float]:
    """Longitude and latitude (deg) of a 3-vector."""
    lon = math.degrees(math.atan2(vec[1], vec[0])) % 360.0
    lat = math.degrees(math.asin(vec[2] / np.linalg.norm(vec)))
    return lon, lat


def direction_errors(vec: np.ndarray, cov: np.ndarray) -> tuple[float, float, float]:
    """Linear propagation: sigma(lon) as a coordinate, sigma(lat), correlation."""
    x, y, z = vec
    r2 = x * x + y * y
    r = math.sqrt(r2)
    n2 = r2 + z * z
    j = np.array([[-y / r2, x / r2, 0.0],
                  [-x * z / (n2 * r), -y * z / (n2 * r), r / n2]])
    c = j @ cov @ j.T
    s = np.sqrt(np.diag(c))
    return math.degrees(s[0]), math.degrees(s[1]), float(c[0, 1] / (s[0] * s[1]))


def amplitude(vec: np.ndarray, cov: np.ndarray) -> dict:
    """|g|, its error, and Klioner et al.'s Appendix A bias correction (A.4)."""
    r = float(np.linalg.norm(vec))
    err = float(math.sqrt(vec @ cov @ vec) / r)
    bias = 0.5 * np.trace(cov) / r - 0.5 * (vec @ cov @ vec) / r**3
    return {"value": r, "err": err, "bias": float(bias), "corrected": r - float(bias)}


def angle_between(v1: np.ndarray, v2: np.ndarray) -> float:
    c = float(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


def gc_unit_icrs() -> np.ndarray:
    """Unit vector to Sgr A*, ICRS."""
    ra, dec = math.radians(C.SGRA_RA_DEG), math.radians(C.SGRA_DEC_DEG)
    return np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])


def describe(res: fit.FitResult, cov: str = "scaled", boot: np.ndarray | None = None) -> dict:
    """Glide and rotation of a fit, in both frames, with the chosen errors."""
    out = {"lmax": res.lmax, "kappa": res.kappa, "n_used": res.n_used,
           "n_outliers": res.n_outliers, "uwe": res.uwe, "iterations": res.iterations,
           "errors": cov}
    for kind, name in (("S", "glide"), ("T", "rotation")):
        for frame in ("icrs", "galactic"):
            vec, c = res.vector(kind, frame, cov)
            if c is None:
                continue
            lon, lat = direction(vec)
            slon, slat, rho = direction_errors(vec, c)
            sig = np.sqrt(np.diag(c))
            corr = c / np.outer(sig, sig)
            out[f"{name}_{frame}"] = {
                "vec": vec.tolist(), "err": sig.tolist(),
                "corr": {"xy": corr[0, 1], "xz": corr[0, 2], "yz": corr[1, 2]},
                "lon": lon, "lat": lat, "lon_err": slon, "lat_err": slat, "lonlat_corr": rho,
                "amplitude": amplitude(vec, c),
            }
    g_icrs, _ = res.vector("S", "icrs", cov)
    out["glide_angle_from_sgra_deg"] = angle_between(g_icrs, gc_unit_icrs())
    if boot is not None:
        m = vsh.l1_vector_matrix("S", res.term_list)
        gb = boot @ m.T
        gal = gb @ A.T
        lb = np.array([direction(v) for v in gal])
        lon = (lb[:, 0] + 180.0) % 360.0 - 180.0
        out["bootstrap_direction"] = {
            "l_sd": float(np.std(lon)), "b_sd": float(np.std(lb[:, 1])),
            "angle_from_sgra_median": float(np.median([angle_between(v, gc_unit_icrs()) for v in gb])),
        }
    return out


# --------------------------------------------------------------------------
# 1. The baseline: Klioner et al. (2021)'s analysis, redone
# --------------------------------------------------------------------------
def _boot_worker(args):
    seed, n_rep, kappa = args
    from threadpoolctl import threadpool_limits
    with threadpool_limits(1):
        return fit.bootstrap(_SHARED["design"], n_rep, seed, kappa=kappa)


_SHARED: dict = {}


def parallel_bootstrap(design: fit.Design, n_rep: int, seed: int, kappa: float = 3.0,
                       workers: int = 4) -> np.ndarray:
    """Bootstrap replicates in forked worker processes sharing the cached basis."""
    _SHARED["design"] = design
    per = [n_rep // workers + (1 if i < n_rep % workers else 0) for i in range(workers)]
    ctx = mp.get_context("fork")
    with ctx.Pool(workers) as pool:
        parts = pool.map(_boot_worker, [(seed + 1000 * i, k, kappa) for i, k in enumerate(per) if k])
    return np.vstack(parts)


def baseline(df: pd.DataFrame, n_boot: int = C.N_BOOTSTRAP, workers: int = 4) -> dict:
    """Five-parameter Gaia-CRF3 sources, l_max = 3, kappa = 3, full weights,
    bootstrap errors: Klioner et al. (2021)'s choices exactly."""
    sel = df["five_param"].to_numpy()
    sample = make_sample(df, sel)
    tl = vsh.terms(C.LMAX_BASELINE)
    design = fit.Design(sample, tl)
    t0 = time.time()
    res = fit.fit_design(design, kappa=C.KAPPA_BASELINE)
    t_fit = time.time() - t0
    t0 = time.time()
    boot = parallel_bootstrap(design, n_boot, C.SEED, C.KAPPA_BASELINE, workers)
    res.boot_coef = boot
    res.cov_bootstrap = np.cov(boot.T)
    out = {"n_sample": int(sel.sum()), "fit_seconds": t_fit,
           "bootstrap_seconds": time.time() - t0, "n_bootstrap": len(boot)}
    for kind in ("formal", "scaled", "sandwich", "bootstrap"):
        out[kind] = describe(res, kind, boot if kind == "bootstrap" else None)
    m = vsh.l1_vector_matrix("S", tl)
    out["bootstrap_over_scaled_err"] = (
        np.sqrt(np.diag(m @ res.cov_bootstrap @ m.T)) / np.sqrt(np.diag(m @ res.cov_scaled @ m.T))
    ).tolist()
    out["vsh_coefficients"] = {t.label: float(c) for t, c in zip(tl, res.coef)}
    out["outlier_source_ids"] = df.loc[sel, "source_id"].to_numpy()[res.outlier].tolist()
    # bootstrap glide vectors (galactic) for the direction cloud figure
    np.save(C.CACHE_DIR / "baseline_bootstrap_glide_gal.npy", (boot @ m.T) @ A.T)
    np.save(C.CACHE_DIR / "baseline_x.npy", res.x.astype(np.float32))
    out["grade"] = grade(out["bootstrap"])
    return out


def grade(desc: dict) -> dict:
    """Our numbers against Klioner et al. (2021) Table 2, in their sigmas."""
    k = C.REFERENCE["klioner2021"]
    g = {}
    for frame, key, keyerr in (("icrs", "g_icrs", "g_icrs_err"), ("galactic", "g_gal", "g_gal_err")):
        ours = desc[f"glide_{frame}"]
        for i, comp in enumerate("xyz" if frame == "icrs" else "XYZ"):
            g[f"g{comp}"] = {"ours": ours["vec"][i], "ours_err": ours["err"][i],
                             "published": k[key][i], "published_err": k[keyerr][i],
                             "diff_sigma": (ours["vec"][i] - k[key][i]) / k[keyerr][i]}
    amp = desc["glide_icrs"]["amplitude"]
    g["amplitude"] = {"ours": amp["value"], "ours_err": amp["err"],
                      "published": k["amplitude_uas_yr"][0], "published_err": k["amplitude_uas_yr"][1],
                      "diff_sigma": (amp["value"] - k["amplitude_uas_yr"][0]) / k["amplitude_uas_yr"][1]}
    for name, (lon, lat) in (("ra_dec", ("lon", "lat")),):
        pass
    gi, gg = desc["glide_icrs"], desc["glide_galactic"]
    g["ra"] = {"ours": gi["lon"], "ours_err": gi["lon_err"], "published": k["ra_dec"][0],
               "published_err": k["ra_dec_err"][0]}
    g["dec"] = {"ours": gi["lat"], "ours_err": gi["lat_err"], "published": k["ra_dec"][1],
                "published_err": k["ra_dec_err"][1]}
    g["l"] = {"ours": gg["lon"], "ours_err": gg["lon_err"], "published": k["l_b"][0],
              "published_err": k["l_b_err"][0]}
    g["b"] = {"ours": gg["lat"], "ours_err": gg["lat_err"], "published": k["l_b"][1],
              "published_err": k["l_b_err"][1]}
    return g


# --------------------------------------------------------------------------
# 2. How the glide depends on the analysis choices
# --------------------------------------------------------------------------
def _glide_row(res: fit.FitResult, label, **extra) -> dict:
    row = {"label": label, "n_used": res.n_used, "n_outliers": res.n_outliers, "uwe": res.uwe,
           "lmax": res.lmax, "kappa": res.kappa, **extra}
    for frame in ("icrs", "galactic"):
        for kind, name in (("S", "glide"), ("T", "rotation")):
            vec, cov = res.vector(kind, frame, "scaled")
            if not np.any(vec) and kind == "T":
                continue
            row[f"{name}_{frame}"] = vec.tolist()
            row[f"{name}_{frame}_err"] = np.sqrt(np.diag(cov)).tolist()
            if res.cov_sandwich is not None:
                _, cs = res.vector(kind, frame, "sandwich")
                row[f"{name}_{frame}_err_sandwich"] = np.sqrt(np.diag(cs)).tolist()
    g, cov = res.vector("S", "icrs", "scaled")
    row["amplitude"] = amplitude(g, cov)
    row["l"], row["b"] = direction(A @ g)
    row["angle_from_sgra_deg"] = angle_between(g, gc_unit_icrs())
    return row


def lmax_scan(df: pd.DataFrame, lmaxes=C.LMAX_SCAN) -> list[dict]:
    sample = make_sample(df, df["five_param"].to_numpy())
    rows = []
    alone = [t for t in vsh.terms(1) if t.kind == "S"]
    res = fit.fit_vsh(sample, term_list=alone, kappa=C.KAPPA_BASELINE)
    rows.append(_glide_row(res, "alone"))
    for lmax in lmaxes:
        t0 = time.time()
        res = fit.fit_vsh(sample, lmax=lmax, kappa=C.KAPPA_BASELINE)
        row = _glide_row(res, str(lmax), seconds=time.time() - t0)
        row["power"] = {k: v.tolist() for k, v in vsh.power_per_degree(res.coef, res.term_list).items()}
        noise = {"T": np.zeros(lmax), "S": np.zeros(lmax)}
        for i, t in enumerate(res.term_list):
            noise[t.kind][t.l - 1] += res.cov_scaled[i, i]
        row["noise_power"] = {k: np.sqrt(v / (4 * math.pi)).tolist() for k, v in noise.items()}
        rows.append(row)
        print(f"    l_max {lmax:>2}: |g| = {row['amplitude']['value']:.3f} +- "
              f"{row['amplitude']['err']:.3f}, {row['seconds']:.0f} s", flush=True)
    return rows


def choice_tests(df: pd.DataFrame) -> list[dict]:
    """Clip limit, correlations on/off, and the six-parameter sources."""
    rows = []
    s5 = make_sample(df, df["five_param"].to_numpy())
    d5 = fit.Design(s5, vsh.terms(C.LMAX_BASELINE))
    for kappa in C.KAPPA_SCAN:
        res = fit.fit_design(d5, kappa=kappa)
        rows.append(_glide_row(res, f"kappa={kappa}", group="kappa", kappa_value=kappa))
    res = fit.fit_design(d5, kappa=C.KAPPA_BASELINE, use_corr=False)
    rows.append(_glide_row(res, "no correlations", group="weights"))
    del d5
    for label, mask in (("six-parameter only", ~df["five_param"].to_numpy()),
                        ("five + six parameter", np.ones(len(df), bool))):
        res = fit.fit_vsh(make_sample(df, mask), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="solution type", n_sample=int(mask.sum())))
    return rows


def subset_fits(df: pd.DataFrame) -> list[dict]:
    """The glide must not depend on what the quasars are, only on where we are."""
    rows = []
    five = df["five_param"].to_numpy()
    g = df["phot_g_mean_mag"].to_numpy()
    for (lo, hi), n_pub in zip(C.G_BINS, C.G_BIN_COUNTS_K21):
        m = five & (g > lo) & (g <= hi)
        res = fit.fit_vsh(make_sample(df, m), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, f"G {lo}-{hi}", group="magnitude", n_sample=int(m.sum()),
                               n_published=n_pub, g_median=float(np.median(g[m]))))
    nu = df["nu_eff_used_in_astrometry"].to_numpy()
    q = np.nanquantile(nu[five], [0, 0.25, 0.5, 0.75, 1.0])
    for i in range(4):
        m = five & (nu >= q[i]) & ((nu < q[i + 1]) if i < 3 else (nu <= q[i + 1]))
        res = fit.fit_vsh(make_sample(df, m), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, f"nu_eff {q[i]:.3f}-{q[i + 1]:.3f}", group="colour",
                               n_sample=int(m.sum()), nu_median=float(np.median(nu[m]))))
    # Redshift: Gaia's own all-sky redshifts (QSOC), cleanest flag only.
    zq = df["redshift_qsoc"].to_numpy()
    flag = df["flags_qsoc"].to_numpy()
    ok = five & np.isfinite(zq) & (flag == 0)
    for lo, hi in ((0.0, 0.8), (0.8, 1.3), (1.3, 1.8), (1.8, 2.5), (2.5, 6.0)):
        m = ok & (zq >= lo) & (zq < hi)
        res = fit.fit_vsh(make_sample(df, m), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, f"z_QSOC {lo}-{hi}", group="redshift (Gaia)",
                               n_sample=int(m.sum()), z_median=float(np.median(zq[m]))))
    # Which external catalogue put the source in the Gaia-CRF3.
    cat = df["catalogue_name"].astype(str).to_numpy()
    # SDSS covers a third of the sky: fit glide + rotation only there.
    groups = {"AllWISE": ("AllWISE", 3), "R90": ("R90", 3), "Gaia-unWISE": ("Gaia-unWISE", 3),
              "Milliquas": ("Milliquas", 3), "SDSS DR14Q (l_max=1)": ("SDSS DR14Q", 1)}
    assigned = np.zeros(len(df), bool)
    for label, (key, lmax) in groups.items():
        m = five & np.char.startswith(cat.astype(str), key)
        assigned |= m
        res = fit.fit_vsh(make_sample(df, m), lmax=lmax, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="catalogue", n_sample=int(m.sum())))
    m = five & ~assigned
    res = fit.fit_vsh(make_sample(df, m), lmax=1, kappa=C.KAPPA_BASELINE)
    rows.append(_glide_row(res, "other catalogues (l_max=1)", group="catalogue", n_sample=int(m.sum())))
    return rows


def dr16q_redshift_fits(df: pd.DataFrame, coefficients: dict) -> list[dict]:
    """Spectroscopic redshifts from SDSS cover a third of the sky. There a fit
    cannot tell the glide from Gaia's large-scale patterns, so two versions:
    glide + rotation only (l_max = 1), and the glide alone after subtracting
    the all-sky fit's non-glide terms (its rotation and l = 2-3 patterns).
    A control sample -- every Gaia-CRF3 quasar in the same sky pixels --
    shows what the footprint alone does."""
    import healpy as hp
    rows = []
    five = df["five_param"].to_numpy()
    z = df["z_dr16q"].to_numpy()
    has = five & np.isfinite(z)
    pix = hp.ang2pix(16, df["ra"].to_numpy(), df["dec"].to_numpy(), lonlat=True)
    footprint = np.isin(pix, np.unique(pix[has])) & five
    clean = remove_nonglide(df, five, coefficients)
    alone = [t for t in vsh.terms(1) if t.kind == "S"]
    groups = [("CRF3 quasars in the SDSS footprint", footprint, np.nan),
              ("all DR16Q matches", has, float(np.median(z[has])))]
    for lo, hi in ((0.0, 1.0), (1.0, 1.6), (1.6, 2.2), (2.2, 7.0)):
        m = has & (z >= lo) & (z < hi)
        groups.append((f"z_SDSS {lo}-{hi}", m, float(np.median(z[m]))))
    for label, m, zmed in groups:
        res = fit.fit_vsh(make_sample(df, m), lmax=1, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="redshift (SDSS, l_max=1)", n_sample=int(m.sum()),
                               z_median=zmed))
        res = fit.fit_vsh(make_sample(clean, m), term_list=alone, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="redshift (SDSS, all-sky systematics removed)",
                               n_sample=int(m.sum()), z_median=zmed))
    return rows


# --------------------------------------------------------------------------
# 3. Gaia's error bars on quasars, which should not move at all
# --------------------------------------------------------------------------
def gaussian_width(u: np.ndarray, lim: float = 5.0, bins: int = 200) -> dict:
    """Width of a normalised quantity: robust (MAD) and a Gaussian fitted to
    its histogram within |u| < lim (Poisson likelihood), as in Klioner et al.'s
    Fig. 6."""
    from scipy.optimize import minimize
    u = u[np.isfinite(u)]
    mad = 1.4826 * np.median(np.abs(u - np.median(u)))
    edges = np.linspace(-lim, lim, bins + 1)
    counts, _ = np.histogram(u, edges)
    centres = 0.5 * (edges[1:] + edges[:-1])
    width = edges[1] - edges[0]

    def nll(p):
        amp, mu, sig = np.exp(p[0]), p[1], np.exp(p[2])
        lam = amp * width * np.exp(-0.5 * ((centres - mu) / sig) ** 2) / (sig * np.sqrt(2 * np.pi))
        return np.sum(lam - counts * np.log(lam + 1e-300))

    r = minimize(nll, [np.log(u.size), 0.0, np.log(mad)], method="Nelder-Mead",
                 options={"xatol": 1e-8, "fatol": 1e-8, "maxiter": 4000})
    amp, mu, sig = np.exp(r.x[0]), r.x[1], np.exp(r.x[2])
    expected_tail = 2 * 0.5 * math.erfc(4.0 / (sig * math.sqrt(2))) * u.size
    return {"n": int(u.size), "mad_sigma": float(mad), "fit_sigma": float(sig), "fit_mean": float(mu),
            "std": float(np.std(u)), "frac_beyond_4": float(np.mean(np.abs(u) > 4)),
            "expected_frac_beyond_4": float(expected_tail / u.size),
            "hist_edges": edges.tolist(), "hist_counts": counts.tolist()}


def error_model(df: pd.DataFrame) -> dict:
    five = df["five_param"].to_numpy()
    d = df[five]
    zp = -0.017  # the median parallax of EDR3 quasars (Lindegren et al. 2021)
    out = {"parallax": gaussian_width((d["parallax"] / d["parallax_error"]).to_numpy()),
           "parallax_minus_zp": gaussian_width(((d["parallax"] - zp) / d["parallax_error"]).to_numpy()),
           "pmra": gaussian_width((d["pmra"] / d["pmra_error"]).to_numpy()),
           "pmdec": gaussian_width((d["pmdec"] / d["pmdec_error"]).to_numpy())}
    # The weighted residual X^2 of each quasar against zero proper motion:
    # for honest Gaussian errors it follows chi-square with 2 dof, median 2 ln 2.
    s = make_sample(df, five)
    w_aa, w_ad, w_dd = s.weights()
    x2 = w_aa * s.pmra**2 + 2 * w_ad * s.pmra * s.pmdec + w_dd * s.pmdec**2
    out["k_overall"] = float(np.sqrt(np.median(x2) / (2 * math.log(2))))

    def binned(col: np.ndarray, edges: np.ndarray) -> list[dict]:
        rows = []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (col >= lo) & (col < hi)
            if m.sum() < 300:
                continue
            k = math.sqrt(np.median(x2[m]) / (2 * math.log(2)))
            rows.append({"lo": float(lo), "hi": float(hi), "n": int(m.sum()), "k": k,
                         "k_err": k / math.sqrt(2 * m.sum()) * 1.25,  # median of chi2_2: ~1.25x mean's error
                         "median_pm_error_uas": float(np.median(np.hypot(s.pmra_error[m], s.pmdec_error[m])
                                                                / math.sqrt(2)))})
        return rows

    gmag = d["phot_g_mean_mag"].to_numpy()
    out["k_vs_g"] = binned(gmag, np.arange(14.0, 21.51, 0.25))
    out["k_vs_nu"] = binned(d["nu_eff_used_in_astrometry"].to_numpy(), np.arange(1.30, 2.01, 0.05))
    out["k_vs_ruwe"] = binned(d["ruwe"].to_numpy(), np.array([0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.4, 2, 5]))
    out["k_vs_ecl_lat"] = binned(np.abs(d["ecl_lat"].to_numpy()), np.arange(0, 91, 10.0))
    out["k_vs_aen"] = binned(d["astrometric_excess_noise"].to_numpy(), np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 3]))
    out["k_vs_bprp_excess"] = binned(d["phot_bp_rp_excess_factor"].to_numpy(),
                                     np.array([1.0, 1.15, 1.2, 1.25, 1.3, 1.4, 1.6, 2.0, 3.0]))
    # The magnitude-and-colour error model: k on a G x nu_eff grid, smoothed
    # by fitting k = k0 + k1 (G-20) + k2 (G-20)^2 + k3 (nu-1.59) to the cells.
    nu = d["nu_eff_used_in_astrometry"].to_numpy()
    ge, ne = np.arange(15.0, 21.6, 0.5), np.arange(1.35, 1.86, 0.05)
    cells = []
    for i in range(len(ge) - 1):
        for j in range(len(ne) - 1):
            m = (gmag >= ge[i]) & (gmag < ge[i + 1]) & (nu >= ne[j]) & (nu < ne[j + 1])
            if m.sum() >= 500:
                cells.append((0.5 * (ge[i] + ge[i + 1]), 0.5 * (ne[j] + ne[j + 1]), int(m.sum()),
                              math.sqrt(np.median(x2[m]) / (2 * math.log(2)))))
    cells = np.array(cells)
    dm, dn = cells[:, 0] - 20.0, cells[:, 1] - 1.59
    design = np.column_stack([np.ones_like(dm), dm, dm**2, dn])
    coef, *_ = np.linalg.lstsq(design * np.sqrt(cells[:, 2])[:, None], cells[:, 3] * np.sqrt(cells[:, 2]),
                               rcond=None)
    out["k_model"] = {"form": "k = c0 + c1 (G-20) + c2 (G-20)^2 + c3 (nu_eff-1.59)",
                      "coef": coef.tolist(), "cells": cells.tolist(),
                      "rms_cell_residual": float(np.sqrt(np.average((design @ coef - cells[:, 3]) ** 2,
                                                                     weights=cells[:, 2])))}
    return out


def k_model_factor(df: pd.DataFrame, coef) -> np.ndarray:
    # one five-parameter source has no G magnitude: give it the median's factor
    g = np.clip(np.nan_to_num(df["phot_g_mean_mag"].to_numpy(), nan=20.0), 15.0, 21.3) - 20.0
    nu = np.clip(np.nan_to_num(df["colour"].to_numpy(), nan=1.59), 1.35, 1.85) - 1.59
    return coef[0] + coef[1] * g + coef[2] * g**2 + coef[3] * nu


def error_model_refit(df: pd.DataFrame, k_coef) -> dict:
    """Refit with Gaia's errors scaled by the magnitude-and-colour model."""
    five = df["five_param"].to_numpy()
    k = k_model_factor(df, k_coef)
    res = fit.fit_vsh(make_sample(df, five, error_scale=k), lmax=C.LMAX_BASELINE,
                      kappa=C.KAPPA_BASELINE)
    return _glide_row(res, "errors x k(G, nu_eff)", group="error model",
                      k_range=[float(k[five].min()), float(k[five].max())])


# --------------------------------------------------------------------------
# 4. Stars pretending to be quasars
# --------------------------------------------------------------------------
def contamination(df: pd.DataFrame, baseline_outliers: list[int]) -> dict:
    five = df["five_param"].to_numpy()
    d = df[five]
    out = {}
    plx = ((d["parallax"] + 0.017) / d["parallax_error"]).to_numpy()
    n = plx.size
    # Excess of significantly positive over negative parallaxes: stars have
    # positive parallaxes, quasars scatter symmetrically about the zero point.
    for t in (2.0, 3.0, 4.0):
        pos, neg = int(np.sum(plx > t)), int(np.sum(plx < -t))
        out[f"parallax_excess_{t:.0f}sigma"] = {"positive": pos, "negative": neg,
                                                "excess": pos - neg, "excess_frac": (pos - neg) / n}
    pstar = d["classprob_dsc_combmod_star"].to_numpy()
    out["dsc_star_prob_gt_0.5"] = int(np.sum(pstar > 0.5))
    out["dsc_star_prob_finite"] = int(np.isfinite(pstar).sum())
    out["n_five_param"] = int(n)
    rows = []
    for label, m in (
        ("DSC P(star) > 0.5", five & (df["classprob_dsc_combmod_star"].to_numpy() > 0.5)),
        ("DSC P(quasar) > 0.5", five & (df["classprob_dsc_combmod_quasar"].to_numpy() > 0.5)),
        ("parallax > 2 sigma", five & (((df["parallax"] + 0.017) / df["parallax_error"]).to_numpy() > 2)),
        ("parallax < -2 sigma", five & (((df["parallax"] + 0.017) / df["parallax_error"]).to_numpy() < -2)),
    ):
        if m.sum() < 2000:
            continue
        res = fit.fit_vsh(make_sample(df, m), lmax=1, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="contamination", n_sample=int(m.sum())))
    # The headline fit without the sources Gaia's own classifier calls stars.
    pstar = df["classprob_dsc_combmod_star"].to_numpy()
    for label, m in (("baseline without DSC P(star) > 0.5", five & ~(pstar > 0.5)),
                     ("DSC P(star) > 0.5 only (l_max=3)", five & (pstar > 0.5))):
        res = fit.fit_vsh(make_sample(df, m), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="contamination", n_sample=int(m.sum())))
    # The quasars the clipping threw out: what do they do on their own?
    out_ids = np.isin(df["source_id"].to_numpy(), np.asarray(baseline_outliers))
    res = fit.fit_vsh(make_sample(df, out_ids), lmax=1, kappa=math.inf)
    rows.append(_glide_row(res, "baseline outliers (no clipping)", group="contamination",
                           n_sample=int(out_ids.sum())))
    out["fits"] = rows
    return out


# --------------------------------------------------------------------------
# 5. How the signal emerges from the noise
# --------------------------------------------------------------------------
def emergence(df: pd.DataFrame, sizes=(1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000),
              draws: int = 5, lmax: int = C.LMAX_BASELINE) -> list[dict]:
    rng = np.random.default_rng(C.SEED + 7)
    idx5 = np.flatnonzero(df["five_param"].to_numpy())
    rows = []
    for n in sizes:
        for k in range(draws):
            pick = np.zeros(len(df), bool)
            pick[rng.choice(idx5, size=n, replace=False)] = True
            res = fit.fit_vsh(make_sample(df, pick), lmax=lmax, kappa=C.KAPPA_BASELINE, sandwich=False)
            row = _glide_row(res, f"N={n}", n_sample=n, draw=k)
            rows.append(row)
    return rows


# --------------------------------------------------------------------------
# 6. The picture: proper motions averaged over patches of sky
# --------------------------------------------------------------------------
def unit_vectors(ra_deg, dec_deg) -> np.ndarray:
    ra, dec = np.radians(ra_deg), np.radians(dec_deg)
    return np.column_stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def tangent_basis(u: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Local east (increasing longitude) and north unit vectors for unit vectors u
    in whatever frame u is expressed in."""
    lon = np.arctan2(u[:, 1], u[:, 0])
    lat = np.arcsin(np.clip(u[:, 2], -1, 1))
    e_lon = np.column_stack([-np.sin(lon), np.cos(lon), np.zeros_like(lon)])
    e_lat = np.column_stack([-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat)])
    return e_lon, e_lat


def cap_averages(df: pd.DataFrame, mask, centres_gal: np.ndarray, radius_deg: float,
                 model_glide_icrs: np.ndarray | None = None, exclude=None) -> dict:
    """Weighted mean proper-motion vector of the quasars within ``radius_deg``
    of each centre (galactic unit vectors), as (east, north) components in
    galactic coordinates at the centre. Each quasar enters with its full 2 x 2
    weight; the vectors are moved to the centre by projecting their 3D form
    onto the centre's tangent plane. The same average of the fitted glide
    field is returned alongside, so data and model are compared like for like."""
    from scipy.spatial import cKDTree
    m = np.asarray(mask) if exclude is None else np.asarray(mask) & ~np.asarray(exclude)
    d = df[m]
    s = make_sample(df, m)
    u_icrs = unit_vectors(d["ra"].to_numpy(), d["dec"].to_numpy())
    e_ra, e_dec = tangent_basis(u_icrs)
    w_aa, w_ad, w_dd = s.weights()
    # Per-quasar 3D weight tensor W3 = E W E^T and weighted vector E W mu, in galactic axes.
    eg_ra, eg_dec = e_ra @ A.T, e_dec @ A.T
    wv_ra = w_aa * s.pmra + w_ad * s.pmdec
    wv_dec = w_ad * s.pmra + w_dd * s.pmdec
    vec = eg_ra * wv_ra[:, None] + eg_dec * wv_dec[:, None]
    if model_glide_icrs is not None:
        gm_ra, gm_dec = vsh.glide_field(s.ra, s.dec, model_glide_icrs)
        mv_ra = w_aa * gm_ra + w_ad * gm_dec
        mv_dec = w_ad * gm_ra + w_dd * gm_dec
        mvec = eg_ra * mv_ra[:, None] + eg_dec * mv_dec[:, None]
    u_gal = u_icrs @ A.T
    tree = cKDTree(u_gal)
    chord = 2 * math.sin(math.radians(radius_deg) / 2)
    c_e, c_n = tangent_basis(centres_gal)
    res = {"east": [], "north": [], "cov": [], "n": [], "model_east": [], "model_north": []}
    for i, c in enumerate(centres_gal):
        idx = tree.query_ball_point(c, chord)
        if len(idx) < 50:
            for k in res:
                res[k].append(np.nan if k != "cov" else [np.nan] * 3)
            continue
        idx = np.asarray(idx)
        E = np.column_stack([c_e[i], c_n[i]])  # 3 x 2
        # W at the centre = sum_k E^T (E_k W_k E_k^T) E
        Pa = eg_ra[idx] @ E  # n x 2: projections of the quasar's e_ra on the centre basis
        Pd = eg_dec[idx] @ E
        Wc = (Pa.T * w_aa[idx]) @ Pa + (Pa.T * w_ad[idx]) @ Pd + (Pd.T * w_ad[idx]) @ Pa \
            + (Pd.T * w_dd[idx]) @ Pd
        bc = E.T @ vec[idx].sum(axis=0)
        cov = np.linalg.inv(Wc)
        v = cov @ bc
        res["east"].append(float(v[0]))
        res["north"].append(float(v[1]))
        res["cov"].append([float(cov[0, 0]), float(cov[0, 1]), float(cov[1, 1])])
        res["n"].append(int(idx.size))
        if model_glide_icrs is not None:
            mvc = cov @ (E.T @ mvec[idx].sum(axis=0))
            res["model_east"].append(float(mvc[0]))
            res["model_north"].append(float(mvc[1]))
    return res


def drift_profile(df: pd.DataFrame, mask, glide_icrs: np.ndarray, rotation_icrs=None,
                  bin_deg: float = 10.0, exclude=None) -> dict:
    """Mean proper motion along the great circle towards the Galactic Centre
    (and across it) in rings of angular distance from the centre. A glide
    towards the centre predicts |g| sin(theta) along and zero across.
    The fitted frame rotation, if given, is removed first."""
    m = np.asarray(mask) if exclude is None else np.asarray(mask) & ~np.asarray(exclude)
    d = df[m]
    s = make_sample(df, m)
    pmra, pmdec = s.pmra.copy(), s.pmdec.copy()
    if rotation_icrs is not None:
        r_ra, r_dec = vsh.rotation_field(s.ra, s.dec, rotation_icrs)
        pmra -= r_ra
        pmdec -= r_dec
    u = unit_vectors(d["ra"].to_numpy(), d["dec"].to_numpy())
    gc = gc_unit_icrs()
    cos_t = np.clip(u @ gc, -1, 1)
    theta = np.degrees(np.arccos(cos_t))
    # unit tangent towards the centre at each quasar
    t = gc[None, :] - cos_t[:, None] * u
    t /= np.linalg.norm(t, axis=1)[:, None]
    e_ra, e_dec = tangent_basis(u)
    ta, td = np.einsum("ij,ij->i", t, e_ra), np.einsum("ij,ij->i", t, e_dec)
    # across = u x t (perpendicular, in the tangent plane)
    xa, xd = -td, ta
    w_aa, w_ad, w_dd = s.weights()
    rows = []
    edges = np.arange(0.0, 180.0 + bin_deg / 2, bin_deg)
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (theta >= lo) & (theta < hi)
        if k.sum() < 30:  # the Galactic Centre itself has no quasars
            continue
        # weighted LSQ for a constant (along, across) pair within the ring
        Pa = np.column_stack([ta[k], xa[k]])
        Pd = np.column_stack([td[k], xd[k]])
        W = (Pa.T * w_aa[k]) @ Pa + (Pa.T * w_ad[k]) @ Pd + (Pd.T * w_ad[k]) @ Pa + (Pd.T * w_dd[k]) @ Pd
        b = Pa.T @ (w_aa[k] * pmra[k] + w_ad[k] * pmdec[k]) + Pd.T @ (w_ad[k] * pmra[k] + w_dd[k] * pmdec[k])
        cov = np.linalg.inv(W)
        v = cov @ b
        rows.append({"lo": float(lo), "hi": float(hi), "n": int(k.sum()),
                     "theta_mean": float(np.average(theta[k])),
                     "along": float(v[0]), "along_err": float(math.sqrt(cov[0, 0])),
                     "across": float(v[1]), "across_err": float(math.sqrt(cov[1, 1]))})
    g_par = float(glide_icrs @ gc)
    return {"rings": rows, "glide_towards_gc": g_par,
            "glide_perp_to_gc": float(np.linalg.norm(glide_icrs - g_par * gc))}


# --------------------------------------------------------------------------
# 7. Using what we know: a quasar's true parallax is zero
# --------------------------------------------------------------------------
def lindegren_zero_point(df: pd.DataFrame) -> np.ndarray:
    """Gaia EDR3 parallax zero point of every source (Lindegren et al. 2021,
    the gaiadr3_zeropoint package), clipped to its calibrated ranges
    (G 6-21, nu_eff 1.1-1.9, pseudocolour 1.24-1.72). Quasars fainter than
    G = 21 get the G = 21 value."""
    from zero_point import zpt
    zpt.load_tables()
    # The package returns NaN if ANY input is NaN (even the colour column a
    # solution type does not use) or sits on a range limit: clip just inside,
    # and give each solution type a dummy for the colour it ignores.
    g = np.clip(np.nan_to_num(df["phot_g_mean_mag"].to_numpy(float), nan=20.0), 6.001, 20.999)
    nu = np.clip(np.nan_to_num(df["nu_eff_used_in_astrometry"].to_numpy(float), nan=1.5), 1.1001, 1.8999)
    pc = np.clip(np.nan_to_num(df["pseudocolour"].to_numpy(float), nan=1.5), 1.2401, 1.7199)
    beta = df["ecl_lat"].to_numpy(float)
    sol = df["astrometric_params_solved"].to_numpy()
    out = np.full(len(df), np.nan)
    for code in (31, 95):
        m = sol == code
        out[m] = zpt.get_zpt(g[m], nu[m], pc[m], beta[m], sol[m], _warnings=False)
    return out


def conditioned_sample(df: pd.DataFrame, mask, zp_mas) -> fit.Sample:
    """Proper motions conditioned on the parallax being its zero point: the
    part of each quasar's proper-motion error that is correlated with its
    parallax error is predicted from (parallax - zero point) and removed.
    For Gaussian errors: mu|plx = mu - C_mu,plx / s_plx^2 (plx - zp), with
    the conditional covariance C_mu - C_mu,plx C_plx,mu / s_plx^2."""
    m = np.asarray(mask)
    d = df[m]
    zp = np.broadcast_to(np.asarray(zp_mas, float), (len(df),))[m]
    resid = d["parallax"].to_numpy() - zp
    sp = d["parallax_error"].to_numpy(float)
    ra = d["parallax_pmra_corr"].to_numpy(float)
    rd = d["parallax_pmdec_corr"].to_numpy(float)
    rad = d["pmra_pmdec_corr"].to_numpy(float)
    sa = d["pmra_error"].to_numpy(float)
    sd = d["pmdec_error"].to_numpy(float)
    return fit.Sample(
        ra=np.radians(d["ra"].to_numpy()), dec=np.radians(d["dec"].to_numpy()),
        pmra=(d["pmra"].to_numpy() - ra * sa / sp * resid) * MAS_TO_UAS,
        pmdec=(d["pmdec"].to_numpy() - rd * sd / sp * resid) * MAS_TO_UAS,
        pmra_error=sa * np.sqrt(1 - ra**2) * MAS_TO_UAS,
        pmdec_error=sd * np.sqrt(1 - rd**2) * MAS_TO_UAS,
        pmra_pmdec_corr=(rad - ra * rd) / np.sqrt((1 - ra**2) * (1 - rd**2)),
    )


def conditioning_tests(df: pd.DataFrame, zp_l21: np.ndarray) -> list[dict]:
    five = df["five_param"].to_numpy()
    rows = []
    for label, zp in (("conditioned, zero point 0 (naive)", 0.0),
                      ("conditioned, zero point -17 uas", -0.017),
                      ("conditioned, Lindegren+21 zero point", zp_l21)):
        res = fit.fit_vsh(conditioned_sample(df, five, zp), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
        rows.append(_glide_row(res, label, group="parallax conditioning"))
    s = conditioned_sample(df, five, zp_l21)
    gain = np.median(np.hypot(s.pmra_error, s.pmdec_error) /
                     (np.hypot(df.loc[five, "pmra_error"], df.loc[five, "pmdec_error"]).to_numpy() * MAS_TO_UAS))
    rows.append({"label": "median error ratio, conditioned / raw", "value": float(gain)})
    return rows


# --------------------------------------------------------------------------
# 8. Does the drift depend on distance? (Quaia redshifts)
# --------------------------------------------------------------------------
QUAIA_Z_EDGES = (0.0, 0.8, 1.25, 1.6, 2.0, 2.5, 5.0)


def attach_quaia(df: pd.DataFrame) -> pd.DataFrame:
    from astropy.io import fits
    with fits.open(C.RAW_DIR / "quaia_G20.5.fits", memmap=True) as h:
        d = h[1].data
        q = pd.DataFrame({"source_id": np.asarray(d["source_id"], np.int64),
                          "z_quaia": np.asarray(d["redshift_quaia"], float),
                          "z_quaia_err": np.asarray(d["redshift_quaia_err"], float)})
    out = df.merge(q, on="source_id", how="left")
    assert len(out) == len(df)
    return out


def _consistency(rows: list[dict], comp: int = 0) -> dict:
    v = np.array([r["glide_galactic"][comp] for r in rows])
    e = np.array([r["glide_galactic_err"][comp] for r in rows])
    w = 1 / e**2
    mean = float(np.sum(w * v) / np.sum(w))
    chi2 = float(np.sum(((v - mean) / e) ** 2))
    from scipy.stats import chi2 as chi2d
    return {"weighted_mean": mean, "chi2": chi2, "dof": len(v) - 1,
            "p_value": float(chi2d.sf(chi2, len(v) - 1))}


def redshift_tests(dfq: pd.DataFrame, zp_l21: np.ndarray) -> dict:
    five = dfq["five_param"].to_numpy()
    z = dfq["z_quaia"].to_numpy()
    has = five & np.isfinite(z)
    out = {"n_with_quaia_z": int(has.sum())}
    for variant in ("raw", "conditioned_l21"):
        rows = []
        for lo, hi in zip(QUAIA_Z_EDGES[:-1], QUAIA_Z_EDGES[1:]):
            m = has & (z >= lo) & (z < hi)
            s = make_sample(dfq, m) if variant == "raw" else conditioned_sample(dfq, m, zp_l21)
            res = fit.fit_vsh(s, lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
            rows.append(_glide_row(res, f"z {lo}-{hi}", group=f"redshift (Quaia, {variant})",
                                   n_sample=int(m.sum()), z_median=float(np.median(z[m])),
                                   g_median=float(np.median(dfq.loc[m, "phot_g_mean_mag"])),
                                   nu_median=float(np.median(dfq.loc[m, "nu_eff_used_in_astrometry"]))))
        out[variant] = {"bins": rows, "consistency_gX": _consistency(rows, 0)}
        split = []
        for lo, hi in ((0.0, 1.25), (1.25, 5.0)):
            m = has & (z >= lo) & (z < hi)
            s = make_sample(dfq, m) if variant == "raw" else conditioned_sample(dfq, m, zp_l21)
            res = fit.fit_vsh(s, lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
            split.append(_glide_row(res, f"z {lo}-{hi}", n_sample=int(m.sum())))
        out[variant]["split_1.25"] = split
        out[variant]["split_consistency_gX"] = _consistency(split, 0)
    # The confounder check: magnitude bins within the low- and high-z halves.
    conf = []
    g = dfq["phot_g_mean_mag"].to_numpy()
    for zlo, zhi in ((0.0, 1.25), (1.25, 5.0)):
        for glo, ghi in ((0, 19.5), (19.5, 20.0), (20.0, 20.5)):
            m = has & (z >= zlo) & (z < zhi) & (g >= glo) & (g < ghi)
            res = fit.fit_vsh(make_sample(dfq, m), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE)
            conf.append(_glide_row(res, f"z {zlo}-{zhi}, G {glo}-{ghi}", n_sample=int(m.sum())))
    out["z_by_magnitude"] = conf
    return out


# --------------------------------------------------------------------------
# 9. The l = 1 rotation: a check on Gaia's reference frame
# --------------------------------------------------------------------------
def frame_rotation_tests(df: pd.DataFrame) -> dict:
    """Gaia set its frame's spin to zero using the proper motions of 428,034
    quasar-like sources (gaiadr3.frame_rotator_source). Fitted on exactly those,
    the rotation should vanish; on the other Gaia-CRF3 sources it need not."""
    from astropy.io import fits
    with fits.open(C.RAW_DIR / "frame_rotator_source.fits") as h:
        d = h[1].data
        used = np.asarray(d["source_id"], np.int64)[np.asarray(d["used_for_reference_frame_spin"]).astype(bool)]
    in_fr = np.isin(df["source_id"].to_numpy(), used)
    five = df["five_param"].to_numpy()
    out = {"n_frame_rotator_spin": int(used.size), "n_in_crf3": int(in_fr.sum()),
           "n_in_crf3_5p": int((in_fr & five).sum()), "n_in_crf3_6p": int((in_fr & ~five).sum())}
    rows = []
    for label, m in (("frame-rotator sources (5p)", in_fr & five),
                     ("other five-parameter sources", ~in_fr & five),
                     ("all five-parameter sources", five),
                     ("frame-rotator sources (5p + 6p)", in_fr),
                     ("all Gaia-CRF3 sources (5p + 6p)", np.ones(len(df), bool))):
        for lmax in (1, 3):
            res = fit.fit_vsh(make_sample(df, m), lmax=lmax, kappa=C.KAPPA_BASELINE)
            rows.append(_glide_row(res, label, n_sample=int(m.sum())))
    out["fits"] = rows
    # The Gaia-CRF3 paper's own check (Gaia Collaboration, Klioner et al. 2022,
    # Sect. 3.4): the frame rotator's algorithm -- spin only, no glide -- on all
    # sources, and on the frame rotator's own sources.
    spin_only = [t for t in vsh.terms(1) if t.kind == "T"]
    check = {}
    for label, m in (("all Gaia-CRF3 sources", np.ones(len(df), bool)), ("frame-rotator sources", in_fr)):
        res = fit.fit_vsh(make_sample(df, m), term_list=spin_only, kappa=C.KAPPA_BASELINE, sandwich=False)
        w, c = res.vector("T", "icrs", "formal")
        check[label] = {"n_accepted": res.n_used, "n_total": int(m.sum()), "omega": w.tolist(),
                        "omega_err": np.sqrt(np.diag(c)).tolist()}
    check["published_all"] = {"n_accepted": 1_607_289, "omega": [-3.44, 1.57, -1.24],
                              "omega_err": [0.30, 0.28, 0.32]}
    check["published_frame_rotator_abs"] = 0.8
    out["crf3_check"] = check
    g = df["phot_g_mean_mag"].to_numpy()
    mag = []
    for lo, hi in ((0, 18), (18, 19), (19, 19.5), (19.5, 20), (20, 20.5), (20.5, 22)):
        for name, m in (("frame rotator", in_fr & five), ("others", ~in_fr & five)):
            mm = m & (g >= lo) & (g < hi)
            if mm.sum() < 5000:
                continue
            res = fit.fit_vsh(make_sample(df, mm), lmax=1, kappa=C.KAPPA_BASELINE)
            mag.append(_glide_row(res, f"G {lo}-{hi}", group=name, n_sample=int(mm.sum()),
                                  g_median=float(np.median(g[mm]))))
    out["by_magnitude"] = mag
    return out


# --------------------------------------------------------------------------
# 10. Why the sample must be pure: stars that pass the same cuts
# --------------------------------------------------------------------------
def lookalike_tests(df: pd.DataFrame, fractions=(0.001, 0.003, 0.01, 0.03)) -> dict:
    """Glide of a 0.5% random sample of Gaia DR3 sources that pass quasar-like
    cuts (zero parallax and proper motion within 5 sigma) but are not in the
    Gaia-CRF3 -- mostly stars; then mix a known fraction of them into the
    quasars and watch the fitted acceleration move."""
    from astropy.table import Table
    t = Table.read(C.RAW_DIR / "lookalikes_halfpct.fits").to_pandas()
    t["five_param"] = True
    out = {"n": len(t), "g_median": float(np.median(t["phot_g_mean_mag"]))}
    s = make_sample(t)
    for lmax in (1, 3):
        res = fit.fit_vsh(s, lmax=lmax, kappa=C.KAPPA_BASELINE)
        out[f"fit_lmax{lmax}"] = _glide_row(res, f"look-alikes, l_max={lmax}", n_sample=len(t))
    five = df["five_param"].to_numpy()
    q = df.loc[five, ["ra", "dec", "pmra", "pmdec", "pmra_error", "pmdec_error", "pmra_pmdec_corr"]]
    rng = np.random.default_rng(C.SEED + 3)
    mix = []
    for f in fractions:
        k = int(round(f * len(q)))
        pick = t.iloc[rng.choice(len(t), size=k, replace=k > len(t))][q.columns]
        both = pd.concat([q, pick], ignore_index=True)
        both["five_param"] = True
        res = fit.fit_vsh(make_sample(both), lmax=C.LMAX_BASELINE, kappa=C.KAPPA_BASELINE, sandwich=False)
        row = _glide_row(res, f"{f:.1%} look-alikes mixed in", fraction=f, n_injected=k)
        row["injected_rejected_frac"] = float(res.outlier[len(q):].mean())
        mix.append(row)
    out["mixing"] = mix
    return out


def remove_nonglide(df: pd.DataFrame, mask, coefficients: dict, lmax: int = C.LMAX_BASELINE) -> pd.DataFrame:
    """A copy of df whose proper motions (rows in mask) have the fit's
    non-glide terms subtracted: the l = 1 rotation and every l = 2..lmax
    pattern -- Gaia's own large-scale systematics as the simultaneous fit
    sees them. What is left is the glide plus noise plus smaller-scale
    systematics. Used only for pictures; every number comes from the fits."""
    tl = vsh.terms(lmax)
    coef = np.array([coefficients[t.label] for t in tl])
    keep = np.array([0.0 if (t.l == 1 and t.kind == "S") else 1.0 for t in tl])
    m = np.asarray(mask)
    d = df[m]
    out = df.copy()
    pm_ra = np.empty(m.sum())
    pm_dec = np.empty(m.sum())
    ra, dec = np.radians(d["ra"].to_numpy()), np.radians(d["dec"].to_numpy())
    for sl in fit._chunks(m.sum()):
        b_ra, b_dec = vsh.basis(ra[sl], dec[sl], lmax, tl)
        pm_ra[sl] = b_ra @ (coef * keep)
        pm_dec[sl] = b_dec @ (coef * keep)
    out.loc[m, "pmra"] = d["pmra"].to_numpy() - pm_ra / MAS_TO_UAS
    out.loc[m, "pmdec"] = d["pmdec"].to_numpy() - pm_dec / MAS_TO_UAS
    return out


def nonglide_rms(df: pd.DataFrame, mask, coefficients: dict, lmax: int = C.LMAX_BASELINE) -> dict:
    """RMS over the quasars of the fitted non-glide field and of the glide."""
    tl = vsh.terms(lmax)
    coef = np.array([coefficients[t.label] for t in tl])
    d = df[np.asarray(mask)]
    b_ra, b_dec = vsh.basis(np.radians(d["ra"].to_numpy()), np.radians(d["dec"].to_numpy()), lmax, tl)
    out = {}
    for name, sel in (("glide", lambda t: t.l == 1 and t.kind == "S"),
                      ("rotation", lambda t: t.l == 1 and t.kind == "T"),
                      ("l2_l3", lambda t: t.l >= 2),
                      ("all_nonglide", lambda t: not (t.l == 1 and t.kind == "S"))):
        k = np.array([1.0 if sel(t) else 0.0 for t in tl])
        out[name] = float(np.sqrt(np.mean((b_ra @ (coef * k)) ** 2 + (b_dec @ (coef * k)) ** 2)))
    return out


def redshift_null_tests(dfq: pd.DataFrame, zp_l21: np.ndarray, n_shuffle: int = 20, split: float = 1.25) -> dict:
    """Is the low/high-redshift difference tied to redshift, or to where on the
    sky each bin's quasars happen to be? Shuffle the redshifts among quasars in
    the same ~7 degree sky pixel (HEALPix nside 8): the shuffled bins keep each
    bin's sky coverage but lose its redshifts. Also: the same split with more
    VSH terms (l_max = 6) and with parallax-conditioned proper motions."""
    import healpy as hp
    five = dfq["five_param"].to_numpy()
    z = dfq["z_quaia"].to_numpy()
    has = five & np.isfinite(z)

    def split_diff(zz, lmax=C.LMAX_BASELINE, conditioned=False):
        out = []
        for m in (has & (zz < split), has & (zz >= split)):
            s = conditioned_sample(dfq, m, zp_l21) if conditioned else make_sample(dfq, m)
            res = fit.fit_vsh(s, lmax=lmax, kappa=C.KAPPA_BASELINE, sandwich=False)
            v, c = res.vector("S", "galactic", "scaled")
            out.append((v[0], math.sqrt(c[0, 0])))
        (lo, elo), (hi, ehi) = out
        return {"low": lo, "low_err": elo, "high": hi, "high_err": ehi,
                "diff": hi - lo, "diff_err": math.hypot(elo, ehi)}

    res = {"real": split_diff(z), "real_lmax6": split_diff(z, lmax=6),
           "real_conditioned_l21": split_diff(z, conditioned=True)}
    pix = hp.ang2pix(8, dfq["ra"].to_numpy(), dfq["dec"].to_numpy(), lonlat=True)
    rng = np.random.default_rng(C.SEED + 11)
    idx = np.flatnonzero(has)
    order = idx[np.argsort(pix[idx], kind="stable")]
    bounds = np.flatnonzero(np.diff(pix[order])) + 1
    groups = np.split(order, bounds)
    shuffled = []
    for k in range(n_shuffle):
        zz = z.copy()
        for g in groups:
            zz[g] = z[rng.permutation(g)]
        shuffled.append(split_diff(zz))
    d = np.array([s["diff"] for s in shuffled])
    res["shuffled"] = shuffled
    res["shuffled_diff_mean"] = float(d.mean())
    res["shuffled_diff_sd"] = float(d.std(ddof=1))
    res["shuffled_frac_ge_real"] = float(np.mean(np.abs(d) >= abs(res["real"]["diff"])))
    res["real_over_shuffled_sd"] = float(res["real"]["diff"] / d.std(ddof=1))
    return res


def extendedness_tests(dfq: pd.DataFrame, coefficients: dict) -> dict:
    """Nearby quasars are more often partly resolved (their host galaxies show):
    more astrometric excess noise and more BP/RP flux excess. Is the redshift
    trend really an extendedness trend? Fit compact and extended quasars
    separately, and split each by redshift."""
    five = dfq["five_param"].to_numpy()
    aen = dfq["astrometric_excess_noise"].to_numpy()
    bx = dfq["phot_bp_rp_excess_factor"].to_numpy()
    zq = dfq["z_quaia"].to_numpy()
    zs = dfq["z_dr16q"].to_numpy()
    compact = five & (aen == 0) & (bx < 1.25)
    extended = five & ((aen > 0) | (bx >= 1.25))
    out = {"n_compact": int(compact.sum()), "n_extended": int(extended.sum()),
           "definition": "compact: astrometric_excess_noise = 0 and phot_bp_rp_excess_factor < 1.25"}
    rows = []

    def add(label, m, group, lmax=C.LMAX_BASELINE, **extra):
        res = fit.fit_vsh(make_sample(dfq, m), lmax=lmax, kappa=C.KAPPA_BASELINE, sandwich=False)
        rows.append(_glide_row(res, label, group=group, n_sample=int(m.sum()), **extra))

    add("compact", compact, "extendedness")
    add("extended", extended, "extendedness")
    q = np.nanquantile(bx[five], [0, 1 / 3, 2 / 3, 1])
    for i in range(3):
        m = five & (bx >= q[i]) & (bx <= q[i + 1] if i == 2 else bx < q[i + 1])
        add(f"BP/RP excess {q[i]:.2f}-{q[i + 1]:.2f}", m, "BP/RP excess tercile")
    add("excess noise = 0", five & (aen == 0), "astrometric excess noise")
    add("excess noise > 0", five & (aen > 0), "astrometric excess noise")
    for name, sel in (("compact", compact), ("extended", extended)):
        for lo, hi in ((0.0, 1.25), (1.25, 5.0)):
            m = sel & np.isfinite(zq) & (zq >= lo) & (zq < hi)
            add(f"{name}, Quaia z {lo}-{hi}", m, "redshift split by extendedness",
                z_median=float(np.median(zq[m])))
    # Spectroscopic redshifts, SDSS footprint, all-sky patterns removed
    clean = remove_nonglide(dfq, five, coefficients)
    alone = [t for t in vsh.terms(1) if t.kind == "S"]
    for name, sel in (("compact", compact), ("extended", extended)):
        for lo, hi in ((0.0, 1.0), (1.0, 7.0)):
            m = sel & np.isfinite(zs) & (zs >= lo) & (zs < hi)
            res = fit.fit_vsh(make_sample(clean, m), term_list=alone, kappa=C.KAPPA_BASELINE, sandwich=False)
            rows.append(_glide_row(res, f"{name}, SDSS z {lo}-{hi}", group="SDSS redshift split by extendedness",
                                   n_sample=int(m.sum()), z_median=float(np.median(zs[m]))))
    out["fits"] = rows
    # what fraction of each redshift bin is extended
    frac = []
    for lo, hi in ((0.0, 0.5), (0.5, 1.0), (1.0, 1.25), (1.25, 1.6), (1.6, 2.0), (2.0, 2.5), (2.5, 5.0)):
        m = five & np.isfinite(zq) & (zq >= lo) & (zq < hi)
        frac.append({"lo": lo, "hi": hi, "n": int(m.sum()), "extended_frac": float(extended[m].mean())})
    out["extended_fraction_vs_z"] = frac
    return out
