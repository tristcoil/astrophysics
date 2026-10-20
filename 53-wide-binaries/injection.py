"""Injection tests: does each method recover a gravity law we put in?

Mock versions of Chae's (2023) 26,615-pair sample are built from our forward
model: every real pair keeps its separation s, separation direction, masses,
distances, proper-motion errors and systemic motion, but its relative proper
motion is replaced by that of a model binary of the same s and total mass,
drawn at chosen population parameters (eccentricities, hidden close
companions, gravity law). We then add

* Gaia-like measurement noise from the pair's own proper-motion errors;
* the pair's real perspective effect: the systemic velocity (mean proper
  motion plus radial velocity) projected along two slightly different lines
  of sight, which Chae's v_p = |Delta mu| d does not remove;
* Hwang et al.-style eccentricity estimates from the mock v-r angle
  (hwang.py), since Chae's method draws each pair's eccentricity from them.

Running Chae's full pipeline (chae.py) on a Newtonian mock shows whether it
reports an anomaly when there is none; on a MOND mock, whether it finds one.
"""

from __future__ import annotations

import numpy as np

import chae as CH
import forward as FW
import hwang as HW
import kinematics as KIN

K_PM = 4.740470463533348


def _cells(log_r, log_m, dr: float, dm: float):
    return (np.floor(log_r / dr).astype(np.int64) * 100003 + np.floor(log_m / dm).astype(np.int64))


def draw_matches(pool_logr, pool_logm, target_logr, target_logm, rng, dr: float = 0.02,
                 dm: float = 0.04):
    """For each target, the index of a random pool sample in the same
    (log r_sky, log M) cell (widening the cell where it is empty)."""
    out = np.full(len(target_logr), -1, np.int64)
    todo = np.arange(len(target_logr))
    for widen in (1, 2, 4, 8, 16):
        cp = _cells(pool_logr, pool_logm, dr * widen, dm * widen)
        order = np.argsort(cp, kind="stable")
        cs = cp[order]
        ct = _cells(target_logr[todo], target_logm[todo], dr * widen, dm * widen)
        lo = np.searchsorted(cs, ct, "left")
        hi = np.searchsorted(cs, ct, "right")
        ok = hi > lo
        pick = lo[ok] + (rng.random(ok.sum()) * (hi[ok] - lo[ok])).astype(np.int64)
        out[todo[ok]] = order[pick]
        todo = todo[~ok]
        if not len(todo):
            break
    return out


def perspective_dmu(d, rv_fill_sigma: float = 30.0, rng=None):
    """Spurious relative proper motion (mas/yr; star 2 minus star 1) from the
    pair's systemic motion seen along two lines of sight."""
    n1, ea1, ed1 = KIN.unit(d["RA1"], d["DEC1"])
    n2, ea2, ed2 = KIN.unit(d["RA2"], d["DEC2"])
    dist = d["dM"]
    mua = 0.5 * (d["mu1ra"] + d["mu2ra"])
    mud = 0.5 * (d["mu1dec"] + d["mu2dec"])
    nm = n1 + n2
    nm /= np.linalg.norm(nm, axis=1, keepdims=True)
    # systemic velocity: tangential part from the mean proper motion at the
    # midpoint (using star 1's local basis), radial part from the RVs
    rv1 = np.where(d["RV1"] > -1000, d["RV1"], np.nan)
    rv2 = np.where(d["RV2"] > -1000, d["RV2"], np.nan)
    e1 = np.where(d["RV1"] > -1000, d["RV1_err"], np.nan)
    e2 = np.where(d["RV2"] > -1000, d["RV2_err"], np.nan)
    rv, _ = KIN.systemic_rv(rv1, e1, rv2, e2)
    if rng is not None:
        rv = np.where(np.isfinite(rv), rv, rv_fill_sigma * rng.standard_normal(len(rv)))
    else:
        rv = np.nan_to_num(rv)
    V = K_PM * (dist / 1000.0)[:, None] * (mua[:, None] * ea1 + mud[:, None] * ed1) + rv[:, None] * n1
    def pm(n, ea, ed):
        vt = V - np.sum(V * n, axis=1, keepdims=True) * n
        return (np.sum(vt * ea, axis=1) / (K_PM * dist / 1000.0), np.sum(vt * ed, axis=1) / (K_PM * dist / 1000.0))
    a1, d1 = pm(n1, ea1, ed1)
    a2, d2 = pm(n2, ea2, ed2)
    return a2 - a1, d2 - d1


def r_acceptance_table():
    """P(R_chance_align < 0.01) as a function of log s and v~ for pairs that
    pass all of Chae's other cuts, measured from El-Badry et al.'s catalogue:
    the chance-alignment ratio uses the proper-motion difference, so this cut
    preferentially removes fast pairs at wide separations (cached)."""
    from config import DERIVED_DIR
    path = DERIVED_DIR / "chae_R_acceptance.npz"
    if path.exists():
        z = np.load(path)
        return z["ls_edges"], z["vt_edges"], z["acc"]
    import elbadry
    import masses as ML
    eb = elbadry.load("near")
    g = lambda c: np.asarray(eb[c])
    bt = np.char.strip(g("binary_type").astype(str))
    d1, d2 = 1000 / g("parallax1"), 1000 / g("parallax2")
    s1, s2 = d1 * g("parallax_error1") / g("parallax1"), d2 * g("parallax_error2") / g("parallax2")
    rel = np.maximum.reduce([np.abs(g(f"{c}_error{w}") / g(f"{c}{w}")) for c in ("pmra", "pmdec") for w in (1, 2)])
    dM = KIN.weighted_distance(g("parallax1"), g("parallax_error1"), g("parallax2"), g("parallax_error2"))
    MG1 = g("phot_g_mean_mag1") - 5 * np.log10(dM / 10)
    MG2 = g("phot_g_mean_mag2") - 5 * np.log10(dM / 10)
    sk = g("sep_AU") / 1000
    base = ((bt == "MSMS") & (np.abs(d1 - d2) < 3 * np.hypot(s1, s2)) & (rel < 0.01) & (sk > 0.2) & (sk < 30)
            & (dM < 200) & (MG1 > 4) & (MG1 < 14) & (MG2 > 4) & (MG2 < 14))
    dmu = np.hypot(g("pmra1") - g("pmra2"), g("pmdec1") - g("pmdec2"))
    mt = ML.chae_mass(MG1) + ML.chae_mass(MG2)
    vt = K_PM * dmu * dM / 1000 / np.sqrt(KIN.GM_KMS2_AU * mt / (sk * 1000))
    ls_edges = np.log10([0.2, 1.0, 2.0, 3.5, 6.0, 10.0, 16.0, 30.0])
    vt_edges = np.array([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0])
    R = g("R_chance_align")
    num, _, _ = np.histogram2d(np.log10(sk[base & (R < 0.01)]), vt[base & (R < 0.01)], [ls_edges, vt_edges])
    den, _, _ = np.histogram2d(np.log10(sk[base]), vt[base], [ls_edges, vt_edges])
    acc = np.where(den > 10, num / np.maximum(den, 1), np.nan)
    # fill sparse cells from the nearest filled cell in v~ (same separation row)
    for i in range(acc.shape[0]):
        row = acc[i]
        ok = np.isfinite(row)
        if ok.any():
            acc[i] = np.interp(np.arange(len(row)), np.flatnonzero(ok), row[ok])
    np.savez(path, ls_edges=ls_edges, vt_edges=vt_edges, acc=acc)
    return ls_edges, vt_edges, acc


def elbadry_pm_cut(d, dmx, dmy):
    """El-Badry et al. (2021) bound-orbit cut: Delta mu < Delta mu_orbit(5 Msun)
    + 2 sigma, with their error propagation (find_binaries_edr3.py)."""
    theta = np.degrees(np.arccos(np.clip(np.sum(KIN.unit(d["RA1"], d["DEC1"])[0] * KIN.unit(d["RA2"], d["DEC2"])[0],
                                                axis=1), -1, 1))) * 3600.0
    plx1 = 1000.0 / d["d1"]
    dmu = np.hypot(dmx, dmy)
    sig = np.sqrt((dmx / dmu) ** 2 * (d["mu1ra_err"] ** 2 + d["mu2ra_err"] ** 2)
                  + (dmy / dmu) ** 2 * (d["mu1dec_err"] ** 2 + d["mu2dec_err"] ** 2))
    return dmu < 0.44428 * plx1**1.5 * theta**-0.5 + 2.0 * sig


def mock_sample(d: dict, theta7, rng, *, n_pool: int = 4_000_000, cb_mode: str = "relative",
                perspective: bool = True, noise: bool = True, alpha=None, selection: bool = True,
                e_alpha_fn=None, max_rounds: int = 40, err_scale: float = 1.0) -> dict:
    """A mock of Chae's sample ``d`` (chae.chae_sample output) at the forward
    model's parameters theta7. With ``selection``, each mock pair must also
    pass the cuts that shaped the real catalogue and that depend on the
    relative motion: El-Badry et al.'s proper-motion cut and the R < 0.01
    chance-alignment filter (acceptance measured from the catalogue);
    rejected pairs are redrawn."""
    mA, mB = d["M1"], d["M2"]
    pool = FW.simulate(theta7, mA, mB, n_pool, rng, a_lo=0.05, a_hi=100.0, cb_mode=cb_mode, alpha=alpha,
                       e_alpha_fn=e_alpha_fn)
    ok = pool["valid"] & (pool["r_sky_kau"] > 0)
    idx = np.flatnonzero(ok)
    n = len(d["s"])
    j = np.full(n, -1, np.int64)
    dmx = np.zeros(n)
    dmy = np.zeros(n)
    todo = np.arange(n)
    ls_e, vt_e, acc = r_acceptance_table() if selection else (None, None, None)
    if perspective:
        px_all, py_all = perspective_dmu(d, rng=rng)
    for rnd in range(max_rounds):
        sub = {k: (v[todo] if isinstance(v, np.ndarray) and v.shape[:1] == (n,) else v) for k, v in d.items()}
        jj = draw_matches(np.log10(pool["r_sky_kau"][idx]), np.log10(pool["m_est"][idx]),
                          np.log10(sub["s"]), np.log10(sub["M1"] + sub["M2"]), rng)
        jj = idx[jj]
        x, y = _mock_dmu(sub, pool, jj, rng, noise, err_scale)
        if perspective:
            x, y = x + px_all[todo], y + py_all[todo]
        if selection and rnd < max_rounds - 1:
            keep = elbadry_pm_cut(sub, x, y)
            vt = K_PM * np.hypot(x, y) * sub["dM"] / 1000 / np.sqrt(KIN.GM_KMS2_AU * (sub["M1"] + sub["M2"])
                                                                     / (sub["s"] * 1000))
            il = np.clip(np.searchsorted(ls_e, np.log10(sub["s"])) - 1, 0, len(ls_e) - 2)
            iv = np.clip(np.searchsorted(vt_e, vt) - 1, 0, len(vt_e) - 2)
            keep &= rng.random(len(todo)) < acc[il, iv]
        else:
            keep = np.ones(len(todo), bool)
        j[todo[keep]] = jj[keep]
        dmx[todo[keep]], dmy[todo[keep]] = x[keep], y[keep]
        todo = todo[~keep]
        if not len(todo):
            break
    m = dict(d)
    mua = 0.5 * (d["mu1ra"] + d["mu2ra"])
    mud = 0.5 * (d["mu1dec"] + d["mu2dec"])
    m["mu1ra"], m["mu2ra"] = mua - 0.5 * dmx, mua + 0.5 * dmx
    m["mu1dec"], m["mu2dec"] = mud - 0.5 * dmy, mud + 0.5 * dmy
    # Hwang-style eccentricities from the mock v-r angle
    g, ge, _ = HW.vr_angle(d["RA1"], d["DEC1"], d["RA2"], d["DEC2"], m["mu1ra"], m["mu1dec"], m["mu2ra"],
                           m["mu2dec"], d["mu1ra_err"], d["mu1dec_err"], d["mu2ra_err"], d["mu2dec_err"])
    e, e0, e1 = HW.posterior_summary(g, ge, d["s"] * 1000.0)
    m["e"], m["e0"], m["e1"] = e, e0, e1
    m["mock_true_e"] = pool["e"][j]
    m["mock_has_cb"] = pool["hasA"][j] | pool["hasB"][j]
    m["mock_rounds"] = rnd + 1
    return m


def _mock_dmu(d, pool, j, rng, noise, err_scale: float = 1.0):
    """Relative proper motion (mas/yr, star 2 minus star 1) of model binaries
    j placed at each real pair's separation, direction, distance and mass."""
    dra = (d["RA2"] - d["RA1"]) * np.cos(np.radians(d["DEC1"]))
    ddec = d["DEC2"] - d["DEC1"]
    phi_real = np.arctan2(ddec, dra)
    sd = pool["sep_dir"][j]
    phi_mod = np.arctan2(sd[:, 1], sd[:, 0])
    flip = rng.random(len(j)) < 0.5             # the handedness of the view is arbitrary
    vt = pool["vt_vec"][j].copy()
    vt[flip, 1] *= -1
    phi_mod = np.where(flip, -phi_mod, phi_mod)
    rot = phi_real - phi_mod
    vx = vt[:, 0] * np.cos(rot) - vt[:, 1] * np.sin(rot)
    vy = vt[:, 0] * np.sin(rot) + vt[:, 1] * np.cos(rot)
    vc = np.sqrt(KIN.GM_KMS2_AU * (d["M1"] + d["M2"]) / (d["s"] * 1000.0))
    dist = d["dM"]
    dmx = vx * vc / (K_PM * dist / 1000.0)
    dmy = vy * vc / (K_PM * dist / 1000.0)
    if noise:   # err_scale > 1: Gaia's quoted errors too small (the cuts still use the quoted ones)
        dmx = dmx + err_scale * np.hypot(d["mu1ra_err"], d["mu2ra_err"]) * rng.standard_normal(len(dmx))
        dmy = dmy + err_scale * np.hypot(d["mu1dec_err"], d["mu2dec_err"]) * rng.standard_normal(len(dmy))
    return dmx, dmy


def chae_on(d: dict, *, n_mc: int = 200, seed: int = 11, grid=None, echoice: str = "i", **kw):
    """Chae's calibration + 200-set Monte Carlo on any sample in his format
    (kw: mock_noise, mock_persp for the corrected method)."""
    f, gr, top = CH.calibrate_fmulti(d, grid=np.linspace(0.0, 1.0, 11) if grid is None else grid,
                                     n_mc=30, seed=seed, echoice=echoice, **kw)
    dl, mo, mn = CH.run(d, f, n_mc=n_mc, seed=seed + 1, echoice=echoice, **kw)
    med, sd = np.median(dl, axis=0), np.std(dl, axis=0)
    return dict(f_multi=f, calib_grid=gr, calib_top=top, delta=med, delta_sd=sd,
                gamma_g=10 ** (np.sqrt(2) * med), obs_medians=np.median(mo, 0), newton_medians=np.median(mn, 0),
                bin_counts=_bin_counts(d, f, seed))


def _bin_counts(d, f, seed):
    rng = np.random.default_rng(seed)
    m = CH.companion_masses(d["MA"], d["MB"], f, rng)
    mt = m["A1"] + m["A2"] + m["B1"] + m["B2"]
    gN, g = CH.deproject(d["s"], CH.observed_vp(d), mt, CH._draw_e(d["e"], d["e0"], d["e1"], rng), rng)
    x0 = 0.5 * (gN + g)
    return np.histogram(x0, bins=CH.X0_BINS)[0]


def vtilde_profile(d, edges=(0.2, 0.5, 1, 2, 3.5, 6, 12, 30), vmax: float = 2.0):
    """Median v~ (v~ < vmax) in separation bins, for all pairs and for the
    precisely measured ones (sigma_v~ < 0.05), plus the fraction above 1.5."""
    vp = CH.observed_vp(d)
    vc = np.sqrt(KIN.GM_KMS2_AU * (d["M1"] + d["M2"]) / (d["s"] * 1000))
    vt = vp / vc
    sig = (K_PM * np.hypot(np.hypot(d["mu1ra_err"], d["mu2ra_err"]), np.hypot(d["mu1dec_err"], d["mu2dec_err"]))
           / np.sqrt(2) * d["dM"] / 1000 / vc)
    out = dict(edges=list(edges), median=[], median_precise=[], n=[], n_precise=[], frac_gt15=[])
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (d["s"] > lo) & (d["s"] < hi)
        kv = k & (vt < vmax)
        kp = kv & (sig < 0.05)
        out["median"].append(float(np.median(vt[kv])))
        out["median_precise"].append(float(np.median(vt[kp])) if kp.sum() > 20 else None)
        out["n"].append(int(kv.sum()))
        out["n_precise"].append(int(kp.sum()))
        out["frac_gt15"].append(float(np.mean(vt[k] > 1.5)))
    return out
