"""Does gravity explain the flows? The 2MRS pull against the measured motions.

Two tests, both comparing where the mass is (2MRS) with how galaxies move (the
CMB dipole and Cosmicflows-4):

1. **The flux dipole** (Erdogdu et al. 2006). Gravity and light both fall off
   as 1/r^2, so if mass follows near-infrared light, the vector sum of the
   galaxies' fluxes points along the gravitational pull on the Local Group --
   no distances needed. Redshifts only sort the galaxies into spheres, showing
   how far one must look before the pull stops changing (convergence). In
   linear theory v = (H0 f / (b_L rho_L)) sum w_i S_i rhat_i, so comparing the
   amplitude with the Local Group's 620 km/s gives beta_L = f / b_L.

2. **Predicted against measured peculiar velocities** (Carrick et al. 2015).
   The reconstructed 2MRS density field predicts a velocity everywhere,
   beta* times the field v1 computed for beta* = 1. Each CF4 group then gets a
   likelihood for its observed redshift given its measured distance modulus,
   integrating over its unknown true distance along the line of sight (with
   the 2MRS density as the prior, which corrects inhomogeneous Malmquist bias).
   Free parameters: beta*, an external bulk flow V_ext from mass beyond the
   2MRS volume, a small-scale velocity scatter sigma_v, and the Hubble
   constant that matches CF4's distance scale.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import optimize

import catalogs
import config
import cosmology as cg
import density as dn
from config import C_KMS

OM = config.PLANCK18.Om0

#: Local Group members in 2MRS, excluded from the dipole (we compute the pull
#: *on* the Local Group). Named as in the 2MRS CAT column.
LOCAL_GROUP = {"MESSIER_031", "MESSIER_032", "MESSIER_033", "MESSIER_110", "NGC_0147", "NGC_0185", "IC_0010",
               "NGC_6822", "NGC_0205", "IC_1613"}


# --------------------------------------------------------------------------
# 1. The flux dipole
# --------------------------------------------------------------------------

def flux_dipole(lf: dict, *, frame: str = "lg", radii=None) -> dict:
    """Cumulative flux-weighted dipole v(<R)/beta_L in km/s, in Galactic Cartesian components.

    Every 2MRS galaxy with a redshift counts except Local Group members (nearby
    galaxies are kept: their fluxes are what they are; only their sphere is
    uncertain). The Milky Way's shadow is filled by cloning. Distances are
    redshift distances in the chosen frame ("lg" or "cmb"), floored at 1 h^-1 Mpc.
    """
    allg = catalogs.twomrs_with_groups()
    g = allg[allg["has_cz"] & (allg["Kcmag"] <= config.TWOMRS_KLIM)].copy()
    g = g[~g["CAT"].str.strip().isin(LOCAL_GROUP)]
    cz = g["cz_lg"] if frame == "lg" else g["cz_cmb"]
    z = np.maximum(cz.to_numpy(), 100.0) / C_KMS
    r = cg.comoving_distance(z, OM)
    g["r"] = r
    g["completeness_w"] = dn.completeness_weights(g)
    clones = dn.zoa_clones(g)
    r_all = np.concatenate([r, r[clones["source_index"].to_numpy()]])
    z_all = np.concatenate([z, z[clones["source_index"].to_numpy()]])
    n_all = np.vstack([g[["gx", "gy", "gz"]].to_numpy(), clones[["gx", "gy", "gz"]].to_numpy()])
    K_all = np.concatenate([g["Kcmag"].to_numpy(), clones["Kcmag"].to_numpy()])
    cw_all = np.concatenate([g["completeness_w"].to_numpy(), clones["completeness_w"].to_numpy()])
    # Flux in Lsun/Mpc^2 (independent of h): L / (4 pi D_L^2).
    S = 10 ** (-0.4 * (K_all - config.KCORR_COEFF * np.log10(1 + z_all) - dn.M_SUN_K - 25.0)) / (4 * np.pi)
    wl = dn.luminosity_weight(dn.limiting_magnitude(r_all, z_all), lf)
    terms = (wl * cw_all * S)[:, None] * n_all
    rho_L = mean_luminosity_density(lf)
    radii = np.arange(5.0, 201.0, 5.0) if radii is None else np.asarray(radii)
    order = np.argsort(r_all)
    cum = np.cumsum(terms[order], axis=0)
    cum2 = np.cumsum((terms[order] ** 2).sum(axis=1))
    rs = r_all[order]
    idx = np.searchsorted(rs, radii, side="right") - 1
    pref = 100.0 / rho_L  # H0 / rho_L with H0 = 100 h, rho_L in h Lsun Mpc^-3
    vec = pref * cum[np.maximum(idx, 0)]
    shot = pref * np.sqrt(cum2[np.maximum(idx, 0)])
    lg = catalogs.local_group_motion()
    ang = catalogs.angle_between(vec, lg["vector_gal"])
    l, b = catalogs.lonlat(vec)
    return {"R": radii.tolist(), "v_over_beta": vec.tolist(), "amp": np.linalg.norm(vec, axis=1).tolist(),
            "shot_noise_3d": shot.tolist(), "l": l.tolist(), "b": b.tolist(), "misalignment": ang.tolist(),
            "rho_L": rho_L, "frame": frame, "n_galaxies": int(len(g)), "n_clones": int(len(clones)),
            "lg_cmb_speed": lg["speed"]}


def mean_luminosity_density(lf: dict, rmin: float = 20.0, rmax: float = 150.0) -> float:
    """Mean Ks luminosity density (h Lsun Mpc^-3, isophotal magnitudes) in a shell
    rmin-rmax, from luminosity-weighted, ZoA-filled 2MRS galaxies (LG frame)."""
    g = dn.galaxies()
    z = g["cz_group"].to_numpy() / C_KMS
    r = cg.comoving_distance(z, OM)
    clones = dn.zoa_clones(g)
    src = clones["source_index"].to_numpy()
    r_all = np.concatenate([r, r[src]])
    K_all = np.concatenate([g["Kcmag"].to_numpy(), clones["Kcmag"].to_numpy()])
    cw = np.concatenate([g["completeness_w"].to_numpy(), clones["completeness_w"].to_numpy()])
    zc = cg.z_from_comoving(r_all, OM)
    L = 10 ** (-0.4 * (dn.absolute_magnitude(K_all, r_all, zc) - dn.M_SUN_K))
    wl = dn.luminosity_weight(dn.limiting_magnitude(r_all, zc), lf)
    m = (r_all > rmin) & (r_all < rmax)
    return float(np.sum((L * wl * cw)[m]) / (4 / 3 * np.pi * (rmax**3 - rmin**3)))


def beta_from_dipole(fd: dict, R_conv: float = 60.0) -> dict:
    """beta_L = |v_LG| / |v(<R)/beta_L| at a convergence radius."""
    i = int(np.argmin(np.abs(np.array(fd["R"]) - R_conv)))
    amp = fd["amp"][i]
    beta = fd["lg_cmb_speed"] / amp
    beta_err = beta * fd["shot_noise_3d"][i] / np.sqrt(3) / amp
    return {"R": fd["R"][i], "beta_L": beta, "beta_L_shot_err": beta_err, "misalignment": fd["misalignment"][i]}


# --------------------------------------------------------------------------
# 2. Velocity-velocity comparison
# --------------------------------------------------------------------------

R_GRID = np.linspace(0.5, 220.0, 300)  # full sightline, for the distance-only normalisation
LOCAL = np.linspace(-30.0, 30.0, 121)  # around the redshift distance, where the redshift term lives
H_GRID = np.linspace(0.66, 0.84, 73)


def cf4_sample(version: str = "corrected") -> pd.DataFrame:
    """CF4 groups for the comparison: inside VCOMP_RMAX in redshift space, no catastrophic outliers."""
    import flows

    g = flows.group_velocities(version)
    cz = g["V3k"]
    sel = g["map_ok"] & (cz > config.VCOMP_CZ_MIN) & (g["r_z"] < config.VCOMP_RMAX)
    return g[sel].reset_index(drop=True)


def _sample_along(n, rr, field3=None, field1=None):
    pos = n[:, None, :] * rr[..., None]
    flat = pos.reshape(-1, 3)
    inside = np.all(np.abs(flat) < config.RECON_HALF - 2, axis=1)
    out = np.zeros(len(flat))
    if field3 is not None:
        vv = dn.sample_field(field3, flat[inside])
        out[inside] = np.sum(vv * np.repeat(n, rr.shape[1], axis=0)[inside], axis=1)
    else:
        out[inside] = dn.sample_field(field1, flat[inside])
    return out.reshape(rr.shape)


def sightline_fields(sample: pd.DataFrame, delta: np.ndarray, v1: np.ndarray, data: dict) -> dict:
    """The radial velocity field (per unit beta*) and density along each sightline.

    The velocity (CMB frame, from linear theory) is needed only near the
    redshift distance; the density also along the whole sightline, for the
    normalisation of the distance prior, which is tabulated against h.
    """
    n = sample[["gx", "gy", "gz"]].to_numpy()
    rloc = np.clip(sample["r_z"].to_numpy()[:, None] + LOCAL[None, :], 0.5, None)
    nsel = data["radial_prior"]
    vr = _sample_along(n, rloc, field3=v1)
    dloc = _sample_along(n, rloc, field1=delta)
    dfull = _sample_along(n, np.broadcast_to(R_GRID, (len(n), len(R_GRID))).copy(), field1=delta)
    prior_full = nsel(R_GRID)[None, :] * np.clip(1 + dfull, 0.0, None)
    dm, e = data["DM"][:, None], data["e_DM"][:, None]
    lognorm = np.empty((len(n), len(H_GRID)))
    for j, h in enumerate(H_GRID):
        mu = data["mu_grid"][None, :] - 5 * np.log10(h)
        lognorm[:, j] = np.log(np.maximum((prior_full * np.exp(-0.5 * ((dm - mu) / e) ** 2)).sum(axis=1)
                                          * (R_GRID[1] - R_GRID[0]), 1e-300))
    zc = cg.z_from_comoving(rloc, OM)
    return {"vr": vr, "prior_loc": nsel(rloc) * np.clip(1 + dloc, 0.0, None), "rloc": rloc, "n": n,
            "zc_loc": zc, "mu_loc": 5 * np.log10((1 + zc) * rloc) + 25.0, "lognorm": lognorm}


def _loglike(params, data, fields, beta_fixed=None, return_terms=False):
    if beta_fixed is None:
        beta, vx, vy, vz, sig_v, h = params
    else:
        vx, vy, vz, sig_v, h = params
        beta = beta_fixed
    if sig_v <= 10 or h <= H_GRID[0] or h >= H_GRID[-1] or beta < 0:
        return -1e30
    n = fields["n"]
    u = beta * fields["vr"] + (n @ np.array([vx, vy, vz]))[:, None]
    czp = C_KMS * ((1 + fields["zc_loc"]) * (1 + u / C_KMS) - 1)
    mu = fields["mu_loc"] - 5 * np.log10(h)
    dm, e = data["DM"][:, None], data["e_DM"][:, None]
    w = fields["prior_loc"] * np.exp(-0.5 * ((dm - mu) / e) ** 2)
    lc = np.exp(-0.5 * ((data["cz"][:, None] - czp) / sig_v) ** 2) / (np.sqrt(2 * np.pi) * sig_v)
    num = (w * lc).sum(axis=1) * (LOCAL[1] - LOCAL[0])
    j = np.interp(h, H_GRID, np.arange(len(H_GRID)))
    j0 = int(np.floor(j))
    fr = j - j0
    lognorm = fields["lognorm"][:, j0] * (1 - fr) + fields["lognorm"][:, min(j0 + 1, len(H_GRID) - 1)] * fr
    ll = np.log(np.maximum(num, 1e-300)) - lognorm
    if return_terms:
        wl = w * lc
        r_mean = (wl * fields["rloc"]).sum(axis=1) / np.maximum(wl.sum(axis=1), 1e-300)
        u_at = (wl * u).sum(axis=1) / np.maximum(wl.sum(axis=1), 1e-300)
        return ll, r_mean, u_at
    return float(ll.sum())


def radial_prior(version: str = "corrected", bandwidth: float = 4.0):
    """How CF4 groups are distributed along a typical line of sight: the smoothed
    histogram of redshift distances of *all* CF4 groups (so it includes r^2,
    the mean density and CF4's selection, which an r^2 prior would ignore and
    so over-correct Malmquist bias). Returned as a function of r."""
    import flows

    r = flows.group_velocities(version).query("ok")["r_z"].to_numpy()
    grid = np.linspace(0.0, 400.0, 801)
    h, edges = np.histogram(r, bins=np.append(grid - 0.25, grid[-1] + 0.25))
    from scipy.ndimage import gaussian_filter1d

    hs = gaussian_filter1d(h.astype(float), bandwidth / 0.5) + 1e-3
    return lambda x: np.interp(x, grid, hs)


def prepare(sample: pd.DataFrame, version: str = "corrected") -> dict:
    z_cos = cg.z_from_comoving(R_GRID, OM)
    dl = (1 + z_cos) * R_GRID
    return {"cz": sample["V3k"].to_numpy(float), "DM": sample["DMzp"].to_numpy(float),
            "e_DM": sample["e_DMzp"].to_numpy(float), "mu_grid": 5 * np.log10(dl) + 25.0,
            "radial_prior": radial_prior(version)}


def fit_at(data, fields, beta_fixed=None, x0=None):
    """Maximum-likelihood parameters (beta_fixed=None: beta free)."""
    if x0 is None:
        x0 = [0.0, 0.0, 0.0, 250.0, 0.75] if beta_fixed is not None else [0.4, 0.0, 0.0, 0.0, 250.0, 0.75]
    f = lambda p: -_loglike(p, data, fields, beta_fixed)
    sol = optimize.minimize(f, x0, method="Powell", options={"xtol": 1e-3, "ftol": 1e-5, "maxfev": 4000})
    return sol


def velocity_comparison(rec: dict, version: str = "corrected", betas=None, verbose: bool = True) -> dict:
    """Profile likelihood in beta*, using the reconstruction made at each beta*."""
    sample = cf4_sample(version)
    data = prepare(sample, version)
    betas = [b for b in rec["betas"] if 0.2 <= b <= 0.7] if betas is None else betas
    _, cell = dn.grid_axes()
    prof = []
    x0 = None
    best = None
    for b in betas:
        delta = dn.field_for(rec, b)
        v1 = dn.velocity_field(delta, 1.0, cell)
        fields = sightline_fields(sample, delta, v1, data)
        sol = fit_at(data, fields, beta_fixed=b, x0=x0)
        x0 = sol.x
        prof.append({"beta": float(b), "lnL": float(-sol.fun), "Vext_gal": sol.x[:3].tolist(),
                     "sigma_v": float(sol.x[3]), "h": float(sol.x[4])})
        if verbose:
            print(f"    beta* = {b:.2f}: lnL {-sol.fun:.2f}  |Vext| {np.linalg.norm(sol.x[:3]):5.0f}  sigma_v {sol.x[3]:5.0f}  "
                  f"H0 {100 * sol.x[4]:.2f}", flush=True)
        if best is None or -sol.fun > best[0]:
            best = (-sol.fun, b, sol, fields, delta, v1)
    lnL = np.array([p["lnL"] for p in prof])
    bb = np.array([p["beta"] for p in prof])
    # Parabola through the profile near its peak for the 1-sigma interval.
    k = int(np.argmax(lnL))
    lo, hi = max(0, k - 3), min(len(bb), k + 4)
    c = np.polyfit(bb[lo:hi], lnL[lo:hi], 2)
    beta_hat = -c[1] / (2 * c[0])
    beta_err = 1 / np.sqrt(-2 * c[0])
    _, b_best, sol, fields, delta, v1 = best
    # Free-amplitude fit on the best reconstruction: is beta self-consistent?
    free = fit_at(data, fields, beta_fixed=None, x0=[b_best, *sol.x])
    vext = np.array(sol.x[:3])
    l, b_ = catalogs.lonlat(vext)
    # Predicted motion of the Local Group: internal field at the origin plus V_ext.
    v_int0 = b_best * dn.sample_field(v1, np.zeros((1, 3)))[0]
    v_lg = v_int0 + vext
    lg = catalogs.local_group_motion()
    ll, r_mean, u_pred = _loglike(sol.x, data, fields, b_best, return_terms=True)
    out = {
        "version": version, "n_groups": int(len(sample)), "profile": prof,
        "beta": float(beta_hat), "beta_err": float(beta_err), "beta_grid_best": float(b_best),
        "beta_free_fit": float(free.x[0]),
        "Vext": float(np.linalg.norm(vext)), "Vext_l": float(l), "Vext_b": float(b_), "Vext_gal": vext.tolist(),
        "sigma_v": float(sol.x[3]), "H0": float(100 * sol.x[4]),
        "v_lg_internal": v_int0.tolist(), "v_lg_pred": v_lg.tolist(), "v_lg_pred_speed": float(np.linalg.norm(v_lg)),
        "v_lg_pred_lb": [float(a) for a in catalogs.lonlat(v_lg)],
        "v_lg_misalignment": float(catalogs.angle_between(v_lg, lg["vector_gal"])),
        "v_lg_internal_speed": float(np.linalg.norm(v_int0)),
        "v_lg_internal_lb": [float(a) for a in catalogs.lonlat(v_int0)],
    }
    # Per-object comparison at the most likely distance: observed WF velocity vs prediction.
    out["_per_object"] = pd.DataFrame({"u_pred": u_pred, "v_obs": sample["v_wf"].to_numpy(),
                                       "e_v": sample["e_v"].to_numpy(), "r_z": sample["r_z"].to_numpy(),
                                       "sdss": sample["sdss_corrected"].to_numpy() if "sdss_corrected" in sample else False})
    out["_best"] = {"delta": delta, "v1": v1, "beta": b_best, "vext": vext}
    return out


def shell_decomposition(rec: dict, beta: float, vext: np.ndarray, shells=((0, 10), (10, 30), (30, 60), (60, 100), (100, 160))) -> list[dict]:
    """How much of the Local Group's predicted motion each shell of mass supplies."""
    delta = dn.field_for(rec, beta)
    R = dn._radius_grid()
    _, cell = dn.grid_axes()
    out = []
    for lo, hi in shells:
        part = np.where((R >= lo) & (R < hi), delta, 0.0)
        v = beta * dn.sample_field(dn.velocity_field(part, 1.0, cell), np.zeros((1, 3)))[0]
        l, b = catalogs.lonlat(v)
        out.append({"shell": [lo, hi], "v": v.tolist(), "speed": float(np.linalg.norm(v)), "l": float(l), "b": float(b)})
    l, b = catalogs.lonlat(vext)
    out.append({"shell": "external", "v": list(map(float, vext)), "speed": float(np.linalg.norm(vext)),
                "l": float(l), "b": float(b)})
    return out


def predicted_at_redshift_position(sample: pd.DataFrame, v1: np.ndarray, beta: float, vext: np.ndarray) -> np.ndarray:
    """Predicted radial peculiar velocity at each object's redshift-space position (km/s)."""
    n = sample[["gx", "gy", "gz"]].to_numpy()
    pos = n * sample["r_z"].to_numpy()[:, None]
    inside = np.all(np.abs(pos) < config.RECON_HALF - 2, axis=1)
    u = np.zeros(len(n))
    u[inside] = beta * np.sum(dn.sample_field(v1, pos[inside]) * n[inside], axis=1)
    return u + n @ np.asarray(vext)


def sdss_sign_test(rec: dict, beta: float, vext: np.ndarray) -> dict:
    """Which SDSS distances does 2MRS gravity agree with: CF4's as published, or sign-corrected?

    2MRS knows nothing of CF4's tables, so this is an independent test. For each
    version the full likelihood is refitted (V_ext, sigma_v, H0 free; beta* fixed)
    and the log-likelihood of the groups whose distances involve SDSS galaxies is
    compared; and, model-light, the correlation between each group's measured
    velocity and the 2MRS prediction at its position.
    """
    delta = dn.field_for(rec, beta)
    _, cell = dn.grid_axes()
    v1 = dn.velocity_field(delta, 1.0, cell)
    import cf4fix

    sdss_ids = set(cf4fix.corrected_groups().query("sdss_corrected")["1PGC"])
    out = {}
    for version in ("corrected", "published"):
        sample = cf4_sample(version)
        data = prepare(sample, version)
        fields = sightline_fields(sample, delta, v1, data)
        sol = fit_at(data, fields, beta_fixed=beta)
        ll, _, _ = _loglike(sol.x, data, fields, beta, return_terms=True)
        is_sdss = sample["1PGC"].isin(sdss_ids).to_numpy()
        upred = predicted_at_redshift_position(sample, v1, beta, sol.x[:3])
        vobs = sample["v_wf"].to_numpy()
        ok = np.isfinite(vobs)
        rs = np.corrcoef(vobs[ok & is_sdss], upred[ok & is_sdss])[0, 1]
        rn = np.corrcoef(vobs[ok & ~is_sdss], upred[ok & ~is_sdss])[0, 1]
        # Slope of measured on predicted (inverse-variance weighted), SDSS groups only.
        w = 1 / (sample["e_v"].to_numpy() ** 2 + sol.x[3] ** 2)
        m = ok & is_sdss
        slope = np.sum(w[m] * vobs[m] * upred[m]) / np.sum(w[m] * upred[m] ** 2)
        slope_err = 1 / np.sqrt(np.sum(w[m] * upred[m] ** 2))
        out[version] = {"lnL_total": float(ll.sum()), "lnL_sdss": float(ll[is_sdss].sum()),
                        "n_sdss": int(is_sdss.sum()), "n_other": int((~is_sdss).sum()),
                        "corr_sdss": float(rs), "corr_other": float(rn), "corr_err": float(1 / np.sqrt(m.sum())),
                        "slope_sdss": float(slope), "slope_sdss_err": float(slope_err),
                        "Vext": float(np.linalg.norm(sol.x[:3])), "sigma_v": float(sol.x[3]), "H0": float(100 * sol.x[4])}
    out["delta_lnL_sdss"] = out["corrected"]["lnL_sdss"] - out["published"]["lnL_sdss"]
    out["delta_lnL_total"] = out["corrected"]["lnL_total"] - out["published"]["lnL_total"]
    return out
