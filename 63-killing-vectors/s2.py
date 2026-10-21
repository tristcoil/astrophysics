"""S2 as a test: are its energy and angular momentum conserved?

Schwarzschild spacetime has a time-translation Killing vector and three
rotations, so a star falling freely around Sgr A* keeps a constant energy E and
a constant angular-momentum vector L. Two tests with the real measurements
H1 fitted (Keck 1995-2018, NTT/VLT 1992-2016):

1. **Kepler's second law, season by season.** In each observing season we fit
   that season's positions alone (a quadratic in time) to get S2's position and
   velocity on the sky, relative to the black hole as located by H1's fit
   (including each team's reference-frame offset and drift). The sky-projected
   angular momentum h = X dY/dt - Y dX/dt (twice the areal velocity) must be
   the same in every season, although S2's speed changes 16-fold. This uses no
   orbit model, only where the mass is. Times are light-travel corrected
   (emission times); the uncorrected version is kept to show the effect.

2. **The energy, season by season.** E = v^2/2 - GM/r needs the full 3D speed
   and distance: the sky velocity from the same fits, the line-of-sight velocity
   from that season's spectra (minus H1's zero point and the relativistic
   redshift), and the line-of-sight position from the orbital plane. The plane
   itself is fixed because L's direction is conserved.

Plus the posterior of E and L over H1's 80,000 samples, in Newtonian and in
general-relativistic form, and the pericentre shift computed from E and L.

**What "general-relativistic form" means here.** The Killing energy and angular
momentum, E = (1 - 2GM/rc^2) c^2 dt/dtau and L = r^2 dphi/dtau, are exactly
conserved along a Schwarzschild geodesic. The Newtonian expressions computed
along the same orbit are not quite: near pericentre they drift by parts in
10^4 (L) to a few per cent of the binding energy (E). ``newton_vs_gr_along_orbit``
computes the size of that drift for S2, to set against the measurement errors.
"""

from __future__ import annotations

import numpy as np

import h1link
import orbits as ob

# Seasons need enough epochs and time span for a quadratic fit.
MIN_EPOCHS = 4
MIN_SPAN_YR = 0.15
#: Seasons this close to a pericentre are flagged: S2 turns by up to 180 degrees in
#: a few months there, which a one-season quadratic cannot follow.
PERICENTRE_FLAG_YR = 0.6


def _consts():
    op = h1link.modules().op
    return op, op.GM_SUN_AU3YR2, op.C_AUYR, op.KMS_PER_AUYR


def _h_model(gm, a, e, inc, r0, gm6):
    """Model sky-projected angular momentum X dY/dt - Y dX/dt in mas^2/yr."""
    op = h1link.modules().op
    k = op.precession_per_orbit(gm6, a, e) / (2 * np.pi)
    return -np.sqrt(gm * a * (1 - e**2)) * np.cos(np.radians(inc)) * (1 + k) / r0**2


def posterior_EL(thin: int = 1) -> dict:
    """E and L of S2 (per unit mass) for every sample of H1's joint GR chain."""
    op, gm_sun, c_auyr, kms = _consts()
    ch = h1link.chain()
    sl = slice(None, None, thin)
    gm6, period, e, inc = ch["gm6"][sl], ch["period"][sl], ch["e"][sl], ch["inc"][sl]
    r0 = ch["r0"][sl]
    gm = gm6 * 1e6 * gm_sun                       # AU^3/yr^2
    a = op.semi_major_axis_au(gm6, period)        # AU
    rg = gm / c_auyr**2                           # GM/c^2 in AU
    p = a * (1 - e**2) / rg                       # semi-latus rectum in GM/c^2
    E_gr, L_gr = ob.bound_orbit_EL(p, e)
    out = {
        "a_au": a, "rg_au": rg, "p": p, "e": e, "r0_kpc": r0, "inc": inc,
        "r_peri_au": a * (1 - e), "r_apo_au": a * (1 + e),
        # Newtonian: E = -GM/2a (km^2/s^2), L = sqrt(GM a (1-e^2)) (AU km/s)
        "E_newton_kms2": -gm / (2 * a) * kms**2,
        "L_newton_aukms": np.sqrt(gm * a * (1 - e**2)) * kms,
        # the same in geometric units: binding energy per c^2, L in GM/c
        "binding_newton": (1 - e**2) / (2 * p),
        "L_newton_geom": np.sqrt(p),
        # general relativity (exact Schwarzschild, Cutler, Kennefick & Poisson 1994)
        "E_gr": E_gr, "binding_gr": 1 - E_gr, "L_gr_geom": L_gr,
        "precession_exact_arcmin": np.degrees(ob.precession_exact(p, e)) * 60,
        "precession_6pi_arcmin": np.degrees(6 * np.pi / p) * 60,
        # sky-projected angular momentum, mas^2/yr: h = -|L| cos(i) (1 + k) / R0^2 in H1's
        # (East, North, away) frame (sign checked against H1's orbit_state), where
        # 1 + k is the areal-speed factor of H1's closed-form precession (k = 5e-4)
        "h_model_mas2yr": _h_model(gm, a, e, inc, r0, gm6),
    }
    return out


def summarise(x) -> dict:
    q = np.percentile(x, [16, 50, 84])
    return {"median": float(q[1]), "lo68": float(q[0]), "hi68": float(q[2]),
            "mean": float(np.mean(x)), "sd": float(np.std(x))}


def _orbit_plane_normal(inc_deg, Omega_deg, omega_deg):
    """Unit normal of the orbital plane in H1's (East, North, away) frame."""
    op = h1link.modules().op
    # positions a quarter orbit apart, from H1's own projection code
    vecs = []
    for nu in (0.0, 90.0):
        X, Y, Z, *_ = op.orbit_state(np.array([2000.0]), 4.0, 16.0, 2000.0, 1e-9, inc_deg,
                                     omega_deg + nu, Omega_deg)
        vecs.append(np.array([X[0], Y[0], Z[0]]))
    n = np.cross(vecs[0], vecs[1])
    return n / np.linalg.norm(n)


def _season_fit(t, X, Y, ex, ey, tc):
    """Weighted quadratic fits X(t), Y(t) about tc: position, velocity and their covariance."""
    A = np.vstack([np.ones_like(t), t - tc, 0.5 * (t - tc) ** 2]).T
    out = {}
    for name, v_, e_ in (("X", X, ex), ("Y", Y, ey)):
        w = 1 / e_**2
        cov = np.linalg.inv(A.T @ (A * w[:, None]))
        beta = cov @ (A.T @ (w * v_))
        chi2 = float(np.sum(w * (v_ - A @ beta) ** 2))
        dof = max(len(t) - 3, 1)
        scale = max(chi2 / dof, 1.0)  # inflate errors if the season scatters more than quoted
        out[name] = (beta, cov * scale, chi2, dof)
    return out


def seasonal_tests(n_draws: int = 400, seed: int = 4) -> dict:
    """Areal velocity (and energy) per observing season, over H1 posterior draws."""
    h = h1link.modules()
    op, gm_sun, c_auyr, kms = _consts()
    m = h1link.model()
    ch = h1link.chain()
    rng = np.random.default_rng(seed)
    idx = rng.choice(ch["_n"], size=n_draws, replace=False)
    draws = [ch["_best"]] + [{n_: ch[n_][i] for n_ in ch["_names"]} for i in idx]

    # all positions (both teams) with team labels
    teams = []
    for fr in m.frames:
        teams.append((fr.key, fr.t_ref, fr.data))
    rv = m.rv
    seasons = sorted({int(np.floor(tt)) for _, _, d in teams for tt in d.t})

    records = {s: {"h": [], "h_obs": [], "E": [], "speed": [], "tc": None, "n": 0, "teams": set(),
                   "h_err": None, "E_err": None, "n_rv": 0} for s in seasons}
    model_h, model_E = [], []
    for k, P in enumerate(draws):
        el = (P["gm6"], P["period"], P["t0"], P["e"], P["inc"], P["omega"], P["Omega"])
        gm = P["gm6"] * 1e6 * gm_sun
        a_ = op.semi_major_axis_au(P["gm6"], P["period"])
        model_E.append(-gm / (2 * a_) * kms**2)
        model_h.append(_h_model(gm, a_, P["e"], P["inc"], P["r0"], P["gm6"]))
        nvec = _orbit_plane_normal(P["inc"], P["Omega"], P["omega"])
        T, Xs, Ys, EX, EY, TEAM = [], [], [], [], [], []
        for key, tref, d in teams:
            dt = d.t - tref
            Xs.append(d.x - (P[f"x0_{key}"] + P[f"vx0_{key}"] * dt))
            Ys.append(d.y - (P[f"y0_{key}"] + P[f"vy0_{key}"] * dt))
            T.append(d.t)
            EX.append(d.x_err)
            EY.append(d.y_err)
            TEAM += [key] * d.n
        T, Xs, Ys, EX, EY = (np.concatenate(v_) for v_ in (T, Xs, Ys, EX, EY))
        TEAM = np.array(TEAM)
        T_em = op.roemer_emission_time(T, *el, f_sp=1.0)
        # velocities: subtract zero point, instrument offset and the GR redshift term
        rv_t_em = op.roemer_emission_time(rv.t, *el, f_sp=1.0)
        _, _, _, vxm, vym, vzm, rm = op.orbit_state(rv_t_em, *el, f_sp=1.0)
        redshift = op.relativistic_redshift_kms(P["gm6"], rm, np.sqrt(vxm**2 + vym**2 + vzm**2))
        v_kin = rv.v - P["vz0"] + np.where(rv.instrument == "NIRC2", P["nirc2_offset"], 0.0) - redshift
        for s in seasons:
            sel = np.floor(T).astype(int) == s
            if sel.sum() < MIN_EPOCHS or np.ptp(T[sel]) < MIN_SPAN_YR:
                continue
            w = 1 / (EX[sel] ** 2 + EY[sel] ** 2)
            tc = float(np.sum(w * T_em[sel]) / np.sum(w))
            fit = _season_fit(T_em[sel], Xs[sel], Ys[sel], EX[sel], EY[sel], tc)
            (bx, cx, _, _), (by, cy, _, _) = fit["X"], fit["Y"]
            X0, VX, Y0, VY = bx[0], bx[1], by[0], by[1]
            hval = X0 * VY - Y0 * VX
            # the same with observed (not light-travel corrected) times
            tco = float(np.sum(w * T[sel]) / np.sum(w))
            fo = _season_fit(T[sel], Xs[sel], Ys[sel], EX[sel], EY[sel], tco)
            hobs = fo["X"][0][0] * fo["Y"][0][1] - fo["Y"][0][0] * fo["X"][0][1]
            rec = records[s]
            rec["h"].append(hval)
            rec["h_obs"].append(hobs)
            # sky speed (km/s) for context
            rec["speed"].append(np.hypot(VX, VY) * P["r0"] * kms)
            if k == 0:  # statistical errors from the best-fit draw
                # h = X0 VY - Y0 VX; (X0, VX) and (Y0, VY) come from independent fits
                var = VY**2 * cx[0, 0] + Y0**2 * cx[1, 1] - 2 * VY * Y0 * cx[0, 1] \
                    + VX**2 * cy[0, 0] + X0**2 * cy[1, 1] - 2 * VX * X0 * cy[0, 1]
                rec["h_err"] = float(np.sqrt(var))
                rec["tc"] = tc
                rec["n"] = int(sel.sum())
                rec["teams"] = sorted(set(TEAM[sel]))
                rec["near_pericentre"] = bool(min(abs(tc - P["t0"]), abs(tc - (P["t0"] - P["period"])))
                                              < PERICENTRE_FLAG_YR)
            # energy: needs velocities in the same season
            rsel = np.floor(rv.t).astype(int) == s
            if rsel.sum() >= 2:
                tr = rv_t_em[rsel]
                wr = 1 / rv.v_err[rsel] ** 2
                Ar = np.vstack([np.ones_like(tr), tr - tc]).T
                covr = np.linalg.inv(Ar.T @ (Ar * wr[:, None]))
                br = covr @ (Ar.T @ (wr * v_kin[rsel]))
                vz = br[0]
                pos = np.array([X0, Y0]) * P["r0"]  # AU
                z_au = -(nvec[0] * pos[0] + nvec[1] * pos[1]) / nvec[2]
                rr = np.sqrt(pos @ pos + z_au**2)
                v2 = ((VX * P["r0"] * kms) ** 2 + (VY * P["r0"] * kms) ** 2 + vz**2)
                Eval = 0.5 * v2 - gm / rr * kms**2
                rec["E"].append(Eval)
                if k == 0:
                    # error: propagate the velocity uncertainties (positions matter less)
                    vxk, vyk = VX * P["r0"] * kms, VY * P["r0"] * kms
                    sig2 = (vxk * np.sqrt(cx[1, 1]) * P["r0"] * kms) ** 2 + \
                           (vyk * np.sqrt(cy[1, 1]) * P["r0"] * kms) ** 2 + (vz * np.sqrt(covr[0, 0])) ** 2
                    # distance error from the positional covariance
                    dr = np.sqrt(cx[0, 0] + cy[0, 0]) * P["r0"]
                    sig2 += (gm / rr**2 * kms**2 * dr) ** 2
                    rec["E_err"] = float(np.sqrt(sig2))
                    rec["n_rv"] = int(rsel.sum())

    model_h, model_E = np.array(model_h), np.array(model_E)
    table = []
    for s in seasons:
        rec = records[s]
        if not rec["h"]:
            continue
        hs = np.array(rec["h"])
        row = {"season": s, "t_center": rec["tc"], "n_epochs": rec["n"], "teams": rec["teams"],
               "near_pericentre": rec["near_pericentre"],
               "h_best": float(hs[0]), "h_stat_err": rec["h_err"],
               "h_frame_sd": float(np.std(hs[1:])), "h_obs_best": float(rec["h_obs"][0]),
               "sky_speed_kms": float(rec["speed"][0])}
        if rec["E"]:
            Es = np.array(rec["E"])
            row.update({"E_best": float(Es[0]), "E_stat_err": rec["E_err"], "E_frame_sd": float(np.std(Es[1:])),
                        "n_rv": rec["n_rv"]})
        table.append(row)

    def consistency(key, err_key, model):
        rows = [r_ for r_ in table if key in r_ and not r_["near_pericentre"]]
        vals = np.array([r_[key] for r_ in rows])
        errs = np.array([np.hypot(r_[err_key], r_[key.replace("_best", "_frame_sd")]) for r_ in rows])
        mu = float(np.median(model))
        chi2 = float(np.sum(((vals - mu) / errs) ** 2))
        wmean = float(np.sum(vals / errs**2) / np.sum(1 / errs**2))
        werr = float(1 / np.sqrt(np.sum(1 / errs**2)))
        return {"n": len(rows), "chi2_vs_model": chi2, "weighted_mean": wmean, "weighted_mean_err": werr,
                "model_median": mu, "model_sd": float(np.std(model)),
                "max_dev_sigma": float(np.max(np.abs(vals - mu) / errs)) if len(rows) else None}

    return {"seasons": table, "model_h": summarise(model_h), "model_E": summarise(model_E),
            "h_consistency": consistency("h_best", "h_stat_err", model_h),
            "E_consistency": consistency("E_best", "E_stat_err", model_E), "n_draws": n_draws}


def newton_vs_gr_along_orbit(n_points: int = 20001) -> dict:
    """Integrate S2's orbit as an exact Schwarzschild geodesic (H1's best-fit p, e)
    and compute along it the Killing quantities (constant) and the Newtonian
    expressions E_N = v^2/2 - GM/r, L_N = r^2 dphi/dt with coordinate velocities."""
    post = posterior_EL()
    ch = h1link.chain()
    best = ch["_best"]
    op, gm_sun, c_auyr, kms = _consts()
    gm = best["gm6"] * 1e6 * gm_sun
    a_ = op.semi_major_axis_au(best["gm6"], best["period"])
    rg = gm / c_auyr**2
    p, e = a_ * (1 - best["e"] ** 2) / rg, best["e"]
    E, L = ob.bound_orbit_EL(p, e)
    # integrate in the orbital plane with chi (relativistic anomaly): r = p / (1 + e cos chi)
    chi = np.linspace(-np.pi, np.pi, n_points)
    r = p / (1 + e * np.cos(chi))
    # dphi/dchi and dt/dchi for Schwarzschild bound orbits (Cutler et al. 1994)
    dphi = np.sqrt(p / (p - 6 - 2 * e * np.cos(chi)))
    dt = (p**2 / ((p - 2 - 2 * e * np.cos(chi)) * (1 + e * np.cos(chi)) ** 2)) * \
        np.sqrt((p - 2 - 2 * e) * (p - 2 + 2 * e) / (p - 6 - 2 * e * np.cos(chi)))
    t = np.concatenate([[0], np.cumsum(0.5 * (dt[1:] + dt[:-1]) * np.diff(chi))])
    phi = np.concatenate([[0], np.cumsum(0.5 * (dphi[1:] + dphi[:-1]) * np.diff(chi))])
    drdchi = p * e * np.sin(chi) / (1 + e * np.cos(chi)) ** 2
    rdot, phidot = drdchi / dt, dphi / dt               # d/dt (coordinate time), c = 1
    v2 = rdot**2 + (r * phidot) ** 2
    E_newton = 0.5 * v2 - 1 / r
    L_newton = r**2 * phidot
    # exact Killing quantities recomputed from the state (should be constant)
    f = 1 - 2 / r
    gamma = 1 / np.sqrt(f - rdot**2 / f - (r * phidot) ** 2)  # dt/dtau
    E_killing = f * gamma
    L_killing = r**2 * phidot * gamma
    t_yr = (t - t[n_points // 2]) * rg / c_auyr               # years from pericentre
    binding = (1 - e**2) / (2 * p)
    return {
        "t_yr": t_yr, "r_au": r * rg, "E_newton_rel": (E_newton + binding) / binding,
        "L_newton_rel": L_newton / np.sqrt(p) - 1,
        "E_killing_spread": float(np.ptp(E_killing)), "L_killing_spread": float(np.ptp(L_killing)),
        "E_killing": float(np.mean(E_killing)), "E_formula": float(E),
        "L_killing": float(np.mean(L_killing)), "L_formula": float(L),
        "E_newton_swing_frac": float(np.ptp(E_newton) / binding),
        "L_newton_swing_frac": float(np.ptp(L_newton) / np.mean(L_newton)),
        "p": float(p), "e": float(e), "rg_au": float(rg), "binding": float(binding),
        "period_check_yr": float(t[-1] * rg / c_auyr),
    }
