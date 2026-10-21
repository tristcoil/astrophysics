"""TIC 290061484 as a four-body system: Kostov et al. (2024)'s solution
integrated with REBOUND, its eclipse times compared with ours, and a fit of
the outer orbits to our eclipse times alone.

Kostov et al. give osculating Jacobi elements at t0 = BJD 2458711.0 for
Aa-Ab (1.79 d), A-B (24.6 d) and AB-C (3,205 d), with nodes relative to the
inner orbit's node. We place them in the sky frame (z toward the observer).
Their omega is the secondary's argument of periastron in the
radial-velocity convention; the sign convention is checked against their
published times of inferior conjunction (a test in this module).

Eclipse times: near each predicted eclipse we integrate to the predicted
time and take one Newton step on the sky-projected separation of Aa and Ab
(the relative motion is nearly uniform over minutes). The light-travel time
across the system is added: an eclipse seen when the pair's barycentre is
displaced by z toward us arrives z/c earlier.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import rebound

import config as C
import nbody as NB

SYS = C.SYSTEMS["TIC 290061484"]
DEG = np.pi / 180


#: REBOUND omega = paper omega - shift, per orbit (see module doc and
#: check_conventions; the inner orbit's shift is fixed by its published time
#: of inferior conjunction, the outer ones by our eclipse times)
OMEGA_SHIFT = {"Aa_Ab": 180.0, "A_B": 0.0, "AB_C": 180.0}


def build(params=None, omega_shift=None) -> rebound.Simulation:
    """Four-body system from Kostov et al.'s Table 3 (optionally modified).
    omega_shift: number or dict per orbit; REBOUND omega = omega_paper - shift."""
    if omega_shift is None:
        omega_shift = OMEGA_SHIFT
    if not isinstance(omega_shift, dict):
        omega_shift = {k: float(omega_shift) for k in ("Aa_Ab", "A_B", "AB_C")}
    p = {"Aa_Ab": dict(SYS["Aa_Ab"]), "A_B": dict(SYS["A_B"]), "AB_C": dict(SYS["AB_C"]),
         "masses": dict(SYS["masses"])}
    if params:
        for k, v in params.items():
            orb, key = k.split(".")
            p[orb][key] = v
    t0 = SYS["t0_bjd"]
    m = p["masses"]
    sim = NB.new_sim("day")
    sim.t = t0
    sim.add(m=m["Aa"])
    for orb, mass in (("Aa_Ab", m["Ab"]), ("A_B", m["B"]), ("AB_C", m["C"])):
        o = p[orb]
        M = 2 * np.pi * (t0 - o["tau"]) / o["P"]
        sim.add(m=mass, P=o["P"], e=o["e"], inc=o["i"] * DEG, Omega=o["Omega"] * DEG,
                omega=(o["omega"] - omega_shift[orb]) * DEG, M=M)
    sim.move_to_com()
    return sim


def _setup(sim):
    """Integrator for eclipse timing: IAS15 (adaptive, error at machine
    precision). Not WHFast: with safe_mode=0, every integrate() call that ends
    on a shortened step (exact_finish_time=1, as when sampling an eclipse)
    leaves a small error that accumulates - tens of minutes of eclipse timing
    by 2022 here - and at dt = P/80 even safe WHFast drifts ~6 min in 3 yr."""
    if sim.integrator != "ias15":
        sim.integrator = "ias15"


def _rel(sim, i, j):
    a, b = sim.particles[i], sim.particles[j]
    return (np.array([a.x - b.x, a.y - b.y, a.z - b.z]), np.array([a.vx - b.vx, a.vy - b.vy, a.vz - b.vz]))


def _bary(sim, idx):
    ps = sim.particles
    m = sum(ps[i].m for i in idx)
    r = sum(ps[i].m * np.array([ps[i].x, ps[i].y, ps[i].z]) for i in idx) / m
    return r


def sky_track(sim, front, back, t_start, t_end, step):
    """Sky separation and line-of-sight order of two groups of bodies
    (barycentres of `front` and `back` index lists) on a time grid."""
    ts = np.arange(t_start, t_end, step)
    sep, dz = [], []
    for t in ts:
        sim.integrate(t, exact_finish_time=1)
        d = _bary(sim, front) - _bary(sim, back)
        sep.append(np.hypot(d[0], d[1]))
        dz.append(d[2])
    return ts, np.array(sep), np.array(dz)


def inferior_conjunctions(ts, sep, dz):
    """Local minima of the sky separation with the 'front' group nearer us."""
    k = np.where((sep[1:-1] < sep[:-2]) & (sep[1:-1] < sep[2:]) & (dz[1:-1] > 0))[0] + 1
    return ts[k]


def eclipse_times(sim, t_near, nstep=160):
    """Model times of the inner eclipses nearest each time in t_near.

    The pair is followed through +-0.55 inner periods around each time; the
    minimum sky separation of Ab and Aa (either one in front) nearest t_near
    is refined with a parabola through the three lowest samples, and the
    light-travel time across the system is added (an eclipse seen when the
    pair's barycentre is z nearer to us arrives z/c earlier).

    Returns (times, kind) with kind = +1 when Ab is in front (primary)."""
    P = SYS["Aa_Ab"]["P"]
    _setup(sim)
    order = np.argsort(t_near)
    out = np.full(len(t_near), np.nan)
    kind = np.zeros(len(t_near))
    grid = np.linspace(-0.55, 0.55, nstep) * P
    for k in order:
        ts = t_near[k] + grid
        if ts[0] < sim.t:
            ts = ts[ts > sim.t]
        s2, zz, zc = [], [], []
        for t in ts:
            sim.integrate(t, exact_finish_time=1)
            r, _ = _rel(sim, 1, 0)
            s2.append(r[0] ** 2 + r[1] ** 2)
            zz.append(r[2])
            zc.append(_bary(sim, [0, 1])[2])
        s2, zz, zc = np.array(s2), np.array(zz), np.array(zc)
        mins = np.where((s2[1:-1] < s2[:-2]) & (s2[1:-1] < s2[2:]))[0] + 1
        if len(mins) == 0:
            continue
        j = mins[np.argmin(np.abs(ts[mins] - t_near[k]))]
        # parabola through three samples
        y0, y1, y2 = s2[j - 1], s2[j], s2[j + 1]
        h = ts[1] - ts[0]
        den = y0 - 2 * y1 + y2
        tmin = ts[j] + (0.5 * h * (y0 - y2) / den if den > 0 else 0.0)
        out[k] = tmin - zc[j] / C.C_AU_DAY
        kind[k] = 1.0 if zz[j] > 0 else -1.0
    return out, kind


def check_conventions():
    """Which omega convention puts Kostov et al.'s inferior conjunctions
    (Ab in front of Aa; B in front of the pair A) at their published times?"""
    res = {}
    for shift in (0.0, 180.0):
        sim = build(omega_shift=shift)
        sim.integrator = "ias15"
        ts, sep, dz = sky_track(sim, [1], [0], SYS["t0_bjd"], SYS["t0_bjd"] + 4.0, 0.002)
        t_in = inferior_conjunctions(ts, sep, dz)
        sim = build(omega_shift=shift)
        sim.integrator = "ias15"
        ts, sep, dz = sky_track(sim, [2], [0, 1], SYS["t0_bjd"], SYS["t0_bjd"] + 26.0, 0.01)
        t_out = inferior_conjunctions(ts, sep, dz)
        res[shift] = dict(inner=[float(x) for x in t_in], outer=[float(x) for x in t_out])
    return res


def model_oc(sim, obs: pd.DataFrame):
    """Model eclipse times at the observed eclipses; returns model times and
    the eclipse kind (+1 Ab in front)."""
    t, kind = eclipse_times(sim, obs["t_obs"].to_numpy())
    return t, kind


def compare(obs: pd.DataFrame, t_model: np.ndarray) -> dict:
    """Residuals after letting the linear ephemeris (T0, P) float: the model
    and our times may use slightly different mean periods."""
    cyc = obs["cycle"].to_numpy()
    d = obs["t_obs"].to_numpy() - t_model
    w = 1 / np.hypot(obs["err"].to_numpy(), 30 / 86400)
    A = np.vstack([np.ones_like(cyc), cyc]).T
    coef, *_ = np.linalg.lstsq(A * w[:, None], d * w, rcond=None)
    r = d - A @ coef
    return dict(rms_min=float(np.sqrt(np.mean(r**2)) * 1440), chi2=float(np.sum((r * w) ** 2)),
                n=int(len(r)), dT0_min=float(coef[0] * 1440), dP_s=float(coef[1] * 86400), resid=r)


def all_eclipses(sim, t_end, per_orbit=120):
    """Every inner conjunction from sim.t to t_end, in order: times (with the
    light-travel-time correction), kind (+1 Ab in front) and sky separation in
    units of the semimajor axis."""
    P = SYS["Aa_Ab"]["P"]
    _setup(sim)
    h = P / per_orbit
    ts = np.arange(sim.t + h, t_end, h)
    s2 = np.empty(len(ts))
    zz = np.empty(len(ts))
    zc = np.empty(len(ts))
    for k, t in enumerate(ts):
        sim.integrate(t, exact_finish_time=1)
        r, _ = _rel(sim, 1, 0)
        s2[k] = r[0] ** 2 + r[1] ** 2
        zz[k] = r[2]
        zc[k] = _bary(sim, [0, 1])[2]
    j = np.where((s2[1:-1] < s2[:-2]) & (s2[1:-1] < s2[2:]))[0] + 1
    y0, y1, y2 = s2[j - 1], s2[j], s2[j + 1]
    den = y0 - 2 * y1 + y2
    tmin = ts[j] + np.where(den > 0, 0.5 * h * (y0 - y2) / np.where(den > 0, den, 1.0), 0.0)
    a = SYS["Aa_Ab"]["a_rsun"] * C.R_SUN_AU
    return tmin - zc[j] / C.C_AU_DAY, np.sign(zz[j]), np.sqrt(y1) / a


def match_by_cycle(obs: pd.DataFrame, t_model, kind_model, T0_obs=None):
    """Model times at the observed cycle numbers: model primaries (Ab in
    front) are counted from the one nearest the observed cycle-0 primary;
    secondaries sit at half-integer cycles."""
    pri = t_model[kind_model > 0]
    sec = t_model[kind_model < 0]
    T0 = obs.loc[obs["cycle"] == obs["cycle"].min(), "t_obs"].iloc[0] - obs["cycle"].min() * 1.7927307 \
        if T0_obs is None else T0_obs
    k0 = int(np.argmin(np.abs(pri - T0)))
    out = np.full(len(obs), np.nan)
    for r, (c, kd) in enumerate(zip(obs["cycle"].to_numpy(), obs["kind"].to_numpy())):
        if kd == "pri":
            k = k0 + int(round(c))
            if 0 <= k < len(pri):
                out[r] = pri[k]
        else:
            # the secondary following primary floor(c)
            k = k0 + int(np.floor(c))
            if 0 <= k < len(pri):
                cand = sec[(sec > pri[k]) & (sec < pri[k] + 1.5)]
                if len(cand):
                    out[r] = cand[0]
    return out


def tertiary_events(sim, t_end, step=0.01):
    """Times when star B overlaps Aa or Ab on the sky (third-body eclipses),
    with who is in front. Radii from Kostov et al. (Rsun -> au)."""
    R = {k: v * C.R_SUN_AU for k, v in SYS["radii"].items()}
    _setup(sim)
    ts = np.arange(sim.t + step, t_end, step)
    rows = []
    prev = None
    for t in ts:
        sim.integrate(t, exact_finish_time=1)
        ps = sim.particles
        B = np.array([ps[2].x, ps[2].y, ps[2].z])
        over = []
        for i, nm in ((0, "Aa"), (1, "Ab")):
            S = np.array([ps[i].x, ps[i].y, ps[i].z])
            d = np.hypot(B[0] - S[0], B[1] - S[1])
            if d < R["B"] + R[nm]:
                over.append((nm, "B in front" if B[2] > S[2] else f"{nm} in front", d / (R["B"] + R[nm])))
        if over:
            rows.append((t, over))
    # group into events
    events = []
    for t, over in rows:
        if events and t - events[-1]["t_end"] < 0.3:
            events[-1]["t_end"] = t
            events[-1]["who"] |= {o[0] + ": " + o[1] for o in over}
            events[-1]["min_sep"] = min(events[-1]["min_sep"], min(o[2] for o in over))
        else:
            events.append(dict(t_start=t, t_end=t, who={o[0] + ": " + o[1] for o in over},
                               min_sep=min(o[2] for o in over)))
    for e in events:
        e["t_mid"] = 0.5 * (e["t_start"] + e["t_end"])
        e["who"] = sorted(e["who"])
    return events


def eclipse_times_fast(sim, t_guess):
    """Model eclipse times near guesses that are already within ~0.03 d:
    approach from 0.03 d before, two Newton steps on the sky separation, and
    the light-travel-time correction."""
    P = SYS["Aa_Ab"]["P"]
    _setup(sim)
    order = np.argsort(t_guess)
    out = np.full(len(t_guess), np.nan)
    kind = np.zeros(len(t_guess))
    for k in order:
        t = max(sim.t, t_guess[k] - 0.03)
        sim.integrate(t, exact_finish_time=1)
        r, v = _rel(sim, 1, 0)
        t1 = t - (r[0] * v[0] + r[1] * v[1]) / (v[0] ** 2 + v[1] ** 2)
        if t1 > sim.t:
            sim.integrate(t1, exact_finish_time=1)
            r, v = _rel(sim, 1, 0)
            t1 = sim.t - (r[0] * v[0] + r[1] * v[1]) / (v[0] ** 2 + v[1] ** 2)
        out[k] = t1 - _bary(sim, [0, 1])[2] / C.C_AU_DAY
        kind[k] = 1.0 if r[2] > 0 else -1.0
    return out, kind


#: parameters fitted to our eclipse times (Kostov et al.'s values start them)
FIT_PARAMS = ["Aa_Ab.P", "Aa_Ab.tau", "Aa_Ab.e", "Aa_Ab.omega",
              "A_B.P", "A_B.e", "A_B.omega", "A_B.tau", "A_B.Omega", "A_B.i",
              "masses.B",
              "AB_C.P", "AB_C.e", "AB_C.omega", "AB_C.tau", "masses.C"]
#: typical parameter scales (the inner period and phase are tightly pinned:
#: 1e-6 d of period moves the 2024 eclipses by ~1.5 min)
FIT_SCALES = {"Aa_Ab.P": 1e-6, "Aa_Ab.tau": 1e-4, "Aa_Ab.e": 0.001, "Aa_Ab.omega": 10.0,
              "A_B.P": 0.001, "A_B.e": 0.005, "A_B.omega": 1.0, "A_B.tau": 0.01, "A_B.Omega": 1.0,
              "A_B.i": 0.5, "masses.B": 0.1, "AB_C.P": 100.0, "AB_C.e": 0.05, "AB_C.omega": 10.0,
              "AB_C.tau": 50.0, "masses.C": 1.0}


def start_values():
    out = {}
    for k in FIT_PARAMS:
        orb, key = k.split(".")
        out[k] = float(SYS[orb][key]) if orb != "masses" else float(SYS["masses"][key])
    return out


def _mean_primary_period(params, t_end):
    sim = build(params)
    tm, kd, _ = all_eclipses(sim, t_end, per_orbit=60)
    pri = tm[kd > 0]
    return np.polyfit(np.arange(len(pri)), pri, 1)[0]


def calibrate_inner_period(obs, params, P_obs):
    """Adjust the inner orbit's osculating period and phase (tau) so that the
    integrated eclipses match ours in mean period and phase, all other
    elements held. Codes differ in how osculating elements map to mean
    motions, and the mapping here depends steeply on the inner phase relative
    to star B (short-period perturbations), so both are solved together: a
    two-parameter least-squares fit to the cycle-matched eclipse times,
    started from a proportional period rescaling."""
    from scipy.optimize import least_squares
    params = dict(params)
    Pm = _mean_primary_period(params, obs["t_obs"].max() + 2)
    params["Aa_Ab.P"] *= P_obs / Pm
    names = ["Aa_Ab.P", "Aa_Ab.tau"]

    def f(v):
        return residuals(v, obs, params, names=names)[:-1]
    x0 = np.array([params[k] for k in names])
    res = least_squares(f, x0, x_scale=[1e-5, 0.01], diff_step=None, method="lm",
                        jac=lambda v: np.array([(f(v + d) - f(v)) / d[d != 0][0]
                                                for d in (np.array([1e-7, 0.0]), np.array([0.0, 1e-4]))]).T)
    params.update(dict(zip(names, res.x.tolist())))
    return params


#: the variant with the fourth star held at Kostov et al.'s orbit: five years
#: of timing cannot constrain a 3,200-day orbit (the 16-parameter fit leaves
#: its mass at 19 +- 38 Msun)
FIT_PARAMS_TRIPLE = [k for k in FIT_PARAMS if not (k.startswith("AB_C") or k == "masses.C")]


def residuals(vec, obs, base, floor_days=20 / 86400, names=None):
    names = names or FIT_PARAMS
    params = dict(base)
    params.update({k: v for k, v in zip(names, vec)})
    try:
        sim = build(params)
        tm = model_times_cycles(sim, obs)
    except Exception:
        return np.full(len(obs) + 1, 1e3)
    bad = ~np.isfinite(tm)
    tm[bad] = obs["t_obs"].to_numpy()[bad] + 1.0
    pen = np.sum(bad) * 10.0
    err = np.hypot(obs["err"].to_numpy(), floor_days)
    return np.concatenate([(obs["t_obs"].to_numpy() - tm) / err, [pen]])


#: absolute finite-difference steps for the Jacobian (a relative step on a
#: Julian date would be hundreds of days)
FIT_STEPS = {k: 0.01 * v for k, v in FIT_SCALES.items()}


def _resid_task(args):
    vec, obs, base, names = args
    return residuals(vec, obs, base, names=names)


def fit_etv(obs: pd.DataFrame, P_obs: float, max_nfev=400, names=None, nproc=4, omega_shift=None):
    """Least-squares fit of the four-body model to the eclipse times, from
    Kostov et al.'s solution with the inner period calibrated; forward-
    difference Jacobian with absolute steps, evaluated in parallel."""
    import multiprocessing as mp
    from scipy.optimize import least_squares
    names = names or FIT_PARAMS
    base = calibrate_inner_period(obs, start_values(), P_obs)
    x0 = np.array([base[k] for k in names])
    scales = np.array([FIT_SCALES[k] for k in names])
    steps = np.array([FIT_STEPS[k] for k in names])
    r0 = residuals(x0, obs, base, names=names)
    pool = mp.get_context("fork").Pool(nproc)

    def jac(x, obs_, base_, names=None):
        tasks = [(x, obs_, base_, names)]
        for k in range(len(x)):
            xk = np.array(x, float)
            xk[k] += steps[k]
            tasks.append((xk, obs_, base_, names))
        f = pool.map(_resid_task, tasks)
        return np.array([(f[k + 1] - f[0]) / steps[k] for k in range(len(x))]).T

    try:
        res = least_squares(residuals, x0, jac=jac, x_scale=scales, args=(obs, base), kwargs=dict(names=names),
                            max_nfev=max_nfev, method="trf", verbose=2)
    finally:
        pool.close()
    J = res.jac
    dof = max(1, len(obs) - len(x0))
    chi2 = float(np.sum(res.fun[:-1] ** 2))
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, chi2 / dof)
        err = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        cov = None
        err = np.full(len(x0), np.nan)
    params = dict(base)
    params.update(dict(zip(names, res.x.tolist())))
    t_model = model_times_cycles(build(params), obs)
    return dict(start=dict(zip(names, x0.tolist())), best=dict(zip(names, res.x.tolist())),
                err=dict(zip(names, err.tolist())), chi2_start=float(np.sum(r0[:-1] ** 2)),
                chi2=chi2, dof=dof, nfev=int(res.nfev), njev=int(res.njev or 0), success=bool(res.success),
                message=str(res.message), base=base, names=names,
                cov=None if cov is None else cov.tolist(),
                t_obs=obs["t_obs"].tolist(), t_model=t_model.tolist(), integrator="ias15")


def model_times_cycles(sim, obs: pd.DataFrame, fine=120, coarse=12, pad=1.0):
    """Model time for every observed eclipse, matched by cycle number.

    The pair is followed continuously from t0: sampled every P/coarse between
    TESS sectors (enough to count every conjunction: the sky separation
    squared varies like sin^2 of the orbital phase) and every P/fine inside
    them (minimum refined with a parabola, light-travel time added). Model
    primaries (Ab in front) are numbered in order; the model primary nearest
    our first observed primary takes its cycle number, and every observed
    eclipse is paired with the model eclipse of the same cycle (secondaries:
    the one after the primary of the same integer cycle). Unlike pairing with
    the nearest model minimum, this stays smooth however far the model drifts."""
    _setup(sim)
    P = SYS["Aa_Ab"]["P"]
    wins = []
    for _, g in obs.groupby("sector"):
        wins.append([max(g["t_obs"].min() - pad, sim.t + P / fine), g["t_obs"].max() + pad])
    wins.sort()
    merged = [wins[0]]
    for w in wins[1:]:
        if w[0] <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], w[1])
        else:
            merged.append(w)
    grid, is_fine = [], []
    t = sim.t
    for t1, t2 in merged:
        tc = np.arange(t + P / coarse, t1 - P / coarse, P / coarse)
        tf = np.arange(t1, t2, P / fine)
        grid += [tc, tf]
        is_fine += [np.zeros(len(tc), bool), np.ones(len(tf), bool)]
        t = tf[-1]
    ts = np.concatenate(grid)
    fine_mask = np.concatenate(is_fine)
    s2 = np.empty(len(ts)); zz = np.empty(len(ts)); zc = np.empty(len(ts))
    for k, tt in enumerate(ts):
        sim.integrate(tt, exact_finish_time=1)
        r, _ = _rel(sim, 1, 0)
        s2[k] = r[0] ** 2 + r[1] ** 2
        zz[k] = r[2]
        zc[k] = _bary(sim, [0, 1])[2]
    j = np.where((s2[1:-1] < s2[:-2]) & (s2[1:-1] < s2[2:]))[0] + 1
    kind = np.sign(zz[j])
    tm = ts[j].copy()
    ok = fine_mask[j - 1] & fine_mask[j] & fine_mask[j + 1]
    y0, y1, y2 = s2[j - 1], s2[j], s2[j + 1]
    den = y0 - 2 * y1 + y2
    h = P / fine
    tm[ok] += np.where(den[ok] > 0, 0.5 * h * (y0[ok] - y2[ok]) / np.where(den[ok] > 0, den[ok], 1.0), 0.0)
    tm -= zc[j] / C.C_AU_DAY
    # number the model eclipses
    pri_idx = np.cumsum(kind > 0) - 1                 # index of the latest primary
    first = obs.loc[obs["kind"] == "pri"].iloc[0]
    pri_t = tm[kind > 0]
    k_first = int(np.argmin(np.abs(pri_t - first["t_obs"])))
    offset = k_first - int(round(first["cycle"]))
    out = np.full(len(obs), np.nan)
    pri_map = {i: tt for i, tt in enumerate(pri_t)}
    sec_map = {}
    for tt, kd, pi in zip(tm, kind, pri_idx):
        if kd < 0 and pi >= 0 and pi not in sec_map:
            sec_map[pi] = tt
    for r_i, (c, kd) in enumerate(zip(obs["cycle"].to_numpy(), obs["kind"].to_numpy())):
        k = int(np.floor(c)) + offset
        tt = pri_map.get(k) if kd == "pri" else sec_map.get(k)
        if tt is not None and (fine_mask[np.searchsorted(ts, tt) - 1] if 0 < np.searchsorted(ts, tt) <= len(ts) else False):
            out[r_i] = tt
    return out


def model_times_sectorwise(sim, obs: pd.DataFrame, per_orbit=120, pad=1.0):
    """Model eclipse time for every observed eclipse: each TESS sector is
    sampled continuously, every model minimum found (with which star is in
    front), and each observed eclipse paired with the nearest model minimum
    of the same kind. Robust to slow drifts between model and data."""
    P = SYS["Aa_Ab"]["P"]
    _setup(sim)
    h = P / per_orbit
    out = np.full(len(obs), np.nan)
    a = SYS["Aa_Ab"]["a_rsun"] * C.R_SUN_AU
    for sec, g in obs.groupby("sector", sort=False):
        pass
    order = obs.groupby("sector")["t_obs"].min().sort_values().index
    for sec in order:
        g = obs[obs["sector"] == sec]
        t1, t2 = g["t_obs"].min() - pad, g["t_obs"].max() + pad
        if t1 > sim.t:
            sim.integrate(t1, exact_finish_time=1)
        ts = np.arange(max(t1, sim.t) + h, t2, h)
        s2 = np.empty(len(ts)); zz = np.empty(len(ts)); zc = np.empty(len(ts))
        for k, t in enumerate(ts):
            sim.integrate(t, exact_finish_time=1)
            r, _ = _rel(sim, 1, 0)
            s2[k] = r[0] ** 2 + r[1] ** 2
            zz[k] = r[2]
            zc[k] = _bary(sim, [0, 1])[2]
        j = np.where((s2[1:-1] < s2[:-2]) & (s2[1:-1] < s2[2:]))[0] + 1
        y0, y1, y2 = s2[j - 1], s2[j], s2[j + 1]
        den = y0 - 2 * y1 + y2
        tm = ts[j] + np.where(den > 0, 0.5 * h * (y0 - y2) / np.where(den > 0, den, 1.0), 0.0) - zc[j] / C.C_AU_DAY
        km = np.sign(zz[j])
        for r_i, row in g.iterrows():
            want = 1.0 if row["kind"] == "pri" else -1.0
            cand = tm[km == want]
            if len(cand):
                out[obs.index.get_loc(r_i)] = cand[np.argmin(np.abs(cand - row["t_obs"]))]
    return out


def run_fit(max_nfev=600, variant="all") -> dict:
    """Fit the four-body model to our eclipse times (cached). variant 'all'
    frees the fourth star too; 'triple' holds it at Kostov et al.'s orbit."""
    out = C.CACHE / ("tic290061484_etv_fit.json" if variant == "all" else f"tic290061484_etv_fit_{variant}.json")
    if out.exists():
        return json.loads(out.read_text())
    obs = pd.read_csv(C.DERIVED / "tic290061484" / "eclipse_times.csv")
    obs = obs[obs["good"]].sort_values("t_obs").reset_index(drop=True)
    eph = json.loads((C.DERIVED / "tic290061484" / "ephemeris.json").read_text())
    res = fit_etv(obs, eph["P"], max_nfev=max_nfev,
                  names=FIT_PARAMS if variant == "all" else FIT_PARAMS_TRIPLE)
    res["variant"] = variant
    out.write_text(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import sys
    import time
    if len(sys.argv) > 1 and sys.argv[1] == "fit":
        t0 = time.time()
        r = run_fit(variant=sys.argv[2] if len(sys.argv) > 2 else "all")
        print(json.dumps({k: r[k] for k in ("chi2_start", "chi2", "dof", "nfev", "success")}, indent=1))
        print(f"{(time.time() - t0) / 60:.1f} min")


# ---------------------------------------------------------------------------
# Results for the README and figures
# ---------------------------------------------------------------------------

#: TESS seasons (pairs of consecutive sectors)
SEASONS = {"2019": (15, 16), "2022": (55, 56), "2024a": (75, 76), "2024b": (82, 83)}


def best_params(variant="all") -> dict:
    return dict(run_fit(variant=variant)["best"])


def observed_times() -> pd.DataFrame:
    obs = pd.read_csv(C.DERIVED / "tic290061484" / "eclipse_times.csv")
    return obs[obs["good"]].sort_values("t_obs").reset_index(drop=True)


def model_minima(params, windows, per_orbit=120):
    """Every model inner eclipse inside each (t1, t2) window (ascending):
    times (light-time corrected) and kind (+1 Ab in front)."""
    sim = build(params)
    out_t, out_k = [], []
    for t1, t2 in windows:
        _setup(sim)
        if t1 > sim.t:
            sim.integrate(t1, exact_finish_time=1)
        t, k, _ = all_eclipses(sim, t2, per_orbit=per_orbit)
        out_t.append(t)
        out_k.append(k)
    return np.concatenate(out_t), np.concatenate(out_k)


def outer_mean_periods(params=None, t_end=2460650.0, step=0.02) -> dict:
    """Mean outer periods of a model over the TESS years: the anomalistic
    period (B's periastron passages about the A pair) and the period of B's
    conjunctions with the pair as seen from Earth (linear fits to the event
    times). The osculating period at t0 differs because the inner pair's
    quadrupole speeds up the outer orbit's apsidal and nodal motion."""
    sim = build(params)
    _setup(sim)
    ts = np.arange(sim.t + step, t_end, step)
    r3 = np.empty(len(ts)); sep = np.empty(len(ts)); dz = np.empty(len(ts))
    for k, t in enumerate(ts):
        sim.integrate(t, exact_finish_time=1)
        d = np.array([sim.particles[2].x, sim.particles[2].y, sim.particles[2].z]) - _bary(sim, [0, 1])
        r3[k] = np.linalg.norm(d)
        sep[k] = np.hypot(d[0], d[1])
        dz[k] = d[2]

    def minima(y, mask=None):
        j = np.where((y[1:-1] < y[:-2]) & (y[1:-1] < y[2:]))[0] + 1
        if mask is not None:
            j = j[mask[j]]
        y0, y1, y2 = y[j - 1], y[j], y[j + 1]
        return ts[j] + 0.5 * step * (y0 - y2) / (y0 - 2 * y1 + y2)

    out = {}
    for key, t in (("anomalistic", minima(r3)), ("conjunction_B_in_front", minima(sep, dz > 0)),
                   ("conjunction_B_behind", minima(sep, dz < 0))):
        n = np.arange(len(t))
        A = np.vstack([np.ones_like(n), n]).T
        coef, *_ = np.linalg.lstsq(A, t, rcond=None)
        out[key] = dict(P=float(coef[1]), n_events=int(len(t)),
                        scatter_d=float(np.std(t - A @ coef)))
    out["osculating_at_t0"] = float((params or {}).get("A_B.P", SYS["A_B"]["P"]))
    return out


def etv_period_scan(obs=None, periods=np.linspace(12.0, 45.0, 1321), nharm=2) -> dict:
    """The outer period from the timing variations alone: in each season fit
    O-C = linear + nharm harmonics of period P (primaries and secondaries
    together), sum chi-square over seasons (phases free per season, so no
    cycle counting across the gaps)."""
    obs = observed_times() if obs is None else obs
    eph = json.loads((C.DERIVED / "tic290061484" / "ephemeris.json").read_text())
    oc = obs["t_obs"].to_numpy() - (eph["T0"] + obs["cycle"].to_numpy() * eph["P"])
    w = 1 / np.hypot(obs["err"].to_numpy(), 30 / 86400)
    chi = np.zeros(len(periods))
    for _, secs in SEASONS.items():
        sel = obs["sector"].isin(secs).to_numpy()
        t = obs["t_obs"].to_numpy()[sel]
        tt = t - t.mean()
        for k, P in enumerate(periods):
            cols = [np.ones_like(tt), tt]
            for h in range(1, nharm + 1):
                cols += [np.cos(2 * np.pi * h * tt / P), np.sin(2 * np.pi * h * tt / P)]
            A = np.vstack(cols).T
            coef, *_ = np.linalg.lstsq(A * w[sel, None], oc[sel] * w[sel], rcond=None)
            chi[k] += np.sum(((oc[sel] - A @ coef) * w[sel]) ** 2)
    k = int(np.argmin(chi))
    # error: delta chi2 = 1 after rescaling chi2 to chi2/dof = 1
    n_par = len(SEASONS) * (2 + 2 * nharm) + 1
    scale = chi[k] / (len(obs) - n_par)
    ok = (chi - chi[k]) / scale < 1
    return dict(periods=periods, chi2=chi, P_best=float(periods[k]),
                P_lo=float(periods[ok].min()), P_hi=float(periods[ok].max()), chi2_scale=float(scale))
