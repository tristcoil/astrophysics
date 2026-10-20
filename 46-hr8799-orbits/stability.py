"""N-body stability test of posterior orbits (REBOUND) - a labelled simulation.

    python stability.py W18/near1248_coplanar 100000     # draws from a fit's posterior

Set-up of Wang et al. (2018, Sec. 4.1), reproduced exactly:

* start on MJD 56609 (the first GPI epoch) with the fitted Jacobi elements;
  planets added e, d, c, b; star mass from the fit;
* masses: c, d and e equal, M_cde ~ U(4, 11) MJup; M_b ~ U(3 MJup, M_cde);
* velocities reversed, so integrating forwards runs the system *backwards*
  in time, for 40 Myr (the age of HR 8799), WHFast, fixed step = 1% of e's
  initial period;
* unstable if any two bodies come closer than the initial mutual Hill radius
  of d and e, R = a_e ((M_e + M_d) / 3 M*)^(1/3), or a planet goes beyond 500 au.

Every draw's outcome and survival time is appended to a CSV as it finishes,
so an interrupted run resumes where it stopped.
"""

from __future__ import annotations

import csv
import multiprocessing as mp
import os
import sys
import time

import numpy as np
import pandas as pd

import config as C
import kepler as K

FIELDS = ["draw", "sample", "m_cde", "m_b", "outcome", "t_end_yr", "cpu_s"]


def draw_masses(rng, n, prior: str = "wang"):
    """Planet masses (MJup).

    "wang": Wang et al.'s priors, M_cde ~ U(4, 11), M_b ~ U(3, M_cde).
    "hotstart": the luminosity (hot-start) masses Wang et al. used to re-weight
    their stable orbits, M_b ~ N(6.0, 0.7), M_cde ~ N(8.7, 1.0), truncated to
    the same support. Drawing from them directly samples the combined
    dynamics + luminosity posterior, instead of re-weighting a handful of
    high-mass survivors.
    """
    if prior == "wang":
        m_cde = rng.uniform(*C.MASS_CDE_RANGE, n)
        m_b = rng.uniform(C.MASS_B_MIN, m_cde)
        return m_cde, m_b
    if prior != "hotstart":
        raise KeyError(prior)
    m_cde = np.empty(n)
    m_b = np.empty(n)
    k = 0
    while k < n:
        c = rng.normal(*C.HOTSTART_MASS["cde"], 2 * n)
        b = rng.normal(*C.HOTSTART_MASS["b"], 2 * n)
        ok = (c >= C.MASS_CDE_RANGE[0]) & (c <= C.MASS_CDE_RANGE[1]) & (b >= C.MASS_B_MIN) & (b <= c)
        take = min(int(ok.sum()), n - k)
        m_cde[k:k + take] = c[ok][:take]
        m_b[k:k + take] = b[ok][:take]
        k += take
    return m_cde, m_b


def build_sim(elements, mstar, masses_mjup, reverse=True, t0=C.NBODY_EPOCH_MJD):
    """REBOUND simulation of the star + 4 planets at MJD t0 (sky frame, au, yr)."""
    import rebound
    masses = {p: masses_mjup[p] * C.MJUP_MSUN for p in C.INNER_TO_OUTER}
    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")
    sim.G = C.G_AU_YR
    sign = -1.0 if reverse else 1.0
    for m, x, y, z, vx, vy, vz in K.state_vectors(t0, elements, mstar, masses):
        sim.add(m=m, x=x, y=y, z=z, vx=sign * vx, vy=sign * vy, vz=sign * vz)
    sim.move_to_com()
    return sim


def hill_radius_de(elements, mstar, masses_mjup):
    a_e = elements["e"][0]
    mu = (masses_mjup["e"] + masses_mjup["d"]) * C.MJUP_MSUN / (3.0 * mstar)
    return a_e * mu ** (1.0 / 3.0)


def integrate(elements, mstar, masses_mjup, years=C.STABILITY_YEARS, reverse=True):
    """Run one system. Returns (outcome, time reached in yr)."""
    import rebound
    sim = build_sim(elements, mstar, masses_mjup, reverse=reverse)
    P_e = sim.particles[1].orbit(primary=sim.particles[0]).P
    sim.integrator = "whfast"
    sim.integrator.safe_mode = 0
    sim.dt = C.DT_FRACTION_OF_PE * P_e
    sim.exit_min_distance = hill_radius_de(elements, mstar, masses_mjup)
    sim.exit_max_distance = C.EJECTION_AU
    try:
        sim.integrate(years, exact_finish_time=0)
    except rebound.Encounter:
        return "encounter", sim.t
    except rebound.Escape:
        return "ejection", sim.t
    return "stable", sim.t


# ---------------------------------------------------------------------------
# Batch driver
# ---------------------------------------------------------------------------

_JOB = {}


def _elements_from_row(row: dict) -> dict:
    return {p: (row[f"a_{p}"], row[f"e_{p}"], np.radians(row[f"i_{p}"]), np.radians(row[f"omega_{p}"]),
                np.radians(row[f"Omega_{p}"]), row[f"tau_{p}"]) for p in C.INNER_TO_OUTER}


def _worker(task):
    draw, sample, m_cde, m_b = task
    row = _JOB["rows"][sample]
    masses = {"e": m_cde, "d": m_cde, "c": m_cde, "b": m_b}
    t0 = time.process_time()
    outcome, t_end = integrate(_elements_from_row(row), row["mstar"], masses, years=_JOB["years"])
    return dict(draw=draw, sample=sample, m_cde=m_cde, m_b=m_b, outcome=outcome,
                t_end_yr=t_end, cpu_s=time.process_time() - t0)


def _init(rows, years):
    _JOB["rows"] = rows
    _JOB["years"] = years
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def results_path(key: str, tag: str = ""):
    return C.CACHE / f"stability_{key.replace('/', '_')}{tag}.csv"


def read_results(path) -> pd.DataFrame:
    """Read a results CSV exactly. pandas' default float parser is off by one unit
    in the last place for ~10% of values; replaying a chaotic system with a mass
    that differs in the 16th digit is a different experiment (it broke up 15 of
    the 31 stable systems it touched). Hence round-trip parsing."""
    return pd.read_csv(path, float_precision="round_trip")


def run(key: str, n_draws: int, nproc: int = C.N_PROCESSES, years: float = C.STABILITY_YEARS,
        tag: str = "", seed_offset: int = 0, mass_prior: str = "wang") -> pd.DataFrame:
    """Integrate ``n_draws`` posterior draws of fit ``key``; resumable."""
    import fits
    fit = fits.load(key)
    table = fit.table()
    rows = table.to_dict("records")
    rng = np.random.default_rng([C.RANDOM_SEED, 4, seed_offset])
    sample = rng.integers(0, len(rows), n_draws)     # posterior draws (with replacement)
    m_cde, m_b = draw_masses(rng, n_draws, mass_prior)
    path = results_path(key, tag)
    done = set()
    if path.exists():
        done = set(pd.read_csv(path, usecols=["draw"]).draw)
    todo = [(i, int(sample[i]), float(m_cde[i]), float(m_b[i])) for i in range(n_draws) if i not in done]
    print(f"[{key}] {len(done)} done, {len(todo)} to go ({years / 1e6:.0f} Myr, {nproc} processes)")
    if todo:
        new = not path.exists()
        ctx = mp.get_context("fork")
        t0 = time.time()
        with open(path, "a", newline="") as f, ctx.Pool(nproc, initializer=_init, initargs=(rows, years)) as pool:
            w = csv.DictWriter(f, FIELDS)
            if new:
                w.writeheader()
            n_stable = 0
            for k, res in enumerate(pool.imap_unordered(_worker, todo, chunksize=4), start=1):
                w.writerow(res)
                n_stable += res["outcome"] == "stable"
                if k % 200 == 0 or k == len(todo):
                    f.flush()
                    os.fsync(f.fileno())
                    el = time.time() - t0
                    print(f"  {k}/{len(todo)}  stable so far {n_stable}  {el / 60:.1f} min "
                          f"(~{el / k * (len(todo) - k) / 60:.0f} min left)", flush=True)
    return read_results(path)


def pilot(key: str, n: int = 240, nproc: int = 8) -> pd.DataFrame:
    """A small first batch: how long does a draw take, and how many survive?"""
    t0 = time.time()
    df = run(key, n, nproc=nproc, tag="_pilot")
    wall = time.time() - t0
    st = df.outcome.eq("stable")
    print(f"[{key}] pilot: {st.sum()}/{len(df)} stable; median survival {df.t_end_yr.median():.3g} yr; "
          f"mean CPU per draw {df.cpu_s.mean():.2f} s (unstable {df.cpu_s[~st].mean():.2f} s); wall {wall / 60:.1f} min")
    return df


# ---------------------------------------------------------------------------
# Stable systems, second pass: resonant angles and the future
# ---------------------------------------------------------------------------

#: resonant angles (Wang et al. 2018, Sec. 4.3): for outer i and inner j,
#: theta_{i:j,k} = lambda_j - 2 lambda_i + varpi_k (k = i or j), the 1:2:4
#: Laplace angle lambda_e - 3 lambda_d + 2 lambda_c, and the four-body
#: lambda_e - 2 lambda_d - lambda_c + 2 lambda_b.
ANGLES = ["b:c,b", "b:c,c", "c:d,c", "c:d,d", "d:e,d", "d:e,e", "c:d:e", "b:c:d:e"]
WINDOW_YR = 100e3          # libration is judged in 100 kyr windows (500 snapshots)
GAP_DEG = 30.0             # librating if the angle leaves a gap > 30 deg on the circle


def _angles(lam, pom):
    """lam, pom: (n, 4) in the order e, d, c, b. Returns dict name -> (n,) radians."""
    e, d, c, b = 0, 1, 2, 3
    out = {}
    for name, (i, j) in (("b:c", (b, c)), ("c:d", (c, d)), ("d:e", (d, e))):
        for k, lab in ((i, name.split(":")[0]), (j, name.split(":")[1])):
            out[f"{name},{lab}"] = lam[:, j] - 2 * lam[:, i] + pom[:, k]
    out["c:d:e"] = lam[:, e] - 3 * lam[:, d] + 2 * lam[:, c]
    out["b:c:d:e"] = lam[:, e] - 2 * lam[:, d] - lam[:, c] + 2 * lam[:, b]
    return {k: np.mod(v, 2 * np.pi) for k, v in out.items()}


def libration_fraction(theta, n_window):
    """Fraction of windows in which the angle stays inside an arc < 360 - GAP_DEG."""
    n = len(theta) // n_window
    lib = 0
    for w in range(n):
        s = np.sort(theta[w * n_window:(w + 1) * n_window])
        gaps = np.diff(np.concatenate([s, [s[0] + 2 * np.pi]]))
        lib += np.degrees(gaps.max()) > GAP_DEG
    return lib / max(n, 1)


def _replay(row_draw, rows, reverse=True, years=C.STABILITY_YEARS, snapshot=None, nudge=0.0):
    """Rebuild one draw's system; integrate, optionally recording (lambda, varpi)."""
    import rebound
    draw = row_draw
    row = rows[int(draw["sample"])]
    masses = {"e": draw["m_cde"], "d": draw["m_cde"], "c": draw["m_cde"], "b": draw["m_b"]}
    el = _elements_from_row(row)
    sim = build_sim(el, row["mstar"], masses, reverse=reverse)
    P_e = sim.particles[1].orbit(primary=sim.particles[0]).P
    # No rotation into the invariable plane: the replay must start bit-for-bit as
    # in the first pass. The resonant angles are d'Alembert combinations (their
    # coefficients sum to zero), so for these nearly coplanar orbits they come out
    # the same with longitudes measured in the sky frame.
    sim.integrator = "whfast"
    sim.integrator.safe_mode = 0
    # reading the state for a snapshot must not touch the integration: without
    # this flag a 200-yr snapshot cadence alone moves planet b by ~70 au within
    # 0.2 Myr (the system is chaotic); with it the replay is bit-for-bit the
    # first pass
    sim.integrator.keep_unsynchronized = 1
    sim.dt = C.DT_FRACTION_OF_PE * P_e
    sim.exit_min_distance = hill_radius_de(el, row["mstar"], masses)
    sim.exit_max_distance = C.EJECTION_AU
    if nudge:
        sim.particles[4].x += nudge          # a perturbed twin: planet b moved by `nudge` au
    lam, pom = [], []
    try:
        if snapshot:
            for t in np.arange(snapshot, years + 1, snapshot):
                sim.integrate(t, exact_finish_time=0)
                o = sim.orbits()
                lam.append([o[k].l for k in range(4)])
                pom.append([o[k].pomega for k in range(4)])
        else:
            sim.integrate(years, exact_finish_time=0)
        outcome = "stable"
    except rebound.Encounter:
        outcome = "encounter"
    except rebound.Escape:
        outcome = "ejection"
    return outcome, sim.t, np.array(lam), np.array(pom)


def _resonance_worker(task):
    draw = task
    outcome, t_end, lam, pom = _replay(draw, _JOB["rows"], snapshot=C.SNAPSHOT_YEARS)
    res = dict(draw=int(draw["draw"]), outcome=outcome, t_end_yr=t_end, m_cde=draw["m_cde"], m_b=draw["m_b"])
    if outcome == "stable":
        th = _angles(lam, pom)
        nw = int(WINDOW_YR / C.SNAPSHOT_YEARS)
        for name, v in th.items():
            res[f"lib_{name}"] = libration_fraction(v, nw)
        # keep a thinned record (every 2,000 yr) of the angles for plots
        res["_series"] = {k: v[::10].astype(np.float32) for k, v in th.items()}
    return res


def _twin_worker(task):
    draw, nudge = task
    outcome, t_end, _, _ = _replay(draw, _JOB["rows"], nudge=nudge)
    return dict(draw=int(draw["draw"]), outcome=outcome, t_end_yr=t_end)


def twins(key: str, nudge_au: float = 1e-12, nproc: int = C.N_PROCESSES, tag: str = "") -> pd.DataFrame:
    """Is 'stable for 40 Myr' a property of the orbit or of one trajectory? Each
    stable system is run again with planet b moved by 1e-12 au."""
    import fits
    out = C.CACHE / f"twins_{key.replace('/', '_')}{tag}.csv"
    stable = stable_systems(key, tag)
    done = set(pd.read_csv(out).draw) if out.exists() else set()
    todo = [(r, nudge_au) for r in stable.to_dict("records") if int(r["draw"]) not in done]
    if todo:
        rows = fits.load(key).table().to_dict("records")
        ctx = mp.get_context("fork")
        new = not out.exists()
        with open(out, "a", newline="") as f, ctx.Pool(nproc, initializer=_init, initargs=(rows, C.STABILITY_YEARS)) as pool:
            w = csv.DictWriter(f, ["draw", "outcome", "t_end_yr"])
            if new:
                w.writeheader()
            for res in pool.imap_unordered(_twin_worker, todo):
                w.writerow(res)
                f.flush()
    return pd.read_csv(out).merge(stable[["draw", "m_cde", "m_b"]], on="draw")


def _future_worker(task):
    draw, years = task
    outcome, t_end, _, _ = _replay(draw, _JOB["rows"], reverse=False, years=years)
    return dict(draw=int(draw["draw"]), outcome=outcome, t_end_yr=t_end)


def stable_systems(key: str, tag: str = "") -> pd.DataFrame:
    df = read_results(results_path(key, tag))
    return df[df.outcome == "stable"].reset_index(drop=True)


def resonances(key: str, nproc: int = C.N_PROCESSES, tag: str = "") -> pd.DataFrame:
    """Libration fractions of the resonant angles of every stable system (cached)."""
    import fits
    out = C.CACHE / f"resonances_{key.replace('/', '_')}{tag}.npz"
    stable = stable_systems(key, tag)
    done = {}
    if out.exists():
        z = np.load(out, allow_pickle=True)
        done = {int(r["draw"]): r for r in z["records"]}
    todo = [r for r in stable.to_dict("records") if int(r["draw"]) not in done]
    if todo:
        rows = fits.load(key).table().to_dict("records")
        ctx = mp.get_context("fork")
        with ctx.Pool(nproc, initializer=_init, initargs=(rows, C.STABILITY_YEARS)) as pool:
            for k, res in enumerate(pool.imap_unordered(_resonance_worker, todo), start=1):
                done[res["draw"]] = res
                if k % 10 == 0 or k == len(todo):
                    np.savez_compressed(out, records=np.array(list(done.values()), dtype=object))
                    print(f"  resonances {k}/{len(todo)}", flush=True)
    recs = list(done.values())
    return pd.DataFrame([{k: v for k, v in r.items() if k != "_series"} for r in recs])


def future(key: str, years: float = 100e6, nproc: int = C.N_PROCESSES, tag: str = "") -> pd.DataFrame:
    """Integrate each stable system forwards from today for up to ``years``."""
    import fits
    out = C.CACHE / f"future_{key.replace('/', '_')}{tag}_{years / 1e6:.0f}Myr.csv"
    stable = stable_systems(key, tag)
    done = set(pd.read_csv(out).draw) if out.exists() else set()
    todo = [(r, years) for r in stable.to_dict("records") if int(r["draw"]) not in done]
    if todo:
        rows = fits.load(key).table().to_dict("records")
        ctx = mp.get_context("fork")
        new = not out.exists()
        with open(out, "a", newline="") as f, ctx.Pool(nproc, initializer=_init, initargs=(rows, years)) as pool:
            w = csv.DictWriter(f, ["draw", "outcome", "t_end_yr"])
            if new:
                w.writeheader()
            for k, res in enumerate(pool.imap_unordered(_future_worker, todo), start=1):
                w.writerow(res)
                f.flush()
    return pd.read_csv(out).merge(stable[["draw", "m_cde", "m_b"]], on="draw")


if __name__ == "__main__":
    key = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    run(key, n)
