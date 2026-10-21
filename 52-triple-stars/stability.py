"""When does a triple fall apart? Chaos maps and a triple pushed past the edge.

Labelled simulations with the masses of the real system TIC 290061484
(Aa 6.85, Ab 6.11, B 7.90 Msun; Kostov et al. 2024):

* MEGNO maps (Cincotta & Simo 2000; REBOUND's variational equations): the
  Mean Exponential Growth of Nearby Orbits tends to 2 for regular motion and
  grows for chaotic motion. Grids in outer-to-inner semimajor-axis ratio and
  outer eccentricity (orbits nearly coplanar, as in TIC 290061484), and in
  ratio and mutual inclination. Mardling & Aarseth's (2001) criterion is
  drawn on top.
* Pushed past the limit: the same three stars with the outer orbit shrunk to
  a_out/a_in = 2.5 (below the limit), 400 random orbital phases: who gets
  ejected, how fast, and what is left behind.

Scale-free units: G = 1, inner semimajor axis 1, masses in Msun.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import time

import numpy as np
import rebound

import config as C

OUT = C.CACHE / "stability"
OUT.mkdir(parents=True, exist_ok=True)
MASSES = C.SYSTEMS["TIC 290061484"]["masses"]
M1, M2, M3 = MASSES["Aa"], MASSES["Ab"], MASSES["B"]


def build(ratio, e_out, i_mut_deg, seed, m=(M1, M2, M3), e_in=0.0):
    rng = np.random.default_rng(seed)
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.add(m=m[0])
    sim.add(m=m[1], a=1.0, e=e_in, inc=0.0, Omega=0.0, omega=rng.uniform(0, 2 * np.pi),
            M=rng.uniform(0, 2 * np.pi))
    sim.add(m=m[2], a=ratio, e=e_out, inc=np.radians(i_mut_deg), Omega=rng.uniform(0, 2 * np.pi),
            omega=rng.uniform(0, 2 * np.pi), M=rng.uniform(0, 2 * np.pi))
    sim.move_to_com()
    return sim


def megno_point(args):
    ratio, e_out, i_mut, seed, n_out = args
    sim = build(ratio, e_out, i_mut, seed)
    sim.integrator = "ias15"
    sim.init_megno()
    P_out = sim.particles[2].P
    t_end = n_out * P_out
    sim.exit_max_distance = 30.0 * ratio * (1 + e_out)
    try:
        sim.integrate(t_end, exact_finish_time=0)
        meg = sim.megno()
        t_surv = t_end
        escaped = False
    except rebound.Escape:
        meg = np.nan
        t_surv = sim.t
        escaped = True
    # did the inner pair keep its identity and size?
    o = sim.particles[1].orbit(primary=sim.particles[0])
    broken = escaped or not (0.5 < o.a < 2.0) or o.e > 0.99
    return float(meg) if np.isfinite(meg) else 99.0, float(t_surv / P_out), bool(broken)


def _init_worker():
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except ImportError:
        pass


def megno_map(kind="ratio_e", n_out=300, nproc=4, seed=C.RANDOM_SEED):
    """kind 'ratio_e': ratio 2..9 x e_out 0..0.8 at TIC 290061484's 9.2 deg;
    kind 'ratio_i': ratio 2..9 x mutual inclination 0..180 at e_out = 0.2."""
    path = OUT / f"megno_{kind}.npz"
    if path.exists():
        return dict(np.load(path))
    ratios = np.linspace(2.0, 9.0, 36)
    if kind == "ratio_e":
        ys = np.linspace(0.0, 0.8, 33)
        tasks = [(r, y, 9.17, seed + k, n_out) for k, (y, r) in enumerate((y, r) for y in ys for r in ratios)]
    else:
        ys = np.linspace(0.0, 180.0, 25)
        tasks = [(r, 0.2011, y, seed + k, n_out) for k, (y, r) in enumerate((y, r) for y in ys for r in ratios)]
    t0 = time.time()
    with mp.get_context("fork").Pool(nproc, initializer=_init_worker) as pool:
        res = pool.map(megno_point, tasks, chunksize=4)
    meg = np.array([r[0] for r in res]).reshape(len(ys), len(ratios))
    surv = np.array([r[1] for r in res]).reshape(len(ys), len(ratios))
    broken = np.array([r[2] for r in res]).reshape(len(ys), len(ratios))
    np.savez(path, ratios=ratios, ys=ys, megno=meg, survival=surv, broken=broken, n_out=n_out,
             seconds=time.time() - t0)
    return dict(np.load(path))


# ---------------------------------------------------------------------------
# Past the limit
# ---------------------------------------------------------------------------

def find_escaper(sim, ratio):
    """The star that has left, if any: far away (beyond 20 x the initial
    outer orbit), unbound to the other two, and those two bound to each other
    in a pair much tighter than the escaper's distance. Without the last two
    tests a member of the surviving binary, measured from the centre of mass
    of its partner and the true escaper, can pass for the escaper. Returns
    (index, energy per unit mass, orbit of the pair) or None."""
    ps = sim.particles
    cands = []
    for k in range(3):
        a_, b_ = [j for j in range(3) if j != k]
        o = ps[b_].orbit(primary=ps[a_])
        if not (0 < o.a and o.e < 1):
            continue
        mo = ps[a_].m + ps[b_].m
        rc = np.array([(ps[a_].m * getattr(ps[a_], c) + ps[b_].m * getattr(ps[b_], c)) / mo for c in ("x", "y", "z")])
        vc = np.array([(ps[a_].m * getattr(ps[a_], c) + ps[b_].m * getattr(ps[b_], c)) / mo
                       for c in ("vx", "vy", "vz")])
        r = np.array([ps[k].x, ps[k].y, ps[k].z]) - rc
        v = np.array([ps[k].vx, ps[k].vy, ps[k].vz]) - vc
        energy = 0.5 * v @ v - sim.G * (mo + ps[k].m) / np.linalg.norm(r)
        if energy > 0 and np.linalg.norm(r) > 20 * ratio and o.a < 0.1 * np.linalg.norm(r):
            cands.append((o.a, k, energy, o))
    if not cands:
        return None
    _, k, energy, o = min(cands, key=lambda c: c[0])     # the hardest pair stays
    return k, energy, o


def disrupt_one(args):
    ratio, e_out, i_mut, seed, t_max_out = args
    sim = build(ratio, e_out, i_mut, seed)
    sim.integrator = "ias15"
    P_out0 = sim.particles[2].P
    names = ["Aa", "Ab", "B"]
    t_end = t_max_out * P_out0
    dt_chk = 2.0 * P_out0
    t = 0.0
    while t < t_end:
        t += dt_chk
        sim.integrate(t, exact_finish_time=0)
        found = find_escaper(sim, ratio)
        if found:
            k, energy, o = found
            return dict(seed=seed, ejected=names[k], t_out_orbits=float(t / P_out0),
                        v_inf=float(np.sqrt(2 * energy)), a_bin=float(o.a), e_bin=float(o.e))
    return dict(seed=seed, ejected=None, t_out_orbits=float(t_end / P_out0))


def disruption_ensemble(ratio=2.5, e_out=0.2011, i_mut=9.17, n=400, t_max_out=2000, nproc=4,
                        seed=C.RANDOM_SEED + 7):
    path = OUT / f"disrupt_ratio{ratio:.2f}.json"
    if path.exists():
        return json.loads(path.read_text())
    tasks = [(ratio, e_out, i_mut, seed + k, t_max_out) for k in range(n)]
    t0 = time.time()
    with mp.get_context("fork").Pool(nproc, initializer=_init_worker) as pool:
        res = pool.map(disrupt_one, tasks, chunksize=2)
    out = dict(ratio=ratio, e_out=e_out, i_mut=i_mut, n=n, t_max_out=t_max_out, results=res,
               seconds=time.time() - t0)
    path.write_text(json.dumps(out))
    return out


def ejection_theory(m=(M1, M2, M3)):
    """Escape probabilities proportional to m^-3 (Valtonen & Karttunen 2006,
    'The Three-Body Problem', for an ergodic, equal-energy breakup)."""
    w = np.array(m, float) ** -3.0
    return dict(zip(["Aa", "Ab", "B"], (w / w.sum()).tolist()))


def track(ratio=2.5, e_out=0.2011, i_mut=9.17, seed=C.RANDOM_SEED + 7, n_out=200, per_out=200):
    """Positions of one disruption for the animation (labelled simulation)."""
    sim = build(ratio, e_out, i_mut, seed)
    sim.integrator = "ias15"
    P_out = sim.particles[2].P
    ts = np.linspace(0, n_out * P_out, n_out * per_out)
    xyz = np.empty((len(ts), 3, 3))
    for k, t in enumerate(ts):
        sim.integrate(t, exact_finish_time=1)
        for j in range(3):
            p = sim.particles[j]
            xyz[k, j] = (p.x, p.y, p.z)
    return ts / P_out, xyz


if __name__ == "__main__":
    import sys
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("all", "maps"):
        for kind in ("ratio_e", "ratio_i"):
            r = megno_map(kind)
            print(kind, "done", float(r["seconds"]), "s")
    if what in ("all", "disrupt"):
        r = disruption_ensemble()
        ej = [x["ejected"] for x in r["results"]]
        print({k: ej.count(k) for k in ("Aa", "Ab", "B", None)}, "theory", ejection_theory())
