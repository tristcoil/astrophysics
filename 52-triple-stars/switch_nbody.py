"""N-body checks of the eclipse-switching fits (labelled simulations).

From a fitted parameter set we build the full three-body system with REBOUND
in the sky frame (z toward the observer):

* inner pair: the position and velocity of the secondary relative to the
  primary at the reference epoch, taken from the fitted (secular) model,
  scaled to the Keplerian semimajor axis;
* third star: the outer orbit's period and eccentricity (spectroscopy), its
  orientation on the opposite side of the precession cone (from the fit), its
  argument of periastron from the outer radial-velocity orbit and its mean
  anomaly from the outer periastron epoch (where published; otherwise a free
  phase).

We then integrate 1880-2300 (forwards and, with the velocities reversed,
backwards) and compare the inner orbit's tilt to the sky with the secular
model's, including the outer-period wiggles the secular model leaves out.
Optional REBOUNDx physics (GR, tidal and rotational bulges) for the
Kozai-Lidov question.
"""

from __future__ import annotations

import numpy as np
import rebound

import config as C
import nbody as NB
import switching as S

DEG = np.pi / 180
#: outer periastron epochs (JD) from the spectroscopic orbits (config.SYSTEMS)
T_PERI_OUT = {"SS Lac": C.SYSTEMS["SS Lac"]["T_peri_out"], "HS Hya": C.SYSTEMS["HS Hya"]["T_peri_out"],
              "V907 Sco": None}


def build(name, theta, M_out=None, omega_out=None):
    """REBOUND simulation (units au, days, Msun) at the fit's reference epoch."""
    model = S.Model(name)
    sp = model.spec
    p = model.unpack(theta)
    g = model.geometry(p)
    t0 = sp.t_ref
    h = 1e-3
    xyz = model.positions(g, np.array([t0 - h, t0, t0 + h]))
    a_in = g["a_in"]
    r = a_in * xyz[1]
    v = a_in * (xyz[2] - xyz[0]) / (2 * h)
    # rescale the velocity to the Keplerian speed for this separation and
    # mass (the secular model's mean motion differs by ~1e-5)
    mu = C.G_AU_DAY * (g["m1"] + g["m2"])
    vk2 = mu * (2 / np.linalg.norm(r) - 1 / a_in)
    v *= np.sqrt(vk2) / np.linalg.norm(v)
    sim = NB.new_sim("day")
    sim.t = t0
    NB.add_inner_from_state(sim, g["m1"], g["m2"], r, v)
    # outer orbit: opposite side of the precession cone
    hh = g["h0"] + np.pi
    eo = g["eps_out"]
    n_inv = np.array([np.sin(eo) * np.sin(hh), -np.sin(eo) * np.cos(hh), np.cos(eo)])
    I = g["I"]
    n_out = np.array([n_inv[0], np.cos(I) * n_inv[1] - np.sin(I) * n_inv[2],
                      np.sin(I) * n_inv[1] + np.cos(I) * n_inv[2]])
    inc, Om = NB.inc_node_from_normal(n_out)
    if omega_out is None:
        omega_out = sp.outer_rv.get("omega", 0.0)
    if M_out is None:
        tp = T_PERI_OUT.get(name)
        M_out = 0.0 if tp is None else 2 * np.pi * ((t0 - tp) / sp.P_out % 1.0)
    NB.add_third(sim, g["m3"], sp.P_out, sp.e_out, inc, Om, omega_out * DEG, M_out)
    sim.move_to_com()
    return sim, g


#: apsidal-motion constant k2 and squared gyration radius of the inner stars
#: (main-sequence values from Claret's tables, approximate)
BULGES = {"SS Lac": (0.004, 0.06), "HS Hya": (0.004, 0.07), "V907 Sco": (0.004, 0.06)}


def _reverse(sim):
    for p in sim.particles:
        p.vx, p.vy, p.vz = -p.vx, -p.vy, -p.vz


def tilt_history(name, theta, years=(1880.0, 2300.0), step_days=20.0, M_out=None, physics=None):
    """Inclination of the inner orbit to the sky (deg), mutual inclination
    and inner eccentricity on a time grid, by N-body."""
    out = {}
    for direction in (-1, 1):
        sim, g = build(name, theta, M_out=M_out)
        t0 = sim.t
        t_end = C.year_to_jd(years[0] if direction < 0 else years[1])
        if direction < 0:
            _reverse(sim)
        P_in = 2 * np.pi / g["n_lam"]
        if physics:
            import reboundx
            rebx = reboundx.Extras(sim)
            if "gr" in physics:
                f = rebx.load_force("gr_full")
                rebx.add_force(f)
                f.params["c"] = C.C_AU_DAY
            if "tides" in physics:
                # conservative tidal + rotational bulges of both inner stars,
                # spins synchronised with the orbit (REBOUNDx 'tides_spin')
                k2a, rg2 = BULGES[name]
                sp_ = S.SPECS[name]
                radii = (sp_.R1, sp_.R2) if sp_.R1 else (g["r1"] * g["a_in"] / C.R_SUN_AU,
                                                        g["r2"] * g["a_in"] / C.R_SUN_AU)
                f = rebx.load_force("tides_spin")
                rebx.add_force(f)
                hn = NB.inner_normal(sim)
                n_orb = np.sqrt(C.G_AU_DAY * (g["m1"] + g["m2"]) / g["a_in"] ** 3)
                for k in (0, 1):
                    pk = sim.particles[k]
                    pk.r = radii[k] * C.R_SUN_AU
                    pk.params["k2"] = 2 * k2a
                    pk.params["tau"] = 0.0
                    pk.params["I"] = rg2 * pk.m * pk.r**2
                    pk.params["Omega"] = rebound.Vec3d(*(direction * n_orb * hn))
                rebx.initialize_spin_ode(f)
            sim.integrator = "ias15"
        else:
            NB.setup_whfast(sim, P_in / 40.0)
        n = int(abs(t_end - t0) / step_days)
        ts = t0 + np.arange(1, n + 1) * step_days
        inc, J, e_in = [], [], []
        for t in ts:
            sim.integrate(t, exact_finish_time=0)
            hn = NB.inner_normal(sim)
            ho = NB.outer_normal(sim)
            # for the backwards run (velocities reversed) the normals flip
            if direction < 0:
                hn, ho = -hn, -ho
            inc.append(np.degrees(np.arccos(np.clip(hn[2], -1, 1))))
            J.append(np.degrees(np.arccos(np.clip(hn @ ho, -1, 1))))
            e_in.append(sim.particles[1].orbit(primary=sim.particles[0]).e)
        tt = t0 + direction * (ts - t0)
        out[direction] = (tt, np.array(inc), np.array(J), np.array(e_in))
    t = np.concatenate([out[-1][0][::-1], out[1][0]])
    o = [np.concatenate([out[-1][k][::-1], out[1][k]]) for k in (1, 2, 3)]
    return dict(jd=t, year=C.jd_to_year(t), i_sky=o[0], J=o[1], e_in=o[2])


def secular_tilt(name, theta, jd):
    model = S.Model(name)
    g = model.geometry(model.unpack(theta))
    return np.degrees(model.inclination(g, jd))


# ---------------------------------------------------------------------------
# Is a fitted geometry dynamically quiet? (Kozai-Lidov check on posterior draws)
# ---------------------------------------------------------------------------

#: a draw is "quiet" if the inner eccentricity rises by less than this above
#: its observed value over the century of data (Kozai-Lidov would pump it)
E_QUIET = 0.05


def is_quiet(name, check):
    return check["e_max"] < S.SPECS[name].e + E_QUIET


def _check_one(args):
    name, theta, M_out, years, physics = args
    try:
        h = tilt_history(name, theta, years=years, step_days=30.0, M_out=M_out, physics=physics)
    except Exception as exc:          # pragma: no cover - unusual geometries
        return dict(ok=False, error=str(exc))
    sec = secular_tilt(name, theta, h["jd"])
    return dict(ok=True, e_max=float(np.max(h["e_in"])), J_range=float(np.ptp(h["J"])),
                di_rms=float(np.sqrt(np.mean((h["i_sky"] - sec) ** 2))),
                di_max=float(np.max(np.abs(h["i_sky"] - sec))))


def classify(name, thetas, years=(1885.0, 2030.0), nproc=6, seed=C.RANDOM_SEED, physics=None):
    """N-body check of many parameter draws (outer phase random where no
    periastron epoch is published)."""
    import multiprocessing as mp
    rng = np.random.default_rng(seed)
    tasks = []
    for th in thetas:
        M_out = None if T_PERI_OUT.get(name) is not None else float(rng.uniform(0, 2 * np.pi))
        tasks.append((name, np.asarray(th, float), M_out, years, physics))
    with mp.get_context("fork").Pool(nproc) as pool:
        res = pool.map(_check_one, tasks, chunksize=2)
    for r, t in zip(res, tasks):
        r["M_out"] = t[2]
    return res


def run_checks(name, n=300, nproc=4, seed=C.RANDOM_SEED, tag=""):
    """N-body check of n posterior draws (cached): Kozai-quiet or not, and how
    far the N-body tilt history departs from the secular one."""
    import json
    import switch_fit as F
    out = F.CACHE / f"nbody_check_{F._slug(name)}{'_' + tag if tag else ''}.json"
    if out.exists():
        return json.loads(out.read_text())
    post = np.load(F.posterior_path(name, tag), allow_pickle=True)["posterior"]
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(post), min(n, len(post)), replace=False)
    res = classify(name, post[idx], nproc=nproc, seed=seed)
    for r, k in zip(res, idx):
        r["index"] = int(k)
    out.write_text(json.dumps(res))
    return res


if __name__ == "__main__":
    import sys
    r = run_checks(sys.argv[1], tag=sys.argv[2] if len(sys.argv) > 2 else "")
    ok = [x for x in r if x.get("ok")]
    quiet = [x for x in ok if is_quiet(sys.argv[1], x)]
    print(sys.argv[1], f"{len(quiet)}/{len(ok)} quiet; di_rms median "
          f"{np.median([x['di_rms'] for x in quiet]):.3f} deg, max {np.max([x['di_max'] for x in quiet]):.3f}")
