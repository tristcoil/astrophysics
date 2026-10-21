"""The keeper animation: Mercury with general relativity turned up a million times.

A labelled **simulation**. Mercury's real relativistic precession is 0.1036
arcseconds per orbit, far too small to see. Here GR's extra pull is multiplied
by 10^6 (the speed of light slowed 1000 times, to 300 km/s) and the orbit
turns into a rosette.

At that strength GR's first-post-Newtonian equations are no longer a small
correction (Mercury would move at 16% of the slowed light speed), and
REBOUNDx's iterative 1PN force stops converging. We therefore use the
effective-potential form of GR's correction, V = -3 (GM)^2 / (c^2 r^2)
(Nobili & Roxburgh 1986; REBOUNDx ``gr_potential``), which gives exactly
GR's precession at real strength (our runs: 42.98068"/cy) and keeps the
orbit an exact turning ellipse at any strength: Newton's "revolving orbits"
(Principia, Book I, Prop. 44). The starting speed is chosen so that the
rosette keeps Mercury's real perihelion and aphelion distances.

    python animate.py         # mp4 + gif + final-frame png in figures/
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, PillowWriter

import ephemeris as eph
import orbit_physics as op
import plotstyle as PS
from config import ARCSEC_PER_RAD, FIGURE_DIR, J2000_JD, ROSETTE_BOOST, ROSETTE_ORBITS

PS.use_style()
FRAMES_PER_ORBIT = 30


def simulate(boost: float = ROSETTE_BOOST, n_orbits: int = ROSETTE_ORBITS,
             samples_per_orbit: int = 120) -> dict:
    import rebound
    import reboundx

    r, v = eph.heliocentric_state("mercury", J2000_JD)
    mu = eph.gm("sun") + eph.gm("mercury")
    el = op.elements(r, v, mu)
    a, e = float(el["a"][0]), float(el["e"][0])
    q, Q = a * (1 - e), a * (1 + e)
    c = eph.c_au_per_day() / np.sqrt(boost)
    beta = 3.0 * eph.gm("sun") ** 2 / c**2
    h = np.sqrt(2 * beta + 2 * mu * q * Q / (q + Q))
    sim = rebound.Simulation()
    sim.G = 1.0
    sim.add(m=eph.gm("sun"))
    sim.add(m=eph.gm("mercury"), x=q, y=0.0, z=0.0, vx=0.0, vy=h / q, vz=0.0)
    sim.move_to_com()
    sim.integrator = "ias15"
    rebx = reboundx.Extras(sim)
    force = rebx.load_force("gr_potential")
    rebx.add_force(force)
    force.params["c"] = c
    period = 2 * np.pi * np.sqrt(a**3 / mu)
    t = np.arange(0, n_orbits * period * 1.25, period / samples_per_orbit)
    xy = np.empty((len(t), 2))
    buf = np.zeros((2, 6))
    for i, ti in enumerate(t):
        sim.integrate(ti)
        sim.serialize_particle_data(xyzvxvyvz=buf)
        xy[i] = buf[1, :2] - buf[0, :2]
    rr = np.hypot(xy[:, 0], xy[:, 1])
    peri = np.flatnonzero((rr[1:-1] < rr[:-2]) & (rr[1:-1] <= rr[2:])) + 1
    keep = t <= t[peri[min(n_orbits, len(peri) - 1)]]
    peri = peri[peri < keep.sum()]
    ang = np.degrees(np.unwrap(np.arctan2(xy[peri, 1], xy[peri, 0])))
    # For V = -beta/r^2 the orbit is r = p / (1 + e cos(k phi)), k^2 = 1 - 2 beta / h^2; with h
    # chosen to keep Mercury's q and Q, 1/k = sqrt(1 + 6 GM / (c^2 p)) exactly.
    x = 6 * mu / (c**2 * a * (1 - e**2))
    return {"t": t[keep], "xy": xy[keep], "peri": peri, "turn_per_orbit_deg": float(np.diff(ang).mean()),
            "exact_deg": float(np.degrees(2 * np.pi * (np.sqrt(1 + x) - 1))),
            "linear_deg": float(np.degrees(np.pi * x)), "a": a, "e": e, "q": q, "Q": Q,
            "real_per_orbit_arcsec": float(op.schwarzschild_per_orbit(mu, a, e, eph.c_au_per_day())
                                           * ARCSEC_PER_RAD)}


def _figure(sim):
    km = eph.au_km() / 1e6
    fig = plt.figure(figsize=(8.4, 8.4), dpi=100)
    ax = fig.add_axes([0.06, 0.12, 0.88, 0.76])
    lim = sim["Q"] * km * 1.12
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.grid(False)
    # Newton: the closed ellipse with the same perihelion and aphelion.
    a, e = sim["a"], sim["e"]
    E = np.linspace(0, 2 * np.pi, 400)
    ax.plot(a * (np.cos(E) - e) * km, a * np.sqrt(1 - e**2) * np.sin(E) * km, color=PS.C_NEWTON,
            lw=1.2, alpha=0.8, zorder=1)
    ax.plot(0, 0, "o", color=PS.YELLOW, ms=13, zorder=5)
    fig.text(0.06, 0.955, "Mercury with general relativity turned up a million times",
             fontsize=14, fontweight="bold", color=PS.INK)
    fig.text(0.06, 0.92, f"Simulation. Real orbit: perihelion turns {sim['real_per_orbit_arcsec']:.4f}\" per lap. "
             f"Here: about {sim['turn_per_orbit_deg']:.0f} degrees per lap.", fontsize=9.5,
             color=PS.INK_SECONDARY)
    fig.text(0.06, 0.035, "Grey: Newton's closed ellipse. At real strength one full turn takes 3.0 million "
             "years (2,250 centuries with the planets' pull).\nGR's correction as the potential "
             "-3(GM)^2/(c^2 r^2) (Nobili & Roxburgh 1986), x 10^6; perihelion and aphelion kept at "
             "Mercury's real 46.0 and 69.8 million km.", fontsize=7.5, color=PS.MUTED)
    trail_old, = ax.plot([], [], color=PS.C_GR, lw=0.9, alpha=0.35, zorder=2)
    trail_new, = ax.plot([], [], color=PS.C_GR, lw=1.8, zorder=3)
    peris = ax.scatter([], [], s=14, color=PS.INK, zorder=4)
    planet, = ax.plot([], [], "o", color=PS.INK, ms=7, zorder=6)
    counter = ax.text(0.02, 0.97, "", transform=ax.transAxes, fontsize=10, color=PS.INK, va="top")
    return fig, (trail_old, trail_new, peris, planet, counter)


def make(stem: str = "mercury_rosette"):
    sim = simulate()
    km = eph.au_km() / 1e6
    xy = sim["xy"] * km
    step = 120 // FRAMES_PER_ORBIT
    frames = np.arange(0, len(xy), step)
    peri = sim["peri"]
    fig, (trail_old, trail_new, peris, planet, counter) = _figure(sim)
    lap_len = 120

    def draw(i):
        trail_old.set_data(xy[:i + 1, 0], xy[:i + 1, 1])
        j0 = max(0, i - lap_len)
        trail_new.set_data(xy[j0:i + 1, 0], xy[j0:i + 1, 1])
        done = peri[peri <= i]
        peris.set_offsets(xy[done] if len(done) else np.empty((0, 2)))
        planet.set_data([xy[i, 0]], [xy[i, 1]])
        laps = max(len(done) - 1, 0)
        counter.set_text(f"lap {min(laps + 1, ROSETTE_ORBITS)}   ·   perihelion has turned {laps * sim['turn_per_orbit_deg']:.0f} degrees")
        return trail_old, trail_new, peris, planet, counter

    mp4 = FIGURE_DIR / f"{stem}.mp4"
    w = FFMpegWriter(fps=30, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "27"])
    with w.saving(fig, str(mp4), dpi=100):
        for i in frames:
            draw(i)
            w.grab_frame()
    # Final frame as a still.
    draw(frames[-1])
    fig.savefig(FIGURE_DIR / f"{stem}.png", dpi=150)
    # A lighter GIF: the first full turn of the rosette, at 15 frames per second.
    n_turn = int(np.ceil(360 / sim["turn_per_orbit_deg"])) + 1
    gif_frames = frames[frames <= peri[min(n_turn, len(peri) - 1)]][::3]
    gif = FIGURE_DIR / f"{stem}.gif"
    gw = PillowWriter(fps=10)
    with gw.saving(fig, str(gif), dpi=50):
        for i in gif_frames:
            draw(i)
            gw.grab_frame()
    plt.close(fig)
    return {"mp4": str(mp4), "gif": str(gif), "frames": int(len(frames)),
            "turn_per_orbit_deg": sim["turn_per_orbit_deg"], "exact_deg": sim["exact_deg"],
            "linear_deg": sim["linear_deg"]}


if __name__ == "__main__":
    print(make())
