"""Black-hole shadows from conserved quantities, against the EHT rings measured in I5.

A light ray near a black hole has a conserved energy E and angular momentum L
(the two Killing vectors) and, in Kerr, Carter's constant Q (the Killing
tensor). Only their ratios matter for light: the impact parameter b = L/E and
eta = Q/E^2. Rays that come in with exactly the right (b, eta) wind forever
around unstable spherical photon orbits; slightly more and they escape,
slightly less and they fall in. Seen from far away, those critical rays draw a
closed curve on the sky -- the edge of the shadow, the "critical curve".

* Schwarzschild: the potential for light, (1 - 2M/r) L^2 / r^2, peaks at the
  photon sphere r = 3M, so the critical impact parameter is b = sqrt(27) M and
  the shadow is a circle 2 sqrt(27) = 10.39 GM/c^2 across.
* Kerr (Bardeen 1973): spherical photon orbits exist between the prograde and
  retrograde equatorial ones. An orbit at radius r has

      xi(r)  = L/E   = -[r^3 - 3 M r^2 + a^2 r + a^2 M] / [a (r - M)],
      eta(r) = Q/E^2 = r^3 [4 a^2 M - r (r - 3M)^2] / [a^2 (r - M)^2],

  and an observer at inclination i sees it at
  alpha = -xi / sin i,  beta = +-sqrt(eta + a^2 cos^2 i - xi^2 cot^2 i).

Spin flattens one side of the shadow and shifts it, but changes its mean
diameter by at most a few per cent, which is why a shadow measures M/D.
The EHT's bright ring is not the critical curve itself: it is emission from
gas just outside it, so the EHT calibrates ring = alpha x theta_g with
simulations (alpha ~ 10.8-11.6 instead of 10.39); I5 discussed this.
"""

from __future__ import annotations

import json

import numpy as np
from scipy import optimize

import config as C
import h1link
import literature as lit


def critical_curve(a: float, inc_deg: float, n: int = 4000) -> tuple:
    """Bardeen's critical curve (alpha, beta) in units of GM/c^2, closed and ordered."""
    inc = np.radians(inc_deg)
    if a < 1e-6:
        th = np.linspace(0, 2 * np.pi, n)
        rb = np.sqrt(27.0)
        return rb * np.cos(th), rb * np.sin(th)
    r_pro = 2 * (1 + np.cos(2 / 3 * np.arccos(-a)))
    r_ret = 2 * (1 + np.cos(2 / 3 * np.arccos(a)))

    def parts(r):
        xi = -(r**3 - 3 * r**2 + a**2 * r + a**2) / (a * (r - 1))
        eta = r**3 * (4 * a**2 - r * (r - 3) ** 2) / (a**2 * (r - 1) ** 2)
        return xi, eta + a**2 * np.cos(inc) ** 2 - xi**2 / np.tan(inc) ** 2

    # The visible orbits are those with beta^2 >= 0; find where beta^2 = 0 exactly,
    # so that the curve closes on the alpha axis without a seam.
    grid = np.linspace(r_pro, r_ret, 20001)
    b2g = parts(grid)[1]
    lo, hi = r_pro, r_ret
    if b2g[0] < 0:
        i = int(np.argmax(b2g >= 0))
        lo = optimize.brentq(lambda q: parts(q)[1], grid[i - 1], grid[i], xtol=1e-14)
    if b2g[-1] < 0:
        i = len(grid) - 1 - int(np.argmax(b2g[::-1] >= 0))
        hi = optimize.brentq(lambda q: parts(q)[1], grid[i], grid[i + 1], xtol=1e-14)
    u = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n)))   # points clustered at both ends
    r = lo + (hi - lo) * u
    xi, b2 = parts(r)
    alpha = -xi / np.sin(inc)
    beta = np.sqrt(np.maximum(b2, 0.0))
    al = np.concatenate([alpha, alpha[::-1]])
    be = np.concatenate([beta, -beta[::-1]])
    return al, be


def mean_diameter(alpha, beta) -> dict:
    """Average diameter about the geometric centre (Johannsen & Psaltis 2010) and
    the diameter of the circle with the same area; also the centre's offset."""
    ac = 0.5 * (alpha.max() + alpha.min())
    x, y = alpha - ac, beta
    ang = np.arctan2(y, x)
    rad = np.hypot(x, y)
    o = np.argsort(ang)
    grid = np.linspace(-np.pi, np.pi, 2001)
    R = np.interp(grid, ang[o], rad[o], period=2 * np.pi)
    area = 0.5 * np.abs(np.sum(alpha * np.roll(beta, -1) - np.roll(alpha, -1) * beta))
    return {"mean_diameter": float(2 * np.mean(R)), "area_diameter": float(2 * np.sqrt(area / np.pi)),
            "offset": float(ac), "width": float(alpha.max() - alpha.min()),
            "height": float(beta.max() - beta.min())}


def diameter_grid(spins=None, incs=None) -> dict:
    spins = np.array([0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.94, 0.98, 0.998]) if spins is None else spins
    incs = np.array([1, 17, 30, 45, 60, 75, 90]) if incs is None else incs
    D = np.zeros((len(spins), len(incs)))
    for i, s in enumerate(spins):
        for j, ic in enumerate(incs):
            D[i, j] = mean_diameter(*critical_curve(float(s), float(ic)))["mean_diameter"]
    return {"spins": spins.tolist(), "incs": incs.tolist(), "mean_diameter_M": D.tolist(),
            "min": float(D.min()), "max": float(D.max())}


def theta_g_sgra_samples(thin: int = 4) -> np.ndarray:
    """theta_g = GM/(c^2 D) for Sgr A* in micro-arcseconds, over H1's posterior."""
    op = h1link.modules().op
    ch = h1link.chain()
    return op.gravitational_angle_uas(ch["gm6"][::thin], ch["r0"][::thin])


def analyse() -> dict:
    grid = diameter_grid()
    i5 = json.loads(C.I5_RESULTS.read_text()) if C.I5_RESULTS.exists() else None
    if i5:
        sg_ring = (i5["sgra_test"]["amplitudes"]["d_peak"]["value"], i5["sgra_test"]["amplitudes"]["d_peak"]["err"])
        m87 = i5["m87_2017_ring"]["d"]
        m87_ring = (m87["mean"], m87["error_of_mean"], m87["rms_scatter"])
        pub = i5["published"]
        m87_theta = {"stars, Gebhardt+11": tuple(pub["theta_dyn_stars_uas"]),
                     "stars, Liepold+23": (i5["m87_mass"]["predictions"]["stars (Liepold+23)"]["theta_g"],
                                           i5["m87_mass"]["predictions"]["stars (Liepold+23)"]["theta_up"],
                                           i5["m87_mass"]["predictions"]["stars (Liepold+23)"]["theta_lo"]),
                     "stars, Simon+24": (i5["m87_mass"]["predictions"]["stars (Simon+24)"]["theta_g"],
                                         i5["m87_mass"]["predictions"]["stars (Simon+24)"]["theta_up"],
                                         i5["m87_mass"]["predictions"]["stars (Simon+24)"]["theta_lo"]),
                     "gas, Walsh+13": tuple(pub["theta_dyn_gas_uas"])}
        src = "36-m87-shadow/data/results.json"
    else:
        sg_ring = C.FALLBACK["sgra_ring_uas"]
        m87_ring = C.FALLBACK["m87_ring_uas"] + (1.1,)
        m87_theta = {"stars, Gebhardt+11": (3.62, 0.60, 0.34)}
        src = "fallback values (I5 results missing)"

    tg = theta_g_sgra_samples()
    sg_inc = lit.SHADOW["sgra_inclination_deg_max"]
    sg = {"theta_g_uas": [float(np.median(tg)), float(np.std(tg))],
          "schwarzschild_shadow_uas": [float(np.median(2 * np.sqrt(27) * tg)), float(np.std(2 * np.sqrt(27) * tg))],
          "kerr_range_uas_all_incl": [grid["min"] * float(np.median(tg)), grid["max"] * float(np.median(tg))],
          "ring_uas": list(sg_ring),
          "ring_over_theta_g": [sg_ring[0] / float(np.median(tg)),
                                sg_ring[0] / float(np.median(tg)) * np.hypot(sg_ring[1] / sg_ring[0],
                                                                             np.std(tg) / np.median(tg))],
          "inclination_assumed_deg": sg_inc}
    sg["ring_over_shadow"] = [sg["ring_over_theta_g"][0] / (2 * np.sqrt(27)), sg["ring_over_theta_g"][1] / (2 * np.sqrt(27))]
    sg["spin_spread_at_incl"] = {f"a={s}": mean_diameter(*critical_curve(s, sg_inc))["mean_diameter"] * sg["theta_g_uas"][0]
                                 for s in (0.0, 0.5, 0.94, 0.998)}
    m87_inc = lit.SHADOW["m87_inclination_deg"]
    m = {"ring_uas": list(m87_ring), "inclination_assumed_deg": m87_inc, "by_mass": {}}
    for k, (th, up, lo) in m87_theta.items():
        m["by_mass"][k] = {"theta_g_uas": [th, up, lo], "schwarzschild_shadow_uas": 2 * np.sqrt(27) * th,
                           "ring_over_theta_g": m87_ring[0] / th,
                           "ring_over_theta_g_up": m87_ring[0] / (th - lo) - m87_ring[0] / th,
                           "ring_over_theta_g_lo": m87_ring[0] / th - m87_ring[0] / (th + up)}
    m["spin_spread_at_incl_M"] = {f"a={s}": mean_diameter(*critical_curve(s, m87_inc))["mean_diameter"]
                                  for s in (0.0, 0.5, 0.94, 0.998)}
    return {"grid": grid, "sgra": sg, "m87": m, "source": src,
            "published_range_M": lit.SHADOW["kerr_diameter_range_M"]}
