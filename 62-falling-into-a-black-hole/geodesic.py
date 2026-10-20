"""Brute-force light rays: an independent numerical check of raymap.py.

Instead of conserved quantities and elliptic integrals, integrate the geodesic
equation  d p^a / d lambda = -Gam^a_bc p^b p^c  directly, with the Christoffel
symbols that SymPy derived from the GP metric (symbolic.py), in the equatorial
plane, coordinates (T, r, phi), M = 1.

The ray starts at the observer with the momentum the observer measures:
for a raindrop at radius r looking at angle chi from the inward direction, the
photon propagates with p = e_0 + cos(chi) e_1 - sin(chi) e_phi, where
e_0 = d_T - beta d_r, e_1 = d_r and e_phi = (1/r) d_phi (its frame, built from
the metric, not from the formulas raymap.py uses). Reversing p follows the ray
back in time. An adaptive Runge-Kutta (DOP853) step, controlled on the relative
error, takes care of rays that skim the photon sphere.

The integration stops at r = r_far (escape: the remaining angle out to infinity
is added from a separate quadrature of the weak-field tail), or when the ray
has crept back to the past horizon (captured: the ray belongs to the shadow).
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate

from symbolic import christoffel_functions


def _rhs_factory():
    gam_flat = christoffel_functions()

    def rhs(lam, y):
        # y = (T, r, phi, pT, pr, pphi)
        r = y[1]
        G = gam_flat(r)
        p = y[3:]
        dp = [0.0, 0.0, 0.0]
        for a in range(3):
            s = 0.0
            base = 9 * a
            for b in range(3):
                for c in range(3):
                    s += G[base + 3 * b + c] * p[b] * p[c]
            dp[a] = -s
        return [p[0], p[1], p[2], dp[0], dp[1], dp[2]]

    return rhs


_RHS = None


def trace(r_obs: float, chi: float, r_far: float = 2000.0, max_lambda: float = 1e6):
    """Follow one ray back from a raindrop at r_obs looking at angle chi.

    Returns (escaped, phi_swept, n_steps). phi_swept is the total angle around the
    black hole between the observer and infinity (tail included) if it escaped.
    """
    global _RHS
    if _RHS is None:
        _RHS = _rhs_factory()
    beta = math.sqrt(2.0 / r_obs)
    # forward momentum from the observer's frame, then reversed to go back in time
    pT, pr, pphi = 1.0, -beta + math.cos(chi), -math.sin(chi) / r_obs
    y0 = [0.0, r_obs, 0.0, -pT, -pr, -pphi]

    def escape(lam, y):
        return y[1] - r_far
    escape.terminal = True
    escape.direction = 1

    def stuck(lam, y):
        # backward ray hugging the horizon from outside or crossing r < 0.05
        return y[1] - 0.05
    stuck.terminal = True

    def horizon_past(lam, y):
        # going back in time a ray can only approach the (past) horizon asymptotically;
        # T running to -1e4 near r ~ 2 means it never came from the sky
        return y[0] + 1e4
    horizon_past.terminal = True

    sol = integrate.solve_ivp(_RHS, (0.0, max_lambda), y0, method="DOP853", rtol=1e-11, atol=1e-12,
                              events=(escape, stuck, horizon_past), dense_output=False)
    if sol.status == 1 and sol.t_events[0].size:
        y_end = sol.y_events[0][0]
        phi_end = abs(y_end[2])
        # tail from r_far to infinity: the backward ray is outgoing there
        L = r_obs * math.sin(chi)
        E = 1.0 - beta * math.cos(chi)
        k = (E / L) ** 2 if L > 0 else math.inf
        if not math.isfinite(k):
            return True, phi_end, sol.t.size
        u_end = 1.0 / y_end[1]
        tail, _ = integrate.quad(lambda u: 1.0 / math.sqrt(2 * u**3 - u**2 + k), 0.0, u_end,
                                 epsabs=1e-14, epsrel=1e-12)
        return True, phi_end + tail, sol.t.size
    return False, math.nan, sol.t.size
