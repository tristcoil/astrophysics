"""Where does each ray come from? The exact light-ray map of a Schwarzschild black hole.

An observer at radius r looks in some direction. Following that light ray
backward in time, does it come from the sky (infinity), and from which point of
the sky? Spherical symmetry makes this a one-dimensional question: everything
depends only on the angle chi between the viewing direction and the direction
to the black hole, measured in the observer's own frame. The ray then lies in
the plane through the black hole, the observer and that direction, and the
answer is one number, the total angle phi it sweeps around the black hole
between infinity and the observer.

Units: G = c = M = 1 (lengths in units of GM/c^2; the horizon is at r = 2).

**The reference frame.** We work in the frame of a "raindrop": an observer
falling radially from rest at infinity. In Gullstrand-Painleve (GP) coordinates

    ds^2 = -dT^2 + (dr + beta dT)^2 + r^2 dOmega^2,   beta = sqrt(2/r),

the raindrop's proper time is T and its orthonormal frame is simply
e_0 = d_T - beta d_r, e_1 = d_r, e_2, e_3 = (1/r) angular. This frame is
regular at and inside the horizon. Any other radially moving observer (one
hovering at fixed r, or one dropped from rest at finite r) is a Lorentz boost
along e_1 from it (``boost_to_rain``).

A photon seen by the raindrop with local energy 1, arriving from apparent
direction at angle chi from the inward direction (propagating with outward
component n_r = cos chi), has conserved Killing energy and angular momentum

    E = 1 - beta cos chi,    L = r sin chi,    b = L / E   (impact parameter),

its frequency is shifted, relative to light emitted by a source at rest far
away, by g = 1 / E, and it moves outward forward in time iff cos chi > beta.

**The orbit integral.** With u = 1/r, (du/dphi)^2 = P(u) = 2u^3 - u^2 + 1/b^2.
Backward rays that escape sweep

    phi = int_0^u_obs du / sqrt(P)                (moving out the whole way), or
    phi = int_0^u_p du/sqrt(P) + int_u_obs^u_p du/sqrt(P)   (via periapsis u_p),

elliptic integrals that we evaluate exactly with Carlson's symmetric R_F
(DLMF 19.29.4), for the three-real-root case (b > sqrt 27) and the
one-real-root case (b < sqrt 27; complex-conjugate roots enter the same
formula). ``deflection_quad`` and geodesic.py check the result independently.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import integrate
from scipy.special import elliprf

B_CRIT = math.sqrt(27.0)
U_PHOTON = 1.0 / 3.0


# ---------------------------------------------------------------------------
# Frames
# ---------------------------------------------------------------------------

def beta(r):
    """Inward speed of the raindrop (the 'river') relative to hovering observers: sqrt(2M/r)."""
    return np.sqrt(2.0 / np.asarray(r, float))


def boost_to_rain(cos_chi_obs, v_rel):
    """Convert an observer's viewing angle to the raindrop frame.

    ``v_rel`` is the observer's radial velocity relative to the raindrop at the
    same event (positive outward; a hovering observer has v_rel = +beta).
    ``cos_chi_obs`` is the cosine of the angle between the viewing direction and
    the inward radial direction in the observer's frame. Returns the raindrop-frame
    cos chi and the ratio (photon energy seen by the observer) / (seen by the rain).
    """
    c = np.asarray(cos_chi_obs, float)
    v = float(v_rel)
    gamma = 1.0 / math.sqrt(1.0 - v * v)
    # The photon propagates with outward component n_r = cos chi. Boosting from the
    # observer (moving at +v outward relative to the rain) back to the rain frame:
    c_rain = (c + v) / (1.0 + v * c)
    energy_ratio = 1.0 / (gamma * (1.0 + v * c))   # eps_obs / eps_rain
    return c_rain, energy_ratio


def velocity_rel_rain(r, kind: str, r0: float | None = None) -> float:
    """Radial velocity (outward positive) of an observer relative to the raindrop at r.

    kind = 'rain' (from rest at infinity), 'static' (hovering, r > 2 only), or
    'drop' (released from rest at r0 > r).
    """
    b = math.sqrt(2.0 / r)
    if kind == "rain":
        return 0.0
    if kind == "static":
        if r <= 2.0:
            raise ValueError("no hovering observers at or inside the horizon")
        return b
    if kind == "drop":
        # A body released from rest at r0 has energy per unit mass E0 = sqrt(1 - 2/r0).
        # Write its 4-velocity in the rain frame as u = gamma (e_0 + v e_1). With GP
        # components, gamma = -u.e_0 = u^T and w = gamma v = u.e_1 = u_r = beta u^T + u^r,
        # while E0 = -u_T = u^T - beta w and u^r = -sqrt(E0^2 - 1 + beta^2). Eliminating
        # gives the horizon-regular form w = (1 - E0^2) / (beta E0 + sqrt(E0^2 - 1 + beta^2)).
        if r0 is None or r0 <= r:
            raise ValueError("a dropped observer needs r0 > r")
        e0 = math.sqrt(1.0 - 2.0 / r0)
        root = math.sqrt(max(e0 * e0 - 1.0 + b * b, 0.0))
        w = (1.0 - e0 * e0) / (b * e0 + root)
        return w / (e0 + b * w)
    raise ValueError(kind)


# ---------------------------------------------------------------------------
# Conserved quantities of a photon seen in the raindrop frame
# ---------------------------------------------------------------------------

def photon_constants(r, cos_chi):
    """Killing energy E, angular momentum L, impact parameter b and shift g for unit local energy."""
    bt = math.sqrt(2.0 / r)
    c = np.asarray(cos_chi, float)
    s = np.sqrt(np.clip(1.0 - c * c, 0.0, None))
    energy = 1.0 - bt * c
    ang = r * s
    with np.errstate(divide="ignore", invalid="ignore"):
        b = np.where(energy > 0, ang / energy, np.inf)
        g = np.where(energy > 0, 1.0 / energy, np.inf)
    return energy, ang, b, g


def shadow_edge(r: float) -> float:
    """Angular radius chi_sh of the dark region (raindrop frame), radians.

    The edge is where b = sqrt(27): r sin chi + sqrt(27) beta cos chi = sqrt(27).
    Outside the photon sphere (r > 3) it lies on the branch where b grows with chi
    (cos chi > beta); at and inside r = 3 on the branch where b falls.
    """
    bt = math.sqrt(2.0 / r)
    amp = math.hypot(r, B_CRIT * bt)
    delta = math.atan2(B_CRIT * bt, r)
    base = math.asin(min(1.0, B_CRIT / amp))
    small, large = base - delta, math.pi - base - delta
    return small if r > 3.0 else large


# ---------------------------------------------------------------------------
# Roots of P(u) = 2u^3 - u^2 + k,  k = 1/b^2
# ---------------------------------------------------------------------------

def _roots_real(k):
    """Three real roots u1 < u2 <= u3 for 0 < k <= 1/27 (b >= sqrt 27)."""
    th = np.arccos(np.clip(1.0 - 54.0 * k, -1.0, 1.0))
    u3 = 1.0 / 6.0 + np.cos(th / 3.0) / 3.0
    u2 = 1.0 / 6.0 + np.cos(th / 3.0 - 2.0 * np.pi / 3.0) / 3.0
    u1 = 1.0 / 6.0 + np.cos(th / 3.0 - 4.0 * np.pi / 3.0) / 3.0
    return u1, u2, u3


def _roots_complex(k):
    """One real root u1 < 0 and a complex-conjugate pair (u2, u3) for k > 1/27 (b < sqrt 27)."""
    x = np.maximum(54.0 * k - 1.0, 1.0)
    u1 = 1.0 / 6.0 - np.cosh(np.arccosh(x) / 3.0) / 3.0
    p = 0.5 * (0.5 - u1)                      # real part: roots sum to 1/2
    q2 = -k / (2.0 * u1) - p * p              # |u2|^2 = -k/(2 u1)
    q = np.sqrt(np.maximum(q2, 0.0))
    return u1, p + 1j * q, p - 1j * q


def _carlson_segment(x, y, u1, u2, u3):
    """int_y^x du / sqrt(2 (u - u1)(u2 - u)(u3 - u)) for y < x, all factors positive on (y, x).

    DLMF 19.29.4: int_y^x dt / sqrt(prod (a_j + b_j t)) = 2 R_F(U12^2, U13^2, U23^2)
    with X_j = sqrt(a_j + b_j x), Y_j = sqrt(a_j + b_j y) and
    U12 = (X1 X2 Y3 + Y1 Y2 X3)/(x - y), etc. Here (a, b) = (-u1, 1), (u2, -1), (u3, -1);
    u2, u3 may be a complex-conjugate pair, in which case the result is real.
    """
    x1 = np.sqrt(x - u1 + 0j)
    x2 = np.sqrt(u2 - x + 0j)
    x3 = np.sqrt(u3 - x + 0j)
    y1 = np.sqrt(y - u1 + 0j)
    y2 = np.sqrt(u2 - y + 0j)
    y3 = np.sqrt(u3 - y + 0j)
    d = x - y
    u12 = (x1 * x2 * y3 + y1 * y2 * x3) / d
    u13 = (x1 * x3 * y2 + y1 * y3 * x2) / d
    u23 = (x2 * x3 * y1 + y2 * y3 * x1) / d
    val = 2.0 * elliprf(u12 * u12, u13 * u13, u23 * u23) / math.sqrt(2.0)
    return np.real(val)


def deflection(r: float, cos_chi):
    """Total angle phi swept by each backward ray from the raindrop observer at r to infinity.

    Returns (phi, g): phi is NaN for rays that do not come from the sky (they end on
    the horizon in the past: the black hole's shadow); g is the frequency shift
    (observed / emitted by a source at rest at infinity).
    """
    c = np.atleast_1d(np.asarray(cos_chi, float))
    u_obs = 1.0 / r
    bt = math.sqrt(2.0 / r)
    energy, _, b, g = photon_constants(r, c)
    phi = np.full(c.shape, np.nan)

    out_back = bt - c > 0          # backward ray moves outward (forward: inward)
    sky_possible = energy > 0
    small_b = b < B_CRIT

    if r >= 3.0:
        esc_out = out_back & sky_possible
        esc_in = (~out_back) & (b > B_CRIT)
    else:
        esc_out = out_back & sky_possible & small_b
        esc_in = np.zeros_like(out_back)

    # --- moving outward the whole way: int_0^u_obs ---
    for mask, branch in ((esc_out & small_b, "complex"), (esc_out & ~small_b, "real")):
        if not mask.any():
            continue
        bb = np.maximum(b[mask], 1e-12)
        k = 1.0 / (bb * bb)
        if branch == "complex":
            u1, u2, u3 = _roots_complex(k)
        else:
            u1, u2, u3 = _roots_real(k)
            u2 = np.maximum(u2, u_obs)       # guard round-off at the tangent point
        phi[mask] = _carlson_segment(np.full(k.shape, u_obs), np.zeros_like(k), u1, u2, u3)

    # --- via periapsis: int_0^u_p + int_u_obs^u_p ---
    if esc_in.any():
        bb = b[esc_in]
        k = 1.0 / (bb * bb)
        u1, u2, u3 = _roots_real(k)
        up = u2
        y = np.minimum(u_obs, up)
        inner = np.where(up - y > 1e-15, _carlson_segment(up, y, u1, up, u3), 0.0)
        outer = _carlson_segment(up, np.zeros_like(up), u1, up, u3)
        phi[esc_in] = outer + inner

    return phi, g


def deflection_quad(r: float, cos_chi: float) -> float:
    """Slow, independent check of ``deflection`` by adaptive quadrature (scalar)."""
    u_obs = 1.0 / r
    bt = math.sqrt(2.0 / r)
    energy, _, b, _ = photon_constants(r, np.array([cos_chi]))
    energy, b = float(energy[0]), float(b[0])
    if energy <= 0:
        return math.nan
    if b == 0.0:
        return 0.0 if bt - cos_chi > 0 else math.nan    # radial ray from straight overhead
    k = 1.0 / b**2

    def f(u):
        return 1.0 / math.sqrt(max(2 * u**3 - u**2 + k, 1e-300))

    out_back = bt - cos_chi > 0
    if out_back:
        if r < 3.0 and b > B_CRIT:
            return math.nan
        # substitution u = u_obs (1 - t^2) removes nothing singular here; plain quad
        val, _ = integrate.quad(f, 0.0, u_obs, limit=400, epsabs=1e-13, epsrel=1e-12)
        return val
    if r < 3.0 or b < B_CRIT:
        return math.nan
    _, up, _ = (np.asarray(v).item() for v in _roots_real(k))

    # remove the 1/sqrt singularity at the turning point: u = up - s^2
    def g_(s):
        u = up - s * s
        return 2.0 * s * f(u)

    a1, _ = integrate.quad(g_, 0.0, math.sqrt(up), limit=400, epsabs=1e-13, epsrel=1e-12)
    a2, _ = integrate.quad(g_, 0.0, math.sqrt(max(up - u_obs, 0.0)), limit=400, epsabs=1e-13, epsrel=1e-12)
    return a1 + a2


# ---------------------------------------------------------------------------
# Per-frame lookup table
# ---------------------------------------------------------------------------

class RayTable:
    """phi(chi) and g(chi) for one observer position and velocity, tabulated for fast lookup.

    Near the shadow edge phi diverges like -ln(chi - chi_sh) (rays that wind
    around the photon sphere), so the table is built on s = ln(chi - chi_sh),
    where phi is smooth. All angles are in the *observer's* frame.
    """

    def __init__(self, r: float, v_rel: float = 0.0, n: int = 4000, s_min: float = -30.0):
        self.r = r
        self.v_rel = v_rel
        # Edge in the observer frame: boost the raindrop-frame edge.
        cs_rain = math.cos(shadow_edge(r))
        # inverse aberration: observer cos from rain cos (boost by -v)
        cs_obs = (cs_rain - v_rel) / (1.0 - v_rel * cs_rain)
        self.chi_sh = math.acos(max(-1.0, min(1.0, cs_obs)))
        span = math.pi - self.chi_sh
        s_max = math.log(span)
        self.s = np.linspace(s_min, s_max, n)
        chi = self.chi_sh + np.exp(self.s)
        chi[-1] = math.pi
        self.chi = chi
        c_rain, eratio = boost_to_rain(np.cos(chi), v_rel)
        phi, g_rain = deflection(r, c_rain)
        self.phi = phi
        self.g = g_rain * eratio
        bad = ~np.isfinite(phi)
        if bad.any():
            # Only the very first nodes (round-off at the edge) may fail; fill by extrapolation.
            idx = np.flatnonzero(~bad)
            self.phi[bad] = np.interp(self.s[bad], self.s[idx], self.phi[idx])
        # derivative dphi/dchi (for magnification and texture filtering)
        dphids = np.gradient(self.phi, self.s)
        self.dphi_dchi = dphids / np.exp(self.s)

    def lookup(self, chi):
        """phi, g, dphi/dchi at viewing angles chi (radians); NaN inside the shadow."""
        chi = np.asarray(chi, float)
        inside = chi <= self.chi_sh
        s = np.log(np.maximum(chi - self.chi_sh, 1e-300))
        # uniform grid in s: direct index arithmetic instead of a binary search
        t = np.clip((s - self.s[0]) / (self.s[1] - self.s[0]), 0.0, len(self.s) - 1.000001)
        i = t.astype(np.int64)
        f = t - i
        phi = self.phi[i] * (1 - f) + self.phi[i + 1] * f
        g = self.g[i] * (1 - f) + self.g[i + 1] * f
        dphi = self.dphi_dchi[i] * (1 - f) + self.dphi_dchi[i + 1] * f
        phi[inside] = np.nan
        return phi, g, dphi

    def invert(self, phi_target):
        """Viewing angle chi at which the sky angle phi_target is seen (phi decreases with chi)."""
        phi_target = np.asarray(phi_target, float)
        order = np.argsort(self.phi)
        s = np.interp(phi_target, self.phi[order], self.s[order], left=np.nan, right=np.nan)
        return self.chi_sh + np.exp(s)
