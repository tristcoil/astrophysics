"""Orbits from conserved quantities: effective potentials, special radii, precession.

Units: G = c = M = 1 (lengths in GM/c^2, energies per unit rest mass in c^2,
angular momenta per unit mass in GM/c).

**Why symmetries make orbits solvable.** Along a geodesic with four-velocity
u, every Killing vector xi gives a constant g(xi, u). In Schwarzschild the
time translation gives the energy E = (1 - 2/r) dt/dtau and the rotation about
the orbit's axis gives the angular momentum L = r^2 dphi/dtau. Put them into
the normalisation g(u, u) = -1 and the four coupled geodesic equations
collapse into one equation for r alone:

    (dr/dtau)^2 = E^2 - V(r)^2,    V(r)^2 = (1 - 2/r)(1 + L^2/r^2),

a particle in a one-dimensional potential. Every special radius follows from
this curve: circular orbits sit where V has a minimum (stable) or maximum
(unstable); the innermost stable one (ISCO) is where minimum and maximum merge;
light (no rest mass) has V^2 = (1 - 2/r) L^2/r^2, whose peak at r = 3 is the
photon sphere and gives the shadow.

For Kerr the same works in the equatorial plane with E, L from the two
Killing vectors; off the plane the third constant, Carter's Q from the Killing
tensor, is needed (``kerr_hamiltonian``).

Everything is computed twice where possible: once from the potential with a
root finder ("our derivation") and once from the closed formulas of Bardeen,
Press & Teukolsky (1972), which the tests compare.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy import integrate, optimize, special

# ---------------------------------------------------------------------------
# Schwarzschild
# ---------------------------------------------------------------------------


def V2_schwarzschild(r, L):
    """Effective potential squared for a massive particle: (1 - 2/r)(1 + L^2/r^2)."""
    return (1 - 2 / r) * (1 + L**2 / r**2)


def V2_photon_schwarzschild(r):
    """Light: (1 - 2/r)/r^2 per unit b^2 (b = L/E, the impact parameter)."""
    return (1 - 2 / r) / r**2


def circular_schwarzschild(r):
    """E, L of the circular orbit at radius r (exists for r > 3, stable for r > 6)."""
    r = np.asarray(r, float)
    return (r - 2) / np.sqrt(r * (r - 3)), r / np.sqrt(r - 3)


def schwarzschild_special_radii() -> dict:
    """Derive the photon sphere, ISCO and marginally bound orbit from V(r)."""
    # Photon sphere: maximum of (1 - 2/r)/r^2.
    rp = optimize.minimize_scalar(lambda q: -V2_photon_schwarzschild(q), bounds=(2.01, 20),
                                  method="bounded", options={"xatol": 1e-12}).x
    b = 1 / np.sqrt(V2_photon_schwarzschild(rp))
    # ISCO: the L at which dV/dr = d2V/dr2 = 0 coincide; dE/dr = 0 along circular orbits.
    risco = optimize.minimize_scalar(lambda q: circular_schwarzschild(q)[1], bounds=(3.5, 20),
                                     method="bounded", options={"xatol": 1e-12}).x
    E_isco, L_isco = circular_schwarzschild(risco)
    # Marginally bound circular orbit: E = 1.
    rmb = optimize.brentq(lambda q: circular_schwarzschild(q)[0] - 1, 3.0001, 10, xtol=1e-14)
    return {"r_photon": float(rp), "b_crit": float(b), "r_isco": float(risco), "E_isco": float(E_isco),
            "L_isco": float(L_isco), "efficiency": float(1 - E_isco), "r_mb": float(rmb)}


def bound_orbit_EL(p, e):
    """Exact Schwarzschild E, L for an orbit with semi-latus rectum p and eccentricity e
    (r_peri = p/(1+e), r_apo = p/(1-e)); Cutler, Kennefick & Poisson (1994)."""
    p, e = np.asarray(p, float), np.asarray(e, float)
    E = np.sqrt((p - 2 - 2 * e) * (p - 2 + 2 * e) / (p * (p - 3 - e**2)))
    L = p / np.sqrt(p - 3 - e**2)
    return E, L


def precession_exact(p, e):
    """Pericentre advance per radial period (radians), exact in Schwarzschild:
    Delta phi = 4 sqrt(p/(p - 6 + 2e)) K(4e/(p - 6 + 2e)) - 2 pi."""
    p, e = np.asarray(p, float), np.asarray(e, float)
    m = 4 * e / (p - 6 + 2 * e)
    return 4 * np.sqrt(p / (p - 6 + 2 * e)) * special.ellipk(m) - 2 * np.pi


def precession_by_quadrature(E, L):
    """The same advance from E and L alone: integrate dphi/dr between the turning points.

    This is the 'pencil and paper' route: L and E fix the potential, the
    turning points are where E^2 = V^2, and phi advances by
    2 * int L/r^2 / sqrt(E^2 - V^2) dr per radial period.
    """
    # turning points: roots of E^2 - V^2 = 0 -> cubic in u = 1/r: 2 L^2 u^3 - L^2 u^2 + 2u - (1 - E^2) = 0
    roots = np.sort(np.roots([2 * L**2, -L**2, 2, -(1 - E**2)]).real)
    u1, u2, u3 = roots  # u1 < u2 < u3: apocentre 1/u1, pericentre 1/u2
    # With u = u1 + (u2 - u1) sin^2(chi), the integrand is smooth (no endpoint singularity).
    def integrand(chi):
        uu = u1 + (u2 - u1) * np.sin(chi) ** 2
        return 2 / np.sqrt(2 * (u3 - uu))  # from L^2 (u-u1)(u2-u)(u3-u)*2 factorisation
    val, _ = integrate.quad(integrand, 0, np.pi / 2, epsabs=1e-14, epsrel=1e-14, limit=200)
    return 2 * val - 2 * np.pi


# ---------------------------------------------------------------------------
# Kerr, equatorial plane (prograde a > 0, retrograde a < 0 in the formulas below)
# ---------------------------------------------------------------------------


def R_equatorial(r, E, L, a):
    """r^4 (dr/dtau)^2 for an equatorial Kerr geodesic of unit mass."""
    return (E**2 - 1) * r**4 + 2 * r**3 + (a**2 * (E**2 - 1) - L**2) * r**2 + 2 * (L - a * E) ** 2 * r


def V_eff_kerr(r, L, a):
    """Equatorial effective potential: the energy at which dr/dtau = 0 (upper root)."""
    A_ = r**4 + a**2 * r**2 + 2 * a**2 * r
    B_ = 2 * a * L * r
    D_ = r**2 - 2 * r + a**2
    C_ = a**2 * L**2 - D_ * (r**2 + L**2)
    return (B_ + np.sqrt(np.maximum(B_**2 - A_ * C_, 0))) / A_


def dR_equatorial(r, E, L, a):
    return 4 * (E**2 - 1) * r**3 + 6 * r**2 + 2 * (a**2 * (E**2 - 1) - L**2) * r + 2 * (L - a * E) ** 2


def d2R_equatorial(r, E, L, a):
    return 12 * (E**2 - 1) * r**2 + 12 * r + 2 * (a**2 * (E**2 - 1) - L**2)


# Sign convention for the Kerr functions: the hole's spin is |a| and the orbit
# has L > 0; a > 0 means the orbit goes round the same way as the hole spins
# (prograde), a < 0 the opposite way (retrograde). R(r) depends only on
# a^2, L^2 and (L - aE)^2, so flipping the signs of both a and L changes nothing.

def isco_kerr_derived(a: float) -> dict:
    """ISCO from the radial equation alone: solve R = R' = R'' = 0 for (r, E, L).

    Equivalently R/r is a cubic with a triple root at r_isco; matching its
    coefficients gives E^2 = 1 - 2/(3r), L^2 = 2r - 2a^2/(3r), (L - aE)^2 = r^2/3.
    We solve the full system numerically (a root finder, no closed formula), then
    check it against those relations in the tests.
    """
    if a == 0:
        guess = (6.0, np.sqrt(8 / 9), np.sqrt(12))
    else:
        rb = float(isco_bpt(a))  # only a starting point for the root finder
        guess = (rb, np.sqrt(1 - 2 / (3 * rb)), np.sqrt(max(2 * rb - 2 * a**2 / (3 * rb), 1e-6)))
    if abs(a) >= 1 - 1e-12 and a > 0:
        return {"r_isco": 1.0, "E": 1 / np.sqrt(3), "L": 2 / np.sqrt(3), "efficiency": 1 - 1 / np.sqrt(3)}

    def F(v):
        rr, E, L = v
        return [R_equatorial(rr, E, L, a) / rr**3, dR_equatorial(rr, E, L, a) / rr**2,
                d2R_equatorial(rr, E, L, a) / rr]

    sol = optimize.root(F, guess, method="hybr", tol=1e-15)
    rr, E, L = sol.x
    return {"r_isco": float(rr), "E": float(E), "L": float(L), "efficiency": float(1 - E),
            "converged": bool(sol.success)}


def isco_bpt(a):
    """Bardeen, Press & Teukolsky (1972) closed form; a < 0 = retrograde."""
    a = np.asarray(a, float)
    aa, s = np.abs(a), np.where(a >= 0, 1.0, -1.0)
    Z1 = 1 + np.cbrt(1 - aa**2) * (np.cbrt(1 + aa) + np.cbrt(1 - aa))
    Z2 = np.sqrt(3 * aa**2 + Z1**2)
    return 3 + Z2 - s * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))


def circular_kerr(r, a):
    """E, L of prograde (a > 0) or retrograde (a < 0) circular orbits at radius r (BPT 1972)."""
    r, a = np.asarray(r, float), np.asarray(a, float)
    s, aa = np.where(a >= 0, 1.0, -1.0), np.abs(a)
    den = r**0.75 * np.sqrt(r**1.5 - 3 * r**0.5 + 2 * s * aa)
    E = (r**1.5 - 2 * r**0.5 + s * aa) / den
    L = s * (r**2 - 2 * s * aa * r**0.5 + aa**2) / den
    return E, L


def photon_orbit_bpt(a):
    a = np.asarray(a, float)
    return 2 * (1 + np.cos(2 / 3 * np.arccos(-a)))  # a<0 -> retrograde: arccos(|a|)


def photon_orbit_derived(a: float) -> dict:
    """Equatorial circular photon orbit from the null radial equation.

    For light R(r)/E^2 = (r^2 + a^2 - a b)^2 - Delta (b - a)^2 with b = L/E; a
    circular orbit needs R = dR/dr = 0: two equations for (r, b), solved
    numerically. (Eliminating b by hand gives r (r - 3)^2 = 4 a^2.)
    """
    if abs(a) >= 1 - 1e-12 and a > 0:
        return {"r": 1.0, "b": 2.0}
    rg = float(photon_orbit_bpt(a)) if a != 0 else 3.0

    def F(v):
        rr, b = v
        D_ = rr**2 - 2 * rr + a**2
        return [(rr**2 + a**2 - a * b) ** 2 - D_ * (b - a) ** 2,
                4 * rr * (rr**2 + a**2 - a * b) - (2 * rr - 2) * (b - a) ** 2]

    bg = -(rg**3 - 3 * rg**2 + a**2 * rg + a**2) / (a * (rg - 1)) if a != 0 else np.sqrt(27)
    sol = optimize.root(F, (rg, bg), method="hybr", tol=1e-15)
    return {"r": float(sol.x[0]), "b": float(sol.x[1])}


def marginally_bound_bpt(a):
    a = np.asarray(a, float)
    return 2 - a + 2 * np.sqrt(1 - a)


def marginally_bound_derived(a: float) -> float:
    """Circular orbit with E = 1: solve R = R' = 0 for (r, L) with E fixed at 1."""
    rg = float(marginally_bound_bpt(a))

    def F(v):
        rr, L = v
        return [R_equatorial(rr, 1.0, L, a) / rr**3, dR_equatorial(rr, 1.0, L, a) / rr**2]

    sol = optimize.root(F, (rg, 2 * np.sqrt(rg) * 0.99 + 0.02), method="hybr", tol=1e-15)
    return float(sol.x[0])


def circular_from_R(r: float, a: float) -> tuple:
    """E, L of the circular orbit at r: solve R(r) = R'(r) = 0 for (E, L)."""
    E0, L0 = circular_kerr(r, a)

    def F(v):
        E, L = v
        return [R_equatorial(r, E, L, a) / r**3, dR_equatorial(r, E, L, a) / r**2]

    sol = optimize.root(F, (float(E0) * 1.001, float(L0) * 0.999), method="hybr", tol=1e-15)
    return float(sol.x[0]), float(sol.x[1])


def horizon(a):
    a = np.asarray(a, float)
    return 1 + np.sqrt(1 - a**2)


def special_radii_table(spins) -> dict:
    """All radii and efficiencies against spin, derived and closed-form."""
    out = {k: [] for k in ("a", "r_plus", "r_isco_pro", "r_isco_retro", "r_ph_pro", "r_ph_retro",
                           "r_mb_pro", "r_mb_retro", "eff_pro", "eff_retro", "max_dev_from_bpt")}
    for s in spins:
        pro = isco_kerr_derived(s)
        retro = isco_kerr_derived(-s) if s > 0 else pro
        ph_pro, ph_retro = photon_orbit_derived(s), photon_orbit_derived(-s)
        mb_pro = marginally_bound_derived(s) if s < 1 else 1.0
        mb_retro = marginally_bound_derived(-s)
        dev = max(abs(pro["r_isco"] - float(isco_bpt(s))), abs(retro["r_isco"] - float(isco_bpt(-s))),
                  abs(ph_pro["r"] - float(photon_orbit_bpt(s))), abs(ph_retro["r"] - float(photon_orbit_bpt(-s))),
                  abs(mb_pro - float(marginally_bound_bpt(s))), abs(mb_retro - float(marginally_bound_bpt(-s))))
        out["a"].append(float(s))
        out["r_plus"].append(float(horizon(s)))
        out["r_isco_pro"].append(pro["r_isco"])
        out["r_isco_retro"].append(retro["r_isco"])
        out["r_ph_pro"].append(ph_pro["r"])
        out["r_ph_retro"].append(ph_retro["r"])
        out["r_mb_pro"].append(mb_pro)
        out["r_mb_retro"].append(mb_retro)
        out["eff_pro"].append(pro["efficiency"])
        out["eff_retro"].append(retro["efficiency"])
        out["max_dev_from_bpt"].append(float(dev))
    return out


def efficiency_from_spin(a):
    """Accretion efficiency 1 - E_isco for prograde discs, vectorised (BPT ISCO)."""
    r = isco_bpt(a)
    return 1 - np.sqrt(1 - 2 / (3 * r))


# ---------------------------------------------------------------------------
# Black hole horizon quantities (Kerr, G = c = M = 1)
# ---------------------------------------------------------------------------

def horizon_quantities(chi):
    """r+, surface gravity, horizon angular velocity, area, irreducible mass for spin chi."""
    chi = np.asarray(chi, float)
    rp = 1 + np.sqrt(1 - chi**2)
    rm = 1 - np.sqrt(1 - chi**2)
    kappa = (rp - rm) / (2 * (rp**2 + chi**2))
    omega = chi / (rp**2 + chi**2)
    area = 4 * np.pi * (rp**2 + chi**2)
    m_irr = np.sqrt(area / (16 * np.pi))
    return {"r_plus": rp, "kappa": kappa, "Omega_H": omega, "area": area, "M_irr": m_irr}


def ergosphere_radius(theta, a):
    """Static limit r_E(theta) = 1 + sqrt(1 - a^2 cos^2 theta), where g_tt = 0."""
    return 1 + np.sqrt(1 - a**2 * np.cos(theta) ** 2)


# ---------------------------------------------------------------------------
# Numerical geodesics from the Christoffel symbols (for Noether checks)
# ---------------------------------------------------------------------------

def christoffel_function(metric, subs):
    """Return f(x) -> Gamma[a, b, c] as a numpy array, lambdified from the symbolic metric."""
    G_ = metric.christoffel
    n = metric.dim
    exprs = [sp.cancel(sp.together(sp.sympify(G_[i][j][k]).subs(subs)))
             for i in range(n) for j in range(n) for k in range(n)]
    fn = sp.lambdify(metric.coords, exprs, "math", cse=True)
    return lambda X: np.array(fn(*X), dtype=float).reshape(n, n, n)


def metric_function(metric, subs):
    fn = sp.lambdify(metric.coords, metric.g.subs(subs), "numpy")
    return lambda X: np.array(fn(*X), dtype=float)


def integrate_geodesic(metric, subs, x0, u0, tau_max, n_out=4000, rtol=1e-12, atol=1e-12):
    """Integrate the geodesic equation d2x/dtau2 = -Gamma u u from (x0, u0)."""
    Gam = christoffel_function(metric, subs)

    def rhs(_, s):
        X, U = s[:4], s[4:]
        G_ = Gam(X)
        return np.concatenate([U, -np.einsum("abc,b,c->a", G_, U, U)])

    taus = np.linspace(0, tau_max, n_out)
    sol = integrate.solve_ivp(rhs, (0, tau_max), np.concatenate([x0, u0]), method="DOP853",
                              t_eval=taus, rtol=rtol, atol=atol)
    return sol.t, sol.y[:4].T, sol.y[4:].T


def conserved_along(metric, subs, X, U, fields: dict, tensors: dict | None = None) -> dict:
    """g(xi, u) for each Killing vector and K(u, u) for each tensor, along an orbit."""
    gfn = metric_function(metric, subs)
    out = {}
    vec_fns = {n_: sp.lambdify(metric.coords, [sp.sympify(c_).subs(subs) for c_ in f_], "numpy")
               for n_, f_ in fields.items()}
    ten_fns = {n_: sp.lambdify(metric.coords, T.subs(subs), "numpy") for n_, T in (tensors or {}).items()}
    for n_, fn in vec_fns.items():
        vals = []
        for Xi, Ui in zip(X, U):
            xi = np.array(fn(*Xi), dtype=float) * np.ones(4)
            vals.append(xi @ gfn(Xi) @ Ui)
        out[n_] = np.array(vals)
    for n_, fn in ten_fns.items():
        out[n_] = np.array([Ui @ np.array(fn(*Xi), dtype=float) @ Ui for Xi, Ui in zip(X, U)])
    out["norm"] = np.array([Ui @ gfn(Xi) @ Ui for Xi, Ui in zip(X, U)])
    return out


# ---------------------------------------------------------------------------
# Hamiltonian geodesics in stationary axisymmetric spacetimes (Poincare sections)
# ---------------------------------------------------------------------------

class AxisymmetricHamiltonian:
    """H = 1/2 g^ab p_a p_b with p_t = -E and p_phi = L conserved (the two Killing
    vectors). What is left is motion in the meridional plane (q1, q2) = (r, cos
    theta) or (X, Y): two degrees of freedom with one constant, H = -1/2. A third
    constant (Carter's, from a Killing tensor) makes the motion regular; without
    one it may be chaotic. Coordinates: (t, q1, q2, phi).
    """

    def __init__(self, metric, subs=None):
        subs = subs or {}
        gi = metric.ginv.subs(subs)
        q1, q2 = metric.coords[1], metric.coords[2]
        self.coords = metric.coords
        E, L, p1, p2 = sp.symbols("E L p1 p2")
        p = [-E, p1, p2, L]
        Ham = sp.Rational(1, 2) * sum(gi[i, j] * p[i] * p[j] for i in range(4) for j in range(4))
        args = (q1, q2, p1, p2, E, L)
        self.H = sp.lambdify(args, Ham, "math", cse=True)
        flows = [sp.diff(Ham, p1), sp.diff(Ham, p2), -sp.diff(Ham, q1), -sp.diff(Ham, q2)]
        flows = [sp.cancel(sp.together(f_)) for f_ in flows]
        self.flow = sp.lambdify(args, flows, "math", cse=True)
        # dt/dlambda and dphi/dlambda (for drawing the orbit in 3D)
        self.tdot = sp.lambdify(args, -sp.diff(Ham, E), "math", cse=True)
        self.phidot = sp.lambdify(args, sp.diff(Ham, L), "math", cse=True)
        # H = A(q) p2^2 + (terms without p2), since the meridional block of g^ab is diagonal
        self._A = sp.lambdify((q1, q2), sp.Rational(1, 2) * gi[2, 2], "math")
        self._rest = sp.lambdify(args, Ham.subs(p2, 0), "math", cse=True)
        assert sp.simplify(gi[1, 2]) == 0, "meridional block must be diagonal"
        self._gtt = sp.lambdify((q1, q2), gi[0, 0], "math")
        self._gtp = sp.lambdify((q1, q2), gi[0, 3], "math")
        self._gpp = sp.lambdify((q1, q2), gi[3, 3], "math")
        self.metric = metric
        self.subs = subs

    def energy_at_rest(self, q1v, q2v, L):
        """Effective potential: the energy of a body momentarily at rest in the meridional
        plane (p1 = p2 = 0) with angular momentum L, i.e. the future-directed root of
        g^tt E^2 - 2 g^tphi E L + g^phiphi L^2 = -1 (valid also inside an ergosphere)."""
        gtt, gtp, gpp = self._gtt(q1v, q2v), self._gtp(q1v, q2v), self._gpp(q1v, q2v)
        disc = (gtp * L) ** 2 - gtt * (gpp * L**2 + 1)
        if disc < 0:
            return float("nan")
        return (gtp * L - np.sqrt(disc)) / gtt

    def p2_for(self, q1v, q2v, p1v, E, L, sign=1.0):
        """The p2 that puts the orbit on the mass shell H = -1/2 (None if forbidden)."""
        val = (-0.5 - self._rest(q1v, q2v, p1v, 0.0, E, L)) / self._A(q1v, q2v)
        return None if val < 0 else sign * np.sqrt(val)

    def integrate(self, y0, E, L, lam_max, section_q2=0.0, rtol=1e-11, atol=1e-12, dense=False):
        def rhs(_, s):
            return np.array(self.flow(s[0], s[1], s[2], s[3], E, L), dtype=float)

        def cross(_, s):
            return s[1] - section_q2
        cross.direction = 1.0  # upward crossings of the equatorial plane only
        sol = integrate.solve_ivp(rhs, (0, lam_max), y0, method="DOP853", rtol=rtol, atol=atol,
                                  events=cross, dense_output=dense)
        return sol


# ---------------------------------------------------------------------------
# Generic Kerr orbits from their turning points (also what the explorer page does)
# ---------------------------------------------------------------------------

def R_kerr(r, a, E, L, Q):
    """Radial function of a unit-mass Kerr geodesic: Sigma^2 (dr/dtau)^2."""
    D_ = r**2 - 2 * r + a**2
    P = E * (r**2 + a**2) - a * L
    return P**2 - D_ * (r**2 + (L - a * E) ** 2 + Q)


def kerr_constants(a, rp, ra, inc_deg):
    """(E, L, Q) of the bound Kerr orbit with turning points rp, ra and inclination
    angle iota defined by cos(iota) = L / sqrt(L^2 + Q) (the usual definition for Kerr
    orbits; exactly the tilt of the orbital plane when a = 0). iota > 90 deg is
    retrograde. Unknowns (E, Lt) with L = Lt cos(iota), Q = Lt^2 sin^2(iota), solved
    with a root finder from the Schwarzschild values. Returns ok = False if the
    turning points do not bound a genuine orbit."""
    ci, si = np.cos(np.radians(inc_deg)), np.sin(np.radians(inc_deg))
    p, e = 2 * rp * ra / (rp + ra), (ra - rp) / (ra + rp)
    if p - 3 - e**2 > 0.5:
        E0, Lt0 = bound_orbit_EL(p, e)
    else:  # close to a spinning hole: start from the circular orbit in between
        E0, L0 = circular_kerr(0.5 * (rp + ra), a if inc_deg <= 90 else -a)
        if not (np.isfinite(E0) and np.isfinite(L0)):
            return {"ok": False}
        Lt0 = abs(L0)

    def F(v):
        E, Lt = v
        L, Q = Lt * ci, (Lt * si) ** 2
        return [R_kerr(rp, a, E, L, Q) / rp**4, R_kerr(ra, a, E, L, Q) / ra**4]

    sol = optimize.root(F, (float(E0), float(Lt0)), method="hybr", tol=1e-14)
    E, Lt = sol.x
    L, Q = Lt * ci, (Lt * si) ** 2
    rs = np.linspace(rp, ra, 41)[1:-1]
    resid = float(np.max(np.abs(F(sol.x))))
    ok = bool(resid < 1e-10 and 0 < E < 1 and Lt > 0 and np.all(R_kerr(rs, a, E, L, Q) > 0))
    return {"E": float(E), "L": float(L), "Q": float(Q), "ok": ok, "residual": resid}


def kerr_initial_state(a, rp, ra, inc_deg):
    """Position and four-velocity (t, r, x = cos theta, phi) at apocentre on the equator."""
    c = kerr_constants(a, rp, ra, inc_deg)
    E, L, Q = c["E"], c["L"], c["Q"]
    r = ra
    Sig = r**2
    D_ = r**2 - 2 * r + a**2
    P = E * (r**2 + a**2) - a * L
    ut = ((r**2 + a**2) / D_ * P - a * (a * E - L)) / Sig
    uph = (a * P / D_ - a * E + L) / Sig
    uth = np.sqrt(max(Q, 0.0)) / Sig          # moving towards the south (theta increasing)
    ux = -uth                                 # dx/dtau = -sin(theta) dtheta/dtau, sin = 1 at the equator
    return np.array([0.0, r, 0.0, 0.0]), np.array([ut, 0.0, ux, uph]), c
