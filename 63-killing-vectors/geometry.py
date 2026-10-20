"""Spacetime geometry with computer algebra (SymPy): metrics, Christoffel symbols,
curvature, and Killing's equation.

A **metric** g_ab(x) gives the interval between nearby events,
ds^2 = g_ab dx^a dx^b. Everything else follows from it:

* the **Christoffel symbols** Gamma^a_bc = 1/2 g^ad (d_b g_dc + d_c g_db - d_d g_bc)
  say how coordinate directions twist from point to point; free fall obeys
  d^2x^a/dtau^2 = -Gamma^a_bc (dx^b/dtau)(dx^c/dtau) (the geodesic equation);
* the **Riemann tensor** R^a_bcd measures real curvature (tidal forces), which
  no choice of coordinates can remove;
* a vector field xi is a **Killing vector** if moving every point a little
  along it leaves the metric unchanged. In equations, its Lie derivative
  vanishes:

      (L_xi g)_ab = xi^c d_c g_ab + g_cb d_a xi^c + g_ac d_b xi^c = 0,

  which is the same thing as Killing's equation nabla_a xi_b + nabla_b xi_a = 0.
  Both forms are coded here, and the tests check that they agree.

Units are geometric (G = c = 1): masses are lengths, M stands for GM/c^2.

Coordinates. Spherical metrics are written with x = cos(theta) instead of theta
wherever possible. That turns every metric component into a ratio of
polynomials (sin^2 theta = 1 - x^2), so exact Taylor expansions about a point
with rational coordinates have rational coefficients and the Killing-vector
count in ``jets.py`` can be done in exact arithmetic.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from functools import cached_property

import sympy as sp

t, r, x, y, z, phi, theta, chi, u, v = sp.symbols("t r x y z phi theta chi u v", real=True)
M, a, H, k = sp.symbols("M a H k", positive=True)


def simplify_rational(expr):
    """Cheap, reliable simplification for rational functions: cancel common factors."""
    return sp.cancel(sp.together(expr))


@dataclass
class Metric:
    """A metric tensor in a coordinate chart.

    ``g`` is a symmetric sympy Matrix whose entries are functions of ``coords``
    and of parameters (M, a, ...). ``simplifier`` is applied to every derived
    component; the default suits metrics whose components are rational
    functions of the coordinates.
    """

    name: str
    coords: tuple
    g: sp.Matrix
    description: str = ""
    simplifier: object = field(default=simplify_rational, repr=False)

    @property
    def dim(self) -> int:
        return len(self.coords)

    @cached_property
    def ginv(self) -> sp.Matrix:
        return self.g.inv(method="LU").applyfunc(self.simplifier)

    @cached_property
    def dg(self) -> list:
        """dg[c][a][b] = d_c g_ab."""
        return [[[sp.diff(self.g[a_, b_], c_) for b_ in range(self.dim)] for a_ in range(self.dim)]
                for c_ in self.coords]

    @cached_property
    def christoffel(self) -> list:
        """Gamma[a][b][c] = Gamma^a_bc (second kind)."""
        n = self.dim
        dg, gi = self.dg, self.ginv
        first = [[[sp.Rational(1, 2) * (dg[b_][d_][c_] + dg[c_][d_][b_] - dg[d_][b_][c_])
                   for c_ in range(n)] for b_ in range(n)] for d_ in range(n)]
        out = [[[0] * n for _ in range(n)] for _ in range(n)]
        for a_ in range(n):
            for b_ in range(n):
                for c_ in range(b_, n):
                    val = self.simplifier(sum(gi[a_, d_] * first[d_][b_][c_] for d_ in range(n)))
                    out[a_][b_][c_] = out[a_][c_][b_] = val
        return out

    @cached_property
    def riemann(self) -> list:
        """R[a][b][c][d] = R^a_bcd = d_c Gamma^a_db - d_d Gamma^a_cb + Gamma^a_ce Gamma^e_db - ..."""
        n, G_ = self.dim, self.christoffel
        X = self.coords
        R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for a_, b_, c_, d_ in itertools.product(range(n), repeat=4):
            if c_ >= d_:
                continue
            val = sp.diff(G_[a_][d_][b_], X[c_]) - sp.diff(G_[a_][c_][b_], X[d_])
            val += sum(G_[a_][c_][e] * G_[e][d_][b_] - G_[a_][d_][e] * G_[e][c_][b_] for e in range(n))
            val = self.simplifier(val)
            R[a_][b_][c_][d_] = val
            R[a_][b_][d_][c_] = -val
        return R

    @cached_property
    def ricci(self) -> sp.Matrix:
        n, R = self.dim, self.riemann
        return sp.Matrix(n, n, lambda b_, d_: self.simplifier(sum(R[a_][b_][a_][d_] for a_ in range(n))))

    @cached_property
    def ricci_scalar(self):
        n = self.dim
        return self.simplifier(sum(self.ginv[a_, b_] * self.ricci[a_, b_] for a_ in range(n) for b_ in range(n)))

    def kretschmann(self):
        """R_abcd R^abcd: a curvature invariant (48 M^2 / r^6 for Schwarzschild)."""
        n, R, g, gi = self.dim, self.riemann, self.g, self.ginv
        low = [[[[self.simplifier(sum(g[a_, e] * R[e][b_][c_][d_] for e in range(n)))
                  for d_ in range(n)] for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
        up = [[[[self.simplifier(sum(gi[b_, f] * gi[c_, h] * gi[d_, l] * R[a_][f][h][l]
                                     for f in range(n) for h in range(n) for l in range(n)))
                 for d_ in range(n)] for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
        return self.simplifier(sum(low[a_][b_][c_][d_] * up[a_][b_][c_][d_]
                                   for a_, b_, c_, d_ in itertools.product(range(n), repeat=4)))

    # -- symmetries --------------------------------------------------------

    def lie_derivative(self, xi) -> sp.Matrix:
        """(L_xi g)_ab for a vector field with contravariant components ``xi``."""
        n, X, g = self.dim, self.coords, self.g
        xi = [sp.sympify(c) for c in xi]
        out = sp.zeros(n, n)
        for a_ in range(n):
            for b_ in range(a_, n):
                val = sum(xi[c_] * self.dg[c_][a_][b_] for c_ in range(n))
                val += sum(g[c_, b_] * sp.diff(xi[c_], X[a_]) + g[a_, c_] * sp.diff(xi[c_], X[b_])
                           for c_ in range(n))
                out[a_, b_] = out[b_, a_] = val
        return out

    def killing_operator(self, xi) -> sp.Matrix:
        """nabla_a xi_b + nabla_b xi_a, computed with covariant derivatives."""
        n, X, G_ = self.dim, self.coords, self.christoffel
        low = [sum(self.g[b_, c_] * sp.sympify(xi[c_]) for c_ in range(n)) for b_ in range(n)]
        out = sp.zeros(n, n)
        for a_ in range(n):
            for b_ in range(a_, n):
                val = sp.diff(low[b_], X[a_]) + sp.diff(low[a_], X[b_])
                val -= 2 * sum(G_[c_][a_][b_] * low[c_] for c_ in range(n))
                out[a_, b_] = out[b_, a_] = val
        return out

    def is_killing(self, xi, simplify=None) -> bool:
        s = simplify or self.simplifier
        L = self.lie_derivative(xi)
        return all(s(L[a_, b_]) == 0 for a_ in range(self.dim) for b_ in range(a_, self.dim))

    def conformal_factor(self, xi):
        """If L_xi g = 2 psi g, return psi (else None): xi is a conformal Killing vector."""
        L = self.lie_derivative(xi)
        n = self.dim
        psi = None
        for a_ in range(n):
            for b_ in range(n):
                if self.g[a_, b_] == 0:
                    if self.simplifier(L[a_, b_]) != 0:
                        return None
                    continue
                ratio = self.simplifier(L[a_, b_] / (2 * self.g[a_, b_]))
                if psi is None:
                    psi = ratio
                elif self.simplifier(ratio - psi) != 0:
                    return None
        return psi

    def lower(self, xi) -> list:
        return [self.simplifier(sum(self.g[b_, c_] * sp.sympify(xi[c_]) for c_ in range(self.dim)))
                for b_ in range(self.dim)]

    def killing_tensor_operator(self, K) -> dict:
        """nabla_(a K_bc) for a symmetric covariant tensor K (sympy Matrix).

        Returns {(a, b, c): value} for a <= b <= c, each the full symmetrisation
        (sum over the distinct permutations, unnormalised).
        """
        n, X, G_ = self.dim, self.coords, self.christoffel

        def nabla(a_, b_, c_):
            val = sp.diff(K[b_, c_], X[a_])
            val -= sum(G_[d_][a_][b_] * K[d_, c_] + G_[d_][a_][c_] * K[b_, d_] for d_ in range(n))
            return val

        out = {}
        for a_, b_, c_ in itertools.combinations_with_replacement(range(n), 3):
            perms = set(itertools.permutations((a_, b_, c_)))
            out[(a_, b_, c_)] = sum(nabla(*p) for p in perms)
        return out


def lie_bracket(X_, Y_, coords) -> list:
    """[X, Y]^a = X^b d_b Y^a - Y^b d_b X^a."""
    n = len(coords)
    return [sp.simplify(sum(X_[b_] * sp.diff(Y_[a_], coords[b_]) - Y_[b_] * sp.diff(X_[a_], coords[b_])
                            for b_ in range(n))) for a_ in range(n)]


# ---------------------------------------------------------------------------
# The metrics
# ---------------------------------------------------------------------------

def minkowski() -> Metric:
    return Metric("Minkowski", (t, x, y, z), sp.diag(-1, 1, 1, 1),
                  "flat spacetime of special relativity, Cartesian coordinates")


def minkowski_spherical() -> Metric:
    """Flat spacetime in (t, r, x = cos theta, phi)."""
    g = sp.diag(-1, 1, r**2 / (1 - x**2), r**2 * (1 - x**2))
    return Metric("Minkowski (spherical coordinates)", (t, r, x, phi), g,
                  "flat spacetime in spherical coordinates, x = cos(theta)")


def schwarzschild_angles() -> Metric:
    """Schwarzschild (1916) in its usual coordinates (t, r, theta, phi)."""
    f = 1 - 2 * M / r
    g = sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(theta)**2)
    return Metric("Schwarzschild", (t, r, theta, phi), g, "outside a non-spinning mass M",
                  simplifier=lambda e: sp.simplify(sp.cancel(sp.together(e))))


def schwarzschild() -> Metric:
    """Schwarzschild in (t, r, x = cos theta, phi): rational components."""
    f = 1 - 2 * M / r
    g = sp.diag(-f, 1 / f, r**2 / (1 - x**2), r**2 * (1 - x**2))
    return Metric("Schwarzschild", (t, r, x, phi), g, "outside a non-spinning mass M; x = cos(theta)")


def schwarzschild_cartesian() -> Metric:
    """Schwarzschild in Cartesian-like coordinates: x = r sin(th) cos(ph), etc.

    g = -f dt^2 + [delta_ij + (1/f - 1) n_i n_j] dx^i dx^j with n = x/r.
    """
    rr = sp.sqrt(x**2 + y**2 + z**2)
    f = 1 - 2 * M / rr
    X = (x, y, z)
    g = sp.zeros(4)
    g[0, 0] = -f
    for i in range(3):
        for j in range(3):
            g[i + 1, j + 1] = sp.KroneckerDelta(i, j) + (1 / f - 1) * X[i] * X[j] / rr**2
    return Metric("Schwarzschild (Cartesian)", (t, x, y, z), g,
                  "Schwarzschild in Cartesian-like coordinates",
                  simplifier=lambda e: sp.simplify(e))


def kerr() -> Metric:
    """Kerr (1963) in Boyer-Lindquist coordinates (t, r, x = cos theta, phi).

    Sigma = r^2 + a^2 cos^2(theta), Delta = r^2 - 2 M r + a^2; a = J/M is the
    spin per unit mass (a = M is a maximally spinning hole).
    """
    Sig = r**2 + a**2 * x**2
    Dl = r**2 - 2 * M * r + a**2
    s2 = 1 - x**2
    g = sp.zeros(4)
    g[0, 0] = -(1 - 2 * M * r / Sig)
    g[0, 3] = g[3, 0] = -2 * M * a * r * s2 / Sig
    g[1, 1] = Sig / Dl
    g[2, 2] = Sig / s2
    g[3, 3] = (r**2 + a**2 + 2 * M * a**2 * r * s2 / Sig) * s2
    return Metric("Kerr", (t, r, x, phi), g, "outside a spinning mass M with spin a; x = cos(theta)")


def kerr_angles() -> Metric:
    """Kerr in Boyer-Lindquist (t, r, theta, phi), for display and plotting."""
    Sig = r**2 + a**2 * sp.cos(theta)**2
    Dl = r**2 - 2 * M * r + a**2
    s2 = sp.sin(theta)**2
    g = sp.zeros(4)
    g[0, 0] = -(1 - 2 * M * r / Sig)
    g[0, 3] = g[3, 0] = -2 * M * a * r * s2 / Sig
    g[1, 1] = Sig / Dl
    g[2, 2] = Sig
    g[3, 3] = (r**2 + a**2 + 2 * M * a**2 * r * s2 / Sig) * s2
    return Metric("Kerr", (t, r, theta, phi), g, "Boyer-Lindquist coordinates",
                  simplifier=lambda e: sp.simplify(e))


def flrw(scale=None) -> Metric:
    """Spatially flat FLRW universe, ds^2 = -dt^2 + a(t)^2 (dx^2 + dy^2 + dz^2).

    ``scale`` is an expression in t; by default an arbitrary function a(t).
    """
    A = sp.Function("a")(t) if scale is None else scale
    g = sp.diag(-1, A**2, A**2, A**2)
    return Metric("FLRW (flat)", (t, x, y, z), g, "expanding universe, comoving coordinates",
                  simplifier=lambda e: sp.simplify(e))


def goedel(scale=1, placeholder=False) -> Metric:
    """Goedel's (1949) rotating universe: ds^2 = A^2[-(dt + e^x dz)^2 + dx^2 + dy^2 + e^(2x) dz^2 / 2].

    With ``placeholder`` the exponential is an undefined function E(x), so that
    the jet count can substitute its Taylor series (SymPy would otherwise turn
    exp(x)**2 into exp(2x)).
    """
    E = sp.Function("E")(x) if placeholder else sp.exp(x)
    g = sp.zeros(4)
    g[0, 0] = -1
    g[0, 3] = g[3, 0] = -E
    g[1, 1] = 1
    g[2, 2] = 1
    g[3, 3] = -E**2 + E**2 / 2
    return Metric("Goedel universe", (t, x, y, z), scale**2 * g, "Goedel's rotating dust universe",
                  simplifier=lambda e: sp.simplify(e))


def einstein_static() -> Metric:
    """Einstein's static universe R x S^3: ds^2 = -dt^2 + R^2 [d chi^2 + sin^2 chi dOmega^2].

    Written with w = cos(chi) and x = cos(theta) so the components are rational:
    d chi^2 = dw^2 / (1 - w^2).
    """
    w = sp.Symbol("w", real=True)
    g = sp.diag(-1, 1 / (1 - w**2), (1 - w**2) / (1 - x**2), (1 - w**2) * (1 - x**2))
    return Metric("Einstein static universe", (t, w, x, phi), g, "a closed, static universe (radius 1)")


def plane_wave(A_of_u=None, B_of_u=None) -> Metric:
    """An exact plane gravitational wave in Brinkmann coordinates (u, v, x, y):

        ds^2 = 2 du dv + [A(u) (x^2 - y^2) + 2 B(u) x y] du^2 + dx^2 + dy^2.

    A(u) and B(u) are the two polarisation profiles (plus and cross).
    """
    A_ = sp.Function("A")(u) if A_of_u is None else A_of_u
    B_ = sp.Function("B")(u) if B_of_u is None else B_of_u
    g = sp.zeros(4)
    g[0, 0] = A_ * (x**2 - y**2) + 2 * B_ * x * y
    g[0, 1] = g[1, 0] = 1
    g[2, 2] = 1
    g[3, 3] = 1
    return Metric("plane gravitational wave", (u, v, x, y), g, "exact plane wave (Brinkmann form)",
                  simplifier=lambda e: sp.simplify(e))


def zipoy_voorhees(delta=2) -> Metric:
    """Zipoy-Voorhees (gamma-) metric: a static, axially symmetric vacuum solution
    of Einstein's equations, in prolate spheroidal coordinates (t, X, Y, phi) with
    sigma = 1 (mass delta):

        ds^2 = -F dt^2 + F^-1 [ ((X^2-1)/(X^2-Y^2))^(delta^2) (X^2-Y^2)
               (dX^2/(X^2-1) + dY^2/(1-Y^2)) + (X^2-1)(1-Y^2) dphi^2 ],
        F = ((X-1)/(X+1))^delta.

    delta = 1 is Schwarzschild (X = r/M - 1, Y = cos theta). Other values give a
    body with the same two symmetries as Kerr (time and rotation) but a
    different shape (quadrupole moment). Integer delta keeps every component
    rational. Used as the "no hidden symmetry" contrast to Kerr.
    """
    X_, Y_ = x, y
    F = ((X_ - 1) / (X_ + 1)) ** delta
    W = ((X_**2 - 1) / (X_**2 - Y_**2)) ** (delta**2) * (X_**2 - Y_**2)
    g = sp.diag(-F, W / ((X_**2 - 1) * F), W / ((1 - Y_**2) * F), (X_**2 - 1) * (1 - Y_**2) / F)
    return Metric(f"Zipoy-Voorhees (delta={delta})", (t, X_, Y_, phi), g,
                  "static axisymmetric vacuum spacetime with a non-Kerr quadrupole")
