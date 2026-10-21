"""Find the Killing vectors themselves, check them, and find Kerr's hidden symmetry.

``jets.py`` proves how many Killing vectors a spacetime can have at most. Here
we find them in closed form and check that each really is one.

**Finding them.** We look for every vector field whose components are
polynomials of degree <= 2 in Cartesian-like coordinates (t, x, y, z): 4
components x 15 monomials = 60 unknown coefficients. Killing's equation,
evaluated at a few dozen points with exact rational coordinates, gives a
linear system for them, solved exactly. For flat spacetime this search is
complete by itself (Killing's equation in flat space forces second
derivatives to vanish, so the components are at most linear). For the curved
spacetimes the jet bound closes the gap: if the search finds as many as the
bound allows, there are no others, polynomial or not.

**Checking them.** Each solution is converted to the usual spherical (or
Boyer-Lindquist) coordinates and substituted back into L_xi g = 0 with the
full symbolic metric.

**Their algebra.** The commutator [xi_1, xi_2] of two Killing vectors is again
a Killing vector. The table of commutators names the symmetry group: the
Poincare group for flat spacetime, rotations x time translations for
Schwarzschild, two commuting symmetries for Kerr, rotations and translations
of space (the Euclidean group) for an expanding universe.

**Kerr's hidden symmetry.** Kerr has only two Killing vectors, yet its orbits
have a third constant of motion (Carter 1968). It comes from a Killing tensor
K_ab, which we verify in the form of Penrose and Floyd: K = Y Y^T, built from
the Killing-Yano 2-form Y (nabla_(a Y_b)c = 0). We also show that K is not
built from the Killing vectors and the metric, and that for a = 0 it becomes
the total angular momentum squared, which in Schwarzschild *is* built from the
three rotations.
"""

from __future__ import annotations

import itertools
import random
from dataclasses import dataclass, field

import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

import geometry as geo
from geometry import M, a, t, x, y, z, r, theta, phi, H

R_ = sp.Symbol("R", positive=True)  # auxiliary radius in the Cartesian charts

# ---------------------------------------------------------------------------
# Charts with exact evaluation at rational points
# ---------------------------------------------------------------------------


@dataclass
class Chart:
    """A metric in Cartesian-like coordinates, evaluable exactly at sample points.

    ``g`` may contain an auxiliary symbol R_ (a radius defined implicitly by
    the coordinates); ``aux_grad`` gives dR/dx^i as expressions, and each
    sample point fixes R_ consistently.
    """

    name: str
    coords: tuple
    g: sp.Matrix
    points: list
    aux_grad: tuple | None = None
    _dg: list = field(default=None, repr=False)

    def dg(self):
        if self._dg is None:
            n = len(self.coords)
            out = []
            for i, c_ in enumerate(self.coords):
                Dg = sp.zeros(n, n)
                for p_ in range(n):
                    for q_ in range(n):
                        e = sp.diff(self.g[p_, q_], c_)
                        if self.aux_grad is not None:
                            e += sp.diff(self.g[p_, q_], R_) * self.aux_grad[i]
                        Dg[p_, q_] = e
                out.append(Dg)
            self._dg = out
        return self._dg


def pythagorean_points(n_points: int, seed: int = 3, avoid=(2,)):
    """Points (t, x, y, z) with rational r = sqrt(x^2 + y^2 + z^2) (and r not in ``avoid``)."""
    rng = random.Random(seed)
    pts = []
    while len(pts) < n_points:
        m_, n_, p_, q_ = (rng.randint(-4, 4) for _ in range(4))
        rr = m_**2 + n_**2 + p_**2 + q_**2
        if rr == 0:
            continue
        X = (m_**2 + n_**2 - p_**2 - q_**2, 2 * (m_ * q_ + n_ * p_), 2 * (n_ * q_ - m_ * p_))
        s = sp.Rational(rng.randint(2, 9), rng.randint(1, 3))  # rescale: vary r
        Xs = [s * c_ for c_ in X]
        rs = s * rr
        if rs in avoid or any(abs(c_) == rs for c_ in Xs):
            continue
        tt = sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))
        pts.append({t: tt, x: Xs[0], y: Xs[1], z: Xs[2], R_: rs})
    return pts


def chart_minkowski():
    rng = random.Random(5)
    pts = [{t: sp.Rational(rng.randint(-9, 9), 3), x: sp.Rational(rng.randint(-9, 9), 2),
            y: sp.Rational(rng.randint(-9, 9), 5), z: sp.Rational(rng.randint(-9, 9), 7)}
           for _ in range(30)]
    return Chart("Minkowski", (t, x, y, z), sp.diag(-1, 1, 1, 1), pts)


def chart_schwarzschild(mass=1):
    """Schwarzschild in Cartesian-like coordinates, R_ = sqrt(x^2 + y^2 + z^2)."""
    X = (x, y, z)
    g = sp.zeros(4)
    g[0, 0] = -(1 - 2 * mass / R_)
    for i in range(3):
        for j in range(3):
            g[i + 1, j + 1] = sp.KroneckerDelta(i, j) + 2 * mass / (R_ - 2 * mass) * X[i] * X[j] / R_**2
    grad = (0, x / R_, y / R_, z / R_)
    return Chart("Schwarzschild", (t, x, y, z), g, pythagorean_points(30, avoid=(2 * mass,)), grad)


def chart_kerr_schild(mass=1, spin=sp.Rational(3, 5), n_points=36, seed=11):
    """Kerr in Kerr-Schild coordinates: g = eta + H l l, with R_ defined by
    (x^2 + y^2)/(R^2 + a^2) + z^2/R^2 = 1.
    """
    a_ = spin
    Hf = 2 * mass * R_**3 / (R_**4 + a_**2 * z**2)
    l_ = [1, (R_ * x + a_ * y) / (R_**2 + a_**2), (R_ * y - a_ * x) / (R_**2 + a_**2), z / R_]
    g = sp.diag(-1, 1, 1, 1) + Hf * sp.Matrix(4, 4, lambda i, j: l_[i] * l_[j])
    F = R_**4 - (x**2 + y**2 + z**2 - a_**2) * R_**2 - a_**2 * z**2
    FR = sp.diff(F, R_)
    grad = (0, -sp.diff(F, x) / FR, -sp.diff(F, y) / FR, -sp.diff(F, z) / FR)
    rng = random.Random(seed)
    pts = []
    while len(pts) < n_points:
        rr = sp.Rational(rng.randint(5, 40), rng.randint(1, 3))
        s1, s2 = sp.Rational(rng.randint(1, 9), rng.randint(1, 9)), sp.Rational(rng.randint(-9, 9), rng.randint(1, 9))
        st, ct = 2 * s1 / (1 + s1**2), (1 - s1**2) / (1 + s1**2)
        sph, cph = 2 * s2 / (1 + s2**2), (1 - s2**2) / (1 + s2**2)
        X = (rr * cph - a_ * sph) * st
        Y = (rr * sph + a_ * cph) * st
        Z = rr * ct
        pts.append({t: sp.Rational(rng.randint(-9, 9), 4), x: X, y: Y, z: Z, R_: rr})
    return Chart("Kerr (Kerr-Schild)", (t, x, y, z), g, pts, grad)


def chart_flrw_matter():
    """Flat FLRW in the matter era, a = t^(2/3), at t = s^3 so that a = s^2 is rational."""
    s = sp.Symbol("s", positive=True)
    rng = random.Random(7)
    pts = []
    for _ in range(30):
        sv = sp.Rational(rng.randint(2, 9), rng.randint(1, 4))
        pts.append({t: sv**3, x: sp.Rational(rng.randint(-9, 9), 2), y: sp.Rational(rng.randint(-9, 9), 3),
                    z: sp.Rational(rng.randint(-9, 9), 5)})
    g = sp.diag(-1, t**sp.Rational(4, 3), t**sp.Rational(4, 3), t**sp.Rational(4, 3))
    return Chart("FLRW, matter era", (t, x, y, z), g, pts)


def chart_de_sitter_conformal(hubble=1):
    """De Sitter space in conformal time eta < 0: ds^2 = (-d eta^2 + dx^2) / (H eta)^2.

    It is the flat-slicing FLRW universe with a = e^(Ht), eta = -e^(-Ht)/H:
    the far future of a universe dominated by a cosmological constant.
    """
    eta = t  # reuse the symbol; it plays the role of conformal time here
    rng = random.Random(9)
    pts = [{t: -sp.Rational(rng.randint(1, 9), rng.randint(1, 4)), x: sp.Rational(rng.randint(-9, 9), 2),
            y: sp.Rational(rng.randint(-9, 9), 3), z: sp.Rational(rng.randint(-9, 9), 5)} for _ in range(30)]
    g = sp.diag(-1, 1, 1, 1) / (hubble * eta) ** 2
    return Chart("de Sitter (conformal time)", (t, x, y, z), g, pts)


# ---------------------------------------------------------------------------
# The polynomial search
# ---------------------------------------------------------------------------


def _monomials(coords, degree):
    out = [sp.Integer(1)]
    for d_ in range(1, degree + 1):
        for comb in itertools.combinations_with_replacement(coords, d_):
            out.append(sp.Mul(*comb))
    return out


def solve_polynomial(chart: Chart, degree: int = 2) -> dict:
    """All Killing vectors of ``chart`` with polynomial components of degree <= ``degree``.

    Returns the basis (list of component lists) and the size of the search.
    """
    n = len(chart.coords)
    monos = _monomials(chart.coords, degree)
    unknowns = [(rho, k_) for rho in range(n) for k_ in range(len(monos))]
    dG = chart.dg()
    rows = []
    for P in chart.points:
        gP = chart.g.subs(P)
        dgP = [D.subs(P) for D in dG]
        mv = [m_.subs(P) for m_ in monos]
        mg = [[sp.diff(m_, c_).subs(P) for c_ in chart.coords] for m_ in monos]
        for a_ in range(n):
            for b_ in range(a_, n):
                row = []
                for rho, k_ in unknowns:
                    val = mv[k_] * dgP[rho][a_, b_] + gP[rho, b_] * mg[k_][a_] + gP[a_, rho] * mg[k_][b_]
                    row.append(QQ(sp.Rational(val).p, sp.Rational(val).q))
                rows.append(row)
    A = DomainMatrix(rows, (len(rows), len(unknowns)), QQ)
    null = A.nullspace().to_Matrix()
    basis = []
    for vec in null.tolist():
        comps = [sp.Integer(0)] * n
        for (rho, k_), c_ in zip(unknowns, vec):
            if c_ != 0:
                comps[rho] += sp.Rational(c_) * monos[k_]
        basis.append([sp.expand(c_) for c_ in comps])
    return {"basis": basis, "n_unknowns": len(unknowns), "n_equations": len(rows),
            "rank": len(unknowns) - len(basis), "n_points": len(chart.points)}


def coefficient_vector(field_, coords, degree=2):
    monos = _monomials(coords, degree)
    vec = []
    for comp in field_:
        P = sp.Poly(sp.expand(comp), *coords)
        d = P.as_dict()
        for m_ in monos:
            mono = sp.Poly(m_, *coords).as_dict()
            (exps,) = mono.keys()
            vec.append(sp.Rational(d.get(exps, 0)))
    return vec


def same_span(found: list, named: dict, coords, degree=2) -> dict:
    """Do the ``named`` fields span exactly the space the solver found?"""
    F = sp.Matrix([coefficient_vector(f_, coords, degree) for f_ in found])
    Nm = sp.Matrix([coefficient_vector(f_, coords, degree) for f_ in named.values()])
    rF, rN = F.rank(), Nm.rank()
    rFN = F.col_join(Nm).rank()
    return {"dim_found": rF, "dim_named": rN, "dim_union": rFN,
            "equal": rF == rN == rFN == len(named)}


# ---------------------------------------------------------------------------
# Named Killing vectors (contravariant components in the given coordinates)
# ---------------------------------------------------------------------------

def minkowski_named() -> dict:
    P = {f"P_{c_}": [int(i == k_) for i in range(4)] for k_, c_ in enumerate("txyz")}
    J = {"J_x": [0, 0, -z, y], "J_y": [0, z, 0, -x], "J_z": [0, -y, x, 0]}
    K = {"K_x": [x, t, 0, 0], "K_y": [y, 0, t, 0], "K_z": [z, 0, 0, t]}
    return {**P, **J, **K}


def rotations_cartesian() -> dict:
    return {"J_x": [0, 0, -z, y], "J_y": [0, z, 0, -x], "J_z": [0, -y, x, 0]}


def rotations_spherical() -> dict:
    """The three rotations in (t, r, theta, phi)."""
    return {"J_x": [0, 0, -sp.sin(phi), -sp.cot(theta) * sp.cos(phi)],
            "J_y": [0, 0, sp.cos(phi), -sp.cot(theta) * sp.sin(phi)],
            "J_z": [0, 0, 0, 1]}


def de_sitter_named() -> dict:
    """The 10 Killing vectors of de Sitter space in conformal coordinates (eta, x, y, z)."""
    X = (x, y, z)
    named = {f"P_{c_}": [0] + [int(i == k_) for i in range(3)] for k_, c_ in enumerate("xyz")}
    named.update({k_: v_ for k_, v_ in rotations_cartesian().items()})
    named["D"] = [t, x, y, z]  # dilatation: eta -> lam eta, x -> lam x (cosmic-time translation)
    rho2 = x**2 + y**2 + z**2 - t**2
    for k_, c_ in enumerate("xyz"):  # "special conformal" combinations that preserve eta = 0
        comp = [2 * X[k_] * t] + [2 * X[k_] * X[i] - (rho2 if i == k_ else 0) for i in range(3)]
        named[f"S_{c_}"] = comp
    return named


# ---------------------------------------------------------------------------
# Checking by substitution, and the algebra
# ---------------------------------------------------------------------------

def verify(metric: geo.Metric, fields: dict, simplify=None) -> dict:
    """Substitute each field into L_xi g; True if every component simplifies to 0."""
    s = simplify or (lambda e: sp.simplify(e))
    out = {}
    for name, f_ in fields.items():
        L = metric.lie_derivative(f_)
        out[name] = all(s(L[i, j]) == 0 for i in range(metric.dim) for j in range(i, metric.dim))
    return out


def structure_constants(fields: dict, coords) -> dict:
    """[X_i, X_j] = sum_k c_ij^k X_k: returns the non-zero commutators as strings."""
    names = list(fields)
    vecs = [fields[n_] for n_ in names]
    cs = sp.symbols(f"c0:{len(names)}")
    table = {}
    for i, j in itertools.combinations(range(len(names)), 2):
        br = geo.lie_bracket(vecs[i], vecs[j], coords)
        if all(sp.simplify(c_) == 0 for c_ in br):
            continue
        comb = [sum(cs[k_] * vecs[k_][m_] for k_ in range(len(names))) for m_ in range(len(coords))]
        eqs = []
        for m_ in range(len(coords)):
            diff = sp.expand(sp.simplify(comb[m_] - br[m_]))
            if diff == 0:
                continue
            free = [s_ for s_ in coords if diff.has(s_)]
            eqs.extend(sp.Poly(diff, *free).coeffs() if free else [diff])
        sol = sp.solve(eqs, cs, dict=True)
        if not sol:
            table[f"[{names[i]}, {names[j]}]"] = "not in span"
            continue
        sol = sol[0]
        terms = []
        for k_, n_ in enumerate(names):
            c_ = sp.simplify(sol.get(cs[k_], 0))
            if c_ != 0:
                terms.append(f"{'+' if c_ > 0 else '-'}{'' if abs(c_) == 1 else abs(c_)}{n_}")
        table[f"[{names[i]}, {names[j]}]"] = " ".join(terms).lstrip("+")
    return table


def norms(metric: geo.Metric, fields: dict) -> dict:
    """g(xi, xi) for each field: negative = timelike, zero = null, positive = spacelike."""
    out = {}
    for name, f_ in fields.items():
        out[name] = sp.simplify(sum(metric.g[i, j] * f_[i] * f_[j] for i in range(metric.dim)
                                    for j in range(metric.dim)))
    return out


# ---------------------------------------------------------------------------
# FLRW: why there is no time-translation symmetry
# ---------------------------------------------------------------------------

def flrw_killing_equations() -> dict:
    """Killing's equation for a general vector field in flat FLRW with arbitrary a(t).

    The tt component gives d_t xi^t = 0; the ij components give
    d_i xi^j + d_j xi^i = -2 (a'/a) xi^t delta_ij. Differentiating the latter
    in time and using the ti components shows that a nonzero xi^t needs
    a^2 d/dt(a'/a) to be constant, which fails for any matter or radiation
    content. ``flrw_condition`` returns that combination for given a(t).
    """
    met = geo.flrw()
    fs = [sp.Function(f"xi{i}")(t, x, y, z) for i in range(4)]
    L = met.lie_derivative(fs)
    return {"metric": met, "fields": fs, "equations": {f"({i}{j})": sp.simplify(L[i, j])
                                                       for i in range(4) for j in range(i, 4)}}


def flrw_condition(scale) -> sp.Expr:
    """a^2 d/dt (a'/a): constant iff a timelike-ish Killing vector component can exist."""
    return sp.simplify(scale**2 * sp.diff(sp.diff(scale, t) / scale, t))


# ---------------------------------------------------------------------------
# Kerr: the Killing-Yano tensor and Carter's Killing tensor
# ---------------------------------------------------------------------------

def kerr_killing_yano() -> sp.Matrix:
    """Penrose-Floyd Killing-Yano 2-form of Kerr in (t, r, x = cos theta, phi):

    Y = a cos(th) dr ^ (dt - a sin^2(th) dphi) + r sin(th) dth ^ ((r^2 + a^2) dphi - a dt),

    and with r sin(th) dth = -r dx:  Y = a x dr^(dt - a(1-x^2)dphi) - r dx^((r^2+a^2)dphi - a dt).
    """
    Y = sp.zeros(4)

    def wedge(alpha, beta):  # (alpha ^ beta)_ij = alpha_i beta_j - alpha_j beta_i
        return sp.Matrix(4, 4, lambda i, j: alpha[i] * beta[j] - alpha[j] * beta[i])

    dt_, dr_, dx_, dph = ([int(i == k_) for i in range(4)] for k_ in range(4))
    one = [dt_[i] - a * (1 - x**2) * dph[i] for i in range(4)]
    two = [(r**2 + a**2) * dph[i] - a * dt_[i] for i in range(4)]
    Y += a * x * wedge(dr_, one)
    Y += -r * wedge(dx_, two)
    return Y


def killing_yano_check(metric: geo.Metric, Y: sp.Matrix) -> bool:
    """nabla_a Y_bc + nabla_b Y_ac = 0 for all a, b, c."""
    n, X, G_ = metric.dim, metric.coords, metric.christoffel

    def nabla(a_, b_, c_):
        return (sp.diff(Y[b_, c_], X[a_]) - sum(G_[d_][a_][b_] * Y[d_, c_] + G_[d_][a_][c_] * Y[b_, d_]
                                               for d_ in range(n)))
    for a_, b_, c_ in itertools.product(range(n), repeat=3):
        if a_ > b_:
            continue
        if sp.cancel(sp.together(nabla(a_, b_, c_) + nabla(b_, a_, c_))) != 0:
            return False
    return True


def carter_tensor(metric: geo.Metric, Y: sp.Matrix) -> sp.Matrix:
    """K_ab = Y_ac Y_bd g^cd (Penrose & Floyd 1973)."""
    return (Y * metric.ginv * Y.T).applyfunc(lambda e: sp.cancel(sp.together(e)))


def killing_tensor_check(metric: geo.Metric, K: sp.Matrix) -> bool:
    ops = metric.killing_tensor_operator(K)
    return all(sp.cancel(sp.together(v_)) == 0 for v_ in ops.values())


def reducibility(metric: geo.Metric, K: sp.Matrix, killing: dict, sample_points: list) -> dict:
    """Is K a constant combination of g and symmetrised products of Killing vectors?

    Compares the components at ``sample_points``: if adding K raises the rank of
    {g, xi_i xi_j}, it is irreducible -- a genuinely new, hidden symmetry.
    """
    lows = {n_: metric.lower(f_) for n_, f_ in killing.items()}
    names = list(lows)
    basis = {"g": metric.g}
    for i, j in itertools.combinations_with_replacement(range(len(names)), 2):
        li, lj = lows[names[i]], lows[names[j]]
        basis[f"{names[i]}*{names[j]}"] = sp.Matrix(4, 4, lambda p_, q_: (li[p_] * lj[q_] + li[q_] * lj[p_]) / 2)

    def vec(T):
        return [T[p_, q_].subs(P) for P in sample_points for p_ in range(4) for q_ in range(p_, 4)]
    B = sp.Matrix([vec(T) for T in basis.values()])
    rB = B.rank()
    rBK = B.col_join(sp.Matrix([vec(K)])).rank()
    return {"n_products": len(basis), "rank_products": rB, "rank_with_K": rBK,
            "irreducible": rBK > rB}
