"""How many symmetries can a spacetime have? An exact count from Taylor series.

The idea, in four steps:

1. A Killing vector field is completely fixed by its value and its first
   derivatives at a single point p: 4 + 16 = 20 numbers in four dimensions.
   (Killing's equation forces all higher derivatives to follow from these.)
   Ten of those combinations are already ruled out at p itself, because the
   symmetric part of the derivative must vanish, so no spacetime has more
   than 10 Killing vectors.

2. Expand everything in a Taylor series about p: the metric (known), and the
   unknown vector field xi^a = sum_alpha c^a_alpha y^alpha, with y = x - p,
   up to degree N. Killing's equation, expanded to degree N - 1, becomes a
   set of linear equations for the coefficients c, with exact rational
   numbers in them when p and the metric's parameters are rational.

3. Every true Killing vector satisfies these truncated equations. So the
   number of independent *first-order data* (the 20 numbers of step 1) that
   the truncated equations allow is an UPPER BOUND on the number of Killing
   vectors, for every N. Increasing N tightens it: the equations at higher
   order feel the curvature, then how the curvature changes, and so on.

4. Once the bound stops dropping, and we have exhibited that many independent
   Killing vectors explicitly (``killing.py``), the count is exact.

The same machinery counts **Killing tensors** (``killing_tensor_bound``):
symmetric K_ab with nabla_(a K_bc) = 0. Each gives a conserved quantity
K_ab u^a u^b that is quadratic in the velocity, like Carter's constant in
Kerr. A Killing 2-tensor is fixed by its derivatives up to second order at a
point, so the projection there is onto the 2-jet.

All arithmetic is exact (Python Fractions; ranks over the rationals with
SymPy's DomainMatrix), so the bounds are theorems, not floating-point guesses.
"""

from __future__ import annotations

import itertools
import time
from fractions import Fraction

import sympy as sp
from sympy import QQ
from sympy.polys.matrices import DomainMatrix

# ---------------------------------------------------------------------------
# Truncated multivariate power series: dict {exponent tuple: Fraction}
# ---------------------------------------------------------------------------


def monomials(n: int, N: int, active=None) -> list[tuple]:
    """Exponent tuples of total degree <= N in n variables (only ``active`` ones used)."""
    act = list(range(n)) if active is None else list(active)
    out = []
    for d in range(N + 1):
        for comb in itertools.combinations_with_replacement(act, d):
            e = [0] * n
            for i in comb:
                e[i] += 1
            out.append(tuple(e))
    return out


def s_mul(A: dict, B: dict, N: int) -> dict:
    out: dict = {}
    for ka, va in A.items():
        da = sum(ka)
        for kb, vb in B.items():
            if da + sum(kb) > N:
                continue
            kk = tuple(p + q for p, q in zip(ka, kb))
            out[kk] = out.get(kk, 0) + va * vb
    return {kk: vv for kk, vv in out.items() if vv != 0}


def s_add(*terms) -> dict:
    out: dict = {}
    for T in terms:
        for kk, vv in T.items():
            out[kk] = out.get(kk, 0) + vv
    return {kk: vv for kk, vv in out.items() if vv != 0}


def s_scale(A: dict, c) -> dict:
    return {kk: vv * c for kk, vv in A.items() if vv * c != 0}


def s_deriv(A: dict, i: int) -> dict:
    out: dict = {}
    for kk, vv in A.items():
        if kk[i] == 0:
            continue
        k2 = list(kk)
        k2[i] -= 1
        out[tuple(k2)] = out.get(tuple(k2), 0) + vv * kk[i]
    return out


def _to_fraction(q) -> Fraction:
    return Fraction(int(q.numerator), int(q.denominator))


def taylor_series(expr, coords, point, N: int, functions=None) -> dict:
    """Exact Taylor series to total degree N of ``expr`` about ``point``.

    ``expr`` must be a rational function of the coordinates once every entry
    of ``functions`` has been replaced: {sympy function application: list of
    Taylor coefficients about the point} (e.g. {a(t): [a0, a1, a2, ...]} for
    a cosmic scale factor, or {exp(x): [1, 1, 1/2, ...]}).
    """
    n = len(coords)
    ys = sp.symbols(f"_y0:{n}")
    e = sp.sympify(expr)
    if functions:
        # Replace each function by its Taylor polynomial in the shifted variable
        # first: afterwards the coordinates can be shifted everywhere safely.
        reps = {}
        for fapp, coeffs in functions.items():
            var = [i for i, c_ in enumerate(coords) if fapp.has(c_)]
            assert len(var) == 1, f"{fapp}: functions of one coordinate only"
            yv = ys[var[0]]
            reps[fapp] = sum(sp.Rational(cf) * yv**j for j, cf in enumerate(coeffs[:N + 1]))
        e = e.subs(reps)
        assert not any(e.has(fapp) for fapp in functions), "function left unsubstituted"
    e = sp.together(e.subs({c_: p_ + y_ for c_, p_, y_ in zip(coords, point, ys)}))
    num, den = sp.fraction(e)
    P = sp.Poly(sp.expand(num), *ys, domain=QQ).as_dict()
    Q = sp.Poly(sp.expand(den), *ys, domain=QQ).as_dict()
    P = {kk: _to_fraction(vv) for kk, vv in P.items() if sum(kk) <= N}
    Q = {kk: _to_fraction(vv) for kk, vv in Q.items() if sum(kk) <= N}
    zero = tuple([0] * n)
    q0 = Q.get(zero, Fraction(0))
    if q0 == 0:
        raise ZeroDivisionError(f"{expr} is singular at {point}")
    rest = {kk: -vv / q0 for kk, vv in Q.items() if kk != zero}
    inv, term = {zero: Fraction(1)}, {zero: Fraction(1)}
    for _ in range(N):  # 1/Q = (1/q0) * sum_k rest^k
        term = s_mul(term, rest, N)
        inv = s_add(inv, term)
    return s_mul(P, s_scale(inv, 1 / q0), N)


# ---------------------------------------------------------------------------
# Exact ranks
# ---------------------------------------------------------------------------


def _rank(entries: dict, nrow: int, cols: list) -> int:
    if not cols:
        return 0
    cmap = {j: i for i, j in enumerate(cols)}
    rows: dict = {}
    for (rr, j), vv in entries.items():
        if j in cmap:
            rows.setdefault(rr, {})[cmap[j]] = QQ(vv.numerator, vv.denominator)
    return DomainMatrix(rows, (nrow, len(cols)), QQ).rank()


def _projected_dimension(entries, nrow, ncol, low_cols) -> dict:
    """dim of the solution space of A c = 0, projected onto the ``low_cols`` coordinates.

    projection = nullity(A) - nullity(A restricted to the other columns).
    """
    high = [j for j in range(ncol) if j not in set(low_cols)]
    rA = _rank(entries, nrow, list(range(ncol)))
    rH = _rank(entries, nrow, high)
    return {"bound": len(low_cols) - (rA - rH), "rank": rA, "rank_high": rH,
            "rows": nrow, "cols": ncol}


# ---------------------------------------------------------------------------
# Killing vectors
# ---------------------------------------------------------------------------


def metric_series(metric, point, N, subs=None, functions=None) -> list:
    g = metric.g.subs(subs or {})
    n = metric.dim
    out = [[None] * n for _ in range(n)]
    for a_ in range(n):
        for b_ in range(a_, n):
            out[a_][b_] = out[b_][a_] = taylor_series(g[a_, b_], metric.coords, point, N, functions)
    return out


def killing_vector_bound(metric, point, N, subs=None, functions=None) -> dict:
    """Upper bound on the number of Killing vectors from the degree-N jet at ``point``."""
    t0 = time.time()
    n = metric.dim
    G = metric_series(metric, point, N, subs, functions)
    dG = [[[s_deriv(G[a_][b_], c_) for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    cols = [(c_, al) for al in monomials(n, N) for c_ in range(n)]
    out_monos = monomials(n, N - 1)
    row_of = {}
    for a_ in range(n):
        for b_ in range(a_, n):
            for gm in out_monos:
                row_of[(a_, b_, gm)] = len(row_of)
    entries = {}
    for j, (c_, al) in enumerate(cols):
        mono = {al: Fraction(1)}
        dmono = [s_deriv(mono, i) for i in range(n)]
        for a_ in range(n):
            for b_ in range(a_, n):
                tot = s_add(s_mul(mono, dG[a_][b_][c_], N - 1),
                            s_mul(G[c_][b_], dmono[a_], N - 1),
                            s_mul(G[a_][c_], dmono[b_], N - 1))
                for kk, vv in tot.items():
                    entries[(row_of[(a_, b_, kk)], j)] = vv
    low = [j for j, (c_, al) in enumerate(cols) if sum(al) <= 1]
    res = _projected_dimension(entries, len(row_of), len(cols), low)
    res.update({"N": N, "seconds": time.time() - t0})
    return res


def killing_vector_sequence(metric, point, N_max, subs=None, functions=None, stop_after=2) -> list:
    """Bounds for N = 1, 2, ... until the bound has been flat ``stop_after`` times."""
    seq, flat = [], 0
    for N in range(1, N_max + 1):
        res = killing_vector_bound(metric, point, N, subs, functions)
        if seq and res["bound"] == seq[-1]["bound"]:
            flat += 1
        else:
            flat = 0
        seq.append(res)
        if N >= 3 and flat >= stop_after:
            break
    return seq


# ---------------------------------------------------------------------------
# Killing tensors (valence 2), invariant under the obvious symmetries
# ---------------------------------------------------------------------------


def killing_tensor_bound(metric, point, N, subs=None, active=(1, 2)) -> dict:
    """Upper bound on the number of Killing 2-tensors that depend only on the
    ``active`` coordinates (for Kerr: r and cos theta, i.e. tensors unchanged
    by time translation and rotation about the axis).

    Every Killing tensor of a stationary, axisymmetric spacetime built from its
    Killing vectors and metric, and Carter's tensor in Kerr, is of this kind.
    """
    t0 = time.time()
    n = metric.dim
    s = subs or {}
    Gam = [[[taylor_series(sp.sympify(metric.christoffel[a_][b_][c_]).subs(s), metric.coords,
                           point, N) for c_ in range(n)] for b_ in range(n)] for a_ in range(n)]
    for a_ in range(n):  # the metric must not depend on the inactive coordinates
        for kk in metric_series(metric, point, 1, subs)[a_][a_]:
            assert all(kk[i] == 0 for i in range(n) if i not in active), "metric not invariant"
    pairs = [(b_, c_) for b_ in range(n) for c_ in range(b_, n)]
    cols = [(pc, al) for al in monomials(n, N, active) for pc in pairs]
    triples = list(itertools.combinations_with_replacement(range(n), 3))
    row_of = {}
    for tr in triples:
        for gm in monomials(n, N - 1, active):
            row_of[(tr, gm)] = len(row_of)
    entries = {}
    for j, ((p_, q_), al) in enumerate(cols):
        mono = {al: Fraction(1)}
        Kc = {}  # nonzero components of the basis tensor
        Kc[(p_, q_)] = Kc[(q_, p_)] = mono

        def comp(b_, c_):
            return Kc.get((b_, c_), {})

        def nabla(a_, b_, c_):
            terms = [s_deriv(comp(b_, c_), a_)]
            for d_ in range(n):
                if comp(d_, c_):
                    terms.append(s_scale(s_mul(Gam[d_][a_][b_], comp(d_, c_), N - 1), -1))
                if comp(b_, d_):
                    terms.append(s_scale(s_mul(Gam[d_][a_][c_], comp(b_, d_), N - 1), -1))
            return s_add(*terms)

        for tr in triples:
            a_, b_, c_ = tr
            tot = s_add(nabla(a_, b_, c_), nabla(b_, c_, a_), nabla(c_, a_, b_))
            for kk, vv in tot.items():
                if sum(kk) <= N - 1:
                    entries[(row_of[(tr, kk)], j)] = vv
    low = [j for j, (pc, al) in enumerate(cols) if sum(al) <= 2]
    res = _projected_dimension(entries, len(row_of), len(cols), low)
    res.update({"N": N, "seconds": time.time() - t0})
    return res


def killing_tensor_sequence(metric, point, N_min, N_max, subs=None, active=(1, 2), stop_after=2):
    seq, flat = [], 0
    for N in range(N_min, N_max + 1):
        res = killing_tensor_bound(metric, point, N, subs, active)
        if seq and res["bound"] == seq[-1]["bound"]:
            flat += 1
        else:
            flat = 0
        seq.append(res)
        if flat >= stop_after:
            break
    return seq
