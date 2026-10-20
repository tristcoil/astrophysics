"""The geometry, derived with computer algebra (SymPy) rather than quoted.

What this module establishes, each as an exact symbolic identity:

1. The Gullstrand-Painleve (GP) metric is the Schwarzschild metric in other
   coordinates: substituting T = t + 4M [sqrt(r/2M) + 1/2 ln|(sqrt(r/2M) - 1)/(sqrt(r/2M) + 1)|]
   turns one into the other.
2. In GP coordinates space (T = const) is flat: dr^2 + r^2 dOmega^2. All the
   curvature is in how that flat space flows inward ("the river").
3. The raindrop u = d_T - beta d_r, beta = sqrt(2M/r), is a geodesic: a body
   falling from rest at infinity, whose own clock reads T.
4. The Kretschmann scalar R_abcd R^abcd = 48 M^2 / r^6: finite at the horizon
   (r = 2M), infinite only at r = 0. The horizon is not a place where anything
   local happens; the singularity is.
5. The tidal tensor felt by any radially falling observer: R^i_0j0 =
   diag(-2M/r^3, M/r^3, M/r^3): stretching along the fall at 2GM/r^3 per metre,
   squeezing sideways at half that. The same for every radial faller (boosts
   along r leave it unchanged), which we verify for a body dropped from rest.

``christoffel_functions`` lambdifies the Christoffel symbols of the equatorial
GP metric for geodesic.py, which integrates light rays numerically to check
raymap.py independently.
"""

from __future__ import annotations

from functools import lru_cache

import sympy as sp

T, t, th, ph = sp.symbols("T t theta phi", real=True)
r = sp.symbols("r", positive=True)
M = sp.symbols("M", positive=True)
COORDS = (T, r, th, ph)


def gp_metric(mass=M) -> sp.Matrix:
    b = sp.sqrt(2 * mass / r)
    return sp.Matrix([[-(1 - b**2), b, 0, 0],
                      [b, 1, 0, 0],
                      [0, 0, r**2, 0],
                      [0, 0, 0, r**2 * sp.sin(th) ** 2]])


def schwarzschild_metric(mass=M) -> sp.Matrix:
    f = 1 - 2 * mass / r
    return sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2)


def christoffel(g: sp.Matrix, coords=COORDS):
    n = len(coords)
    ginv = sp.simplify(g.inv())
    gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                expr = sum(ginv[a, d] * (sp.diff(g[d, b], coords[c]) + sp.diff(g[d, c], coords[b])
                                         - sp.diff(g[b, c], coords[d])) for d in range(n)) / 2
                expr = sp.simplify(expr)
                gam[a][b][c] = gam[a][c][b] = expr
    return gam


def riemann(gam, coords=COORDS):
    """R^a_bcd = d_c Gam^a_db - d_d Gam^a_cb + Gam^a_ce Gam^e_db - Gam^a_de Gam^e_cb."""
    n = len(coords)
    R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(c + 1, n):
                    expr = (sp.diff(gam[a][d][b], coords[c]) - sp.diff(gam[a][c][b], coords[d])
                            + sum(gam[a][c][e] * gam[e][d][b] - gam[a][d][e] * gam[e][c][b] for e in range(n)))
                    expr = sp.simplify(expr)
                    R[a][b][c][d] = expr
                    R[a][b][d][c] = -expr
    return R


def kretschmann(g, R, coords=COORDS):
    n = len(coords)
    ginv = g.inv()
    # lower the first index, raise the other three
    low = [[[[sum(g[a, e] * R[e][b][c][d] for e in range(n)) for d in range(n)] for c in range(n)]
            for b in range(n)] for a in range(n)]
    total = 0
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    if low[a][b][c][d] == 0:
                        continue
                    up = sum(ginv[a, p] * ginv[b, q] * ginv[c, s] * ginv[d, w] * low[p][q][s][w]
                             for p in range(n) for q in range(n) for s in range(n) for w in range(n)
                             if low[p][q][s][w] != 0)
                    total += low[a][b][c][d] * up
    return sp.simplify(total)


def tidal_tensor(g, R, u, basis):
    """E_ij = R^i_{0 j 0} in an orthonormal frame: e_0 = u, e_i = basis[i] (vectors in coordinates)."""
    n = 4
    frame = [u] + list(basis)
    out = sp.zeros(3, 3)
    # R(e_i, u, e_j, u) with the first index lowered: R_abcd e_i^a u^b e_j^c u^d
    low = [[[[sum(g[a, e] * R[e][b][c][d] for e in range(n)) for d in range(n)] for c in range(n)]
            for b in range(n)] for a in range(n)]
    for i in range(3):
        for j in range(3):
            ei, ej = frame[i + 1], frame[j + 1]
            val = sum(low[a][b][c][d] * ei[a] * u[b] * ej[c] * u[d]
                      for a in range(n) for b in range(n) for c in range(n) for d in range(n)
                      if low[a][b][c][d] != 0)
            out[i, j] = sp.simplify(val)
    return out


def derive() -> dict:
    """Run every symbolic check; returns printable strings and booleans."""
    res = {}
    g = gp_metric()
    gs = schwarzschild_metric()

    # 1. GP <- Schwarzschild: dT = dt + (dF/dr) dr with F as in the docstring.
    x = sp.sqrt(r / (2 * M))
    F = 4 * M * (x + sp.log((x - 1) / (x + 1)) / 2)
    dF = sp.simplify(sp.diff(F, r))
    # Jacobian: T = t + F(r). Pull the GP metric back to (t, r): g_s = J^T g J.
    J = sp.Matrix([[1, dF, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    pulled = sp.simplify(J.T * g * J)
    res["gp_equals_schwarzschild"] = bool(sp.simplify(pulled - gs) == sp.zeros(4, 4))
    res["dT_dt_term"] = str(dF)

    # 2. Flat space at T = const.
    res["spatial_metric_flat"] = bool(sp.simplify(g[1:, 1:] - sp.diag(1, r**2, r**2 * sp.sin(th) ** 2))
                                      == sp.zeros(3, 3))

    gam = christoffel(g)
    # 3. The raindrop is a geodesic: a^mu = u^nu d_nu u^mu + Gam^mu_ab u^a u^b = 0.
    b = sp.sqrt(2 * M / r)
    u = [1, -b, 0, 0]
    acc = [sp.simplify(sum(u[n] * sp.diff(u[m], COORDS[n]) for n in range(4))
                       + sum(gam[m][a][c] * u[a] * u[c] for a in range(4) for c in range(4))) for m in range(4)]
    res["raindrop_is_geodesic"] = all(a == 0 for a in acc)

    # 4. Curvature.
    R = riemann(gam)
    K = kretschmann(g, R)
    res["kretschmann"] = str(K)
    res["kretschmann_ok"] = bool(sp.simplify(K - 48 * M**2 / r**6) == 0)

    # 5. Tidal tensor for the raindrop (orthonormal: e_r = d_r, e_th = d_th / r, e_ph = d_ph / (r sin th)).
    basis = ([0, 1, 0, 0], [0, 0, 1 / r, 0], [0, 0, 0, 1 / (r * sp.sin(th))])
    E = tidal_tensor(g, R, u, basis)
    res["tidal_rain"] = str(E)
    res["tidal_rain_ok"] = bool(sp.simplify(E - sp.diag(-2 * M / r**3, M / r**3, M / r**3)) == sp.zeros(3, 3))

    # ... and for a body dropped from rest at r0 (energy E0 = sqrt(1 - 2M/r0)): its
    # 4-velocity u_d and the radial unit vector orthogonal to it.
    r0 = sp.symbols("r_0", positive=True)
    E0 = sp.sqrt(1 - 2 * M / r0)
    ur = -sp.sqrt(E0**2 - 1 + 2 * M / r)
    w = (1 - E0**2) / (b * E0 - ur)          # gamma v relative to the rain (horizon-regular)
    uT = E0 + b * w
    ud = [uT, ur, 0, 0]
    # radial unit vector orthogonal to u_d: n = gamma (v e_0 + e_1) in rain-frame terms
    gam_rel, gv = uT, w
    n_vec = [gv, -b * gv + gam_rel, 0, 0]      # gv e_0 + gamma e_1 with e_0 = (1, -b), e_1 = (0, 1)
    basis_d = (n_vec, [0, 0, 1 / r, 0], [0, 0, 0, 1 / (r * sp.sin(th))])
    Ed = tidal_tensor(g, R, ud, basis_d)
    res["tidal_drop_ok"] = bool(sp.simplify(Ed - sp.diag(-2 * M / r**3, M / r**3, M / r**3)) == sp.zeros(3, 3))
    return res


@lru_cache(maxsize=1)
def christoffel_functions():
    """Numerical Christoffel symbols of the equatorial GP metric (M = 1), coordinates (T, r, phi).

    Returns f(r) -> 3x3x3 nested list Gam[a][b][c] for a, b, c in (T, r, phi).
    """
    g = gp_metric(1).subs(th, sp.pi / 2)
    g3 = sp.Matrix([[g[0, 0], g[0, 1], 0], [g[1, 0], g[1, 1], 0], [0, 0, r**2]])
    coords = (T, r, ph)
    gam = christoffel(g3, coords)
    flat = [gam[a][b][c] for a in range(3) for b in range(3) for c in range(3)]
    f = sp.lambdify(r, flat, "math")
    return f


if __name__ == "__main__":
    for key, val in derive().items():
        print(f"{key:28s} {val}")
