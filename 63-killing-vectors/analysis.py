"""Every computation of project N6, each cached in data/cache/ (delete or --force to redo).

Sections:
  symbolic    Killing vectors: exact jet bounds, explicit solutions, checks, algebra;
              Killing tensors and Kerr's hidden symmetry
  noether     numerical geodesics: the Killing quantities stay constant
  orbits      effective potentials, special radii against spin, precession
  s2          S2's energy and angular momentum (H1's data and posterior)
  shadows     critical curves against I5's rings
  thermo      horizon thermodynamics for the atlas's black holes, GW250114
  cosmology   photon redshift and peculiar velocities; D7's T(z)
  chaos       Poincare sections: Kerr against Zipoy-Voorhees
"""

from __future__ import annotations

import pickle
import time
from fractions import Fraction

import numpy as np
import sympy as sp

import config as C

R = sp.Rational


def cached(name: str, fn, force: bool = False):
    path = C.CACHE_DIR / f"{name}.pkl"
    if path.exists() and not force:
        return pickle.loads(path.read_bytes())
    t0 = time.time()
    out = fn()
    if isinstance(out, dict):
        out["_seconds"] = time.time() - t0
    path.write_bytes(pickle.dumps(out))
    return out


# ---------------------------------------------------------------------------
# Symbolic: Killing vectors and tensors
# ---------------------------------------------------------------------------

def _binomial_coeffs(p, n):
    out, c = [], Fraction(1)
    for k in range(n):
        out.append(c)
        c = c * (p - k) / (k + 1)
    return out


def lcdm_taylor(n: int = 9) -> list:
    """Taylor coefficients of a(t) about today for Planck 2018 LCDM, time in 1/H0.

    From (a da/dt)^2 = Omega_r + Omega_m a + Omega_L a^4 (H0 = 1), solved as a
    power series with exact rational arithmetic in the (rounded) parameters.
    """
    Om = Fraction(str(C.PLANCK18["Omega_m"]))
    Or = Fraction(str(C.PLANCK18["Omega_r"]))
    OL = 1 - Om - Or
    # recursion for a = sum c_k t^k with c0 = 1: differentiate a' = F(a) repeatedly via sympy
    A = sp.Symbol("A")
    F = sp.sqrt(sp.Rational(Or.numerator, Or.denominator) / A**2 + sp.Rational(Om.numerator, Om.denominator) / A
                + sp.Rational(OL.numerator, OL.denominator) * A**2)
    derivs = [A, F]
    for _ in range(n - 2):
        derivs.append(sp.diff(derivs[-1], A) * F)
    coeffs = []
    for k, d_ in enumerate(derivs):
        val = sp.N(d_.subs(A, 1), 40) / sp.factorial(k)
        coeffs.append(Fraction(str(val)).limit_denominator(10**18))
    return coeffs


def symbolic() -> dict:
    import geometry as geo
    import jets
    import killing as kl
    from math import factorial

    out = {"vector_bounds": {}, "explicit": {}, "checks": {}, "tensor_bounds": {}}
    NMAX = 6
    E_ = sp.Function("E")(geo.x)
    A_ = sp.Function("a")(geo.t)
    Au, Bu = sp.Function("A")(geo.u), sp.Function("B")(geo.u)
    flrw_pt = (0, R(1, 3), R(2, 7), R(-1, 5))
    cases = [
        ("Minkowski", geo.minkowski(), (R(1, 3), R(2, 7), R(-1, 5), R(3, 11)), {}, None),
        ("Minkowski (spherical coordinates)", geo.minkowski_spherical(), (0, 3, R(1, 3), 0), {}, None),
        ("Schwarzschild", geo.schwarzschild(), (0, R(7, 2), R(1, 3), 0), {geo.M: 1}, None),
        ("Kerr", geo.kerr(), (0, R(7, 2), R(1, 3), 0), {geo.M: 1, geo.a: R(3, 5)}, None),
        ("FLRW, generic a(t)", geo.flrw(), flrw_pt, {},
         {A_: [1, R(7, 10), R(-1, 3), R(2, 9), R(-1, 7), R(3, 11), R(5, 13), R(-2, 17), R(1, 19)]}),
        ("FLRW, matter era", geo.flrw(), flrw_pt, {}, {A_: _binomial_coeffs(Fraction(2, 3), NMAX + 3)}),
        ("FLRW, radiation era", geo.flrw(), flrw_pt, {}, {A_: _binomial_coeffs(Fraction(1, 2), NMAX + 3)}),
        ("FLRW, Planck 2018 LCDM", geo.flrw(), flrw_pt, {}, {A_: lcdm_taylor(NMAX + 3)}),
        ("de Sitter (FLRW form)", geo.flrw(), flrw_pt, {}, {A_: [Fraction(1, factorial(k_)) for k_ in range(NMAX + 3)]}),
        ("Goedel universe", geo.goedel(placeholder=True), (0, 0, 0, 0), {},
         {E_: [Fraction(1, factorial(k_)) for k_ in range(NMAX + 3)]}),
        ("Einstein static universe", geo.einstein_static(), (0, R(1, 3), R(1, 5), 0), {}, None),
        ("plane gravitational wave", geo.plane_wave(), (0, 0, R(1, 3), R(1, 5)), {},
         {Au: [R(1, 2), R(1, 3), R(-1, 5), R(2, 7), R(1, 11), R(-3, 13), R(1, 17), R(2, 19), R(1, 23)],
          Bu: [R(1, 3), R(-1, 7), R(1, 4), R(-2, 9), R(1, 5), R(1, 23), R(-1, 29), R(3, 31), R(1, 37)]}),
        ("Zipoy-Voorhees (delta = 2)", geo.zipoy_voorhees(2), (0, 3, R(1, 3), 0), {}, None),
    ]
    for name, met, pt, subs, fn in cases:
        t0 = time.time()
        seq = jets.killing_vector_sequence(met, pt, NMAX, subs, fn, stop_after=1 if name == "Kerr" else 2)
        out["vector_bounds"][name] = {"point": [str(p_) for p_ in pt], "sequence": [s_["bound"] for s_ in seq],
                                      "sizes": [(s_["rows"], s_["cols"]) for s_ in seq],
                                      "final": seq[-1]["bound"], "known": C.KILLING_COUNTS.get(name),
                                      "seconds": time.time() - t0}
        print(f"  jets {name}: {out['vector_bounds'][name]['sequence']} ({time.time() - t0:.0f} s)", flush=True)

    # Explicit Killing vectors by the polynomial search, and checks by substitution
    rot_c = kl.rotations_cartesian()
    charts = [
        ("Minkowski", kl.chart_minkowski(), kl.minkowski_named()),
        ("Schwarzschild", kl.chart_schwarzschild(), {"P_t": [1, 0, 0, 0], **rot_c}),
        ("Kerr", kl.chart_kerr_schild(), {"P_t": [1, 0, 0, 0], "J_z": [0, -geo.y, geo.x, 0]}),
        ("FLRW, matter era", kl.chart_flrw_matter(),
         {**{f"P_{c_}": [0] + [int(i == k_) for i in range(3)] for k_, c_ in enumerate("xyz")}, **rot_c}),
        ("de Sitter (conformal time)", kl.chart_de_sitter_conformal(), kl.de_sitter_named()),
    ]
    for name, ch, named in charts:
        sol = kl.solve_polynomial(ch, 2)
        span = kl.same_span(sol["basis"], named, ch.coords)
        out["explicit"][name] = {"n_found": len(sol["basis"]), "n_unknowns": sol["n_unknowns"],
                                 "n_equations": sol["n_equations"], "n_points": sol["n_points"], "span": span,
                                 "named": {k_: [str(c_) for c_ in v_] for k_, v_ in named.items()},
                                 "algebra": kl.structure_constants(named, ch.coords)}
        print(f"  explicit {name}: {len(sol['basis'])} found; same span as named: {span['equal']}", flush=True)

    # Checks by substitution in the usual coordinates
    chk = out["checks"]
    S = geo.schwarzschild_angles()
    fields = {"d_t": [1, 0, 0, 0], **kl.rotations_spherical()}
    chk["schwarzschild_spherical"] = kl.verify(S, fields)
    chk["schwarzschild_norms"] = {k_: str(v_) for k_, v_ in kl.norms(S, {"d_t": [1, 0, 0, 0]}).items()}
    chk["schwarzschild_not_killing"] = {"r d_r": kl.verify(S, {"r d_r": [0, geo.r, 0, 0]})["r d_r"],
                                        "boost-like t d_r + r d_t": kl.verify(S, {"b": [geo.r, geo.t, 0, 0]})["b"]}
    Kb = geo.kerr_angles()
    chk["kerr_bl"] = kl.verify(Kb, {"d_t": [1, 0, 0, 0], "d_phi": [0, 0, 0, 1]})
    Lx = Kb.lie_derivative(kl.rotations_spherical()["J_x"])
    chk["kerr_rotation_J_x_fails"] = {
        "nonzero_components": int(sum(1 for i in range(4) for j in range(i, 4) if sp.simplify(Lx[i, j]) != 0)),
        "all_vanish_when_a_is_0": all(sp.simplify(Lx[i, j].subs(geo.a, 0)) == 0 for i in range(4) for j in range(i, 4)),
        "all_vanish_when_M_is_0": all(sp.simplify(Lx[i, j].subs(geo.M, 0)) == 0 for i in range(4) for j in range(i, 4)),
        "example_component_t_theta": str(sp.factor(sp.simplify(Lx[0, 2])))}
    Kn = kl.norms(Kb, {"d_t": [1, 0, 0, 0]})["d_t"]
    chk["kerr_dt_norm"] = str(sp.simplify(Kn))
    # FLRW with arbitrary a(t): translations and rotations are Killing, d_t is not, a d_t is conformal
    F = geo.flrw()
    flfields = {"P_x": [0, 1, 0, 0], "J_z": [0, -geo.y, geo.x, 0], "J_x": [0, 0, -geo.z, geo.y]}
    chk["flrw_spatial"] = kl.verify(F, flfields)
    a_t = sp.Function("a")(geo.t)
    chk["flrw_d_t_killing"] = kl.verify(F, {"d_t": [1, 0, 0, 0]})["d_t"]
    chk["flrw_conformal_factor_of_a_dt"] = str(F.conformal_factor([a_t, 0, 0, 0]))
    chk["flrw_condition"] = {"matter": str(kl.flrw_condition(geo.t ** R(2, 3))),
                             "radiation": str(kl.flrw_condition(sp.sqrt(geo.t))),
                             "de Sitter": str(kl.flrw_condition(sp.exp(geo.H * geo.t)))}
    # curvature facts
    Sx = geo.schwarzschild()
    chk["schwarzschild_ricci_flat"] = Sx.ricci == sp.zeros(4)
    chk["schwarzschild_kretschmann"] = str(Sx.kretschmann())
    Kx = geo.kerr()
    chk["kerr_ricci_flat"] = all(sp.cancel(sp.together(Kx.ricci[i, j])) == 0 for i in range(4) for j in range(4))
    zv = geo.zipoy_voorhees(2)
    chk["zv_ricci_flat"] = all(sp.cancel(sp.together(zv.ricci[i, j])) == 0 for i in range(4) for j in range(4))
    # Kerr's Killing-Yano and Carter tensors
    Y = kl.kerr_killing_yano()
    chk["kerr_killing_yano"] = kl.killing_yano_check(Kx, Y)
    Kc = kl.carter_tensor(Kx, Y)
    chk["kerr_carter_is_killing_tensor"] = kl.killing_tensor_check(Kx, Kc)
    pts = [{geo.M: 1, geo.a: R(3, 5), geo.r: R(p_, q_), geo.x: R(s_, w_)}
           for p_, q_, s_, w_ in [(7, 2, 1, 3), (5, 1, -2, 7), (9, 4, 3, 5), (11, 3, 1, 9), (13, 2, -1, 4)]]
    chk["kerr_carter_irreducible"] = kl.reducibility(Kx, Kc, {"T": [1, 0, 0, 0], "Phi": [0, 0, 0, 1]}, pts)
    sq = sp.sqrt(1 - geo.x**2)
    rots_x = {"J_x": [0, 0, sq * sp.sin(geo.phi), -geo.x / sq * sp.cos(geo.phi)],
              "J_y": [0, 0, -sq * sp.cos(geo.phi), -geo.x / sq * sp.sin(geo.phi)], "J_z": [0, 0, 0, 1]}
    L2 = sp.zeros(4)
    for f_ in rots_x.values():
        lo = Sx.lower(f_)
        L2 += sp.Matrix(4, 4, lambda i, j: lo[i] * lo[j])
    chk["carter_at_a0_equals_L_squared"] = (Kc.subs(geo.a, 0).applyfunc(sp.cancel) - L2).applyfunc(sp.simplify) == sp.zeros(4)
    out["carter_tensor"] = {f"K_{i}{j}": str(sp.factor(Kc[i, j])) for i in range(4) for j in range(i, 4) if Kc[i, j] != 0}
    print("  checks done", flush=True)

    # Killing tensors (invariant under time translation and the axial rotation)
    for name, met, pt, subs in [("Schwarzschild", geo.schwarzschild(), (0, R(7, 2), R(1, 3), 0), {geo.M: 1}),
                                ("Kerr", geo.kerr(), (0, R(7, 2), R(1, 3), 0), {geo.M: 1, geo.a: R(3, 5)}),
                                ("Zipoy-Voorhees (delta = 2)", geo.zipoy_voorhees(2), (0, 3, R(1, 3), 0), {})]:
        m2 = geo.Metric(met.name, met.coords, met.g.subs(subs)) if subs else met
        t0 = time.time()
        seq = jets.killing_tensor_sequence(m2, pt, 3, 5)
        out["tensor_bounds"][name] = {"sequence": [s_["bound"] for s_ in seq], "N": [s_["N"] for s_ in seq],
                                      "seconds": time.time() - t0}
        print(f"  Killing tensors {name}: {out['tensor_bounds'][name]['sequence']} ({time.time() - t0:.0f} s)", flush=True)
    out["tensor_explicit"] = {
        "Schwarzschild": ["g", "T T", "T Phi", "Phi Phi", "J_x J_x + J_y J_y + J_z J_z (total L^2)"],
        "Kerr": ["g", "T T", "T Phi", "Phi Phi", "Carter's K = Y Y (irreducible)"],
        "Zipoy-Voorhees (delta = 2)": ["g", "T T", "T Phi", "Phi Phi"]}
    return out


# ---------------------------------------------------------------------------
# Numerical geodesics: Noether's theorem at work
# ---------------------------------------------------------------------------

def noether() -> dict:
    import geometry as geo
    import killing as kl
    import orbits as ob

    out = {}
    S = geo.schwarzschild()
    p, e = 12.0, 0.5
    E, L = ob.bound_orbit_EL(p, e)
    r0 = p / (1 - e)
    th0 = np.arccos(0.3)
    w = L / r0**2
    u0 = np.array([E / (1 - 2 / r0), 0.0, -np.sin(th0) * 0.6 * w, np.sqrt(w**2 - (0.6 * w) ** 2) / np.sin(th0)])
    tau, X, U = ob.integrate_geodesic(S, {geo.M: 1}, np.array([0, r0, 0.3, 0]), u0, 3000.0, n_out=500)
    sq = sp.sqrt(1 - geo.x**2)
    fields = {"E (d_t)": [1, 0, 0, 0], "L_z (d_phi)": [0, 0, 0, 1],
              "L_x (J_x)": [0, 0, sq * sp.sin(geo.phi), -geo.x / sq * sp.cos(geo.phi)],
              "L_y (J_y)": [0, 0, -sq * sp.cos(geo.phi), -geo.x / sq * sp.sin(geo.phi)],
              "r d_r (not Killing)": [0, geo.r, 0, 0]}
    c = ob.conserved_along(S, {geo.M: 1}, X, U, fields)
    out["schwarzschild"] = {"tau": tau.tolist(), "r": X[:, 1].tolist(),
                            "values": {k_: v_.tolist() for k_, v_ in c.items()},
                            "rel_spread": {k_: float(np.ptp(v_) / max(np.max(np.abs(v_)), 1e-30)) for k_, v_ in c.items()}}
    K_ = geo.kerr()
    sub = {geo.M: 1, geo.a: R(9, 10)}
    Kc = kl.carter_tensor(K_, kl.kerr_killing_yano())
    # Schwarzschild's conserved "total angular momentum squared", J_x^2 + J_y^2 + J_z^2, built
    # with the same vector fields but the Kerr metric: NOT conserved in Kerr (the J_x, J_y
    # rotations are not symmetries there). Its wobble shows that Carter's K is something new.
    L2 = sp.zeros(4)
    for f_ in (fields["L_x (J_x)"], fields["L_y (J_y)"], fields["L_z (d_phi)"]):
        lo = K_.lower(f_)
        L2 += sp.Matrix(4, 4, lambda i, j: lo[i] * lo[j])
    x0, u0, consts = ob.kerr_initial_state(0.9, 4.0, 9.0, 40.0)
    tau, X, U = ob.integrate_geodesic(K_, sub, x0, u0, 3000.0, n_out=500)
    c = ob.conserved_along(K_, sub, X, U, {"E (d_t)": [1, 0, 0, 0], "L_z (d_phi)": [0, 0, 0, 1]},
                           {"Carter K (Killing tensor)": Kc, "total L^2 (not a symmetry of Kerr)": L2})
    out["kerr"] = {"tau": tau.tolist(), "r": X[:, 1].tolist(), "x": X[:, 2].tolist(), "phi": X[:, 3].tolist(),
                   "constants": consts, "orbit": {"a": 0.9, "r_p": 4.0, "r_a": 9.0, "inclination_deg": 40.0},
                   "values": {k_: v_.tolist() for k_, v_ in c.items()},
                   "rel_spread": {k_: float(np.ptp(v_) / max(np.max(np.abs(v_)), 1e-30)) for k_, v_ in c.items()}}
    return out


# ---------------------------------------------------------------------------
# Orbits, S2, shadows, thermodynamics, cosmology, chaos
# ---------------------------------------------------------------------------

def orbits_section() -> dict:
    import orbits as ob
    spins = np.concatenate([np.linspace(0, 0.9, 19), [0.92, 0.94, 0.96, 0.98, 0.99, 0.995, 0.998, 0.999, 0.9999, 1.0]])
    tab = ob.special_radii_table(spins)
    sch = ob.schwarzschild_special_radii()
    prec = []
    for p, e in [(20.0, 0.3), (100.0, 0.6), (5541.0, 0.8859)]:
        E, L = ob.bound_orbit_EL(p, e)
        prec.append({"p": p, "e": e, "exact": float(ob.precession_exact(p, e)),
                     "quadrature": float(ob.precession_by_quadrature(E, L)), "six_pi_over_p": 6 * np.pi / p})
    isco_checks = []
    for s in (0.3, 0.7, 0.95):
        d = ob.isco_kerr_derived(s)
        isco_checks.append({"a": s, "E2_vs_1-2/3r": d["E"] ** 2 - (1 - 2 / (3 * d["r_isco"])),
                            "L-aE_vs_r/sqrt3": (d["L"] - s * d["E"]) - d["r_isco"] / np.sqrt(3)})
    return {"table": tab, "schwarzschild": sch, "precession": prec, "isco_relations": isco_checks}


def s2_section() -> dict:
    import s2
    post = s2.posterior_EL()
    keys = ["a_au", "p", "r_peri_au", "r_apo_au", "E_newton_kms2", "L_newton_aukms", "binding_newton", "binding_gr",
            "L_newton_geom", "L_gr_geom", "E_gr", "precession_exact_arcmin", "precession_6pi_arcmin", "h_model_mas2yr"]
    summ = {k_: s2.summarise(post[k_]) for k_ in keys}
    summ["binding_gr_over_newton_minus_1"] = s2.summarise(post["binding_gr"] / post["binding_newton"] - 1)
    summ["L_gr_over_newton_minus_1"] = s2.summarise(post["L_gr_geom"] / post["L_newton_geom"] - 1)
    thin = slice(None, None, 20)
    samples = {k_: post[k_][thin].tolist() for k_ in ("E_newton_kms2", "L_newton_aukms", "binding_gr", "L_gr_geom",
                                                     "binding_newton", "L_newton_geom", "precession_exact_arcmin")}
    seasons = s2.seasonal_tests(n_draws=400)
    geod = s2.newton_vs_gr_along_orbit()
    geod_small = {k_: (v_[::10].tolist() if isinstance(v_, np.ndarray) else v_) for k_, v_ in geod.items()}
    return {"summary": summ, "samples": samples, "seasons": seasons, "geodesic": geod_small}


def shadows_section() -> dict:
    import shadows as sh
    res = sh.analyse()
    curves = {}
    for tgt, inc in (("sgra", res["sgra"]["inclination_assumed_deg"]), ("m87", res["m87"]["inclination_assumed_deg"])):
        curves[tgt] = {}
        for s in (0.0, 0.5, 0.94, 0.998):
            al, be = sh.critical_curve(s, inc, n=600)
            curves[tgt][str(s)] = {"alpha": al.tolist(), "beta": be.tolist()}
    edge = {}
    for s in (0.5, 0.94, 0.998):
        al, be = sh.critical_curve(s, 90.0, n=600)
        edge[str(s)] = {"alpha": al.tolist(), "beta": be.tolist()}
    res["curves"] = curves
    res["curves_edge_on"] = edge
    return res


def thermo_section() -> dict:
    import thermo as th
    cen = th.census()
    rows = {r_["name"]: r_ for r_ in cen["rows"]}
    ff = th.far_future({"1 solar mass": 1.0, "GW250114 remnant": rows["GW250114 remnant"]["M_msun"],
                        "Sgr A*": rows["Sgr A*"]["M_msun"], "M87*": rows["M87*"]["M_msun"]})
    one = {k_: float(v_) for k_, v_ in th.kerr_thermo(1.0, 0.0).items()}
    life = {"page_photons_gravitons_yr": float(th.evaporation_time_yr(1.0)),
            "page_with_neutrinos_yr": float(th.evaporation_time_yr(1.0, 2.011e-4)),
            "blackbody_photons_yr": float(th.blackbody_evaporation_time_yr(1.0))}
    return {"census": cen, "far_future": ff, "one_solar_mass": one, "lifetimes_1msun": life,
            "mass_T_equals_CMB_kg": float(th.mass_at_temperature_msun(C.T_CMB) * C.M_SUN)}


def cosmology_section() -> dict:
    import cosmology as cs
    return {"tcmb": cs.tcmb_fit(), "geodesics": cs.geodesics_demo()}


def chaos_section(force=False) -> dict:
    import chaos
    return chaos.run(force=force)


SECTIONS = {"symbolic": symbolic, "noether": noether, "orbits": orbits_section, "s2": s2_section,
            "shadows": shadows_section, "thermo": thermo_section, "cosmology": cosmology_section}


def run_all(force: bool = False, only=None) -> dict:
    out = {}
    for name, fn in SECTIONS.items():
        if only and name not in only:
            continue
        t0 = time.time()
        print(f"== {name}", flush=True)
        out[name] = cached(name, fn, force=force)
        print(f"   {name}: {time.time() - t0:.0f} s", flush=True)
    if not only or "chaos" in only:
        print("== chaos", flush=True)
        t0 = time.time()
        out["chaos"] = chaos_section(force=force)
        print(f"   chaos: {time.time() - t0:.0f} s", flush=True)
    return out
