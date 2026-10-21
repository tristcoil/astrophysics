"""Regenerate every result, figure and the orbit explorer of project N6.

    python fetch_data.py   # once: 27 MB (GW250114 samples) + fingerprints of earlier projects' files
    python run.py          # ~15-25 min the first time (symbolic algebra, orbits), seconds afterwards
    python run.py --force  # recompute everything

Writes data/results.json, figures/ and figures/kerr_orbit_explorer.html, and
prints each headline number next to the value it is checked against.
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np

import config as C


def _clean(x):
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items() if not str(k).startswith("_")}
    if isinstance(x, (list, tuple)):
        if len(x) > 60 and all(isinstance(v, (int, float, np.floating)) for v in x[:5]):
            return f"[{len(x)} values; see data/cache]"
        return [_clean(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    if isinstance(x, np.ndarray):
        return f"[array {x.shape}; see data/cache]"
    return x


def headline(res) -> list:
    sym, orb, s2r, shd, th, cs, ch = (res[k] for k in ("symbolic", "orbits", "s2", "shadows", "thermo", "cosmology",
                                                         "chaos"))
    L = []
    L.append("== Killing vectors: jet upper bound -> explicit solutions found (textbook)")
    for name, d in sym["vector_bounds"].items():
        ex = sym["explicit"].get(name, {}).get("n_found", "-")
        L.append(f"   {name:36s} bound {d['final']:2d}  explicit {ex!s:>2}  ({d['known'] or C.KILLING_COUNTS.get(name)})")
    L.append(f"   de Sitter (conformal coordinates): explicit {sym['explicit']['de Sitter (conformal time)']['n_found']} (10)")
    L.append("== Killing tensors invariant under t and phi: " +
             ", ".join(f"{k} {v['sequence'][-1]}" for k, v in sym["tensor_bounds"].items()) + "  (5, 5, 4)")
    ck = sym["checks"]
    L.append(f"   Kerr Killing-Yano {ck['kerr_killing_yano']}, Carter tensor Killing {ck['kerr_carter_is_killing_tensor']}, "
             f"irreducible {ck['kerr_carter_irreducible']['irreducible']}, = L^2 at a=0 {ck['carter_at_a0_equals_L_squared']}")
    sc = orb["schwarzschild"]
    L.append(f"== Schwarzschild: ISCO {sc['r_isco']:.6f} (6), photon sphere {sc['r_photon']:.6f} (3), "
             f"b {sc['b_crit']:.6f} (5.196152), efficiency {100 * sc['efficiency']:.3f}% (5.719%)")
    tab = orb["table"]
    i998 = tab["a"].index(0.998)
    L.append(f"== Kerr: efficiency a=0.998 {100 * tab['eff_pro'][i998]:.2f}% (32.1%), a=1 {100 * tab['eff_pro'][-1]:.2f}% "
             f"(42.26%), retrograde a=1 {100 * tab['eff_retro'][-1]:.3f}% (3.775%); max |derived - BPT| "
             f"{max(tab['max_dev_from_bpt']):.1e}")
    sm = s2r["summary"]
    L.append(f"== S2: E = 1 - {sm['binding_gr']['median']:.4e} c^2, L = {sm['L_gr_geom']['median']:.2f} GM/c; "
             f"precession {sm['precession_exact_arcmin']['median']:.2f}' +- {sm['precession_exact_arcmin']['sd']:.2f}' "
             "(H1 11.7'; GRAVITY 12.1' x f_SP)")
    hc, ec = s2r["seasons"]["h_consistency"], s2r["seasons"]["E_consistency"]
    L.append(f"   Kepler II: {hc['n']} seasons, chi2 {hc['chi2_vs_model']:.1f}, mean {hc['weighted_mean']:.0f} +- "
             f"{hc['weighted_mean_err']:.0f} vs fit {hc['model_median']:.0f} mas^2/yr; energy: {ec['n']} seasons, "
             f"chi2 {ec['chi2_vs_model']:.1f}")
    sg = shd["sgra"]
    L.append(f"== Shadows: Kerr mean diameter {shd['grid']['min']:.2f}-{shd['grid']['max']:.2f} GM/c^2 (9.6-10.4); "
             f"Sgr A* {sg['schwarzschild_shadow_uas'][0]:.1f} +- {sg['schwarzschild_shadow_uas'][1]:.1f} uas vs ring "
             f"{sg['ring_uas'][0]:.1f} +- {sg['ring_uas'][1]:.1f} (ring/theta_g {sg['ring_over_theta_g'][0]:.2f})")
    one = th["one_solar_mass"]
    L.append(f"== Hawking: 1 Msun T = {one['T_K']:.4e} K (6.1703e-8), S = {one['S_kB']:.4e} k_B (1.0489e77), "
             f"lifetime {th['lifetimes_1msun']['page_with_neutrinos_yr']:.3e} yr (Page: 2.16e66)")
    gw = th["census"]["gw250114"]
    L.append(f"   GW250114: area {gw['area_initial_km2']['median']:.0f} -> {gw['area_final_km2']['median']:.0f} km^2, "
             f"Delta S = {gw['delta_S_kB']['median']:.2e} k_B; remnant T = {gw['T_final_K']['median']:.3e} K")
    tc = cs["tcmb"]
    L.append(f"== CMB: beta = {tc['beta']:.4f} +{tc['plus']:.4f}/-{tc['minus']:.4f} (D7: 0.0040 +0.0108/-0.0107); "
             f"photon a*E spread {cs['geodesics']['photon']['aE_spread']:.1e}")
    L.append("== Chaos test: separation growth of orbits 1e-8 apart: " +
             "; ".join(f"{n}: " + ", ".join(f"{d['final_over_initial']:.3g}" for d in v) for n, v in ch["divergence"].items()))
    return L


def main(force: bool = False) -> dict:
    t0 = time.time()
    import fetch_data
    changed = fetch_data.check_inputs()
    import analysis
    res = analysis.run_all(force=force)
    print("== figures")
    import figures
    for p in figures.make_all():
        if p.endswith(".png"):
            print("  ", p.split("/")[-1])
    import viewer
    print("== explorer:", viewer.make())
    lines = headline(res)
    print("\n".join(lines))
    out = {k: _clean(v) for k, v in res.items()}
    out["inputs_changed_since_fetch"] = changed
    out["headline"] = lines
    out["runtime_s"] = time.time() - t0
    C.RESULTS_PATH.write_text(json.dumps(out, indent=1, default=str))
    print(f"done in {time.time() - t0:.0f} s -> {C.RESULTS_PATH}")
    return out


if __name__ == "__main__":
    main(force="--force" in sys.argv)
