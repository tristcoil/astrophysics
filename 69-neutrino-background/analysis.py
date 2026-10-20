"""Every measurement in the project, in order, with the expensive steps cached.

    python analysis.py           # use cached steps where present
    python analysis.py --force   # recompute everything (CAMB-heavy steps: ~1 h)

Steps:

1. **Relic neutrinos today**: densities, temperature, speeds, the PTOLEMY rate.
2. **Decoupling**: our own N_eff from the temperature evolution through
   neutrino decoupling and e± annihilation.
3. **Masses**: what oscillations require (NuFIT 6.0), in Σm_ν, m_β and m_ββ.
4. **The CMB and BAO fits** (CAMB): ΛCDM, then N_eff free, then Σm_ν free,
   for a ladder of data combinations, with profile likelihoods.
5. **DESI's own chains**: the published limits recomputed from DESI's MCMC.
6. **Illustrations**: spectra at fixed acoustic scale, matter power.

Results go to ``data/results.json``.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import config as cfg  # noqa: E402


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None  # strict JSON: no NaN/Infinity
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return obj


def _cached_json(name: str, fn, force: bool):
    path = cfg.DERIVED_DIR / name
    if path.exists() and not force:
        return json.loads(path.read_text())
    out = _jsonable(fn())
    path.write_text(json.dumps(out, indent=1))
    return out


# ---------------------------------------------------------------------------
# 1. Relic neutrinos today
# ---------------------------------------------------------------------------


def step_relic(neff_ours: float, heating=None) -> dict:
    import masses as ms
    import relic

    o = cfg.NUFIT["NO"]
    m1, m2, m3 = ms.masses(0.0, o["dm21"][0], o["dm3l"][0], "NO")
    states = {"m1 (lightest, if massless)": 0.0, "m2 (normal ordering, m1 = 0)": float(m2),
              "m3 (normal ordering, m1 = 0)": float(m3), "0.1 eV (illustration)": 0.1}
    out = relic.relic_summary(cfg.T_CMB_K, neff_ours, states, cfg.OMEGA_M_REF,
                              cfg.TRITIUM_SIGMA_V_CM2, heating=heating)
    out["ratio_n_nu_to_n_gamma"] = out["n_nu_total_cm3"] / out["n_gamma_cm3"]
    return out


# ---------------------------------------------------------------------------
# 2. Decoupling
# ---------------------------------------------------------------------------


def step_decoupling() -> dict:
    import decoupling as dc

    plasma = dc.Plasma()
    table = dc.load_rate_table()
    runs = {}
    for label, rates, qed in cfg.DECOUPLING_RUNS:
        r = dc.solve(rates=rates, qed_order=qed, plasma=plasma, table=table)
        row = {k: v for k, v in r.items() if not isinstance(v, np.ndarray)}
        if rates != "off":
            row["T_dec_e_MeV"] = dc.decoupling_temperature(r, "e")
            row["T_dec_mu_MeV"] = dc.decoupling_temperature(r, "mu")
        row["escudero_2020"] = cfg.ESCUDERO_2020_TABLE1.get(label)
        step = max(1, len(r["T"]) // 300)
        row["track"] = {k: r[k][::step].tolist() for k in
                        ("lna", "T", "Tnue", "Tnumu", "neff_track", "gamma_e_over_H",
                         "gamma_mu_over_H")}
        runs[label] = row
    head = runs[cfg.DECOUPLING_HEADLINE]
    out = {"runs": runs, "headline": cfg.DECOUPLING_HEADLINE,
           "check_A12": {k: (head[k], v) for k, v in cfg.ESCUDERO_2020_A12.items()}}
    try:  # project P2's independent implementation of the same calculation
        p2 = json.loads(cfg.P2_RESULTS.read_text())["neutrinos"]["cases"]
        # same set-up: Fermi-Dirac rates with the electron mass, QED to O(e^3),
        # separate nu_e and nu_mu,tau temperatures (P2 starts no-interaction runs
        # at 10 MeV rather than our 300 MeV, so only this case is like for like)
        out["p2_crosscheck"] = {"ours": head["neff"], "p2": p2["fdme_qed3_separate"]["N_eff"]}
    except (OSError, KeyError, json.JSONDecodeError):
        out["p2_crosscheck"] = None
    return out


# ---------------------------------------------------------------------------
# 3. Masses
# ---------------------------------------------------------------------------


def nufit_ranges(profiles: dict) -> dict:
    """Best fit and 1σ / 3σ ranges from NuFIT's 1D Δχ² profiles, per ordering.

    NuFIT's files measure Δχ² from the global minimum over both orderings, so
    each ordering's own minimum is subtracted first. Ranges are where the
    profile lies below 1 and 9 (interpolated on NuFIT's grid).
    """
    conv = {"T12": ("s12sq", lambda x: x), "T13": ("s13sq", lambda x: x),
            "DMS": ("dm21", lambda x: 10.0**x), "DMA": ("dm3l", lambda x: x * 1e-3)}
    out = {}
    for order, prof in profiles.items():
        y_any = prof["DMA"][1]
        offset = float(np.min(y_any))
        row = {"dchi2_vs_global_min": offset}
        for tag, (name, f) in conv.items():
            x, y = prof[tag]
            y = y - np.min(y)
            fine = np.linspace(x.min(), x.max(), 200001)
            yf = np.interp(fine, x, y)
            # the grid quantises the minimum: take the vertex of a parabola
            # through the grid points around it (Δχ² < 1)
            near = np.argsort(np.abs(x - x[np.argmin(y)]))[:7]
            near = near[y[near] < 1.5] if np.sum(y[near] < 1.5) >= 3 else near[:3]
            a2, a1, _ = np.polyfit(x[near], y[near], 2)
            best = float(f(-a1 / (2 * a2))) if a2 > 0 else float(f(fine[np.argmin(yf)]))
            r1 = fine[yf <= 1.0]
            r3 = fine[yf <= 9.0]
            row[name] = {"best": best, "lo1": float(f(r1.min())), "hi1": float(f(r1.max())),
                         "lo3": float(f(r3.min())), "hi3": float(f(r3.max()))}
        out[order] = row
    return out


def step_masses() -> dict:
    import fetch_data as fd
    import masses as ms

    out = {"nufit": {}}
    for analysis in ("TByes", "TBoff"):
        out["nufit"][analysis] = nufit_ranges(fd.load_nufit(analysis))
    nf = out["nufit"]["TByes"]
    out["nufit_check"] = {o: {k: {"ours": nf[o][k]["best"], "published": cfg.NUFIT[o][k][0]}
                              for k in ("s12sq", "s13sq", "dm21", "dm3l")} for o in ("NO", "IO")}
    for order in ("NO", "IO"):
        b = {k: nf[order][k]["best"] for k in ("dm21", "dm3l", "s12sq", "s13sq")}
        smin = ms.minimum_sum(b["dm21"], b["dm3l"], order)
        corners = [ms.minimum_sum(d21, d3, order)
                   for d21 in (nf[order]["dm21"]["lo3"], nf[order]["dm21"]["hi3"])
                   for d3 in (nf[order]["dm3l"]["lo3"], nf[order]["dm3l"]["hi3"])]
        m = ms.masses(0.0, b["dm21"], b["dm3l"], order)
        mb0 = float(ms.m_beta(0.0, b["dm21"], b["dm3l"], order, b["s12sq"], b["s13sq"]))
        mbb = ms.m_betabeta_range(0.0, b["dm21"], b["dm3l"], order, b["s12sq"], b["s13sq"])
        grid = np.geomspace(1e-5, 3.0, 40001)
        mbeta = ms.m_beta(grid, b["dm21"], b["dm3l"], order, b["s12sq"], b["s13sq"])
        ml_katrin = float(np.interp(cfg.KATRIN_MBETA_LIMIT, mbeta, grid))
        out[order] = {"min_sum_eV": smin, "min_sum_3sigma_eV": [min(corners), max(corners)],
                      "masses_at_min_eV": [float(x) for x in m], "m_beta_at_min_eV": mb0,
                      "m_betabeta_at_min_eV": [float(mbb[0]), float(mbb[1])],
                      "katrin_lightest_limit_eV": ml_katrin,
                      "katrin_sum_limit_eV": float(ms.sum_masses(ml_katrin, b["dm21"], b["dm3l"],
                                                                 order))}
    return out


# ---------------------------------------------------------------------------
# 6. Illustrations
# ---------------------------------------------------------------------------


def step_illustrations(force: bool) -> dict:
    import illustrations as il

    s = il.neff_spectra(force=force)
    m = il.matter_power(force=force)
    out = {"neff_spectra": {}, "matter_power": {}}
    ell = s["ell"]
    band = (ell >= 1500) & (ell <= 2500)
    for tag in ("theta", "zeq", "thetad"):
        row = {}
        for i, n in enumerate(s["neffs"]):
            row[f"{n:.3f}"] = {q: float(s[f"{tag}_{q}_{i}"]) for q in
                               ("H0", "YHe", "zeq", "thetad", "thetastar", "rstar")}
            ref = s[f"{tag}_tt_1"]
            row[f"{n:.3f}"]["tt_ratio_1500_2500"] = float(np.mean(s[f"{tag}_tt_{i}"][band]
                                                                 / ref[band]))
        out["neff_spectra"][tag] = row
    kh = m["kh"]
    for fam in ("fixed_h", "fixed_theta"):
        rows = []
        for i, mass in enumerate(m["masses"]):
            j = int(np.argmin(np.abs(kh - 1.0)))
            fnu = float(m[f"{fam}_fnu_{i}"])
            lin = float(m[f"{fam}_lin_{i}"][j] / m[f"{fam}_lin_0"][j] - 1)
            nl = float(m[f"{fam}_nl_{i}"][j] / m[f"{fam}_nl_0"][j] - 1)
            rows.append({"mnu": float(mass), "f_nu": fnu, "dP_lin_k1": lin, "dP_nl_k1": nl,
                         "lin_over_fnu": lin / fnu if fnu else None,
                         "nl_over_fnu": nl / fnu if fnu else None,
                         "h": float(m[f"{fam}_h_{i}"]), "sigma8": float(m[f"{fam}_sigma8_{i}"])})
        out["matter_power"][fam] = rows
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def run(force: bool = False, *, skip_fits: bool = False) -> dict:
    t0 = time.monotonic()
    results = {}
    print("1-2. Decoupling: our own N_eff")
    results["decoupling"] = _cached_json("decoupling.json", step_decoupling, force)
    head = results["decoupling"]["runs"][cfg.DECOUPLING_HEADLINE]
    neff_ours = head["neff"]
    print(f"     N_eff = {neff_ours:.4f}")
    print("1. Relic neutrinos today")
    results["relic"] = _jsonable(step_relic(
        neff_ours, (head["ratio_vs_instantaneous_e"], head["ratio_vs_instantaneous_mu"])))
    print("3. Masses from oscillations")
    results["masses"] = _jsonable(step_masses())
    if not skip_fits:
        import cmb_fits

        print("4. CMB and BAO fits (CAMB)")
        results["fits"] = cmb_fits.run(force=force)
        print("5. DESI's chains")
        results["desi"] = cmb_fits.desi_answer_key()
    print("6. Illustrations")
    results["illustrations"] = _cached_json("illustrations.json",
                                            lambda: step_illustrations(force), force)
    results["seconds"] = time.monotonic() - t0
    cfg.RESULTS_PATH.write_text(json.dumps(_jsonable(results), indent=1))
    print(f"   results -> {cfg.RESULTS_PATH.relative_to(cfg.PROJECT_DIR)} "
          f"({results['seconds']:.0f} s)")
    return results


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--skip-fits", action="store_true", help="skip the CAMB fits")
    args = ap.parse_args(argv)
    run(force=args.force, skip_fits=args.skip_fits)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
