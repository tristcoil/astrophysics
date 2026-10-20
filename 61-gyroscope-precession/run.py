"""Regenerate every result and figure of project N4 from the cached data.

    python fetch_data.py     # once: ~0.9 GB of weekly orbits (8 analysis centres) + small files
    python campaign.py       # once: ~9,000 weekly orbit fits, ~4-5 h on 8 cores (cached, resumable)
    python run.py            # minutes: the frame-dragging fits, GP-B, pulsars, figures
    python run.py --force    # recompute the cached fits and predictions too

Writes data/results.json and figures/.
"""

from __future__ import annotations

import json
import sys
import time

import numpy as np

import config as C


def main(force: bool = False) -> dict:
    t0 = time.time()
    import campaign
    import framedrag as FD
    import nodes

    print("== frame dragging: every analysis centre")
    acs = [C.PRIMARY_AC] + [a for a in C.ACS if a != C.PRIMARY_AC]
    per_ac = {}
    for ac in acs:
        if not (nodes.ARC_DIR / f"lares2{C.BASE_TAG}" / ac).exists():
            print(f"   {ac}: no fits cached (run campaign.py)")
            continue
        variants = FD.VARIANTS if ac == C.PRIMARY_AC else {"baseline": FD.VARIANTS["baseline"]}
        per_ac[ac] = FD.analyse(ac, C.BASE_TAG, variants=variants, force=force)

    print("== systematics on ASI's orbits")
    syst = {}
    for name, tag in (("injection (model with frame dragging)", "_c20_cr_mu1"),
                      ("C_R fitted every week", "_c20"),
                      ("C_R fitted every week, injection", "_c20_mu1"),
                      ("C_R at the nominal ILRS values", "_c20_crfix"), ("CSR's C20 for J2", "_csr_cr"),
                      ("GGM05S gravity field", "_c20_cr_ggm05s"), ("no ocean tides", "_c20_cr_noocean"),
                      ("no de Sitter in the model", "_c20_cr_nods")):
        if (nodes.ARC_DIR / f"lares2{tag}" / C.PRIMARY_AC).exists():
            r = FD.analyse(C.PRIMARY_AC, tag, variants={"baseline": FD.VARIANTS["baseline"]}, force=force)
            syst[name] = {"tag": tag, **r["variants"]["baseline"], "raw_sum_rate": r["raw_sum_rate"]}
    spans = FD.span_test(tag=C.BASE_TAG)
    incl = FD.inclination_fit(tag=C.BASE_TAG)

    print("== J2 from the pair")
    have = lambda tag: (nodes.ARC_DIR / f"lares2{tag}" / C.PRIMARY_AC).exists()
    j2_tag = "_cr" if have("_cr") else ""                 # pass 1 with C_R fixed, if run
    j2 = FD.j2_comparison(tag=j2_tag)
    j2["tag"] = j2_tag
    om1_tag = "_om1_cr" if have("_om1_cr") else ("_om1" if have("_om1") else None)
    j2_om1 = FD.j2_comparison(tag=om1_tag) if om1_tag else None

    print("== Gravity Probe B")
    import gpb
    gp = gpb.run(force=force)

    print("== pulsars and the ladder")
    import pulsars
    base = per_ac[C.PRIMARY_AC]["variants"]["baseline"]
    summary = summarise(per_ac, syst, spans)

    print("== figures")
    import figures
    import tides
    figures.gpb_diagram()
    figures.gpb_result(gp)
    we = {s: nodes.weekly_elements(s, C.PRIMARY_AC, C.BASE_TAG) for s in ("lageos1", "lares2")}
    rates = node_rates(we)
    figures.node_drift(we, rates)
    figures.frame_dragging(FD.curves(C.PRIMARY_AC, C.BASE_TAG, "baseline"), summary=summary)
    figures.j2_figure(j2, nodal_terms(j2, j2_om1))
    figures.ladder(pulsars.ladder(slr_result={"mu": summary["mu"], "err": summary["total_err"],
                                              "rate_mas_yr": summary["lt_sum"]}))
    figures.systematics(summary, systematics_rows(per_ac, syst, spans, summary))

    results = {"frame_dragging": summary, "per_ac": {k: {kk: vv for kk, vv in v.items() if kk != "variants"}
                                                     | {"baseline": v["variants"]["baseline"]}
                                                     for k, v in per_ac.items()},
               "asi_variants": per_ac[C.PRIMARY_AC]["variants"], "systematics": syst, "spans": spans,
               "inclination_rates_mas_yr": incl, "node_drift": rates,
               "j2": {k: v for k, v in j2.items() if np.isscalar(v)},
               "j2_om1": ({k: v for k, v in j2_om1.items() if np.isscalar(v)} | {"tag": om1_tag})
               if j2_om1 else None, "gpb": gp,
               "pulsar_predictions": pulsars.predictions(), "runtime_s": time.time() - t0}
    C.RESULTS_PATH.write_text(json.dumps(results, indent=1, default=float))
    print(f"done in {time.time() - t0:.0f} s -> {C.RESULTS_PATH}")
    return results


def summarise(per_ac: dict, syst: dict, spans: dict) -> dict:
    """The headline number: ASI's baseline, statistical error, and the systematic spread."""
    base = per_ac[C.PRIMARY_AC]["variants"]["baseline"]
    mus = np.array([v["variants"]["baseline"]["mu"] for v in per_ac.values()])
    errs = np.array([v["variants"]["baseline"]["mu_err"] for v in per_ac.values()])
    # Centres whose weekly orbits carry large week-to-week rotations about the Earth's axis
    # (NSGF: common node jumps of 10-18 mas in both satellites) are consistent but imprecise;
    # they are kept in the table and left out of the spread and the combined mean.
    good = errs < 0.1
    good_acs = [ac for ac, g in zip(per_ac, good) if g]
    asi_v = per_ac[C.PRIMARY_AC]["variants"]
    nuis = np.array([asi_v[k]["mu"] for k in ("K1+S2", "tides", "tides+srp", "baseline") if k in asi_v])
    model = [s["mu"] for k, s in syst.items()
             if k in ("C_R at the nominal ILRS values", "CSR's C20 for J2", "GGM05S gravity field")]
    parts = {"analysis_centres": float(np.std(mus[good], ddof=1)) if good.sum() > 1 else 0.0,
             "nuisance_terms": float(np.std(nuis, ddof=1)) if nuis.size > 1 else 0.0,
             "force_model": float(np.max(np.abs(np.array(model) - base["mu"]))) if model else 0.0}
    syst_err = float(np.sqrt(sum(v ** 2 for v in parts.values())))
    lt = sum(per_ac[C.PRIMARY_AC]["lt_rate"].values())
    import framedrag as FD
    ser = FD.series(C.PRIMARY_AC, C.BASE_TAG)
    ds = FD.de_sitter_node_rate(ser)
    comb = None
    if len(good_acs) > 2:
        sc = FD.series_combined(good_acs, C.BASE_TAG)
        f = FD.fit(sc, FD.NUISANCE_SETS["full"], noise="ar1")
        comb = {"mu": f["mu"], "mu_err": f["mu_err"], "n_weeks": int(sc["mjd"].size), "noise": f["noise"]}
    return {"mu": base["mu"], "stat_err": base["mu_err"], "syst_err": syst_err, "syst_parts": parts,
            "de_sitter_node_rate": ds, "ds_over_lt": ds["sum"] / lt,
            **({"mean of centres": comb} if comb else {}),
            "total_err": float(np.hypot(base["mu_err"], syst_err)), "lt_sum": lt,
            "rate_mas_yr": base["mu"] * lt, "ac_mus": dict(zip(per_ac, mus.tolist())),
            "ac_errs": dict(zip(per_ac, errs.tolist())),
            "ac_weighted_mean": float(np.sum(mus / errs ** 2) / np.sum(1 / errs ** 2)),
            "ac_good": good_acs, "ac_excluded": [ac for ac in per_ac if ac not in good_acs],
            "ac_range_good": [float(mus[good].min()), float(mus[good].max())],
            "n_weeks": per_ac[C.PRIMARY_AC]["n_weeks"],
            "raw_sum_rate": per_ac[C.PRIMARY_AC]["raw_sum_rate"], "spans": spans}


def systematics_rows(per_ac: dict, syst: dict, spans: dict, summary: dict) -> list:
    """Rows of the forest plot (figures.systematics)."""
    import framedrag as FD
    rows = [{"header": True, "label": "Analysis choices (ASI's orbits)"}]
    v = per_ac[C.PRIMARY_AC]["variants"]
    names = {"baseline": "adopted: K1, K2, S2, P1, S1, annual terms", "none": "no periodic terms",
             "K1": "K1 only", "K1+S2": "K1 + S2", "tides": "K1, K2, S2", "tides+srp": "K1, K2, S2, S1",
             "baseline, white jumps": "adopted terms, white jump noise only",
             "baseline, per satellite": "adopted terms, one noise model per satellite"}
    for k, lab in names.items():
        if k in v:
            rows.append({"label": lab, "mu": v[k]["mu"], "err": v[k]["mu_err"], "primary": k == "baseline"})
    rows.append({"header": True, "label": "Data span"})
    for k, lab in (("first", "first 1,050 days"), ("last", "last 1,050 days")):
        rows.append({"label": lab, "mu": spans[k]["mu"], "err": spans[k]["mu_err"]})
    rows.append({"header": True, "label": "Analysis centre (whose weekly orbits)"})
    import figures
    for ac, r in per_ac.items():
        b = r["variants"]["baseline"]
        row = {"label": figures.AC_NAMES.get(ac, ac), "mu": b["mu"], "err": b["mu_err"],
               "primary": ac == C.PRIMARY_AC}
        if ac in summary.get("ac_excluded", []):
            row.update(test=True, note=f"{b['mu']:.2f} ± {b['mu_err']:.2f}: weekly frame rotations")
        rows.append(row)
    if "mean of centres" in summary:
        m = summary["mean of centres"]
        rows.append({"label": "weekly mean of the centres (without NSGF)", "mu": m["mu"], "err": m["mu_err"]})
    rows.append({"header": True, "label": "Force model (ASI's orbits)"})
    labels = {"C_R at the nominal ILRS values": "C_R at the nominal ILRS values (1.13, 1.08)",
              "C_R fitted every week": "C_R fitted every week (absorbs ~3%)",
              "CSR's C20 for J2": "CSR's monthly C20 for J2", "GGM05S gravity field": "GGM05S gravity field",
              "no ocean tides": "no ocean tides"}
    for k, lab in labels.items():
        if k in syst:
            rows.append({"label": lab, "mu": syst[k]["mu"], "err": syst[k]["mu_err"]})
    rows.append({"header": True, "label": "Tests (shifted by the expected amount)"})
    if "injection (model with frame dragging)" in syst:
        r = syst["injection (model with frame dragging)"]
        rows.append({"label": "frame dragging in the model: μ + 1", "mu": r["mu"] + 1, "err": r["mu_err"],
                     "test": True, "note": f"fit {r['mu']:+.3f}"})
    if "C_R fitted every week, injection" in syst:
        r = syst["C_R fitted every week, injection"]
        rows.append({"label": "the same with C_R fitted weekly: μ + 1", "mu": r["mu"] + 1, "err": r["mu_err"],
                     "test": True, "note": f"fit {r['mu']:+.3f}"})
    if "no de Sitter in the model" in syst:
        r = syst["no de Sitter in the model"]
        ds = summary.get("ds_over_lt", 0.574)
        rows.append({"label": f"no de Sitter in the model: μ − {ds:.3f}", "mu": r["mu"] - ds, "err": r["mu_err"],
                     "test": True, "note": f"fit {r['mu']:.3f}"})
    return rows


def node_rates(we: dict) -> dict:
    import frames as F
    rates = {}
    for s, w in we.items():
        yr = F.mjd_to_year(w["mjd"])
        rates[s] = float(np.polyfit(yr, np.degrees(np.unwrap(w["node_cip"])), 1)[0])
    rates["lt_sum"] = 61.35
    rates["ds_sum"] = float(2 * 19.2 * np.cos(np.radians(23.44)))
    rates["inc_sum"] = float(np.degrees(we["lageos1"]["inc_cip"].mean() + we["lares2"]["inc_cip"].mean()))
    rates["da_km"] = float((we["lageos1"]["a"].mean() - we["lares2"]["a"].mean()) / 1e3)
    return rates


def nodal_terms(j: dict, j_om1: dict | None = None) -> dict:
    """The 18.6-year tide terms for the J2 figure and the phase-fix check.

    With j_om1 (a full re-run with the ocean term in phase) its pair series and its offset
    from CSR are used; otherwise the analytic correction 2 x 6.58e-11 cos N'.
    """
    import tides
    to_mjd = lambda y: (np.asarray(y) - 2000) * 365.25 + 51_544.5
    mjd = np.arange(56_000.0, 64_500.0, 5.0)
    beta = tides.doodson_arguments(mjd, mjd)
    # amplitudes of the 18.6-yr terms, fitted to our own tide model over 19 years (test_physics)
    solid_amp, ocean_amp = 3.891e-10, -6.584e-11
    cN = np.cos(beta[:, 4])
    bp = tides.doodson_arguments(to_mjd(j["year"]), to_mjd(j["year"]))
    pm = j["pair_monthly"]
    ok = np.isfinite(pm)
    bm = tides.doodson_arguments(to_mjd(j["csr_year"][ok]), to_mjd(j["csr_year"][ok]))
    d = pm[ok] - j["csr"][ok] + 2 * ocean_amp * np.cos(bm[:, 4])
    out = {"year": 2000 + (mjd - 51_544.5) / 365.25, "solid": solid_amp * cN,
           "ocean_iers": ocean_amp * cN, "ocean_flip": -ocean_amp * cN,
           "fix_on_pair": -2 * ocean_amp * np.cos(bp[:, 4]),
           "corrected_offset": float(d.mean()), "corrected_err": float(d.std() / np.sqrt(d.size)),
           "source": "analytic"}
    if j_om1 is not None:
        out["fix_on_pair"] = j["pair_smooth"] - np.interp(j["year"], j_om1["year"], j_om1["pair_smooth"])
        out["corrected_offset"] = j_om1["offset"]
        out["corrected_err"] = j_om1["offset_err"]
        out["source"] = "re-run"
    return out


if __name__ == "__main__":
    main(force="--force" in sys.argv)
