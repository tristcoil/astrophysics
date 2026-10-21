"""Run the whole project from cached data: python run.py [--force]

1. lunar laser ranging: our range model for 3,846 APOLLO normal points with two
   lunar ephemerides (DE440, EPM2021), then the residual fits (llr_fit);
2. the two labelled simulations: the Nordtvedt effect on the Moon (nordtvedt) and
   an equivalence-principle violation in the triple pulsar (triple_pulsar);
3. GR's predictions for the classic tests (gr_predictions) and the scorecard of
   published and atlas results (scorecard, literature);
4. every figure (figures), the interactive scorecard page (explorer) and
   data/results.json.

The model runs and simulations are cached in data/derived; --force recomputes
them (6.5 min here: ~3 min triple pulsar, ~1 min Nordtvedt, ~1 min lunar models
and the EPM2021 parse, ~1 min fits and bootstraps).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import config as C
import explorer
import figures
import gr_predictions
import literature as L
import llr_fit
import moon_compare
import nordtvedt
import scorecard
import triple_pulsar


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return o


def main(force: bool = False) -> dict:
    t0 = time.time()
    if force:
        for p in C.DERIVED_DIR.glob("*"):
            if p.name.startswith(("llr_model_", "llr_results", "llr_fit_tables", "nordtvedt_sim", "triple_sim", "moon_compare", "epm2021_records")):
                p.unlink()
    print("1. lunar laser ranging fits")
    llr = llr_fit.run_all(force=force)
    tables = pd.read_pickle(C.DERIVED_DIR / "llr_fit_tables.pkl")
    print("2. simulations")
    nord = nordtvedt.amplitude()
    tri = triple_pulsar.signature()
    tri_wn = triple_pulsar.white_noise_sensitivity()
    moon = moon_compare.compute()
    print("3. scorecard")
    table = scorecard.table(llr)
    table.to_csv(C.DERIVED_DIR / "scorecard_table.csv", index=False)
    print("4. figures")
    figs = []
    figs += figures.fig_scorecard(table)
    figs += figures.fig_report_card(table)
    figs += figures.fig_equivalence_principle()
    if L.gamma_table() is not None:
        figs += figures.fig_gamma_history()
    p = llr["EPM2021"]["params"]["gamma_minus_1"]
    figs += figures.fig_gamma_beta((p["value"], p["boot"]))
    figs += figures.fig_binary_pulsars()
    figs += figures.fig_triple_pulsar(tri, tri_wn)
    figs += figures.fig_lunar_ranging(tables, llr, nord["cosD"])
    figs += figures.fig_moon_along_track(moon_compare.curves(), llr)
    figs += figures.fig_nordtvedt(nord)
    figs.append(explorer.write(table))
    results = {
        "llr": llr,
        "nordtvedt_simulation": nord,
        "triple_pulsar_simulation": {**tri, "white_noise": tri_wn},
        "moon_along_track_vs_DE440": {k: {"tau_ms_2006_2010_2014_2018": [float(np.interp(y, v["year"], v["tau_ms"]))
                                                                          for y in (2006.5, 2010.5, 2014.5, 2018.5)],
                                         "max_remainder_m": v["rest_max_m"]} for k, v in moon.items()},
        "gr_predictions": gr_predictions.table(),
        "scorecard": table.drop(columns=["precision"]).to_dict(orient="records"),
        "scorecard_summary": {"n": int(len(table)), "within_1sigma": int((table["pull"].abs() <= 1).sum()),
                              "within_2sigma": int((table["pull"].abs() <= 2).sum()),
                              "max_pull": {"test": table.loc[table["pull"].abs().idxmax(), "label"],
                                           "pull": float(table.loc[table["pull"].abs().idxmax(), "pull"])}},
        "figures": [str(Path(f).resolve().relative_to(C.PROJECT_DIR)) for f in figs],
        "runtime_s": time.time() - t0,
    }
    C.RESULTS_PATH.write_text(json.dumps(_clean(results), indent=1))
    print(f"done in {time.time() - t0:.0f} s; results in {C.RESULTS_PATH.relative_to(C.PROJECT_DIR)}")
    return results


if __name__ == "__main__":
    main(force="--force" in sys.argv)
