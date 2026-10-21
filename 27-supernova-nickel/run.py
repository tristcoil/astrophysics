"""Run the whole project end to end.

    python run.py                  # download what is missing, analyse, draw, summarise
    python run.py --figures-only   # redraw figures from the last analysis run

The first run downloads ~50 MB (and ~700 small ZTF light curves, ~20 min).
Once the cache is warm everything runs offline in about two minutes.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
import config as C
import figures
import interactive
import literature as LIT
from fetch_data import fetch_all

CACHE = C.CACHE_DIR / "analysis.pkl"


def _v(d, nd=3):
    if d["err"] is None:
        return f"{d['value']:.{nd}f}"
    return f"{d['value']:.{nd}f} ± {d['err']:.{nd}f}"


def print_summary(s: dict) -> None:
    print("\n" + "=" * 78)
    print("RESULTS: radioactive nickel-56 in four supernovae and the ZTF sample")
    print("=" * 78)
    for sn in ("SN2011fe", "SN2014J"):
        x = s[sn]
        print(f"\n{C.SUPERNOVAE[sn]['label']} (Type Ia)")
        print(f"  peak luminosity        {_v(x['L_peak_1e43'])} x 1e43 erg/s on day {x['t_peak_days']}")
        print(f"  Arnett's rule          M(Ni) = {_v(x['m_ni_arnett_rule'])} Msun")
        print(f"  Arnett's model fit     M(Ni) = {_v(x['m_ni_arnett_fit'])} Msun, "
              f"T0 = {_v(x['t0_days'], 1)} d, diffusion time {_v(x['tau_m_days'], 1)} d")
        print(f"  tail (T0 fixed)        M(Ni) = {_v(x['m_ni_tail_t0_fixed'])} Msun")
        print(f"  Katz integral          M(Ni) = {x['m_ni_katz']:.3f} Msun")
        if sn == "SN2014J":
            print(f"  other dust law         Arnett's rule M(Ni) = {x['m_ni_arnett_rule_alt_dust']:.3f} Msun "
                  f"(peak {x['L_peak_alt_dust_1e43']:.2f} x 1e43 erg/s)")
        for lit in LIT.NICKEL[sn]:
            e = f" ± {lit['err']}" if lit["err"] else ""
            print(f"  published: {lit['ref']}: {lit['m_ni']}{e} ({lit['method']})")
    g = s["SN2014J_gamma"]
    print("\nSN 2014J from INTEGRAL's cobalt-56 gamma-ray lines")
    for k, lab in (("transparent", "ejecta transparent (lower limit)"), ("uniform", "uniform ejecta"),
                   ("exponential", "W7-like ejecta")):
        extra = ""
        if "escape_75d" in g[k]:
            extra = f"   escape at day 75: {g[k]['escape_75d']['847']:.2f} (847 keV)"
        print(f"  {lab:<34} M(Ni) = {_v(g[k]['m_ni'])} Msun{extra}")
    print(f"  check (Churazov et al. 2014): {LIT.CHECKS['SN2014J_gamma'][0]} ± {LIT.CHECKS['SN2014J_gamma'][1]} Msun")
    for sn, keys in (("SN2023ixf", ("tail_measured_nir", "tail_full_year")),
                     ("SN1987A", ("tail_full_trapping", "tail_with_escape"))):
        print(f"\n{C.SUPERNOVAE[sn]['label']} (Type II)")
        for k in keys:
            x = s[sn][k]
            t0 = f", T0 = {_v(x['t0_days'], 1)} d" if x["t0_days"] else " (all gamma rays trapped)"
            print(f"  days {x['window'][0]:.0f}-{x['window'][1]:.0f}: M(Ni) = {_v(x['m_ni'], 4)} Msun{t0}")
        for lit in LIT.NICKEL[sn][:2]:
            e = f" ± {lit['err']}" if lit["err"] else ""
            print(f"  published: {lit['ref']}: {lit['m_ni']}{e}")
    z = s["ztf"]
    print(f"\nZTF Type Ia, Arnett's rule: {z['n_measured']} of {z['n_candidates']} measured, "
          f"{z['n_low_dust']} low-dust (peak g-r <= {z['blue_cut_gr']}), {z['n_outliers']} outliers")
    print(f"  median M(Ni) = {z['m_ni_median']:.2f} Msun (68%: {z['m_ni_p16']:.2f}-{z['m_ni_p84']:.2f}); "
          f"median rise to bolometric peak {z['rise_bol_median']:.1f} d")
    print(f"  brightness-decline (Phillips) slope {z['phillips_slope_mag_per_mag']:.2f} mag/mag, scatter "
          f"{z['phillips_scatter_mag']:.2f} mag; corr(Delta m15, log M_Ni) = {z['corr_dm15_logmni']:.2f}")
    print()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from the cached analysis")
    args = ap.parse_args()
    if args.figures_only:
        if not CACHE.exists():
            print("No cached analysis yet: run `python run.py` first.")
            return 1
        res = pickle.loads(CACHE.read_bytes())
    else:
        fetch_all()
        res = analysis.run()
        for sn in C.SUPERNOVAE:  # the GP objects are bulky and not needed for figures
            res[sn]["bolo"].pop("bol", None)
        CACHE.write_bytes(pickle.dumps(res))
        C.RESULTS_PATH.write_text(json.dumps(res["summary"], indent=1, default=str))
    paths = figures.make_all(res)
    paths.append(interactive.write(res))
    print_summary(res["summary"])
    print("Figures:")
    for p in paths:
        if str(p).endswith((".png", ".html")):
            print("  ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
