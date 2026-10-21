"""Run the whole project end to end.

    python run.py                  # download what is missing, measure, infer H0, draw, summarise
    python run.py --figures-only   # redraw figures and the page from cached results
    python run.py --workers 2      # fewer CPU cores for the delay measurement

The first run downloads ~190 MB (mostly lens-model posteriors) and measures
13 delays with ~3,600 simulated light curves (about an hour on 4 cores; every
result is cached in results/ and reused). After that everything runs offline.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import warnings

warnings.filterwarnings("ignore")

import config as C

PY = sys.executable


def _v(t, nd=1):
    return f"{t[0]:.{nd}f} (-{t[1]:.{nd}f}/+{t[2]:.{nd}f})"


def summary():
    import figures as F

    print("\n" + "=" * 78)
    print("RESULTS: H0 from lensed-quasar time delays")
    print("=" * 78)
    print("\nDelays t(second) - t(first) in days: ours (joint spline) vs published")
    for key, a, b in F.DELAY_ROWS:
        comb = F._j(f"{key}_combined")
        v = F.pair_from(comb["taus"], comb["ref"], a, b)
        e = F.pair_err(comb, a, b)
        pv, pem, pep, src = F.published_delay(key, a, b)
        print(f"  {C.LENSES[key]['name']:<16} {a}->{b}: {v:+8.2f} +- {e:5.2f}   published {pv:+7.1f} "
              f"(-{pem:.1f}/+{pep:.1f})  [{src}]")
    res = json.loads((C.RESULTS / "h0.json").read_text())
    print("\nH0 (km/s/Mpc), flat LCDM: published lens models with published delays / with our delays")
    for key in C.H0LICOW_SIX + ["ALL"]:
        p, o = res["published"][key], res["ours"][key]
        w = C.LENSES[key]["wong20_h0"] if key != "ALL" else C.H0_LIT["H0LiCOW 2020"]
        print(f"  {key:<8} {_v(p)}  (Wong+2020: {_v(w)})   ours {_v(o)}")
    td = res["tdcosmo2025"]
    print(f"\nTDCOSMO 2025 chains: LCDM1a {_v(td['LambdaCDM1a']['h0'])}, lambda_int "
          f"{td['LambdaCDM1a']['lambda_int'][0]:.3f}; with our delays (estimate) x{td['our_delay_shift_factor']:.4f}"
          f" -> {td['LambdaCDM1a_with_our_delays_estimate'][0]:.1f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    if not args.figures_only:
        subprocess.run([PY, "fetch_data.py"], check=True)
        subprocess.run([PY, "measure_delays.py", "--workers", str(args.workers)], check=True)
        subprocess.run([PY, "analysis_h0.py"], check=True)
    subprocess.run([PY, "figures.py"], check=True)
    subprocess.run([PY, "interactive.py"], check=True)
    summary()


if __name__ == "__main__":
    main()
