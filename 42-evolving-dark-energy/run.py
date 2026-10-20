"""Run the whole project end to end.

    python run.py                 # fetch (if needed), every fit, every figure, summary
    python run.py --figures-only  # redraw figures from the cached results
    python run.py --processes 4   # parallel MCMC workers (default 6)

Once ``fetch_data.py`` has filled ``data/`` this runs offline. The first full run
takes about an hour, almost all of it MCMC; the chains are cached afterwards.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings

warnings.filterwarnings("ignore")

from config import DESI_DCHI2  # noqa: E402


def _row(label, ours, sig, desi=None, dsig=None):
    s = f"  {label:<26} Δχ² = {ours:6.2f}  ({sig:.2f}σ)"
    if desi is not None:
        s += f"    DESI: {desi:5.2f} ({dsig:.2f}σ)"
    return s


def print_summary(r: dict) -> None:
    print("\n" + "=" * 78)
    print("RESULTS — is dark energy changing? (DESI DR2 BAO + CMB + supernovae)")
    print("=" * 78)
    v = r["validation"]
    d_un = [abs(x["chi2"] - x["desi_chi2"]) for x in v["rows"] if "chi2_cmb" not in x]
    d_cmb = [abs(x["chi2"] - x["desi_chi2"]) for x in v["rows"] if "chi2_cmb" in x]
    print(f"\nOur chi^2 at DESI's best-fit points: max |diff| {max(d_un):.4f} without the CMB, "
          f"{max(d_cmb):.3f} with it (theta_* offset {v['theta_offset_rel']:.1e}).")
    print("\nPreference for evolving dark energy (w0waCDM vs LCDM, MAP):")
    for row in r["desi_comparison"]:
        print(_row(row["label"], row["dchi2"], row["sigma"], row.get("desi_files_dchi2"),
                   row.get("desi_files_sigma")))
    for key, rec in r["maps"].items():
        if "dchi2" in rec and "desi_label" not in rec:
            print(_row(rec["label"], rec["dchi2"], rec["sigma"]))
    print("\nFull-CMB values DESI quotes (not reproducible with a compressed prior):")
    for k, (d, s) in DESI_DCHI2.items():
        if "full" in k:
            print(f"  {k:<26} Δχ² = {-d:5.1f}  ({s}σ)")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--figures-only", action="store_true")
    ap.add_argument("--processes", type=int, default=6)
    ap.add_argument("--force-mcmc", action="store_true")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)

    import analysis

    if args.figures_only:
        results = json.loads(analysis.RESULTS_PATH.read_text())
    else:
        from fetch_data import fetch_all

        fetch_all()
        results = analysis.run_core(processes=args.processes, force_mcmc=args.force_mcmc)
    print_summary(results)
    if not args.no_figures:
        import figures

        figures.make_all(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
