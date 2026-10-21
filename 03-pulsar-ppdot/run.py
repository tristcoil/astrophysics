"""Run the whole project end to end.

    python run.py

Downloads the catalogue if it is not already cached, sorts every pulsar into its
family, measures the populations, writes the figures, and prints a summary of
every number the README quotes. Once the cache is warm this runs offline.
"""

from __future__ import annotations

import argparse
import sys
import warnings

warnings.filterwarnings("ignore")

import analysis
from config import FIGURE_DIR
from fetch_data import fetch_catalogue, load_manifest


def _pct(p: dict, digits: int = 0) -> str:
    lo, hi = p["ci95"]
    return (f"{p['fraction']:.{digits}%} ({p['k']}/{p['n']}; 95% CI "
            f"{lo:.{digits}%}-{hi:.{digits}%})")


def print_summary(r: dict) -> None:
    m = load_manifest()
    c = r["census"]
    print("\n" + "=" * 78)
    print(f"ATNF Pulsar Catalogue v{r['catalogue_version']} "
          f"(retrieved {m['retrieved_utc'][:10]})")
    print("=" * 78)
    print(f"  {c['n_total']:,} pulsars; {c['n_with_period']:,} with a period; "
          f"{c['n_pdot_positive']:,} with Pdot > 0, {c['n_pdot_negative']} with Pdot < 0, "
          f"{c['n_period_no_pdot']:,} with no Pdot yet")
    print(f"  P-Pdot sample (Pdot > 0, outside globular clusters): {c['n_in_ppdot']:,}")
    for k, v in c["by_class_ppdot"].items():
        print(f"    {analysis.CLASS_NAMES[k]:<34}{v:>6,}")

    print("\nRECYCLING -- binaries among pulsars with a timing solution")
    bf = r["binary_fractions"]
    print(f"  millisecond pulsars:  {_pct(bf['msp'])}")
    print(f"  ordinary pulsars:     {_pct(bf['non_recycled'], 1)}")
    print(f"  (MSPs incl. untimed:  {_pct(bf['msp_including_untimed'])})")
    print("  median spin period by companion type:")
    for code, v in r["companions"].items():
        print(f"    {v['name']:<38}{v['median_P_ms']:>8.1f} ms   (n = {v['n']})")

    print("\nDEATH VALLEY -- radio pulsars relative to the death lines")
    dv = r["death_valley"]
    print(f"  above the dipole line: {_pct(dv['above_dipole_line'])}")
    print(f"  inside the valley:     {_pct(dv['in_valley'])}")
    print(f"  beyond the valley:     {_pct(dv['beyond_valley'], 1)}")
    print(f"    slowest of those: {', '.join(dv['beyond_valley_slowest'][:5])}")

    print("\nRECORD HOLDERS")
    ex = r["extremes"]
    f = ex["fastest"]
    print(f"  fastest spin:  {f['jname']} in {f.get('cluster', '?')}, "
          f"{ex['fastest_spin_hz']:.2f} Hz")
    print(f"     -> mean density > {ex['fastest_min_density_g_cm3']:.2e} g/cm3; "
          f"radius < {ex['fastest_max_radius_km_1p4']:.1f} km for 1.4 Msun; "
          f"equator at {ex['fastest_equator_fraction_c_10km']:.2f} c for R = 10 km")
    fo = ex["fastest_outside_clusters"]
    print(f"  fastest outside a cluster: {fo['jname']}, {1 / fo['P_s']:.1f} Hz")
    print(f"  slowest:       {ex['slowest']['jname']}, P = {ex['slowest']['P_s']:.2f} s")
    print(f"  strongest B:   {ex['strongest_field']['jname']}, "
          f"{ex['strongest_field']['B_gauss']:.2e} G")
    print(f"  youngest (rotation-powered): {ex['youngest_rotation_powered']['jname']}, "
          f"tau = {ex['youngest_rotation_powered']['tau_yr']:.0f} yr")
    print(f"  most powerful: {ex['most_powerful']['jname']}, "
          f"Edot = {ex['most_powerful']['edot']:.2e} erg/s")

    print("\nTHE CRAB -- characteristic age against history")
    cr = r["crab"]
    print(f"  at the timing epoch {cr['epoch_year']:.1f}: tau = {cr['characteristic_age_yr']:.0f} yr, "
          f"true age {cr['true_age_yr']:.0f} yr (+{cr['overestimate_percent']:.0f}%)")
    print(f"  braking index from f0, f1, f2: n = {cr['braking_index']:.3f} "
          f"(Lyne et al. 2015: 2.342)")
    print(f"  implied birth period: {cr['birth_period_ms_n3']:.1f} ms (n = 3), "
          f"{cr['birth_period_ms_measured_n']:.1f} ms (measured n)")

    print("\nWHERE THE NUMBERS STOP MEANING WHAT THEY SAY")
    ma = r["msp_ages"]
    sk = r["shklovskii_example"]
    print(f"  MSPs with tau > age of the Universe: {_pct(ma['older_than_universe'])}")
    print(f"  Shklovskii: median {ma['median_kinematic_fraction_of_pdot']:.0%} of an MSP's "
          f"observed Pdot is motion; {sk['jname']}: {sk['kinematic_fraction']:.0%}")
    gc = r["globular_clusters"]
    print(f"  globular clusters: {gc['n']} pulsars in {gc['n_clusters']} clusters; "
          f"negative Pdot in {_pct(gc['negative_pdot'])}")
    bi = r["braking_indices"]
    print(f"  braking indices: |n| > 10 for {_pct(bi['abs_n_above_10'])}")

    print("\nTHE GALAXY")
    h = r["heights"]
    for row in h["by_age"]:
        print(f"  tau {row['tau_lo_yr']:>8.0e}-{row['tau_hi_yr']:<8.0e} yr: median |z| = "
              f"{row['median_z_pc']:5.0f} pc  (n = {row['n']})")
    lat = h["latitude"]
    print(f"  median |b|: ordinary {lat['non_recycled']['median_abs_b_deg']:.1f} deg, "
          f"MSPs {lat['msp']['median_abs_b_deg']:.1f} deg "
          f"({lat['msp_not_fast_gpps']['median_abs_b_deg']:.1f} deg without FAST GPPS "
          f"discoveries, which are {h['msp_found_by_fast_gpps']['fraction']:.0%} of field MSPs)")

    print("\nCONSISTENCY -- our formulas against psrqpy's derived columns")
    for k, v in r["consistency"].items():
        print(f"  {k:<11} max relative difference {v['max_rel_diff']:.1e} over {v['n']:,} pulsars")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    fetch_catalogue(quiet=True)
    table, results = analysis.analyse()
    analysis.save_results(results)

    if not args.no_figures:
        import figures

        written = figures.make_all(table, results)
        print(f"wrote {len(written)} figure files to {FIGURE_DIR}")

    print_summary(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
