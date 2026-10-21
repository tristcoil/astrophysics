"""Run the whole project end to end.

    python run.py

Downloads whatever is not already cached (~265 MB the first time), analyses the
hour of data around GW150914, writes the figures and audio, and prints our
measurements beside the published ones. Once the cache is warm this runs
offline, in a few minutes.
"""

from __future__ import annotations

import argparse
import warnings

warnings.filterwarnings("ignore")

import analysis
import figures
import sonify
from config import AUDIO_DIR, FIGURE_DIR, PRL
from fetch_data import fetch_all, load_event


def _row(label, ours, published, note=""):
    print(f"  {label:<34}{ours:>24}   {published:<26}{note}")


def print_summary(r: dict) -> None:
    arr, ch, cp, sig, pub = r["arrival"], r["chirp"], r["compactness"], r["significance"], r["published"]
    h1, l1 = r["detectors"]["H1"], r["detectors"]["L1"]
    lag_pub = pub["predicted_lag_ms"]
    mc_pub = pub["chirp_mass_detector_frame"]

    print("\n" + "=" * 96)
    print(f"{r['event']} — our measurements vs published values")
    print("=" * 96)
    _row("quantity", "ours", "published")
    print("  " + "-" * 92)

    print("  The signal")
    _row("matched-filter SNR, H1 / L1", f"{h1['matched_filter']['snr']:.1f} / {l1['matched_filter']['snr']:.1f}", "—")
    snrs = pub["network_snr_by_search"]
    _row("network SNR", f"{r['network_snr']:.1f}",
         f"{min(snrs.values())}–{max(snrs.values())} (GWTC-1 searches)")
    _row("peak strain, filtered (H1)", f"{h1['peak_strain_filtered']:.2e}", f"{PRL['peak_strain']:.1e} (PRL)")
    _row("chi^2/dof at the event, H1 / L1",
         f"{h1['matched_filter']['chisq'] / h1['matched_filter']['chisq_dof']:.2f} / "
         f"{l1['matched_filter']['chisq'] / l1['matched_filter']['chisq_dof']:.2f}", "~1 expected for a signal")
    for name, d in (("H1", h1), ("L1", l1)):
        m = d["paper_figure1"]["observed"]
        _row(f"PRL Fig. 1 reproduced, {name}", f"r = {m['correlation']:.3f}, amp ×{m['rms_ratio_ours_over_paper']:.3f}",
             "1, ×1 if identical")

    print("  Arrival")
    _row("H1 after L1, matched filter (ms)", f"{arr['lag_mf_ms']:.2f} ± {arr['lag_mf_err_ms']:.2f}",
         f"{PRL['arrival_lag_ms']} +{PRL['arrival_lag_ms_upper']}/{PRL['arrival_lag_ms_lower']} (PRL)")
    _row("H1 after L1, cross-correlation (ms)", f"{arr['lag_xcorr_ms']:.2f} ± {arr['lag_xcorr_err_ms']:.2f}",
         f"{lag_pub['median']:.2f} [{lag_pub['p5']:.2f}, {lag_pub['p95']:.2f}] (GWTC-1 sky)")
    _row("relative sign (H1 vs L1)", f"{'inverted' if arr['xcorr_peak'] < 0 else 'same'} "
         f"(r = {arr['xcorr_peak']:+.2f})", "inverted (PRL)")
    _row("light travel time (ms)", f"{r['geometry']['light_travel_ms']:.2f}", f"{PRL['light_travel_ms']} (PRL)")

    print("  Chirp mass (detector frame, Msun)")
    _row("Newtonian zero-crossing fit", f"{ch['raw']:.1f}", "~37 (same method, Ann. Phys. 2017)")
    _row("... calibrated by injections", f"{ch['calibrated']:.1f} [{ch['calibrated_68'][0]:.1f}, {ch['calibrated_68'][1]:.1f}]",
         f"{mc_pub['median']:.1f} [{mc_pub['p5']:.1f}, {mc_pub['p95']:.1f}] (GWTC-1)", "  68% vs 90%")
    _row("method bias on known waveform", f"×{ch['bias_factor']:.2f}", "—")

    print("  Compactness (conservative inputs)")
    _row("minimum total mass (Msun)", f"{cp['min_total_mass']:.1f}", f">~{PRL['total_mass_bound_msun']:.0f} (PRL)")
    _row("sum of Schwarzschild radii (km)", f"{cp['schwarzschild_sum_km']:.0f}", f">~{PRL['schwarzschild_sum_km']:.0f} (PRL)")
    _row("separation at peak (km)", f"{cp['separation_km']:.0f}", f"~{PRL['separation_km']:.0f} (PRL)")
    _row("orbital speed at peak (v/c)", f"{cp['orbital_speed_c']:.2f}", "—")

    print("  Significance")
    _row("ranking statistic rho_c", f"{sig['event_statistic']:.1f}", f"{PRL['detection_statistic']} (PRL)")
    _row("loudest of all time slides", f"{sig['loudest_background']:.1f} ({sig['slides']:,} slides)", "—")
    _row("false-alarm rate (per year)", f"< {sig['far_per_year_bound']:.1f}", f"< 1/{PRL['far_bound_years']:,} (PRL)")
    _row("significance", f">= {sig['sigma_bound']:.1f} sigma", f"> {PRL['significance_sigma']} sigma (PRL)",
         "  bound set by data length")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--no-figures", action="store_true", help="skip figures")
    parser.add_argument("--no-audio", action="store_true", help="skip the WAV files")
    args = parser.parse_args(argv)

    print("Data (downloaded only if missing)")
    fetch_all(quiet=True)

    print("Analysis")
    results, art = analysis.analyse(load_event())
    art = analysis.compact_artefacts(art)

    if not args.no_figures:
        written, extra = figures.make_all(results, art)
        results["figure_checks"] = extra
        print(f"  wrote {len(written)} figure files to {FIGURE_DIR.name}/")
    if not args.no_audio:
        written = sonify.make_all(results, art)
        print(f"  wrote {len(written)} audio files to {AUDIO_DIR.name}/")

    analysis.save_results(results)
    print_summary(results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
