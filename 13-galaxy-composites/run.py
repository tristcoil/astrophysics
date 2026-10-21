"""Run the whole project end to end.

    python run.py

Fetches whatever is not already cached, validates the data, builds the
pictures of all twelve galaxies, measures what their colours mean, writes every
figure, and prints a summary. Once the cache is warm this runs offline, and
galaxies are processed one at a time, so it needs well under 2 GB of memory.

    python run.py --only m51 --only m104   # just these galaxies (no gallery figure)
    python run.py --no-figures
"""

from __future__ import annotations

import argparse
import gc
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np

import analysis
import color_physics as cp
import imaging
from config import FIGURE_DIR, FLAGSHIP, GALAXIES, LUPTON_Q, LUPTON_SATURATION, galaxy_by_key
from fetch_data import fetch_all, missing_files

WHITE_POINT_GALAXY = "m81"


def _flagship_extras(g, res, art, ref) -> dict:
    """The extra measurements the explanatory figures need for the flagship."""
    fr = art["frame"]
    keep = ~art["star_mask"]
    for m in art["companions"].values():
        keep &= ~m
    art["sky"] = res.sky
    art["peak_mu_r"] = res.luminance["peak_mu_r"]
    return {
        "display_profile": analysis.display_color_profile(fr, art["rc3"], keep,
                                                          ref["white_balance"]["ours"]["weights"],
                                                          res.sky["r"]["std"]),
        "color_map": analysis.color_index_map(fr, art["star_mask"], res.ebv_sfd, sky=res.sky),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", action="append", default=[], help="analyse only this galaxy key (repeatable)")
    parser.add_argument("--no-figures", action="store_true", help="skip figure generation")
    args = parser.parse_args(argv)

    if missing_files():
        print("Some data are not cached yet; fetching")
        fetch_all()

    galaxies = [galaxy_by_key(k) for k in args.only] if args.only else GALAXIES
    print("Reference points: the Sun, the stars, the eye")
    ref = analysis.reference()
    eye = analysis.eye_matrix()
    ref["sun_white_xyz"] = cp.sun_white_xyz()
    ref["y_per_intensity"] = analysis.reference_xyz_per_intensity(eye)

    results, small, star_rows, jpeg_summaries, jpeg_pixels, noise_r = [], {}, [], [], [], []
    flag_art = flag_extra = None
    for g in galaxies:
        print(f"  {g.name:<9} ({g.common})")
        res, art = analysis.analyse_galaxy(g, eye, ref)
        comp = analysis.make_composites(art["frame"], ref, eye)
        res.composite = {"true_color_" + k: v for k, v in comp["true_diag"].items()}
        art["composite"] = comp
        if not args.no_figures:
            import figures
            figures.save_gallery_pngs({g.key: art})
        star_rows += analysis.star_photometry_rows(art["frame"])
        noise_r.append(res.sky["r"]["std"])
        summ, pix = analysis.validate_viewer_jpeg(g)
        jpeg_summaries.append(summ)
        jpeg_pixels.append(pix)
        small[g.key] = {"composite": {k: np.asarray(cm_small, dtype=np.float32) for k, cm_small in
                                      (("standard", _small(comp["standard"])), ("true", _small(comp["true"])))},
                        "frame": _FrameStub(g)}
        if g.key == FLAGSHIP:
            flag_extra = _flagship_extras(g, res, art, ref)
            flag_art = art
        else:
            del art
        results.append(res)
        gc.collect()

    print("Validation")
    stars = analysis.summarize_star_photometry(star_rows, float(np.median(noise_r)))
    stitch = [analysis.validate_stitching(k) for k in ("m81", "m82") if k in [g.key for g in galaxies]]
    brown_rows = [{"name": name, **{k: v for k, v in b.items() if k in ("mag_diff", "our_colors", "published_colors",
                                                                      "our_mags", "published_mags", "aperture")}}
                  for r in results for name, b in r.brown.items()]
    validation = {"jpeg": jpeg_summaries, "stitch": stitch, "stars": stars, "brown": brown_rows}
    context = analysis.context_luminances()

    written = []
    if not args.no_figures:
        import figures

        print("Figures")
        if flag_art is not None:
            written += figures.figure_how_made(flag_art, ref)
            written += figures.figure_stretch(flag_art, ref, flag_extra["display_profile"])
            written += figures.figure_true_color(flag_art, results, ref, eye.summary)
            written += figures.figure_too_faint(flag_art, results, context, analysis._dark_sky_luminance())
            written += figures.figure_colour_meaning(flag_art, results, ref, flag_extra["color_map"])
        if flag_art is not None and (analysis.DATA_DIR / "multiwavelength" / "m51_lsdr9north.fits").exists():
            caveat = analysis.pipeline_caveat(ref)
            from PIL import Image as _Image
            jpg = np.asarray(_Image.open(analysis.DATA_DIR / "multiwavelength" / "m51_lsdr9north_viewer.jpg")
                             .convert("RGB"))
            written += figures.figure_pipeline_caveat(caveat, flag_art["composite"]["standard"], jpg,
                                                      flag_art["frame"].galaxy.fov_arcmin)
            validation["pipeline_caveat"] = {k: v for k, v in caveat.items() if k != "images"}
        brown = cp.load_brown()
        written += figures.figure_filters_vs_eye(ref, analysis.eye_visibility_per_band(),
                                                 brown.spectra["NGC_5194"][:2])
        written += figures.figure_validation({**validation, "jpeg_pixels": jpeg_pixels}, ref)
        written += figures.save_swatches(results)
        if not args.only:
            written += figures.figure_gallery(results, small)
        if WHITE_POINT_GALAXY in [g.key for g in galaxies]:
            fr = imaging.load_frame(galaxy_by_key(WHITE_POINT_GALAXY))
            written += figures.figure_white_point({"frame": fr}, ref)
        try:
            import multiwave

            written += multiwave.make_figures(ref, eye)
            validation["alignment"] = getattr(multiwave.figure_multiwavelength, "alignment", {})
            validation["hubble"] = getattr(multiwave.figure_hubble, "alignment", {})
        except FileNotFoundError as err:
            print(f"  (multi-wavelength figures skipped: {err})")
        import interactive

        keys = [k for k in interactive.EXPLORER_GALAXIES if k in [g.key for g in galaxies]]
        if keys:
            frames = {k: imaging.load_frame(galaxy_by_key(k)) for k in keys}
            written.append(interactive.write_explorer(ref, eye, frames))
            del frames
        print(f"Wrote {len(written)} figure files to {FIGURE_DIR}")

    analysis.save_results(reference=ref, eye_matrix=eye.summary, validation={k: v for k, v in validation.items()},
                          galaxies=results, context_luminances=context,
                          stretch={"Q": LUPTON_Q, "saturation_nmgy_arcsec2": LUPTON_SATURATION})
    print_summary(ref, eye, validation, results, context)
    return 0


def _small(rgb, max_px=520):
    import figures
    return figures._downsample(rgb, max_px)


class _FrameStub:
    """Just enough of a Frame for the gallery figure (its galaxy)."""

    def __init__(self, g):
        self.galaxy = g


def print_summary(ref, eye, val, results, context) -> None:
    print("\n" + "=" * 96)
    print("G1 -- Colour composites of real galaxies")
    print("=" * 96)
    sun = ref["sun"]
    print("\nThe Sun through the SDSS filters (defines 'white' in every picture)")
    print("  absolute AB mag: " + ", ".join(f"{b} {sun['abs_mag_ab'][b]:.3f} ({sun['abs_mag_diff'][b]:+.3f} vs "
                                          f"Willmer)" for b in sun["abs_mag_diff"]))
    print(f"  colours g-r {sun['colors']['g-r']:.3f}, r-i {sun['colors']['r-i']:.3f}; colour temperature "
          f"{sun['eye']['cct_K']:.0f} K; illuminance at 1 AU {sun['illuminance_1au_lux']:,.0f} lux")
    wb = ref["white_balance"]
    print(f"  the survey viewer's weights make a {wb['viewer']['neutral_star']} star white "
          f"(g-r {wb['viewer']['neutral_color']['g-r']:.2f}); equal weights make an "
          f"{wb['flat']['neutral_star']} star white")

    s = eye.summary
    print(f"\nFrom three filters to the eye: 3x3 matrix fitted to {s['n_training_spectra']} galaxy spectra")
    print(f"  leave-one-out colour error du'v': median {s['loo_duv_median']:.4f}, 90% {s['loo_duv_p90']:.4f}, "
          f"worst {s['loo_duv_max']:.4f} ({s['loo_duv_worst_galaxy']}); one JND ~ 0.004")
    print(f"  on stars (never seen in training): median {s['stars_duv_median']:.4f}")

    print("\nValidation")
    j = val["jpeg"]
    print(f"  viewer JPEG, after the same compression: median |diff| "
          f"{np.median([x['compressed']['median_abs_diff'] for x in j]):.1f}, values within +-2: "
          f"{min(x['compressed']['frac_within_2'] for x in j) * 100:.3f}% (worst galaxy)")
    for cat in ("sdss", "gaia"):
        st = val["stars"][cat]
        print(f"  stars vs {cat.upper():<5}: " + ", ".join(
            f"{b} {st[b]['median']:+.3f} (sigma {st[b]['robust_sigma']:.3f}, n={st[b]['n']})" for b in "gri"))
    for st in val["stitch"]:
        print(f"  stitching {st['galaxy']}: shift " + ", ".join(
            f"{b} ({st[b]['shift_px'][0]:+.3f}, {st[b]['shift_px'][1]:+.3f}) px" for b in "gri"))
        print("      typical pixel change / noise " + ", ".join(
            f"{b} {st[b]['typical_diff_over_noise']:.2f}" for b in "gri") + "; bright pixels median "
            + ", ".join(f"{st[b]['bright_pixels_median_abs_frac'] * 100:.2f}%" for b in "gri")
            + "; flux ratio " + ", ".join(f"{st[b]['total_flux_ratio']:.5f}" for b in "gri"))
    if val.get("alignment"):
        print("  alignment with Gaia (slow foreground stars): " + ", ".join(
            f"{k} {v['median_total_arcsec']:.2f}\" (n={v['n_stars']})" for k, v in val["alignment"].items()
            if "median_total_arcsec" in v))
    if val.get("hubble"):
        h = val["hubble"]
        reg = h["registration_to_sdss"]
        print(f"  Hubble registered to SDSS: shifted {reg['east_arcsec']:.2f}\" east, {reg['north_arcsec']:.2f}\" north; "
              "brightness vs SDSS / predicted: " + ", ".join(
                  f"{k} {v['measured_over_predicted']:.3f}" for k, v in h["photometry_vs_sdss"].items()))
    if val.get("pipeline_caveat"):
        r = val["pipeline_caveat"]["retained"]
        print(f"  Legacy Surveys DR9-north keeps {r['g'] * 100:.0f}% (g), {r['r'] * 100:.0f}% (r), "
              f"{r['z'] * 100:.0f}% (z) of M51's light")
    if val["brown"]:
        d = np.array([[r["mag_diff"][b] for b in "gri"] for r in val["brown"]])
        c = np.array([[r["our_colors"][k] - r["published_colors"][k] for k in ("g-r", "r-i")] for r in val["brown"]])
        print(f"  Brown et al. apertures ({len(val['brown'])}): mag diff median {np.median(d):+.3f}, "
              f"rms {np.sqrt(np.mean(d**2)):.3f}; colour diff median |{np.median(np.abs(c)):.3f}|")

    print(f"\n  {'galaxy':<9}{'D25':>6}{'g-r':>7}{'r-i':>7}{'dered g-r':>11}{'eye colour':>12}{'CCT':>8}"
          f"{'vs spec':>9}{'peak cd/m2':>12}{'>0.005':>8}")
    for r in results:
        e = r.eye_from_image
        vs = f"{e['duv_vs_spectrum']:.4f}" if "duv_vs_spectrum" in e else "--"
        print(f"  {r.name:<9}{r.d25_arcmin:>5.1f}'{r.colors_d25['g-r']:>7.3f}{r.colors_d25['r-i']:>7.3f}"
              f"{r.colors_d25_dereddened['g-r']:>11.3f}{e['hex']:>12}{e['cct_K']:>7.0f}K{vs:>9}"
              f"{r.luminance['peak_cd_m2']:>12.4f}{r.luminance['frac_d25_above_mesopic_lower'] * 100:>7.2f}%")
    print("  (eye colour: the galaxy's total light as an eye adapted to sunlight would see it;"
          "\n   >0.005: share of the galaxy's area bright enough for cones to begin responding)")


if __name__ == "__main__":
    sys.exit(main())
