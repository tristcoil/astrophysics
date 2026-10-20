"""Run the whole project: measurements, images, figures, and the comparison table.

    python run.py                  # everything except the neural network
    python run.py --cnn            # also train and test the neural network (slow first time)
    python run.py --no-figures     # measurements and results.json only

The first run downloads the few hundred galaxy images the figures need (and,
with --cnn, the training set); after that everything runs offline.
"""

from __future__ import annotations

import argparse
import time
import warnings

import numpy as np

import analysis as an
import fetch_data as fd
import literature as lit
import morphology as mo
import physics as ph
from config import FIGURE_DIR

# astropy-healpix warns about NaN neighbours of pixels at the poles; harmless.
warnings.filterwarnings("ignore", message="invalid value encountered in neighbours_nested")


def say(msg: str) -> None:
    print(msg, flush=True)


def build(quiet: bool = False) -> dict:
    """All measurements. Returns results plus the per-galaxy arrays figures need."""
    t0 = time.time()
    c, meta = fd.load_catalog()
    mgs, mgs_meta = fd.load_mgs()
    na10 = fd.load_nair_abraham()
    if not quiet:
        say(f"catalogue: {len(c['dr7objid']):,d} GZ2 galaxies, {mgs_meta['n_galaxies']:,d} tracers")

    d = an.derive(c)
    kz = an.red_kcorrection(d)
    vl = an.volume_limit(d, kz)
    mass_info = an.mass_completeness(d, kz)
    cls = mo.classes(c)

    m = __import__("environment").derive_mgs(mgs)
    divider = an.colour_divider(m["M"], m["gr0"], m["z"] < 0.2)
    if not quiet:
        say(f"  colour divider (D2's method, D2's data): g-r = {divider['slope']:.5f} M_r "
            f"{divider['intercept']:+.5f}")
    dens = an.densities(c, d, mgs)
    sigma = dens["sigma"]

    res = {
        "catalog": meta, "tracers": mgs_meta,
        "sample": {
            "n_good": int(d["good"].sum()),
            "volume_limited": {k: v for k, v in vl.items() if k != "sel"},
            "mass_bins": [{k: v for k, v in b.items() if k != "sel"} for b in mass_info],
            "classes_all": {k: int(v[d["good"]].sum()) for k, v in cls.items()},
            "classes_vl": {k: int(v[vl["sel"]].sum()) for k, v in cls.items()},
        },
        "colour_divider": {"slope": divider["slope"], "intercept": divider["intercept"]},
        "density": {k: v for k, v in dens.items() if not k.startswith("_")
                    and k not in ("sigma", "is_tracer")},
    }
    if not quiet:
        say("  redshift bias")
    res["redshift_bias"] = an.redshift_bias(c, d, kz, vl)
    if not quiet:
        say("  colours by shape")
    res["wrong_colours"] = an.wrong_colours(c, d, cls, divider, vl, kz)
    if not quiet:
        say("  bars")
    res["bars"] = an.bars(c, d, cls, divider, vl, mass_info)
    res["bars"]["masters11_n20"] = _masters11_votes(c, d)
    if not quiet:
        say("  spiral arms")
    res["arms"] = an.spiral_arms(c, d, cls, sigma)
    if not quiet:
        say("  environment")
    res["environment"] = an.environment_trends(c, d, cls, divider, mass_info, sigma)
    res["environment"]["robustness"] = _environment_robustness(c, d, cls, divider, mass_info,
                                                               sigma)
    if not quiet:
        say("  volunteers versus experts")
    res["experts"] = an.experts(c, d, na10)
    if not quiet:
        say(f"  measurements done in {time.time() - t0:.0f}s")
    return {"results": res, "c": c, "d": d, "cls": cls, "vl": vl, "mass_info": mass_info,
            "sigma": sigma, "divider": divider, "kz": kz}


def _masters11_votes(c: dict, d: dict) -> list:
    """Masters et al. (2011)'s bar fraction against the vote-count threshold.

    Their sample used GZ2 votes collected up to July 2009, about half of the
    final total, and required >= 10 votes on the bar question. With the final
    catalogue the same threshold admits galaxies that only reached 10 votes
    because twice as many volunteers saw them, so the threshold is varied.
    """
    nq = c["n_bar"] + c["n_nobar"]
    base = (d["good"] & (d["z"] > 0.01) & (d["z"] < 0.06) & (d["M70"] < -19.38)
            & (d["log_ab"] < 0.3))
    rows = []
    for nmin in (10, 15, 20, 25, 30):
        s = base & (nq >= nmin)
        rows.append({"min_votes": nmin, **an.frac((c["n_bar"] >= c["n_nobar"])[s])})
    return rows


def _environment_robustness(c, d, cls, divider, mass_info, sigma) -> dict:
    """Does the environment decomposition survive the obvious objections?

    * Dust: edge-on disks are red because of dust, whatever their stars do.
      Drop them.
    * Redshift: each mass bin spans its own redshift range, and the debiased
      votes drift slightly with redshift. Use one common, narrow range.
    """
    out = {}
    not_edgeon = ~cls["edgeon"]
    mi = [dict(b, sel=b["sel"] & not_edgeon) for b in mass_info]
    out["no_edgeon"] = an.environment_trends(c, d, cls, divider, mi, sigma, n_boot=100)
    common = [dict(b, sel=b["sel"] & (d["z"] < 0.037)) for b in mass_info]
    out["z_below_0.037"] = an.environment_trends(c, d, cls, divider, common, sigma, n_boot=100)
    return {k: [{"mass_range": b["mass_range"], "n": b["n"],
                 "decomposition": b["decomposition"],
                 "decomposition_err": b["decomposition_err"]} for b in v["bins"]]
            for k, v in out.items()}


# ---------------------------------------------------------------------------
# The comparison table
# ---------------------------------------------------------------------------

def table(res: dict) -> list[str]:
    rows = []

    def row(what, ours, published):
        rows.append(f"  {what:<58} {ours:>24}   {published}")

    wc, b, a, ex = res["wrong_colours"], res["bars"], res["arms"], res["experts"]
    pct = lambda f: f"{100 * f['f']:.1f} +/- {100 * f['err']:.1f}%"  # noqa: E731
    rows.append("Colours by shape")
    m10 = wc["masters10_their_cut"]
    row("Red face-on disky spirals (Masters+10 cuts)", pct(m10["red"]),
        f"{100 * lit.MASTERS10['red_fraction']:.0f}% (294 of 5,433)")
    row("  bar fraction, red vs blue", f"{100 * m10['bar_fraction_red']['f']:.0f} vs "
        f"{100 * m10['bar_fraction_blue']['f']:.0f}%", "67-70 vs 27% (by eye, two authors)")
    s09 = wc["schawinski09"]
    row("Blue early types (Schawinski+09 cuts)", pct(s09["blue"]), "5.7 +/- 0.4%")
    rows.append("Bars")
    m11 = {r["min_votes"]: r for r in b["masters11_n20"]}
    row("Masters+11 cuts, >= 10 bar votes", pct(m11[10]), "29.4 +/- 0.5% (13,665 disks)")
    row("Masters+11 cuts, >= 20 bar votes (vote count of their snapshot)", pct(m11[20]),
        f"n = {m11[20]['n']:,d}")
    row("Our bar fraction, volume-limited disks", pct(b["vl"]["barred"]), "")
    rows.append("Spiral arms (Hart+16 samples)")
    row("Two-armed, luminosity-limited", pct(a["m2_lum"]), "62.1 +/- 0.4%")
    row("Two-armed, mass-limited", pct(a["m2_mass"]), "64.2 +/- 0.4%")
    top = a["vs_mass"][-1]
    row(f"Two-armed, highest mass bin (log M = {top['mass']:.2f})",
        f"{100 * top['f2']:.1f} +/- {100 * top['e2']:.1f}%", "53 +/- 2% (11.2 +/- 0.1)")
    row(f"Five or more arms, highest mass bin", f"{100 * top['f5+']:.1f} +/- {100 * top['e5+']:.1f}%",
        "15 +/- 2%")
    lo, hi = a["vs_colour"][0], a["vs_colour"][-1]
    row(f"Two-armed, bluest / reddest colour bin", f"{100 * lo['f2']:.0f} / {100 * hi['f2']:.0f}%",
        "32 / 84%")
    rows.append("Volunteers versus experts (Nair & Abraham 2010)")
    raw, w13 = ex["raw"], ex["willett13"]
    row("Ellipticals: median p_smooth (raw votes)",
        f"{raw['E_median_psmooth']:.2f} +/- {raw['E_psmooth_sd']:.2f}", "0.86 +/- 0.07")
    row("Ellipticals: highest p_features", f"{raw['E_max_pfeatures']:.3f}", "0.741")
    row("Experts' barred disks with GZ2 p_bar < 0.5", f"{100 * w13['bars']['barred_below_05']:.1f}%",
        "39%")
    row("Strong / intermediate / weak bars with p_bar > 0.8",
        f"{100 * w13['bars']['strong_above_08']:.0f} / {100 * w13['bars']['intermediate_above_08']:.0f}"
        f" / {100 * w13['bars']['weak_above_08']:.0f}%", "78 / 40 / 9%")
    rows.append("Environment (D2's Sigma_5; sparsest to densest fifth): rise in the red fraction")
    for bb in res["environment"]["bins"]:
        dd = bb["decomposition"]
        row(f"  log M* {bb['mass_range'][0]:.1f}-{bb['mass_range'][1]:.1f}: total / of it at fixed shape",
            f"{100 * dd['delta_red']:+.1f} / {100 * dd['red_at_fixed_morphology']:+.1f} pts",
            "Bamford+09: mostly at fixed shape")
    if "cnn" in res:
        cn = res["cnn"]
        rows.append("Neural network (test galaxies never seen in training)")
        row("Smooth or features, clear majorities: our CNN", f"{100 * cn['cnn']['features']['accuracy_confident']:.1f}%",
            "")
        if cn.get("zoobot", {}).get("features"):
            row("  same galaxies: Zoobot (Walmsley+23, DESI images)",
                f"{100 * cn['zoobot']['features']['accuracy_confident']:.1f}%", "")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--cnn", action="store_true", help="train and test the neural network")
    ap.add_argument("--offline", action="store_true", help="never touch the network")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    state = build(quiet=args.quiet)
    res = state["results"]
    if not args.no_figures:
        import figures

        figures.make_all(state, offline=args.offline, quiet=args.quiet)
    if args.cnn:
        import cnn

        res["cnn"] = cnn.run(state, offline=args.offline, quiet=args.quiet)
    if not args.no_figures:
        import interactive

        interactive.build(state, offline=args.offline, quiet=args.quiet)
    an.save_results(res)
    print("\nComparison with the literature")
    for line in table(res):
        print(line)
    print(f"\nFigures in {FIGURE_DIR}")


if __name__ == "__main__":
    main()
