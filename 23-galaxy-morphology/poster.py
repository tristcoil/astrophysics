"""Hubble's tuning fork, drawn with real galaxies chosen by the volunteers.

Hubble (1926, 1936) arranged galaxies in a "tuning fork": ellipticals along
the handle, ordered by how flattened they look (E0 round to E7 elongated),
lenticulars (S0) at the junction, and spirals along two prongs, unbarred (S)
and barred (SB), ordered from big bulges and tight arms (a) to small bulges
and open arms (c).

Galaxy Zoo 2 never asked "which Hubble type is this?" but its questions map
onto exactly the features Hubble used, so every class here is a *product of
vote fractions along a path through the decision tree*. No galaxy is picked
by eye: within each class the highest-scoring candidates are shown, among
galaxies that are big on the sky, nearby, and well classified.

* Ellipticals: smooth (p_smooth), with Hubble's E number from the measured
  axis ratio b/a of the SDSS de Vaucouleurs fit, En with n = 10 (1 - b/a).
* S0: a disk, not edge-on, no spiral arms, no bar, a prominent bulge.
* Sa / Sb / Sc: a face-on spiral without a bar, with a big bulge and tight
  arms / a medium bulge and medium arms / a small bulge and loose arms.
  SBa / SBb / SBc: the same, with a bar.
"""

from __future__ import annotations

import numpy as np

from config import (
    POSTER_MIN_R90_ARCSEC, POSTER_PER_CLASS, POSTER_Z_RANGE,
)

#: Elliptical classes: Hubble's En, with the axis-ratio window each is drawn from.
E_CLASSES = {"E0": (0.95, 1.00), "E2": (0.77, 0.83), "E4": (0.57, 0.63), "E6": (0.37, 0.43)}

FORK_UPPER = ("Sa", "Sb", "Sc")
FORK_LOWER = ("SBa", "SBb", "SBc")
EXTRAS = ("edge-on", "ring", "dust lane", "merger", "irregular")


def class_scores(c: dict, d: dict) -> tuple[dict, dict]:
    """Score (a path probability) and an extra eligibility mask for each class."""
    g = lambda k: np.nan_to_num(c["p_" + k].astype(float))  # noqa: E731
    ba = c["deVAB_r"].astype(float)
    faceon = g("features") * g("notedgeon")
    spiral = faceon * g("spiral")
    big_bulge = g("bulge_obvious") + g("bulge_dominant")
    mid_bulge = g("bulge_noticeable") + g("bulge_obvious")
    small_bulge = g("bulge_none") + g("bulge_noticeable")
    flat_enough = c["expAB_r"].astype(float) > 0.55  # arms and bars need a view from above
    scores, extra = {}, {}
    for name, (lo, hi) in E_CLASSES.items():
        scores[name] = g("smooth")
        extra[name] = (ba >= lo) & (ba <= hi)
    scores["S0"] = faceon * g("nospiral") * g("nobar") * big_bulge
    extra["S0"] = flat_enough
    for bar_key, names in (("nobar", FORK_UPPER), ("bar", FORK_LOWER)):
        base = spiral * g(bar_key)
        scores[names[0]] = base * big_bulge * g("arms_tight")
        scores[names[1]] = base * mid_bulge * g("arms_medium")
        scores[names[2]] = base * small_bulge * g("arms_loose")
        for n in names:
            extra[n] = flat_enough
    scores["edge-on"] = g("features") * g("edgeon")
    scores["ring"] = g("odd") * g("ring")
    scores["dust lane"] = g("odd") * g("dustlane")
    scores["merger"] = g("odd") * g("merger")
    scores["irregular"] = g("odd") * g("irregular")
    for n in EXTRAS:
        extra[n] = np.ones(len(ba), bool)
    return scores, extra


#: Votes needed on the questions that define a class, so that its score is
#: not the opinion of a handful of volunteers.
MIN_QUESTION_VOTES = 15
_QUESTIONS = {
    "round": ("round", "inbetween", "cigar"),
    "bulge": ("bulge_none", "bulge_noticeable", "bulge_obvious", "bulge_dominant"),
    "winding": ("arms_tight", "arms_medium", "arms_loose"),
    "bar": ("bar", "nobar"),
    "odd": ("ring", "lens", "disturbed", "irregular", "other", "merger", "dustlane"),
}
_NEEDS = {**{n: ("round",) for n in E_CLASSES}, "S0": ("bulge", "bar"),
          **{n: ("bulge", "winding", "bar") for n in FORK_UPPER + FORK_LOWER},
          "edge-on": (), **{n: ("odd",) for n in ("ring", "dust lane", "merger", "irregular")}}


def _enough_votes(c: dict, name: str) -> np.ndarray:
    ok = np.ones(len(c["dr7objid"]), bool)
    for q in _NEEDS[name]:
        ok &= sum(c["n_" + a] for a in _QUESTIONS[q]) >= MIN_QUESTION_VOTES
    return ok


def select(c: dict, d: dict, per_class: int = POSTER_PER_CLASS) -> dict:
    """The best-scoring galaxies of each class, by rule.

    Candidates are nearby (z in POSTER_Z_RANGE), big on the sky (petroR90 at
    least POSTER_MIN_R90_ARCSEC, so structure is resolved), classified by at
    least 30 volunteers, and not suspected of being an artefact. A galaxy is
    used once, in the class where it scores highest.
    """
    scores, extra = class_scores(c, d)
    base = (d["good"] & (d["z"] >= POSTER_Z_RANGE[0]) & (d["z"] < POSTER_Z_RANGE[1])
            & (d["R90_arcsec"] >= POSTER_MIN_R90_ARCSEC) & (c["n_class"] >= 30)
            & (np.nan_to_num(c["wf_artifact"]) < 0.1))
    names = list(scores)
    stack = np.vstack([np.where(base & extra[n] & _enough_votes(c, n), scores[n], -1.0)
                       for n in names])
    best_class = np.argmax(stack, axis=0)
    out = {}
    for k, n in enumerate(names):
        s = np.where(best_class == k, stack[k], -1.0)
        order = np.argsort(-s)[:per_class]
        order = order[s[order] > 0]
        out[n] = [{"index": int(i), "dr7objid": int(c["dr7objid"][i]), "score": float(s[i]),
                   "z": float(d["z"][i]), "r90": float(d["R90_arcsec"][i]),
                   "ba": float(c["deVAB_r"][i]), "n_class": int(c["n_class"][i])}
                  for i in order]
    return out


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

POSTER_BG = "#0d0d0d"
POSTER_INK = "#f2f1ec"
POSTER_INK2 = "#c3c2b7"
POSTER_MUTED = "#898781"
POSTER_LINE = "#52514e"


def _show(ax, path, label: str, sub: str, *, label_size: float = 13,
          sub_size: float = 8) -> None:
    from PIL import Image

    img = np.asarray(Image.open(path).convert("RGB"))
    # GZ2 images span 8.5 x the galaxy's 90%-light radius; show the middle 60%.
    n = img.shape[0]
    cut = int(round(n * 0.2))
    ax.imshow(img[cut:n - cut, cut:n - cut], interpolation="lanczos")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0.5, -0.06, label, transform=ax.transAxes, ha="center", va="top", color=POSTER_INK,
            fontsize=label_size, fontweight="bold")
    ax.text(0.5, -0.06 - 0.016 * label_size, sub, transform=ax.transAxes, ha="center", va="top",
            color=POSTER_INK2, fontsize=sub_size, linespacing=1.3)


def _save_dark(fig, stem, dpi: int) -> list:
    out = f"{stem}.png"
    fig.savefig(out, dpi=dpi, facecolor=POSTER_BG)
    import matplotlib.pyplot as plt

    plt.close(fig)
    return [out]


def draw(selection: dict, image_path, out_stem, *, n_volunteers: int, n_galaxies: int) -> list:
    """The poster: the fork itself, then the rest of the zoo underneath."""
    import matplotlib.pyplot as plt

    W, H = 16.0, 11.0
    fig = plt.figure(figsize=(W, H), facecolor=POSTER_BG)
    tile = 1.55  # inches

    def place(xc, yc, entry, label, sub, size=tile, label_size=13, sub_size=8):
        ax = fig.add_axes([(xc - size / 2) / W, (yc - size / 2) / H, size / W, size / H])
        ax.set_facecolor("black")
        _show(ax, image_path(entry["dr7objid"]), label, sub, label_size=label_size,
              sub_size=sub_size)
        return ax

    y_up, y_mid, y_lo = 8.2, 5.97, 3.78
    xs_e = [1.15, 2.95, 4.75, 6.55]
    x_s0 = 8.55
    xs_prong = [10.75, 12.55, 14.35]

    line = fig.add_axes([0, 0, 1, 1], facecolor="none")
    line.set_xlim(0, W)
    line.set_ylim(0, H)
    line.axis("off")
    kw = {"color": POSTER_LINE, "lw": 2.0, "solid_capstyle": "round", "zorder": 0}
    line.plot([xs_e[0], x_s0], [y_mid, y_mid], **kw)
    for y in (y_up, y_lo):
        line.plot([x_s0 + tile / 2, xs_prong[0] - tile / 2 - 0.3, xs_prong[-1]], [y_mid, y, y], **kw)
    line.text(xs_e[0] - tile / 2, y_mid + tile / 2 + 0.22, "ELLIPTICALS  ·  rounder to flatter",
              color=POSTER_MUTED, fontsize=8.5, va="bottom")
    line.text(xs_prong[0] - tile / 2, y_up + tile / 2 + 0.22,
              "SPIRALS  ·  bigger bulge, tighter arms → smaller bulge, looser arms",
              color=POSTER_MUTED, fontsize=8.5, va="bottom")
    line.text(xs_prong[0] - tile / 2, y_lo + tile / 2 + 0.22, "BARRED SPIRALS", color=POSTER_MUTED,
              fontsize=8.5, va="bottom")
    line.text(x_s0 - tile / 2, y_mid + tile / 2 + 0.22, "LENTICULAR", color=POSTER_MUTED,
              fontsize=8.5, va="bottom")

    for x, name in zip(xs_e, E_CLASSES):
        e = selection[name][0]
        place(x, y_mid, e, name, f"{100 * e['score']:.0f}% smooth · b/a {e['ba']:.2f}\n"
                                 f"z = {e['z']:.3f}")
    e = selection["S0"][0]
    place(x_s0, y_mid, e, "S0", f"disk, no arms, no bar\nscore {e['score']:.2f} · z = {e['z']:.3f}")
    for row, names, word in ((y_up, FORK_UPPER, "spiral"), (y_lo, FORK_LOWER, "barred spiral")):
        for x, name in zip(xs_prong, names):
            e = selection[name][0]
            place(x, row, e, name, f"{word} · score {e['score']:.2f}\nz = {e['z']:.3f}")

    line.text(0.6, 2.3, "…and the rest of the zoo", color=POSTER_INK2, fontsize=12,
              fontstyle="italic")
    xs_x = np.linspace(1.7, W - 1.7, len(EXTRAS))
    labels = {"edge-on": "edge-on disk", "ring": "ring", "dust lane": "dust lane",
              "merger": "merger", "irregular": "irregular"}
    for x, name in zip(xs_x, EXTRAS):
        e = selection[name][0]
        place(x, 1.45, e, labels[name], f"score {e['score']:.2f} · z = {e['z']:.3f}", size=1.15,
              label_size=10.5, sub_size=7.5)

    fig.text(0.6 / W, 10.55 / H, "Hubble's tuning fork, drawn by volunteers", color=POSTER_INK,
             fontsize=24, fontweight="bold", va="top")
    fig.text(0.6 / W, 9.98 / H,
             f"Real galaxies from the Sloan Digital Sky Survey, chosen by rule from the "
             f"{n_galaxies:,d} classified by {n_volunteers:,d} Galaxy Zoo 2 volunteers. Each class "
             "is a path through the volunteers' decision tree;\nshown is the galaxy with the "
             "highest product of vote fractions along it (ellipticals: the measured axis ratio b/a "
             "sets Hubble's E number, n = 10 (1 - b/a)). None was picked by eye.",
             color=POSTER_INK2, fontsize=10.5, va="top", linespacing=1.45)
    fig.text(0.6 / W, 0.1 / H,
             "Images: the gri colour composites the volunteers classified (SDSS; Galaxy Zoo 2, "
             "Willett et al. 2013; Zenodo record 3565489), central 60%. Votes debiased by Hart et "
             "al. (2016). Score: product of the debiased vote fractions along the class's path.",
             color=POSTER_MUTED, fontsize=8)
    return _save_dark(fig, out_stem, dpi=150)


def draw_gallery(selection: dict, image_path, out_stem, classes=None) -> list:
    """Every selected candidate, three per class: the choice was a ranking, not luck."""
    import matplotlib.pyplot as plt

    classes = classes or [*E_CLASSES, "S0", *FORK_UPPER, *FORK_LOWER, *EXTRAS]
    per = max(len(selection[n]) for n in classes)
    half = (len(classes) + 1) // 2
    cols = 2 * per + 1  # two groups of classes, a spacer column between them
    fig, axes = plt.subplots(half, cols, figsize=(cols * 1.75, half * 2.05), facecolor=POSTER_BG,
                             gridspec_kw={"width_ratios": [1] * per + [0.25] + [1] * per})
    for ax in axes.ravel():
        ax.set_facecolor(POSTER_BG)
        ax.axis("off")
    for k, name in enumerate(classes):
        r, c0 = k % half, (0 if k < half else per + 1)
        for j, e in enumerate(selection[name]):
            ax = axes[r, c0 + j]
            ax.axis("on")
            _show(ax, image_path(e["dr7objid"]), f"{name} #{j + 1}",
                  f"score {e['score']:.2f} · z {e['z']:.3f}", label_size=9, sub_size=7)
    fig.suptitle("The three best-scoring candidates of each class", color=POSTER_INK, fontsize=13,
                 fontweight="bold", y=0.995)
    fig.subplots_adjust(hspace=0.62, wspace=0.06, top=0.965, bottom=0.03, left=0.02, right=0.98)
    return _save_dark(fig, out_stem, dpi=110)
