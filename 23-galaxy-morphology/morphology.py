"""From volunteer votes to galaxy shapes.

Every GZ2 galaxy was classified by ~40 volunteers, each answering a tree of
questions (Willett et al. 2013, Fig. 1). For each answer the catalogue gives
the number of votes, the fraction of the volunteers *asked that question* who
gave it, and that fraction corrected for classification bias.

Two things make vote fractions subtle:

* **They are conditional.** Only volunteers who answered "features or disk"
  were asked "is it edge-on?", and only those who then answered "no" were asked
  about bars and spiral arms. So p_bar = 0.9 means "90% of those who saw a
  face-on disk saw a bar", not "90% sure it is a barred galaxy". The
  probability of the whole path is the product down the tree, e.g.
  p_features x p_not-edge-on x p_spiral for a face-on spiral (Hart et al. 2016).

* **They depend on distance.** A galaxy seen farther away is smaller and
  fainter in its image, so its arms and bars are harder to see and more
  volunteers call it smooth. The *debiased* fractions correct for this
  statistically: they assume that galaxies of the same luminosity and physical
  size have, on average, the same mix of shapes at every redshift, and rescale
  the votes until they do. We use Hart et al. (2016)'s debiasing, which the
  Galaxy Zoo team recommends; Willett et al. (2013)'s is kept for comparison.

A vote fraction is not a probability. It measures how many people saw a
feature, which mixes how obvious the feature is with how sure anyone could be.
Thresholds (> 0.8 for "clean" samples) turn fractions into classes at the cost
of leaving the ambiguous galaxies out; averages of fractions keep everyone but
assign no individual a class. Results here are checked both ways.
"""

from __future__ import annotations

import numpy as np

from config import (
    BAR_MIN_AXIS_RATIO, BAR_MIN_P_FEATURES, BAR_MIN_P_NOTEDGEON, BAR_MIN_VOTES,
    BAR_MIN_VOTES_DEBIASED, BAR_P, CLEAN_P,
    HART_MIN_ARM_VOTES, HART_SPIRAL_P,
)

ARM_LABELS = ("1", "2", "3", "4", "5+")
ARM_KEYS = ("arms_1", "arms_2", "arms_3", "arms_4", "arms_5plus")


def path_probabilities(c: dict, prefix: str = "p_") -> dict:
    """Probabilities of whole paths through the decision tree.

    ``prefix`` selects debiased (``p_``) or raw consistency-weighted (``wf_``)
    fractions. Missing answers (a question nobody was asked) count as zero.
    """
    f = {k: np.nan_to_num(c[prefix + k].astype(float)) for k in
         ("smooth", "features", "edgeon", "notedgeon", "spiral", "bar")}
    disk = f["features"]
    faceon = disk * f["notedgeon"]
    return {
        "smooth": f["smooth"],
        "features": disk,
        "edgeon_disk": disk * f["edgeon"],
        "faceon_disk": faceon,
        "spiral": faceon * f["spiral"],
        "barred_disk": faceon * f["bar"],
    }


def classes(c: dict) -> dict:
    """Boolean samples used throughout, from Hart et al. (2016) debiased fractions.

    * ``smooth`` / ``features``: Willett et al. (2013)'s clean samples, a
      debiased fraction above 0.8 for the first question. Smooth galaxies are
      mostly ellipticals; the disk-dominated lenticulars (S0) split between the
      two (Willett et al. 2013, section 5.2.1).
    * ``spiral``: Hart et al. (2016)'s spiral sample, p_features x
      p_not-edge-on x p_spiral > 0.5; ``spiral_clean`` the same above 0.8.
    * ``edgeon``: p_features x p_edge-on > 0.5.
    * ``disk``: disks seen face-on enough to show a bar: p_features x
      p_not-edge-on > 0.5, axis ratio b/a > 0.5 (inclination below ~60 deg, as
      in Masters et al. 2011), and at least ``BAR_MIN_VOTES_DEBIASED`` votes on
      the bar question. ``barred``: such a disk with debiased p_bar >= 0.5, the
      majority rule. Eligibility uses debiased fractions, because raw ones
      admit fewer disks the farther away they are; read this way the bar
      fraction of the volume-limited sample stays within 4 points from
      z = 0.03 to 0.085 (README, "Reading the votes").
    * ``oblique`` / ``barred_majority``: Willett et al. (2013)'s oblique disks
      (Table 3 thresholds on the raw fractions: p_features > 0.430,
      p_not-edge-on > 0.715, >= 20 bar votes) and the debiased majority rule
      p_bar >= 0.5 -- the definitions their expert comparison used.
    """
    pp = path_probabilities(c)
    n_bar_question = c["n_bar"] + c["n_nobar"]
    oblique = ((np.nan_to_num(c["wf_features"]) > BAR_MIN_P_FEATURES)
               & (np.nan_to_num(c["wf_notedgeon"]) > BAR_MIN_P_NOTEDGEON)
               & (n_bar_question >= BAR_MIN_VOTES))
    disk = ((pp["faceon_disk"] > 0.5) & (n_bar_question >= BAR_MIN_VOTES_DEBIASED)
            & (np.nan_to_num(c["expAB_r"].astype(float)) > BAR_MIN_AXIS_RATIO))
    return {
        "smooth": np.nan_to_num(c["p_smooth"]) > CLEAN_P,
        "features": np.nan_to_num(c["p_features"]) > CLEAN_P,
        "spiral": pp["spiral"] > HART_SPIRAL_P,
        "spiral_clean": pp["spiral"] > CLEAN_P,
        "edgeon": pp["edgeon_disk"] > 0.5,
        "disk": disk,
        "barred": disk & (np.nan_to_num(c["p_bar"]) >= BAR_P),
        "oblique": oblique,
        "barred_majority": oblique & (np.nan_to_num(c["p_bar"]) >= BAR_P),
    }


def arm_number(c: dict) -> tuple[np.ndarray, np.ndarray]:
    """Spiral arm number, Hart et al. (2016)'s way.

    Among galaxies with at least 5 votes naming a definite number of arms,
    the arm number is the answer with the highest debiased fraction ("can't
    tell" excluded). Returns (index into ARM_LABELS or -1, that fraction).
    """
    p = np.column_stack([np.nan_to_num(c["p_" + k]) for k in ARM_KEYS])
    n_definite = sum(c["n_" + k] for k in ARM_KEYS)
    enough = n_definite >= HART_MIN_ARM_VOTES
    modal = np.where(enough & (p.max(axis=1) > 0), p.argmax(axis=1), -1)
    return modal, p.max(axis=1)


def weighted_fraction(values, weights) -> tuple[float, float]:
    """Mean of a vote fraction over a set of galaxies, with its standard error.

    Bamford et al. (2009)'s "likelihood-weighted" fractions: each galaxy counts
    in proportion to its vote fraction instead of being forced into one class.
    The error is the standard error of the mean.
    """
    v = np.asarray(values, float)
    w = np.asarray(weights, float) if weights is not None else np.ones_like(v)
    ok = np.isfinite(v) & (w > 0)
    v, w = v[ok], w[ok]
    if w.sum() == 0:
        return float("nan"), float("nan")
    mean = float(np.sum(w * v) / np.sum(w))
    n_eff = w.sum() ** 2 / np.sum(w**2)
    var = np.sum(w * (v - mean) ** 2) / np.sum(w)
    return mean, float(np.sqrt(var / max(n_eff - 1.0, 1.0)))


def binomial_fraction(flags) -> tuple[float, float, int]:
    """Fraction of True values, its binomial standard error, and the count."""
    f = np.asarray(flags, bool)
    n = f.size
    if n == 0:
        return float("nan"), float("nan"), 0
    p = f.mean()
    return float(p), float(np.sqrt(max(p * (1 - p), 1.0 / n) / n)), int(n)


# ---------------------------------------------------------------------------
# Hubble types for the tuning-fork poster
# ---------------------------------------------------------------------------
#
# GZ2 never asked "what Hubble type is this?" but its questions map onto the
# features Hubble used. Ellipticals E0..E7 are ordered by apparent flattening
# (question 7: completely round / in between / cigar-shaped). Spirals Sa..Sc
# are ordered by bulge size and how tightly the arms wind (questions 5 and 10;
# Willett et al. 2013 show both track expert T-types). The bar question splits
# the fork's two prongs. Lenticulars (S0) are disks with no spiral arms.

def hubble_scores(c: dict) -> dict:
    """A score in [0, 1] for each tuning-fork class: the probability of the path.

    Each score multiplies the debiased fractions along the decision-tree path
    that defines the class; the poster shows each class's highest scorers.
    """
    g = lambda k: np.nan_to_num(c["p_" + k].astype(float))  # noqa: E731
    faceon = g("features") * g("notedgeon")
    spiral = faceon * g("spiral")
    unbarred, barred = spiral * g("nobar"), spiral * g("bar")
    early = g("bulge_dominant") + g("bulge_obvious")  # big bulge
    late = g("bulge_noticeable") + g("bulge_none")    # small bulge
    sa = early * g("arms_tight")
    sb = g("bulge_obvious") * g("arms_medium") + g("bulge_noticeable") * g("arms_medium")
    sc = late * g("arms_loose")
    s0 = faceon * g("nospiral")
    return {
        "E0": g("smooth") * g("round"),
        "E3": g("smooth") * g("inbetween"),
        "E7": g("smooth") * g("cigar"),
        "S0": s0 * g("nobar") * (g("bulge_obvious") + g("bulge_dominant")),
        "SB0": s0 * g("bar"),
        "Sa": unbarred * sa,
        "Sb": unbarred * sb,
        "Sc": unbarred * sc,
        "SBa": barred * sa,
        "SBb": barred * sb,
        "SBc": barred * sc,
        "edge-on": g("features") * g("edgeon"),
        "merger": g("odd") * g("merger"),
        "ring": g("odd") * g("ring"),
        "irregular": g("odd") * g("irregular"),
        "dust lane": g("odd") * g("dustlane"),
    }
