"""How much of each element came promptly, how much late? The two-process model.

Weinberg et al. (2019, 2022) and Griffith et al. (2019, 2022) describe every
disc star as the sum of two ingredients: a *prompt* one (core-collapse
supernovae of massive stars, within ~40 Myr) and a *delayed* one (Type Ia
supernovae, ~0.1 Gyr to many Gyr later). Each element X has a fixed recipe
(q_cc, q_Ia) at a given metallicity; each star has its own amounts (A_cc,
A_Ia). Normalised so a solar star has A_cc = A_Ia = 1,

    10^[X/Mg] = q_cc^X(z) + q_Ia^X(z) * A_Ia / A_cc,     z = 10^[Mg/H].

Magnesium is taken as purely prompt. Stars on the high-alpha plateau have no
Type Ia iron at all, so [Fe/Mg] there fixes iron's prompt share:
q_cc^Fe = 10^(-[Mg/Fe]_plateau). For any star, A_Ia/A_cc then follows from its
[Fe/Mg]. Comparing the median [X/Mg] of the two disc sequences at the same
[Mg/H] (they differ only in A_Ia/A_cc) gives two equations for q_cc^X and
q_Ia^X, and the prompt fraction of X in a solar-like star:

    f_cc^X = q_cc / (q_cc + q_Ia * A_Ia/A_cc),  on the high-Ia sequence at [Mg/H] = 0.

The plateau is measured here (the knee fit), not assumed. Everything is done
separately for dwarfs and giants, with each class's [X/Fe] scale zero-pointed
so its low-alpha sequence is solar at [Mg/H] = 0.

The "delayed" ingredient is whatever tracks Type Ia iron in time. For the
neutron-capture elements that is mostly AGB stars (the s-process), and for
europium possibly neutron-star mergers; the model cannot tell those apart,
only prompt from delayed.
"""

from __future__ import annotations

import numpy as np

#: [Mg/H] bins (Griffith et al. 2022 used 0.1 dex bins in -0.5 < [Mg/H] < 0.5).
MGH_EDGES = np.round(np.arange(-0.5, 0.51, 0.1), 2)
MIN_PER_BIN = 30
#: A weak line is "detected" (flag 0) mostly in stars where the element is
#: enhanced, so where many stars are upper limits the median of the
#: detections is biased high -- by different amounts on the two sequences.
#: A bin is used only if both sequences have this fraction detected.
MIN_DETECTED = 0.7


def zero_points(d, sel, high_ia, elements, *, window=0.05):
    """Offsets that put the high-Ia (low-alpha) sequence at solar ratios at
    [Mg/H] = 0. Returns offsets for [Mg/Fe] and for each [X/Mg]."""
    near = sel & high_ia & (np.abs(d["mg_h"]) < window)
    zp = {"mg_fe": float(np.nanmedian(d["mg_fe"][near]))}
    for el in elements:
        x = d[f"{el}_fe"] - d["mg_fe"]
        zp[el] = float(np.nanmedian(x[near & np.isfinite(x)])) if np.isfinite(x[near]).sum() > 20 else np.nan
    return zp


def ratio_amplitude(fe_mg, q_cc_fe):
    """A_Ia / A_cc from [Fe/Mg] (iron's process vector is metallicity independent)."""
    return (10 ** fe_mg - q_cc_fe) / (1.0 - q_cc_fe)


def process_vectors(mg_h, x_mg, fe_mg, low_ia, high_ia, q_cc_fe, *, edges=MGH_EDGES):
    """q_cc and q_Ia of one element, bin by bin in [Mg/H], from the medians of
    the two sequences."""
    rows = []
    det = np.isfinite(x_mg)
    for lo, hi in zip(edges[:-1], edges[1:]):
        inbin = (mg_h >= lo) & (mg_h < hi)
        b = inbin & det
        lo_s, hi_s = b & low_ia, b & high_ia
        if lo_s.sum() < MIN_PER_BIN or hi_s.sum() < MIN_PER_BIN:
            continue
        if (lo_s.sum() / max((inbin & low_ia).sum(), 1) < MIN_DETECTED
                or hi_s.sum() / max((inbin & high_ia).sum(), 1) < MIN_DETECTED):
            continue
        xl, xh = np.median(x_mg[lo_s]), np.median(x_mg[hi_s])
        rl = ratio_amplitude(np.median(fe_mg[lo_s]), q_cc_fe)
        rh = ratio_amplitude(np.median(fe_mg[hi_s]), q_cc_fe)
        # 10^xl = qcc + qIa rl ; 10^xh = qcc + qIa rh
        q_ia = (10 ** xh - 10 ** xl) / (rh - rl)
        q_cc = 10 ** xl - q_ia * rl
        rows.append(dict(mg_h=0.5 * (lo + hi), x_low=float(xl), x_high=float(xh), r_low=float(rl),
                         r_high=float(rh), q_cc=float(q_cc), q_ia=float(q_ia),
                         f_cc_high=float(q_cc / (q_cc + q_ia * rh)), n_low=int(lo_s.sum()),
                         n_high=int(hi_s.sum())))
    return rows


def f_cc_solar(rows):
    """Prompt fraction on the high-Ia sequence at [Mg/H] = 0 (the two bins
    straddling zero, averaged)."""
    near = [r for r in rows if abs(r["mg_h"]) < 0.1]
    if not near:
        return np.nan
    return float(np.mean([r["f_cc_high"] for r in near]))


def decompose(d, sel, low_ia_mask, elements, plateau, *, n_boot=100, seed=0, extra=None):
    """Prompt fractions of every element for one stellar class.

    ``sel`` picks the class (e.g. disc dwarfs); ``low_ia_mask`` is True for
    high-alpha stars. ``plateau`` is the class's [Mg/Fe] plateau on the same
    scale as d["mg_fe"] (before zero-pointing). Errors: bootstrap over stars.
    """
    high_ia = sel & ~low_ia_mask
    low_ia = sel & low_ia_mask
    zp = zero_points(d, sel, high_ia, elements)
    mg_fe = d["mg_fe"] - zp["mg_fe"]
    plat = plateau - zp["mg_fe"]
    q_cc_fe = 10 ** (-plat)
    fe_mg = -mg_fe
    mg_h = d["mg_h"]
    rng = np.random.default_rng(seed)
    out = {"_zero_points": zp, "_plateau_zp": float(plat), "_q_cc_fe": float(q_cc_fe)}
    # Extra elements come as (values, coverage): coverage marks the stars the
    # catalogue looked at, the denominator of the detection fraction (the
    # europium VAC covers only high-S/N spectra, a selection unrelated to Eu).
    extra = extra or {}
    for k, (v, cov) in extra.items():
        near = sel & high_ia & (np.abs(d["mg_h"]) < 0.05) & np.isfinite(v)
        zp[k] = float(np.median((v - d["mg_fe"])[near])) if near.sum() > 20 else np.nan
    for el in ["fe"] + list(elements) + list(extra):
        cover = np.ones_like(sel)
        if el == "fe":
            x_mg = fe_mg
        elif el in extra:
            x_mg = extra[el][0] - d["mg_fe"] - zp[el]
            cover = extra[el][1]
        else:
            x_mg = d[f"{el}_fe"] - d["mg_fe"] - zp[el]
        s_el = sel & cover
        idx_el = np.flatnonzero(s_el)
        rows = process_vectors(mg_h[s_el], x_mg[s_el], fe_mg[s_el], low_ia[s_el], high_ia[s_el],
                               q_cc_fe)
        f = f_cc_solar(rows)
        boots = []
        for _ in range(n_boot if rows else 0):
            ii = idx_el[rng.integers(0, idx_el.size, idx_el.size)]
            r = process_vectors(mg_h[ii], x_mg[ii], fe_mg[ii], low_ia[ii], high_ia[ii], q_cc_fe)
            boots.append(f_cc_solar(r))
        boots = np.array(boots)
        n_x = int(np.isfinite(x_mg[s_el]).sum())
        out[el] = dict(f_cc=f, f_cc_err=float(np.nanstd(boots)) if np.isfinite(boots).sum() > 5 else np.nan,
                       n=n_x, rows=rows)
    return out
