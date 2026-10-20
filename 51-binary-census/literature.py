"""Published results we compare against, each with its source.

Every number was read from the paper's own arXiv source (LaTeX, or the PDF for
Raghavan et al. 2010, which has no TeX source), with section, table or equation
given. Nothing here is an input to a measurement unless a function says it
adopts it (the Moe & Di Stefano distributions are used as the *assumed*
intrinsic mass-ratio and eccentricity distributions in the completeness model,
and that choice is varied as a systematic).
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Raghavan et al. 2010, ApJS 190, 1: 454 F6-K3 stars within 25 pc (Hipparcos)
# ---------------------------------------------------------------------------

#: Sec. 5.3.3 and Fig. 13: a Gaussian fitted to the histogram of the 259
#: confirmed companions' periods (all hierarchy levels pooled; periods of
#: common-proper-motion pairs estimated from separations), P in days. Not
#: corrected for selection. Median 252 yr.
RAGHAVAN_LOGP_MEAN = 5.03
RAGHAVAN_LOGP_SIGMA = 2.28
RAGHAVAN_N_PRIMARIES = 454
RAGHAVAN_N_COMPANIONS_FIG13 = 259

#: Table 16: observed fractions of single/double/triple/higher-order systems;
#: with candidate companions; and the completeness-adjusted single fraction
#: (the only completeness-adjusted value given). 46% "of stars having
#: companions" (Sec. 5.3.9) is 100 - 54.
RAGHAVAN_FRACTIONS_OBSERVED = {"single": 0.56, "binary": 0.33, "triple": 0.08, "higher": 0.03}
RAGHAVAN_FRACTIONS_WITH_CANDIDATES = {"single": 0.54, "binary": 0.34, "triple": 0.09, "higher": 0.03}
RAGHAVAN_SINGLE_CORRECTED = (0.54, 0.02)
RAGHAVAN_MULTIPLICITY = (0.46, 0.02)

#: Sec. 5.3.2: companions by colour, F6-G2 (B-V 0.5-0.625, N=131) vs G2-K3 (N=323).
RAGHAVAN_MF_BLUE = (0.50, 0.04)
RAGHAVAN_MF_RED = (0.41, 0.03)

#: Sec. 5.3.4: "Pairs with orbital periods longer than the circularization
#: limit of about 12 days have a roughly flat eccentricity distribution out to
#: e ~ 0.6"; Sec. 5.3.5: mass ratios "roughly flat" for q = 0.2-0.95, and of 27
#: twins (q > 0.95) only six have P > ~200 yr, none > ~1000 yr (~115 AU).
RAGHAVAN_PCIRC_DAYS = 12.0

#: Raghavan's conversion for pairs without orbits (from Duquennoy & Mayor 1991):
#: log a" = log rho" + 0.13, i.e. a = 1.35 x projected separation.
RAGHAVAN_A_OVER_S = 10 ** 0.13

# ---------------------------------------------------------------------------
# Tokovinin 2014, AJ 147, 86/87: 4,847 F5-G6 dwarfs within 67 pc ("FG-67")
# ---------------------------------------------------------------------------

#: Paper II, Table 2 ("all pairs"): maximum-likelihood, selection-aware
#: log-normal in x = log10(P/day); epsilon = companions per primary.
TOKOVININ_ALL = {"x0": 4.54, "x0_err": 0.06, "sigma": 2.40, "sigma_err": 0.06, "eps": 0.571}
TOKOVININ_OUTER = {"x0": 4.93, "x0_err": 0.06, "sigma": 2.34, "sigma_err": 0.06, "eps": 0.464}
TOKOVININ_MULTIPLICITY = (0.46, 0.01)
#: Paper II, L1043: fractions of systems with n = 1, 2, 3, 4, 5+ components
#: (from simulations). This, not Raghavan et al., is the source of the often
#: quoted "54:33:8:4" split.
TOKOVININ_FN = (0.543, 0.323, 0.076, 0.043, 0.013)

# ---------------------------------------------------------------------------
# Moe & Di Stefano 2017, ApJS 230, 15
# ---------------------------------------------------------------------------

#: Table 11: solar-type f_logP;q>0.3 (companions with q > 0.3 per decade of P)
#: at log P = 0.5 ... 7.5.
MDS_SOLAR_TABLE11_LOGP = np.arange(0.5, 8.0, 1.0)
MDS_SOLAR_TABLE11_F = np.array([0.017, 0.020, 0.025, 0.045, 0.067, 0.077, 0.057, 0.052])

#: Table 13, f_logP;q>0.1 in four decades, by primary class (value, error).
MDS_TABLE13_CLASSES = {
    # name: (M1 range, representative M1)
    "solar": ((0.8, 1.2), 1.0), "A/late-B": ((2, 5), 3.5), "mid-B": ((5, 9), 7.0),
    "early-B": ((9, 16), 12.0), "O": ((16, 60), 28.0),
}
MDS_TABLE13_FLOGP_Q01 = {
    # log P centre: values per class in the order above
    1.0: [(0.027, 0.009), (0.07, 0.02), (0.14, 0.04), (0.19, 0.06), (0.29, 0.08)],
    3.0: [(0.057, 0.016), (0.12, 0.04), (0.22, 0.07), (0.26, 0.09), (0.32, 0.11)],
    5.0: [(0.095, 0.018), (0.13, 0.03), (0.20, 0.06), (0.23, 0.07), (0.30, 0.09)],
    7.0: [(0.075, 0.015), (0.09, 0.02), (0.11, 0.03), (0.13, 0.04), (0.18, 0.05)],
}
#: Table 13: multiplicity frequency (companions per primary, q > 0.1, 0.2 < log P < 8)
#: and the single / binary / higher-order fractions.
MDS_TABLE13_FMULT = [(0.50, 0.04), (0.84, 0.11), (1.3, 0.2), (1.6, 0.2), (2.1, 0.3)]
MDS_TABLE13_CLOSE = [(0.15, 0.03), (0.37, 0.08), (0.63, 0.13), (0.8, 0.2), (1.0, 0.2)]  # log P < 3.7
MDS_TABLE13_SINGLE = [(0.60, 0.04), (0.41, 0.08), (0.24, 0.08), (0.16, 0.09), (0.06, 0.06)]


def mds_flogp_q03(logp, m1):
    """Moe & Di Stefano (2017) Eqs. 20-23: companions with q > 0.3 per decade of
    period, for a primary of mass m1 (Msun), valid for 0.2 <= log P < 8."""
    logp = np.asarray(logp, dtype=float)
    lm = np.log10(m1)
    f1 = 0.020 + 0.04 * lm + 0.07 * lm ** 2          # Eq. 20, log P < 1
    f27 = 0.039 + 0.07 * lm + 0.01 * lm ** 2         # Eq. 21, log P = 2.7
    f55 = 0.078 - 0.05 * lm + 0.04 * lm ** 2         # Eq. 22, log P = 5.5
    alpha, dlp = 0.018, 0.7
    out = np.full(logp.shape, np.nan)
    a = (logp >= 0.2) & (logp < 1.0)
    out[a] = f1
    b = (logp >= 1.0) & (logp < 2.7 - dlp)
    out[b] = f1 + (logp[b] - 1) / (1.7 - dlp) * (f27 - f1 - alpha * dlp)
    c = (logp >= 2.7 - dlp) & (logp < 2.7 + dlp)
    out[c] = f27 + alpha * (logp[c] - 2.7)
    d = (logp >= 2.7 + dlp) & (logp < 5.5)
    out[d] = f27 + alpha * dlp + (logp[d] - 2.7 - dlp) / (2.8 - dlp) * (f55 - f27 - alpha * dlp)
    e = (logp >= 5.5) & (logp < 8.0)
    out[e] = f55 * np.exp(-0.3 * (logp[e] - 5.5))
    return out


def mds_gamma_largeq(logp, m1):
    """Eq. 9 (solar), 10 (3.5 Msun), 11 (> 6 Msun), interpolated in mass."""
    logp = np.asarray(logp, dtype=float)
    g_sol = np.where(logp < 5.0, -0.5, -0.5 - 0.3 * (logp - 5))
    g_35 = np.select([logp < 1.0, logp < 4.5, logp < 6.5],
                     [-0.5, -0.5 - 0.2 * (logp - 1.0), -1.2 - 0.4 * (logp - 4.5)], -2.0)
    g_6 = np.select([logp < 1.0, logp < 2.0, logp < 4.0],
                    [-0.5, -0.5 - 0.9 * (logp - 1), -1.4 - 0.3 * (logp - 2)], -2.0)
    return _interp_mass(m1, g_sol, g_35, g_6)


def mds_gamma_smallq(logp, m1):
    """Eq. 13 (solar), 14 (3.5 Msun), 15 (> 6 Msun), interpolated in mass."""
    logp = np.asarray(logp, dtype=float)
    g_sol = np.full(logp.shape, 0.3)
    g_35 = np.select([logp < 2.5, logp < 5.5], [0.2, 0.2 - 0.3 * (logp - 2.5)], -0.7 - 0.2 * (logp - 5.5))
    g_6 = np.select([logp < 1.0, logp < 3.0, logp < 5.6],
                    [0.1, 0.1 - 0.15 * (logp - 1), -0.2 - 0.50 * (logp - 3)], -1.5)
    return _interp_mass(m1, g_sol, g_35, g_6)


def mds_ftwin(logp, m1):
    """Eqs. 5-7: excess fraction of twins (q > 0.95) relative to the power law."""
    logp = np.asarray(logp, dtype=float)
    f_short = 0.30 - 0.15 * np.log10(m1)
    logp_twin = 8.0 - m1 if m1 <= 6.5 else 1.5
    out = np.where(logp < 1, f_short, f_short * (1 - (logp - 1) / (logp_twin - 1)))
    return np.clip(np.where(logp >= logp_twin, 0.0, out), 0, None)


def mds_eta(logp, m1):
    """Eqs. 17-18: eccentricity power-law index p(e) ~ e^eta (0 < e < e_max)."""
    logp = np.maximum(np.asarray(logp, dtype=float), 0.51)
    late = 0.6 - 0.7 / (logp - 0.5)
    early = 0.9 - 0.2 / (logp - 0.5)
    w = np.clip((m1 - 3.0) / 4.0, 0, 1)
    return (1 - w) * late + w * early


def _interp_mass(m1, g_sol, g_35, g_6):
    if m1 <= 1.2:
        return g_sol
    if m1 <= 3.5:
        w = (m1 - 1.2) / (3.5 - 1.2)
        return (1 - w) * g_sol + w * g_35
    if m1 <= 6.0:
        w = (m1 - 3.5) / (6.0 - 3.5)
        return (1 - w) * g_35 + w * g_6
    return g_6


# ---------------------------------------------------------------------------
# Multiplicity against primary mass
# ---------------------------------------------------------------------------

#: Offner et al. 2023 (PPVII), Table 1: multiplicity fraction MF and companion
#: frequency CF, corrected for incompleteness. (M1 range in Msun, MF, err, CF, err, source)
OFFNER23_MF = [
    ((0.075, 0.15), 0.19, 0.03, 0.21, 0.03, "Winters+2019"),
    ((0.15, 0.30), 0.23, 0.02, 0.27, 0.03, "Winters+2019"),
    ((0.30, 0.60), 0.30, 0.02, 0.38, 0.03, "Winters+2019"),
    ((0.75, 1.25), 0.46, 0.03, 0.60, 0.04, "Raghavan+2010"),
    ((0.85, 1.5), 0.47, 0.03, 0.62, 0.04, "Tokovinin 2014"),
    ((1.6, 2.4), 0.68, 0.07, 0.99, 0.13, "Moe & Kratter 2021"),
    ((3, 5), 0.81, 0.06, 1.28, 0.17, "Moe & Di Stefano 2017"),
    ((5, 8), 0.89, 0.05, 1.55, 0.24, "Moe & Di Stefano 2017"),
    ((8, 17), 0.93, 0.04, 1.8, 0.3, "Moe & Di Stefano 2017"),
    ((17, 50), 0.96, 0.04, 2.1, 0.3, "Sana et al."),
]
#: Offner et al. 2023, Table 1, wide-binary fraction WBF (a companion at
#: 100-10,000 au), where given.
OFFNER23_WBF = [
    ((0.15, 0.30), 0.048, 0.012, "Winters+2019"),
    ((0.30, 0.60), 0.12, 0.02, "Winters+2019"),
    ((0.75, 1.25), 0.22, 0.02, "Raghavan+2010"),
    ((0.85, 1.5), 0.19, 0.02, "Tokovinin 2014"),
    ((1.6, 2.4), 0.28, 0.04, "De Rosa+2014"),
    ((3, 5), 0.40, 0.06, "Moe & Di Stefano 2017"),
    ((5, 8), 0.49, 0.07, "Moe & Di Stefano 2017"),
    ((8, 17), 0.55, 0.08, "Moe & Di Stefano 2017"),
]

#: Winters et al. 2019 (M dwarfs within 25 pc): corrected multiplicity and
#: companion rates; separation distribution peaks at 20 au (sigma_log a 1.16).
WINTERS19_MR = (0.268, 0.014)
WINTERS19_CR = (0.324, 0.014)

# ---------------------------------------------------------------------------
# Twins among wide binaries: El-Badry et al. 2019, MNRAS 489, 5822
# ---------------------------------------------------------------------------

#: Table G1: F_twin (excess of q > ~0.95 relative to the power law, as a
#: fraction of all binaries with q > 0.3), median and 2-sigma errors, by
#: separation bin (AU) and primary mass bin (Msun).
ELBADRY19_SEP_BINS = [(50, 350), (350, 600), (600, 1000), (1000, 2500), (2500, 5000), (5000, 15000),
                      (15000, 50000)]
ELBADRY19_MASS_BINS = [(0.1, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.2), (1.2, 2.5)]
ELBADRY19_FTWIN = np.array([  # rows: separation bins; columns: mass bins; (median, +2sig, -2sig)
    [(0.189, 0.042, 0.049), (0.112, 0.034, 0.038), (0.059, 0.019, 0.021), (0.101, 0.027, 0.029), (0.086, 0.058, 0.075)],
    [(0.105, 0.036, 0.035), (0.062, 0.018, 0.020), (0.043, 0.016, 0.019), (0.054, 0.014, 0.018), (0.088, 0.037, 0.047)],
    [(0.049, 0.032, 0.039), (0.054, 0.017, 0.020), (0.032, 0.016, 0.019), (0.035, 0.015, 0.016), (0.015, 0.024, 0.034)],
    [(0.023, 0.025, 0.033), (0.049, 0.018, 0.022), (0.024, 0.016, 0.016), (0.013, 0.010, 0.011), (0.002, 0.015, 0.020)],
    [(-0.001, 0.048, 0.039), (0.047, 0.025, 0.029), (0.018, 0.020, 0.019), (0.023, 0.017, 0.021), (-0.003, 0.023, 0.023)],
    [(0.005, 0.040, 0.043), (0.043, 0.030, 0.031), (0.014, 0.018, 0.021), (0.007, 0.015, 0.015), (0.003, 0.023, 0.035)],
    [(0.024, 0.094, 0.120), (-0.009, 0.039, 0.043), (-0.012, 0.028, 0.032), (0.009, 0.024, 0.032), (0.001, 0.037, 0.057)],
])

# ---------------------------------------------------------------------------
# El-Badry, Rix & Heintz 2021, MNRAS 506, 2269
# ---------------------------------------------------------------------------

#: Table 1: 1,817,594 candidates, 1,256,400 with R < 0.1 (~0.4% of which are
#: expected chance alignments). Sec. 4.1: contrast limit ~4 mag at 1" and ~7 mag
#: at 2"; below ~0.7" only near-equal pairs.
ELBADRY21_N_CANDIDATES = 1_817_594
ELBADRY21_N_R01 = 1_256_400
ELBADRY21_CONTRAST = [(1.0, 4.0), (2.0, 7.0)]  # (theta arcsec, Delta G mag)
