"""The galaxy census: masses, star-formation rates, and who has stopped.

Per-galaxy quantities come from the MPA-JHU catalogue:

* **Stellar mass** from fitting stellar-population models to the five-band
  SDSS photometry (Kauffmann et al. 2003; Salim et al. 2007).
* **Star-formation rate (SFR)**. Inside the 3" fibre, from the H-alpha line
  (corrected for dust with the H-alpha/H-beta ratio) for galaxies whose lines
  are powered by young stars; from the 4000 A break for galaxies whose lines
  are weak or powered by an active nucleus. Outside the fibre, from fitting
  the photometry of the rest of the galaxy (Brinchmann et al. 2004, revised
  for DR7). The specific SFR (sSFR = SFR / M) is the rate at which a galaxy
  would double its stars: log sSFR = -10 means every 10 billion years.

The functions here build the samples, classify galaxies by what powers
their emission lines, find the star-forming main sequence and the green
valley, and measure the fraction of galaxies that have stopped forming stars.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

import cosmo_physics as cp
from config import (
    BPT_MIN_SN,
    LINE_ERROR_SCALE,
    MS_BIN_DEX,
    MS_FIT_LOGM,
)

# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------


def derive(g: dict) -> dict:
    """Per-galaxy arrays: magnitudes, distances, MPA-JHU values with validity.

    MPA-JHU marks missing values with -9999 (and SkyServer's LEFT JOIN leaves
    galaxies without a galSpec entry at -9999 too). Salim et al. (2016) found
    that when the mass is missing but a total SFR is present, the "total" SFR
    is really the fibre SFR, so an SFR counts as valid only with a valid mass.
    """
    z = g["z"].astype(float)
    r = (g["petroMag_r"] - g["extinction_r"]).astype(float)
    gr = (g["dered_g"] - g["dered_r"]).astype(float)
    # 453 galaxies have broken photometry (observed g - r beyond -0.5 to 2.5;
    # some reach 5 to 15), where the K-correction polynomials diverge. Real red
    # galaxies near z = 0.3 do reach g - r ~ 2, so only those beyond the range
    # are treated: their K-correction uses a typical colour (0.9) and their
    # colour is not used.
    color_ok = (gr > -0.5) & (gr < 2.5)
    gr_k = np.where(color_ok, gr, 0.9)
    kr = cp.k_correction("r", z, gr_k)
    kg = cp.k_correction("g", z, gr_k)
    logm = g["lgm_tot_p50"].astype(float)
    logsfr = g["sfr_tot_p50"].astype(float)
    has_mpa = g["reliable"] > -9000
    mass_ok = has_mpa & (logm > 5.0) & (logm < 13.0)
    # log sSFR > -8 would double a galaxy's stars in under 100 Myr. The 67
    # objects claiming it are mostly star-forming knots of big galaxies that
    # the photometry split off (stellar masses 10^6-10^7), plus one with a
    # claimed 7,500 M_sun/yr -- alone 3.6% of the local SFR density. Dropped.
    sfr_ok = mass_ok & (logsfr > -9.0) & (logsfr < 4.0) & (logsfr - logm < -8.0)
    out = {
        "ra": g["ra"].astype(float), "dec": g["dec"].astype(float), "z": z,
        "r": r, "gr_obs": gr_k, "gr_raw": gr, "color_ok": color_ok,
        "gr0": np.where(color_ok, gr - kg + kr, np.nan),
        "M_r": cp.absolute_magnitude(r, z, kr),
        "d_c": cp.comoving_distance(z),
        "logm": logm, "logm_lo": g["lgm_tot_p16"].astype(float),
        "logm_hi": g["lgm_tot_p84"].astype(float),
        "logsfr": logsfr, "logsfr_lo": g["sfr_tot_p16"].astype(float),
        "logsfr_hi": g["sfr_tot_p84"].astype(float),
        "logsfr_p2p5": g["sfr_tot_p2p5"].astype(float),
        "logsfr_p97p5": g["sfr_tot_p97p5"].astype(float),
        "logsfr_fib": g["sfr_fib_p50"].astype(float),
        "logm_fib": g["lgm_fib_p50"].astype(float),
        "logssfr": logsfr - logm,
        "d4000": g["d4000_n"].astype(float),
        "bptclass_mpa": g["bptclass"].astype(int),
        "reliable": g["reliable"].astype(int),
        "has_mpa": has_mpa, "mass_ok": mass_ok, "sfr_ok": sfr_ok,
        "objID": g["objID"], "specObjID": g["specObjID"],
        "plate": g["plate"], "mjd": g["mjd"], "fiberID": g["fiberID"],
    }
    out.update(bpt_classify(g))
    return out


# ---------------------------------------------------------------------------
# What powers the emission lines: the BPT diagram
# ---------------------------------------------------------------------------
#
# Baldwin, Phillips & Terlevich (1981) noticed that two line ratios separate
# gas lit by young stars from gas lit by an accreting black hole, whose harder
# radiation boosts [O III] and [N II] relative to the hydrogen lines:
#
#   x = log([N II] 6584 / H-alpha),  y = log([O III] 5007 / H-beta).
#
# Each pair of lines is close in wavelength, so dust barely changes the
# ratios. Boundaries:
#   Kauffmann et al. (2003): y = 0.61 / (x - 0.05) + 1.30  -- the empirical
#     upper edge of the star-forming sequence;
#   Kewley et al. (2001):    y = 0.61 / (x - 0.47) + 1.19  -- the most a
#     population of young stars can produce (photoionization models);
#   Schawinski et al. (2007): y = 1.05 x + 0.45 -- Seyferts above, LINERs below.
# Between the first two lie "composites": star formation plus an AGN.

BPT_NONE, BPT_SF, BPT_COMPOSITE, BPT_SEYFERT, BPT_LINER, BPT_WEAK = 0, 1, 2, 3, 4, 5
BPT_NAMES = {BPT_NONE: "no lines", BPT_SF: "star-forming", BPT_COMPOSITE: "composite",
             BPT_SEYFERT: "Seyfert", BPT_LINER: "LINER", BPT_WEAK: "weak lines"}


def kauffmann03(x):
    return 0.61 / (np.asarray(x) - 0.05) + 1.30


def kewley01(x):
    return 0.61 / (np.asarray(x) - 0.47) + 1.19


def schawinski07(x):
    return 1.05 * np.asarray(x) + 0.45


def bpt_classify(g: dict) -> dict:
    """Our own BPT classes from the MPA-JHU line fluxes.

    Flux errors are first scaled up by the factors the MPA-JHU team derived
    from repeat observations (config.LINE_ERROR_SCALE). A galaxy is classified
    on the full diagram if all four lines have S/N >= 3. Otherwise, if H-alpha
    is detected at S/N >= 3 it is "weak lines", and if not, "no lines" --
    in both cases its SFR in MPA-JHU came from the 4000 A break, not from lines.
    """
    def sn(line):
        f = g[f"{line}_flux"].astype(float)
        e = g[f"{line}_flux_err"].astype(float) * LINE_ERROR_SCALE[line]
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where((e > 0) & (f > -9000), f / e, 0.0), f

    sn_ha, ha = sn("h_alpha")
    sn_hb, hb = sn("h_beta")
    sn_o3, o3 = sn("oiii_5007")
    sn_n2, n2 = sn("nii_6584")
    full = (sn_ha >= BPT_MIN_SN) & (sn_hb >= BPT_MIN_SN) & (sn_o3 >= BPT_MIN_SN) \
        & (sn_n2 >= BPT_MIN_SN)
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.where(full, np.log10(n2 / ha), np.nan)
        y = np.where(full, np.log10(o3 / hb), np.nan)
    cls = np.full(len(ha), BPT_NONE, dtype=np.int8)
    cls[sn_ha >= BPT_MIN_SN] = BPT_WEAK
    below_k03 = (x < 0.05) & (y < kauffmann03(x))
    below_k01 = (x < 0.47) & (y < kewley01(x))
    cls[full & below_k03] = BPT_SF
    cls[full & ~below_k03 & below_k01] = BPT_COMPOSITE
    agn = full & ~below_k01
    cls[agn & (y >= schawinski07(x))] = BPT_SEYFERT
    cls[agn & (y < schawinski07(x))] = BPT_LINER
    return {"bpt": cls, "bpt_x": x, "bpt_y": y, "ha_sn": sn_ha}


def is_agn(bpt) -> np.ndarray:
    """Composite, Seyfert or LINER: lines at least partly powered by an AGN."""
    bpt = np.asarray(bpt)
    return (bpt == BPT_COMPOSITE) | (bpt == BPT_SEYFERT) | (bpt == BPT_LINER)


# ---------------------------------------------------------------------------
# The main sequence: Renzini & Peng's ridge line
# ---------------------------------------------------------------------------


def weighted_hist2d(x, y, w, x_edges, y_edges) -> np.ndarray:
    h, _, _ = np.histogram2d(x, y, bins=[x_edges, y_edges], weights=w)
    return h


def ridge_line(logm, logsfr, w, *, bin_dex: float = MS_BIN_DEX,
               logm_range=MS_FIT_LOGM, sfr_floor: float | None = None) -> dict:
    """The star-forming ridge: mode of log SFR in each mass column.

    Following Renzini & Peng (2015): weighted counts in 0.2 x 0.2 dex cells of
    the SFR-mass plane, no pre-selection of star-forming galaxies. In each mass
    column the star-forming peak is the highest local maximum on the
    high-SFR side; its position is refined by a parabola through the peak cell
    and its two neighbours. A straight line fitted to the peaks over
    ``logm_range`` is the main sequence.

    ``sfr_floor``: in the most massive columns the quenched peak can outweigh
    the star-forming one; maxima below this line (a function of mass, applied
    as log SFR > sfr_floor(logm)) are ignored when locating the SF peak.
    """
    x_edges = np.arange(logm_range[0], logm_range[1] + bin_dex / 2, bin_dex)
    y_edges = np.arange(-4.0, 2.6 + bin_dex / 2, bin_dex)
    h = weighted_hist2d(logm, logsfr, w, x_edges, y_edges)
    yc = 0.5 * (y_edges[1:] + y_edges[:-1])
    xc = 0.5 * (x_edges[1:] + x_edges[:-1])
    peaks = np.full(len(xc), np.nan)
    for i, col in enumerate(h):
        ok = np.ones(len(yc), bool) if sfr_floor is None else yc > sfr_floor(xc[i])
        cand = [j for j in range(1, len(yc) - 1)
                if ok[j] and col[j] >= col[j - 1] and col[j] >= col[j + 1] and col[j] > 0]
        if not cand:
            continue
        j = max(cand, key=lambda k: col[k])
        a, b, c = col[j - 1], col[j], col[j + 1]
        denom = a - 2 * b + c
        shift = 0.5 * (a - c) / denom if denom != 0 else 0.0
        peaks[i] = yc[j] + np.clip(shift, -0.5, 0.5) * bin_dex
    good = np.isfinite(peaks)
    coef, cov = np.polyfit(xc[good], peaks[good], 1, cov=True)
    return {"logm": xc, "peak": peaks, "slope": float(coef[0]),
            "intercept": float(coef[1]), "slope_err": float(np.sqrt(cov[0, 0])),
            "intercept_err": float(np.sqrt(cov[1, 1])), "hist": h,
            "x_edges": x_edges, "y_edges": y_edges}


def ms_sfr(logm, ms: dict) -> np.ndarray:
    """log SFR of the main sequence at log M."""
    return ms["slope"] * np.asarray(logm) + ms["intercept"]


def gaussian(x, a, mu, s):
    return a * np.exp(-0.5 * ((x - mu) / s) ** 2)


def sf_peak_width(dms, w, *, fit_above: float = -0.3) -> dict:
    """Width of the star-forming peak: a Gaussian fitted to its upper side.

    The low-SFR side of the peak is contaminated by galaxies on their way to
    quenching, so the Gaussian is fitted to bins above ``fit_above`` dex below
    the ridge, where the peak is clean (Renzini & Peng's Fig. 5 shows the same
    asymmetry).
    """
    edges = np.arange(-3.0, 1.6, 0.05)
    hst, _ = np.histogram(dms, bins=edges, weights=w)
    c = 0.5 * (edges[1:] + edges[:-1])
    sel = c > fit_above
    p, cov = curve_fit(gaussian, c[sel], hst[sel], p0=[hst.max(), 0.0, 0.3])
    return {"amp": float(p[0]), "mu": float(p[1]), "sigma": float(abs(p[2])),
            "sigma_err": float(np.sqrt(cov[2, 2]))}


def valley(dms, w, *, smooth_dex: float = 0.1, sf_side: float = -0.5,
           q_side: float = -1.0) -> dict:
    """The green valley: the minimum between the two peaks in Delta_MS.

    Delta_MS = log SFR - log SFR_MS(M), the distance below the main sequence.
    The weighted histogram is smoothed by a Gaussian of ``smooth_dex``. The
    star-forming peak is its highest local maximum above ``sf_side``; the
    quenched peak its highest local maximum below ``q_side``. Only if both exist
    is there a valley: the minimum between them. Otherwise ``found`` is False
    (a single peak with a tail is not bimodal).
    """
    edges = np.arange(-3.5, 1.5 + 0.025, 0.025)
    hst, _ = np.histogram(dms, bins=edges, weights=w)
    sm = gaussian_filter1d(hst.astype(float), smooth_dex / 0.025)
    c = 0.5 * (edges[1:] + edges[:-1])
    peaks = [j for j in range(1, len(sm) - 1) if sm[j] > sm[j - 1] and sm[j] >= sm[j + 1]]
    sf = [j for j in peaks if c[j] > sf_side]
    qq = [j for j in peaks if c[j] < q_side]
    out = {"centres": c, "smoothed": sm, "found": False, "position": np.nan,
           "depth_vs_sf_peak": np.nan, "sf_peak": np.nan, "q_peak": np.nan}
    if not sf or not qq:
        return out
    js = max(sf, key=lambda j: sm[j])
    jq = max(qq, key=lambda j: sm[j])
    jv = jq + int(np.argmin(sm[jq:js + 1]))
    # a real valley must dip below both peaks by at least 10%
    if sm[jv] > 0.9 * min(sm[js], sm[jq]):
        return out
    out.update(found=True, position=float(c[jv]), depth_vs_sf_peak=float(sm[jv] / sm[js]),
               sf_peak=float(c[js]), q_peak=float(c[jq]), q_over_sf=float(sm[jq] / sm[js]))
    return out


def likelihood_mean_sfr(p2p5, p16, p50, p84, p97p5, *, chunk: int = 50_000) -> np.ndarray:
    """Mean SFR (linear) from five percentiles of the log-SFR likelihood.

    MPA-JHU gives each galaxy's log SFR as a likelihood distribution, quoted by
    its 2.5, 16, 50, 84 and 97.5 percentiles. A sum of SFRs over a volume needs
    the mean of SFR, not the median of log SFR, and the likelihoods are
    skewed. We take log SFR to be piecewise linear in the Gaussian quantile
    z = Phi^-1(q) through the five percentiles (z = -1.96, -1, 0, 1, 1.96),
    extended beyond them with the outer slopes, and integrate 10^x(z) phi(z) dz
    numerically. This is exact when the likelihood is Gaussian in log SFR and
    follows its skewness otherwise. Where a percentile is missing the median
    is returned.
    """
    from scipy.stats import norm

    knots = np.array([-1.959964, -1.0, 0.0, 1.0, 1.959964])
    zg = np.linspace(-6.0, 6.0, 481)
    wz = norm.pdf(zg) * (zg[1] - zg[0])
    x = np.vstack([p2p5, p16, p50, p84, p97p5]).astype(float)
    good = np.all(x > -50, axis=0) & np.all(np.diff(x, axis=0) >= -1e-6, axis=0)
    out = 10 ** np.asarray(p50, float)
    idx = np.flatnonzero(good)
    # interpolation weights on the z grid, shared by every galaxy
    seg = np.clip(np.searchsorted(knots, zg) - 1, 0, 3)
    t = (zg - knots[seg]) / (knots[seg + 1] - knots[seg])  # extrapolates outside
    for lo in range(0, len(idx), chunk):
        ii = idx[lo:lo + chunk]
        xa = x[seg][:, ii]          # (nz, n) lower knot values
        xb = x[seg + 1][:, ii]
        xz = xa + t[:, None] * (xb - xa)
        out[ii] = np.sum(wz[:, None] * 10 ** np.clip(xz, -10, 6), axis=0)
    return out


def weighted_fraction(flag, w) -> tuple[float, float]:
    """Weighted fraction and its binomial error with the effective sample size."""
    flag = np.asarray(flag, bool)
    w = np.asarray(w, float)
    if w.sum() <= 0:
        return np.nan, np.nan
    f = float(np.sum(w * flag) / w.sum())
    n_eff = w.sum() ** 2 / np.sum(w**2)
    return f, float(np.sqrt(max(f * (1 - f), 1e-12) / n_eff))
