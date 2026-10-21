"""Would Gaia DR3 have published a spectroscopic orbit for this binary?

Gaia's spectrometer (RVS) measured each bright star's radial velocity on about
20 transits over 34 months. DR3's single-lined pipeline (Gosset et al. 2025,
A&A 693, A124) only looked at stars that

* had an RV template temperature of 3,875-8,125 K, a G_RVS of 5.5-12 mag and at
  least 10 transits, and
* varied: rv_renormalised_gof > 4 (an F2 statistic of the velocities, N(0,1)
  for constant stars; Katz et al. 2023) with rv_chisq_pvalue <= 0.01;

and it published an orbit only if, among other cuts, the semi-amplitude was
significant (K/sigma_K >= 5), the orbit's phases were covered (no gap > 0.3),
the period was found with high confidence, and the period was below 1.5 times
the time span (longer ones became "trend" solutions).

We simulate exactly those steps for a mock binary at a real star: its RVS
epochs drawn from Gaia's scan times at its position (gaiamock's scan-law
files), its own number of transits, and the per-transit precision of constant
stars of its brightness. Two cuts are not modelled -- the fit's goodness of
fit (F2 <= 3) and the "efficiency" >= 0.1 of its parameter correlations -- so
the model is calibrated on SB9 binaries of known orbit (``validate_on_sb9``).

What the spectrum shows depends on the companion's light. If the two stars'
lines overlap (|v1 - v2| below ~1.5 resolution elements), a single-lined
analysis measures the light-weighted mean velocity, whose amplitude is
K1 (1 - beta/q)/(1 + beta) for an RVS-band flux ratio beta: zero for twins.
When the lines separate and the companion is bright enough (beta > 0.3),
Gaia's double-lined (SB2) pipeline can catch it instead.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
from astropy.table import Table

import binaries as B
import config as C

TEFF_RANGE = (3875.0, 8125.0)
GRVS_RANGE = (5.5, 12.0)
N_MIN = 10
GOF_MIN = 4.0
SIGNIFICANCE_MIN = 5.0
PHASE_GAP_MAX = 0.3
ERROR_FLOOR = 0.113  # km/s (Katz et al. 2023, Eq. 1)
F_MAX_PER_DAY = 4.0  # highest frequency searched (Gosset et al. 2025)
#: Line blending: the RVS resolution element is ~26 km/s (R = 11,500); we take
#: two lines as separable when they are > 40 km/s apart.
DV_RESOLVED = 40.0
SB2_FLUX_MIN = 0.3
SB2_EPOCHS_MIN = 3
DR3_JD = (2456891.5, 2457902.0)  # window used by gaiamock for DR3 scan times


def epoch_sigma_from_error(rv_error, n_transits):
    """Per-transit RV scatter from DR3's radial_velocity_error (Katz et al.
    2023, Eq. 1: eps^2 = (pi/(2N)) sigma^2 + 0.113^2)."""
    eps = np.asarray(rv_error, dtype=float)
    n = np.asarray(n_transits, dtype=float)
    return np.sqrt(2 * n / np.pi) * np.sqrt(np.clip(eps ** 2 - ERROR_FLOOR ** 2, 0, None))


class NoiseModel:
    """Median per-transit scatter of RV-constant stars, as a function of G_RVS
    (binned by 0.25 mag), measured from a sample's own DR3 data."""

    def __init__(self, grvs, rv_error, n_transits, gof, pvalue, teff):
        grvs = np.asarray(grvs, dtype=float)
        sig = epoch_sigma_from_error(rv_error, n_transits)
        const = ((np.asarray(gof) < 2) & (np.asarray(pvalue) > 0.01) & (np.asarray(n_transits) >= N_MIN)
                 & np.isfinite(sig) & np.isfinite(grvs) & (np.asarray(teff) > TEFF_RANGE[0])
                 & (np.asarray(teff) < TEFF_RANGE[1]))
        self.edges = np.arange(2.0, 12.51, 0.25)
        c = 0.5 * (self.edges[1:] + self.edges[:-1])
        med = np.full(len(c), np.nan)
        for k in range(len(c)):
            m = const & (grvs >= self.edges[k]) & (grvs < self.edges[k + 1])
            if m.sum() >= 20:
                med[k] = np.median(sig[m])
        ok = np.isfinite(med)
        self.grvs, self.sigma = c[ok], med[ok]

    def __call__(self, grvs):
        return np.interp(np.asarray(grvs, dtype=float), self.grvs, self.sigma)


@lru_cache(maxsize=60000)
def _scan_times(pix: int) -> np.ndarray:
    t = Table.read(C.GAIAMOCK_DIR / "healpix_scans" / f"healpix_64_{pix}.fits")
    jd = np.asarray(t["ObservationTimeAtBarycentre[BarycentricJulianDateInTCB]"], dtype=float)
    return np.sort(jd[(jd > DR3_JD[0]) & (jd < DR3_JD[1])])


def rvs_epochs(ra_deg, dec_deg, n_transits, rng) -> np.ndarray:
    """RVS epochs (days, barycentric JD - 2457389, i.e. from J2016.0): a random
    subset of the field-of-view transits at this position. RVS covers about
    half of the focal plane's rows, so not every transit yields a spectrum."""
    import healpy as hp
    pix = int(hp.ang2pix(64, np.radians(90.0 - dec_deg), np.radians(ra_deg), nest=False))
    jd = _scan_times(pix)
    n = int(min(max(n_transits, 0), len(jd)))
    pick = np.sort(rng.choice(len(jd), size=n, replace=False)) if n > 0 else np.array([], int)
    return jd[pick] - 2457389.0


def f2_statistic(chi2, nu):
    """Wilson-Hilferty transform of a chi-square to an approximately N(0,1)
    variable (the 'F2' of Gaia's goodness-of-fit columns)."""
    return np.sqrt(9 * nu / 2) * ((chi2 / nu) ** (1 / 3) + 2 / (9 * nu) - 1)


def _shape(t, period, t_peri, ecc, omega):
    """Keplerian RV curve for K = 1."""
    return B.radial_velocity(t, period, t_peri, ecc, omega, 1.0)


def max_phase_gap(t, period) -> float:
    ph = np.sort(np.mod(t / period, 1.0))
    if len(ph) < 2:
        return 1.0
    gaps = np.diff(np.concatenate([ph, [ph[0] + 1.0]]))
    return float(gaps.max())


def period_confidence_threshold(period):
    """Gosset et al. 2025: 0.95 above 10 d, 0.995 below 1 d, log-linear between."""
    p = float(period)
    if p >= 10:
        return 0.95
    if p <= 1:
        return 0.995
    return 0.995 - 0.045 * np.log10(p)


def simulate(times, sigma, period, t_peri, ecc, omega, k1, q, beta_rvs, rng) -> str:
    """One mock binary through the DR3 spectroscopic selection.

    Returns 'orbit' (SB1/SB2 orbit published), 'trend' (variable, period too
    long: a trend solution), 'variable' (flagged variable, no solution) or ''.
    """
    n = len(times)
    if n < N_MIN or not np.isfinite(sigma) or sigma <= 0:
        return ""
    shape = _shape(times, period, t_peri, ecc, omega)
    k2 = k1 / q
    sep = np.abs(k1 * shape + k2 * shape)  # |v1 - v2|
    resolved = sep > DV_RESOLVED
    # SB2 route: companion bright in the RVS band and the lines apart often enough
    sb2 = (beta_rvs >= SB2_FLUX_MIN) and (resolved.sum() >= SB2_EPOCHS_MIN)
    # single-lined velocity: light-weighted mean where lines blend, primary where they separate
    k_blend = B.light_weighted_k(k1, q, beta_rvs)
    v = np.where(resolved & (beta_rvs < 1.0), k1 * shape, k_blend * shape)
    v = v + rng.normal(0.0, sigma, n)
    chi2 = np.sum((v - v.mean()) ** 2) / sigma ** 2
    if f2_statistic(chi2, n - 1) <= GOF_MIN and not sb2:
        return ""
    span = times.max() - times.min()
    if period > 1.5 * span:
        return "trend"
    if period < 1.0 / F_MAX_PER_DAY:
        return "variable"
    # significance of the semi-amplitude at the true period (linear in gamma, K)
    x = np.column_stack([np.ones(n), shape])
    cov = np.linalg.inv(x.T @ x) * sigma ** 2
    amp = k1 if (sb2 or np.mean(resolved) > 0.5) else k_blend
    signif = amp / np.sqrt(cov[1, 1])
    if signif < SIGNIFICANCE_MIN:
        return "variable"
    if max_phase_gap(times, period) > PHASE_GAP_MAX:
        return "variable"
    # period confidence: false-alarm probability of the best peak among
    # ~f_max * span independent frequencies, for a chi-square drop of signif^2
    n_freq = max(F_MAX_PER_DAY * span, 1.0)
    fap = min(1.0, n_freq * np.exp(-0.5 * signif ** 2))
    if 1.0 - fap < period_confidence_threshold(period):
        return "variable"
    return "orbit_sb2" if sb2 else "orbit"


# ---------------------------------------------------------------------------
# Calibration on SB9: binaries whose orbits we already know
# ---------------------------------------------------------------------------

def sb9_with_gaia(sb9_table) -> Table:
    """SB9 systems inside DR3's SB1 input box, with their DR3 RV statistics."""
    x = Table.read(C.VIZIER_DIR / "sb9_gaia_xmatch.fits")
    f = lambda c: np.ma.filled(np.ma.asarray(x[c], dtype=float), np.nan)
    src = np.asarray(x["Source"]).astype(np.int64)
    info = {int(s): i for i, s in enumerate(src)}
    t = sb9_table
    rows = np.array([info.get(int(s), -1) for s in np.asarray(t["gaia_source_id"])])
    has = rows >= 0
    out = t[has]
    r = rows[has]
    out["grvs"] = f("GRVSmag")[r]
    out["teff_rv"] = f("RVtempTeff")[r]
    out["n_rv"] = f("o_RV")[r]
    out["rv_err"] = f("e_RV")[r]
    out["rv_gof"] = f("RVgof")[r]
    out["rv_pval"] = f("RVchi2")[r]
    out["ra"] = f("RAdeg")[r]
    out["dec"] = f("DEdeg")[r]
    return out


def validate_on_sb9(sb9_table, orbits_table, noise: NoiseModel, *, n_sim: int = 40, seed: int = 7) -> dict:
    """Predicted vs actual DR3 spectroscopic orbits for SB9 binaries.

    For each SB9 orbit (P, e, omega, K1, and K2 for double-lined systems) we
    simulate Gaia's RVS epochs at the star's position many times, with random
    orbital phase, and record the fraction that DR3's selection would publish.
    """
    import sb9 as sb9mod
    rng = np.random.default_rng(seed)
    s = sb9_with_gaia(sb9_table)
    teff, grvs, n = np.asarray(s["teff_rv"]), np.asarray(s["grvs"]), np.asarray(s["n_rv"])
    inbox = ((teff > TEFF_RANGE[0]) & (teff < TEFF_RANGE[1]) & (grvs >= GRVS_RANGE[0]) & (grvs <= GRVS_RANGE[1])
             & (n >= N_MIN) & np.isfinite(np.asarray(s["period"])) & np.isfinite(np.asarray(s["k1"])))
    s = s[inbox]
    types = np.char.strip(np.asarray(orbits_table["nss_solution_type"]).astype(str))
    spec_ids = set(np.asarray(orbits_table["source_id"])[np.isin(types, ("SB1", "SB1C", "SB2", "SB2C",
                                                                         "AstroSpectroSB1", "EclipsingSpectro"))]
                   .tolist())
    pred, found, sig = [], [], []
    for row in s:
        p, e = float(row["period"]), float(row["ecc"]) if np.isfinite(row["ecc"]) else 0.0
        k1 = float(row["k1"])
        k2 = float(row["k2"]) if np.isfinite(row["k2"]) else np.nan
        q = k1 / k2 if np.isfinite(k2) and k2 > 0 else 0.3
        beta = 0.0
        if np.isfinite(k2):  # double-lined: companion visible; light ratio from the mass ratio (L ~ M^4)
            beta = min(q ** 4, 1.0)
        sigma = float(noise(row["grvs"]))
        hits = 0
        for _ in range(n_sim):
            t = rvs_epochs(float(row["ra"]), float(row["dec"]), int(row["n_rv"]), rng)
            res = simulate(t, sigma, p, rng.uniform(0, p), min(e, 0.95), rng.uniform(0, 2 * np.pi), k1, q, beta, rng)
            hits += res.startswith("orbit")
        pred.append(hits / n_sim)
        found.append(int(row["gaia_source_id"]) in spec_ids)
        sig.append(k1 / sigma)
    s["p_model"] = np.array(pred)
    s["found"] = np.array(found)
    s["k_over_sigma"] = np.array(sig)
    return {"table": s, "n": len(s), "predicted": float(np.sum(pred)), "found": int(np.sum(found))}


#: Bins of the empirical acceptance efficiency (DR3 published / modelled).
EFF_LOGP_EDGES = np.array([-1.0, 0.5, 1.0, 2.0, 3.0, 3.5])
EFF_E_EDGES = np.array([0.0, 0.5, 1.0])


def efficiency_from_validation(val: dict, dwarf_mask) -> dict:
    """DR3's orbit acceptance efficiency, calibrated on SB9.

    For single-lined SB9 dwarfs, the ratio (number with a DR3 orbit) /
    (number our idealised selection predicts) in bins of log P and e; for
    double-lined systems a single overall ratio (there are too few to bin).
    These ratios absorb what the model leaves out: the fit-quality and
    'efficiency' cuts, line broadening in tidally locked short-period stars,
    and the troubles of double-lined spectra.
    """
    t = val["table"]
    lp = np.log10(np.asarray(t["period"], dtype=float))
    ecc = np.nan_to_num(np.asarray(t["ecc"], dtype=float))
    dl = np.asarray(t["double_lined"], dtype=bool)
    found = np.asarray(t["found"], dtype=float)
    pm = np.asarray(t["p_model"], dtype=float)
    dw = np.asarray(dwarf_mask, dtype=bool)
    ne, nl = len(EFF_E_EDGES) - 1, len(EFF_LOGP_EDGES) - 1
    eff = np.full((nl, ne), np.nan)
    err = np.full((nl, ne), np.nan)
    nobj = np.zeros((nl, ne), int)
    for i in range(nl):
        for j in range(ne):
            m = (~dl & dw & (lp >= EFF_LOGP_EDGES[i]) & (lp < EFF_LOGP_EDGES[i + 1])
                 & (ecc >= EFF_E_EDGES[j]) & (ecc < EFF_E_EDGES[j + 1]))
            nobj[i, j] = int(m.sum())
            if pm[m].sum() >= 3:
                eff[i, j] = min(found[m].sum() / pm[m].sum(), 1.0)
                err[i, j] = np.sqrt(max(found[m].sum(), 1.0)) / pm[m].sum()
    # fill empty cells from the same period bin's other eccentricity range
    for i in range(nl):
        if np.isnan(eff[i]).all():
            eff[i] = np.nanmean(eff)
        eff[i] = np.where(np.isnan(eff[i]), np.nanmean(eff[i]), eff[i])
    msb2 = dl & dw & (lp < 3)
    eff_sb2 = float(found[msb2].sum() / max(pm[msb2].sum(), 1e-9))
    return {"eff": eff, "err": err, "n": nobj, "logp_edges": EFF_LOGP_EDGES, "e_edges": EFF_E_EDGES,
            "eff_sb2": min(eff_sb2, 1.0), "n_sb2": int(msb2.sum()),
            "n_single_dwarfs": int((~dl & dw & (lp < 3.5)).sum())}


def efficiency(logp, ecc, eff_table: dict, sb2: bool = False):
    if sb2:
        return eff_table["eff_sb2"]
    i = int(np.clip(np.searchsorted(eff_table["logp_edges"], logp, side="right") - 1, 0,
                    len(eff_table["logp_edges"]) - 2))
    j = int(np.clip(np.searchsorted(eff_table["e_edges"], ecc, side="right") - 1, 0,
                    len(eff_table["e_edges"]) - 2))
    return float(eff_table["eff"][i, j])
