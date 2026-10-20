"""Multiplicity against primary mass: wide companions from M dwarfs to B stars.

The cleanest number Gaia gives across the whole main sequence is the
*wide-binary fraction*: the fraction of primaries with a resolved companion
at a projected separation of 100-10,000 au (the quantity Offner et al. 2023
tabulate from imaging surveys, their "WBF"). Within 100 pc such companions sit
1" or more from the primary, where Gaia resolves all but very faint
companions; for the rarer A and B stars we go out to 500 pc and correct for
the companions lost to Gaia's resolution limit with the contrast curve
measured in ``wide.py``.

Two bookkeeping rules make the fraction a property of *systems*:

* a star whose resolved companion is brighter is a secondary, not a primary,
  so it leaves the sample (its system is counted from its brighter partner);
* only companions with mass ratio q > 0.1 count, so that faint M-dwarf and
  brown-dwarf companions, which Gaia sees around nearby stars but not around
  distant ones, do not make nearby primaries look more multiple.
"""

from __future__ import annotations

import numpy as np

import binaries as B
import literature as L
import mocks
import samples
import stellar
import wide

SEP_RANGE_AU = (100.0, 10_000.0)
Q_MIN = 0.1

#: Primary-mass bins (Msun) and the volume used for each (pc): nearby for the
#: common low-mass stars, wider for the rare massive ones.
MASS_BINS = [((0.15, 0.30), 50.0), ((0.30, 0.60), 100.0), ((0.60, 0.80), 100.0), ((0.80, 1.20), 100.0),
             ((1.20, 1.60), 100.0), ((1.60, 2.40), 250.0), ((2.40, 5.0), 500.0)]


def _pairs_for(primary_ids, pairs):
    """Rows of the pair table whose brighter star is one of ``primary_ids``."""
    return pairs[np.isin(np.asarray(pairs["source_id1"]), primary_ids)]


def resolved_completeness(rng, dist_pc, m1, curve, n_draw=4000):
    """Average probability that a companion with s in SEP_RANGE_AU (log-flat,
    i.e. Opik's law) and q from Moe & Di Stefano's distribution is in
    El-Badry et al.'s catalogue, for primaries at these distances."""
    k = rng.integers(0, len(dist_pc), n_draw)
    d, m = np.asarray(dist_pc)[k], np.asarray(m1)[k]
    log_s = rng.uniform(np.log10(SEP_RANGE_AU[0]), np.log10(SEP_RANGE_AU[1]), n_draw)
    s = 10 ** log_s
    # log P for the q draw (a ~ 1.2 s on average; only sets the MDS q-shape)
    logp = np.log10(B.period_days_from_a(1.2 * s, 1.5 * m))
    q = mocks.sample_q(rng, logp, m)
    q = np.maximum(q, Q_MIN)
    dg = stellar.mg_from_mass(q * m) - stellar.mg_from_mass(m)
    theta = s / d  # arcsec
    p = wide.detection_probability(theta, np.where(np.isfinite(dg), dg, 12.0), curve)
    # companions too faint for a parallax with sigma < plx/5 (G > ~20.3) are lost
    g2 = stellar.mg_from_mass(q * m) + 5 * np.log10(d / 10.0)
    p = np.where(np.isfinite(g2) & (g2 < 20.3), p, 0.0)
    return float(np.mean(p))


def _f(col):
    return np.ma.filled(np.ma.asarray(col, dtype=float), np.nan)


def primaries_100pc(max_dist_pc: float):
    """Main-sequence stars within ``max_dist_pc`` that are not the fainter member
    of a resolved pair, with masses from their colours."""
    t = samples.load_100pc()
    keep = (t["dist_pc"] < max_dist_pc) & t["good_astrometry"] & t["ms_band"] & (t["wide_role"] != "secondary")
    return t[keep]


def primaries_bright(max_dist_pc: float):
    """B and A stars within ``max_dist_pc`` (from the 500 pc bright sample).

    Their masses come from Gaia's GSP-Phot temperatures, which fit extinction
    along with the spectrum (BP-RP alone would be reddened by dust at a few
    hundred pc), on Mamajek's temperature-mass sequence.
    """
    t = samples.gaia_table("bright_500pc").copy()
    samples.add_photometry(t)
    teff = _f(t["teff_gspphot"])
    t["mass_teff"] = np.interp(teff, *_teff_mass())
    t["mass_teff"] = np.where(np.isfinite(teff), t["mass_teff"], np.nan)
    pairs = _bright_pairs()
    secondaries = np.asarray(pairs["source_id2"])
    keep = (t["dist_pc"] < max_dist_pc) & ~np.isin(np.asarray(t["source_id"]), secondaries)
    return t[keep]


def _teff_mass():
    tab = stellar.mamajek_table()
    teff, mass = tab["Teff"], tab["Msun"]
    ok = np.isfinite(teff) & np.isfinite(mass)
    order = np.argsort(teff[ok])
    return teff[ok][order], mass[ok][order]


def _bright_pairs():
    from astropy.table import Table
    import config as C
    e = Table.read(C.VIZIER_DIR / "elbadry21_bright.fits")
    e = e[_f(e["R"]) < C.R_CHANCE_MAX]
    out = Table()
    out["source_id1"] = np.asarray(e["Source1"], dtype=np.int64)
    out["source_id2"] = np.asarray(e["Source2"], dtype=np.int64)
    out["sep_au"] = _f(e["sepAU"])
    out["theta_arcsec"] = _f(e["theta"]) * 3600
    out["dg"] = _f(e["Gmag2"]) - _f(e["Gmag1"])
    out["g2"] = _f(e["Gmag2"])
    return out


def wide_binary_fraction(rng, curve, mass_bins=MASS_BINS) -> list[dict]:
    """Observed and completeness-corrected fraction of primaries with a
    resolved companion at 100-10,000 au (q > 0.1), by primary mass."""
    near_pairs = samples.wide_pairs()
    bright_pairs = _bright_pairs()
    out = []
    for (m_lo, m_hi), dmax in mass_bins:
        if dmax <= 100:
            t = primaries_100pc(dmax)
            m1 = np.asarray(t["mass_colour"])
            pairs = near_pairs
            g1 = _f(t["phot_g_mean_mag"])
        else:
            t = primaries_bright(dmax)
            m1 = np.asarray(t["mass_teff"])
            pairs = bright_pairs
            g1 = _f(t["phot_g_mean_mag"])
        sel = (m1 >= m_lo) & (m1 < m_hi)
        t, m1, g1 = t[sel], m1[sel], g1[sel]
        ids = np.asarray(t["source_id"], dtype=np.int64)
        pp = pairs[np.isin(np.asarray(pairs["source_id1"]), ids)]
        in_range = (np.asarray(pp["sep_au"]) >= SEP_RANGE_AU[0]) & (np.asarray(pp["sep_au"]) < SEP_RANGE_AU[1])
        pp = pp[in_range]
        # mass ratio from the brightness difference (same distance, same dust)
        idx = {s: i for i, s in enumerate(ids)}
        rows = np.array([idx[s] for s in np.asarray(pp["source_id1"])], dtype=int)
        mg1 = stellar.mg_from_mass(m1[rows])
        m2 = stellar.mass_from_mg(mg1 + np.asarray(pp["dg"]))
        q = m2 / m1[rows]
        q = np.where(np.isfinite(q), q, 0.0)  # companions fainter than the table (brown dwarfs) have q < 0.1
        has = np.zeros(len(t), bool)
        has[rows[q > Q_MIN]] = True
        n, k = len(t), int(has.sum())
        comp = resolved_completeness(rng, np.asarray(t["dist_pc"]), m1, curve)
        f_obs = k / n if n else np.nan
        err = np.sqrt(f_obs * (1 - f_obs) / n) if n else np.nan
        out.append({"m_lo": m_lo, "m_hi": m_hi, "m_med": float(np.median(m1)) if n else np.nan,
                    "dist_max": dmax, "n_primaries": n, "n_with_companion": k, "wbf_obs": f_obs,
                    "wbf_obs_err": err, "completeness": comp, "wbf": f_obs / comp, "wbf_err": err / comp})
    return out
