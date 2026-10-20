"""How common is each kind of supernova? Luminosity functions, rates, fractions.

A magnitude-limited survey sees bright supernovae out to large distances and
faint ones only nearby. The 1/V_max method undoes this: each supernova counts
as 1/V_max, where V_max is the volume in which a supernova of its peak
luminosity would still have peaked brighter than the survey limit. Divided by
the survey time and the fraction of sky and time the survey really watched,
the sum is a volumetric rate (Schmidt 1968; Perley et al. 2020, their Eq. 2):

    R = 1/T * sum_i 1 / (V_max,i * f_sky * f_ext * f_rec * f_cl,i)

f_sky: the average sky area under active watch; f_ext: the volume lost to
Galactic dust; f_rec: the probability that a supernova peaking inside the
watched area ends up in the quality-cut sample; f_cl: the probability that it
was classified (depends on how bright it got).

Two versions are computed:

* **Perley et al.'s, replicated**: their window (June 2018 - mid July 2020),
  their f-factors, their Euclidean V_max with heliocentric BTS magnitudes.
  This checks our implementation against their published rates.
* **ours**, over 2018-06 to 2025-10: CMB-frame and Cosmicflows-4 distances, a
  cosmological V_max that counts time dilation, f_cl measured from the
  catalogue itself, f_ext from the SFD dust map, f_rec re-measured with their
  recipe for each season, and the survey's down time removed from T.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd

import config as C
from bts import fcl_model, fcl_perley, classification_completeness, statistical_sample
from distances import cosmology, distance_modulus

GPC3_PER_MPC3 = 1e-9


# ---------------------------------------------------------------------------
# V_max
# ---------------------------------------------------------------------------

def vmax_euclidean(M, m_lim: float = C.M_LIM) -> np.ndarray:
    """Perley et al.'s V_max: (4 pi / 3) D_max^3, D_max from the distance
    modulus with no cosmology (Mpc^3, full sky)."""
    d_max = 10 ** (0.2 * (m_lim - np.asarray(M, float)) - 5.0)  # Mpc
    return 4.0 / 3.0 * np.pi * d_max ** 3


@lru_cache(maxsize=1)
def _veff_table():
    """Cumulative time-dilation-weighted comoving volume, full sky:
    V_eff(z) = int_0^z dV_c/dz / (1 + z') dz'  (Mpc^3)."""
    z = np.concatenate([np.geomspace(1e-5, 0.02, 300), np.linspace(0.02, 2.0, 4000)[1:]])
    dvdz = 4 * np.pi * cosmology().differential_comoving_volume(z).value  # Mpc^3 per unit z
    integrand = dvdz / (1 + z)
    veff = np.concatenate([[0.0], np.cumsum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(z))])
    # the apparent magnitude a source of absolute magnitude 0 would have at z
    m0 = distance_modulus(z) - 2.5 * np.log10(1 + z)
    return z, veff, m0


def zmax_cosmological(M, m_lim: float = C.M_LIM) -> np.ndarray:
    """Redshift at which a peak of absolute magnitude M (BTS K-convention)
    reaches m_lim: M + mu(z) - 2.5 log10(1+z) = m_lim."""
    z, _, m0 = _veff_table()
    return np.interp(m_lim - np.asarray(M, float), m0, z, left=np.nan, right=z[-1])


def vmax_cosmological(M, m_lim: float = C.M_LIM) -> np.ndarray:
    """Time-dilation-weighted comoving V_max (Mpc^3, full sky): a supernova at
    redshift z is seen at a (1+z)-times lower rate in observer time."""
    z, veff, _ = _veff_table()
    return np.interp(zmax_cosmological(M, m_lim), z, veff)


# ---------------------------------------------------------------------------
# Survey factors
# ---------------------------------------------------------------------------

def f_rec_perley_recipe(df: pd.DataFrame, start, end, mag_col: str = "M_bts") -> dict:
    """Perley et al.'s recovery fraction: of the SNe Ia brighter than -18.5
    absolute and 18.0 apparent that peaked in the window (easy to find and to
    classify), the fraction that passes the quality cut."""
    t = df["peak_date"]
    win = (t >= pd.Timestamp(start)) & (t <= pd.Timestamp(end) + pd.Timedelta(days=1))
    easy = win & (df["cls"] == "Ia") & (df[mag_col] < C.EASY_IA_MABS) & (df["peakmag"] < C.EASY_IA_MAPP)
    n, k = int(easy.sum()), int((easy & df["quality"]).sum())
    return {"n_easy": n, "n_pass": k, "f_rec": k / n if n else np.nan}


def f_ext_from_dust(nside: int = 128) -> dict:
    """Average fraction of 1/V_max volume left after Galactic dust, over the
    BTS footprint (|b| > 7, dec > -26, A_V < 1), for peaks in r and in g."""
    import healpy as hp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    # LAMBDA's SFD map is in Galactic coordinates; read_map returns RING order.
    ebv = hp.ud_grade(hp.read_map(str(C.DUST_DIR / "lambda_sfd_ebv.fits")), nside)
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix))
    gal = SkyCoord(l=np.degrees(phi) * u.deg, b=(90 - np.degrees(theta)) * u.deg, frame="galactic")
    b, dec = gal.b.deg, gal.icrs.dec.deg
    a_v = 2.742 * ebv   # Schlafly & Finkbeiner (2011) recalibration of SFD, R_V = 3.1
    foot = (np.abs(b) > C.FOOTPRINT_MIN_ABS_B) & (dec > C.FOOTPRINT_MIN_DEC)
    keep = foot & (a_v < C.FOOTPRINT_MAX_AV)
    out = {"area_deg2_footprint": float(foot.sum() * hp.nside2pixarea(nside, degrees=True)),
           "area_deg2_low_extinction": float(keep.sum() * hp.nside2pixarea(nside, degrees=True))}
    for band, k in (("r", C.A_R_PER_AV), ("g", C.A_G_PER_AV)):
        out[f"f_ext_{band}"] = float(np.mean(10 ** (-0.6 * k * a_v[keep])))
    return out


def survey_downtime(df: pd.DataFrame, start, end, min_gap_days: float = 21.0) -> list[tuple]:
    """Stretches with no saved candidate peaking at all (the survey was off).
    Any gap longer than three weeks in the dates of *all* candidates' peaks."""
    t = np.sort(df.loc[(df["peak_date"] >= start) & (df["peak_date"] <= end), "peak_date"].values)
    gaps = []
    for a, b in zip(t[:-1], t[1:]):
        if (b - a) / np.timedelta64(1, "D") > min_gap_days:
            gaps.append((pd.Timestamp(a), pd.Timestamp(b)))
    return gaps


def effective_years(df, start, end) -> tuple[float, list]:
    span = (pd.Timestamp(end) + pd.Timedelta(days=1) - pd.Timestamp(start)).days
    gaps = survey_downtime(df, pd.Timestamp(start), pd.Timestamp(end) + pd.Timedelta(days=1))
    lost = sum((b - a).days for a, b in gaps)
    return (span - lost) / 365.25, gaps


# ---------------------------------------------------------------------------
# Rates
# ---------------------------------------------------------------------------

def rate_weights(sample: pd.DataFrame, *, M: np.ndarray, vmax: np.ndarray, fcl: np.ndarray,
                 T_yr: float, f_sky: float, f_ext: float, f_rec: float) -> np.ndarray:
    """Each supernova's contribution to the volumetric rate (Gpc^-3 yr^-1)."""
    return 1.0 / (T_yr * vmax * GPC3_PER_MPC3 * f_sky * f_ext * f_rec * fcl)


def summed_rate(w: np.ndarray, rng: np.random.Generator, n_boot: int = C.N_BOOT,
                ci: float = 0.95) -> dict:
    """Sum of weights with a bootstrap interval (resampling the supernovae)."""
    w = np.asarray(w, float)
    n = len(w)
    if n == 0:
        return {"rate": 0.0, "lo": 0.0, "hi": np.nan, "n": 0}
    boots = np.array([w[rng.integers(0, n, n)].sum() for _ in range(n_boot)])
    a = (1 - ci) / 2
    return {"rate": float(w.sum()), "lo": float(np.quantile(boots, a)),
            "hi": float(np.quantile(boots, 1 - a)), "n": n,
            "err68": [float(np.quantile(boots, 0.16)), float(np.quantile(boots, 0.84))]}


def perley_replication(df: pd.DataFrame, rng) -> dict:
    """Perley et al. (2020) Section 4.3, with today's classifications."""
    mask = statistical_sample(df, C.PERLEY_PEAK_START, C.PERLEY_PEAK_END)
    s = df[mask & df["cls"].notna()].copy()
    M = s["M_bts"].to_numpy()
    w = rate_weights(s, M=M, vmax=vmax_euclidean(M), fcl=fcl_perley(s["peakmag"]),
                     T_yr=C.PERLEY_T_YR, f_sky=C.PERLEY_F_SKY, f_ext=C.PERLEY_F_EXT,
                     f_rec=C.PERLEY_F_REC)
    s["w"] = w
    out = {"n_sample": int(mask.sum()), "n_classified": int(len(s)),
           "f_rec_today": f_rec_perley_recipe(df, C.PERLEY_PEAK_START, C.PERLEY_PEAK_END)}
    sel = {"Ia": (s["cls"] == "Ia") & (s["M_bts"] < C.RATE_MLIM_IA),
           "CC": s["is_cc"] & (s["M_bts"] < C.RATE_MLIM_CC),
           "SLSN": s["cls"].isin(["SLSN-I", "SLSN-II"]) & (s["M_bts"] < C.RATE_MLIM_SLSN)}
    for key, m in sel.items():
        out[key] = summed_rate(s.loc[m, "w"].to_numpy(), rng)
    return out


def seasons(start=C.STAT_PEAK_START, end=C.STAT_PEAK_END) -> list[tuple]:
    """Twelve-month survey seasons starting in June (the last one truncated)."""
    out = []
    for a in pd.date_range(start, end, freq="12MS"):
        b = min(a + pd.DateOffset(years=1) - pd.Timedelta(days=1), pd.Timestamp(end))
        out.append((a, b))
    return out


def _in(df, a, b) -> pd.Series:
    return (df["peak_date"] >= a) & (df["peak_date"] < b + pd.Timedelta(days=1))


def our_sample(df: pd.DataFrame, start=C.STAT_PEAK_START, end=C.STAT_PEAK_END,
               m_lim: float = C.M_LIM) -> tuple[pd.DataFrame, dict]:
    """The statistical sample with our distances, and survey factors measured
    season by season: classification efficiency vs magnitude, the recovery
    fraction (Perley's recipe) and the time the survey was actually running.
    The exposure that enters the rate is sum_s T_s * f_rec,s."""
    mask = statistical_sample(df, start, end, m_lim)
    s = df[mask & df["cls"].notna() & df["M_peak"].notna()].copy()
    s["fcl"] = np.nan
    s["season"] = ""
    rows = []
    for a, b in seasons(start, end):
        comp = classification_completeness(df, statistical_sample(df, a, b, 19.0))
        sel = _in(s, a, b)
        s.loc[sel, "fcl"] = fcl_model(s.loc[sel, "peakmag"].to_numpy(), comp)
        s.loc[sel, "season"] = f"{a:%Y-%m}"
        T, gaps = effective_years(df, a, b)
        fr = f_rec_perley_recipe(df, a, b, mag_col="M_peak")
        rows.append({"season": f"{a:%Y-%m}", "start": a.date().isoformat(),
                     "end": b.date().isoformat(), "T_eff_yr": T,
                     "downtime": [(x.date().isoformat(), y.date().isoformat()) for x, y in gaps],
                     "f_rec": fr["f_rec"], "n_easy": fr["n_easy"],
                     "f_cl_18p5": float(fcl_model(np.array([18.375]), comp)[0])})
    s["vmax"] = vmax_cosmological(s["M_peak"].to_numpy(), m_lim)
    s["zmax"] = zmax_cosmological(s["M_peak"].to_numpy(), m_lim)
    tab = pd.DataFrame(rows)
    info = {"n_candidates": int(mask.sum()), "n_used": int(len(s)),
            "T_eff_yr": float(tab["T_eff_yr"].sum()),
            "exposure_yr": float((tab["T_eff_yr"] * tab["f_rec"]).sum()),
            "seasons": tab}
    return s, info


def apply_weights(s: pd.DataFrame, *, exposure_yr, f_sky, f_ext) -> pd.DataFrame:
    """Weights for a sample whose recovery-weighted survey time is
    exposure_yr = sum_s T_s f_rec,s (so f_rec is already inside it).
    f_ext is a number, or the dict from f_ext_from_dust (then each supernova
    gets the dust loss of the band its peak was measured in)."""
    s = s.copy()
    if isinstance(f_ext, dict):
        f_ext = np.where(s["peak_filter"] == "g", f_ext["f_ext_g"], f_ext["f_ext_r"])
    s["w"] = rate_weights(s, M=s["M_peak"].to_numpy(), vmax=s["vmax"].to_numpy(),
                          fcl=s["fcl"].to_numpy(), T_yr=exposure_yr, f_sky=f_sky, f_ext=f_ext,
                          f_rec=1.0)
    return s


def rates_by_class(s: pd.DataFrame, rng) -> dict:
    out = {}
    groups = {
        "Ia": (s["cls"] == "Ia") & (s["M_peak"] < C.RATE_MLIM_IA),
        "Ia_all": s["cls"] == "Ia",
        "CC": s["is_cc"] & (s["M_peak"] < C.RATE_MLIM_CC),
        "II": (s["cls"] == "II") & (s["M_peak"] < C.RATE_MLIM_CC),
        "IIb": (s["cls"] == "IIb") & (s["M_peak"] < C.RATE_MLIM_CC),
        "IIn": (s["cls"] == "IIn") & (s["M_peak"] < C.RATE_MLIM_CC),
        "Ibc": (s["cls"] == "Ibc") & (s["M_peak"] < C.RATE_MLIM_CC),
        "IcBL": (s["cls"] == "IcBL") & (s["M_peak"] < C.RATE_MLIM_CC),
        "Ibn/Icn": (s["cls"] == "Ibn/Icn"),
        "SLSN": s["cls"].isin(["SLSN-I", "SLSN-II"]) & (s["M_peak"] < C.RATE_MLIM_SLSN),
        "SLSN-I": (s["cls"] == "SLSN-I") & (s["M_peak"] < C.RATE_MLIM_SLSN),
        "TDE": s["cls"] == "TDE",
    }
    for key, m in groups.items():
        out[key] = summed_rate(s.loc[m, "w"].to_numpy(), rng)
    return out


def rate_ratio(s: pd.DataFrame, num: pd.Series, den: pd.Series, rng,
               n_boot: int = C.N_BOOT) -> dict:
    """Ratio of two rates from the same sample, with a joint bootstrap. The
    survey's sky coverage, running time and recovery cancel: only the
    luminosities, distances and classification efficiencies remain."""
    w = s["w"].to_numpy()
    a, b = num.to_numpy(), den.to_numpy()
    r = w[a].sum() / w[b].sum()
    n = len(w)
    boots = []
    for _ in range(n_boot):
        k = rng.integers(0, n, n)
        boots.append(w[k][a[k]].sum() / w[k][b[k]].sum())
    return {"ratio": float(r), "lo68": float(np.quantile(boots, 0.16)),
            "hi68": float(np.quantile(boots, 0.84)), "lo95": float(np.quantile(boots, 0.025)),
            "hi95": float(np.quantile(boots, 0.975))}


def luminosity_function(s: pd.DataFrame, mask, bins=np.arange(-23.0, -11.99, 0.5),
                        rng=None, n_boot: int = 500) -> pd.DataFrame:
    """Rate per unit magnitude in bins of peak absolute magnitude, with a
    bootstrap 68% band and the raw count per bin."""
    sub = s[mask]
    M, w = sub["M_peak"].to_numpy(), sub["w"].to_numpy()
    idx = np.digitize(M, bins) - 1
    nb = len(bins) - 1
    phi = np.array([w[idx == i].sum() for i in range(nb)]) / np.diff(bins)
    cnt = np.array([(idx == i).sum() for i in range(nb)])
    lo = hi = np.full(nb, np.nan)
    if rng is not None and len(w):
        boots = np.empty((n_boot, nb))
        for b in range(n_boot):
            k = rng.integers(0, len(w), len(w))
            boots[b] = np.bincount(idx[k].clip(0, nb - 1), weights=w[k] * ((idx[k] >= 0) & (idx[k] < nb)),
                                   minlength=nb)[:nb] / np.diff(bins)
        lo, hi = np.quantile(boots, 0.16, axis=0), np.quantile(boots, 0.84, axis=0)
    return pd.DataFrame({"M_lo": bins[:-1], "M_hi": bins[1:], "phi": phi, "phi_lo": lo,
                         "phi_hi": hi, "n": cnt})


def volumetric_fractions(s: pd.DataFrame, rng, M_faint: float = C.RATE_MLIM_CC) -> dict:
    """Ia / Ibc / II shares of all supernovae brighter than M_faint, by volume,
    in Li et al.'s grouping (II includes IIb and IIn; Ibc includes Ic-BL, Ibn)."""
    li_group = {"Ia": ["Ia"], "Ibc": ["Ibc", "IcBL", "Ibn/Icn"], "II": ["II", "IIb", "IIn"]}
    use = s[(s["M_peak"] < M_faint) & s["cls"].isin(sum(li_group.values(), []))]
    lab = np.full(len(use), "", dtype=object)
    for g, members in li_group.items():
        lab[use["cls"].isin(members).to_numpy()] = g
    w = use["w"].to_numpy()
    n = len(w)
    keys = list(li_group)
    frac = {g: float(w[lab == g].sum() / w.sum()) for g in keys}
    boots = {g: [] for g in keys}
    for _ in range(C.N_BOOT):
        k = rng.integers(0, n, n)
        tot = w[k].sum()
        for g in keys:
            boots[g].append(w[k][lab[k] == g].sum() / tot)
    out = {g: {"frac": frac[g], "lo68": float(np.quantile(boots[g], 0.16)),
               "hi68": float(np.quantile(boots[g], 0.84)), "n": int((lab == g).sum())}
           for g in keys}
    # within Type II: II-P/II-L are not split by TNS for most events, so report
    # the IIb and IIn shares of the hydrogen-rich class
    ii = use[use["cls"].isin(li_group["II"])]
    wi = ii["w"].to_numpy()
    for c in ("IIb", "IIn"):
        out[f"{c}_of_II"] = float(wi[ii["cls"].to_numpy() == c].sum() / wi.sum())
    return out


def magnitude_limited_fractions(df: pd.DataFrame, m_grid=np.arange(14.0, 19.01, 0.25),
                                start=C.STAT_PEAK_START, end=C.STAT_PEAK_END) -> pd.DataFrame:
    """Share of each family among classified supernovae peaking brighter than
    m (quality + purity cut; each weighted by 1/f_cl to stand for the
    unclassified ones at its magnitude)."""
    mask19 = statistical_sample(df, start, end, 19.0)
    comp = classification_completeness(df, mask19)
    s = df[mask19 & df["family"].notna()].copy()
    s["wcl"] = 1.0 / fcl_model(s["peakmag"].to_numpy(), comp)
    rows = []
    for m in m_grid:
        sub = s[s["peakmag"] <= m]
        tot = sub["wcl"].sum()
        row = {"m_lim": m, "n": len(sub)}
        for fam in C.FAMILY_ORDER:
            row[fam] = sub.loc[sub["family"] == fam, "wcl"].sum() / tot if tot else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def seasonal_check(s_all: pd.DataFrame, info: dict, rng, f_sky: float, f_ext: float) -> pd.DataFrame:
    """Is the survey stable? The SN Ia and core-collapse rates season by season,
    each season with its own f_cl, f_rec and running time, f_sky fixed. A true
    rate does not change over seven years, so any trend is the survey's."""
    rows = []
    for _, ss in info["seasons"].iterrows():
        sub = s_all[s_all["season"] == ss["season"]]
        sub = apply_weights(sub, exposure_yr=ss["T_eff_yr"] * ss["f_rec"], f_sky=f_sky, f_ext=f_ext)
        ia = summed_rate(sub.loc[(sub["cls"] == "Ia") & (sub["M_peak"] < C.RATE_MLIM_IA), "w"].to_numpy(),
                         rng, n_boot=500)
        cc = summed_rate(sub.loc[sub["is_cc"] & (sub["M_peak"] < C.RATE_MLIM_CC), "w"].to_numpy(),
                         rng, n_boot=500)
        rows.append({**ss.drop(labels=["downtime"]).to_dict(), "n_Ia": ia["n"], "rate_Ia": ia["rate"],
                     "Ia_lo68": ia["err68"][0], "Ia_hi68": ia["err68"][1], "n_CC": cc["n"],
                     "rate_CC": cc["rate"], "CC_lo68": cc["err68"][0], "CC_hi68": cc["err68"][1]})
    return pd.DataFrame(rows)


def malmquist_table(s: pd.DataFrame, zbins=(0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2)) -> pd.DataFrame:
    """Median observed peak absolute magnitude by redshift: a magnitude-limited
    survey sees only the brightest members of each class far away."""
    rows = []
    for fam in ("Ia", "II", "Ibc", "IIn"):
        f = s[s["family"] == fam]
        for a, b in zip(zbins[:-1], zbins[1:]):
            g = f[(f["z_cmb"] >= a) & (f["z_cmb"] < b)]
            if len(g) >= 5:
                rows.append({"family": fam, "z_lo": a, "z_hi": b, "n": len(g),
                             "M_median": float(g["M_peak"].median())})
    return pd.DataFrame(rows)
