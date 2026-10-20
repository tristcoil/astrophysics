"""Published results we compare against, and parsers for the papers' own tables.

Values are typed in from the papers (each with its source) or parsed from the
arXiv source files cached by ``fetch_data.py lit``. Nothing here enters a fit
unless a module says so explicitly, e.g. when we deliberately reproduce a
published fit, or use another paper's eclipse times for epochs we have no
photometry for.
"""

from __future__ import annotations

import re
import tarfile

import numpy as np

import config as C

# ---------------------------------------------------------------------------
# V1309 Sco
# ---------------------------------------------------------------------------

#: Tylenda et al. (2011, A&A 528, A114), Eq. 1 and Section 5.
TYLENDA2011 = {
    "p0": 1.4456, "c": 15.29, "t0": 2455233.5,   # P = P0 exp(c / (t - t0)), t in JD
    "period_drop_percent": 1.2,                   # "decreased by 1.2% during the pre-outburst observations"
    "onset": "March 2008",                        # "the violent phase of the merger ... began in March 2008"
    "local_max": "April 2007",
    "n_pre_discovery": 1340,
    "discovery_jd": 2454712.0,                    # 2008 Sep 2.5 (Nakano 2008)
    "peak_i": 6.8,
}

#: Molnar et al. (2017, ApJ 840, 1), Section 2: for V1309 Sco the onset of the
#: outburst preceded t0 by 1.8 years.
MOLNAR2017_V1309 = {"onset_before_t0_yr": 1.8}

# ---------------------------------------------------------------------------
# KIC 9832227
# ---------------------------------------------------------------------------

#: Molnar et al. (2017), Table 4 and Section 4.2. P(t) = P0 exp(tau / (t - t0)).
MOLNAR2017 = {
    "p0": 0.4579615, "p0_err": 0.0000009, "tau_d": 0.114, "tau_err": 0.014,
    "t0_year": 2022.2, "t0_year_err": 0.7,
    "me15": {"p0": 0.4579603, "tau_d": 0.089, "t0_year": 2020.8, "t0_year_err": 1.5},
    "q": 0.228, "q_err": 0.003, "incl": 53.19, "incl_err": 0.10, "fillout": 0.430, "fillout_err": 0.015,
    "t1": 5800.0, "t2": 5920.0, "t_ratio": 1.021,
    "m1": 1.395, "m2": 0.318, "a_rsun": 2.992, "distance_pc": 565.0,
    "bjd0_mjd": 55688.49913,  # their MJD = JD - 2400000.0
    "third_body": {"semiamp_min": 1.69, "period_d": 590.0},
    "pdot_threshold": 1e-8,
}

#: Socia et al. (2018, ApJL 864, L32): the NSVS eclipse time used by Molnar
#: et al. carried a +0.5-day (MJD convention) error and a one-cycle miscount,
#: which together put it ~1.01 h early. Their reference period for O - C plots.
SOCIA2018 = {"p_ref": 0.4579515, "nsvs_error_hours": 1.01, "nsvs_oc_min": 42.0}

#: Kovacs, Hartman & Bakos (2019, A&A 631, A126): the reference ephemeris of
#: their Table 1 and the singularity date of Molnar's model in their units.
KOVACS2019 = {"epoch_bjd": 2454953.949183, "period": 0.4579515, "molnar_singularity_hjd": 2459651.0}

# ---------------------------------------------------------------------------
# Contact-binary census and period-luminosity relations
# ---------------------------------------------------------------------------

#: Rucinski (2007, MNRAS 382, 393): ASAS contact binaries.
RUCINSKI2007 = {"cutoff_lo": 0.215, "cutoff_hi": 0.22, "peak": 0.27, "record_field": 0.2178,
                "cc_com": 0.2207, "v34_47tuc": 0.2155}

#: Chen, Wang & Deng (2020, ApJS 249, 18): class counts of the ZTF catalogue.
CHEN2020 = {"n_total": 781602, "n_ew": 369707, "n_ea": 49943, "amp_min": 0.08}

#: Chen et al. (2018, ApJ 859, 140), Table 2: M = a log P + b, 183 W UMa stars
#: with TGAS parallaxes, -0.575 < log P < -0.25.
CHEN2018_PL = {
    "W1_max": (-5.67, 0.17, -0.15, 0.08, 0.16), "W1_mean": (-5.83, 0.18, -0.05, 0.08, 0.17),
    "W2_mean": (-5.87, 0.19, -0.04, 0.09, 0.18),
    "G_max": (-8.53, 0.33, 0.05, 0.15, 0.31), "G_mean": (-8.58, 0.35, 0.20, 0.16, 0.32),
    "V_max": (-8.98, 0.38, 0.05, 0.17, 0.35),
    "Ks_mean": (-5.95, 0.21, 0.00, 0.09, 0.19),
    # period-luminosity-colour: M_G(mean) = -5.02 log P + 1.28 (G - W1)_0 - 0.11, sigma 0.16
    "PLC_G_mean": (-5.02, 1.28, -0.11, 0.16),
    "logp_range": (-0.575, -0.25),
}

#: Rucinski & Duerbeck (1997, PASP 109, 1340): M_V = -4.44 log P + 3.02 (B-V)_0 + 0.12.
RUCINSKI_DUERBECK1997 = {"a": -4.44, "b": 3.02, "c": 0.12}

# ---------------------------------------------------------------------------
# T CrB
# ---------------------------------------------------------------------------

#: Schaefer (2023, MNRAS 524, 3146).
SCHAEFER2023 = {
    "peak_jd": {1866: 2402734.5, 1946: 2431860.0},
    "peak_v": 2.0, "t2": 3.0, "t3": 5.0, "t6": 12.0,
    "secondary_peak_after_d": 160.0, "secondary_duration_d": 100.0,
    "high_state": (-10.0, 9.0),             # years around eruption
    "high_state_start": {1946: 1936.0, "next": 2015.0},
    "dip_lead_yr": 1.0,                     # eruption ~1.0 yr after the dip starts
    "after_dip_bottom_months": 3.0,
    "prediction": (2025.5, 1.3),            # 2024.2 - 2026.8
    "template_date": 2023.8,
    "p_orb_fekel": 227.5687, "p_orb_fekel_err": 0.0099, "t_vmax_fekel": 2447918.62,
    "delta_p_1946": (0.185, 0.056), "pdot_1947_2022": (-8.9e-6, 1.6e-6),
    "n_magnitudes": 213730,
}

#: Pei et al. (2026, arXiv:2607.15245): three effective recurrence intervals
#: from the eruption epochs of Schaefer (2024, RNAAS 8, 272) - AD 1217 to 1787
#: counts as seven cycles - and their survival-conditioned probabilities as of
#: 2026 July 11.
PEI2026 = {"intervals_yr": (81.44, 78.39, 79.75), "mean": 79.86, "sd": 1.53, "asof": "2026-07-11",
           "p_rest_of_2026": 0.302, "p_following_12_months": 0.569, "dip_analogue_window": "2026 December",
           "accretion_deficit_earliest": "2029 May",
           "phase_windows": ["2026 Aug 5-8", "2027 Feb 8-10", "2027 Mar 21-23", "2027 Sep 24-25", "2027 Nov 3-6"]}

#: Merc et al. (2025): the 2023 dimming is not the 1945 pre-eruption dip.
MERC2025 = {"verdict": "the 2023 decline deviates from the deep minimum before 1946 and does not "
                       "reliably predict the eruption timing"}

# ---------------------------------------------------------------------------
# Parsers for the arXiv sources
# ---------------------------------------------------------------------------


def _member(aid: str, suffix: str) -> bytes:
    path = C.LIT_DIR / f"arXiv_{aid}.tar.gz"
    with tarfile.open(path) as tf:
        for m in tf.getmembers():
            if m.name.endswith(suffix):
                return tf.extractfile(m).read()
    raise FileNotFoundError(f"{suffix} not in {path}")


def tylenda_periods() -> dict:
    """Tylenda et al. (2011) Fig. 2: their 17 period measurements, read from the
    Supermongo vector file (period.eps) of the arXiv source.

    Axis calibration from the figure's own tick marks: x = JD - 2450000 from
    2200 (device 292) to 4500 (device 2261); y major ticks 1.42, 1.43, 1.44
    at device 670, 1428, 2185. One device unit is 1.17 d in time and
    1.3e-5 d in period. Error bars are the segments drawn from each marker.
    """
    eps = _member("1012.0163", "period.eps").decode("latin-1").splitlines()
    marks = [tuple(map(int, ln.split()[:2])) for ln in eps if re.fullmatch(r"\d+ \d+ M B", ln.strip())]
    segs = []
    for ln in eps:
        m = re.fullmatch(r"(\d+) (\d+) M (\d+) (\d+) L", ln.strip())
        if m:
            segs.append(tuple(map(int, m.groups())))

    def x2jd(x):
        return 2452200.0 + (x - 292) * 2300.0 / (2261 - 292)

    def y2p(y):
        return 1.42 + (y - 670) * 0.01 / (1428 - 670)

    out = {"t": [], "p": [], "t_lo": [], "t_hi": [], "p_lo": [], "p_hi": []}
    for x, y in marks:
        from_mark = [s for s in segs if s[0] == x and s[1] == y]
        xs = [s[2] for s in from_mark if s[3] == y and s[2] != x]
        ys = [s[3] for s in from_mark if s[2] == x and s[3] != y]
        out["t"].append(x2jd(x))
        out["p"].append(y2p(y))
        out["t_lo"].append(x2jd(min(xs)) if xs else np.nan)
        out["t_hi"].append(x2jd(max(xs)) if xs else np.nan)
        out["p_lo"].append(y2p(min(ys)) if ys else np.nan)
        out["p_hi"].append(y2p(max(ys)) if ys else np.nan)
    return {k: np.array(v) for k, v in out.items()}


def socia_eclipse_times() -> list[dict]:
    """Socia et al. (2018) Table 1: primary eclipse times (BJD) with uncertainties."""
    tex = _member("1809.02771", "ms.tex").decode("utf-8", "replace")
    rows = []
    for m in re.finditer(r"\n\s*(\\textit\{Kepler\}|[A-Za-z]+) \(([\d-]+)\)\s*&\s*([\d.]+)\s*&\s*\$\\pm\$([\d.]+)", tex):
        name = m.group(1).replace("\\textit{Kepler}", "Kepler")
        rows.append({"source": name, "years": m.group(2), "bjd": float(m.group(3)), "err": float(m.group(4))})
    return rows


def kovacs_table1() -> list[dict]:
    """Kovacs, Hartman & Bakos (2019) Table 1: 37 primary-minimum times of
    KIC 9832227 from NSVS (1999) to AAVSO (2018), with their O - C values
    (minutes) for five estimators. T_obs is BJD - 2450000."""
    tex = _member("1909.00255", "red_nova_v3_aph.tex").decode("utf-8", "replace")
    body = tex[tex.index("Times of the primary minima"):]
    body = body[:body.index("\\end{tabular}")]
    rows = []
    for ln in body.splitlines():
        if "\\em" not in ln:
            continue
        parts = [p.strip() for p in ln.split("&")]
        nums = [float(re.sub(r"[$\\ ]", "", p)) for p in parts[:13]]
        src = re.search(r"\{\\em ([^}]+)\}", parts[13]).group(1)
        rows.append({"bjd": 2450000.0 + nums[0], "t_mid": 2455000.0 + nums[1], "half_span": nums[2],
                     "oc": nums[3:13:2], "eoc": nums[4:13:2], "source": src})
    return rows


def molnar_v1309_oc() -> dict | None:
    """Placeholder for Molnar et al.'s V1309 Sco O - C figure (not needed: we
    have the OGLE photometry itself)."""
    return None
