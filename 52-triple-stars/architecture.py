"""How triples are built: periods, stability, alignment, and who has a third star.

* The MSC's inner and outer periods against the Mardling-Aarseth stability
  limit (which real triples sit closest to it?).
* Orbit alignment: Tokovinin (2017) reproduced from his own tables (sign
  correlation C, its mass and separation dependence, and the two-orbit
  sample), then repeated on the 2026 MSC.
* Tertiary companions of close binaries: Tokovinin et al. (2006) raw
  fractions reproduced from their table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import catalogs as K
import config as C

# ---------------------------------------------------------------------------
# Periods and stability
# ---------------------------------------------------------------------------


def msc_pairs():
    comp, sys_, orb, notes = K.read_msc()
    pairs = K.inner_outer_pairs(sys_, orb, comp)
    return comp, sys_, orb, notes, pairs


def stability_table(pairs: pd.DataFrame) -> pd.DataFrame:
    """Add the Mardling-Aarseth critical period ratio to every couple with
    known masses (e_out = 0 where no outer orbit is known; coplanar)."""
    p = pairs.copy()
    ok = np.isfinite(p["P_in"]) & np.isfinite(p["P_out"]) & np.isfinite(p["m_in"]) & np.isfinite(p["m3"])
    p = p[ok].copy()
    e_out = p["e_out"].fillna(0.0).clip(0, 0.95)
    p["e_out_used"] = e_out
    p["ratio_crit"] = C.ma01_period_ratio(p["m_in"], p["m3"], e_out)
    p["margin"] = p["ratio"] / p["ratio_crit"]
    p["both_orbits"] = p["inner_from_orbit"] & p["outer_from_orbit"]
    return p


def stability_summary(p: pd.DataFrame) -> dict:
    both = p[p["both_orbits"]]
    viol = both[both["margin"] < 1]
    compact = both[both["P_out"] < 1000].sort_values("P_out")
    out = {
        "n_couples": int(len(p)),
        "n_both_from_orbits": int(len(both)),
        "n_both_orbits_below_ma01": int(len(viol)),
        "below_ma01": viol[["wds", "inner", "outer", "P_in", "P_out", "ratio", "ratio_crit", "e_out", "ids"]]
        .round(3).to_dict("records"),
        "median_log_ratio_all": float(np.median(np.log10(p["ratio"]))),
        "median_log_ratio_both_orbits": float(np.median(np.log10(both["ratio"]))),
        "frac_ratio_lt_10_both_orbits": float(np.mean(both["ratio"] < 10)),
        "shortest_outer_periods": compact.head(8)[["wds", "inner", "outer", "P_in", "P_out", "ratio",
                                                    "ratio_crit", "ids"]].round(3).to_dict("records"),
        "tightest_both_orbits": both.sort_values("margin").head(8)[
            ["wds", "inner", "outer", "P_in", "P_out", "ratio", "ratio_crit", "margin", "e_out", "ids"]]
        .round(3).to_dict("records"),
        "min_margin_both_orbits": float(both["margin"].min()),
    }
    return out


# ---------------------------------------------------------------------------
# Orbit alignment: Tokovinin (2017)
# ---------------------------------------------------------------------------

def sign_correlation(signs) -> tuple[float, float, int]:
    """C = (n+ - n-)/(n+ + n-) with its binomial error (Tokovinin 2017 eq. 1)."""
    s = np.asarray(signs)
    s = s[s != 0]
    n = len(s)
    if n == 0:
        return np.nan, np.nan, 0
    c = (np.sum(s > 0) - np.sum(s < 0)) / n
    return float(c), float(np.sqrt((1 + c) * (1 - c) / n)), n


def mean_phi_from_c(c: float) -> float:
    """<Phi> = (pi/2)(1 - C) for random viewing (Worley 1967), in degrees."""
    return float(np.degrees(np.pi / 2 * (1 - c)))


def tokovinin2017_reproduction() -> dict:
    t1, t3 = K.read_tokovinin2017()
    known = t1[t1["revdir"] != 0].copy()
    # recompute the sign from the inner inclination and the outer direction
    inner_sense = np.where(known["i1"] < 90, 1, -1)
    known["sign_ours"] = inner_sense * known["revdir"]
    agree = float(np.mean(known["sign_ours"] == known["sign"]))
    c_all, c_err, n = sign_correlation(known["sign"])
    known["s_au"] = known["rho"] / (known["plx"] / 1000.0)
    by_mass = []
    for lo, hi, lab in [(0, 1, "<1"), (1, 2, "1 to 2"), (2, 99, ">2"), (0, 99, "All")]:
        g = known[(known["M1"] >= lo) & (known["M1"] < hi)] if lab != "All" else known
        c, ce, nn = sign_correlation(g["sign"])
        e_plus = g.loc[g["sign"] > 0, "e1"].mean()
        e_minus = g.loc[g["sign"] < 0, "e1"].mean()
        by_mass.append(dict(mass=lab, N=nn, median_M1=float(g["M1"].median()), C=c, C_err=ce,
                            median_s_au=float(g["s_au"].median()), e_plus=float(e_plus),
                            e_minus=float(e_minus)))
    # running mean in groups of 24 sorted by separation, and the linear fit
    ks = known.sort_values("s_au").reset_index(drop=True)
    run = []
    for k in range(0, len(ks) - 23):
        g = ks.iloc[k:k + 24]
        c, ce, _ = sign_correlation(g["sign"])
        run.append((float(np.median(g["s_au"])), c, ce))
    run = np.array(run)
    # least squares on the raw signs: C(s) = a + b log s
    A = np.vstack([np.ones(len(ks)), np.log10(ks["s_au"])]).T
    coef, *_ = np.linalg.lstsq(A, ks["sign"].astype(float), rcond=None)
    tight = ks[ks["s_au"] < 50]
    c_tight, ce_tight, n_tight = sign_correlation(tight["sign"])
    wide = ks[ks["s_au"] > 1000]
    c_wide, ce_wide, n_wide = sign_correlation(wide["sign"])
    # Table 3: systems with two visual orbits. Rows come in (I, O) pairs; three
    # systems list two inner orbits, so pair the rows in order, not by WDS.
    t3 = t3.reset_index(drop=True)
    assert (t3["io"].iloc[0::2] == "I").all() and (t3["io"].iloc[1::2] == "O").all()
    t3i = t3.iloc[0::2].reset_index(drop=True)
    t3o = t3.iloc[1::2].reset_index(drop=True)
    # The VizieR file has one defect: the I row of 20002-5522 is followed by
    # the O row of 20181+4044 (the two rows in between are missing). Dropping
    # that mismatched pair leaves the paper's 54 triples.
    same = (t3i["wds"] == t3o["wds"]).to_numpy()
    mismatched = [f"{a} I / {b} O" for a, b in zip(t3i["wds"][~same], t3o["wds"][~same])]
    both = t3i[same].join(t3o[same].drop(columns="wds"), lsuffix="_in", rsuffix="_out")
    phi1, phi2 = K.mutual_inclinations(both["i_in"], both["Omega_in"], both["i_out"], both["Omega_out"])
    both["phi1_ours"], both["phi2_ours"] = phi1, phi2
    dphi = np.minimum(np.abs(phi1 - both["phi_in"]) + np.abs(phi2 - both["phi_out"]),
                      np.abs(phi2 - both["phi_in"]) + np.abs(phi1 - both["phi_out"]))
    # one row (23393+4543) prints the inner inclination in both Phi columns
    phi_bad = both.loc[dphi > 1.0, "wds"].tolist()
    both = both.assign(phi_in=np.where(dphi > 1.0, phi1, both["phi_in"]),
                       phi_out=np.where(dphi > 1.0, phi2, both["phi_out"]))
    sense_in = np.where(both["i_in"] < 90, 1, -1)
    sense_out = np.where(both["i_out"] < 90, 1, -1)
    c3, c3e, n3 = sign_correlation(sense_in * sense_out)
    return dict(
        n_table1=int(len(t1)), n_known_direction=int(n), sign_recomputed_agreement=agree,
        C_all=c_all, C_all_err=c_err, mean_phi_all=mean_phi_from_c(c_all),
        by_mass=by_mass, fit_a=float(coef[0]), fit_b=float(coef[1]),
        C_s_lt_50=c_tight, C_s_lt_50_err=ce_tight, N_s_lt_50=n_tight,
        mean_phi_s_lt_50=mean_phi_from_c(c_tight),
        C_s_gt_1000=c_wide, C_s_gt_1000_err=ce_wide, N_s_gt_1000=n_wide,
        running=run.tolist(),
        n_two_orbits=int(len(both)), table3_mismatched=mismatched, C_two_orbits=c3, C_two_orbits_err=c3e,
        mean_phi_two_orbits=mean_phi_from_c(c3), phi_n_match=int(np.sum(dphi <= 1.0)),
        phi_table_errors=phi_bad,
        two_orbit_phis=np.vstack([both["phi_in"], both["phi_out"]]).T.tolist(),
    )


# ---------------------------------------------------------------------------
# The "f aligned within Phi0, the rest random" model (Tokovinin 2017 sec. 4)
# ---------------------------------------------------------------------------

def _random_unit(n, rng):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1)[:, None]


def _tilt(n1, phi, rng):
    """Unit vectors at angle phi from n1, random azimuth."""
    a = np.cross(n1, np.array([0.0, 0.0, 1.0]))
    bad = np.linalg.norm(a, axis=1) < 1e-8
    a[bad] = np.cross(n1[bad], np.array([1.0, 0.0, 0.0]))
    a /= np.linalg.norm(a, axis=1)[:, None]
    b = np.cross(n1, a)
    psi = rng.uniform(0, 2 * np.pi, len(n1))
    perp = a * np.cos(psi)[:, None] + b * np.sin(psi)[:, None]
    return n1 * np.cos(phi)[:, None] + perp * np.sin(phi)[:, None]


def _visual_elements(n):
    """Inclination and node (mod 180) a visual orbit gives for spin vector n,
    the observer looking along -z."""
    i = np.degrees(np.arccos(np.clip(n[:, 2], -1, 1)))
    node = np.degrees(np.arctan2(n[:, 0], -n[:, 1])) % 180.0   # node line is z x n
    return i, node


def simulate_phi_pairs(f: float, phi0: float, n: int, rng) -> np.ndarray:
    """(Phi_a, Phi_b) pairs an observer would derive from two visual orbits
    (nodes known only mod 180) for the model: fraction f with true Phi uniform
    in (0, phi0), the rest isotropic. Returned sorted (min, max)."""
    n1 = _random_unit(n, rng)
    aligned = rng.uniform(size=n) < f
    cosphi = np.where(aligned, np.cos(np.radians(rng.uniform(0, phi0, n))), rng.uniform(-1, 1, n))
    n2 = _tilt(n1, np.arccos(cosphi), rng)
    i1, O1 = _visual_elements(n1)
    i2, O2 = _visual_elements(n2)
    p1, p2 = K.mutual_inclinations(i1, O1, i2, O2)
    return np.sort(np.vstack([p1, p2]).T, axis=1)


def fit_alignment_model(phi_pairs: np.ndarray, *, nsim=200_000, seed=C.RANDOM_SEED,
                        f_grid=None, phi0_grid=None, bins=12) -> dict:
    """Likelihood of the observed unordered (Phi1, Phi2) pairs on a grid of
    (f, Phi0), using 2-D histograms of simulated pairs (min, max)."""
    rng = np.random.default_rng(seed)
    f_grid = np.linspace(0.0, 1.0, 21) if f_grid is None else f_grid
    phi0_grid = np.arange(10, 181, 10.0) if phi0_grid is None else phi0_grid
    obs = np.sort(np.asarray(phi_pairs, float), axis=1)
    edges = np.linspace(0, 180, bins + 1)
    io = np.clip(np.digitize(obs[:, 0], edges) - 1, 0, bins - 1)
    jo = np.clip(np.digitize(obs[:, 1], edges) - 1, 0, bins - 1)
    lnl = np.full((len(f_grid), len(phi0_grid)), np.nan)
    for b, phi0 in enumerate(phi0_grid):
        aligned = simulate_phi_pairs(1.0, phi0, nsim, rng)
        rand = simulate_phi_pairs(0.0, 180.0, nsim, rng)
        ha, _, _ = np.histogram2d(aligned[:, 0], aligned[:, 1], bins=[edges, edges])
        hr, _, _ = np.histogram2d(rand[:, 0], rand[:, 1], bins=[edges, edges])
        ha /= ha.sum()
        hr /= hr.sum()
        for a, f in enumerate(f_grid):
            prob = f * ha + (1 - f) * hr + 1e-6
            lnl[a, b] = float(np.sum(np.log(prob[io, jo])))
    k = np.unravel_index(np.nanargmax(lnl), lnl.shape)
    post = np.exp(lnl - np.nanmax(lnl))
    post /= post.sum()
    f_mean = float(np.sum(post.sum(axis=1) * f_grid))
    lnl_t17 = float(lnl[np.argmin(np.abs(f_grid - 0.8)), np.argmin(np.abs(phi0_grid - 70))])
    lnl_random = float(lnl[np.argmin(np.abs(f_grid - 0.0)), 0])
    return dict(f_best=float(f_grid[k[0]]), phi0_best=float(phi0_grid[k[1]]), f_mean=f_mean,
                lnL_best=float(np.nanmax(lnl)), lnL_tokovinin=lnl_t17, lnL_random=lnl_random,
                f_grid=f_grid.tolist(), phi0_grid=phi0_grid.tolist(), lnL=lnl.tolist(), n=int(len(obs)))


def msc2026_two_orbits(pairs: pd.DataFrame, orb: pd.DataFrame, p_out_max_yr: float = 1000.0) -> pd.DataFrame:
    """Couples with two visual orbits in the 2026 MSC, with Tokovinin's cuts:
    outer period below 1000 yr, both nodes and inclinations given."""
    v = K.visual_orbit_pairs(pairs, orb)
    v = v[(v["P_out_yr"] < p_out_max_yr) & (v["Omega_in"] > 0) & (v["Omega_out"] > 0)]
    v = v.drop_duplicates(["wds", "inner", "outer"])
    v["s_au"] = v["a_out"] / (v["plx"] / 1000.0)
    return v.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Tertiary companions of close binaries: Tokovinin et al. (2006)
# ---------------------------------------------------------------------------

def tokovinin2006_reproduction() -> dict:
    t2, t5 = K.read_tokovinin2006()
    sb = t5[t5["type"].str.contains("S") & (t5["P_unit"] == "d") & t5["P"].between(0.9, 31)].copy()
    nsub = t5.groupby("hip").size()
    sb["has_tertiary"] = sb["hip"].map(nsub) > 1
    # uncertain companions: every other subsystem of the system flagged '?'
    other = t5.merge(sb[["hip", "comp"]], on="hip", suffixes=("", "_sb"))
    other = other[other["comp"] != other["comp_sb"]]
    certain = other.groupby("hip").apply(
        lambda g: bool(((g["m2_code"] != "?") & (~g["remark"].str.startswith("?"))).any()))
    sb["certain_tertiary"] = sb["hip"].map(certain).fillna(False).astype(bool)
    sb["uncertain_only"] = sb["has_tertiary"] & ~sb["certain_tertiary"]
    groups = []
    for lab, sel in [("P1 < 7 d", sb["P"] <= 7.0), ("P1 > 7 d", sb["P"] > 7.0), ("All", sb["P"] > 0),
                     ("P1 < 2.9 d", sb["P"] < 2.9)]:
        g = sb[sel]
        g2 = g[~g["uncertain_only"]]
        groups.append(dict(sample=lab, N=int(len(g)), f_raw=float(g["has_tertiary"].mean()),
                           N_certain=int(len(g2)), f_raw_certain=float(g2["has_tertiary"].mean())))
    # quartiles in period, as in their Fig. 14
    qs = np.quantile(sb["P"], [0, 0.25, 0.5, 0.75, 1.0])
    quart = []
    for a, b in zip(qs[:-1], qs[1:]):
        g = sb[(sb["P"] >= a) & (sb["P"] <= b)]
        k, n = int(g["has_tertiary"].sum()), len(g)
        quart.append(dict(P_lo=float(a), P_hi=float(b), N=n, k=k, f_raw=k / n,
                          f_err=float(np.sqrt(k / n * (1 - k / n) / n))))
    boundary = sb[(sb["P"] > 6.9) & (sb["P"] < 7.1)][["hip", "comp", "P"]].to_dict("records")
    return dict(n_sb=int(len(sb)), n_systems=int(sb["hip"].nunique()), groups=groups, quartiles=quart,
                near_7d=boundary, periods=sb["P"].tolist(), has_tertiary=sb["has_tertiary"].tolist())


def run_all() -> dict:
    comp, sys_, orb, notes, pairs = msc_pairs()
    st = stability_table(pairs)
    summ = stability_summary(st)
    t17 = tokovinin2017_reproduction()
    model_t17 = fit_alignment_model(np.array(t17["two_orbit_phis"]))
    v26 = msc2026_two_orbits(pairs, orb)
    sense = np.where(v26["i_in"] < 90, 1, -1) * np.where(v26["i_out"] < 90, 1, -1)
    c26, c26e, n26 = sign_correlation(sense)
    model_26 = fit_alignment_model(v26[["phi1", "phi2"]].to_numpy())
    tight = v26[v26["s_au"] < 50]
    c26t, c26te, n26t = sign_correlation(np.where(tight["i_in"] < 90, 1, -1) * np.where(tight["i_out"] < 90, 1, -1))
    t06 = tokovinin2006_reproduction()
    st.to_csv(C.DERIVED / "msc_couples.csv", index=False)
    v26.to_csv(C.DERIVED / "msc2026_two_visual_orbits.csv", index=False)
    return dict(
        msc=dict(n_systems=int(comp[comp["sep"].fillna(0) == 0]["wds"].nunique()),
                 n_subsystems=int(len(sys_)), n_orbits=int(len(orb)), **summ),
        tokovinin2017=t17, tokovinin2017_model=model_t17,
        msc2026_alignment=dict(n=int(len(v26)), C=c26, C_err=c26e, mean_phi=mean_phi_from_c(c26),
                               N_s_lt_50=n26t, C_s_lt_50=c26t, C_s_lt_50_err=c26te, model=model_26),
        tokovinin2006=t06,
    )


if __name__ == "__main__":
    import json
    res = run_all()
    slim = json.loads(json.dumps(res, default=float))
    for k in ("running", "two_orbit_phis"):
        slim["tokovinin2017"].pop(k)
    for m in ("tokovinin2017_model",):
        slim[m].pop("lnL")
    slim["msc2026_alignment"]["model"].pop("lnL")
    slim["tokovinin2006"].pop("periods")
    slim["tokovinin2006"].pop("has_tertiary")
    print(json.dumps(slim, indent=1)[:9000])
