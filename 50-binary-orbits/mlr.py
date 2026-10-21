"""The mass-luminosity relation from visual binary orbits alone.

For every grade 1-2 orbit in ORB6 we find the stars in Gaia DR3, choose the
best parallax (both components' Gaia parallaxes where Gaia resolves the pair,
the photocentre's otherwise, Hipparcos as a fallback), and weigh the pair with
Kepler's third law: M1 + M2 = (a/parallax)^3 / P^2.

A visual orbit gives only the *total* mass, but ORB6 also lists each star's V
magnitude. If every main-sequence star obeys one mass-luminosity relation
log M = f(M_V), the total must equal 10^f(M_V1) + 10^f(M_V2). Fitting f to a
few hundred totals recovers the relation itself, star by star: the method of
Hertzsprung (1923) and of every visual-binary calibration since. Parallax
errors move the mass (as parallax^-3) and the luminosities (as parallax^-2)
together; the fit propagates that correlation.
"""

from __future__ import annotations

import gzip
import re

import numpy as np
import pandas as pd
from scipy import optimize

import config as C
import literature as L
import orb6
import orbits

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def read_hip2(path=C.HIP_DIR / "hip2.dat.gz") -> pd.DataFrame:
    """Hipparcos new reduction (van Leeuwen 2007): HIP, solution type, parallax, errors, Hp."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            # byte layout from I/311 ReadMe (1-indexed): HIP 1-6, Sn 8-10, So 12, Nc 14,
            # Plx 44-50, e_Plx 84-89, pmRA 52-59, pmDE 61-68, Hpmag 130-136
            rows.append((int(line[0:6]), int(line[7:10]), int(line[11]), float(line[43:50]),
                         float(line[83:89]), float(line[51:59]), float(line[60:68]), float(line[129:136])))
    return pd.DataFrame(rows, columns=["hip", "sn", "so", "plx", "e_plx", "pmra", "pmdec", "hpmag"])


def mamajek_table() -> pd.DataFrame:
    """Pecaut & Mamajek dwarf sequence: SpT, Teff, BCv, logL, Mv, B-V, Bp-Rp, M_G, Msun."""
    path = C.LIT_DIR / "mamajek" / "EEM_dwarf_UBVIJHK_colors_Teff.txt"
    lines = path.read_text(encoding="latin-1").splitlines()
    head = None
    rows = []
    for ln in lines:
        if ln.startswith("#SpT"):
            if head is not None and rows:
                break          # a second table follows the main one
            head = ln[1:].split()
            continue
        if head is None or not ln.strip() or ln.startswith("#"):
            continue
        parts = ln.split()
        if len(parts) < len(head) - 1:
            continue
        rows.append(parts[:len(head)])
    df = pd.DataFrame(rows, columns=head[:len(rows[0])])
    out = pd.DataFrame({"spt": df["SpT"]})
    for col in ("Teff", "BCv", "logL", "Mv", "B-V", "Bp-Rp", "M_G", "Msun", "R_Rsun"):
        out[col] = pd.to_numeric(df[col].str.replace("...", "", regex=False).str.replace(":", "", regex=False),
                                 errors="coerce")
    return out


def bc_v(teff):
    """Bolometric correction BC_V as a function of Teff (interpolated in Mamajek's table)."""
    m = mamajek_table().dropna(subset=["Teff", "BCv"]).sort_values("Teff")
    return np.interp(np.log10(teff), np.log10(m["Teff"].to_numpy()), m["BCv"].to_numpy())


def eker_curve(mass):
    """Eker et al. 2018's relations as M_V(mass): MLR -> L, MRR/MTR -> Teff, Mamajek BC_V."""
    mass = np.asarray(mass, float)
    logm = np.log10(mass)
    logL = np.full_like(mass, np.nan)
    for lo, hi, alpha, b, _ in L.EKER18_MLR:
        sel = (mass > lo) & (mass <= hi)
        logL[sel] = alpha * logm[sel] + b
    a2, a1, a0 = L.EKER18_MRR
    R = a2 * mass**2 + a1 * mass + a0
    teff_lo = C.TEFF_SUN * (10**logL / R**2) ** 0.25
    c2, c1, c0 = L.EKER18_MTR
    teff_hi = 10 ** (c2 * logm**2 + c1 * logm + c0)
    teff = np.where(mass <= 1.5, teff_lo, teff_hi)
    mbol = C.MBOL_SUN - 2.5 * logL
    return {"logL": logL, "teff": teff, "MV": mbol - bc_v(teff), "Mbol": mbol}


def eker_stars() -> pd.DataFrame:
    """Eker et al. 2018 Table 1 (detached eclipsing binaries): M, R, Teff -> L, M_V."""
    path = C.LIT_DIR / "eker2018" / "table1.dat"
    rows = []
    for ln in path.read_text(encoding="latin-1").splitlines():
        ln = ln.ljust(196)
        def num(a, b):
            s = ln[a - 1:b].strip()
            return float(s) if s else np.nan
        rows.append({"name": ln[4:28].strip(), "comp": ln[30], "M": num(94, 100), "e_M": num(102, 107),
                     "R": num(109, 115), "e_R": num(117, 122), "Teff": num(160, 164), "remark": ln[193:196].strip()})
    df = pd.DataFrame(rows)
    df["selected"] = df["remark"].str.startswith("*")      # the 509 main-sequence stars Eker et al. used
    df["logL"] = np.log10(df["R"] ** 2 * (df["Teff"] / C.TEFF_SUN) ** 4)
    df["Mbol"] = C.MBOL_SUN - 2.5 * df["logL"]
    df["MV"] = df["Mbol"] - bc_v(df["Teff"].to_numpy())
    return df


# ---------------------------------------------------------------------------
# Gaia matching
# ---------------------------------------------------------------------------


def load_gaia():
    xm = pd.read_csv(C.GAIA_DIR / "orb6_xmatch_gaiadr3.csv")
    det = pd.read_csv(C.GAIA_DIR / "orb6_gaia_details.csv")
    det.columns = [c.lower() for c in det.columns]
    nss = pd.read_csv(C.GAIA_DIR / "orb6_gaia_nss.csv")
    hbn = pd.read_csv(C.GAIA_DIR / "orb6_hip_best_neighbour.csv")
    hip1 = pd.read_csv(C.GAIA_DIR / "orb6_hipparcos1.csv")
    return xm, det, nss, hbn, hip1


def zero_point(det: pd.DataFrame) -> np.ndarray:
    """Lindegren et al. 2021 parallax zero point (valid 6 < G < 21); 0 outside, with a flag."""
    from zero_point import zpt
    zpt.load_tables()
    g = det["phot_g_mean_mag"].to_numpy()
    sol = det["astrometric_params_solved"].to_numpy()
    ok = (g > 6.0) & (g < 21.0) & np.isin(sol, (31, 95))
    z = np.zeros(len(det))
    if ok.any():
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            z[ok] = zpt.get_zpt(g[ok], det["nu_eff_used_in_astrometry"].to_numpy()[ok],
                                det["pseudocolour"].to_numpy()[ok], det["ecl_lat"].to_numpy()[ok], sol[ok])
    return z


def _hip_spclass(sp: str) -> str:
    """Luminosity class from a Hipparcos spectral type ('' when not given)."""
    sp = str(sp) if isinstance(sp, str) else ""
    m = re.search(r"(VI|V|IV|III|II|Ia|Ib|I)(?![a-z])", sp.replace("IV-V", "V").replace("V-IV", "IV"))
    return m.group(1) if m else ""


def build_sample() -> pd.DataFrame:
    """One row per grade 1-2 ORB6 orbit: Gaia identification, parallax choice, dynamical mass, M_V."""
    o = orb6.read_orb6()
    o["oid"] = o.index
    # the shortest period of any *other* orbit in the same WDS system
    min_other_P = {}
    for wds, grp in o.groupby("wds"):
        for k in grp.index:
            others = grp.drop(index=k)["P_yr"]
            min_other_P[k] = others.min() if len(others) else np.inf
    g = o[o["grade"].isin(C.MLR_GRADES)].copy()
    xm, det, nss, hbn, hip1 = load_gaia()
    det = det.drop_duplicates("source_id").set_index("source_id")
    det["zp"] = zero_point(det.reset_index())
    hip2 = read_hip2().set_index("hip")
    hip1 = hip1.drop_duplicates("hip").set_index("hip")
    nss_by = nss.sort_values("significance", ascending=False).drop_duplicates("source_id").set_index("source_id")

    # default fractional errors where a paper gave none: the grade's 84th percentile
    frac_a = {gr: np.nanpercentile((o.loc[o.grade == gr, "a_err_arcsec"] / o.loc[o.grade == gr, "a_arcsec"]), 84)
              for gr in C.MLR_GRADES}
    frac_P = {gr: np.nanpercentile((o.loc[o.grade == gr, "P_err_yr"] / o.loc[o.grade == gr, "P_yr"]), 84)
              for gr in C.MLR_GRADES}

    rows = []
    for _, r in g.iterrows():
        rec = {k: r[k] for k in ("oid", "wds", "disc", "grade", "hip", "v1", "v2", "v1_flag", "v2_flag",
                                 "P_yr", "a_arcsec", "e", "i", "ref", "last_obs")}
        rec["a_err"] = r["a_err_arcsec"] if np.isfinite(r["a_err_arcsec"]) and r["a_err_arcsec"] > 0 \
            else frac_a[r["grade"]] * r["a_arcsec"]
        rec["P_err"] = r["P_err_yr"] if np.isfinite(r["P_err_yr"]) and r["P_err_yr"] > 0 \
            else frac_P[r["grade"]] * r["P_yr"]
        rec["err_assumed"] = not (np.isfinite(r["a_err_arcsec"]) and np.isfinite(r["P_err_yr"]))
        # an outer orbit of a hierarchy: another orbit in the system is shorter, so
        # this one's "component" is itself a pair (three or more stars in the total)
        rec["multiple"] = bool(min_other_P.get(r["oid"], np.inf) < r["P_yr"])
        comp = re.sub(r"^[A-Z]+\s*\d*", "", r["disc"]).replace(" ", "")
        rec["components"] = comp
        # Hipparcos colour/type describe star A: usable for AB, Aa,Ab or unlabelled pairs only
        rec["hip_is_pair"] = comp in ("", "AB", "Aa,Ab", "Aa1,Aa2", "AaAb") or comp.startswith("Aa")
        # predicted separation/PA at Gaia's epoch
        rho16, th16 = orbits.rho_theta(2016.0, *orb6.elements(r))
        rec["rho2016"], rec["pa2016"] = float(rho16), float(th16)
        vcomb = r["v1"] if not np.isfinite(r["v2"]) else -2.5 * np.log10(10 ** (-0.4 * r["v1"]) + 10 ** (-0.4 * r["v2"]))
        rec["v_comb"] = vcomb
        # --- Gaia candidates near the WDS position (J2000)
        c = xm[(xm["oid"] == r["oid"]) & (xm["Gmag"] < vcomb + 3.0)].copy()
        prim = sec = None
        if len(c):
            ra_j = np.where(np.isfinite(c["RAJ2000"]), c["RAJ2000"], c["RAdeg"])
            de_j = np.where(np.isfinite(c["DEJ2000"]), c["DEJ2000"], c["DEdeg"])
            dra = (ra_j - r["ra_deg"]) * np.cos(np.radians(r["dec_deg"])) * 3600
            dde = (de_j - r["dec_deg"]) * 3600
            c["d2000"] = np.hypot(dra, dde)
            near = c[(c["d2000"] < C.GAIA_MATCH_TOL_ARCSEC + rho16) & (np.abs(c["Gmag"] - vcomb) < 2.5)]
            if r["hip"] > 0 and (hbn["original_ext_source_id"] == r["hip"]).any():
                sid = int(hbn.loc[hbn["original_ext_source_id"] == r["hip"], "source_id"].iloc[0])
                if sid in set(c["Source"]):
                    near = pd.concat([c[c["Source"] == sid], near[near["Source"] != sid]])
            if len(near):
                prim = near.sort_values(["Gmag"]).iloc[0] if not (r["hip"] > 0 and len(near) and
                                                                  near.iloc[0]["Source"] in set(hbn["source_id"])) \
                    else near.iloc[0]
                # brightest close source is the primary unless the HIP best neighbour says otherwise
                bright = near.sort_values("Gmag").iloc[0]
                if bright["Gmag"] < prim["Gmag"] - 0.5:
                    prim = bright
                # secondary: a source at the predicted separation (either quadrant)
                cos_d = np.cos(np.radians(prim["DEdeg"]))
                sx = (c["RAdeg"] - prim["RAdeg"]) * cos_d * 3600
                sy = (c["DEdeg"] - prim["DEdeg"]) * 3600
                c["sep_p"] = np.hypot(sx, sy)
                c["pa_p"] = np.degrees(np.arctan2(sx, sy)) % 360
                tol = max(0.15, 0.15 * rho16)
                dpa = np.minimum(np.abs(orbits.wrap_angle_deg(c["pa_p"] - th16)),
                                 np.abs(orbits.wrap_angle_deg(c["pa_p"] - th16 - 180)))
                cand = c[(c["Source"] != prim["Source"]) & (np.abs(c["sep_p"] - rho16) < tol) & (dpa < 15.0)]
                if len(cand):
                    sec = cand.sort_values("Gmag").iloc[0]
        rec["resolved"] = sec is not None
        for lab, s in (("p", prim), ("s", sec)):
            if s is None:
                continue
            sid = int(s["Source"])
            rec[f"gaia_{lab}"] = sid
            rec[f"G_{lab}"] = s["Gmag"]
            rec[f"bp_rp_{lab}"] = s["BP-RP"]
            rec[f"plx_{lab}"] = s["Plx"]
            rec[f"e_plx_{lab}"] = s["e_Plx"]
            rec[f"ruwe_{lab}"] = s["RUWE"]
            rec[f"solved_{lab}"] = s["Solved"]
            rec[f"zp_{lab}"] = det["zp"].get(sid, 0.0) if sid in det.index else 0.0
            rec[f"nss_{lab}"] = nss_by["nss_solution_type"].get(sid, "") if sid in nss_by.index else ""
        rows.append(rec)
    df = pd.DataFrame(rows)
    return _choose_parallax(df, hip2, hip1, nss_by)


def _good(df, lab):
    return ((df.get(f"solved_{lab}", np.nan).isin([31, 95])) & (df.get(f"ruwe_{lab}", np.inf) < C.GAIA_RUWE_MAX)
            & (df.get(f"plx_{lab}", np.nan) / df.get(f"e_plx_{lab}", np.nan) > C.GAIA_PLX_SNR_MIN))


def _choose_parallax(df, hip2, hip1, nss_by) -> pd.DataFrame:
    """Gaia (both components if resolved and consistent) > Gaia NSS orbit > Hipparcos-2."""
    for col in ("plx_p", "e_plx_p", "ruwe_p", "solved_p", "zp_p", "plx_s", "e_plx_s", "ruwe_s", "solved_s",
                "zp_s", "gaia_p", "gaia_s", "nss_p"):
        if col not in df:
            df[col] = np.nan
    gp, gs = _good(df, "p"), _good(df, "s")
    plx_p = df["plx_p"] - df["zp_p"].fillna(0.0)     # zero point: plx_true = plx - zp
    plx_s = df["plx_s"] - df["zp_s"].fillna(0.0)
    plx = np.full(len(df), np.nan)
    eplx = np.full(len(df), np.nan)
    src = np.full(len(df), "", dtype=object)
    both = gp & gs
    w1, w2 = 1 / df["e_plx_p"] ** 2, 1 / df["e_plx_s"] ** 2
    consistent = (np.abs(plx_p - plx_s) < 3 * np.hypot(df["e_plx_p"], df["e_plx_s"]))
    m = both & consistent
    plx[m] = ((plx_p * w1 + plx_s * w2) / (w1 + w2))[m]
    eplx[m] = (1 / np.sqrt(w1 + w2))[m]
    src[m] = "Gaia A+B"
    m = (gp & ~(both & consistent)) | (both & ~consistent & (df["ruwe_p"] <= df["ruwe_s"]))
    plx[m], eplx[m], src[m] = plx_p[m], df["e_plx_p"][m], "Gaia"
    m = gs & ~gp & (src == "")
    plx[m], eplx[m], src[m] = plx_s[m], df["e_plx_s"][m], "Gaia (B)"
    df["plx_gaia_inconsistent"] = both & ~consistent
    # Gaia's own orbital solution for the photocentre (parallax corrected for the wobble)
    for k in np.where(src == "")[0]:
        sid = df.at[k, "gaia_p"]
        if np.isfinite(sid) and int(sid) in nss_by.index and np.isfinite(nss_by.at[int(sid), "parallax"]):
            p_, e_ = nss_by.at[int(sid), "parallax"], nss_by.at[int(sid), "parallax_error"]
            if p_ / e_ > C.GAIA_PLX_SNR_MIN:
                plx[k], eplx[k], src[k] = p_ - (df.at[k, "zp_p"] if np.isfinite(df.at[k, "zp_p"]) else 0), e_, \
                    f"Gaia NSS ({nss_by.at[int(sid), 'nss_solution_type']})"
    # Hipparcos new reduction as the fallback (skip stochastic solutions, s = 1)
    df["plx_hip2"] = df["hip"].map(hip2["plx"])
    df["e_plx_hip2"] = df["hip"].map(hip2["e_plx"])
    df["sn_hip2"] = df["hip"].map(hip2["sn"])
    for k in np.where(src == "")[0]:
        p_, e_, sn = df.at[k, "plx_hip2"], df.at[k, "e_plx_hip2"], df.at[k, "sn_hip2"]
        if np.isfinite(p_) and e_ > 0 and p_ / e_ > C.GAIA_PLX_SNR_MIN and int(sn) % 10 != 1:
            plx[k], eplx[k], src[k] = p_, e_, "Hipparcos"
    df["plx"], df["e_plx"], df["plx_source"] = plx, eplx, src
    # dynamical mass
    a, P = df["a_arcsec"], df["P_yr"]
    df["Mtot"] = orbits.total_mass(a, df["plx"] / 1e3, P)
    df["Mtot_frac_err"] = orbits.total_mass_frac_err(a, df["a_err"], df["plx"], df["e_plx"], P, df["P_err"])
    df["dist_pc"] = 1e3 / df["plx"]
    df["a_au"] = a / (df["plx"] / 1e3)
    # absolute magnitudes (V); only plain V magnitudes for both stars are used
    vok = (df["v1_flag"] == "") & (df["v2_flag"] == "") & np.isfinite(df["v2"])
    df["v_ok"] = vok
    mu = 5 * np.log10(df["plx"] / 1e3) + 5
    df["MV1"], df["MV2"] = df["v1"] + mu, df["v2"] + mu
    df["MV_comb"] = df["v_comb"] + mu
    df["dm"] = df["v2"] - df["v1"]
    # Hipparcos colour and spectral type, for the evolved-star test
    df["b_v"] = df["hip"].map(hip1["b_v"]).where(df["hip_is_pair"])
    df["sptype"] = df["hip"].map(hip1["sptype"]).where(df["hip_is_pair"]).fillna("")
    df["lumclass"] = df["sptype"].map(_hip_spclass)
    ms = mamajek_table().dropna(subset=["B-V", "Mv"]).sort_values("B-V")
    mv_ms = np.interp(df["b_v"], ms["B-V"], ms["Mv"], left=np.nan, right=np.nan)
    df["dMV_ms"] = df["MV_comb"] - mv_ms          # negative = brighter than a single MS star
    df["evolved"] = (df["dMV_ms"] < -C.EVOLVED_DMAG) | df["lumclass"].isin(["III", "II", "I", "Ia", "Ib"])
    return df


# ---------------------------------------------------------------------------
# The fit: log M = f(M_V), from total masses
# ---------------------------------------------------------------------------

#: Knots of the piecewise-linear relation (absolute V magnitude).
KNOTS = np.array([-1.0, 1.0, 2.5, 4.0, 5.5, 7.0, 8.5, 10.0, 12.0, 14.0])
#: Error assumed for ORB6's (heterogeneous) component V magnitudes.
V_ERR = 0.05
#: Width of the outlier component, dex.
OUT_WIDTH_DEX = 0.3


def fit_selection(s: pd.DataFrame, parallax="any") -> pd.Series:
    """Systems used in the fit: a parallax, < 100 pc, both V magnitudes, not evolved, not an outer orbit."""
    m = (s["plx_source"] != "") & (s["dist_pc"] < C.MLR_MAX_DIST_PC) & s["v_ok"] & ~s["evolved"] \
        & ~s["multiple"] & (s["Mtot_frac_err"] < C.MLR_MAX_MASS_FRAC_ERR) & s["MV2"].between(-5, 18)
    if parallax == "gaia":
        m &= s["plx_source"].str.startswith("Gaia")
    return m


class MLRLikelihood:
    def __init__(self, d: pd.DataFrame, knots=KNOTS):
        self.knots = knots
        self.logM = np.log10(d["Mtot"].to_numpy())
        self.MV1, self.MV2 = d["MV1"].to_numpy(), d["MV2"].to_numpy()
        self.slogplx = (d["e_plx"] / d["plx"]).to_numpy() / np.log(10)
        self.sloga = (d["a_err"] / d["a_arcsec"]).to_numpy() / np.log(10)
        self.slogP = (d["P_err"] / d["P_yr"]).to_numpy() / np.log(10)

    def f(self, theta, MV):
        nodes = self.nodes(theta)
        # linear interpolation, linear extrapolation beyond the end knots
        k = self.knots
        y = np.interp(MV, k, nodes)
        lo, hi = MV < k[0], MV > k[-1]
        y = np.where(lo, nodes[0] + (MV - k[0]) * (nodes[1] - nodes[0]) / (k[1] - k[0]), y)
        y = np.where(hi, nodes[-1] + (MV - k[-1]) * (nodes[-1] - nodes[-2]) / (k[-1] - k[-2]), y)
        return y

    def slope(self, theta, MV, h=0.01):
        return (self.f(theta, MV + h) - self.f(theta, MV - h)) / (2 * h)

    def nodes(self, theta):
        """Node values: log M at the brightest knot, then strictly decreasing steps."""
        n = len(self.knots)
        steps = np.exp(theta[1:n])
        return theta[0] - np.concatenate([[0.0], np.cumsum(steps)])

    def residuals(self, theta):
        f1, f2 = self.f(theta, self.MV1), self.f(theta, self.MV2)
        m1, m2 = 10**f1, 10**f2
        pred = np.log10(m1 + m2)
        w1, w2 = m1 / (m1 + m2), m2 / (m1 + m2)
        s1, s2 = self.slope(theta, self.MV1), self.slope(theta, self.MV2)
        # parallax moves log M_obs by -3 dlogplx and log M_pred by 5 (w1 s1 + w2 s2) dlogplx
        dr_dplx = -3.0 - 5.0 * (w1 * s1 + w2 * s2)
        var = (dr_dplx * self.slogplx) ** 2 + (3 * self.sloga) ** 2 + (2 * self.slogP) ** 2 \
            + (w1 * s1 * V_ERR) ** 2 + (w2 * s2 * V_ERR) ** 2
        return self.logM - pred, var

    def __call__(self, theta):
        n = len(self.knots)
        ln_sint, lf = theta[n], theta[n + 1]
        if not (-8 < ln_sint < 0 and -8 < lf < 0) or np.any(np.abs(theta[1:n]) > 8):
            return -np.inf
        r, var = self.residuals(theta)
        s2 = var + np.exp(2 * ln_sint)
        f = 1 / (1 + np.exp(-lf))
        lg = -0.5 * r * r / s2 - 0.5 * np.log(2 * np.pi * s2)
        so2 = s2 + OUT_WIDTH_DEX**2
        lo = -0.5 * r * r / so2 - 0.5 * np.log(2 * np.pi * so2)
        return float(np.sum(np.logaddexp(np.log1p(-f) + lg, np.log(f) + lo)))

    def outlier_prob(self, theta):
        n = len(self.knots)
        r, var = self.residuals(theta)
        s2 = var + np.exp(2 * theta[n])
        f = 1 / (1 + np.exp(-theta[n + 1]))
        lg = np.log1p(-f) - 0.5 * r * r / s2 - 0.5 * np.log(s2)
        so2 = s2 + OUT_WIDTH_DEX**2
        lo = np.log(f) - 0.5 * r * r / so2 - 0.5 * np.log(so2)
        return np.exp(lo - np.logaddexp(lg, lo))

    def start(self):
        """Start from Mamajek's dwarf sequence."""
        m = mamajek_table().dropna(subset=["Mv", "Msun"]).sort_values("Mv")
        nodes = np.log10(np.interp(self.knots, m["Mv"], m["Msun"]))
        steps = np.maximum(-np.diff(nodes), 1e-3)
        return np.concatenate([[nodes[0]], np.log(steps), [np.log(0.03), np.log(0.05)]])


def fit_mlr(d: pd.DataFrame, nwalkers=48, nsteps=6000, burn=2000, seed=C.SEED, cache=None):
    import emcee
    like = MLRLikelihood(d)
    if cache is not None and cache.exists():
        z = np.load(cache)
        return like, z["chain"], z["best"]
    th0 = like.start()
    res = optimize.minimize(lambda t: -like(t), th0, method="Nelder-Mead",
                            options={"maxiter": 40000, "maxfev": 40000, "xatol": 1e-7, "fatol": 1e-9})
    res = optimize.minimize(lambda t: -like(t), res.x, method="Powell", options={"maxiter": 40000, "xtol": 1e-8})
    best = res.x
    rng = np.random.default_rng(seed)
    p0 = best + 1e-3 * rng.standard_normal((nwalkers, len(best)))
    sampler = emcee.EnsembleSampler(nwalkers, len(best), like)
    sampler.random_state = np.random.RandomState(seed).get_state()
    sampler.run_mcmc(p0, nsteps, progress=False)
    chain = sampler.get_chain(discard=burn, thin=10, flat=True)
    if cache is not None:
        np.savez_compressed(cache, chain=chain, best=best)
    return like, chain, best


def relation_band(like, chain, MV=np.linspace(-1.5, 15.0, 200), n=2000, seed=C.SEED):
    rng = np.random.default_rng(seed)
    pick = chain[rng.choice(len(chain), min(n, len(chain)), replace=False)]
    curves = np.array([like.f(t, MV) for t in pick])
    return MV, np.percentile(curves, [2.5, 16, 50, 84, 97.5], axis=0)
