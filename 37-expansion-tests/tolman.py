"""Test 3: are distant galaxies too dim for a static universe?

Surface brightness (light per square arcsecond) is independent of distance in
ordinary static space: a galaxy twice as far is four times fainter but also
four times smaller in area. In an expanding universe it drops as (1 + z)^4
(Tolman 1930): one (1 + z) because each photon arrives with less energy, one
because photons arrive less often (time dilation), and two because the galaxy
was closer when the light left, so it looks bigger than its distance today
suggests (aberration). Tired light keeps only the first factor: (1 + z)^1.

Lubin & Sandage (2001, Papers I-IV) measured the mean surface brightness
<SB> of 34 early-type galaxies in three clusters at z = 0.76-0.92 with the
Hubble Space Telescope, and compared it with nearby elliptical galaxies *of
the same physical size*, measured inside the same Petrosian radius (a radius
defined by the light profile itself, eta = <SB> inside minus SB at the
radius). The depression, written as 2.5 log10 (1 + z)^n magnitudes, gives n.

The catch: converting an angular radius to a physical one needs a distance,
which depends on the cosmological model. So the test is one of consistency:
compute radii *in each model*, measure n, and ask whether n matches that
model's prediction once galaxies' evolution is allowed for (they were younger,
hence brighter, in the past, which makes n come out *below* the true value).

This module re-reads their tables from the arXiv TeX sources, refits the local
relations from the 100-odd nearby galaxies, and redoes the test for several
cosmologies, including tired light.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from config import (C_KM_S, ETAS, ETAS_FINAL, H0_TOLMAN, LUBIN_SANDAGE01, PAPER_DIR,
                    R_MINUS_I_LOCAL, TOLMAN_CLUSTERS)

ARCSEC = np.pi / (180 * 3600)
PAPER1 = PAPER_DIR / "sandage_lubin_2001_paper1" / "Sandage.tex"
PAPER3 = PAPER_DIR / "lubin_sandage_2001_paper3" / "tolman3.tex"
LOCAL_ETAS = (1.0, 1.3, 1.5, 1.7, 2.0, 2.5)

# ---------------------------------------------------------------------------
# Reading the published tables
# ---------------------------------------------------------------------------


def _tables(path):
    """The \\startdata ... \\enddata blocks of a TeX file, in order."""
    return re.findall(r"\\startdata(.*?)\\enddata", path.read_text(), flags=re.S)


def _cells(line):
    line = re.sub(r"\\tablenotemark\{[a-z]\}", "", line)
    line = line.replace("\\\\", "").replace("$", "")
    return [c.strip() for c in line.split("&")]


def _value(cell):
    """'(4.969)' -> (4.969, True); '' or \\nodata -> (nan, False)."""
    cell = cell.strip()
    if not cell or "nodata" in cell or cell == "-":
        return np.nan, False
    uncertain = cell.startswith("(")
    return float(cell.strip("()").split("\\pm")[0]), uncertain


def local_galaxies() -> pd.DataFrame:
    """Paper I Table 1: first-ranked cluster ellipticals (Postman & Lauer 1995), R band.

    One row per galaxy and Petrosian eta: redshift, the distance factor f
    (log R[pc] = f + log r''), M, log R and <SB>, and whether the value was
    extrapolated (parentheses in the table).
    """
    rows = [l for l in _tables(PAPER1)[0].splitlines() if "&" in l]
    out = []
    i = 0
    while i < len(rows):
        cells = _cells(rows[i])
        if cells[0] and cells[0][0].isalpha():
            name = cells[0]
            zc, mm, fc = (_cells(rows[i + k]) for k in (1, 2, 3))
            z = float(zc[0])
            f = float(fc[0])
            for j, eta in enumerate(LOCAL_ETAS, start=1):
                M, u1 = _value(zc[j])
                logR, u2 = _value(mm[j])
                sb, u3 = _value(fc[j])
                if np.isfinite(logR) and np.isfinite(sb):
                    out.append({"cluster": name, "z": z, "dist_mod": float(mm[0]), "f": f,
                                "eta": eta, "M": M, "logR": logR, "SB": sb,
                                "uncertain": u1 or u2 or u3})
            i += 4
        else:
            i += 1
    return pd.DataFrame(out)


def local_fits_published() -> pd.DataFrame:
    """Paper I Table 2: <SB> = a log R + b at <z> = 0.037 and at z = 0."""
    rows = [_cells(l) for l in _tables(PAPER1)[1].splitlines() if "&" in l]
    out = []
    for c in rows:
        a, ae = (float(v) for v in c[1].split("\\pm"))
        b, be = (float(v) for v in c[2].split("\\pm"))
        out.append({"eta": float(c[0]), "slope": a, "slope_err": ae, "zp": b, "zp_err": be,
                    "zp_z0": float(c[3]), "valid_min": float(c[4].replace(">", ""))})
    return pd.DataFrame(out)


def nonlinearity_correction() -> pd.DataFrame:
    """Paper I Table 3: correction to the local lines at small radii (from
    Sandage & Perelmuter 1991). The lines are made fainter by |Delta|."""
    rows = [_cells(l) for l in _tables(PAPER1)[2].splitlines() if "&" in l]
    return pd.DataFrame({"logR": [float(c[0]) for c in rows],
                         "delta": [abs(float(c[1])) for c in rows],
                         "sigma": [float(c[2]) for c in rows]}).sort_values("logR")


def distant_galaxies() -> pd.DataFrame:
    """Paper III Tables 5-7: K-corrected photometry of the 34 cluster galaxies."""
    tables = _tables(PAPER3)[4:7]          # Tables 1-3 are single galaxies, 4 is K terms
    names = list(TOLMAN_CLUSTERS)
    out = []
    for name, tab in zip(names, tables):
        for line in tab.splitlines():
            if "&" not in line:
                continue
            photz = "tablenotemark{c}" in line
            c = _cells(line)
            gal = c[0]
            for k, eta in enumerate(ETAS):
                mag_cell, sb_cell, r_cell = c[1 + 3 * k: 4 + 3 * k]
                if "nodata" in mag_cell:
                    continue
                mag, merr = (float(v) for v in mag_cell.split("\\pm"))
                out.append({"cluster": name, "gal": gal, "eta": eta, "mag": mag,
                            "mag_err": merr, "SB": float(sb_cell), "logr_arcsec": float(r_cell),
                            "photz_member": photz})
    df = pd.DataFrame(out)
    df["z"] = df["cluster"].map({k: v["z"] for k, v in TOLMAN_CLUSTERS.items()})
    df["band"] = df["cluster"].map({k: v["band"] for k, v in TOLMAN_CLUSTERS.items()})
    return df


# ---------------------------------------------------------------------------
# Distances in each world model
# ---------------------------------------------------------------------------


def angular_diameter_distance_pc(z, model, h0=H0_TOLMAN):
    """Distance that converts angle to physical size, in parsecs.

    * ``q0=x``: matter-only Friedmann models (Mattig 1958), as used by Lubin &
      Sandage;
    * ``lcdm``: flat, Omega_m = 0.3, Omega_Lambda = 0.7;
    * ``tired``: static Euclidean space with D = (c/H0) ln(1 + z) (Lubin &
      Sandage's choice, from Sandage 1995). Nothing expands, so the distance
      when the light left is the distance now.
    """
    z = np.asarray(z, float)
    dh = C_KM_S / h0 * 1e6
    if model == "tired":
        return dh * np.log1p(z)
    if model == "lcdm":
        from astropy.cosmology import FlatLambdaCDM
        return FlatLambdaCDM(H0=h0, Om0=0.3, Tcmb0=0).angular_diameter_distance(z).value * 1e6
    q0 = float(model.split("=")[1])
    if q0 == 0:
        return dh * z * (1 + z / 2) / (1 + z) ** 2
    dl = dh / q0 ** 2 * (z * q0 + (q0 - 1) * (np.sqrt(2 * q0 * z + 1) - 1))
    return dl / (1 + z) ** 2


def log_a(z, model):
    """log10 of parsecs per arcsecond."""
    return np.log10(angular_diameter_distance_pc(z, model) * ARCSEC)


MODELS = {
    "q0=0.5": {"label": "matter only, q0 = 1/2 (Lubin & Sandage)", "n_true": 4},
    "q0=0": {"label": "empty, q0 = 0", "n_true": 4},
    "q0=1": {"label": "closed, q0 = 1", "n_true": 4},
    "lcdm": {"label": "flat LambdaCDM, Omega_m = 0.3 (today's model)", "n_true": 4},
    "tired": {"label": "static universe with tired light", "n_true": 1},
}


# ---------------------------------------------------------------------------
# The test
# ---------------------------------------------------------------------------


#: Where Paper I's small-radius correction takes over from the linear fit.
LOGR_ANCHOR = 4.5
#: Mean redshift of Paper I's local galaxies.
Z_LOCAL = 0.037


def bisector_fit(x, y):
    """OLS-bisector line (Isobe et al. 1990), symmetric in x and y."""
    a1, _ = np.polyfit(x, y, 1)
    a2 = 1.0 / np.polyfit(y, x, 1)[0]
    a = (a1 * a2 - 1 + np.sqrt((1 + a1 ** 2) * (1 + a2 ** 2))) / (a1 + a2)
    return a, float(np.mean(y) - a * np.mean(x))


def refit_local(local):
    """Our own fits of the local relations (q0 = 1/2 radii, <z> = 0.037).

    Paper I does not say which least-squares line it fitted; the OLS bisector
    through all galaxies reproduces its eta = 2.0 line, and through the
    non-extrapolated values its eta = 1.7 line, to 0.01-0.02 mag. (eta = 1.0
    is not from these galaxies at all but from Sandage & Perelmuter 1991.)
    """
    out = {}
    for eta in ETAS:
        if eta == 1.0:
            continue
        d = local[local["eta"] == eta]
        a, b = bisector_fit(d["logR"].to_numpy(), d["SB"].to_numpy())
        rng = np.random.default_rng(1)
        boot = []
        for _ in range(500):
            k = rng.integers(0, len(d), len(d))
            aa, bb = bisector_fit(d["logR"].to_numpy()[k], d["SB"].to_numpy()[k])
            boot.append(aa * LOGR_ANCHOR + bb)
        out[eta] = {"slope": float(a), "zp": float(b), "n": int(len(d)),
                    "line_anchor": float(a * LOGR_ANCHOR + b),
                    "line_anchor_err": float(np.std(boot))}
    return out


def fit_local(local, model):
    """Zero-redshift <SB>-log R lines in one world model.

    Paper I's published lines (Table 2) are used, because Paper I's
    small-radius corrections (Table 3) were measured relative to exactly those
    lines. They are moved into the model's frame: local radii shift by the
    model's angular-diameter distance relative to q0 = 1/2 at each local
    galaxy's redshift, and the surface brightness is brightened to z = 0 by
    the model's own dimming at <z> = 0.037 (0.16 mag for n = 4, 0.04 mag for
    tired light's n = 1).
    """
    pub = local_fits_published().set_index("eta")
    shift = float(np.mean(log_a(local["z"], model) - log_a(local["z"], "q0=0.5")))
    dim = 2.5 * MODELS[model]["n_true"] * np.log10(1 + Z_LOCAL)
    fits = {}
    for eta in ETAS:
        a, b = pub.loc[eta, "slope"], pub.loc[eta, "zp"]
        fits[eta] = {"slope": float(a), "zp": float(b - a * shift - dim),
                     "line_err_anchor": float(pub.loc[eta, "zp_err"]), "radius_shift": shift}
    return fits, None


def local_line(fit, logr, band, nl):
    """Local <SB> at log R (z = 0), with the small-radius correction and colour."""
    corr = np.interp(logr, nl["logR"], nl["delta"], left=nl["delta"].iloc[0], right=0.0)
    corr = np.where(logr > nl["logR"].max(), 0.0, corr)
    sig = np.interp(logr, nl["logR"], nl["sigma"], left=nl["sigma"].iloc[0], right=0.0)
    sb = fit["slope"] * logr + fit["zp"] + corr
    if band == "I":
        sb = sb - R_MINUS_I_LOCAL
    return sb, sig


def tolman_signal(local, distant, model):
    """Exponent n of the surface-brightness depression, per cluster and eta.

    Error per cluster: the scatter of its galaxies about their mean, divided
    by sqrt(N), combined with errors common to all of them (the local line's
    zero point and small-radius correction, and the K-correction).
    """
    fits, _ = fit_local(local, model)
    nl = nonlinearity_correction()
    d = distant.copy()
    d["logR_model"] = d["logr_arcsec"] + log_a(d["z"], model)
    rows = []
    for (cl, eta), g in d.groupby(["cluster", "eta"]):
        z = TOLMAN_CLUSTERS[cl]["z"]
        band = TOLMAN_CLUSTERS[cl]["band"]
        # Paper I's eta = 1.0 line was fitted to small galaxies directly
        # (Sandage & Perelmuter 1991), so it gets no small-radius correction.
        line, sig_nl = local_line(fits[eta], g["logR_model"].to_numpy(), band,
                                  nl.assign(delta=0.0) if eta == 1.0 else nl)
        dsb = g["SB"].to_numpy() - line
        mean = float(np.mean(dsb))
        stat = float(np.std(dsb, ddof=1) / np.sqrt(len(dsb)))
        sys = float(np.sqrt(fits[eta]["line_err_anchor"] ** 2 + np.mean(sig_nl) ** 2
                            + TOLMAN_CLUSTERS[cl]["K"][1] ** 2))
        err = np.hypot(stat, sys)
        unit = 2.5 * np.log10(1 + z)
        rows.append({"model": model, "cluster": cl, "band": band, "z": z, "eta": eta,
                     "n_gal": int(len(g)), "dSB": mean, "dSB_err": float(err),
                     "n": mean / unit, "n_err": float(err / unit),
                     "mean_logR": float(g["logR_model"].mean())})
    return pd.DataFrame(rows), fits


def combine(signal):
    """Final n in R (Cl 1604+4321) and I (the other two), over eta = 1.7 and 2.0."""
    s = signal[signal["eta"].isin(ETAS_FINAL)]
    out = {}
    for band in ("R", "I"):
        b = s[s["band"] == band]
        w = 1.0 / b["n_err"] ** 2
        out[f"n_{band}"] = float(np.sum(w * b["n"]) / np.sum(w))
        # The eta = 1.7 and 2.0 values share galaxies, so the error is not
        # reduced by combining them: quote the mean single-eta error / sqrt(n clusters).
        out[f"n_{band}_err"] = float(np.mean(b.groupby("eta")["n_err"].apply(
            lambda e: 1 / np.sqrt(np.sum(1 / e ** 2)))))
    return out


def lookback_gyr(z, model, h0=70.0):
    """Look-back time to redshift z (Gyr), for the evolution discussion."""
    if model == "tired":
        return 977.79 / h0 * float(np.log1p(z))      # 1/H0 in Gyr, times ln(1 + z)
    from astropy.cosmology import FlatLambdaCDM
    return float(FlatLambdaCDM(H0=h0, Om0=0.3).lookback_time(z).value)


def run():
    local = local_galaxies()
    distant = distant_galaxies()
    out = {"n_local": int(local["cluster"].nunique()),
           "n_distant": int(distant.groupby(["cluster", "gal"]).ngroups)}
    # How well Paper I's published lines describe the local galaxies.
    ours = refit_local(local)
    pub = local_fits_published().set_index("eta")
    out["local_fit_check"] = {eta: {**v, "slope_pub": float(pub.loc[eta, "slope"]),
                                    "line_anchor_pub": float(pub.loc[eta, "slope"] * LOGR_ANCHOR
                                                             + pub.loc[eta, "zp"])}
                              for eta, v in ours.items()}
    signals, results = [], {}
    for model in MODELS:
        sig, _ = tolman_signal(local, distant, model)
        signals.append(sig)
        results[model] = {**combine(sig), "n_true": MODELS[model]["n_true"],
                          "label": MODELS[model]["label"]}
        for band in ("R", "I"):
            n, e = results[model][f"n_{band}"], results[model][f"n_{band}_err"]
            # Evolution needed to reach the model's prediction: p = n_true - n.
            results[model][f"p_{band}"] = MODELS[model]["n_true"] - n
            results[model][f"p_{band}_err"] = e
    out["signal"] = pd.concat(signals, ignore_index=True)
    out["results"] = results
    out["published"] = LUBIN_SANDAGE01
    out["lookback_gyr"] = {m: {cl: lookback_gyr(v["z"], "tired" if m == "tired" else "lcdm")
                               for cl, v in TOLMAN_CLUSTERS.items()} for m in ("lcdm", "tired")}
    out["local"] = local
    out["distant"] = distant
    return out
