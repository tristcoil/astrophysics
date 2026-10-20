"""Eclipsing binaries whose orbits precess mostly because of general relativity.

In an eccentric eclipsing binary the secondary eclipse does not fall halfway
between primary eclipses: how far it is displaced depends on e cos(omega),
where omega is the argument of periastron. As the orbit precesses, omega turns
and the spacing ``T2 - T1`` between primary and secondary eclipse drifts.
Sterne's (1939) formula, for an edge-on orbit:

    T2 - T1 = P/2 + (P/pi) [ arctan(e cos w / sqrt(1-e^2))
                              + e cos w sqrt(1-e^2) / (1 - e^2 sin^2 w) ]

so the drift rate d(T2-T1)/dN (per orbital cycle N) divided by
d(T2-T1)/dw gives the apsidal motion rate dw/dN.

Two physical effects turn the orbit: **classical** ones (each star is
distorted by tides and by its own spin, and a non-spherical star does not pull
like a point mass) and **general relativity**. Baroch et al. (2021, A&A 649,
A64) chose binaries of small, centrally concentrated Sun-like stars in wide
orbits, where GR does 60-99% of the turning. We read their tables from the
paper's LaTeX source, re-derive the rates from the eclipse timings they
tabulate, and take their classical terms (which need stellar-structure models)
as published.
"""

from __future__ import annotations

import io
import re
import tarfile

import numpy as np

from config import AU_M, C_SI, GM_SUN_SI, PUBLISHED, RAW_DIR

ARXIV_DIR = RAW_DIR / "arxiv"
BAROCH_TARBALL = ARXIV_DIR / "2103.03140v2.tar.gz"
BAROCH_TEX_NAME = "aaMAIN.tex"
DAY_S = 86400.0

#: The systems Baroch et al. used for their test of GR (Sect. 6); V459 Cas and
#: V501 Mon were left out for characterisation reasons, V501 Her as a 2-sigma
#: detection.
BAROCH_TEST_SYSTEMS = ("KX Cnc", "AL Dor", "V541 Cyg", "KW Hya", "GG Ori", "EY Cep")


# --------------------------------------------------------------------------
# Reading the LaTeX
# --------------------------------------------------------------------------


def baroch_tex() -> str:
    if not BAROCH_TARBALL.exists():
        raise FileNotFoundError(f"{BAROCH_TARBALL} missing -- run `python fetch_data.py` first")
    with tarfile.open(BAROCH_TARBALL) as tf:
        member = next(m for m in tf.getmembers() if m.name.endswith(BAROCH_TEX_NAME))
        return tf.extractfile(member).read().decode("utf-8", "replace")


def _clean(cell: str) -> str:
    cell = re.sub(r"\\multirow\{[^}]*\}\{[^}]*\}\{([^}]*)\}", r"\1", cell)
    cell = re.sub(r"\\multicolumn\{[^}]*\}\{[^}]*\}\{([^}]*)\}", r"\1", cell)
    cell = cell.replace("$-$", "-").replace("\\,", " ").replace("$", "").replace("\\pm", "+-")
    return re.sub(r"\s+", " ", cell).strip()


def tabular_rows(tex: str, start: int) -> list[list[str]]:
    """Rows of the first tabular after ``start``, comments dropped, cells cleaned."""
    b = tex.index("\\begin{tabular}", start)
    e = tex.index("\\end{tabular}", b)
    body = tex[b:e].split("\n", 1)[1]
    lines = [ln for ln in body.splitlines() if not ln.lstrip().startswith("%")]
    body = "\n".join(lines)
    for junk in ("\\hline", "\\noalign{\\smallskip}"):
        body = body.replace(junk, "")
    rows = []
    for raw in body.split("\\\\"):
        if not raw.strip():
            continue
        cells = [_clean(c.replace("\\&", "&")) for c in re.split(r"(?<!\\)&", raw)]
        rows.append(cells)
    return rows


def table_by_label(tex: str, label: str) -> list[list[str]]:
    return tabular_rows(tex, tex.index(f"\\label{{{label}}}"))


def value(cell: str) -> tuple[float, float]:
    """'0.000163(6)' -> (0.000163, 0.000006); '6.4(1.0)' -> (6.4, 1.0); 'x+-y' -> (x, y)."""
    cell = cell.replace(" ", "")
    if "+-" in cell:
        v, err = cell.split("+-")
        return float(v), float(err)
    m = re.fullmatch(r"(-?\d*\.?\d*)\(([\d.]+)\)", cell)
    if not m:
        return float(cell), np.nan
    v, e = m.group(1), m.group(2)
    if "." in e:  # error written out in full, e.g. 6.4(1.0)
        return float(v), float(e)
    decimals = len(v.split(".")[1]) if "." in v else 0
    return float(v), float(e) * 10.0 ** (-decimals)


def _name(cell: str) -> str:
    return cell.replace("  ", " ").strip()


def baroch_tables() -> dict:
    tex = baroch_tex()
    out = {"systems": {}, "timings": {}}
    for row in table_by_label(tex, "tab:systemsprop")[1:]:
        if len(row) < 9:
            continue
        name = _name(row[0])
        P, eP = value(row[1])
        e, ee = value(row[2])
        out["systems"][name] = {"P": P, "P_err": eP, "e": e, "e_err": ee,
                                "m1": value(row[3])[0], "m2": value(row[4])[0],
                                "m1_err": value(row[3])[1], "m2_err": value(row[4])[1],
                                "r1": value(row[5])[0], "r2": value(row[6])[0],
                                "gr_share_predicted_pct": float(row[7])}
    for row in table_by_label(tex, "tab:otherprop")[2:]:
        if len(row) < 8:
            continue
        name = _name(row[0])
        out["systems"][name].update({"T0": float(row[1]), "omega0_deg": float(row[2]),
                                     "inc_deg": value(row[7])[0]})
    for row in table_by_label(tex, "tab:results")[2:]:
        if len(row) < 9:
            continue
        name = _name(row[0])
        keys = ["obs", "theory", "classical", "gr", "gr_measured", "ratio"]
        for k, cell in zip(keys, row[3:9]):
            v, err = value(cell)
            out["systems"][name][f"wdot_{k}"] = v
            out["systems"][name][f"wdot_{k}_err"] = err
    # Table 4 (TESS T2 - T1) and its continuation.
    i = tex.index("\\label{tab:a1}")
    rows = tabular_rows(tex, i) + tabular_rows(tex, tex.index("\\caption{Continued.}", i))
    current = None
    for row in rows[1:]:
        if len(row) < 3:
            continue
        first = row[0]
        if first and not re.fullmatch(r"[\d.]+", first) and not first.startswith("System"):
            current = _name(first)
        if not row[1] or not row[2] or not re.fullmatch(r"-?\d+", row[1]):
            continue  # header rows, and the T0 line under each system name
        n = float(row[1])
        d, derr = value(row[2])
        out["timings"].setdefault(current, []).append((n, d, derr, "TESS"))
    for label, name in (("tab:V541Cyg", "V541 Cyg"), ("tab:V459Cas", "V459 Cas"),
                        ("tab:EYCep", "EY Cep")):
        for row in table_by_label(tex, label)[1:]:
            if len(row) < 3 or not row[0]:
                continue
            d, derr = value(row[1])
            out["timings"].setdefault(name, []).append((float(row[0]), d, derr, row[2]))
    return out


# --------------------------------------------------------------------------
# Eclipse spacing and apsidal motion
# --------------------------------------------------------------------------


def sterne_spacing(P: float, e: float, omega: np.ndarray) -> np.ndarray:
    """T2 - T1 (same units as P) for an edge-on orbit (Sterne 1939)."""
    s = np.sqrt(1.0 - e**2)
    ec = e * np.cos(omega)
    return P / 2.0 + P / np.pi * (np.arctan(ec / s) + ec * s / (1.0 - e**2 * np.sin(omega) ** 2))


def spacing_derivative(P: float, e: float, omega: float) -> float:
    """d(T2 - T1)/d omega (days per radian), by central difference."""
    h = 1e-6
    return (sterne_spacing(P, e, omega + h) - sterne_spacing(P, e, omega - h)) / (2 * h)


def fit_apsidal_rate(name: str, sysinfo: dict, timings: list) -> dict:
    """Weighted straight line through T2 - T1 against cycle -> d omega / d N (deg/cycle)."""
    n = np.array([t[0] for t in timings])
    d = np.array([t[1] for t in timings])
    err = np.array([t[2] for t in timings])
    w = 1.0 / err**2
    A = np.vstack([np.ones_like(n), n]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    coef = cov @ (A.T @ (w * d))
    resid = d - A @ coef
    chi2 = float((w * resid**2).sum())
    dof = max(len(n) - 2, 1)
    slope_err = np.sqrt(cov[1, 1])
    scale = max(1.0, np.sqrt(chi2 / dof))  # inflate if the scatter exceeds the errors
    deriv = spacing_derivative(sysinfo["P"], sysinfo["e"], np.deg2rad(sysinfo["omega0_deg"]))
    rate = np.degrees(coef[1] / deriv)
    rate_err = abs(np.degrees(slope_err * scale / deriv))
    return {"system": name, "n_points": int(len(n)), "slope_days_per_cycle": float(coef[1]),
            "slope_err": float(slope_err * scale), "chi2": chi2, "dof": dof,
            "wdot_deg_per_cycle": float(rate), "wdot_err": float(rate_err),
            "span_cycles": float(n.max() - n.min())}


def compactness(m_total: float, P_days: float, e: float) -> tuple[float, float]:
    """(GM/(c^2 a (1-e^2)), a in au) for a binary of total mass m_total (Msun)."""
    gm = GM_SUN_SI * m_total
    a = (gm * (P_days * DAY_S) ** 2 / (4 * np.pi**2)) ** (1 / 3)
    return gm / (C_SI**2 * a * (1 - e**2)), a / AU_M


def gr_rate_deg_per_cycle(m_total: float, P_days: float, e: float) -> float:
    """6 pi GM / (c^2 a (1-e^2)), in degrees per orbit (Baroch et al.'s Eq. 1)."""
    x, _ = compactness(m_total, P_days, e)
    return float(np.degrees(6 * np.pi * x))


def analyse() -> dict:
    tab = baroch_tables()
    systems = tab["systems"]
    out = {"systems": {}, "refits": {}}
    for name, s in systems.items():
        if "wdot_obs" not in s:
            continue
        x, a_au = compactness(s["m1"] + s["m2"], s["P"], s["e"])
        gr = gr_rate_deg_per_cycle(s["m1"] + s["m2"], s["P"], s["e"])
        meas = s["wdot_obs"] - s["wdot_classical"]
        meas_err = np.hypot(s["wdot_obs_err"], s["wdot_classical_err"])
        out["systems"][name] = {
            **s, "compactness": x, "a_au": a_au, "gr_ours": gr,
            "gr_measured_ours": meas, "gr_measured_err": meas_err,
            "ratio_ours": meas / gr, "ratio_err": meas_err / gr,
            "gr_share_pct": 100 * gr / s["wdot_theory"],
            "used_in_test": name in BAROCH_TEST_SYSTEMS,
        }
    for name in ("AL Dor", "V541 Cyg", "V459 Cas", "EY Cep"):
        fit = fit_apsidal_rate(name, systems[name], tab["timings"][name])
        s = out["systems"][name]
        fit["published"] = s["wdot_obs"]
        fit["published_err"] = s["wdot_obs_err"]
        fit["difference_sigma"] = (fit["wdot_deg_per_cycle"] - s["wdot_obs"]) / np.hypot(
            fit["wdot_err"], s["wdot_obs_err"])
        out["refits"][name] = fit
    # KX Cnc: the paper quotes the change of T2-T1 between Sowell et al.'s data and TESS.
    kx = systems["KX Cnc"]
    slope, slope_err = -0.00343 / 177.0, 0.00015 / 177.0
    deriv = spacing_derivative(kx["P"], kx["e"], np.deg2rad(kx["omega0_deg"]))
    out["refits"]["KX Cnc"] = {
        "system": "KX Cnc", "slope_days_per_cycle": slope, "slope_err": slope_err,
        "wdot_deg_per_cycle": float(np.degrees(slope / deriv)),
        "wdot_err": float(abs(np.degrees(slope_err / deriv))),
        "published": kx["wdot_obs"], "published_err": kx["wdot_obs_err"],
        "note": "from the quoted change of T2-T1 (-0.00343 +- 0.00015 d over 177 cycles)"}
    # The combined test: slope A through the origin of measured vs predicted GR rates.
    sel = [out["systems"][n] for n in BAROCH_TEST_SYSTEMS]
    xg = np.array([s["wdot_gr"] for s in sel])
    yg = np.array([s["gr_measured_ours"] for s in sel])
    sg = np.array([s["gr_measured_err"] for s in sel])
    w = 1 / sg**2
    A = float((w * xg * yg).sum() / (w * xg**2).sum())
    A_err = float(1 / np.sqrt((w * xg**2).sum()))
    out["combined_A"] = {"A": A, "err": A_err, "systems": list(BAROCH_TEST_SYSTEMS),
                         "published": PUBLISHED["baroch2021_A"]}
    gamma_ppn = PUBLISHED["cassini_gamma"]
    out["combined_A"]["beta_minus_1"] = 2 * gamma_ppn["gamma_minus_1"] - 3 * (A - 1)
    out["combined_A"]["beta_minus_1_err"] = 3 * A_err
    out["tess_timings"] = {k: [list(t) for t in v] for k, v in tab["timings"].items()}
    return out


def di_her() -> dict:
    """DI Herculis: the 'anomaly', and how tilted spins resolve it (published values)."""
    d = PUBLISHED["di_her"]
    gr = gr_rate_deg_per_cycle(d["m1"] + d["m2"], d["p_days"], d["e"])
    x, a_au = compactness(d["m1"] + d["m2"], d["p_days"], d["e"])
    latest = d["history"][3]  # Claret et al. 2021 (TESS + photoelectric)
    classical_mis = d["theory_misaligned_claret21"] - d["gr_claret10"]
    meas_gr = latest["obs"] - classical_mis
    meas_err = np.hypot(latest["err"], d["theory_misaligned_claret21_err"])
    cycles_per_century = 36525.0 / d["p_days"]
    return {"gr_ours": gr, "compactness": x, "a_au": a_au,
            "cycles_per_century": cycles_per_century,
            "gr_deg_per_century": gr * cycles_per_century,
            "aligned_prediction": d["theory_aligned_gm85"],
            "observed_latest": latest["obs"], "observed_latest_err": latest["err"],
            "gr_measured_misaligned": meas_gr, "gr_measured_err": meas_err,
            "ratio": meas_gr / gr, "ratio_err": meas_err / gr,
            "ratio_if_aligned": (latest["obs"] - d["classical_aligned_gm85"]) / gr}
