"""Build the interactive page "Reduce the plates yourself": one self-contained HTML file.

    python explorer.py      # writes figures/plate_lab.html and data/derived/plate_lab_artifact.html

The page carries the measures of 1919 and 1922 as transcribed here, and refits
them in the browser while you choose which plates, stars and coordinates to
use, or how to treat the plate scale: the Sobral 4-inch (joint least squares
over 14 plates), the Sobral astrographic (16 plates, five stars, declinations),
Principe (four comparisons, both coordinates) and Wallal 1922 (radial
displacements against distance). The same fits in Python
(eclipse1919.py, principe.py, wallal.py) are embedded as reference values, and
the page's own numbers must match them (checked in a headless browser by
tools/check_explorer.js). So are the jackknife errors of the 4-inch and
astrographic fits with every plate selected, which the result cards quote.

The artifact copy is the same page as a fragment for a viewer that supplies
the document skeleton and the light/dark theme (no theme buttons).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # works under python -I too

import numpy as np  # noqa: E402

import config as C  # noqa: E402

TEMPLATE = C.PROJECT_DIR / "explorer_template.html"
OUT = C.FIGURE_DIR / "plate_lab.html"
DERIVED_DIR = C.DATA_DIR / "derived"
OUT_ARTIFACT = DERIVED_DIR / "plate_lab_artifact.html"
DATA_JSON = DERIVED_DIR / "plate_lab_data.json"
THEME_BUTTONS = """      <div class="seg themebtn" role="group" aria-label="Theme">
        <button type="button" data-theme-set="light" aria-pressed="false">Light</button>
        <button type="button" data-theme-set="dark" aria-pressed="false">Dark</button>
      </div>"""


def _r(v, n=5):
    if isinstance(v, (list, tuple, np.ndarray)):
        return [_r(x, n) for x in v]
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    return round(float(v), n)


def field() -> dict:
    import stars as S

    t1 = S.table1().set_index("no")
    f = S.fit_table1_frame()
    out = {"stars": [], "sun": {}}
    for no, r in t1.iterrows():
        out["stars"].append({"no": int(no), "name_1920": r["name"], "name": C.TABLE1_GAIA[no][1],
                             "mag": float(r["photog_mag"]), "x": float(r["x"]), "y": float(r["y"])})
    for site in ("sobral", "principe"):
        pat = S.predicted_pattern(f, site)
        out["sun"][site] = {"x": _r(pat["sun_x"].iloc[0]), "y": _r(pat["sun_y"].iloc[0]),
                            "radius_units": _r(pat["sun_radius_arcmin"].iloc[0] / C.UNIT_TABLE1_ARCMIN)}
        # GR's exact displacement (arcsec) of each star at that station
        out["sun"][site]["gr_dx"] = _r(pat["dx"].tolist(), 4)
        out["sun"][site]["gr_dy"] = _r(pat["dy"].tolist(), 4)
    return out


def inch4() -> dict:
    import eclipse1919 as E

    d = E.load_4inch()
    dp = E.load_4inch(as_printed=True)
    co = d["coeff"]
    plates = E.ECLIPSE_4IN + E.COMPARISON_4IN
    ref = {}
    for key, coords in (("both", ("RA", "Dec")), ("RA", ("RA",)), ("Dec", ("Dec",))):
        f = E.joint_fit_4inch(d, coords=coords)
        ref[key] = {"alpha50": _r(f["alpha"], 6), "alpha50_err": _r(f["alpha_err"], 6),
                    "sigma_by_group": {k: _r(v, 6) for k, v in f["sigma_by_group"].items()}}
    f = E.joint_fit_4inch(dp)
    ref["both_as_printed"] = {"alpha50": _r(f["alpha"], 6), "alpha50_err": _r(f["alpha_err"], 6)}
    for p in ("I", "15_1"):
        f = E.joint_fit_4inch(d, drop_plates=(p,))
        ref[f"drop_plate_{p}"] = {"alpha50": _r(f["alpha"], 6), "alpha50_err": _r(f["alpha_err"], 6)}
    for s in (11, 2):
        f = E.joint_fit_4inch(d, drop_stars=(s,))
        ref[f"drop_star_{s}"] = {"alpha50": _r(f["alpha"], 6), "alpha50_err": _r(f["alpha_err"], 6)}
    # jackknife of the full selection, one plate or one star left out at a time (as analysis.step_sobral_4inch);
    # the write-up adopts the largest of the formal error and these two
    jk = E.jackknife_4inch(d)
    jk_err = {}
    for what in ("plate", "star"):
        a = jk.loc[jk["left_out"].str.startswith(what), "alpha"].to_numpy()
        jk_err[what] = float(np.sqrt((len(a) - 1) / len(a) * np.sum((a - a.mean()) ** 2)))
    rep = E.replicate_4inch(d)
    return {"stars": E.STARS_4IN, "x": _r(co["x"].tolist()), "y": _r(co["y"].tolist()),
            "Ex": _r(co["Ex"].tolist()), "Ey": _r(co["Ey"].tolist()),
            "plates": plates, "eclipse": E.ECLIPSE_4IN,
            "Dx": [_r(d["Dx_used"][p].tolist(), 4) for p in plates],          # rev, per plate, per star
            "Dy": [_r(d["Dy_used"][p].tolist(), 4) for p in plates],
            "Dx_as_printed": [_r(dp["Dx_used"][p].tolist(), 4) for p in plates],
            "rev_arcsec": C.REV_4INCH,
            "per_plate_1920": rep["solutions"][["coord", "plate", "alpha"]].assign(alpha=lambda t: t["alpha"].round(4))
            .to_dict("records"),
            "published": {"limb": 1.98, "pe": 0.12, "ra": 2.06, "dec": 1.94}, "reference": ref,
            "jackknife": {"alpha50_err_plates": _r(jk_err["plate"], 6), "alpha50_err_stars": _r(jk_err["star"], 6)}}


def astro() -> dict:
    import eclipse1919 as E

    d = E.load_astro()
    co = d["coeff"]
    ref = {}
    for key, kw in (("common", {"scale": "common"}), ("per_plate", {"scale": "per_plate"}),
                    ("fixed_0.020", {"scale": "fixed", "e_fixed": 0.020}), ("fixed_0.082", {"scale": "fixed", "e_fixed": 0.082})):
        f = E.fit_astro(d, **kw)
        ref[key] = {"alpha50": _r(f["alpha"], 6), "alpha50_err": _r(f["alpha_err"], 6)}
        if "e_rev" in f:
            ref[key].update({"e_rev": _r(f["e_rev"], 6), "e_err_rev": _r(f["e_err_rev"], 6)})
    # jackknife over the 16 plates of the fit with one common e (as analysis.step_astrographic); the write-up
    # adopts the larger of it and the formal error
    a = np.array([E.fit_astro({**d, "Dy": d["Dy"].drop(index=p)}, scale="common")["alpha"] for p in d["Dy"].index])
    jk_plates = float(np.sqrt((len(a) - 1) / len(a) * np.sum((a - a.mean()) ** 2)))
    t9 = d["t9"]
    return {"plates": list(d["Dy"].index), "stars": E.STARS_ASTRO, "x": _r(co["x"].tolist()), "y": _r(co["y"].tolist()),
            "Ey": _r(co["Ey"].tolist()), "Dy": [_r(row.tolist(), 4) for _, row in d["Dy"].iterrows()],
            "rev_arcsec": C.REV_ASTRO,
            "table9": [{"plate": r.eclipse_plate, "comparison": r.comparison_plate, "n_stars": int(r.n_stars),
                        "e_rev": _r(r.e_rev), "limb": _r(r.alpha_limb_arcsec)} for r in t9.itertuples()],
            "published": {"limb_plates_scale": 0.93, "limb_refraction_scale": 1.52, "e_plates": 0.082, "e_refraction": 0.020,
                          "harvey": [1.55, 0.34]},
            "reference": ref, "jackknife": {"alpha50_err_plates": _r(jk_plates, 6)}}


def principe() -> dict:
    import principe as PR

    d = PR.load()
    comps = {}
    for c in PR.COMPARISONS:
        t = PR._comparison(d, c)
        gr = PR.gr_pattern(t["star"].tolist())
        comps[c] = {"plate": PR.PLATE_OF[c], "stars": [int(s) for s in t["star"]], "x": _r(t["x"].tolist()),
                    "y": _r(t["y"].tolist()), "d1x": _r(t["d1x"].tolist(), 4), "d1y": _r(t["d1y"].tolist(), 4),
                    "Ex": _r(gr["Ex"].tolist(), 4), "Ey": _r(gr["Ey"].tolist(), 4)}
    ref = {}
    for key, kw in (("check_scale", {}), ("free_scale", {"free_scale": True})):
        f = PR.fit(d, **kw)
        ref[key] = {"limb": _r(f["limb"], 6), "limb_err": _r(f["limb_err"], 6),
                    "limb_err_measures": _r(f["limb_err_measures"], 6),
                    "per_comparison": {c: _r(v["limb"], 6) for c, v in f["per_comparison"].items()}}
    return {"comparisons": comps, "unit_arcsec": PR.UNIT, "gr_limb": C.DEFLECTION_LIMB_GR,
            "scale_err_limb_1sigma": PR.SCALE_PE_LIMB / 0.6745,
            "published": {"limb": 1.61, "pe": 0.30, "per_comparison": {"X-G1": 1.94, "X-H1": 1.44, "W-D1": 1.55, "W-I2": 1.67}},
            "reference": ref}


def wallal() -> dict | None:
    if not (C.TRANSCRIBED_DIR / "wallal1928_table6_radial.csv").exists():
        return None
    import wallal as W

    t6 = W.table6()
    pr = W.pairs(t6)
    fit = W.fit_5ft(t6)
    pw = W.fit_5ft(t6, n_free=True)
    rep = W.replicate(t6)
    return {"model": "D_r = beta (R/r)^n - alpha_p r for each plate pair p; R in degrees; weights w; "
                     "one-sigma of a measure of weight 1 = sigma_unit; n = 1 unless fitted",
            "R_deg": W.R_SUN_DEG, "sigma_unit": W.PE_UNIT_5FT / W.PE,
            "r_deg": _r(t6["distance_deg"].tolist(), 4), "Dr_mean": _r(t6["Dr_mean"].tolist(), 3),
            "wt_total": _r(t6["wt_total"].tolist(), 3), "pairs": pr,
            "Dr": {p: _r(t6[f"Dr_{p}"].tolist(), 3) for p in pr}, "wt": {p: _r(t6[f"wt_{p}"].tolist(), 3) for p in pr},
            "alpha_mean": _r(rep["alpha"], 6),
            "published": {"five_foot": [1.82, 0.15], "fifteen_foot": [1.72, 0.11], "combined": [1.75, 0.09],
                          "per_pair": {"M15_M7": 2.38, "M16_M11": 1.66, "M17_M9": 1.36, "N18_N11": 2.00,
                                       "N19_N13": 2.54, "N20_N9": 1.01}, "error_type": "probable error"},
            "reference": {"fit": {"beta": _r(fit["beta"], 6), "beta_err": _r(fit["beta_err"], 6),
                                  "chi2": _r(fit["chi2"], 3), "dof": fit["dof"]},
                          "power": {"beta": _r(pw["beta"], 6), "n": _r(pw["n"], 6), "n_err": _r(pw["n_err"], 6)},
                          "per_pair_beta": {r["pair"]: _r(r["beta"], 6) for r in rep["per_pair"].to_dict("records")}}}


def data() -> dict:
    return {"meta": {"gr_limb": C.DEFLECTION_LIMB_GR, "limb_1919": C.DEFLECTION_LIMB_1919,
                     "sun_radius_true_arcmin": C.SUN_RADIUS_TRUE_ARCMIN, "unit_arcmin": C.UNIT_TABLE1_ARCMIN},
            "field": field(), "inch4": inch4(), "astro": astro(), "principe": principe(), "wallal": wallal()}


def build() -> list[Path]:
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    d = data()
    DATA_JSON.write_text(json.dumps(d, indent=1))
    if not TEMPLATE.exists():
        return [DATA_JSON]
    tpl = TEMPLATE.read_text()
    payload = json.dumps(d, separators=(",", ":"))
    fragment = tpl.replace("/*__DATA__*/null", payload)
    # standalone page: a full document (standards mode), with theme buttons
    page = fragment.replace("<!--__THEME_BUTTONS__-->", THEME_BUTTONS)
    head, body = page.split("<main>", 1)
    OUT.write_text('<!doctype html>\n<html lang="en">\n<head>\n' + head + "</head>\n<body>\n<main>" + body
                   + "</body>\n</html>\n")
    # artifact: the fragment, for a viewer that supplies the document and the theme
    OUT_ARTIFACT.write_text(fragment.replace("<!--__THEME_BUTTONS__-->\n", "").replace("<!--__THEME_BUTTONS__-->", ""))
    return [OUT, OUT_ARTIFACT, DATA_JSON]


def main() -> int:
    for p in build():
        print(f"  {p.name} ({p.stat().st_size / 1e3:.0f} kB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
