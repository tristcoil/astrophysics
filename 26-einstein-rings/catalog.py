"""Read the SLACS catalogues into tidy tables.

Three CDS tables carry the whole survey:

* **SLACS V** (Bolton et al. 2008a), table 4: all 131 candidates Hubble
  observed -- SDSS plate-MJD-fibre, the foreground (lens) and background
  (source) redshifts, de Vaucouleurs F814W photometry, the SDSS velocity
  dispersion, morphology and the lens grade (A = definite, B = probable,
  X = not a lens or not detected).
* **SLACS V**, table 5: singular-isothermal-ellipsoid models of 63 grade-A lenses.
* **SLACS IX** (Auger et al. 2009), tables 3 and 4: all 85 grade-A lenses,
  15 of them new, with multi-band photometry, Einstein radii in kpc and
  stellar masses for two choices of the stellar initial mass function.

All names are normalised to the short form ``J0946+1006``.
"""

from __future__ import annotations

import re
import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
from astropy.io import ascii

from config import TABLE_DIR

V_DIR = TABLE_DIR / "J_ApJ_682_964"
IX_DIR = TABLE_DIR / "J_ApJ_705_1099"


def _read(folder, table: str) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tab = ascii.read(folder / table, format="cds", readme=folder / "ReadMe")
    df = tab.to_pandas()
    return df


def short_name(name: str) -> str:
    """'SDSSJ0946+1006', 'J0946+1006' or '094656.68+100652.8' -> 'J0946+1006'."""
    name = name.strip().replace("SDSS", "").replace("−", "-")
    m = re.fullmatch(r"J?(\d{4})\d*(?:\.\d+)?([+-])(\d{4})\d*(?:\.\d+)?", name)
    if not m:
        raise ValueError(f"unrecognised SLACS name {name!r}")
    return f"J{m.group(1)}{m.group(2)}{m.group(3)}"


def _radec_from_sdss_name(name: str) -> tuple[float, float]:
    """'094656.68+100652.8' -> (146.73617, 10.11467) degrees."""
    m = re.fullmatch(r"(\d\d)(\d\d)(\d\d\.\d+)([+-])(\d\d)(\d\d)(\d\d\.\d+)", name.strip())
    if not m:
        raise ValueError(f"unrecognised SDSS name {name!r}")
    hh, mm, ss, sign, dd, am, asec = m.groups()
    ra = 15.0 * (int(hh) + int(mm) / 60 + float(ss) / 3600)
    dec = int(dd) + int(am) / 60 + float(asec) / 3600
    return ra, (-dec if sign == "-" else dec)


def _clean(series: pd.Series) -> pd.Series:
    """CDS masked values arrive as NaN or as masked strings; make them NaN."""
    return pd.to_numeric(series, errors="coerce").astype("float64")


@lru_cache(maxsize=None)
def read_candidates() -> pd.DataFrame:
    """All 131 SLACS candidates with HST imaging (SLACS V, table 4)."""
    raw = _read(V_DIR, "table4.dat")
    out = pd.DataFrame({
        "name": [short_name(n) for n in raw["Name"]],
        "sdss_name": raw["SDSS"].astype(str),
        "plate": raw["Plate"].astype(int),
        "mjd": raw["MJD"].astype(int),
        "fiber": raw["Fiber"].astype(int),
        "z_lens": raw["zFG"].astype(float),
        "z_source": raw["zBG"].astype(float),
        "i814": _clean(raw["Imag"]),
        "i814_origin": raw["n_Imag"].astype(str),  # s = 420-s snapshot, d = full orbit
        "a_i": _clean(raw["AI"]),
        "lum_v_1e9": _clean(raw["L(V555)"]),
        "re_arcsec": _clean(raw["Re"]),  # de Vaucouleurs, intermediate axis
        "light_frac_re2": _clean(raw["Le2/LdV"]),
        "q_light": _clean(raw["b/a"]),
        "pa_light": _clean(raw["PA"]),
        "sigma": _clean(raw["sigma"]),  # SDSS, within the 3" fibre
        "e_sigma": _clean(raw["e_sigma"]),
        "morphology": raw["Mph"].astype(str),  # E early, L late, U unclassified
        "multiplicity": raw["Mul"].astype(str),  # S single, M multiple
        "grade": raw["Lens"].astype(str),  # A, B or X
        "not_modelled_note": raw["n_Lens"].astype(str).eq("*"),
    })
    radec = [_radec_from_sdss_name(n) for n in out["sdss_name"]]
    out["ra"] = [r for r, _ in radec]
    out["dec"] = [d for _, d in radec]
    return out


@lru_cache(maxsize=None)
def read_bolton_models() -> pd.DataFrame:
    """SIE and light-traces-mass lens models of 63 grade-A lenses (SLACS V, table 5)."""
    raw = _read(V_DIR, "table5.dat")
    return pd.DataFrame({
        "name": [short_name(n) for n in raw["Name"]],
        "b_sie": raw["bSIE"].astype(float),  # Einstein radius, intermediate-axis normalisation
        "q_sie": raw["qSIE"].astype(float),
        "pa_sie": raw["PASIE"].astype(float),  # degrees east of north
        "light_frac_bsie": raw["LSIE/LdV"].astype(float),
        "b_ltm": raw["bLTM"].astype(float),
        "shear_ltm": raw["gammaLTM"].astype(float),
        "n_source_components": raw["Nsrc"].astype(int),
        "ring_subset": raw["Ring?"].astype(str).str.strip().eq("Yes"),
        "good_sigma": raw["Good?"].astype(str).str.strip().eq("Yes"),
    })


@lru_cache(maxsize=None)
def read_notes() -> pd.DataFrame:
    raw = _read(V_DIR, "notes.dat")
    return pd.DataFrame({"name": [short_name(n) for n in raw["Name"]],
                         "note": raw["Note"].astype(str)})


def _radec_hms(row) -> tuple[float, float]:
    ra = 15.0 * (row["RAh"] + row["RAm"] / 60 + row["RAs"] / 3600)
    dec = row["DEd"] + row["DEm"] / 60 + row["DEs"] / 3600
    return ra, (-dec if str(row["DE-"]).strip() == "-" else dec)


@lru_cache(maxsize=None)
def read_auger() -> pd.DataFrame:
    """All 85 grade-A lenses with stellar masses (SLACS IX, tables 3 and 4)."""
    t3 = _read(IX_DIR, "table3.dat")
    t4 = _read(IX_DIR, "table4.dat")
    radec = [_radec_hms(r) for _, r in t3.iterrows()]
    a = pd.DataFrame({
        "name": [short_name(n) for n in t3["SDSS"]],
        "ra": [r for r, _ in radec],
        "dec": [d for _, d in radec],
        "z_lens": t3["zlens"].astype(float),
        "z_source": t3["zsrc"].astype(float),
        "sigma": _clean(t3["sigma"]),
        "e_sigma": _clean(t3["e_sigma"]),
        "morph_type": t3["MType"].astype(str).str.strip(),
        "companion": t3["f_MType"].astype(str).str.strip().eq("c"),
        "b_mag": _clean(t3["Bmag"]), "re_b": _clean(t3["Re(B)"]),
        "v_mag": _clean(t3["Vmag"]), "v_filter": t3["f_Vmag"].astype(str), "re_v": _clean(t3["Re(V)"]),
        "i_mag": _clean(t3["Imag"]), "re_i": _clean(t3["Re(I)"]),
        "h_mag": _clean(t3["Hmag"]),
    })
    b = pd.DataFrame({
        "name": [short_name(n) for n in t4["SDSS"]],
        "r_ein_kpc": _clean(t4["RE"]),
        "log_m_ein": _clean(t4["Mass"]),
        "fstar_chab": _clean(t4["Fc"]), "e_fstar_chab": _clean(t4["e_Fc"]),
        "fstar_salp": _clean(t4["Fs"]), "e_fstar_salp": _clean(t4["e_Fs"]),
        "fstar_ratio_ein_eff": _clean(t4["F"]),
        "log_mstar_chab": _clean(t4["logMc"]), "e_log_mstar_chab": _clean(t4["e_logMc"]),
        "log_mstar_salp": _clean(t4["logMs"]), "e_log_mstar_salp": _clean(t4["e_logMs"]),
        "abs_b": _clean(t4["BMag"]), "abs_v": _clean(t4["VMag"]),
        "abs_v_z0": _clean(t4["VMag0"]),
    })
    return a.merge(b, on="name", how="outer", validate="one_to_one")


@lru_cache(maxsize=None)
def read_grade_a() -> pd.DataFrame:
    """The 85 definite lenses, one row each, merging SLACS V and IX.

    Positions, redshifts and dispersions come from SLACS IX (which covers all
    85); SLACS V adds the SDSS fibre, the de Vaucouleurs fit and, for 63, the
    published SIE lens model. Columns from SLACS V are NaN for the 15 lenses
    confirmed only in SLACS IX.
    """
    aug = read_auger()
    cand = read_candidates().drop(columns=["ra", "dec", "z_lens", "z_source", "sigma", "e_sigma"])
    cand = cand.rename(columns={"grade": "grade_slacs_v"})
    models = read_bolton_models()
    out = aug.merge(cand, on="name", how="left", validate="one_to_one")
    out = out.merge(models, on="name", how="left", validate="one_to_one")
    out["in_slacs_v"] = out["plate"].notna()
    out["has_sie_model"] = out["b_sie"].notna()
    return out.sort_values("name").reset_index(drop=True)


def summary() -> dict:
    cand = read_candidates()
    grade_a = read_grade_a()
    return {
        "n_candidates": len(cand),
        "n_grade": cand["grade"].value_counts().to_dict(),
        "n_grade_a_slacs_ix": len(grade_a),
        "n_sie_models_slacs_v": int(grade_a["has_sie_model"].sum()),
        "n_with_einstein_radius_slacs_ix": int(grade_a["r_ein_kpc"].notna().sum()),
        "n_with_stellar_mass": int(grade_a["log_mstar_chab"].notna().sum()),
    }


if __name__ == "__main__":
    print(summary())
    g = read_grade_a()
    print(g[["name", "z_lens", "z_source", "sigma", "b_sie", "r_ein_kpc", "grade_slacs_v",
             "morph_type"]].to_string())
    new = g[~g["in_slacs_v"] | (g["grade_slacs_v"] != "A")]
    print("not grade A in SLACS V:", list(new["name"]), list(new["grade_slacs_v"]))
    print(np.nanmax(np.abs(g["z_lens"] - read_candidates().set_index("name").reindex(g["name"])["z_lens"].values)))
