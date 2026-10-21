"""The Sun's gravitational redshift, measured in its own iron lines.

Light leaving the Sun's surface climbs out of its potential: GR predicts its
frequency drops by GM/(R c^2), 636.3 m/s in velocity units at the surface and
633.1 m/s as received on Earth (after the small fall into the potential at
1 AU and the Earth's own; Gonzalez Hernandez et al. 2020).

The complication is convection. Hot, bright gas rises in granules and cool, dim
gas sinks between them, so a disc-integrated line is *blueshifted* by a few
hundred m/s, and by different amounts for different lines: weak lines form deep
where the granulation is vigorous, the cores of the strongest lines form high
in the photosphere where it has died away. Plotting each line's core shift
against its depth gives the "third signature" of granulation; the strongest
lines should sit at the bare gravitational redshift.

Data: the IAG solar flux atlas (Reiners et al. 2016, Fourier-transform
spectrometer in Goettingen, absolute wavelength scale good to +-10 m/s,
corrected for the observer's motion), laboratory iron wavelengths from Nave et
al. (1994, Fourier-transform hollow-cathode spectra), and the list of 1,251
clean Fe I lines with depths and equivalent widths that Reiners et al. used
(Table A1, after Allende Prieto & Garcia Lopez 1998).
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd

import config as C

CDS = "https://cdsarc.cds.unistra.fr/ftp/"
FILES = {
    "spvis.dat.gz": CDS + "J/A+A/587/A65/spvis.dat.gz",      # IAG atlas, 405-1065 nm (65 MB)
    "iag_tablea1.dat": CDS + "J/A+A/587/A65/tablea1.dat",    # 1,251 Fe I lines (depth, EW, blueshift)
    "iag_ReadMe": CDS + "J/A+A/587/A65/ReadMe",
    "nave1994_table2.gz": CDS + "J/ApJS/94/221/table2.gz",   # Fe I multiplet table, 9,759 lines
    "nave1994_ReadMe": CDS + "J/ApJS/94/221/ReadMe",
    "gh2020_tablea1.dat": CDS + "J/A+A/643/A146/tablea1.dat",  # HARPS-LFC Moon core shifts, 188 lines
    "gh2020_tablea2.dat": CDS + "J/A+A/643/A146/tablea2.dat",
    "gh2020_ReadMe": CDS + "J/A+A/643/A146/ReadMe",
}

#: Nave et al.'s wavenumber scale was tied to Norlen's (1973) Ar II standards; Whaling
#: et al. (1995) remeasured them, and NIST (Kramida et al. 2019) raises Fe I
#: wavenumbers by 6.7 +- 0.8 parts in 1e8: laboratory wavelengths get shorter and
#: every solar line shift grows by 6.7e-8 c = 20 m/s (Gonzalez Hernandez et al. 2020).
LAB_SCALE = 6.7e-8
LAB_SCALE_ERR = 0.8e-8
CORE_HALF_WINDOW_KMS = 1.6     # five HARPS pixels around the core (Gonzalez Hernandez et al. 2020)

V_GRS_SURFACE = C.GM_SUN / (C.R_SUN * C.C_LIGHT)                     # 636.3 m/s
V_GRS_EARTH = V_GRS_SURFACE - C.GM_SUN / (C.AU * C.C_LIGHT) - 0.21   # 633.1 m/s (Earth's own: 0.21)


def fetch() -> None:
    import fetch_data
    for name, url in FILES.items():
        fetch_data.fetch_file(url, C.SOLAR_DIR / name, note="CDS VizieR")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def load_atlas() -> tuple[np.ndarray, np.ndarray]:
    """IAG visible atlas: vacuum wavenumber (cm^-1, increasing) and normalised flux."""
    cache = C.CACHE_DIR / "iag_vis.npz"
    if cache.exists():
        with np.load(cache) as z:
            return z["sigma"], z["flux"]
    with gzip.open(C.SOLAR_DIR / "spvis.dat.gz", "rt") as fh:
        arr = pd.read_csv(fh, sep=r"\s+", header=None, usecols=[0, 1], names=["sigma", "flux"],
                          dtype=np.float64, engine="c").to_numpy()
    order = np.argsort(arr[:, 0])
    sigma, flux = arr[order, 0], arr[order, 1]
    np.savez(cache, sigma=sigma, flux=flux)
    return sigma, flux


def load_nave() -> pd.DataFrame:
    rows = []
    with gzip.open(C.SOLAR_DIR / "nave1994_table2.gz", "rt") as fh:
        for line in fh:
            lvac = line[37:47].strip()
            if not lvac or "." not in lvac:
                continue
            decimals = len(lvac.split(".")[1])
            rows.append({"lvac": float(lvac), "lair": float(line[48:58]) if line[48:58].strip() else np.nan,
                         "sigma": float(line[65:75]) if line[65:75].strip() else np.nan,
                         "quality": line[76:77].strip(), "ft": decimals == 4,
                         "intensity": line[30:36].strip(), "elow_ev": float(line[84:89] or "nan"),
                         "loggf": float(line[116:122]) if line[116:122].strip() else np.nan,
                         "blend": line[126:131].strip()})
    return pd.DataFrame(rows)


def _num(s: str) -> float:
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return np.nan


def load_iag_table() -> pd.DataFrame:
    """Reiners et al.'s Table A1: wavelength (nm, rounded to 1 mA), EW (mA), depth, blueshift (m/s)."""
    rows = []
    with open(C.SOLAR_DIR / "iag_tablea1.dat") as fh:
        for line in fh:
            if line.strip():
                rows.append({"lam_nm": float(line[0:8]), "ew_mA": _num(line[11:15]),
                             "depth": _num(line[18:22]), "convbs": _num(line[25:32])})
    return pd.DataFrame(rows)


def match_lines(iag: pd.DataFrame, nave: pd.DataFrame) -> pd.DataFrame:
    """Attach Nave et al.'s vacuum wavelength to each IAG line (the IAG list rounds them to 1 mA).

    The IAG list gives *vacuum* wavelengths (all 1,251 lie within 1.5 mA of a Nave
    vacuum wavelength, median 0.3 mA). A line is kept only if exactly one
    FT-measured, unblended Nave line of quality A or B lies within 1.5 mA.
    """
    good = nave[nave["ft"] & (nave["blend"] == "") & nave["quality"].isin(["A", "B"])].sort_values("lvac")
    lv = good["lvac"].to_numpy()
    out = []
    for r in iag.itertuples(index=False):
        target = r.lam_nm * 10
        j = np.searchsorted(lv, target)
        cand = [k for k in (j - 2, j - 1, j, j + 1) if 0 <= k < len(lv) and abs(lv[k] - target) < 0.0015]
        if len(cand) != 1:
            continue
        g = good.iloc[cand[0]]
        out.append({**r._asdict(), "lvac": g["lvac"], "quality": g["quality"],
                    "loggf": g["loggf"], "elow_ev": g["elow_ev"]})
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Measuring line cores
# ---------------------------------------------------------------------------

def core_shift(sigma: np.ndarray, flux: np.ndarray, lvac_A: float, *,
               half_window_kms: float = CORE_HALF_WINDOW_KMS, search_kms: float = 3.0,
               iterations: int = 3, lab_scale: float = LAB_SCALE) -> dict | None:
    """Velocity (m/s) of a line's core: the minimum of a parabola fitted to the central
    +-half_window (default +-1.6 km/s, like the five-pixel HARPS fit of Gonzalez
    Hernandez et al. 2020; Reiners et al. 2016 used a 4 km/s = 50-90 mA wide window),
    re-centred on the fitted minimum until it settles. The laboratory wavenumber is
    put on the NIST scale (lab_scale; 0 keeps Nave et al.'s original)."""
    c = C.C_LIGHT
    s_lab = 1e8 / lvac_A * (1 + lab_scale)
    # observed wavenumber = s_lab / (1 + v/c); search around the redshifted position
    s_exp = s_lab / (1 + 400.0 / c)
    lo, hi = s_exp * (1 - search_kms * 1e3 / c), s_exp * (1 + search_kms * 1e3 / c)
    i0, i1 = np.searchsorted(sigma, [lo, hi])
    if i1 - i0 < 5:
        return None
    k = i0 + int(np.argmin(flux[i0:i1]))
    if k in (i0, i1 - 1):          # minimum at the edge of the search box: blended or misidentified
        return None
    w = s_lab * half_window_kms * 1e3 / c
    centre = sigma[k]
    for _ in range(iterations):
        j0, j1 = np.searchsorted(sigma, [centre - w, centre + w])
        x = sigma[j0:j1] - centre
        y = flux[j0:j1]
        if len(x) < 7:
            return None
        p = np.polyfit(x, y, 2)
        if p[0] <= 0:
            return None
        x_min = -p[1] / (2 * p[0])
        if abs(x_min) > w:
            return None
        centre = centre + x_min
    v = c * (s_lab - centre) / centre
    return {"v_core": float(v), "depth_meas": float(1 - np.polyval(p, x_min)), "n_pts": int(len(x))}


def measure(lines: pd.DataFrame, sigma: np.ndarray, flux: np.ndarray, **kw) -> pd.DataFrame:
    """core_shift for every line (keyword arguments passed on)."""
    rows = []
    for r in lines.itertuples(index=False):
        m = core_shift(sigma, flux, r.lvac, **kw)
        if m is None:
            continue
        rows.append({**r._asdict(), **m})
    return pd.DataFrame(rows)


def strong_line_estimate(df: pd.DataFrame, ew_lo: float = 150.0, ew_hi: float = 550.0) -> dict:
    """Mean core shift of the strongest lines, whose cores form above the convection
    (Gonzalez Hernandez et al. 2020 used 150 < EW < 550 mA)."""
    s = df[(df["ew_mA"] > ew_lo) & (df["ew_mA"] < ew_hi)]["v_core"].to_numpy()
    med = np.median(s)
    rs = 1.4826 * np.median(np.abs(s - med))
    keep = np.abs(s - med) < 4 * max(rs, 20.0)
    s = s[keep]
    return {"n_lines": int(len(s)), "mean": float(np.mean(s)), "sem": float(np.std(s, ddof=1) / np.sqrt(len(s))),
            "sd": float(np.std(s, ddof=1)), "median": float(np.median(s)), "rejected": int((~keep).sum())}


def binned(df: pd.DataFrame, col: str = "depth", width: float = 0.05) -> pd.DataFrame:
    edges = np.arange(0.0, 1.0 + width, width)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = df[(df[col] >= lo) & (df[col] < hi)]["v_core"]
        if len(s) < 5:
            continue
        rows.append({"center": (lo + hi) / 2, "median": float(s.median()), "n": int(len(s)),
                     "sem": float(1.2533 * s.std(ddof=1) / np.sqrt(len(s)))})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Like-for-like replication of Gonzalez Hernandez et al. (2020)
# ---------------------------------------------------------------------------

HARPS_R = 115_000
HARPS_PIXEL_KMS = 0.82


def load_gh2020() -> pd.DataFrame:
    """Their Table A.1: per-line core shifts on HARPS-LFC Moon spectra (vs Nave's original wavelengths)."""
    rows = []
    with open(C.SOLAR_DIR / "gh2020_tablea1.dat") as fh:
        for line in fh:
            if not line.strip():
                continue
            rows.append({"lab": float(line[0:9]), "nist": _num(line[10:19]), "el": line[38:40].strip(),
                         "ion": int(line[41:42]), "ew_gh": _num(line[54:59]), "v_gh": _num(line[65:70]),
                         "e_v_gh": _num(line[71:75])})
    return pd.DataFrame(rows)


def harps_like_core(sigma: np.ndarray, flux: np.ndarray, lvac_A: float, *, lab_scale: float = 0.0,
                    phases: int = 8) -> float | None:
    """Core shift (m/s) after blurring the atlas to HARPS resolution and fitting a parabola
    to five 0.82 km/s pixels around the minimum, averaged over pixel phases."""
    c = C.C_LIGHT / 1e3                            # km/s
    s_lab = 1e8 / lvac_A * (1 + lab_scale)
    u = np.arange(-15.0, 15.0001, 0.02)            # km/s, velocity grid; v > 0 = redshift
    s_grid = s_lab / (1 + u / c)
    order = np.argsort(s_grid)
    i0, i1 = np.searchsorted(sigma, [s_grid.min() - 0.05, s_grid.max() + 0.05])
    if i1 - i0 < 50:
        return None
    f = np.interp(s_grid[order], sigma[i0:i1], flux[i0:i1])
    fu = np.empty_like(f)
    fu[order] = f
    sd = c / HARPS_R / 2.3548                      # Gaussian sigma of the instrumental profile
    k = np.arange(-int(4 * sd / 0.02), int(4 * sd / 0.02) + 1) * 0.02
    kern = np.exp(-0.5 * (k / sd) ** 2)
    blurred = np.convolve(fu, kern / kern.sum(), mode="same")
    centres = []
    for ph in np.arange(phases) / phases * HARPS_PIXEL_KMS:
        px = np.arange(-12.0, 12.0, HARPS_PIXEL_KMS) + ph
        fp = np.interp(px, u, blurred)
        j = int(np.argmin(np.where(np.abs(px) < 3.0, fp, np.inf)))
        if j < 2 or j > len(px) - 3:
            return None
        p = np.polyfit(px[j - 2:j + 3], fp[j - 2:j + 3], 2)
        if p[0] <= 0:
            return None
        centres.append(-p[1] / (2 * p[0]))
    return float(np.mean(centres) * 1e3)


def compare_with_gh2020(sigma: np.ndarray, flux: np.ndarray) -> pd.DataFrame:
    gh = load_gh2020()
    gh = gh[(gh["el"] == "Fe") & (gh["ion"] == 1) & np.isfinite(gh["v_gh"])]
    nave = load_nave()
    lv = nave.sort_values("lvac")
    out = []
    for r in gh.itertuples(index=False):
        # their lambda_lab is Nave's air or vacuum wavelength? match on either
        for col in ("lair", "lvac"):
            arr = lv[col].to_numpy()
            j = int(np.nanargmin(np.abs(arr - r.lab)))
            if abs(arr[j] - r.lab) < 0.00015:
                lvac = lv["lvac"].to_numpy()[j]
                break
        else:
            continue
        v_h = harps_like_core(sigma, flux, lvac)
        full = core_shift(sigma, flux, lvac, lab_scale=0.0)
        if v_h is None or full is None:
            continue
        out.append({**r._asdict(), "lvac": lvac, "v_iag_harpslike": v_h, "v_iag_full": full["v_core"]})
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Everything at once (cached)
# ---------------------------------------------------------------------------

WINDOWS_KMS = (1.0, 1.3, 1.6, 2.0)
IAG_SCALE_ERR = 10.0       # m/s, Reiners et al. 2016's stated accuracy of the IAG wavelength scale


def analyse(force: bool = False) -> dict:
    import json
    path = C.CACHE_DIR / "solar_results.json"
    if path.exists() and not force:
        return json.loads(path.read_text())
    nave, iag = load_nave(), load_iag_table()
    lines = match_lines(iag, nave)
    sigma, flux = load_atlas()
    per_window = []
    for hw in WINDOWS_KMS:
        m = measure(lines, sigma, flux, half_window_kms=hw)
        if hw == CORE_HALF_WINDOW_KMS:
            m.to_parquet(C.CACHE_DIR / "solar_lines.parquet")
            primary = m
        s = strong_line_estimate(m)
        g = strong_line_estimate(m[(m["lvac"] > 4760) & (m["lvac"] < 5850)])
        per_window.append({"half_window_kms": hw, "mean": s["mean"], "sem": s["sem"], "n": s["n_lines"],
                           "mean_476_585": g["mean"], "sem_476_585": g["sem"], "n_476_585": g["n_lines"]})
    pw = pd.DataFrame(per_window)
    pw.to_csv(C.CACHE_DIR / "solar_windows.csv", index=False)
    cmp = compare_with_gh2020(sigma, flux)
    cmp.to_parquet(C.CACHE_DIR / "solar_gh2020_compare.parquet")
    strong_gh = cmp[(cmp["ew_gh"] > 150) & (cmp["ew_gh"] < 550)]
    shift_nist = LAB_SCALE * C.C_LIGHT
    prim = strong_line_estimate(primary)
    spread = 0.5 * (pw["mean"].max() - pw["mean"].min())
    # check against Reiners et al.'s own numbers with their recipe (4 km/s wide window, Nave's scale)
    rm = measure(lines, sigma, flux, half_window_kms=2.0, lab_scale=0.0).dropna(subset=["convbs"])
    diff = rm["v_core"] - (rm["convbs"] + 636.0)
    out = {
        "theory_earth": V_GRS_EARTH, "theory_surface": V_GRS_SURFACE,
        "n_lines_measured": int(len(primary)), "n_strong": prim["n_lines"],
        "strong_mean": prim["mean"], "strong_sem": prim["sem"], "strong_sd": prim["sd"],
        "window_spread": float(spread), "iag_scale_err": IAG_SCALE_ERR,
        "lab_scale_err": LAB_SCALE_ERR * C.C_LIGHT,
        "strong_total_err": float(np.sqrt(prim["sem"] ** 2 + spread ** 2 + IAG_SCALE_ERR ** 2
                                          + (LAB_SCALE_ERR * C.C_LIGHT) ** 2)),
        "windows": per_window,
        "gh_compare_n": int(len(cmp)),
        "iag_minus_harps_fullres": float((cmp["v_iag_full"] - cmp["v_gh"]).mean()),
        "iag_minus_harps_fullres_sem": float((cmp["v_iag_full"] - cmp["v_gh"]).std() / np.sqrt(len(cmp))),
        "iag_minus_harps_harpslike": float((cmp["v_iag_harpslike"] - cmp["v_gh"]).mean()),
        "iag_minus_harps_harpslike_sem": float((cmp["v_iag_harpslike"] - cmp["v_gh"]).std() / np.sqrt(len(cmp))),
        "gh_lines_n": int(len(strong_gh)),
        "gh_lines_harps": float(strong_gh["v_gh"].mean() + shift_nist),
        "gh_lines_harps_err": float(strong_gh["v_gh"].std() / np.sqrt(len(strong_gh))),
        "gh_lines_iag": float(strong_gh["v_iag_full"].mean() + shift_nist),
        "gh_lines_iag_err": float(strong_gh["v_iag_full"].std() / np.sqrt(len(strong_gh))),
        "reiners_check_median": float(diff.median()),
        "reiners_check_robust_sd": float(1.4826 * np.median(np.abs(diff - diff.median()))),
        "binned": binned(primary).to_dict("records"),
    }
    path.write_text(json.dumps(out, indent=1, default=float))
    return out
