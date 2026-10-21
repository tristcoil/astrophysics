"""Neutron-capture elements: barium (mostly slow, s-process) and europium
(mostly rapid, r-process).

* Trends of [Ba/Fe], [Eu/Fe] and [Eu/Mg] with [Fe/H] in the disc and in
  Gaia-Sausage-Enceladus, from GALAH DR4 giants (Eu from Kane et al.'s VAC).
* [Ba/Eu]: how pure the r-process was. Solar r-process matter alone has
  [Ba/Eu] = log10[(1 - s_Ba)/(1 - s_Eu)] = -0.8 (Bisterzo et al. 2014's s
  fractions); AGB stars push it up toward solar as they start to contribute.
* Scatter: rare, prolific sources (one merger enriching a whole gas cloud)
  leave star-to-star scatter wherever few events have been mixed in -- at low
  metallicity -- while common sources (core-collapse supernovae for Mg)
  average out early. GALAH reaches [Fe/H] ~ -1.7 for Eu; below that we use
  JINAbase, a literature compilation of metal-poor stars.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from config import LIT_DIR

#: Asplund et al. (2009) photospheric abundances, the scale JINAbase uses.
SOLAR_JINA = {"Fe": 7.50, "Mg": 7.60, "Ba": 2.18, "Eu": 0.52}


def robust_sigma(y):
    y = y[np.isfinite(y)]
    if y.size < 5:
        return np.nan
    return float(1.4826 * np.median(np.abs(y - np.median(y))))


def binned(x, y, edges, *, err=None, min_n=15):
    """Median, 16/84th percentiles, robust scatter and intrinsic scatter
    (robust scatter with the median measurement error removed) per bin."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (x >= lo) & (x < hi) & np.isfinite(y)
        if s.sum() < min_n:
            continue
        rs = robust_sigma(y[s])
        e = float(np.median(err[s])) if err is not None else 0.0
        rows.append(dict(x=float(np.median(x[s])), lo=float(lo), hi=float(hi), n=int(s.sum()),
                         median=float(np.median(y[s])), p16=float(np.percentile(y[s], 16)),
                         p84=float(np.percentile(y[s], 84)), robust_sigma=rs, median_err=e,
                         intrinsic=float(np.sqrt(max(rs ** 2 - e ** 2, 0.0)))))
    return rows


def offset_at_fixed_feh(feh, y, ref_feh, ref_y, edges, *, n_boot=300, seed=0):
    """Median of a population minus the reference population's median in the
    same [Fe/H] bins, star-number weighted, with a bootstrap error."""
    rng = np.random.default_rng(seed)

    def one(f, v, rf, rv):
        diffs, ws = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            a = (f >= lo) & (f < hi) & np.isfinite(v)
            b = (rf >= lo) & (rf < hi) & np.isfinite(rv)
            if a.sum() >= 5 and b.sum() >= 20:
                diffs.append(np.median(v[a]) - np.median(rv[b]))
                ws.append(a.sum())
        return np.average(diffs, weights=ws) if diffs else np.nan

    best = one(feh, y, ref_feh, ref_y)
    boots = []
    for _ in range(n_boot):
        i = rng.integers(0, feh.size, feh.size)
        j = rng.integers(0, ref_feh.size, ref_feh.size)
        boots.append(one(feh[i], y[i], ref_feh[j], ref_y[j]))
    return float(best), float(np.nanstd(boots)), int(np.isfinite(y).sum())


# ---------------------------------------------------------------------------
# JINAbase
# ---------------------------------------------------------------------------


def load_jinabase() -> dict:
    """One row per star (the compilation's preferred entry where a star has
    several), Milky Way field stars only (no dwarf galaxies or bulge), and no
    CEMP-s stars (their barium came from a binary companion, not the gas the
    star formed from)."""
    d = pd.read_csv(LIT_DIR / "JINAbase_2021.csv", low_memory=False)
    d = d[~d["Loc"].isin(["UF", "DW", "BU"])]
    d = d[d["Sci_key"] != "S"]
    d = d.sort_values("Priority", na_position="last").drop_duplicates("Name", keep="first")
    out = {"name": d["Name"].to_numpy(), "feh": pd.to_numeric(d["Fe/H"], errors="coerce").to_numpy(),
           "sci_key": d["Sci_key"].to_numpy()}
    for el in ("Mg", "Ba", "Eu"):
        raw = d[el].astype(str).str.strip()
        upper = raw.str.startswith("<").to_numpy()
        val = pd.to_numeric(raw.str.lstrip("<"), errors="coerce").to_numpy()
        out[f"{el.lower()}_fe"] = val - SOLAR_JINA[el] - out["feh"]
        out[f"{el.lower()}_upper"] = upper & np.isfinite(val)
    return out


def bisterzo_s_fraction(element: str) -> float:
    """Solar s-process fraction of an element (Bisterzo et al. 2014, 'TW+' column)."""
    t = pd.read_csv(LIT_DIR / "bisterzo2014_sprocess.csv")
    row = t[t["Isotope"].astype(str).str.strip() == element]
    return float(row["TW+"].iloc[0]) / 100.0


def pure_r_ba_eu() -> float:
    """[Ba/Eu] of pure solar r-process matter."""
    return float(np.log10((1 - bisterzo_s_fraction("Ba")) / (1 - bisterzo_s_fraction("Eu"))))


def pure_s_ba_eu() -> float:
    return float(np.log10(bisterzo_s_fraction("Ba") / bisterzo_s_fraction("Eu")))
