"""Published numbers we compare against, read from the papers' LaTeX where
possible (so they cannot be mistyped) and otherwise quoted with their source.
"""

from __future__ import annotations

import re

import config as C

_VAL = re.compile(r"\$?\s*(-?[\d.]+)\s*(?:\^\{\+?([\d.]+)\}_\{-?([\d.]+)\}|\\pm\s*([\d.]+))\s*\$?")
_NAMES = {r"$a_": "a_", r"$\tau_": "tau_", r"$\omega_": "omega_", r"$\Omega_": "Omega_",
          r"$i_": "i_", r"$e_": "e_"}


def _value(cell: str):
    m = _VAL.search(cell)
    if not m:
        return None
    med = float(m.group(1))
    if m.group(4):
        return med, float(m.group(4)), float(m.group(4))
    return med, float(m.group(2)), float(m.group(3))


def _param(cell: str, planet: str | None):
    for tex, name in _NAMES.items():
        if tex in cell:
            return name + cell.split(tex)[1][0]
    if "Parallax" in cell:
        return "plx"
    if "M_{\\star}" in cell:
        return "mstar"
    if "chi^2" in cell:
        return "chi2nu"
    if "BIC" in cell:
        return "dBIC"
    return None


def _table(label: str, ncols: int) -> list[dict]:
    tex = (C.ARXIV["Wang+2018"].folder / "main_clean.tex").read_text()
    i = tex.index(r"\label{" + label + "}")
    body = tex[tex.index(r"\startdata", i):tex.index(r"\enddata", i)]
    cols = [dict() for _ in range(ncols)]
    for line in body.split(r"\\"):
        cells = [c.strip() for c in line.replace(r"\hline", "").split("&")]
        if len(cells) < 2 + ncols:
            continue
        name = _param(cells[1], None)
        if name is None:
            continue
        for k in range(ncols):
            v = _value(cells[2 + k])
            if v is not None:
                cols[k][name] = v
    return cols


def wang2018_models() -> dict:
    """Table 2 of Wang et al. 2018: model -> {param: (median, +, -)}."""
    cols = _table("table:mcmcfits", 4)
    out = dict(zip(C.MODELS, cols))
    # the stable-orbit counts (first 10^6 draws) are in the last row
    tex = (C.ARXIV["Wang+2018"].folder / "main_clean.tex").read_text()
    row = re.search(r"Stable Orbits.*?&\s*(\d+)\s*&\s*(\d+)\s*&\s*(\d+)\s*&\s*(\d+)", tex, re.S)
    for model, n in zip(C.MODELS, row.groups()):
        out[model]["stable_per_million"] = int(n)
    return out


def wang2018_stable() -> dict:
    """Table 3: the 9792 stable coplanar orbits."""
    return _table("table:stablefits", 1)[0]


#: Other published values used as checks (quoted from the papers' text).
WANG2018_TEXT = {
    "n_stable_coplanar": 9792, "n_sims_coplanar": 22_000_000,
    "frac_mcde_lt6": 0.956, "frac_mcde_lt5": 0.730,
    "mcde_999": 7.6, "mb_999": 6.3,
    "Pb/Pc": (2.11, 0.06, 0.06), "Pc/Pd": (2.10, 0.04, 0.05), "Pd/Pe": (2.06, 0.03, 0.04),
    "omega_e-omega_d": (94, 11, 9),
    "mass_b_lum": (5.8, 0.5, 0.5), "mass_cde_lum": (7.2, 0.6, 0.7),
    # Sec. 4.3-4.5: only 1.2% of stable orbits have theta_c:d:e librating >= 90% of
    # the time; theta_d:e,e locked for most stable orbits above ~6 MJup, theta_c:d,d
    # above ~7; no angle involving b ever reaches 100%. With the hot-start prior,
    # 95% of stable systems have M_b > 4.9 and M_cde > 6.1. Forwards: stable for
    # 1 Gyr only up to ~6 MJup; above 7 MJup none lasts 100 Myr.
    "theta_cde_lib90_frac": 0.012, "dee_locked_above_mcde": 6.0, "cdd_locked_above_mcde": 7.0,
    "hotstart_95_lower_b": 4.9, "hotstart_95_lower_cde": 6.1,
    "future_none_100Myr_above_mcde": 7.0, "future_1Gyr_up_to_mcde": 6.0,
}
#: GRAVITY Collaboration 2019, Sec. 3: Wang et al.'s 9792 stable coplanar
#: orbits predicted e at (-352.6 +3.1/-2.6, 157.9 +-1.8) mas on 2018 Aug 28
#: (the paper prints -157.9; the sign must be + - e was north of the star),
#: 5 mas from the measurement in both coordinates.
GRAVITY2019_PREDICTION = {"ra": (-352.6, 3.1, 2.6), "dec": (157.9, 1.8, 1.8)}
BRANDT2021_MASS_E = (9.6, 1.9, 1.8)                 # MJup, HGCA + orbits
ZURLO2022_PARALLAX = (24.49, 0.07)
SEPULVEDA2022_MSTAR = (1.47, 0.12, 0.17)            # Msun, dynamical (uniform prior)
DEROSA2020_NORTH_CHANGE = 0.36                      # deg, weighted mean revision of GPI north


if __name__ == "__main__":
    import json
    print(json.dumps(wang2018_models(), indent=1)[:1500])
    print(wang2018_stable())
