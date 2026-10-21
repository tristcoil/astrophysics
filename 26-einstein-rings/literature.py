"""Published per-lens results, used only to grade ours.

Two SLACS papers list per-lens numbers that CDS does not hold, so we read them
from the LaTeX sources the authors posted to arXiv:

* Auger et al. (2010, SLACS X), Table 1: for 59 early-type lenses, the
  velocity dispersion inside r_e/2, stellar masses, the total mass inside
  r_e/2, the dark-matter fraction there for Chabrier and Salpeter IMFs, and
  the slope gamma' of a power-law total density rho ~ r^-gamma' fitted to
  lensing plus stellar dynamics.
* Treu et al. (2010), Table 1: for 56 lenses, sigma_SIE and the "IMF mismatch
  parameter" log(alpha) = log(M_star from lensing + dynamics / M_star from a
  Salpeter-IMF stellar population model).

The analysis never imports this module; only the grading at the end does.
"""

from __future__ import annotations

import re
import tarfile
from functools import lru_cache

import pandas as pd

from catalog import short_name
from config import LIT_DIR

AUGER10_TAR = LIT_DIR / "auger2010_1007.2880.tar.gz"
TREU10_TAR = LIT_DIR / "treu2010_0911.3392.tar.gz"

_PM = r"\$?\s*([-+]?\d+\.?\d*)\s*\\pm\s*(?:\\phantom\{0\})?\s*(\d+\.?\d*)\s*\$?"


def _member_text(path, member_name: str) -> str:
    with tarfile.open(path, "r:gz") as tar:
        member = next(m for m in tar.getmembers() if m.name.endswith(member_name))
        return tar.extractfile(member).read().decode("latin-1")


def _clean_tex(s: str) -> str:
    s = s.replace("$-$", "-").replace("$+$", "+").replace("\\phantom{0}", "").replace("\\phantom{-}", "")
    return s


@lru_cache(maxsize=None)
def auger2010() -> pd.DataFrame:
    """Auger et al. (2010), Table 1."""
    text = _member_text(AUGER10_TAR, "T_lens_data.tex")
    rows = []
    for line in text.splitlines():
        if not line.strip().startswith("SDSSJ"):
            continue
        cells = [_clean_tex(c).strip() for c in line.split("&")]
        name = short_name(cells[0].replace("$\\dagger$", "").strip())
        vals = []
        for c in cells[2:]:
            m = re.search(_PM, c.replace("\\\\", ""))
            vals.append((float(m.group(1)), float(m.group(2))))
        rows.append({
            "name": name,
            "flagged": "dagger" in cells[0],
            "re_v_kpc": float(cells[1]),
            "sigma_e2": vals[0][0], "e_sigma_e2": vals[0][1],
            "log_mstar_chab": vals[1][0], "log_mstar_salp": vals[2][0],
            "log_m_re2": vals[3][0], "e_log_m_re2": vals[3][1],
            "fdm_chab": vals[4][0], "e_fdm_chab": vals[4][1],
            "fdm_salp": vals[5][0], "e_fdm_salp": vals[5][1],
            "gamma": vals[6][0], "e_gamma": vals[6][1],
        })
    return pd.DataFrame(rows)


@lru_cache(maxsize=None)
def treu2010() -> pd.DataFrame:
    """Treu et al. (2010), Table 1."""
    text = _member_text(TREU10_TAR, "ms.tex")
    rows = []
    pat = re.compile(
        r"^(J\d{4}[+-]\d{4})\s*&\s*([\d.]+)\s*&\s*(\d+)\$\\pm\$(\d+)\s*&\s*([\d.]+)\$\\pm\$([\d.]+)"
        r"\s*&\s*([\d.]+)\$\\pm\$([\d.]+)\s*&\s*([-+]?[\d.]+)\$\^\{\+([\d.]+)\}_\{-([\d.]+)\}\$")
    for line in text.splitlines():
        m = pat.match(line.strip())
        if not m:
            continue
        g = m.groups()
        rows.append({
            "name": g[0], "z_lens": float(g[1]), "sigma": float(g[2]), "e_sigma": float(g[3]),
            "sigma_sie": float(g[4]), "e_sigma_sie": float(g[5]),
            "lum_v": float(g[6]), "log_alpha": float(g[8]),
            "e_log_alpha_hi": float(g[9]), "e_log_alpha_lo": float(g[10]),
        })
    return pd.DataFrame(rows)


#: Ensemble results quoted in the papers' abstracts and text.
PUBLISHED = {
    # Bolton et al. 2008b (SLACS VII), section 4: 53 lenses
    "f_sdss_over_sie": (0.948, 0.008),
    "f_e8_over_sie": (1.019, 0.008),
    "f_rms_sdss": 0.09, "f_rms_e8": 0.10,
    "f_intrinsic_scatter": 0.075,
    # Auger et al. 2009 (SLACS IX): stellar fraction inside the Einstein radius
    "fstar_ein_chab": (0.4, 0.1), "fstar_ein_salp": (0.7, 0.2),
    # Koopmans et al. 2009: 58 lenses
    "gamma_koopmans09": (2.085, 0.025, 0.018), "gamma_scatter_koopmans09": 0.20,
    # Auger et al. 2010 (SLACS X): 73 lenses
    "gamma_auger10": (2.078, 0.027), "gamma_scatter_auger10": (0.16, 0.02),
    # Treu et al. 2010: 56 lenses
    "log_alpha_salp": (0.00, 0.03), "log_alpha_chab": (0.25, 0.03),
    # Bolton et al. 2008a: adopted empirical error on Einstein radii
    "einstein_radius_error": 0.02,
}


if __name__ == "__main__":
    a, t = auger2010(), treu2010()
    print(len(a), "lenses in Auger+2010 Table 1;", len(t), "in Treu+2010 Table 1")
    print(a.head().to_string())
    print(t.head().to_string())
    print("mean gamma (Auger+2010 table):", a["gamma"].mean().round(3))
    print("mean log alpha (Treu+2010 table):", t["log_alpha"].mean().round(3))
