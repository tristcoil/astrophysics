"""The Galaxy's supernova remnants: Green's catalogue, and how blast waves grow.

A remnant's blast wave goes through stages:

1. **Free expansion.** For the first few hundred years the ejecta outweigh the
   gas they have swept up and coast: R = v t, with v ~ 5,000-10,000 km/s.
2. **Sedov-Taylor.** Once the swept-up mass exceeds the ejecta's, the blast is
   a point explosion of energy E in gas of density rho, and dimensional
   analysis alone fixes its size: R = xi (E t^2 / rho)^(1/5), xi = 1.15 for a
   monatomic gas (Sedov 1959; Taylor 1950). It decelerates: R ~ t^0.4.
3. Later, radiative cooling (not reached by any remnant here).

The expansion parameter m = (dlnR/dlnt) = v t / R marks the stage: 1 for free
expansion, 0.4 for Sedov-Taylor.
"""

from __future__ import annotations

import html
import re

import numpy as np
import pandas as pd

import config
import literature as L

PC_CM = 3.0857e18
M_H = 1.6726e-24
YEAR_S = 3.15576e7
KM_CM = 1e5


def sedov_radius_pc(t_yr, E_erg: float = 1e51, n_cm3: float = 1.0, mu: float = 1.4, xi: float = 1.15):
    """Sedov-Taylor blast-wave radius (pc): R = xi (E t^2 / rho)^(1/5),
    rho = mu m_H n (mu = 1.4 for 10% helium by number)."""
    rho = mu * M_H * n_cm3
    return xi * (E_erg * (np.asarray(t_yr) * YEAR_S) ** 2 / rho) ** 0.2 / PC_CM


def free_radius_pc(t_yr, v_kms: float):
    """Radius (pc) of debris coasting at v_kms for t_yr."""
    return v_kms * KM_CM * np.asarray(t_yr) * YEAR_S / PC_CM


def load_green() -> pd.DataFrame:
    """Parse the summary table of Green's catalogue (2024 October version).

    Each row: l, b, RA (h m s), Dec (d m), size (arcmin, 'a' or 'axb'),
    type (S shell, C composite, F filled-centre; '?' = uncertain), 1-GHz flux,
    spectral index, other names. Sizes 'axb' become the geometric mean."""
    text = (config.CATALOG_DIR / "green_snrs_2024oct.html").read_text(errors="replace")
    pre = text[text.index("<PRE>") : text.index("</PRE>")]
    rows = []
    for line in pre.splitlines():
        if "HREF" not in line:
            continue
        m = re.search(r'HREF="snrs\.(G[^"]+)\.html">\s*([\d.]+)\s+([+-][\d.]+)</A>(.*)', line)
        if not m:
            continue
        name, l, b, rest = m.group(1), float(m.group(2)), float(m.group(3)), html.unescape(m.group(4))
        parts = rest.split()
        rah, ram, ras, dd, dm = parts[:5]
        size_txt, type_txt = parts[5], parts[6]
        flux_txt = parts[7] if len(parts) > 7 else ""
        names = " ".join(parts[9:]) if len(parts) > 9 else ""
        ra = 15.0 * (float(rah) + float(ram) / 60.0 + float(ras) / 3600.0)
        sign = -1.0 if dd.startswith("-") else 1.0
        dec = sign * (abs(float(dd)) + float(dm) / 60.0)
        dims = [float(x) for x in re.findall(r"[\d.]+", size_txt)]
        size = float(np.sqrt(np.prod(dims))) if len(dims) == 2 else (dims[0] if dims else np.nan)
        fl = re.findall(r"[\d.]+", flux_txt)
        rows.append({"name": name, "l": l, "b": b, "ra": ra, "dec": dec, "size_arcmin": size,
                     "size_text": size_txt, "type": type_txt.rstrip("?"), "type_uncertain": "?" in type_txt,
                     "flux_1ghz": float(fl[0]) if fl else np.nan, "names": names})
    t = pd.DataFrame(rows)
    if len(t) != 310:
        raise RuntimeError(f"expected 310 remnants in Green's 2024 catalogue, parsed {len(t)}")
    return t


def historical(green: pd.DataFrame, now: float = 2025.0) -> pd.DataFrame:
    """The remnants of supernovae seen (or dated) in history, with physical sizes."""
    rows = []
    for event, remnant, gname, year, d, d_rng, ref, note in L.HISTORICAL:
        g = green[green.name == gname]
        if g.empty:
            raise KeyError(gname)
        size = float(g.iloc[0].size_arcmin)
        r_pc = 0.5 * size / 60.0 * np.pi / 180.0 * d * 1000.0
        rows.append({"event": event, "remnant": remnant, "green": gname, "year": year, "age": now - year,
                     "size_arcmin": size, "size_text": g.iloc[0].size_text, "type": g.iloc[0].type,
                     "d_kpc": d, "d_lo": d_rng[0], "d_hi": d_rng[1], "radius_pc": r_pc,
                     "radius_lo": r_pc * d_rng[0] / d, "radius_hi": r_pc * d_rng[1] / d,
                     "mean_speed_kms": r_pc * PC_CM / ((now - year) * YEAR_S) / KM_CM,
                     "disputed": "?" in event, "reference": ref, "note": note})
    return pd.DataFrame(rows)


def analyse() -> dict:
    green = load_green()
    hist = historical(green)
    counts = green["type"].value_counts().to_dict()
    return {"n_remnants": int(len(green)), "type_counts": {str(k): int(v) for k, v in counts.items()},
            "size_median_arcmin": float(np.nanmedian(green.size_arcmin)),
            "historical": hist.to_dict("records")}
