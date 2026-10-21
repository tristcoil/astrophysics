"""The stellar graveyard: neutron stars and black holes weighed by light and by
gravitational waves.

* Neutron stars from radio pulsar timing: Paulo Freire's table of precise
  masses (the first of its three tables, "systems with individual NS mass
  measurements"). Values in red there are neutron stars; blue entries may be
  neutron stars or massive white dwarfs and are left out.
* Neutron stars weighed in X-ray binaries (Ozel & Freire 2016).
* Gravitational-wave events from GWOSC, with approximate selection weights.
"""

from __future__ import annotations

import html
import json
import re

import numpy as np
import pandas as pd

from config import GW_PASTRO_MIN, GW_SNR_THRESHOLD
from fetch_data import FREIRE_PATH, GWOSC_PATH

# ---------------------------------------------------------------------------
# Neutron stars from pulsar timing
# ---------------------------------------------------------------------------

_MASS = re.compile(r"^\s*(?P<v>\d+\.\d+)\+?\((?P<a>[^)]*)\)")


def parse_freire_value(text: str) -> tuple[float, float, float] | None:
    """'1.4380(+12/-11)' -> (1.4380, 0.0012, 0.0011); '2.01(4)' -> (2.01, 0.04, 0.04).
    The digits in brackets are the uncertainty in the last digits."""
    text = html.unescape(text).replace("−", "-").replace("–", "-").replace(" ", "")
    m = _MASS.match(text)
    if not m:
        return None
    v, a = m.group("v"), m.group("a")
    scale = 10.0 ** -len(v.split(".")[1])
    if "/" in a:
        up, down = a.split("/")
        return float(v), abs(float(up)) * scale, abs(float(down)) * scale
    return float(v), float(a) * scale, float(a) * scale


def load_freire_ns() -> pd.DataFrame:
    """One row per neutron star with a precise mass. Where a cell lists more
    than one measurement of the same star, the last (most recent) is used."""
    s = FREIRE_PATH.read_text(encoding="utf-8", errors="replace")
    start = s.find('name="part1"')
    end = s.find('name="part2"', start)  # parts 2 and 3 repeat systems
    if start < 0 or end < 0:
        raise RuntimeError("Freire's page layout changed: part1/part2 anchors not found")
    block = s[start:end]
    rows = []
    for tr in re.findall(r"<tr[^>]*>(.*?)</tr>", block, re.S):
        cells = re.findall(r"<td[^>]*>(.*?)</td>", tr, re.S)
        if len(cells) < 8:
            continue
        name = html.unescape(re.sub(r"<[^>]+>", "", cells[0])).strip()
        if not re.match(r"^[JB]\d{4}", name):
            continue
        for role, cell in (("pulsar", cells[6]), ("companion", cells[7])):
            reds = re.findall(r'<font color="?#?ff0000"?>(.*?)</font>', cell, re.I | re.S)
            vals = [parse_freire_value(r) for r in reds]
            vals = [v for v in vals if v is not None]
            if vals:
                v, up, down = vals[-1]
                rows.append({"name": name + ("" if role == "pulsar" else " companion"),
                             "role": role, "mass": v, "plus": up, "minus": down,
                             "n_measurements": len(vals)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Gravitational-wave events
# ---------------------------------------------------------------------------

def load_gw_events() -> pd.DataFrame:
    """GWOSC cumulative list -> events with source-frame masses.

    Kept: p_astro >= 0.5, or no p_astro listed (the GWTC-1 'confident' events
    predate it). Each component is labelled NS (< 2.5 M_sun), gap (2.5-5) or
    BH (> 5) from its median mass: labels for plotting, not classifications."""
    ev = json.loads(GWOSC_PATH.read_text())["events"]
    rows = []
    for key, e in ev.items():
        if e.get("mass_1_source") is None or e.get("mass_2_source") is None:
            continue
        pa = e.get("p_astro")
        if pa is not None and pa < GW_PASTRO_MIN:
            continue
        rows.append({
            "name": e["commonName"], "catalog": e["catalog.shortName"],
            "m1": e["mass_1_source"], "m1_lo": e["mass_1_source_lower"], "m1_hi": e["mass_1_source_upper"],
            "m2": e["mass_2_source"], "m2_lo": e["mass_2_source_lower"], "m2_hi": e["mass_2_source_upper"],
            "mf": e.get("final_mass_source"), "mchirp": e.get("chirp_mass_source"),
            "mchirp_det": e.get("chirp_mass"), "dl_mpc": e.get("luminosity_distance"),
            "z": e.get("redshift"), "snr": e.get("network_matched_filter_snr"),
            "chi_eff": e.get("chi_eff"), "p_astro": pa,
        })
    df = pd.DataFrame(rows).sort_values("name").reset_index(drop=True)

    def label(m):
        return np.where(m < 2.5, "NS", np.where(m < 5.0, "gap", "BH"))

    df["kind1"], df["kind2"] = label(df["m1"].values), label(df["m2"].values)
    return df


def gw_selection_weights(df: pd.DataFrame) -> pd.DataFrame:
    """Two crude corrections for 'heavier mergers are louder'.

    * chirp: an inspiral is heard to a distance proportional to Mc^(5/6), so
      the surveyed volume grows as Mc^(5/2); weight = Mc^(-5/2) (source frame).
    * vmax: each event could have been seen out to d_max = d_L x SNR / 8
      (SNR scales as 1/d_L); weight = 1 / V_max, the comoving volume to
      d_max (Planck 2018 cosmology). This uses each event's own SNR, so it
      folds in spin, mass ratio and the merger, but it ignores that the
      redshifted signal changes shape with distance.
    Both treat every observing run alike. Weights are normalised to sum to 1.
    """
    from astropy.cosmology import Planck18
    import astropy.units as u

    out = df.copy()
    w = out["mchirp"].values ** -2.5
    out["w_chirp"] = w / np.nansum(w)
    dmax = out["dl_mpc"].values * out["snr"].values / GW_SNR_THRESHOLD
    zgrid = np.geomspace(1e-4, 20, 4000)
    dl_grid = Planck18.luminosity_distance(zgrid).to_value(u.Mpc)
    zmax = np.interp(dmax, dl_grid, zgrid)
    vmax = Planck18.comoving_volume(zmax).to_value(u.Gpc**3)
    w = 1.0 / vmax
    out["vmax_gpc3"] = vmax
    out["w_vmax"] = w / np.nansum(w)
    return out


def weighted_quantile(x, w, q):
    x, w = np.asarray(x, float), np.asarray(w, float)
    ok = np.isfinite(x) & np.isfinite(w)
    x, w = x[ok], w[ok]
    order = np.argsort(x)
    cw = np.cumsum(w[order])
    cw /= cw[-1]
    return np.interp(q, cw, x[order])
