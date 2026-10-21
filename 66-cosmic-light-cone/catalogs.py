"""Load the real objects we place on the light cone, all into one simple form.

Each population is a redshift array plus labels and a citation. Loaders read
only from ``data/`` (run ``fetch_data.py`` first) and never touch the network.

Redshift is what places an object on our past light cone: a source we now see
at redshift z sent out its light when the universe was 1/(1+z) of its present
size, and that one number fixes when (cosmic time), from how far (proper
distance then) and how fast it was receding then and now, once the expansion
history is known.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field

import numpy as np

import config as C
from fetch_data import parse_vizier_tsv


@dataclass
class Population:
    key: str
    label: str
    z: np.ndarray
    source: str
    names: list = field(default_factory=list)
    z_lo: np.ndarray | None = None   # 90% bounds where the catalogue gives them (GW events)
    z_hi: np.ndarray | None = None
    ra: np.ndarray | None = None
    dec: np.ndarray | None = None
    note: str = ""

    @property
    def n(self) -> int:
        return int(self.z.size)


def _require(path):
    if not path.exists():
        raise FileNotFoundError(f"{path} is not cached; run `python fetch_data.py` first")
    return path


# --- SDSS ------------------------------------------------------------------------------


def sdss_galaxies() -> Population:
    d = np.load(_require(C.DERIVED_DIR / "sdss_mgs_z.npz"))
    return Population("sdss_galaxies", "SDSS galaxies (D2)", d["z"],
                      "SDSS DR17 Main Galaxy Sample (Strauss et al. 2002), D2's selection",
                      ra=d["ra"], dec=d["dec"])


def dr16q_quasars() -> Population:
    d = np.load(_require(C.DERIVED_DIR / "dr16q_z.npz"))
    z = d["z"]
    ok = z > 0  # DR16Q keeps a few objects with z <= 0 (stars and failures flagged by -1)
    return Population("quasars", "SDSS quasars (I3)", z[ok],
                      "SDSS DR16Q (Lyke et al. 2020, ApJS 250, 8), VizieR VII/289",
                      ra=d["ra"][ok], dec=d["dec"][ok],
                      note=f"{int((~ok).sum())} entries with z <= 0 dropped")


def highz_quasars() -> Population:
    rows = parse_vizier_tsv(_require(C.RAW_DIR / "highz_quasars" / "banados2016_highzqso.tsv")
                            .read_text(encoding="utf-8"))
    rows = [r for r in rows if r["z"]]
    return Population("highz_quasars", "Quasars at z > 5.6", np.array([float(r["z"]) for r in rows]),
                      "Banados et al. 2016, ApJS 227, 11 (VizieR J/ApJS/227/11)",
                      names=[r["Name"] for r in rows],
                      ra=np.array([float(r["RAJ2000"]) for r in rows]),
                      dec=np.array([float(r["DEJ2000"]) for r in rows]))


# --- Supernovae -----------------------------------------------------------------------------


def pantheon_table():
    """The Pantheon+ table as a dict of columns (all 1,701 light curves)."""
    path = _require(C.RAW_DIR / "sn" / "Pantheon+SH0ES.dat")
    lines = path.read_text().split("\n")
    header = lines[0].split()
    cols = {h: [] for h in header}
    for ln in lines[1:]:
        if ln.strip():
            for h, v in zip(header, ln.split()):
                cols[h].append(v)
    out = {}
    for h, v in cols.items():
        try:
            out[h] = np.array(v, dtype=float)
        except ValueError:
            out[h] = np.array(v)
    return out


def pantheon_cov() -> np.ndarray:
    """The 1701 x 1701 STAT+SYS covariance (first line holds the size)."""
    path = _require(C.RAW_DIR / "sn" / "Pantheon+SH0ES_STAT+SYS.cov")
    cache = C.DERIVED_DIR / "pantheon_cov.npy"
    if cache.exists() and cache.stat().st_mtime > path.stat().st_mtime:
        return np.load(cache)
    flat = np.loadtxt(path)
    n = int(flat[0])
    cov = flat[1:].reshape(n, n)
    C.DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    np.save(cache, cov)
    return cov


def supernovae() -> Population:
    """One entry per supernova (several surveys observed some: average their zHD)."""
    t = pantheon_table()
    cid, z = t["CID"], t["zHD"]
    names, inv = np.unique(cid, return_inverse=True)
    zz = np.bincount(inv, weights=z) / np.bincount(inv)
    ra = np.bincount(inv, weights=t["RA"]) / np.bincount(inv)
    dec = np.bincount(inv, weights=t["DEC"]) / np.bincount(inv)
    return Population("supernovae", "Type Ia supernovae (D1)", zz,
                      "Pantheon+ (Scolnic et al. 2022, ApJ 938, 113; Brout et al. 2022)",
                      names=list(names), ra=ra, dec=dec,
                      note=f"{t['CID'].size} light curves of {names.size} supernovae")


# --- Gravitational waves ---------------------------------------------------------------------


def gw_events(universe=None) -> Population:
    """Confident GWOSC mergers with a distance. The catalogue redshifts assume the LVK
    collaboration's cosmology; with ``universe`` given, z is recomputed from the
    measured luminosity distance (median and 90% bounds) in that universe."""
    ev = json.loads(_require(C.RAW_DIR / "gwosc" / "gwosc_GWTC.json").read_text())["events"]
    rows = []
    for v in ev.values():
        if v.get("luminosity_distance") is None or v.get("redshift") is None:
            continue
        p = v.get("p_astro")
        if p is not None and p < 0.5:
            continue
        dl = v["luminosity_distance"]
        rows.append((v["commonName"], v["catalog.shortName"], v["redshift"],
                     v["redshift"] + (v.get("redshift_lower") or 0.0),
                     v["redshift"] + (v.get("redshift_upper") or 0.0),
                     dl, dl + (v.get("luminosity_distance_lower") or 0.0),
                     dl + (v.get("luminosity_distance_upper") or 0.0)))
    rows.sort(key=lambda r: r[2])
    names = [r[0] for r in rows]
    z_cat = np.array([r[2] for r in rows])
    if universe is None:
        z, lo, hi = z_cat, np.array([r[3] for r in rows]), np.array([r[4] for r in rows])
        note = "catalogue redshifts (LVK cosmology)"
    else:
        dl = np.array([[r[5], r[6], r[7]] for r in rows]) * C.GLY_PER_MPC
        zgrid = np.geomspace(1e-4, 20.0, 20000)
        dgrid = universe.dl_of_z(zgrid)
        z, lo, hi = (np.interp(dl[:, k], dgrid, zgrid) for k in range(3))
        note = "redshifts recomputed from the measured luminosity distances in our cosmology"
    pop = Population("gw_events", "Gravitational-wave mergers", z,
                     "GWOSC cumulative catalogue GWTC-1 to GWTC-5.0 (LIGO/Virgo/KAGRA)",
                     names=names, z_lo=lo, z_hi=hi, note=note)
    pop.z_catalogue = z_cat
    return pop


# --- JWST: JADES DR4 -------------------------------------------------------------------------


def jades_galaxies(max_sep_arcsec: float = 0.3) -> Population:
    """JADES DR4 robust redshifts (flags A, B, C), repeat observations merged by position."""
    from astropy.io import fits

    with fits.open(_require(C.RAW_DIR / "jades" / "jades_dr4_obs_info.fits"), memmap=False) as h:
        d = h[1].data
        flag = np.char.strip(np.asarray(d["z_Spec_flag"]).astype(str))
        z = np.asarray(d["z_Spec"], float)
        ra = np.asarray(d["RA_TARG"], float)
        dec = np.asarray(d["Dec_TARG"], float)
        pname = np.char.strip(np.asarray(d["z_paper_name"]).astype(str))
        pref = np.char.strip(np.asarray(d["z_paper_ref"]).astype(str))
    robust = np.isin(flag, ["A", "B", "C"]) & (z > 0)
    n_rows = int(robust.sum())
    z, ra, dec, pname, pref, flag = (v[robust] for v in (z, ra, dec, pname, pref, flag))
    # Merge repeats: same target within max_sep (sorted by dec, small neighbourhoods).
    keep = np.ones(z.size, bool)
    order = np.argsort(dec)
    sep_deg = max_sep_arcsec / 3600.0
    for ii, i in enumerate(order):
        if not keep[i]:
            continue
        for j in order[ii + 1:]:
            if dec[j] - dec[i] > sep_deg:
                break
            if keep[j] and abs((ra[j] - ra[i]) * np.cos(np.radians(dec[i]))) < sep_deg:
                if not pname[i] and pname[j]:
                    pname[i], pref[i] = pname[j], pref[j]
                keep[j] = False
    pop = Population("jades", "JWST galaxies (JADES)", z[keep],
                     "JADES DR4 (Curtis-Lake et al. 2025; Scholtz, Carniani et al. 2025)",
                     names=[n for n in pname[keep]], ra=ra[keep], dec=dec[keep],
                     note=f"{n_rows} robust redshifts (flags A/B/C) of {int(keep.sum())} "
                          f"distinct galaxies (repeats within {max_sep_arcsec}\" merged)")
    pop.refs = list(pref[keep])
    pop.flags = list(flag[keep])
    pop.n_rows = n_rows
    return pop


# --- DESI ---------------------------------------------------------------------------------

DESI_TRACERS = {0.295: "BGS", 0.51: "LRG", 0.706: "LRG", 0.934: "LRG+ELG", 1.321: "ELG",
                1.484: "QSO", 2.33: "Ly-alpha"}


def desi_bao() -> dict:
    """DESI DR2 BAO: redshifts, values, observable names and covariance (13 numbers)."""
    rows = []
    for ln in _require(C.RAW_DIR / "bao" / "desi_dr2_mean.txt").read_text().splitlines():
        if ln.strip() and not ln.startswith("#"):
            a, b, c = ln.split()
            rows.append((float(a), float(b), c))
    cov = np.loadtxt(_require(C.RAW_DIR / "bao" / "desi_dr2_cov.txt"))
    z = np.array([r[0] for r in rows])
    if cov.shape != (z.size, z.size):
        raise ValueError("DESI covariance does not match the data vector")
    return {"z": z, "value": np.array([r[1] for r in rows]), "kind": [r[2] for r in rows],
            "cov": cov, "tracer": [DESI_TRACERS[round(v, 3)] for v in z]}


# --- Landmarks --------------------------------------------------------------------------------


def landmarks(jades: Population | None = None) -> list[dict]:
    """Literature record-holders from config, plus JADES's named galaxies at z > 10."""
    out = [{"key": L.key, "name": L.name, "kind": L.kind, "z": L.z, "z_err": L.z_err,
            "ref": L.ref, "note": L.note} for L in C.LANDMARKS]
    if jades is not None:
        seen = set()
        for n, zz, ref in sorted(zip(jades.names, jades.z, jades.refs), key=lambda r: -r[1]):
            if n and zz > 10 and n not in seen:
                seen.add(n)
                out.append({"key": "jades_" + re.sub(r"[^A-Za-z0-9]+", "_", n).lower(),
                            "name": ("JADES-" + n) if n.startswith("GS-") else n,
                            "kind": "galaxy", "z": float(zz), "z_err": np.nan, "ref": ref + " (via JADES DR4)",
                            "note": "named in JADES DR4"})
    return sorted(out, key=lambda r: -r["z"])


def all_populations(universe=None) -> dict:
    jades = jades_galaxies()
    return {
        "sdss_galaxies": sdss_galaxies(),
        "supernovae": supernovae(),
        "gw_events": gw_events(universe),
        "quasars": dr16q_quasars(),
        "highz_quasars": highz_quasars(),
        "jades": jades,
    }
