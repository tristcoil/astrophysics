"""Download and cache every input this project uses (about 50 MB in all).

    python fetch_data.py

1. **OGLE-IV classical Cepheids of the Magellanic Clouds** (Soszynski et al.
   2015, extended to 2023): periods, modes, mean V and I magnitudes, and the
   full V and I light curves of ~9,600 stars.
2. **OGLE's list of all known Milky Way classical Cepheids** (Pietrukowicz et
   al. 2021, updated March 2026), with Gaia DR3 source_ids.
3. **Gaia DR3** (ESA archive): the whole ``vari_cepheid`` table joined to
   ``gaia_source``, and the astrometry and photometry of every Milky Way
   Cepheid in OGLE's list.
4. **CDS xMatch** of the OGLE Magellanic Cepheids against Gaia DR3 (positions),
   which also counts each star's Gaia neighbours (a crowding diagnostic).
5. **Literature tables** used for checks (see ``literature.py``).

Files land in ``data/raw/`` and are never re-downloaded once present, so after
one run the analysis works offline. ``data/manifest.json`` records each file's
source, size, SHA-256 and download date.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import io
import json
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
import uuid
from pathlib import Path

import numpy as np

import config as C

UA_TOOL = "astro-atlas-fetch/1.0 (python-urllib; hobby research project)"


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _fetch(url: str, *, data: bytes | None = None, headers: dict | None = None,
           timeout: float = 180.0, retries: int = 4) -> bytes:
    """GET (or POST, with ``data``) ``url`` and return the body, retrying transient failures."""
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": UA_TOOL, **(headers or {})})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError, ConnectionError) as exc:
            last = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code in (400, 403, 404, 410):
                body = exc.read()[:500]
                raise RuntimeError(f"HTTP {exc.code} for {url}: {body!r}") from exc
            time.sleep(3.0 * (attempt + 1))
    raise RuntimeError(f"could not fetch {url}: {last}")


def _save(dest: Path, payload: bytes) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    tmp.write_bytes(payload)
    tmp.replace(dest)
    return dest


def _load_manifest() -> dict:
    if C.MANIFEST_PATH.exists():
        return json.loads(C.MANIFEST_PATH.read_text())
    return {"files": {}}


def _record(manifest: dict, dest: Path, url: str, note: str) -> None:
    blob = dest.read_bytes()
    sha = hashlib.sha256(blob).hexdigest()
    key = str(dest.relative_to(C.DATA_DIR))
    old = manifest["files"].get(key, {})
    manifest["files"][key] = {
        "url": url,
        "note": note,
        "bytes": len(blob),
        "sha256": sha,
        "fetched": old["fetched"] if old.get("sha256") == sha
        else dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }


def download(url: str, dest: Path, manifest: dict, note: str, **kw) -> Path:
    """Fetch ``url`` to ``dest`` unless it is already cached."""
    if not dest.exists():
        print(f"  downloading {url}")
        _save(dest, _fetch(url, **kw))
    _record(manifest, dest, url, note)
    return dest


# ---------------------------------------------------------------------------
# 1-2. OGLE
# ---------------------------------------------------------------------------

def fetch_ogle(manifest: dict) -> None:
    print("OGLE-IV classical Cepheids (LMC, SMC) ...")
    for gal in C.GALAXIES:
        gdir = C.OGLE_DIR / gal
        listing = _fetch(f"{C.OGLE_BASE}{gal}/cep/").decode("latin-1")
        for name in C.OGLE_CEP_FILES:
            if f'href="{name}"' not in listing:
                continue  # not every galaxy has every multi-mode file
            download(f"{C.OGLE_BASE}{gal}/cep/{name}", gdir / name, manifest,
                     f"OGLE OCVS {gal.upper()} classical Cepheids: {name}")
        phot = gdir / "phot"
        if not (phot / "I").is_dir():
            print(f"  unpacking {gal} light curves")
            with tarfile.open(gdir / "phot.tar.gz") as tar:
                tar.extractall(gdir, filter="data")
    download(C.OGLE_BASE + C.OGLE_GALACTIC_LIST, C.OGLE_DIR / C.OGLE_GALACTIC_LIST, manifest,
             "OGLE list of all known Galactic classical Cepheids (Pietrukowicz et al. 2021, updated 2026)")


# ---------------------------------------------------------------------------
# 3. Gaia DR3
# ---------------------------------------------------------------------------

#: gaia_source columns needed for the parallax zero point (Lindegren et al.
#: 2021 needs G, nu_eff or pseudocolour, ecliptic latitude and solution type),
#: quality cuts and photometry.
GAIA_SOURCE_COLS = (
    "source_id", "ra", "dec", "l", "b", "ecl_lat", "parallax", "parallax_error",
    "pmra", "pmra_error", "pmdec", "pmdec_error", "ruwe", "astrometric_params_solved",
    "nu_eff_used_in_astrometry", "pseudocolour", "astrometric_excess_noise",
    "astrometric_excess_noise_sig", "visibility_periods_used", "ipd_frac_multi_peak",
    "ipd_gof_harmonic_amplitude", "duplicated_source", "phot_g_mean_mag", "phot_bp_mean_mag",
    "phot_rp_mean_mag", "phot_g_mean_flux_over_error", "phot_bp_mean_flux_over_error",
    "phot_rp_mean_flux_over_error", "phot_bp_rp_excess_factor", "phot_variable_flag",
    "radial_velocity", "radial_velocity_error",
)
VARI_CEP_COLS = (
    "pf", "pf_error", "p1_o", "p1_o_error", "p2_o", "p2_o_error", "epoch_g",
    "int_average_g", "int_average_g_error", "int_average_bp", "int_average_bp_error",
    "int_average_rp", "int_average_rp_error", "peak_to_peak_g", "peak_to_peak_bp",
    "peak_to_peak_rp", "metallicity", "metallicity_error", "r21_g", "phi21_g", "r31_g", "phi31_g",
    "num_clean_epochs_g", "num_clean_epochs_bp", "num_clean_epochs_rp",
    "type_best_classification", "type2_best_sub_classification",
    "mode_best_classification", "multi_mode_best_classification",
)


def gaia_query(adql: str, fmt: str = "votable") -> bytes:
    """One synchronous ADQL query against the ESA Gaia archive."""
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": fmt,
                                   "QUERY": adql}).encode()
    body = _fetch(C.GAIA_TAP_SYNC, data=data, timeout=300)
    if b'value="OVERFLOW"' in body:
        raise RuntimeError("Gaia query result truncated (OVERFLOW)")
    if b'value="ERROR"' in body:
        raise RuntimeError("Gaia query error: " + body[:2000].decode("latin-1"))
    return body


def _votable_to_table(blob: bytes):
    from astropy.io.votable import parse_single_table
    tab = parse_single_table(io.BytesIO(blob)).to_table(use_names_over_ids=True)
    for name in list(tab.colnames):  # the archive's column-name case varies
        if name != name.lower():
            tab.rename_column(name, name.lower())
    return tab


def _write_fits(tab, dest: Path) -> None:
    from astropy.table import Table
    tab = Table(tab, masked=True)
    for name in tab.colnames:  # object-dtype strings cannot go to FITS
        if tab[name].dtype.kind == "O":
            tab[name] = np.asarray([str(x) for x in np.ma.filled(tab[name], "")])
    dest.parent.mkdir(parents=True, exist_ok=True)
    tab.write(dest, overwrite=True)


def fetch_gaia_vari_cepheid(manifest: dict) -> None:
    """Every Gaia DR3 Cepheid (classical, type II, anomalous), joined to gaia_source."""
    dest = C.GAIA_DIR / "vari_cepheid_x_gaia_source.fits"
    if not dest.exists():
        print("Gaia DR3 vari_cepheid x gaia_source ...")
        cols = ", ".join([f"g.{c}" for c in GAIA_SOURCE_COLS] + [f"v.{c}" for c in VARI_CEP_COLS])
        adql = (f"SELECT {cols} FROM gaiadr3.vari_cepheid AS v "
                f"JOIN gaiadr3.gaia_source AS g ON g.source_id = v.source_id")
        tab = _votable_to_table(gaia_query(adql))
        n_server = int(gaia_query("SELECT COUNT(*) AS n FROM gaiadr3.vari_cepheid", fmt="csv").split()[-1])
        if len(tab) != n_server:
            raise RuntimeError(f"vari_cepheid: got {len(tab)} rows, server has {n_server}")
        _write_fits(tab, dest)
    _record(manifest, dest, C.GAIA_TAP_SYNC, "Gaia DR3 vari_cepheid JOIN gaia_source (all rows)")


def galactic_list_ids() -> np.ndarray:
    """Gaia DR3 source_ids in OGLE's Galactic Cepheid list (0 where none)."""
    ids = []
    for line in (C.OGLE_DIR / C.OGLE_GALACTIC_LIST).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        tok = line.split()
        ids.append(int(tok[-2]) if tok[-2].isdigit() else 0)
    return np.array(ids, dtype=np.int64)


def fetch_gaia_galactic(manifest: dict, chunk: int = 400) -> None:
    """gaia_source (+ vari_cepheid where present) for OGLE's Milky Way Cepheids."""
    dest = C.GAIA_DIR / "galactic_cepheids_gaia.fits"
    if not dest.exists():
        from astropy.table import vstack
        ids = np.unique(galactic_list_ids())
        ids = ids[ids > 0]
        print(f"Gaia DR3 astrometry for {len(ids)} Milky Way Cepheids ...")
        cols = ", ".join([f"g.{c}" for c in GAIA_SOURCE_COLS] + [f"v.{c}" for c in VARI_CEP_COLS])
        parts = []
        for k in range(0, len(ids), chunk):
            id_list = ",".join(str(i) for i in ids[k:k + chunk])
            adql = (f"SELECT {cols} FROM gaiadr3.gaia_source AS g "
                    f"LEFT OUTER JOIN gaiadr3.vari_cepheid AS v ON v.source_id = g.source_id "
                    f"WHERE g.source_id IN ({id_list})")
            parts.append(_votable_to_table(gaia_query(adql)))
        tab = vstack(parts)
        missing = set(ids) - set(np.asarray(tab["source_id"]))
        if missing:
            print(f"  note: {len(missing)} listed source_ids not found in gaia_source")
        _write_fits(tab, dest)
    _record(manifest, dest, C.GAIA_TAP_SYNC,
            "Gaia DR3 gaia_source LEFT JOIN vari_cepheid for OGLE's Galactic Cepheid list")


# ---------------------------------------------------------------------------
# 4. CDS xMatch: OGLE Magellanic Cepheids -> Gaia DR3
# ---------------------------------------------------------------------------

def _xmatch(csv_bytes: bytes) -> bytes:
    boundary = uuid.uuid4().hex
    fields = {"request": "xmatch", "distMaxArcsec": str(C.XMATCH_RADIUS_ARCSEC),
              "RESPONSEFORMAT": "csv", "colRA1": "ra", "colDec1": "dec", "cat2": C.VIZIER_GAIA_DR3}
    out = io.BytesIO()
    for k, v in fields.items():
        out.write(f"--{boundary}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n{v}\r\n".encode())
    out.write(f"--{boundary}\r\nContent-Disposition: form-data; name=\"cat1\"; filename=\"ogle.csv\"\r\n"
              f"Content-Type: text/csv\r\n\r\n".encode())
    out.write(csv_bytes)
    out.write(f"\r\n--{boundary}--\r\n".encode())
    return _fetch(C.CDS_XMATCH, data=out.getvalue(), timeout=600,
                  headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})


def fetch_xmatch(manifest: dict) -> None:
    import ogle
    dest = C.GAIA_DIR / "xmatch_ogle_mc_gaiadr3.csv"
    if not dest.exists():
        print("CDS xMatch: OGLE Magellanic Cepheids vs Gaia DR3 ...")
        lines = ["id,ra,dec"]
        for gal in C.GALAXIES:
            ident = ogle.read_ident(gal)
            lines += [f"{i},{r:.6f},{d:.6f}" for i, r, d in zip(ident["id"], ident["ra"], ident["dec"])]
        _save(dest, _xmatch(("\n".join(lines) + "\n").encode()))
    _record(manifest, dest, C.CDS_XMATCH,
            f"CDS xMatch of OGLE LMC+SMC Cepheids with Gaia DR3 (I/355), radius {C.XMATCH_RADIUS_ARCSEC}\"")


# ---------------------------------------------------------------------------

def fetch_all() -> dict:
    manifest = _load_manifest()
    try:
        fetch_ogle(manifest)
        fetch_gaia_vari_cepheid(manifest)
        fetch_gaia_galactic(manifest)
        fetch_xmatch(manifest)
        import literature_data
        literature_data.fetch(manifest, download)
    finally:
        C.MANIFEST_PATH.write_text(json.dumps(manifest, indent=1, sort_keys=True))
    return manifest


if __name__ == "__main__":
    m = fetch_all()
    total = sum(f["bytes"] for f in m["files"].values())
    print(f"{len(m['files'])} files cached, {total / 1e6:.1f} MB")
    sys.exit(0)
