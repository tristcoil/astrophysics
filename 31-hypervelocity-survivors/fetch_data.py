"""Download and cache everything the analysis needs.

All sources are free and public:

1. **Gaia DR3** (ESA), through full mirrors of the catalogue (ARI Heidelberg,
   NOIRLab Astro Data Lab; ESA's own archive as a fallback):
   - El-Badry et al. (2023)'s search query, verbatim (25 sources published);
   - our own broad search: every source with a face-value tangential velocity
     above 600 km/s and a parallax measured at > 3 sigma, *and* its mirror
     image, sources whose parallax is significantly negative (the cut runs on
     ARI's narrow gaia_source_lite table, patch by patch over the sky; full
     records are then fetched by source_id);
   - the full astrometric record of every known runaway (literature.py);
   - Shen et al. (2018)'s seven Gaia DR2 candidates, matched to DR3 through
     the archive's DR2 -> DR3 neighbourhood table, plus their DR2 records;
   - stars within 100 pc for a background HR diagram;
   - a random all-sky subsample, to measure parallax precision vs magnitude;
   - the Rybizki et al. (2022) spurious-astrometry classifier for every
     candidate (ARI only).
2. **Green (2025)'s catalogue of 310 Galactic supernova remnants** (VizieR
   VII/297).
3. **Bedard et al. (2020) white dwarf cooling sequences**, 0.20-1.30 Msun, from
   the Montreal group's website (radius vs temperature for each mass).
4. **SFD dust reddening** (Schlafly & Finkbeiner 2011 recalibration) toward
   every star we place in the HR diagram, from IRSA's DUST service.
5. **Finkbeiner (2003)'s all-sky H-alpha map** (NASA LAMBDA), in which old
   supernova remnants show up as faint filaments.

A manifest (data/manifest.json) records each file's URL, date, size and
SHA-256; every Gaia query also gets a JSON sidecar with its exact ADQL and the
service that answered. Nothing cached is fetched again, so after one run the
project works offline.

    python fetch_data.py               # fill the cache
    python fetch_data.py --only gaia   # one group: gaia, green, bedard, dust, halpha
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

import gaia_tap
from binary import BEDARD_MASSES, bedard_path
from config import (AUX_DIR, BACKGROUND_PLX_MIN, BACKGROUND_PLX_SNR, DATA_DIR, ELBADRY_QUERY,
                    G_MAX, GAIA_DIR, GAIA_SERVICES, MANIFEST_PATH, PLX_SNR_MIN, PM_MIN,
                    PRECISION_RANDOM_INDEX_MAX, V_TAN_MIN)
from literature import CENSUS, ELBADRY_TABLE1, SHEN18_TABLE1

USER_AGENT = "astrophysics-atlas-L5/1.0 (hobbyist research; hypervelocity white dwarfs)"
TIMEOUT_S = 180
RETRIES = 5

GREEN_URL = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=VII/297/snrs"
             "&-out.max=1000&-out.all&-out.add=_Glon,_Glat")
GREEN_PATH = AUX_DIR / "green2025_snrs.tsv"
BEDARD_URL = "https://www.astro.umontreal.ca/~bergeron/CoolingModels/CoolingModels/{name}"
DUST_URL = "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust"
HALPHA_URL = "https://lambda.gsfc.nasa.gov/data/foregrounds/halpha/lambda_halpha_fwhm06_0512.fits"
HALPHA_PATH = AUX_DIR / "lambda_halpha_fwhm06_0512.fits"
DUST_PATH = AUX_DIR / "irsa_dust.json"

#: Columns kept from gaia_source for the search samples (a full row is ~150).
SEARCH_COLS = ", ".join([
    "source_id", "ra", "dec", "l", "b", "ecl_lat", "parallax", "parallax_error",
    "parallax_over_error", "pm", "pmra", "pmdec", "pmra_error", "pmdec_error",
    "parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr",
    "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "bp_rp",
    "phot_bp_rp_excess_factor", "ruwe", "astrometric_excess_noise",
    "astrometric_excess_noise_sig", "astrometric_params_solved", "astrometric_sigma5d_max",
    "visibility_periods_used", "ipd_frac_multi_peak", "ipd_gof_harmonic_amplitude",
    "ipd_frac_odd_win", "duplicated_source", "nu_eff_used_in_astrometry", "pseudocolour",
    "radial_velocity", "radial_velocity_error", "random_index"])

#: Both signs at once: |parallax| significant and pm > 600 km/s * |parallax| / 4.74.
SEARCH_QUERY = (f"select {{cols}} from {{table}} where pm > {PM_MIN} and phot_g_mean_mag < {G_MAX} "
                f"and abs(parallax_over_error) > {PLX_SNR_MIN} "
                f"and 4.74*pm > {V_TAN_MIN}*abs(parallax)")

#: The same cut in the columns gaia_source_lite has (no `pm`, so squared):
#: pm^2 > PM_MIN^2 and (4.74 pm)^2 > (V_TAN_MIN parallax)^2.
SEARCH_LITE_WHERE = (f"pmra*pmra + pmdec*pmdec > {PM_MIN ** 2} and phot_g_mean_mag < {G_MAX} "
                     f"and abs(parallax_over_error) > {PLX_SNR_MIN} "
                     f"and 22.4676*(pmra*pmra + pmdec*pmdec) > {V_TAN_MIN ** 2}*parallax*parallax")

#: El-Badry et al.'s query in gaia_source_lite's columns (pm from pmra, pmdec),
#: identical cuts; used only if the verbatim query cannot be run.
_PM = "sqrt(pmra*pmra + pmdec*pmdec)"
ELBADRY_LITE_WHERE = (f"phot_g_mean_mag < 20 and ruwe < 1.4 and pmra*pmra + pmdec*pmdec > 2500 "
                      f"and (4.74*{_PM}/(parallax + parallax_error) > 400 or (parallax + parallax_error) < 0) "
                      f"and (4.74*{_PM}/parallax > 600 or parallax < 0) "
                      f"and bp_rp < 0.5 and parallax_over_error < 5")

BACKGROUND_QUERY = (f"select source_id, phot_g_mean_mag, bp_rp, parallax, parallax_error, ruwe "
                    f"from {{table}} where parallax > {BACKGROUND_PLX_MIN} "
                    f"and parallax_over_error > {BACKGROUND_PLX_SNR}")

PRECISION_QUERY = (f"select source_id, l, b, phot_g_mean_mag, bp_rp, parallax_error, pmra_error, "
                   f"pmdec_error, ruwe, astrometric_params_solved, visibility_periods_used "
                   f"from {{table}} where random_index < {PRECISION_RANDOM_INDEX_MAX}")

ARI_ONLY = [s for s in GAIA_SERVICES if s[0].startswith("ARI")]
ARI_OR_ESA = [s for s in GAIA_SERVICES if s[0].startswith(("ARI", "ESA"))]


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------

def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _load_manifest() -> dict:
    if MANIFEST_PATH.exists():
        return json.loads(MANIFEST_PATH.read_text())
    return {"files": {}}


def record(path: Path, url: str, source: str, **extra) -> None:
    manifest = _load_manifest()
    data = path.read_bytes()
    manifest["files"][str(path.relative_to(DATA_DIR))] = {
        "url": url, "source": source, "retrieved_utc": _now_utc(), "bytes": len(data),
        "sha256": hashlib.sha256(data).hexdigest(), **extra}
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))


def http_get(url: str, *, params: dict | None = None) -> bytes:
    """GET with retries (DNS on this machine fails intermittently)."""
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code < 500:
                raise
            last = err
        except (urllib.error.URLError, TimeoutError, ConnectionError) as err:
            last = err
        time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"giving up on {url}: {last}")


def _gaia_record(name: str) -> None:
    path = GAIA_DIR / f"{name}.fits"
    prov = gaia_tap.provenance(name)
    if path.exists() and prov:
        record(path, prov["url"], f"Gaia DR3 via {prov['service']}", adql=prov["adql"],
               rows=prov["rows"])


# ---------------------------------------------------------------------------
# Gaia
# ---------------------------------------------------------------------------

def census_ids() -> list[int]:
    return [r.source_id for r in CENSUS]


def fetch_gaia(quiet: bool = False) -> None:
    q = gaia_tap.cached_query
    # 1. the known runaways and El-Badry et al.'s 25 candidates: full rows
    ids = sorted(set(census_ids()) | {i for i, _ in ELBADRY_TABLE1})
    q("known", f"select * from {{table}} where source_id in ({gaia_tap.id_list(ids)})",
      quiet=quiet)
    _gaia_record("known")
    # 2. Shen et al. (2018)'s DR2 candidates -> DR3, and their DR2 records
    dr2 = gaia_tap.id_list([s[1] for s in SHEN18_TABLE1])
    q("shen18_dr3", "select n.dr2_source_id, n.dr3_source_id, n.angular_distance, "
      "n.magnitude_difference, n.proper_motion_propagation, g.* "
      "from gaiadr3.dr2_neighbourhood as n join gaiadr3.gaia_source as g "
      f"on g.source_id = n.dr3_source_id where n.dr2_source_id in ({dr2})",
      services=ARI_OR_ESA, quiet=quiet)
    _gaia_record("shen18_dr3")
    q("shen18_dr2", f"select * from gaiadr2.gaia_source where source_id in ({dr2})",
      services=ARI_OR_ESA, quiet=quiet)
    _gaia_record("shen18_dr2")
    # 3. HR-diagram background and the parallax-precision sample
    q("background_100pc", BACKGROUND_QUERY, quiet=quiet)
    _gaia_record("background_100pc")
    q("precision", PRECISION_QUERY, quiet=quiet)
    _gaia_record("precision")
    # 4. El-Badry et al. (2023)'s query, verbatim (a whole-catalogue scan)
    try:
        q("elbadry_query", ELBADRY_QUERY, cols="*", quiet=quiet, scan=True)
    except RuntimeError as exc:
        if not quiet:
            print(f"    verbatim query failed ({exc}); same cuts on gaia_source_lite", flush=True)
        gaia_tap.lite_then_full("elbadry_query", ELBADRY_LITE_WHERE, cols="*", quiet=quiet)
    _gaia_record("elbadry_query")
    # 5. our broad search and its negative-parallax mirror (one scan of ARI's
    #    narrow gaia_source_lite table, then full rows by source_id)
    gaia_tap.lite_then_full("search", SEARCH_LITE_WHERE, cols=SEARCH_COLS, quiet=quiet)
    _gaia_record("search")
    # 6. Rybizki et al. (2022) fidelity for every candidate (ARI mirror only)
    fetch_fidelity(quiet)


def fetch_fidelity(quiet: bool = False) -> None:
    path = GAIA_DIR / "fidelity.fits"
    if path.exists():
        return
    ids = set()
    for name in ("search", "elbadry_query", "known"):
        t = Table.read(GAIA_DIR / f"{name}.fits")
        ids |= {int(i) for i in t["source_id"]}
    ids = sorted(ids)
    parts = []
    for k in range(0, len(ids), 400):
        chunk = ids[k:k + 400]
        t, prov = gaia_tap.run_query(
            "select source_id, fidelity_v2, dist_nearest_neighbor_at_least_m2_brighter, "
            "dist_nearest_neighbor_at_least_0_brighter, dist_nearest_neighbor_at_least_2_brighter, "
            "dist_nearest_neighbor_at_least_4_brighter, norm_dg "
            f"from gaiaedr3_complements.spurious_astrometry where source_id in ({gaia_tap.id_list(chunk)})",
            services=ARI_ONLY, quiet=quiet)
        parts.append(t)
    table = vstack(parts)
    table.write(path, format="fits", overwrite=True)
    (GAIA_DIR / "fidelity.json").write_text(json.dumps(dict(
        service=prov["service"], url=prov["url"], table="gaiaedr3_complements.spurious_astrometry",
        reference="Rybizki et al. 2022, MNRAS 510, 2597", rows=len(table), requested=len(ids),
        retrieved=_now_utc()), indent=1))
    record(path, prov["url"], "Rybizki et al. 2022 fidelity, via ARI", rows=len(table))
    if not quiet:
        print(f"  fidelity: {len(table):,} of {len(ids):,} sources", flush=True)


# ---------------------------------------------------------------------------
# Remnants, white dwarf models, dust
# ---------------------------------------------------------------------------

def fetch_green(quiet: bool = False) -> None:
    if GREEN_PATH.exists():
        return
    GREEN_PATH.write_bytes(http_get(GREEN_URL))
    record(GREEN_PATH, GREEN_URL, "Green 2025, J. Astrophys. Astron. 46, 14 (VizieR VII/297)")
    if not quiet:
        print("  Green (2025) remnant catalogue", flush=True)


def fetch_bedard(quiet: bool = False) -> None:
    for m in BEDARD_MASSES:
        path = bedard_path(m, "thick")
        if path.exists():
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        url = BEDARD_URL.format(name=path.name)
        path.write_bytes(http_get(url))
        record(path, url, "Bedard et al. 2020, ApJ 901, 93 (Montreal cooling sequences)")
        time.sleep(0.3)
    if not quiet:
        print(f"  Bedard et al. (2020) sequences: {len(BEDARD_MASSES)} masses", flush=True)


def _parse_dust(xml: str) -> dict:
    out = {}
    for key in ("refPixelValueSandF", "refPixelValueSFD", "meanValueSandF", "meanValueSFD"):
        m = re.search(rf"<{key}>\s*([-\d.]+)", xml)
        out[key] = float(m.group(1)) if m else float("nan")
    return out


def fetch_dust(quiet: bool = False) -> None:
    """E(B-V) toward every star placed in the HR diagram (the census, El-Badry
    et al.'s candidates and Shen et al.'s DR2 candidates)."""
    cache = json.loads(DUST_PATH.read_text()) if DUST_PATH.exists() else {}
    rows = Table.read(GAIA_DIR / "known.fits")
    extra = Table.read(GAIA_DIR / "shen18_dr3.fits")
    targets = [(int(r["source_id"]), float(r["ra"]), float(r["dec"])) for r in rows]
    targets += [(int(r["source_id"]), float(r["ra"]), float(r["dec"])) for r in extra]
    new = 0
    for sid, ra, dec in targets:
        if str(sid) in cache:
            continue
        xml = http_get(DUST_URL, params={"locstr": f"{ra:.5f} {dec:.5f} equ j2000",
                                         "regSize": "2.0"}).decode("utf-8", "replace")
        cache[str(sid)] = dict(ra=ra, dec=dec, **_parse_dust(xml))
        new += 1
        time.sleep(0.3)
    DUST_PATH.write_text(json.dumps(cache, indent=1, sort_keys=True))
    if new:
        record(DUST_PATH, DUST_URL, "IRSA DUST service: SFD 1998 and Schlafly & Finkbeiner 2011",
               positions=len(cache))
    if not quiet:
        print(f"  dust: {len(cache)} positions ({new} new)", flush=True)


def fetch_halpha(quiet: bool = False) -> None:
    """Finkbeiner (2003) all-sky H-alpha map (HEALPix nside 512, 6' FWHM),
    from NASA's LAMBDA archive: the background for D6-2's remnant."""
    if HALPHA_PATH.exists():
        return
    HALPHA_PATH.write_bytes(http_get(HALPHA_URL))
    record(HALPHA_PATH, HALPHA_URL, "Finkbeiner 2003, ApJS 146, 407 (LAMBDA)")
    if not quiet:
        print("  Finkbeiner (2003) H-alpha map", flush=True)


GROUPS = {"gaia": fetch_gaia, "green": fetch_green, "bedard": fetch_bedard, "dust": fetch_dust,
          "halpha": fetch_halpha}


def ensure_data(quiet: bool = False, only=None) -> None:
    for name, fn in GROUPS.items():
        if only and name not in only:
            continue
        fn(quiet=quiet)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", action="append", choices=list(GROUPS))
    args = ap.parse_args(argv)
    ensure_data(only=args.only)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
