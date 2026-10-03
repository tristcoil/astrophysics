"""Download and cache everything the analysis needs.

Two kinds of data are fetched, both free and public:

1. **Kepler light curves** from the Mikulski Archive for Space Telescopes (MAST),
   operated by STScI. These are the actual calibrated brightness measurements the
   spacecraft made. We use the PDCSAP flux, which is the pipeline's
   systematics-corrected photometry.

2. **Published system parameters** from the NASA Exoplanet Archive (Caltech/IPAC).
   These are used for two things only: the star's radius, mass and temperature,
   which the light curve alone cannot give us, and the published planet
   parameters we grade our own measurements against at the end.

Everything is cached under ``data/``. Re-running is cheap and works offline once
the cache is warm; nothing is re-downloaded that is already present.

Run directly to populate the cache::

    python fetch_data.py
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from datetime import datetime, timezone

import numpy as np

from config import CACHE_DIR, STELLAR_PARAMS_PATH, TARGETS, Target

warnings.filterwarnings("ignore", category=UserWarning)

#: Columns pulled from the NASA Exoplanet Archive's `pscomppars` table. That
#: table holds one composite row per planet, assembled from the literature, which
#: is what we want for a single set of reference values per system.
ARCHIVE_COLUMNS = [
    "pl_name",       # planet name
    "hostname",      # host star name
    "pl_orbper",     # orbital period, days
    "pl_orbpererr1", "pl_orbpererr2",
    "pl_rade",       # planet radius, Earth radii
    "pl_radeerr1", "pl_radeerr2",
    "pl_radj",       # planet radius, Jupiter radii
    "pl_trandep",    # transit depth, PERCENT in this table
    "pl_trandeperr1", "pl_trandeperr2",
    "pl_trandur",    # transit duration, hours
    "pl_orbsmax",    # semi-major axis, AU
    "pl_imppar",     # impact parameter
    "pl_eqt",        # equilibrium temperature, K
    "st_rad",        # stellar radius, solar radii
    "st_raderr1", "st_raderr2",
    "st_mass",       # stellar mass, solar masses
    "st_masserr1", "st_masserr2",
    "st_teff",       # stellar effective temperature, K
    "st_tefferr1", "st_tefferr2",
    "st_logg",       # surface gravity
    "sy_kepmag",     # Kepler-band magnitude
    "disc_year",     # year of discovery
    "disc_facility",
]


#: Cadences flagged by ``quality_bitmask='default'`` are those Kepler itself
#: marks as corrupted -- attitude tweaks, safe-mode recovery, thruster firings,
#: detector anomalies. The stricter 'hard'/'hardest' masks additionally discard
#: cadences flagged as impulsive outliers and rolling-band artefacts, which can
#: clip genuine short-duration transit signal, so 'default' is the right choice
#: for transit work.
QUALITY_BITMASK = "default"


def _cached_fits(target: Target) -> list[str]:
    """Find already-downloaded FITS files for this target, if any.

    Kepler filenames embed the KIC number zero-padded to nine digits, e.g.
    ``kplr011904151-2009350155506_llc.fits``, which is what we match on.
    """
    import glob
    import os

    digits = "".join(c for c in target.kic if c.isdigit())
    if not digits:
        return []
    prefix = f"kplr{int(digits):09d}"
    pattern = os.path.join(str(CACHE_DIR), "mastDownload", "Kepler", "**", "*_llc.fits")
    return sorted(f for f in glob.glob(pattern, recursive=True)
                  if os.path.basename(f).startswith(prefix))


def _stitch(collection, *, quiet: bool = False):
    """Normalise each quarter and concatenate, in time order.

    Two things matter here. First, ``stitch`` divides each quarter by its own
    median before concatenating: Kepler rolled 90 degrees every three months, so
    the star landed on a different detector each quarter and the raw flux level
    jumps by several percent between them. Left uncorrected those steps swamp any
    transit.

    Second, ``stitch`` concatenates in *collection order*, not time order, so the
    collection is sorted by quarter first. Out-of-order data would not break the
    fold -- which only cares about phase -- but it does break the gap detection
    used by the detrending filter.
    """
    import lightkurve as lk

    ordered = sorted(collection, key=lambda x: x.meta.get("QUARTER", 0))
    lc = lk.LightCurveCollection(ordered).stitch().remove_nans()
    return lc[np.argsort(lc.time.value)] if _out_of_order(lc) else lc


def _out_of_order(lc) -> bool:
    t = lc.time.value
    return bool(len(t) > 1 and (t[1:] < t[:-1]).any())


def fetch_light_curves(target: Target, *, quiet: bool = False, refresh: bool = False):
    """Load every available Kepler long-cadence light curve for a target.

    Returns a stitched, normalised ``LightCurve``.

    If the FITS files are already in ``data/mast/`` they are read straight from
    disk and **no network access happens at all**. This matters more than it
    might seem: ``search_lightcurve`` contacts MAST every time it is called, even
    when every file it would return is already cached, so a download-then-search
    pattern is not actually offline-capable. Reading the cached files directly is.

    A note on the two flux columns Kepler provides: SAP is raw aperture
    photometry; PDCSAP has had common-mode instrumental systematics removed by
    fitting trends shared across many stars on the same detector. PDCSAP is the
    right choice for transit work and is what lightkurve returns by default.
    """
    import lightkurve as lk

    cached = [] if refresh else _cached_fits(target)
    if cached:
        if not quiet:
            print(f"  reading {len(cached)} cached quarters (offline)")
        collection = lk.LightCurveCollection(
            [lk.read(f, quality_bitmask=QUALITY_BITMASK) for f in cached]
        )
    else:
        if not quiet:
            print(f"  searching MAST for {target.name} ({target.kic})...")
        search = lk.search_lightcurve(
            target.name, mission="Kepler", cadence="long", author="Kepler",
        )
        if len(search) == 0:
            raise RuntimeError(f"no Kepler light curves found for {target.name}")

        if target.quarters is not None:
            keep = [
                i for i, m in enumerate(search.mission)
                if any(f"Quarter {q:02d}" in m for q in target.quarters)
            ]
            search = search[keep]

        if not quiet:
            print(f"  {len(search)} quarters available; downloading to data/mast/")
        collection = search.download_all(
            download_dir=str(CACHE_DIR),
            quality_bitmask=QUALITY_BITMASK,
            flux_column="pdcsap_flux",
        )
        if collection is None or len(collection) == 0:
            raise RuntimeError(f"download returned nothing for {target.name}")

    lc = _stitch(collection, quiet=quiet)

    if not quiet:
        span = float(lc.time.value.max() - lc.time.value.min())
        print(f"  stitched: {len(lc):,} cadences over {span:.0f} days "
              f"({len(collection)} quarters)")
    return lc


def fetch_stellar_parameters(targets: list[Target], *, quiet: bool = False) -> dict:
    """Pull published system parameters from the NASA Exoplanet Archive.

    Cached to ``data/stellar_params.json`` with the retrieval date, so results
    stay reproducible even as the archive is updated.
    """
    from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

    out = {
        "_retrieved_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "_source": "NASA Exoplanet Archive, pscomppars table",
        "_note": (
            "pscomppars is the composite-parameters table: one row per planet, "
            "assembled from the literature. Values may come from different papers."
        ),
        "planets": {},
    }

    for target in targets:
        if not quiet:
            print(f"  querying archive for {target.planet}...")
        table = NasaExoplanetArchive.query_criteria(
            table="pscomppars",
            select=",".join(ARCHIVE_COLUMNS),
            where=f"pl_name='{target.planet}'",
        )
        if len(table) == 0:
            raise RuntimeError(f"archive returned no row for {target.planet}")

        row = table[0]
        record = {}
        for col in ARCHIVE_COLUMNS:
            if col not in table.colnames:
                continue
            value = row[col]
            try:
                # Strip astropy units to plain floats so the cache is plain JSON.
                value = value.value if hasattr(value, "value") else value
                if hasattr(value, "item"):
                    value = value.item()
            except Exception:
                pass
            if isinstance(value, float) and value != value:  # NaN
                value = None
            record[col] = value

        # The archive stores transit depth in percent; convert once, here, so no
        # downstream code has to remember the unit. This is a classic place to
        # introduce a factor-of-100 error.
        if record.get("pl_trandep") is not None:
            record["depth_fraction"] = record["pl_trandep"] / 100.0
            record["depth_ppm"] = record["pl_trandep"] * 1e4
        else:
            record["depth_fraction"] = None
            record["depth_ppm"] = None

        out["planets"][target.planet] = record

    STELLAR_PARAMS_PATH.write_text(json.dumps(out, indent=2))
    if not quiet:
        print(f"  wrote {STELLAR_PARAMS_PATH.relative_to(STELLAR_PARAMS_PATH.parent.parent)}")
    return out


def load_stellar_parameters() -> dict:
    """Read the cached archive parameters, erroring clearly if not fetched yet."""
    if not STELLAR_PARAMS_PATH.exists():
        raise FileNotFoundError(
            f"{STELLAR_PARAMS_PATH} not found -- run `python fetch_data.py` first"
        )
    return json.loads(STELLAR_PARAMS_PATH.read_text())


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "--target", action="append",
        help="fetch only this target (repeatable); default is all",
    )
    parser.add_argument(
        "--skip-archive", action="store_true",
        help="skip the NASA Exoplanet Archive query (light curves only)",
    )
    args = parser.parse_args(argv)

    targets = TARGETS
    if args.target:
        wanted = {t.lower() for t in args.target}
        targets = [t for t in TARGETS if t.name.lower() in wanted]
        if not targets:
            print(f"no targets matched {args.target}", file=sys.stderr)
            return 1

    print("Fetching Kepler light curves from MAST")
    for target in targets:
        print(f"\n{target.planet} -- {target.difficulty} case")
        lc = fetch_light_curves(target)
        del lc  # only warming the cache here; analysis reloads from disk

    if not args.skip_archive:
        print("\nFetching published parameters from the NASA Exoplanet Archive")
        fetch_stellar_parameters(targets)

    print("\nDone. Cache is warm; analysis can now run offline.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
