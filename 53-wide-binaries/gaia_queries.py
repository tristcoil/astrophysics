"""The Gaia DR3 queries of this project, each cached under data/gaia/.

1. ``nearby_lean``: every Gaia DR3 source with parallax > 10/3 mas (closer
   than ~300 pc) and G < 20, all sky, with only the columns needed to find
   pairs. This is the parent sample of Pittordis & Sutherland (2023), whose
   pairs Banik et al. (2024) analysed (P&S: parallax > 10/3 mas, G < 17,
   |b| > 15 deg for the binary members; G < 20 for the "faint star" sample
   used to veto triples). The magnitude and latitude cuts are applied later,
   so they can be varied.
2. ``pair_stars_full``: all gaia_source columns we use, for the stars of the
   candidate pairs only (fetched by source_id, which is indexed).
3. ``pair_stars_flame``: Gaia DR3 FLAME masses for the same stars.
"""

from __future__ import annotations

import numpy as np

from config import GAIA_SERVICES
from gaia_tap import cached_query, chunked_scan, id_list

PLX_MIN_MAS = 10.0 / 3.0
G_MAX_FAINT = 20.0

LEAN_COLS = ("source_id, ra, dec, l, b, parallax, parallax_error, pmra, pmra_error, pmdec, "
             "pmdec_error, phot_g_mean_mag, ruwe")

FULL_COLS = (
    "source_id, ref_epoch, ra, ra_error, dec, dec_error, l, b, ecl_lat, parallax, parallax_error, "
    "pmra, pmra_error, pmdec, pmdec_error, ra_dec_corr, ra_parallax_corr, ra_pmra_corr, "
    "ra_pmdec_corr, dec_parallax_corr, dec_pmra_corr, dec_pmdec_corr, parallax_pmra_corr, "
    "parallax_pmdec_corr, pmra_pmdec_corr, astrometric_n_obs_al, astrometric_n_good_obs_al, "
    "astrometric_chi2_al, astrometric_excess_noise, astrometric_excess_noise_sig, "
    "astrometric_params_solved, nu_eff_used_in_astrometry, pseudocolour, visibility_periods_used, "
    "ruwe, ipd_gof_harmonic_amplitude, ipd_frac_multi_peak, ipd_frac_odd_win, duplicated_source, "
    "phot_g_mean_mag, phot_g_mean_flux_over_error, phot_bp_mean_mag, phot_rp_mean_mag, "
    "phot_bp_mean_flux_over_error, phot_rp_mean_flux_over_error, phot_bp_rp_excess_factor, "
    "phot_bp_n_blended_transits, phot_rp_n_blended_transits, bp_rp, radial_velocity, "
    "radial_velocity_error, rv_nb_transits, rv_template_teff, non_single_star"
)


def nearby_lean(quiet: bool = False):
    """Step 1 (5.8 M rows, ~10 s per patch on ESA on 9 Oct 2026)."""
    esa = next(s for s in GAIA_SERVICES if s[0].startswith("ESA"))
    where = (f"parallax > {PLX_MIN_MAS!r} and phot_g_mean_mag < {G_MAX_FAINT!r} "
             "and source_id >= {lo} and source_id < {hi}")
    return chunked_scan("nearby_lean", "select {cols} from {table} where " + where,
                        cols=LEAN_COLS, service=esa, base_level=1, max_level=4, workers=3,
                        quiet=quiet)


def full_rows(name: str, source_ids, *, cols: str = FULL_COLS, batch: int = 2000,
              quiet: bool = False):
    """Full gaia_source rows for ``source_ids`` (cached as data/gaia/<name>.fits)."""
    from astropy.table import Table, vstack
    from config import GAIA_DIR
    import json
    path = GAIA_DIR / f"{name}.fits"
    if path.exists() and path.with_suffix(".json").exists():
        return Table.read(path)
    sids = np.unique(np.asarray(source_ids, dtype=np.int64))
    parts = []
    for k in range(0, len(sids), batch):
        sub = f"{name}_part{k // batch:04d}"
        q = f"select {{cols}} from {{table}} where source_id in ({id_list(sids[k:k + batch])})"
        parts.append(cached_query(sub, q, cols=cols, quiet=True))
        if not quiet and (k // batch) % 10 == 0:
            print(f"    {name}: {min(k + batch, len(sids)):,} / {len(sids):,}", flush=True)
    full = vstack(parts, metadata_conflicts="silent")
    if len(full) != len(sids):
        raise RuntimeError(f"{name}: asked for {len(sids)} sources, got {len(full)}")
    full.write(path, format="fits", overwrite=True)
    path.with_suffix(".json").write_text(json.dumps(dict(
        rows=len(full), batches=len(parts), cols=cols,
        note="fetched by source_id in batches; each batch's provenance is in its _partNNNN.json"),
        indent=1))
    for k in range(len(parts)):  # the stacked file replaces the batch files
        for ext in (".fits",):
            (GAIA_DIR / f"{name}_part{k:04d}{ext}").unlink(missing_ok=True)
    return full


if __name__ == "__main__":
    t = nearby_lean()
    print(len(t))
