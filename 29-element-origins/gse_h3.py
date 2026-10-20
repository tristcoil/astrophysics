"""Project H3's Gaia-Sausage-Enceladus stars, for cross-matching into GALAH.

H3 (folder ``*-galactic-archaeology``) picked GSE out of 3.2 million Gaia DR3
stars with full 3D motions using two published boxes, each in the potential
it was defined in:

* Helmi et al. (2018): -1500 < L_z < 150 kpc km/s, E > -1.8e5 km^2/s^2 in the
  Helmi et al. (2017) potential;
* Koppelman et al. (2019b): near-radial orbits, |circularity| < 0.2 and
  -1.5 < E/10^5 < -1.1 in the McMillan (2017) potential;

both among H3's halo stars (within 5 kpc, |v - v_LSR| > 210 km/s, RUWE < 1.4,
clean BP/RP photometry).

H3 did not save member lists, so we rerun H3's own code on H3's cached Gaia
file and keep the source IDs. H3's modules share names with ours (config,
fetch_data), so the work runs in a separate Python process rooted in H3's
folder. The result is cached in data/cache/h3_gse.npz.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np

# Our config is imported lazily: the worker process must import H3's instead.


def _cache() -> Path:
    from config import CACHE_DIR

    return CACHE_DIR / "h3_gse.npz"


def h3_folder() -> Path:
    from config import H3_GLOB, HERE

    hits = sorted(HERE.parent.glob(H3_GLOB))
    if not hits:
        raise FileNotFoundError(f"project H3 ({H3_GLOB}) not found next to {HERE}")
    return hits[0]


def load() -> dict:
    cache = _cache()
    if not cache.exists():
        h3 = h3_folder()
        cache.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(h3),
                        str(cache)], cwd=h3, check=True)
    with np.load(cache) as z:
        return {k: z[k] for k in z.files}


def _worker(h3: str, out: str) -> None:
    sys.path = [h3] + [p for p in sys.path if Path(p or ".").resolve() != Path(__file__).parent]
    import halo  # H3's module

    fast, _mp, _h17, _m17 = halo.load_samples()
    masks = halo.substructure_masks(fast)
    # H3's "halo5" set, exactly as its counts were made (halo.analyse).
    halo5 = (fast["good"] & (fast["distance"] < 5.0) & (fast["speed_lsr"] > halo.HALO_SPEED_KMS)
             & (fast["phot_bp_rp_excess_factor"] < 1.27))
    helmi = masks["GSE (Helmi+18)"] & halo5
    kopp = masks["GSE (Koppelman+19b)"] & halo5
    either = helmi | kopp
    np.savez(out, source_id=fast["source_id"][either], helmi2018=helmi[either],
             koppelman2019b=kopp[either], L_z=fast["L_z"][either], E_m17=fast["E_m17"][either],
             circ_m17=fast["circ_m17"][either], distance=fast["distance"][either],
             n_helmi=int(helmi.sum()), n_koppelman=int(kopp.sum()))
    print(f"H3 GSE: {helmi.sum():,} stars in Helmi et al.'s box, {kopp.sum():,} in "
          f"Koppelman et al.'s (H3's halo set: 5 kpc, |v - v_LSR| > 210 km/s)")


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--worker":
        _worker(sys.argv[2], sys.argv[3])
    else:
        d = load()
        print({k: (v.shape if v.ndim else v) for k, v in d.items()})
