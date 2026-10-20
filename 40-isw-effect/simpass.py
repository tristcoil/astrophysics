"""Push the simulated CMB skies through every linear statistic and cache the results.

Each statistic in this project -- a stacked filter value, a template fit, a
galaxy-CMB cross-spectrum -- is linear in the CMB map, so its null distribution
comes from evaluating it on simulated skies that know nothing about galaxies.
Sky i is regenerated from seed RANDOM_SEED + i (cmb.sim_alm), so separate passes
see identical skies and their results can be combined with full covariance.

    python simpass.py bank        # every stacking filter (~20 min for 1000 skies)
"""

from __future__ import annotations

import sys
import time

import numpy as np

from config import DERIVED_DIR, N_SIMS

BANK_SIMS = DERIVED_DIR / "sims_bank.npy"


def bank_sims(n_sims: int = N_SIMS, *, force: bool = False) -> np.ndarray:
    if BANK_SIMS.exists() and not force:
        v = np.load(BANK_SIMS)
        if v.shape[0] >= n_sims:
            return v[:n_sims]
    import stacking as S
    import superstructures as SS

    layout = SS.build_layout()
    t0 = time.monotonic()
    v = S.run_sims(layout.bank, n_sims, log_every=50)
    np.save(BANK_SIMS, v)
    print(f"    {n_sims} skies through {len(layout.bank)} filters in {time.monotonic() - t0:.0f} s")
    return v


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "bank"
    if what == "bank":
        bank_sims(force=True)
