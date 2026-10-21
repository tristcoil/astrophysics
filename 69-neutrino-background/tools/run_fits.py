"""Run (part of) the CMB/BAO fit campaign and print a summary.

    python tools/run_fits.py                 # every rung
    python tools/run_fits.py tt plik         # only these rungs

Results are cached in data/derived/fits/; rerunning costs nothing.
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import cmb_fits as cf  # noqa: E402

rungs = tuple(sys.argv[1:]) or cf.RUNGS
t0 = time.monotonic()
out = cf.run(rungs=rungs)
for key, s in out["neff"].items():
    lp, iv = s["laplace"], s["interval68"]
    print(f"N_eff {key:18s} Laplace {lp['mean']:.3f} ± {lp['sigma']:.3f}; profile "
          f"{iv['best']:.3f} -{iv['minus']:.3f} +{iv['plus']:.3f}; vs {s['compare']}")
for key, s in out["mnu"].items():
    print(f"mnu   {key:18s} 95% < {s['limit95_flat_ge_0']} eV (gauss {s['limit95_flat_ge_0_gauss']:.3f});"
          f" parabola mu {s['parabola']['mu']:+.3f} sigma {s['parabola']['sigma']:.3f}; vs {s['compare']}")
print(f"done in {(time.monotonic() - t0) / 60:.1f} min")
