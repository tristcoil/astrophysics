"""Check data/gaia/nearby_lean.fits against server-side COUNT(*)s, patch by
patch (each chunk file's source_id range is in its name), as the guard
against a silently truncated download."""
import json, re, warnings
import pyvo
from astropy.table import Table
from config import GAIA_DIR
from gaia_queries import PLX_MIN_MAS, G_MAX_FAINT
warnings.simplefilter("ignore")
svc = pyvo.dal.TAPService("https://gea.esac.esa.int/tap-server/tap")
out = {}
total_local = total_server = 0
for f in sorted((GAIA_DIR / "chunks" / "nearby_lean").glob("L*_*.fits")):
    lev, pix = map(int, re.match(r"L(\d+)_(\d+)", f.stem).groups())
    width = 2**35 * 4 ** (12 - lev)
    lo, hi = pix * width, (pix + 1) * width
    n_local = len(Table.read(f))
    q = (f"select count(*) from gaiadr3.gaia_source where parallax > {PLX_MIN_MAS!r} and "
         f"phot_g_mean_mag < {G_MAX_FAINT!r} and source_id >= {lo} and source_id < {hi}")
    for attempt in range(4):
        try:
            n_srv = int(svc.run_sync(q, timeout=600).to_table()[0][0])
            break
        except Exception as e:
            print(f.stem, "retry", type(e).__name__, flush=True)
    out[f.stem] = (n_local, n_srv)
    total_local += n_local; total_server += n_srv
    print(f.stem, n_local, n_srv, "OK" if n_local == n_srv else "MISMATCH", flush=True)
res = dict(total_local=total_local, total_server=total_server, per_patch=out,
           all_match=all(a == b for a, b in out.values()))
(GAIA_DIR / "nearby_lean_counts.json").write_text(json.dumps(res, indent=1))
print(res["total_local"], res["total_server"], res["all_match"])
