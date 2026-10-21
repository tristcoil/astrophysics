"""Choose the red giants we analyse.

The sample is drawn at random (fixed seed) so that it is not cherry-picked,
but stratified so that it spans the whole range red giants occupy:

* the parent list is APOKASC-3's Gold sample (numax and Dnu agreed by 5+ of
  their 10 pipelines, an APOGEE spectrum, and an evolutionary state from
  seismology or spectroscopy) that Yu et al. (2018) also measured, so that
  every star can be checked against both catalogues;
* at least ``MIN_QUARTERS`` quarters of Kepler long-cadence data (~3.5-4
  years) and Kepler magnitude ``Kp <= MAX_KEPMAG``;
* ten bins in Yu et al.'s numax from 5 to 230 uHz, ``STARS_PER_BIN`` stars
  per bin; where red-clump stars exist, half of each bin's draw is red clump
  and half red giant branch.

Catalogue numax is used *only* to stratify the draw; our pipeline never sees
it. The four literature showcase stars in ``config.SHOWCASE_STARS`` are added
by hand.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import catalogs as K
import config as C


def eligible() -> pd.DataFrame:
    a = K.apokasc3()
    y = K.yu2018()
    j = a.join(y[["numax", "dnu", "kepmag", "n_quarters", "state"]], rsuffix="_yu", how="inner")
    ok = ((j["category"] == "Gold") & j["state"].isin(["RGB", "RC"])
          & (j["n_quarters_yu"] >= C.MIN_QUARTERS) & (j["n_quarters"] >= C.MIN_QUARTERS)
          & (j["kepmag"] <= C.MAX_KEPMAG))
    return j[ok]


def build_sample() -> pd.DataFrame:
    rng = np.random.default_rng(C.SAMPLE_SEED)
    pool = eligible()
    edges = np.asarray(C.NUMAX_BIN_EDGES)
    rows = []
    for i in range(len(edges) - 1):
        in_bin = pool[(pool["numax_yu"] >= edges[i]) & (pool["numax_yu"] < edges[i + 1])]
        rc = in_bin[in_bin["state"] == "RC"]
        rgb = in_bin[in_bin["state"] == "RGB"]
        n_rc = C.STARS_PER_BIN // 2 if len(rc) >= C.STARS_PER_BIN // 2 else 0
        n_rgb = C.STARS_PER_BIN - n_rc
        for grp, n in ((rgb, n_rgb), (rc, n_rc)):
            if n == 0:
                continue
            pick = rng.choice(np.sort(grp.index.to_numpy()), size=min(n, len(grp)), replace=False)
            for k in pick:
                rows.append({"kic": int(k), "numax_bin": i, "selection": "random"})
    chosen = {r["kic"] for r in rows}
    for k in C.SHOWCASE_STARS:
        if k not in chosen:
            rows.append({"kic": int(k), "numax_bin": -1, "selection": "showcase"})
        else:
            for r in rows:
                if r["kic"] == k:
                    r["selection"] = "random+showcase"
    s = pd.DataFrame(rows).set_index("kic")
    a = K.apokasc3()
    y = K.yu2018()
    s["state_apokasc"] = a.loc[s.index, "state"]
    s["state_yu"] = y.loc[s.index, "state"]
    s["kepmag"] = y.loc[s.index, "kepmag"]
    s["n_quarters"] = a.loc[s.index, "n_quarters"]
    s["numax_yu"] = y.loc[s.index, "numax"]
    s["teff_apogee"] = a.loc[s.index, "teff"]
    s["feh_apogee"] = a.loc[s.index, "feh"]
    s["gaia_dr3"] = a.loc[s.index, "gaia_dr3"].astype("Int64")
    s["tmass"] = a.loc[s.index, "tmass"]
    return s.sort_values("numax_yu")


def load_or_build_sample(*, quiet: bool = False) -> pd.DataFrame:
    if C.SAMPLE_PATH.exists():
        s = pd.read_csv(C.SAMPLE_PATH, dtype={"gaia_dr3": "Int64"})
    else:
        s = build_sample().reset_index()
        s.to_csv(C.SAMPLE_PATH, index=False)
        if not quiet:
            n = s["selection"].value_counts().to_dict()
            print(f"  sample: {len(s)} stars ({n}); "
                  f"{(s['state_apokasc'] == 'RC').sum()} red clump, "
                  f"{(s['state_apokasc'] == 'RGB').sum()} red giant branch")
    return s


if __name__ == "__main__":
    s = load_or_build_sample()
    print(s.to_string())
