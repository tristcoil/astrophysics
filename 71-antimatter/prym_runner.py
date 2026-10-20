"""Run PRyMordial for a list of baryon-to-photon ratios; print one JSON line each.

Meant to be launched by ``bbn.py`` as an isolated subprocess:

    python -I prym_runner.py <PRyMordial dir> <primat|nacreii> <small|full> grid eta1 eta2 ...
    python -I prym_runner.py <PRyMordial dir> <primat|nacreii> <small|full> mc eta n seed

In ``mc`` mode the 12 key nuclear rates are drawn from PRyMordial's own
log-normal uncertainties (rate = median x expsigma^p, p ~ N(0, 1)) and the
neutron lifetime from 878.4 +- 0.5 s (PDG), to give the theory spread of D/H.

PRyMordial reads its rate tables relative to the current directory and is
configured through module-level flags, so each rate set runs in its own
process started inside the PRyMordial folder. Nothing is written to disk.
Output per eta: Neff, Y_p (BBN mass fraction), D/H, 3He/H, 7Li/H, run time.
"""

from __future__ import annotations

import json
import os
import sys
import time


KEY_REACTIONS = ("npdg", "dpHe3g", "ddHe3n", "ddtp", "tpag", "tdan", "taLi7g", "He3ntp",
                 "He3dap", "He3aBe7g", "Be7nLi7p", "Li7paa")


def main(argv: list[str]) -> int:
    prym_dir, rates, net, mode = argv[0], argv[1], argv[2], argv[3]
    if mode == "grid":
        jobs = [(float(v), {}) for v in argv[4:]]
    else:
        import numpy as np

        eta, n, seed = float(argv[4]), int(argv[5]), int(argv[6])
        rng = np.random.default_rng(seed)
        jobs = []
        for _ in range(n):
            draw = {f"p_{k}": float(rng.standard_normal()) for k in KEY_REACTIONS}
            draw["tau_n"] = 878.4 + 0.5 * float(rng.standard_normal())
            jobs.append((eta, draw))
    os.chdir(prym_dir)
    sys.path.insert(0, prym_dir)
    import PRyM.PRyM_init as ini

    ini.verbose_flag = False
    ini.numba_flag = True
    ini.aTid_flag = True
    ini.compute_bckg_flag = True
    ini.compute_nTOp_flag = True
    ini.compute_nTOp_thermal_flag = False
    ini.save_bckg_flag = False
    ini.save_nTOp_flag = False
    ini.save_nTOp_thermal_flag = False
    ini.julia_flag = False
    ini.smallnet_flag = (net == "small")
    ini.nacreii_flag = (rates == "nacreii")
    ini.ReloadKeyRates()
    import PRyM.PRyM_main as main_mod

    for eta, draw in jobs:
        t0 = time.time()
        ini.eta0b = eta
        ini.Omegabh2 = eta / ini.Omegabh2_to_eta0b
        for k, v in draw.items():
            setattr(ini, k, v)
        res = main_mod.PRyMclass().PRyMresults()
        print(json.dumps({"eta": eta, "rates": rates, "net": net, "draw": draw,
                          "Neff": float(res[0]),
                          "Yp_cmb": float(res[3]), "Yp_bbn": float(res[4]),
                          "DoH_e5": float(res[5]), "He3oH_e5": float(res[6]),
                          "Li7oH_e10": float(res[7]), "seconds": time.time() - t0,
                          "tau_n": float(ini.tau_n)}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
