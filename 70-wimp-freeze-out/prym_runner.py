"""Run PRyMordial once with given switches; print one JSON line with its abundances.

    python -I prym_runner.py <PRyMordial dir> '<json config>'

Config keys (all optional): eta, tau_n, born (Born weak rates), small (12-reaction network),
aTid (non-instantaneous decoupling in a(T)), save_bckg (write its background table into the
PRyMordial folder's PRyMrates/thermo/). PRyMordial reads its tables relative to the current
directory and is configured through module globals, so each run is a fresh subprocess started
inside our own copy of the code (data/external). Launched by ``neutrons.prymordial``.
"""

from __future__ import annotations

import json
import os
import sys
import time


def main(argv: list[str]) -> int:
    prym_dir, cfg = argv[0], json.loads(argv[1])
    os.chdir(prym_dir)
    sys.path.insert(0, prym_dir)
    import PRyM.PRyM_init as ini

    ini.verbose_flag = False
    ini.numba_flag = True
    ini.julia_flag = False
    ini.aTid_flag = bool(cfg.get("aTid", True))
    ini.compute_bckg_flag = True
    ini.save_bckg_flag = bool(cfg.get("save_bckg", False))
    ini.compute_nTOp_flag = True
    ini.nTOpBorn_flag = bool(cfg.get("born", False))
    ini.compute_nTOp_thermal_flag = False
    ini.save_nTOp_flag = False
    ini.save_nTOp_thermal_flag = False
    ini.smallnet_flag = bool(cfg.get("small", False))
    ini.nacreii_flag = False
    ini.ReloadKeyRates()
    if "tau_n" in cfg:
        ini.tau_n = float(cfg["tau_n"]) * ini.second
    if "eta" in cfg:
        ini.eta0b = float(cfg["eta"])
        ini.Omegabh2 = ini.eta0b / ini.Omegabh2_to_eta0b
    import PRyM.PRyM_main as main_mod

    t0 = time.time()
    res = main_mod.PRyMclass().PRyMresults()
    print(json.dumps({"config": cfg, "Neff": float(res[0]), "Yp_cmb": float(res[3]),
                      "Yp_bbn": float(res[4]), "DoH": float(res[5]) * 1e-5,
                      "He3oH": float(res[6]) * 1e-5, "Li7oH": float(res[7]) * 1e-10,
                      "seconds": time.time() - t0}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
