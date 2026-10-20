"""Run the whole project: data, analyses, theory, figures.

    python run.py                 # fetch what is missing, analyse, draw (first run ~40 min + download)
    python run.py --offline       # everything from the caches in data/ (~2 min)
    python run.py --refresh       # recompute every cached analysis step
    python run.py --toys 5000     # fewer pseudo-experiments (default 20,000)

Writes data/results.json (every number quoted in the README) and figures/.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time

import numpy as np

import config as cfg

FULL = cfg.CACHE_DIR / "full_results.json"
DIMUON_FITS = cfg.CACHE_DIR / "dimuon_fits.json"
BULKY = {"edges", "counts", "expected", "profile_mass", "profile_mu", "q0_bg_hist",
         "max_q0_bg_hist", "scan_masses", "scan_q0", "data_events", "same_sign_events",
         "signal_shapes"}


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def _trim(o):
    if isinstance(o, dict):
        return {k: _trim(v) for k, v in o.items() if k not in BULKY}
    if isinstance(o, list):
        return [_trim(v) for v in o]
    return o


def dimuon_fits(refresh: bool = False) -> dict:
    import dimuon
    if DIMUON_FITS.exists() and not refresh:
        return json.loads(DIMUON_FITS.read_text())
    d = dimuon.load()
    fits = dimuon.fit_all(d, "all")
    cats = dimuon.categories(d)
    out = {"fits": fits, "z": {}}
    for reg in ("barrel", "endcap", "all"):
        m = d["mass"][cats[reg]].astype(float)
        m = m[(m > 80) & (m < 102)]
        out["z"][reg] = dimuon.fit_z(m)
        out["z"][f"{reg}_width_free"] = dimuon.fit_z(m, gamma_free=True)
    out["cutflow"] = {k: int(d[k]) for k in ("n_events", "n_two_muons", "n_opposite_charge")}
    DIMUON_FITS.write_text(json.dumps(_jsonable(out), indent=1))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--offline", action="store_true", help="do not touch the network")
    ap.add_argument("--refresh", action="store_true", help="recompute all cached analyses")
    ap.add_argument("--toys", type=int, default=20000)
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()

    if not args.offline:
        import fetch_data
        fetch_data.fetch_all()

    import electroweak as ew
    import events
    import higgs_search
    import stability as st

    print("Event passes (cached after the first run):")
    events.process_all(refresh=args.refresh)
    events.dimuon_pass(refresh=args.refresh)
    events.dielectron_pass(refresh=args.refresh)

    print("Dimuon resonances ...")
    dm = dimuon_fits(args.refresh)
    print("Four-lepton Higgs search ...")
    higgs = higgs_search.run(n_toys=args.toys, refresh=args.refresh)
    print("Running couplings and vacuum stability ...")
    stab = st.summary(refresh=args.refresh)
    st.grids(refresh=args.refresh)
    print("Electroweak calendar ...")
    elw = {"potential": ew.higgs_potential(), "yukawas": ew.yukawas(),
           "atlas": ew.atlas_couplings(), "calendar": ew.calendar()}

    R = {"meta": {"project": "P7 The Higgs field switches on: Higgs bosons in CMS open data",
                  "python": platform.python_version(), "inputs_pdg2026": {
                      "M_h": cfg.M_HIGGS, "M_t": cfg.M_TOP, "M_W": cfg.M_W, "M_Z": cfg.M_Z,
                      "alpha_s": cfg.ALPHA_S_MZ, "G_F": cfg.G_FERMI},
                  "luminosity_fb": cfg.LUMINOSITY_PB / 1000},
         "dimuon": {"cutflow": dm["cutflow"],
                    "peaks": {k: v for fit in dm["fits"].values() for k, v in fit["peaks"].items()},
                    "z": {k: {kk: v[kk] for kk in ("mass", "mass_err", "ratio_minus_1", "width",
                                                    "width_err", "resolution", "chi2", "ndof")}
                          for k, v in dm["z"].items()}},
         "dimuon_fits": dm["fits"], "dimuon_z_barrel": dm["z"]["barrel"],
         "higgs": higgs, "stability": stab, "electroweak": elw}
    R = _jsonable(R)

    if not args.no_figures:
        print("Figures ...")
        import figures
        written = figures.make_all(R)
        print(f"  {len(written)} files in figures/")
    R["meta"]["runtime_s"] = round(time.time() - t0, 1)
    FULL.write_text(json.dumps(R, indent=1))
    cfg.RESULTS_PATH.write_text(json.dumps(_trim(R), indent=1))

    v = R["higgs"]["variants"]["ours_calibrated_with_zx"]
    s = R["stability"]
    h = s["results"][s["headline"]]
    print("\nHeadline numbers")
    print(f"  replica: every published bin reproduced = "
          f"{all(r['data_identical'] for r in R['higgs']['replica']['as_published'].values())}")
    print(f"  2e2mu collisions counted twice in the outreach analysis: "
          f"{R['higgs']['duplicates']['duplicated_events']}")
    print(f"  m_H = {v['mh']:.2f} (+{v['mh_hi'] - v['mh']:.2f} -{v['mh'] - v['mh_lo']:.2f}) "
          f"± {R['higgs']['mass_systematic']['total']:.2f} GeV; mu = {v['mu']:.2f}; "
          f"local {v['z_local_known_mass']:.2f} sigma at 125.13 GeV")
    print(f"  lambda < 0 above {h['Lambda_lambda']:.2e} GeV; stability needs M_t < "
          f"{h['stability_bound_mt']:.2f} GeV ({s['stability_gap']['n_sigma_exp_th']:.1f} sigma "
          f"away); log10 P(decayed in past light cone) = "
          f"{h['log10_p_decayed_in_past_light_cone']:.0f}")
    print(f"Done in {R['meta']['runtime_s']:.0f} s.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
