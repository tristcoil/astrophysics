"""Stage 2: H0 from the published lens models, first with their delays, then with ours.

    python analysis_h0.py

1. Reproduce H0LiCOW (Wong et al. 2020): each lens's published distance
   posterior -> H0 in flat LCDM, and the six combined.
2. Swap in our delays (``measure_delays.py``): for each lens we re-measured,
   the published distance posterior is rescaled by r = (our delays / the
   delays the model was fitted to), and the delay part of its width is
   replaced by ours (``h0.delay_scale``, ``h0.swap_delays``).
3. TDCOSMO 2025: the released chains, the mass-sheet parameter lambda_int,
   and an estimate of how our delays would move its H0.

Writes ``results/h0.json`` and caches the grid log-likelihoods in
``results/h0_grid_*.npz``.
"""

from __future__ import annotations

import json

import h5py
import numpy as np

import config as C
import h0 as H

OUT = C.RESULTS

#: Independent delay pairs used to compare our delays with the model's.
#: (Quads publish all pairs, but only N-1 are independent; we keep the most
#: precise independent set. WFI 2033's A1 and A2 are merged into A in the
#: Euler/SMARTS data, so the model's A1B and A2B delays are averaged.)
RATIO_PAIRS = {
    "RXJ1131": {("B", "A"): None, ("B", "C"): None, ("B", "D"): None},
    "HE0435": {("A", "B"): None, ("A", "C"): None, ("A", "D"): None},
    "PG1115": {("A", "C"): None, ("B", "C"): None},
    "WFI2033": {("B", "A"): (36.75, 2.25), ("B", "C"): None},
    "J1206": {("A", "B"): None},
    "DES0408": {("A", "B"): None, ("A", "D"): None},
}


def published_pairs(key: str):
    """(pairs, values, symmetric errors) of the delays the published model used."""
    L = C.LENSES[key]
    vals, errs = {}, {}
    for pair, override in RATIO_PAIRS[key].items():
        if override is not None:
            vals[pair], errs[pair] = override
        else:
            v, em, ep = L["model_delays"][pair]
            vals[pair], errs[pair] = v, 0.5 * (em + ep)
    return vals, errs


def our_pairs(key: str):
    """Our delays for the same pairs, with their covariance (from the combined joint-spline result)."""
    res = json.loads((OUT / "delays" / f"{key}_combined.json").read_text())
    ref, others = res["ref"], res["others"]
    cov = np.array(res["cov"])
    tau = np.array([res["taus"][o] for o in others])
    pairs = list(RATIO_PAIRS[key])
    T = np.zeros((len(pairs), len(others)))
    for i, (x, y) in enumerate(pairs):
        if y != ref:
            T[i, others.index(y)] += 1.0
        if x != ref:
            T[i, others.index(x)] -= 1.0
    d = T @ tau
    return dict(zip(pairs, d.tolist())), T @ cov @ T.T


def tdcosmo_lens(name: str) -> dict:
    """One lens's likelihood input from the TDCOSMO 2025 release (a pickled dict).

    Objects of classes we do not have installed (hierArc helpers) are replaced
    by inert stubs; only the plain arrays are used.
    """
    import pickle

    class _Stub:
        def __init__(self, *a, **k):
            pass

        def __setstate__(self, state):
            self.__dict__["state"] = state

    class _U(pickle.Unpickler):
        def find_class(self, mod, cls):
            try:
                return super().find_class(mod, cls)
            except Exception:
                return type(cls, (_Stub,), {})

    with open(C.DATA / "tdcosmo2025" / "TDCOSMO_sample" / f"{name}_const_processed.pkl", "rb") as f:
        return _U(f).load()


def _wquantile(x, w, qs):
    o = np.argsort(x)
    cw = np.cumsum(w[o])
    cw /= cw[-1]
    return np.interp(qs, cw, x[o])


def grid_cache(name: str, builder):
    p = OUT / f"h0_grid_{name}.npz"
    if p.exists():
        return np.load(p)["logl"]
    logl = builder()
    OUT.mkdir(exist_ok=True)
    np.savez_compressed(p, logl=logl)
    return logl


def summ(post: H.H0Posterior) -> list:
    return [round(v, 3) for v in post.summary()]


def main():
    out = {"published": {}, "ours": {}, "swap": {}}
    logl_pub, logl_ours = {}, {}
    liks = {}
    for key in C.H0LICOW_SIX:
        def build(k=key):
            liks[k] = H.load_h0licow(k)
            return H.grid_logl(liks[k])
        logl_pub[key] = grid_cache(f"pub_{key}", build)
        out["published"][key] = summ(H.posterior([logl_pub[key]]))
        print(f"{key:8s} published H0 = {out['published'][key]}   Wong+2020 {C.LENSES[key]['wong20_h0']}")
    allpub = H.posterior(list(logl_pub.values()))
    out["published"]["ALL"] = summ(allpub)
    print(f"ALL six  published H0 = {out['published']['ALL']}   Wong+2020 {C.H0_LIT['H0LiCOW 2020']}")

    # --- our delays ---
    for key in C.MEASURED:
        pub, pub_err = published_pairs(key)
        ours, cov = our_pairs(key)
        sc = H.delay_scale(ours, cov, pub, pub_err)
        frac_pub = C.LENSES[key]["ddt_delay_frac"] or sc["frac_pub"]
        sc.update(frac_pub_used=frac_pub, ours={f"{a}{b}": v for (a, b), v in ours.items()},
                  ours_err={f"{a}{b}": float(np.sqrt(cov[i, i])) for i, (a, b) in enumerate(ours)},
                  pub={f"{a}{b}": v for (a, b), v in pub.items()}, pub_err={f"{a}{b}": v for (a, b), v in pub_err.items()})
        sc["pairs"] = [f"{a}{b}" for a, b in sc["pairs"]]
        out["swap"][key] = sc
        print(f"{key:8s} delay scale r = {sc['r']:.4f} +- {sc['frac_ours'] * sc['r']:.4f}  "
              f"(delay share of D_dt error: ours {100 * sc['frac_ours']:.1f}%, published {100 * frac_pub:.1f}%)")
        if key not in C.H0LICOW_SIX:
            continue

        def build(k=key, s=sc, fp=frac_pub):
            lik = liks.get(k) or H.load_h0licow(k)
            return H.grid_logl(H.swap_delays(lik, s["r"], s["frac_ours"], fp))
        logl_ours[key] = grid_cache(f"ours_{key}_{sc['r']:.5f}_{sc['frac_ours']:.5f}", build)
        out["ours"][key] = summ(H.posterior([logl_ours[key]]))
        print(f"{key:8s} H0 with our delays = {out['ours'][key]}")
    logl_ours["B1608"] = logl_pub["B1608"]  # radio delays, not re-measured
    out["ours"]["B1608"] = out["published"]["B1608"]
    allours = H.posterior([logl_ours[k] for k in C.H0LICOW_SIX])
    out["ours"]["ALL"] = summ(allours)
    five = H.posterior([logl_ours[k] for k in C.H0LICOW_SIX if k != "B1608"])
    out["ours"]["FIVE_OPTICAL"] = summ(five)
    print(f"ALL six with our delays (B1608 as published) H0 = {out['ours']['ALL']}")
    np.savez_compressed(OUT / "h0_posteriors.npz", h0=H.H0_GRID,
                        **{f"pub_{k}": H.posterior([v]).p_h0 for k, v in logl_pub.items()},
                        **{f"ours_{k}": H.posterior([v]).p_h0 for k, v in logl_ours.items()},
                        pub_ALL=allpub.p_h0, ours_ALL=allours.p_h0, ours_FIVE=five.p_h0)

    # --- TDCOSMO 2025 ---
    td = {}
    for model in ("ULambdaCDM1", "LambdaCDM1a", "LambdaCDM1d", "LambdaCDM3a"):
        with h5py.File(C.DATA / "tdcosmo2025" / "chains_export" / f"{model}.h5") as f:
            s = f["samples"][:]
            names = [x.decode() for x in f["parameters"][:]]
            ds = f.attrs["dataset"]
        h = s[:, names.index("h0")]
        lam = s[:, names.index("lambda_mst")]
        q = np.percentile(h, [15.865525, 50, 84.134475])
        ql = np.percentile(lam, [15.865525, 50, 84.134475])
        td[model] = dict(dataset=str(ds), h0=[q[1], q[1] - q[0], q[2] - q[1]], lambda_int=[ql[1], ql[1] - ql[0], ql[2] - ql[1]],
                         corr_h0_lambda=float(np.corrcoef(h, lam)[0, 1]),
                         h0_over_lambda=list(np.percentile(h / lam, [15.865525, 50, 84.134475])))
        print(f"TDCOSMO 2025 {model} ({ds}): H0 {q[1]:.1f} -{q[1] - q[0]:.1f} +{q[2] - q[1]:.1f}, "
              f"lambda_int {ql[1]:.3f}, corr(H0, lambda) {td[model]['corr_h0_lambda']:.2f}")
    # Our delays' effect on TDCOSMO 2025: each lens's D_dt moves by r; H0 moves by
    # the inverse-variance-weighted mean of 1/r over the eight lenses (unmeasured
    # lenses: r = 1), weights from each lens's published D_dt precision.
    # The weights are each lens's D_dt precision in TDCOSMO 2025's own
    # likelihood inputs (power-law models, before kappa_ext and kinematics).
    w, invr, fracs = [], [], {}
    for k in ["B1608", "DES0408", "HE0435", "PG1115", "RXJ1131", "J1206", "WFI2033", "WGD2038"]:
        lik = tdcosmo_lens(C.LENSES[k]["tdcosmo"])
        x = np.asarray(lik["ddt_samples"], float)
        wt = lik.get("ddt_weights")
        wt = np.ones_like(x) if wt is None else np.asarray(wt, float)
        q16, q50, q84 = _wquantile(x, wt, [0.15865525, 0.5, 0.84134475])
        fracs[k] = 0.5 * (q84 - q16) / q50
        r = out["swap"][k]["r"] if k in out["swap"] else 1.0
        w.append(1 / fracs[k] ** 2)
        invr.append(1 / r)
    shift = float(np.average(invr, weights=w))
    td["ddt_frac_tdcosmo_inputs"] = fracs
    td["our_delay_shift_factor"] = shift
    td["LambdaCDM1a_with_our_delays_estimate"] = [td["LambdaCDM1a"]["h0"][0] * shift, td["LambdaCDM1a"]["h0"][1] * shift,
                                                  td["LambdaCDM1a"]["h0"][2] * shift]
    print(f"TDCOSMO 2025 (LCDM1a) H0 with our delays, estimated: x{shift:.4f} -> {td['LambdaCDM1a_with_our_delays_estimate'][0]:.1f}")
    out["tdcosmo2025"] = td
    (OUT / "h0.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
