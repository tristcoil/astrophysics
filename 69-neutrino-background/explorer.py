"""Build the interactive "Weigh a Neutrino" page from the cached results.

    python explorer.py

Writes ``figures/neutrino_mass_explorer.html``: a single self-contained file
(fonts from Google Fonts, everything else inline) that works offline apart from
the fonts, which fall back to system faces. The template is a body fragment;
the same fragment with the data goes to
``data/derived/neutrino_mass_explorer_artifact.html`` for publishing as an
Artifact. The page embeds:

* NuFIT 6.0's best fits and 3σ ranges (from NuFIT's own Δχ² tables, as read
  by analysis.py), for both mass orderings;
* our profile likelihoods of Σm_ν (CAMB fits; see cmb_fits.py) and the
  published limits;
* a grid of CAMB linear matter power spectra for Σm_ν = 0 … 3 eV, so the page
  can show the suppression for any mass by interpolation;
* the Fermi–Dirac momentum quantiles and the relic temperature, for speeds.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402

import config as cfg  # noqa: E402

TEMPLATE = cfg.PROJECT_DIR / "explorer_template.html"
OUT = cfg.FIGURE_DIR / "neutrino_mass_explorer.html"
#: The same page as a bare fragment (no doctype or <head>), for publishing as an Artifact.
OUT_ARTIFACT = cfg.DERIVED_DIR / "neutrino_mass_explorer_artifact.html"
PK_GRID = cfg.DERIVED_DIR / "explorer_pk_grid.npz"
PK_MASSES = (0.0, 0.03, 0.06, 0.09, 0.12, 0.18, 0.24, 0.36, 0.48, 0.72, 1.0, 1.5, 2.0, 3.0)
#: Momentum quantiles of a relativistic Fermi-Dirac gas, x = p/kT (16, 50, 84%).
FD_QUANTILES = (1.5206, 2.83907, 4.77692)


def pk_grid(force: bool = False) -> dict:
    """Linear P(k) at z = 0 for three degenerate neutrinos of total mass Σ,
    at fixed Ω_m, ω_b, h and A_s (as figure 07), on a fine mass grid."""
    if PK_GRID.exists() and not force:
        with np.load(PK_GRID) as z:
            return {k: z[k] for k in z.files}
    import camb

    import cmb_model as cm
    import illustrations as il

    ref = cm.camb_params(il.REFERENCE, "one", want_cls=False)
    h = ref.H0 / 100.0
    omh2 = il.REFERENCE["ombh2"] + il.REFERENCE["omch2"] + ref.omnuh2
    kh = np.geomspace(1e-3, 5.0, 160)
    out = {"kh": kh, "masses": np.array(PK_MASSES)}
    pk0 = None
    for m in PK_MASSES:
        p = {k: v for k, v in il.REFERENCE.items() if k != "theta100"}
        p.update(H0=100 * h, mnu=m)
        pars = cm.camb_params(p, "degenerate", want_cls=False)
        pars.set_cosmology(ombh2=p["ombh2"], omch2=omh2 - p["ombh2"] - pars.omnuh2, mnu=m,
                           nnu=cfg.NEFF_SM, num_massive_neutrinos=3 if m > 0 else 0,
                           neutrino_hierarchy="degenerate", tau=p["tau"], H0=100 * h)
        pars.set_matter_power(redshifts=[0.0], kmax=10.0, nonlinear=False)
        pars.NonLinear = camb.model.NonLinear_none
        res = camb.get_results(pars)
        k, _, pk = res.get_matter_power_spectrum(minkh=1e-3, maxkh=5.0, npoints=160)
        pk = np.interp(kh, k, pk[0])
        if pk0 is None:
            pk0 = pk
        out[f"ratio_{m:g}"] = pk / pk0
        out[f"sigma8_{m:g}"] = np.array(res.get_sigma8_0())
        out[f"fnu_{m:g}"] = np.array(pars.omnuh2 / omh2)
    np.savez(PK_GRID, **out)
    return out


def _omega_m_h2() -> float:
    """ω_m of the reference universe used for the P(k) grid (Planck 2018's best fit)."""
    import cmb_model as cm
    import illustrations as il

    ref = cm.camb_params(il.REFERENCE, "one", want_cls=False)
    return float(il.REFERENCE["ombh2"] + il.REFERENCE["omch2"] + ref.omnuh2)


def _profile(fits: dict, key: str, nf: dict):
    import masses as ms

    s = fits.get("mnu", {}).get(key)
    if not s:
        return None
    x = np.array(s["profile"]["values"])
    y = np.array(s["profile"]["dchi2"])
    o = np.argsort(x)
    out = {"label": s["label"], "values": x[o].round(4).tolist(), "dchi2": y[o].round(3).tolist(),
           "limit95": s.get("limit95_flat_ge_0"), "parabola": s.get("parabola")}
    # Each ordering on its own (prior starting at its minimum sum), and the lightest mass
    # that limit allows.
    for order in ("NO", "IO"):
        lim = s.get(f"limit95_flat_ge_{order}")
        out[f"limit95_{order}"] = lim
        out[f"lightest95_{order}"] = (ms.lightest_from_sum(lim, nf[order]["dm21"]["best"],
                                                           nf[order]["dm3l"]["best"], order)
                                      if lim else None)
    return out


def _finite(obj):
    """Replace NaN and infinities by None: JavaScript's JSON.parse rejects them."""
    if isinstance(obj, dict):
        return {k: _finite(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_finite(v) for v in obj]
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    return obj


def build() -> Path:
    r = json.loads(cfg.RESULTS_PATH.read_text())
    nf = r["masses"]["nufit"]["TByes"]
    g = pk_grid()
    fits = r.get("fits", {})
    data = {
        "nufit": {o: {k: {q: nf[o][k][q] for q in ("best", "lo3", "hi3")}
                      for k in ("dm21", "dm3l", "s12sq", "s13sq")} for o in ("NO", "IO")},
        "nufit_dchi2_io": nf["IO"]["dchi2_vs_global_min"] - nf["NO"]["dchi2_vs_global_min"],
        "profiles": {k: _profile(fits, k, nf) for k in ("cmb_desi", "planck18", "tt_bao")},
        "limits": {
            "DESI DR2 + CMB": cfg.MNU_PUBLISHED["DESI DR2 + CMB (CamSpec)"],
            "Planck 2018 + BAO": cfg.MNU_PUBLISHED["Planck TT,TE,EE+lowE+lensing+BAO"],
        },
        "katrin": cfg.KATRIN_MBETA_LIMIT,
        "pk": {"kh": g["kh"].round(5).tolist(), "masses": g["masses"].tolist(),
               "ratio": [g[f"ratio_{m:g}"].round(5).tolist() for m in g["masses"]],
               "sigma8": [float(g[f"sigma8_{m:g}"]) for m in g["masses"]]},
        "relic": {"kT_eV": r["relic"]["T_nu_eV"], "mean_p_eV": r["relic"]["mean_p_eV"],
                  "x_quantiles": FD_QUANTILES, "n_per_flavour": r["relic"]["n_nu_per_flavour_cm3"],
                  "mass_to_density_eV": r["relic"].get("mass_to_density_eV_dec", 93.03),
                  "omega_m_h2": _omega_m_h2()},
        "neff": r["decoupling"]["runs"][cfg.DECOUPLING_HEADLINE]["neff"],
    }
    blob = json.dumps(_finite(data), separators=(",", ":"), allow_nan=False)
    frag = TEMPLATE.read_text().replace("/*__DATA__*/{}", blob)
    OUT_ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
    OUT_ARTIFACT.write_text(frag)
    # The standalone copy for opening from disk: a full document (standards mode, UTF-8).
    head, body = frag.split("<!--BODY-->", 1)
    doc = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
           "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
           + head + "</head>\n<body>\n" + body + "\n</body>\n</html>\n")
    OUT.write_text(doc)
    print(f"  explorer -> {OUT.relative_to(cfg.PROJECT_DIR)} ({OUT.stat().st_size / 1e3:.0f} kB)")
    return OUT


def main() -> int:
    build()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
