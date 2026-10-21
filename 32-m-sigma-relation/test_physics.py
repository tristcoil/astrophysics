"""Check the physics, the fitting code and the table readers against answers known independently.

    python test_physics.py

Each check prints PASS or FAIL with the numbers compared; the script exits
non-zero if any check fails. Checks that need the downloaded tables are
skipped (with a note) when data/raw/ is empty.
"""

from __future__ import annotations

import sys

import numpy as np

import catalogs as K
import config as C
import fitting as F
import physics as P

RESULTS: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, bool(ok), detail))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))


def close(a, b, rtol=1e-6, atol=0.0):
    return bool(np.all(np.isclose(a, b, rtol=rtol, atol=atol)))


# --------------------------------------------------------------------------

def test_scales():
    print("Scales")
    check("G = 4.301e-3 pc (km/s)^2 / Msun (G Msun / pc)", close(C.G_PC, 1.32712440018e20 / 3.0856775814913673e16 / 1e6,
                                                                   1e-5), f"{C.G_PC:.6e}")
    check("sphere of influence of 1e8 Msun at 200 km/s is 10.75 pc", close(P.r_influence_pc(1e8, 200.0), 10.752, 1e-3),
          f"{P.r_influence_pc(1e8, 200.0):.3f} pc")
    rs = P.schwarzschild_radius_pc(1.0) * 3.0856775814913673e13  # km
    check("Schwarzschild radius of the Sun is 2.953 km", close(rs, 2.953, 1e-3), f"{rs:.4f} km")
    check("1 pc at 1 Mpc subtends 0.206 arcsec", close(P.angular_size_arcsec(1.0, 1.0), 0.2062648, 1e-6))
    m = P.min_measurable_mass(150.0, 18.0, 0.1)
    check("the smallest measurable mass has a sphere of influence of exactly the resolution",
          close(P.angular_size_arcsec(P.r_influence_pc(m, 150.0), 18.0), 0.1, 1e-10))
    # Sgr A*: 4e6 Msun at sigma ~ 105 km/s has r_infl ~ 1.5-2 pc (the classic number).
    r = P.r_influence_pc(4.0e6, 105.0)
    check("Sgr A*'s sphere of influence is ~1.6 pc", 1.4 < r < 1.8, f"{r:.2f} pc")


def test_feedback():
    print("Feedback and energy")
    check("electron scattering: kappa = sigma_T / m_p = 0.398 cm^2/g for pure hydrogen",
          close(P.kappa_es(1.0), 0.3977, 1e-3), f"{P.kappa_es(1.0):.4f}")
    check("kappa = 0.338 cm^2/g for X = 0.70", close(P.kappa_es(0.7), 0.338, 2e-3), f"{P.kappa_es(0.7):.4f}")
    led = P.eddington_luminosity(1.0)
    check("Eddington luminosity of 1 Msun is 1.26e38 erg/s", close(led, 1.257e38, 2e-3), f"{led:.4e}")
    s = np.array([100.0, 200.0, 400.0])
    m = P.king_msigma(s)
    check("King's M_sigma scales exactly as sigma^4", close(np.diff(np.log10(m)) / np.diff(np.log10(s)), 4.0, 1e-12))
    # The same number by hand, in cgs: f kappa sigma^4 / (pi G^2) / Msun.
    hand = 0.16 * 0.3977 * (2e7) ** 4 / (np.pi * (6.6743e-8) ** 2) / 1.98841e33
    check("King's M_sigma at 200 km/s (f = 0.16, pure H) = 3.66e8 Msun by hand",
          close(P.king_msigma(200.0, f_gas=0.16, kappa=P.kappa_es(1.0)), hand, 2e-3),
          f"{P.king_msigma(200.0, f_gas=0.16, kappa=P.kappa_es(1.0)):.3e} vs {hand:.3e}")
    er = P.energy_ratio(1e8, 2e10, 200.0)
    check("energy ratio 0.1 M c^2 / (M_b sigma^2) = 1124 for M/M_b = 0.5%, sigma = 200",
          close(er, 0.1 * 0.005 * (299792.458 / 200.0) ** 2, 1e-9), f"{er:.0f}")


def test_fitting():
    print("Line fitting (a simulation with known answers)")
    rng = np.random.default_rng(42)
    x = rng.uniform(-0.6, 0.35, 400)
    y = 8.3 + 5.0 * x
    f = F.fit_fitexy(x, np.full_like(x, 1e-6), y, np.full_like(x, 0.01))
    check("noiseless line: least squares recovers a, b exactly, eps = 0",
          close([f.a, f.b], [8.3, 5.0], 1e-6) and f.eps < 1e-6, f"a={f.a:.6f} b={f.b:.6f} eps={f.eps:.2g}")
    # A realistic case: scatter 0.4 dex, x errors 0.03 dex, y errors 0.15 dex.
    xt = rng.normal(-0.1, 0.2, 300)
    yt = 8.3 + 5.0 * xt + 0.4 * rng.standard_normal(300)
    xo = xt + 0.03 * rng.standard_normal(300)
    yo = yt + 0.15 * rng.standard_normal(300)
    ex, ey = np.full(300, 0.03), np.full(300, 0.15)
    fs = F.fit_bayes(xo, ex, yo, ey, steps=3000, burn=1000)
    ok = abs(fs.b - 5.0) < 3 * fs.eb and abs(fs.a - 8.3) < 3 * fs.ea and abs(fs.eps - 0.4) < 3 * fs.eeps
    check("Bayesian fit recovers slope 5, zero point 8.3, scatter 0.4 (within 3 sigma)", ok,
          f"b={fs.b:.2f}±{fs.eb:.2f} a={fs.a:.2f}±{fs.ea:.2f} eps={fs.eps:.2f}±{fs.eeps:.2f}")
    # Censoring: hide every true value below a threshold behind an upper limit.
    lim = 7.6
    up = yo < lim
    yc = np.where(up, lim, yo)
    fc = F.fit_bayes(xo, ex, yc, np.where(up, 0.0, ey), up, steps=3000, burn=1000)
    fd = F.fit_fitexy(xo[~up], ex[~up], yo[~up], ey[~up])
    check(f"with {up.sum()} values turned into upper limits, the censored fit still finds slope 5",
          abs(fc.b - 5.0) < 3 * fc.eb, f"b={fc.b:.2f}±{fc.eb:.2f}")
    check("dropping those limits instead biases the slope low (the selection effect we simulate)",
          fd.b < fc.b - 0.3, f"detections only b={fd.b:.2f}")


def test_estimator_bias():
    print("Estimator bias over 20 simulated samples (true slope 5, x errors 0.03 dex, x spread 0.2 dex)")
    from scipy import optimize
    out = {"structural": [], "conditional": [], "fitexy": []}
    for seed in range(20):
        rng = np.random.default_rng(1000 + seed)
        xt = rng.normal(-0.1, 0.2, 300)
        xo = xt + 0.03 * rng.standard_normal(300)
        yo = 8.3 + 5.0 * xt + 0.4 * rng.standard_normal(300) + 0.15 * rng.standard_normal(300)
        ex, ey, w, up = np.full(300, 0.03), np.full(300, 0.15), np.ones(300), np.zeros(300, bool)
        span = (xo.min() - 1, xo.max() + 1)
        for model, ng in (("structural", 1), ("conditional", 0)):
            p0 = [8.3, 5.0, 0.4] + ([np.mean(xo), np.log(np.std(xo))] if ng else [])
            nll = lambda t: -float(F.log_posterior(np.array(t), xo, ex, yo, ey, up, w, model, ng, span)[0])  # noqa: E731
            out[model].append(optimize.minimize(nll, p0, method="Nelder-Mead",
                                                options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 20000}).x[1])
        out["fitexy"].append(F.fit_fitexy(xo, ex, yo, ey).b)
    m = {k: (np.mean(v), np.std(v, ddof=1) / np.sqrt(len(v))) for k, v in out.items()}
    check("the structural likelihood (our default) is unbiased", abs(m["structural"][0] - 5.0) < 3 * m["structural"][1],
          f"mean slope {m['structural'][0]:.3f} ± {m['structural'][1]:.3f}")
    check("so is least squares with tuned scatter (FITEXY, as Kormendy & Ho and McConnell & Ma used)",
          abs(m["fitexy"][0] - 5.0) < 3 * m["fitexy"][1], f"{m['fitexy'][0]:.3f} ± {m['fitexy'][1]:.3f}")
    dil = 5.0 * 0.2**2 / (0.2**2 + 0.03**2)
    check("taking measured sigma as given dilutes the slope by tau^2/(tau^2 + ex^2), as expected",
          abs(m["conditional"][0] - dil) < 3 * m["conditional"][1],
          f"{m['conditional'][0]:.3f} ± {m['conditional'][1]:.3f}, expected {dil:.3f}")


def test_names():
    print("Galaxy names")
    cases = {"N4486 (M87)": "NGC 4486", "M{\\ts}87": "NGC 4486", "NGC0221": "NGC 221", "M 31": "NGC 224",
             "A1836-BCG": "A1836 BCG", "A1836 BCG": "A1836 BCG", "Galaxy": "Milky Way", "MW": "Milky Way",
             "Milky Way $^a$": "Milky Way", "U3789 $^c$": "UGC 3789", "UGC03789": "UGC 3789", "Mrk0335": "Mrk 335",
             "NGC 4486B": "NGC 4486B", "N4486A $^x$": "NGC 4486A", "IC 1459": "IC 1459", "Arp151": "Mrk 40",
             "N5128 (Cen A) $^{h,w}$": "NGC 5128", "Cygnus A": "Cygnus A"}
    bad = {k: K.canonical(k) for k, v in cases.items() if K.canonical(k) != v}
    check(f"{len(cases)} spellings map to one name each", not bad, str(bad) if bad else "")


def test_tables():
    print("Reading the published tables")
    if not (C.RAW_DIR / "kh13" / "annrev4K-S.ps").exists():
        print("  (skipped: run fetch_data.py first)")
        return
    kh = K.load_kh13()
    fit = kh[kh["kh13_fit"]]
    check("Kormendy & Ho: 45 ellipticals, 21 classical-bulge rows, 22 pseudobulges",
          ((kh["table"] == 2).sum(), (kh["bulge"] == "classical").sum(), (kh["bulge"] == "pseudo").sum()) == (45, 21, 22))
    check("Kormendy & Ho's fit sample: 45 galaxies, 29 below 270 km/s (their Section 6.7)",
          (len(fit), int((fit["sigma"] < 270).sum())) == (45, 29))
    check("their mean log M error over the fitted galaxies is 0.117 dex (their Section 6.6)",
          abs(fit["elogM"].mean() - 0.117) < 0.0005, f"{fit['elogM'].mean():.4f}")
    check("colour flags read off the rendered table: NGC 1277 and NGC 4342 masses cyan, NGC 4596 row cyan",
          set(kh.loc[(kh["table"] == 3) & (kh["mass_color"] == "cyan"), "name"]) ==
          {"NGC 1277", "NGC 4342", "NGC 4596", "NGC 3998"})
    mm = K.load_mm13()
    u = mm.drop_duplicates("name")
    check("McConnell & Ma: 72 galaxies, 53 early and 19 late types",
          (len(u), int(u["early_type"].sum()), int((~u["early_type"]).sum())) == (72, 53, 19))
    v = K.load_vdb16()
    vs = v[v["selected"] == 1]
    check("van den Bosch: 230 fitted galaxies, 49 upper limits", (len(vs), int(vs["upper"].sum())) == (230, 49))
    m33 = v[v["name"] == "NGC 598"].iloc[0]
    check("M33's upper limit decodes to < 1,400 Msun (Gebhardt et al. 2001: < 1,500)",
          1300 < 10 ** m33["logM_limit"] < 1500, f"{10 ** m33['logM_limit']:.0f}")
    ag = K.load_agn("virial_products")
    ag43 = K.load_agn("grier13")
    d = (ag43.set_index("object")["logM"] - ag.set_index("object")["logM"]).dropna()
    check("AGN database: masses with <f> = 4.3 are exactly 4.3 x the virial products",
          close(d.to_numpy(), np.log10(4.3), 1e-9), f"{len(d)} objects")
    print("Reproducing published fits (least squares, fast)")
    f = F.fit_msigma(fit, method="fitexy", equal_logm_err=C.KH13_EQUAL_LOGM_ERR)
    check("Kormendy & Ho eq. 3: slope 4.377 ± 0.290, zero point -0.510", abs(f.b - 4.377) < 0.1 * 0.29
          and abs(f.a - 9 + 0.510) < 0.1 * 0.049, f"{f.b:.3f}, {f.a - 9:.3f}")
    f = F.fit_msigma(mm, method="fitexy")
    check("McConnell & Ma: slope 5.64 ± 0.32, zero point 8.32, scatter 0.38",
          abs(f.b - 5.64) < 0.1 * 0.32 and abs(f.a - 8.32) < 0.15 * 0.05 and abs(f.eps - 0.38) < 0.01,
          f"{f.b:.3f}, {f.a:.3f}, {f.eps:.3f}")


def main() -> int:
    for t in (test_scales, test_feedback, test_fitting, test_estimator_bias, test_names, test_tables):
        t()
    n_fail = sum(not ok for _, ok, _ in RESULTS)
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks passed")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
