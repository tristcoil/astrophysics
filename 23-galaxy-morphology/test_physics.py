"""Self-checks for the physics, the vote handling and the estimators.

Run with ``python test_physics.py``. No test framework and no network needed:
each function is checked against a case with a known answer, not against the
data. A slip in a K-correction, a decision-tree product or the density
estimator would give plausible-looking but wrong fractions, so each gets an
independent check.
"""

from __future__ import annotations

import importlib.util
import io
import sys
import zipfile
from pathlib import Path

import numpy as np
from astropy import units as u

import analysis as an
import environment as env
import fetch_data as fd
import morphology as mo
import physics as ph
from config import COMPLETENESS_R, DATA_DIR, PROJECT_DIR

FAILURES: list[str] = []


def check(label: str, got: float, expected: float, tol: float, unit: str = "",
          *, absolute: bool = False) -> None:
    """Assert ``got`` is within ``tol`` of ``expected`` (fractional unless ``absolute``)."""
    err = abs(got - expected) if (absolute or expected == 0) else abs(got - expected) / abs(expected)
    ok = err <= tol
    if not ok:
        FAILURES.append(label)
    kind = "abs" if (absolute or expected == 0) else "rel"
    print(f"  [{'ok  ' if ok else 'FAIL'}] {label:<56} {got:>12.6g} vs {expected:>12.6g} "
          f"{unit:<5} ({kind} err {err:.2g}, tol {tol:g})")


def test_distances() -> None:
    print("Distances")
    for z in (0.02, 0.085, 0.3):
        exact = ph.COSMO.comoving_distance(z).to_value(u.Mpc)
        check(f"interpolated D_C(z={z}) vs astropy", float(ph.comoving_distance(z)), exact, 1e-7, "Mpc")
    check("distance modulus z=0.05 vs astropy", float(ph.distance_modulus(0.05)),
          ph.COSMO.distmod(0.05).value, 1e-6, "mag")
    exact = (ph.COSMO.angular_diameter_distance(0.05) * (1 * u.arcsec).to(u.rad).value).to_value(u.kpc)
    check("kpc per arcsec at z=0.05 vs astropy", float(ph.kpc_per_arcsec(0.05)), exact, 1e-6, "kpc")
    check("magnitude shift to H0 = 70 (5 log10(0.70/0.677))", ph.h_magnitude_shift(0.70),
          5 * np.log10(0.70 / ph.H_LITTLE), 1e-12, "mag")
    check("LOS depth of 1000 km/s at z=0.001", float(ph.los_window(1000.0, 0.001)), 10.0, 3e-3,
          "Mpc/h")


def test_k_corrections() -> None:
    print("\nK-corrections (Chilingarian et al. 2010, 2012)")
    # Values from running the authors' calc_kcor.py (kcor.sai.msu.ru) itself.
    reference = [("u", 0.05, "ur", 2.6, 0.19133420305000004), ("r", 0.05, "ur", 2.6, 0.05548049980375001),
                 ("u", 0.08, "ur", 1.6, 0.11064947691520041), ("r", 0.08, "ur", 1.6, 0.03334781986816004),
                 ("g", 0.06, "gi", 1.2, 0.13032334916256), ("i", 0.06, "gi", 1.2, 0.03726332023391996),
                 ("g", 0.03, "gi", 0.7, 0.03692265732333), ("i", 0.03, "gi", 0.7, -0.0037589987286900176),
                 ("g", 0.1, "gr", 0.9, 0.23910890019999978), ("r", 0.1, "gr", 0.9, 0.08405887050000001),
                 ("g", 0.2, "gr", 1.1, 0.5209713975999992)]
    for band, z, col, val, expected in reference:
        check(f"K_{band}(z={z}, {col[0]}-{col[1]}={val}) vs authors' code",
              float(ph.k_correction(band, z, val, col)), expected, 1e-12, "mag", absolute=True)
    # If the authors' file is cached, compare on a whole grid too.
    path = DATA_DIR / "calc_kcor.py"
    if path.exists():
        spec = importlib.util.spec_from_file_location("calc_kcor_ref", path)
        ref = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(ref)
        worst = 0.0
        for band, col in (("u", "ur"), ("r", "ur"), ("g", "gi"), ("i", "gi"), ("g", "gr"), ("r", "gr")):
            for z in np.linspace(0.0, 0.3, 7):
                for cv in np.linspace(0.2, 3.0, 8):
                    worst = max(worst, abs(float(ph.k_correction(band, z, cv, col))
                                           - ref.calc_kcor(band, z, f"{col[0]} - {col[1]}", cv)))
        check("all six fits on a 7 x 8 grid vs the authors' file", worst, 0.0, 1e-11, "mag",
              absolute=True)
    # Rest-frame colour = observed - K_blue + K_red.
    obs = 0.9
    expect = obs - ph.k_correction("g", 0.1, obs, "gr") + ph.k_correction("r", 0.1, obs, "gr")
    check("rest-frame g-r definition", float(ph.rest_frame_colour("g", "r", obs, 0.1)), float(expect),
          1e-14, "mag", absolute=True)


def test_same_as_d2() -> None:
    """The physics must be D2's, so luminosities and densities mean the same thing."""
    print("\nSame formulas as project D2")
    d2 = sorted(PROJECT_DIR.parent.glob("*-cosmic-web/cosmo_physics.py"))
    if not d2:
        print("  (project D2 not found next to this one; skipped)")
        return
    spec = importlib.util.spec_from_file_location("d2_cosmo_physics", d2[0])
    cp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cp)
    rng = np.random.default_rng(1)
    z = rng.uniform(0.01, 0.25, 1000)
    gr = rng.uniform(0.2, 1.1, 1000)
    m = rng.uniform(14, 17.8, 1000)
    worst_k = np.max(np.abs(ph.k_correction("r", z, gr, "gr") - cp.k_correction("r", z, gr)))
    worst_c = np.max(np.abs(ph.rest_frame_colour("g", "r", gr, z) - cp.rest_frame_color(gr, z)))
    worst_m = np.max(np.abs(ph.absolute_magnitude(m, z, 0.1) - cp.absolute_magnitude(m, z, 0.1)))
    check("K_r against D2's code, 1000 random galaxies", float(worst_k), 0.0, 1e-13, "mag", absolute=True)
    check("rest-frame g-r against D2's code", float(worst_c), 0.0, 1e-13, "mag", absolute=True)
    check("absolute magnitude against D2's code", float(worst_m), 0.0, 1e-12, "mag", absolute=True)


def _toy_catalog(rows: list[dict]) -> dict:
    keys = set()
    for r in rows:
        keys |= set(r)
    names = list(mo.ARM_KEYS) + ["smooth", "features", "artifact", "edgeon", "notedgeon", "spiral",
                                 "nospiral", "bar", "nobar", "arms_cant_tell"]
    c = {}
    for k in names:
        c["p_" + k] = np.array([r.get("p_" + k, 0.0) for r in rows], np.float32)
        c["wf_" + k] = np.array([r.get("wf_" + k, r.get("p_" + k, 0.0)) for r in rows], np.float32)
        c["n_" + k] = np.array([r.get("n_" + k, 0) for r in rows], np.int32)
    c["expAB_r"] = np.array([r.get("ba", 0.9) for r in rows], float)
    return c


def test_votes() -> None:
    print("\nReading the votes")
    c = _toy_catalog([
        # A face-on barred spiral: 90% features, 80% of those not edge-on, ...
        {"p_features": 0.9, "p_smooth": 0.1, "p_notedgeon": 0.8, "p_edgeon": 0.2, "p_spiral": 0.75,
         "p_bar": 0.6, "n_bar": 12, "n_nobar": 8, "ba": 0.8},
        # The same votes but seen nearly edge-on (b/a 0.3): bars cannot be judged.
        {"p_features": 0.9, "p_smooth": 0.1, "p_notedgeon": 0.8, "p_spiral": 0.75, "p_bar": 0.6,
         "n_bar": 12, "n_nobar": 8, "ba": 0.3},
        # A clean elliptical; nobody was asked about bars.
        {"p_smooth": 0.95, "p_features": 0.05},
    ])
    pp = mo.path_probabilities(c)
    check("face-on spiral = p_features x p_notedgeon x p_spiral", float(pp["spiral"][0]),
          0.9 * 0.8 * 0.75, 1e-6)
    cls = mo.classes(c)
    check("spiral class (product 0.54 > 0.5)", float(cls["spiral"][0]), 1.0, 0, absolute=True)
    check("not a clean spiral (0.54 < 0.8)", float(cls["spiral_clean"][0]), 0.0, 0, absolute=True)
    check("barred: disk, b/a 0.8, p_bar 0.6 >= 0.5", float(cls["barred"][0]), 1.0, 0, absolute=True)
    check("not judged for bars when b/a = 0.3", float(cls["disk"][1]), 0.0, 0, absolute=True)
    check("clean elliptical", float(cls["smooth"][2]), 1.0, 0, absolute=True)
    check("elliptical is not a disk", float(cls["disk"][2]), 0.0, 0, absolute=True)

    # Arm number: the modal debiased answer, if >= 5 votes name a number.
    c = _toy_catalog([
        {"p_arms_2": 0.6, "p_arms_3": 0.3, "n_arms_2": 6, "n_arms_3": 3},
        {"p_arms_5plus": 0.5, "p_arms_2": 0.4, "n_arms_5plus": 2, "n_arms_2": 2, "n_arms_cant_tell": 9},
    ])
    arm, pmax = mo.arm_number(c)
    check("modal arm number (two arms)", float(arm[0]), 1.0, 0, absolute=True)
    check("too few definite arm votes (4 < 5) -> unclassified", float(arm[1]), -1.0, 0, absolute=True)

    v, e = mo.weighted_fraction([1, 0, 1, 1], [1, 1, 1, 1])
    check("unweighted mean of 1,0,1,1", v, 0.75, 1e-12)
    v, _ = mo.weighted_fraction([1, 0], [3, 1])
    check("weighted mean (weights 3:1)", v, 0.75, 1e-12)
    f, e, n = mo.binomial_fraction([True] * 30 + [False] * 70)
    check("binomial fraction 30/100", f, 0.30, 1e-12)
    check("binomial error sqrt(0.3*0.7/100)", e, np.sqrt(0.21 / 100), 1e-12)


def test_decomposition() -> None:
    """The colour-density split must add up exactly, whatever the cells."""
    print("\nEnvironment decomposition")
    rng = np.random.default_rng(3)
    worst_a, worst_b = 0.0, 0.0
    for _ in range(200):
        cells = []
        for _ in range(2):
            n = 400
            ws = rng.uniform(0, 1, n)
            red = rng.uniform(0, 1, n) < rng.uniform(0.2, 0.9)
            cells.append(an._cell_stats(red, ws, ws > 0.8, ws < 0.2))
        dd = an._decompose(*cells)
        worst_a = max(worst_a, abs(dd["red_from_morphology_mix"] + dd["red_at_fixed_morphology"]
                                   - dd["delta_red"]))
        worst_b = max(worst_b, abs(dd["smooth_from_colour_mix"] + dd["smooth_at_fixed_colour"]
                                   - dd["delta_smooth"]))
    check("red: morphology-mix + fixed-morphology = total (200 cases)", worst_a, 0.0, 1e-12,
          absolute=True)
    check("smooth: colour-mix + fixed-colour = total (200 cases)", worst_b, 0.0, 1e-12, absolute=True)
    # A world where colour depends only on morphology: all of the change is morphology mix.
    n = 20000
    lo_ws = (rng.uniform(0, 1, n) < 0.3).astype(float)
    hi_ws = (rng.uniform(0, 1, n) < 0.6).astype(float)
    lo_red = np.where(lo_ws > 0.5, rng.uniform(0, 1, n) < 0.9, rng.uniform(0, 1, n) < 0.3)
    hi_red = np.where(hi_ws > 0.5, rng.uniform(0, 1, n) < 0.9, rng.uniform(0, 1, n) < 0.3)
    dd = an._decompose(an._cell_stats(lo_red, lo_ws, lo_ws > 0.5, lo_ws < 0.5),
                       an._cell_stats(hi_red, hi_ws, hi_ws > 0.5, hi_ws < 0.5))
    check("pure morphology-density world: fixed-morphology term ~ 0",
          dd["red_at_fixed_morphology"], 0.0, 0.015, absolute=True)
    check("  ...and the mix term = 0.3 x (0.9 - 0.3) = 0.18", dd["red_from_morphology_mix"], 0.18, 0.08)


def test_density() -> None:
    print("\nLocal density (D2's Sigma_5)")
    # A target at 100 h^-1 Mpc along x; tracers on a circle across the line of
    # sight at projected radius 1, 2, ..., 8 h^-1 Mpc, plus decoys far behind.
    d0 = 100.0
    z0 = float(np.interp(d0, ph.comoving_distance_h(np.linspace(0, 0.1, 2001)), np.linspace(0, 0.1, 2001)))
    target = np.array([[d0, 0.0, 0.0]])
    ang = np.linspace(0, 2 * np.pi, 9)[:-1]
    radii = np.arange(1.0, 9.0)
    tracers = np.column_stack([np.full(8, d0), radii * np.cos(ang), radii * np.sin(ang)])
    window = float(ph.los_window(1000.0, z0))
    decoys = np.array([[d0 + 3 * window, 0.3, 0.0], [d0 - 3 * window, 0.0, 0.4]])  # outside +/-1000 km/s
    own = np.array([[d0, 0.0, 0.0]])  # the target itself, also a tracer
    allt = np.vstack([tracers, decoys, own])
    ids = np.arange(len(allt))
    res = env.local_density(target, np.array([z0]), np.array([ids[-1]]), allt, ids, k=10)
    # The pair line of sight is the midpoint direction, so r_p is the projected
    # radius to within ~(r/d)^2 -- a few parts in 10^4 here.
    check("5th neighbour: the tracer at 5 h^-1 Mpc", float(res["d_n"][0]), 5.0, 2e-3, "Mpc/h")
    check("Sigma_5 = 5 / (pi d5^2)", float(res["sigma"][0]), 5 / (np.pi * 25), 4e-3, "h2/Mpc2")


def test_mask() -> None:
    print("\nFootprint mask")
    # A survey covering a 20-degree-radius cap around the north pole.
    rng = np.random.default_rng(5)
    n = 400_000
    ra = rng.uniform(0, 360, n)
    dec = np.degrees(np.arcsin(rng.uniform(np.sin(np.radians(70)), 1, n)))
    m = env.build_mask(ra, dec)
    edge = float(m.edge_distance_deg(np.array([0.0]), np.array([90.0]))[0])
    check("pole of a 20-deg cap is ~20 deg from the edge", edge, 20.0, 0.06, "deg")
    check("a point at dec 60 is outside", float(m.contains([10.0], [60.0], "surveyed")[0]), 0.0, 0,
          absolute=True)


def test_completeness() -> None:
    print("\nVolume limits")
    kz = (np.array([0.0, 0.2]), np.array([0.0, 0.2]))  # K = z for the reddest galaxies (toy)
    z = 0.085
    m = float(an.faint_limit(z, kz))
    r = m + float(ph.distance_modulus(z)) + 0.085
    check("M_lim(z) puts the reddest galaxy exactly on the flux limit", r, COMPLETENESS_R, 1e-12, "mag")


def test_colour_divider() -> None:
    print("\nColour divider (D2's method) on synthetic galaxies")
    rng = np.random.default_rng(11)
    M = rng.uniform(-23, -19, 40000)
    blue = rng.uniform(0, 1, 40000) < 0.5
    # Two populations, each with a colour that drifts with luminosity.
    gr0 = np.where(blue, 0.45 - 0.03 * (M + 21) + rng.normal(0, 0.06, 40000),
                   0.78 - 0.02 * (M + 21) + rng.normal(0, 0.04, 40000))
    div = an.colour_divider(M, gr0, np.ones_like(M, bool))
    # Equal-probability point at M = -21 for equal weights: where the two
    # weighted Gaussians cross, computed directly.
    x = np.linspace(0.45, 0.78, 100001)
    g1 = np.exp(-0.5 * ((x - 0.45) / 0.06) ** 2) / 0.06
    g2 = np.exp(-0.5 * ((x - 0.78) / 0.04) ** 2) / 0.04
    cross = x[np.argmin(np.abs(g1 - g2))]
    check("divider at M = -21 on the analytic crossing", div["slope"] * -21 + div["intercept"],
          float(cross), 0.01, "mag", absolute=True)


def test_zip_reading() -> None:
    """The byte-range zip reader, on a zip built in memory."""
    print("\nReading single images out of a zip")
    buf = io.BytesIO()
    payloads = {"images/7.jpg": b"\xff\xd8" + bytes(range(256)) * 40,
                "__MACOSX/images/._7.jpg": b"resource fork",
                "images/12.jpg": b"\xff\xd8" + b"galaxy" * 500}
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for name, data in payloads.items():
            z.writestr(name, data)
    raw = buf.getvalue()
    i = raw.rfind(b"PK\x05\x06")
    cd_size = int.from_bytes(raw[i + 12:i + 16], "little")
    cd_off = int.from_bytes(raw[i + 16:i + 20], "little")
    index = fd.parse_central_directory(raw[cd_off:cd_off + cd_size])
    check("two images indexed, the resource fork skipped", float(len(index["asset_id"])), 2.0, 0,
          absolute=True)
    ok = True
    for a, off, cs, meth in zip(index["asset_id"], index["offset"], index["csize"], index["method"]):
        blob = raw[off:off + 30 + 512 + cs]
        ok &= fd.extract_member(blob, int(cs), int(meth)) == payloads[f"images/{a}.jpg"]
    check("each image's bytes recovered exactly from its range", float(ok), 1.0, 0, absolute=True)


def test_loss() -> None:
    try:
        import torch
    except ImportError:
        print("\n(PyTorch not installed: neural network loss not checked)")
        return
    import cnn

    print("\nNeural network loss")
    counts = torch.tensor([[30.0, 10.0, 0.0, 4.0, 6.0, 3.0, 7.0, 5.0, 5.0]])
    # Logits giving exactly the observed fractions: the loss is the votes' entropy.
    probs = [counts[0, :3] / 40, counts[0, 3:5] / 10, counts[0, 5:7] / 10, counts[0, 7:9] / 10]
    logits = torch.log(torch.cat(probs).clamp(min=1e-12)).unsqueeze(0).requires_grad_(True)
    loss = cnn.multinomial_nll(logits, counts)
    ent = sum(-(cn * torch.log(p.clamp(min=1e-12))).sum() for cn, p in
              zip((counts[0, :3], counts[0, 3:5], counts[0, 5:7], counts[0, 7:9]), probs))
    check("loss at the observed fractions = entropy per vote", float(loss), float(ent / counts.sum()),
          1e-6)
    loss.backward()
    check("and it is a minimum (gradient ~ 0)", float(logits.grad.abs().max()), 0.0, 1e-6,
          absolute=True)


def main() -> int:
    test_distances()
    test_k_corrections()
    test_same_as_d2()
    test_votes()
    test_decomposition()
    test_density()
    test_mask()
    test_completeness()
    test_colour_divider()
    test_zip_reading()
    test_loss()
    print()
    if FAILURES:
        print(f"{len(FAILURES)} check(s) FAILED: {FAILURES}")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
