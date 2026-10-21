"""The redshift ladder: one law from a 33-cm laboratory step to a neutron star's surface.

Every rung is a measured frequency shift between two places whose gravitational
potentials differ by dPhi; GR predicts a fractional shift dPhi/c^2 (GM/(r c^2)
for light leaving a body of mass M and radius r). Our own measurements (this
project, and S2's from H1) sit beside published ones.

Published values live in LITERATURE with their sources; see the README's ladder
table. "inferred" rungs are not frequency-shift measurements (neutron-star surface
redshifts come from pulse-profile models of mass and radius) and are drawn hollow.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import config as C
import plotstyle as S

S.use_style()
C2 = C.C_LIGHT**2
G_SURF = 9.80
KMS_C = C.C_LIGHT / 1e3


def _height(h_m: float) -> float:
    return G_SURF * h_m / C2


#: Published measurements, checked against the primary papers (October 2026). x: dPhi/c^2 of the
#: experiment; ratio: measured / GR with its 1-sigma error; kind: clock | mossbauer | line | orbit | inferred.
MSUN_KM = 1.4766          # GM_sun / c^2 in km


def _ns(m, r_km):
    return m * MSUN_KM / r_km


LITERATURE = [
    # Mossbauer gamma rays up the Harvard tower (74 ft = 22.6 m)
    dict(name="Pound & Rebka 1960", group="tower", x=_height(22.6), ratio=1.05, err=0.10, kind="mossbauer",
         ref="PRL 4, 337: (Δν)exp/(Δν)theor = +1.05 ± 0.10"),
    dict(name="Pound & Snider 1965", group="tower", x=_height(22.6), ratio=0.9990, err=0.0076, kind="mossbauer",
         ref="Phys. Rev. 140, B788: 0.9990 ± 0.0076 of 2gh/c² (systematic limit 0.010)"),
    # optical clocks in the laboratory and up a tower
    dict(name="Chou et al. 2010", group="33 cm", x=_height(0.33), ratio=4.1 / 3.6, err=1.6 / 3.6, kind="clock",
         ref="Science 329, 1630: (4.1 ± 1.6)e-17 for 33 cm"),
    dict(name="Zheng et al. 2023", group="1 cm", x=_height(0.01), ratio=1.13, err=0.23, kind="clock",
         ref="Nat. Commun. 14, 4886: α = 0.13 ± 0.23"),
    dict(name="Bothwell et al. 2022", group="1 mm", x=1.09e-19, ratio=9.8 / 10.9, err=2.3 / 10.9, kind="clock",
         ref="Nature 602, 420: −9.8(2.3)e-20 /mm vs −1.09e-19 /mm"),
    dict(name="Takamoto et al. 2020", group="450 m", x=_height(452.6), ratio=1 + 1.4e-5, err=9.1e-5, kind="clock",
         ref="Nat. Photon. 14, 411: α = (1.4 ± 9.1)e-5, Tokyo Skytree"),
    # clocks in space
    dict(name="Gravity Probe A 1976", group="GP-A", x=4.28e-10, ratio=1 + 5e-6, err=1.4e-4, kind="clock",
         ref="Vessot & Levine 1979 final report: ε = (5 ± 140)e-6 (Vessot et al. 1980, PRL 45, 2081)"),
    dict(name="Delva et al. 2018", group="Galileo", x=None, ratio=1 + 0.19e-5, err=2.48e-5, kind="clock",
         ref="PRL 121, 231101: α = (0.19 ± 2.48)e-5"),
    dict(name="Herrmann et al. 2018", group="Galileo", x=None, ratio=1 + 4.5e-5, err=3.1e-5, kind="clock",
         ref="PRL 121, 231102: α_rs = (4.5 ± 3.1)e-5"),
    dict(name="Nunes et al. 2023 (RadioAstron)", group="RadioAstron", x=6.8e-10, ratio=1 + 2.1e-4,
         err=np.hypot(3.3e-4, 0.5e-4), kind="clock", ref="CQG 40, 175005: ε = (2.1 ± 3.3 ± 0.5)e-4"),
    # the Sun
    dict(name="LoPresto et al. 1991", group="Sun", x=None, ratio=627 / 633.1, err=10 / 633.1, kind="line",
         ref="ApJ 376, 757: 627 ± 10 m/s (K I / O I limb lines)"),
    dict(name="González Hernández et al. 2020", group="Sun", x=None, ratio=638 / 633.1, err=6 / 633.1,
         kind="line", ref="A&A 643, A146: 638 ± 6 m/s (3D-corrected); 639 ± 14 (strong lines)"),
    # white dwarfs
    dict(name="40 Eri B (Reid 1996 vs Bond et al. 2017)", group="40 Eri B", x=27.82 / KMS_C,
         ratio=25.8 / 27.82, err=(25.8 / 27.82) * np.hypot(1.4 / 25.8, 0.97 / 27.82), kind="line",
         ref="AJ 111, 2000 (25.8 ± 1.4 km/s); ApJ 848, 16 (predicted 27.82 ± 0.97)"),
    # S2 around Sgr A*
    dict(name="GRAVITY 2018", group="S2", x=None, ratio=0.90, err=np.hypot(0.09, 0.15), kind="orbit",
         ref="A&A 615, L15: f = 0.90 ± 0.09 ± 0.15"),
    dict(name="GRAVITY 2019", group="S2", x=None, ratio=1.04, err=0.05, kind="orbit", ref="A&A 625, L10: f = 1.04 ± 0.05"),
    dict(name="Do et al. 2019", group="S2", x=None, ratio=0.88, err=0.17, kind="orbit",
         ref="Science 365, 664: Υ = 0.88 ± 0.17"),
    # neutron-star surfaces: inferred from pulse-profile models that assume GR (not measured shifts)
    dict(name="PSR J0030+0451 (NICER, Riley et al. 2019)", group="neutron stars", x=_ns(1.34, 12.71), ratio=1.0,
         err=0.0, kind="inferred", ref="ApJL 887, L21: 1.34 M☉, 12.71 km"),
    dict(name="PSR J0740+6620 (NICER, Riley et al. 2021)", group="neutron stars", x=_ns(2.072, 12.39), ratio=1.0,
         err=0.0, kind="inferred", ref="ApJL 918, L27: 2.072 M☉, 12.39 km"),
]


def schwarzschild_z(x):
    """Surface redshift of light leaving a non-rotating body with compactness x = GM/(Rc^2):
    1/sqrt(1 - 2x) - 1, which is x for x << 1."""
    x = np.asarray(x, float)
    return np.expm1(-0.5 * np.log1p(-2 * x))     # stable for tiny x (1/sqrt(1-2x) - 1 underflows)


def galileo_x() -> float:
    """Potential difference between E14's perigee and apogee (the modulated part)."""
    a, e = 27_977.8e3, 0.1625
    return C.GM_EARTH / C2 * (1 / (a * (1 - e)) - 1 / (a * (1 + e)))


def s2_x() -> float:
    """GM/(r c^2) at S2's pericentre (Sgr A* 4.3e6 Msun, r_p = 120 au)."""
    return C.GM_SUN * 4.3e6 / (120 * C.AU) / C2


def table(gal: dict, sun: dict, wd: dict, gps: dict) -> pd.DataFrame:
    rows = []
    xg, xs, xsun = galileo_x(), s2_x(), C.GM_SUN / (C.R_SUN * C2)
    c = gal["combined"]
    rows.append(dict(name="Galileo E14 + E18 (this work)", group="Galileo", x=xg, ratio=1 + c["alpha"],
                     err=c["total"], kind="clock", ours=True, ref="section 3"))
    rows.append(dict(name="GPS clocks (this work)", group="GPS", x=gps["gps_circular"]["gravity"] / 86_400e6,
                     ratio=gps["implied_gain_median"] / gps["predicted_gain_median"],
                     err=gps["published_drift_robust_sd"] / gps["predicted_gain_median"], kind="clock", ours=True,
                     ref="section 1 (consistency check: operators adjust the clocks)"))
    rows.append(dict(name="the Sun (this work)", group="Sun", x=xsun, ratio=sun["strong_mean"] / sun["theory_earth"],
                     err=sun["strong_total_err"] / sun["theory_earth"], kind="line", ours=True, ref="section 7"))
    rows.append(dict(name="Sirius B (Joyce et al. 2018 vs M1's prediction)", group="Sirius B",
                     x=wd["potential_over_c2"], ratio=wd["ratio"], err=wd["ratio_err"], kind="line", ours=True,
                     ref="section 8"))
    h1 = C.PROJECT_DIR.parent / "17-galactic-centre" / "data" / "results.json"
    if h1.exists():
        f = json.loads(h1.read_text())["fits"]["keck_redshift"]["parameters"]["f"]
        rows.append(dict(name="S2 (H1, this atlas)", group="S2", x=xs, ratio=f["median"], err=f["sd"],
                         kind="orbit", ours=True, ref="17-galactic-centre: Keck positions + all velocities"))
    fill = {"Galileo": xg, "Sun": xsun, "S2": xs}
    for r in LITERATURE:
        r = dict(r)
        if r["x"] is None:
            r["x"] = fill[r["group"]]
        rows.append({**r, "ours": False})
    return pd.DataFrame(rows)


KIND_STYLE = {
    "mossbauer": dict(marker="^", label="gamma rays up a tower"),
    "clock": dict(marker="o", label="atomic clocks"),
    "line": dict(marker="s", label="spectral lines of the Sun and white dwarfs"),
    "orbit": dict(marker="D", label="a star passing the Milky Way's black hole"),
    "inferred": dict(marker="o", label="neutron-star surfaces (inferred from models)"),
}

#: group -> (label, offset in points, horizontal alignment) for the left panel
LEFT_LABELS = {
    "1 mm": ("1 mm across an atom cloud (2022)", (10, 0), "left"),
    "1 cm": ("1 cm (2023)", (10, 0), "left"),
    "33 cm": ("33 cm between two clocks (2010)", (10, 0), "left"),
    "tower": ("a 22.6 m tower (1960, 1965)", (10, 0), "left"),
    "450 m": ("Tokyo Skytree, 450 m (2020)", (10, 0), "left"),
    "Galileo": ("Galileo E14/E18 (2018; this work)", (10, 0), "left"),
    "GP-A": ("Gravity Probe A, GPS, RadioAstron", (10, -2), "left"),
    "Sun": ("the Sun", (10, 0), "left"),
    "40 Eri B": ("white dwarfs", (-10, 6), "right"),
    "S2": ("S2 at Sgr A*", (10, 0), "left"),
    "neutron stars": ("neutron stars (NICER)", (-10, 4), "right"),
}

#: group -> (label, y in data units) for the ratio panel
RIGHT_LABELS = {
    "1 mm": ("1 mm", 1.15), "1 cm": ("1 cm", 1.42), "33 cm": ("33 cm", 1.62), "tower": ("tower", 1.19),
    "450 m": ("Skytree", 0.93), "Galileo": ("Galileo", 1.06), "GP-A": ("GP-A, GPS,\nRadioAstron", 0.92),
    "Sun": ("Sun", 1.06), "40 Eri B": ("40 Eri B", 0.84), "Sirius B": ("Sirius B", 1.05),
    "S2": ("S2", 1.17),
}
#: small horizontal spreads (decades) so that rungs at the same depth do not hide each other
JITTER = {"Galileo": [-0.25, 0.0, 0.25], "Sun": [-0.25, 0.0, 0.25], "S2": [0.3, 0.5, 0.7, 0.9],
          "tower": [-0.15, 0.15], "GP-A": [0.0]}
LABEL_SHIFT = {"S2": 0.6}       # decades, where the jitter moved a group


def figure(df: pd.DataFrame) -> list:
    """Left: each rung's shift against the potential (everything on one curve). Right: measured/predicted
    with its 1-sigma error, i.e. how well each rung confirms general relativity."""
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 6.0), gridspec_kw={"width_ratios": [1, 1.2]})
    xs = np.logspace(-20, np.log10(0.42), 300)
    ax.plot(xs, schwarzschild_z(xs), color=S.INK, lw=1.0, zorder=1)
    ax.text(1e-19, 1e-4, "the line: general relativity,\nshift = ΔΦ/c²\n(1/√(1−2ΔΦ/c²) − 1 near\na neutron star)",
            fontsize=8, color=S.INK_SECONDARY, va="top")
    for _, r in df.iterrows():
        st = KIND_STYLE[r["kind"]]
        col = S.BLUE if r["ours"] else S.MUTED
        hollow = r["kind"] == "inferred"
        kw = dict(marker=st["marker"], ms=7, ls="none", color=col, mec=col if hollow else S.SURFACE,
                  mfc="none" if hollow else col, mew=1.3 if hollow else 0.6, zorder=4 if r["ours"] else 3)
        ax.plot([r["x"]], [float(schwarzschild_z(r["x"])) * r["ratio"]], **kw)
    for group, (text, off, ha) in LEFT_LABELS.items():
        g = df[df["group"] == group]
        if g.empty:
            continue
        x = float(g["x"].median())
        ax.annotate(text, (x, float(schwarzschild_z(x))), xytext=off, textcoords="offset points", ha=ha,
                    va="center", fontsize=7.5, color=S.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(3e-20, 3)
    ax.set_ylim(3e-20, 3)
    ax.set_xlabel("depth of the potential difference, ΔΦ/c²")
    ax.set_ylabel("measured fractional frequency shift")
    ax.set_title("One law over twenty decades of gravity")

    bx.axhline(1, color=S.INK, lw=1.1)
    for group, g in df[df["kind"] != "inferred"].groupby("group", sort=False):
        jit = JITTER.get(group, [0.0] * len(g))
        if len(jit) < len(g):
            jit = list(np.linspace(-0.2, 0.2, len(g)))
        for (_, r), dj in zip(g.iterrows(), jit):
            st = KIND_STYLE[r["kind"]]
            col = S.BLUE if r["ours"] else S.MUTED
            bx.errorbar([r["x"] * 10**dj], [r["ratio"]], yerr=[r["err"]], color=col, marker=st["marker"], ms=6.5,
                        capsize=2.5, lw=1.3, ls="none", mec=S.SURFACE, mew=0.6, zorder=4 if r["ours"] else 3)
    for group, (text, y) in RIGHT_LABELS.items():
        g = df[df["group"] == group]
        if g.empty:
            continue
        xg = float(g["x"].median()) * 10 ** LABEL_SHIFT.get(group, 0.0)
        bx.text(xg, y, text, ha="center", va="center", fontsize=7.3, color=S.INK_SECONDARY)
    bx.set_xscale("log")
    bx.set_yscale("log")
    bx.set_xlim(3e-20, 3e-1)
    bx.set_ylim(0.6, 1.75)
    ticks = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6]
    bx.set_yticks(ticks)
    bx.set_yticklabels([f"{t:.1f}" for t in ticks])
    bx.minorticks_off()
    bx.set_xlabel("depth of the potential difference, ΔΦ/c²")
    bx.set_ylabel("measured ÷ general relativity's prediction")
    bx.set_title("How well each rung agrees")
    handles = [plt.Line2D([], [], color=S.MUTED, marker=st["marker"], ls="none", ms=7,
                          mfc="none" if k == "inferred" else S.MUTED, label=st["label"])
               for k, st in KIND_STYLE.items() if (df["kind"] == k).any()]
    handles.append(plt.Line2D([], [], color=S.BLUE, marker="o", ls="none", ms=7,
                              label="measured or predicted in this atlas (N3, M1, H1)"))
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.06))
    S.caption(fig, "ΔΦ/c² is gh/c² for a height h and GM/(rc²) for light leaving a star (for Galileo, the perigee-apogee "
              "difference). Neutron-star points use NICER's masses and radii, which come from models that assume GR: "
              "they show where the ladder ends, not a test. Many bars are smaller than their markers. Sources in the "
              "README.", y=-0.1)
    return S.save(fig, C.FIGURE_DIR / "redshift_ladder")
