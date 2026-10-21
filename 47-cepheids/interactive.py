"""The interactive 3D map: every Cepheid placed by its own period and brightness.

A self-contained HTML page (plotly.js and the data are embedded, so it works
offline and can be shared as one file). Each Cepheid's distance comes from the
Leavitt law calibrated in this project: Magellanic Cepheids from OGLE's W_I
with our LMC distance, Milky Way Cepheids from Gaia's W_G with the same
calibration (and their spectroscopic metallicity where known). It follows the
viewer's light/dark setting, with colours re-stepped for each surface.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import plotstyle as ps

THEMES = {
    "light": dict(surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED,
                  lmc=ps.BLUE, smc=ps.ORANGE, mw=ps.AQUA, menu_bg="#f0efec", menu_border="#c3c2b7"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781",
                 lmc="#3987e5", smc="#d95926", mw="#199e70", menu_bg="#2c2c2a", menu_border="#383835"),
}


def _galactic_xyz(ra, dec, d_kpc):
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg, frame="icrs").galactic
    l, b = c.l.radian, c.b.radian
    return d_kpc * np.cos(b) * np.cos(l), d_kpc * np.cos(b) * np.sin(l), d_kpc * np.sin(b)


def star_distances(res) -> dict:
    """Distance (kpc) of every Cepheid from the calibrated Leavitt law."""
    import analysis
    s = analysis.summarize(res)
    mu_lmc = s["lmc_distance"]["mu"]
    out = {}
    for gal in ("lmc", "smc"):
        t = res["tabs"][gal]
        mu = np.full(len(t["id"]), np.nan)
        for mode in ("F", "1O"):
            f = res["laws"]["lmc"][(mode, "W_I")]
            m = (t["mode"] == mode) & np.isfinite(t["W_I"])
            mu[m] = mu_lmc + t["W_I"][m] - (f["a"] + f["b"] * (t["logP"][m] - 1))
        # only the stars kept by the 3-sigma-clipped fits of the geometry section
        # (blends, which look too bright and so too near, are dropped; so are SMC
        # fundamental-mode stars shorter than 2.5 d, where the law bends)
        ok = np.isfinite(mu) & np.isfinite(res["geometry"][gal]["dmu"])
        d = 10 ** (0.2 * mu + 1) / 1e3
        out[gal] = {"ok": ok, "d": d, "ra": t["ra"], "dec": t["dec"], "id": t["id"], "P": t["P"], "mode": t["mode"]}
    # Milky Way: Gaia W_G with the fitted R, the LMC's law shifted by mu_LMC, and [Fe/H] where known
    import gaia_cal as G
    g = res["gaia"]
    b = g["baseline"]
    p = {n: b[n]["median"] for n in ("mu", "gamma", "R", "d1O") if n in b}
    # use the bias-corrected calibration (see analysis.injection_test)
    p["mu"] = s["lmc_distance"]["mu"]
    p["gamma"] = s["lmc_distance"]["gamma_corrected"]
    p["R"] = s["lmc_distance"]["R_corrected"]
    mw = G.milky_way()
    feh, _ = analysis.mw_metallicities(mw)
    feh0 = np.where(np.isfinite(feh), feh, 0.0)
    R = p["R"]
    w = mw["int_average_g"] - R * (mw["int_average_bp"] - mw["int_average_rp"])
    shp = g["lmc_shapes"]
    is1o = mw["mode"] == "1O"
    isF = np.char.startswith(mw["mode"].astype(str), "F")
    a = np.where(is1o, shp["1O"]["a_g"] - R * shp["1O"]["a_c"], shp["F"]["a_g"] - R * shp["F"]["a_c"])
    bb = np.where(is1o, shp["1O"]["b_g"] - R * shp["1O"]["b_c"], shp["F"]["b_g"] - R * shp["F"]["b_c"])
    M = a + bb * (mw["logP"] - 1) - p["mu"] + p["gamma"] * (feh0 + 0.409) + np.where(is1o, p.get("d1O", 0.0), 0.0)
    mu = w - M
    ok = (is1o | isF) & np.isfinite(mu)
    out["mw"] = {"ok": ok, "d": 10 ** (0.2 * mu + 1) / 1e3, "ra": mw["ra"], "dec": mw["dec"], "id": mw["name"],
                 "P": mw["P"], "mode": mw["mode"], "has_feh": np.isfinite(feh)}
    return out


def write(res) -> str:
    import plotly.graph_objects as go
    L, D = THEMES["light"], THEMES["dark"]
    dist = star_distances(res)
    fig = go.Figure()
    styles = []
    names = {"lmc": "Large Magellanic Cloud", "smc": "Small Magellanic Cloud", "mw": "Milky Way (Gaia DR3)"}
    for key in ("mw", "lmc", "smc"):
        d = dist[key]
        ok = d["ok"]
        x, y, z = _galactic_xyz(d["ra"][ok], d["dec"][ok], d["d"][ok])
        custom = np.column_stack([d["id"][ok], np.round(d["P"][ok], 3), d["mode"][ok], np.round(d["d"][ok], 2)])
        fig.add_trace(go.Scatter3d(
            x=np.round(x, 3), y=np.round(y, 3), z=np.round(z, 3), mode="markers", name=f"{names[key]} ({ok.sum():,})",
            marker=dict(size=2.2 if key != "mw" else 2.6, color=L[key], opacity=0.75),
            customdata=custom,
            hovertemplate="%{customdata[0]}<br>period %{customdata[1]} d (%{customdata[2]})"
                          "<br>distance %{customdata[3]} kpc<extra></extra>"))
        styles.append({"light": {"marker.color": L[key]}, "dark": {"marker.color": D[key]}})
    # the Sun and the Galactic centre
    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers+text", text=["Sun"], textposition="bottom center",
                               marker=dict(size=5, color=L["ink"]), name="Sun", showlegend=False,
                               textfont=dict(color=L["ink"]), hovertemplate="the Sun<extra></extra>"))
    styles.append({"light": {"marker.color": L["ink"], "textfont.color": L["ink"]},
                   "dark": {"marker.color": D["ink"], "textfont.color": D["ink"]}})
    fig.add_trace(go.Scatter3d(x=[8.2], y=[0], z=[0], mode="markers+text", text=["Galactic centre"],
                               textposition="top center", marker=dict(size=4, color=L["muted"], symbol="diamond"),
                               name="Galactic centre", showlegend=False, textfont=dict(color=L["ink2"]),
                               hovertemplate="Galactic centre (8.2 kpc)<extra></extra>"))
    styles.append({"light": {"textfont.color": L["ink2"]}, "dark": {"textfont.color": D["ink2"]}})

    axis = dict(title="", showbackground=False, showgrid=True, gridcolor="rgba(137,135,129,0.25)",
                zeroline=False, color=L["muted"])
    scenes = {
        "everything": dict(xaxis=dict(range=[-25, 25]), yaxis=dict(range=[-60, 25]), zaxis=dict(range=[-60, 15])),
    }
    cams = {
        "the whole map": dict(eye=dict(x=1.5, y=-1.2, z=0.7)),
        "the Clouds up close": dict(eye=dict(x=0.35, y=-0.75, z=-0.35)),
        "the Milky Way from above": dict(eye=dict(x=0.0, y=0.0, z=2.2), up=dict(x=1, y=0, z=0)),
    }
    buttons = [dict(label=f"view: {name}", method="relayout", args=[{"scene.camera": cam}]) for name, cam in cams.items()]
    fig.update_layout(
        height=820, margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor=L["surface"],
        font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif", size=12, color=L["ink"]),
        scene=dict(xaxis={**axis, "title": "x (kpc), towards the Galactic centre"},
                   yaxis={**axis, "title": "y (kpc), direction of rotation"},
                   zaxis={**axis, "title": "z (kpc), Galactic north"}, aspectmode="data", bgcolor=L["surface"],
                   camera=cams["the whole map"]),
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(0,0,0,0)", font=dict(size=12), itemsizing="constant"),
        hoverlabel=dict(bgcolor=L["surface"], font=dict(color=L["ink"], size=12)),
        updatemenus=[dict(buttons=buttons, direction="down", x=0.99, xanchor="right", y=0.98, yanchor="top",
                          bgcolor=L["menu_bg"], bordercolor=L["menu_border"], font=dict(color=L["ink"], size=12),
                          showactive=True)])
    theme_layout = {name: {"paper_bgcolor": th["surface"], "scene.bgcolor": th["surface"], "font.color": th["ink"],
                           "hoverlabel.bgcolor": th["surface"], "hoverlabel.font.color": th["ink"],
                           "updatemenus[0].bgcolor": th["menu_bg"], "updatemenus[0].bordercolor": th["menu_border"],
                           "updatemenus[0].font.color": th["ink"]} for name, th in THEMES.items()}
    themes_js = json.dumps({name: {"layout": theme_layout[name], "traces": [s[name] for s in styles]}
                            for name in THEMES})
    plot_html = fig.to_html(full_html=False, include_plotlyjs=True, div_id="cepheids3d",
                            config={"displaylogo": False, "responsive": True})
    import analysis
    s = analysis.summarize(res)
    dl = s["lmc_distance"]
    n_mc = int(dist["lmc"]["ok"].sum() + dist["smc"]["ok"].sum())
    n_mw = int(dist["mw"]["ok"].sum())
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Cepheids in 3D</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']}; --ink2: {L['ink2']}; --muted: {L['muted']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']}; --ink2: {D['ink2']}; --muted: {D['muted']}; }}
  }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 80ch; }}
  .note {{ font-size: 0.85rem; color: var(--muted); }}
</style>
</head>
<body>
<main>
<h1>Every Cepheid, placed by its own heartbeat</h1>
<p>{n_mc + n_mw:,} pulsating stars, each at the distance its pulsation period and reddening-free brightness
imply: {n_mw:,} in our Milky Way and {n_mc:,} in the two Magellanic Clouds. The period tells us how
luminous a Cepheid really is (the Leavitt law); comparing that with how bright it looks gives its
distance. The law's zero point comes from Gaia parallaxes of Milky Way Cepheids, which put the Large
Magellanic Cloud at {dl['d_kpc']:.1f} ± {dl['d_err_kpc']:.1f} kpc (about {round(dl['d_kpc'] * 3261.56, -3):,.0f} light-years).</p>
<p>Drag to rotate, scroll to zoom, hover any star for its name, period and distance; the menu at top
right offers three views. Up close, the Large Cloud is a tilted disk (its north-east side nearer to us),
and the Small Cloud is stretched along our line of sight; the Milky Way's Cepheids trace the young disk
and its spiral arms, with a warp at large radii.</p>
{plot_html}
<p class="note">Data: OGLE-IV Collection of Variable Stars (Magellanic Clouds; Soszyński et al. 2015-2023) and
OGLE's list of Galactic classical Cepheids (Pietrukowicz et al. 2021, 2026 update); Gaia DR3 photometry.
Single-star distances carry the law's ~0.07 mag (3%) scatter, which smears each Cloud along the line of
sight by about 1.5 kpc; stars more than 3σ off the law (mostly blends) and SMC fundamental-mode stars
shorter than 2.5 days (where the SMC's law bends) are left out. The SMC uses the LMC's law with no metallicity correction (its metal-poor Cepheids
are a few per cent fainter, so it appears ~2% farther). Milky Way distances use each star's spectroscopic
metallicity where known (else solar) and are least certain for the most reddened stars. Coordinates:
heliocentric Galactic, kpc. Made by run.py in 47-cepheids.</p>
</main>
<script>
(function () {{
  var gd = document.getElementById("cepheids3d");
  var THEMES = {themes_js};
  function apply(dark) {{
    var th = THEMES[dark ? "dark" : "light"];
    th.traces.forEach(function (s, i) {{ Plotly.restyle(gd, s, [i]); }});
    Plotly.relayout(gd, th.layout);
  }}
  var mq = window.matchMedia("(prefers-color-scheme: dark)");
  if (mq.matches) apply(true);
  mq.addEventListener("change", function (e) {{ apply(e.matches); }});
}})();
</script>
</body>
</html>
"""
    out = C.FIGURE_DIR / "cepheids_3d.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
