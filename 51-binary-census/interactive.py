"""An interactive page: explore the binary census yourself.

A self-contained HTML file (plotly.js and the data are embedded, so it works
offline and can be shared as one file) with two linked ideas:

1. the period-eccentricity plane of Gaia DR3's spectroscopic orbits and SB9's
   binaries, where tides have erased the eccentricity of the closest orbits
   (hover an SB9 point for the star's name and orbit);
2. how many companions Sun-like stars have at each orbital period: Gaia DR3
   corrected for what it cannot see, against Raghavan et al. 2010, Tokovinin
   2014 and Moe & Di Stefano 2017 (click legend entries to toggle them).

It follows the viewer's light/dark setting, with colours re-stepped for each
surface.
"""

from __future__ import annotations

import json

import numpy as np

import config as C
import literature as L
import plotstyle as ps

THEMES = {
    "light": dict(surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED, grid="#e1e0d9",
                  blue=ps.BLUE, orange=ps.ORANGE, aqua=ps.AQUA, band="rgba(42,120,214,0.18)",
                  density=["#e8f1fd", "#9ec5f4", "#3987e5", "#184f95", "#0d366b"]),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#383835",
                 blue="#3987e5", orange="#d95926", aqua="#199e70", band="rgba(57,135,229,0.25)",
                 density=["#1f2a3a", "#1c4a80", "#2a78d6", "#6da7ec", "#cde2fb"]),
}


def _colorscale(steps):
    return [[i / (len(steps) - 1), c] for i, c in enumerate(steps)]


def build(res) -> str:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    Lt, Dk = THEMES["light"], THEMES["dark"]
    g = res["gaia_circ"]
    s9 = res["sb9"]
    sc = res["solar"]

    fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], horizontal_spacing=0.09,
                        subplot_titles=("Tides circularize the closest orbits",
                                        "Companions of Sun-like stars, per decade of period"))
    styles = []  # per-trace restyle for each theme, in trace order

    # 1. period-eccentricity: Gaia as a density, SB9 as hoverable points
    lp = np.log10(g["P"])
    e = np.where(g["circ"], 0.0, g["e"])
    xe, ye = np.arange(-0.6, 3.6001, 0.05), np.arange(-0.02, 1.0001, 0.0125)
    h, _, _ = np.histogram2d(lp, e, bins=[xe, ye])
    z = np.where(h > 0, np.log10(np.maximum(h, 1)), np.nan).T  # log colour scale: counts span 1-1000
    fig.add_trace(go.Heatmap(x=0.5 * (xe[1:] + xe[:-1]), y=0.5 * (ye[1:] + ye[:-1]), z=z, customdata=h.T,
                             colorscale=_colorscale(Lt["density"]), showscale=False, name="Gaia DR3 orbits",
                             showlegend=True, hoverongaps=False,
                             hovertemplate="log P %{x:.2f}, e %{y:.2f}: %{customdata:.0f} Gaia orbits<extra></extra>"),
                  row=1, col=1)
    styles.append({"light": {"colorscale": [_colorscale(Lt["density"])]},
                   "dark": {"colorscale": [_colorscale(Dk["density"])]}})
    t = s9["table"]
    sel = s9["fgk_dwarf"]
    names = [str(n) for n in np.asarray(t["name"])[sel]]
    sp = [str(x) for x in np.asarray(t["sp1"])[sel]]
    P9, E9 = np.asarray(t["period"])[sel], np.asarray(t["ecc"])[sel]
    text = [f"{n} ({s})<br>P = {p:.3f} d, e = {ee:.3f}" for n, s, p, ee in zip(names, sp, P9, E9)]
    fig.add_trace(go.Scatter(x=np.log10(P9), y=E9, mode="markers", name="SB9 binaries (F, G, K dwarfs)",
                               marker=dict(size=6, color=Lt["orange"], line=dict(width=0.5, color=Lt["surface"])),
                               text=text, hovertemplate="%{text}<extra></extra>"), row=1, col=1)
    styles.append({"light": {"marker.color": Lt["orange"], "marker.line.color": Lt["surface"]},
                   "dark": {"marker.color": Dk["orange"], "marker.line.color": Dk["surface"]}})
    import binaries as B
    pp = np.logspace(-0.6, 3.6, 300)
    fig.add_trace(go.Scatter(x=np.log10(pp), y=B.e_max_envelope(pp), mode="lines", name="e_max envelope",
                             line=dict(color=Lt["muted"], dash="dash", width=1.5), hoverinfo="skip"), row=1, col=1)
    styles.append({"light": {"line.color": Lt["muted"]}, "dark": {"line.color": Dk["muted"]}})
    p_prime = g["all"]["p_circ"] - B.meibom_mathieu_pcirc(0.0)
    fig.add_trace(go.Scatter(x=np.log10(pp), y=B.meibom_mathieu_e(pp, p_prime), mode="lines",
                             name=f"circularization fit (P_circ = {g['all']['p_circ']:.1f} d)",
                             line=dict(color=Lt["ink"], width=2), hoverinfo="skip"), row=1, col=1)
    styles.append({"light": {"line.color": Lt["ink"]}, "dark": {"line.color": Dk["ink"]}})

    # 2. the period distribution
    comp = sc["completeness"]
    x = comp["centres"]
    f, fe, fi = sc["f_corr"], sc["f_err"], sc["f_ideal"]
    ok = np.isfinite(f)
    lo = np.fmin(f, fi)
    hi = np.fmax(f, fi)
    short = ok & (x < 3.5)  # the two acceptance assumptions differ only where spectroscopy matters
    xb = np.concatenate([x[short], x[short][::-1]])
    yb = np.concatenate([hi[short], lo[short][::-1]])
    fig.add_trace(go.Scatter(x=xb, y=yb, fill="toself", fillcolor=Lt["band"], line=dict(width=0), mode="lines",
                             name="range from DR3's spectroscopic acceptance", hoverinfo="skip"), row=1, col=2)
    styles.append({"light": {"fillcolor": Lt["band"]}, "dark": {"fillcolor": Dk["band"]}})
    fig.add_trace(go.Scatter(x=x[ok], y=f[ok], mode="markers", name="Gaia DR3, corrected (this work)",
                             error_y=dict(type="data", array=fe[ok], color=Lt["blue"], thickness=1.2, width=0),
                             marker=dict(size=7, color=Lt["blue"]),
                             customdata=np.column_stack([comp["period_known"][ok], (sc["h_orbit"] + sc["h_wide"])[ok]]),
                             hovertemplate="log P %{x:.2f}: %{y:.3f} companions per star per dex<br>"
                                           "%{customdata[1]:.0f} seen, Gaia catches %{customdata[0]:.0%}<extra></extra>"),
                  row=1, col=2)
    styles.append({"light": {"marker.color": Lt["blue"], "error_y.color": Lt["blue"]},
                   "dark": {"marker.color": Dk["blue"], "error_y.color": Dk["blue"]}})
    r10 = res["raghavan"]
    ed = r10["hist_edges"]
    xs = np.repeat(ed, 2)[1:-1]
    ys = np.repeat(r10["per_primary_per_dex"], 2)
    fig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", name="Raghavan et al. 2010 (25 pc)",
                             line=dict(color=Lt["ink"], width=2), hoverinfo="skip"), row=1, col=2)
    styles.append({"light": {"line.color": Lt["ink"]}, "dark": {"line.color": Dk["ink"]}})
    xx = np.linspace(-0.5, 9.0, 300)
    gauss = lambda m, s_: np.exp(-0.5 * ((xx - m) / s_) ** 2) / (s_ * np.sqrt(2 * np.pi))
    for name, y, dash, visible in (
            ("Raghavan's log-normal", 259 / 454 * gauss(L.RAGHAVAN_LOGP_MEAN, L.RAGHAVAN_LOGP_SIGMA), "dash", True),
            ("Tokovinin 2014 (67 pc)", L.TOKOVININ_ALL["eps"] * gauss(L.TOKOVININ_ALL["x0"], L.TOKOVININ_ALL["sigma"]),
             "dot", "legendonly")):
        fig.add_trace(go.Scatter(x=xx, y=y, mode="lines", name=name, visible=visible,
                                 line=dict(color=Lt["ink"], dash=dash, width=1.5), hoverinfo="skip"), row=1, col=2)
        styles.append({"light": {"line.color": Lt["ink"]}, "dark": {"line.color": Dk["ink"]}})
    mds = np.where((x >= 0.2) & (x < 8), sc["models"]["mds"], np.nan)
    fig.add_trace(go.Scatter(x=x, y=mds, mode="lines", name="Moe & Di Stefano 2017", visible="legendonly",
                             line=dict(color=Lt["muted"], dash="dashdot", width=1.5), hoverinfo="skip"), row=1, col=2)
    styles.append({"light": {"line.color": Lt["muted"]}, "dark": {"line.color": Dk["muted"]}})

    axis = dict(gridcolor=Lt["grid"], zeroline=False, color=Lt["muted"], linecolor=Lt["grid"])
    fig.update_xaxes(title_text="orbital period, log10(P / day)", **axis)
    fig.update_yaxes(**axis)
    fig.update_yaxes(title_text="eccentricity", range=[-0.03, 1.0], row=1, col=1)
    fig.update_xaxes(range=[-0.6, 3.6], row=1, col=1)
    fig.update_yaxes(title_text="companions per star per decade", range=[0, 0.13], row=1, col=2)
    fig.update_xaxes(range=[-0.5, 9.0], row=1, col=2)
    fig.update_layout(height=620, margin=dict(l=60, r=20, t=60, b=50), paper_bgcolor=Lt["surface"],
                      plot_bgcolor=Lt["surface"],
                      font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif", size=12, color=Lt["ink"]),
                      legend=dict(orientation="h", y=-0.18, x=0, bgcolor="rgba(0,0,0,0)"),
                      hoverlabel=dict(bgcolor=Lt["surface"], font=dict(color=Lt["ink"], size=12)))
    theme_layout = {}
    for name, th in THEMES.items():
        lay = {"paper_bgcolor": th["surface"], "plot_bgcolor": th["surface"], "font.color": th["ink"],
               "hoverlabel.bgcolor": th["surface"], "hoverlabel.font.color": th["ink"]}
        for ax_name in ("xaxis", "yaxis", "xaxis2", "yaxis2"):
            lay[f"{ax_name}.gridcolor"] = th["grid"]
            lay[f"{ax_name}.color"] = th["muted"]
            lay[f"{ax_name}.linecolor"] = th["grid"]
        for k in range(2):
            lay[f"annotations[{k}].font.color"] = th["ink"]
        theme_layout[name] = lay
    themes_js = json.dumps({name: {"layout": theme_layout[name], "traces": [s[name] for s in styles]}
                            for name in THEMES})
    plot_html = fig.to_html(full_html=False, include_plotlyjs=True, div_id="census",
                            config={"displaylogo": False, "responsive": True})

    obs = sc["observed"]
    n = sc["n_primaries"]
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Binary census explorer</title>
<style>
  :root {{ color-scheme: light; --surface: {Lt['surface']}; --ink: {Lt['ink']}; --ink2: {Lt['ink2']}; --muted: {Lt['muted']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {Dk['surface']}; --ink: {Dk['ink']}; --ink2: {Dk['ink2']}; --muted: {Dk['muted']}; }}
  }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1240px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 85ch; }}
  .note {{ font-size: 0.85rem; color: var(--muted); }}
</style>
</head>
<body>
<main>
<h1>Is the lonely Sun typical? A census of binary stars</h1>
<p>Of the {n:,} Sun-like stars within 100 parsecs, Gaia DR3 shows some sign of a companion around
{100 * obs['any_flag'] / n:.0f}%: an orbit measured from radial velocities or from the star's wobble on the
sky, a curving path, a poor single-star fit, or a resolved partner. Each method sees only part of the
range of orbits, so the counts are corrected with a simulation of what Gaia could have detected.</p>
<p><b>Left:</b> orbits closer than about 10 days are circular, because tides between the two stars have
drained their eccentricity; wider orbits keep the eccentricities they were born with. <b>Right:</b>
companions per star per decade of orbital period. Hover for details; click legend entries to switch the
published comparisons on and off.</p>
{plot_html}
<p class="note">Data: Gaia DR3 non-single-star tables (Gaia Collaboration, Arenou et al. 2023); SB9 (Pourbaix et al.
2004); El-Badry, Rix &amp; Heintz (2021) wide binaries; Raghavan et al. (2010) Table 13. Gaia's orbits are
dwarfs within 500 pc with significance &gt; 20. Gaia cannot measure periods between ~3,000 and ~100,000 days
(the gap): there it sees companions only as accelerations and excess RUWE. Made by run.py in
51-binary-census.</p>
</main>
<script>
(function () {{
  var gd = document.getElementById("census");
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
    out = C.FIGURE_DIR / "binary_census.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
