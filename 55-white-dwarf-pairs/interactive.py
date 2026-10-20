"""The interactive page: LISA's sensitivity with every known binary, for 3 months to 4 years.

A self-contained HTML file (plotly.js and the data embedded; works offline). Buttons
switch the observing time: each binary's height grows as sqrt(T) while the Galactic
foreground fades as more sources are resolved, so the height above the curve is
always the sky-averaged signal-to-noise ratio. Hover any binary for its parameters;
the table below the plot holds every number. It follows the viewer's light/dark
setting, with colours re-stepped for the dark surface.
"""

from __future__ import annotations

import html
import json

import numpy as np

import config as C
import gw
import plotstyle as ps

THEMES = {
    "light": dict(surface=ps.SURFACE, ink=ps.INK, ink2=ps.INK_SECONDARY, muted=ps.MUTED, grid=ps.GRID,
                  axis=ps.AXIS, blue=ps.BLUE, orange=ps.ORANGE, aqua=ps.AQUA, violet=ps.VIOLET,
                  context=ps.CONTEXT, menu_bg="#f0efec", menu_border="#c3c2b7", menu_active="#d9d8d1"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#2c2c2a",
                 axis="#383835", blue="#3987e5", orange="#d95926", aqua="#199e70", violet="#9085e9",
                 context="#52514e", menu_bg="#2c2c2a", menu_border="#383835", menu_active="#52514e"),
}
TYPES = [("AMCVn", "blue", "circle", "AM CVn (transferring mass)"),
         ("DWD", "orange", "square", "detached white-dwarf pair"),
         ("sdB", "aqua", "triangle-up", "hot subdwarf + white dwarf"),
         ("UCXB", "ink2", "pentagon", "ultracompact X-ray binary")]
TIMES = [("3 months", 0.25), ("1 year", 1.0), ("4 years", 4.0)]


def _fac(r):
    return 0.8 if r["incl_kind"] == "assumed" else float(gw.inclination_factor(np.radians(r["incl"])))


def build(res, art) -> str:
    import plotly.graph_objects as go

    b = art["binaries"]
    L, D = THEMES["light"], THEMES["dark"]
    f = np.logspace(-4, -1, 500)
    popd = art["pop_default"]
    rng = np.random.default_rng(5)
    k = rng.choice(popd["resolved_f"].size, size=min(4000, popd["resolved_f"].size), replace=False)
    pf, pp = popd["resolved_f"][k], popd["resolved_pw"][k]

    def heights(t_yr):
        T = t_yr * C.YEAR
        out = {}
        for typ, *_ in TYPES:
            for status in ("verification", "detectable"):
                s = b[(b["type"] == typ) & (b["status"] == status)]
                out[(typ, status)] = [float(np.sqrt(4 * _fac(r) * r["f_Hz"] * T) * r["A"]) for r in s.to_dict("records")]
        j = b[b["name"] == "ZTFJ1539"].iloc[0]
        out["j1539"] = [float(np.sqrt(4 * gw.inclination_factor(np.radians(84.15)) * j["f_Hz"] * T) * j["A"])]
        out["pop"] = list(np.sqrt(4 * pf * T * pp))
        out["noise"] = list(gw.characteristic_strain_noise(f, t_yr, confusion=True))
        return out

    fig = go.Figure()
    styles = {"light": [], "dark": []}
    # trace 0: resolved model population (context)
    fig.add_trace(go.Scatter(x=pf, y=heights(4.0)["pop"], mode="markers", name="model Galaxy: pairs LISA resolves",
                               marker=dict(size=3, color=L["context"]), hoverinfo="skip"))
    for th in ("light", "dark"):
        styles[th].append({"marker.color": THEMES[th]["context"]})
    # trace 1: instrument noise; trace 2: instrument + foreground (time dependent)
    fig.add_trace(go.Scatter(x=f, y=gw.characteristic_strain_noise(f, 4.0, confusion=False), mode="lines",
                             name="LISA instrument noise", line=dict(color=L["ink"], width=2), hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=f, y=heights(4.0)["noise"], mode="lines", name="+ Galactic foreground (Robson et al.)",
                             line=dict(color=L["violet"], width=2), hoverinfo="skip"))
    for th in ("light", "dark"):
        styles[th].append({"line.color": THEMES[th]["ink"]})
        styles[th].append({"line.color": THEMES[th]["violet"]})
    order = []
    for typ, col, sym, label in TYPES:
        for status in ("verification", "detectable"):
            s = b[(b["type"] == typ) & (b["status"] == status)]
            if s.empty:
                continue
            order.append((typ, status))
            filled = status == "verification"
            text = []
            for r in s.to_dict("records"):
                mk = lambda kind: " (assumed)" if kind == "assumed" else ""
                text.append(
                    f"<b>{html.escape(r['name'])}</b><br>{label}, {status}<br>"
                    f"orbit {r['P_s'] / 60:.2f} min, f = {r['f_Hz'] * 1e3:.3f} mHz<br>"
                    f"masses {r['m1']:.2f}{mk(r['m1_kind'])} + {r['m2']:.3f}{mk(r['m2_kind'])} M☉<br>"
                    f"distance {r['d_pc']:,.0f} pc, inclination {r['incl']:.0f}°{mk(r['incl_kind'])}<br>"
                    f"SNR: {r['snr_3mo']:.0f} (3 mo), {r['snr_1yr']:.0f} (1 yr), {r['snr_4yr']:.0f} (4 yr)")
            fig.add_trace(go.Scatter(
                x=s["f_Hz"], y=heights(4.0)[(typ, status)], mode="markers",
                name=f"{label}, {'verification' if filled else 'detectable'}",
                marker=dict(symbol=sym if filled else f"{sym}-open", size=11, color=L[col],
                            line=dict(color=L["surface"] if filled else L[col], width=1 if filled else 2)),
                text=text, hovertemplate="%{text}<extra></extra>"))
            for th in ("light", "dark"):
                c = THEMES[th][col]
                styles[th].append({"marker.color": c, "marker.line.color": THEMES[th]["surface"] if filled else c})
    j = b[b["name"] == "ZTFJ1539"].iloc[0]
    fig.add_trace(go.Scatter(x=[j["f_Hz"]], y=heights(4.0)["j1539"], mode="markers+text", name="ZTF J1539+5027",
                             text=["ZTF J1539+5027"], textposition="top center",
                             textfont=dict(color=L["ink"], size=12),
                             marker=dict(symbol="star", size=20, color=L["ink"], line=dict(color=L["surface"], width=1)),
                             hovertemplate=(f"<b>ZTF J1539+5027</b><br>6.91-minute eclipsing pair<br>"
                                            f"our dP/dt = {res['joint_ztf_ptf']['Pdot'] * 1e11:.4f}e-11 s/s<br>"
                                            f"SNR {j['snr_4yr']:.0f} in 4 yr<extra></extra>")))
    for th in ("light", "dark"):
        styles[th].append({"marker.color": THEMES[th]["ink"], "marker.line.color": THEMES[th]["surface"],
                           "textfont.color": THEMES[th]["ink"]})
    # buttons: observing time
    buttons = []
    for label, t_yr in TIMES:
        h = heights(t_yr)
        ys = [h["pop"], list(gw.characteristic_strain_noise(f, 4.0, confusion=False)), h["noise"]]
        ys += [h[key] for key in order] + [h["j1539"]]
        buttons.append(dict(label=f"LISA observes for {label}", method="update",
                            args=[{"y": ys}, {"title.text": f"LISA's view after {label}"}]))
    axis = dict(type="log", gridcolor=L["grid"], linecolor=L["axis"], zeroline=False, ticks="outside",
                tickcolor=L["muted"], tickfont=dict(color=L["muted"]), dtick=1, exponentformat="power")
    fig.update_layout(
        height=700, margin=dict(l=70, r=20, t=90, b=60), paper_bgcolor=L["surface"], plot_bgcolor=L["surface"],
        font=dict(family="system-ui, -apple-system, 'Segoe UI', sans-serif", size=12, color=L["ink"]),
        xaxis={**axis, "title": "gravitational-wave frequency (Hz)", "range": [-4, -1]},
        yaxis={**axis, "title": "characteristic strain", "range": [-21.3, -17.7]},
        legend=dict(x=0.99, xanchor="right", y=0.99, bgcolor="rgba(0,0,0,0)", font=dict(size=11)),
        hoverlabel=dict(bgcolor=L["surface"], font=dict(color=L["ink"], size=12)),
        updatemenus=[dict(type="buttons", direction="right", buttons=buttons, active=2, x=0.33, xanchor="left", y=1.08,
                          yanchor="bottom", bgcolor=L["menu_bg"], bordercolor=L["menu_border"],
                          font=dict(color=L["ink"], size=12), showactive=False)],
        title=dict(text="LISA's view after 4 years", x=0.03, xanchor="left", y=0.965, font=dict(size=15)))
    theme_layout = {name: {"paper_bgcolor": t["surface"], "plot_bgcolor": t["surface"], "font.color": t["ink"],
                           "xaxis.gridcolor": t["grid"], "yaxis.gridcolor": t["grid"], "xaxis.linecolor": t["axis"],
                           "yaxis.linecolor": t["axis"], "hoverlabel.bgcolor": t["surface"],
                           "hoverlabel.font.color": t["ink"], "updatemenus[0].bgcolor": t["menu_bg"],
                           "updatemenus[0].bordercolor": t["menu_border"], "updatemenus[0].font.color": t["ink"]}
                    for name, t in THEMES.items()}
    themes_js = json.dumps({name: {"layout": theme_layout[name], "traces": styles[name]} for name in THEMES})
    plot_html = fig.to_html(full_html=False, include_plotlyjs=True, div_id="lisa",
                            config={"displaylogo": False, "responsive": True})
    rows = []
    for r in b.sort_values("f_Hz", ascending=False).to_dict("records"):
        rows.append(f"<tr><td>{html.escape(r['name'])}</td><td>{html.escape(r['type'])}</td><td>{r['status']}</td>"
                    f"<td>{r['P_s'] / 60:.2f}</td><td>{r['f_Hz'] * 1e3:.3f}</td><td>{r['Mc']:.3f}</td>"
                    f"<td>{r['d_pc']:,.0f}</td><td>{r['snr_3mo']:.1f}</td><td>{r['snr_1yr']:.1f}</td>"
                    f"<td>{r['snr_4yr']:.1f}</td></tr>")
    e = res["joint_ztf_ptf"]
    ls = res["lisa_summary"]
    p = res["population"]["default"]
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LISA's guaranteed sources</title>
<style>
  :root {{ color-scheme: light; --surface: {L['surface']}; --ink: {L['ink']}; --ink2: {L['ink2']}; --muted: {L['muted']}; --grid: {L['grid']}; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ color-scheme: dark; --surface: {D['surface']}; --ink: {D['ink']}; --ink2: {D['ink2']}; --muted: {D['muted']}; --grid: {D['grid']}; }}
  }}
  body {{ margin: 0; background: var(--surface); color: var(--ink);
         font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 24px 16px 40px; }}
  h1 {{ font-size: 1.35rem; margin: 0 0 6px; }}
  h2 {{ font-size: 1.05rem; margin: 24px 0 6px; }}
  p {{ margin: 0 0 10px; color: var(--ink2); max-width: 80ch; }}
  .note {{ font-size: 0.85rem; color: var(--muted); }}
  .tablewrap {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; font-size: 0.85rem; font-variant-numeric: tabular-nums; }}
  th, td {{ padding: 4px 10px; border-bottom: 1px solid var(--grid); text-align: right; white-space: nowrap; }}
  th:nth-child(-n+3), td:nth-child(-n+3) {{ text-align: left; }}
  th {{ color: var(--ink2); font-weight: 600; }}
</style>
</head>
<body>
<main>
<h1>The binaries LISA is guaranteed to hear</h1>
<p>Pairs of white dwarfs orbiting each other in minutes radiate gravitational waves at twice their orbital
frequency, in the millihertz band of LISA, the space observatory ESA plans to launch around 2035. Telescopes
have already found {ls['n_binaries']} such binaries strong enough to detect (Kupfer et al. 2024). Each is drawn
here at a height chosen so that its distance above the noise curve is its signal-to-noise ratio. Switch
the observing time and watch the sources rise out of the noise. Within three months, between {ls['n_snr3mo_gt7']}
(with the Galactic foreground expected early in the mission) and {ls['n_snr3mo_gt7_instrument_only']} (instrument noise
alone) pass a signal-to-noise ratio of 7 in our estimate.</p>
<p>The star is ZTF J1539+5027, whose 6.9-minute orbit we timed from public survey data: it is shrinking by
{abs(e['Pdot']) * C.YEAR * 1e6:.0f} microseconds a year, exactly as gravitational-wave emission requires. The grey cloud is
the {p['n_resolved_snr7_4yr']:,} pairs (of {p['n_in_band'] / 1e6:.0f} million in the band) that our observationally driven
model of the Galaxy predicts LISA will pick out individually after four years.</p>
{plot_html}
<p class="note">Sensitivity: the analytic sky-averaged curve of Robson, Cornish &amp; Liu (2019), with their fit to the
Galactic foreground for each observing time (their 6-month curve stands in for 3 months). Binary parameters
(masses, inclinations, Gaia DR3 distances) from Kupfer et al. (2024); where they list a value in brackets
(assumed) we use it as given, and unknown inclinations are averaged over. Sky-averaged SNRs ignore where each
binary sits relative to LISA's antenna pattern, which shifts individual values by tens of per cent.
Made by run.py in 55-white-dwarf-pairs.</p>
<h2>Table view</h2>
<div class="tablewrap"><table>
<thead><tr><th>binary</th><th>type</th><th>Kupfer et al.</th><th>P (min)</th><th>f (mHz)</th><th>chirp mass (M☉)</th>
<th>distance (pc)</th><th>SNR 3 mo</th><th>SNR 1 yr</th><th>SNR 4 yr</th></tr></thead>
<tbody>{''.join(rows)}</tbody></table></div>
</main>
<script>
(function () {{
  var gd = document.getElementById("lisa");
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
    out = C.FIGURE_DIR / "lisa_explorer.html"
    out.write_text(page, encoding="utf-8")
    return str(out)
