"""The wide-binary explorer: choose a sample and its cuts, watch the test respond.

A self-contained HTML page (plotly.js and the data embedded; works offline).
For the chosen sample and cuts it shows

* the median scaled velocity v~ in bins of r_sky / r_M (the model-light test
  of Banik et al. 2024 and Cookson et al. 2026), with the Newtonian
  expectation (flat, normalised to the two innermost bins) and MOND's (the
  same, times the rise our forward model predicts with the Galactic field);
* the v~ distributions of tight (r_sky < 0.3 r_M) and wide (r_sky > r_M) pairs;
* the wide/tight ratio of medians, which MOND predicts near 1.2 and Newton near 1.

Samples: Chae (2023)'s 26,615 pairs (his public file), our rebuild of Banik
et al. (2024)'s 8,228, and our Cookson et al. (2026)-style 1,222.
"""

from __future__ import annotations

import json

import numpy as np

import plotstyle as PS
from config import DERIVED_DIR, FIG_DIR

THEMES = {
    "light": dict(surface=PS.SURFACE, ink=PS.INK, ink2=PS.INK_SECONDARY, muted=PS.MUTED, grid=PS.GRID,
                  newton=PS.BLUE, mond=PS.ORANGE, third=PS.AQUA, panel="#f0efec", border="#c3c2b7"),
    "dark": dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", muted="#898781", grid="#383835",
                 newton="#3987e5", mond="#d95926", third="#199e70", panel="#2c2c2a", border="#383835"),
}


def _round(a, nd=4):
    return [None if not np.isfinite(x) else round(float(x), nd) for x in np.asarray(a, float)]


def payload():
    import analysis as A
    import chae as CH
    rm = lambda m: A.RM1_KAU * np.sqrt(m)
    out = {}
    d = CH.chae_sample(CH.load_chae_base())
    sc = A.chae_for_forward(d)
    both = (d["RV1"] > -1000) & (d["RV2"] > -1000)
    out["chae2023"] = dict(
        label="Chae (2023): 26,615 pairs from El-Badry et al., within 200 pc, 0.2–30 kAU (no perspective correction)",
        x=_round(sc["rsky_kau"] / rm(sc["m1"] + sc["m2"])), vt=_round(sc["vt"]), sig=_round(sc["sig_vt"]),
        d=_round(d["dM"], 1), ruwe=_round(np.maximum(d["ruwe1"], d["ruwe2"]), 3), rv2=both.astype(int).tolist())
    for name, label in (("banik2024", "our rebuild of Banik et al. (2024): 8,228 pairs within 250 pc, 2–30 kAU"),
                        ("cookson2026", "Cookson et al. (2026)-style: 1,222 pairs within 130 pc, 1–30 kAU")):
        sub, _ = A.sample(name)
        both = np.isfinite(sub["radial_velocity1"]) & np.isfinite(sub["radial_velocity2"])
        out[name] = dict(label=label, x=_round(sub["rsky_kau"] / rm(sub["mtot"])), vt=_round(sub["vt"]),
                         sig=_round(sub["vt_err"]), d=_round(sub["dist"], 1),
                         ruwe=_round(np.maximum(sub["ruwe1"], sub["ruwe2"]), 3), rv2=both.astype(int).tolist())
    # MOND/Newton ratio of medians vs x from the forward model (v~ < 2), if computed
    ratio = None
    p = DERIVED_DIR / "results" / "median_profiles.json"
    if p.exists():
        mp = json.loads(p.read_text())["banik2024"]["cuts"]["2.0"]
        xs = np.array(mp["newton"]["x"])
        rr = np.array(mp["mond"]["median"]) / np.array(mp["newton"]["median"])
        ok = np.isfinite(rr)
        ratio = dict(x=_round(xs[ok]), r=_round(rr[ok]))
    if ratio is None:   # fall back to the angle-averaged boost: sqrt(eta(1.2 x))
        import qumond as Q
        gf = Q.galactic_field()
        f = Q.cached_field(gf["gNe_a0"])
        xs = np.logspace(-2, 1.5, 40)
        ratio = dict(x=_round(xs), r=_round(np.sqrt(np.interp(np.log(1.2 * xs), f.lnr, f.eta))))
    return out, ratio


PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Wide binary explorer</title>
<style>
:root{--surface:%(l_surface)s;--ink:%(l_ink)s;--ink2:%(l_ink2)s;--muted:%(l_muted)s;--grid:%(l_grid)s;
--newton:%(l_newton)s;--mond:%(l_mond)s;--third:%(l_third)s;--panel:%(l_panel)s;--border:%(l_border)s}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--surface:%(d_surface)s;--ink:%(d_ink)s;
--ink2:%(d_ink2)s;--muted:%(d_muted)s;--grid:%(d_grid)s;--newton:%(d_newton)s;--mond:%(d_mond)s;--third:%(d_third)s;
--panel:%(d_panel)s;--border:%(d_border)s}}
:root[data-theme="dark"]{--surface:%(d_surface)s;--ink:%(d_ink)s;--ink2:%(d_ink2)s;--muted:%(d_muted)s;
--grid:%(d_grid)s;--newton:%(d_newton)s;--mond:%(d_mond)s;--third:%(d_third)s;--panel:%(d_panel)s;--border:%(d_border)s}
body{margin:0;background:var(--surface);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1180px;margin:0 auto;padding:16px}
h1{font-size:21px;margin:4px 0 2px} p.lead{color:var(--ink2);margin:0 0 12px;max-width:900px}
.controls{display:flex;flex-wrap:wrap;gap:10px 18px;align-items:center;background:var(--panel);
border:1px solid var(--border);border-radius:8px;padding:10px 12px;margin-bottom:10px}
.controls label{color:var(--ink2);font-size:13px;display:flex;align-items:center;gap:6px}
select,input[type=range]{accent-color:var(--newton)} output{color:var(--ink);font-variant-numeric:tabular-nums;min-width:3.2em;display:inline-block}
.stats{display:flex;flex-wrap:wrap;gap:12px;margin:8px 0 6px}
.tile{background:var(--panel);border:1px solid var(--border);border-radius:8px;padding:8px 12px;min-width:150px}
.tile .v{font-size:20px;font-weight:600;font-variant-numeric:tabular-nums}.tile .k{color:var(--ink2);font-size:12px}
.plots{display:grid;grid-template-columns:1.25fr 1fr;gap:10px}
@media (max-width:820px){.plots{grid-template-columns:1fr}}
.plot{height:420px;border:1px solid var(--border);border-radius:8px}
footer{color:var(--muted);font-size:12px;margin-top:10px}
</style></head><body><main>
<h1>Wide binaries: does gravity strengthen at low acceleration?</h1>
<p class="lead">MOND predicts that binaries wider than about one MOND radius (r_M = 7 kAU × √(M/M☉)) orbit ~19%% faster
than Newton allows, even after the Galaxy's own pull damps the effect. Pick a published sample and tighten its cuts:
the wide/tight ratio of the median scaled velocity v~ should sit near 1.2 under MOND, and near 1 under Newton if
eccentricities, hidden companions and measurement errors do not change with separation (they do, a little).</p>
<div class="controls">
<label>Sample <select id="sample"></select></label>
<label>max distance <input type="range" id="dmax" min="50" max="250" step="10" value="250"><output id="dmax_o"></output>pc</label>
<label>max σ(v~) <input type="range" id="smax" min="0.02" max="0.3" step="0.01" value="0.3"><output id="smax_o"></output></label>
<label>max RUWE <input type="range" id="rmax" min="1.0" max="3.0" step="0.05" value="3.0"><output id="rmax_o"></output></label>
<label>v~ below <select id="vcut"><option>1.5</option><option selected>2</option><option>2.5</option><option>5</option></select></label>
<label><input type="checkbox" id="rv2"> both stars have Gaia RVs</label>
</div>
<div class="stats" id="stats"></div>
<div class="plots"><div class="plot" id="p1"></div><div class="plot" id="p2"></div></div>
<footer>Data: Gaia DR3 via El-Badry, Rix &amp; Heintz (2021), Chae (2023, Zenodo 8416435) and our own pair search (project M3).
The MOND curve uses our QUMOND forward model with the Galactic external field; both curves are normalised to the data's two
innermost bins. σ(v~) is each pair's propagated Gaia error. Error bars: 1.25 σ/√N. Built offline; no network needed.</footer>
</main>
<script>%(plotly)s</script>
<script>
const DATA = %(data)s; const RATIO = %(ratio)s;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const sel = document.getElementById('sample');
for (const [k, v] of Object.entries(DATA)) { const o = document.createElement('option'); o.value = k; o.textContent = v.label; sel.appendChild(o); }
const ids = ['dmax','smax','rmax','vcut','rv2'];
function median(a){ if(!a.length) return NaN; const b=[...a].sort((x,y)=>x-y); const n=b.length; return n%%2?b[(n-1)/2]:0.5*(b[n/2-1]+b[n/2]); }
function interpRatio(x){ const xs=RATIO.x, rs=RATIO.r; if(x<=xs[0]) return rs[0]; if(x>=xs[xs.length-1]) return rs[rs.length-1];
  for(let i=1;i<xs.length;i++){ if(x<xs[i]){ const t=(Math.log(x)-Math.log(xs[i-1]))/(Math.log(xs[i])-Math.log(xs[i-1])); return rs[i-1]+t*(rs[i]-rs[i-1]); } } }
function update(){
  const S = DATA[sel.value]; const dmax=+dmax_i.value, smax=+smax_i.value, rmax=+rmax_i.value, vc=+vcut_i.value, need=rv2_i.checked;
  dmax_o.textContent=dmax; smax_o.textContent=smax.toFixed(2); rmax_o.textContent=rmax.toFixed(2);
  const x=[], v=[];
  for(let i=0;i<S.x.length;i++){ if(S.vt[i]===null||S.x[i]===null) continue; if(S.d[i]>dmax||S.sig[i]>smax||(S.ruwe[i]!==null&&S.ruwe[i]>rmax)||S.vt[i]>=vc||(need&&!S.rv2[i])) continue; x.push(S.x[i]); v.push(S.vt[i]); }
  const order = x.map((_,i)=>i).sort((a,b)=>x[a]-x[b]); const nb = Math.max(3, Math.min(10, Math.floor(x.length/150)));
  const bx=[], bm=[], be=[];
  for(let k=0;k<nb;k++){ const lo=Math.floor(k*order.length/nb), hi=Math.floor((k+1)*order.length/nb); const idx=order.slice(lo,hi); if(idx.length<10) continue;
    const vv=idx.map(i=>v[i]); const xx=idx.map(i=>x[i]); const m=median(vv); const sd=Math.sqrt(vv.reduce((s,a)=>s+(a-m)*(a-m),0)/vv.length);
    bx.push(Math.exp(xx.reduce((s,a)=>s+Math.log(a),0)/xx.length)); bm.push(m); be.push(1.2533*sd/Math.sqrt(vv.length)); }
  const norm = bm.length>=2 ? 0.5*(bm[0]+bm[1])/(0.5*(interpRatio(bx[0])+interpRatio(bx[1]))) : NaN;
  const gx=[]; for(let i=0;i<=60;i++) gx.push(Math.pow(10,-1.7+i*2.6/60));
  const tight=[], wide=[]; for(let i=0;i<x.length;i++){ if(x[i]<0.3) tight.push(v[i]); else if(x[i]>1) wide.push(v[i]); }
  const mt=median(tight), mw=median(wide); const et=1.2533*Math.sqrt(tight.reduce((s,a)=>s+(a-mt)*(a-mt),0)/Math.max(tight.length,1))/Math.sqrt(Math.max(tight.length,1));
  const ew=1.2533*Math.sqrt(wide.reduce((s,a)=>s+(a-mw)*(a-mw),0)/Math.max(wide.length,1))/Math.sqrt(Math.max(wide.length,1));
  const r=mw/mt, er=r*Math.sqrt((et/mt)**2+(ew/mw)**2);
  const rM = interpRatio(3)/interpRatio(0.15);
  stats.innerHTML = `<div class="tile"><div class="v">${x.length.toLocaleString()}</div><div class="k">pairs pass the cuts</div></div>`+
   `<div class="tile"><div class="v">${isFinite(mt)?mt.toFixed(3):'–'}</div><div class="k">median v~, tight (r_sky &lt; 0.3 r_M; ${tight.length})</div></div>`+
   `<div class="tile"><div class="v">${isFinite(mw)?mw.toFixed(3):'–'}</div><div class="k">median v~, wide (r_sky &gt; r_M; ${wide.length})</div></div>`+
   `<div class="tile"><div class="v">${isFinite(r)?r.toFixed(3)+' ± '+er.toFixed(3):'–'}</div><div class="k">wide / tight (Newton ≈ 1, MOND ≈ ${rM.toFixed(2)})</div></div>`;
  const ink=css('--ink'), ink2=css('--ink2'), grid=css('--grid'), bg=css('--surface');
  const lay = (t,xa,ya,xlog)=>({title:{text:t,font:{size:14,color:ink},x:0.02}, paper_bgcolor:bg, plot_bgcolor:bg, margin:{l:56,r:14,t:44,b:48},
     font:{color:ink2,size:12}, xaxis:Object.assign({title:{text:xa},type:xlog?'log':'linear',gridcolor:grid,zeroline:false},
       xlog?{tickvals:[0.02,0.05,0.1,0.2,0.5,1,2,5,10],ticktext:['0.02','0.05','0.1','0.2','0.5','1','2','5','10']}:{}),
     yaxis:{title:{text:ya},gridcolor:grid,zeroline:false},
     legend:{orientation:'h',y:-0.22,x:0}, hovermode:'closest'});
  Plotly.react('p1', [
     {x:gx, y:gx.map(_=>norm*interpRatio(bx.length?bx[0]:0.1)), mode:'lines', name:'Newton (flat)', line:{color:css('--newton'),width:2}, hoverinfo:'skip'},
     {x:gx, y:gx.map(g=>norm*interpRatio(g)), mode:'lines', name:'MOND with the Galactic field', line:{color:css('--mond'),width:2,dash:'dash'}, hoverinfo:'skip'},
     {x:bx, y:bm, error_y:{type:'data',array:be,color:ink,thickness:1.2,width:0}, mode:'markers', name:'data: median v~ per bin', marker:{color:ink,size:9,line:{color:bg,width:2}},
      hovertemplate:'r_sky/r_M %%{x:.2f}<br>median v~ %%{y:.3f}<extra></extra>'}],
     lay('Median scaled velocity vs separation','r_sky / r_M (projected separation / MOND radius)','median v~',true), {displayModeBar:false, responsive:true});
  const step = Math.min(0.08, vc/20), nbh = Math.round(vc/step);
  function hist(a){ const h=new Array(nbh).fill(0); for(const t of a){ const i=Math.floor(t/step); if(i>=0&&i<nbh) h[i]++; }
    const n=Math.max(a.length,1); return h.map(c=>c/n/step); }
  const edges=[...Array(nbh+1).keys()].map(i=>i*step);
  const ht=hist(tight), hw=hist(wide);
  Plotly.react('p2', [
     {x:edges, y:[...ht, ht[ht.length-1]], mode:'lines', line:{shape:'hv',color:css('--newton'),width:2}, name:`tight: r_sky < 0.3 r_M (${tight.length})`,
      hovertemplate:'v~ %%{x:.2f}: %%{y:.2f}<extra>tight</extra>'},
     {x:edges, y:[...hw, hw[hw.length-1]], mode:'lines', line:{shape:'hv',color:css('--mond'),width:2}, name:`wide: r_sky > r_M (${wide.length})`,
      hovertemplate:'v~ %%{x:.2f}: %%{y:.2f}<extra>wide</extra>'}],
     lay('v~ distributions: tight vs wide pairs','v~','fraction per unit v~',false), {displayModeBar:false, responsive:true});
}
const dmax_i=document.getElementById('dmax'), smax_i=document.getElementById('smax'), rmax_i=document.getElementById('rmax'), vcut_i=document.getElementById('vcut'), rv2_i=document.getElementById('rv2');
const dmax_o=document.getElementById('dmax_o'), smax_o=document.getElementById('smax_o'), rmax_o=document.getElementById('rmax_o'), stats=document.getElementById('stats');
for (const el of [sel, dmax_i, smax_i, rmax_i, vcut_i, rv2_i]) el.addEventListener('input', update);
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', update);
update();
</script></body></html>
"""


def build(out=None):
    import plotly.offline
    data, ratio = payload()
    js = plotly.offline.get_plotlyjs()
    fill = {f"l_{k}": v for k, v in THEMES["light"].items()}
    fill.update({f"d_{k}": v for k, v in THEMES["dark"].items()})
    html = PAGE % dict(plotly=js, data=json.dumps(data, separators=(",", ":")),
                       ratio=json.dumps(ratio, separators=(",", ":")), **fill)
    out = FIG_DIR / "wide_binary_explorer.html" if out is None else out
    out.write_text(html)
    return out


if __name__ == "__main__":
    print(build())
