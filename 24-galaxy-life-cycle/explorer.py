"""An offline page for exploring the galaxy life cycle: figures/life_cycle_explorer.html.

Three panels, all drawn from numbers this project measured:

* **The SFR-mass plane.** Galaxies per unit volume (1/V_max weighted), with
  a switch between all environments, the emptiest quarter and the most
  crowded quarter, and an overlay of the fraction of galaxies whose gas is lit
  by an AGN. A mass slider shows the distribution of distance below the main
  sequence at that mass and the fraction quenched.
* **The cosmic clock.** Drag through cosmic time: the star-formation rate of
  the universe then (Madau & Dickinson's points, our fit), and how much of
  today's stellar mass already existed.
* **The gallery.** The rule-picked SDSS images of Fig. 09 (typical colour per cell).

The page is self-contained (no network, no libraries): data embedded as JSON,
drawn as SVG, light and dark themes.
"""

from __future__ import annotations

import base64
import json

import numpy as np

import galaxies as gx
import sfh
from config import FIGURE_DIR, IMF_TO_SALPETER_SFR

X_EDGES = np.round(np.arange(8.5, 12.01, 0.1), 2)
Y_EDGES = np.round(np.arange(-3.5, 1.81, 0.1), 2)


def _grid(logm, logsfr, w):
    h, _, _ = np.histogram2d(logm, logsfr, bins=[X_EDGES, Y_EDGES], weights=w)
    return h


def _norm(h):
    return np.round(h / h.max(), 5)


def build_payload(res: dict) -> dict:
    g, w, ms = res["_g"], res["_w_local"], res["ms"]
    env = res["env"]
    ok = ms["ok"]
    opd = env["one_plus_delta"]
    usable = env["usable"] & ok
    q = np.nanpercentile(np.log10(opd[usable]), [25, 75])
    low = usable & (np.log10(opd) <= q[0])
    high = usable & (np.log10(opd) >= q[1])
    agn = gx.is_agn(g["bpt"])
    grids = {}
    for name, sel in (("all", ok), ("low", low), ("high", high)):
        h = _grid(g["logm"][sel], g["logsfr"][sel], w[sel])
        ha = _grid(g["logm"][sel & agn], g["logsfr"][sel & agn], w[sel & agn])
        n = _grid(g["logm"][sel], g["logsfr"][sel], np.ones(sel.sum()))
        frac = np.where(n >= 20, ha / np.where(h > 0, h, 1), -1)
        grids[name] = {"dens": _norm(h).tolist(), "agn": np.round(frac, 3).tolist()}
    # Delta_MS histograms per 0.2-dex mass bin, by environment
    dms = ms["dms"]
    mb = np.round(np.arange(9.0, 11.81, 0.2), 2)
    de = np.round(np.arange(-3.5, 1.51, 0.1), 2)
    hists = {}
    for name, sel in (("all", ok), ("low", low), ("high", high)):
        rows = []
        for a, b in zip(mb[:-1], mb[1:]):
            m = sel & (g["logm"] >= a) & (g["logm"] < b)
            hh, _ = np.histogram(dms[m], bins=de, weights=w[m])
            tot = hh.sum()
            fq = float(np.sum(w[m] * ms["quenched"][m]) / max(np.sum(w[m]), 1e-30))
            rows.append({"h": np.round(hh / max(tot, 1e-30), 5).tolist(), "fq": round(fq, 4),
                         "n": int(m.sum())})
        hists[name] = rows
    r = ms["ridge"]
    h = res["sfh"]
    pts = h["points"]
    ours, pub = h["ours"], h["published"]
    tgrid = np.round(np.linspace(0.45, float(sfh.cosmic_time_gyr(0.0)), 300), 4)
    zt = np.array([sfh.find_redshift_at_time(t) for t in tgrid])
    ho = sfh.history(ours)
    f_formed = np.interp(zt, ho["z"], ho["formed_before"] / ho["total_formed"])
    c = res["census"]
    return {
        "xEdges": X_EDGES.tolist(), "yEdges": Y_EDGES.tolist(), "grids": grids,
        "ms": {"slope": r["slope"], "intercept": r["intercept"],
               "valley": ms["valley"]["position"]},
        "rp": {"slope": 0.76, "intercept": -7.64},
        "massBins": mb.tolist(), "dmsEdges": de.tolist(), "hists": hists,
        "envQuartiles": [float(q[0]), float(q[1])],
        "clock": {"t": tgrid.tolist(), "z": np.round(zt, 4).tolist(),
                  "psi": np.round(sfh.psi_md(zt, ours["A"], ours["B"], ours["C"], ours["D"]),
                                  5).tolist(),
                  "psiPub": np.round(sfh.psi_md(zt, pub["A"], pub["B"], pub["C"], pub["D"]),
                                     5).tolist(),
                  "formed": np.round(f_formed, 4).tolist()},
        "points": [{"t": float(sfh.cosmic_time_gyr(z)), "z": float(z), "psi": float(10 ** y),
                    "band": str(b), "ref": str(ref)}
                   for z, y, b, ref in zip(pts["z"], pts["logpsi"], pts["band"],
                                           pts["reference"])],
        "sdss": {"t": float(sfh.cosmic_time_gyr(0.1)),
                 "psi": c["rho_sfr_salpeter"],
                 "err": c["rho_sfr_err"] * IMF_TO_SALPETER_SFR["kroupa"]},
        "n": int(ms["n"]),
    }


def gallery_payload(res: dict) -> list:
    import figures

    out = []
    for p in figures.gallery_selection(res):
        if p is None:
            out.append(None)
            continue
        figures.fetch_cutout(p)
        from config import CUTOUT_DIR

        path = CUTOUT_DIR / f"{p['objID']}_{int(figures.GALLERY_FRAME_KPC)}kpc.jpg"
        out.append({**{k: p[k] for k in ("row", "col", "z", "logm", "logsfr", "dms", "bpt",
                                           "ra", "dec")},
                    "img": "data:image/jpeg;base64," + base64.b64encode(path.read_bytes())
                    .decode()})
    return out


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Galaxy Life Cycle</title>
<style>
:root{--bg:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;
--blue:#2a78d6;--red:#e34948;--amber:#eda100;--amberText:#8a5d00;--violet:#4a3aa7;--orange:#eb6834;
--r0:#fcfcfb;--r1:#cde2fb;--r2:#9ec5f4;--r3:#6da7ec;--r4:#3987e5;--r5:#2a78d6;--r6:#1c5cab;--r7:#104281;--r8:#0d366b}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#141413;--panel:#1c1c1a;--ink:#f2f1ec;--ink2:#c3c2b7;--muted:#9a998f;--grid:#2e2d2a;--axis:#4a4944;
--blue:#5b9bea;--red:#f06e6c;--amber:#f2b632;--amberText:#f2b632;--violet:#8f82e0;--orange:#f08a5d;
--r0:#1c1c1a;--r1:#0d2a52;--r2:#123c73;--r3:#1a4f95;--r4:#2463b5;--r5:#3a7fd1;--r6:#68a1e3;--r7:#9cc3f0;--r8:#d4e6fb}}
:root[data-theme="dark"]{--bg:#141413;--panel:#1c1c1a;--ink:#f2f1ec;--ink2:#c3c2b7;--muted:#9a998f;--grid:#2e2d2a;--axis:#4a4944;
--blue:#5b9bea;--red:#f06e6c;--amber:#f2b632;--amberText:#f2b632;--violet:#8f82e0;--orange:#f08a5d;
--r0:#1c1c1a;--r1:#0d2a52;--r2:#123c73;--r3:#1a4f95;--r4:#2463b5;--r5:#3a7fd1;--r6:#68a1e3;--r7:#9cc3f0;--r8:#d4e6fb}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1100px;margin:0 auto;padding:20px 16px 48px}
h1{font-size:22px;margin:4px 0 2px}h2{font-size:16px;margin:0 0 4px}.sub{color:var(--ink2);margin:0 0 18px}
.card{background:var(--panel);border:1px solid var(--grid);border-radius:10px;padding:16px;margin:0 0 18px}
.row{display:flex;flex-wrap:wrap;gap:16px;align-items:flex-start}.grow{flex:1 1 380px;min-width:0}
.controls{display:flex;flex-wrap:wrap;gap:8px 16px;align-items:center;margin:6px 0 10px;color:var(--ink2);font-size:13px}
.seg button{border:1px solid var(--axis);background:transparent;color:var(--ink2);padding:4px 10px;cursor:pointer;font:inherit}
.seg button:first-child{border-radius:6px 0 0 6px}.seg button:last-child{border-radius:0 6px 6px 0}
.seg button+button{border-left:none}.seg button[aria-pressed="true"]{background:var(--blue);color:#fff;border-color:var(--blue)}
svg{display:block;width:100%;height:auto}text{fill:var(--ink2);font-size:11px}.tick{fill:var(--muted);font-size:10px}
.note{color:var(--muted);font-size:12px;margin:6px 0 0}.big{font-size:22px;font-weight:600;color:var(--ink)}
.stat{min-width:150px}.stat .lab{color:var(--ink2);font-size:12px}
input[type=range]{width:100%;accent-color:var(--blue)}
#tip{position:fixed;pointer-events:none;background:var(--panel);border:1px solid var(--axis);border-radius:6px;padding:6px 8px;font-size:12px;color:var(--ink);display:none;z-index:9;box-shadow:0 2px 8px rgba(0,0,0,.15)}
.gal{display:grid;grid-template-columns:120px repeat(5,minmax(0,1fr));gap:6px;align-items:center}
.gal img{width:100%;border-radius:6px;display:block}.gal .lab{font-size:12px;font-weight:600}
.gal .cap{font-size:10px;color:var(--muted);margin-top:2px}.gal .colh{font-size:12px;text-align:center;color:var(--ink2)}
.theme{position:absolute;top:14px;right:16px}
table{border-collapse:collapse;font-size:12px}td,th{padding:2px 8px;text-align:right;border-bottom:1px solid var(--grid)}
@media (max-width:640px){.gal{grid-template-columns:70px repeat(5,minmax(0,1fr))}.gal .cap{display:none}}
</style></head><body><div id="tip"></div>
<main style="position:relative">
<div class="theme seg"><button id="thL">light</button><button id="thD">dark</button></div>
<h1>The galaxy life cycle</h1>
<p class="sub">__N__ SDSS galaxies (0.02 &lt; z &lt; 0.085) with MPA-JHU masses and star-formation rates, and the universe's
star-formation history from Madau &amp; Dickinson (2014). Everything here was measured in project J3.</p>

<section class="card"><h2>The SFR–mass plane</h2>
<div class="controls"><span>environment</span><span class="seg" id="envSeg">
<button data-v="all" aria-pressed="true">all</button><button data-v="low">emptiest quarter</button><button data-v="high">most crowded quarter</button></span>
<span>colour by</span><span class="seg" id="modeSeg"><button data-v="dens" aria-pressed="true">galaxies per volume</button><button data-v="agn">AGN fraction</button></span></div>
<div class="row"><div class="grow"><svg id="plane" viewBox="0 0 560 420" role="img" aria-label="SFR versus stellar mass"></svg></div>
<div class="grow"><div class="controls" style="margin-top:0"><span>stellar mass</span><span class="big" id="massLab"></span></div>
<input type="range" id="mass" min="0" max="13" step="1" value="7" aria-label="stellar mass bin">
<svg id="hist" viewBox="0 0 440 250" role="img" aria-label="distance from the main sequence"></svg>
<div class="row" style="gap:24px;margin-top:6px"><div class="stat"><div class="lab">quenched at this mass</div><div class="big" id="fq"></div></div>
<div class="stat"><div class="lab">galaxies in the bin</div><div class="big" id="nb"></div></div></div>
<p class="note">Distance from the main sequence = log SFR − (our ridge line). The amber band is the green valley.
Switch environment to see crowded places shift galaxies from the blue peak to the red one at every mass.</p></div></div>
</section>

<section class="card"><h2>The cosmic clock</h2>
<div class="controls"><span>age of the universe</span><span class="big" id="tLab"></span><span id="zLab"></span></div>
<input type="range" id="clock" min="0" max="299" value="120" aria-label="cosmic time">
<div class="row"><div class="grow"><svg id="sfr" viewBox="0 0 560 300" role="img" aria-label="cosmic star formation rate"></svg></div>
<div class="stat" style="flex:0 0 220px"><div class="lab">star formation then, compared with today</div><div class="big" id="ratio"></div>
<div class="lab" style="margin-top:12px">today's stars already formed</div><div class="big" id="formed"></div>
<p class="note">Points: 51 survey measurements (UV in blue, infrared in red). Orange: our fit; dashed violet: Madau &amp; Dickinson's curve.
Black star: this project's SDSS measurement.</p></div></div>
</section>

<section class="card"><h2>What it looks like</h2>
<p class="note" style="margin-top:0">In each cell of the plane, the galaxy whose colour is closest to the cell's median: picked by rule, not by eye. Every frame spans 50 kpc. Hover for details.</p>
<div class="gal" id="gal"></div></section>
<p class="note">Data: SDSS DR17 / MPA-JHU (Brinchmann et al. 2004; Kauffmann et al. 2003), GSWLC-2 (Salim et al. 2016, 2018),
Madau &amp; Dickinson (2014). Built by explorer.py; AI-generated, see the project README.</p>
</main>
<script>
const D = __DATA__; const G = __GALLERY__;
const $ = s => document.querySelector(s), NS = "http://www.w3.org/2000/svg";
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
function el(t, a, p){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);if(p)p.appendChild(e);return e;}
const tip=$("#tip");function showTip(ev,html){tip.innerHTML=html;tip.style.display="block";tip.style.left=(ev.clientX+14)+"px";tip.style.top=(ev.clientY+14)+"px";}
function hideTip(){tip.style.display="none";}
let env="all", mode="dens";
function seg(id, cb){document.querySelectorAll(id+" button").forEach(b=>b.onclick=()=>{document.querySelectorAll(id+" button").forEach(x=>x.setAttribute("aria-pressed","false"));b.setAttribute("aria-pressed","true");cb(b.dataset.v);});}
seg("#envSeg", v=>{env=v;drawPlane();drawHist();}); seg("#modeSeg", v=>{mode=v;drawPlane();});
function ramp(v){const st=["--r1","--r2","--r3","--r4","--r5","--r6","--r7","--r8"];if(v<=0)return null;const l=Math.max(0,Math.min(1,(Math.log10(v)+3)/3));return css(st[Math.min(7,Math.floor(l*8))]);}
function agnCol(f){if(f<0)return null;const a=Math.min(1,f/0.6);return `color-mix(in oklab, var(--violet) ${Math.round(15+85*a)}%, var(--panel))`;}
function axes(svg, x0,x1,y0,y1, X,Y, xt,yt, xl,yl){
 for(const v of xt){const x=X(v);el("line",{x1:x,x2:x,y1:Y(y0),y2:Y(y1),stroke:css("--grid")},svg);el("text",{x:x,y:Y(y0)+14,"text-anchor":"middle",class:"tick"},svg).textContent=v;}
 for(const v of yt){const y=Y(v);el("line",{x1:X(x0),x2:X(x1),y1:y,y2:y,stroke:css("--grid")},svg);el("text",{x:X(x0)-6,y:y+3,"text-anchor":"end",class:"tick"},svg).textContent=v;}
 el("text",{x:(X(x0)+X(x1))/2,y:Y(y0)+30,"text-anchor":"middle"},svg).textContent=xl;
 const t=el("text",{x:14,y:(Y(y0)+Y(y1))/2,"text-anchor":"middle",transform:`rotate(-90 14 ${(Y(y0)+Y(y1))/2})`},svg);t.textContent=yl;}
function drawPlane(){
 const svg=$("#plane");svg.innerHTML="";const xe=D.xEdges,ye=D.yEdges,L=50,R=548,T=10,B=380;
 const X=v=>L+(v-xe[0])/(xe[xe.length-1]-xe[0])*(R-L), Y=v=>B-(v-ye[0])/(ye[ye.length-1]-ye[0])*(B-T);
 axes(svg,xe[0],xe[xe.length-1],ye[0],ye[ye.length-1],X,Y,[9,10,11,12],[-3,-2,-1,0,1],"log stellar mass (M☉)","log SFR (M☉/yr)");
 const g=D.grids[env];const dw=X(xe[1])-X(xe[0]),dh=Y(ye[0])-Y(ye[1]);
 for(let i=0;i<xe.length-1;i++)for(let j=0;j<ye.length-1;j++){
  const v=g.dens[i][j];const c=mode==="dens"?ramp(v):(v>1e-3?agnCol(g.agn[i][j]):null);if(!c)continue;
  const r=el("rect",{x:X(xe[i]),y:Y(ye[j+1]),width:dw+0.4,height:dh+0.4,fill:c},svg);
  r.addEventListener("mousemove",ev=>showTip(ev,`log M ${xe[i]}–${xe[i+1]}, log SFR ${ye[j]}–${ye[j+1]}<br>`+(mode==="dens"?`relative density ${v.toPrecision(2)}`:(g.agn[i][j]>=0?`AGN fraction ${(100*g.agn[i][j]).toFixed(0)}%`:"too few galaxies"))));
  r.addEventListener("mouseleave",hideTip);}
 const m=D.ms, xs=[8.6,11.95];
 const line=(s,b,off,col,dash,w)=>el("line",{x1:X(xs[0]),y1:Y(s*xs[0]+b+off),x2:X(xs[1]),y2:Y(s*xs[1]+b+off),stroke:col,"stroke-width":w,"stroke-dasharray":dash||""},svg);
 line(m.slope,m.intercept,0,css("--blue"),"",2);line(m.slope,m.intercept,m.valley,css("--amber"),"5 4",1.5);
 line(D.rp.slope,D.rp.intercept,0,css("--violet"),"3 3",1.2);
 el("text",{x:X(9.0),y:Y(m.slope*9.0+m.intercept)-8,fill:css("--blue")},svg).textContent="main sequence (ours)";
 el("text",{x:X(9.0),y:Y(m.slope*9.0+m.intercept+m.valley)+14,fill:css("--amberText")},svg).textContent="valley floor";
 const mb=D.massBins,k=+$("#mass").value;
 el("rect",{x:X(mb[k]),y:T,width:X(mb[k+1])-X(mb[k]),height:B-T,fill:"none",stroke:css("--ink2"),"stroke-dasharray":"2 2"},svg);}
function drawHist(){
 const svg=$("#hist");svg.innerHTML="";const k=+$("#mass").value,row=D.hists[env][k],e=D.dmsEdges;
 $("#massLab").textContent=`10^${D.massBins[k].toFixed(1)}–10^${D.massBins[k+1].toFixed(1)} M☉`;
 $("#fq").textContent=(100*row.fq).toFixed(0)+"%";$("#nb").textContent=row.n.toLocaleString();
 const L=40,R=430,T=10,B=210,ym=Math.max(...D.hists.all.map(r=>Math.max(...r.h)),...D.hists.high.map(r=>Math.max(...r.h)),...D.hists.low.map(r=>Math.max(...r.h)));
 const X=v=>L+(v-e[0])/(e[e.length-1]-e[0])*(R-L),Y=v=>B-v/ym*(B-T);
 axes(svg,e[0],e[e.length-1],0,ym,X,Y,[-3,-2,-1,0,1],[],"distance from the main sequence (dex)","share");
 el("rect",{x:X(D.ms.valley-0.25),y:T,width:X(D.ms.valley+0.25)-X(D.ms.valley-0.25),height:B-T,fill:css("--amber"),opacity:0.22},svg);
 for(let i=0;i<row.h.length;i++){const c=(e[i]+e[i+1])/2;const col=c<D.ms.valley-0.25?css("--red"):(c<D.ms.valley+0.25?css("--amber"):css("--blue"));
  const r=el("rect",{x:X(e[i])+0.5,y:Y(row.h[i]),width:Math.max(0,X(e[i+1])-X(e[i])-1),height:B-Y(row.h[i]),fill:col,rx:1},svg);
  r.addEventListener("mousemove",ev=>showTip(ev,`${e[i].toFixed(1)} to ${e[i+1].toFixed(1)} dex: ${(100*row.h[i]).toFixed(1)}% of galaxies`));r.addEventListener("mouseleave",hideTip);}
 drawPlane();}
$("#mass").oninput=drawHist;
function drawClock(){
 const svg=$("#sfr");svg.innerHTML="";const C=D.clock,k=+$("#clock").value,L=50,R=550,T=10,B=260;
 const ymax=0.2,X=v=>L+v/14*(R-L),Y=v=>B-v/ymax*(B-T);
 axes(svg,0,14,0,ymax,X,Y,[0,2,4,6,8,10,12,14],[0,0.05,0.1,0.15,0.2],"age of the universe (billion years)","M☉ / yr / Mpc³ (Salpeter)");
 const path=(arr,col,dash,w)=>{let d="";arr.forEach((v,i)=>{d+=(i?"L":"M")+X(C.t[i]).toFixed(1)+","+Y(v).toFixed(1);});el("path",{d:d,fill:"none",stroke:col,"stroke-width":w,"stroke-dasharray":dash||""},svg);};
 let a=`M${X(C.t[0])},${Y(0)}`;for(let i=0;i<=k;i++)a+=`L${X(C.t[i]).toFixed(1)},${Y(C.psi[i]).toFixed(1)}`;a+=`L${X(C.t[k])},${Y(0)}Z`;
 el("path",{d:a,fill:css("--orange"),opacity:0.15},svg);
 path(C.psiPub,css("--violet"),"5 4",1.4);path(C.psi,css("--orange"),"",2);
 for(const p of D.points){const c=el("circle",{cx:X(p.t),cy:Y(p.psi),r:3.5,fill:p.band==="UV"?css("--blue"):css("--red"),stroke:css("--panel"),"stroke-width":1},svg);
  c.addEventListener("mousemove",ev=>showTip(ev,`${p.ref}<br>z = ${p.z.toFixed(2)} (${p.t.toFixed(2)} Gyr), ${p.band}<br>ψ = ${p.psi.toFixed(3)}`));c.addEventListener("mouseleave",hideTip);}
 const s=el("text",{x:X(D.sdss.t),y:Y(D.sdss.psi)+5,"text-anchor":"middle",style:"font-size:16px",fill:css("--ink")},svg);s.textContent="★";
 el("line",{x1:X(C.t[k]),x2:X(C.t[k]),y1:T,y2:B,stroke:css("--ink2"),"stroke-dasharray":"3 3"},svg);
 el("circle",{cx:X(C.t[k]),cy:Y(C.psi[k]),r:5,fill:css("--orange"),stroke:css("--panel"),"stroke-width":2},svg);
 $("#tLab").textContent=C.t[k].toFixed(2)+" billion years";$("#zLab").textContent=`(redshift ${C.z[k].toFixed(2)}, light emitted ${(C.t[C.t.length-1]-C.t[k]).toFixed(2)} Gyr ago)`;
 $("#ratio").textContent=(C.psi[k]/C.psi[C.psi.length-1]).toFixed(1)+"×";$("#formed").textContent=(100*C.formed[k]).toFixed(1)+"%";}
$("#clock").oninput=drawClock;
function drawGallery(){
 const g=$("#gal");g.innerHTML="";const rows=["above the sequence","on the main sequence","green valley","quenched"],rc=["--blue","--blue","--amberText","--red"],masses=[9.6,10.1,10.6,11.1,11.5];
 g.appendChild(document.createElement("div"));for(const m of masses){const d=document.createElement("div");d.className="colh";d.textContent=`10^${m} M☉`;g.appendChild(d);}
 for(let r=0;r<4;r++){const l=document.createElement("div");l.className="lab";l.style.color=`var(${rc[r]})`;l.textContent=rows[r];g.appendChild(l);
  for(let c=0;c<5;c++){const p=G[r*5+c],d=document.createElement("div");if(p){const im=document.createElement("img");im.src=p.img;im.alt=`SDSS image, ${rows[r]}, log M ${p.logm.toFixed(2)}`;
   im.addEventListener("mousemove",ev=>showTip(ev,`RA ${p.ra.toFixed(4)}, Dec ${p.dec.toFixed(4)}<br>z = ${p.z.toFixed(4)}, log M = ${p.logm.toFixed(2)}<br>SFR = ${Math.pow(10,p.logsfr).toPrecision(2)} M☉/yr (${p.dms>=0?"+":""}${p.dms.toFixed(2)} dex)<br>lines: ${p.bpt}`));
   im.addEventListener("mouseleave",hideTip);d.appendChild(im);const cap=document.createElement("div");cap.className="cap";cap.textContent=`z=${p.z.toFixed(3)} · SFR ${Math.pow(10,p.logsfr).toPrecision(2)}`;d.appendChild(cap);}
   else{d.className="cap";d.textContent="none in cell";}g.appendChild(d);}}}
function theme(t){document.documentElement.setAttribute("data-theme",t);try{localStorage.setItem("j3theme",t);}catch(e){}drawAll();}
$("#thL").onclick=()=>theme("light");$("#thD").onclick=()=>theme("dark");
function drawAll(){drawHist();drawClock();}
try{const t=localStorage.getItem("j3theme");if(t)document.documentElement.setAttribute("data-theme",t);}catch(e){}
drawGallery();drawAll();
</script></body></html>
"""


def write(res: dict) -> str:
    data = build_payload(res)
    gal = gallery_payload(res)
    html = (PAGE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
            .replace("__GALLERY__", json.dumps(gal, separators=(",", ":")))
            .replace("__N__", f"{data['n']:,}"))
    path = FIGURE_DIR / "life_cycle_explorer.html"
    path.write_text(html)
    return str(path)
