"""An offline page for exploring the supernova census: figures/supernova_explorer.html.

* **The zoo.** Every supernova of the statistical sample in the plane of peak
  luminosity and timescale, by family (click a family to hide it). Click a
  point for its ZTF light curve (where we fetched one), its host galaxy's mass
  and star-formation rate, and, for Type II, the fitted end of the plateau.
* **Survey simulator.** Drag the survey's magnitude limit and watch Type Ia
  supernovae take over, next to what really explodes per unit volume.
* **Missing stars.** Smartt's counting argument as a slider: how many of the
  24 progenitors should have been born above 18 solar masses, if red
  supergiants up to a given mass exploded.

Self-contained: no network, no libraries; the data are embedded as JSON and
drawn on canvas and SVG; light and dark themes follow the system setting (or
the toggle).
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import bts
import config as C
import lightcurves as LC
import progenitors as P

OUT = C.FIGURE_DIR / "supernova_explorer.html"
FAM_CODE = {f: i for i, f in enumerate(C.FAMILY_ORDER)}


def _r(x, nd):
    return None if x is None or not np.isfinite(x) else round(float(x), nd)


def build_payload(res: dict, df: pd.DataFrame) -> dict:
    stat = df[bts.statistical_sample(df) & df["family"].notna() & df["M_peak"].notna()].copy()
    hosts = pd.read_csv(C.DERIVED_DIR / "hosts.csv", dtype={"objID": str, "specObjID": str})
    hosts = hosts.set_index("ztf_id")
    lct = res["lightcurves"]["table"].set_index("ztf_id")
    have_lc = {p.stem for p in C.ALERCE_DIR.glob("*.json")}
    pts, lcs = [], {}
    for _, r in stat.iterrows():
        oid = r["ztf_id"]
        dur = r["duration"] / (1 + r["z_cmb"]) if not r["duration_is_limit"] else np.nan
        h = hosts.loc[oid] if oid in hosts.index else None
        p = {"id": oid, "n": r["iau_id"] if r["iau_id"] != "-" else oid, "t": r["tns_type"],
             "f": FAM_CODE[r["family"]], "M": _r(r["M_peak"], 2), "d": _r(dur, 1),
             "z": _r(r["redshift"], 4), "m": _r(r["peakmag"], 2),
             "ds": "CF4" if r["mu_source"] == "Cosmicflows-4 group" else "z",
             "y": str(r["peak_date"].date())}
        if h is not None and bool(h["has_mpa"]):
            p["hm"], p["hs"] = _r(h["lgm_tot_p50"], 2), _r(h["sfr_tot_p50"], 2)
        if oid in lct.index and bool(lct.loc[oid, "plateau_ok"]):
            row = lct.loc[oid]
            p["tp"] = _r(row["fit_t_pt"], 1)
            p["E"], p["Mej"] = _r(row["E51"], 2), _r(row["M_ej"], 1)
        if oid in have_lc:
            got = LC.load_lc(oid)
            if got is not None:
                nl = LC.nightly(got[0])
                z = r["z_cmb"]
                ph = (nl["mjd"] + 2400000.5 - r["peak_jd"]) / (1 + z)
                M = nl["mag"] - r["mu"] - np.where(nl["fid"] == 1, C.A_G_PER_AV, C.A_R_PER_AV) * r["A_V"] \
                    + 2.5 * np.log10(1 + z)
                keep = (ph > -60) & (ph < 400) & np.isfinite(M)
                lcs[oid] = [[int(f), round(float(a), 1), round(float(b), 2)]
                            for f, a, b in zip(nl["fid"][keep], ph[keep], M[keep])]
                p["lc"] = 1
        pts.append(p)
    # survey simulator: magnitude-limited shares vs depth (already computed)
    ml = res["rates"]["mag_limited"]
    sim = {"m": ml["m_lim"].round(2).tolist(),
           "share": {f: ml[f].round(4).tolist() for f in C.FAMILY_ORDER}}
    fr = res["rates"]["main"]["fractions"]["M<-14.0"]
    vol = {g: round(fr[g]["frac"], 4) for g in ("Ia", "Ibc", "II")}
    # missing-stars slider
    pr = res["progenitors"]
    n_prog = int(pr["p_none_above18"]["N"])
    mm = np.arange(18.0, 60.01, 0.5)
    expect = [round(n_prog * P.salpeter_number_fraction(18.0, m, 8.0, m), 3) for m in mm]
    rates = res["rates"]["main"]["rates"]
    table = [{"family": C.FAMILY_LABEL[f], "n": int((stat["family"] == f).sum())}
             for f in C.FAMILY_ORDER]
    summary = {
        "rate_Ia": rates["Ia"]["rate"], "rate_Ia_lo": rates["Ia"]["lo"], "rate_Ia_hi": rates["Ia"]["hi"],
        "rate_CC": rates["CC"]["rate"], "rate_CC_lo": rates["CC"]["lo"], "rate_CC_hi": rates["CC"]["hi"],
        "ratio": res["rates"]["main"]["ratio_cc_ia"]["ratio"],
        "p_nothing_brighter": pr["nothing_brighter"]["p"],
        "counting_expected": pr["counting"]["n_expected_above"],
    }
    return {"pts": pts, "lcs": lcs, "sim": sim, "vol": vol, "li_vol": {"Ia": 0.241, "Ibc": 0.187, "II": 0.572},
            "fam": [{"key": f, "label": C.FAMILY_LABEL[f]} for f in C.FAMILY_ORDER],
            "miss": {"N": n_prog, "M": mm.round(1).tolist(), "expect": expect},
            "table": table, "summary": summary}


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Supernova census explorer</title>
<style>
:root{color-scheme:light;--bg:#fcfcfb;--panel:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;
--grid:#e1e0d9;--axis:#c3c2b7;--cloud:#d6d4cc;--sel:#0b0b0b;
--f0:#2a78d6;--f1:#eb6834;--f2:#1baf7a;--f3:#eda100;--f4:#e87ba4;--f5:#008300;--hi:#e34948}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#1a1a19;
--panel:#222220;--ink:#ffffff;--ink2:#c3c2b7;--muted:#9a988f;--grid:#33332f;--axis:#4a4944;--cloud:#45443f;
--sel:#ffffff;--f0:#3987e5;--f1:#d95926;--f2:#199e70;--f3:#c98500;--f4:#d55181;--f5:#008300;--hi:#e66767}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#1a1a19;--panel:#222220;--ink:#ffffff;--ink2:#c3c2b7;
--muted:#9a988f;--grid:#33332f;--axis:#4a4944;--cloud:#45443f;--sel:#ffffff;--f0:#3987e5;--f1:#d95926;
--f2:#199e70;--f3:#c98500;--f4:#d55181;--f5:#008300;--hi:#e66767}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1180px;margin:0 auto;padding:16px}
h1{font-size:22px;margin:4px 0 2px}h2{font-size:16px;margin:0 0 6px}
.sub{color:var(--ink2);margin:0 0 14px}.grid{display:grid;grid-template-columns:1.25fr 1fr;gap:14px}
.card{background:var(--panel);border:1px solid var(--grid);border-radius:10px;padding:12px;min-width:0}
.chips{display:flex;flex-wrap:wrap;gap:6px;margin:4px 0 8px}
.chip{border:1px solid var(--axis);background:none;color:var(--ink);border-radius:999px;padding:3px 10px;
font:inherit;font-size:12px;cursor:pointer;display:flex;align-items:center;gap:6px}
.chip[aria-pressed="false"]{opacity:.45}.dot{width:10px;height:10px;border-radius:50%;display:inline-block}
canvas{width:100%;display:block;touch-action:manipulation}.note{color:var(--muted);font-size:12px}
.kv{display:grid;grid-template-columns:auto 1fr;gap:2px 10px;font-size:13px}.kv b{font-weight:600}
table{border-collapse:collapse;width:100%;font-size:13px}td,th{padding:4px 6px;border-bottom:1px solid var(--grid);
text-align:left}th{color:var(--ink2);font-weight:600}input[type=range]{width:100%}
.big{font-size:22px;font-weight:700}#tip{position:fixed;pointer-events:none;background:var(--panel);
border:1px solid var(--axis);border-radius:6px;padding:4px 8px;font-size:12px;display:none;z-index:5}
.top{display:flex;justify-content:space-between;align-items:flex-start;gap:12px}
.toggle{border:1px solid var(--axis);background:none;color:var(--ink);border-radius:8px;padding:4px 10px;
font:inherit;font-size:12px;cursor:pointer}
@media (max-width:860px){.grid{grid-template-columns:1fr}}
</style></head><body><main>
<div class="top"><div><h1>Supernova census explorer</h1>
<p class="sub">Seven years of the ZTF Bright Transient Survey: every classified supernova that peaked brighter
than 18.5 mag. Project L2 of the Astrophysics Project Atlas (all data public; see the project README).</p></div>
<button class="toggle" id="theme" type="button">Theme</button></div>
<div class="grid">
<section class="card"><h2>The zoo: how bright, and for how long</h2>
<div class="chips" id="chips"></div>
<canvas id="zoo" height="420" aria-label="Peak absolute magnitude against duration for every supernova"></canvas>
<p class="note">Duration = time above half of peak brightness (rest frame). Click a point to inspect it;
click a family to show or hide it. Absolute magnitudes are not corrected for dust in the host galaxy.</p></section>
<section class="card"><h2 id="selname">Click a supernova</h2>
<canvas id="lc" height="260" aria-label="Light curve of the selected supernova"></canvas>
<div class="kv" id="info"><span class="note">Light curves (ZTF g and r alert photometry via ALeRCE) are
available for every Type II and the best-observed few of each other class.</span></div></section>
<section class="card"><h2>Survey simulator: magnitude-limited vs per volume</h2>
<label for="mlim">Survey depth: peak brighter than <b id="mlv">18.5</b> mag</label>
<input type="range" id="mlim" min="0" max="1" step="1" value="0">
<svg id="bars" viewBox="0 0 520 170" width="100%" role="img" aria-label="Type shares"></svg>
<p class="note">Top: what the survey finds at that depth (classified events, weighted for classification
completeness). Below: per unit volume (1/V<sub>max</sub>, brighter than M = −14), and the Lick Observatory
Supernova Search's volume-limited census (Li et al. 2011). Ib/c includes Ic-BL, Ibn, Icn; II includes IIb, IIn.</p></section>
<section class="card"><h2>Missing stars: the red supergiant problem</h2>
<label for="mmax">If red supergiants up to <b id="mmv">25</b> M☉ exploded as Type II-P...</label>
<input type="range" id="mmax" min="0" max="84" step="1" value="14">
<p><span class="big" id="nexp">–</span> of the <span id="nprog"></span> progenitors in the sample would have been born
above 18 M☉, and the chance of finding none would be <span class="big" id="pnone">–</span>.</p>
<p class="note">Salpeter IMF between 8 M☉ and the chosen mass; the sample is Davies &amp; Beasor (2018)'s 24
Type II-P/L supernovae with progenitor detections or limits. None was more massive than ~18–20 M☉. With
measurement errors and luminosities instead of masses the tension is milder: a sample like this shows nothing
brighter than the brightest seen <b id="pnb"></b> of the time (see Fig. 10).</p>
<h2 style="margin-top:14px">The numbers</h2><table id="tab"></table></section>
</div></main><div id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('data').textContent);
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const famColor=i=>css('--f'+i);const shown=D.fam.map(()=>true);let sel=null;
const root=document.documentElement;
document.getElementById('theme').onclick=()=>{const dark=matchMedia('(prefers-color-scheme: dark)').matches;
const cur=root.dataset.theme||(dark?'dark':'light');root.dataset.theme=cur==='dark'?'light':'dark';drawAll();};
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',drawAll);
const chips=document.getElementById('chips');
D.fam.forEach((f,i)=>{const b=document.createElement('button');b.className='chip';b.type='button';
b.setAttribute('aria-pressed','true');b.innerHTML=`<span class="dot" style="background:var(--f${i})"></span>${f.label}`;
b.onclick=()=>{shown[i]=!shown[i];b.setAttribute('aria-pressed',shown[i]);drawZoo();};chips.appendChild(b);});
function setup(cv){const r=cv.getBoundingClientRect();const dpr=window.devicePixelRatio||1;
const h=+cv.getAttribute('height');cv.width=r.width*dpr;cv.height=h*dpr;cv.style.height=h+'px';
const c=cv.getContext('2d');c.setTransform(dpr,0,0,dpr,0,0);return [c,r.width,h];}
// ---- zoo ----
const zoo=document.getElementById('zoo');let zmap=null;
const X0=Math.log10(4),X1=Math.log10(400),Y0=-23.2,Y1=-12.5;  // brighter at the top
function drawZoo(){const [c,w,h]=setup(zoo);const L=46,R=10,T=8,B=34;
const sx=x=>L+(Math.log10(x)-X0)/(X1-X0)*(w-L-R),sy=y=>T+(y-Y0)/(Y1-Y0)*(h-T-B);
c.fillStyle=css('--panel');c.fillRect(0,0,w,h);c.strokeStyle=css('--grid');c.lineWidth=1;
c.fillStyle=css('--ink2');c.font='11px system-ui';
[5,10,20,50,100,200].forEach(v=>{const x=sx(v);c.beginPath();c.moveTo(x,T);c.lineTo(x,h-B);c.stroke();
c.textAlign='center';c.fillText(v,x,h-B+14);});
for(let m=-13;m>=-23;m-=1){const y=sy(m);c.beginPath();c.moveTo(L,y);c.lineTo(w-R,y);c.stroke();
c.textAlign='right';c.fillText(m,L-6,y+4);}
c.textAlign='center';c.fillText('days above half brightness (rest frame)',L+(w-L-R)/2,h-4);
c.save();c.translate(12,T+(h-T-B)/2);c.rotate(-Math.PI/2);c.fillText('peak absolute magnitude',0,0);c.restore();
zmap=[];const order=[...D.pts.keys()].sort((a,b)=>D.pts[b].f-D.pts[a].f);
for(const k of order){const p=D.pts[k];if(!shown[p.f]||p.d==null||p.M==null)continue;
const x=sx(p.d),y=sy(p.M);if(x<L||x>w-R||y<T||y>h-B)continue;
c.fillStyle=famColor(p.f);c.beginPath();c.arc(x,y,p.lc?2.6:2.0,0,6.283);c.fill();zmap.push([x,y,k]);}
if(sel!=null){const p=D.pts[sel];if(p.d!=null){c.strokeStyle=css('--sel');c.lineWidth=2;c.beginPath();
c.arc(sx(p.d),sy(p.M),7,0,6.283);c.stroke();}}}
function nearest(ev){const r=zoo.getBoundingClientRect();const x=ev.clientX-r.left,y=ev.clientY-r.top;
let best=null,bd=144;for(const [px,py,k] of zmap){const d=(px-x)**2+(py-y)**2;if(d<bd){bd=d;best=k;}}return best;}
const tip=document.getElementById('tip');
zoo.addEventListener('mousemove',ev=>{const k=nearest(ev);if(k==null){tip.style.display='none';return;}
const p=D.pts[k];tip.innerHTML=`<b>${p.n}</b> ${p.t}<br>M = ${p.M}, z = ${p.z}${p.lc?' · light curve':''}`;
tip.style.left=(ev.clientX+12)+'px';tip.style.top=(ev.clientY+12)+'px';tip.style.display='block';});
zoo.addEventListener('mouseleave',()=>tip.style.display='none');
zoo.addEventListener('click',ev=>{const k=nearest(ev);if(k!=null){sel=k;drawZoo();drawLC();}});
// ---- light curve ----
function drawLC(){const cv=document.getElementById('lc');const [c,w,h]=setup(cv);c.fillStyle=css('--panel');
c.fillRect(0,0,w,h);const info=document.getElementById('info');if(sel==null)return;const p=D.pts[sel];
document.getElementById('selname').textContent=`${p.n} — ${p.t}`;
const rows=[['ZTF name',p.id],['peak',`${p.m} mag on ${p.y}`],['redshift',p.z],
['peak absolute mag',`${p.M} (${p.ds==='CF4'?'Cosmicflows-4 distance':'redshift distance'})`]];
if(p.hm!=null)rows.push(['host galaxy',`log M★ = ${p.hm}, log SFR = ${p.hs} (MPA-JHU)`]);
if(p.tp!=null)rows.push(['plateau',`ends at day ${p.tp}; E ≈ ${p.E}×10⁵¹ erg, M_ej ≈ ${p.Mej} M☉ (R★ = 500 R☉)`]);
info.innerHTML=rows.map(r=>`<b>${r[0]}</b><span>${r[1]}</span>`).join('');
const L=46,R=10,T=10,B=30;const pts=D.lcs[p.id];c.fillStyle=css('--ink2');c.font='11px system-ui';
if(!pts){c.textAlign='center';c.fillText('no light curve fetched for this event',w/2,h/2);return;}
let t0=Math.min(...pts.map(q=>q[1]),-20),t1=Math.max(...pts.map(q=>q[1]),60);t1=Math.min(t1,400);
let m0=Math.max(...pts.map(q=>q[2]))+0.3,m1=Math.min(...pts.map(q=>q[2]))-0.3;
const sx=x=>L+(x-t0)/(t1-t0)*(w-L-R),sy=y=>T+(y-m1)/(m0-m1)*(h-T-B);  // brighter at the top
c.strokeStyle=css('--grid');const step=(t1-t0)>200?100:25;
for(let t=Math.ceil(t0/step)*step;t<=t1;t+=step){c.beginPath();c.moveTo(sx(t),T);c.lineTo(sx(t),h-B);c.stroke();
c.textAlign='center';c.fillText(t,sx(t),h-B+14);}
for(let m=Math.ceil(m1);m<=m0;m+=1){c.beginPath();c.moveTo(L,sy(m));c.lineTo(w-R,sy(m));c.stroke();
c.textAlign='right';c.fillText(m,L-6,sy(m)+4);}
c.textAlign='center';c.fillText('rest-frame days from peak',L+(w-L-R)/2,h-2);
for(const [f,t,m] of pts){c.fillStyle=f===1?css('--f2'):css('--hi');c.beginPath();c.arc(sx(t),sy(m),2.6,0,6.283);c.fill();}
c.textAlign='right';c.fillStyle=css('--f2');c.fillText('● g',w-R-34,T+12);c.fillStyle=css('--hi');c.fillText('● r',w-R-6,T+12);}
// ---- survey simulator ----
const famShort=['Ia','II','Ib/c','IIn','SLSN','TDE'];
function drawBars(){const i=+document.getElementById('mlim').value;const m=D.sim.m[i];
document.getElementById('mlv').textContent=m.toFixed(2);const svg=document.getElementById('bars');
let s='';const W=500,x0=10;let x=x0;
D.fam.forEach((f,k)=>{const v=D.sim.share[f.key][i]||0;const wd=v*W;
s+=`<rect x="${x+1}" y="20" width="${Math.max(wd-2,0)}" height="36" fill="var(--f${k})"><title>${f.label}: ${(100*v).toFixed(1)}%</title></rect>`;
if(wd>34)s+=`<text x="${x+wd/2}" y="43" text-anchor="middle" font-size="12" fill="${k===0?'#fff':'#0b0b0b'}">${famShort[k]} ${(100*v).toFixed(0)}%</text>`;x+=wd;});
s+=`<text x="${x0}" y="14" font-size="12" fill="var(--ink2)">what the survey finds</text>`;
const rows=[['per unit volume (this work)',D.vol],['LOSS volume-limited (Li et al. 2011)',D.li_vol]];
rows.forEach((r,j)=>{let y=86+j*46;x=x0;s+=`<text x="${x0}" y="${y-6}" font-size="12" fill="var(--ink2)">${r[0]}</text>`;
[['Ia',0],['Ibc',2],['II',1]].forEach(([g,k])=>{const v=r[1][g];const wd=v*W;
s+=`<rect x="${x+1}" y="${y}" width="${wd-2}" height="26" fill="var(--f${k})"><title>${g}: ${(100*v).toFixed(1)}%</title></rect>`;
if(wd>40)s+=`<text x="${x+wd/2}" y="${y+18}" text-anchor="middle" font-size="12" fill="${k===0?'#fff':'#0b0b0b'}">${g==='Ibc'?'Ib/c':g} ${(100*v).toFixed(0)}%</text>`;x+=wd;});});
svg.innerHTML=s;}
(()=>{const el=document.getElementById('mlim');el.max=D.sim.m.length-1;
let k=D.sim.m.findIndex(v=>Math.abs(v-18.5)<1e-6);el.value=k<0?D.sim.m.length-1:k;})();
document.getElementById('mlim').addEventListener('input',drawBars);
// ---- missing stars ----
document.getElementById('nprog').textContent=D.miss.N;
document.getElementById('pnb').textContent=(100*D.summary.p_nothing_brighter).toFixed(0)+'%';
function drawMiss(){const i=+document.getElementById('mmax').value;const M=D.miss.M[i],e=D.miss.expect[i];
document.getElementById('mmv').textContent=M.toFixed(1);document.getElementById('nexp').textContent=e.toFixed(1);
const p=Math.pow(1-e/D.miss.N,D.miss.N);document.getElementById('pnone').textContent=p<0.001?p.toExponential(1):(100*p).toFixed(1)+'%';}
document.getElementById('mmax').addEventListener('input',drawMiss);
function fmt(v){return (v/1e4).toFixed(2)+'×10⁴';}
const S=D.summary;document.getElementById('tab').innerHTML=
`<tr><th>Quantity</th><th>Value</th></tr>`+
`<tr><td>Type Ia rate (M &lt; −16.5)</td><td>${fmt(S.rate_Ia)} Gpc⁻³ yr⁻¹ (95%: ${fmt(S.rate_Ia_lo)}–${fmt(S.rate_Ia_hi)})</td></tr>`+
`<tr><td>Core-collapse rate (M &lt; −14)</td><td>${fmt(S.rate_CC)} Gpc⁻³ yr⁻¹ (95%: ${fmt(S.rate_CC_lo)}–${fmt(S.rate_CC_hi)})</td></tr>`+
`<tr><td>Core collapse per Type Ia (same volume)</td><td>${S.ratio.toFixed(2)}</td></tr>`+
`<tr><td>Smartt's "expected above 18 M☉"</td><td>${S.counting_expected.toFixed(1)} (found: at most 1)</td></tr>`+
D.table.map(t=>`<tr><td>${t.family}</td><td>${t.n.toLocaleString()} in the zoo</td></tr>`).join('');
function drawAll(){drawZoo();drawLC();drawBars();drawMiss();}
window.addEventListener('resize',()=>{drawZoo();drawLC();});drawAll();
</script></body></html>
"""


def write(res: dict, df: pd.DataFrame) -> str:
    payload = build_payload(res, df)
    text = HTML.replace("__DATA__", json.dumps(payload, separators=(",", ":")).replace("</", "<\\/"))
    OUT.write_text(text)
    return str(OUT)
