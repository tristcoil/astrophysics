"""An offline page for exploring the fresh supernovae: figures/fresh_supernova_explorer.html.

* **The Hubble residuals** of every ZTF DR2 supernova we fitted, against D1's
  universe, with D1's own Pantheon+ supernovae behind them and the newest
  classifications on top.
* **Standardise them yourself.** Switch the corrections on one by one (none,
  stretch, stretch + colour, D1's full recipe) and watch the scatter fall.
* **Click a supernova** for its ZTF light curve and SALT3 fit, its redshift
  source, sub-type and host galaxy.

Self-contained: no network, no libraries; data embedded as JSON, drawn on
canvas; light and dark themes follow the system setting (or the toggle).
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

import config
import saltfit

OUT = config.FIGURE_DIR / "fresh_supernova_explorer.html"
BANDS = ("ztfg", "ztfr", "ztfi")
#: SALT3 is defined for -20..+50 rest-frame days (zero outside); stay inside.
MODEL_PHASES = np.arange(-15, 47.6, 2.5)


def _r(x, nd):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return round(x, nd) if np.isfinite(x) else None


def _lc_payload(lc: pd.DataFrame, row, kind="salt3") -> dict:
    """Light-curve points and model curves, normalised to the model's peak."""
    z, t0 = float(row["redshift"]), float(row["t0"])
    ph = (lc["mjd"] - t0) / (1 + z)
    lc = lc[lc["good"] & (ph > -20) & (ph < 60)]
    ph = ph[lc.index]
    model = saltfit.get_model(kind)
    model.set(z=z, t0=t0, x0=row["x0"], x1=row["x1"], c=row["c"], mwebv=row["mwebv"])
    times = t0 + MODEL_PHASES * (1 + z)
    curves, peak = {}, 0.0
    for b in BANDS:
        if (lc["band"] == b).any():
            curves[b] = model.bandflux(b, times, zp=30, zpsys="ab")
            peak = max(peak, float(curves[b].max()))
    if peak <= 0:
        return {}
    pts = []
    for (p, b, f, e) in zip(ph, lc["band"], lc["flux"], lc["fluxerr"]):
        pts += [int(round(p * 10)), BANDS.index(b), int(round(f / peak * 1000)), int(round(e / peak * 1000))]
    return {"p": pts, "m": {BANDS.index(b): [int(round(v / peak * 1000)) for v in c]
                            for b, c in curves.items()}}


def build_payload(R: dict) -> dict:
    from lightcurves import read_dr2_lightcurve

    t = pd.read_parquet(config.CACHE_DIR / "hubble_table.parquet")
    p = pd.read_parquet(config.CACHE_DIR / "pantheon_table.parquet")
    zf = R["ztf_standardisation"]["with_mass_step"]
    a, b = zf["params"]["alpha"], zf["params"]["beta"]
    keep = t["ok"] & t["good_sampling"].fillna(False).astype(bool) & (t["c_err"] < 0.15) \
        & (t["x1_err"] < 1.5) & t["fresh"] & t["redshift"].between(0.005, 0.3)
    t = t[keep]
    cos = t["cosmo"]
    # Residuals at each stage, each centred on the cosmology sample's median.
    raw = t["resid_raw"]
    stages = {"raw": raw, "x1": raw + a * t["x1"], "x1c": raw + a * t["x1"] - b * t["c"]}
    stages = {k: v - np.median(v[cos]) for k, v in stages.items()}
    sn, lcs = [], {}
    for name, r in t.iterrows():
        sn.append({
            "id": name, "n": r["iau_name"] if isinstance(r["iau_name"], str) else "",
            "z": _r(r["redshift"], 4), "x1": _r(r["x1"], 2), "c": _r(r["c"], 3), "mB": _r(r["mB"], 3),
            "r0": _r(stages["raw"][name], 3), "r1": _r(stages["x1"][name], 3),
            "r2": _r(stages["x1c"][name], 3), "rz": _r(r["resid_ztf"], 3),
            "rd": _r(r["resid_d1_full"], 3) if r["d1_test"] else None,
            "cos": bool(r["cosmo"]), "d1": bool(r["d1_test"]),
            "st": r["sub_type"] if isinstance(r["sub_type"], str) else "",
            "pec": bool(r["peculiar"]), "zs": r["z_group"],
            "hm": _r(r["host_mass"], 2), "chi": _r(r["chisq"] / max(r["ndof"], 1), 2),
            "np": int(r["lc_n_phases"]) if np.isfinite(r["lc_n_phases"]) else None,
        })
        lc, _ = read_dr2_lightcurve(name)
        lcs[name] = _lc_payload(lc, r)
    pm = (p["zHD"] > 0.01) & (p["IS_CALIBRATOR"] == 0) & (p["zHD"] < 0.3)
    pan = [[_r(z, 4), _r(rr, 3)] for z, rr in zip(p.loc[pm, "zHD"], p.loc[pm, "resid"])]

    live = []
    lv = R.get("live", {})
    if lv.get("available"):
        lt = pd.read_parquet(config.CACHE_DIR / "live_table.parquet")
        llc = pd.read_parquet(config.CACHE_DIR / "live_lightcurves.parquet")
        for name, r in lt[lt["d1_test"]].iterrows():
            live.append({"id": name, "oid": r["oid"], "z": _r(r["redshift"], 4), "x1": _r(r["x1"], 2),
                         "c": _r(r["c"], 3), "rd": _r(r["resid_d1_full"], 3),
                         "disc": str(r["discovery_date"])[:10], "age": _r(r["days_since_peak"], 0),
                         "zs": r["z_group"]})
            lcs[name] = _lc_payload(llc[llc["tns_name"] == name], r)
    off = R["landing"]["full"]["offset_hubble_flow"]
    meta = {"d1": R["d1_universe"], "alpha": a, "beta": b, "d1_alpha": R["d1_recipe"]["alpha"],
            "d1_beta": R["d1_recipe"]["beta"], "offset": off["mean"], "offset_err": off["err"],
            "snapshot": lv.get("meta", {}).get("snapshot"), "phases": MODEL_PHASES.tolist(),
            "n_cos": int(cos.sum()), "n_d1": int(t["d1_test"].sum())}
    return {"sn": sn, "pan": pan, "live": live, "lc": lcs, "meta": meta}


HTML = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Fresh Supernovae Explorer</title>
<style>
:root{color-scheme:light;--bg:#fcfcfb;--panel:#f5f4f0;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;
--grid:#e1e0d9;--axis:#c3c2b7;--band:#f0efec;--ztf:#2a78d6;--live:#eb6834;--pan:#a8a69f;--pec:#eb6834;
--g:#008300;--r:#e34948;--i:#eda100;--model:#0b0b0b}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;--bg:#1a1a19;--panel:#222220;
--ink:#ffffff;--ink2:#c3c2b7;--muted:#8f8d86;--grid:#2e2e2b;--axis:#4a4945;--band:#262624;--ztf:#3987e5;--live:#d95926;
--pan:#6b6a65;--pec:#d95926;--g:#2fa52f;--r:#e66767;--i:#c98500;--model:#ffffff}}
:root[data-theme="dark"]{color-scheme:dark;--bg:#1a1a19;--panel:#222220;--ink:#ffffff;--ink2:#c3c2b7;--muted:#8f8d86;
--grid:#2e2e2b;--axis:#4a4945;--band:#262624;--ztf:#3987e5;--live:#d95926;--pan:#6b6a65;--pec:#d95926;--g:#2fa52f;
--r:#e66767;--i:#c98500;--model:#ffffff}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1240px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:4px 0 2px}
.sub{color:var(--ink2);margin:0 0 12px;max-width:80ch}
.row{display:flex;gap:16px;flex-wrap:wrap}
.card{background:var(--panel);border-radius:10px;padding:12px;min-width:0}
.main{flex:1 1 640px}.side{flex:1 1 360px;max-width:460px}
.controls{display:flex;flex-wrap:wrap;gap:8px;align-items:center;margin-bottom:8px}
.seg{display:inline-flex;border:1px solid var(--axis);border-radius:8px;overflow:hidden}
.seg button{border:0;background:transparent;color:var(--ink2);padding:5px 10px;font:inherit;cursor:pointer}
.seg button.on{background:var(--ztf);color:#fff}
.chip{border:1px solid var(--axis);background:transparent;color:var(--ink);border-radius:999px;padding:3px 10px;
font:inherit;font-size:12.5px;cursor:pointer;display:inline-flex;gap:6px;align-items:center}
.chip.off{opacity:.45}.dot{width:9px;height:9px;border-radius:50%;display:inline-block}
canvas{width:100%;display:block;touch-action:none}
#hd{height:440px}#lcc{height:260px}
.stat{color:var(--ink2);font-size:13px;margin-top:6px}
.stat b{color:var(--ink);font-weight:600}
#tip{position:fixed;pointer-events:none;background:var(--bg);border:1px solid var(--axis);border-radius:6px;
padding:5px 8px;font-size:12px;display:none;z-index:5;max-width:260px}
table{border-collapse:collapse;width:100%;font-size:13px}td{padding:2px 4px;border-bottom:1px solid var(--grid)}
td:first-child{color:var(--ink2)}
.legend{display:flex;gap:12px;font-size:12px;color:var(--ink2);margin:4px 0}
#theme{margin-left:auto}
.note{font-size:12px;color:var(--muted);margin-top:10px}
</style></head><body><main>
<div class="controls"><h1>Fresh supernovae on D1's Hubble diagram</h1>
<button class="chip" id="theme" aria-label="Toggle light or dark theme">◐ theme</button></div>
<p class="sub">Each blue dot is a Type Ia supernova that ZTF found in 2018–2020 and that D1's sample never contained.
Its height is how far its brightness-based distance sits from the prediction of the universe D1 fitted.
Switch the corrections on one at a time, then click a supernova to see its light curve.</p>
<div class="row">
<section class="card main">
<div class="controls">
<span style="color:var(--ink2);font-size:13px">Correct for:</span>
<div class="seg" id="stage">
<button data-k="r0">nothing</button><button data-k="r1">stretch</button><button data-k="r2">stretch + colour</button><button data-k="rd" class="on">D1's full recipe</button>
</div>
</div>
<div class="controls" id="chips">
<button class="chip" data-s="cos"><span class="dot" style="background:var(--ztf)"></span>ZTF DR2, normal</button>
<button class="chip off" data-s="pec"><span class="dot" style="background:var(--pec)"></span>peculiar sub-types</button>
<button class="chip" data-s="pan"><span class="dot" style="background:var(--pan)"></span>D1's supernovae</button>
<button class="chip" data-s="live"><span class="dot" style="background:var(--live)"></span>last few weeks</button>
</div>
<canvas id="hd" aria-label="Hubble residuals against redshift"></canvas>
<div class="stat" id="stat"></div>
</section>
<section class="card side">
<div id="lchead" style="font-weight:600">Click a supernova</div>
<div class="legend"><span><span class="dot" style="background:var(--g)"></span> g (circles)</span>
<span><span class="dot" style="background:var(--r)"></span> r (squares)</span>
<span><span class="dot" style="background:var(--i)"></span> i (diamonds)</span></div>
<canvas id="lcc" aria-label="Light curve of the selected supernova"></canvas>
<table id="info"></table>
<p class="note">Light curves: ZTF forced photometry (DR2 for 2018–2020; alert packets via ALeRCE for the newest),
divided by the fitted peak; lines are our SALT3 fit; time runs in the supernova's own frame.</p>
</section></div>
<p class="note" id="foot"></p>
</main><div id="tip"></div>
<script id="data" type="application/json">__DATA__</script>
<script>
const D=JSON.parse(document.getElementById('data').textContent);
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
let stage='rd', show={cos:true,pec:false,pan:true,live:true}, sel=null, hover=null;
const X0=Math.log10(0.008),X1=Math.log10(0.3),Y0=-1.6,Y1=1.6;
function setup(cv){const r=cv.getBoundingClientRect(),dp=devicePixelRatio||1;cv.width=r.width*dp;cv.height=r.height*dp;
 const c=cv.getContext('2d');c.setTransform(dp,0,0,dp,0,0);return [c,r.width,r.height]}
const PAD={l:52,r:12,t:10,b:38};
function geom(w,h){return {x:z=>PAD.l+(Math.log10(z)-X0)/(X1-X0)*(w-PAD.l-PAD.r),
 y:v=>PAD.t+(v-Y0)/(Y1-Y0)*(h-PAD.t-PAD.b)}}
function nmad(a){if(!a.length)return NaN;const s=[...a].sort((p,q)=>p-q),m=s[s.length>>1];
 const d=s.map(v=>Math.abs(v-m)).sort((p,q)=>p-q);return 1.4826*d[d.length>>1]}
function val(s){if(stage!=='rd')return s[stage];return s.pec?s.rz:(s.rd??null)}
function pts(){const out=[];
 if(show.pan&&stage==='rd')for(const [z,r] of D.pan)out.push({k:'pan',z,v:r});
 for(const s of D.sn){if(s.pec?!show.pec:!(show.cos&&s.cos))continue;
  const v=val(s);if(v==null)continue;out.push({k:s.pec?'pec':'cos',z:s.z,v,s})}
 if(show.live&&stage==='rd')for(const s of D.live)out.push({k:'live',z:s.z,v:s.rd,s});return out}
let P=[];
function drawHD(){const cv=document.getElementById('hd');const [c,w,h]=setup(cv);const g=geom(w,h);
 c.clearRect(0,0,w,h);c.fillStyle=css('--band');c.fillRect(g.x(0.008),PAD.t,g.x(0.06)-g.x(0.008),h-PAD.t-PAD.b);
 c.font='11px system-ui';c.fillStyle=css('--muted');c.textAlign='right';c.fillText('volume-limited: no selection bias',g.x(0.06)-4,h-PAD.b-6);
 c.strokeStyle=css('--grid');c.lineWidth=1;c.fillStyle=css('--muted');c.textAlign='right';
 for(let v=-1.5;v<=1.51;v+=0.5){const y=g.y(-v);c.beginPath();c.moveTo(PAD.l,y);c.lineTo(w-PAD.r,y);c.stroke();c.fillText((v>0?'+':'')+v.toFixed(1),PAD.l-6,y+4)}
 c.textAlign='center';for(const z of [0.01,0.02,0.03,0.05,0.1,0.2,0.3]){const x=g.x(z);c.beginPath();c.moveTo(x,PAD.t);c.lineTo(x,h-PAD.b);c.stroke();c.fillText(z,x,h-PAD.b+15)}
 c.fillStyle=css('--ink2');c.fillText('redshift z',(PAD.l+w-PAD.r)/2,h-6);
 c.save();c.translate(13,(PAD.t+h-PAD.b)/2);c.rotate(-Math.PI/2);c.fillText('residual from D1 [mag]  (up = fainter, farther)',0,0);c.restore();
 c.strokeStyle=css('--model');c.lineWidth=1.4;c.beginPath();c.moveTo(PAD.l,g.y(0));c.lineTo(w-PAD.r,g.y(0));c.stroke();
 P=pts();const col={pan:css('--pan'),cos:css('--ztf'),pec:css('--pec'),live:css('--live')};
 for(const k of ['pan','cos','pec']){c.fillStyle=col[k];c.globalAlpha=k==='cos'?0.45:0.7;
  for(const p of P){if(p.k!==k)continue;p.px=g.x(p.z);p.py=g.y(-p.v);c.beginPath();c.arc(p.px,Math.max(PAD.t,Math.min(h-PAD.b,p.py)),k==='cos'?2:2.4,0,7);c.fill()}}
 c.globalAlpha=1;for(const p of P){if(p.k!=='live')continue;p.px=g.x(p.z);p.py=g.y(-p.v);star(c,p.px,p.py,7,col.live)}
 if(sel){const p=P.find(q=>q.s&&q.s.id===sel);if(p){c.strokeStyle=css('--ink');c.lineWidth=2;c.beginPath();c.arc(p.px,p.py,6,0,7);c.stroke()}}
 const vs=P.filter(p=>(p.k==='cos')&&p.z>=0.03&&p.z<=0.1).map(p=>p.v);
 const names={r0:'no correction',r1:'stretch corrected',r2:'stretch and colour corrected',rd:"D1's full recipe (D1's model, α, β, M)"};
 document.getElementById('stat').innerHTML=`Scatter of normal ZTF supernovae at 0.03 &lt; z &lt; 0.1, <b>${names[stage]}</b>: nMAD <b>${nmad(vs).toFixed(3)} mag</b> `+
  `(${vs.length} supernovae; 0.15 mag = 7% in distance).`+(stage==='rd'?` ZTF sits <b>${D.meta.offset>=0?'+':''}${D.meta.offset.toFixed(3)} mag</b> from D1's curve on average: DR2's flux scale is blinded, so this offset is not cosmology.`:' Each stage is centred on its median.')}
function star(c,x,y,r,f){c.beginPath();for(let i=0;i<10;i++){const a=Math.PI/2+i*Math.PI/5,rr=i%2?r*0.45:r;c.lineTo(x+rr*Math.cos(a),y-rr*Math.sin(a))}
 c.closePath();c.fillStyle=f;c.fill();c.strokeStyle=css('--bg');c.lineWidth=1;c.stroke()}
function nearest(mx,my){let b=null,bd=100;for(const p of P){if(!p.s)continue;const d=(p.px-mx)**2+(p.py-my)**2;if(d<bd){bd=d;b=p}}return b}
const hd=document.getElementById('hd'),tip=document.getElementById('tip');
hd.addEventListener('pointermove',e=>{const r=hd.getBoundingClientRect();const p=nearest(e.clientX-r.left,e.clientY-r.top);
 if(!p){tip.style.display='none';hd.style.cursor='default';return}hd.style.cursor='pointer';const s=p.s;
 tip.innerHTML=`<b>${s.n?'SN '+s.n:s.id}</b>${p.k==='live'?' (new)':''}<br>z = ${s.z}, x1 = ${s.x1}, c = ${s.c}<br>residual ${p.v>=0?'+':''}${p.v.toFixed(3)} mag`+(s.st?`<br>sub-type: ${s.st}`:'');
 tip.style.display='block';tip.style.left=Math.min(e.clientX+14,innerWidth-270)+'px';tip.style.top=(e.clientY+12)+'px'});
hd.addEventListener('pointerleave',()=>tip.style.display='none');
hd.addEventListener('click',e=>{const r=hd.getBoundingClientRect();const p=nearest(e.clientX-r.left,e.clientY-r.top);if(p){sel=p.s.id;drawHD();drawLC(p)}});
function drawLC(p){const s=p.s,L=D.lc[s.id]||{};const cv=document.getElementById('lcc');const [c,w,h]=setup(cv);c.clearRect(0,0,w,h);
 const pl=44,pr=10,pt=8,pb=30,x=v=>pl+(v+20)/80*(w-pl-pr),y=v=>pt+(1.4-v)/1.6*(h-pt-pb);
 c.strokeStyle=css('--grid');c.fillStyle=css('--muted');c.font='11px system-ui';c.textAlign='center';
 for(const v of [-20,0,20,40,60]){c.beginPath();c.moveTo(x(v),pt);c.lineTo(x(v),h-pb);c.stroke();c.fillText(v,x(v),h-pb+14)}
 c.textAlign='right';for(const v of [0,0.5,1]){c.beginPath();c.moveTo(pl,y(v));c.lineTo(w-pr,y(v));c.stroke();c.fillText(v,pl-5,y(v)+4)}
 c.fillStyle=css('--ink2');c.textAlign='center';c.fillText('days from peak (rest frame)',(pl+w-pr)/2,h-3);
 const bc=[css('--g'),css('--r'),css('--i')];
 if(L.m){for(const b in L.m){c.strokeStyle=bc[b];c.lineWidth=1.6;c.beginPath();L.m[b].forEach((v,i)=>{const X=x(D.meta.phases[i]),Y=y(v/1000);i?c.lineTo(X,Y):c.moveTo(X,Y)});c.stroke()}}
 const q=L.p||[];for(let i=0;i<q.length;i+=4){const ph=q[i]/10,b=q[i+1],f=q[i+2]/1000,e=q[i+3]/1000;if(ph<-20||ph>60)continue;
  const X=x(ph),Y=y(Math.max(-0.2,Math.min(1.4,f)));c.strokeStyle=bc[b];c.globalAlpha=0.5;c.lineWidth=1;c.beginPath();c.moveTo(X,y(f-e));c.lineTo(X,y(f+e));c.stroke();
  c.globalAlpha=0.9;c.fillStyle=bc[b];c.beginPath();if(b===0)c.arc(X,Y,2.6,0,7);else if(b===1)c.rect(X-2.3,Y-2.3,4.6,4.6);else{c.moveTo(X,Y-3);c.lineTo(X+3,Y);c.lineTo(X,Y+3);c.lineTo(X-3,Y)}c.fill()}
 c.globalAlpha=1;
 document.getElementById('lchead').textContent=(s.n?'SN '+s.n+' · ':'')+s.id+(p.k==='live'?' (classified in the last weeks)':'');
 const row=(k,v)=>v==null||v===''?'':`<tr><td>${k}</td><td>${v}</td></tr>`;
 document.getElementById('info').innerHTML=row('redshift',s.z+(s.zs?' ('+s.zs+')':''))+row('stretch x1',s.x1)+row('colour c',s.c)+
  row('peak mB (blinded scale)',s.mB)+row('sub-type',s.st)+row('host stellar mass',s.hm!=null?'10^'+s.hm+' M☉':null)+
  row('detected phases',s.np)+row('fit χ²/dof',s.chi)+row('discovered',s.disc)+row('days past peak at last point',s.age)+
  row("residual, D1's recipe",s.rd!=null?(s.rd>=0?'+':'')+s.rd.toFixed(3)+' mag':'outside D1\'s range/cuts')}
document.querySelectorAll('#stage button').forEach(b=>b.onclick=()=>{stage=b.dataset.k;document.querySelectorAll('#stage button').forEach(x=>x.classList.toggle('on',x===b));drawHD()});
document.querySelectorAll('#chips .chip').forEach(b=>b.onclick=()=>{show[b.dataset.s]=!show[b.dataset.s];b.classList.toggle('off',!show[b.dataset.s]);drawHD()});
document.getElementById('theme').onclick=()=>{const dark=matchMedia('(prefers-color-scheme: dark)').matches;const cur=document.documentElement.dataset.theme||(dark?'dark':'light');
 document.documentElement.dataset.theme=cur==='dark'?'light':'dark';drawAll()};
function drawAll(){drawHD();if(sel){const p=P.find(q=>q.s&&q.s.id===sel);if(p)drawLC(p)}}
matchMedia('(prefers-color-scheme: dark)').addEventListener('change',drawAll);addEventListener('resize',drawAll);
document.getElementById('foot').textContent=`D1's universe: flat ΛCDM, Ωm = ${D.meta.d1.Om.toFixed(3)}, H0 = ${D.meta.d1.H0.toFixed(2)} km/s/Mpc, M = ${D.meta.d1.M.toFixed(3)}. `+
 `Stretch/colour stages use α = ${D.meta.alpha.toFixed(3)}, β = ${D.meta.beta.toFixed(2)} fitted on ZTF; D1's recipe uses α = ${D.meta.d1_alpha.toFixed(3)}, β = ${D.meta.d1_beta.toFixed(2)}. `+
 `Data: ZTF SN Ia DR2 (Rigault et al. 2025), Pantheon+ (Scolnic et al. 2022), TNS + ALeRCE${D.meta.snapshot?' (snapshot '+D.meta.snapshot+')':''}.`;
drawHD();const first=P.find(p=>p.k==='live')||P.find(p=>p.k==='cos');if(first){sel=first.s.id;drawHD();drawLC(first)}
</script></body></html>
"""


def build(R: dict) -> str:
    payload = build_payload(R)
    text = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    OUT.write_text(HTML.replace("__DATA__", text))
    return str(OUT)
