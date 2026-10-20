"""figures/element_explorer.html: click an element, see where it came from.

A self-contained page (no network, no libraries): a periodic table coloured
by the measured prompt share; for the chosen element, the density of GALAH
disc giants in [X/Fe] against [Fe/H], the Gaia-Sausage-Enceladus giants on
top, and the numbers behind its colour.
"""

from __future__ import annotations

import json

import numpy as np

from config import ELEMENTS, FIGURES, REFERENCE

XR, YR = (-2.0, 0.6), (-1.0, 1.4)
NX, NY = 78, 72

#: One-line summaries of where each element is thought to come from (the
#: textbook view, for context next to the measurement).
ORIGIN_NOTES = {
    "C": "Massive-star winds and AGB stars both contribute.",
    "N": "Hydrogen burning (CNO cycle) in intermediate-mass stars; in giants it is also changed by their own mixing.",
    "O": "Core-collapse supernovae almost entirely.",
    "Na": "Core-collapse supernovae, with yields that grow with metallicity, plus AGB stars.",
    "Mg": "Core-collapse supernovae: the reference 'prompt' element here (100% by definition).",
    "Al": "Core-collapse supernovae (metallicity-dependent yields).",
    "Si": "Core-collapse supernovae, with a share from Type Ia.",
    "K": "Core-collapse supernovae; models under-produce it.",
    "Ca": "Core-collapse supernovae, with a Type Ia share.",
    "Sc": "Core-collapse supernovae; models under-produce it.",
    "Ti": "Core-collapse supernovae; models under-produce it.",
    "V": "Core-collapse and Type Ia supernovae.",
    "Cr": "Type Ia and core-collapse supernovae.",
    "Mn": "Mostly Type Ia supernovae (especially near-Chandrasekhar-mass ones).",
    "Fe": "About half from Type Ia, half from core-collapse supernovae.",
    "Co": "Core-collapse and Type Ia supernovae.",
    "Ni": "Type Ia and core-collapse supernovae.",
    "Cu": "Weak s-process in massive stars, plus supernovae.",
    "Zn": "Core-collapse supernovae (energetic 'hypernovae' help).",
    "Rb": "s-process in AGB stars plus r-process.",
    "Sr": "s-process in AGB stars, weak s-process in massive stars, light r-process.",
    "Y": "Mostly the s-process in AGB stars.",
    "Zr": "Mostly the s-process in AGB stars.",
    "Mo": "Mixed s- and r-process (some isotopes from the p-process).",
    "Ru": "Mostly the r-process.",
    "Ba": "Mostly the s-process in AGB stars (~85%).",
    "La": "Mostly the s-process in AGB stars.",
    "Ce": "Mostly the s-process in AGB stars.",
    "Nd": "Roughly half s-process, half r-process.",
    "Sm": "Mostly the r-process.",
    "Eu": "~94% r-process; the site (neutron-star mergers? rare supernovae?) is debated.",
}


def _element_payload(res, ctx):
    import ncapture as nc
    import twoprocess  # noqa: F401  (documentation link)

    d, eu = ctx.d, ctx.eu
    g = ctx.giants
    disc = ctx.disc & g
    gse = ctx.gse["GALAH orbits (Feuillet+21)"] & g
    pr = res["process"]
    s_frac = {}
    try:
        import budget

        s_frac = budget.s_fractions().to_dict()
    except Exception:
        pass
    out = {}
    for el in ELEMENTS:
        if el == "Fe":
            y = d["mg_fe"] * 0.0      # [Fe/Fe] is trivially zero: show Mg instead? keep empty
            continue
        y = eu["eu_best"] if el == "Eu" else d[f"{el.lower()}_fe"]
        s = disc & np.isfinite(y)
        H, _, _ = np.histogram2d(d["fe_h"][s], y[s], bins=(NX, NY), range=[XR, YR])
        Hq = np.where(H > 0, np.round(255 * np.log10(H + 1) / max(np.log10(H.max() + 1), 1e-9)), 0).astype(int)
        gs = gse & np.isfinite(y)
        rng = np.random.default_rng(0)
        idx = np.flatnonzero(gs)
        if idx.size > 900:
            idx = rng.choice(idx, 900, replace=False)
        med = nc.binned(d["fe_h"][s], y[s], np.round(np.arange(-1.6, 0.51, 0.1), 2))
        key = "eu" if el == "Eu" else el.lower()
        fcc = {k: (pr[k].get(key) or {}).get("f_cc") if el != "Mg" else 1.0
               for k in ("giants_ours", "dwarfs_griffith_recipe", "dwarfs_ours")}
        fp = res["ncapture"]["fingerprint"].get(el, {})
        out[el] = dict(
            hist=Hq.T.flatten().tolist(),
            gse=[[round(float(a), 2), round(float(b), 2)] for a, b in zip(d["fe_h"][idx], y[idx])],
            med=[[round(r["x"], 3), round(r["median"], 3)] for r in med],
            n_disc=int(s.sum()), n_gse=int(gs.sum()), fcc=fcc,
            griffith=REFERENCE["griffith2022_fcc"].get(el),
            s_frac=s_frac.get(el), gse_offset=fp.get("offset"), gse_offset_err=fp.get("err"),
            note=ORIGIN_NOTES.get(el, ""))
    # Iron: no [Fe/Fe] plot; give its numbers only.
    out["Fe"] = dict(hist=[], gse=[], med=[], n_disc=int(disc.sum()), n_gse=int(gse.sum()),
                     fcc={k: pr[k]["fe"]["f_cc"] for k in ("giants_ours", "dwarfs_griffith_recipe", "dwarfs_ours")},
                     griffith=REFERENCE["griffith2022_fcc"]["Fe"], s_frac=None, gse_offset=None,
                     gse_offset_err=None, note=ORIGIN_NOTES["Fe"])
    return out


PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Element Origins Explorer</title>
<style>
:root{--bg:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--grey:#e1e0d9;
--prompt:#eb6834;--delayed:#2a78d6;--bb:#4a3aa7;--gse:#e87ba4;--line:#0b0b0b;--panel:#f3f2ee}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#151514;--ink:#f0efec;--ink2:#c3c2b7;
--muted:#8f8d86;--grid:#33322f;--grey:#2c2b29;--line:#f0efec;--panel:#1e1d1b}}
:root[data-theme="dark"]{--bg:#151514;--ink:#f0efec;--ink2:#c3c2b7;--muted:#8f8d86;--grid:#33322f;--grey:#2c2b29;
--line:#f0efec;--panel:#1e1d1b}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1180px;margin:0 auto;padding:16px}
h1{font-size:20px;margin:4px 0 2px}p.sub{color:var(--ink2);margin:0 0 12px;font-size:13px}
.wrap{display:grid;grid-template-columns:minmax(0,1.15fr) minmax(0,1fr);gap:18px;align-items:start}
@media (max-width:860px){.wrap{grid-template-columns:1fr}}
#table{display:grid;grid-template-columns:repeat(18,minmax(0,1fr));gap:2px}
.cell{aspect-ratio:1/1.05;border-radius:3px;position:relative;font-size:clamp(8px,1.25vw,13px);font-weight:700;
display:flex;align-items:center;justify-content:center;color:var(--muted);background:var(--grey);overflow:hidden;
border:0;padding:0;font-family:inherit}
.cell.m{color:#fff;cursor:pointer}.cell.m:hover,.cell.sel{outline:2px solid var(--line);outline-offset:1px}
.cell .z{position:absolute;top:1px;left:2px;font-size:7px;font-weight:400;opacity:.8}
.cell span.s{position:relative;z-index:1}
.gap{height:8px;grid-column:1/-1}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:var(--ink2);margin:10px 0}
.legend i{display:inline-block;width:12px;height:12px;border-radius:2px;margin-right:5px;vertical-align:-2px}
.panel{background:var(--panel);border-radius:8px;padding:12px 14px}
.panel h2{margin:0 0 4px;font-size:18px}.note{color:var(--ink2);font-size:13px;margin:0 0 8px}
canvas{width:100%;height:auto;display:block}
table.kv{width:100%;border-collapse:collapse;font-size:13px;margin-top:8px}
table.kv td{padding:3px 4px;border-top:1px solid var(--grid)}table.kv td:last-child{text-align:right;font-variant-numeric:tabular-nums}
.bar{height:14px;border-radius:3px;overflow:hidden;display:flex;margin:6px 0 2px}
.bar div{height:100%}
footer{color:var(--muted);font-size:12px;margin-top:16px}
</style></head><body><main>
<h1>Where your atoms were made</h1>
<p class="sub">Click an element. Colours: the share of each element that arrived <b style="color:var(--prompt)">promptly</b>
(massive stars exploding within ~40 Myr) versus <b style="color:var(--delayed)">late</b> (Type Ia supernovae, AGB stars),
measured in GALAH DR4 disc giants with the two-process model.</p>
<div class="wrap"><div>
<div id="table"></div>
<div class="legend"><span><i style="background:var(--prompt)"></i>prompt</span><span><i style="background:var(--delayed)"></i>delayed</span>
<span><i style="background:var(--bb)"></i>Big Bang</span><span><i style="background:var(--grey)"></i>not measured here</span></div>
</div>
<div class="panel" id="panel"><h2 id="title"></h2><p class="note" id="note"></p>
<div class="bar" id="bar"></div><div id="barlab" style="font-size:12px;color:var(--ink2)"></div>
<canvas id="plot" width="560" height="430" aria-label="abundance trend"></canvas>
<table class="kv" id="kv"></table></div></div>
<footer>Data: GALAH DR4 (Buder et al. 2025), europium from Kane et al. (2026); comparison values: Griffith et al. (2022),
Bisterzo et al. (2014). Grey density: disc giants; pink dots: Gaia-Sausage-Enceladus giants; line: disc median.
Made in project L4 of the Astrophysics Project Atlas (AI-generated).</footer>
</main>
<script>
const DATA = __DATA__;
const SYM = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U".split(" ");
function pos(z){if(z==1)return[0,0];if(z==2)return[0,17];if(z<=10)return[1,z<=4?z-3:z-5+12];if(z<=18)return[2,z<=12?z-11:z-13+12];
if(z<=36)return[3,z-19];if(z<=54)return[4,z-37];if(z>=57&&z<=71)return[8,z-57+2];if(z>=89&&z<=103)return[9,z-89+2];
if(z<=56)return[5,z-55];if(z<=86)return[5,z-86+17];return[6,z-87];}
const css=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const T=document.getElementById('table');const cells={};
for(let r=0;r<10;r++){for(let c=0;c<18;c++){const d=document.createElement('div');d.style.gridRow=(r+1);d.style.gridColumn=(c+1);d.dataset.rc=r+','+c;T.appendChild(d);}}
T.querySelectorAll('div').forEach(d=>{if(d.dataset.rc.startsWith('7,'))d.className='gap';});
SYM.forEach((s,i)=>{const z=i+1;const [r,c]=pos(z);const slot=T.querySelector(`[data-rc="${r},${c}"]`);
 const b=document.createElement('button');b.className='cell';b.innerHTML=`<span class="z">${z}</span><span class="s">${s}</span>`;
 b.style.gridRow=slot.style.gridRow;b.style.gridColumn=slot.style.gridColumn;slot.replaceWith(b);cells[s]=b;
 const e=DATA.el[s];
 if(s==='H'||s==='He'){const f=s==='H'?1:0.91;b.style.background=`linear-gradient(to top,var(--bb) ${100*f}%,var(--grey) ${100*f}%)`;b.style.color='#fff';}
 else if(e){const f=Math.max(0,Math.min(1,e.fcc.giants_ours??NaN));if(isFinite(f)){b.classList.add('m');
  b.style.background=`linear-gradient(to top,var(--prompt) ${100*f}%,var(--delayed) ${100*f}%)`;b.onclick=()=>show(s);
  b.setAttribute('aria-label',`${s}: ${Math.round(100*f)}% prompt`);}}
});
T.querySelectorAll('div[data-rc]').forEach(d=>{if(d.className!=='gap')d.style.visibility='hidden';});
const fmt=(v,d=2)=>v==null||!isFinite(v)?'—':(+v).toFixed(d);
const pct=v=>v==null||!isFinite(v)?'—':Math.round(100*v)+'%';
function show(s){Object.values(cells).forEach(c=>c.classList.remove('sel'));cells[s].classList.add('sel');
 const e=DATA.el[s];document.getElementById('title').textContent=s+' — '+(DATA.names[s]||'');
 document.getElementById('note').textContent=e.note;
 const f=Math.max(0,Math.min(1,e.fcc.giants_ours));
 document.getElementById('bar').innerHTML=`<div style="width:${100*f}%;background:var(--prompt)"></div><div style="width:${100*(1-f)}%;background:var(--delayed)"></div>`;
 document.getElementById('barlab').textContent=`prompt ${pct(e.fcc.giants_ours)} · delayed ${pct(1-e.fcc.giants_ours)} (disc giants)`+(e.fcc.giants_ours>1?' — above 100%: no measurable delayed share':'');
 const rows=[['Prompt share, GALAH giants (this work)',pct(e.fcc.giants_ours)],['… GALAH dwarfs, Griffith et al.’s recipe',pct(e.fcc.dwarfs_griffith_recipe)],
 ['Griffith et al. 2022 (GALAH DR3)',pct(e.griffith)],['Solar s-process share (Bisterzo et al. 2014)',pct(e.s_frac)],
 ['GSE − disc at fixed [Fe/H]',e.gse_offset==null?'—':`${e.gse_offset>0?'+':''}${fmt(e.gse_offset)} ± ${fmt(e.gse_offset_err)} dex`],
 ['Disc giants plotted',e.n_disc.toLocaleString()],['GSE giants',e.n_gse.toLocaleString()]];
 document.getElementById('kv').innerHTML=rows.map(r=>`<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('');
 draw(s,e);}
function draw(s,e){const cv=document.getElementById('plot');const ctx=cv.getContext('2d');const W=cv.width,H=cv.height;
 const L=52,R=12,Tp=12,B=40;ctx.clearRect(0,0,W,H);ctx.font='12px system-ui';
 const X=[DATA.xr[0],DATA.xr[1]],Y=[DATA.yr[0],DATA.yr[1]];const px=x=>L+(x-X[0])/(X[1]-X[0])*(W-L-R);const py=y=>H-B-(y-Y[0])/(Y[1]-Y[0])*(H-Tp-B);
 ctx.strokeStyle=css('--grid');ctx.fillStyle=css('--muted');ctx.lineWidth=1;
 for(let x=-2;x<=0.5;x+=0.5){ctx.beginPath();ctx.moveTo(px(x),Tp);ctx.lineTo(px(x),H-B);ctx.stroke();ctx.fillText(x.toFixed(1),px(x)-9,H-B+16);}
 ctx.textAlign='right';for(let y=-1;y<=1.4;y+=0.5){ctx.beginPath();ctx.moveTo(L,py(y));ctx.lineTo(W-R,py(y));ctx.stroke();ctx.fillText(y.toFixed(1),L-6,py(y)+4);}ctx.textAlign='left';
 ctx.fillStyle=css('--ink2');ctx.fillText('[Fe/H]',W/2-16,H-6);ctx.save();ctx.translate(12,H/2+20);ctx.rotate(-Math.PI/2);ctx.fillText(`[${s}/Fe]`,0,0);ctx.restore();
 if(!e.hist.length){ctx.fillStyle=css('--ink2');ctx.fillText('Iron is the reference: [Fe/Fe] = 0 by definition.',L+20,H/2);return;}
 const nx=DATA.nx,ny=DATA.ny,dx=(W-L-R)/nx,dy=(H-Tp-B)/ny;const dark=css('--bg')!=='#fcfcfb';
 for(let j=0;j<ny;j++)for(let i=0;i<nx;i++){const v=e.hist[j*nx+i];if(!v)continue;const a=0.12+0.88*v/255;
  ctx.fillStyle=dark?`rgba(200,198,190,${a})`:`rgba(60,58,54,${a})`;ctx.fillRect(L+i*dx,H-B-(j+1)*dy,dx+0.5,dy+0.5);}
 ctx.fillStyle=css('--gse');e.gse.forEach(p=>{ctx.beginPath();ctx.arc(px(p[0]),py(p[1]),2.2,0,7);ctx.fill();});
 ctx.strokeStyle=css('--line');ctx.lineWidth=2;ctx.beginPath();e.med.forEach((p,i)=>i?ctx.lineTo(px(p[0]),py(p[1])):ctx.moveTo(px(p[0]),py(p[1])));ctx.stroke();}
show('Eu');
</script></body></html>
"""

NAMES = {"C": "carbon", "N": "nitrogen", "O": "oxygen", "Na": "sodium", "Mg": "magnesium", "Al": "aluminium",
         "Si": "silicon", "K": "potassium", "Ca": "calcium", "Sc": "scandium", "Ti": "titanium",
         "V": "vanadium", "Cr": "chromium", "Mn": "manganese", "Fe": "iron", "Co": "cobalt", "Ni": "nickel",
         "Cu": "copper", "Zn": "zinc", "Rb": "rubidium", "Sr": "strontium", "Y": "yttrium",
         "Zr": "zirconium", "Mo": "molybdenum", "Ru": "ruthenium", "Ba": "barium", "La": "lanthanum",
         "Ce": "cerium", "Nd": "neodymium", "Sm": "samarium", "Eu": "europium"}


def write_explorer(res, art):
    from analysis import _clean

    els = _clean(_element_payload(res, art["ctx"]))          # NaN -> null
    payload = dict(el={k: v for k, v in els.items() if v["fcc"].get("giants_ours") is not None},
                   names=NAMES, xr=XR, yr=YR, nx=NX, ny=NY)
    path = FIGURES / "element_explorer.html"
    path.write_text(PAGE.replace("__DATA__", json.dumps(payload, separators=(",", ":"), allow_nan=False)))
    return path
