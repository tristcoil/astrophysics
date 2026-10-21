"""A self-contained, offline HTML page: the look-elsewhere playground.

    python interactive.py   ->  figures/isw_explorer.html

Pick a filter radius and how many of Granett et al.'s structures to stack;
the page shows the stacked Planck image, the filtered signal and its S/N.
It also keeps your best S/N so far and tells you how often pure-noise skies,
searched over the same choices, beat it -- the cost of looking around.

All numbers are precomputed by run.py (the per-structure filter values on the
real sky, the per-choice scatter from 1,000 simulated skies, and each noise
sky's best S/N on the grid); the page only averages them.
"""

from __future__ import annotations

import json

import numpy as np

import superstructures as SS
from config import DERIVED_DIR, FIGURE_DIR, RESULTS_PATH, SCAN_N, SCAN_RADII_DEG, STACK_IMAGE_DEG


def build() -> str:
    import simpass

    r = json.loads(RESULTS_PATH.read_text())
    layout = SS.build_layout()
    data = np.load(DERIVED_DIR / "real_bank_smica.npy")
    sims = simpass.bank_sims()
    im = np.load(DERIVED_DIR / "stack_images.npz")
    payload = {"radii": list(SCAN_RADII_DEG), "ns": list(SCAN_N), "patch_deg": STACK_IMAGE_DEG}
    for kind, key in (("voids", "granett_voids"), ("clusters", "granett_clusters")):
        vals = np.array([[data[layout.rows(key, ("R", rad))[i]] for rad in SCAN_RADII_DEG]
                         for i in range(50)])
        payload[f"{kind}_values"] = np.round(np.nan_to_num(vals, nan=-9999), 2).tolist()
        payload[f"{kind}_valid"] = np.isfinite(vals).tolist()
        p = im[f"{key}__patches"]
        payload[f"{kind}_patches"] = np.round(np.nan_to_num(p, nan=0.0)).astype(int).tolist()
    # simulated scatter of every (statistic, radius, N) cell
    sig = {}
    for stat in ("voids", "clusters", "combined"):
        g = np.zeros((len(SCAN_RADII_DEG), len(SCAN_N)))
        for i, rad in enumerate(SCAN_RADII_DEG):
            for j, n in enumerate(SCAN_N):
                if stat == "combined":
                    s = SS.granett_combined(sims, layout, rad, n)
                else:
                    s = SS.stack(sims, layout, f"granett_{stat}", ("R", rad), n)
                g[i, j] = np.std(s, ddof=1)
        sig[stat] = np.round(g, 3).tolist()
    payload["sigma"] = sig
    ms = np.load(DERIVED_DIR / "lookelsewhere_max_snr_sims.npy")
    payload["max_snr_sims"] = np.round(np.sort(ms), 3).tolist()
    payload["pred_combined"] = r["stacking"]["scans"]["granett_combined"]["pred"]
    return HTML.replace("__DATA__", json.dumps(payload, separators=(",", ":")))


HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Look-Elsewhere Playground</title>
<style>
:root{--bg:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--grid:#e1e0d9;--card:#f4f3ef;
--blue:#2a78d6;--orange:#eb6834;--violet:#4a3aa7;--accent:#2a78d6}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#171716;--ink:#f2f1ec;--ink2:#c3c2b7;
--muted:#9a998f;--grid:#33332f;--card:#212120;--blue:#5598e7;--orange:#f08a5d;--violet:#8b7fd9;--accent:#5598e7}}
:root[data-theme="dark"]{--bg:#171716;--ink:#f2f1ec;--ink2:#c3c2b7;--muted:#9a998f;--grid:#33332f;--card:#212120;
--blue:#5598e7;--orange:#f08a5d;--violet:#8b7fd9;--accent:#5598e7}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:980px;margin:0 auto;padding:24px 16px 48px}
h1{font-size:24px;margin:0 0 4px}
p.lede{color:var(--ink2);margin:0 0 20px;max-width:70ch}
.grid{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:20px}
@media (max-width:760px){.grid{grid-template-columns:1fr}}
.card{background:var(--card);border-radius:10px;padding:16px}
label{display:block;font-size:13px;color:var(--ink2);margin:10px 0 4px}
input[type=range]{width:100%}
.seg{display:flex;gap:6px;flex-wrap:wrap}
.seg button{border:1px solid var(--grid);background:var(--bg);color:var(--ink);border-radius:999px;
padding:4px 12px;font:inherit;font-size:13px;cursor:pointer}
.seg button[aria-pressed=true]{background:var(--accent);border-color:var(--accent);color:#fff}
canvas{width:100%;height:auto;border-radius:8px;display:block;image-rendering:pixelated}
.big{font-size:34px;font-weight:700;font-variant-numeric:tabular-nums}
.row{display:flex;gap:20px;flex-wrap:wrap;margin-top:8px}
.stat .k{font-size:12px;color:var(--muted)}.stat .v{font-size:18px;font-weight:600;font-variant-numeric:tabular-nums}
.note{font-size:13px;color:var(--ink2)}
.bar{height:10px;border-radius:5px;background:var(--grid);overflow:hidden;margin-top:6px}
.bar>div{height:100%;background:var(--orange)}
.legend{display:flex;justify-content:space-between;font-size:11px;color:var(--muted);margin-top:4px}
button.reset{margin-top:10px;border:1px solid var(--grid);background:var(--bg);color:var(--ink);border-radius:6px;
padding:4px 10px;font:inherit;font-size:13px;cursor:pointer}
</style>
</head>
<body>
<main>
<h1>Look-elsewhere playground</h1>
<p class="lede">In 2008 Granett, Neyrinck &amp; Szapudi stacked the microwave sky on 50 supervoids and 50 superclusters
and found a 9.6&nbsp;µK signal at 4.4σ. Here are the same structures on Planck's map. Choose the filter radius and how many
of the most significant structures to stack. Try to make the signal as strong as you can, then see what that search costs.</p>
<div class="grid">
<section class="card">
<label>What to stack</label>
<div class="seg" id="stat">
<button data-v="combined" aria-pressed="true">Combined</button>
<button data-v="voids" aria-pressed="false">Supervoids (cold)</button>
<button data-v="clusters" aria-pressed="false">Superclusters (hot)</button>
</div>
<label for="rad">Filter radius: <b id="radv"></b></label>
<input id="rad" type="range" min="0" max="18" step="1" value="6">
<label for="num">Number of structures (most significant first): <b id="numv"></b></label>
<input id="num" type="range" min="0" max="4" step="1" value="4">
<canvas id="img" width="196" height="196" aria-label="Stacked CMB image"></canvas>
<div class="legend"><span>−20 µK</span><span>0</span><span>+20 µK</span></div>
<p class="note">24° × 24° patches of Planck's SMICA map (multipoles below 10 removed). Solid circle: the filter
radius; dashed: the equal-area ring it subtracts.</p>
</section>
<section class="card">
<div class="k note">Filtered signal</div>
<div class="big" id="sig"></div>
<div class="row">
<div class="stat"><div class="k">Chance scatter (1,000 simulated skies)</div><div class="v" id="sd"></div></div>
<div class="stat"><div class="k">S/N</div><div class="v" id="snr"></div></div>
<div class="stat"><div class="k">ΛCDM predicts (BOSS)</div><div class="v" id="pred"></div></div>
</div>
<hr style="border:none;border-top:1px solid var(--grid);margin:16px 0">
<div class="note">Your best S/N so far: <b id="best"></b> (<span id="bestset"></span>)</div>
<div class="note" style="margin-top:8px">If you search all 19 radii × 5 sample sizes, how often does a sky with
<i>no signal at all</i> give a best S/N at least this big?</div>
<div class="big" id="pglob" style="font-size:26px"></div>
<div class="bar"><div id="pbar"></div></div>
<p class="note" id="explain"></p>
<button class="reset" id="reset">Reset my best</button>
</section>
</div>
</main>
<script>
const D=__DATA__;
const stops=[[13,54,107],[57,135,229],[158,197,244],[240,239,236],[245,163,138],[227,73,72],[143,32,32]];
function col(v){const t=Math.max(0,Math.min(1,(v+20)/40))*(stops.length-1);const i=Math.min(stops.length-2,Math.floor(t)),f=t-i;
return stops[i].map((a,k)=>Math.round(a+(stops[i+1][k]-a)*f));}
let stat="combined",ri=6,ni=4,best=-Infinity,bestset="";
function mean(arr){let s=0,n=0;for(const x of arr){if(x>-9000){s+=x;n++}}return n?s/n:NaN}
function value(st,r,n){const take=k=>D[k+"_values"].slice(0,n).map(row=>row[r]);
if(st==="voids")return mean(take("voids"));if(st==="clusters")return mean(take("clusters"));
return (mean(take("clusters"))-mean(take("voids")))/2}
function draw(){const n=D.ns[ni],R=D.radii[ri];const c=document.getElementById("img"),g=c.getContext("2d");
const P=49,img=g.createImageData(P,P);const sum=new Float32Array(P*P);
for(let i=0;i<n;i++){const vp=D.voids_patches[i],cp=D.clusters_patches[i];
for(let y=0;y<P;y++)for(let x=0;x<P;x++){const vv=vp[y][x],cc=cp[y][x];let s;
if(stat==="voids")s=vv;else if(stat==="clusters")s=cc;else s=(cc-vv)/2;sum[(P-1-y)*P+(P-1-x)]+=s/n;}}
// light smoothing for display
const sm=new Float32Array(P*P);for(let y=0;y<P;y++)for(let x=0;x<P;x++){let s=0,w=0;
for(let dy=-2;dy<=2;dy++)for(let dx=-2;dx<=2;dx++){const yy=y+dy,xx=x+dx;if(yy<0||xx<0||yy>=P||xx>=P)continue;
const ww=Math.exp(-(dx*dx+dy*dy)/2.5);s+=sum[yy*P+xx]*ww;w+=ww;}sm[y*P+x]=s/w;}
for(let p=0;p<P*P;p++){const [r,gg,b]=col(sm[p]);img.data.set([r,gg,b,255],4*p);}
const off=document.createElement("canvas");off.width=P;off.height=P;off.getContext("2d").putImageData(img,0,0);
g.imageSmoothingEnabled=false;g.drawImage(off,0,0,c.width,c.height);
const px=c.width/D.patch_deg,cx=c.width/2;g.strokeStyle="#0b0b0b";g.lineWidth=1.5;
g.beginPath();g.arc(cx,cx,R*px,0,2*Math.PI);g.stroke();g.setLineDash([4,3]);g.beginPath();g.arc(cx,cx,R*Math.SQRT2*px,0,2*Math.PI);g.stroke();g.setLineDash([]);}
function pglobal(s){const a=D.max_snr_sims;let k=0;for(const x of a)if(x>=s)k++;return (k+1)/(a.length+1)}
function update(){const n=D.ns[ni],R=D.radii[ri];document.getElementById("radv").textContent=R.toFixed(1)+"°";
document.getElementById("numv").textContent=n+(stat==="combined"?" voids + "+n+" clusters":"");
const v=value(stat,ri,n),sd=D.sigma[stat][ri][ni],sgn=stat==="voids"?-1:1,snr=sgn*v/sd;
document.getElementById("sig").textContent=(v>=0?"+":"")+v.toFixed(1)+" µK";
document.getElementById("sd").textContent="± "+sd.toFixed(1)+" µK";
document.getElementById("snr").textContent=snr.toFixed(2)+"σ";
document.getElementById("pred").textContent=stat==="combined"?((D.pred_combined[ri]>=0?"+":"")+D.pred_combined[ri].toFixed(2)+" µK"):"≈ 0 µK";
if(snr>best){best=snr;bestset=`${stat}, ${R.toFixed(1)}°, N=${n}`;}
document.getElementById("best").textContent=best.toFixed(2)+"σ";document.getElementById("bestset").textContent=bestset;
const p=pglobal(best);document.getElementById("pglob").textContent=(100*p).toFixed(1)+"% of noise skies";
document.getElementById("pbar").style.width=Math.min(100,100*p*5)+"%";
document.getElementById("explain").textContent=`A single pre-chosen test reaching ${best.toFixed(2)}σ would be a ${(100*0.5*erfc(best/Math.SQRT2)).toFixed(2)}% fluke. `+
`Searched over radius and N, pure noise beats it ${(100*p).toFixed(1)}% of the time. (The noise skies were searched with the combined statistic; flipping to voids-only or clusters-only is a further choice, so the true cost of your search is higher still.)`;
draw();}
function erfc(x){const t=1/(1+0.5*Math.abs(x));const y=t*Math.exp(-x*x-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))));return x>=0?y:2-y}
document.getElementById("rad").addEventListener("input",e=>{ri=+e.target.value;update()});
document.getElementById("num").addEventListener("input",e=>{ni=+e.target.value;update()});
document.querySelectorAll("#stat button").forEach(b=>b.addEventListener("click",()=>{stat=b.dataset.v;
document.querySelectorAll("#stat button").forEach(x=>x.setAttribute("aria-pressed",x===b));update()}));
document.getElementById("reset").addEventListener("click",()=>{best=-Infinity;update()});
update();
</script>
</body>
</html>
"""


def main():
    out = FIGURE_DIR / "isw_explorer.html"
    out.write_text(build())
    print(f"  wrote {out.name} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
