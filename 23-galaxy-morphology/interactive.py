"""An offline page: classify galaxies yourself, then see what everyone else said.

``figures/galaxy_zoo_explorer.html`` is self-contained (images embedded, no
network needed) and follows the viewer's light or dark setting.

* **Classify**: a galaxy appears exactly as the volunteers saw it (the
  central part of its GZ2 image). Answer the first Galaxy Zoo question, then
  see how the ~40 volunteers split, the debiased fraction, our neural
  network's prediction and Zoobot's.
* **Browse**: every galaxy on the page, filterable by class and sortable.

The galaxies are a fixed random draw from the neural network's test set --
galaxies the network never saw while learning -- so its predictions shown
here are honest out-of-sample guesses.
"""

from __future__ import annotations

import base64
import io
import json

import numpy as np

import fetch_data as fd
import morphology as mo
from config import FIGURE_DIR, SEED

N_GALAXIES = 360
THUMB = 160
OUT_PATH = FIGURE_DIR / "galaxy_zoo_explorer.html"


def _thumb(path) -> str | None:
    from PIL import Image

    if not path.exists():
        return None
    img = Image.open(path).convert("RGB")
    n = img.size[0]
    k = int(round(n * 0.2))
    img = img.crop((k, k, n - k, n - k)).resize((THUMB, THUMB), Image.Resampling.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format="JPEG", quality=82, optimize=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def _choose(state: dict) -> np.ndarray:
    """Indices into the catalogue: a seeded draw from the network's test galaxies."""
    import cnn_sample

    c = state["c"]
    pred = None
    try:
        import cnn

        if cnn.PRED_PATH.exists():
            with np.load(cnn.PRED_PATH) as f:
                pred = set(int(x) for x in f["dr7objid"])
    except ImportError:
        pass
    has = np.array([fd.gz2_image_path(o).exists() for o in c["dr7objid"]], bool)
    import cnn

    s = cnn_sample.choose(c, cnn.has_archived_image(c))
    test = s["index"][s["split"] == "test"]
    test = test[has[test] & state["d"]["good"][test]]
    if pred:
        test = test[np.isin(c["dr7objid"][test], list(pred))]
    rng = np.random.default_rng(SEED)
    return np.sort(rng.choice(test, size=min(N_GALAXIES, len(test)), replace=False))


def _records(state: dict, idx: np.ndarray) -> list[dict]:
    c, d, cls = state["c"], state["d"], state["cls"]
    arm, _ = mo.arm_number(c)
    preds, zoo = {}, {}
    try:
        import cnn

        if cnn.PRED_PATH.exists():
            with np.load(cnn.PRED_PATH) as f:
                ids = f["dr7objid"]
                for k, o in enumerate(ids):
                    preds[int(o)] = {a: float(f[f"pred_{a}"][k]) for a in
                                     ("smooth", "features", "artifact", "edgeon", "bar", "spiral")}
    except ImportError:
        pass
    if fd.GZDESI_SUBSET_PATH.exists():
        import pandas as pd

        z = pd.read_parquet(fd.GZDESI_SUBSET_PATH).drop_duplicates("dr7objid").set_index("dr7objid")
        for o in c["dr7objid"][idx]:
            if int(o) in z.index:
                row = z.loc[int(o)]
                zoo[int(o)] = {
                    "smooth": float(row["smooth-or-featured_smooth_fraction"]),
                    "features": float(row["smooth-or-featured_featured-or-disk_fraction"]),
                    "artifact": float(row["smooth-or-featured_artifact_fraction"])}

    def f(x):
        return None if not np.isfinite(x) else round(float(x), 3)

    red = (d["gr0"] > state["divider"]["slope"] * d["M"] + state["divider"]["intercept"])
    from physics import COSMO

    lookback = COSMO.lookback_time(np.clip(d["z"][idx], 0, None)).to_value("Myr")
    lookback = dict(zip(idx.tolist(), lookback.tolist()))
    out = []
    for i in idx:
        o = int(c["dr7objid"][i])
        img = _thumb(fd.gz2_image_path(o))
        if img is None:
            continue
        tags = [k for k in ("smooth", "features", "spiral", "edgeon", "barred") if cls[k][i]]
        if np.nan_to_num(c["p_odd"][i]) * np.nan_to_num(c["p_merger"][i]) > 0.4:
            tags.append("merger")
        out.append({
            "id": str(o), "img": img, "z": f(d["z"][i]), "M": f(d["M"][i]),
            "lookback_myr": int(round(lookback[int(i)], -1)),
            "mass": f(d["mass"][i]), "gr": f(d["gr0"][i]), "red": bool(red[i]),
            "n": int(c["n_class"][i]), "tags": tags,
            "raw": {a: f(c["wf_" + a][i]) for a in ("smooth", "features", "artifact", "edgeon",
                                                     "notedgeon", "bar", "nobar", "spiral",
                                                     "nospiral")},
            "deb": {a: f(c["p_" + a][i]) for a in ("smooth", "features", "artifact")},
            "nq": {"q1": int(c["n_smooth"][i] + c["n_features"][i] + c["n_artifact"][i]),
                   "edge": int(c["n_edgeon"][i] + c["n_notedgeon"][i]),
                   "bar": int(c["n_bar"][i] + c["n_nobar"][i]),
                   "spiral": int(c["n_spiral"][i] + c["n_nospiral"][i])},
            "arms": mo.ARM_LABELS[arm[i]] if arm[i] >= 0 and cls["spiral"][i] else None,
            "cnn": preds.get(o), "zoo": zoo.get(o),
        })
    return out


def build(state: dict, *, offline: bool = False, quiet: bool = False) -> str:
    idx = _choose(state)
    recs = _records(state, idx)
    html = TEMPLATE.replace("__DATA__", json.dumps(recs, separators=(",", ":")))
    OUT_PATH.write_text(html, encoding="utf-8")
    if not quiet:
        print(f"  explorer: {len(recs)} galaxies, {OUT_PATH.stat().st_size / 1e6:.1f} MB -> {OUT_PATH.name}")
    return str(OUT_PATH)


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Galaxy Zoo Explorer</title>
<style>
:root {
  color-scheme: light;
  --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --axis: #c3c2b7; --ring: rgba(11,11,11,0.10);
  --smooth: #eb6834; --features: #2a78d6; --artifact: #898781; --accent: #2a78d6;
  --good: #006300; --chip: #f0efec;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --axis: #383835; --ring: rgba(255,255,255,0.10);
    --smooth: #d95926; --features: #3987e5; --artifact: #898781; --accent: #3987e5;
    --good: #0ca30c; --chip: #2c2c2a;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --axis: #383835; --ring: rgba(255,255,255,0.10);
  --smooth: #d95926; --features: #3987e5; --artifact: #898781; --accent: #3987e5;
  --good: #0ca30c; --chip: #2c2c2a;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
       font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1040px; margin: 0 auto; padding: 20px 16px 48px; }
h1 { font-size: 24px; margin: 4px 0 2px; }
.sub { color: var(--ink2); margin: 0 0 16px; max-width: 760px; }
nav { display: flex; gap: 8px; margin-bottom: 16px; }
nav button, .chip { border: 1px solid var(--ring); background: var(--surface); color: var(--ink);
       border-radius: 999px; padding: 6px 14px; font: inherit; cursor: pointer; }
nav button[aria-pressed="true"], .chip[aria-pressed="true"] { background: var(--ink); color: var(--surface); }
.card { background: var(--surface); border: 1px solid var(--ring); border-radius: 12px; padding: 16px; }
.play { display: grid; grid-template-columns: minmax(0, 340px) minmax(0, 1fr); gap: 20px; align-items: start; }
@media (max-width: 720px) { .play { grid-template-columns: 1fr; } }
.galaxy { width: 100%; aspect-ratio: 1; border-radius: 10px; background: #000; display: block;
          image-rendering: auto; }
.q { font-weight: 600; margin: 4px 0 10px; }
.answers { display: flex; flex-wrap: wrap; gap: 8px; }
.answers button { flex: 1 1 140px; padding: 10px 12px; border-radius: 10px; border: 1px solid var(--ring);
                  background: var(--chip); color: var(--ink); font: inherit; cursor: pointer; }
.answers button:hover { border-color: var(--accent); }
.answers button:disabled { cursor: default; opacity: 0.6; }
.answers button.picked { outline: 2px solid var(--accent); opacity: 1; }
.reveal h3 { font-size: 13px; text-transform: uppercase; letter-spacing: .04em; color: var(--muted);
             margin: 16px 0 6px; }
.bar { display: flex; height: 22px; border-radius: 6px; overflow: hidden; gap: 2px; background: var(--surface); }
.bar span { display: block; height: 100%; }
.legend { display: flex; gap: 14px; font-size: 13px; color: var(--ink2); margin-top: 4px; flex-wrap: wrap; }
.legend i { display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 5px; }
.row { display: grid; grid-template-columns: 150px 1fr 48px; gap: 10px; align-items: center; margin: 6px 0;
       font-size: 13px; }
.row .lab { color: var(--ink2); }
.row .val { text-align: right; font-variant-numeric: tabular-nums; }
.verdict { font-size: 15px; margin: 6px 0 0; }
.meta { font-size: 13px; color: var(--ink2); margin-top: 10px; }
.score { font-size: 13px; color: var(--ink2); margin-top: 10px; }
.next { margin-top: 14px; padding: 9px 16px; border-radius: 10px; border: 0; background: var(--ink);
        color: var(--surface); font: inherit; cursor: pointer; }
.filters { display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 12px; align-items: center; }
.filters select { font: inherit; padding: 5px 8px; border-radius: 8px; border: 1px solid var(--ring);
                  background: var(--surface); color: var(--ink); }
.grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(120px, 1fr)); gap: 10px; }
.thumb { border: 0; padding: 0; background: none; cursor: pointer; text-align: left; color: var(--ink); }
.thumb img { width: 100%; aspect-ratio: 1; border-radius: 8px; display: block; background: #000; }
.thumb small { display: block; font-size: 12px; color: var(--ink2); margin-top: 3px; }
.hidden { display: none !important; }
footer { color: var(--muted); font-size: 12.5px; margin-top: 28px; }
</style>
</head>
<body>
<main>
  <h1>Classify like a volunteer</h1>
  <p class="sub">Between 2009 and 2010, 83,943 people answered a tree of questions about 304,122 galaxies
  for Galaxy Zoo 2. Here are some of those galaxies, as they saw them. Make your call, then see how the
  volunteers split, and what two neural networks would have said.</p>
  <nav>
    <button id="tab-play" aria-pressed="true">Classify</button>
    <button id="tab-browse" aria-pressed="false">Browse</button>
  </nav>

  <section id="play" class="card">
    <div class="play">
      <div><img id="g-img" class="galaxy" alt="SDSS image of a galaxy"></div>
      <div>
        <p class="q">Is the galaxy simply smooth and rounded, with no sign of a disk?</p>
        <div class="answers" id="answers">
          <button data-a="smooth">Smooth</button>
          <button data-a="features">Features or disk</button>
          <button data-a="artifact">Star or artifact</button>
        </div>
        <div id="reveal" class="reveal hidden"></div>
        <p class="score" id="score"></p>
      </div>
    </div>
  </section>

  <section id="browse" class="hidden">
    <div class="filters" id="filters"></div>
    <div class="grid" id="grid"></div>
    <div id="detail" class="card hidden" style="margin-top:14px"></div>
  </section>

  <footer>Images: SDSS gri composites exactly as shown to Galaxy Zoo 2 volunteers (Willett et al. 2013;
  Zenodo record 3565489). "Volunteers": consistency-weighted vote fractions; "debiased": corrected for
  distance by Hart et al. (2016). "Our network": a small CNN trained in this project on 9,000 other GZ2
  galaxies; these galaxies were held out of its training. "Zoobot": Walmsley et al. (2023), trained on
  deeper DESI Legacy Surveys images, which show more disks. Part of the Astrophysics Project Atlas, J2.</footer>
</main>
<script>
const DATA = __DATA__;
const $ = (id) => document.getElementById(id);
const pct = (x) => x == null ? "–" : Math.round(100 * x) + "%";
const ANS = {smooth: "Smooth", features: "Features or disk", artifact: "Star or artifact"};
let order = [], pos = 0, agree = 0, played = 0;

function store(k, v) { try { localStorage.setItem("gz2-" + k, JSON.stringify(v)); } catch (e) {} }
function load(k, d) { try { const v = localStorage.getItem("gz2-" + k); return v ? JSON.parse(v) : d; } catch (e) { return d; } }

function shuffle(n) {
  const a = [...Array(n).keys()];
  for (let i = n - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; }
  return a;
}
function bar3(o) {
  // Debiased fractions are corrected answer by answer and need not sum to one,
  // so every bar is normalised to the total of its three answers.
  let s = o.smooth || 0, f = o.features || 0, a = o.artifact == null ? Math.max(0, 1 - s - f) : o.artifact;
  const t = s + f + a || 1; s /= t; f /= t; a /= t;
  return `<div class="bar" role="img" aria-label="smooth ${pct(s)}, features ${pct(f)}">
    <span style="width:${100*s}%;background:var(--smooth)"></span>
    <span style="width:${100*f}%;background:var(--features)"></span>
    <span style="width:${100*a}%;background:var(--artifact)"></span></div>`;
}
function row(label, o) {
  if (!o) return `<div class="row"><span class="lab">${label}</span><span class="lab">not available</span><span></span></div>`;
  return `<div class="row"><span class="lab">${label}</span>${bar3(o)}<span class="val">${pct(o.features)}</span></div>`;
}
function ago(myr) { return myr >= 1000 ? (myr / 1000).toFixed(1) + "\u202fbillion years" : myr + "\u202fmillion years"; }
function majority(o) {
  const k = ["smooth", "features", "artifact"].reduce((a, b) => ((o[b] || 0) > (o[a] || 0) ? b : a));
  return k;
}
function details(g) {
  const extra = [];
  if (g.raw.edgeon != null && g.nq.edge >= 5) extra.push(`edge-on: ${pct(g.raw.edgeon)} of ${g.nq.edge} asked`);
  if (g.raw.bar != null && g.nq.bar >= 5) extra.push(`bar: ${pct(g.raw.bar)} of ${g.nq.bar}`);
  if (g.raw.spiral != null && g.nq.spiral >= 5) extra.push(`spiral arms: ${pct(g.raw.spiral)} of ${g.nq.spiral}`);
  if (g.arms) extra.push(`most-voted arm count: ${g.arms}`);
  return `<h3>How it was answered</h3>
    ${row("Volunteers", g.raw)}${row("Debiased", g.deb)}${row("Our network", g.cnn)}${row("Zoobot", g.zoo)}
    <div class="legend"><span><i style="background:var(--smooth)"></i>smooth</span>
      <span><i style="background:var(--features)"></i>features or disk</span>
      <span><i style="background:var(--artifact)"></i>star or artifact</span>
      <span>(right: % saying features)</span></div>
    <h3>Further down the tree</h3><div class="meta">${extra.join(" · ") || "–"}</div>
    <h3>The galaxy</h3><div class="meta">redshift ${g.z} · its light set out ${ago(g.lookback_myr)} ago ·
      ${g.mass ? "stellar mass 10<sup>" + g.mass.toFixed(1) + "</sup> suns · " : ""}
      ${g.red ? "red" : "blue"} (g−r ${g.gr}) · ${g.n} volunteers · SDSS DR7 ${g.id}</div>`;
}
function showPlay() {
  const g = DATA[order[pos]];
  $("g-img").src = g.img;
  $("reveal").classList.add("hidden");
  document.querySelectorAll("#answers button").forEach(b => { b.disabled = false; b.classList.remove("picked"); });
  $("score").textContent = played ? `You matched the volunteers' majority on ${agree} of ${played}.` : "";
}
function answer(a) {
  const g = DATA[order[pos]];
  document.querySelectorAll("#answers button").forEach(b => { b.disabled = true; if (b.dataset.a === a) b.classList.add("picked"); });
  const m = majority(g.raw);
  played += 1; if (m === a) agree += 1;
  store("score", {agree, played});
  const share = g.raw[a] == null ? "–" : pct(g.raw[a]);
  $("reveal").innerHTML = `<p class="verdict">${share} of the volunteers agreed with you; most said <b>${ANS[m]}</b>.</p>
    ${details(g)}<button class="next" id="next">Next galaxy</button>`;
  $("reveal").classList.remove("hidden");
  $("next").onclick = () => { pos = (pos + 1) % order.length; showPlay(); };
  $("score").textContent = `You matched the volunteers' majority on ${agree} of ${played}.`;
}

const FILTERS = [["all", "All"], ["smooth", "Smooth"], ["features", "Disks"], ["spiral", "Spirals"],
                 ["barred", "Barred"], ["edgeon", "Edge-on"], ["merger", "Mergers"]];
let filt = "all", sortKey = "features";
function drawGrid() {
  let list = DATA.filter(g => filt === "all" || g.tags.includes(filt));
  const key = {features: g => -(g.raw.features || 0), smooth: g => -(g.raw.smooth || 0), z: g => g.z,
               colour: g => -(g.gr || 0), disagree: g => -Math.abs((g.cnn ? g.cnn.features : g.raw.features) - g.raw.features)}[sortKey];
  list = list.slice().sort((a, b) => key(a) - key(b));
  $("grid").innerHTML = list.map(g => `<button class="thumb" data-id="${g.id}">
      <img src="${g.img}" alt="galaxy ${g.id}" loading="lazy">
      <small>features: ${pct(g.raw.features)} raw, ${pct(g.deb.features)} debiased · z ${g.z}</small></button>`).join("");
  document.querySelectorAll(".thumb").forEach(b => b.onclick = () => {
    const g = DATA.find(x => x.id === b.dataset.id);
    $("detail").innerHTML = `<div class="play"><div><img class="galaxy" src="${g.img}" alt="galaxy ${g.id}"></div><div>${details(g)}</div></div>`;
    $("detail").classList.remove("hidden");
    $("detail").scrollIntoView({behavior: "smooth", block: "nearest"});
  });
}
function initBrowse() {
  $("filters").innerHTML = FILTERS.map(([k, l]) =>
      `<button class="chip" data-f="${k}" aria-pressed="${k === filt}">${l} (${k === "all" ? DATA.length : DATA.filter(g => g.tags.includes(k)).length})</button>`).join("")
    + `<label style="margin-left:auto;font-size:13px;color:var(--ink2)">sort
       <select id="sort"><option value="features">most 'features'</option><option value="smooth">most 'smooth'</option>
       <option value="z">nearest first</option><option value="colour">reddest first</option>
       <option value="disagree">network disagrees most</option></select></label>`;
  document.querySelectorAll(".chip").forEach(b => b.onclick = () => {
    filt = b.dataset.f; store("filter", filt);
    document.querySelectorAll(".chip").forEach(c => c.setAttribute("aria-pressed", c.dataset.f === filt));
    drawGrid();
  });
  $("sort").value = sortKey;
  $("sort").onchange = (e) => { sortKey = e.target.value; drawGrid(); };
  drawGrid();
}
function tab(which) {
  $("play").classList.toggle("hidden", which !== "play");
  $("browse").classList.toggle("hidden", which !== "browse");
  $("tab-play").setAttribute("aria-pressed", which === "play");
  $("tab-browse").setAttribute("aria-pressed", which === "browse");
  store("tab", which);
}
document.querySelectorAll("#answers button").forEach(b => b.onclick = () => answer(b.dataset.a));
$("tab-play").onclick = () => tab("play");
$("tab-browse").onclick = () => tab("browse");
const saved = load("score", null); if (saved) { agree = saved.agree; played = saved.played; }
filt = load("filter", "all");
order = shuffle(DATA.length);
showPlay(); initBrowse(); tab(load("tab", "play"));
</script>
</body>
</html>
"""
