// Headless check of the "1919 Plate Lab" page (figures/plate_lab.html), built by explorer.py.
//
//   node tools/check_explorer.js                 # the real page
//   node tools/check_explorer.js --page <file>   # another build of it (e.g. a scratch test page)
//
// Opens the page in light, dark and phone (390 px) contexts and fails on any console or page error. Drives the
// controls and compares the numbers the page computes (data-value attributes and window.__plateLab) with the
// Python reference values embedded in the page's data (D.*.reference): Sobral 4-inch, astrographic, Principe,
// and Wallal when its data is present; and checks that the jackknife notes quote D.*.jackknife correctly. Also
// checks that nothing scrolls sideways at phone width, that chart labels do not overlap, and runs the artifact
// fragment inside a viewer-like skeleton. Screenshots go to data/scratch/explorer/. Exit code 1 on any failure.
"use strict";
const path = require("path");
const fs = require("fs");
// Playwright: PLAYWRIGHT_PATH if set, else a normal install, else the npx cache on this machine.
function loadPlaywright() {
  const tries = [process.env.PLAYWRIGHT_PATH, "playwright",
                 path.join(process.env.HOME || "", ".npm/_npx/e41f203b7505f1fb/node_modules/playwright")];
  for (const t of tries) {
    if (!t) continue;
    try { return require(t); } catch (e) { /* try the next */ }
  }
  throw new Error("Playwright not found: npm install playwright, or set PLAYWRIGHT_PATH");
}
const {chromium} = loadPlaywright();

const PROJECT = path.resolve(__dirname, "..");
const args = process.argv.slice(2);
const pageArg = args.includes("--page") ? path.resolve(args[args.indexOf("--page") + 1]) : null;
const PAGE = pageArg || path.join(PROJECT, "figures", "plate_lab.html");
const ARTIFACT = path.join(PROJECT, "data", "derived", "plate_lab_artifact.html");
const OUT = path.join(PROJECT, "data", "scratch", "explorer");
const TAG = pageArg ? path.basename(PAGE, ".html") + "_" : "";
fs.mkdirSync(OUT, {recursive: true});

const failures = [];
const lines = [];
function log(s) { lines.push(s); console.log(s); }
function fail(s) { failures.push(s); log("FAIL " + s); }
function compare(name, got, ref, tol) {
  const d = Math.abs(got - ref);
  const ok = Number.isFinite(got) && d <= tol;
  const s = `${name.padEnd(42)} page ${Number.isFinite(got) ? got.toFixed(7) : String(got)}  reference ${ref}  |diff| ${Number.isFinite(d) ? d.toExponential(2) : "n/a"}  (tol ${tol})`;
  if (ok) log("ok   " + s); else fail(s);
}

async function value(page, sel) {
  const v = await page.$eval(sel, el => el.dataset.value);
  return v === "" || v === undefined ? NaN : parseFloat(v);
}
async function click(page, sel) { await page.click(sel); await page.waitForTimeout(30); }
async function setRange(page, sel, v) {
  await page.$eval(sel, (el, val) => { el.value = val; el.dispatchEvent(new Event("input", {bubbles: true})); }, String(v));
  await page.waitForTimeout(30);
}
async function setChecked(page, sel, on) {
  const cur = await page.$eval(sel, el => el.checked);
  if (cur !== on) await click(page, sel);
}
async function tab(page, name) { await click(page, `#tab-${name}`); await page.waitForTimeout(80); }

// Overlapping text in the visible chart, and text outside the chart's box (screen coordinates).
async function labelAudit(page, where) {
  const res = await page.evaluate(() => {
    const out = [];
    document.querySelectorAll("svg.chart").forEach(svg => {
      if (!svg.getClientRects().length || svg.closest("[hidden]")) return;
      const sb = svg.getBoundingClientRect();
      if (svg.querySelectorAll("*").length < 20) out.push(`${svg.id}: almost empty (${svg.querySelectorAll("*").length} elements)`);
      for (const el of svg.querySelectorAll("*")) {
        for (const a of el.attributes) {
          if (/NaN|undefined|Infinity/.test(a.value)) { out.push(`${svg.id}: <${el.tagName} ${a.name}="${a.value.slice(0, 60)}">`); break; }
        }
      }
      const boxes = [...svg.querySelectorAll("text")].filter(t => t.textContent.trim()).map(t => {
        const b = t.getBoundingClientRect();
        return {t: t.textContent, x0: b.left + 1, y0: b.top + 1.5, x1: b.right - 1, y1: b.bottom - 1.5};
      });
      for (const b of boxes) {
        if (b.x0 < sb.left - 1 || b.x1 > sb.right + 1 || b.y0 < sb.top - 1 || b.y1 > sb.bottom + 1) out.push(`${svg.id}: "${b.t}" outside the chart`);
      }
      for (let i = 0; i < boxes.length; i++) {
        for (let j = i + 1; j < boxes.length; j++) {
          const a = boxes[i], b = boxes[j];
          const w = Math.min(a.x1, b.x1) - Math.max(a.x0, b.x0), h = Math.min(a.y1, b.y1) - Math.max(a.y0, b.y0);
          if (w > 0.5 && h > 0.5) out.push(`${svg.id}: "${a.t}" overlaps "${b.t}"`);
        }
      }
      // text over small marks (data points, markers, arrowheads)
      const marks = [...svg.querySelectorAll("circle, rect, path")].filter(el => !el.closest("defs, clipPath"))
        .map(el => ({el, b: el.getBoundingClientRect()})).filter(({b}) => b.width > 0 && b.width <= 20 && b.height <= 20);
      for (const a of boxes) {
        for (const {el, b} of marks) {
          const w = Math.min(a.x1, b.right - 1) - Math.max(a.x0, b.left + 1), h = Math.min(a.y1, b.bottom - 1) - Math.max(a.y0, b.top + 1);
          if (w > 1 && h > 1) { out.push(`${svg.id}: "${a.t}" covers a <${el.tagName} class="${el.getAttribute("class")}">`); break; }
        }
      }
    });
    return out;
  });
  for (const r of res) fail(`${where}: ${r}`);
  return res.length;
}
// The page's current Wallal solution must be a minimum of its weighted sum of squares over the measures of the
// selected pairs: no small step in beta, n or any pair's scale term alpha may lower it.
async function localMinimum(page, where) {
  const r = await page.evaluate(() => {
    const L = window.__plateLab, r = L.results.wallal;
    if (!r.ok) return {ok: false, why: r.why};
    const ms = L.wallalMeasures(r.options), su2 = L.D.wallal.sigma_unit ** 2;
    const ssr = q => ms.reduce((s, m) => s + m.w / su2 * (m.y - q.beta * Math.pow(m.x, q.n) + (r.scale ? q.alpha[m.k] * m.rDeg : 0)) ** 2, 0);
    const base = {beta: r.beta, n: r.n, alpha: r.alpha.slice()}, s0 = ssr(base);
    let worst = 0;
    const steps = [["beta", null]].concat(r.power ? [["n", null]] : [], r.scale ? r.alpha.map((_, k) => ["alpha", k]) : []);
    for (const [key, k] of steps) {
      for (const sgn of [-1, 1]) {
        const q = {beta: base.beta, n: base.n, alpha: base.alpha.slice()};
        if (k === null) q[key] += sgn * 1e-4 * Math.max(Math.abs(base[key]), 1);
        else q.alpha[k] += sgn * 1e-5;
        worst = Math.min(worst, (ssr(q) - s0) / s0);
      }
    }
    return {ok: true, worst, beta: r.beta, n: r.n, s0, chi2: r.chi2};
  });
  if (!r.ok) fail(`${where}: no fit (${r.why})`);
  else if (r.worst < -1e-12 || !Number.isFinite(r.beta)) fail(`${where}: not at a minimum (chi2 drops by a fraction ${r.worst})`);
  else if (Math.abs(r.s0 - r.chi2) > 1e-6 * r.chi2) fail(`${where}: chi2 ${r.chi2} differs from the recomputed ${r.s0}`);
  else log(`ok   ${where}: beta ${r.beta.toFixed(5)}, n ${r.n.toFixed(5)}: a minimum of chi2 (${r.s0.toFixed(3)})`);
}
// The reference power fit (scipy, default tolerances) gives beta and n only: with the scale terms profiled out
// (best alpha for each pair at that beta and n), the page's chi2 must be no larger than the reference's.
async function profileSSR(page, where, ref) {
  const res = await page.evaluate(([bRef, nRef]) => {
    const L = window.__plateLab, r = L.results.wallal, ms = L.wallalMeasures(r.options), su2 = L.D.wallal.sigma_unit ** 2;
    const prof = (beta, n) => {
      const P = r.options.pairs.length, swr = new Array(P).fill(0), swrr = new Array(P).fill(0);
      for (const m of ms) { const e = m.y - beta * Math.pow(m.x, n); swr[m.k] += m.w * m.rDeg * e; swrr[m.k] += m.w * m.rDeg * m.rDeg; }
      const alpha = swr.map((v, k) => -v / swrr[k]);
      return ms.reduce((s, m) => s + m.w / su2 * (m.y - beta * Math.pow(m.x, n) + alpha[m.k] * m.rDeg) ** 2, 0);
    };
    return {page: prof(r.beta, r.n), ref: prof(bRef, nRef)};
  }, [ref.beta, ref.n]);
  if (res.page <= res.ref * (1 + 1e-12)) log(`ok   ${where}: chi2 ${res.page.toFixed(6)} <= ${res.ref.toFixed(6)} at the reference beta and n`);
  else fail(`${where}: chi2 ${res.page} is larger than at the reference beta and n (${res.ref})`);
}
async function noSideScroll(page, where) {
  const r = await page.evaluate(() => ({sw: document.documentElement.scrollWidth, cw: document.documentElement.clientWidth,
    bw: document.body.scrollWidth}));
  if (r.sw > r.cw + 1 || r.bw > r.cw + 1) fail(`${where}: page scrolls sideways (scrollWidth ${r.sw}/${r.bw} > ${r.cw})`);
  else log(`ok   ${where}: no sideways scroll (${r.sw} <= ${r.cw})`);
}

async function numbers(page, D) {
  const I = D.inch4, R4 = I.reference, tol4 = 2e-5;
  // ---- Sobral 4-inch ----
  await tab(page, "inch4");
  const read4 = async (key, withSigma) => {
    compare(`4-inch ${key}: alpha50`, await value(page, "#i4-alpha50"), R4[key].alpha50, tol4);
    compare(`4-inch ${key}: alpha50_err`, await value(page, "#i4-alpha50-err"), R4[key].alpha50_err, tol4);
    if (withSigma && R4[key].sigma_by_group) {
      const sig = await page.evaluate(() => window.__plateLab.results.inch4.sigma);
      for (const g of Object.keys(R4[key].sigma_by_group)) compare(`4-inch ${key}: sigma ${g}`, sig[g], R4[key].sigma_by_group[g], tol4);
    }
    const limb = await value(page, "#i4-limb"), a = await value(page, "#i4-alpha50");
    const k = await page.evaluate(() => window.__plateLab.TO_LIMB);
    if (Math.abs(limb - a * D.meta.unit_arcmin / D.meta.sun_radius_true_arcmin) > 1e-9 || Math.abs(k - 50 / D.meta.sun_radius_true_arcmin) > 1e-12) {
      fail(`4-inch ${key}: limb ${limb} is not alpha50 x 50/${D.meta.sun_radius_true_arcmin}`);
    }
  };
  await read4("both", true);
  await click(page, '#i4-coords button[data-v="RA"]'); await read4("RA", true);
  await click(page, '#i4-coords button[data-v="Dec"]'); await read4("Dec", true);
  await click(page, '#i4-coords button[data-v="both"]');
  await setChecked(page, "#i4-printed", true); await read4("both_as_printed");
  await setChecked(page, "#i4-printed", false);
  for (const p of ["I", "15_1"]) {
    await setChecked(page, `#i4-plate-${p}`, false); await read4(`drop_plate_${p}`); await setChecked(page, `#i4-plate-${p}`, true);
  }
  for (const s of [11, 2]) {
    await setChecked(page, `#i4-star-${s}`, false); await read4(`drop_star_${s}`); await setChecked(page, `#i4-star-${s}`, true);
  }
  await read4("both");
  // guard rail: no comparison plates -> a message, no numbers, no error
  for (const p of I.plates.filter(p => !I.eclipse.includes(p))) await setChecked(page, `#i4-plate-${p}`, false);
  const msg = await page.$eval("#i4-msg", el => (el.hidden ? "" : el.textContent));
  if (!msg) fail("4-inch: no message when every comparison plate is left out"); else log("ok   4-inch: message shown with no comparison plates");
  await click(page, "#i4-reset");
  await read4("both");

  // ---- Sobral astrographic ----
  const A = D.astro.reference, tolA = 2e-5;
  await tab(page, "astro");
  const readA = async key => {
    compare(`astrographic ${key}: alpha50`, await value(page, "#ia-alpha50"), A[key].alpha50, tolA);
    compare(`astrographic ${key}: alpha50_err`, await value(page, "#ia-alpha50-err"), A[key].alpha50_err, tolA);
    if (A[key].e_rev !== undefined) {
      compare(`astrographic ${key}: e_rev`, await value(page, "#ia-efit"), A[key].e_rev, tolA);
      compare(`astrographic ${key}: e_err_rev`, await value(page, "#ia-efit-err"), A[key].e_err_rev, tolA);
    }
  };
  await readA("common");
  await click(page, '#ia-mode button[data-v="per_plate"]'); await readA("per_plate");
  await setRange(page, "#ia-e", "0.02");
  const mode = await page.$eval('#ia-mode button[aria-pressed="true"]', b => b.dataset.v);
  if (mode !== "fixed") fail(`astrographic: moving the slider left the mode at ${mode}`);
  await readA("fixed_0.020");
  await click(page, '#panel-astro .ticks button[data-e="0.082"]'); await readA("fixed_0.082");
  await click(page, "#ia-reset"); await readA("common");

  // ---- Principe ----
  const P = D.principe.reference, tolP = 1e-4;
  await tab(page, "principe");
  const readP = async key => {
    compare(`Principe ${key}: limb`, await value(page, "#ip-limb"), P[key].limb, tolP);
    compare(`Principe ${key}: limb_err`, await value(page, "#ip-limb-err"), P[key].limb_err, tolP);
    compare(`Principe ${key}: limb_err_measures`, await value(page, "#ip-err-meas"), P[key].limb_err_measures, tolP);
    for (const c of Object.keys(P[key].per_comparison)) {
      compare(`Principe ${key}: ${c}`, await value(page, `#ip-c-${c}`), P[key].per_comparison[c], tolP);
    }
  };
  await readP("check_scale");
  await click(page, '#ip-scale button[data-v="free"]'); await readP("free_scale");
  await click(page, "#ip-reset"); await readP("check_scale");

  // ---- Wallal ----
  await tab(page, "wallal");
  if (!D.wallal) {
    const shown = await page.$eval("#iw-empty", el => !el.hidden && el.getClientRects().length > 0);
    if (!shown) fail("Wallal: D.wallal is null but the 'coming soon' card is not shown"); else log("ok   Wallal: data not yet available, 'coming soon' card shown");
    return;
  }
  const Wr = D.wallal.reference, tolW = 1e-4;
  const readW = async (name, ref) => {
    compare(`Wallal ${name}: beta`, await value(page, "#iw-k"), ref.beta, tolW);
    if (ref.beta_err !== undefined) compare(`Wallal ${name}: beta_err`, await value(page, "#iw-k-err"), ref.beta_err, tolW);
    if (ref.n !== undefined) compare(`Wallal ${name}: n`, await value(page, "#iw-n"), ref.n, tolW);
    if (ref.n_err !== undefined) compare(`Wallal ${name}: n_err`, await value(page, "#iw-n-err"), ref.n_err, tolW);
    if (ref.chi2 !== undefined) {
      const res = await page.evaluate(() => ({chi2: window.__plateLab.results.wallal.chi2, dof: window.__plateLab.results.wallal.dof}));
      compare(`Wallal ${name}: chi2`, res.chi2, ref.chi2, 2e-3);
      compare(`Wallal ${name}: dof`, res.dof, ref.dof, 0);
    }
  };
  await readW("all six pairs", Wr.fit);
  await localMinimum(page, "Wallal fit");
  await setChecked(page, "#iw-power", true);
  await readW("power fit", Wr.power);
  await profileSSR(page, "Wallal power fit", Wr.power);
  await localMinimum(page, "Wallal power fit");
  await setChecked(page, "#iw-noscale", true);
  await localMinimum(page, "Wallal power fit without scale terms");
  await setChecked(page, "#iw-power", false);
  await localMinimum(page, "Wallal fit without scale terms");
  await click(page, "#iw-reset");
  // each plate pair alone (the 1928 Table 9 solutions redone): the main fit and the card's per-pair table
  for (let i = 0; i < D.wallal.pairs.length; i++) {
    const p = D.wallal.pairs[i];
    compare(`Wallal ${p} alone (card table)`, await value(page, `#iw-p-${i}`), Wr.per_pair_beta[p], tolW);
    for (let j = 0; j < D.wallal.pairs.length; j++) await setChecked(page, `#iw-pair-${j}`, j === i);
    compare(`Wallal ${p} alone (main fit)`, await value(page, "#iw-k"), Wr.per_pair_beta[p], tolW);
  }
  await click(page, "#iw-reset");
  await readW("after reset", Wr.fit);
}

// The jackknife notes of the 4-inch and astrographic result cards: each must quote the formal and jackknife errors of
// the full selection (D.*.reference, D.*.jackknife), and the distances in sigma from Einstein's deflection and the
// Newtonian half that the largest of them (the error the write-up adopts) gives.
async function jackknifeNotes(page, D) {
  const k = D.meta.unit_arcmin / D.meta.sun_radius_true_arcmin, gr = D.meta.gr_limb;
  const f = (v, d) => v.toFixed(d).replace("-", "−");
  for (const [sel, inst, ref, jk] of [["#i4-jk", "4-inch", D.inch4.reference.both, D.inch4.jackknife],
                                       ["#ia-jk", "astrographic", D.astro.reference.common, D.astro.jackknife]]) {
    if (!jk) { fail(`${inst}: no jackknife errors in the page's data`); continue; }
    const errs = [jk.alpha50_err_plates, jk.alpha50_err_stars].filter(v => v !== undefined).map(v => v * k);
    const limb = ref.alpha50 * k, formal = ref.alpha50_err * k, err = Math.max(formal, ...errs);
    const want = [`${f(limb, 2)}″ ± ${f(formal, 2)}″`, ...errs.map(e => `± ${f(e, 2)}″`),
                  `${f(Math.abs(limb - gr) / err, 1)}σ ${limb >= gr ? "above" : "below"} Einstein's`,
                  `${f(Math.abs(limb - gr / 2) / err, 1)}σ ${limb >= gr / 2 ? "above" : "below"} the Newtonian 0.87″`];
    const text = await page.$eval(sel, el => (el.hidden ? "" : el.textContent.replace(/\u00a0/g, " ")));
    const missing = want.filter(w => !text.includes(w));
    if (!text || missing.length) fail(`${inst}: jackknife note "${text}" lacks ${missing.map(w => `"${w}"`).join(", ")}`);
    else log(`ok   ${inst}: jackknife note quotes ${want.map(w => `"${w}"`).join(", ")}`);
  }
}

// Selections too small to fit must give a message and no numbers (and no error); other states get a label audit.
async function edgeCases(page, D) {
  const expectMsg = async (panel, msgSel, what) => {
    const m = await page.$eval(msgSel, el => (el.hidden ? "" : el.textContent));
    if (!m) fail(`${panel}: no message for ${what}`); else log(`ok   ${panel}: message for ${what}`);
  };
  await tab(page, "inch4");
  for (const s of [11, 5, 4, 3]) await setChecked(page, `#i4-star-${s}`, false);
  await expectMsg("4-inch", "#i4-msg", "three stars");
  await labelAudit(page, "light/inch4 three stars");
  await click(page, "#i4-reset");
  await click(page, '#i4-coords button[data-v="RA"]');
  for (const p of D.inch4.eclipse) await setChecked(page, `#i4-plate-${p}`, false);
  await expectMsg("4-inch", "#i4-msg", "no eclipse plates");
  await click(page, "#i4-reset");
  await tab(page, "astro");
  for (const p of D.astro.plates) await setChecked(page, `#ia-plate-${p}`, false);
  await expectMsg("astrographic", "#ia-msg", "no plates");
  await setChecked(page, "#ia-plate-15", true);
  await click(page, '#ia-mode button[data-v="per_plate"]');
  await expectMsg("astrographic", "#ia-msg", "plate 15 alone with e free (4 measures, 4 unknowns)");
  await click(page, "#ia-reset");
  await click(page, '#ia-mode button[data-v="per_plate"]');
  await labelAudit(page, "light/astro e per plate");
  await click(page, "#ia-reset");
  await tab(page, "principe");
  for (const c of Object.keys(D.principe.comparisons)) await setChecked(page, `#ip-comp-${c}`, false);
  await expectMsg("Principe", "#ip-msg", "no comparisons");
  await setChecked(page, "#ip-comp-W-I2", true);
  await click(page, '#ip-scale button[data-v="free"]');
  const one = await value(page, "#ip-limb");
  if (!Number.isFinite(one)) fail("Principe: no result from W-I2 alone with the scale free"); else log(`ok   Principe: W-I2 alone, scale free: ${one.toFixed(3)}″`);
  await labelAudit(page, "light/principe one comparison");
  await click(page, "#ip-reset");
  if (D.wallal) {
    await tab(page, "wallal");
    for (let i = 0; i < D.wallal.pairs.length; i++) await setChecked(page, `#iw-pair-${i}`, false);
    await expectMsg("Wallal", "#iw-msg", "no plate pairs");
    await click(page, "#iw-reset");
    await setRange(page, "#iw-rmin", await page.$eval("#iw-rmin", el => el.max));
    await setChecked(page, "#iw-power", true);
    await setChecked(page, "#iw-noscale", true);
    const k = await value(page, "#iw-k");
    const why = await page.$eval("#iw-msg", el => (el.hidden ? "" : el.textContent));
    if (Number.isFinite(k)) log(`ok   Wallal: largest cut, power free, no scale terms: beta = ${k.toFixed(3)}`);
    else if (why) log(`ok   Wallal: largest cut, power free, no scale terms: message "${why.slice(0, 60)}..."`);
    else fail("Wallal: largest cut with the power free: neither a result nor a message");
    await setChecked(page, "#iw-noscale", false);
    const k2 = await value(page, "#iw-k"), why2 = await page.$eval("#iw-msg", el => (el.hidden ? "" : el.textContent));
    if (!Number.isFinite(k2) && !why2) fail("Wallal: largest cut, power free, with scale terms: neither a result nor a message");
    await setChecked(page, "#iw-power", false);
    const k3 = await value(page, "#iw-k");
    if (!Number.isFinite(k3)) fail("Wallal: no 1/r result at the largest distance cut");
    else log(`ok   Wallal: largest cut, n = 1: beta = ${k3.toFixed(3)} ± ${(await value(page, "#iw-k-err")).toFixed(3)}`);
    await labelAudit(page, "light/wallal largest cut");
    await page.screenshot({path: path.join(OUT, `${TAG}light_wallal_largest_cut.png`)});
    await click(page, "#iw-reset");
  }
}

async function shots(page, name, full) {
  for (const t of ["inch4", "astro", "principe", "wallal"]) {
    await tab(page, t);
    await page.screenshot({path: path.join(OUT, `${TAG}${name}_${t}.png`), fullPage: full});
  }
  await tab(page, "inch4");
}

async function run() {
  const browser = await chromium.launch();
  const url = "file://" + PAGE;
  const D = JSON.parse(await (async () => {
    const html = fs.readFileSync(PAGE, "utf8");
    const m = html.match(/const D = (\{.*?\});\n/s);
    if (!m) throw new Error("no data block in the page");
    return m[1];
  })());
  log(`page: ${PAGE}`);
  for (const [name, opts] of [["light", {colorScheme: "light", viewport: {width: 1280, height: 1000}}],
                              ["dark", {colorScheme: "dark", viewport: {width: 1280, height: 1000}}],
                              ["phone", {colorScheme: "light", viewport: {width: 390, height: 844}, deviceScaleFactor: 2, isMobile: true, hasTouch: true}]]) {
    const ctx = await browser.newContext(opts);
    const page = await ctx.newPage();
    const errors = [];
    page.on("pageerror", e => errors.push("pageerror: " + e.message));
    page.on("console", m => { if (m.type() === "error" || m.type() === "warning") errors.push(`console ${m.type()}: ${m.text()}`); });
    await page.goto(url);
    await page.waitForFunction(() => window.__plateLab && window.__plateLab.ready, null, {timeout: 10000});
    await page.waitForTimeout(150);
    log(`--- ${name} (${await page.evaluate(() => document.compatMode)}) ---`);
    if (name === "light") {
      await numbers(page, D);
      await jackknifeNotes(page, D);
      await edgeCases(page, D);
      await tab(page, "inch4");
      await click(page, "#i4-reset");
    }
    await shots(page, name, true);
    for (const t of ["inch4", "astro", "principe", "wallal"]) {
      await tab(page, t);
      await labelAudit(page, `${name}/${t}`);
      if (name === "phone") await noSideScroll(page, `${name}/${t}`);
    }
    if (name === "light") {
      // hover a star on the 4-inch chart and the astrographic line: tooltips appear
      await tab(page, "inch4");
      const box = await page.locator("#i4-chart").boundingBox();
      const hit = await page.evaluate(() => {
        const svg = document.getElementById("i4-chart"), c = svg.querySelector("circle.inst-fill");
        const b = c.getBoundingClientRect();
        return {x: b.left + b.width / 2, y: b.top + b.height / 2};
      });
      await page.mouse.move(hit.x, hit.y);
      await page.waitForTimeout(120);
      const tipShown = await page.$eval("#panel-inch4 .tooltip", el => !el.hidden && el.textContent.length > 0);
      if (!tipShown) fail("light: no tooltip over a 4-inch star"); else log("ok   light: tooltip over a 4-inch star");
      await page.screenshot({path: path.join(OUT, `${TAG}light_inch4_tooltip.png`), clip: {x: box.x - 10, y: box.y - 10, width: box.width + 20, height: box.height + 20}});
      await page.mouse.move(5, 5);
      // theme buttons on the standalone page
      if (await page.$("[data-theme-set=dark]")) {
        await click(page, "[data-theme-set=dark]");
        const bg = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
        if (bg !== "rgb(13, 13, 13)") fail(`light: the Dark button did not switch the theme (body ${bg})`); else log("ok   light: Dark button switches the theme");
        await page.screenshot({path: path.join(OUT, `${TAG}light_dark_button.png`)});
        await click(page, "[data-theme-set=dark]");
        const bg2 = await page.evaluate(() => getComputedStyle(document.body).backgroundColor);
        if (bg2 !== "rgb(249, 249, 247)") fail(`light: pressing Dark again did not return to the system theme (body ${bg2})`);
      }
      // keyboard: arrow keys move between tabs
      await page.focus("#tab-inch4");
      await page.keyboard.press("ArrowRight");
      const sel = await page.$eval('[role="tab"][aria-selected="true"]', b => b.id);
      if (sel !== "tab-astro") fail(`light: ArrowRight from the first tab selected ${sel}`); else log("ok   light: arrow keys move between tabs");
      // a non-default astrographic state for the record
      await setRange(page, "#ia-e", "0.02");
      await page.screenshot({path: path.join(OUT, `${TAG}light_astro_fixed_0.020.png`), fullPage: true});
    }
    if (errors.length) for (const e of errors) fail(`${name}: ${e}`); else log(`ok   ${name}: no console or page errors`);
    await ctx.close();
  }

  // the artifact fragment inside a viewer-like skeleton (standards mode, theme from data-theme or the OS)
  if (!pageArg && fs.existsSync(ARTIFACT)) {
    const wrapped = path.join(OUT, "artifact_wrapped.html");
    fs.writeFileSync(wrapped, '<!doctype html>\n<html lang="en"><head><meta charset="utf-8">' +
      '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">' +
      "<style>:root{color-scheme:light;padding:env(safe-area-inset-top) 0 env(safe-area-inset-bottom)}" +
      "body{margin:0;font:14px system-ui,sans-serif}img{max-width:100%}[hidden]{display:none!important}</style></head><body>\n" +
      fs.readFileSync(ARTIFACT, "utf8") + "\n</body></html>\n");
    for (const [name, opts, theme] of [["artifact_light", {colorScheme: "light", viewport: {width: 1100, height: 900}}, null],
                                       ["artifact_dark_attr", {colorScheme: "light", viewport: {width: 1100, height: 900}}, "dark"],
                                       ["artifact_phone_dark", {colorScheme: "dark", viewport: {width: 390, height: 844}, deviceScaleFactor: 2, isMobile: true}, null]]) {
      const ctx = await browser.newContext(opts);
      const page = await ctx.newPage();
      const errors = [];
      page.on("pageerror", e => errors.push("pageerror: " + e.message));
      page.on("console", m => { if (m.type() === "error" || m.type() === "warning") errors.push(`console ${m.type()}: ${m.text()}`); });
      await page.goto("file://" + wrapped);
      await page.waitForFunction(() => window.__plateLab && window.__plateLab.ready, null, {timeout: 10000});
      if (theme) await page.evaluate(t => document.documentElement.setAttribute("data-theme", t), theme);
      await page.waitForTimeout(120);
      const info = await page.evaluate(() => ({mode: document.compatMode, buttons: document.querySelectorAll("[data-theme-set]").length,
        bg: getComputedStyle(document.body).backgroundColor, a: document.getElementById("i4-alpha50").dataset.value}));
      log(`--- ${name} (${info.mode}) ---`);
      if (info.buttons) fail(`${name}: theme buttons present in the artifact fragment`);
      const wantDark = theme === "dark" || name.includes("dark");
      if (wantDark !== (info.bg === "rgb(13, 13, 13)")) fail(`${name}: body background ${info.bg} does not match the theme`);
      compare(`${name}: 4-inch alpha50`, parseFloat(info.a), D.inch4.reference.both.alpha50, 2e-5);
      await page.screenshot({path: path.join(OUT, `${name}.png`), fullPage: true});
      for (const t of ["inch4", "astro", "principe", "wallal"]) {
        await tab(page, t);
        await labelAudit(page, `${name}/${t}`);
        if (name.includes("phone")) await noSideScroll(page, `${name}/${t}`);
      }
      if (errors.length) for (const e of errors) fail(`${name}: ${e}`); else log(`ok   ${name}: no console or page errors`);
      await ctx.close();
    }
  }
  await browser.close();
  log(failures.length ? `\n${failures.length} FAILURE(S)` : "\nALL CHECKS PASSED");
  fs.writeFileSync(path.join(OUT, `${TAG}check_log.txt`), lines.join("\n") + "\n");
  process.exit(failures.length ? 1 : 0);
}
run().catch(e => { console.error(e); process.exit(2); });
