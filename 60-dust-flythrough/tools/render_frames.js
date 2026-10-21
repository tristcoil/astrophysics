// Render tour frames of the viewer in headless Chromium on the GPU (ANGLE + Vulkan).
//   node tools/render_frames.js <page.html> <outdir> --mode 0 --times 3,30,60 [--size 1280x720]
//   node tools/render_frames.js <page.html> <outdir> --mode 0 --fps 30 [--from 0 --to 100]
// Frames are written as <outdir>/frame_00000.png (movie) or <outdir>/t<time>_m<mode>.png (stills).
const fs = require("fs");
const path = require("path");
const PW = process.env.PLAYWRIGHT_DIR || path.join(process.env.HOME, ".npm/_npx/e41f203b7505f1fb/node_modules/playwright");
const { chromium } = require(PW);

function arg(name, def) {
  const i = process.argv.indexOf("--" + name);
  return i > 0 ? process.argv[i + 1] : def;
}

(async () => {
  const [pageHtml, outDir] = process.argv.slice(2);
  const mode = +arg("mode", "0");
  const [w, h] = arg("size", "1280x720").split("x").map(Number);
  const step = +arg("step", "1");
  const layers = arg("layers", null);
  fs.mkdirSync(outDir, { recursive: true });
  const browser = await chromium.launch({
    headless: true,
    args: ["--use-angle=vulkan", "--enable-features=Vulkan", "--ignore-gpu-blocklist", "--enable-gpu"],
  });
  const page = await browser.newPage({ viewport: { width: w, height: h }, deviceScaleFactor: 1,
                                       colorScheme: mode === 1 ? "light" : "dark" });
  page.on("pageerror", e => console.log("pageerror", e.message));
  page.on("console", m => { if (m.type() === "error") console.log("console error", m.text()); });
  await page.goto("file://" + path.resolve(pageHtml));
  await page.waitForFunction(() => window.__ready === true, null, { timeout: 300000 });
  await page.evaluate(m => window.__setMode(m), mode);
  await page.evaluate(e => window.__setExposure(e), +arg("exposure", "0"));
  const renderer = await page.evaluate(() => {
    const gl = document.getElementById("gl").getContext("webgl2");
    const ext = gl.getExtension("WEBGL_debug_renderer_info");
    return ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : "unknown";
  });
  console.log("renderer:", renderer);
  if (!/AMD|Radeon|NVIDIA|Intel/.test(renderer) && !process.argv.includes("--allow-software")) {
    throw new Error("no GPU renderer (" + renderer + "); pass --allow-software to render anyway (slow)");
  }
  let times;
  if (arg("times", null)) times = arg("times").split(",").map(Number);
  else {
    const fps = +arg("fps", "30"), t0 = +arg("from", "0"), t1 = +arg("to", "100");
    times = []; for (let i = 0; t0 + i / fps <= t1 + 1e-9; i++) times.push(t0 + i / fps);
  }
  const movie = !arg("times", null);
  const t0 = Date.now();
  for (let i = 0; i < times.length; i++) {
    const t = times[i];
    await page.evaluate(([t, w, h, step, layers]) => {
      window.__renderTourFrame(t, w, h, step);
      if (layers) window.__setLayers(JSON.parse(layers));
      const gl = document.getElementById("gl").getContext("webgl2");
      const px = new Uint8Array(4); gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
    }, [t, w, h, step, layers]);
    const jpeg = process.argv.includes("--jpeg");
    const ext = jpeg ? "jpg" : "png";
    const name = movie ? `frame_${String(i).padStart(5, "0")}.${ext}` : `t${t}_m${mode}.${ext}`;
    await page.screenshot(jpeg ? { path: path.join(outDir, name), type: "jpeg", quality: 95 }
                               : { path: path.join(outDir, name) });
    if (i % 100 === 0) console.log(`frame ${i + 1}/${times.length} (t=${t.toFixed(2)} s) ${((Date.now() - t0) / 1000).toFixed(0)} s`);
  }
  await browser.close();
  console.log("done", times.length, "frames in", ((Date.now() - t0) / 1000).toFixed(0), "s");
})().catch(e => { console.error("FAILED", e.message); process.exit(1); });
