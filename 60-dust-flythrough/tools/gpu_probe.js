// Which WebGL renderer does headless Chromium give us, and how fast is one tour frame?
//   node tools/gpu_probe.js <page.html> <out.png> [flagset]
const path = require("path");
const PW = process.env.PLAYWRIGHT_DIR || path.join(process.env.HOME, ".npm/_npx/e41f203b7505f1fb/node_modules/playwright");
const { chromium } = require(PW);

const FLAGSETS = {
  default: [],
  vulkan: ["--use-angle=vulkan", "--enable-features=Vulkan", "--ignore-gpu-blocklist", "--enable-gpu"],
  egl: ["--use-gl=angle", "--use-angle=gl-egl", "--ignore-gpu-blocklist", "--enable-gpu"],
  swiftshader: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader"],
};

(async () => {
  const [page_html, out_png, flagset = "default"] = process.argv.slice(2);
  const browser = await chromium.launch({ headless: true, args: FLAGSETS[flagset] });
  const page = await browser.newPage({ viewport: { width: 1280, height: 720 } });
  page.on("pageerror", e => console.log("pageerror", e.message));
  page.on("console", m => { if (m.type() === "error") console.log("console error", m.text()); });
  await page.goto("file://" + path.resolve(page_html));
  await page.waitForFunction(() => window.__ready === true, null, { timeout: 180000 });
  const info = await page.evaluate(() => {
    const gl = document.getElementById("gl").getContext("webgl2");
    const ext = gl.getExtension("WEBGL_debug_renderer_info");
    return ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
  });
  console.log("flagset", flagset, "renderer:", info);
  for (const t of [3, 30]) {
    const ms = await page.evaluate(t => {
      const t0 = performance.now();
      window.__renderTourFrame(t, 1280, 720, 1);
      const gl = document.getElementById("gl").getContext("webgl2");
      const px = new Uint8Array(4); gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
      return performance.now() - t0;
    }, t);
    console.log(`tour t=${t}s frame ${ms.toFixed(0)} ms`);
  }
  await page.screenshot({ path: out_png });
  await browser.close();
})().catch(e => { console.error("FAILED", e.message); process.exit(1); });
