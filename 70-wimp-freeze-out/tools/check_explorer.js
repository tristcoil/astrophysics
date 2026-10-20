// Drive the offline explorer page in headless Chromium: run its in-browser solver on the test
// cases given on stdin (JSON list of [m, sigv, g, wave]), report Omega h^2 for each, any console
// or page errors, and whether anything overflows the page width; optionally take screenshots.
//
//   node tools/check_explorer.js <page.html> [screenshot-dir] < cases.json
//
// Playwright is resolved from PLAYWRIGHT_PATH (default: the npx cache on this machine).
const path = require("path");
const fs = require("fs");
const pw = require(process.env.PLAYWRIGHT_PATH || path.join(process.env.HOME, ".npm/_npx/e41f203b7505f1fb/node_modules/playwright"));

(async () => {
  const page_file = path.resolve(process.argv[2]);
  const shots = process.argv[3] ? path.resolve(process.argv[3]) : null;
  const cases = JSON.parse(fs.readFileSync(0, "utf8"));
  const browser = await pw.chromium.launch();
  const out = { errors: [], omegas: [], overflow: {}, thermal_100: null };
  for (const [name, width, scheme] of [["desktop-light", 1280, "light"], ["phone-dark", 400, "dark"],
                                         ["desktop-dark", 1280, "dark"], ["phone-light", 400, "light"]]) {
    const ctx = await browser.newContext({ viewport: { width, height: 1000 }, colorScheme: scheme });
    const page = await ctx.newPage();
    page.on("pageerror", (e) => out.errors.push(`${name}: ${e.message}`));
    page.on("console", (m) => { if (m.type() === "error") out.errors.push(`${name}: console: ${m.text()}`); });
    await page.goto("file://" + page_file, { waitUntil: "load" });
    await page.waitForTimeout(400);
    out.overflow[name] = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
    if (name === "desktop-light") {
      out.omegas = await page.evaluate((cs) => cs.map(([m, sv, g, w]) => window.__freezeout.relic(m, sv, g, w, false).omega), cases);
      out.thermal_100 = await page.evaluate(() => window.__freezeout.thermalValue(100, 2, "s"));
      out.readout = await page.evaluate(() => ({ omega: document.getElementById("omega").textContent,
        verdict: document.getElementById("verdict").textContent, checks: document.getElementById("checks").textContent }));
    }
    if (shots) await page.screenshot({ path: path.join(shots, `explorer_${name}.png`), fullPage: true });
    await ctx.close();
  }
  await browser.close();
  process.stdout.write(JSON.stringify(out));
})();
