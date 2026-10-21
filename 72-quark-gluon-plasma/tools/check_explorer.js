// Drive figures/freezeout_thermometer.html in headless Chromium:
// JS model vs Python yields, the Fit button, console errors, screenshots.
// usage: node tools/check_explorer.js <playwright module dir> <outdir>
const path = require("path");
const { chromium } = require(process.argv[2]);
const outdir = process.argv[3];
const page_url = "file://" + path.resolve(__dirname, "..", "figures", "freezeout_thermometer.html");
(async () => {
  const browser = await chromium.launch();
  const report = { errors: [], datasets: [] };
  for (const [name, opts] of [["light", { viewport: { width: 1280, height: 1500 }, colorScheme: "light" }],
                              ["dark", { viewport: { width: 1280, height: 1500 }, colorScheme: "dark" }],
                              ["phone", { viewport: { width: 400, height: 1900 }, colorScheme: "light" }]]) {
    const ctx = await browser.newContext(opts);
    const page = await ctx.newPage();
    page.on("pageerror", e => report.errors.push(name + ": " + e.message));
    page.on("console", m => { if (m.type() === "error" && !/fonts\.g/.test(m.text())) report.errors.push(name + " console: " + m.text()); });
    await page.goto(page_url);
    await page.waitForTimeout(1500);
    if (name === "light") {
      report.datasets = await page.evaluate(() => {
        const P = JSON.parse(document.getElementById("payload").textContent);
        return P.datasets.map((d, i) => {
          const js = window.__model(i, d.best.T, d.best.muB, d.best.V);
          let worst = 0;
          js.forEach((v, o) => { worst = Math.max(worst, Math.abs(v / d.python_model[o] - 1)); });
          return { name: d.name, worst };
        });
      });
      await page.selectOption("#ds", "5");
      await page.click("#fit");
      await page.waitForFunction(() => /Best fit/.test(document.getElementById("busy").textContent), null, { timeout: 60000 });
      report.fit_star77 = await page.textContent("#busy");
      await page.selectOption("#ds", "0");
      await page.fill("#T", "140");
      await page.dispatchEvent("#T", "input");
      await page.click("#fit");
      await page.waitForFunction(() => /Best fit/.test(document.getElementById("busy").textContent) && !document.getElementById("fit").disabled, null, { timeout: 60000 });
      report.fit_alice_from_140 = await page.textContent("#busy");
      await page.click("#reset");
      report.overflow = await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth);
    } else {
      report["overflow_" + name] = await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth);
    }
    await page.screenshot({ path: path.join(outdir, `explorer_${name}.png`), fullPage: true });
    await ctx.close();
  }
  await browser.close();
  console.log(JSON.stringify(report, null, 1));
})();
