// Render the explorer at desktop and phone widths in light and dark, report JS errors
// and horizontal overflow, click the controls once, and save screenshots.
const path = require("path");
// Playwright: from PLAYWRIGHT_DIR, a normal install (npm i playwright), or the npx cache.
function loadPlaywright() {
  const tries = [process.env.PLAYWRIGHT_DIR, "playwright",
    path.join(process.env.HOME || "", ".npm/_npx/e41f203b7505f1fb/node_modules/playwright")];
  for (const p of tries) {
    if (!p) continue;
    try { return require(p); } catch (e) { /* try the next */ }
  }
  throw new Error("Playwright not found: run `npm i playwright` or set PLAYWRIGHT_DIR");
}
const { chromium } = loadPlaywright();
(async () => {
  const file = "file://" + path.resolve(__dirname, "../figures/neutrino_mass_explorer.html");
  const out = path.resolve(__dirname, "../data/logs/shots");
  const browser = await chromium.launch();
  for (const [w, h, scheme] of [[1200, 1500, "light"], [400, 2400, "dark"]]) {
    const ctx = await browser.newContext({ viewport: { width: w, height: h }, colorScheme: scheme });
    const page = await ctx.newPage();
    const errors = [];
    page.on("pageerror", (e) => errors.push(String(e)));
    page.on("console", (m) => { if (m.type() === "error") errors.push(m.text()); });
    await page.goto(file);
    await page.waitForTimeout(800);
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
    await page.screenshot({ path: `${out}/explorer_${w}_${scheme}.png`, fullPage: true });
    // The fit's verdict in its three states: normal ordering at 0, 15 and 40 meV, inverted at 0.
    const verdicts = [];
    for (const [ord, v] of [["#ordNO", "0"], ["#ordNO", "544"], ["#ordNO", "650"], ["#ordIO", "0"]]) {
      await page.click(ord);
      await page.fill("#mslider", v);
      await page.dispatchEvent("#mslider", "input");
      await page.waitForTimeout(100);
      verdicts.push(await page.evaluate(() => [document.getElementById("vsum").textContent,
        document.getElementById("tcos").className, document.getElementById("scos").textContent]));
    }
    await page.click("#ordIO");
    await page.fill("#mslider", "700");
    await page.dispatchEvent("#mslider", "input");
    await page.waitForTimeout(300);
    const tiles = await page.evaluate(() => ["vsum", "vcos", "vkat", "vbb"].map((id) => document.getElementById(id).textContent));
    await page.click("#pkat");
    const kat = await page.evaluate(() => [document.getElementById("mread").textContent, document.getElementById("vkat").textContent]);
    await page.screenshot({ path: `${out}/explorer_${w}_${scheme}_io.png`, fullPage: true });
    console.log(JSON.stringify({ w, scheme, errors, overflow, tiles, kat, verdicts }, null, 1));
    await ctx.close();
  }
  await browser.close();
})();
