// Load the viewer like a person would, click its controls, report errors and take screenshots.
//   node tools/ui_check.js <page.html> <outdir>
const path = require("path");
const fs = require("fs");
const PW = process.env.PLAYWRIGHT_DIR || path.join(process.env.HOME, ".npm/_npx/e41f203b7505f1fb/node_modules/playwright");
const { chromium } = require(PW);

(async () => {
  const [pageHtml, outDir] = process.argv.slice(2);
  fs.mkdirSync(outDir, { recursive: true });
  const browser = await chromium.launch({ headless: true,
    args: ["--use-angle=vulkan", "--enable-features=Vulkan", "--ignore-gpu-blocklist", "--enable-gpu"] });
  const report = [];
  for (const [name, vp, scheme] of [["desktop_dark", { width: 1366, height: 820 }, "dark"],
                                    ["desktop_light", { width: 1366, height: 820 }, "light"],
                                    ["phone_dark", { width: 400, height: 820 }, "dark"]]) {
    const ctx = await browser.newContext({ viewport: vp, colorScheme: scheme, deviceScaleFactor: 1 });
    const page = await ctx.newPage();
    const errors = [];
    page.on("pageerror", e => errors.push("pageerror: " + e.message));
    page.on("console", m => { if (m.type() === "error") errors.push("console: " + m.text()); });
    await page.goto("file://" + path.resolve(pageHtml));
    await page.waitForFunction(() => window.__ready === true, null, { timeout: 300000 });
    await page.waitForTimeout(800);
    await page.screenshot({ path: path.join(outDir, `${name}_start.png`) });
    const res = { name };
    res.mode = await page.evaluate(() => [...document.querySelectorAll(".seg button")].find(b => b.getAttribute("aria-pressed") === "true")?.textContent);
    if (name.startsWith("phone")) {
      await page.click("#panelToggle");
      await page.waitForTimeout(300);
      res.panelOpen = await page.evaluate(() => !document.getElementById("panel").classList.contains("closed"));
      res.bodyScrollWidth = await page.evaluate(() => document.documentElement.scrollWidth);
    }
    await page.click("#modeGlow");
    res.glowPressed = await page.evaluate(() => document.getElementById("modeGlow").getAttribute("aria-pressed"));
    await page.click("[data-goto=wave]");
    await page.waitForTimeout(2200);
    res.readoutAfterWave = await page.evaluate(() => document.getElementById("readout").textContent.trim());
    await page.screenshot({ path: path.join(outDir, `${name}_wave.png`) });
    await page.click("#tourBtn");
    await page.waitForTimeout(2500);
    res.tourTime = await page.evaluate(() => document.getElementById("tourTime").textContent);
    res.caption = await page.evaluate(() => document.getElementById("caption").textContent.trim().slice(0, 60));
    await page.click("#tourBtn");
    await page.click("#lyrShell");
    res.shellChecked = await page.evaluate(() => document.getElementById("lyrShell").checked);
    // drag to orbit
    const box = await page.locator("#gl").boundingBox();
    await page.mouse.move(box.x + box.width / 2, box.y + box.height / 2);
    await page.mouse.down(); await page.mouse.move(box.x + box.width / 2 + 120, box.y + box.height / 2 + 40, { steps: 8 }); await page.mouse.up();
    await page.mouse.wheel(0, -400);
    await page.waitForTimeout(600);
    await page.screenshot({ path: path.join(outDir, `${name}_after.png`) });
    res.errors = errors;
    report.push(res);
    await ctx.close();
  }
  console.log(JSON.stringify(report, null, 1));
  await browser.close();
})().catch(e => { console.error("FAILED", e.message); process.exit(1); });
