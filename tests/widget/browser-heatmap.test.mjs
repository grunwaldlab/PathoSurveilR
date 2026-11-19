import test from "node:test";
import assert from "node:assert/strict";
import { execFileSync } from "node:child_process";
import { existsSync, readFileSync, unlinkSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";

const root = resolve(new URL("../..", import.meta.url).pathname);
const bundlePath = join(root, "inst/htmlwidgets/psHeatmap.js");
const stylesheetPath = join(root, "inst/htmlwidgets/psHeatmap.css");
const browserPath = [
  process.env.CHROME_BIN,
  "/usr/bin/google-chrome",
  "/usr/bin/chromium",
  "/usr/bin/chromium-browser",
].find((path) => path && existsSync(path));
const rscriptPath = [process.env.RSCRIPT_BIN, "/usr/bin/Rscript", "/usr/local/bin/Rscript"]
  .find((path) => path && existsSync(path));

function fixturePayload(size = 93) {
  const matrix = Array.from({ length: size }, () => Array(size).fill(null));
  for (let row = 0; row < size; row += 1) {
    matrix[row][row] = 100;
    for (let column = row + 1; column < size; column += 1) {
      const low = Math.min(row, column);
      const high = Math.max(row, column);
      const value = 45 + ((low * 131 + high * 73 + low * high * 17) % 55);
      matrix[row][column] = value;
      matrix[column][row] = value;
    }
  }
  const labels = Array.from({ length: size }, (_, index) => `label-${String(index).padStart(3, "0")}`);
  return {
    labels,
    display_labels: labels,
    full_labels: labels,
    matrix,
    sample_ids: [],
    metric: "ani",
    metric_label: "Estimated ANI similarity (%)",
  };
}

function reportedPayload() {
  const fixturePath = resolve(
    root,
    "../pathogensurveillance/test_profiles_out/test_bacteria/report_group_data/all_inputs",
  );
  if (!rscriptPath || !existsSync(fixturePath)) {
    return { source: "synthetic", payload: fixturePayload() };
  }
  const rRoot = JSON.stringify(root);
  const rFixture = JSON.stringify(fixturePath);
  const expression = `
    suppressPackageStartupMessages(devtools::load_all(${rRoot}, quiet = TRUE))
    widget <- estimated_ani_heatmap(${rFixture}, interactive = TRUE)
    cat(as.character(jsonlite::toJSON(widget$x, auto_unbox = TRUE, null = "null", na = "null")))
  `;
  const output = execFileSync(rscriptPath, ["--vanilla", "-e", expression], {
    encoding: "utf8",
    timeout: 120000,
    maxBuffer: 20 * 1024 * 1024,
  });
  return { source: "reported-93x93", payload: JSON.parse(output) };
}

function browserProbe() {
  return `
    try {
    const require = (condition, message) => {
      if (!condition) throw new Error(message);
    };
    const close = (actual, expected, message) => {
      require(Math.abs(actual - expected) < 0.02, \`\${message}: \${actual} != \${expected}\`);
    };
    const matrixTransform = (svg) => {
      const grid = svg.querySelector(".ani-heatmap-grid");
      return grid.transform.baseVal.consolidate().matrix;
    };
    const rowCenters = (svg) => new Map(
      [...svg.querySelectorAll(".ani-heatmap-row-labels text")].map((label) => [
        Number(label.getAttribute("data-row-index")),
        Number(label.getAttribute("y")) - Number(label.getAttribute("font-size")) * 0.35,
      ]),
    );
    const columnCenters = (svg) => new Map(
      [...svg.querySelectorAll(".ani-heatmap-column-labels text")].map((label) => {
        const match = label.getAttribute("transform").match(/^translate\\(([-+\\d.eE]+),/);
        return [Number(label.getAttribute("data-column-index")), Number(match[1])];
      }),
    );
    const leafYs = (path, treeWidth) => {
      const number = "([-+]?(?:\\\\d*\\\\.)?\\\\d+(?:[eE][-+]?\\\\d+)?)";
      const matches = path.getAttribute("d").matchAll(new RegExp(\`M\${number},\${number}H\${number}\`, "g"));
      return [...matches]
        .filter((match) => Math.abs(Number(match[3]) - treeWidth) < 0.02)
        .map((match) => Number(match[2]));
    };
    const leafXs = (path) => {
      const number = "([-+]?(?:\\\\d*\\\\.)?\\\\d+(?:[eE][-+]?\\\\d+)?)";
      return [...path.getAttribute("d").matchAll(new RegExp(\`M\${number},\${number}V\${number}\`, "g"))]
        .map((match) => Number(match[1]));
    };
    const assertAlignment = (svg, stage) => {
      const plotLeft = Number(svg.getAttribute("data-plot-left"));
      const plotTop = Number(svg.getAttribute("data-plot-top"));
      const plotSide = Number(svg.getAttribute("data-plot-side"));
      const rowTreeWidth = Number(svg.getAttribute("data-row-tree-space"));
      const grid = matrixTransform(svg);
      const rows = rowCenters(svg);
      const columns = columnCenters(svg);
      const rowBranches = leafYs(svg.querySelector(".ani-heatmap-row-dendrogram path"), rowTreeWidth);
      const columnBranches = leafXs(svg.querySelector(".ani-heatmap-column-dendrogram path"));

      for (const [index, labelCenter] of rows) {
        const cell = svg.querySelector(\`rect.ani-heatmap-cell[data-row="\${index}"][data-column="0"]\`);
        const matrixCenter = grid.d * (Number(cell.getAttribute("y")) + Number(cell.getAttribute("height")) / 2) + grid.f;
        close(labelCenter, matrixCenter, \`\${stage}: row label/cell alignment for leaf \${index}\`);
        require(rowBranches.some((center) => Math.abs(center - matrixCenter) < 0.02),
          \`\${stage}: row dendrogram leaf \${index} is not aligned\`);
      }
      for (const [index, labelCenter] of columns) {
        const cell = svg.querySelector(\`rect.ani-heatmap-cell[data-row="0"][data-column="\${index}"]\`);
        const matrixCenter = grid.a * (Number(cell.getAttribute("x")) + Number(cell.getAttribute("width")) / 2) + grid.e;
        close(labelCenter, matrixCenter, \`\${stage}: column label/cell alignment for leaf \${index}\`);
        require(columnBranches.some((center) => Math.abs(center - matrixCenter) < 0.02),
          \`\${stage}: column dendrogram leaf \${index} is not aligned\`);
      }
      return { visibleRows: rows.size, visibleColumns: columns.size };
    };
    const focus = (svg) => {
      const matrix = matrixTransform(svg);
      const plotLeft = Number(svg.getAttribute("data-plot-left"));
      const plotTop = Number(svg.getAttribute("data-plot-top"));
      const plotSide = Number(svg.getAttribute("data-plot-side"));
      return {
        scale: matrix.a,
        x: ((plotLeft + plotSide / 2 - matrix.e) / matrix.a - plotLeft) / plotSide,
        y: ((plotTop + plotSide / 2 - matrix.f) / matrix.d - plotTop) / plotSide,
      };
    };
    const firstHost = document.querySelector("#first-host");
    const firstElement = firstHost.querySelector(".psHeatmap");
    const binding = window.__heatmapDefinition.factory(firstElement);
    binding.renderValue(window.__payload);
    let svg = firstElement.querySelector("svg");
    require(svg?.getAttribute("data-dendrograms") === "true", "large fixture should show both trees");
    require(svg.querySelector(".ani-heatmap-row-dendrogram path")?.getAttribute("d").length > 0,
      "row dendrogram should render");
    require(svg.querySelector(".ani-heatmap-column-dendrogram path")?.getAttribute("d").length > 0,
      "column dendrogram should render");
    const initialAlignment = assertAlignment(svg, "initial");

    const rowPath = svg.querySelector(".ani-heatmap-row-dendrogram path");
    const columnPath = svg.querySelector(".ani-heatmap-column-dendrogram path");
    const rowBeforeZoom = rowPath.getAttribute("d");
    const columnBeforeZoom = columnPath.getAttribute("d");
    const slider = firstHost.querySelector(".ani-heatmap-zoom-slider");
    slider.value = "350";
    slider.dispatchEvent(new Event("input", { bubbles: true }));
    svg = firstElement.querySelector("svg");
    require(rowPath.getAttribute("d") !== rowBeforeZoom, "row tree should track zoom");
    require(columnPath.getAttribute("d") !== columnBeforeZoom, "column tree should track zoom");
    const zoomAlignment = assertAlignment(svg, "zoomed");

    const rowBeforePan = rowPath.getAttribute("d");
    const columnBeforePan = columnPath.getAttribute("d");
    const focusBeforePan = focus(svg);
    const interaction = svg.querySelector(".ani-heatmap-interaction");
    const bounds = interaction.getBoundingClientRect();
    const startX = bounds.left + bounds.width / 2;
    const startY = bounds.top + bounds.height / 2;
    interaction.dispatchEvent(new MouseEvent("mousedown", {
      bubbles: true, view: window, clientX: startX, clientY: startY, button: 0, buttons: 1,
    }));
    window.dispatchEvent(new MouseEvent("mousemove", {
      bubbles: true, view: window, clientX: startX + 35, clientY: startY + 27, button: 0, buttons: 1,
    }));
    window.dispatchEvent(new MouseEvent("mouseup", {
      bubbles: true, view: window, clientX: startX + 35, clientY: startY + 27, button: 0, buttons: 0,
    }));
    svg = firstElement.querySelector("svg");
    require(rowPath.getAttribute("d") !== rowBeforePan, "row tree should track pan");
    require(columnPath.getAttribute("d") !== columnBeforePan, "column tree should track pan");
    const focusAfterPan = focus(svg);
    require(
      Math.abs(focusAfterPan.x - focusBeforePan.x) > 0.001
        || Math.abs(focusAfterPan.y - focusBeforePan.y) > 0.001,
      "synthetic drag should pan the matrix",
    );
    const panAlignment = assertAlignment(svg, "panned");
    const beforeResize = focus(svg);

    firstHost.style.width = "1200px";
    firstHost.style.height = "900px";
    firstElement.style.width = "100%";
    firstElement.style.height = "100%";
    binding.resize(1200, 900);
    svg = firstElement.querySelector("svg");
    const afterResize = focus(svg);
    close(afterResize.scale, beforeResize.scale, "resize should preserve zoom scale");
    close(afterResize.x, beforeResize.x, "resize should preserve horizontal zoom focus");
    close(afterResize.y, beforeResize.y, "resize should preserve vertical zoom focus");
    assertAlignment(svg, "resized");

    const smallHost = document.querySelector("#small-host");
    const smallElement = smallHost.querySelector(".psHeatmap");
    const smallBinding = window.__heatmapDefinition.factory(smallElement);
    smallBinding.renderValue(window.__payload);
    const smallSvg = smallElement.querySelector("svg");
    require(smallSvg?.getAttribute("data-dendrograms") === "false",
      \`small view should hide both trees (widget \${smallElement.clientWidth}x\${smallElement.clientHeight}, plot \${smallSvg?.getAttribute("data-plot-side")})\`);

    document.body.setAttribute("data-browser-result", btoa(JSON.stringify({
      ok: true,
      initialAlignment,
      zoomAlignment,
      panAlignment,
      zoomScale: afterResize.scale,
      resizeFocusPreserved: true,
      smallTreesHidden: true,
    })));
    } catch (error) {
      document.body.setAttribute("data-browser-result", btoa(JSON.stringify({
        ok: false,
        error: String(error?.stack ?? error),
      })));
    }
  `;
}

test("browser keeps dendrogram leaves aligned through zoom, pan, and resize", {
  skip: !browserPath && "Chrome/Chromium is unavailable",
}, () => {
  const bundle = readFileSync(bundlePath, "utf8").replaceAll("</script", "<\\/script");
  const stylesheet = readFileSync(stylesheetPath, "utf8");
  const { source, payload: fixture } = reportedPayload();
  assert.equal(fixture.labels.length, 93, "browser alignment fixture must be 93x93");
  const payload = JSON.stringify(fixture);
  const html = `<!doctype html>
    <html><head><meta charset="utf-8"><style>
      html, body { margin: 0; width: 100%; height: 100%; }
      .host { position: relative; }
      #first-host { width: 1300px; height: 1000px; }
      #small-host { width: 480px; height: 360px; }
      .psHeatmap { width: 100%; height: 100%; }
      ${stylesheet}
    </style></head><body>
      <div id="first-host" class="host"><div class="psHeatmap"></div></div>
      <div id="small-host" class="host"><div class="psHeatmap"></div></div>
      <script>window.HTMLWidgets = { widget: (definition) => { window.__heatmapDefinition = definition; } };</script>
      <script>${bundle}</script>
      <script>window.__payload = ${payload}; ${browserProbe()}</script>
    </body></html>`;
  const htmlPath = join(tmpdir(), `pathosurv-heatmap-browser-${process.pid}.html`);
  writeFileSync(htmlPath, html);
  try {
    const dom = execFileSync(browserPath, [
      "--headless",
      "--no-sandbox",
      "--disable-gpu",
      "--allow-file-access-from-files",
      "--hide-scrollbars",
      "--window-size=1600,1200",
      "--virtual-time-budget=1500",
      "--dump-dom",
      `file://${htmlPath}`,
    ], { encoding: "utf8", timeout: 30000, maxBuffer: 20 * 1024 * 1024 });
    const encodedResult = dom.match(/data-browser-result="([^"]+)"/)?.[1];
    assert.ok(encodedResult, "browser did not report the dendrogram assertions");
    const result = JSON.parse(Buffer.from(encodedResult, "base64").toString("utf8"));
    assert.equal(result.ok, true, result.error ?? "browser assertions failed");
    assert.ok(result.initialAlignment.visibleRows > 0);
    assert.ok(result.zoomAlignment.visibleRows > 0);
    assert.ok(result.panAlignment.visibleRows > 0);
    assert.ok(result.resizeFocusPreserved);
    assert.ok(result.smallTreesHidden);
    console.log(`browser alignment check (${source}): ${JSON.stringify(result)}`);
  } finally {
    unlinkSync(htmlPath);
  }
});
