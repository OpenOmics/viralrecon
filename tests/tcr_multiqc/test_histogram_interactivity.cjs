/* Controller checks with a mocked DOM/Plotly API; not a browser rendering test. */
"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const script = fs.readFileSync(path.join(__dirname, "../../assets/tcr_pod5/multiqc/histogram_interactivity.js"), "utf8");
const tick = () => new Promise(resolve => setImmediate(resolve));
const plain = value => JSON.parse(JSON.stringify(value));

function setup(mode = "ok", readyState = "complete") {
  const spec = {points: [[0, 3], [500, 2], [1000, 0], [10000, 1]], bin_width: 500,
    title: "Read length", xlab: "Read length (nt)", xmin: 0, xmax: 10500,
    focus_upper: 10000, filename: "test_histogram"};
  const chart = {style: {}, clientWidth: 800};
  const fallback = {hidden: false};
  const controls = {hidden: true};
  const status = {textContent: ""};
  const nodes = {".tcr-hist-chart": chart, ".tcr-hist-fallback": fallback,
    ".tcr-hist-controls": controls, ".tcr-hist-status": status,
    ".tcr-hist-data": {textContent: JSON.stringify(spec)}};
  const actions = {};
  const buttons = ["focus", "full", "reset", "linear", "log", "svg"].map(action => ({
    dataset: {tcrAction: action}, addEventListener: (event, fn) => { actions[action] = fn; }
  }));
  const root = {dataset: {}, querySelector: sel => nodes[sel], querySelectorAll: () => buttons};
  const calls = {newPlot: [], relayout: [], download: []};
  const Plotly = {
    newPlot: (...args) => {
      calls.newPlot.push(args);
      if (mode === "throw") throw Error("unavailable");
      return mode === "reject" ? Promise.reject(Error("unavailable")) : Promise.resolve();
    },
    relayout: (...args) => { calls.relayout.push(args); return Promise.resolve(); },
    downloadImage: (...args) => { calls.download.push(args); return Promise.resolve(); }
  };
  const listeners = {};
  const document = {currentScript: {closest: () => root}, readyState,
    addEventListener: (event, fn) => { listeners[event] = fn; }};
  vm.runInNewContext(script, {document, window: mode === "absent" ? {} : {Plotly}, Promise});
  return {chart, fallback, controls, status, calls, actions, listeners};
}

(async function () {
  const c = setup();
  await tick();
  assert.equal(c.fallback.hidden, true);
  assert.equal(c.controls.hidden, false);
  const [, traces, layout, config] = c.calls.newPlot[0];
  assert.deepEqual(plain(traces[0].x), [250, 750, 1250, 10250]);
  assert.deepEqual(plain(traces[0].y), [3, 2, 0, 1]);
  assert.equal(traces[0].width, 500);
  assert.deepEqual(plain(traces[0].customdata[3]), [10000, 10500]);
  assert.match(traces[0].hovertemplate, /Reads:/);
  assert.deepEqual(plain(layout.xaxis.range), [0, 10000]);
  assert.equal(config.toImageButtonOptions.format, "svg");
  c.actions.full();
  assert.deepEqual(plain(c.calls.relayout.at(-1)[1]["xaxis.range"]), [0, 10500]);
  c.actions.log();
  assert.equal(c.calls.relayout.at(-1)[1]["yaxis.type"], "log");
  c.actions.linear();
  assert.equal(c.calls.relayout.at(-1)[1]["yaxis.type"], "linear");
  c.actions.focus();
  assert.deepEqual(plain(c.calls.relayout.at(-1)[1]["xaxis.range"]), [0, 10000]);
  c.actions.reset();
  assert.equal(c.calls.relayout.at(-1)[1]["yaxis.type"], "linear");
  c.actions.svg();
  assert.equal(c.calls.download[0][1].format, "svg");
  for (const mode of ["absent", "reject", "throw"]) {
    const f = setup(mode);
    await tick();
    assert.equal(f.fallback.hidden, false);
    assert.equal(f.controls.hidden, true);
    assert.equal(f.chart.style.display, "none");
    assert.match(f.status.textContent, /SVG preview/);
  }
  const deferred = setup("ok", "loading");
  assert.equal(deferred.calls.newPlot.length, 0);
  deferred.listeners.DOMContentLoaded();
  await tick();
  assert.equal(deferred.calls.newPlot.length, 1);
  deferred.listeners.DOMContentLoaded();
  assert.equal(deferred.calls.newPlot.length, 1);
  process.stdout.write("Interactive histogram controller checks passed (mocked DOM/Plotly).\n");
})().catch(error => { console.error(error); process.exitCode = 1; });
