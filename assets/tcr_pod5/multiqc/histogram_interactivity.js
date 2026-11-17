/* Embedded into MultiQC custom HTML. Uses the report's bundled Plotly runtime. */
(function () {
  "use strict";
  const root = document.currentScript.closest(".tcr-histogram");
  if (!root) return;

  function init() {
    if (root.dataset.initialized) return;
    root.dataset.initialized = "true";
    const chart = root.querySelector(".tcr-hist-chart");
    const fallback = root.querySelector(".tcr-hist-fallback");
    const controls = root.querySelector(".tcr-hist-controls");
    const status = root.querySelector(".tcr-hist-status");
    const spec = JSON.parse(root.querySelector(".tcr-hist-data").textContent);
    const plotly = window.Plotly;
    function fail() {
      chart.style.display = "none";
      fallback.hidden = false;
      controls.hidden = true;
      status.textContent = "Interactive plotting is unavailable in this report template; showing the SVG preview.";
    }
    if (!plotly || typeof plotly.newPlot !== "function") {
      fail();
      return;
    }
    const fullRange = [spec.xmin, spec.xmax];
    const initialRange = [spec.xmin, spec.focus_upper === null ? spec.xmax : Math.min(spec.xmax, spec.focus_upper)];
    const trace = {
      type: "bar",
      x: spec.points.map(p => p[0] + spec.bin_width / 2),
      y: spec.points.map(p => p[1]),
      width: spec.bin_width,
      customdata: spec.points.map(p => [p[0], p[0] + spec.bin_width]),
      marker: {color: "#397fb8", line: {color: "white", width: 0.4}},
      hovertemplate: "Bin: [%{customdata[0]:,.0f}, %{customdata[1]:,.0f})<br>Reads: %{y:,.0f}<extra></extra>",
      name: "Reads"
    };
    const layout = {
      title: {text: spec.title, x: 0.06, xanchor: "left", font: {size: 17}},
      height: 430, autosize: true, showlegend: false, bargap: 0,
      dragmode: "zoom", hovermode: "closest",
      paper_bgcolor: "white", plot_bgcolor: "white",
      font: {family: "Arial, sans-serif", size: 12, color: "#364152"},
      margin: {l: 72, r: 25, t: 65, b: 64},
      xaxis: {type: "linear", title: {text: spec.xlab}, range: initialRange,
              autorange: false, tickformat: ",.0f", gridcolor: "#e6ebf0", zeroline: false},
      yaxis: {type: "linear", title: {text: "Reads"}, rangemode: "tozero", autorange: true,
              tickformat: ",.0f", gridcolor: "#e6ebf0", zeroline: false}
    };
    const options = {
      responsive: true, displayModeBar: true, displaylogo: false, scrollZoom: false,
      modeBarButtonsToRemove: ["select2d", "lasso2d"],
      toImageButtonOptions: {format: "svg", filename: spec.filename, width: 1200, height: 450}
    };
    chart.style.display = "block";
    let pending;
    try {
      pending = plotly.newPlot(chart, [trace], layout, options);
    } catch (error) {
      fail();
      return;
    }
    Promise.resolve(pending).then(function () {
      fallback.hidden = true;
      controls.hidden = false;
      status.textContent = "Hover for bin counts. Drag to zoom; use the plot toolbar to pan. All bins are retained in every view.";
      function update(settings) {
        Promise.resolve(plotly.relayout(chart, settings)).catch(function () {
          status.textContent = "That view change failed. The SVG preview and exact bin counts remain available in the report output directory.";
        });
      }
      root.querySelectorAll("[data-tcr-action]").forEach(function (button) {
        button.addEventListener("click", function () {
          const action = button.dataset.tcrAction;
          if (action === "full" || action === "focus" || action === "reset") {
            const change = {"xaxis.range": action === "full" ? fullRange : initialRange,
                            "xaxis.autorange": false, "yaxis.range": null, "yaxis.autorange": true};
            if (action === "reset") {
              change["yaxis.type"] = "linear";
              change["yaxis.title.text"] = "Reads";
            }
            update(change);
          } else if (action === "linear" || action === "log") {
            update({"yaxis.type": action, "yaxis.range": null, "yaxis.autorange": true,
                    "yaxis.title.text": action === "log" ? "Reads (log scale)" : "Reads"});
          } else if (action === "svg") {
            Promise.resolve(plotly.downloadImage(chart, {format: "svg", filename: spec.filename,
                                                         width: 1200, height: 450})).catch(function () {
              status.textContent = "Image export was blocked by the browser. SVG copies are also saved in the report plots directory.";
            });
          }
        });
      });
      // The native report can resize its content when the navigation/toolbox opens.
      if (window.ResizeObserver && plotly.Plots && plotly.Plots.resize) {
        let previousWidth = chart.clientWidth;
        const observer = new ResizeObserver(function () {
          if (chart.clientWidth > 0 && chart.clientWidth !== previousWidth) {
            previousWidth = chart.clientWidth;
            plotly.Plots.resize(chart);
          }
        });
        observer.observe(chart);
      }
    }).catch(fail);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init, {once: true});
  } else {
    init();
  }
})();
