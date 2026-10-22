#!/usr/bin/env python3
"""Build a self-contained MultiQC config from viralrecon TCR POD5 outputs.

Python 3.9+ and Matplotlib. Reads published JSON/TSV/FASTA, never BAM/POD5.
This is reporting only: no reclassification, trimming, filtering, or consensus calls.
"""
import argparse
import base64
import csv
import hashlib
import html
import io
import json
import math
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

VERSION = "1.0.3"
GROUPS = {
    "00": "Overview across steps",
    "01": "POD5 input audit",
    "02": "Primer, barcode and reference preparation",
    "03": "Dorado basecalling",
    "04": "Two-ended demultiplexing",
    "05": "Primer trimming",
    "06": "Alignment and competitive construct assignment",
    "07": "Batch target inference",
    "08": "Per-well mixture and consensus screening",
    "09": "Dorado polishing and final validation",
    "10": "Diagnostic consensus and delivered sequences",
    "11": "Recorded settings and report completeness",
    "12": "Exploratory: assignment comparison",
    "13": "Exploratory: barcode margin comparison",
    "14": "Exploratory: read-end audit",
    "15": "Exploratory: parent-tracked splitting comparison",
    "16": "Exploratory: internal primer junctions",
}
LABELS = {
    "raw_assigned_reads": "Raw assigned reads", "trim_pass": "Trimmed reads",
    "usable_reads": "Assigned to panel (pre-polish)", "selected_reads": "Selected reads (pre-polish)",
    "plus": "+ strand reads", "minus": "− strand reads", "callable_pct": "Callable % (pre-polish)",
    "status": "Final status", "diagnostic_kind": "Diagnostic kind", "diagnostic_pct": "Diagnostic callable %",
    "dominant_construct": "Dominant construct", "polish_selected_reads": "Selected reads (post-polish)",
    "post_callable_pct": "Unmasked % (post-polish)", "masked_bases": "Masked bases",
    "diagnostic_consensuses": "All diagnostic consensuses", "review_consensuses": "Polished diagnostic subset",
    "assigned_reads": "Panel-assigned reads", "unresolved_for_qc": "Unresolved for QC",
    "ambiguous_reference_reads": "Ambiguous reference reads", "both_ends": "Both ends required",
    "raw_identity": "Raw alignment identity", "deletion_aware_identity": "Deletion-aware identity",
    "resolved_reasons": "Resolved review reasons", "variation_reassessment": "Original variation reassessment",
    "split_alignment_reads": "Split reads", "split_alignment_pct": "Split reads (%)",
    "split_alignment_total_reads": "Reads evaluated for split fraction",
}
HELP = {
    "raw_assigned_reads": "Reads assigned to this well by production Dorado demultiplexing, before primer/QC filters.",
    "usable_reads": "Pre-polish reads assigned to any reference; not necessarily selected for the dominant target.",
    "selected_reads": "Dominant-construct target-supported reads before polishing; not an independent template count.",
    "callable_pct": "Pre-polish callable fraction × 100, even when the displayed final status comes from polishing.",
    "post_callable_pct": "100 × (1 − masked_bases / length) from the final polishing decision.",
    "diagnostic_pct": "Callable percentage of the diagnostic sequence actually selected by the pipeline report.",
    "review_consensuses": "Subset of all diagnostic sequences; do not add these two counts together.",
    "unresolved_for_qc": "Unresolved reads used by the mixture/QC gate; ambiguous-reference reads are treated separately.",
    "resolved_reasons": "Historical review reasons resolved by recorded post-polish evidence; other acceptance checks still apply.",
    "variation_reassessment": "RESOLVED, RETAINED or NOT_APPLICABLE for the original variation flag. Blank means an older report did not record this check.",
    "split_alignment_reads": "Read IDs classified as split alignments at this stage; excluded from consensus even below the review thresholds.",
    "split_alignment_pct": "100 × split reads / all read IDs evaluated at this stage, including unassigned reads. Blank in older reports without split_alignment_qc.",
    "split_alignment_total_reads": "Stage-specific denominator: all trimmed read IDs at initial QC; the selected input cohort realigned to the polished sequence at final QC. Supplementary records do not add reads.",
}
FILTER_NAMES = {
    "primer_errors": "maximum full-primer edits", "primer_window": "primer search window (nt)",
    "max_overtrim": "maximum overtrim (nt)", "max_trim_fail_fraction": "trim-failure warning fraction",
    "min_read_q": "minimum mean read Q", "min_identity": "minimum deletion-aware identity",
    "min_query_coverage": "minimum query coverage", "min_mapq": "minimum mapping quality",
    "score_gap": "minimum alternative-reference alignment-score gap", "deletion_min_length": "minimum exempt internal deletion length (nt)",
    "deletion_min_anchor": "minimum aligned anchor on each deletion side (nt)",
    "min_reads": "minimum selected reads for acceptance", "min_strand_reads": "minimum reads on each strand",
    "min_depth": "minimum per-position depth", "min_baseq": "minimum base Q",
    "consensus_fraction": "minimum base-support fraction", "min_callable_fraction": "minimum callable target fraction",
    "minor_reads": "minimum minor-construct reads", "minor_fraction": "minimum minor-construct fraction",
    "minor_each_strand": "minimum minor-construct reads on each strand", "alert_reads": "minor-construct / split-read alert count",
    "alert_fraction": "minor-construct / split-read alert fraction", "phase_spacing": "minimum phased-site spacing (nt)",
    "max_unresolved_fraction": "maximum unresolved fraction for QC", "min_target_coverage": "minimum target coverage",
    "polish_min_reads": "minimum selected reads to attempt polishing", "polish_min_callable_fraction": "minimum callable fraction to attempt polishing",
    "diagnostic_min_reads": "minimum reads for an unpolished diagnostic", "diagnostic_min_strand_reads": "minimum diagnostic reads on each strand",
    "diagnostic_min_depth": "minimum diagnostic per-position depth", "diagnostic_min_callable_fraction": "minimum diagnostic callable fraction",
    "target_min_alignment_length": "minimum aligned length for target inference (nt)", "target_min_query_coverage": "minimum query coverage for target inference",
    "target_min_reads": "minimum target-inference supporting reads", "target_min_reads_per_well": "minimum supporting reads per eligible well",
    "target_min_wells": "minimum independent supporting wells", "target_min_each_strand": "minimum target-inference reads on each strand",
    "target_min_cluster_fraction": "minimum well-balanced endpoint cluster fraction", "target_boundary_tolerance": "endpoint tolerance (nt)",
}


def scalar(v):
    if v is None or v == "":
        return None
    if not isinstance(v, str):
        return v
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    try:
        f = float(v)
        return (int(f) if f.is_integer() else f) if math.isfinite(f) else None
    except ValueError:
        return v


def pct(v):
    return None if v is None or v == "" else round(100 * float(v), 3)


def ratio(a, b):
    return None if a is None or b in (None, 0) else round(100 * a / b, 3)


def natural(s):
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", str(s))]


def yaml_lines(value, depth=0):
    """Small YAML 1.2 emitter: JSON-quoted keys/scalars avoid implicit type traps."""
    pad = " " * depth
    if isinstance(value, dict):
        for key, val in value.items():
            prefix = pad + json.dumps(str(key), ensure_ascii=False) + ":"
            if isinstance(val, (dict, list)) and val:
                yield prefix
                yield from yaml_lines(val, depth + 2)
            else:
                yield prefix + " " + json.dumps(val, ensure_ascii=False, allow_nan=False)
    elif isinstance(value, list):
        for val in value:
            if isinstance(val, (dict, list)) and val:
                yield pad + "-"
                yield from yaml_lines(val, depth + 2)
            else:
                yield pad + "- " + json.dumps(val, ensure_ascii=False, allow_nan=False)


def clean(row):
    return {k: (html.escape(v) if isinstance(v, str) else v) for k, v in row.items() if v is not None}


def flat(obj, prefix=""):
    for k, v in obj.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            yield from flat(v, key)
        elif isinstance(v, list):
            if not any(isinstance(x, (dict, list)) for x in v):
                yield key, "; ".join(map(str, v))
        else:
            yield key, v


def histogram_points(counts, width):
    """Numeric order, explicit zero-count bins, and exact original bin boundaries."""
    ordered = sorted((float(x), int(n)) for x, n in counts.items())
    if not ordered:
        return []
    if width <= 0 or any(n < 0 or not math.isfinite(x) for x, n in ordered):
        raise ValueError("Invalid histogram bin or count")
    return ordered


def histogram_svg(points, width, title, xlab, upper=None, log_counts=False):
    """Render numeric bins directly, avoiding text-key ordering in line-plot parsers."""
    try:
        import matplotlib
        from matplotlib.figure import Figure
        from matplotlib.ticker import FuncFormatter, MaxNLocator
    except ImportError as e:
        raise ValueError("Histogram rendering requires Matplotlib. Install it in this Python environment with: python3 -m pip install matplotlib") from e
    visible = [(x, n) for x, n in points if upper is None or x < upper]
    xmin = min(0, points[0][0])
    xmax = points[-1][0] + width
    if upper is not None:
        xmax = min(xmax, upper)
    with matplotlib.rc_context({"svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 11,
                               "axes.titleweight": "bold", "axes.labelcolor": "#364152", "text.color": "#1f2937",
                               "xtick.color": "#526174", "ytick.color": "#526174", "svg.hashsalt": "tcr-multiqc-v1.0.1"}):
        fig = Figure(figsize=(11.8, 4.0), facecolor="white", layout="constrained")
        ax = fig.subplots()
        ax.bar([x for x, n in visible], [n for x, n in visible], width=width, align="edge",
               color="#397fb8", edgecolor="white", linewidth=0.5, zorder=3)
        ax.set_xlim(xmin, xmax)
        ax.set_xlabel(xlab, labelpad=10)
        ax.set_ylabel("Reads (log scale)" if log_counts else "Reads", labelpad=8)
        ax.set_title(title, loc="left", pad=14, fontsize=13)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=10, integer=True))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:,.0f}"))
        ax.grid(axis="y", color="#e6ebf0", linewidth=.8, zorder=0)
        ax.set_axisbelow(True)
        if log_counts and any(n > 0 for x, n in visible):
            ax.set_yscale("log")
            ax.set_ylim(bottom=.8)
        else:
            ax.set_ylim(bottom=0)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:,.0f}" if x >= 1 or x == 0 else ""))
        ax.spines[["top", "right"]].set_visible(False)
        for side in ("bottom", "left"):
            ax.spines[side].set_color("#c6ced7")
        buf = io.StringIO()
        fig.savefig(buf, format="svg", metadata={"Date": None})
        return buf.getvalue()


class Report:
    def __init__(self, args):
        self.args = args
        self.root = Path(args.results_dir).expanduser().resolve()
        if (self.root / "tcr_pod5").is_dir():
            self.root /= "tcr_pod5"
        if not self.root.is_dir():
            raise ValueError(f"Results directory not found: {self.root}")
        self.sources, self.availability, self.warnings = {}, {}, []
        self.data, self.order, self.tables, self.group_descriptions = {}, [], {}, {}
        self.settings, self.wells, self.batches = {}, {}, {}
        self.plot_assets = {}

    def track(self, path):
        path = Path(path).resolve()
        if str(path) not in self.sources:
            h = hashlib.sha256()
            with path.open("rb") as f:
                for b in iter(lambda: f.read(1024 * 1024), b""):
                    h.update(b)
            self.sources[str(path)] = {"bytes": path.stat().st_size, "sha256": h.hexdigest()}

    def read_json(self, path, required=False):
        path = Path(path)
        if not path.is_file():
            if required:
                raise ValueError(f"Missing required file: {path}")
            return None
        self.track(path)
        try:
            return json.loads(path.read_text())
        except Exception as e:
            raise ValueError(f"Invalid JSON in {path}: {e}") from e

    def tsv(self, path):
        path = Path(path)
        if not path.is_file():
            return []
        self.track(path)
        with path.open(newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def add(self, group, key, title, rows, description="", kind="table", headers=None, **pconfig):
        sid = f"tcr_{group}_{key}"
        parent = f"tcr_{group}"
        if parent not in self.order:
            self.order.append(parent)
        if group not in self.group_descriptions and description:
            self.group_descriptions[group] = description
        if kind != "html":
            rows = {k: v for k, v in rows.items() if any(x is not None for x in v.values())}
        if not rows:
            rows = {"Not available": {"availability": "No matching output was found."}}
            kind = "table"
        pc = {"id": sid + "_plot", "title": GROUPS[group] + ": " + title}
        if kind == "table":
            pc.update(col1_header="Batch / well / item", sort_rows=False, only_defined_headers=True, no_violin=True)
            keys = list(dict.fromkeys(k for row in rows.values() for k in row))
            hs = {k: {"title": LABELS.get(k, k.replace("_", " ").capitalize()),
                      "description": HELP.get(k, LABELS.get(k, k.replace("_", " "))),
                      "scale": False} for k in keys}
            for k in keys:
                hs[k]["placement"] = keys.index(k)
                values = [r[k] for r in rows.values() if k in r and r[k] is not None]
                if values and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in values):
                    hs[k]["format"] = "{:,.2f}" if any(isinstance(v, float) for v in values) else "{:,.0f}"
            for k, val in (headers or {}).items():
                hs.setdefault(k, {}).update(val)
            self.tables[sid] = rows
        elif kind == "bargraph":
            pc.update(ylab="Reads", cpswitch=True, sort_samples=False)
        pc.update(pconfig)
        section = {"id": sid, "parent_id": parent, "parent_name": group + " · " + GROUPS[group],
                   "parent_description": self.group_descriptions.get(group, ""),
                   "section_name": title, "description": description, "plot_type": kind, "pconfig": pc,
                   "data": rows if kind == "html" else {html.escape(str(k)): clean(v) for k, v in rows.items()}}
        if kind == "table":
            section["headers"] = hs
        elif kind == "bargraph":
            colors = {"CONSENSUS_PASS": "#168477", "REVIEW": "#d69b26", "REJECT": "#c74747", "NO_READS": "#939ca6",
                      "Assigned": "#168477", "Unclassified": "#939ca6", "Passed": "#168477", "Failed": "#c74747",
                      "Callable": "#168477", "Not callable": "#ccd2d9", "+": "#286ca8", "−": "#8f6ac8"}
            section["categories"] = {k: {"name": k, **({"color": colors[k]} if k in colors else {})}
                                      for k in dict.fromkeys(k for row in rows.values() for k in row)}
        self.data[sid] = section

    def add_histogram(self, key, title, series, width, xlab, description, focus_upper=None):
        chunks, bin_rows = [], {}
        script_path = Path(__file__).with_name("histogram_interactivity.js")
        if not script_path.is_file():
            raise ValueError("Missing histogram_interactivity.js. Extract the complete collector bundle.")
        script = script_path.read_text()
        for batch, counts in series.items():
            points = histogram_points(counts, width)
            if not points:
                continue
            total = sum(n for x, n in points)
            tail = sum(n for x, n in points if focus_upper is not None and x >= focus_upper)
            tag = re.sub(r"[^A-Za-z0-9_-]", "_", batch) + "_" + hashlib.sha256(batch.encode()).hexdigest()[:8]
            def picture(svg, suffix, alt):
                name = f"tcr_03_{key}_{tag}_{suffix}.svg"
                self.plot_assets[name] = svg
                encoded = base64.b64encode(svg.encode()).decode("ascii")
                return f'<img alt="{html.escape(alt, quote=True)}" src="data:image/svg+xml;base64,{encoded}" style="width:100%;height:auto;display:block;background:white;" />'
            chunks.append(f'<div class="tcr-histogram" style="margin:0 0 24px;"><p><strong>{html.escape(batch)}</strong> · {total:,} reads</p>')
            spec = {"points": points, "bin_width": width, "title": title, "xlab": xlab,
                    "xmin": min(0, points[0][0]), "xmax": points[-1][0] + width,
                    "focus_upper": focus_upper, "filename": f"tcr_03_{key}_{tag}"}
            payload = json.dumps(spec, ensure_ascii=True).replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
            actions = [("reset", "Reset view"), ("full", "Full range"), ("linear", "Linear counts"), ("log", "Log counts"), ("svg", "Download SVG")]
            if focus_upper is not None:
                actions.insert(0, ("focus", "0–10 kb"))
            chunks.append('<div class="tcr-hist-controls" hidden><div style="display:flex;flex-wrap:wrap;gap:6px;margin-bottom:8px;">')
            for action, label in actions:
                chunks.append(f'<button type="button" class="btn btn-sm btn-default" data-tcr-action="{action}" '
                              f'style="border:1px solid #c6ced7;background:white;color:#364152;">{label}</button>')
            tsv = "bin_start_inclusive\tbin_end_exclusive\treads\n" + "".join(f"{x:g}\t{x+width:g}\t{n}\n" for x, n in points)
            encoded_tsv = base64.b64encode(tsv.encode()).decode("ascii")
            chunks.append(f'<a class="btn btn-sm btn-default" download="{spec["filename"]}_bins.tsv" '
                          f'href="data:text/tab-separated-values;base64,{encoded_tsv}" '
                          'style="border:1px solid #c6ced7;background:white;color:#364152;">Download bin counts</a></div></div>')
            chunks.append('<p class="tcr-hist-status" role="status" style="font-size:12px;color:#526174;">'
                          'SVG preview; interactive controls load when the report opens.</p>')
            chunks.append('<div class="tcr-hist-chart" style="display:none;width:100%;height:430px;"></div>')
            chunks.append('<div class="tcr-hist-fallback">')
            svg = histogram_svg(points, width, title, xlab, upper=focus_upper)
            chunks.append(picture(svg, "main", title + " for " + batch))
            if tail:
                chunks.append('<details style="margin-top:12px;"><summary style="cursor:pointer;">Full read-length range — all reads, logarithmic count axis</summary>')
                svg = histogram_svg(points, width, title + " · full range", xlab, log_counts=True)
                chunks.append(picture(svg, "full_range", "Full read-length histogram for " + batch))
                chunks.append('</details>')
            chunks.append('</div>')
            if tail:
                chunks.append(f'<p>Initial 0–10 kb view: {total-tail:,} reads below {focus_upper:,} nt ({100*(total-tail)/total:.2f}%); '
                              f'{tail:,} reads at or above {focus_upper:,} nt are retained. Use Full range to inspect them.</p>')
            chunks.append(f'<script type="application/json" class="tcr-hist-data">{payload}</script>')
            chunks.append('<script>' + script + '</script></div>')
            for x, n in points:
                start = int(x) if x.is_integer() else x
                bin_rows[f"{batch} / [{start}, {start+width})"] = {
                    "batch": batch, "bin_start_inclusive": start, "bin_end_exclusive": start+width, "reads": n}
        self.add("03", key, title, "".join(chunks), description, kind="html")
        self.tables[f"tcr_03_{key}_bins"] = bin_rows

    def filters(self, keys):
        lines = []
        for batch, (settings, source) in self.settings.items():
            s = "; ".join(f"{FILTER_NAMES.get(k, k)}: {settings.get(k, 'not recorded')}" for k in keys)
            lines.append(f"<p><strong>{html.escape(batch)}</strong>: {html.escape(s)}.</p>")
        return "".join(lines)

    def check_equal(self, label, a, b):
        if a is not None and b is not None and a != b:
            self.warnings.append(f"{label}: {a} != {b}. Check that these outputs are from the same completed run.")

    def load(self):
        self.run = self.read_json(self.root / "summary/run_summary.json") or {}
        summary = self.tsv(self.root / "summary/well_summary.tsv")
        summary_by = {}
        for r in summary:
            k = (r["batch"], r["well"])
            if k in summary_by:
                raise ValueError(f"Duplicate well in well_summary.tsv: {k}")
            summary_by[k] = {k: scalar(v) for k, v in r.items()}
        discovered = {p.name for p in self.root.iterdir() if p.is_dir() and
                      any((p / x).exists() for x in ("basecalling", "demultiplexing", "wells", "scheme", "target_inference"))}
        discovered |= {r[0] for r in summary_by}
        discovered |= set(self.run.get("batch_targets", {}))
        batches = self.args.batch or sorted(discovered)
        if not batches:
            raise ValueError("No TCR batches found. Pass the results directory or its tcr_pod5 subdirectory.")
        for batch in batches:
            if batch not in discovered:
                raise ValueError(f"Batch not found: {batch}")
            bdir = self.root / batch
            manifest = self.read_json(bdir / "demultiplexing/demux/manifest.json")
            scheme = self.read_json(bdir / "scheme/scheme/scheme.json") or {}
            target_file = bdir / "target_inference/targets/inferred_targets.json"
            target = self.read_json(target_file)
            if target is None:
                target = self.run.get("batch_targets", {}).get(batch)
            if target and target.get("batch", batch) != batch:
                raise ValueError(f"Target report batch does not match directory: {batch}")
            settings = (target or {}).get("settings", {})
            source = str(target_file) if target_file.is_file() else "summary/run_summary.json → batch_targets"
            if not settings:
                source = "NOT AVAILABLE — defaults were not substituted"
            self.settings[batch] = (settings, source)
            raw = {}
            for item in manifest or []:
                if item["well"] in raw:
                    raise ValueError(f"Duplicate demux manifest well: {batch}/{item['well']}")
                if item.get("batch", batch) != batch:
                    raise ValueError(f"Wrong batch in {bdir}/demultiplexing/demux/manifest.json")
                raw[item["well"]] = item["reads"]
            wells = set(raw) - {"unclassified"}
            wells |= {x["well"] for x in scheme.get("pairs", [])}
            wells |= {w for b, w in summary_by if b == batch}
            if (bdir / "wells").is_dir():
                wells |= {p.name for p in (bdir / "wells").iterdir() if p.is_dir()}
            b = {"dir": bdir, "manifest": manifest, "raw": raw, "scheme": scheme, "target": target,
                 "base": self.read_json(bdir / "basecalling/basecall_provenance.json"),
                 "demux_settings": self.read_json(bdir / "demultiplexing/demux/settings.json")}
            self.batches[batch] = b
            for well in sorted(wells, key=natural):
                wd = bdir / "wells" / well
                q = self.read_json(wd / "heterogeneity/qc/decision.json")
                pol = self.read_json(wd / "consensus/consensus_result/decision.json")
                trim = self.read_json(wd / "trimming/trim_stats.json")
                s = summary_by.get((batch, well), {})
                if trim is None and q and q.get("trim_stats"):
                    trim = q["trim_stats"]
                raw_n = raw.get(well, s.get("raw_assigned_reads"))
                if raw_n is not None:
                    self.check_equal(f"{batch}/{well} raw reads", raw_n, s.get("raw_assigned_reads"))
                if trim:
                    self.check_equal(f"{batch}/{well} trim input", raw_n, trim.get("input_reads"))
                    if all(k in trim for k in ("input_reads", "passed_reads", "failed_reads")):
                        self.check_equal(f"{batch}/{well} trim accounting", trim["input_reads"], trim["passed_reads"] + trim["failed_reads"])
                row = dict(s)
                row["raw_assigned_reads"] = raw_n
                row["trim_pass"] = (trim or {}).get("passed_reads")
                if q:
                    for dest, src in (("usable_reads", "assigned_reads"), ("selected_reads", "selected_reads"),
                                      ("dominant_construct", "dominant_construct")):
                        if dest in s:
                            self.check_equal(f"{batch}/{well} {dest}", q.get(src), s[dest])
                        row[dest] = q.get(src)
                    row.update(plus=q["strand_counts"].get("+", 0) if "strand_counts" in q else None,
                               minus=q["strand_counts"].get("-", 0) if "strand_counts" in q else None)
                    row["callable_pct"] = pct(q.get("callable_fraction"))
                else:
                    row["callable_pct"] = pct(s.get("callable_fraction"))
                if not row.get("status"):
                    row["status"] = (pol or {}).get("status", "NOT REPORTED" if raw_n != 0 else "NO_READS")
                row["diagnostic_pct"] = pct(s.get("diagnostic_callable_fraction"))
                key = batch + " / " + well
                self.wells[key] = {"batch": batch, "well": well, "dir": wd, "summary": s,
                                   "row": row, "q": q, "pol": pol, "trim": trim}
                self.availability[key] = {
                    "raw_reads": raw_n,
                    "trim": "PRESENT" if trim is not None else "NOT EXPECTED" if raw_n == 0 else "MISSING",
                    "qc": "PRESENT" if q is not None else "NOT EXPECTED" if raw_n == 0 else "MISSING",
                    "polish": "PRESENT" if pol else "MISSING" if (q or {}).get("polish_eligible") else "NOT ELIGIBLE" if q or raw_n == 0 else "UNKNOWN",
                    "final_summary": "PRESENT" if s else "MISSING",
                }

    def overview(self):
        rows = {}
        audit = self.read_json(self.root / "audit/pod5_audit.json") or {}
        for batch, b in self.batches.items():
            wells = [w for w in self.wells.values() if w["batch"] == batch]
            r = {"wells": len(wells), "pod5_reads": audit.get("batch_reads", {}).get(batch), "basecalled_reads": (b["base"] or {}).get("reads")}
            if b["manifest"] is not None:
                r.update(raw_assigned_reads=sum(n for w, n in b["raw"].items() if w != "unclassified"),
                         unclassified_reads=b["raw"].get("unclassified"),
                         assigned_pct=ratio(sum(n for w, n in b["raw"].items() if w != "unclassified"), sum(b["raw"].values())))
                self.check_equal(batch + " basecall/demux accounting", r["basecalled_reads"], sum(b["raw"].values()))
            for col, source, metric in (("trimmed_reads", "trim", "passed_reads"), ("panel_assigned_reads", "q", "assigned_reads"),
                                         ("selected_reads", "q", "selected_reads")):
                relevant = [w for w in wells if w["row"].get("raw_assigned_reads") != 0]
                complete = all(w[source] is not None and metric in w[source] for w in relevant)
                r[col] = sum(w[source][metric] for w in relevant) if complete and wells else None
            complete_summary = wells and all(w["summary"] for w in wells)
            if complete_summary:
                statuses = Counter(w["summary"]["status"] for w in wells)
                r.update(accepted_consensuses=statuses.get("CONSENSUS_PASS", 0), review_wells=statuses.get("REVIEW", 0),
                         rejected_wells=statuses.get("REJECT", 0), no_read_wells=statuses.get("NO_READS", 0),
                         diagnostic_consensuses=sum(bool(w["summary"].get("diagnostic_available")) for w in wells))
            r["summary_complete"] = bool(complete_summary)
            rows[batch] = r
        self.add("00", "batch", "Batch overview", rows,
                 "Production counts only. Downstream batch totals are blank unless all applicable wells have that metric. "
                 "Unclassified reads have no well assignment. Diagnostic sequences are separate from accepted sequences.")
        cols = ["status", "raw_assigned_reads", "trim_pass", "usable_reads", "selected_reads", "plus", "minus", "callable_pct", "dominant_construct", "diagnostic_kind", "diagnostic_pct"]
        self.add("00", "well", "Across-step metrics by well", {k: {c: w["row"].get(c) for c in cols} for k, w in self.wells.items()},
                 "Read left to right to locate losses. Selected reads and callable percentage are pre-polish; status is the final reported decision. "
                 "Counts are reads, not depth at every base or independent starting templates. Empty cells are unavailable.")
        statuses = {}
        for b in self.batches:
            ws = [w for w in self.wells.values() if w["batch"] == b]
            if ws and all(w["summary"] for w in ws):
                statuses[b] = dict(Counter(w["row"]["status"] for w in ws))
        if statuses:
            self.add("00", "outcomes", "Final well outcomes", statuses, kind="bargraph", ylab="Wells")

    def inputs_and_basecalls(self):
        audit = self.read_json(self.root / "audit/pod5_audit.json")
        rows, file_rows = {}, {}
        if audit is not None:
            for b in self.batches:
                files = [f for f in audit.get("files", []) if f.get("batch") == b]
                rows[b] = {"pod5_reads": audit.get("batch_reads", {}).get(b), "files": len(files),
                           "empty_files": sum(x.get("batch") == b for x in audit.get("empty_files", [])),
                           "duplicate_read_entries": sum(b in (x.get("first_batch"), x.get("second_batch")) for x in audit.get("duplicate_reads", [])),
                           "pod5_version": audit.get("pod5_version")}
                for i, f in enumerate(files):
                    file_rows[f"{b} / {i+1} / {Path(f['source']).name}"] = {"reads": f.get("reads"), "bytes": f.get("bytes")}
                self.check_equal(b + " POD5/basecall count", rows[b]["pod5_reads"], (self.batches[b]["base"] or {}).get("reads"))
        self.add("01", "audit", "Input validation", rows,
                 "Checks POD5 readability, read UUID uniqueness, sequencing-kit agreement and a consistent signal sample rate within each batch. "
                 "Zero-read batches fail; empty-file handling is explicit. This table reports the recorded audit, not a new scan of raw signals.")
        if file_rows:
            self.add("01", "files", "Input files", file_rows)
        schemes, refs, scoring = {}, {}, {}
        base_rows, length_hist, q_hist = {}, {}, {}
        for batch, b in self.batches.items():
            pairs = b["scheme"].get("pairs", [])
            if pairs:
                schemes[batch] = {"well_pairs": len(pairs), "distinct_forward_indexes": len({x['forward']['index'] for x in pairs}),
                                  "distinct_reverse_indexes": len({x['reverse']['index'] for x in pairs})}
            fp = b["dir"] / "scheme/scheme/reference.fasta"
            if fp.is_file():
                for name, seq in self.fasta(fp).items():
                    refs[batch + " / " + name] = {"reference_length": len(seq),
                        "expected_construct": name in b["target"].get("expected_constructs", []) if b["target"] else None}
            toml = b["dir"] / "scheme/scheme/barcodes.toml"
            if toml.is_file():
                self.track(toml)
                in_scoring = False
                for line in toml.read_text().splitlines():
                    line = line.strip()
                    if line.startswith("["):
                        in_scoring = line == "[scoring]"
                    elif in_scoring and "=" in line:
                        k, v = line.split("=", 1)
                        scoring[batch + " / " + k.strip()] = {"value": scalar(v.split("#", 1)[0].strip().strip('"'))}
            base = b["base"]
            if base:
                base_rows[batch] = {k: base.get(k) for k in ("reads", "model", "read_splitting", "basecall_quality_filter")}
                base_rows[batch]["dorado_version"] = base.get("dorado_version", "").splitlines()[0]
                base_rows[batch]["command"] = " ".join(base.get("command", []))
            summary_path = b["dir"] / "demultiplexing/demux/dorado_sequencing_summary.txt"
            if summary_path.is_file():
                records = self.tsv(summary_path)
                lengths, qualities = [], []
                seen = set()
                for r in records:
                    rid = r.get("read_id")
                    if rid in seen:
                        raise ValueError(f"Duplicate read in {summary_path}: {rid}")
                    seen.add(rid)
                    ln, q = scalar(r.get("sequence_length_template")), scalar(r.get("mean_qscore_template"))
                    if isinstance(ln, (int, float)):
                        lengths.append(ln)
                    if isinstance(q, (int, float)):
                        qualities.append(q)
                self.check_equal(batch + " basecall/summary count", (base or {}).get("reads"), len(records))
                row = base_rows.setdefault(batch, {})
                row.update(summary_reads=len(records), median_length=statistics.median(lengths) if lengths else None,
                           median_q=statistics.median(qualities) if qualities else None)
                if lengths:
                    counts = Counter(int(x // 500) * 500 for x in lengths)
                    length_hist[batch] = {k: counts[k] for k in range(0, max(counts) + 500, 500)}
                if qualities:
                    counts = Counter(math.floor(x) for x in qualities)
                    q_hist[batch] = {k: counts[k] for k in range(min(counts), max(counts) + 1)}
        self.add("02", "scheme", "Assay definitions", schemes,
                 "Builds the forward/reverse barcode-pair scheme, full-primer definitions and competitive reference panel. "
                 "The same well must be supported by both barcode ends. Reference IDs define possible construct assignments.")
        self.add("02", "refs", "Reference panel", refs)
        self.add("02", "scoring", "Recorded barcode scoring parameters", scoring,
                 "Values come from the published barcodes.toml. These are Dorado scoring parameters; independent exploratory audits use different scoring.")
        self.add("03", "base", "Basecalling and read characteristics", base_rows,
                 "Converts POD5 signal to reads. Production preserves barcode/primer context with --no-trim. "
                 "The recorded command and provenance show the model, splitting setting and basecall quality cutoff. "
                 "Length/Q summaries below include all reads in the production demultiplexing sequencing summary.")
        if length_hist:
            self.add_histogram("length", "Read-length distribution", length_hist, 500, "Read length (nt)",
                     "Histogram with 500-nt bins: [0, 500), [500, 1000), and so on. All reads, including unclassified reads. "
                     "The initial view focuses on 0–10,000 nt; use Full range to inspect all lengths. Hover for exact bin counts, "
                     "drag to zoom, and switch between linear and logarithmic counts. SVG previews are retained as a fallback. "
                     "Exact counts are also exported in tables/tcr_03_length_bins.tsv.", focus_upper=10000)
        if q_hist:
            self.add_histogram("quality", "Mean read quality distribution", q_hist, 1, "Mean read Q score",
                     "Histogram with one-Q-score bins: [Q, Q+1). Bar height is the number of reads in that bin. "
                     "All reads are included; mean read Q is not consensus accuracy. "
                     "Hover for exact bin counts, drag to zoom, or use the toolbar to pan and reset. "
                     "Exact bin counts are exported in tables/tcr_03_quality_bins.tsv.")

    def demux_trim_assignment(self):
        demux, demux_plot, trim_rows, trim_plot, trim_reasons = {}, {}, {}, {}, {}
        assignment_rows, assignment_plot, deletion_rows = {}, {}, {}
        for batch, b in self.batches.items():
            if b["manifest"] is not None:
                assigned = sum(n for w, n in b["raw"].items() if w != "unclassified")
                demux[batch] = {"input_reads": sum(b["raw"].values()), "assigned_reads": assigned,
                               "unclassified_reads": b["raw"].get("unclassified"),
                               "assigned_pct": ratio(assigned, sum(b["raw"].values())),
                               "wells_with_reads": sum(n > 0 for w, n in b["raw"].items() if w != "unclassified"),
                               "both_ends": (b["demux_settings"] or {}).get("both_ends"),
                               "trim_during_demux": (b["demux_settings"] or {}).get("trim")}
                if "unclassified" in b["raw"]:
                    demux_plot[batch] = {"Assigned": assigned, "Unclassified": b["raw"]["unclassified"]}
        for key, w in self.wells.items():
            t, q = w["trim"], w["q"]
            if t is not None:
                trim_rows[key] = {k: t.get(k) for k in ("input_reads", "passed_reads", "failed_reads")}
                trim_rows[key]["failed_pct"] = ratio(t.get("failed_reads"), t.get("input_reads"))
                trim_plot[key] = {"Passed": t.get("passed_reads", 0), "Failed": t.get("failed_reads", 0)}
                trim_reasons[key] = t.get("reasons", {})
            if q:
                assignment_rows[key] = {k: q.get(k) for k in ("total_reads", "assigned_reads", "selected_reads", "ambiguous_reference_reads", "unresolved_for_qc")}
                counts = {**q.get("construct_counts", {}), **q.get("unresolved_counts", {})}
                if counts:
                    assignment_plot[key] = counts
            fp = w["dir"] / "heterogeneity/qc/read_assignments.tsv"
            if fp.is_file():
                records = self.tsv(fp)
                counts = Counter(r["assignment"] for r in records)
                if key in assignment_plot:
                    self.check_equal(key + " assignment accounting", sum(assignment_plot[key].values()), len(records))
                assignment_plot[key] = dict(counts)
                affected = [r for r in records if (scalar(r.get("exempt_deleted_bases")) or 0) > 0]
                deletion_rows[key] = {"evaluated_reads": len(records), "reads_with_exempt_deletions": len(affected),
                                      "exempt_deleted_bases": sum(scalar(r.get("exempt_deleted_bases")) or 0 for r in records),
                                      "median_raw_identity_affected": statistics.median(float(r["raw_identity"]) for r in affected) if affected else None,
                                      "median_deletion_aware_identity_affected": statistics.median(float(r["deletion_aware_identity"]) for r in affected) if affected else None}
        self.add("04", "counts", "Production barcode assignment", demux,
                 "Dorado assigns a well using the custom paired barcode scheme with both ends required. Unclassified reads remain batch-level. "
                 "A score of −1 does not identify the reason for failing assignment; read-end audits are separate diagnostic evidence.")
        if demux_plot:
            self.add("04", "yield", "Assigned versus unclassified", demux_plot, kind="bargraph")
        self.add("04", "well", "Raw reads per well", {k: {"raw_assigned_reads": w["row"].get("raw_assigned_reads")} for k, w in self.wells.items()})
        self.add("05", "counts", "Trimming outcomes", trim_rows,
                 "Confirms full primers at both ends, checks that trimming retains the expected interior sequence, and preserves base qualities. "
                 "Failed reads do not proceed to alignment. Exceeding max_trim_fail_fraction raises a warning; it does not itself reject the whole well."
                 + self.filters(["primer_errors", "primer_window", "max_overtrim", "max_trim_fail_fraction"]))
        if trim_plot:
            self.add("05", "yield", "Reads retained after trimming", trim_plot, kind="bargraph")
        if trim_reasons:
            self.add("05", "reasons", "Recorded trimming reasons", trim_reasons)
        self.add("06", "counts", "Read assignment and selection", assignment_rows,
                 "Aligns clean reads competitively to the reference panel. Requires one primary alignment, no supplementary/SA/spliced alignment, "
                 "sufficient read quality, query coverage, deletion-aware identity, MAPQ and separation from alternative reference scores. "
                 "Long internal deletions with sufficient aligned anchors are excluded from the identity denominator; this is not permission to ignore all deletions. "
                 "Reads assigned to another construct and unresolved reads remain visible in the counts."
                 + self.filters(["min_read_q", "min_identity", "min_query_coverage", "min_mapq", "score_gap", "deletion_min_length", "deletion_min_anchor"]))
        if assignment_plot:
            self.add("06", "outcomes_table", "Assignment categories by well", assignment_plot,
                     "Mutually exclusive per-read assignments/first-failing categories. A read may have additional problems not represented by its first-failing reason.")
            self.add("06", "outcomes", "Assignment outcomes", assignment_plot, kind="bargraph")
        self.add("06", "deletions", "Deletion-aware identity audit", deletion_rows,
                 "Summarizes the published per-read identity audit. 'Affected' means at least one exempt deleted base; medians refer only to those reads. "
                 "This table does not infer how many would change assignments under another method; see the optional assignment comparison.")

    def target_and_qc(self):
        target_rows, target_batch, target_excluded = {}, {}, {}
        for batch, b in self.batches.items():
            t = b["target"]
            if not t:
                continue
            target_batch[batch] = {"input_reads": t.get("input_reads"), "input_wells": t.get("input_wells"), "method": t.get("method")}
            target_excluded[batch] = t.get("excluded_reads", {})
            for ref, r in t.get("references", {}).items():
                row = {k: r.get(k) for k in ("status", "consensus_eligible", "start", "end", "length", "reference_length", "candidate_reads", "candidate_wells", "supporting_reads", "supporting_wells", "reason")}
                row.update(cluster_pct=pct(r.get("balanced_cluster_fraction")), plus=r.get("strand_counts", {}).get("+", 0),
                           minus=r.get("strand_counts", {}).get("-", 0), expected_construct=ref in t.get("expected_constructs", []))
                target_rows[batch + " / " + ref] = row
        self.add("07", "targets", "Inferred target intervals", target_rows,
                 "Infers a shared amplicon interval using long, high-quality reads from multiple wells. Joint start/end clusters are balanced by well, "
                 "with support from both strands. Coordinates are original-reference, 0-based start-inclusive/end-exclusive. "
                 "Non-target panel references may legitimately have no support. These empirical bounds do not prove a complete PCR product; shared truncation/deletion can be missed."
                 + self.filters(["target_min_alignment_length", "target_min_query_coverage", "target_min_reads", "target_min_reads_per_well", "target_min_wells", "target_min_each_strand", "target_min_cluster_fraction", "target_boundary_tolerance"]))
        if target_batch:
            self.add("07", "input", "Target-inference input", target_batch)
        if target_excluded:
            self.add("07", "excluded", "Target-inference exclusions", target_excluded,
                     "These counts belong to the pooled target-inference stage, not the per-well final consensus gate.")
        qc_rows, strands, callable_rows, reasons, blockers = {}, {}, {}, {}, {}
        for key, w in self.wells.items():
            q = w["q"]
            if not q:
                continue
            qc_rows[key] = {"screen_status": q.get("status"), **{k: q.get(k) for k in ("dominant_construct", "selected_reads", "candidate_sites", "strong_candidate_sites", "blocking_variation_sites", "target_partial_reads", "target_boundary_outlier_reads", "polish_eligible")},
                            "split_alignment_reads": q.get("split_alignment_qc", {}).get("reads"),
                            "split_alignment_pct": pct(q.get("split_alignment_qc", {}).get("fraction")),
                            "split_alignment_total_reads": q.get("split_alignment_qc", {}).get("total_reads"),
                            "linked_haplotype_groups": len(q["linked_haplotypes"]) if isinstance(q.get("linked_haplotypes"), list) else None,
                            "callable_pct": pct(q.get("callable_fraction")), "reasons": "; ".join(q.get("reasons", [])), "warnings": "; ".join(q.get("warnings", []))}
            if "strand_counts" in q:
                strands[key] = {"+": q["strand_counts"].get("+", 0), "−": q["strand_counts"].get("-", 0)}
            if q.get("callable_fraction") is not None:
                callable_rows[key] = {"Callable": pct(q["callable_fraction"]), "Not callable": pct(1 - q["callable_fraction"])}
            blockers[key] = {"polish_eligible": q.get("polish_eligible"), "polish_blockers": "; ".join(q.get("polish_blockers", []))}
            for reason in set(q.get("reasons", [])):
                row = reasons.setdefault(w["batch"], Counter())
                row[reason] += 1
        self.add("08", "qc", "Pre-polish screening", qc_rows,
                 "Screens construct mixtures, strand-supported minor alleles, linked variation, target support and consensus callability. "
                 "A strong minor construct needs the minimum count, fraction and support on each strand; smaller signals can trigger REVIEW. "
                 "An ambiguous reference assignment is excluded from consensus and reported separately from unresolved_for_qc. "
                 "From v0.4.1, split-read REVIEW requires both alert_reads and alert_fraction of all read IDs evaluated at this stage. "
                 "All split reads remain excluded from consensus and included in unresolved_for_qc. "
                 "Screening cannot establish that a PCR well started from one physical template."
                 + self.filters(["minor_reads", "minor_fraction", "minor_each_strand", "alert_reads", "alert_fraction", "phase_spacing", "max_unresolved_fraction", "min_target_coverage"]))
        self.add("08", "polish_eligibility", "Polishing eligibility", blockers,
                 "Polishing can be attempted for eligible REVIEW wells, but mixture/structural/target blockers prevent promotion. "
                 "Eligibility does not mean a final consensus pass."
                 + self.filters(["polish_min_reads", "min_strand_reads", "polish_min_callable_fraction"]))
        if strands:
            self.add("08", "strands", "Selected evidence by strand", strands, kind="bargraph")
        if callable_rows:
            self.add("08", "callability", "Pre-polish callability", callable_rows,
                     "Percentage of the inferred target meeting the pre-polish callability rules. This is not a per-base depth trace."
                     + self.filters(["min_depth", "min_baseq", "consensus_fraction", "min_callable_fraction"]), kind="bargraph", ylab="Target positions (%)", cpswitch=False)
        if reasons:
            self.add("08", "reasons", "Wells with each screening reason", reasons,
                     "Reasons are nonexclusive: one well may contribute to several columns. Do not add reason counts to estimate affected wells.")

    def fasta(self, path):
        self.track(path)
        seqs, name, pieces = {}, None, []
        for line in Path(path).read_text().splitlines():
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(pieces)
                name = line[1:].split()[0]
                if name in seqs:
                    raise ValueError(f"Duplicate FASTA identifier {name} in {path}")
                pieces = []
            elif line.strip():
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                pieces.append(line.strip().upper())
        if name is not None:
            if name in seqs:
                raise ValueError(f"Duplicate FASTA identifier {name} in {path}")
            seqs[name] = "".join(pieces)
        return seqs

    def consensus(self):
        polished, diagnostic, final = {}, {}, {}
        for key, w in self.wells.items():
            p, q, s = w["pol"], w["q"], w["summary"]
            if p:
                polished[key] = {k: p.get(k) for k in ("status", "length", "masked_bases", "candidate_sites", "strong_candidate_sites", "blocking_variation_sites", "ambiguous_reference_reads", "unresolved_for_qc")}
                polished[key].update(polish_selected_reads=p.get("selected_reads"),
                                     split_alignment_reads=p.get("split_alignment_qc", {}).get("reads"),
                                     split_alignment_pct=pct(p.get("split_alignment_qc", {}).get("fraction")),
                                     split_alignment_total_reads=p.get("split_alignment_qc", {}).get("total_reads"),
                                     post_callable_pct=ratio(p["length"] - p["masked_bases"], p["length"]) if p.get("length") and p.get("masked_bases") is not None else None,
                                     plus=p["strand_counts"].get("+", 0) if "strand_counts" in p else None,
                                     minus=p["strand_counts"].get("-", 0) if "strand_counts" in p else None,
                                     residual_fixed_indels=len(p["residual_fixed_indels"]) if isinstance(p.get("residual_fixed_indels"), list) else None,
                                     resolved_reasons="; ".join(p.get("resolved_reasons", [])),
                                     variation_reassessment=p.get("variation_reassessment", {}).get("status", ""),
                                     reasons="; ".join(p.get("reasons", [])), warnings="; ".join(p.get("warnings", [])))
            if s or q:
                d = (q or {}).get("diagnostic", {})
                diagnostic[key] = {"status": s.get("status"), "diagnostic_available": s.get("diagnostic_available"),
                                   "diagnostic_kind": s.get("diagnostic_kind"), "diagnostic_pct": pct(s.get("diagnostic_callable_fraction")),
                                   "unpolished_draft_eligible": d.get("eligible"), "unpolished_draft_callable_pct": pct(d.get("callable_fraction")),
                                   "blockers": s.get("diagnostic_blockers", "; ".join(d.get("blockers", [])))}
            if s:
                final[key] = {k: s.get(k) for k in ("status", "reasons", "warnings", "diagnostic_kind")}
        self.add("09", "polished", "Post-polish validation", polished,
                 "Realigns selected reads to a draft, runs Dorado polish, realigns again and repeats the QC checks. "
                 "Low-depth, low-support and problematic positions are masked. Persistent biological blockers remain relevant. "
                 "An earlier SNP variation flag can resolve only when the original reads and supporting bases remain evaluable, "
                 "every affected region is uniquely mapped and unmasked with sufficient support, and blocking variation is absent. "
                 "From v0.4.1, the post-polish split-read flag also requires both alert thresholds, using all read IDs evaluated in this realignment. "
                 "Split reads remain excluded from consensus and included in unresolved_for_qc. "
                 "The table uses the final polishing decision; per-well overview callability remains explicitly pre-polish."
                 + self.filters(["min_reads", "min_strand_reads", "min_depth", "min_baseq", "consensus_fraction", "min_callable_fraction", "alert_reads", "alert_fraction"]))
        self.add("10", "diagnostics", "Masked diagnostic sequences", diagnostic,
                 "Non-rejected REVIEW wells may supply a clearly labelled masked sequence for investigation. "
                 "The final report prefers a polished diagnostic over an eligible unpolished draft and exports at most one diagnostic per well. "
                 "REJECT wells do not contribute. Diagnostic eligibility is not an accepted assay call."
                 + self.filters(["diagnostic_min_reads", "diagnostic_min_strand_reads", "diagnostic_min_depth", "min_baseq", "consensus_fraction", "diagnostic_min_callable_fraction"]))
        self.add("10", "final", "Final decisions and reasons", final)
        run_rows = {}
        if self.run:
            run_rows["Whole run (all batches)"] = {k: self.run.get(k) for k in ("wells", "accepted_consensuses", "diagnostic_consensuses", "review_consensuses", "implementation_version")}
            for status, n in self.run.get("statuses", {}).items():
                run_rows["Whole run (all batches)"][status] = n
        self.add("10", "totals", "Reported output totals", run_rows,
                 "run_summary.json totals cover the whole run, even when --batch selects a subset for detailed reporting. "
                 "review_consensus.fasta is the polished subset of diagnostic_consensus.fasta: their counts must not be added.")
        frows, seqrows, sequences = {}, {}, {}
        for name, kind, metric in (("accepted_consensus.fasta", "ACCEPTED", "accepted_consensuses"),
                                   ("diagnostic_consensus.fasta", "DIAGNOSTIC", "diagnostic_consensuses"),
                                   ("review_consensus.fasta", "POLISHED DIAGNOSTIC SUBSET", "review_consensuses")):
            fp = self.root / "summary" / name
            if not fp.is_file():
                frows[name] = {"availability": "MISSING", "kind": kind}
                continue
            seqs = self.fasta(fp)
            sequences[name] = seqs
            frows[name] = {"availability": "PRESENT", "kind": kind, "sequences": len(seqs), "total_bases": sum(map(len, seqs.values()))}
            self.check_equal(name + " versus run summary", len(seqs), self.run.get(metric))
            for sample, seq in seqs.items():
                if self.args.batch and not any(sample.startswith(b + "__") for b in self.args.batch):
                    continue
                seqrows[kind + " / " + sample] = {"length": len(seq), "n_bases": seq.count("N"), "n_pct": ratio(seq.count("N"), len(seq))}
        acc, dia, rev = (sequences.get(n) for n in ("accepted_consensus.fasta", "diagnostic_consensus.fasta", "review_consensus.fasta"))
        if acc is not None and dia is not None and set(acc) & set(dia):
            self.warnings.append("Accepted and diagnostic FASTAs overlap in sample IDs.")
        if rev is not None and dia is not None and any(s not in dia or seq != dia[s] for s, seq in rev.items()):
            self.warnings.append("Polished diagnostic FASTA is not an identical-sequence subset of the diagnostic FASTA.")
        self.add("10", "files", "Delivered FASTA inventory", frows,
                 "Counts and N masking are measured directly from published FASTAs when present. File absence is not a zero-sequence result.")
        if seqrows:
            self.add("10", "sequence_metrics", "Delivered sequence lengths and masking", seqrows,
                     "Polished diagnostic subset rows also appear in the all-diagnostic set; these rows are deliberately labelled and are not additive.")

    def audits(self):
        paths = {Path(p).expanduser().resolve() for p in self.args.audit_dir}
        if not self.args.no_discover_audits:
            for b in self.batches.values():
                for pattern in ("demux_margin1_check_*", "demux_end_audit_*", "read_splitting_check_*", "read_splitting_check_*/internal_junction_audit_*"):
                    paths.update(p.resolve() for p in b["dir"].glob(pattern) if p.is_dir())
        for path in sorted(paths):
            provenance = self.read_json(path / "run_provenance.json") or {}
            if provenance.get("status") not in (None, "COMPLETE"):
                self.warnings.append(f"Skipped incomplete exploratory audit: {path} ({provenance['status']}).")
                continue
            if provenance.get("batch") and provenance["batch"] not in self.batches:
                self.warnings.append(f"Skipped audit for an unselected batch: {path}")
                continue
            files = [p for p in (path / "comparison_summary.json", path / "audit_summary.json") if p.is_file()]
            if len(files) != 1:
                self.warnings.append(f"No unique supported audit summary in {path}; skipped.")
                continue
            d = self.read_json(files[0], True)
            if "setting_change" in d:
                group, description = "13", "Compares a barcode scoring change while keeping the recorded two-ended scheme. Counts describe assignments only; downstream consensus QC was not rerun. Length/Q cohorts overlap and must not be added."
                scalar_fields = ("reads", "baseline_assigned", "trial_assigned", "baseline_unclassified", "trial_unclassified", "net_assigned_change")
                chart = {c: {k: v.get(k, 0) for k in ("unchanged_assigned", "newly_assigned", "unchanged_unclassified", "lost_assignment", "switched_well")} for c, v in d.get("cohorts", {}).items()}
            elif "totals" in d and "split_output_records" in d["totals"]:
                group, description = "15", "Selected parent reads are basecalled with splitting disabled/enabled and compared with unchanged two-ended barcoding. Child records are not independent starting molecules. Target-like recovery is not accepted consensus evidence."
                scalar_fields = ()
                chart = {c: {"split_parents": v["split_event"], "unsplit_parents": v["parents"] - v["split_event"]} for c, v in d.get("cohorts", {}).items()}
            elif "pair_failure_reasons" in d:
                group, description = "16", "Searches for adjacent inward-facing full-primer/barcode pairs and target-like internal segments. This is a stratified sample, not a batch-wide recovery estimate. Provisional segments require further validation; multiple-well parents require review. Pair failure reasons are nonexclusive."
                scalar_fields = ()
                chart = {c: {"parents_with_proposals": v["parents_with_provisional_segments"], "parents_without_proposals": v["parents"] - v["parents_with_provisional_segments"]} for c, v in d.get("cohorts", {}).items()}
            elif "selected_reads" in d and "cohorts" in d:
                group, description = "14", "Independently inspects primer/barcode ends and reference-alignment structure in a length/Q-selected subset. Diagnostic categories are not Dorado's actual rejection reasons, and target-like span does not prove a correct complete amplicon."
                scalar_fields = ("input_reads", "selected_reads", "excluded_by_length_or_q")
                chart = {c: v.get("end_categories", {}) for c, v in d.get("cohorts", {}).items()}
            else:
                self.warnings.append(f"Unknown audit schema: {files[0]}; skipped.")
                continue
            ident = hashlib.sha256(str(path).encode()).hexdigest()[:10]
            label = (provenance.get("batch") or "Batch not recorded") + " / " + path.name
            rows = {c: {k: v for k, v in x.items() if not isinstance(v, (dict, list))} for c, x in d.get("cohorts", {}).items()}
            intro = description + "<p><strong>EXPLORATORY — excluded from production totals.</strong> " + html.escape(label) + "</p>"
            if scalar_fields:
                self.add(group, ident + "_totals", label + " — totals", {label: {k: d.get(k) for k in scalar_fields}}, intro)
            self.add(group, ident + "_cohorts", label + " — cohort metrics", rows, intro)
            if chart:
                self.add(group, ident + "_plot", label + " — outcomes", chart, description, kind="bargraph", ylab="Parents" if group in ("15", "16") else "Reads")
            for nested in ("end_categories", "target_categories", "agreement_with_recorded_well", "transitions"):
                nested_rows = {c: x[nested] for c, x in d.get("cohorts", {}).items() if nested in x}
                if nested_rows:
                    self.add(group, ident + "_" + nested, label + " — " + nested.replace("_", " "), nested_rows)
            settings = {"summary": {k: d[k] for k in ("selection", "setting_change", "barcode_settings", "pair_failure_reasons") if k in d},
                        "provenance": {k: provenance[k] for k in ("settings", "parameters") if k in provenance}}
            self.add(group, ident + "_settings", label + " — settings and sampling", {k: {"value": v} for k, v in flat(settings)},
                     "Recorded settings and sample sizes. Repeated audit runs are shown separately, never summed. Pair failure counts, when present, are nonexclusive.")
        assignment_files = {Path(p).expanduser().resolve() for p in self.args.assignment_evaluation}
        if not self.args.no_discover_audits:
            assignment_files.update(self.root.rglob("assignment_evaluation.json"))
        for p in sorted(assignment_files):
            d = self.read_json(p, True)
            if d.get("batch") and d["batch"] not in self.batches:
                continue
            ident = hashlib.sha256(str(p).encode()).hexdigest()[:10]
            label = str(d.get("batch", "Batch not recorded")) + " / " + p.parent.name
            rows, outcomes = {}, {}
            for well, v in d.get("wells", {}).items():
                rows[well] = {k: v.get(k) for k in ("deletion_changed_reads", "informative_changed_reads")}
                for method, counts in v.get("counts", {}).items():
                    outcomes[well + " / " + method] = counts
            desc = "Comparison only: evaluates legacy, deletion-aware and informative-position assignments. Informative scoring uses covered canonical variable reference columns; gaps/deletions are not nucleotide identity evidence. The alternative omits MAPQ/score-gap and whole-read identity gates. Changed assignments are not automatically recovered production reads."
            self.add("12", ident + "_summary", label + " — comparison", {label: {"reads": d.get("reads"), "wells": len(d.get("wells", {})),
                      "informative_columns": len(d.get("informative_columns_0based", [])), "excluded_variable_columns": len(d.get("excluded_variable_columns_0based", []))}}, desc)
            self.add("12", ident + "_wells", label + " — changes by well", rows, desc)
            self.add("12", ident + "_outcomes", label + " — assignments by method", outcomes)
            self.add("12", ident + "_transitions", label + " — transitions", {str(i+1): r for i, r in enumerate(d.get("transitions", []))})
            self.add("12", ident + "_settings", label + " — recorded thresholds", {k: {"value": v} for k, v in d.get("settings", {}).items()})

    def completeness(self):
        self.add("11", "settings", "Effective QC settings", {b + " / " + k: {"value": v, "source": source}
                 for b, (cfg, source) in self.settings.items() for k, v in cfg.items()},
                 "Values are read from the production target report or run_summary.json. Defaults are never substituted for missing provenance. "
                 "These are recorded pipeline thresholds, not recommendations or validated assay specifications.")
        for batch, (cfg, source) in self.settings.items():
            if not cfg:
                self.warnings.append(batch + ": no recorded QC settings were found.")
        self.add("11", "coverage", "Per-well output availability", self.availability,
                 "Missing files produce blank metrics. Zero-read wells do not require trimming/QC outputs; noneligible wells do not require polishing. "
                 "Published files alone cannot determine whether absent outputs indicate an incomplete copy or an unfinished run.")
        self.add("11", "warnings", "Consistency checks", {str(i+1): {"finding": s} for i, s in enumerate(self.warnings)} or
                 {"Checks performed": {"finding": "No inconsistencies found among comparable available outputs; missing outputs are still listed separately."}},
                 "Checks include available read-count accounting, duplicate keys, summary/FASTA counts and the polished-diagnostic subset relationship. "
                 "This is a reporting check, not pipeline validation.")

    def write(self):
        out = Path(self.args.outdir).expanduser().resolve()
        out.mkdir(parents=True, exist_ok=False)
        (out / "tables").mkdir()
        if self.plot_assets:
            (out / "plots").mkdir()
            for name, svg in self.plot_assets.items():
                (out / "plots" / name).write_text(svg)
        (out / "empty_scan").mkdir()
        base = Path(self.args.base_config or Path(__file__).with_name("multiqc_base_config.yaml"))
        if not base.is_file():
            raise ValueError(f"Missing base config: {base}")
        # MultiQC otherwise sorts grouped subsections alphabetically.
        config = {"custom_content": {"order": sorted(self.order)},
                  "report_section_order": {sid: {"order": 10000 - i} for i, sid in enumerate(self.data)},
                  "custom_data": self.data}
        cfg = out / "multiqc_config.yaml"
        cfg.write_text(base.read_text().rstrip() + "\n\n# Generated measured content\n" + "\n".join(yaml_lines(config)) + "\n")
        (out / "report_data.json").write_text(json.dumps(config, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
        for sid, rows in self.tables.items():
            cols = list(dict.fromkeys(k for r in rows.values() for k in r))
            with (out / "tables" / (sid + ".tsv")).open("w", newline="") as f:
                writer = csv.DictWriter(f, fieldnames=["item"] + cols, delimiter="\t")
                writer.writeheader()
                for key, row in rows.items():
                    writer.writerow({"item": key, **row})
        manifest = {"collector_version": VERSION, "created_utc": datetime.now(timezone.utc).isoformat(),
                    "results_root": str(self.root), "batches": list(self.batches), "warnings": self.warnings,
                    "sources": self.sources, "settings_sources": {b: src for b, (_, src) in self.settings.items()},
                    "report_sections": list(self.data), "render_status": "NOT_REQUESTED"}
        mp = out / "collector_manifest.json"
        mp.write_text(json.dumps(manifest, indent=2) + "\n")
        cmd = [self.args.multiqc, str(out / "empty_scan"), "--config", str(cfg), "--outdir", str(out), "--filename", "tcr_multiqc.html", "--strict"]
        print(f"Generated {cfg}\nTables: {out / 'tables'}\nSections: {len(self.data)}; consistency findings: {len(self.warnings)}")
        if self.args.render:
            if not shutil.which(self.args.multiqc):
                raise ValueError("Configuration and tables were generated, but MultiQC is not installed/on PATH. Install MultiQC and rerun the render command in README.md.")
            manifest["multiqc_command"] = cmd
            with (out / "multiqc.log").open("w") as log:
                result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=False)
            manifest["render_status"] = "COMPLETE" if result.returncode == 0 and (out / "tcr_multiqc.html").is_file() else "FAILED"
            manifest["multiqc_exit_code"] = result.returncode
            mp.write_text(json.dumps(manifest, indent=2) + "\n")
            if manifest["render_status"] != "COMPLETE":
                raise ValueError(f"MultiQC failed; see {out / 'multiqc.log'}. Collected outputs are preserved.")
            print(f"Report: {out / 'tcr_multiqc.html'}")
        else:
            import shlex
            print("Render with:\n" + shlex.join(cmd))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--results-dir", required=True, help="Nextflow results root, or its tcr_pod5 directory")
    p.add_argument("--outdir", required=True, help="New report directory; existing directories are not overwritten")
    p.add_argument("--batch", action="append", help="Limit detailed production tables to this batch; repeatable")
    p.add_argument("--audit-dir", action="append", default=[], help="Completed comparison/audit directory; repeatable")
    p.add_argument("--assignment-evaluation", action="append", default=[], help="Path to assignment_evaluation.json; repeatable")
    p.add_argument("--no-discover-audits", action="store_true", help="Include only explicitly supplied exploratory audits")
    p.add_argument("--base-config", help="Alternative styling config; must not define custom_data or custom_content")
    p.add_argument("--render", action="store_true", help="Run MultiQC after collection")
    p.add_argument("--multiqc", default="multiqc", help="MultiQC command or executable path")
    p.add_argument("--version", action="version", version=VERSION)
    args = p.parse_args()
    try:
        if Path(args.outdir).expanduser().exists():
            raise ValueError("Output directory already exists; choose a new --outdir to avoid overwriting a report.")
        r = Report(args)
        r.load()
        r.overview()
        r.inputs_and_basecalls()
        r.demux_trim_assignment()
        r.target_and_qc()
        r.consensus()
        r.audits()
        r.completeness()
        r.write()
    except (ValueError, OSError, KeyError, TypeError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
