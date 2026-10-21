# Production TCR MultiQC report

The final `TCR_MULTIQC` task creates one combined report for all batches after
`TCR_REPORT` and its upstream inputs finish. It runs by default in TCR POD5 mode;
the existing `--skip_multiqc` parameter disables it. Stock viralrecon workflows
are unchanged.

## Production sections

| Sections | Included by default | Purpose |
|---|---|---|
| 00–09 | Yes | Overview, input audit, preparation, basecalling, demultiplexing, trimming, assignment, target inference, well QC, and polishing. |
| 10 | Yes | Final accepted and diagnostic consensus outputs, final decisions, FASTA counts, sequence lengths, and N masking. Polished diagnostics are a subset of all diagnostics. |
| 11 | Yes | Recorded QC settings, per-well availability, and consistency findings. |
| 12–16 | No | Exploratory assignment, barcode-margin, read-end, splitting, and internal-junction comparisons. These remain available through the standalone collector. |

The module passes `--no-discover-audits` and supplies no audit paths. It never
runs those exploratory analyses. Production assignment/deletion statistics
remain in section 06. Diagnostic sequences in section 10 remain explicitly
separate from accepted consensus sequences.

## Inputs, execution, and outputs

The task receives selected JSON, TSV, TOML, and small FASTA files directly from
Nextflow output channels. It does not scan the published output directory.
Unique staging directories preserve repeated filenames across batches and wells.
The mapping is sorted before task creation for deterministic input order.
Empty polishing channels and zero accepted consensus sequences are supported.
It does not stage BAM/POD5 files, rerun analyses, or alter QC thresholds.

The module reconstructs a private input layout, runs the collector and MultiQC,
and publishes:

```text
<outdir>/tcr_pod5/multiqc/tcr_multiqc.html
```

The same directory contains `multiqc_config.yaml`, exported tables and SVGs,
`report_data.json`, the collector manifest, the source-to-destination input
manifest, reporting software versions, and `multiqc.log`. Interactive quality
and read-length histograms retain SVG fallbacks. Their controls are local to the
panels and do not participate in MultiQC's global sample toolbox/bulk export.

MultiQC runs with `--strict`. A failed report task fails visibly; it can be
resumed after fixing the reporting environment. Consistency findings are shown
in section 11; they do not reclassify wells. The report is intended for completed
workflow outputs, not as a partial-run monitor.

## Software

The `tcr_conda` profile uses a separate reporting environment pinned to Python
3.11, MultiQC 1.30, and Matplotlib 3.10.3. Existing Dorado/QC environments are
unchanged. With Conda disabled, Python, MultiQC, and Matplotlib must be available
on the node running `TCR_MULTIQC`. No custom MultiQC plugin or Python Plotly
package is required. This integration does not provide a new container image.

The TCR report uses `assets/tcr_pod5/multiqc/multiqc_base_config.yaml`. Generic
stock-workflow styling parameters such as `--multiqc_config` and
`--multiqc_title` are not wired into this custom TCR report; edit its base config
to change its appearance. That directory is a declared task input, so edits to
the collector, histogram controller, or styling invalidate the report task cache.

After applying the integration, rerun the original Nextflow command with
`-resume`. Retain the original work directory and cache. Existing computational
process definitions, scripts, and their environments are unchanged; Nextflow
should reuse valid cached tasks and run the new report task. Cache reuse still
depends on unchanged inputs and execution settings. Simply placing old published
results under `--outdir` does not provide a Nextflow cache.

For completed results without their work/cache, use the bundled standalone
collector:

```bash
python3 assets/tcr_pod5/multiqc/collect_tcr_multiqc.py \
  --results-dir /path/to/completed/results \
  --outdir /path/to/new/production_report \
  --no-discover-audits --render
```

For exploratory reporting, omit `--no-discover-audits` or explicitly pass
`--audit-dir` / `--assignment-evaluation` to the standalone collector. Such
reports do not change the production output totals.

## Validation

Run from the repository root:

```bash
python3 -m unittest discover -s tests/tcr_multiqc -v
node tests/tcr_multiqc/test_histogram_interactivity.cjs
```

Python checks cover existing accounting, missing-versus-zero outputs, multiple
batches, post-polish distinctions, numeric histogram bins, duplicate filenames,
staging failures, and a zero-assigned-read run without polishing. The JavaScript
test checks controller behavior using a mocked DOM/Plotly interface. These are
not Nextflow scheduling or browser rendering tests. The integration still needs
a resumed run on the target system: Nextflow, MultiQC, and a browser runtime were
unavailable in the authoring environment, and the pinned Conda environment was
not solved there.
