# Optional TCR consensus phylogeny

Enable with `--tcr_phylogeny`. The default is **off**. The optional
`TCR_POD5:TCR_PHYLOGENY` task runs once per batch, after final consensus reporting.
It compares sequences with MAFFT and RAxML-NG and never changes well decisions,
consensus bases, or consensus masks. `--skip_multiqc` suppresses the report but
does not suppress a requested phylogeny.

## Run or resume

Add these arguments to the same complete Nextflow command used for the run:

```bash
--tcr_phylogeny --tcr_phylogeny_bootstraps 200 -resume
```

Keep the original `-work-dir`, Nextflow launch/cache directory, input paths and
execution settings for cache reuse. The new helper is a separate task input;
this feature does not modify upstream QC scripts. MultiQC reruns to include the
new results. Published results alone are not a Nextflow cache.

The `tcr_conda` profile provides a separate environment with Python 3.11,
MAFFT 7.526, RAxML-NG 1.2.2, Biopython 1.85 and Matplotlib 3.10.3. Without Conda,
install these on the execution node. External installations can be selected with
`--tcr_mafft_bin /path/to/mafft` and `--tcr_raxml_ng_bin /path/to/raxml-ng`.
This feature does not supply a container image. Executables and dependencies
must be visible inside any container used by the site's execution profile.

| Parameter | Default | Meaning |
|---|---|---|
| `tcr_phylogeny` | `false` | Enable the optional per-batch analysis. |
| `tcr_phylogeny_model` | `GTR+G` | Nucleotide substitution model; no automatic model selection. |
| `tcr_phylogeny_bootstraps` | `200` | Fixed number of Felsenstein bootstrap replicates. |
| `tcr_phylogeny_seed` | `17` | RAxML-NG random seed. |
| `tcr_phylogeny_min_target_overlap` | `0.85` | Minimum fraction of each donor's inferred span retained in the common crop. |

Resources default to 4 CPUs, 8 GB and 12 hours. RAxML-NG uses `--threads auto{N}`
to choose an appropriate thread count up to the allocation. Override resources
through `withName: 'TCR_POD5:TCR_PHYLOGENY'` in a custom Nextflow config. A fixed
seed does not guarantee identical output across software versions or settings;
versions, commands and input hashes are recorded.

## Sequences and reference cropping

Each batch includes:

* Every `CONSENSUS_PASS` record from `accepted_consensus.fasta`, labelled `PASS`.
* Every available final `REVIEW` record from `diagnostic_consensus.fasta`, labelled
  `REVIEW_POLISHED` or `REVIEW_UNPOLISHED`. The pipeline has already chosen one
  sequence per well; `review_consensus.fasta` is a subset and is not appended again.
* Every reference in that batch's complete competitive panel, including
  contamination sentinel references, labelled `REFERENCE`.

REJECT and NO_READS wells do not contribute consensus sequences. REVIEW wells
without an available diagnostic sequence are listed in `excluded_wells.tsv`.
Missing FASTA records promised by the final summary are errors, not exclusions.
No new minimum length, coverage or quality filter is applied to available
consensuses. Original REVIEW reasons remain in the sequence metadata.

The full reference panel is first aligned **without the consensus sequences**.
Consensus-eligible references with supported `INFERRED` boundaries define a
common overlapping span in reference-alignment coordinates. That span is
projected back into **every** original reference. Each donor must retain at
least the configured fraction of its own inferred span. This handles coordinate
shifts caused by insertions/deletions; it does not reuse APC coordinates directly
on other constructs. Crop coordinates, donor references, original inference
status and the reference-only alignment are saved. The target report must match
the batch and reference digest.

These projections are homology-based comparison intervals, not additional
experimental evidence of primer boundaries. There is no fallback to untrimmed
full-vector references when inference fails. Incompatible/nonoverlapping target
spans, unsupported boundaries or an empty projected reference produce a clear
`SKIPPED` report.

## Alignment and inference

MAFFT uses `--auto --thread N --threadit 0 --inputorder` for both reference and
final sequence alignment. Iterative-refinement threading is disabled, as
recommended in the [MAFFT threading documentation](https://mafft.cbrc.jp/alignment/software/multithreading.html).
The code verifies every ID and ungapped sequence against its input; it never
reorients a sequence, replaces Ns with reference bases, or drops identical tips.

The original final alignment is preserved. Only columns consisting entirely of
N/gaps are removed from the RAxML-NG input; **invariant columns are retained**.
All included tips remain present. Gaps and fully ambiguous bases are missing
information for the nucleotide model, not a separate indel state or evidence
that a consensus matches the reference. Per-tip gap/ambiguity percentages,
unambiguous variable/informative site counts and identical tip groups are saved.

RAxML-NG runs `--all --model GTR+G --bs-trees 200 --bs-metric fbp` by default,
with the recorded seed and thread cap. See the
[RAxML-NG tutorial](https://github.com/amkozlov/raxml-ng/wiki/Tutorial).
The optional reduced alignment produced by RAxML-NG is **not used**, because it
can remove duplicate sequences. Output trees are checked against all input IDs,
including bootstrap replicates. The requested bootstrap count must be present.

The alignment is retained but no tree is estimated if there are fewer than four
distinguishable sequences, no parsimony-informative nucleotide sites, or a
requested tip has no unambiguous aligned bases. No requested tip is silently
dropped to make inference succeed. Executable failures, invalid provenance,
missing expected files and lost/renamed tips fail the task visibly.

The tree is **exploratory and unrooted**. Its rectangular drawing uses an arbitrary
display orientation; none of the references is assumed to be an outgroup.
Engineered APC/PIK3CA differences and shared vector sequence can dominate the
result. Masked sequences, nearly identical sequences, homopolymer artifacts and
limited within-construct signal can make some placements uncertain. A tree or
high bootstrap value does not establish one physical template, convert REVIEW
to PASS, or by itself establish natural evolutionary ancestry. All available
REVIEW consensuses are included for inspection, even when heavily masked.

## Outputs and MultiQC

Outputs are under `<outdir>/tcr_pod5/<batch>/phylogeny/`:

| File | Contents |
|---|---|
| `phylogeny.json` | COMPLETE/SKIPPED state, reason, counts, settings, versions, input hashes and exact commands. |
| `sequence_metadata.tsv` | Unique tip IDs, original sample/well, status, polishing kind, construct, sequence composition and original QC reasons. |
| `excluded_wells.tsv` | Wells without an included consensus and why. |
| `reference_crops.tsv` | Crop intervals (0-based, end-exclusive), projection provenance and donor overlap. |
| `references_aligned.fasta` | Full reference-only MAFFT alignment used to project boundaries. |
| `cropped_references.fasta` | Cropped, ungapped reference sequences. |
| `sequences.fasta`, `sequences_aligned.fasta` | All selected sequences before and after final MAFFT alignment. |
| `tree_input.fasta`, `alignment_metrics.tsv` | Alignment used for inference and per-tip N/ambiguity/gap metrics. |
| `best_tree.nwk`, `supported_tree.nwk` | Maximum-likelihood tree and tree labelled with bootstrap support. |
| `bootstrap_trees.nwk`, `fitted_model.txt` | Bootstrap replicates and fitted model. |
| `tree.svg` | Tree with coloured status labels and bootstrap values. |
| `raxml/`, `*.log` | Native RAxML-NG outputs and execution logs. |

Skipped analyses contain available tables and the explicit skip reason, but no
invented tree. Early crop skips do not have a final alignment.

MultiQC adds a summary in section 00, **section 17** for selection/crops/alignment
(with counts and per-tip sequence composition plots), and **section 18** for
tree settings, bootstrap summaries, versions and the embedded tree. Alignment
and Newick downloads are embedded in the report. Sections 12–16 remain reserved
for the earlier offline audits. Optional phylogeny sections appear only when
their outputs were supplied. They also appear in a standalone report rebuilt
from a completed results directory containing `phylogeny/` outputs.

For completed results without a Nextflow cache, the helper can run directly:

```bash
python3 bin/tcr_pod5_phylogeny.py \
  --batch YOUR_BATCH \
  --references RESULTS/tcr_pod5/YOUR_BATCH/scheme/scheme/reference.fasta \
  --targets RESULTS/tcr_pod5/YOUR_BATCH/target_inference/targets/inferred_targets.json \
  --summary RESULTS/tcr_pod5/summary/well_summary.tsv \
  --accepted RESULTS/tcr_pod5/summary/accepted_consensus.fasta \
  --review RESULTS/tcr_pod5/summary/diagnostic_consensus.fasta \
  --outdir RESULTS/tcr_pod5/YOUR_BATCH/phylogeny --threads 4
```

Use actual paths if your results root is already `tcr_pod5`; the output directory
must be new. Then rebuild MultiQC using the standalone collector documented in
[tcr_multiqc.md](tcr_multiqc.md).

## Validation

```bash
python3 -m unittest discover -s tests/tcr_phylogeny -v
python3 -m unittest discover -s tests/tcr_multiqc -v
```

The tests cover final-status selection, batch separation, reference crop
projection with insertions, provenance mismatches, missing/duplicate inputs,
mask preservation, identical tips, invariant columns, skips, command failures,
report embedding and display aliases. A real MAFFT/RAxML-NG smoke test and SVG
rendering test run when the optional environment is installed; otherwise they
are explicitly skipped. Mocked command tests check file interfaces only and do
not validate scientific inference. The integration also requires a resumed
Nextflow run on the target cluster; the authoring environment did not provide
Nextflow, MAFFT, RAxML-NG, Biopython or MultiQC, or allow solving the Conda environment.
