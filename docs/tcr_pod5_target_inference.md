# BED-free target inference: migration and validation

## Changes

The `tcr-pod5` entry point no longer accepts a nonempty `targets_bed` manifest value. Mapping now uses the full reference panel. One additional `TCR_INFER_TARGETS` task pools primer-validated, competitively assigned reads across each batch's wells, estimates per-reference boundaries, and supplies those fixed intervals to well QC. The stock workflow entry point is not modified by this update.

The six required manifest fields are:

```csv
batch,pod5_dir,primer_table,sequencing_kit,basecall_model,reference_fasta
```

For a panel containing possible contaminants, also use the optional `expected_constructs` field. Separate FASTA IDs with semicolons, not commas. For the PID201588 experiments this is `PID201588_APC;PID201588_PIK3CA`, assuming those are the headers in the full panel. All other panel members remain alignment/mixture candidates but cannot yield an accepted consensus. Without an expected list, all panel members are consensus-eligible.

Construct names now remain the original FASTA IDs throughout. Previously the fourth BED column could rename a construct (for example `PID201588_APC` to `APC`); those aliases are no longer applied. Update downstream report consumers accordingly. The provided schema table is `assets/tcr_pod5/batches.example.csv`.

## Outputs and statuses

Each batch publishes:

```text
tcr_pod5/<batch>/target_inference/targets/
  inferred_targets.bed       # only INFERRED rows; original-reference BED4
  inferred_targets.json      # every reference, digest, method, settings and support
  target_summary.tsv         # every reference, including absent sentinels
  boundary_evidence.tsv      # qualifying read IDs, wells, endpoints and membership
```

Target statuses are `INFERRED`, `INSUFFICIENT_SUPPORT`, and `UNSTABLE_BOUNDARIES`. Failing inference is not an error that fabricates a full-length interval: affected wells are `REVIEW_TARGET_NOT_INFERRED`. Mixtures can still be rejected based on full-panel assignments even when a minor reference has no inferred target. Sentinel-dominated wells are `REVIEW_UNEXPECTED_CONSTRUCT` when an expected list is supplied.

Coverage/draft evidence uses the frozen batch interval. Original reference coordinates are retained in heterogeneity reports, while final masked-position coordinates refer to the polished draft. Target-partial reads remain in reference-assignment counts. Excess partial reads or displaced endpoints trigger review instead of shrinking the per-well target.

## Initial inference settings

| JSON key | Default |
|---|---:|
| `target_min_reads` | 50 |
| `target_min_wells` | 2 |
| `target_min_reads_per_well` | 5 |
| `target_min_each_strand` | 5 |
| `target_min_alignment_length` | 3000 |
| `target_min_query_coverage` | 0.95 |
| `target_boundary_tolerance` | 75 |
| `target_min_cluster_fraction` | 0.80 |

These supplement the existing identity, quality and alignment-score separation filters. Full-reference coverage is deliberately NOT a candidate-assignment filter. MAPQ 255, supplementary/SA-tagged alignments and hard-clipped alignments are excluded. Each qualifying well has equal weight in endpoint-cluster selection, so a single very deep well cannot overwhelm multiple other wells. Supporting-well medians are combined into the shared interval, rounded outward to integer boundaries.

A one-well pilot can audit/basecall/demultiplex and map, but cannot pass the default two-well inference gate. Prefer a pilot containing multiple wells; do not weaken thresholds merely to obtain a consensus. Uniform truncation across all wells can evade an empirical-boundary screen, and identical templates remain indistinguishable. Controls and assay knowledge remain necessary.

## Validation performed for this update

```bash
python3 -m unittest discover -s tests/tcr_pod5 -v
python3 -m compileall -q bin/tcr_pod5*.py tests/tcr_pod5
git diff --check
```

Result in the development environment: **92 tests run, 91 passed, 1 skipped**. The skipped test requires real Dorado and samtools for synthetic barcode demultiplexing/primer trimming. POD5 audit tests use mocked POD5 readers. Real signal decoding, GPU basecalling, polishing, Nextflow channel execution and the cluster scheduler were **not** run here.

The suite includes:

- Empty batches and absent references, duplicate UUID rejection and provenance mismatches.
- Multiple-well interval recovery, depth imbalance, outliers, bimodal endpoints and strand requirements.
- Full-panel assignment before target inference, missing MAPQ, SA/supplementary/hard-clipped evidence.
- Inferred-interval coverage, target-only drafts, original-reference variant coordinates, shifted/partial wells.
- Low-support sentinels, pure unexpected constructs, mixed-construct rejection and consensus gating.
- CLI-level preparation without BED and synthetic SAM -> inference -> well QC/reporting.
- Existing primer, trimming, mixture, indel and read-group regression tests brought into the repository.

Next required validation is a small multi-well POD5 run on the target HPC/GPU environment. Inspect `target_summary.tsv`, the boundary evidence and alignment plots before accepting intervals or running the full dataset. These tests establish software behavior, not assay sensitivity or specificity.
