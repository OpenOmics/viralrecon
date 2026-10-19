# POD5-first TCR analysis extension for nf-core/viralrecon

**Status:** development prototype, not an assay-validated pipeline. The branch includes synthetic SAM-to-inference-to-QC tests plus an optional Dorado/samtools smoke test; it still requires full Nextflow/Dorado/GPU validation on the target system. See [validation and migration notes](tcr_pod5_target_inference.md).

This is the opt-in, POD5-only TCR path in the `OpenOmics/viralrecon` `tcr-pod5` branch, originally based on upstream commit `fa23078485cb75e96add952045b2b897aab61b42`. Stock workflows remain separate. No stock Nextflow regression run was performed for this update.

## 1. Implemented workflow

```
POD5 files grouped by supplied batch
    -> read-ID / metadata audit
    -> Dorado simplex basecalling, without trimming
    -> Dorado custom double-ended index classification and demultiplexing
       (both ends required; unclassified reads retained)
    -> per-well Dorado custom primer trimming and boundary validation
    -> per-well BAM + FASTQ
    -> competitive alignment to ALL full TCR references
    -> pool qualifying alignments across wells within each batch
    -> infer shared, reference-specific amplicon intervals
       (inferred_targets.bed is an OUTPUT, not an input)
    -> within-well heterogeneity / reconstruction screen
         hard mixture or structural blocker -> no polish; retain evidence and reads
         PASS_SCREEN or polish-eligible REVIEW
                         -> read-supported draft
                         -> Dorado polishing
                         -> realignment / coverage and variation recheck
                         -> CONSENSUS_PASS: accepted_consensus.fasta
                         -> REVIEW: diagnostic review_consensus.fasta
```

For PID201588, the expected constructs are **APC and PIK3CA**. The competitive panel can also include PID201403 and PID202478 as contamination sentinels. Use `expected_constructs` to declare which FASTA IDs are eligible for consensus; all records remain available for mapping. Each batch/well is a separate analysis unit, for example `preinfusion_201588_run1__A1`. Acquisitions are not equated to biological wells. Batch boundaries remain explicit even if acquisition IDs appear in more than one directory.

This branch bypasses ARTIC, Guppyplex, viral lineage analysis, and the old Cutadapt demultiplexer. It does not claim that a consensus or tree proves a well began with one physical molecule. Its interpretation is **no detectable mixture above the configured thresholds**.

## 2. Installation and software

Use the updated branch directly (do not reapply the old package installer):

```bash
git clone --branch tcr-pod5 https://github.com/OpenOmics/viralrecon.git
cd viralrecon
git rev-parse HEAD
```

For an existing clean checkout on `tcr-pod5`, use `git pull --ff-only origin tcr-pod5`. Pin and record the resulting commit for reproducibility.

Runtime requirements:

| Component | Implementation target |
|---|---|
| Nextflow | 25.04 or later, matching the pinned pipeline's requirements; not locally executed |
| Dorado | **2.1.2**, checked at runtime; official ONT installation provided separately |
| CPU dependencies | Python 3.11, samtools 1.21, POD5 Python package 0.3.28 |
| Basecalling | Compatible NVIDIA GPU allocation; default `cuda:0` |
| Polishing | CPU by default; compatible pre-populated Dorado model cache |

The `tcr_conda` profile supplies the CPU dependencies, **not Dorado or GPU drivers**. Alternatively, supply an approved environment on the execution nodes. This package does not supply tested Docker/Singularity images. Model caches must be accessible on those nodes. Dorado may try downloading a missing polishing model; pre-populate the cache on restricted/offline clusters and record the selected model from the log.

## 3. Inputs

### Batch manifest

Use `--tcr_input`, not the standard viralrecon FASTQ samplesheet. Required CSV columns are:

```
batch,pod5_dir,primer_table,sequencing_kit,basecall_model,reference_fasta
```

`assets/tcr_pod5/batches.example.csv` shows two illustrative batches with an optional `expected_constructs` column. Edit its placeholders before use. Paths are resolved relative to the CSV's directory unless absolute. The current implementation supports local/shared-filesystem paths, not remote object-store URLs. There must be exactly one row per batch; `.pod5` files are discovered recursively under that row's directory.

* `sequencing_kit` is the actual ONT library kit, **not** the custom TCR barcode scheme name. The example uses the previously reported `SQK-LSK114`; the audit checks it against nonblank POD5 metadata.
* `basecall_model` is a pre-downloaded, versioned simplex model directory containing `config.toml`, not the moving alias `sup`. Select a model compatible with the acquisition and Dorado release. Do not assume the old FASTQ's v4.1.0 model is supported by a new Dorado release.
* `reference_fasta` is the full competitive panel, with unique, safe FASTA IDs. All four supplied constructs can be included; no reference is removed just because its target cannot be inferred.
* Optional `expected_constructs` is a **semicolon-separated** list of FASTA IDs eligible for consensus, e.g. `PID201588_APC;PID201588_PIK3CA`. Other references are sentinels. A sentinel-dominated well is `REVIEW_UNEXPECTED_CONSTRUCT`, never an accepted consensus. If omitted or empty, all panel members are eligible; the pipeline does not guess biological roles from names.
* Remove the old `targets_bed` column. A nonempty legacy value is rejected with a migration message rather than silently ignored.

### Full primer table

Supply the complete TSV for each batch, including:

```
Shorthand  Direction  Sequence  I-start  I-end  Pair_with
```

The parser expects `fA1` paired with `rA1`, and so on. It interprets `I-start` / `I-end` as **1-based inclusive** index coordinates; `3` and `9` extract seven bases. It validates shared flanking sequences, equal index lengths, unique indexes, and complete pairing. The default expects 48 pairs. A complete table is required; unmatched or missing reverse primers fail validation.

From the table, the pipeline generates the custom Dorado TOML, barcode FASTA, well mapping, and a full-primer FASTA for each well. Rear-primer sequences are supplied in their original oligo orientation; Dorado handles reverse-complement recognition. `Frag-start` / `Frag-end` are not interpreted: their previous software-specific semantics are not established. This workflow intentionally trims the entire indexed PCR primer rather than reproducing an undocumented index-only trimming rule.

### Automatically inferred amplicon intervals

You no longer supply a BED. Full-panel alignment happens **before** target inference, without a full-vector target-coverage filter. High-quality, unique primary alignments with no split/SA/hard-clipped evidence are eligible. Inference additionally requires >=3,000 query-aligned bases and >=95% query coverage by default; missing MAPQ (255) is not treated as high confidence.

Each batch is handled separately. The algorithm finds a joint start/end cluster using per-well median endpoints as candidate centers, weights wells equally rather than by read count, and uses the median of the supporting per-well medians as its interval. The interval is fixed for all wells in that batch/reference. It does not infer a separate, easier target for each failing well.

Default acceptance requires 50 supporting reads, at least 2 wells with at least 5 supporting reads each, at least 5 reads on each alignment strand, and 80% well-balanced support within +/-75 bp of both endpoints. All thresholds live in `assets/tcr_pod5/qc_defaults.json` and can be overridden with `--tcr_qc_config`. These are exploratory assay-specific defaults.

Unsupported or unstable references receive explicit statuses and no BED interval. Assignment still counts reads supporting those references, so a low-support contaminant cannot disappear from mixture screening. A dominant construct lacking an inferred interval is held for review. There is **no full-vector fallback**.

The inferred BED/JSON/TSV and per-read boundary evidence are published under `<batch>/target_inference/targets/`. The JSON ties each interval to its batch, reference sequence digest, settings, support and strand counts. A batch with no demultiplexed reads still gets a target report for every panel member.

**Limits:** estimates are not independent proof of the complete PCR product. Systematic truncation, shared deletions, or a biased set of input wells can affect them. The default assumes one linear amplicon population per construct and is not designed for tiled or circular-origin-spanning amplicons. Inspect the inference report and validate with assay documentation/controls before interpreting production results. Historical PAF endpoints computed without base-level alignment are comparison material, not hard-coded boundaries.

## 4. Running

After editing the batch CSV, installing Dorado, and provisioning the model directories:

```bash
nextflow run /path/to/viralrecon-pod5 \
  --tcr_pod5 \
  --platform nanopore \
  --tcr_input /data/project/config/batches.csv \
  --tcr_dorado_bin /opt/dorado-2.1.2-linux-x64/bin/dorado \
  --tcr_polish_models /data/models/dorado-polish \
  --outdir /data/project/results \
  -profile tcr_conda \
  -c /data/project/config/site-slurm.config \
  -resume
```

Replace the paths and site configuration. Provide a site-specific Slurm configuration; GPU allocation remains an explicit site responsibility. The `tcr_gpu` label needs a real GPU resource request. Without a scheduler configuration, Nextflow runs locally.

TCR mode rejects `--input` and `--fastq_dir`. FASTQs are outputs only. Keep the existing FASTQ-based patch in a separate checkout for comparison, not layered beneath this one.

Useful options:

| Parameter | Default / purpose |
|---|---|
| `--tcr_expected_pairs` | `48`; lower only for intentionally reduced schemes/tests |
| `--tcr_barcode_errors` | `1`; barcode edit penalty, not whole-primer matching tolerance |
| `--tcr_qc_config` | Optional JSON overriding settings in `assets/tcr_pod5/qc_defaults.json` |
| `--tcr_allow_empty_pod5` | `false`; explicit opt-in to ignore audited empty files |
| `--tcr_save_basecalls` | `true`; use `false` to omit publishing the large basecalled BAM |
| `--tcr_basecall_device` | `cuda:0` |
| `--tcr_polish_device` | `cpu` |
| `--tcr_samtools_bin` | `samtools`; an absolute path is also accepted |

Zero-byte or zero-read POD5s fail by default. Repeated raw read UUIDs within or across batches always fail rather than silently deduplicating. The audit reads POD5 containers and metadata, not every signal sample; full signal decoding occurs during basecalling. On audit failure, inspect the task work directory and logs, since failed-task outputs may not be published.

## 5. Demultiplexing, trimming, and metadata

Basecalling uses `--no-trim --emit-moves --min-qscore 0 --disable-read-splitting`. Keeping barcodes intact allows subsequent classification. Disabling read splitting makes the initial raw-read accounting explicit; chimeric/split alignment evidence is inspected downstream rather than silently creating extra basecalled records.

A separate `dorado demux` call uses the generated arrangement, `--kit-name TCR_WELLS`, `--barcode-both-ends`, `--no-trim`, and `--emit-summary`. It performs classification and splitting together. Assignment uses the BAM `BC` tag, not chunk filenames. Read counts and uniqueness are checked, and unclassified reads are retained.

Dorado then trims each assigned well using its full forward/reverse primers. An independent boundary check requires both full primers in the original read-end windows and verifies that Dorado removed them without excessive extra trimming. Failures are recorded, not silently rescued or deleted.

**BAM is authoritative** throughout. Per-well FASTQs are also written, but read-group/model information and move tables remain available in BAM for polishing. Reads are not manually reverse-complemented into canonical FASTQ orientation; alignment handles both strands. This avoids rewriting sequences while leaving move metadata inconsistent.

The seven-base barcode scoring settings are deliberately conservative starting values. Dorado's standard kit defaults are not an assay validation. The generated TOML records the exact settings, including a strict separation requirement and disabled high-penalty fallback for this short-index scheme.

## 6. Heterogeneity and consensus gate

Competitive Dorado/minimap2 alignment against the full reference panel retains secondary hits so shared vector or constant-region sequence is not automatically treated as an unambiguous construct assignment. Ambiguous, partial-query, low-quality, and split alignments are accounted for separately. A target-coverage filter is NOT used to assign references. Once a target has been inferred, target coverage is assessed for draft eligibility without erasing assignment counts. Excess target-incomplete reads or endpoints departing from the batch interval cause review.

The custom Python screen examines both:

1. **Between-construct mixtures:** substantial, strand-supported populations assigned to different panel members, including unexpected sentinel constructs.
2. **Within-construct variation:** candidate SNPs and primitive indels, including linked alleles on the same reads. Two coherent, strand-supported patterns at separated sites trigger rejection. Isolated or weaker variation triggers review rather than a claim of two templates.

Examples of the **exploratory defaults**, not established sensitivity/specificity:

| Check | Default |
|---|---|
| Minimum support for final accepted consensus | 50 reads |
| Minimum support / callable fraction to attempt polishing | 20 reads / 80% |
| Minimum read / base quality | Q10 / Q15 |
| Competitive alignment identity / query coverage | 90% / 85%; no full-vector coverage requirement |
| Internal deletion exemption from base-error identity | D >=50 bp, with >=20 contiguous aligned bases on each side |
| Diagnostic draft support / per-position depth | 10 dominant reads / 10; at least 2 reads per strand |
| Diagnostic draft callable fraction | >=10%; unsupported positions remain N |
| Dominant-read target coverage for draft eligibility | 85% of the batch-inferred interval |
| Mapping quality; alignment-score separation | 20; 20 |
| Strong minor population | At least 5 reads and 10%, with at least 2 per strand |
| Weaker variation alert | At least 3 reads and 3% |
| Minimum consensus depth / base support | 20 / 80% |
| Required callable amplicon fraction | 98% |

The minority fraction is calculated among relevant assigned or informative reads, not raw POD5 counts. Candidate sites remain sensitive at >=3 reads and >=3%, but a site blocks consensus only when the minor allele has >=5 reads, >=10% support, >=2 reads on each strand, and is not homopolymer-associated. Strong homopolymer-associated, weaker, or strand-biased sites remain explicit warnings/evidence. The same-construct linked test still requires two separated informative sites and remains a hard rejection criterion. Fixed differences shared by all reads are not themselves evidence of a mixture.

`ambiguous_reference` reads remain excluded from construct evidence and consensus, but are reported separately rather than counted in the unresolved-fraction numerator because shared vector sequence can make otherwise good reads non-discriminating. A high primer-trim failure fraction is likewise retained as a warning; only successfully trimmed reads can enter downstream evidence.

Polishing and acceptance are separate decisions. Expected constructs with an inferred interval may enter Dorado polish when at least 20 complete dominant reads remain, the draft is at least 80% callable, strand support is adequate, and there is no hard mixture/structural blocker. Final `CONSENSUS_PASS` still requires the configured 50-read support, depth/callable thresholds, coherent post-polish alignments, no fixed-indel disagreement, and no strong bidirectional non-homopolymer variation. A polished sequence that remains `REVIEW` is retained in `review_consensus.fasta` for inspection but is never added to `accepted_consensus.fasta`.

### Deletions and assignment evaluation (v0.4.0)

Long internal deletions are evaluated separately from base-error identity. For an alignment spanning `S` columns with `NM` edit bases and `D` qualifying deleted bases, the assignment identity is `(S - NM) / (S - D)`. Only CIGAR D operations at least `deletion_min_length=50` bp long with at least `deletion_min_anchor=20` contiguous aligned bases on **each** side qualify. Short or terminal deletions, insertions and substitutions remain penalties. These length/anchor defaults are exploratory. Read quality, query coverage, MAPQ, competitive score separation, and split-alignment checks still apply; a missing NM tag cannot imply a perfect match.

This exemption retains candidates for evaluation; it does not certify the deletion. Deletion events still contribute to the well's heterogeneity and linked-haplotype evidence. A consistent, supported deletion is represented in the draft; intact/deleted subpopulations remain visible. A long deletion mixture is not downgraded to a homopolymer warning merely because its breakpoint is near a homopolymer. Equivalent repeat deletions are projected onto the normalized event coordinates so normalization does not introduce spurious masked bases. `qc/read_assignments.tsv` reports both identity calculations and the exempted deletion coordinates/lengths. Batch target inference still requires 3,000 query-aligned bases by default: a strongly deleted read can contribute to QC using an already inferred target without necessarily contributing endpoints to target inference.

`evaluate-assignment` provides an **offline comparison**, without changing production assignments or acceptance rules. It compares the previous NM-based method, deletion-aware production assignment, and informative-position assignment. Use the full panel, including contamination sentinels. Supply a reference-only multiple alignment with exactly the original IDs and ungapped sequences; a mismatch fails validation.

For example, align the reference panel once using MAFFT (an optional audit dependency, not required for the production workflow), and name-sort the existing per-well full-panel BAMs:

```bash
mafft --auto --thread 4 reference.fasta > reference_panel.aligned.fasta
samtools sort -n -o A1.name_sorted.bam /results/tcr_pod5/Run1/wells/A1/alignment/aligned.bam
```

Create a JSON manifest listing **all wells being compared**, with paths accessible from the working directory:

```json
{
  "batch": "Run1",
  "alignments": [
    {"well": "A1", "bam": "/absolute/path/A1.name_sorted.bam"},
    {"well": "A2", "bam": "/absolute/path/A2.name_sorted.bam"}
  ]
}
```

```bash
python3 bin/tcr_pod5.py evaluate-assignment \
  --manifest assignment_inputs.json \
  --references reference.fasta \
  --reference-alignment reference_panel.aligned.fasta \
  --settings assets/tcr_pod5/qc_defaults.json \
  --outdir assignment_evaluation
```

The audit writes `read_assignment_comparison.tsv` and `assignment_evaluation.json`, including per-well counts, transitions, reference/alignment digests, and exact informative columns. Default informative criteria are identity >=0.75, a best-versus-second-best margin >=0.15, and at least 50 quality-filtered informative positions. Every candidate is scored over the **same observed positions**. Only variable columns containing A/C/G/T in every reference are used in this first implementation; gapped or N-containing variable columns are listed separately. Deleted read positions are missing nucleotide evidence, not mismatches or support for an alternative construct. Thus this audit need not reproduce the collaborator's count of 798 positions.

Read bases are projected through the primary BAM alignment into the reference alignment; original reads and insertions are preserved. The alternative scoring omits MAPQ/AS and whole-read identity gates so it can examine reads rejected as ambiguous by the production method. Quality, query-coverage and split-alignment exclusions still apply. Primary placement in repeats can bias these scores, and a deleted insert may remove the ability to distinguish constructs. Inspect changed reads, sentinel assignments and known controls before adopting this method. Neither extra assignments nor fewer unresolved reads demonstrate improved accuracy. This is an evaluation of the informative-position approach, not an exact implementation of MAFFT read placement or an automatic rescue.

### Diagnostic drafts (v0.4.0)

Wells held for review can now produce a separate **unpolished, masked diagnostic draft**, even when they cannot enter Dorado polishing. Defaults require 10 complete, assigned dominant-construct reads, at least 2 on each strand, an inferred expected-construct target, and at least 10% callable target positions. At each called position, depth must be at least `diagnostic_min_depth=10`, base quality at least Q15, and support at least 80%. Strong candidate variation is additionally masked, including strong homopolymer-associated variation. No reference bases fill coverage gaps. The draft includes supported insertions/deletions; conflicting fixed indels block its output.

These support levels are for inspecting the available evidence and are **not validated consensus accuracy thresholds**. Unresolved-read, endpoint, partial-read, split-read or weak second-construct flags remain in the decision; the draft represents only the selected dominant subset and does not resolve those concerns. Rejected mixtures, unexpected constructs, missing targets, inadequate support and conflicting indels do not receive a diagnostic draft. A zero-assignment well now has a null dominant construct instead of an arbitrary panel member.

`qc/decision.json` records diagnostic eligibility, blockers, support, masking thresholds and original reasons. Its per-well draft is a pre-polish artifact. Final aggregation includes only wells whose final status is `REVIEW`, prefers the polished diagnostic when available, and otherwise uses the unpolished draft. It excludes later post-polish rejections even if a pre-polish draft exists. FASTA descriptions explicitly say `DIAGNOSTIC_UNPOLISHED` or `DIAGNOSTIC_POLISHED` and `accepted=false`. The TSV also records the kind, availability and callable fraction; this fraction can differ from the stricter pre-polish `callable_fraction` column. `usable_reads` remains the count assigned across all constructs; `selected_reads` identifies the complete dominant subset used for reconstruction.

The accepted thresholds remain 50 reads, at least 5 per strand, depth 20, 80% base support and 98% callability. Strand support is rechecked after polishing, and post-polish linked haplotypes are a hard rejection. Strong candidate variation is masked in polished output too. Pre-polish biological/structural review reasons persist: losing the minor population during remapping is not evidence that a mixture was resolved. Only numerical support/coverage checks are reassessed. Diagnostic availability never changes a REVIEW well to `CONSENSUS_PASS`.

**Limits:** identical templates are indistinguishable; low-abundance, closely related, or poorly sequenced mixtures can be missed; systematic ONT/PCR errors can cause review or false mixture signals. Complex rearrangements and repeat-associated indels require manual examination. This is a conservative candidate screen, not a validated haplotype caller, a phylogenetic proof of single-template origin, or a clinical pipeline. Do not use it to estimate the number of original molecules from sequencing read counts.

## 7. Outputs

Outputs are published beneath `outdir/tcr_pod5/`:

```
audit/                              POD5 inventory and raw UUID provenance
<batch>/scheme/scheme/              generated scheme and full reference panel
<batch>/basecalling/                calls.bam, summary, model/version provenance
<batch>/demultiplexing/demux/        assignment TSV, per-well BAMs, unclassified reads
<batch>/wells/<well>/trimming/       clean.bam, <well>.fastq.gz, failed reads, audit
<batch>/wells/<well>/alignment/      full-panel coordinate-sorted BAM and index
<batch>/target_inference/targets/    inferred BED, JSON, summary TSV, read evidence
<batch>/wells/<well>/heterogeneity/qc/
                                    decision, assignment, allele/haplotype evidence,
                                    optional masked diagnostic draft and mask BED
<batch>/wells/<well>/consensus/consensus_result/
                                    post-polish decision and accepted/review sequence
summary/                            well_summary.tsv, run_summary.json,
                                    accepted_consensus.fasta, review_consensus.fasta,
                                    diagnostic_consensus.fasta
```

Every expected well appears in the final summary, including `NO_READS`. Mixed and review wells are excluded from `accepted_consensus.fasta` without deleting their evidence. `diagnostic_consensus.fasta` contains one sequence per eligible final REVIEW well. `review_consensus.fasta` is its polished subset, retained for compatibility; do not concatenate the two files. An all-rejected run legitimately produces empty accepted and diagnostic FASTAs. IDs include batch and well to prevent collisions. Variant positions in heterogeneity QC refer to the **original full reference**, restricted to the inferred target. Draft sequences span that target. Each masked-position BED refers to its associated **output sequence after indels**, whose coordinates can differ from the reference. `inferred_targets.json` records original-reference boundaries.

This branch produces TSV/JSON reports, not a customized MultiQC report or a final multiple-sequence alignment/tree. Those are separate extensions beyond the requested consensus gate.

## 8. Validation and interpretation

Run the repository tests:

```bash
python3 -m unittest discover -s tests/tcr_pod5 -v
```

With Dorado 2.1.2 and samtools available, enable the optional executable smoke test:

```bash
TCR_TEST_DORADO=/opt/dorado-2.1.2-linux-x64/bin/dorado \
TCR_TEST_SAMTOOLS=/path/to/samtools \
python3 -m unittest discover -s tests/tcr_pod5 -v
```

That test exercises **synthetic barcode classification and trimming only**, not POD5 basecalling or polishing. A successful result does not validate the full workflow or biological thresholds.

The reconstruction regressions additionally cover large deletion retention, ordinary error penalties, repeat normalization, informative-site coverage/margins, low-depth masking, diagnostic/accepted output separation, and post-polish strand/linked-haplotype gates. Post-polish driver tests mock external commands; they do not exercise the Dorado model. The Run1 summary files alone cannot establish how many reads or wells these changes recover: rerun QC on the BAMs and compare the assignment audit before reporting a revised yield.

Before broader use, perform a small POD5 end-to-end run on the cluster and challenge the QC with known clean controls and deliberate mixtures. The historical Run1 FASTA and well counts are comparison material, not established ground truth. A new basecaller and stricter index/primer matching can change counts. Historical renamed FASTA IDs do not automatically provide a UUID crosswalk for read-by-read comparison.

## References consulted

- [nf-core/viralrecon pinned source](https://github.com/nf-core/viralrecon/tree/fa23078485cb75e96add952045b2b897aab61b42)
- [Dorado 2.1.2 release](https://github.com/nanoporetech/dorado/releases/tag/v2.1.2)
- [Dorado custom barcode arrangements](https://software-docs.nanoporetech.com/dorado/latest/barcoding/custom_barcodes/)
- [Dorado barcode classification and demultiplexing](https://software-docs.nanoporetech.com/dorado/latest/barcoding/barcoding/)
- [Dorado custom primers](https://software-docs.nanoporetech.com/dorado/latest/barcoding/custom_primers/)
- [Dorado read trimming](https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_trimming/)
- [Dorado alignment](https://software-docs.nanoporetech.com/dorado/latest/basecaller/alignment/)
- [Dorado polishing requirements](https://software-docs.nanoporetech.com/dorado/latest/secondary/polish/)
- [Minimap2 base-level alignment and PAF semantics](https://github.com/lh3/minimap2)
- [POD5 tools and Python API](https://pod5-file-format.readthedocs.io/en/latest/)

The mixture rules and thresholds are custom prototype choices, not recommendations established by those tool manuals.
