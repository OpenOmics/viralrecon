# POD5-first TCR analysis extension for nf-core/viralrecon

**Status:** development prototype, not an assay-validated pipeline. The BED-free update has **91 passing Python tests and 1 skipped external-tool test**. These include synthetic SAM-to-inference-to-QC integration tests, not a Nextflow/Dorado/GPU end-to-end run. See [validation and migration notes](tcr_pod5_target_inference.md).

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
    -> within-well heterogeneity screen
         REJECT / REVIEW -> no consensus; retain evidence and reads
         PASS_SCREEN     -> read-supported draft
                             -> Dorado polishing
                             -> realignment / coverage and variation recheck
                             -> masked consensus only if final checks pass
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
| Minimum usable support for dominant construct | 50 reads |
| Minimum read / base quality | Q10 / Q15 |
| Competitive alignment identity / query coverage | 90% / 85%; no full-vector coverage requirement |
| Dominant-read target coverage for draft eligibility | 85% of the batch-inferred interval |
| Mapping quality; alignment-score separation | 20; 20 |
| Strong minor population | At least 5 reads and 10%, with at least 2 per strand |
| Weaker variation alert | At least 3 reads and 3% |
| Minimum consensus depth / base support | 20 / 80% |
| Required callable amplicon fraction | 98% |

The minority fraction is calculated among relevant assigned or informative reads, not raw POD5 counts. The same-construct linked test additionally requires two separated informative sites. Homopolymer-associated evidence is routed to review rather than treated as strong linked-haplotype evidence. Fixed differences shared by all reads are not themselves evidence of a mixture.

Only expected constructs with an inferred interval and `PASS_SCREEN` wells produce a read-supported draft and enter Dorado polish. Polishing uses an aligned, sorted/indexed BAM preserving Dorado metadata. Combining read groups is allowed only after verifying that all read-group headers specify the same basecalling model. A post-polish realignment rechecks support and residual variation, masks unsupported bases with `N`, and withholds consensus on failure. Final accepted status is `CONSENSUS_PASS`.

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
                                    decision, assignment, allele/haplotype evidence
<batch>/wells/<well>/consensus/consensus_result/
                                    post-polish decision; consensus only on pass
summary/                            well_summary.tsv, run_summary.json,
                                    accepted_consensus.fasta
```

Every expected well appears in the final summary, including `NO_READS`. Mixed and review wells are excluded from `accepted_consensus.fasta` without deleting their evidence. An all-rejected run legitimately produces an empty accepted FASTA. IDs include batch and well to prevent collisions. Variant positions in heterogeneity QC refer to the **original full reference**, restricted to the inferred target. Draft sequences span that target. Masked-position BEDs refer to the **polished amplicon**, whose coordinates can differ after indels. `inferred_targets.json` records original-reference boundaries.

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
