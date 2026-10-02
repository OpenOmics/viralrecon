# POD5-first TCR analysis extension for nf-core/viralrecon

**Status:** implemented development prototype; not assay-validated or end-to-end runtime-tested. The supplied tests report **50 passed and 1 skipped**. Dorado, samtools, GPU basecalling, and Nextflow execution were unavailable in the build environment. Read `validation/VALIDATION.md` before running patient-derived data.

This package installs a separate, opt-in **POD5-only TCR workflow** into a clean checkout of `nf-core/viralrecon` pinned to:

```
fa23078485cb75e96add952045b2b897aab61b42
```

It supersedes the earlier FASTQ/Cutadapt integration. It is a local extension, not an official nf-core release or a submitted GitHub change. The installer retains the original viralrecon entry source alongside the new dispatcher. The stock branch has not been regression-tested after installation.

## 1. Implemented workflow

```
POD5 files grouped by supplied batch
    -> read-ID / metadata audit
    -> Dorado simplex basecalling, without trimming
    -> Dorado custom double-ended index classification and demultiplexing
       (both ends required; unclassified reads retained)
    -> per-well Dorado custom primer trimming and boundary validation
    -> per-well BAM + FASTQ
    -> competitive alignment to expected TCR amplicon references
    -> within-well heterogeneity screen
         REJECT / REVIEW -> no consensus; retain evidence and reads
         PASS_SCREEN     -> read-supported draft
                             -> Dorado polishing
                             -> realignment / coverage and variation recheck
                             -> masked consensus only if final checks pass
```

For PID201588, the expected construct panel should contain both **APC and PIK3CA**. Each batch/well is a separate analysis unit, for example `preinfusion_201588_run1__A1`. Acquisitions are not equated to biological wells. Batch boundaries remain explicit even if acquisition IDs appear in more than one directory.

This branch bypasses ARTIC, Guppyplex, viral lineage analysis, and the old Cutadapt demultiplexer. It does not claim that a consensus or tree proves a well began with one physical molecule. Its interpretation is **no detectable mixture above the configured thresholds**.

## 2. Installation and software

Use a fresh checkout, not the previously FASTQ-patched working tree:

```bash
git clone https://github.com/nf-core/viralrecon.git viralrecon-pod5
cd viralrecon-pod5
git checkout fa23078485cb75e96add952045b2b897aab61b42
cd ..

python3 viralrecon_tcr_pod5/apply.py viralrecon-pod5 --check-only
python3 viralrecon_tcr_pod5/apply.py viralrecon-pod5

git -C viralrecon-pod5 diff --stat
git -C viralrecon-pod5 status --short
```

The installer checks the commit, a clean working tree, the original `main.nf` blob, and expected schema structure. It refuses to overwrite an existing extension. It changes three root files and adds the workflow, processes, Python helpers, configuration, and documentation. It does not create commits or contact GitHub.

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
batch,pod5_dir,primer_table,sequencing_kit,basecall_model,reference_fasta,targets_bed
```

`examples/batches.csv` shows the two supplied batches. Edit its placeholders before use. Paths are resolved relative to the CSV's directory unless absolute. The current implementation supports local/shared-filesystem paths, not remote object-store URLs. There must be exactly one row per batch; `.pod5` files are discovered recursively under that row's directory.

* `sequencing_kit` is the actual ONT library kit, **not** the custom TCR barcode scheme name. The example uses the previously reported `SQK-LSK114`; the audit checks it against nonblank POD5 metadata.
* `basecall_model` is a pre-downloaded, versioned simplex model directory containing `config.toml`, not the moving alias `sup`. Select a model compatible with the acquisition and Dorado release. Do not assume the old FASTQ's v4.1.0 model is supported by a new Dorado release.
* `reference_fasta` contains the expected constructs for that batch. For this patient, combine the APC and PIK3CA reference records with distinct FASTA IDs.

### Full primer table

Supply the complete TSV for each batch, including:

```
Shorthand  Direction  Sequence  I-start  I-end  Pair_with
```

The parser expects `fA1` paired with `rA1`, and so on. It interprets `I-start` / `I-end` as **1-based inclusive** index coordinates; `3` and `9` extract seven bases. It validates shared flanking sequences, equal index lengths, unique indexes, and complete pairing. The default expects 48 pairs. The partial reverse-primer list pasted in the conversation is insufficient to run the whole plate.

From the table, the pipeline generates the custom Dorado TOML, barcode FASTA, well mapping, and a full-primer FASTA for each well. Rear-primer sequences are supplied in their original oligo orientation; Dorado handles reverse-complement recognition. `Frag-start` / `Frag-end` are not interpreted: their previous software-specific semantics are not established. This workflow intentionally trims the entire indexed PCR primer rather than reproducing an undocumented index-only trimming rule.

### Expected amplicon intervals

`targets_bed` must have exactly four **tab-separated** columns:

```
reference_FASTA_ID    start_0based_inclusive    end_0based_exclusive    construct_ID
```

Provide one interval per construct, with construct IDs such as `APC` and `PIK3CA`. The interval must be the **expected primer-trimmed amplicon interior**, including any vector sequence genuinely covered by the assay, not the entire supplied vector reference and not merely the TCR coding sequence. This prevents expected unsequenced vector regions from incorrectly failing coverage QC.

The pipeline crops a local reference panel from these intervals. Verify the primer binding positions and amplicon boundaries against the actual reference records. `examples/targets.template.bed` contains format guidance only; no biological coordinates have been invented.

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

Replace the paths and site configuration. `examples/slurm.config` is a template with explicit partition placeholders; it is not a verified NIH cluster configuration. The `tcr_gpu` label needs a real GPU resource request. Without a scheduler configuration, Nextflow runs locally.

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

Competitive Dorado/minimap2 alignment retains secondary hits so shared vector or constant-region sequence is not automatically treated as an unambiguous construct assignment. Ambiguous, partial, low-quality, and split alignments are accounted for separately.

The custom Python screen examines both:

1. **Between-construct mixtures:** substantial, strand-supported APC and PIK3CA read populations in the same well.
2. **Within-construct variation:** candidate SNPs and primitive indels, including linked alleles on the same reads. Two coherent, strand-supported patterns at separated sites trigger rejection. Isolated or weaker variation triggers review rather than a claim of two templates.

Examples of the **exploratory defaults**, not established sensitivity/specificity:

| Check | Default |
|---|---|
| Minimum usable support for dominant construct | 50 reads |
| Minimum read / base quality | Q10 / Q15 |
| Alignment identity; query and target coverage | 90%; 85% and 85% |
| Mapping quality; alignment-score separation | 20; 20 |
| Strong minor population | At least 5 reads and 10%, with at least 2 per strand |
| Weaker variation alert | At least 3 reads and 3% |
| Minimum consensus depth / base support | 20 / 80% |
| Required callable amplicon fraction | 98% |

The minority fraction is calculated among relevant assigned or informative reads, not raw POD5 counts. The same-construct linked test additionally requires two separated informative sites. Homopolymer-associated evidence is routed to review rather than treated as strong linked-haplotype evidence. Fixed differences shared by all reads are not themselves evidence of a mixture.

Only `PASS_SCREEN` wells produce a read-supported draft and enter Dorado polish. Polishing uses an aligned, sorted/indexed BAM preserving Dorado metadata. Combining read groups is allowed only after verifying that all read-group headers specify the same basecalling model. A post-polish realignment rechecks support and residual variation, masks unsupported bases with `N`, and withholds consensus on failure. Final accepted status is `CONSENSUS_PASS`.

**Limits:** identical templates are indistinguishable; low-abundance, closely related, or poorly sequenced mixtures can be missed; systematic ONT/PCR errors can cause review or false mixture signals. Complex rearrangements and repeat-associated indels require manual examination. This is a conservative candidate screen, not a validated haplotype caller, a phylogenetic proof of single-template origin, or a clinical pipeline. Do not use it to estimate the number of original molecules from sequencing read counts.

## 7. Outputs

Outputs are published beneath `outdir/tcr_pod5/`:

```
audit/                              POD5 inventory and raw UUID provenance
<batch>/scheme/                     generated scheme and cropped references
<batch>/basecalling/                calls.bam, summary, model/version provenance
<batch>/demultiplexing/demux/        assignment TSV, per-well BAMs, unclassified reads
<batch>/wells/<well>/trimming/       clean.bam, <well>.fastq.gz, failed reads, audit
<batch>/wells/<well>/alignment/      coordinate-sorted BAM and index
<batch>/wells/<well>/heterogeneity/qc/
                                    decision, assignment, allele/haplotype evidence
<batch>/wells/<well>/consensus/consensus_result/
                                    post-polish decision; consensus only on pass
summary/                            well_summary.tsv, run_summary.json,
                                    accepted_consensus.fasta
```

Every expected well appears in the final summary, including `NO_READS`. Mixed and review wells are excluded from `accepted_consensus.fasta` without deleting their evidence. An all-rejected run legitimately produces an empty accepted FASTA. IDs include batch and well to prevent collisions. Variant positions in QC refer to the **cropped amplicon**; masked-position BEDs refer to the **polished amplicon**. `scheme.json` records original-reference coordinates.

This branch produces TSV/JSON reports, not a customized MultiQC report or a final multiple-sequence alignment/tree. Those are separate extensions beyond the requested consensus gate.

## 8. Validation and interpretation

Run included tests from the unpacked package:

```bash
python3 -m unittest discover -s viralrecon_tcr_pod5/tests -v
```

With Dorado 2.1.2 and samtools available, enable the optional executable smoke test:

```bash
TCR_TEST_DORADO=/opt/dorado-2.1.2-linux-x64/bin/dorado \
TCR_TEST_SAMTOOLS=/path/to/samtools \
python3 -m unittest discover -s viralrecon_tcr_pod5/tests -v
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
- [POD5 tools and Python API](https://pod5-file-format.readthedocs.io/en/latest/)

The mixture rules and thresholds are custom prototype choices, not recommendations established by those tool manuals.
