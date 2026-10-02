// POD5-only, metadata-preserving TCR path. All dependencies travel through channels.
import groovy.json.JsonOutput

def tq(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }
def toolsArgs() {
    "--dorado ${tq(params.tcr_dorado_bin)} --samtools ${tq(params.tcr_samtools_bin)} --version ${tq(params.tcr_dorado_version)}"
}

process TCR_POD5_AUDIT {
    label 'tcr_cpu'
    publishDir "${params.outdir}/tcr_pod5/audit", mode: 'copy', pattern: 'pod5_*'
    input:
    tuple val(entries), path(pod5s, stageAs: 'pod5_input??/*')
    path tooling
    output:
    path 'pod5_audit.json', emit: report
    path 'audit.ok', emit: gate
    path 'pod5_read_ids.tsv', emit: read_ids
    script:
    def paths = pod5s instanceof List ? pod5s : [pod5s]
    def staged = entries.withIndex().collect { e, i -> e + [staged: paths[i].toString()] }
    """
    cat > inputs.json <<'TCR_PAYLOAD'
    ${JsonOutput.toJson(staged)}
    TCR_PAYLOAD
    python3 tcr_pod5.py audit --manifest inputs.json ${params.tcr_allow_empty_pod5 ? '--allow-empty' : ''}
    """
}

process TCR_PREPARE {
    tag "${meta.id}"
    label 'tcr_cpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.id}/scheme" }, mode: 'copy', pattern: 'scheme'
    input:
    tuple val(meta), path(primers), path(refs)
    path tooling
    output:
    tuple val(meta), path('scheme'), emit: scheme
    script:
    """
    python3 tcr_pod5.py prepare \\
      --primers ${tq(primers)} --references ${tq(refs)} \\
      --expected-constructs ${tq(meta.expected_constructs.join(';'))} \\
      --expected-pairs ${params.tcr_expected_pairs} --barcode-errors ${params.tcr_barcode_errors} \\
      --outdir scheme
    """
}

process TCR_DORADO_BASECALL {
    tag "${meta.id}"
    label 'tcr_gpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.id}/basecalling" }, mode: 'copy',
        saveAs: { name -> name == 'calls.bam' && !params.tcr_save_basecalls ? null : name }
    input:
    tuple val(meta), path(pod5s, stageAs: 'pod5_input??/*'), path(model)
    path audit
    path tooling
    output:
    tuple val(meta), path('calls.bam'), emit: bam
    path 'basecall_provenance.json', emit: provenance
    path 'basecall_summary.tsv', emit: summary
    path 'dorado*.log', emit: logs
    script:
    def paths = pod5s instanceof List ? pod5s : [pod5s]
    """
    cat > pod5_files.json <<'TCR_PAYLOAD'
    ${JsonOutput.toJson(paths.collect { it.toString() })}
    TCR_PAYLOAD
    python3 tcr_pod5.py basecall --files pod5_files.json \\
      --model ${tq(model)} --batch ${tq(meta.id)} --audit ${tq(audit)} \\
      --device ${tq(params.tcr_basecall_device)} ${toolsArgs()}
    """
}

process TCR_DORADO_DEMUX {
    tag "${meta.id}"
    label 'tcr_cpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.id}/demultiplexing" }, mode: 'copy', pattern: 'demux'
    input:
    tuple val(meta), path(scheme), path(bam)
    path tooling
    output:
    tuple val(meta), path('demux'), path(scheme), emit: result
    script:
    """
    python3 tcr_pod5.py demux --bam ${tq(bam)} --scheme ${tq(scheme)} \\
      --batch ${tq(meta.id)} --threads ${task.cpus} ${toolsArgs()}
    """
}

process TCR_DORADO_TRIM {
    tag "${meta.id}"
    label 'tcr_cpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.batch}/wells/${meta.well}/trimming" }, mode: 'copy',
        pattern: '*.{bam,gz,json,tsv,log}'
    input:
    tuple val(meta), path(bam), path(primers), path(refs)
    path settings
    path tooling
    output:
    tuple val(meta), path('clean.bam'), path('trim_stats.json'), path(refs), emit: reads
    path '*.fastq.gz', emit: fastq
    path 'trim_audit.tsv', emit: audit
    path 'trim_failed.bam', emit: failed
    path 'dorado_trim.log', emit: log
    script:
    """
    python3 tcr_pod5.py trim --bam ${tq(bam)} --primers ${tq(primers)} \\
      --well ${tq(meta.well)} --kit ${tq(meta.kit)} --settings ${tq(settings)} \\
      --threads ${task.cpus} ${toolsArgs()}
    """
}

process TCR_ALIGN {
    tag "${meta.id}"
    label 'tcr_cpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.batch}/wells/${meta.well}/alignment" }, mode: 'copy',
        pattern: 'aligned*'
    input:
    tuple val(meta), path(clean), path(trim_stats), path(refs)
    path tooling
    output:
    tuple val(meta), path('name_sorted.bam'), path(refs), path(trim_stats), path(clean), emit: for_qc
    tuple val(meta), path('aligned.bam'), path('aligned.bam.bai'), emit: bam
    path '*.log', emit: log, optional: true
    script:
    """
    python3 tcr_pod5.py align --bam ${tq(clean)} --references ${tq(refs)} \\
      --threads ${task.cpus} ${toolsArgs()}
    """
}


process TCR_INFER_TARGETS {
    tag "${meta.id}"
    label 'tcr_cpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.id}/target_inference" }, mode: 'copy', pattern: 'targets'
    input:
    tuple val(meta), path(scheme), val(wells), path(bams, stageAs: 'well_bams??/*')
    path settings
    path tooling
    output:
    tuple val(meta), path('targets'), emit: result
    script:
    def paths = bams instanceof List ? bams : (bams ? [bams] : [])
    def entries = wells.withIndex().collect { well, i -> [well: well, bam: paths[i].toString()] }
    def payload = [batch: meta.id, expected_constructs: (meta.expected_constructs ?: null), alignments: entries]
    """
    cat > inference_inputs.json <<'TCR_PAYLOAD'
    ${JsonOutput.toJson(payload)}
    TCR_PAYLOAD
    python3 tcr_pod5.py infer-targets --manifest inference_inputs.json \\
      --batch ${tq(meta.id)} --references ${tq(scheme.resolve('reference.fasta'))} \\
      --settings ${tq(settings)} --samtools ${tq(params.tcr_samtools_bin)} --outdir targets
    """
}

process TCR_WELL_QC {
    tag "${meta.id}"
    label 'tcr_cpu'
    publishDir { "${params.outdir}/tcr_pod5/${meta.batch}/wells/${meta.well}/heterogeneity" }, mode: 'copy', pattern: 'qc'
    input:
    tuple val(meta), path(bam), path(refs), path(trim_stats), path(clean), path(targets)
    path settings
    path tooling
    output:
    tuple val(meta), path('qc'), emit: report
    tuple val(meta), path('pass'), path(clean), emit: passed, optional: true
    script:
    """
    python3 tcr_pod5.py qc --bam ${tq(bam)} --references ${tq(refs)} \\
      --settings ${tq(settings)} --trim-stats ${tq(trim_stats)} \\
      --inferred-targets ${tq(targets.resolve('inferred_targets.json'))} \\
      --sample ${tq(meta.id)} --batch ${tq(meta.batch)} --well ${tq(meta.well)} \\
      --samtools ${tq(params.tcr_samtools_bin)} --threads ${task.cpus}
    """
}

process TCR_DORADO_POLISH {
    tag "${meta.id}"
    label 'tcr_polish'
    publishDir { "${params.outdir}/tcr_pod5/${meta.batch}/wells/${meta.well}/consensus" }, mode: 'copy',
        pattern: 'consensus_result'
    input:
    tuple val(meta), path(passed), path(clean)
    path models
    path settings
    path tooling
    output:
    tuple val(meta), path('consensus_result'), emit: report
    script:
    """
    python3 tcr_pod5.py polish --bam ${tq(clean)} --passed ${tq(passed)} \\
      --models ${tq(models)} --settings ${tq(settings)} --sample ${tq(meta.id)} \\
      --device ${tq(params.tcr_polish_device)} --threads ${task.cpus} ${toolsArgs()}
    """
}

process TCR_REPORT {
    label 'tcr_cpu'
    publishDir "${params.outdir}/tcr_pod5/summary", mode: 'copy'
    input:
    path demux_dirs, stageAs: 'demux_input??/*'
    path qc_dirs, stageAs: 'qc_input??/*'
    path polish_dirs, stageAs: 'polish_input??/*'
    path target_dirs, stageAs: 'targets_input??/*'
    path tooling
    output:
    path 'well_summary.tsv', emit: table
    path 'accepted_consensus.fasta', emit: consensus
    path 'run_summary.json', emit: json
    script:
    def asStrings = { v -> (v instanceof List ? v : (v ? [v] : [])).collect { it.toString() } }
    def payload = [demux_dirs: asStrings(demux_dirs), qc_dirs: asStrings(qc_dirs), polish_dirs: asStrings(polish_dirs), target_dirs: asStrings(target_dirs)]
    """
    cat > report_inputs.json <<'TCR_PAYLOAD'
    ${JsonOutput.toJson(payload)}
    TCR_PAYLOAD
    python3 tcr_pod5.py report --manifest report_inputs.json
    """
}
