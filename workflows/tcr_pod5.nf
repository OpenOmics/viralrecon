import groovy.json.JsonSlurper
import java.nio.file.Files
import java.util.stream.Collectors

include {
    TCR_POD5_AUDIT; TCR_PREPARE; TCR_DORADO_BASECALL; TCR_DORADO_DEMUX;
    TCR_DORADO_TRIM; TCR_ALIGN; TCR_WELL_QC; TCR_DORADO_POLISH; TCR_REPORT
} from '../modules/local/tcr_pod5/main'

def localPath(base, value) {
    if (!value || value.toString().contains('\n') || value.toString().contains('\r'))
        error 'Missing path or newline in TCR samplesheet path'
    if (value.toString().contains('://'))
        error 'This TCR extension currently requires local/shared-filesystem paths, not remote URLs'
    def p = java.nio.file.Paths.get(value.toString())
    file(p.isAbsolute() ? p : base.resolve(p).normalize(), checkIfExists: true)
}

workflow TCR_POD5 {
    main:
    if (!params.tcr_input) error 'POD5 mode requires --tcr_input batches.csv'
    if (!params.outdir) error 'POD5 mode requires --outdir'
    if (params.input || params.fastq_dir) error 'TCR POD5 mode does not accept --input or --fastq_dir'
    if (params.platform && params.platform != 'nanopore') error 'TCR POD5 mode requires platform nanopore'
    if (!params.tcr_polish_models) error 'Provide --tcr_polish_models with a pre-populated Dorado polishing model cache'
    if (!(params.tcr_expected_pairs instanceof Number) || params.tcr_expected_pairs < 1)
        error '--tcr_expected_pairs must be a positive integer'

    def input = file(params.tcr_input, checkIfExists: true)
    def root = input.toAbsolutePath().parent
    def tooling = Channel.value([
        file("${projectDir}/bin/tcr_pod5.py", checkIfExists: true),
        file("${projectDir}/bin/tcr_pod5_core.py", checkIfExists: true)
    ])
    def settings = Channel.value(file(params.tcr_qc_config ?: "${projectDir}/assets/tcr_pod5/qc_defaults.json", checkIfExists: true))
    def models = Channel.value(file(params.tcr_polish_models, checkIfExists: true))

    ch_batches = Channel.fromPath(input)
        .splitCsv(header: true, strip: true)
        .map { row ->
            for (key in ['batch','pod5_dir','primer_table','sequencing_kit','basecall_model','reference_fasta','targets_bed'])
                if (!row[key]) error "Missing ${key} in TCR samplesheet"
            if (!(row.batch ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) error "Unsafe batch ID: ${row.batch}"
            def pod5dir = localPath(root, row.pod5_dir)
            if (!Files.isDirectory(pod5dir)) error "Not a POD5 directory: ${pod5dir}"
            def stream = Files.walk(pod5dir)
            def pod5s
            try {
                pod5s = stream.filter { Files.isRegularFile(it) && it.fileName.toString().endsWith('.pod5') }
                    .sorted().collect(Collectors.toList())
            } finally { stream.close() }
            if (!pod5s) error "No .pod5 files found under ${pod5dir}"
            def meta = [id: row.batch, kit: row.sequencing_kit]
            tuple(meta, pod5s, localPath(root, row.primer_table), localPath(root, row.reference_fasta),
                  localPath(root, row.targets_bed), localPath(root, row.basecall_model))
        }
        .collect(flat: false)
        .map { batches ->
            if (!batches) error 'TCR samplesheet is empty'
            def ids = batches.collect { it[0].id }
            if (ids.unique(false).size() != ids.size()) error 'Each batch must have exactly one samplesheet row'
            batches
        }
        .flatMap { it }

    ch_audit = ch_batches.collect(flat: false).map { batches ->
        def entries = []
        def paths = []
        batches.each { meta, pod5s, primers, refs, bed, model ->
            pod5s.each { p ->
                entries << [batch: meta.id, source: p.toString(), kit: meta.kit]
                paths << p
            }
        }
        tuple(entries, paths)
    }
    TCR_POD5_AUDIT(ch_audit, tooling)
    TCR_PREPARE(ch_batches.map { meta, p, primers, refs, bed, model -> tuple(meta, primers, refs, bed) }, tooling)
    TCR_DORADO_BASECALL(
        ch_batches.map { meta, p, primers, refs, bed, model -> tuple(meta, p, model) },
        TCR_POD5_AUDIT.out.report.first(), tooling
    )
    TCR_DORADO_DEMUX(TCR_PREPARE.out.scheme.join(TCR_DORADO_BASECALL.out.bam), tooling)

    ch_wells = TCR_DORADO_DEMUX.out.result.flatMap { meta, demux, scheme ->
        def rows = new JsonSlurper().parseText(demux.resolve('manifest.json').text)
        rows.findAll { it.well != 'unclassified' && it.reads > 0 }.collect { row ->
            def m = [id: "${meta.id}__${row.well}".toString(), batch: meta.id, well: row.well, kit: meta.kit]
            tuple(m, demux.resolve("wells/${row.well}.bam"), scheme.resolve("primers/${row.well}.fasta"),
                  scheme.resolve('amplicons.fasta'))
        }
    }
    TCR_DORADO_TRIM(ch_wells, settings, tooling)
    TCR_ALIGN(TCR_DORADO_TRIM.out.reads, tooling)
    TCR_WELL_QC(TCR_ALIGN.out.for_qc, settings, tooling)
    TCR_DORADO_POLISH(TCR_WELL_QC.out.passed, models, settings, tooling)
    TCR_REPORT(
        TCR_DORADO_DEMUX.out.result.map { meta, demux, scheme -> demux }.collect().ifEmpty([]),
        TCR_WELL_QC.out.report.map { meta, result -> result }.collect().ifEmpty([]),
        TCR_DORADO_POLISH.out.report.map { meta, result -> result }.collect().ifEmpty([]),
        tooling
    )

    emit:
    well_summary = TCR_REPORT.out.table
    consensus = TCR_REPORT.out.consensus
    summary = TCR_REPORT.out.json
}
