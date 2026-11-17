import groovy.json.JsonOutput

process TCR_MULTIQC {
    label 'tcr_report'
    conda "${moduleDir}/environment.yml"
    publishDir "${params.outdir}/tcr_pod5", mode: 'copy', pattern: 'multiqc'

    input:
    tuple val(destinations), path(report_files, stageAs: 'report_inputs????/*')
    path report_tools, stageAs: 'tcr_multiqc_tools'

    output:
    path 'multiqc', emit: report

    script:
    def paths = report_files instanceof List ? report_files : [report_files]
    if (paths.size() != destinations.size()) error 'TCR MultiQC input metadata/file count mismatch'
    def entries = destinations.withIndex().collect { destination, i ->
        [source: paths[i].toString(), destination: destination]
    }
    """
    export MPLCONFIGDIR="\$PWD/.matplotlib"
    export XDG_CACHE_HOME="\$PWD/.cache"
    cat > report_inputs.json <<'TCR_MULTIQC_INPUTS'
    ${JsonOutput.toJson(entries)}
    TCR_MULTIQC_INPUTS
    python3 tcr_multiqc_tools/stage_tcr_multiqc.py \\
      --manifest report_inputs.json --outdir report_input
    python3 tcr_multiqc_tools/collect_tcr_multiqc.py \\
      --results-dir report_input --outdir multiqc --no-discover-audits --render
    cp report_inputs.json multiqc/input_manifest.json
    python3 - <<'TCR_MULTIQC_VERSIONS'
    import importlib.metadata
    import json
    import platform
    from pathlib import Path
    versions = {
        'python': platform.python_version(),
        'multiqc': importlib.metadata.version('multiqc'),
        'matplotlib': importlib.metadata.version('matplotlib'),
        'collector': json.loads(Path('multiqc/collector_manifest.json').read_text())['collector_version'],
    }
    Path('multiqc/report_software_versions.json').write_text(json.dumps(versions, indent=2) + '\\n')
    TCR_MULTIQC_VERSIONS
    """
}
