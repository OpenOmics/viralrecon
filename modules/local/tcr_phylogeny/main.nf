// Optional downstream analysis: never feeds back into QC or consensus decisions.
def phyloQuote(value) { "'" + value.toString().replace("'", "'\"'\"'") + "'" }

process TCR_PHYLOGENY {
    tag "${meta.id}"
    label 'tcr_phylogeny'
    conda "${moduleDir}/environment.yml"
    publishDir { "${params.outdir}/tcr_pod5/${meta.id}" }, mode: 'copy', pattern: 'phylogeny'

    input:
    tuple val(meta), path(refs, stageAs: 'reference.fasta'), path(targets, stageAs: 'inferred_targets.json'),
          path(summary, stageAs: 'well_summary.tsv'), path(accepted, stageAs: 'accepted_consensus.fasta'),
          path(review, stageAs: 'diagnostic_consensus.fasta')
    path tooling

    output:
    tuple val(meta), path('phylogeny'), emit: result

    script:
    """
    python3 tcr_pod5_phylogeny.py \\
      --batch ${phyloQuote(meta.id)} --references reference.fasta --targets inferred_targets.json \\
      --summary well_summary.tsv --accepted accepted_consensus.fasta --review diagnostic_consensus.fasta \\
      --outdir phylogeny --threads ${task.cpus} \\
      --mafft ${phyloQuote(params.tcr_mafft_bin)} --raxml-ng ${phyloQuote(params.tcr_raxml_ng_bin)} \\
      --model ${phyloQuote(params.tcr_phylogeny_model)} --bootstraps ${phyloQuote(params.tcr_phylogeny_bootstraps)} \\
      --seed ${phyloQuote(params.tcr_phylogeny_seed)} \\
      --min-target-overlap ${phyloQuote(params.tcr_phylogeny_min_target_overlap)}
    """
}
