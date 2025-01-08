process IBDLDPruning {
    tag "IBD-LD-pruning"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_6

    output:
    path "merged_indepSNP.prune.in", emit: prune_in
    path "merged_indepSNP.prune.out"

    script:
    """
    plink --bfile target_6 --indep-pairwise ${params.indep_window_size} ${params.indep_step_size} ${params.indep_threshold} --out merged_indepSNP
    """
}
