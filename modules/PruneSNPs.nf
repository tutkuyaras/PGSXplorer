process PruneSNPs {
    tag "PruneSNPs by PLINK"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_9

    output:
    path "target_9_OR.prune.in"
    
    script:
    """
    plink \
        --bfile target_9  \
        --indep-pairwise ${params.indep_window_size} ${params.indep_step_size} ${params.indep_threshold} \
        --out target_9_OR
    """
}
