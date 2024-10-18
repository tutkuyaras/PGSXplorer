process LDpred2grid {
    tag "Calculate PGS with LDpred2 Grid Model"
    debug true
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path pheno_file
    path eigenvec
    path gwas_sumstat
    path target_9

    output:
    path "*"
    
    script:
    """
    Rscript ${params.ldpred_grid_script} ${pheno_file} ${eigenvec} ${gwas_sumstat} ${target_9}
    """
}

