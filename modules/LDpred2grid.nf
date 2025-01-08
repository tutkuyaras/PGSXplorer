process LDpred2grid {
    tag "Calculate PGS with LDpred2 Grid Model"
    debug true
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path pheno_file
    path eigenvec
    path gwas
    path last

    output:
    path "*"
    
    script:
    """
    Rscript ${params.ldpred_grid_script} ${pheno_file} ${eigenvec} ${gwas} ${last}
    """
}

