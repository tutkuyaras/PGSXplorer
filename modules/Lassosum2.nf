process Lassosum2 {
    tag "Lassosum2"
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
    Rscript ${params.lassosum2_script} ${pheno_file} ${eigenvec} ${gwas_sumstat} ${target_9}
    """
}

