process LDpred2auto {
    tag "Calculate PGS with LDpred2 Auto Model"
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
    Rscript ${params.ldpred_auto_script} ${pheno_file} ${eigenvec} ${gwas_sumstat} ${target_9} 
    """
}

