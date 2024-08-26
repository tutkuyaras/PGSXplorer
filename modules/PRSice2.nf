process PRSice2 {
    tag "CalculatePGS with PRSice2"
    debug true
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path target_9
    path gwas_sumstat
    path pheno_file
    path eigenvec

    output:
    path "target_9_PRSice2*"

    script:
    """
    Rscript ${params.prsice_script} \
        --prsice ${params.prsice_executable} \
        --base ${gwas_sumstat} \
        --target target_9 \
        --binary-target T \
        --pheno ${pheno_file} \
        --cov ${eigenvec} \
        --stat OR \
        --or \
        --out target_9_PRSice2
    """
}

