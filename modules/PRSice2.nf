process PRSice2 {
    tag "CalculatePGS with PRSice2"
    debug true
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
    path last
    path gwas
    path pheno_file
    path eigenvec

    output:
    path "last_PRSice2*"

    script:
    """
    Rscript ${params.prsice_script} \
        --prsice ${params.prsice_executable} \
        --base ${gwas} \
        --target last \
        --binary-target T \
        --pheno ${pheno_file} \
        --cov ${eigenvec} \
        --stat OR \
        --or \
        --out last_PRSice2
    """
}

