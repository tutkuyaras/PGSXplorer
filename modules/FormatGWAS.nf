process FormatGWAS {
    tag "Format GWAS for SBayesR-C"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gwas

    output:
    path "formatted_gwas.ma"

    script:
    """
    Rscript ${params.formatgwas} ${gwas} 
    """
}


