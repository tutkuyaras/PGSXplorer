process CreateSNPpvalue {
    tag "CreateSNPpvalue"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path gwas_sumstat

    output:
    path "SNP.pvalue"

    script:
    """
    awk '{print \$3,\$13}' ${gwas_sumstat} > SNP.pvalue
    """
}
