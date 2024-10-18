process CalculatePGS {
    tag "Calculate PGS by PLINK"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_9
    path gwas_sumstat
    path range_list
    path valid_snp
    path snp_pvalue

    output:
    path "target_9_OR.*"

    script:
    """
    plink \
        --bfile target_9 \
        --score ${gwas_sumstat} 3 5 8 header \
        --q-score-range ${range_list} ${snp_pvalue} \
        --extract ${valid_snp} \
        --out target_9_OR
    """
}
