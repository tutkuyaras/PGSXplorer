process CalculateHet {
    tag "Calculate heterozygosity rate"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_7
    path prune_in_file

    output:
    path "het_R_check.het"

    script:
    """
    plink --bfile target_7 --extract ${prune_in_file} --het --out het_R_check
    """
}
