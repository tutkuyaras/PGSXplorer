process InitialHWEFilter {
    tag "Initial HWE filtering"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_4

    output:
    path "target_5*"

    script:
    """
    plink --bfile target_4 --hwe ${params.hwe_case} --make-bed --out target_5
    """
}
