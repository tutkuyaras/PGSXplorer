process SecondHWEFilter {
    tag "Second HWE filtering"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_5

    output:
    path "target_6.*"

    script:
    """
    plink --bfile target_5 --hwe ${params.hwe_ctrl} --hwe-all --make-bed --out target_6
    """
}
