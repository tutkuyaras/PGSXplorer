process CalculateHWE {
    tag "Calculating Hardy-Weinberg Equilibrium"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_4

    output:
    path "HWE.hwe"

    script:
    """
    plink --bfile target_4 --hardy --out HWE
    """
}
