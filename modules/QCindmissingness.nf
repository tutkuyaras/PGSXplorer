process QCindmissingness {
    tag "Removal of missing individuals according to the specified threshold value"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_1

    output:
    path 'target_2.*'

    script:
    """
    plink --bfile target_1 --mind ${params.mind} --make-bed --out target_2
    """
}
