process RemoveLowMAF {
    tag "Removing SNPs with low MAF frequency"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_3

    output:
    path "target_4.*"

    script:
    """
    plink --bfile target_3 --maf ${params.maf} --make-bed --out target_4
    """
}
