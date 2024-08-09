process QCSNPmissingness {
    tag "Removal of missing SNPs according to the specified threshold value"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target

    output:
    path "target_*"

    script:
    """
    plink --bfile $target"/target" --geno ${params.geno} --make-bed --out target_1
    """
}
