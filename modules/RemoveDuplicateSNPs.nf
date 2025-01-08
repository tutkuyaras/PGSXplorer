process RemoveDuplicateSNPs {
    tag "Remove duplicate SNPs"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_8
    path duplicate_snps_list

    output:
    path "target_9.*"

    script:
    """
    plink --bfile target_8 --exclude ${duplicate_snps_list} --make-bed --out target_9
    """
}
