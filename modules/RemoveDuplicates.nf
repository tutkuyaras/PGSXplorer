process RemoveDuplicates {
    tag "Remove duplicates"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path last
    path duplicate_snps_list

    output:
    path "last.*"

    script:
    """
    plink --bfile merged --exclude ${duplicate_snps_list} --make-bed --out last
    """
}
