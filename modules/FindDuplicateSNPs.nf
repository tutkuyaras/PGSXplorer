process FindDuplicateSNPs {
    tag "Find duplicate SNPs"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_8

    output:
    path "target_8_duplicate_snps.list"

    script:
    """
    awk '{print \$2}' target_8.bim | sort | uniq -d > target_8_duplicate_snps.list
    """
}
