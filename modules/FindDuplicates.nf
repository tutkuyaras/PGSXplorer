process FindDuplicates {
    tag "Find duplicate SNPs "
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path merged

    output:
    path "last_duplicate_snps.list"

    script:
    """
    awk '{print \$2}' merged.bim | sort | uniq -d > last_duplicate_snps.list
    """
}
