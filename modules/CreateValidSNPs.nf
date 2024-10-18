process CreateValidSNPs {
    tag "CreateValidSNPs"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path clumped_file

    output:
    path "target_9.valid.snp"

    script:
    """
    awk 'NR!=1{print \$3}' ${clumped_file} > target_9.valid.snp
    """
}
