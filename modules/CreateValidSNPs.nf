process CreateValidSNPs {
    tag "CreateValidSNPs"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path clumped_file

    output:
    path "last.valid.snp"

    script:
    """
    awk 'NR!=1{print \$3}' ${clumped_file} > last.valid.snp
    """
}
