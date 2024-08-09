process CleanHetFailInd {
    tag "Clean het fail individuals file"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path fail_het_qc_file

    output:
    path "het_fail_ind.txt"

    script:
    """
    sed 's/\"//g' ${fail_het_qc_file} | awk '{print \$1, \$2}' > het_fail_ind.txt
    """
}
