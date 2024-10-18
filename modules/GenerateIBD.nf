process GenerateIBD {
    tag "Generate genome-wide IBD estimates"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_6
    path prune_in_file

    output:
    path "ibd_prune.*"

    script:
    """
    plink --bfile target_6 --extract ${prune_in_file} --genome --out ibd_prune
    """
}
