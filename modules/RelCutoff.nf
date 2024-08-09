process RelCutoff {
    tag "Create file with non-relative individuals"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_6
    path prune_in_file

    output:
    path "pihat.rel.id"

    script:
    """
    plink --bfile target_6 --extract ${prune_in_file} --rel-cutoff ${params.pihat} --out pihat
    """
}
