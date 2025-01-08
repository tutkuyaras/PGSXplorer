process CalculatePCA {
    tag "CalculatePCA"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path last
    path prune_in_file

    output:
    tuple path("*.eigenvec"), path("*.eigenval")


    script:
    """
    plink --bfile last \
          --extract ${prune_in_file} \
          --pca ${params.pca} \
          --out last_OR
    """
}
