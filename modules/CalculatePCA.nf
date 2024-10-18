process CalculatePCA {
    tag "CalculatePCA"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_9
    path prune_in_file

    output:
    tuple path("*.eigenvec"), path("*.eigenval")


    script:
    """
    plink --bfile target_9 \
          --extract ${prune_in_file} \
          --pca ${params.pca} \
          --out target_9_OR
    """
}
