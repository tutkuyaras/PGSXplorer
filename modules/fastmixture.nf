process fastmixture {
    tag "Target Ancestry Inference"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path last

    output:
    path "*"

    script:
    """
    fastmixture --bfile last --K ${params.num_of_ancestry} --threads 8 --seed 1 --out ancestry
    """
}
