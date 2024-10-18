process FilterIndividuals {
    tag "Filter individuals based on relatedness"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_6
    path rel_id_file

    output:
    path "target_7.*"

    script:
    """
    plink --bfile target_6 --keep ${rel_id_file} --make-bed --out target_7
    """
}
