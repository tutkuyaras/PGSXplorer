process RemoveHetOutliers {
    tag "Remove heterozygosity rate outliers"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_7
    path het_fail_ind_file

    output:
    path "target_8.*"

    script:
    """
    plink --bfile target_7 --allow-no-sex --remove ${het_fail_ind_file} --make-bed --out target_8
    """
}
