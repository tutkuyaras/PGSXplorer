process FilterMAF {
    tag "Filtering autosomal SNPs and calculating MAF"
    debug true
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path target_2

    output:
    tuple path("target_3.*"), path("MAF_check.frq")

    script:
    """
    awk '{ if (\$1 >= 1 && \$1 <= 22) print \$2 }' target_2.bim > snp_1_22.txt
    plink --bfile target_2 --extract snp_1_22.txt --make-bed --out target_3
    plink --bfile target_3 --freq --out MAF_check
    """
}
